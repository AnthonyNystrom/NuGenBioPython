"""Small helpers for validating request input and normalising error responses.

These are narrow: used by routes that need to harden user-supplied data or
standardise server-side error logging. They intentionally avoid grand
abstractions — each call site can keep its existing try/except shape.
"""
import logging
import os
import re
import socket
import threading
import time
from contextlib import contextmanager

from flask import jsonify, request


log = logging.getLogger(__name__)


def require_json(req=None):
    """Return the JSON body as a dict, or {} if Content-Type was not JSON
    and no body was provided. Prevents the `None.get(...)` crash pattern
    when a client posts with the wrong Content-Type."""
    r = req if req is not None else request
    data = r.get_json(silent=True)
    if data is None:
        return {}
    if not isinstance(data, dict):
        # Some clients may post a JSON array; callers generally expect dicts.
        return {}
    return data


def clamp_int(value, default, *, lo=None, hi=None):
    """Coerce value to int and clamp into [lo, hi]. Returns default on any
    conversion failure."""
    try:
        n = int(value)
    except (TypeError, ValueError):
        return default
    if lo is not None and n < lo:
        return lo
    if hi is not None and n > hi:
        return hi
    return n


def clamp_float(value, default, *, lo=None, hi=None):
    """Coerce value to float and clamp into [lo, hi]."""
    try:
        n = float(value)
    except (TypeError, ValueError):
        return default
    if lo is not None and n < lo:
        return lo
    if hi is not None and n > hi:
        return hi
    return n


# ---------------------------------------------------------------------------
# Outbound rate limiting
# ---------------------------------------------------------------------------
# NCBI publishes hard limits (3 requests/second without an API key, 10 with)
# and enforces them by blocking the offending IP. Nothing in this app limited
# how fast it called out, so a handful of users clicking through the Database
# or BLAST pages could exceed the ceiling and get the whole deployment blocked
# — with no signal beyond sudden failures.
#
# This gates *outbound* traffic per service, which is the direction that gets
# you blocked. It is deliberately not a per-client API limiter: that is a
# different concern (abuse protection) and wants shared state across workers.

class RateGate:
    """Enforce a minimum interval between outbound calls to one service.

    Blocking by design: the caller is a request handler that is about to make
    a network call anyway, and a short sleep is far cheaper than being
    blocked by the upstream service. Thread-safe, though the app is intended
    to run one request per worker (see remote_timeout's threading note).
    """

    def __init__(self, name, min_interval):
        self.name = name
        self.min_interval = float(min_interval)
        self._lock = threading.Lock()
        self._last_call = 0.0

    def wait(self):
        """Block until it is safe to call. Returns seconds actually slept."""
        with self._lock:
            now = time.monotonic()
            earliest = self._last_call + self.min_interval
            delay = earliest - now
            if delay > 0:
                time.sleep(delay)
                now = time.monotonic()
            else:
                delay = 0.0
            self._last_call = now
        if delay > 0:
            log.debug('throttled %s call by %.3fs', self.name, delay)
        return delay


def _ncbi_interval():
    """NCBI allows 10 req/s with an API key, 3 without."""
    return 0.11 if os.environ.get('ENTREZ_API_KEY', '').strip() else 0.34


# One gate per upstream service. NCBI's is computed per call so that setting
# ENTREZ_API_KEY takes effect without a restart.
NCBI_GATE = RateGate('ncbi', 0.34)
KEGG_GATE = RateGate('kegg', 0.34)
BLAST_GATE = RateGate('blast', 3.0)
UNIPROT_GATE = RateGate('uniprot', 0.34)

_GATES = {
    'ncbi': NCBI_GATE,
    'kegg': KEGG_GATE,
    'blast': BLAST_GATE,
    'uniprot': UNIPROT_GATE,
}


def throttle(service):
    """Block until it is safe to call `service`. Unknown services pass through."""
    gate = _GATES.get(service)
    if gate is None:
        return 0.0
    if service == 'ncbi':
        gate.min_interval = _ncbi_interval()
    return gate.wait()


# Serialises the global-socket-timeout window in remote_timeout. See the note
# in that function for why serialising is free in practice.
_SOCKET_TIMEOUT_LOCK = threading.RLock()


@contextmanager
def remote_timeout(seconds, service=None):
    """Temporarily set the default socket timeout for remote calls.

    Pass `service` ('ncbi', 'kegg', 'blast', 'uniprot') to also apply that
    service's outbound rate gate before the call — see RateGate above.

    BioPython's REST/Entrez/NCBIWWW helpers use urllib under the hood and
    don't accept a ``timeout=`` kwarg. We set the default socket timeout
    for the duration of the call so a slow/hung external service aborts
    instead of blocking a Flask worker indefinitely. The original timeout
    is always restored, even on exception.

    THREAD SAFETY: this mutates the process-global default socket timeout,
    so two threads doing it at once would stomp each other's window — one
    request could inherit another's timeout, or restore the wrong value. A
    module-level lock serialises the whole call instead.

    Serialising costs almost nothing here: every caller passes a `service`,
    and the outbound RateGate already spaces those calls (0.34s for NCBI and
    KEGG, 3s for BLAST submissions), so concurrent external calls to the same
    service were never going to overlap anyway. What the lock buys is that a
    threaded server no longer corrupts the timeout window.
    """
    with _SOCKET_TIMEOUT_LOCK:
        if service:
            throttle(service)
        old = socket.getdefaulttimeout()
        socket.setdefaulttimeout(seconds)
        try:
            yield
        finally:
            socket.setdefaulttimeout(old)


# Messages from these are written for the person using the tool — "Sequences
# must all be the same length", "Mismatch, 4 open vs 0 close parentheses" —
# and are the most useful thing we can say back. ValueError also covers
# UploadError, which subclasses it.
_USER_FACING_EXCEPTIONS = (ValueError,)

# Anything raised from inside BioPython is a parse/validation message aimed at
# a biologist, so it is worth surfacing even though it is not a ValueError
# (Bio.Phylo.NewickIO.NewickError subclasses Exception directly).
_USER_FACING_MODULE_PREFIX = 'Bio.'

_GENERIC_ERROR = ('The server could not complete that request. '
                  'The details have been logged.')

# Absolute paths and Windows drive paths, so an upload directory layout does
# not travel to the browser inside an otherwise-safe message.
_PATH_RE = re.compile(r'(?:[A-Za-z]:)?(?:/[\w.\-]+){2,}/?|(?:\\[\w.\-]+){2,}')


def _is_user_facing(exc):
    if isinstance(exc, _USER_FACING_EXCEPTIONS):
        return True
    module = type(exc).__module__ or ''
    return module.startswith(_USER_FACING_MODULE_PREFIX)


def _scrub(message):
    """Remove filesystem paths from a message we are about to return."""
    return _PATH_RE.sub('<path>', message)


def safe_error_message(exc):
    """The most useful thing we can say about `exc` without leaking internals.

    Parse and validation errors carry messages written for the user, so they
    are passed through (with any filesystem path scrubbed). Everything else —
    OSError with an upload path in it, KeyError naming an internal dict key,
    AttributeError exposing object structure — collapses to a generic line,
    because the detail belongs in the log, not the browser.
    """
    if _is_user_facing(exc):
        text = _scrub(str(exc)).strip()
        return text or _GENERIC_ERROR
    return _GENERIC_ERROR



# ---------------------------------------------------------------------------
# Input size limits
# ---------------------------------------------------------------------------
# Bio.Align.PairwiseAligner is O(n*n) in BOTH time and memory. Measured on a
# random DNA pair: 10 kb -> 1.0 s / 243 MB, 20 kb -> 4.0 s / 861 MB,
# 40 kb -> 16.1 s / 4.0 GB, scaling cleanly 4x per doubling. With
# MAX_CONTENT_LENGTH at 16 MB and no length check, a single request could
# carry ~8,000,000 bp — roughly 179 hours of CPU, and an out-of-memory kill
# long before that.
#
# Memory is the sharper edge: extrapolating the measured curve, an 8 GB host
# is exhausted by a ~57 kb sequence. That is not an attack payload, it is one
# assembly contig — an entirely ordinary thing to paste in. Since the app runs
# one request per process (see remote_timeout's threading note), the OOM kill
# takes the whole service down.
#
# Linear work (GC content, translation, reverse complement) is cheap by
# comparison — 5 Mb in 0.18 s / 80 MB — so it gets a much higher ceiling.

# Quadratic algorithms: pairwise alignment and anything built on it.
MAX_ALIGNMENT_LEN = int(os.environ.get('MAX_ALIGNMENT_LEN', 10_000))

# Linear scans: composition, translation, ORF finding, restriction search.
MAX_SEQUENCE_LEN = int(os.environ.get('MAX_SEQUENCE_LEN', 5_000_000))


class SequenceTooLong(ValueError):
    """Raised when a submitted sequence exceeds the limit for an operation.

    Subclasses ValueError so safe_error_message() passes the (deliberately
    actionable) message through to the user unchanged.
    """


def _fmt_bp(n):
    """Human-readable length that still carries the exact figure.

    Rounding alone is misleading at the boundary: 10,001 bp against a 10,000
    limit rendered as "is 10 kb; the limit is 10 kb", which reads as a bug.
    """
    if n >= 1000:
        return f'{n:,} bp'
    return f'{n} bp'


def check_sequence_length(sequence, limit=None, name='sequence', quadratic=False):
    """Validate a sequence's length, returning it unchanged when acceptable.

    Raises SequenceTooLong with a message that names the actual limit, so the
    user knows what to do rather than watching the request hang.
    """
    if sequence is None:
        return sequence
    if limit is None:
        limit = MAX_ALIGNMENT_LEN if quadratic else MAX_SEQUENCE_LEN
    length = len(sequence)
    if length > limit:
        detail = (' This operation scales with the square of the sequence '
                  'length, so the limit is deliberately conservative.'
                  if quadratic else '')
        raise SequenceTooLong(
            f'{name} is {_fmt_bp(length)}; the limit for this operation is '
            f'{_fmt_bp(limit)}.{detail}'
        )
    return sequence


def check_alignment_inputs(*sequences, name='sequence'):
    """Length-check every sequence going into a quadratic alignment."""
    for i, seq in enumerate(sequences, 1):
        label = name if len(sequences) == 1 else f'{name} {i}'
        check_sequence_length(seq, quadratic=True, name=label)
    return sequences

# Errors that mean "the request was not valid", not "the server broke".
# SequenceTooLong is defined above; UploadError is imported lazily to avoid a
# circular import between request_helpers and upload_helpers.
def _input_validation_errors():
    try:
        from utils.upload_helpers import UploadError
        return (SequenceTooLong, UploadError)
    except Exception:
        return (SequenceTooLong,)


_INPUT_VALIDATION_ERRORS = _input_validation_errors()


def error_response(exc, user_msg=None, *, context=None):
    """Log the exception server-side and return a safe JSON response.

    The response shape is unchanged: {'success': False, 'error': <message>}.
    When `user_msg` is not given, the message is chosen by
    safe_error_message() rather than being str(exc) verbatim — 197 route
    handlers used to hand the raw exception text to the browser, which leaked
    upload paths and internal structure.
    """
    label = context or 'route error'
    if isinstance(exc, _INPUT_VALIDATION_ERRORS):
        # Someone pasted a sequence that is too long, or uploaded nothing.
        # That is ordinary use, not a fault: record it at INFO with no
        # traceback so real failures stay visible in the log.
        log.info('%s: %s', label, exc)
    else:
        log.exception('%s: %s', label, exc)
    if user_msg is None:
        user_msg = safe_error_message(exc)
    return jsonify({'success': False, 'error': user_msg})
