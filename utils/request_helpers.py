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

    THREAD-SAFETY CAVEAT: this mutates the process-global default socket
    timeout via ``socket.setdefaulttimeout``. Under a threaded server
    (Flask ``threaded=True`` or gunicorn ``--threads > 1``) concurrent
    requests can stomp each other's timeout window. It is safe under a
    single-threaded / one-request-per-worker model (gunicorn sync workers,
    the default), which is how this app is intended to run. If you move to
    threaded workers, replace this with a per-call timeout instead.
    """
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


def error_response(exc, user_msg=None, *, context=None):
    """Log the exception server-side and return a safe JSON response.

    The response shape is unchanged: {'success': False, 'error': <message>}.
    When `user_msg` is not given, the message is chosen by
    safe_error_message() rather than being str(exc) verbatim — 197 route
    handlers used to hand the raw exception text to the browser, which leaked
    upload paths and internal structure.
    """
    label = context or 'route error'
    log.exception('%s: %s', label, exc)
    if user_msg is None:
        user_msg = safe_error_message(exc)
    return jsonify({'success': False, 'error': user_msg})
