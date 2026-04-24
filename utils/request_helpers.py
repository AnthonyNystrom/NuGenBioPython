"""Small helpers for validating request input and normalising error responses.

These are narrow: used by routes that need to harden user-supplied data or
standardise server-side error logging. They intentionally avoid grand
abstractions — each call site can keep its existing try/except shape.
"""
import logging
import socket
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


@contextmanager
def remote_timeout(seconds):
    """Temporarily set the default socket timeout for remote calls.

    BioPython's REST/Entrez/NCBIWWW helpers use urllib under the hood and
    don't accept a ``timeout=`` kwarg. We set the default socket timeout
    for the duration of the call so a slow/hung external service aborts
    instead of blocking a Flask worker indefinitely. The original timeout
    is always restored, even on exception.
    """
    old = socket.getdefaulttimeout()
    socket.setdefaulttimeout(seconds)
    try:
        yield
    finally:
        socket.setdefaulttimeout(old)


def error_response(exc, user_msg=None, *, context=None):
    """Log the exception server-side and return a safe JSON response.

    For callers that previously returned `{'error': str(e)}`, passing
    `user_msg=None` preserves the same payload while also logging server-side.
    Pass an explicit user_msg to avoid leaking library internals.
    """
    label = context or 'route error'
    log.exception('%s: %s', label, exc)
    if user_msg is None:
        user_msg = str(exc)
    return jsonify({'success': False, 'error': user_msg})
