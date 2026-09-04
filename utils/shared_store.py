"""Bounded key-value stores that can be shared across worker processes.

Four caches (motif, BLAST, HMM, service status) were plain module-level
OrderedDicts. That works for exactly one worker process. Under gunicorn with
`--workers 2` a user's follow-up request lands on a different worker, finds
nothing, and gets "not found" for a result that was created seconds earlier —
a failure that only appears in production and looks random.

This gives them one implementation with a pluggable backend:

* No REDIS_URL  -> in-process LRU, identical to the previous behaviour.
* REDIS_URL set -> Redis, so every worker sees the same data.

The Redis backend is optional; if the package or the server is missing the
store falls back to memory and says so once, rather than failing to start.
"""
import logging
import os
import pickle
import threading
import time
from collections import OrderedDict

log = logging.getLogger(__name__)

REDIS_URL = os.environ.get('REDIS_URL', '').strip()
KEY_PREFIX = os.environ.get('SHARED_STORE_PREFIX', 'nugenbio')

_redis_client = None
_redis_checked = False
_redis_lock = threading.Lock()


def _get_redis():
    """Return a live Redis client, or None. Connects at most once."""
    global _redis_client, _redis_checked
    if _redis_checked:
        return _redis_client
    with _redis_lock:
        if _redis_checked:
            return _redis_client
        _redis_checked = True
        if not REDIS_URL:
            return None
        try:
            import redis  # optional dependency
            client = redis.Redis.from_url(REDIS_URL, socket_timeout=2)
            client.ping()
            _redis_client = client
            log.info('shared store: using Redis at %s', REDIS_URL)
        except Exception as exc:
            log.warning(
                'shared store: REDIS_URL is set but Redis is unavailable (%s). '
                'Falling back to per-process memory — do NOT run more than one '
                'worker in this state, results will be lost between requests.',
                exc,
            )
            _redis_client = None
        return _redis_client


class SharedStore:
    """A bounded, optionally TTL'd store addressed by string keys.

    `serializer` is 'pickle' (default) so arbitrary Python objects — a
    BioPython Motif, a trained HMM — round-trip. Pickle is only ever read back
    from this application's own Redis; if that Redis is shared with untrusted
    writers, use serializer='json' and store plain dicts instead.
    """

    def __init__(self, name, max_entries=64, ttl=None, serializer='pickle'):
        self.name = name
        self.max_entries = int(max_entries)
        self.ttl = ttl
        self.serializer = serializer
        self._local = OrderedDict()
        self._lock = threading.Lock()

    # -- backend ---------------------------------------------------------
    @property
    def backend(self):
        return 'redis' if _get_redis() is not None else 'memory'

    def _key(self, key):
        return f'{KEY_PREFIX}:{self.name}:{key}'

    def _dumps(self, value):
        if self.serializer == 'json':
            import json
            return json.dumps(value).encode()
        return pickle.dumps(value)

    def _loads(self, blob):
        if self.serializer == 'json':
            import json
            return json.loads(blob)
        return pickle.loads(blob)

    # -- api -------------------------------------------------------------
    def set(self, key, value):
        client = _get_redis()
        if client is not None:
            try:
                blob = self._dumps(value)
                if self.ttl:
                    client.setex(self._key(key), int(self.ttl), blob)
                else:
                    client.set(self._key(key), blob)
                return
            except Exception:
                log.warning('shared store %s: Redis write failed, using memory',
                            self.name, exc_info=True)
        with self._lock:
            if key in self._local:
                self._local.move_to_end(key)
            self._local[key] = (time.monotonic(), value)
            while len(self._local) > self.max_entries:
                self._local.popitem(last=False)

    def get(self, key, default=None):
        client = _get_redis()
        if client is not None:
            try:
                blob = client.get(self._key(key))
                return self._loads(blob) if blob is not None else default
            except Exception:
                log.warning('shared store %s: Redis read failed, using memory',
                            self.name, exc_info=True)
        with self._lock:
            entry = self._local.get(key)
            if entry is None:
                return default
            stored_at, value = entry
            if self.ttl and (time.monotonic() - stored_at) >= self.ttl:
                del self._local[key]
                return default
            self._local.move_to_end(key)
            return value

    def delete(self, key):
        client = _get_redis()
        if client is not None:
            try:
                client.delete(self._key(key))
            except Exception:
                pass
        with self._lock:
            self._local.pop(key, None)

    # Mapping-style access so the existing `cache[key]` call sites read
    # naturally. NOTE: with the Redis backend `store[k]` returns a *fresh*
    # deserialised object, so mutating it in place does not persist — do
    # read / modify / store[k] = value instead. See update() below.
    def __getitem__(self, key):
        sentinel = object()
        value = self.get(key, sentinel)
        if value is sentinel:
            raise KeyError(key)
        return value

    def __setitem__(self, key, value):
        self.set(key, value)

    def update_entry(self, key, **fields):
        """Read / modify / write one entry's top-level fields atomically enough.

        Needed because an in-place mutation of a nested dict is invisible to
        the Redis backend: the value it handed back was a fresh object.
        """
        current = self.get(key)
        if current is None:
            return False
        if isinstance(current, dict):
            current.update(fields)
            self.set(key, current)
            return True
        return False

    def __contains__(self, key):
        sentinel = object()
        return self.get(key, sentinel) is not sentinel

    def __len__(self):
        # Only meaningful for the memory backend; Redis entries expire on
        # their own and counting them would mean a keyspace scan.
        with self._lock:
            return len(self._local)

    def clear(self):
        with self._lock:
            self._local.clear()


def reset_backend_for_tests():
    """Force the next call to re-evaluate REDIS_URL."""
    global _redis_client, _redis_checked
    with _redis_lock:
        _redis_client = None
        _redis_checked = False
