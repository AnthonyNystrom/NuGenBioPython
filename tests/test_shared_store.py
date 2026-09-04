"""Cross-worker state.

Four caches (motif, BLAST, HMM, service status) were module-level
OrderedDicts. That is correct for exactly one worker process. Under gunicorn
with --workers 2, a user's follow-up request lands on a different worker,
finds nothing, and is told the result does not exist — seconds after creating
it. The failure only appears in production and looks random.

These tests pin the store that replaced them, including the subtle part: with
a Redis backend a get() returns a *fresh* object, so mutating it in place is
silently discarded.
"""
import time

import pytest

from utils import shared_store
from utils.shared_store import SharedStore


@pytest.fixture(autouse=True)
def _reset_backend(monkeypatch):
    monkeypatch.setattr(shared_store, "REDIS_URL", "")
    shared_store.reset_backend_for_tests()
    yield
    shared_store.reset_backend_for_tests()


# --------------------------------------------------------------------------
# Memory backend — the default, and what runs today
# --------------------------------------------------------------------------

def test_defaults_to_memory_without_redis_url():
    assert SharedStore("t").backend == "memory"


def test_round_trips_a_value():
    s = SharedStore("t")
    s.set("k", {"a": 1})
    assert s.get("k") == {"a": 1}


def test_missing_key_returns_the_default():
    assert SharedStore("t").get("nope", "fallback") == "fallback"


def test_evicts_oldest_at_capacity():
    s = SharedStore("t", max_entries=3)
    for i in range(5):
        s.set(f"k{i}", i)
    assert len(s) == 3
    assert s.get("k0") is None
    assert s.get("k4") == 4


def test_reading_an_entry_makes_it_recently_used():
    s = SharedStore("t", max_entries=3)
    for i in range(3):
        s.set(f"k{i}", i)
    s.get("k0")            # refresh k0
    s.set("k3", 3)         # evicts the true oldest, k1
    assert s.get("k0") == 0
    assert s.get("k1") is None


def test_ttl_expires_entries():
    s = SharedStore("t", ttl=0.2)
    s.set("k", 1)
    assert s.get("k") == 1
    time.sleep(0.3)
    assert s.get("k") is None


def test_contains_and_delete():
    s = SharedStore("t")
    s.set("k", 1)
    assert "k" in s
    s.delete("k")
    assert "k" not in s


def test_mapping_access_matches_the_old_dict_call_sites():
    s = SharedStore("t")
    s["k"] = 5
    assert s["k"] == 5
    with pytest.raises(KeyError):
        s["missing"]


def test_pickle_serializer_round_trips_objects():
    """Motifs and HMMs are not JSON-serialisable."""
    class Thing:
        def __init__(self, v): self.v = v
        def __eq__(self, o): return isinstance(o, Thing) and o.v == self.v
    s = SharedStore("t")
    s.set("k", Thing(7))
    assert s.get("k") == Thing(7)


# --------------------------------------------------------------------------
# Redis backend, exercised through a stub
# --------------------------------------------------------------------------

class _FakeRedis:
    def __init__(self):
        self.data = {}
        self.ttls = {}
    def ping(self): return True
    def set(self, k, v): self.data[k] = v
    def setex(self, k, ttl, v): self.data[k] = v; self.ttls[k] = ttl
    def get(self, k): return self.data.get(k)
    def delete(self, k): self.data.pop(k, None)


@pytest.fixture
def fake_redis(monkeypatch):
    client = _FakeRedis()
    monkeypatch.setattr(shared_store, "_redis_client", client)
    monkeypatch.setattr(shared_store, "_redis_checked", True)
    return client


def test_redis_backend_is_reported(fake_redis):
    assert SharedStore("t").backend == "redis"


def test_values_are_namespaced_by_store_name(fake_redis):
    SharedStore("motifs").set("abc", 1)
    SharedStore("hmm").set("abc", 2)
    keys = sorted(fake_redis.data)
    assert any(k.endswith("motifs:abc") for k in keys)
    assert any(k.endswith("hmm:abc") for k in keys)
    assert SharedStore("motifs").get("abc") == 1
    assert SharedStore("hmm").get("abc") == 2


def test_ttl_is_passed_to_redis(fake_redis):
    SharedStore("t", ttl=60).set("k", 1)
    assert list(fake_redis.ttls.values()) == [60]


def test_json_serializer_stores_plain_json(fake_redis):
    SharedStore("t", serializer="json").set("k", {"a": 1})
    raw = next(iter(fake_redis.data.values()))
    assert raw == b'{"a": 1}'


def test_in_place_mutation_does_not_persist_under_redis(fake_redis):
    """The trap this design has to survive.

    Redis hands back a fresh object, so `store[k]['field'] = v` is discarded.
    update_entry() is the supported read/modify/write.
    """
    s = SharedStore("t")
    s.set("model", {"name": "m1"})

    s["model"]["name"] = "MUTATED"          # in-place: silently lost
    assert s.get("model")["name"] == "m1"

    assert s.update_entry("model", name="UPDATED") is True
    assert s.get("model")["name"] == "UPDATED"


def test_update_entry_reports_a_missing_key():
    assert SharedStore("t").update_entry("nope", a=1) is False


# --------------------------------------------------------------------------
# Fallback
# --------------------------------------------------------------------------

def test_unreachable_redis_falls_back_to_memory_with_a_warning(monkeypatch, caplog):
    """A misconfigured REDIS_URL must not stop the app booting."""
    import logging
    monkeypatch.setattr(shared_store, "REDIS_URL", "redis://127.0.0.1:1/0")
    shared_store.reset_backend_for_tests()
    with caplog.at_level(logging.WARNING, logger="utils.shared_store"):
        store = SharedStore("t")
        assert store.backend == "memory"
        store.set("k", 1)
        assert store.get("k") == 1
    assert "Falling back to per-process memory" in caplog.text


# --------------------------------------------------------------------------
# The routes actually use it
# --------------------------------------------------------------------------

def test_every_route_cache_is_a_shared_store():
    from routes.motif_routes import motif_cache
    from routes.blast_routes import blast_results_cache
    from routes.advanced_routes import _hmm_models
    from routes.status_routes import _CACHE
    for store in (motif_cache, blast_results_cache, _hmm_models, _CACHE):
        assert isinstance(store, SharedStore)


def test_motif_survives_a_round_trip_through_the_store(client):
    """End-to-end: create a motif, then read it back by id."""
    seqs = ["TATAAA", "TATAAT", "TATATA", "TACAAA"]
    created = client.post("/api/motifs/create", json={"sequences": seqs}).get_json()
    assert created["success"] is True
    info = client.get("/api/motifs/info").get_json()
    assert info["success"] is True
