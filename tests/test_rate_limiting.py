"""Outbound rate limiting for external biology services.

NCBI enforces 3 requests/second (10 with an API key) by blocking the
offending IP. Nothing limited how fast this app called out, so a few
concurrent users could get the whole deployment blocked with no signal
beyond sudden failures.
"""
import os
import threading
import time

import pytest

from utils.request_helpers import (
    RateGate, throttle, remote_timeout, _ncbi_interval, _GATES,
)


@pytest.fixture(autouse=True)
def _reset_gates(monkeypatch):
    monkeypatch.delenv("ENTREZ_API_KEY", raising=False)
    for gate in _GATES.values():
        gate._last_call = 0.0
    yield


def test_first_call_is_not_delayed():
    gate = RateGate("t", 0.05)
    assert gate.wait() == 0.0


def test_second_call_is_delayed_by_the_interval():
    gate = RateGate("t", 0.05)
    gate.wait()
    start = time.monotonic()
    gate.wait()
    assert time.monotonic() - start >= 0.045


def test_calls_spaced_beyond_the_interval_are_not_delayed():
    gate = RateGate("t", 0.02)
    gate.wait()
    time.sleep(0.03)
    assert gate.wait() == 0.0


def test_sustained_rate_stays_under_the_ceiling():
    """The property that matters: N calls take at least (N-1)*interval."""
    interval = 0.02
    gate = RateGate("t", interval)
    calls = 6
    start = time.monotonic()
    for _ in range(calls):
        gate.wait()
    elapsed = time.monotonic() - start
    assert elapsed >= (calls - 1) * interval * 0.9
    observed_rate = calls / elapsed
    assert observed_rate <= (1 / interval) * 1.2


def test_gate_is_thread_safe():
    """Concurrent callers must not both slip through the same window."""
    gate = RateGate("t", 0.02)
    timestamps = []
    lock = threading.Lock()

    def worker():
        gate.wait()
        with lock:
            timestamps.append(time.monotonic())

    threads = [threading.Thread(target=worker) for _ in range(5)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    timestamps.sort()
    gaps = [b - a for a, b in zip(timestamps, timestamps[1:])]
    assert all(g >= 0.015 for g in gaps), f"calls bunched together: {gaps}"


def test_ncbi_interval_respects_the_api_key(monkeypatch):
    """NCBI allows 10 req/s with a key, 3 without."""
    monkeypatch.delenv("ENTREZ_API_KEY", raising=False)
    assert _ncbi_interval() == pytest.approx(0.34)
    assert 1 / _ncbi_interval() <= 3.0

    monkeypatch.setenv("ENTREZ_API_KEY", "abc123")
    assert _ncbi_interval() == pytest.approx(0.11)
    assert 1 / _ncbi_interval() <= 10.0


def test_unknown_service_passes_through():
    assert throttle("not-a-service") == 0.0


def test_remote_timeout_applies_the_gate():
    """The gate must fire through the context manager the call sites use."""
    _GATES["kegg"]._last_call = time.monotonic()
    start = time.monotonic()
    with remote_timeout(5, service="kegg"):
        pass
    assert time.monotonic() - start >= 0.3


def test_remote_timeout_without_service_is_not_throttled():
    start = time.monotonic()
    with remote_timeout(5):
        pass
    assert time.monotonic() - start < 0.05


def test_remote_timeout_still_restores_the_socket_timeout():
    import socket
    original = socket.getdefaulttimeout()
    try:
        with remote_timeout(7, service="kegg"):
            assert socket.getdefaulttimeout() == 7
        assert socket.getdefaulttimeout() == original
    finally:
        socket.setdefaulttimeout(original)


def test_socket_timeout_restored_even_on_exception():
    import socket
    original = socket.getdefaulttimeout()
    try:
        with pytest.raises(ValueError):
            with remote_timeout(7, service="kegg"):
                raise ValueError("boom")
        assert socket.getdefaulttimeout() == original
    finally:
        socket.setdefaulttimeout(original)


def test_every_outbound_call_site_declares_a_service(repo_root):
    """Guard: a new ungated remote_timeout() must not slip in."""
    import pathlib
    import re
    offenders = []
    for path in list(pathlib.Path(repo_root, "routes").rglob("*.py")) + \
                list(pathlib.Path(repo_root, "utils").rglob("*.py")):
        if path.name == "request_helpers.py":
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"remote_timeout\(", line) and "service=" not in line:
                offenders.append(f"{path.name}:{i}")
    assert not offenders, f"ungated outbound calls: {offenders}"
