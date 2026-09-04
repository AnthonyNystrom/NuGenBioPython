"""Threaded-worker safety.

Three process-global hazards made threaded workers unsafe:

  1. remote_timeout mutated the default socket timeout, so concurrent calls
     stomped each other's window and could restore the wrong value.
  2. Figures came from pyplot, which registers them in a global dict; app.py
     swept it with plt.close('all') on every teardown — destroying figures
     other in-flight requests were still drawing.
  3. The outbound RateGate (already lock-protected).

These pin the fixes so the app cannot silently regress to sync-only.
"""
import socket
import threading

import pytest


# --------------------------------------------------------------------------
# 1. socket timeout
# --------------------------------------------------------------------------

def test_concurrent_remote_timeout_does_not_corrupt():
    from utils.request_helpers import remote_timeout
    import random
    import time

    errors = []
    original = socket.getdefaulttimeout()

    def worker(value):
        try:
            with remote_timeout(value):
                if socket.getdefaulttimeout() != value:
                    errors.append(f"saw {socket.getdefaulttimeout()}, wanted {value}")
                time.sleep(random.uniform(0.001, 0.01))
                if socket.getdefaulttimeout() != value:
                    errors.append(f"changed under us: wanted {value}")
        except Exception as exc:
            errors.append(str(exc))

    threads = [threading.Thread(target=worker, args=(v,))
               for v in (5, 10, 15, 20, 25) * 4]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, errors
    assert socket.getdefaulttimeout() == original


def test_remote_timeout_is_reentrant():
    """A nested call must not deadlock on its own lock."""
    from utils.request_helpers import remote_timeout
    with remote_timeout(10):
        with remote_timeout(5):
            assert socket.getdefaulttimeout() == 5
        assert socket.getdefaulttimeout() == 10


# --------------------------------------------------------------------------
# 2. figures never enter pyplot's global registry
# --------------------------------------------------------------------------

def test_helper_figures_are_not_registered_with_pyplot():
    import matplotlib.pyplot as plt
    from utils.plot_helpers import subplots

    before = set(plt.get_fignums())
    fig, ax = subplots(figsize=(3, 2))
    ax.plot([0, 1], [1, 0])
    assert set(plt.get_fignums()) == before, \
        "figure leaked into pyplot's global registry"


def test_no_renderer_uses_pyplot_subplots(repo_root):
    """Guard: plt.subplots() reintroduces the leak and the unsafe sweep."""
    import pathlib
    import re
    offenders = []
    for folder in ("routes", "utils"):
        for path in pathlib.Path(repo_root, folder).rglob("*.py"):
            if path.name == "plot_helpers.py":
                continue        # documents the problem in prose
            for i, line in enumerate(path.read_text().splitlines(), 1):
                if line.strip().startswith("#"):
                    continue
                if re.search(r"\bplt\.(subplots|figure)\(", line):
                    offenders.append(f"{path.name}:{i}")
    assert not offenders, (
        "use utils.plot_helpers.subplots() instead of pyplot: %s" % offenders)


def test_app_no_longer_sweeps_all_figures(repo_root):
    """plt.close('all') closed other requests' in-flight figures."""
    import pathlib
    source = pathlib.Path(repo_root, "app.py").read_text()
    code = "\n".join(l for l in source.splitlines()
                     if not l.strip().startswith("#"))
    assert "plt.close" not in code
    assert "teardown_request" not in code


def test_many_renders_leave_the_registry_empty():
    import matplotlib.pyplot as plt
    from utils.plot_helpers import subplots, figure_to_svg_data_url

    before = set(plt.get_fignums())
    for _ in range(20):
        fig, ax = subplots(figsize=(2, 2))
        ax.plot([0, 1], [0, 1])
        figure_to_svg_data_url(fig)
    assert set(plt.get_fignums()) == before


def test_concurrent_renders_are_identical_and_uncorrupted():
    """The scenario plt.close('all') used to break."""
    import base64
    import re
    from utils.plot_helpers import subplots, figure_to_svg_data_url

    results = []
    lock = threading.Lock()

    def render():
        fig, ax = subplots(figsize=(4, 3))
        ax.plot([0, 1, 2], [1, 3, 2])
        ax.set_title("concurrent")
        url = figure_to_svg_data_url(fig, theme="light")
        svg = base64.b64decode(url.split(",", 1)[1])
        svg = re.sub(rb"<dc:date>.*?</dc:date>", b"", svg)
        svg = re.sub(rb'(id="|url\(#|xlink:href="#)[A-Za-z0-9]+', rb"\1ID", svg)
        with lock:
            results.append(svg)

    threads = [threading.Thread(target=render) for _ in range(16)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(results) == 16
    assert len({r for r in results}) == 1, \
        "concurrent renders produced differing output"


# --------------------------------------------------------------------------
# 3. rate gate
# --------------------------------------------------------------------------

def test_rate_gate_is_lock_protected():
    from utils.request_helpers import RateGate
    import time

    gate = RateGate("t", 0.02)
    stamps = []
    lock = threading.Lock()

    def worker():
        gate.wait()
        with lock:
            stamps.append(time.monotonic())

    threads = [threading.Thread(target=worker) for _ in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    stamps.sort()
    gaps = [b - a for a, b in zip(stamps, stamps[1:])]
    assert all(g >= 0.015 for g in gaps), gaps
