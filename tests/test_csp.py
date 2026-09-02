"""Content-Security-Policy hardening.

script-src previously carried 'unsafe-inline', which makes the directive
almost decorative: any injected <script> executes. Dropping it required two
things — every inline <script> block carrying a per-request nonce, and the
~54 inline onclick= handlers becoming data-action dispatches, because inline
event handlers can never be covered by a nonce.
"""
import pathlib
import re

import pytest


def _csp(client, path="/"):
    return client.get(path).headers["Content-Security-Policy"]


def test_script_src_no_longer_allows_unsafe_inline(client):
    csp = _csp(client)
    script_src = re.search(r"script-src ([^;]*)", csp).group(1)
    assert "'unsafe-inline'" not in script_src, script_src


def test_script_src_carries_a_nonce(client):
    script_src = re.search(r"script-src ([^;]*)", _csp(client)).group(1)
    assert re.search(r"'nonce-[A-Za-z0-9_-]{16,}'", script_src), script_src


def test_nonce_is_unique_per_request(client):
    """A reused nonce is as good as no nonce."""
    nonces = set()
    for _ in range(5):
        nonces.add(re.search(r"'nonce-([^']+)'", _csp(client)).group(1))
    assert len(nonces) == 5, "nonce must be minted per request"


def test_page_nonce_matches_the_header(client):
    resp = client.get("/")
    header_nonce = re.search(r"'nonce-([^']+)'",
                             resp.headers["Content-Security-Policy"]).group(1)
    body = resp.data.decode()
    page_nonces = set(re.findall(r'<script nonce="([^"]+)"', body))
    assert page_nonces, "no nonced inline scripts found on the page"
    assert page_nonces == {header_nonce}, \
        "inline scripts must carry the same nonce as the header"


def test_other_directives_are_preserved(client):
    csp = _csp(client)
    for directive in ("default-src 'self'", "object-src 'none'",
                      "base-uri 'self'", "frame-ancestors 'none'",
                      "worker-src 'self' blob:"):
        assert directive in csp, directive


def test_external_script_hosts_still_allowed(client):
    """3Dmol and the CDNs must keep working."""
    script_src = re.search(r"script-src ([^;]*)", _csp(client)).group(1)
    for host in ("https://3Dmol.org", "https://cdn.jsdelivr.net",
                 "https://cdnjs.cloudflare.com"):
        assert host in script_src, host


# --------------------------------------------------------------------------
# Source guards
# --------------------------------------------------------------------------

def _source_files(repo_root):
    for pattern in ("templates/**/*.html", "static/js/*.js"):
        yield from pathlib.Path(repo_root).glob(pattern)


# Matches ANY inline event handler attribute. An earlier version of this test
# enumerated only onclick/onchange/onsubmit/oninput and so missed a live
# onkeyup= search box and an onerror= image fallback, both of which the
# nonce-based CSP silently blocks.
_INLINE_HANDLER = re.compile(r'\son[a-z]+\s*=\s*"')


def test_no_inline_event_handlers_remain(repo_root):
    """Inline handlers cannot be nonced; one reintroduced breaks silently.

    They fail only when the event fires, not on page load, so a page-load
    smoke test will not catch them.
    """
    offenders = []
    for path in _source_files(repo_root):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("//") or stripped.startswith("#"):
                continue
            match = _INLINE_HANDLER.search(line)
            if match:
                offenders.append(f"{path.name}:{i} ({match.group().strip()})")
    assert not offenders, (
        "inline event handlers are blocked by the nonce-based CSP; "
        f"use data-action / data-action-change / data-action-submit / "
        f"data-action-input instead: {offenders}")


def test_every_inline_script_block_is_nonced(repo_root):
    offenders = []
    for path in pathlib.Path(repo_root, "templates").rglob("*.html"):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"<script\s*>", line):
                offenders.append(f"{path.name}:{i}")
    assert not offenders, f"inline <script> without a nonce: {offenders}"
