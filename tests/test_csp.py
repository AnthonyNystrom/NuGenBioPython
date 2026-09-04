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

# --------------------------------------------------------------------------
# style-src
# --------------------------------------------------------------------------

def test_style_src_no_longer_allows_unsafe_inline(client):
    """139 inline style attributes became classes or data-css."""
    style_src = re.search(r"style-src ([^;]*)", _csp(client)).group(1)
    assert "'unsafe-inline'" not in style_src, style_src


def test_style_src_carries_the_nonce(client):
    style_src = re.search(r"style-src ([^;]*)", _csp(client)).group(1)
    assert re.search(r"'nonce-[A-Za-z0-9_-]{16,}'", style_src), style_src


def test_no_style_attributes_remain_in_markup(repo_root):
    """A style= attribute in markup is refused outright under this CSP.

    Verified in Chromium: markup style= and setAttribute('style', ...) are
    blocked, while el.style.prop = x and el.style.cssText are exempt because
    the CSSOM is not covered by style-src. Dynamic values therefore travel as
    data-css and are applied through the CSSOM.
    """
    offenders = []
    for path in _source_files(repo_root):
        text = path.read_text()
        for i, line in enumerate(text.splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith(("//", "*", "#")):
                continue
            # (?<![\w-]) rather than \s: an earlier version of this pattern
            # missed a style attribute written as `style="..."` inside a JS
            # template literal, which would have violated CSP the moment a
            # restriction map rendered.
            if re.search(r'(?<![\w-])style\s*=\s*["\']', line):
                offenders.append(f"{path.name}:{i}  {stripped[:60]}")
    assert not offenders, (
        "inline style attributes are blocked; use a utility class, or "
        f"data-css for a value only known at render time: {offenders}")


def test_no_style_set_via_setattribute(repo_root):
    """setAttribute('style', ...) is blocked too, unlike el.style.prop."""
    offenders = []
    for path in _source_files(repo_root):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip().startswith(("//", "*", "#")):
                continue          # prose describing the rule, not code
            if re.search(r"""setAttribute\(\s*['"]style['"]""", line):
                offenders.append(f"{path.name}:{i}")
    assert not offenders, offenders


def test_every_inline_style_block_is_nonced(repo_root):
    offenders = []
    for path in pathlib.Path(repo_root, "templates").rglob("*.html"):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"<style\s*>", line):
                offenders.append(f"{path.name}:{i}")
    assert not offenders, f"inline <style> without a nonce: {offenders}"


def test_data_css_applier_exists(repo_root):
    """The CSSOM shim dynamic values depend on."""
    utils = pathlib.Path(repo_root, "static", "js", "utils.js").read_text()
    assert "applyDataCss" in utils
    assert "cssText" in utils
    assert "MutationObserver" in utils
