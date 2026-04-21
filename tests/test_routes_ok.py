"""
Smoke test for every page URL in the sidebar.

Locks down the URL surface: any redesign that makes one of these 4xx/5xx will
fail CI. Phase 3 hub consolidation uses this file to prove old URLs stay
reachable even after they become aliases.
"""
import re

import pytest

from tests.fixtures.inventory import PAGE_URLS, HUB_TAB_PRESELECTION


@pytest.mark.parametrize("url", PAGE_URLS)
def test_page_returns_html(client, url):
    """Every known page URL must return 200 with an HTML body."""
    resp = client.get(url)
    assert resp.status_code == 200, f"{url} returned {resp.status_code}"
    body = resp.data
    assert b"<html" in body.lower() or b"<!doctype" in body.lower(), (
        f"{url} response does not look like HTML"
    )


def test_security_headers_applied(client):
    """Security headers from app.py's after_request hook must reach the browser."""
    resp = client.get("/")
    assert resp.headers.get("X-Content-Type-Options") == "nosniff"
    assert resp.headers.get("X-Frame-Options") == "DENY"
    assert resp.headers.get("Referrer-Policy") == "strict-origin-when-cross-origin"


def test_base_template_loads_js_bundle(client):
    """utils.js, formatters.js, and workspace.js must be linked on every
    page via base.html."""
    resp = client.get("/sequence")
    body = resp.data
    assert b"/static/js/utils.js" in body
    assert b"/static/js/formatters.js" in body
    assert b"/static/js/workspace.js" in body


def test_workspace_js_served(client):
    """workspace.js exposes the Phase-4 cross-tool clipboard API."""
    resp = client.get("/static/js/workspace.js")
    assert resp.status_code == 200
    body = resp.data
    assert b"Workspace.add" in body or b"add = add" in body or b"add: add" in body
    assert b"sessionStorage" in body
    assert b"ws-panel" in body


def test_base_template_loads_design_tokens(client):
    """tokens.css and components.css must be linked, and served non-empty.
    Guards Phase 1 of the redesign: anyone removing these breaks the token
    cascade and regresses the sidebar-leak fix."""
    resp = client.get("/")
    body = resp.data
    assert b"/static/css/tokens.css" in body
    assert b"/static/css/components.css" in body

    tokens = client.get("/static/css/tokens.css")
    assert tokens.status_code == 200
    assert b"--slate-" in tokens.data
    assert b"--color-primary" in tokens.data

    components = client.get("/static/css/components.css")
    assert components.status_code == 200
    # The sidebar-nav rules are defined in custom.css (root-cause scoping);
    # components.css only provides the low-specificity fallback for other
    # .nav-link usages. Check that the selector RULE isn't in components.
    assert b".sidebar-nav .nav-link {" not in components.data
    assert b":where(body) .nav-link" in components.data


@pytest.mark.parametrize("url,hub,expected_tab", HUB_TAB_PRESELECTION)
def test_hub_preselects_correct_tab(client, url, hub, expected_tab):
    """Phase-3 invariant: each legacy URL (e.g. /motifs) must render the
    hub view with the correct tab active. Failing this test means someone
    broke the URL-alias → tab-state wiring and a user who bookmarked the
    old URL would land on the wrong sub-tool."""
    resp = client.get(url)
    assert resp.status_code == 200, f"{url} returned {resp.status_code}"
    body = resp.data.decode()
    # Find the active tab id in the hub tab bar
    m = re.search(
        r'id="hub-tab-(\w+)"[^>]*class="[^"]*active|class="nav-link active"[^>]*id="hub-tab-(\w+)"',
        body,
    )
    assert m, f"{url}: no active hub tab found in body"
    active = m.group(1) or m.group(2)
    assert active == expected_tab, (
        f"{url}: expected tab={expected_tab!r} but found {active!r}"
    )
