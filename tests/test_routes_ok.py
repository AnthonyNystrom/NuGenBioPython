"""
Smoke test for every page URL in the sidebar.

Locks down the URL surface: any redesign that makes one of these 4xx/5xx will
fail CI. Phase 3 hub consolidation uses this file to prove old URLs stay
reachable even after they become aliases.
"""
import pytest

from tests.fixtures.inventory import PAGE_URLS


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
    """utils.js and formatters.js must be linked on every page via base.html."""
    resp = client.get("/sequence")
    body = resp.data
    assert b"/static/js/utils.js" in body
    assert b"/static/js/formatters.js" in body


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
