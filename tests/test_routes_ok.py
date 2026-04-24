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
    # Phase-5 polish: Content-Security-Policy + Permissions-Policy
    csp = resp.headers.get("Content-Security-Policy") or ""
    assert "default-src 'self'" in csp
    assert "object-src 'none'" in csp
    assert "frame-ancestors 'none'" in csp
    # Must allowlist the CDN origins we actually load from
    assert "cdn.jsdelivr.net" in csp      # Bootstrap
    assert "cdnjs.cloudflare.com" in csp  # FontAwesome
    assert "3Dmol.org" in csp             # 3Dmol viewer
    # 3Dmol's surface rendering uses blob: Web Workers
    assert "worker-src 'self' blob:" in csp
    assert resp.headers.get("Permissions-Policy")


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


def test_mobile_nav_present(client):
    """Phase-5 polish: every page has the mobile hamburger + backdrop."""
    resp = client.get("/")
    body = resp.data
    assert b'id="mobileNavToggle"' in body
    assert b'id="sidebarBackdrop"' in body
    assert b'id="sidebarContainer"' in body


def test_r3_hub_uses_segmented_control(client):
    """R3: hubs switched from .nav-tabs to .segmented-control at the top
    level so users can visually distinguish hub nav from inner tabs."""
    for url, hub_id in [
        ("/patterns",  "patternsHubTabs"),
        ("/phylogeny", "phylogenyHubTabs"),
        ("/compare",   "compareHubTabs"),
    ]:
        body = client.get(url).data.decode()
        # Must use segmented-control on the hub bar
        assert 'class="segmented-control mb-3"' in body, f"{url} missing segmented-control"
        assert f'id="{hub_id}"' in body, f"{url} missing hub id {hub_id}"
        # And NOT nav-tabs at the hub level any more
        assert f'class="nav nav-tabs mb-3" id="{hub_id}"' not in body
        # Bootstrap tab-switching still attached
        assert 'data-bs-toggle="tab"' in body


def test_overhaul_d1_no_legacy_css(client):
    """D1 demolition: legacy CSS files must not exist and must not be linked."""
    body = client.get("/").data
    assert b"/static/css/custom.css" not in body
    assert b"/static/css/gradient-buttons.css" not in body
    assert b"/static/css/performance-fix.css" not in body
    assert b"/static/css/shell.css" not in body
    assert b"/static/css/components.css" not in body
    # Served files must 404 or return empty
    for css in ("custom.css", "gradient-buttons.css", "performance-fix.css",
                "shell.css", "components.css"):
        r = client.get("/static/css/" + css)
        assert r.status_code == 404, f"/static/css/{css} should be deleted"


def test_overhaul_d2_no_gradient_buttons(client):
    """D2: no template or JS should reference gradient-btn-* classes.
    Walk all served pages + the two JS files that emit button HTML."""
    routes = ["/", "/sequence", "/alignment", "/features", "/structure",
              "/restriction", "/blast", "/phylogeny", "/patterns"]
    for url in routes:
        body = client.get(url).data
        for legacy in (b"gradient-btn-primary", b"gradient-btn-secondary",
                       b"gradient-btn-info", b"gradient-btn-success",
                       b"gradient-btn-warning", b"gradient-btn-outline",
                       b"gradient-btn-small"):
            assert legacy not in body, f"legacy {legacy!r} in {url}"
    for js in ("database.js", "kegg.js"):
        r = client.get("/static/js/" + js)
        assert b"gradient-btn-" not in r.data, f"legacy class in {js}"


def test_overhaul_d3_split_pane_primitive(client):
    """D3: split-pane primitive lives in app.css and /sequence is an exemplar."""
    app_css = client.get("/static/css/app.css").data
    assert b".tool-page" in app_css
    assert b".tool-input-pane" in app_css
    assert b".tool-results-pane" in app_css

    seq = client.get("/sequence").data
    assert b'class="tool-page"' in seq
    assert b'class="tool-input-pane"' in seq
    assert b'class="tool-results-pane"' in seq
    assert b'class="tool-page-mobile-tabs"' in seq


def test_overhaul_d4_app_alerts_toast(client):
    """D4: #app-alerts toast region is on every page; utils.js targets it."""
    body = client.get("/").data
    assert b'id="app-alerts"' in body
    assert b'aria-live="polite"' in body

    utils = client.get("/static/js/utils.js").data
    # showAlert must target app-alerts, no longer query .card-body.p-4
    assert b"app-alerts" in utils
    assert b".card-body.p-4" not in utils


def test_overhaul_d5_no_legacy_tab_chrome(client):
    """D5: no template retains bg-gradient card-header or the inline
    <style>#XxxTabs overrides, or the text-white tab-button class soup."""
    for url in ["/sequence", "/alignment", "/blast", "/structure",
                "/restriction", "/phylogeny", "/patterns"]:
        body = client.get(url).data
        assert b"card-header bg-gradient" not in body, f"{url} has legacy gradient header"
        assert b"text-white border-0 py-3" not in body, f"{url} has legacy tab-button classes"


def test_overhaul_d6_sidebar_sub_tools(client):
    """D6: the sidebar surfaces sub-tools under each hub, and active tool
    highlights the parent hub AND the sub-tool with aria-current=page."""
    # Dashboard: no sub-tool active, but the sidebar contains several sub links
    body = client.get("/").data
    assert b'class="nav-sub"' in body
    assert b'class="nav-sublink' in body
    # Visit /alignment — its sub-link should be aria-current=page
    al = client.get("/alignment").data.decode()
    assert 'aria-current="page"' in al
    # And the Compare parent hub should ALSO be highlighted (active class on parent)
    import re
    # the matching nav-link active line
    hub_active_line = re.search(r'class="nav-link active"[^>]*>\s*<i[^>]*>[^<]*</i>\s*<span>(Compare)</span>', al)
    assert hub_active_line, "Compare hub not highlighted on /alignment"


def test_overhaul_d17_dark_mode(client):
    """D17: tokens.css defines the dark theme and topbar has a toggle."""
    tokens = client.get("/static/css/tokens.css").data
    assert b'html[data-theme="dark"]' in tokens
    body = client.get("/").data
    assert b'id="topbarThemeBtn"' in body
    # Inline FOUC-prevention script that applies saved/preferred theme
    assert b"nug-theme" in body
    assert b'data-theme' in body


def test_overhaul_d18_results_row_marker(client):
    """D18: primary hub tool pages have .results-row on their result row."""
    # Sample one from each of the three hub dirs
    for url, marker in [
        ("/alignment", b"alignmentResults"),
        ("/blast", b"blastResults"),
        ("/motifs", b"searchResults"),
        ("/popgen", b"popgenResults"),
        ("/clustering", b"clusteringCard"),
        ("/phylo", b"treeInfo"),
    ]:
        body = client.get(url).data
        assert marker in body, f"{url}: missing {marker}"
        # the row carrying that id should have the class applied
        assert b"results-row" in body, f"{url}: no .results-row marker"


def test_overhaul_d19_blast_inline_alignment(client):
    """D19: BLAST view-alignment renders inline (no tab-jump)."""
    js = client.get("/compare").data
    assert b"blastInlineAlignment" in js


def test_overhaul_d7_dashboard(client):
    """D7: dashboard is a tool grid + workspace recent + external-data list +
    quickstart. No gradient hero. No 'Core Tools' / 'Comparison & Evolution' section headers."""
    body = client.get("/").data
    assert b'class="dashboard"' in body
    assert b'class="tool-grid"' in body
    assert b'class="tool-tile"' in body
    assert b'class="quickstart"' in body
    assert b'id="dashWorkspace"' in body
    # Old layout markers must be gone
    assert b"display-4" not in body
    assert b"Core Tools</h" not in body
    assert b"Comparison & Evolution" not in body


def test_polish_a11y_skip_link_and_main(client):
    """Polish-2: skip-to-main-content link is the first focusable element,
    and the content region is a real <main> with role=main."""
    body = client.get("/").data
    # Skip link present and points at #main-content
    assert b'class="skip-link"' in body
    assert b'href="#main-content"' in body
    # Main element with id and role
    assert b'<main class="main-content" id="main-content" role="main"' in body


def test_polish_a11y_sidebar_aria_current(client):
    """Polish-2: the active sidebar link has aria-current=page."""
    body = client.get("/").data.decode()
    # On / the Dashboard link should be active
    assert 'aria-current="page"' in body


def test_polish_topbar_search_enabled(client):
    """Bug fix verification: topbar search is no longer disabled."""
    body = client.get("/").data
    assert b'id="topbarSearchInput"' in body
    # Must NOT be disabled
    assert b'<input type="search" id="topbarSearchInput"' in body
    # Old disabled stub should be gone
    topbar_chunk = body[body.index(b'topbar-search'):body.index(b'topbar-search') + 600]
    assert b'disabled' not in topbar_chunk


def test_r11_r17_extended_results_card_mounts(client):
    """R11–R17: structure sub-analyses, features/restriction extras,
    read variants, KEGG secondary, biodata, searchio secondary all
    call ResultsCard.mount().
    """
    structure = client.get("/static/js/structure.js").data
    for key in (b"contactsResults", b"interactionsResults",
                b"dsspResults", b"ramaResults", b"sasaResults",
                b"selectionResults", b"superimposeResults"):
        assert b"ResultsCard.mount('" + key + b"'" in structure, f"structure {key!r} not mounted"

    features = client.get("/static/js/features.js").data
    for key in (b"compoundResults", b"annotateResults"):
        assert b"ResultsCard.mount('" + key + b"'" in features, f"features {key!r} not mounted"

    pathway = client.get("/static/js/pathway.js").data
    assert b"ResultsCard.mount('pathwayVisualization'" in pathway

    restriction = client.get("/static/js/restriction.js").data
    assert b"ResultsCard.mount('compatibleResults'" in restriction

    biodata = client.get("/static/js/biodata.js").data
    assert b"ResultsCard.mount('iupacResults'" in biodata

    swiss = client.get("/static/js/swissprot.js").data
    assert b"ResultsCard.mount('readResults'" in swiss

    searchio = client.get("/static/js/searchio.js").data
    for key in (b"readResults", b"indexResults", b"convertResults",
                b"filterResults", b"writeResults"):
        assert b"ResultsCard.mount('" + key + b"'" in searchio, f"searchio {key!r} not mounted"

    unigene = client.get("/static/js/unigene.js").data
    assert b"ResultsCard.mount('readResults'" in unigene

    kegg = client.get("/static/js/kegg.js").data
    for key in (b"entryModalContent", b"linkResults",
                b"convertResults", b"infoResults"):
        assert b"ResultsCard.mount('" + key + b"'" in kegg or \
               b"ResultsCard.mount(containerId" in kegg, f"kegg {key!r} not mounted"

    database = client.get("/static/js/database.js").data
    for key in (b"globalResults", b"fetchResults", b"linkResults"):
        assert b"ResultsCard.mount('" + key + b"'" in database, f"database {key!r} not mounted"


def test_r10_extended_results_card_mounts(client):
    """R10: pathway network, structure quality, features extract,
    restriction map scripts each call ResultsCard.mount()."""
    pathway = client.get("/static/js/pathway.js")
    assert b"ResultsCard.mount('networkAnalysisResults'" in pathway.data, "pathway network not mounted"

    structure = client.get("/static/js/structure.js")
    assert b"ResultsCard.mount('qualityResults'" in structure.data, "structure quality not mounted"

    features = client.get("/static/js/features.js")
    assert b"ResultsCard.mount('extractResults'" in features.data, "features extract not mounted"

    restriction = client.get("/static/js/restriction.js")
    assert b"ResultsCard.mount('mapResultsCard'" in restriction.data, "restriction map not mounted"
    # Restriction map needs its new mount container in the template
    assert b'id="mapResultsCard"' in client.get("/restriction").data, "missing mapResultsCard container"


def test_r9_extended_results_card_mounts(client):
    """R9: biodata, features parsed, structure geometry, restriction
    advanced scripts each call ResultsCard.mount()."""
    biodata = client.get("/static/js/biodata.js")
    assert b"ResultsCard.mount('codonResults'" in biodata.data, "biodata not mounted"

    features = client.get("/static/js/features.js")
    assert b"ResultsCard.mount('parseResults'" in features.data, "features parse not mounted"

    structure = client.get("/static/js/structure.js")
    assert b"ResultsCard.mount('geometryResults'" in structure.data, "structure geometry not mounted"

    restriction = client.get("/static/js/restriction.js")
    assert b"ResultsCard.mount('advancedResults'" in restriction.data, "restriction advanced not mounted"


def test_r8_extended_results_card_mounts(client):
    """R8: pathway, database (Entrez), hmm, unigene scripts each call
    ResultsCard.mount() to render their primary results.
    """
    pathway = client.get("/static/js/pathway.js")
    assert b"ResultsCard.mount('systemAnalysisResults'" in pathway.data, "pathway not mounted"

    database = client.get("/static/js/database.js")
    assert b"ResultsCard.mount('searchResults'" in database.data, "database not mounted"

    # HMM lives in advanced.js (used by /hmm template)
    advanced = client.get("/static/js/advanced.js")
    assert b"ResultsCard.mount('hmmAnalysis'" in advanced.data, "hmm not mounted"

    unigene = client.get("/static/js/unigene.js")
    assert b"ResultsCard.mount('parseResults'" in unigene.data, "unigene not mounted"


def test_r7_extended_results_card_mounts(client):
    """R7: popgen, swissprot, kegg, searchio scripts each call
    ResultsCard.mount() to render their primary results.
    """
    # popgen lives inside /phylogeny hub
    phylogeny = client.get("/phylogeny").data
    assert b"ResultsCard.mount('popgenResults'" in phylogeny, "popgen not mounted"

    # swissprot, kegg, searchio have their own static JS files
    swiss = client.get("/static/js/swissprot.js")
    assert b"ResultsCard.mount('parseResults'" in swiss.data, "swissprot not mounted"

    kegg = client.get("/static/js/kegg.js")
    assert b"ResultsCard.mount(containerId" in kegg.data, "kegg helper not mounted"

    searchio = client.get("/static/js/searchio.js")
    assert b"ResultsCard.mount('parseResults'" in searchio.data, "searchio not mounted"


def test_r6_extended_results_card_mounts(client):
    """R6: BLAST, Motifs, Clustering, Phylogeny scripts each call
    ResultsCard.mount() to render their primary results.
    """
    # BLAST lives inside the /compare hub
    compare = client.get("/compare").data
    assert b"ResultsCard.mount('blastResults'" in compare, "BLAST not mounted"

    # Motifs lives inside /patterns hub
    patterns = client.get("/patterns").data
    assert b"ResultsCard.mount('searchResults'" in patterns, "Motif search not mounted"

    # Clustering + Phylogeny live inside /phylogeny hub
    phylogeny = client.get("/phylogeny").data
    assert b"ResultsCard.mount('clusteringCard'" in phylogeny, "Clustering not mounted"
    assert b"ResultsCard.mount('treeInfo'" in phylogeny, "Phylo tree not mounted"
    # Clustering needs its mount container
    assert b'id="clusteringCard"' in phylogeny, "Missing clusteringCard container"


def test_r4_structure_viewer_first_layout(client):
    """R4: /structure is viewer-first — 3D viewer card + vertical
    analyses nav with 10 items + content card below."""
    body = client.get("/structure").data
    # Viewer card and its 3Dmol canvas container
    assert b'id="view3dContainer"' in body
    assert b'structure-viewer-canvas' in body
    # Vertical analyses nav (list-group style)
    assert b'structure-analyses-nav' in body
    assert b'id="structureTabs"' in body
    # All 10 analyses still wired
    for analysis in (b"parse", b"geometry", b"quality", b"contacts",
                     b"interactions", b"dssp", b"ramachandran", b"sasa",
                     b"selection", b"superimpose"):
        assert b'id="' + analysis + b'-tab"' in body, f"missing {analysis!r}"
        assert b'id="' + analysis + b'"' in body, f"missing pane {analysis!r}"
    # Legacy element kept for structure.js compat
    assert b'id="advancedAnalysisBtn"' in body


def test_r2_results_card_served(client):
    """R2: results_card.js is served and linked on every page."""
    resp = client.get("/static/js/results_card.js")
    assert resp.status_code == 200
    assert b"ResultsCard.mount" in resp.data or b"mount: mount" in resp.data
    idx = client.get("/")
    assert b"/static/js/results_card.js" in idx.data


def test_r1_shell_structure(client):
    """R1 redesign: every page has the new topbar, Workspace button,
    sidebar partial with 6 hubs, and shell.css."""
    resp = client.get("/")
    body = resp.data
    # Topbar
    assert b'class="app-topbar"' in body
    assert b'id="topbarWorkspaceBtn"' in body
    assert b'id="topbarWorkspaceBadge"' in body
    # Sidebar shows 6 hubs + dashboard, no section headers
    assert b'>Dashboard<' in body
    assert b'>Sequences<' in body
    assert b'>Compare<' in body
    assert b'>Patterns<' in body
    assert b'>Structure<' in body
    assert b'>Phylogeny<' in body
    assert b'>External data<' in body
    assert b'nav-section-header' not in body  # section headers removed
    # app.css is now the single source of truth (D1 demolition)
    assert b'/static/css/app.css' in body


def test_r1_sidebar_activates_correct_hub(client):
    """Visiting /alignment highlights the Compare hub nav link (not the old
    individual-tool link). This guards the R1 6-hub IA."""
    import re
    resp = client.get("/alignment")
    body = resp.data.decode()
    # Find the active nav-link and confirm it's the Compare hub
    m = re.search(r'class="nav-link active"[^>]*>\s*<i[^>]*>[^<]*</i>\s*<span>([^<]+)</span>', body)
    assert m, "no active hub nav link on /alignment"
    assert m.group(1).strip() == "Compare"


def test_structure_3d_viewer_wired(client):
    """Phase-5 polish: /structure has the persistent 3Dmol viewer +
    full control bar (sample loader, style/color selects, surface/spin
    toggles, reset and PNG download)."""
    resp = client.get("/structure")
    body = resp.data
    assert b'id="view3dContainer"' in body
    assert b'3Dmol.org' in body or b'3Dmol-min.js' in body
    # Viewer controls exposed via data-action delegation
    for action in [
        b'data-action="view3dLoadSample"',
        b'data-action="view3dReset"',
        b'data-action="view3dToggleSpin"',
        b'data-action="view3dToggleSurface"',
        b'data-action="view3dDownloadPng"',
    ]:
        assert action in body, f'{action!r} not in /structure'
    # Public viewer API surfaces for per-tab overlays
    assert b'NuGenStructureViewer' in body


def test_base_template_loads_design_tokens(client):
    """tokens.css and app.css must be linked, and served non-empty.
    D1 demolition: app.css is now the single visual source of truth —
    custom.css / gradient-buttons.css / performance-fix.css / shell.css /
    components.css are deleted. Anyone re-introducing them or removing
    tokens.css / app.css breaks the cascade."""
    resp = client.get("/")
    body = resp.data
    assert b"/static/css/tokens.css" in body
    assert b"/static/css/app.css" in body
    # Legacy CSS files must NOT be linked
    assert b"/static/css/custom.css" not in body
    assert b"/static/css/gradient-buttons.css" not in body
    assert b"/static/css/performance-fix.css" not in body
    assert b"/static/css/shell.css" not in body
    assert b"/static/css/components.css" not in body

    tokens = client.get("/static/css/tokens.css")
    assert tokens.status_code == 200
    assert b"--slate-" in tokens.data
    assert b"--color-primary" in tokens.data

    app_css = client.get("/static/css/app.css")
    assert app_css.status_code == 200
    # Core component classes must be present in the consolidated file
    assert b".app-topbar" in app_css.data
    assert b".sidebar-nav" in app_css.data
    assert b".results-card" in app_css.data
    assert b".segmented-control" in app_css.data


@pytest.mark.parametrize("url,hub,expected_tab", HUB_TAB_PRESELECTION)
def test_hub_preselects_correct_tab(client, url, hub, expected_tab):
    """Phase-3 invariant: each legacy URL (e.g. /motifs) must render the
    hub view with the correct tab active. Failing this test means someone
    broke the URL-alias → tab-state wiring and a user who bookmarked the
    old URL would land on the wrong sub-tool."""
    resp = client.get(url)
    assert resp.status_code == 200, f"{url} returned {resp.status_code}"
    body = resp.data.decode()
    # The active hub tab has class="... active" and id="hub-tab-XXX" on the
    # same element. Order of attributes varies (R3 swapped nav-link → segmented-
    # item). Look for id first, then the active class on the same tag.
    active = None
    for m in re.finditer(r'<button[^>]*>', body):
        tag = m.group(0)
        id_m = re.search(r'id="hub-tab-(\w+)"', tag)
        if id_m and re.search(r'class="[^"]*\bactive\b', tag):
            active = id_m.group(1)
            break
    assert active is not None, f"{url}: no active hub tab found in body"
    assert active == expected_tab, (
        f"{url}: expected tab={expected_tab!r} but found {active!r}"
    )
