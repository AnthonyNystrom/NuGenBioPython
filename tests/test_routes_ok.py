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
    # shell.css linked after components.css
    assert b'/static/css/shell.css' in body


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
