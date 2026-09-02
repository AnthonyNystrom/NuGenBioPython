"""3Dmol viewer configuration.

Crystallographic waters are HETATM records. Styling every HETATM as a stick
drew all of them — 350 in 4INS — as red specks around the protein. 3Dmol also
infers bonds by distance, and 4INS contains 49 water pairs closer than 2 A
(alternate, partially-occupied sites; the closest 0.83 A apart), so it drew
bonds *between* waters, rendering them as short red rods that read as a
graphics fault rather than data.
"""
import pathlib

import pytest


@pytest.fixture
def viewer_source(repo_root):
    return pathlib.Path(repo_root, "templates", "structure.html").read_text()


def test_waters_are_excluded_from_the_hetatm_stick_style(viewer_source):
    """The HETATM rule must not catch waters."""
    assert "hetflag: true, resn: WATER_RESN, invert: true" in viewer_source, \
        "the heteroatom stick style must exclude water residues"


def test_bare_hetflag_stick_rule_is_gone(viewer_source):
    """Guard against the original rule coming back."""
    assert "setStyle({ hetflag: true }," not in viewer_source


def test_water_residue_names_cover_the_common_variants(viewer_source):
    for resn in ("HOH", "WAT", "DOD"):
        assert f"'{resn}'" in viewer_source, f"{resn} not in WATER_RESN"


def test_waters_are_hidden_by_default(viewer_source):
    assert "let showWaters = false;" in viewer_source, \
        "waters must default to hidden, as in PyMOL/Chimera/Mol*"


def test_waters_render_as_spheres_not_sticks(viewer_source):
    """A water is one atom; sphere style avoids the inferred-bond artefact."""
    idx = viewer_source.index("if (showWaters)")
    block = viewer_source[idx:idx + 400]
    # Strip // comments first — the explanatory comment mentions "sticks".
    code = "\n".join(line.split("//")[0] for line in block.splitlines())
    assert "sphere" in code
    assert "stick" not in code


def test_a_toggle_exists_and_is_wired(viewer_source):
    assert 'data-action="view3dToggleWaters"' in viewer_source
    assert "window.view3dToggleWaters = function" in viewer_source
    assert 'id="view3dWatersBtn"' in viewer_source


def test_non_water_heteroatoms_are_still_shown(viewer_source):
    """Ligands and ions — the zincs in insulin — must not be hidden too."""
    idx = viewer_source.index("function applyStyle()")
    block = viewer_source[idx:idx + 900]
    assert "hetflag: true" in block and "invert: true" in block


def test_structure_page_still_renders(client):
    resp = client.get("/structure")
    assert resp.status_code == 200
    assert b"view3dWatersBtn" in resp.data
