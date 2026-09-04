"""Ramachandran classification.

The previous in-route classifier could never report an outlier: its final
branch tested `-180 <= phi <= 180 and -180 <= psi <= 180`, true of every
possible pair of angles. Every structure scored zero outliers regardless of
geometry — false confidence in a validation feature.
"""
import pytest

from utils.structure_plots import (
    classify_phi_psi, classify_all, render_ramachandran, _delta,
)

CANONICAL_FAVORED = {
    'right-handed alpha helix': (-63, -43),
    'antiparallel beta sheet': (-139, 135),
    'parallel beta sheet': (-119, 113),
    'polyproline II': (-75, 145),
    '3-10 helix': (-49, -26),
}

IMPOSSIBLE = {
    'steric clash': (120, 120),
    'fully eclipsed': (0, 0),
    'forbidden quadrant': (150, -150),
}


@pytest.mark.parametrize("name,angles", CANONICAL_FAVORED.items())
def test_canonical_conformations_are_favored(name, angles):
    assert classify_phi_psi(*angles, 'ALA') == 'favored', name


@pytest.mark.parametrize("name,angles", IMPOSSIBLE.items())
def test_impossible_geometry_is_an_outlier(name, angles):
    """The regression: outliers must be reachable at all."""
    assert classify_phi_psi(*angles, 'ALA') == 'outlier', name


def test_outliers_are_reachable():
    """Guard against the always-allowed branch coming back."""
    grid = [classify_phi_psi(p, s, 'ALA')
            for p in range(-180, 181, 10)
            for s in range(-180, 181, 10)]
    assert 'outlier' in grid, "no angle anywhere is classified as an outlier"
    # A large share of the plot is genuinely disallowed for non-Gly residues.
    assert grid.count('outlier') / len(grid) > 0.3


def test_glycine_tolerates_positive_phi():
    """Gly has no side chain; positive-phi conformations are normal for it."""
    phi, psi = 106.3, 7.3   # a real GLY from static/sample_protein.pdb
    assert classify_phi_psi(phi, psi, 'ALA') == 'outlier'
    assert classify_phi_psi(phi, psi, 'GLY') in ('favored', 'allowed')


# Angles that are strained even for glycine. Note (150, -150) is NOT one of
# them: it mirrors onto the beta basin, and real glycine Ramachandran plots
# genuinely show density in that mirrored region.
IMPOSSIBLE_EVEN_FOR_GLYCINE = {
    'steric clash': (120, 120),
    'fully eclipsed': (0, 0),
}


@pytest.mark.parametrize("name,angles", IMPOSSIBLE_EVEN_FOR_GLYCINE.items())
def test_glycine_still_cannot_occupy_strained_geometry(name, angles):
    assert classify_phi_psi(*angles, 'GLY') == 'outlier', name


def test_glycine_is_more_permissive_than_other_residues():
    """The biologically meaningful property: Gly's region is broader, but
    it is not unrestricted."""
    grid = [(p, s) for p in range(-180, 181, 15) for s in range(-180, 181, 15)]
    ala = sum(1 for p, s in grid if classify_phi_psi(p, s, 'ALA') == 'outlier')
    gly = sum(1 for p, s in grid if classify_phi_psi(p, s, 'GLY') == 'outlier')
    assert gly < ala, "glycine must tolerate more conformations"
    assert gly / len(grid) > 0.2, "glycine must still exclude strained geometry"


def test_angle_delta_wraps():
    assert _delta(170, -170) == pytest.approx(-20)
    assert _delta(-170, 170) == pytest.approx(20)


def test_classify_all_counts_and_annotates():
    data = [
        {'phi': -63, 'psi': -43, 'resname': 'ALA'},
        {'phi': 120, 'psi': 120, 'resname': 'ALA'},
    ]
    counts, annotated = classify_all(data)
    assert counts['favored'] == 1
    assert counts['outlier'] == 1
    assert [e['region'] for e in annotated] == ['favored', 'outlier']


def test_real_structure_scores_as_well_refined(repo_root):
    """A high-resolution structure must land >95% favored+allowed."""
    import os
    import numpy as np
    from Bio.PDB import PDBParser, PPBuilder

    path = os.path.join(repo_root, 'static', 'sample_protein.pdb')
    structure = PDBParser(QUIET=True).get_structure('s', path)
    data = []
    for chain in structure[0]:
        for pp in PPBuilder().build_peptides(chain):
            residues = list(pp)
            for i, (phi, psi) in enumerate(pp.get_phi_psi_list()):
                if phi and psi:
                    data.append({'phi': np.degrees(phi), 'psi': np.degrees(psi),
                                 'resname': residues[i].get_resname()})

    counts, _ = classify_all(data)
    total = sum(counts.values())
    assert total > 20, "expected a reasonable number of residues"
    good = (counts['favored'] + counts['allowed']) / total
    assert good > 0.95, f"only {good:.1%} favored+allowed"


def test_render_produces_svg():
    _counts, annotated = classify_all([
        {'phi': -63, 'psi': -43, 'resname': 'ALA', 'residue': 'ALA1'},
        {'phi': 120, 'psi': 120, 'resname': 'ALA', 'residue': 'ALA2'},
    ])
    assert render_ramachandran(annotated).lstrip().startswith("<svg")


def test_endpoint_returns_plot_and_reachable_outliers(client, repo_root):
    import io
    import os
    path = os.path.join(repo_root, 'static', 'sample_protein.pdb')
    with open(path, 'rb') as fh:
        payload = fh.read()
    resp = client.post('/api/structure/ramachandran', data={
        'file': (io.BytesIO(payload), 'sample_protein.pdb'),
        'chain_id': 'A',
    }, content_type='multipart/form-data')
    data = resp.get_json()
    assert data['success'] is True
    assert data['plot'].lstrip().startswith("<svg")
    c = data['classification']
    assert c['favored'] + c['allowed'] + c['outliers'] == data['total_residues']
    assert 'outlier_residues' in data
