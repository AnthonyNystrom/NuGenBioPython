"""Regression tests for pairwise-alignment statistics.

These lock down the bug where alignment statistics were derived by scraping
``str(alignment)``. BioPython wraps that text at 60 columns and omits the
trailing coordinate on wrapped lines, so ``line.split()[-2]`` returned the
start coordinate instead of the sequence — every alignment longer than one
block reported 100% identity with alignment_length 1.
"""
import pytest


# 120 bp: long enough that str(alignment) wraps into multiple blocks.
IDENTICAL = "ATGCGTACGT" * 12
# Same length, second half fully divergent -> 50% identity.
HALF_DIVERGENT = "ATGCGTACGT" * 6 + "TTTTTTTTTT" * 6


def _post(client, url, payload):
    resp = client.post(url, json=payload)
    assert resp.status_code == 200, f"{url} -> {resp.status_code}"
    data = resp.get_json()
    assert data.get("success") is True, data.get("error")
    return data


def test_identical_long_sequences_report_full_length(client):
    """A 120 bp self-alignment must report 100% identity over 120 columns."""
    data = _post(client, "/api/alignment/pairwise",
                 {"sequence1": IDENTICAL, "sequence2": IDENTICAL})
    stats = data["statistics"]
    assert stats["alignment_length"] == 120
    assert stats["matches"] == 120
    assert stats["identity_percent"] == 100.0
    assert stats["gaps"] == 0


def test_divergent_long_sequences_are_not_reported_as_identical(client):
    """The headline bug: 50%-divergent sequences reported 100% identity."""
    data = _post(client, "/api/alignment/pairwise",
                 {"sequence1": IDENTICAL, "sequence2": HALF_DIVERGENT})
    stats = data["statistics"]
    assert stats["identity_percent"] < 100.0
    assert stats["alignment_length"] > 1
    # Aligner may introduce gaps; identity should land near 50%, never 100%.
    assert 30.0 <= stats["identity_percent"] <= 70.0


def test_short_sequences_still_correct(client):
    """The unwrapped (<=60 column) path must keep working."""
    seq = "ATGCGTACGT" * 4  # 40 bp, single block
    data = _post(client, "/api/alignment/pairwise",
                 {"sequence1": seq, "sequence2": seq})
    stats = data["statistics"]
    assert stats["alignment_length"] == 40
    assert stats["identity_percent"] == 100.0


def test_gaps_are_counted_across_the_whole_alignment(client):
    """A long deletion must show up in the gap statistics."""
    long_seq = "ATGCGTACGT" * 12   # 120 bp
    short_seq = "ATGCGTACGT" * 8   # 80 bp -> ~40 gap columns
    data = _post(client, "/api/alignment/pairwise",
                 {"sequence1": long_seq, "sequence2": short_seq})
    stats = data["statistics"]
    assert stats["alignment_length"] >= 120
    assert stats["gaps"] > 0, "gaps must be detected beyond the first block"


def test_identity_matrix_diverges_off_diagonal(client):
    """Identity matrix must not read 100% for divergent long sequences."""
    data = _post(client, "/api/alignment/identity_matrix",
                 {"sequences": [IDENTICAL, HALF_DIVERGENT]})
    matrix = data["matrix"]
    assert matrix[0][0] == 100.0
    assert matrix[0][1] < 100.0, "off-diagonal identity must reflect divergence"
    assert matrix[0][1] == matrix[1][0], "matrix must stay symmetric"


def test_detailed_stats_cover_whole_alignment(client):
    """detailed_stats must count matches beyond the first 60-column block."""
    data = _post(client, "/api/alignment/detailed_stats",
                 {"sequence1": IDENTICAL, "sequence2": IDENTICAL})
    stats = data["statistics"]
    assert stats["matches"] == 120, stats


def test_codon_aware_counts_all_codons(client):
    """codon_aware must score all 40 codons of a 120 bp sequence."""
    data = _post(client, "/api/alignment/codon_aware",
                 {"sequence1": IDENTICAL, "sequence2": IDENTICAL})
    cs = data["codon_statistics"]
    assert cs["codon_matches"] == 40, cs
    assert cs["codon_identity_percent"] == 100.0


def test_local_alignment_coordinates_are_the_aligned_region(client):
    """Local alignment must report the matched extent, not the full sequence."""
    target = "TTTTTTTTTT" + "ATGCGTACGT" * 8 + "GGGGGGGGGG"  # motif at 10..90
    query = "ATGCGTACGT" * 8
    resp = client.post("/api/alignment/coordinates",
                       json={"sequence1": target, "sequence2": query,
                             "mode": "local"})
    if resp.status_code == 404:
        pytest.skip("coordinates endpoint not routed at this path")
    data = resp.get_json()
    assert data.get("success") is True, data.get("error")
    c = data["coordinates"]
    assert (c["target_start"], c["target_end"]) == (10, 90)
    assert (c["query_start"], c["query_end"]) == (0, 80)
