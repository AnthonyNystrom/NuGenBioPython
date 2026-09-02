"""GC analysis reporting.

Bio.SeqUtils.GC123 returns a 4-tuple (overall, pos1, pos2, pos3) whose values
are already percentages. The template read indices 0-2 and multiplied each by
100, so it labelled the OVERALL GC as "Pos 1", shifted the other two, dropped
codon position 3 entirely, and printed values like 4984.67%.
"""
import pathlib
import re

import pytest


def test_gc123_is_a_four_tuple_of_percentages(client):
    """Pin the contract the template depends on."""
    seq = "ATGCGTACGT" * 60
    data = client.post("/api/sequence/gc_analysis",
                       json={"sequence": seq, "window": 100}).get_json()
    assert data["success"] is True
    gc123 = data["gc123"]
    assert len(gc123) == 4, "GC123 returns (overall, pos1, pos2, pos3)"
    for value in gc123:
        assert 0 <= value <= 100, f"{value} is not a percentage"


def test_overall_gc123_matches_reported_gc_content(client):
    """gc123[0] is the overall figure — this is what the off-by-one broke."""
    seq = "ATGCGTACGTTAGCCATGGATCC" * 20
    data = client.post("/api/sequence/gc_analysis",
                       json={"sequence": seq, "window": 100}).get_json()
    assert data["gc123"][0] == pytest.approx(data["gc_content"], abs=0.5)


def test_template_does_not_rescale_gc123(repo_root):
    """Guard: the * 100 must not come back."""
    html = pathlib.Path(repo_root, "templates", "sequence.html").read_text()
    offenders = re.findall(r"gc123\[\d\]\s*\*\s*100", html)
    assert not offenders, f"GC123 values are already percentages: {offenders}"


def test_template_reads_all_three_codon_positions(repo_root):
    """Positions live at indices 1, 2, 3 — index 0 is the overall value."""
    html = pathlib.Path(repo_root, "templates", "sequence.html").read_text()
    for index in (1, 2, 3):
        assert f"gc123[{index}]" in html, \
            f"codon position at index {index} is not displayed"
