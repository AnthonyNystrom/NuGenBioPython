"""Sequence length limits.

Bio.Align.PairwiseAligner is O(n*n) in both time and memory. Measured on
random DNA: 10 kb -> 1.0 s / 243 MB, 20 kb -> 4.0 s / 861 MB, 40 kb ->
16.1 s / 4.0 GB, scaling 4x per doubling. With MAX_CONTENT_LENGTH at 16 MB
and no length check, one request could carry ~8,000,000 bp: roughly 179
hours of CPU, and an out-of-memory kill long before that.

Memory is the sharper edge. Extrapolating the curve, an 8 GB host is
exhausted by a ~57 kb sequence — one assembly contig, an entirely ordinary
paste. The app runs one request per process, so that OOM takes the whole
service down. These tests pin the guard that prevents it.
"""
import logging
import time

import pytest

from utils.request_helpers import (
    MAX_ALIGNMENT_LEN, MAX_SEQUENCE_LEN, SequenceTooLong,
    check_sequence_length, check_alignment_inputs,
)


# --------------------------------------------------------------------------
# The guard itself
# --------------------------------------------------------------------------

def test_quadratic_limit_is_far_below_the_oom_threshold():
    """~57 kb OOMs an 8 GB box; the limit must leave real headroom."""
    assert MAX_ALIGNMENT_LEN <= 20_000


def test_linear_limit_is_generous():
    """Linear work is cheap (5 Mb in 0.18 s), so it should not be crippled."""
    assert MAX_SEQUENCE_LEN >= 1_000_000


def test_sequence_at_the_limit_is_accepted():
    seq = "A" * MAX_ALIGNMENT_LEN
    assert check_sequence_length(seq, quadratic=True) == seq


def test_sequence_one_over_the_limit_is_rejected():
    with pytest.raises(SequenceTooLong):
        check_sequence_length("A" * (MAX_ALIGNMENT_LEN + 1), quadratic=True)


def test_message_carries_exact_figures(): 
    """Rounding alone read as 'is 10 kb; the limit is 10 kb' at the boundary."""
    with pytest.raises(SequenceTooLong) as exc:
        check_sequence_length("A" * (MAX_ALIGNMENT_LEN + 1), quadratic=True)
    message = str(exc.value)
    assert f"{MAX_ALIGNMENT_LEN + 1:,}" in message
    assert f"{MAX_ALIGNMENT_LEN:,}" in message


def test_none_passes_through():
    assert check_sequence_length(None) is None


def test_alignment_inputs_names_the_offending_sequence():
    with pytest.raises(SequenceTooLong) as exc:
        check_alignment_inputs("A" * 10, "C" * (MAX_ALIGNMENT_LEN + 1))
    assert "sequence 2" in str(exc.value)


# --------------------------------------------------------------------------
# Endpoints
# --------------------------------------------------------------------------

QUADRATIC_ENDPOINTS = [
    ("/api/alignment/pairwise", {"sequence1": None, "sequence2": None}),
    ("/api/alignment/coordinates", {"sequence1": None, "sequence2": None}),
    ("/api/alignment/detailed_stats", {"sequence1": None, "sequence2": None}),
    ("/api/alignment/codon_aware", {"sequence1": None, "sequence2": None}),
    ("/api/alignment/identity_matrix", {"sequences": None}),
    ("/api/alignment/multiple", {"sequences": None}),
]


@pytest.mark.parametrize("path,template", QUADRATIC_ENDPOINTS)
def test_oversized_input_is_refused_immediately(client, path, template):
    """The DoS: a huge sequence must be rejected, not computed."""
    big = "A" * 200_000
    payload = {}
    for key, _ in template.items():
        payload[key] = [big, big] if key == "sequences" else big

    started = time.time()
    data = client.post(path, json=payload).get_json()
    elapsed = time.time() - started

    assert data["success"] is False, f"{path} accepted a 200 kb sequence"
    assert "limit" in data["error"].lower()
    # Unguarded, 200 kb would take minutes and gigabytes.
    assert elapsed < 5, f"{path} spent {elapsed:.1f}s before refusing"


def test_the_full_content_length_ceiling_is_refused_fast(client):
    """MAX_CONTENT_LENGTH allows ~8 Mb per sequence: ~179 hours unguarded."""
    huge = "A" * 8_000_000
    started = time.time()
    data = client.post("/api/alignment/pairwise",
                       json={"sequence1": huge, "sequence2": huge}).get_json()
    elapsed = time.time() - started
    assert data["success"] is False
    assert elapsed < 10, f"took {elapsed:.1f}s to refuse an 8 Mb sequence"


def test_ordinary_sequences_still_work(client):
    """A normal gene-sized input must be unaffected."""
    seq = "ATGCGTACGT" * 100   # 1 kb
    data = client.post("/api/alignment/pairwise",
                       json={"sequence1": seq, "sequence2": seq}).get_json()
    assert data["success"] is True
    assert data["statistics"]["identity_percent"] == 100.0


def test_a_sequence_at_the_limit_still_aligns(client):
    seq = "ACGT" * (MAX_ALIGNMENT_LEN // 4)
    data = client.post("/api/alignment/pairwise",
                       json={"sequence1": seq, "sequence2": seq}).get_json()
    assert data["success"] is True


# --------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------

def test_oversized_input_does_not_log_a_traceback(client, caplog):
    """A big paste is ordinary use, not a server fault."""
    with caplog.at_level(logging.DEBUG):
        client.post("/api/alignment/pairwise",
                    json={"sequence1": "A" * 50_000, "sequence2": "C" * 50_000})
    assert "Traceback" not in caplog.text
    assert any(r.levelno <= logging.INFO for r in caplog.records)
