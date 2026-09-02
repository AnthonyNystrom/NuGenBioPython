"""SearchIO threshold filters must fail closed.

The filters previously wrapped each comparison in `except Exception: pass`,
leaving `passes_filter` True when a value could not be read as a number. A
hit that was never actually checked against the threshold was therefore
returned as if it had passed — a search filtered at "e-value < 1e-5" could
still contain hits far above it.
"""
import logging

import pytest

from utils.searchio_helpers import _within


class _Hit:
    id = "hit_1"


@pytest.mark.parametrize("value,threshold,expected", [
    (1e-9, 1e-5, True),    # comfortably below the ceiling
    (1e-5, 1e-5, True),    # exactly at the ceiling
    (1e-2, 1e-5, False),   # above the ceiling
])
def test_max_direction(value, threshold, expected):
    assert _within(value, threshold, 'max') is expected


@pytest.mark.parametrize("value,threshold,expected", [
    (200.0, 50.0, True),
    (50.0, 50.0, True),
    (10.0, 50.0, False),
])
def test_min_direction(value, threshold, expected):
    assert _within(value, threshold, 'min') is expected


@pytest.mark.parametrize("bad", [None, "n/a", "", object(), float("nan")])
def test_unreadable_values_are_excluded_not_admitted(bad):
    """The regression: an unevaluable hit must be dropped, never kept."""
    if isinstance(bad, float):
        # NaN compares False either way; it must still not pass a max filter.
        assert _within(bad, 1e-5, 'max') is False
        return
    assert _within(bad, 1e-5, 'max') is False
    assert _within(bad, 50.0, 'min') is False


def test_unreadable_threshold_is_also_excluded():
    assert _within(1e-9, "not-a-number", 'max') is False


def test_exclusion_is_logged(caplog):
    """An operator must be able to see why a hit disappeared."""
    with caplog.at_level(logging.WARNING, logger="utils.searchio_helpers"):
        _within(None, 1e-5, 'max', _Hit(), 'evalue')
    assert caplog.records, "expected a warning when a filter cannot be applied"
    assert "evalue" in caplog.text
    assert "hit_1" in caplog.text
