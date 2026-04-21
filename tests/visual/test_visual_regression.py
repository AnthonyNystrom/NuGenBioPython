"""
Compare a fresh capture against tests/visual/baseline/ and fail on pixel
differences beyond the tolerance.

Skipped cleanly if playwright or pillow/pixelmatch aren't installed, and if
no baseline has been captured yet (Phase 0 ships without running
take_baseline.py — that's deliberate so CI doesn't need a browser).

To enable:
    pip install playwright pytest-playwright pillow
    playwright install chromium
    python tests/visual/take_baseline.py
"""
import os

import pytest


_VISUAL_DIR = os.path.dirname(__file__)
_BASELINE_DIR = os.path.join(_VISUAL_DIR, "baseline")


def _deps_ok():
    try:
        import playwright  # noqa: F401
        from PIL import Image  # noqa: F401
        return True
    except ImportError:
        return False


def _baseline_exists():
    return os.path.isdir(_BASELINE_DIR) and any(
        f.endswith(".png") for f in os.listdir(_BASELINE_DIR)
    )


@pytest.mark.skipif(not _deps_ok(), reason="playwright/pillow not installed")
@pytest.mark.skipif(not _baseline_exists(),
                    reason="no baseline captured; run tests/visual/take_baseline.py")
def test_visual_regression_placeholder():
    """Placeholder. Full implementation lands after Phase 0 baseline is captured.

    The design is:
    1. Capture fresh screenshots into tests/visual/current/
    2. For each baseline PNG, diff against current via pixelmatch
    3. Write diff PNGs to tests/visual/diff/
    4. Fail if any pair differs by more than TOLERANCE_FRACTION pixels
    """
    pytest.skip("visual regression runner deferred until baseline is captured")
