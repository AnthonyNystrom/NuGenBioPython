"""
Bridge test: runs tests/formatters/run.js under Node with jsdom.

Skipped cleanly if Node isn't available. Install jsdom once:
    npm install --prefix tests/formatters jsdom
"""
import os
import shutil
import subprocess

import pytest


_FORMATTERS_DIR = os.path.join(os.path.dirname(__file__), "formatters")


def _node_available():
    return shutil.which("node") is not None


def _jsdom_installed():
    return os.path.isdir(os.path.join(_FORMATTERS_DIR, "node_modules", "jsdom"))


@pytest.mark.skipif(not _node_available(), reason="node is not installed")
@pytest.mark.skipif(
    not _jsdom_installed(),
    reason="jsdom not installed. Run: npm install --prefix tests/formatters jsdom",
)
def test_formatters_js():
    """Run the Node-side formatter assertions and propagate exit status."""
    result = subprocess.run(
        ["node", os.path.join(_FORMATTERS_DIR, "run.js")],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            "Formatter tests failed:\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
