"""
Every sample/example file referenced by a Load Example button must exist.

If one of these disappears (gitignored, renamed, deleted), Playwright tests
in phase 0.6 would fail in a confusing way; this test fails with a clear
'missing file' message up front.
"""
import os

import pytest

from tests.fixtures.inventory import STATIC_SAMPLES, STATIC_AUTO_GENERATED


@pytest.mark.parametrize("path", STATIC_SAMPLES)
def test_static_sample_exists(repo_root, path):
    full = os.path.join(repo_root, path)
    assert os.path.isfile(full), f"missing sample file: {path}"
    size = os.path.getsize(full)
    assert size > 0, f"sample file is empty: {path} ({size} bytes)"


@pytest.mark.parametrize("path", STATIC_AUTO_GENERATED)
def test_server_auto_generates_file(app, repo_root, path):
    """configure_app() creates files in uploads/ on first boot. Booting the
    Flask app (via the `app` fixture) should have produced them."""
    full = os.path.join(repo_root, path)
    assert os.path.isfile(full), f"server did not create {path}"


def test_sample_protein_pdb_is_served(client):
    """/structure/sample downloads static/sample_protein.pdb. This is the
    golden-path smoke for Structure's Load PDB Example flow."""
    resp = client.get("/api/structure/sample")
    assert resp.status_code == 200
    assert resp.data.startswith(b"HEADER") or b"ATOM" in resp.data[:500]
