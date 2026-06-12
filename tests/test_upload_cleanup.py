"""Regression tests for the upload-cleanup contract.

The searchio / swissprot / unigene routes were migrated to
utils.upload_helpers.saved_upload(), which guarantees the temp file is
removed even when the parser raises on malformed input. Before the fix the
os.remove() sat in the normal flow, so a parse error leaked the file into
UPLOAD_FOLDER (disk-fill over time). These tests pin that behavior:

  1. a malformed upload returns JSON {success: false}, never a 500, and
  2. it leaves no new file behind in UPLOAD_FOLDER.
"""
import io
import os

import pytest


# (endpoint, form fields) for each migrated parse route.
_PARSE_ENDPOINTS = [
    ("/api/searchio/parse", {"format": "blast-xml"}),
    ("/api/swissprot/parse", {"max_records": "5"}),
    ("/api/unigene/parse", {}),
]

_GARBAGE = b"this is not a valid bioinformatics file <<< >>> \x00\x01 not parseable"


def _upload_dir(app):
    return app.config["UPLOAD_FOLDER"]


def _snapshot(folder):
    try:
        return set(os.listdir(folder))
    except FileNotFoundError:
        return set()


@pytest.mark.parametrize("endpoint,form", _PARSE_ENDPOINTS)
def test_malformed_upload_is_graceful_and_leaves_no_file(client, app, endpoint, form):
    folder = _upload_dir(app)
    before = _snapshot(folder)

    payload = dict(form)
    payload["file"] = (io.BytesIO(_GARBAGE), "garbage.bin")
    resp = client.post(endpoint, data=payload, content_type="multipart/form-data")

    # Graceful JSON error, not a 500.
    assert resp.status_code == 200, f"{endpoint} → {resp.status_code}"
    body = resp.get_json()
    assert body is not None, f"{endpoint} returned non-JSON"
    assert body.get("success") is False, f"{endpoint} should fail on garbage input"

    # The temp file must have been cleaned up despite the parse error.
    after = _snapshot(folder)
    assert after == before, (
        f"{endpoint} leaked upload file(s): {after - before}"
    )


@pytest.mark.parametrize("endpoint,form", _PARSE_ENDPOINTS)
def test_missing_file_returns_json_error(client, endpoint, form):
    """saved_upload() raises UploadError when no file is provided; the route
    must surface that as JSON {success: false}, not a 500."""
    resp = client.post(endpoint, data=dict(form), content_type="multipart/form-data")
    assert resp.status_code == 200
    body = resp.get_json()
    assert body is not None
    assert body.get("success") is False
    assert body.get("error")
