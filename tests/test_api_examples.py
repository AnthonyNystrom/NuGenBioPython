"""
Regression tests for every API endpoint a Load Example button would hit.

Goal: if a redesign (or accidental backend change) breaks the Load Example
flow for any tool, this file fails with a specific endpoint + payload
difference. External-service endpoints (NCBI, KEGG, BLAST remote) are
deliberately NOT exercised here — they're in the Playwright matrix where
they can be mocked or run on a separate schedule.
"""
import io
import os

import pytest

from tests.fixtures.inventory import API_EXAMPLES, API_UPLOAD_EXAMPLES


def _ids(entries, key):
    return [e[key].split("/")[-1] for e in entries]


@pytest.mark.parametrize(
    "entry",
    API_EXAMPLES,
    ids=_ids(API_EXAMPLES, "endpoint"),
)
def test_json_api_example(client, entry):
    method = entry.get("method", "POST")
    endpoint = entry["endpoint"]
    payload = entry.get("payload", {})

    if method == "GET":
        resp = client.get(endpoint)
    else:
        resp = client.post(endpoint, json=payload)

    assert resp.status_code == 200, (
        f"{method} {endpoint} → {resp.status_code}"
    )
    body = resp.get_json()
    assert body is not None, f"{endpoint} returned non-JSON body"
    assert body.get("success") is True, (
        f"{endpoint} success!=True: error={body.get('error')!r}"
    )
    for key in entry.get("expect_keys", []):
        assert key in body, (
            f"{endpoint} response missing '{key}'; got keys={list(body.keys())}"
        )


@pytest.mark.parametrize(
    "entry",
    API_UPLOAD_EXAMPLES,
    ids=_ids(API_UPLOAD_EXAMPLES, "endpoint"),
)
def test_upload_api_example(client, repo_root, entry):
    endpoint = entry["endpoint"]
    if "fixture_file" in entry:
        full = os.path.join(repo_root, entry["fixture_file"])
        assert os.path.isfile(full), f"fixture missing: {entry['fixture_file']}"
        with open(full, "rb") as fh:
            data = fh.read()
    else:
        data = entry["fixture_bytes"]

    form = dict(entry.get("form", {}))
    form["file"] = (io.BytesIO(data), entry["fixture_name"])

    resp = client.post(endpoint, data=form, content_type="multipart/form-data")
    assert resp.status_code == 200, f"POST {endpoint} → {resp.status_code}"
    body = resp.get_json()
    assert body is not None, f"{endpoint} returned non-JSON"
    assert body.get("success") is True, (
        f"{endpoint} success!=True: error={body.get('error')!r}"
    )
    for key in entry.get("expect_keys", []):
        assert key in body, (
            f"{endpoint} response missing '{key}'; got keys={list(body.keys())}"
        )


def test_missing_file_returns_json_error(client):
    """Invariant: file-upload endpoints must return JSON {success:false, error:...}
    when the file is missing, not 500. This prevents regressions of the
    upload_helpers cleanup contract."""
    resp = client.post(
        "/api/seqio/parse", data={"format": "fasta"}, content_type="multipart/form-data"
    )
    assert resp.status_code == 200
    body = resp.get_json()
    assert body is not None
    assert body.get("success") is False
    assert body.get("error")


def test_request_json_null_safe(client):
    """JSON endpoints must handle wrong Content-Type gracefully."""
    resp = client.post("/api/seqio/convert", data="not-json", content_type="text/plain")
    body = resp.get_json()
    assert body is not None
    assert body.get("success") is False
