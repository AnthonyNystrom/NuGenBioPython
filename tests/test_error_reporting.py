"""Route failures must be observable server-side without changing the client
contract.

197 route handlers returned `{'success': False, 'error': str(e)}` and only two
of them logged anything, so almost every failure mode left no server-side
trace at all. The catch-all handlers now go through
utils.request_helpers.error_response, which logs the traceback and returns the
*same* payload — the migration was deliberately behaviour-preserving.
"""
import logging
import pathlib
import re

import pytest


def test_error_response_preserves_the_legacy_payload(app):
    """With no explicit user message, the body must match the old shape."""
    from utils.request_helpers import error_response
    with app.test_request_context():
        resp = error_response(ValueError("boom"))
        body = resp.get_json()
    assert body == {"success": False, "error": "boom"}


def test_error_response_logs_the_traceback(app, caplog):
    from utils.request_helpers import error_response
    with caplog.at_level(logging.ERROR, logger="utils.request_helpers"):
        with app.test_request_context():
            error_response(ValueError("boom"), context="demo.handler")
    assert caplog.records, "the failure must reach the log"
    assert "demo.handler" in caplog.text
    assert "boom" in caplog.text


def test_error_response_can_hide_internals(app):
    """An explicit user message must not leak the exception text."""
    from utils.request_helpers import error_response
    with app.test_request_context():
        body = error_response(ValueError("/srv/secret/path.py line 42"),
                              user_msg="Could not parse that file").get_json()
    assert body["error"] == "Could not parse that file"
    assert "secret" not in body["error"]


def test_a_failing_route_logs_and_still_answers(client, caplog):
    """End-to-end: a real parse failure logs server-side, client sees JSON."""
    with caplog.at_level(logging.ERROR):
        resp = client.post("/api/phylo/parse",
                           data={"tree_string": "((((broken newick",
                                 "format": "newick"})
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["success"] is False
    assert body["error"]
    assert "phylo_routes.parse_tree" in caplog.text


def test_catch_all_handlers_are_not_left_silent(repo_root):
    """Guard: a new bare `except Exception` catch-all must not slip back in."""
    offenders = []
    for path in sorted(pathlib.Path(repo_root, "routes").rglob("*.py")):
        lines = path.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.strip() != "except Exception as e:":
                continue
            nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
            if nxt == "return jsonify({'success': False, 'error': str(e)})":
                offenders.append(f"{path.name}:{i + 2}")
    assert not offenders, (
        "these catch-alls return str(e) without logging; "
        f"use error_response(e, context=...) instead: {offenders}"
    )


# --------------------------------------------------------------------------
# Information disclosure
# --------------------------------------------------------------------------

def test_parse_errors_reach_the_user(app):
    """Messages written for a biologist are the most useful thing we can say."""
    from Bio.Phylo.NewickIO import NewickError
    from utils.request_helpers import safe_error_message
    assert safe_error_message(
        NewickError("Mismatch, 4 open vs 0 close parentheses.")
    ) == "Mismatch, 4 open vs 0 close parentheses."
    assert safe_error_message(
        ValueError("Sequences must all be the same length")
    ) == "Sequences must all be the same length"


def test_upload_validation_messages_reach_the_user():
    from utils.request_helpers import safe_error_message
    from utils.upload_helpers import UploadError
    assert safe_error_message(UploadError("No file provided")) == "No file provided"


@pytest.mark.parametrize("exc", [
    OSError("[Errno 2] No such file or directory: '/srv/app/uploads/a1_x.pdb'"),
    KeyError("_internal_cache_key"),
    AttributeError("'NoneType' object has no attribute 'seq'"),
    IndexError("list index out of range"),
])
def test_internal_errors_do_not_reach_the_user(exc):
    """Upload paths and object structure must not travel to the browser."""
    from utils.request_helpers import safe_error_message
    message = safe_error_message(exc)
    assert "uploads" not in message
    assert "_internal_cache_key" not in message
    assert "NoneType" not in message
    assert "logged" in message


def test_paths_are_scrubbed_even_from_safe_messages():
    from utils.request_helpers import safe_error_message
    message = safe_error_message(
        ValueError("could not open /srv/app/uploads/secret_file.pdb"))
    assert "/srv/app/uploads" not in message
    assert "<path>" in message


def test_explicit_user_message_still_wins(app):
    from utils.request_helpers import error_response
    with app.test_request_context():
        body = error_response(OSError("/srv/secret"),
                              user_msg="Could not read that file").get_json()
    assert body["error"] == "Could not read that file"


def test_failing_route_no_longer_leaks_internals(client):
    """End-to-end: a route that raises OSError must not echo the path."""
    resp = client.post("/api/structure/ramachandran", data={},
                       content_type="multipart/form-data")
    body = resp.get_json()
    assert body["success"] is False
    assert "uploads" not in body["error"]


def test_no_catch_all_returns_a_raw_exception_string(repo_root):
    """Guard: every catch-all must sanitize before answering the client.

    Covers all the shapes the migration had to handle — a bare
    `{'error': str(e)}`, an f-string wrapping it with context, an
    intermediate `error_msg = str(e)`, and per-item annotations inside a
    larger result dict.
    """
    offenders = []
    for path in sorted(pathlib.Path(repo_root, "routes").rglob("*.py")):
        lines = path.read_text().splitlines()
        for i, line in enumerate(lines):
            if "str(e)" not in line or "safe_error_message" in line:
                continue
            prev = lines[i - 1].strip() if i else ""
            if prev.startswith("except Exception"):
                offenders.append(f"{path.name}:{i + 1}  {line.strip()[:70]}")
    assert not offenders, (
        "raw exception text reaches the client; wrap it in "
        f"safe_error_message(e): {offenders}")
