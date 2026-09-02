"""NCBI identity handling.

NCBI's usage policy asks every Entrez request to carry a real contact
address and offers a higher rate ceiling with an API key. The routes used to
default to the literal placeholder 'user@example.com', which identifies every
deployment as the same fake contact — the fast route to an IP block.
"""
import os

import pytest
from Bio import Entrez

from utils.entrez_helpers import configure_entrez, _looks_like_email


@pytest.fixture(autouse=True)
def _clean_entrez_env(monkeypatch):
    monkeypatch.delenv("ENTREZ_EMAIL", raising=False)
    monkeypatch.delenv("ENTREZ_API_KEY", raising=False)
    yield


@pytest.mark.parametrize("value", [
    "user@example.com", "your.email@example.com", "email@example.com",
    "", None, "not-an-email", "no@domain", "trailing@dot.",
])
def test_placeholders_and_junk_are_rejected(value):
    assert _looks_like_email(value) is False


@pytest.mark.parametrize("value", [
    "researcher@university.edu", "a.b@lab.example.org",
])
def test_real_addresses_are_accepted(value):
    assert _looks_like_email(value) is True


def test_caller_supplied_address_wins():
    assert configure_entrez("scientist@lab.edu") == "scientist@lab.edu"
    assert Entrez.email == "scientist@lab.edu"


def test_placeholder_falls_back_to_environment(monkeypatch):
    monkeypatch.setenv("ENTREZ_EMAIL", "ops@deployment.org")
    assert configure_entrez("user@example.com") == "ops@deployment.org"
    assert Entrez.email == "ops@deployment.org"


def test_missing_address_falls_back_to_environment(monkeypatch):
    monkeypatch.setenv("ENTREZ_EMAIL", "ops@deployment.org")
    assert configure_entrez("") == "ops@deployment.org"


def test_no_identity_anywhere_does_not_raise():
    """Must degrade to a warning, never break the feature outright."""
    assert configure_entrez("") is None


def test_api_key_is_applied_when_set(monkeypatch):
    monkeypatch.setenv("ENTREZ_API_KEY", "abc123")
    configure_entrez("scientist@lab.edu")
    assert Entrez.api_key == "abc123"


def test_tool_is_always_identified():
    configure_entrez("scientist@lab.edu")
    assert Entrez.tool == "NuGenBioPython"


def test_routes_no_longer_ship_a_placeholder_default(repo_root):
    """Guard against the placeholder creeping back into the route layer."""
    import pathlib
    for name in ("routes/database_routes.py", "routes/advanced_routes.py"):
        text = pathlib.Path(repo_root, name).read_text()
        assert "user@example.com" not in text, name
        assert "your.email@example.com" not in text, name
