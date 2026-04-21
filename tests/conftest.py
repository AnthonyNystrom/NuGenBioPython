"""
Shared pytest fixtures.

Keep this file small — complex setup goes in individual test modules or
fixtures/inventory.py.
"""
import os
import sys

import pytest

# Make the app importable regardless of where pytest is invoked from
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


@pytest.fixture(scope="session")
def repo_root():
    return _REPO_ROOT


@pytest.fixture(scope="session")
def app():
    """Boot the Flask app once per test session."""
    from app import app as flask_app

    flask_app.config.update(TESTING=True)
    return flask_app


@pytest.fixture
def client(app):
    """Fresh test client per test. Cheap to construct."""
    return app.test_client()
