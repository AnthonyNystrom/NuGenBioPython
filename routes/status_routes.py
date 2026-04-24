"""External-service health probes.

Each endpoint returns {ok: bool, latency_ms: int, error?: str} after
attempting a lightweight HEAD/GET against the real service. Responses
are cached in-memory for 60s so the probe itself doesn't become a
per-page rate-limit concern against KEGG/NCBI.
"""
import time
from functools import lru_cache

from flask import Blueprint, jsonify

try:
    import requests
except ImportError:  # pragma: no cover
    requests = None  # type: ignore

status_bp = Blueprint('status', __name__, url_prefix='/api/status')

_CACHE: dict[str, tuple[float, dict]] = {}
_CACHE_TTL = 60  # seconds


def _cached(key: str, fn):
    now = time.time()
    if key in _CACHE:
        ts, payload = _CACHE[key]
        if now - ts < _CACHE_TTL:
            return payload
    payload = fn()
    _CACHE[key] = (now, payload)
    return payload


def _probe(url: str, timeout: float = 3.0, method: str = 'HEAD') -> dict:
    if requests is None:
        return {'ok': False, 'error': 'requests library not available'}
    try:
        start = time.time()
        resp = requests.request(method, url, timeout=timeout, allow_redirects=True)
        latency_ms = int((time.time() - start) * 1000)
        return {'ok': resp.status_code < 500, 'latency_ms': latency_ms, 'status': resp.status_code}
    except requests.Timeout:
        return {'ok': False, 'error': 'timeout'}
    except requests.ConnectionError:
        return {'ok': False, 'error': 'connection_error'}
    except Exception as e:  # pragma: no cover
        return {'ok': False, 'error': str(e)[:80]}


@status_bp.route('/kegg')
def kegg_status():
    """KEGG REST endpoint."""
    return jsonify(_cached('kegg', lambda: _probe('https://rest.kegg.jp/info/pathway', method='GET')))


@status_bp.route('/ncbi')
def ncbi_status():
    """NCBI Entrez eutils."""
    return jsonify(_cached('ncbi', lambda: _probe('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi', method='GET')))


@status_bp.route('/ncbi-blast')
def ncbi_blast_status():
    """NCBI BLAST — same host as Entrez, but we probe separately in case
    BLAST-specific endpoints degrade."""
    return jsonify(_cached('ncbi-blast', lambda: _probe('https://www.ncbi.nlm.nih.gov/blast/', method='HEAD')))


@status_bp.route('/uniprot')
def uniprot_status():
    """UniProt REST."""
    return jsonify(_cached('uniprot', lambda: _probe('https://rest.uniprot.org/uniprotkb/P68871', method='HEAD')))


def register(app):
    app.register_blueprint(status_bp)
