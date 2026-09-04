"""
NuGenBioPython - BioPython Web Application
Main application entry point
"""
# Set matplotlib backend before any other imports to prevent NSWindow threading issues
import matplotlib
matplotlib.use('Agg')

import logging
import os
import secrets

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

logging.basicConfig(
    level=os.environ.get('LOG_LEVEL', 'INFO').upper(),
    format='%(asctime)s %(levelname)s %(name)s: %(message)s',
)

from flask import Flask, g
from utils.config import configure_app
from routes import register_blueprints

app = Flask(__name__)
configure_app(app)
register_blueprints(app)


# Cache-busting for static assets: append ?v=<file-mtime> to every
# url_for('static', filename=...) call so browsers fetch fresh CSS/JS
# whenever the source file changes, even if max-age says the old one
# is cacheable.
@app.url_defaults
def _static_cache_bust(endpoint, values):
    if endpoint == 'static' and 'filename' in values and 'v' not in values:
        static_folder = app.static_folder
        if static_folder:
            path = os.path.join(static_folder, values['filename'])
            try:
                values['v'] = int(os.path.getmtime(path))
            except OSError:
                pass


# Content-Security-Policy.
#
# script-src no longer carries 'unsafe-inline'. Every inline <script> block
# carries a per-request nonce (see _csp_nonce below), and the ~54 inline
# onclick= handlers that previously made a nonce impossible — inline event
# handlers can never be nonce-covered — have been migrated to the data-action
# dispatcher in static/js/utils.js.
#
# style-src no longer carries 'unsafe-inline' either. All 139 inline
# style="..." attributes became utility classes (static values) or data-css
# attributes applied through the CSSOM (dynamic values) — el.style.cssText is
# not covered by style-src, unlike a style attribute in markup. The nonce
# covers the inline <style> block in base.html.
#
# Connect-src allowlists the external biology APIs the app talks to.
def _build_csp(nonce):
    return (
        "default-src 'self'; "
        f"script-src 'self' 'nonce-{nonce}' "
            "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com https://3Dmol.org; "
        f"style-src 'self' 'nonce-{nonce}' "
            "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com https://fonts.googleapis.com; "
        "font-src 'self' https://fonts.gstatic.com https://cdnjs.cloudflare.com data:; "
        "img-src 'self' data: https:; "
        "connect-src 'self' "
            "https://eutils.ncbi.nlm.nih.gov https://www.ncbi.nlm.nih.gov "
            "https://rest.kegg.jp https://rest.uniprot.org https://www.uniprot.org "
            # RCSB PDB — /structure fetches .pdb files directly from the browser
            # (see structure.js fetchFromPDB).
            "https://files.rcsb.org "
            # CDN hosts — needed so DevTools can fetch .js.map / .css.map
            # sourcemaps for bootstrap.bundle.min.js and Font Awesome.
            "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com; "
        # 3Dmol.js spawns Web Workers from blob: URLs for surface computation.
        "worker-src 'self' blob:; "
        "object-src 'none'; "
        "base-uri 'self'; "
        "frame-ancestors 'none'"
    )


@app.before_request
def _csp_nonce():
    """Mint a fresh nonce per request for inline <script> blocks."""
    g.csp_nonce = secrets.token_urlsafe(16)


@app.context_processor
def _inject_csp_nonce():
    """Expose the nonce to templates as {{ csp_nonce }}."""
    return {'csp_nonce': getattr(g, 'csp_nonce', '')}


@app.after_request
def _apply_security_headers(response):
    response.headers.setdefault('X-Content-Type-Options', 'nosniff')
    response.headers.setdefault('X-Frame-Options', 'DENY')
    response.headers.setdefault('Referrer-Policy', 'strict-origin-when-cross-origin')
    response.headers.setdefault(
        'Content-Security-Policy', _build_csp(getattr(g, 'csp_nonce', '')))
    response.headers.setdefault('Permissions-Policy', 'camera=(), microphone=(), geolocation=()')
    return response


# Figures no longer need a global sweep. utils.plot_helpers.subplots()
# constructs matplotlib Figure objects directly instead of going through
# pyplot, so nothing is registered in pyplot's process-global figure dict and
# each figure is collected with the request that made it.
#
# The teardown_request hook that used to live here called plt.close('all').
# That reclaimed leaked figures, but it was also the third reason this app
# could not run threaded workers: 'all' meant all, including figures still
# being drawn by other in-flight requests.


if __name__ == '__main__':
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'
    port = int(os.environ.get('PORT', 9000))
    app.run(debug=debug_mode, port=port)
