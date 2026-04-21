"""
NuGenBioPython - BioPython Web Application
Main application entry point
"""
# Set matplotlib backend before any other imports to prevent NSWindow threading issues
import matplotlib
matplotlib.use('Agg')

import logging
import os

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

logging.basicConfig(
    level=os.environ.get('LOG_LEVEL', 'INFO').upper(),
    format='%(asctime)s %(levelname)s %(name)s: %(message)s',
)

from flask import Flask
from utils.config import configure_app
from routes import register_blueprints

app = Flask(__name__)
configure_app(app)
register_blueprints(app)


# Content-Security-Policy. 'unsafe-inline' on script-src is retained for now
# because ~73 inline onclick handlers with arguments + inline <script> blocks
# remain (simple fn() cases are migrated to data-action via utils.js). A
# future pass can externalize those and drop 'unsafe-inline' entirely.
# Connect-src allowlists the external biology APIs the app talks to.
_CSP = (
    "default-src 'self'; "
    "script-src 'self' 'unsafe-inline' "
        "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com https://3Dmol.org; "
    "style-src 'self' 'unsafe-inline' "
        "https://cdn.jsdelivr.net https://cdnjs.cloudflare.com https://fonts.googleapis.com; "
    "font-src 'self' https://fonts.gstatic.com https://cdnjs.cloudflare.com data:; "
    "img-src 'self' data: https:; "
    "connect-src 'self' "
        "https://eutils.ncbi.nlm.nih.gov https://www.ncbi.nlm.nih.gov "
        "https://rest.kegg.jp https://rest.uniprot.org https://www.uniprot.org; "
    # 3Dmol.js spawns Web Workers from blob: URLs for surface computation.
    "worker-src 'self' blob:; "
    "object-src 'none'; "
    "base-uri 'self'; "
    "frame-ancestors 'none'"
)


@app.after_request
def _apply_security_headers(response):
    response.headers.setdefault('X-Content-Type-Options', 'nosniff')
    response.headers.setdefault('X-Frame-Options', 'DENY')
    response.headers.setdefault('Referrer-Policy', 'strict-origin-when-cross-origin')
    response.headers.setdefault('Content-Security-Policy', _CSP)
    response.headers.setdefault('Permissions-Policy', 'camera=(), microphone=(), geolocation=()')
    return response

if __name__ == '__main__':
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'
    port = int(os.environ.get('PORT', 9000))
    app.run(debug=debug_mode, port=port)
