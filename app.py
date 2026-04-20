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


@app.after_request
def _apply_security_headers(response):
    response.headers.setdefault('X-Content-Type-Options', 'nosniff')
    response.headers.setdefault('X-Frame-Options', 'DENY')
    response.headers.setdefault('Referrer-Policy', 'strict-origin-when-cross-origin')
    return response

if __name__ == '__main__':
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'
    port = int(os.environ.get('PORT', 9000))
    app.run(debug=debug_mode, port=port)
