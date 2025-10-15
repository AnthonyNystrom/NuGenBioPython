"""
NuGenBioPython - BioPython Web Application
Main application entry point
"""
# Set matplotlib backend before any other imports to prevent NSWindow threading issues
import matplotlib
matplotlib.use('Agg')

from flask import Flask
from utils.config import configure_app
from routes import register_blueprints

# Create Flask application
app = Flask(__name__)

# Configure application
configure_app(app)

# Register all blueprints
register_blueprints(app)

if __name__ == '__main__':
    import os
    debug_mode = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'
    port = int(os.environ.get('PORT', 9000))
    app.run(debug=debug_mode, port=port)
