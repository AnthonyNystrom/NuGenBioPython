"""
Flask application configuration
"""
import os
import logging
import secrets


log = logging.getLogger(__name__)


def _resolve_secret_key():
    key = os.environ.get('SECRET_KEY')
    if key:
        return key
    log.warning(
        "SECRET_KEY environment variable is not set. "
        "Generating an ephemeral random key for this process. "
        "Sessions will not survive restarts. Set SECRET_KEY in production."
    )
    return secrets.token_hex(32)


def configure_app(app):
    """Configure Flask application with settings"""
    app.config['SECRET_KEY'] = _resolve_secret_key()
    app.config['UPLOAD_FOLDER'] = os.environ.get('UPLOAD_FOLDER', 'uploads')
    app.config['MAX_CONTENT_LENGTH'] = int(os.environ.get('MAX_CONTENT_LENGTH', 16 * 1024 * 1024))
    # D24: cache static assets for 1 hour in dev, 1 day in prod. Browsers
    # still revalidate via Last-Modified, so edits still appear on reload.
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = int(
        os.environ.get('STATIC_MAX_AGE', '86400' if os.environ.get('FLASK_ENV') == 'production' else '3600')
    )

    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

    example_file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'example_popgen.gen')
    if not os.path.exists(example_file_path):
        example_content = """Example Population Genetics Study - Three Populations
Locus1
Locus2
Locus3
POP
Ind1_1, 0101 0202 0101
Ind1_2, 0102 0201 0102
Ind1_3, 0101 0202 0101
Ind1_4, 0102 0201 0102
Ind1_5, 0101 0202 0101
POP
Ind2_1, 0202 0101 0202
Ind2_2, 0201 0102 0201
Ind2_3, 0202 0101 0202
Ind2_4, 0201 0102 0201
Ind2_5, 0202 0101 0202
POP
Ind3_1, 0101 0101 0101
Ind3_2, 0102 0102 0102
Ind3_3, 0101 0101 0101
Ind3_4, 0102 0102 0102
Ind3_5, 0101 0101 0101
"""
        with open(example_file_path, 'w') as f:
            f.write(example_content)
