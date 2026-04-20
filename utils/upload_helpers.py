"""
Shared helpers for uploaded-file handling.

Guarantees cleanup of files written to UPLOAD_FOLDER and avoids collisions
between concurrent requests by prefixing each stored name with a UUID.
"""
import os
import uuid
from contextlib import contextmanager

from flask import current_app
from werkzeug.utils import secure_filename


class UploadError(ValueError):
    """Raised when an uploaded file is missing or invalid."""


def _build_filepath(file_storage):
    original = secure_filename(file_storage.filename or '')
    if not original:
        original = 'upload.bin'
    unique = f"{uuid.uuid4().hex}_{original}"
    return os.path.join(current_app.config['UPLOAD_FOLDER'], unique)


@contextmanager
def saved_upload(file_storage):
    """
    Persist an uploaded FileStorage to disk and yield its filepath.

    The file is always removed on exit, whether the caller returns normally
    or raises. Raises UploadError if file_storage is falsy or has no filename.
    """
    if not file_storage or not file_storage.filename:
        raise UploadError('No file provided')

    filepath = _build_filepath(file_storage)
    file_storage.save(filepath)
    try:
        yield filepath
    finally:
        try:
            if os.path.exists(filepath):
                os.remove(filepath)
        except OSError:
            pass


def require_uploaded_file(request, field='file'):
    """Return the FileStorage at request.files[field] or raise UploadError."""
    if field not in request.files:
        raise UploadError(f"Missing '{field}' in upload")
    fs = request.files[field]
    if not fs or not fs.filename:
        raise UploadError('No file provided')
    return fs


def stage_upload(file_storage):
    """Save an uploaded FileStorage to disk with a UUID-prefixed filename
    and return its filepath. Caller is responsible for cleanup (use with
    try/finally); prefer saved_upload() when the scope is a clean block."""
    if not file_storage or not file_storage.filename:
        raise UploadError('No file provided')
    filepath = _build_filepath(file_storage)
    file_storage.save(filepath)
    return filepath


def cleanup_upload(filepath):
    """Remove a staged upload file; silent no-op if already gone."""
    if not filepath:
        return
    try:
        if os.path.exists(filepath):
            os.remove(filepath)
    except OSError:
        pass
