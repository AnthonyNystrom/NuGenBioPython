# NuGenBioPython — production image.
#
# Runs under gunicorn (see gunicorn.conf.py), not the Flask dev server, which
# prints its own "do not use in a production deployment" warning for good
# reason. Scaling is horizontal: set WEB_CONCURRENCY and REDIS_URL together,
# never WEB_CONCURRENCY alone.
FROM python:3.12-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Build-time only: matplotlib/scipy/reportlab wheels cover most of this, but a
# source build of any of them needs a compiler.
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential curl \
 && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt gunicorn

COPY . .

# Non-root. UPLOAD_FOLDER must be writable by this user.
RUN useradd --create-home --uid 10001 nugen \
 && mkdir -p /app/uploads \
 && chown -R nugen:nugen /app
USER nugen

EXPOSE 9000

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD curl -fsS http://127.0.0.1:9000/ >/dev/null || exit 1

CMD ["gunicorn", "-c", "gunicorn.conf.py", "app:app"]
