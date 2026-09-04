"""Gunicorn configuration.

Worker model
------------
Threaded workers are safe. They were not always: three things in this app used
to be process-global and unsynchronised, and each has been fixed.

  * utils.request_helpers.remote_timeout mutates the default socket timeout
    (BioPython's urllib-based helpers accept no timeout argument). It now holds
    a lock for the duration, which costs nothing because the outbound RateGate
    already spaces calls to each service.
  * matplotlib figures came from pyplot, which registers them in a global dict;
    app.py then swept it with plt.close('all') on every teardown, destroying
    figures other requests were still drawing. utils.plot_helpers.subplots()
    now builds Figure objects directly, so there is no registry and no sweep.
  * the outbound RateGate was already lock-protected.

Verified: 2 workers x 8 threads, 60 concurrent renders, all byte-identical.

`threads` therefore defaults to 4. Threads help here because requests spend
most of their time waiting on NCBI, KEGG and BLAST rather than on CPU.

Scaling out horizontally is still safe, and is what you want for CPU-bound
rendering. Set REDIS_URL whenever workers > 1: motif/BLAST/HMM/status state
lives in a shared store, but without REDIS_URL that store is per-process and a
user's follow-up request can land on a worker that has never seen their result.
"""
import multiprocessing
import os

bind = f"0.0.0.0:{os.environ.get('PORT', '9000')}"

# Long enough for a slow NCBI/KEGG round trip, which utils.request_helpers
# already caps (Entrez 30s, BLAST 180s). This must exceed the longest of those
# or gunicorn kills the worker mid-request.
timeout = int(os.environ.get('GUNICORN_TIMEOUT', '240'))
graceful_timeout = 30
keepalive = 5

workers = int(os.environ.get('WEB_CONCURRENCY',
                             min(multiprocessing.cpu_count(), 4)))
worker_class = 'gthread'
threads = int(os.environ.get('GUNICORN_THREADS', '4'))

max_requests = 500                # recycle workers: matplotlib/BioPython leak
max_requests_jitter = 50          # slowly across very long-lived processes

accesslog = '-'
errorlog = '-'
loglevel = os.environ.get('LOG_LEVEL', 'info').lower()


def on_starting(server):
    """Refuse to scale out silently onto per-process state."""
    if workers > 1 and not os.environ.get('REDIS_URL'):
        server.log.warning(
            'Running %d workers without REDIS_URL. Motif, BLAST and HMM '
            'results are held per-process, so a follow-up request that lands '
            'on a different worker will report the result as missing. Set '
            'REDIS_URL, or run with WEB_CONCURRENCY=1.',
            workers,
        )
    if not os.environ.get('SECRET_KEY'):
        server.log.warning(
            'SECRET_KEY is not set; each worker generates its own ephemeral '
            'key, so sessions break as requests move between workers.'
        )
