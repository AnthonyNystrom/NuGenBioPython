/* External service health probe.
 *
 * On page load, fires a lightweight probe against each external service
 * the current page depends on. If the probe fails or times out, we render
 * a "Service degraded" banner above the tool so users know why their
 * requests may hang or fail.
 *
 * Probes are cached per-session (sessionStorage) so we don't hammer each
 * service on every navigation.
 */
(function () {
    'use strict';

    // Map URL → array of services that page depends on.
    // Hubs (/compare, /data) inherit from their sub-tools since the hub landing
    // can render the BLAST/Entrez sub-panel inline.
    var PAGE_SERVICES = {
        '/kegg':     ['kegg'],
        '/pathway':  ['kegg'],  // pathway uses KEGG fetch in some routes
        '/database': ['ncbi'],
        '/data':     ['ncbi'],  // Data hub wraps database sub-tool
        '/swissprot':['uniprot'],
        '/blast':    ['ncbi-blast'],
        '/compare':  ['ncbi-blast'],  // Compare hub includes BLAST sub-tool
        '/unigene':  ['unigene'],
        '/hmm':      ['ncbi'],  // PubMed literature
    };

    var SERVICE_META = {
        'kegg':       { label: 'KEGG (rest.kegg.jp)',         probe: '/api/status/kegg',       typical: '1-5s' },
        'ncbi':       { label: 'NCBI Entrez',                 probe: '/api/status/ncbi',       typical: '1-10s' },
        'ncbi-blast': { label: 'NCBI BLAST',                  probe: '/api/status/ncbi-blast', typical: '30-90s per query' },
        'uniprot':    { label: 'UniProt (rest.uniprot.org)',  probe: '/api/status/uniprot',    typical: '1-5s' },
        'unigene':    { label: 'NCBI UniGene',                probe: null,                     typical: 'retired',
                        degraded_msg: 'UniGene was retired by NCBI in 2019. External fetch is unavailable — local file parsing still works.' },
    };

    var PROBE_TIMEOUT = 5000;
    var CACHE_TTL = 5 * 60 * 1000;  // 5 min

    function renderBanner(service, status, detail) {
        var host = document.getElementById('app-alerts');
        if (!host) return;
        var existing = host.querySelector('[data-service-status="' + service + '"]');
        if (existing) existing.remove();
        if (status === 'ok') return;

        var meta = SERVICE_META[service] || { label: service };
        var msg;
        if (status === 'retired') {
            msg = meta.degraded_msg || (meta.label + ' is retired.');
        } else if (status === 'timeout') {
            msg = meta.label + ' is slow to respond. Requests may take a while or time out.';
        } else {
            msg = meta.label + ' is unreachable right now. External operations may fail — local parsing still works.';
        }

        var el = document.createElement('div');
        el.className = 'alert alert-warning';
        el.setAttribute('data-service-status', service);
        el.setAttribute('role', 'status');
        el.innerHTML =
            '<i class="fas fa-exclamation-triangle me-2"></i>' +
            '<strong>' + escapeHtml(meta.label) + '</strong> — ' +
            escapeHtml(msg) +
            (detail ? ' <small class="text-muted">(' + escapeHtml(detail) + ')</small>' : '');
        host.appendChild(el);
    }

    function escapeHtml(s) {
        return String(s == null ? '' : s)
            .replace(/&/g, '&amp;').replace(/</g, '&lt;')
            .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }

    function probe(service) {
        var meta = SERVICE_META[service];
        if (!meta) return Promise.resolve({ status: 'ok' });
        if (!meta.probe) {
            return Promise.resolve({ status: 'retired' });
        }
        var cacheKey = 'svc-status:' + service;
        try {
            var cached = sessionStorage.getItem(cacheKey);
            if (cached) {
                var parsed = JSON.parse(cached);
                if (parsed && parsed.t && (Date.now() - parsed.t) < CACHE_TTL) {
                    return Promise.resolve({ status: parsed.status, detail: parsed.detail });
                }
            }
        } catch (_) {}

        var ac = new AbortController();
        var to = setTimeout(function () { ac.abort(); }, PROBE_TIMEOUT);
        return fetch(meta.probe, { signal: ac.signal })
            .then(function (r) { return r.json().catch(function () { return null; }); })
            .then(function (data) {
                clearTimeout(to);
                var status = (data && data.ok) ? 'ok' : 'down';
                var detail = data && data.latency_ms ? data.latency_ms + 'ms' : null;
                try { sessionStorage.setItem(cacheKey, JSON.stringify({ status: status, detail: detail, t: Date.now() })); } catch (_) {}
                return { status: status, detail: detail };
            })
            .catch(function (e) {
                clearTimeout(to);
                var status = (e && e.name === 'AbortError') ? 'timeout' : 'down';
                try { sessionStorage.setItem(cacheKey, JSON.stringify({ status: status, t: Date.now() })); } catch (_) {}
                return { status: status };
            });
    }

    function checkForCurrentPage() {
        var path = window.location.pathname;
        var services = PAGE_SERVICES[path];
        if (!services) return;
        services.forEach(function (svc) {
            probe(svc).then(function (result) {
                renderBanner(svc, result.status, result.detail);
            });
        });
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', checkForCurrentPage);
    } else {
        checkForCurrentPage();
    }

    window.ServiceStatus = { probe: probe, renderBanner: renderBanner };
})();
