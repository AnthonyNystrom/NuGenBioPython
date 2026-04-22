// ResultsCard — one visual component every tool uses for results.
//
// Emits:
//   ┌── title ─ meta ───────────── Save · Copy · Download ▾ ──┐
//   │ [ Summary | Details | Raw ]                              │
//   ├──────────────────────────────────────────────────────────┤
//   │                                                          │
//   │  (selected tab's content)                                │
//   │                                                          │
//   └──────────────────────────────────────────────────────────┘
//
// Usage:
//   ResultsCard.mount('sequenceResults', {
//       title: 'Sequence Analysis',
//       meta:  '<span class="meta">DNA · 300 bp</span>',
//       summary:  '<div>…</div>',            // HTML; Summary tab
//       details:  '<div>…</div>',            // HTML; Details tab
//       raw:      'ACGTACGT…',               // plain text; Raw tab (monospace, scrollable)
//       downloads: [
//           { label: 'FASTA', filename: 'seq.fasta', text: '>seq\nACGT' },
//           { label: 'JSON',  filename: 'seq.json',  text: JSON.stringify(data, null, 2) },
//       ],
//       workspaceItem: {                     // Save-to-workspace button payload
//           type: 'sequence', name: 'Analyzed DNA', data: 'ACGTACGT', meta: {kind:'dna'}
//       },
//       copyText: 'ACGTACGT',                // what Copy writes to clipboard; defaults to raw
//   });
//
// Every section is OPTIONAL — tabs that have no content are hidden.
(function (global) {
    'use strict';

    const esc = (global.NuGenUtils && global.NuGenUtils.escapeHtml) ||
                function (s) {
                    return String(s == null ? '' : s)
                        .replace(/&/g, '&amp;').replace(/</g, '&lt;')
                        .replace(/>/g, '&gt;').replace(/"/g, '&quot;')
                        .replace(/'/g, '&#39;');
                };

    let _cardSeq = 0;

    function buildCard(spec) {
        _cardSeq += 1;
        const prefix = spec.prefix || ('rc' + _cardSeq);

        const tabs = [];
        if (spec.summary) tabs.push({ id: 'summary', title: 'Summary', content: spec.summary, active: true });
        if (spec.details) tabs.push({ id: 'details', title: 'Details', content: spec.details, active: !spec.summary });
        if (spec.raw != null && spec.raw !== '') {
            const rawHtml =
                '<pre class="results-raw">' + esc(spec.raw) + '</pre>';
            tabs.push({ id: 'raw', title: 'Raw', content: rawHtml, active: !spec.summary && !spec.details });
        }
        // Allow custom extra tabs (e.g. Formatted for alignment)
        if (Array.isArray(spec.extraTabs)) {
            spec.extraTabs.forEach(function (t) { tabs.push(t); });
        }

        const tabsHtml = global.ResultsPanel && global.ResultsPanel.tabs
            ? global.ResultsPanel.tabs(tabs, { prefix: prefix })
            : tabs.map(function (t) { return t.content; }).join('');

        const downloads = Array.isArray(spec.downloads) ? spec.downloads : [];

        let actions = '';
        if (spec.workspaceItem) {
            actions += '<button type="button" class="btn-rc-action" data-rc-save aria-label="Save result to workspace">' +
                '<i class="fas fa-layer-group" aria-hidden="true"></i><span>Save</span></button>';
        }
        if (spec.copyText || spec.raw) {
            actions += '<button type="button" class="btn-rc-action" data-rc-copy aria-label="Copy result to clipboard">' +
                '<i class="far fa-copy" aria-hidden="true"></i><span>Copy</span></button>';
        }
        if (downloads.length) {
            const items = downloads.map(function (d, i) {
                return '<li><button type="button" class="dropdown-item" data-rc-download="' + i + '">' +
                    esc(d.label) + '</button></li>';
            }).join('');
            actions += '<div class="dropdown d-inline-block">' +
                '<button type="button" class="btn-rc-action" data-bs-toggle="dropdown" aria-expanded="false" aria-label="Download result">' +
                    '<i class="fas fa-download" aria-hidden="true"></i><span>Download</span>' +
                '</button>' +
                '<ul class="dropdown-menu dropdown-menu-end">' + items + '</ul>' +
            '</div>';
        }

        const header =
            '<header class="results-card-header">' +
                '<div class="results-card-title">' +
                    (spec.title ? '<strong>' + esc(spec.title) + '</strong>' : '') +
                    (spec.meta ? '<span class="results-card-meta">' + spec.meta + '</span>' : '') +
                '</div>' +
                (actions ? '<div class="results-card-actions">' + actions + '</div>' : '') +
            '</header>';

        return '<section class="results-card" data-rc-id="' + prefix + '">' +
                   header +
                   '<div class="results-card-body">' + tabsHtml + '</div>' +
               '</section>';
    }

    function wire(container, spec) {
        const saveBtn = container.querySelector('[data-rc-save]');
        if (saveBtn && spec.workspaceItem) {
            saveBtn.addEventListener('click', function () {
                if (!global.Workspace) return;
                const entry = global.Workspace.add(spec.workspaceItem);
                if (global.showAlert) {
                    global.showAlert('Saved "' + entry.name + '" to workspace', 'success');
                }
            });
        }
        const copyBtn = container.querySelector('[data-rc-copy]');
        if (copyBtn) {
            copyBtn.addEventListener('click', function () {
                const text = spec.copyText || spec.raw || '';
                if (!text) return;
                if (navigator.clipboard && navigator.clipboard.writeText) {
                    navigator.clipboard.writeText(text).then(function () {
                        if (global.showAlert) global.showAlert('Copied to clipboard', 'success');
                    });
                }
            });
        }
        container.querySelectorAll('[data-rc-download]').forEach(function (btn) {
            btn.addEventListener('click', function () {
                const idx = parseInt(btn.dataset.rcDownload, 10);
                const d = spec.downloads && spec.downloads[idx];
                if (!d) return;
                const text = typeof d.text === 'function' ? d.text() : d.text;
                const blob = new Blob([text || ''], { type: d.mime || 'text/plain' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = d.filename || 'download.txt';
                a.click();
                URL.revokeObjectURL(url);
            });
        });
    }

    function mount(containerOrId, spec) {
        const container = typeof containerOrId === 'string'
            ? document.getElementById(containerOrId)
            : containerOrId;
        if (!container) return;
        container.innerHTML = buildCard(spec);
        wire(container, spec);
        return container;
    }

    // Render a consistent loading skeleton in the target container while a
    // long-running analysis fetches results. Tools should call this right
    // before kicking off the request, then call mount() with real data on
    // success — mount() replaces the skeleton in one swap so layout doesn't jump.
    function showLoading(containerOrId, opts) {
        const container = typeof containerOrId === 'string'
            ? document.getElementById(containerOrId)
            : containerOrId;
        if (!container) return;
        const title = (opts && opts.title) || 'Working…';
        container.innerHTML =
            '<section class="results-card results-card-loading" aria-busy="true">' +
                '<header class="results-card-header">' +
                    '<div class="results-card-title">' +
                        '<strong>' + esc(title) + '</strong>' +
                    '</div>' +
                '</header>' +
                '<div class="results-card-body">' +
                    '<div class="rc-skeleton-row"></div>' +
                    '<div class="rc-skeleton-row" style="width:80%"></div>' +
                    '<div class="rc-skeleton-row" style="width:60%"></div>' +
                    '<div class="rc-skeleton-block"></div>' +
                '</div>' +
            '</section>';
    }

    // Render a consistent error state. Tools should call this from their
    // .catch / non-success paths so users see the same recovery UX everywhere.
    function showError(containerOrId, opts) {
        const container = typeof containerOrId === 'string'
            ? document.getElementById(containerOrId)
            : containerOrId;
        if (!container) return;
        const title = (opts && opts.title) || 'Something went wrong';
        const message = (opts && opts.message) || 'The request failed. Please try again.';
        container.innerHTML =
            '<section class="results-card results-card-error" role="alert">' +
                '<header class="results-card-header">' +
                    '<div class="results-card-title">' +
                        '<strong><i class="fas fa-exclamation-triangle"></i> ' + esc(title) + '</strong>' +
                    '</div>' +
                '</header>' +
                '<div class="results-card-body">' +
                    '<p class="rc-error-text">' + esc(message) + '</p>' +
                '</div>' +
            '</section>';
    }

    global.ResultsCard = { mount: mount, build: buildCard, showLoading: showLoading, showError: showError };
})(window);
