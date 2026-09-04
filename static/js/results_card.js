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
        // Auto-scroll the card into view so users never miss a result that
        // rendered below the fold. Skip if the container is already in a
        // split-pane results column (those are always visible side-by-side).
        if (!container.closest('.tool-split-results, .tool-results-pane')) {
            const card = container.querySelector('.results-card');
            if (card && typeof card.scrollIntoView === 'function') {
                // Delay to let layout settle
                setTimeout(function () {
                    card.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }, 50);
            }
        }
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
                    '<div class="rc-skeleton-row u-width80"></div>' +
                    '<div class="rc-skeleton-row u-width60"></div>' +
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

    // Render a consistent empty state — used by tool pages to prime the
    // results pane with a helpful placeholder before the user runs the tool.
    function showEmpty(containerOrId, opts) {
        const container = typeof containerOrId === 'string'
            ? document.getElementById(containerOrId)
            : containerOrId;
        if (!container) return;
        const icon = (opts && opts.icon) || 'fa-chart-bar';
        const title = (opts && opts.title) || 'No results yet';
        const body = (opts && opts.body) || 'Run the tool on the left to see results here.';
        container.innerHTML =
            '<div class="rc-empty">' +
                '<i class="fas ' + esc(icon) + '" aria-hidden="true"></i>' +
                '<div class="rc-empty-title">' + esc(title) + '</div>' +
                '<p class="rc-empty-body">' + esc(body) + '</p>' +
            '</div>';
    }

    global.ResultsCard = { mount: mount, build: buildCard, showLoading: showLoading, showError: showError, showEmpty: showEmpty };

    // D23: prime recognized results containers with a consistent empty state
    // on page load. Containers get a placeholder only if they are currently
    // empty / contain just a legacy "Run an analysis…" stub. Once the tool
    // mounts real content, this is overwritten.
    function primeEmptyStates() {
        // Id → {icon, title, body}
        var EMPTY_DEFAULTS = {
            resultsArea:        { icon: 'fa-chart-bar',       title: 'No results yet',    body: 'Paste a sequence and click Analyze.' },
            alignmentResults:   { icon: 'fa-align-center',    title: 'No alignment yet',  body: 'Enter two sequences and click Align.' },
            blastResults:       { icon: 'fa-rocket',          title: 'No BLAST run yet',  body: 'Paste a query and click Run BLAST Search.' },
            parseResults:       { icon: 'fa-file-alt',        title: 'No records parsed', body: 'Upload a file and click Parse.' },
            searchResults:      { icon: 'fa-search',          title: 'No search run yet', body: 'Enter a query and click Search.' },
            basicResults:       { icon: 'fa-cut',             title: 'No analysis yet',   body: 'Paste a sequence and click Analyze.' },
            advancedResults:    { icon: 'fa-filter',          title: 'No results yet',    body: 'Configure filters and run.' },
            mapResultsCard:     { icon: 'fa-map',             title: 'No map generated',  body: 'Enter a sequence + enzyme list.' },
            compatibleResults:  { icon: 'fa-link',            title: 'No pairs yet',      body: 'Enter 2+ enzymes to find compatible ends.' },
            popgenResults:      { icon: 'fa-users',           title: 'No data parsed',    body: 'Upload or paste a GENEPOP file.' },
            clusteringCard:     { icon: 'fa-circle-nodes',    title: 'No clustering run', body: 'Paste a data matrix and choose a method.' },
            treeInfo:           { icon: 'fa-tree',            title: 'No tree yet',       body: 'Parse or build a phylogenetic tree.' },
            readResults:        { icon: 'fa-book-open',       title: 'No record read',    body: 'Upload a file and click Read.' },
            indexResults:       { icon: 'fa-database',        title: 'No index built',    body: 'Upload a file and click Create Index.' },
            convertResults:     { icon: 'fa-exchange-alt',    title: 'No conversion yet', body: 'Upload a file and run a format conversion.' },
            filterResults:      { icon: 'fa-filter',          title: 'No filter applied', body: 'Upload a file and set thresholds.' },
            writeResults:       { icon: 'fa-file-export',     title: 'No output written', body: 'Upload, pick a format, click Write.' },
            codonResults:       { icon: 'fa-dna',             title: 'No translation',    body: 'Paste DNA and click Translate.' },
            iupacResults:       { icon: 'fa-atom',            title: 'No lookup yet',     body: 'Enter an IUPAC code to look up.' },
            globalResults:      { icon: 'fa-server',          title: 'No global search',  body: 'Enter a term to query every NCBI database.' },
            linkResults:        { icon: 'fa-link',            title: 'No linked records', body: 'Enter an ID and pick source/target.' },
            fetchResults:       { icon: 'fa-download',        title: 'No fetch yet',      body: 'Enter IDs and a format.' },
            infoResults:        { icon: 'fa-info-circle',     title: 'No database info',  body: 'Pick a database to fetch its info.' },
            orfResults:         { icon: 'fa-puzzle-piece',    title: 'No ORFs found yet', body: 'Paste a sequence and click Find ORFs.' },
            extractResults:     { icon: 'fa-cut',             title: 'No feature cut',    body: 'Enter range + strand.' },
            compoundResults:    { icon: 'fa-puzzle-piece',    title: 'No compound',       body: 'Enter sequence + locations.' },
            annotateResults:    { icon: 'fa-file-alt',        title: 'No annotations',   body: 'Define features and generate.' },
            contactsResults:    { icon: 'fa-network-wired',   title: 'No contacts yet',   body: 'Parse a structure first, then calculate contacts.' },
            interactionsResults:{ icon: 'fa-link',            title: 'No interactions',   body: 'Parse a structure first.' },
            dsspResults:        { icon: 'fa-dna',             title: 'No secondary structure', body: 'Upload a PDB and analyze.' },
            ramaResults:        { icon: 'fa-chart-scatter',   title: 'No phi/psi yet',    body: 'Upload a PDB to compute angles.' },
            sasaResults:        { icon: 'fa-water',           title: 'No SASA computed',  body: 'Upload a PDB to compute accessibility.' },
            selectionResults:   { icon: 'fa-crop',            title: 'No selection yet',  body: 'Upload a PDB and define a selection.' },
            superimposeResults: { icon: 'fa-project-diagram', title: 'No superposition',  body: 'Upload two PDBs to align.' },
            geometryResults:    { icon: 'fa-ruler',           title: 'No geometry',       body: 'Parse a structure first.' },
            qualityResults:     { icon: 'fa-check-circle',    title: 'No quality report', body: 'Upload a PDB to analyze quality.' },
            systemAnalysisResults:  { icon: 'fa-sitemap',     title: 'No system analysis', body: 'Add reactions, then Analyze System.' },
            networkAnalysisResults: { icon: 'fa-network-wired', title: 'No network analysis', body: 'Add reactions, then Analyze Network.' },
            pathwayVisualization:   { icon: 'fa-image',       title: 'No graph yet',      body: 'Add reactions, then Generate Visualization.' },
            hmmAnalysis:        { icon: 'fa-circle-nodes',    title: 'No model built',    body: 'Configure states and Build HMM.' },
            trainingResults:    { icon: 'fa-graduation-cap',  title: 'No training yet',   body: 'Build a model first, then train.' },
            decodingResults:    { icon: 'fa-search',          title: 'No decoding yet',   body: 'Build + train a model, then decode.' },
            literatureResults:  { icon: 'fa-book',            title: 'No articles fetched', body: 'Enter a PubMed query.' },
            nexusResults:       { icon: 'fa-file-code',       title: 'No Nexus parsed',   body: 'Upload or paste a Nexus file.' },
            scopResults:        { icon: 'fa-sitemap',         title: 'No SCOP lookup',    body: 'Enter a SCOP ID.' },
        };
        Object.keys(EMPTY_DEFAULTS).forEach(function (id) {
            var el = document.getElementById(id);
            if (!el) return;
            var txt = (el.textContent || '').trim().toLowerCase();
            var hasReal = el.querySelector('.results-card, .rc-stats, table, .accordion');
            if (hasReal) return;
            // Only prime if the container is empty-ish (default stub text)
            if (txt.length > 0 && txt.length > 120) return;
            global.ResultsCard.showEmpty(el, EMPTY_DEFAULTS[id]);
        });
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', primeEmptyStates);
    } else {
        primeEmptyStates();
    }
})(window);
