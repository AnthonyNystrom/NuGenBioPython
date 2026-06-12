/**
 * Restriction Enzymes - Complete JavaScript
 * Handles all tabs and functionality
 */

// Global state
let availableEnzymes = [];
let lastAnalysisResults = null;
let lastSequence = '';
let uploadedFileInfo = {};

// ============================================================================
// FILE UPLOAD HANDLERS
// ============================================================================

function handleBasicFileUpload(event) {
    handleFileUpload(event, 'basicSequence', 'basicFileInfo');
}

function handleAdvancedFileUpload(event) {
    handleFileUpload(event, 'advancedSequence', 'advancedFileInfo');
}

function handleMapFileUpload(event) {
    handleFileUpload(event, 'mapSequence', 'mapFileInfo');
}

function handleFileUpload(event, targetTextareaId, infoId) {
    const file = event.target.files[0];
    if (!file) return;

    const formData = new FormData();
    formData.append('file', file);

    // Show loading
    const infoDiv = document.getElementById(infoId);
    infoDiv.style.display = 'block';
    infoDiv.innerHTML = '<small class="text-muted"><i class="fas fa-spinner fa-spin"></i> Parsing file...</small>';

    fetch('/api/restriction/upload_sequence', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            // Populate textarea
            document.getElementById(targetTextareaId).value = data.sequence;

            // Show file info
            const info = data.info;
            let infoHtml = `<div class="alert alert-success py-1 px-2 small mb-0">`;
            infoHtml += `<i class="fas fa-check-circle"></i> <strong>${data.filename}</strong><br>`;
            infoHtml += `Format: ${info.format} | Length: ${info.length.toLocaleString()} bp`;

            if (info.id) infoHtml += `<br>ID: ${info.id}`;
            if (info.description) infoHtml += `<br>Description: ${info.description}`;
            if (info.features) infoHtml += ` | Features: ${info.features}`;

            infoHtml += `</div>`;
            infoDiv.innerHTML = infoHtml;

            // Store file info
            uploadedFileInfo[targetTextareaId] = data.info;

            showAlert(`File loaded: ${data.filename} (${info.length.toLocaleString()} bp)`, 'success');
        } else {
            infoDiv.innerHTML = `<small class="text-danger"><i class="fas fa-exclamation-triangle"></i> ${escapeHtml(friendlyError(data.error, 'server'))}</small>`;
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        infoDiv.innerHTML = `<small class="text-danger"><i class="fas fa-exclamation-triangle"></i> Upload failed</small>`;
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

// ============================================================================
// INITIALIZATION
// ============================================================================

window.addEventListener('load', function() {
    loadBasicEnzymeList();
    setupEventListeners();

    // Load enzyme browser when tab is shown
    document.getElementById('browser-tab').addEventListener('shown.bs.tab', function() {
        if (document.getElementById('enzymeBrowserResults').children.length === 1) {
            loadEnzymeBrowser();
        }
    });
});

function setupEventListeners() {
    // Basic Analysis Form
    document.getElementById('basicAnalysisForm').addEventListener('submit', handleBasicAnalysis);

    // Advanced Analysis Form
    document.getElementById('advancedAnalysisForm').addEventListener('submit', handleAdvancedAnalysis);

    // Map Form
    document.getElementById('mapForm').addEventListener('submit', handleMapGeneration);

    // Compatible Ends Form
    document.getElementById('compatibleEndsForm').addEventListener('submit', handleCompatibleEnds);
}

// ============================================================================
// BASIC ANALYSIS TAB
// ============================================================================

function loadBasicEnzymeList() {
    fetch('/api/restriction/list_enzymes?filter=common')
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            availableEnzymes = data.enzymes;
            displayBasicEnzymeCheckboxes(data.enzymes);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
            loadFallbackEnzymes();
        }
    })
    .catch(error => {
        showAlert(friendlyError(error, 'server'), 'danger');
        loadFallbackEnzymes();
    });
}

function loadFallbackEnzymes() {
    const fallbackEnzymes = [
        {name: 'EcoRI', site: 'GAATTC', overhang: -4, overhang_type: '5-overhang'},
        {name: 'BamHI', site: 'GGATCC', overhang: -4, overhang_type: '5-overhang'},
        {name: 'HindIII', site: 'AAGCTT', overhang: -4, overhang_type: '5-overhang'},
        {name: 'XbaI', site: 'TCTAGA', overhang: -4, overhang_type: '5-overhang'},
        {name: 'SalI', site: 'GTCGAC', overhang: -4, overhang_type: '5-overhang'},
        {name: 'PstI', site: 'CTGCAG', overhang: -4, overhang_type: '3-overhang'}
    ];
    availableEnzymes = fallbackEnzymes;
    displayBasicEnzymeCheckboxes(fallbackEnzymes);
}

function displayBasicEnzymeCheckboxes(enzymes) {
    const container = document.getElementById('basicEnzymeCheckboxes');
    if (!container) return;

    let html = '';
    enzymes.forEach(enzyme => {
        const overhangBadge = enzyme.is_blunt ?
            '<span class="badge bg-secondary">Blunt</span>' :
            `<span class="badge bg-info">${enzyme.overhang_type}</span>`;

        html += `
            <div class="col-md-4 col-sm-6 mb-1">
                <div class="form-check">
                    <input class="form-check-input" type="checkbox" value="${enzyme.name}" id="basic_enzyme_${enzyme.name}">
                    <label class="form-check-label small" for="basic_enzyme_${enzyme.name}">
                        <strong>${enzyme.name}</strong> ${overhangBadge}<br>
                        <small class="text-muted">${enzyme.site}</small>
                    </label>
                </div>
            </div>
        `;
    });

    container.innerHTML = html;

    // Auto-select common enzymes
    setTimeout(() => {
        selectCommonBasicEnzymes();
    }, 100);
}

function selectAllBasicEnzymes() {
    document.querySelectorAll('#basicEnzymeCheckboxes input[type="checkbox"]').forEach(cb => cb.checked = true);
}

function clearAllBasicEnzymes() {
    document.querySelectorAll('#basicEnzymeCheckboxes input[type="checkbox"]').forEach(cb => cb.checked = false);
}

function selectCommonBasicEnzymes() {
    clearAllBasicEnzymes();
    const commonEnzymes = ['EcoRI', 'BamHI', 'HindIII', 'XbaI', 'SalI', 'PstI'];
    commonEnzymes.forEach(enzyme => {
        const checkbox = document.getElementById(`basic_enzyme_${enzyme}`);
        if (checkbox) checkbox.checked = true;
    });
}

function handleBasicAnalysis(e) {
    e.preventDefault();

    const sequence = document.getElementById('basicSequence').value.trim().toUpperCase();
    const checkedEnzymes = Array.from(document.querySelectorAll('#basicEnzymeCheckboxes input:checked')).map(cb => cb.value);
    const includeSequences = document.getElementById('includeFragmentSequences').checked;

    if (!sequence) {
        showAlert('Please enter a DNA sequence', 'warning');
        return;
    }

    if (checkedEnzymes.length === 0) {
        showAlert('Please select at least one enzyme', 'warning');
        return;
    }

    if (!/^[ATGC]+$/i.test(sequence)) {
        showAlert('Invalid DNA sequence. Only A, T, G, C allowed.', 'danger');
        return;
    }

    showLoading('basicAnalyzeBtn');

    fetch('/api/restriction/analyze', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            enzymes: checkedEnzymes,
            include_sequences: includeSequences
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('basicAnalyzeBtn', '<i class="fas fa-cut me-2"></i>Analyze Restriction Sites');

        if (data.success) {
            lastAnalysisResults = data;
            lastSequence = sequence;
            displayBasicResults(data.analysis, sequence.length);
            document.getElementById('exportBasicBtn').style.display = 'block';
            // D20 flattening: prime the Map tab with the same sequence + enzymes
            // so jumping to Map view is one click (no re-entry). Same for Advanced.
            const mapSeq = document.getElementById('mapSequence');
            const mapEnz = document.getElementById('mapEnzymes');
            if (mapSeq && !mapSeq.value.trim()) mapSeq.value = sequence;
            if (mapEnz && !mapEnz.value.trim()) mapEnz.value = checkedEnzymes.join(', ');
            const advSeq = document.getElementById('advancedSequence');
            if (advSeq && !advSeq.value.trim()) advSeq.value = sequence;
        } else {
            if (typeof ResultsCard !== 'undefined' && ResultsCard.showError) {
                ResultsCard.showError('basicResults', { title: 'Restriction analysis failed', message: data.error });
            }
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('basicAnalyzeBtn', '<i class="fas fa-cut me-2"></i>Analyze Restriction Sites');
        if (typeof ResultsCard !== 'undefined' && ResultsCard.showError) {
            ResultsCard.showError('basicResults', { title: 'Request failed', message: friendlyError(error, 'server') });
        }
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

function displayBasicResults(analysis, sequenceLength) {
    // The sibling containers (#basicSummary, #basicFragmentChart) are cleared
    // so the old split panels collapse; everything renders inside one ResultsCard.
    const summaryEl = document.getElementById('basicSummary');
    const chartEl = document.getElementById('basicFragmentChart');
    if (summaryEl) summaryEl.innerHTML = '<p class="text-muted small mb-0">See results card →</p>';
    if (chartEl) chartEl.innerHTML = '';

    // --- Summary: headline stats + fragment distribution ---
    let totalCuts = 0, enzymesCutting = 0, totalEnzymes = 0;
    const allFragments = [];
    Object.values(analysis).forEach(r => {
        if (r.error) { totalEnzymes++; return; }
        totalEnzymes++;
        if (r.number_of_cuts > 0) { totalCuts += r.number_of_cuts; enzymesCutting++; }
        if (r.fragments) allFragments.push(...r.fragments);
    });
    const summary =
        '<div class="rc-stats">' +
        '<div class="rc-stat"><div class="rc-stat-label">Sequence length</div>' +
            '<div class="rc-stat-value">' + sequenceLength + '</div><div class="rc-stat-sub">bp</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Total cut sites</div>' +
            '<div class="rc-stat-value">' + totalCuts + '</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Enzymes cutting</div>' +
            '<div class="rc-stat-value">' + enzymesCutting + '</div>' +
            '<div class="rc-stat-sub">of ' + totalEnzymes + '</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Fragments</div>' +
            '<div class="rc-stat-value">' + allFragments.length + '</div></div>' +
        '</div>' +
        (allFragments.length ? buildFragmentHistogram(allFragments) : '');

    // --- Details: per-enzyme table ---
    let details = '<div class="table-responsive"><table class="table table-sm table-hover mb-0"><thead><tr>' +
        '<th>Enzyme</th><th>Site</th><th>Type</th><th>Cuts</th><th>Positions</th><th>Fragments</th></tr></thead><tbody>';
    Object.entries(analysis).forEach(([enzyme, r]) => {
        if (r.error) {
            details += '<tr><td>' + enzyme + '</td><td colspan="5" class="text-danger small">' + r.error + '</td></tr>';
            return;
        }
        const badge = r.is_blunt
            ? '<span class="badge bg-secondary">Blunt</span>'
            : '<span class="badge bg-info-subtle text-info-emphasis border border-info-subtle">' + r.overhang_type + '</span>';
        const frags = r.fragments.length ? r.fragments.join(', ') + ' bp' : 'No cuts';
        const positions = r.cut_positions.length ? r.cut_positions.join(', ') : '—';
        details += '<tr>' +
            '<td><strong>' + enzyme + '</strong></td>' +
            '<td><code class="small">' + r.recognition_site + '</code></td>' +
            '<td>' + badge + '</td>' +
            '<td><span class="badge bg-primary-subtle text-primary-emphasis border border-primary-subtle">' + r.number_of_cuts + '</span></td>' +
            '<td class="small text-muted">' + positions + '</td>' +
            '<td class="small">' + frags + '</td>' +
        '</tr>';
        if (r.fragment_sequences && r.fragment_sequences.length) {
            const fragHtml = r.fragment_sequences.map((s, i) =>
                '<div><strong>Fragment ' + (i + 1) + ':</strong> <code>' + s + '</code></div>').join('');
            details += '<tr><td colspan="6" class="bg-light small">' + fragHtml + '</td></tr>';
        }
    });
    details += '</tbody></table></div>';

    // --- Raw: tab-separated CSV of enzyme/cuts ---
    const rawLines = ['Enzyme\tSite\tCuts\tPositions\tFragments\tType'];
    Object.entries(analysis).forEach(([enzyme, r]) => {
        if (r.error) { rawLines.push(enzyme + '\t' + r.error); return; }
        rawLines.push([enzyme, r.recognition_site, r.number_of_cuts,
                       r.cut_positions.join(';'), r.fragments.join(';'), r.overhang_type].join('\t'));
    });
    const raw = rawLines.join('\n');

    ResultsCard.mount('basicResults', {
        title: 'Restriction Analysis',
        meta:  '<span class="meta">' + totalCuts + ' cuts · ' + enzymesCutting + ' enzymes · ' + allFragments.length + ' fragments</span>',
        summary:  summary,
        details:  details,
        raw:      raw,
        copyText: raw,
        downloads: [
            { label: 'Table (TSV)',  filename: 'restriction.tsv',  text: raw, mime: 'text/tab-separated-values' },
            { label: 'Table (JSON)', filename: 'restriction.json',
              text: JSON.stringify({ sequence_length: sequenceLength, analysis: analysis }, null, 2),
              mime: 'application/json' },
        ],
        workspaceItem: {
            type: 'alignment',    // closest existing bucket; UI groups by type
            name: 'Restriction · ' + enzymesCutting + ' enzymes',
            data: raw,
            meta: { source: 'restriction', sequence_length: sequenceLength },
        },
    });
}

function buildFragmentHistogram(fragments) {
    const bins  = [0, 100, 500, 1000, 2000, 5000, 10000, Infinity];
    const labels= ['<100', '100-500', '500-1K', '1K-2K', '2K-5K', '5K-10K', '>10K'];
    const counts = new Array(bins.length - 1).fill(0);
    fragments.forEach(size => {
        for (let i = 0; i < bins.length - 1; i++) {
            if (size >= bins[i] && size < bins[i + 1]) { counts[i]++; break; }
        }
    });
    const max = Math.max(...counts);
    let html = '<div class="mt-3"><div class="small text-muted mb-2">Fragment size distribution</div>' +
               '<div class="d-flex align-items-end gap-2" style="height:80px;">';
    counts.forEach((c, i) => {
        const h = max > 0 ? Math.max((c / max) * 100, 3) : 3;
        html += '<div class="text-center flex-fill">' +
            '<div style="height:' + h + '%; background:var(--color-primary); border-radius:var(--radius-sm); min-height:3px;"></div>' +
            '<div class="rc-stat-sub mt-1">' + c + '</div>' +
            '<div class="rc-stat-sub">' + labels[i] + '</div>' +
        '</div>';
    });
    html += '</div></div>';
    return html;
}


// Superseded in R2 — histogram now rendered inside the ResultsCard Summary tab
// (see buildFragmentHistogram). This stub keeps any legacy caller happy.
function displayBasicFragmentChart() { /* no-op */ }

function loadBasicExample() {
    document.getElementById('basicSequence').value = 'GAATTCAAGCTTATCGATCGAATTCCTGCAGGGATCCAAGCTTTCTAGATGCATGCCTGCAGGAATTC';
    selectCommonBasicEnzymes();

    setTimeout(() => {
        document.getElementById('basicAnalysisForm').dispatchEvent(new Event('submit'));
    }, 500);
}

function exportBasicResults(format) {
    if (!lastAnalysisResults) {
        showAlert('No results to export. Run an analysis first.', 'warning');
        return;
    }

    exportResults(lastAnalysisResults, format);
}

// ============================================================================
// ADVANCED ANALYSIS TAB
// ============================================================================

function handleAdvancedAnalysis(e) {
    e.preventDefault();

    const sequence = document.getElementById('advancedSequence').value.trim().toUpperCase();
    const filterType = document.getElementById('filterType').value;
    const minCuts = parseInt(document.getElementById('minCuts').value) || 0;
    const maxCuts = document.getElementById('maxCuts').value ? parseInt(document.getElementById('maxCuts').value) : null;
    const useAllEnzymes = document.getElementById('useAllEnzymes').checked;

    if (!sequence) {
        showAlert('Please enter a DNA sequence', 'warning');
        return;
    }

    if (!/^[ATGC]+$/i.test(sequence)) {
        showAlert('Invalid DNA sequence. Only A, T, G, C allowed.', 'danger');
        return;
    }

    showLoading('advancedAnalyzeBtn');

    const payload = {
        sequence: sequence,
        filter: filterType,
        min_cuts: minCuts,
        max_cuts: maxCuts
    };

    if (!useAllEnzymes) {
        // Get selected enzymes from basic tab
        const selected = Array.from(document.querySelectorAll('#basicEnzymeCheckboxes input:checked')).map(cb => cb.value);
        if (selected.length === 0) {
            showAlert('Please select some enzymes in Basic Analysis tab or check "Use all enzymes"', 'warning');
            hideLoading('advancedAnalyzeBtn', '<i class="fas fa-filter me-2"></i>Run Advanced Analysis');
            return;
        }
        payload.enzymes = selected;
    }

    fetch('/api/restriction/advanced_analysis', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(payload)
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('advancedAnalyzeBtn', '<i class="fas fa-filter me-2"></i>Run Advanced Analysis');

        if (data.success) {
            displayAdvancedResults(data.result);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('advancedAnalyzeBtn', '<i class="fas fa-filter me-2"></i>Run Advanced Analysis');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

function displayAdvancedResults(result) {
    const tiles = [
        { label: 'Sequence',   value: result.sequence_length + ' bp' },
        { label: 'Tested',     value: result.total_enzymes_tested },
        { label: 'Matched',    value: result.enzyme_count },
        { label: 'Filters',    value: result.filters_applied.length },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const filterNote = `<p class="small text-muted mb-3"><strong>Filters applied:</strong> ${result.filters_applied.join(', ') || 'none'}</p>`;

    let body;
    if (result.enzymes.length === 0) {
        body = '<p class="text-muted text-center">No enzymes match the specified criteria.</p>';
    } else {
        let rows = '';
        result.enzymes.forEach(enzyme => {
            const overhangBadge = enzyme.is_blunt
                ? '<span class="badge bg-secondary">Blunt</span>'
                : `<span class="badge bg-info">${enzyme.overhang_type}</span>`;
            rows += `
                <tr>
                    <td><strong>${enzyme.name}</strong></td>
                    <td><code class="small">${enzyme.site}</code></td>
                    <td>${overhangBadge}</td>
                    <td><span class="badge bg-primary">${enzyme.cuts}</span></td>
                    <td><small class="text-muted">${enzyme.positions.join(', ')}</small></td>
                </tr>`;
        });
        body = `<div class="table-responsive"><table class="table table-hover table-sm">
            <thead class="table-light"><tr><th>Enzyme</th><th>Site</th><th>Type</th><th>Cuts</th><th>Positions</th></tr></thead>
            <tbody>${rows}</tbody></table></div>`;
    }

    const tsv = ['enzyme\tsite\toverhang\tcuts\tpositions'];
    result.enzymes.forEach(e => {
        tsv.push([e.name, e.site, e.is_blunt ? 'blunt' : e.overhang_type, e.cuts, e.positions.join(',')].join('\t'));
    });

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('advancedResults', {
            title: 'Advanced Restriction',
            meta: `${result.enzyme_count} / ${result.total_enzymes_tested} enzymes matched`,
            summary: tilesHtml + filterNote + body,
            raw: JSON.stringify(result, null, 2),
            workspaceItem: { type: 'restriction-advanced', name: `Advanced restriction (${result.enzyme_count} enzymes)`, data: result },
            downloads: [
                { label: 'TSV',  filename: 'restriction-advanced.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'restriction-advanced.json', text: JSON.stringify(result, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('advancedResults').innerHTML = tilesHtml + filterNote + body;
    }
}

function loadAdvancedRestrictionExample() {
    document.getElementById('advancedSequence').value = 'GAATTCAAGCTTATCGATCGAATTCCTGCAGGGATCCAAGCTTTCTAGATGCATGCCTGCAGGAATTC';
    document.getElementById('filterType').value = 'unique';
    document.getElementById('useAllEnzymes').checked = true;

    setTimeout(() => {
        document.getElementById('advancedAnalysisForm').dispatchEvent(new Event('submit'));
    }, 500);
}

// ============================================================================
// ENZYME BROWSER TAB
// ============================================================================

function loadEnzymeBrowser() {
    const filter = document.getElementById('enzymeDatabaseFilter').value;
    const search = document.getElementById('enzymeSearch').value;
    const siteSize = document.getElementById('siteSizeFilter').value;
    const overhang = document.getElementById('overhangFilter').value;

    const resultsDiv = document.getElementById('enzymeBrowserResults');
    resultsDiv.innerHTML = '<div class="text-center"><div class="spinner-border spinner-border-sm text-success"></div><span class="ms-2 small">Loading...</span></div>';

    let url = `/api/restriction/list_enzymes?filter=${filter}`;
    if (search) url += `&search=${encodeURIComponent(search)}`;
    if (siteSize) url += `&site_size=${siteSize}`;
    if (overhang) url += `&overhang=${overhang}`;

    fetch(url)
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            displayEnzymeBrowser(data.enzymes);
        } else {
            resultsDiv.innerHTML = `<p class="text-danger small">${escapeHtml(friendlyError(data.error, 'server'))}</p>`;
        }
    })
    .catch(error => {
        resultsDiv.innerHTML = `<p class="text-danger small">${escapeHtml(friendlyError(error, 'server'))}</p>`;
    });
}

function displayEnzymeBrowser(enzymes) {
    const resultsDiv = document.getElementById('enzymeBrowserResults');

    if (enzymes.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted text-center small">No enzymes match the criteria.</p>';
        return;
    }

    let html = `<p class="small mb-2"><strong>${enzymes.length} enzymes found</strong></p>`;
    html += '<div class="d-flex flex-wrap gap-1">';

    enzymes.forEach(enzyme => {
        const colorClass = enzyme.is_blunt ? 'bg-secondary' : 'bg-success';
        html += `
            <span class="badge ${colorClass} enzyme-badge"
                  onclick="showEnzymeDetails('${enzyme.name}')"
                  title="${enzyme.site} - ${enzyme.overhang_type}">
                ${enzyme.name}
            </span>
        `;
    });

    html += '</div>';
    resultsDiv.innerHTML = html;
}

function showEnzymeDetails(enzymeName) {
    fetch(`/api/restriction/enzyme_details/${enzymeName}`)
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            displayEnzymeDetails(data.enzyme);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

function displayEnzymeDetails(enzyme) {
    const card = document.getElementById('enzymeDetailsCard');
    const placeholder = document.getElementById('enzymeBrowserPlaceholder');
    const content = document.getElementById('enzymeDetailsContent');

    card.style.display = 'block';
    placeholder.style.display = 'none';

    const overhangBadge = enzyme.is_blunt ?
        '<span class="badge bg-secondary">Blunt</span>' :
        `<span class="badge bg-info">${enzyme.overhang_type}</span>`;

    let html = `
        <h6 class="text-success">${enzyme.name}</h6>
        <table class="table table-sm small mb-2">
            <tr><td><strong>Site:</strong></td><td><code>${enzyme.site}</code></td></tr>
            <tr><td><strong>Site Size:</strong></td><td>${enzyme.site_size} bp</td></tr>
            <tr><td><strong>Type:</strong></td><td>${overhangBadge}</td></tr>
            <tr><td><strong>Overhang:</strong></td><td>${enzyme.overhang}</td></tr>
    `;

    if (enzyme.overhang_seq) {
        html += `<tr><td><strong>Overhang Seq:</strong></td><td><code>${enzyme.overhang_seq}</code></td></tr>`;
    }

    html += '</table>';

    if (enzyme.suppliers && enzyme.suppliers.length > 0) {
        html += `<p class="small mb-2"><strong>Suppliers (${enzyme.supplier_count}):</strong></p>`;
        html += '<ul class="list-unstyled small">';
        enzyme.suppliers.slice(0, 5).forEach(supplier => {
            html += `<li><i class="fas fa-building text-primary"></i> ${supplier}</li>`;
        });
        if (enzyme.suppliers.length > 5) {
            html += `<li class="text-muted">... and ${enzyme.suppliers.length - 5} more</li>`;
        }
        html += '</ul>';
    }

    if (enzyme.isoschizomers && enzyme.isoschizomers.length > 0) {
        html += `<p class="small mb-1"><strong>Isoschizomers:</strong></p>`;
        html += `<p class="small text-muted">${enzyme.isoschizomers.join(', ')}</p>`;
    }

    content.innerHTML = html;
}

// ============================================================================
// RESTRICTION MAP TAB
// ============================================================================

function handleMapGeneration(e) {
    e.preventDefault();

    const sequence = document.getElementById('mapSequence').value.trim().toUpperCase();
    const enzymesInput = document.getElementById('mapEnzymes').value.trim();

    if (!sequence) {
        showAlert('Please enter a DNA sequence', 'warning');
        return;
    }

    if (!enzymesInput) {
        showAlert('Please enter enzyme names', 'warning');
        return;
    }

    if (!/^[ATGC]+$/i.test(sequence)) {
        showAlert('Invalid DNA sequence. Only A, T, G, C allowed.', 'danger');
        return;
    }

    const enzymes = enzymesInput.split(',').map(e => e.trim()).filter(e => e);

    showLoading('generateMapBtn');

    fetch('/api/restriction/analyze', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            enzymes: enzymes
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('generateMapBtn', '<i class="fas fa-map me-2"></i>Generate Restriction Map');

        if (data.success) {
            displayRestrictionMap(data.analysis, sequence.length, data.map);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('generateMapBtn', '<i class="fas fa-map me-2"></i>Generate Restriction Map');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

function displayRestrictionMap(analysis, sequenceLength, mapSvg) {
    const mapDiv = document.getElementById('mapVisualization');

    // Prefer the server-rendered vector map (collision-free label packing,
    // proper scale axis). If absent, fall back to the legacy HTML map.
    let html = '<div class="restriction-map">';
    if (mapSvg) {
        html += `<img src="${mapSvg}" alt="Restriction map" class="img-fluid" `
             + `style="width:100%; height:auto; border:1px solid var(--color-border, #e2e8f0); `
             + `border-radius:6px; background:#fff;">`;
    } else {
    html += `<h6 class="small mb-3">Linear Restriction Map (${sequenceLength} bp)</h6>`;

    // Draw scale
    html += '<div style="position: relative; height: 100px; margin-bottom: 20px;">';

    // Main sequence line
    html += '<div style="position: absolute; top: 50px; left: 0; right: 0; height: 3px; background: var(--color-text);"></div>';

    // Start and end markers
    html += '<div style="position: absolute; top: 45px; left: 0; font-size: 10px;">0 bp</div>';
    html += `<div style="position: absolute; top: 45px; right: 0; font-size: 10px;">${sequenceLength} bp</div>`;

    // Plot cut sites
    const colors = ['#dc3545', '#0d6efd', '#198754', '#fd7e14', '#6f42c1', '#0dcaf0'];
    let colorIndex = 0;

    Object.entries(analysis).forEach(([enzyme, result]) => {
        if (!result.error && result.cut_positions && result.cut_positions.length > 0) {
            const color = colors[colorIndex % colors.length];
            colorIndex++;

            result.cut_positions.forEach(pos => {
                const percentage = (pos / sequenceLength) * 100;
                html += `
                    <div style="position: absolute; top: 30px; left: ${percentage}%; width: 2px; height: 40px; background: ${color};">
                        <div style="position: absolute; top: -20px; left: -20px; font-size: 9px; color: ${color}; white-space: nowrap;">
                            ${enzyme}<br>${pos}
                        </div>
                    </div>
                `;
            });
        }
    });

    html += '</div>';

    // Legend
    html += '<div class="mt-4">';
    html += '<p class="small mb-2"><strong>Enzymes:</strong></p>';
    html += '<div class="d-flex flex-wrap gap-2">';

    colorIndex = 0;
    Object.entries(analysis).forEach(([enzyme, result]) => {
        if (!result.error && result.cut_positions && result.cut_positions.length > 0) {
            const color = colors[colorIndex % colors.length];
            colorIndex++;

            html += `
                <span class="badge" style="background-color: ${color};">
                    ${enzyme} (${result.number_of_cuts} cuts)
                </span>
            `;
        }
    });

    html += '</div>';
    html += '</div>';
    }  // end else (legacy HTML map)

    html += '</div>';
    mapDiv.innerHTML = html;

    // Summary card (tiles + per-enzyme table + downloads) next to the visualization
    const enzymes = Object.entries(analysis).filter(([_, r]) => !r.error && r.cut_positions);
    const totalCuts = enzymes.reduce((s, [_, r]) => s + r.cut_positions.length, 0);
    const tiles = [
        { label: 'Sequence', value: sequenceLength + ' bp' },
        { label: 'Enzymes',  value: enzymes.length },
        { label: 'Cuts',     value: totalCuts },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    enzymes.forEach(([enz, r]) => {
        rows += `<tr>
            <td><strong>${enz}</strong></td>
            <td class="text-center"><span class="badge bg-primary">${r.cut_positions.length}</span></td>
            <td class="small text-muted">${r.cut_positions.join(', ')}</td>
        </tr>`;
    });
    const table = `<div class="table-responsive"><table class="table table-sm">
        <thead class="table-light"><tr><th>Enzyme</th><th class="text-center">Cuts</th><th>Positions</th></tr></thead>
        <tbody>${rows}</tbody></table></div>`;

    const tsv = ['enzyme\tcuts\tpositions'];
    enzymes.forEach(([enz, r]) => tsv.push([enz, r.cut_positions.length, r.cut_positions.join(',')].join('\t')));

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('mapResultsCard', {
            title: 'Restriction Map',
            meta: `${enzymes.length} enzyme${enzymes.length === 1 ? '' : 's'} · ${totalCuts} cut${totalCuts === 1 ? '' : 's'}`,
            summary: tilesHtml + table,
            raw: JSON.stringify(analysis, null, 2),
            workspaceItem: { type: 'restriction-map', name: `Map (${enzymes.length} enzymes)`, data: analysis },
            downloads: [
                { label: 'TSV',  filename: 'restriction-map.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'restriction-map.json', text: JSON.stringify(analysis, null, 2), mime: 'application/json' },
            ],
        });
    }
}

function loadMapExample() {
    document.getElementById('mapSequence').value = 'GAATTCAAGCTTATCGATCGAATTCCTGCAGGGATCCAAGCTTTCTAGATGCATGCCTGCAGGAATTC';
    document.getElementById('mapEnzymes').value = 'EcoRI, BamHI, HindIII, PstI';

    setTimeout(() => {
        document.getElementById('mapForm').dispatchEvent(new Event('submit'));
    }, 500);
}

// ============================================================================
// UTILITIES TAB
// ============================================================================

function handleCompatibleEnds(e) {
    e.preventDefault();

    const enzymesInput = document.getElementById('compatibleEnzymes').value.trim();

    if (!enzymesInput) {
        showAlert('Please enter enzyme names', 'warning');
        return;
    }

    const enzymes = enzymesInput.split(',').map(e => e.trim()).filter(e => e);

    if (enzymes.length < 2) {
        showAlert('Please enter at least 2 enzymes', 'warning');
        return;
    }

    fetch('/api/restriction/compatible_ends', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({enzymes: enzymes})
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            displayCompatibleResults(data.compatible_pairs);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

function displayCompatibleResults(pairs) {
    const resultsDiv = document.getElementById('compatibleResults');
    if (pairs.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted small mt-2">No compatible pairs found.</p>';
        return;
    }

    const tiles = [{ label: 'Compatible Pairs', value: pairs.length }];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    pairs.forEach(p => { rows += `<tr><td><strong>${p.enzyme1}</strong></td><td class="text-center">↔</td><td><strong>${p.enzyme2}</strong></td><td><code>${p.overhang_seq}</code></td></tr>`; });
    const table = `<div class="table-responsive"><table class="table table-sm">
        <thead class="table-light"><tr><th>Enzyme 1</th><th></th><th>Enzyme 2</th><th>Overhang</th></tr></thead>
        <tbody>${rows}</tbody></table></div>`;

    const tsv = ['enzyme1\tenzyme2\toverhang'];
    pairs.forEach(p => tsv.push([p.enzyme1, p.enzyme2, p.overhang_seq].join('\t')));

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('compatibleResults', {
            title: 'Compatible Ends',
            meta: `${pairs.length} compatible pair${pairs.length === 1 ? '' : 's'}`,
            summary: tilesHtml + table,
            raw: JSON.stringify(pairs, null, 2),
            workspaceItem: { type: 'restriction-compatible', name: `Compatible ends (${pairs.length} pairs)`, data: pairs },
            downloads: [
                { label: 'TSV',  filename: 'compatible.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'compatible.json', text: JSON.stringify(pairs, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        resultsDiv.innerHTML = tilesHtml + table;
    }
}

function loadCompatibleExample() {
    document.getElementById('compatibleEnzymes').value = 'BamHI, BglII, XbaI, SpeI, EcoRI, MfeI';

    setTimeout(() => {
        document.getElementById('compatibleEndsForm').dispatchEvent(new Event('submit'));
    }, 500);
}

function exportLastResults(format) {
    if (!lastAnalysisResults) {
        showAlert('No results to export. Run an analysis first.', 'warning');
        return;
    }

    exportResults(lastAnalysisResults, format);
}

function exportResults(results, format) {
    fetch('/api/restriction/export', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            results: results,
            format: format
        })
    })
    .then(response => {
        if (response.ok) {
            return response.blob();
        }
        throw new Error('Export failed');
    })
    .then(blob => {
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `restriction_analysis.${format}`;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);
        showAlert(`Results exported as ${format.toUpperCase()}`, 'success');
    })
    .catch(error => {
        showAlert(friendlyError(error, 'server'), 'danger');
    });
}

// showLoading, hideLoading, showAlert are provided globally by static/js/utils.js
