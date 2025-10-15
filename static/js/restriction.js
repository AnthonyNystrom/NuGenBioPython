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
            infoDiv.innerHTML = `<small class="text-danger"><i class="fas fa-exclamation-triangle"></i> ${data.error}</small>`;
            showAlert('Error loading file: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        infoDiv.innerHTML = `<small class="text-danger"><i class="fas fa-exclamation-triangle"></i> Upload failed</small>`;
        showAlert('Upload error: ' + error.message, 'danger');
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
            showAlert('Failed to load enzymes: ' + data.error, 'danger');
            loadFallbackEnzymes();
        }
    })
    .catch(error => {
        console.error('Error loading enzymes:', error);
        showAlert('Error loading restriction enzymes: ' + error.message, 'danger');
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
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('basicAnalyzeBtn', '<i class="fas fa-cut me-2"></i>Analyze Restriction Sites');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

function displayBasicResults(analysis, sequenceLength) {
    displayBasicSummary(analysis, sequenceLength);
    displayBasicDetailedResults(analysis);
    displayBasicFragmentChart(analysis);
}

function displayBasicSummary(analysis, sequenceLength) {
    const summaryDiv = document.getElementById('basicSummary');

    let totalCuts = 0;
    let enzymesCutting = 0;

    Object.values(analysis).forEach(result => {
        if (!result.error && result.number_of_cuts > 0) {
            totalCuts += result.number_of_cuts;
            enzymesCutting++;
        }
    });

    let html = '<div class="row g-2">';
    html += `
        <div class="col-12">
            <div class="border rounded p-2 bg-light text-center">
                <h6 class="text-success mb-1 small">Sequence Length</h6>
                <h5 class="mb-0">${sequenceLength} bp</h5>
            </div>
        </div>
        <div class="col-12">
            <div class="border rounded p-2 bg-light text-center">
                <h6 class="text-success mb-1 small">Total Cut Sites</h6>
                <h5 class="mb-0">${totalCuts}</h5>
            </div>
        </div>
        <div class="col-12">
            <div class="border rounded p-2 bg-light text-center">
                <h6 class="text-success mb-1 small">Enzymes Cutting</h6>
                <h5 class="mb-0">${enzymesCutting}</h5>
            </div>
        </div>
    `;
    html += '</div>';

    summaryDiv.innerHTML = html;
}

function displayBasicDetailedResults(analysis) {
    const resultsDiv = document.getElementById('basicResults');

    let html = '<div class="table-responsive">';
    html += '<table class="table table-hover table-sm"><thead class="table-light">';
    html += '<tr><th>Enzyme</th><th>Site</th><th>Type</th><th>Cuts</th><th>Positions</th><th>Fragments</th></tr>';
    html += '</thead><tbody>';

    Object.entries(analysis).forEach(([enzyme, result]) => {
        if (result.error) {
            html += `<tr><td>${enzyme}</td><td colspan="5" class="text-danger small">${result.error}</td></tr>`;
        } else {
            const overhangBadge = result.is_blunt ?
                '<span class="badge bg-secondary">Blunt</span>' :
                `<span class="badge bg-info">${result.overhang_type}</span>`;

            const fragmentsStr = result.fragments.length > 0 ? result.fragments.join(', ') + ' bp' : 'No cuts';
            const positionsStr = result.cut_positions.length > 0 ? result.cut_positions.join(', ') : 'None';

            html += `
                <tr>
                    <td><strong>${enzyme}</strong></td>
                    <td><code class="small">${result.recognition_site}</code></td>
                    <td>${overhangBadge}</td>
                    <td><span class="badge bg-primary">${result.number_of_cuts}</span></td>
                    <td><small class="text-muted">${positionsStr}</small></td>
                    <td><small>${fragmentsStr}</small></td>
                </tr>
            `;

            // Show fragment sequences if available
            if (result.fragment_sequences && result.fragment_sequences.length > 0) {
                html += `
                    <tr>
                        <td colspan="6" class="bg-light">
                            <small><strong>Fragment Sequences:</strong></small><br>
                            <div class="sequence-display" style="max-height: 100px; overflow-y: auto;">
                                ${result.fragment_sequences.map((seq, i) =>
                                    `<small><strong>Fragment ${i+1}:</strong> ${seq}</small>`
                                ).join('<br>')}
                            </div>
                        </td>
                    </tr>
                `;
            }
        }
    });

    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function displayBasicFragmentChart(analysis) {
    const chartDiv = document.getElementById('basicFragmentChart');

    // Collect all fragments
    const allFragments = [];
    Object.entries(analysis).forEach(([enzyme, result]) => {
        if (!result.error && result.number_of_cuts > 0 && result.fragments.length > 0) {
            allFragments.push(...result.fragments);
        }
    });

    if (allFragments.length === 0) {
        chartDiv.innerHTML = '<p class="text-muted small">No fragments to display</p>';
        return;
    }

    // Create histogram bins
    const bins = [0, 100, 500, 1000, 2000, 5000, 10000, Infinity];
    const binLabels = ['<100', '100-500', '500-1K', '1K-2K', '2K-5K', '5K-10K', '>10K'];
    const counts = new Array(bins.length - 1).fill(0);

    allFragments.forEach(size => {
        for (let i = 0; i < bins.length - 1; i++) {
            if (size >= bins[i] && size < bins[i + 1]) {
                counts[i]++;
                break;
            }
        }
    });

    let html = '<h6 class="small">Fragment Size Distribution</h6>';
    html += '<div class="row">';

    const maxCount = Math.max(...counts);

    counts.forEach((count, i) => {
        const percentage = maxCount > 0 ? (count / maxCount) * 100 : 0;
        html += `
            <div class="col">
                <div class="text-center mb-2">
                    <div class="bg-primary" style="height: ${Math.max(percentage, 5)}px; width: 100%; margin-bottom: 5px;"></div>
                    <small><strong>${count}</strong><br>${binLabels[i]}</small>
                </div>
            </div>
        `;
    });

    html += '</div>';
    chartDiv.innerHTML = html;
}

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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('advancedAnalyzeBtn', '<i class="fas fa-filter me-2"></i>Run Advanced Analysis');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

function displayAdvancedResults(result) {
    const resultsDiv = document.getElementById('advancedResults');

    let html = '<div class="alert alert-info small">';
    html += `<strong>Analysis Summary:</strong><br>`;
    html += `Sequence length: ${result.sequence_length} bp<br>`;
    html += `Enzymes tested: ${result.total_enzymes_tested}<br>`;
    html += `Filters: ${result.filters_applied.join(', ')}<br>`;
    html += `<strong>Results: ${result.enzyme_count} enzymes match criteria</strong>`;
    html += '</div>';

    if (result.enzymes.length === 0) {
        html += '<p class="text-muted text-center">No enzymes match the specified criteria.</p>';
    } else {
        html += '<div class="table-responsive">';
        html += '<table class="table table-hover table-sm"><thead class="table-light">';
        html += '<tr><th>Enzyme</th><th>Site</th><th>Type</th><th>Cuts</th><th>Positions</th></tr>';
        html += '</thead><tbody>';

        result.enzymes.forEach(enzyme => {
            const overhangBadge = enzyme.is_blunt ?
                '<span class="badge bg-secondary">Blunt</span>' :
                `<span class="badge bg-info">${enzyme.overhang_type}</span>`;

            const positionsStr = enzyme.positions.join(', ');

            html += `
                <tr>
                    <td><strong>${enzyme.name}</strong></td>
                    <td><code class="small">${enzyme.site}</code></td>
                    <td>${overhangBadge}</td>
                    <td><span class="badge bg-primary">${enzyme.cuts}</span></td>
                    <td><small class="text-muted">${positionsStr}</small></td>
                </tr>
            `;
        });

        html += '</tbody></table></div>';
    }

    resultsDiv.innerHTML = html;
}

function loadAdvancedExample() {
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
            resultsDiv.innerHTML = `<p class="text-danger small">Error: ${data.error}</p>`;
        }
    })
    .catch(error => {
        resultsDiv.innerHTML = `<p class="text-danger small">Error: ${error.message}</p>`;
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
            showAlert('Error loading enzyme details: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        showAlert('Error: ' + error.message, 'danger');
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
            displayRestrictionMap(data.analysis, sequence.length);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('generateMapBtn', '<i class="fas fa-map me-2"></i>Generate Restriction Map');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

function displayRestrictionMap(analysis, sequenceLength) {
    const mapDiv = document.getElementById('mapVisualization');

    // Create linear restriction map
    let html = '<div class="restriction-map">';
    html += `<h6 class="small mb-3">Linear Restriction Map (${sequenceLength} bp)</h6>`;

    // Draw scale
    html += '<div style="position: relative; height: 100px; margin-bottom: 20px;">';

    // Main sequence line
    html += '<div style="position: absolute; top: 50px; left: 0; right: 0; height: 3px; background: #333;"></div>';

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

    html += '</div>';
    mapDiv.innerHTML = html;
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        showAlert('Network error: ' + error.message, 'danger');
    });
}

function displayCompatibleResults(pairs) {
    const resultsDiv = document.getElementById('compatibleResults');

    if (pairs.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted small mt-2">No compatible pairs found.</p>';
        return;
    }

    let html = '<div class="alert alert-success small mt-2">';
    html += `<strong>${pairs.length} compatible pair(s) found:</strong><br>`;
    pairs.forEach(pair => {
        html += `<i class="fas fa-link"></i> ${pair.enzyme1} ↔ ${pair.enzyme2} (overhang: ${pair.overhang_seq})<br>`;
    });
    html += '</div>';

    resultsDiv.innerHTML = html;
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
        showAlert('Export failed: ' + error.message, 'danger');
    });
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

function showLoading(buttonId) {
    const btn = document.getElementById(buttonId);
    if (btn) {
        btn.disabled = true;
        btn.innerHTML = '<span class="spinner-border spinner-border-sm me-2"></span>Processing...';
    }
}

function hideLoading(buttonId, originalText) {
    const btn = document.getElementById(buttonId);
    if (btn) {
        btn.disabled = false;
        btn.innerHTML = originalText;
    }
}

function showAlert(message, type) {
    // Create and show Bootstrap alert
    const alertDiv = document.createElement('div');
    alertDiv.className = `alert alert-${type} alert-dismissible fade show position-fixed top-0 start-50 translate-middle-x mt-3`;
    alertDiv.style.zIndex = '9999';
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `;
    document.body.appendChild(alertDiv);

    setTimeout(() => {
        alertDiv.remove();
    }, 5000);
}
