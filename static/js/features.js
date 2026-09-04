// Sequence Features JavaScript

// Example Loading Functions
function loadORFExample() {
    document.getElementById('orfSequence').value = "ATGGCTAGCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAA";
    document.getElementById('minLength').value = "50";
    document.getElementById('strand').value = "both";
}

function loadCreateExample() {
    document.getElementById('createSequence').value = "ATGGCTAGCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT";
    document.getElementById('featureType').value = "CDS";
    document.getElementById('createStart').value = "0";
    document.getElementById('createEnd').value = "60";
    document.getElementById('createStrand').value = "1";
}

function loadExtractExample() {
    document.getElementById('extractSequence').value = "ATGGCTAGCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT";
    document.getElementById('extractStart').value = "0";
    document.getElementById('extractEnd').value = "30";
    document.getElementById('extractStrand').value = "1";
    document.getElementById('translateCheck').checked = true;
}

function loadCompoundExample() {
    document.getElementById('compoundSequence').value = "ATGGCTAGCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGA";
    document.getElementById('compoundLocations').value = "0,20,1\n30,50,1\n60,80,1";
}

function loadAnnotateExample() {
    document.getElementById('annotateSequence').value = "ATGGCTAGCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT";
    document.getElementById('annotateId').value = "example_seq";
    document.getElementById('annotateDesc').value = "Example annotated sequence with multiple features";

    // Clear existing features except first one
    const featuresList = document.getElementById('annotateFeaturesList');
    const features = featuresList.querySelectorAll('.feature-entry');
    features.forEach((feat, idx) => {
        if (idx > 0) feat.remove();
    });

    // Fill first feature
    const firstFeature = features[0];
    firstFeature.querySelector('.feat-type').value = "gene";
    firstFeature.querySelector('.feat-start').value = "0";
    firstFeature.querySelector('.feat-end').value = "60";
    firstFeature.querySelector('.feat-strand').value = "1";
    firstFeature.querySelector('.feat-label').value = "example_gene";
}

// ORF Finder Tab
document.getElementById('orfForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const sequence = document.getElementById('orfSequence').value.trim();
    const minLength = document.getElementById('minLength').value;
    const strand = document.getElementById('strand').value;

    if (!sequence) {
        showAlert('Please enter a DNA sequence', 'warning');
        return;
    }

    showLoading('orfBtn');

    fetch('/api/features/orf_find', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            min_length: parseInt(minLength),
            strand: strand
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('orfBtn', '<i class="fas fa-search me-2"></i>Find ORFs');
        if (data.success) {
            displayORFs(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('orfBtn', '<i class="fas fa-search me-2"></i>Find ORFs');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayORFs(data) {
    window.currentORFs = data.orfs;

    // Empty-state shortcut
    if (!data.orfs.length) {
        ResultsCard.mount('orfResults', {
            title: 'ORF Finder',
            meta: '<span class="meta">No ORFs found</span>',
            summary: '<div class="alert alert-warning mb-0">No ORFs found with minimum length ' +
                     data.min_length + ' bp. Try lowering the minimum length.</div>',
        });
        return;
    }

    // --- Summary: headline stats ---
    const longest = data.orfs.reduce((m, o) => Math.max(m, o.length), 0);
    const avg = Math.round(data.orfs.reduce((s, o) => s + o.length, 0) / data.orfs.length);
    const summary =
        '<div class="rc-stats">' +
        '<div class="rc-stat"><div class="rc-stat-label">ORFs found</div>' +
            '<div class="rc-stat-value">' + data.orf_count + '</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Longest</div>' +
            '<div class="rc-stat-value">' + longest + '</div><div class="rc-stat-sub">bp</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Average</div>' +
            '<div class="rc-stat-value">' + avg + '</div><div class="rc-stat-sub">bp</div></div>' +
        '<div class="rc-stat"><div class="rc-stat-label">Min length</div>' +
            '<div class="rc-stat-value">' + data.min_length + '</div><div class="rc-stat-sub">bp</div></div>' +
        '</div>';

    // --- Details: table with View action ---
    let details = '<div class="table-responsive"><table class="table table-sm table-hover mb-0"><thead><tr>';
    details += '<th>Start</th><th>End</th><th>Length</th><th>Frame</th><th>Strand</th><th>Protein length</th><th></th>';
    details += '</tr></thead><tbody>';
    data.orfs.forEach((orf, index) => {
        details += '<tr>' +
            '<td>' + orf.start + '</td>' +
            '<td>' + orf.end + '</td>' +
            '<td>' + orf.length + ' bp</td>' +
            '<td>' + orf.frame + '</td>' +
            '<td>' + escapeHtml(orf.strand) + '</td>' +
            '<td>' + orf.protein_length + ' aa</td>' +
            '<td><button class="btn-app-sm btn-app-secondary" data-action="showORFDetails" data-action-args="[' + index + ']">View</button></td>' +
        '</tr>';
    });
    details += '</tbody></table></div>';

    // --- Raw: FASTA of every ORF's protein sequence ---
    const fasta = data.orfs.map((o, i) =>
        '>ORF_' + (i + 1) + ' ' + o.start + '-' + o.end + ' frame=' + o.frame + ' strand=' + o.strand +
        '\n' + (o.protein || '')
    ).join('\n') + '\n';

    ResultsCard.mount('orfResults', {
        title: 'ORF Finder',
        meta:  '<span class="meta">' + data.orf_count + ' ORFs · min ' + data.min_length + ' bp</span>',
        summary:  summary,
        details:  details,
        raw:      fasta,
        copyText: fasta,
        downloads: [
            { label: 'Proteins (FASTA)', filename: 'orfs_protein.fasta', text: fasta, mime: 'text/plain' },
            { label: 'DNA (FASTA)', filename: 'orfs_dna.fasta',
              text: data.orfs.map((o, i) =>
                  '>ORF_' + (i + 1) + ' ' + o.start + '-' + o.end + ' frame=' + o.frame +
                  '\n' + (o.sequence || '')).join('\n') + '\n',
              mime: 'text/plain' },
            { label: 'Table (JSON)', filename: 'orfs.json',
              text: JSON.stringify(data.orfs, null, 2),
              mime: 'application/json' },
        ],
        workspaceItem: {
            type: 'sequence',
            name: 'ORFs · ' + data.orf_count + ' hits',
            data: fasta,
            meta: { source: 'orf_finder', min_length: data.min_length },
        },
    });
}

function showORFDetails(index) {
    const orf = window.currentORFs[index];
    const modal = `<div class="modal fade" id="orfModal" tabindex="-1">
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header">
                    <h5 class="modal-title">ORF Details (${orf.start}-${orf.end})</h5>
                    <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
                </div>
                <div class="modal-body">
                    <p><strong>Length:</strong> ${orf.length} bp</p>
                    <p><strong>Frame:</strong> ${orf.frame} | <strong>Strand:</strong> ${escapeHtml(orf.strand)}</p>
                    <p><strong>DNA Sequence:</strong></p>
                    <textarea class="form-control sequence-display mb-2" rows="3" readonly aria-label="Sequence">${escapeHtml(orf.sequence)}</textarea>
                    <p><strong>Protein Sequence (${orf.protein_length} aa):</strong></p>
                    <textarea class="form-control sequence-display" rows="3" readonly aria-label="Sequence">${escapeHtml(orf.protein)}</textarea>
                </div>
            </div>
        </div>
    </div>`;

    const existing = document.getElementById('orfModal');
    if (existing) existing.remove();
    document.body.insertAdjacentHTML('beforeend', modal);
    new bootstrap.Modal(document.getElementById('orfModal')).show();
}

// Create Feature Tab
document.getElementById('createFeatureForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const sequence = document.getElementById('createSequence').value.trim();
    const featureType = document.getElementById('featureType').value;
    const start = parseInt(document.getElementById('createStart').value);
    const end = parseInt(document.getElementById('createEnd').value);
    const strand = parseInt(document.getElementById('createStrand').value);

    if (!sequence) {
        showAlert('Please enter a sequence', 'warning');
        return;
    }

    showLoading('createBtn');

    fetch('/api/features/create', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            feature_type: featureType,
            start: start,
            end: end,
            strand: strand
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('createBtn', '<i class="fas fa-plus me-2"></i>Create Feature');
        if (data.success) {
            displayCreatedFeature(data.feature);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('createBtn', '<i class="fas fa-plus me-2"></i>Create Feature');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayCreatedFeature(feature) {
    const tiles = [
        { label: 'Type',   value: escapeHtml(feature.type) },
        { label: 'Length', value: feature.length + ' bp' },
        { label: 'Strand', value: feature.strand },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const summary = tilesHtml +
        `<p class="mb-1"><strong>Location:</strong> ${feature.start}..${feature.end} (${feature.strand} strand)</p>` +
        `<p class="mb-1"><strong>Sequence:</strong></p>` +
        `<textarea class="form-control sequence-display" rows="3" readonly aria-label="Sequence">${escapeHtml(feature.sequence)}</textarea>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('createResults', {
            title: 'Created Feature',
            meta: `${feature.type} · ${feature.length} bp`,
            summary: summary,
            raw: feature.sequence,
            copyText: feature.sequence,
            workspaceItem: { type: 'sequence', name: `${feature.type} ${feature.start}..${feature.end}`, data: feature.sequence, meta: { kind: 'dna' } },
            downloads: [
                { label: 'FASTA', filename: 'feature.fasta', text: `>${feature.type}_${feature.start}_${feature.end}\n${feature.sequence}`, mime: 'text/plain' },
                { label: 'JSON',  filename: 'feature.json',  text: JSON.stringify(feature, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('createResults').innerHTML = summary;
    }
}

// Parse GenBank Tab
document.getElementById('parseForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const fileInput = document.getElementById('parseFile');
    const format = document.getElementById('parseFormat').value;

    if (!fileInput.files[0]) {
        showAlert('Please select a file', 'warning');
        return;
    }

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('format', format);

    showLoading('parseBtn');

    fetch('/api/features/parse_genbank', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Features');
        if (data.success) {
            displayParsedFeatures(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Features');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayParsedFeatures(data) {
    const featureTypes = [...new Set(data.features.map(f => f.type))];
    const tiles = [
        { label: 'Features',     value: data.feature_count },
        { label: 'Feature Types', value: featureTypes.length },
        { label: 'Sequence',     value: data.sequence_length.toLocaleString() + ' bp' },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const header = `<p class="small mb-2"><strong>${escapeHtml(data.record_id)}</strong> — ${escapeHtml(data.record_description)}</p>`;

    let accordion = '<div class="accordion" id="featuresAccordion">';
    data.features.forEach((feat, index) => {
        const isFirst = index === 0;
        accordion += `<div class="accordion-item">
            <h2 class="accordion-header">
                <button class="accordion-button ${isFirst ? '' : 'collapsed'}" type="button" data-bs-toggle="collapse" data-bs-target="#feat${index}">
                    ${escapeHtml(feat.type)} (${feat.start}..${feat.end}, ${escapeHtml(feat.strand_symbol)})
                </button>
            </h2>
            <div id="feat${index}" class="accordion-collapse collapse ${isFirst ? 'show' : ''}" data-bs-parent="#featuresAccordion">
                <div class="accordion-body small">
                    <p><strong>Location:</strong> ${feat.start}..${feat.end} (${feat.length} bp)</p>
                    <p><strong>Strand:</strong> ${escapeHtml(feat.strand_symbol)}</p>`;
        if (Object.keys(feat.qualifiers).length > 0) {
            accordion += '<p><strong>Qualifiers:</strong></p><ul>';
            for (const [key, value] of Object.entries(feat.qualifiers)) {
                accordion += `<li><strong>${escapeHtml(key)}:</strong> ${escapeHtml(value)}</li>`;
            }
            accordion += '</ul>';
        }
        if (feat.sequence) {
            accordion += `<p><strong>Sequence (first 100 bp):</strong></p><code class="small">${escapeHtml(feat.sequence)}</code>`;
        }
        accordion += '</div></div></div>';
    });
    accordion += '</div>';

    const tsv = ['type\tstart\tend\tlength\tstrand'];
    data.features.forEach(f => tsv.push([f.type, f.start, f.end, f.length, f.strand_symbol].join('\t')));

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('parseResults', {
            title: 'Features Parse',
            meta: `${data.record_id} · ${data.feature_count} features`,
            summary: tilesHtml + header + accordion,
            raw: JSON.stringify(data.features, null, 2),
            workspaceItem: { type: 'features', name: `${data.record_id} features (${data.feature_count})`, data: data.features },
            downloads: [
                { label: 'TSV',  filename: 'features.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'features.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('parseResults').innerHTML = tilesHtml + header + accordion;
    }
}

// Extract Feature Tab
document.getElementById('extractForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const sequence = document.getElementById('extractSequence').value.trim();
    const start = parseInt(document.getElementById('extractStart').value);
    const end = parseInt(document.getElementById('extractEnd').value);
    const strand = parseInt(document.getElementById('extractStrand').value);
    const translate = document.getElementById('translateCheck').checked;

    if (!sequence) {
        showAlert('Please enter a sequence', 'warning');
        return;
    }

    showLoading('extractBtn');

    fetch('/api/features/extract', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            start: start,
            end: end,
            strand: strand,
            translate: translate
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('extractBtn', '<i class="fas fa-cut me-2"></i>Extract Feature');
        if (data.success) {
            displayExtractedFeature(data.feature);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('extractBtn', '<i class="fas fa-cut me-2"></i>Extract Feature');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayExtractedFeature(feature) {
    const tiles = [
        { label: 'Start',  value: feature.start },
        { label: 'End',    value: feature.end },
        { label: 'Length', value: feature.length + ' bp' },
        { label: 'Strand', value: feature.strand },
    ];
    if (feature.protein) tiles.push({ label: 'Protein', value: feature.protein_length + ' AA' });
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let summary = tilesHtml +
        '<p class="mb-1"><strong>Extracted Sequence:</strong></p>' +
        `<textarea class="form-control sequence-display mb-2" rows="3" readonly aria-label="Sequence">${escapeHtml(feature.sequence)}</textarea>`;
    if (feature.protein) {
        summary += `<p class="mb-1"><strong>Protein Translation (${feature.protein_length} aa):</strong></p>` +
            `<textarea class="form-control sequence-display" rows="2" readonly aria-label="Sequence">${escapeHtml(feature.protein)}</textarea>`;
    }

    const downloads = [
        { label: 'FASTA (DNA)', filename: 'feature.fasta', text: `>feature_${feature.start}_${feature.end}\n${feature.sequence}`, mime: 'text/plain' },
    ];
    if (feature.protein) {
        downloads.push({ label: 'FASTA (protein)', filename: 'feature-protein.fasta', text: `>feature_protein\n${feature.protein}`, mime: 'text/plain' });
    }
    downloads.push({ label: 'JSON', filename: 'feature.json', text: JSON.stringify(feature, null, 2), mime: 'application/json' });

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('extractResults', {
            title: 'Extracted Feature',
            meta: `${feature.start}..${feature.end} · ${feature.length} bp`,
            summary: summary,
            raw: feature.sequence,
            copyText: feature.sequence,
            workspaceItem: { type: 'sequence', name: `Feature ${feature.start}..${feature.end} (${feature.length} bp)`, data: feature.sequence, meta: { kind: 'dna' } },
            downloads: downloads,
        });
    } else {
        document.getElementById('extractResults').innerHTML = summary;
    }
}

// Compound Location Tab
document.getElementById('compoundForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const sequence = document.getElementById('compoundSequence').value.trim();
    const locationsText = document.getElementById('compoundLocations').value.trim();

    if (!sequence || !locationsText) {
        showAlert('Please enter sequence and locations', 'warning');
        return;
    }

    // Parse locations
    const locations = [];
    const lines = locationsText.split('\n');
    for (const line of lines) {
        const parts = line.trim().split(',');
        if (parts.length >= 2) {
            locations.push({
                start: parseInt(parts[0]),
                end: parseInt(parts[1]),
                strand: parts.length > 2 ? parseInt(parts[2]) : 1
            });
        }
    }

    if (locations.length === 0) {
        showAlert('No valid locations found', 'warning');
        return;
    }

    showLoading('compoundBtn');

    fetch('/api/features/compound_location', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            locations: locations
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('compoundBtn', '<i class="fas fa-layer-group me-2"></i>Create Compound Feature');
        if (data.success) {
            displayCompoundFeature(data.compound_feature);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('compoundBtn', '<i class="fas fa-layer-group me-2"></i>Create Compound Feature');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayCompoundFeature(feature) {
    const tiles = [
        { label: 'Segments',     value: feature.parts },
        { label: 'Total Length', value: feature.total_length + ' bp' },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let locs = '<p class="mb-1"><strong>Locations:</strong></p><ul>';
    feature.locations.forEach(loc => {
        locs += `<li>${loc.start}..${loc.end} (${loc.strand > 0 ? '+' : '-'})</li>`;
    });
    locs += '</ul>';
    const textarea = `<p class="mb-1"><strong>Joined Sequence:</strong></p><textarea class="form-control sequence-display" rows="3" readonly aria-label="Sequence">${escapeHtml(feature.sequence)}</textarea>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('compoundResults', {
            title: 'Compound Feature',
            meta: `${feature.parts} segments · ${feature.total_length} bp`,
            summary: tilesHtml + locs + textarea,
            raw: feature.sequence,
            copyText: feature.sequence,
            workspaceItem: { type: 'sequence', name: `Compound feature (${feature.total_length} bp)`, data: feature.sequence, meta: { kind: 'dna' } },
            downloads: [
                { label: 'FASTA', filename: 'compound.fasta', text: `>compound_${feature.parts}_parts\n${feature.sequence}`, mime: 'text/plain' },
                { label: 'JSON',  filename: 'compound.json',  text: JSON.stringify(feature, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('compoundResults').innerHTML = tilesHtml + locs + textarea;
    }
}

// Annotate Tab
document.getElementById('addFeatureBtn').addEventListener('click', function() {
    const featureEntry = `<div class="feature-entry p-2 border rounded mb-2">
        <div class="row g-2">
            <div class="col-md-3">
                <input type="text" class="form-control form-control-sm feat-type" placeholder="Type (e.g., CDS)" aria-label="Type">
            </div>
            <div class="col-md-2">
                <input type="number" class="form-control form-control-sm feat-start" placeholder="Start" aria-label="Start">
            </div>
            <div class="col-md-2">
                <input type="number" class="form-control form-control-sm feat-end" placeholder="End" aria-label="End">
            </div>
            <div class="col-md-2">
                <select class="form-select form-select-sm feat-strand" aria-label="Feat strand">
                    <option value="1">+</option>
                    <option value="-1">-</option>
                </select>
            </div>
            <div class="col-md-3">
                <input type="text" class="form-control form-control-sm feat-label" placeholder="Label/Gene" aria-label="Label/Gene">
            </div>
        </div>
    </div>`;
    document.getElementById('annotateFeaturesList').insertAdjacentHTML('beforeend', featureEntry);
});

document.getElementById('annotateForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const sequence = document.getElementById('annotateSequence').value.trim();
    const seqId = document.getElementById('annotateId').value.trim();
    const description = document.getElementById('annotateDesc').value.trim();

    if (!sequence) {
        showAlert('Please enter a sequence', 'warning');
        return;
    }

    // Collect features
    const features = [];
    document.querySelectorAll('.feature-entry').forEach(entry => {
        const type = entry.querySelector('.feat-type').value.trim();
        const start = entry.querySelector('.feat-start').value;
        const end = entry.querySelector('.feat-end').value;
        const strand = entry.querySelector('.feat-strand').value;
        const label = entry.querySelector('.feat-label').value.trim();

        if (type && start && end) {
            const qualifiers = {};
            if (label) qualifiers.label = label;

            features.push({
                type: type,
                start: parseInt(start),
                end: parseInt(end),
                strand: parseInt(strand),
                qualifiers: qualifiers
            });
        }
    });

    if (features.length === 0) {
        showAlert('Please add at least one feature', 'warning');
        return;
    }

    showLoading('annotateSubmitBtn');

    fetch('/api/features/annotate', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            seq_id: seqId,
            description: description,
            features: features
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('annotateSubmitBtn', '<i class="fas fa-file-export me-2"></i>Generate GenBank File');
        if (data.success) {
            displayAnnotatedSequence(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('annotateSubmitBtn', '<i class="fas fa-file-export me-2"></i>Generate GenBank File');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayAnnotatedSequence(data) {
    window.currentGenBank = data.genbank;

    const tiles = [
        { label: 'Features', value: data.feature_count },
        { label: 'ID',       value: escapeHtml(data.summary.id) },
        { label: 'Length',   value: data.summary.length + ' bp' },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const summary = tilesHtml +
        `<p class="mb-1"><strong>GenBank Output:</strong></p>` +
        `<textarea class="form-control" rows="10" readonly aria-label="Feature output">${escapeHtml(data.genbank)}</textarea>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('annotateResults', {
            title: 'Annotated Sequence',
            meta: `${data.summary.id} · ${data.feature_count} features`,
            summary: summary,
            raw: data.genbank,
            copyText: data.genbank,
            workspaceItem: { type: 'genbank', name: `${data.summary.id} (${data.feature_count} features)`, data: data.genbank },
            downloads: [
                { label: 'GenBank', filename: 'annotated_sequence.gb', text: data.genbank, mime: 'text/plain' },
                { label: 'JSON',    filename: 'annotated_sequence.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('annotateResults').innerHTML = summary;
    }
}

function downloadGenBank() {
    const blob = new Blob([window.currentGenBank], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'annotated_sequence.gb';
    a.click();
    URL.revokeObjectURL(url);
}

// Helper Functions: showLoading, hideLoading, showAlert are provided globally by static/js/utils.js
