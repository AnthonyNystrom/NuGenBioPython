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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('orfBtn', '<i class="fas fa-search me-2"></i>Find ORFs');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayORFs(data) {
    const resultsDiv = document.getElementById('orfResults');
    let html = `<div class="alert alert-info small">Found ${data.orf_count} ORFs (showing top 50, minimum length: ${data.min_length} bp)</div>`;

    if (data.orfs.length === 0) {
        html += '<div class="alert alert-warning">No ORFs found. Try lowering the minimum length.</div>';
    } else {
        html += '<div class="table-responsive"><table class="table table-sm table-hover"><thead><tr>';
        html += '<th>Start</th><th>End</th><th>Length</th><th>Frame</th><th>Strand</th><th>Protein Length</th><th>Action</th>';
        html += '</tr></thead><tbody>';

        data.orfs.forEach((orf, index) => {
            html += `<tr>
                <td>${orf.start}</td>
                <td>${orf.end}</td>
                <td>${orf.length} bp</td>
                <td>${orf.frame}</td>
                <td>${orf.strand}</td>
                <td>${orf.protein_length} aa</td>
                <td><button class="btn btn-sm btn-outline-primary" onclick="showORFDetails(${index})">View</button></td>
            </tr>`;
        });

        html += '</tbody></table></div>';
    }

    resultsDiv.innerHTML = html;
    window.currentORFs = data.orfs;
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
                    <p><strong>Frame:</strong> ${orf.frame} | <strong>Strand:</strong> ${orf.strand}</p>
                    <p><strong>DNA Sequence:</strong></p>
                    <textarea class="form-control sequence-display mb-2" rows="3" readonly>${orf.sequence}</textarea>
                    <p><strong>Protein Sequence (${orf.protein_length} aa):</strong></p>
                    <textarea class="form-control sequence-display" rows="3" readonly>${orf.protein}</textarea>
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('createBtn', '<i class="fas fa-plus me-2"></i>Create Feature');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayCreatedFeature(feature) {
    const resultsDiv = document.getElementById('createResults');
    let html = '<div class="alert alert-success">Feature created successfully!</div>';
    html += '<div class="card"><div class="card-body">';
    html += `<p><strong>Type:</strong> ${feature.type}</p>`;
    html += `<p><strong>Location:</strong> ${feature.start}..${feature.end} (${feature.strand} strand)</p>`;
    html += `<p><strong>Length:</strong> ${feature.length} bp</p>`;
    html += '<p><strong>Extracted Sequence:</strong></p>';
    html += `<textarea class="form-control sequence-display" rows="3" readonly>${feature.sequence}</textarea>`;
    html += '</div></div>';
    resultsDiv.innerHTML = html;
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Features');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayParsedFeatures(data) {
    const resultsDiv = document.getElementById('parseResults');
    let html = `<div class="alert alert-info small">Parsed ${data.feature_count} features from ${data.record_id}</div>`;
    html += `<p class="small"><strong>Description:</strong> ${data.record_description}</p>`;
    html += `<p class="small"><strong>Sequence Length:</strong> ${data.sequence_length} bp</p>`;

    html += '<div class="accordion" id="featuresAccordion">';
    data.features.forEach((feat, index) => {
        const isFirst = index === 0;
        html += `<div class="accordion-item">
            <h2 class="accordion-header">
                <button class="accordion-button ${isFirst ? '' : 'collapsed'}" type="button" data-bs-toggle="collapse" data-bs-target="#feat${index}">
                    ${feat.type} (${feat.start}..${feat.end}, ${feat.strand_symbol})
                </button>
            </h2>
            <div id="feat${index}" class="accordion-collapse collapse ${isFirst ? 'show' : ''}" data-bs-parent="#featuresAccordion">
                <div class="accordion-body small">
                    <p><strong>Location:</strong> ${feat.start}..${feat.end} (${feat.length} bp)</p>
                    <p><strong>Strand:</strong> ${feat.strand_symbol}</p>`;

        if (Object.keys(feat.qualifiers).length > 0) {
            html += '<p><strong>Qualifiers:</strong></p><ul>';
            for (const [key, value] of Object.entries(feat.qualifiers)) {
                html += `<li><strong>${key}:</strong> ${value}</li>`;
            }
            html += '</ul>';
        }

        if (feat.sequence) {
            html += `<p><strong>Sequence (first 100 bp):</strong></p>`;
            html += `<code class="small">${feat.sequence}</code>`;
        }

        html += '</div></div></div>';
    });
    html += '</div>';

    resultsDiv.innerHTML = html;
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('extractBtn', '<i class="fas fa-cut me-2"></i>Extract Feature');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayExtractedFeature(feature) {
    const resultsDiv = document.getElementById('extractResults');
    let html = '<div class="alert alert-success">Feature extracted successfully!</div>';
    html += '<div class="card"><div class="card-body">';
    html += `<p><strong>Location:</strong> ${feature.start}..${feature.end} (${feature.strand} strand)</p>`;
    html += `<p><strong>Length:</strong> ${feature.length} bp</p>`;
    html += '<p><strong>Extracted Sequence:</strong></p>';
    html += `<textarea class="form-control sequence-display mb-2" rows="3" readonly>${feature.sequence}</textarea>`;

    if (feature.protein) {
        html += `<p><strong>Protein Translation (${feature.protein_length} aa):</strong></p>`;
        html += `<textarea class="form-control sequence-display" rows="2" readonly>${feature.protein}</textarea>`;
    }

    html += '</div></div>';
    resultsDiv.innerHTML = html;
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('compoundBtn', '<i class="fas fa-layer-group me-2"></i>Create Compound Feature');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayCompoundFeature(feature) {
    const resultsDiv = document.getElementById('compoundResults');
    let html = '<div class="alert alert-success">Compound feature created successfully!</div>';
    html += '<div class="card"><div class="card-body">';
    html += `<p><strong>Parts:</strong> ${feature.parts} segments</p>`;
    html += `<p><strong>Total Length:</strong> ${feature.total_length} bp</p>`;
    html += '<p><strong>Locations:</strong></p><ul>';
    feature.locations.forEach(loc => {
        html += `<li>${loc.start}..${loc.end} (${loc.strand > 0 ? '+' : '-'})</li>`;
    });
    html += '</ul>';
    html += '<p><strong>Joined Sequence:</strong></p>';
    html += `<textarea class="form-control sequence-display" rows="3" readonly>${feature.sequence}</textarea>`;
    html += '</div></div>';
    resultsDiv.innerHTML = html;
}

// Annotate Tab
document.getElementById('addFeatureBtn').addEventListener('click', function() {
    const featureEntry = `<div class="feature-entry p-2 border rounded mb-2">
        <div class="row g-2">
            <div class="col-md-3">
                <input type="text" class="form-control form-control-sm feat-type" placeholder="Type (e.g., CDS)">
            </div>
            <div class="col-md-2">
                <input type="number" class="form-control form-control-sm feat-start" placeholder="Start">
            </div>
            <div class="col-md-2">
                <input type="number" class="form-control form-control-sm feat-end" placeholder="End">
            </div>
            <div class="col-md-2">
                <select class="form-select form-select-sm feat-strand">
                    <option value="1">+</option>
                    <option value="-1">-</option>
                </select>
            </div>
            <div class="col-md-3">
                <input type="text" class="form-control form-control-sm feat-label" placeholder="Label/Gene">
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
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('annotateSubmitBtn', '<i class="fas fa-file-export me-2"></i>Generate GenBank File');
        showAlert('Error: ' + error.message, 'danger');
    });
});

function displayAnnotatedSequence(data) {
    const resultsDiv = document.getElementById('annotateResults');
    let html = `<div class="alert alert-success">Generated GenBank file with ${data.feature_count} features!</div>`;
    html += '<div class="card"><div class="card-body">';
    html += `<p><strong>ID:</strong> ${data.summary.id} | <strong>Length:</strong> ${data.summary.length} bp</p>`;
    html += '<p><strong>GenBank Output:</strong></p>';
    html += `<textarea class="form-control" rows="10" readonly>${data.genbank}</textarea>`;
    html += '<button class="btn btn-primary mt-2" onclick="downloadGenBank()">Download GenBank File</button>';
    html += '</div></div>';
    resultsDiv.innerHTML = html;
    window.currentGenBank = data.genbank;
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

// Helper Functions
function showLoading(btnId) {
    const btn = document.getElementById(btnId);
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner-border spinner-border-sm me-2"></span>Processing...';
}

function hideLoading(btnId, originalText) {
    const btn = document.getElementById(btnId);
    btn.disabled = false;
    btn.innerHTML = originalText;
}

function showAlert(message, type = 'danger') {
    const alertDiv = document.createElement('div');
    alertDiv.className = `alert alert-${type} alert-dismissible fade show`;
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    `;
    document.querySelector('.card-body').prepend(alertDiv);
    setTimeout(() => alertDiv.remove(), 5000);
}
