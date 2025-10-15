// Structure Analysis JavaScript
let currentStructureFile = null;

// Parse Tab
document.getElementById('structureForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const fileInput = document.getElementById('structureFile');

    if (!fileInput.files[0] && !window.structureInfo) {
        return;
    }

    if (fileInput.files[0]) {
        currentStructureFile = fileInput.files[0];
        const formData = new FormData();
        formData.append('file', fileInput.files[0]);

        showLoading('parseStructureBtn');

        fetch('/api/structure/parse', {
            method: 'POST',
            body: formData
        })
        .then(response => response.json())
        .then(data => {
            hideLoading('parseStructureBtn', '<i class="fas fa-cube me-2"></i>Parse');
            if (data.success) {
                displayStructureInfo(data.structure_info);
                window.structureInfo = data.structure_info;
                document.getElementById('advancedAnalysisBtn').disabled = false;
            } else {
                showAlert('Error: ' + data.error, 'danger');
            }
        })
        .catch(error => {
            hideLoading('parseStructureBtn', '<i class="fas fa-cube me-2"></i>Parse');
        });
    } else if (window.structureInfo) {
        displayStructureInfo(window.structureInfo);
    }
});

function displayStructureInfo(structureInfo) {
    const infoDiv = document.getElementById('structureInfo');
    let html = '<div class="row g-2">';
    html += `<div class="col-6"><div class="border rounded p-2 bg-light text-center"><h6 class="text-danger mb-1 small">Models</h6><h5 class="mb-0">${structureInfo.models}</h5></div></div>`;
    html += `<div class="col-6"><div class="border rounded p-2 bg-light text-center"><h6 class="text-danger mb-1 small">Chains</h6><h5 class="mb-0">${structureInfo.chains.length}</h5></div></div>`;
    html += `<div class="col-6"><div class="border rounded p-2 bg-light text-center"><h6 class="text-danger mb-1 small">Residues</h6><h5 class="mb-0">${structureInfo.residue_count}</h5></div></div>`;
    html += `<div class="col-6"><div class="border rounded p-2 bg-light text-center"><h6 class="text-danger mb-1 small">Atoms</h6><h5 class="mb-0">${structureInfo.atom_count}</h5></div></div>`;
    html += '</div>';
    infoDiv.innerHTML = html;

    const chainDiv = document.getElementById('chainDetails');
    if (structureInfo.chains.length > 0) {
        let chainHtml = '<div class="table-responsive"><table class="table table-hover table-sm"><thead><tr><th>Chain</th><th>Residues</th><th>Atoms</th><th>Type</th></tr></thead><tbody>';
        structureInfo.chains.forEach(chain => {
            const chainType = chain.residue_count > 50 ? 'Protein' : 'Small molecule/Ion';
            chainHtml += `<tr><td><span class="badge bg-danger">${chain.id}</span></td><td>${chain.residue_count}</td><td>${chain.atoms}</td><td><small>${chainType}</small></td></tr>`;
        });
        chainHtml += '</tbody></table></div>';
        chainDiv.innerHTML = chainHtml;
    }
}

function fetchFromPDB() {
    const pdbId = document.getElementById('pdbId').value.trim().toUpperCase();
    if (!pdbId || pdbId.length !== 4) return;

    const pdbUrl = `https://files.rcsb.org/download/${pdbId}.pdb`;
    fetch(pdbUrl)
        .then(response => {
            if (!response.ok) throw new Error(`PDB ${pdbId} not found`);
            return response.text();
        })
        .then(pdbContent => {
            const formData = new FormData();
            const blob = new Blob([pdbContent], { type: 'text/plain' });
            const file = new File([blob], `${pdbId}.pdb`, { type: 'text/plain' });
            formData.append('file', file);
            currentStructureFile = file;
            return fetch('/api/structure/parse', { method: 'POST', body: formData });
        })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                displayStructureInfo(data.structure_info);
                window.structureInfo = data.structure_info;
                document.getElementById('advancedAnalysisBtn').disabled = false;
            } else {
                showAlert('Error: ' + data.error, 'danger');
            }
        })
        .catch(error => {
            showAlert('Error: ' + error.message, 'danger');
        });
}

function loadExample() {
    const samplePDB = `HEADER    INSULIN SAMPLE\nATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N\nEND`;
    const formData = new FormData();
    const blob = new Blob([samplePDB], { type: 'text/plain' });
    const file = new File([blob], 'sample.pdb', { type: 'text/plain' });
    formData.append('file', file);
    currentStructureFile = file;

    fetch('/api/structure/parse', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                displayStructureInfo(data.structure_info);
                window.structureInfo = data.structure_info;
                document.getElementById('advancedAnalysisBtn').disabled = false;
            }
        });
}

// Superimpose Tab
document.getElementById('superimposeForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const file1 = document.getElementById('superFile1').files[0];
    const file2 = document.getElementById('superFile2').files[0];

    if (!file1 || !file2) return;

    const formData = new FormData();
    formData.append('file1', file1);
    formData.append('file2', file2);

    showLoading('superimposeBtn');

    fetch('/api/structure/superimpose', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('superimposeBtn', '<i class="fas fa-project-diagram me-2"></i>Align Structures');
        if (data.success) {
            displaySuperimposeResults(data);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('superimposeBtn', '<i class="fas fa-project-diagram me-2"></i>Align Structures');
    });
});

function displaySuperimposeResults(data) {
    const resultsDiv = document.getElementById('superimposeResults');
    let html = '<div class="row g-3">';
    html += `<div class="col-md-4"><div class="border rounded p-3 bg-light text-center"><h6 class="text-info mb-1 small">RMSD</h6><h4 class="mb-0">${data.rmsd.toFixed(3)} Å</h4></div></div>`;
    html += `<div class="col-md-4"><div class="border rounded p-3 bg-light text-center"><h6 class="text-success mb-1 small">Atoms Aligned</h6><h4 class="mb-0">${data.atoms_aligned}</h4></div></div>`;
    html += `<div class="col-md-4"><div class="border rounded p-3 bg-light text-center"><h6 class="text-warning mb-1 small">Quality</h6><h4 class="mb-0">${data.rmsd < 2 ? 'Excellent' : data.rmsd < 5 ? 'Good' : 'Fair'}</h4></div></div>`;
    html += '</div>';
    resultsDiv.innerHTML = html;
}

// Geometry Tab
document.getElementById('geometryForm').addEventListener('submit', function(e) {
    e.preventDefault();

    if (!currentStructureFile) return;

    const formData = new FormData();
    formData.append('file', currentStructureFile);
    formData.append('chain_id', document.getElementById('geomChainId').value);

    showLoading('geometryBtn');

    fetch('/api/structure/geometry', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('geometryBtn', '<i class="fas fa-ruler me-2"></i>Calculate Geometry');
        if (data.success) {
            displayGeometry(data.geometry);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('geometryBtn', '<i class="fas fa-ruler me-2"></i>Calculate Geometry');
    });
});

function displayGeometry(geometry) {
    const resultsDiv = document.getElementById('geometryResults');
    let html = '<div class="accordion" id="geometryAccordion">';

    html += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button" data-bs-toggle="collapse" data-bs-target="#distances">Distances (first 10)</button></h2>';
    html += '<div id="distances" class="accordion-collapse collapse show" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    html += '<table class="table table-sm"><thead><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead><tbody>';
    geometry.distances.slice(0, 10).forEach(d => {
        html += `<tr><td>${d.residue1}</td><td>${d.residue2}</td><td>${d.distance}</td></tr>`;
    });
    html += '</tbody></table></div></div></div>';

    html += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#angles">Angles (first 10)</button></h2>';
    html += '<div id="angles" class="accordion-collapse collapse" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    html += '<table class="table table-sm"><thead><tr><th>Residues</th><th>Angle (°)</th></tr></thead><tbody>';
    geometry.angles.slice(0, 10).forEach(a => {
        html += `<tr><td>${a.residues}</td><td>${a.angle_degrees}</td></tr>`;
    });
    html += '</tbody></table></div></div></div>';

    html += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#dihedrals">Dihedrals</button></h2>';
    html += '<div id="dihedrals" class="accordion-collapse collapse" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    html += '<table class="table table-sm"><thead><tr><th>Residues</th><th>Dihedral (°)</th></tr></thead><tbody>';
    geometry.dihedrals.forEach(d => {
        html += `<tr><td>${d.residues}</td><td>${d.dihedral_degrees}</td></tr>`;
    });
    html += '</tbody></table></div></div></div>';

    html += '</div>';
    resultsDiv.innerHTML = html;
}

// Quality Tab
document.getElementById('qualityForm').addEventListener('submit', function(e) {
    e.preventDefault();

    if (!currentStructureFile) return;

    const formData = new FormData();
    formData.append('file', currentStructureFile);

    showLoading('qualityBtn');

    fetch('/api/structure/quality', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('qualityBtn', '<i class="fas fa-check-circle me-2"></i>Analyze Quality');
        if (data.success) {
            displayQuality(data.quality);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('qualityBtn', '<i class="fas fa-check-circle me-2"></i>Analyze Quality');
    });
});

function displayQuality(quality) {
    const resultsDiv = document.getElementById('qualityResults');
    let html = '<div class="row g-3">';

    if (quality.bfactor_stats.mean) {
        html += '<div class="col-md-6"><h6 class="small">B-factor Statistics</h6><table class="table table-sm">';
        html += `<tr><td>Mean</td><td>${quality.bfactor_stats.mean}</td></tr>`;
        html += `<tr><td>Median</td><td>${quality.bfactor_stats.median}</td></tr>`;
        html += `<tr><td>Range</td><td>${quality.bfactor_stats.min} - ${quality.bfactor_stats.max}</td></tr>`;
        html += '</table></div>';
    }

    if (quality.occupancy_stats.mean) {
        html += '<div class="col-md-6"><h6 class="small">Occupancy Statistics</h6><table class="table table-sm">';
        html += `<tr><td>Mean</td><td>${quality.occupancy_stats.mean}</td></tr>`;
        html += `<tr><td>Range</td><td>${quality.occupancy_stats.min} - ${quality.occupancy_stats.max}</td></tr>`;
        html += '</table></div>';
    }

    html += '</div>';

    if (quality.residue_completeness.length > 0) {
        html += '<h6 class="mt-3 small">Incomplete Residues</h6>';
        html += '<div class="alert alert-warning small"><strong>Warning:</strong> ' + quality.residue_completeness.length + ' residues have missing backbone atoms</div>';
    }

    resultsDiv.innerHTML = html;
}

// Contacts Tab
document.getElementById('contactsForm').addEventListener('submit', function(e) {
    e.preventDefault();

    if (!currentStructureFile) return;

    const formData = new FormData();
    formData.append('file', currentStructureFile);
    formData.append('chain_id', document.getElementById('contactChainId').value);
    formData.append('cutoff', document.getElementById('contactCutoff').value);

    showLoading('contactsBtn');

    fetch('/api/structure/contacts', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('contactsBtn', '<i class="fas fa-network-wired me-2"></i>Calculate Contacts');
        if (data.success) {
            displayContacts(data);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('contactsBtn', '<i class="fas fa-network-wired me-2"></i>Calculate Contacts');
    });
});

function displayContacts(data) {
    const resultsDiv = document.getElementById('contactsResults');
    let html = `<p class="small">Found ${data.contact_count} contacts within ${data.cutoff}Å cutoff (showing first 50)</p>`;
    html += '<table class="table table-sm table-hover"><thead><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead><tbody>';
    data.contacts.forEach(contact => {
        html += `<tr><td>${contact.residue1}</td><td>${contact.residue2}</td><td>${contact.distance}</td></tr>`;
    });
    html += '</tbody></table>';
    resultsDiv.innerHTML = html;
}

// Interactions Tab
document.getElementById('interactionsForm').addEventListener('submit', function(e) {
    e.preventDefault();

    if (!currentStructureFile) return;

    const formData = new FormData();
    formData.append('file', currentStructureFile);

    showLoading('interactionsBtn');

    fetch('/api/structure/hbonds', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('interactionsBtn', '<i class="fas fa-link me-2"></i>Detect Interactions');
        if (data.success) {
            displayInteractions(data);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('interactionsBtn', '<i class="fas fa-link me-2"></i>Detect Interactions');
    });
});

function displayInteractions(data) {
    const resultsDiv = document.getElementById('interactionsResults');
    let html = '<div class="accordion" id="interactionsAccordion">';

    html += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button" data-bs-toggle="collapse" data-bs-target="#hbonds">';
    html += `Hydrogen Bonds (${data.hbond_count} found, showing first 50)</button></h2>`;
    html += '<div id="hbonds" class="accordion-collapse collapse show" data-bs-parent="#interactionsAccordion"><div class="accordion-body">';
    html += '<table class="table table-sm"><thead><tr><th>Donor</th><th>Acceptor</th><th>Distance (Å)</th></tr></thead><tbody>';
    data.hbonds.forEach(hb => {
        html += `<tr><td>${hb.donor}</td><td>${hb.acceptor}</td><td>${hb.distance}</td></tr>`;
    });
    html += '</tbody></table></div></div></div>';

    html += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#saltbridges">';
    html += `Salt Bridges (${data.salt_bridge_count} found)</button></h2>`;
    html += '<div id="saltbridges" class="accordion-collapse collapse" data-bs-parent="#interactionsAccordion"><div class="accordion-body">';
    html += '<table class="table table-sm"><thead><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead><tbody>';
    data.salt_bridges.forEach(sb => {
        html += `<tr><td>${sb.residue1}</td><td>${sb.residue2}</td><td>${sb.distance}</td></tr>`;
    });
    html += '</tbody></table></div></div></div>';

    html += '</div>';
    resultsDiv.innerHTML = html;
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

// DSSP Tab
document.getElementById('dsspForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('dsspFile');
    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);

    showLoading('dsspBtn');
    fetch('/api/structure/secondary_structure', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('dsspBtn', '<i class="fas fa-dna me-2"></i>Analyze Secondary Structure');
        if (data.success) {
            displayDsspResults(data.result);
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('dsspBtn', '<i class="fas fa-dna me-2"></i>Analyze Secondary Structure');
        showAlert('Request failed');
    });
});

function displayDsspResults(result) {
    const resultsDiv = document.getElementById('dsspResults');
    let html = '<div class="card border-0 shadow-sm mb-3"><div class="card-body">';

    if (result.dssp_available) {
        const ss = result.secondary_structure;
        html += '<h6 class="mb-2"><i class="fas fa-chart-pie"></i> Secondary Structure Distribution</h6><div class="row g-2">';
        for (const [key, count] of Object.entries(ss.counts)) {
            const name = ss.mapping[key] || key;
            html += `<div class="col-md-4"><div class="border rounded p-2"><strong>${name} (${key}):</strong> ${count}</div></div>`;
        }
        html += '</div>';

        if (result.residue_details.length > 0) {
            html += '<h6 class="mt-3 mb-2"><i class="fas fa-list"></i> Residue Details (first 50)</h6>';
            html += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>SS</th><th>Type</th><th>Accessibility</th><th>Phi</th><th>Psi</th></tr></thead><tbody>';
            result.residue_details.forEach(r => {
                html += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.ss}</td><td class="small">${r.ss_name}</td><td>${r.accessibility}</td><td>${r.phi !== null ? r.phi + '°' : '-'}</td><td>${r.psi !== null ? r.psi + '°' : '-'}</td></tr>`;
            });
            html += '</tbody></table></div>';
        }
    } else {
        html += '<p class="text-warning">DSSP not available or failed. ' + (result.dssp_error || '') + '</p>';
    }
    html += '</div></div>';
    resultsDiv.innerHTML = html;
}

// Ramachandran Tab
document.getElementById('ramachandranForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('ramaFile');
    const chainId = document.getElementById('ramaChainId').value;
    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('chain_id', chainId);

    showLoading('ramaBtn');
    fetch('/api/structure/ramachandran', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('ramaBtn', '<i class="fas fa-chart-scatter me-2"></i>Calculate Phi/Psi Angles');
        if (data.success) {
            displayRamaResults(data);
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('ramaBtn', '<i class="fas fa-chart-scatter me-2"></i>Calculate Phi/Psi Angles');
        showAlert('Request failed');
    });
});

function displayRamaResults(data) {
    const resultsDiv = document.getElementById('ramaResults');
    let html = '<div class="card border-0 shadow-sm mb-3"><div class="card-body">';

    html += `<h6 class="mb-2"><i class="fas fa-chart-pie"></i> Classification (${data.total_residues} residues)</h6>`;
    html += '<div class="row g-2 mb-3">';
    html += `<div class="col-md-4"><div class="border rounded p-2 bg-success text-white text-center"><strong>Favored:</strong> ${data.classification.favored}</div></div>`;
    html += `<div class="col-md-4"><div class="border rounded p-2 bg-warning text-center"><strong>Allowed:</strong> ${data.classification.allowed}</div></div>`;
    html += `<div class="col-md-4"><div class="border rounded p-2 bg-danger text-white text-center"><strong>Outliers:</strong> ${data.classification.outliers}</div></div>`;
    html += '</div>';

    html += '<h6 class="mb-2"><i class="fas fa-table"></i> Phi/Psi Angles (first 100)</h6>';
    html += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Residue</th><th>Phi (°)</th><th>Psi (°)</th></tr></thead><tbody>';
    data.phi_psi_data.forEach(r => {
        html += `<tr><td>${r.residue}</td><td>${r.phi}</td><td>${r.psi}</td></tr>`;
    });
    html += '</tbody></table></div>';

    html += '</div></div>';
    resultsDiv.innerHTML = html;
}

// SASA Tab
document.getElementById('sasaForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('sasaFile');
    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);

    showLoading('sasaBtn');
    fetch('/api/structure/sasa', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('sasaBtn', '<i class="fas fa-water me-2"></i>Calculate SASA');
        if (data.success) {
            displaySasaResults(data);
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('sasaBtn', '<i class="fas fa-water me-2"></i>Calculate SASA');
        showAlert('Request failed');
    });
});

function displaySasaResults(data) {
    const resultsDiv = document.getElementById('sasaResults');
    let html = '<div class="card border-0 shadow-sm mb-3"><div class="card-body">';

    if (data.total_sasa !== undefined) {
        html += `<h6 class="mb-2"><i class="fas fa-water"></i> Total SASA: ${data.total_sasa} Ų</h6>`;
        html += '<h6 class="mt-3 mb-2">Chain SASA:</h6><div class="row g-2 mb-3">';
        for (const [chain, sasa] of Object.entries(data.chain_sasa)) {
            html += `<div class="col-md-3"><div class="border rounded p-2"><strong>Chain ${chain}:</strong> ${sasa} Ų</div></div>`;
        }
        html += '</div>';

        if (data.residue_sasa && data.residue_sasa.length > 0) {
            html += '<h6 class="mb-2">Residue SASA (first 50):</h6>';
            html += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>SASA (Ų)</th></tr></thead><tbody>';
            data.residue_sasa.forEach(r => {
                html += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.sasa}</td></tr>`;
            });
            html += '</tbody></table></div>';
        }
    } else if (data.method) {
        html += `<p class="text-info">Using ${data.method}</p>`;
        html += '<h6 class="mb-2">Chain Accessibility:</h6><div class="row g-2 mb-3">';
        for (const [chain, acc] of Object.entries(data.chain_accessibility)) {
            html += `<div class="col-md-3"><div class="border rounded p-2"><strong>Chain ${chain}:</strong> ${acc}</div></div>`;
        }
        html += '</div>';

        if (data.residue_accessibility) {
            html += '<h6 class="mb-2">Relative Accessibility (first 50):</h6>';
            html += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>Accessibility</th></tr></thead><tbody>';
            data.residue_accessibility.forEach(r => {
                html += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.relative_accessibility}</td></tr>`;
            });
            html += '</tbody></table></div>';
        }
    }

    html += '</div></div>';
    resultsDiv.innerHTML = html;
}

// Selection Tab
document.getElementById('selectionForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('selectionFile');
    const selectionType = document.getElementById('selectionType').value;
    const selectionValue = document.getElementById('selectionValue').value;
    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('selection_type', selectionType);
    formData.append('selection_value', selectionValue);

    showLoading('selectionBtn');
    fetch('/api/structure/extract', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('selectionBtn', '<i class="fas fa-crop me-2"></i>Extract Selection');
        if (data.success) {
            displaySelectionResults(data.extraction);
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('selectionBtn', '<i class="fas fa-crop me-2"></i>Extract Selection');
        showAlert('Request failed');
    });
});

function displaySelectionResults(extraction) {
    const resultsDiv = document.getElementById('selectionResults');
    let html = '<div class="card border-0 shadow-sm mb-3"><div class="card-body">';

    html += `<h6 class="mb-2"><i class="fas fa-info-circle"></i> Selection Type: ${extraction.selection_type}</h6>`;
    html += `<p><strong>Selection Value:</strong> ${extraction.selection_value}</p>`;
    html += `<p><strong>Extracted Count:</strong> ${extraction.extracted_count}</p>`;

    if (extraction.details.length > 0) {
        html += '<h6 class="mt-3 mb-2">Details:</h6>';
        html += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr>';

        const firstItem = extraction.details[0];
        for (const key of Object.keys(firstItem)) {
            html += `<th>${key}</th>`;
        }
        html += '</tr></thead><tbody>';

        extraction.details.forEach(item => {
            html += '<tr>';
            for (const value of Object.values(item)) {
                html += `<td>${value}</td>`;
            }
            html += '</tr>';
        });
        html += '</tbody></table></div>';
    }

    html += '</div></div>';
    resultsDiv.innerHTML = html;
}

// Example button functions
function loadDsspExample() {
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to analyze DSSP', 'info');
}

function loadRamaExample() {
    document.getElementById('ramaChainId').value = 'A';
    showAlert('Click "Parse" tab Example button first to load a structure, then return here for Ramachandran analysis', 'info');
}

function loadSasaExample() {
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to calculate SASA', 'info');
}

function loadSelectionExample() {
    document.getElementById('selectionType').value = 'chain';
    document.getElementById('selectionValue').value = 'A';
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to extract chain A', 'info');
}

function loadSuperimposeExample() {
    showAlert('Upload two PDB files to superimpose. Tip: Use test_structure.pdb from project root for both files', 'info');
}

function loadGeometryExample() {
    document.getElementById('geomChainId').value = 'A';
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to calculate geometry for chain A', 'info');
}

function loadQualityExample() {
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to analyze quality metrics', 'info');
}

function loadContactsExample() {
    document.getElementById('contactChainId').value = 'A';
    document.getElementById('contactCutoff').value = '8.0';
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to calculate contacts for chain A', 'info');
}

function loadInteractionsExample() {
    showAlert('Click "Parse" tab Example button first to load a structure, then return here to detect molecular interactions', 'info');
}
