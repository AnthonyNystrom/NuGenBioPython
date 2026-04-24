// Structure Analysis JavaScript
let currentStructureFile = null;

// Residue-coloring helpers. Each display function can call `applyResidueColors`
// with a chain-prefixed map {"A:10": "#hex", ...}; the 3Dmol viewer recolors
// cartoons/sticks for those residues while leaving the rest in the base style.
function _structureResidueNum(s) {
    const m = String(s).match(/(\d+)/);
    return m ? parseInt(m[1], 10) : null;
}
function _structureSsColor(code) {
    // PyMol-style SS palette. Unknown codes fall back to neutral slate.
    const palette = { H: '#ef4444', G: '#f87171', I: '#dc2626',
                      E: '#f59e0b', B: '#fbbf24',
                      T: '#14b8a6', S: '#94a3b8' };
    return palette[code] || '#cbd5e1';
}
function _structureGradient(t) {
    // t in [0,1]. Buried (low) = slate-blue, exposed (high) = coral.
    t = Math.max(0, Math.min(1, t));
    const c1 = [30, 64, 175];    // #1e40af
    const c2 = [252, 165, 165];  // #fca5a5
    const mix = c1.map((v, i) => Math.round(v + (c2[i] - v) * t));
    return '#' + mix.map(v => v.toString(16).padStart(2, '0')).join('');
}
function applyResidueColors(resiMap) {
    if (!window.NuGenStructureViewer || !resiMap) return;
    const keys = Object.keys(resiMap);
    if (keys.length === 0) return;
    try { window.NuGenStructureViewer.colorResidues(resiMap); } catch (_) {}
}

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
                showAlert(friendlyError(data.error, 'pdb'), 'danger');
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
            // Render the structure in the 3D viewer first so the user
            // sees something immediately while the parse call runs.
            if (window.NuGenStructureViewer && window.NuGenStructureViewer.load) {
                window.NuGenStructureViewer.load(pdbContent, pdbId);
            }
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
                showAlert(friendlyError(data.error, 'pdb'), 'danger');
            }
        })
        .catch(error => {
            showAlert(friendlyError(error, 'pdb'), 'danger');
        });
}

function loadExample() {
    // The bundled sample (Crambin 1CRN, 46 residues, α-helix + β-sheet) is
    // served from /api/structure/sample. Parse it on the server so the
    // sidebar stats populate, and push the raw text into the 3D viewer so
    // the canvas actually shows the structure.
    fetch('/api/structure/sample')
        .then(function (r) {
            if (!r.ok) throw new Error('HTTP ' + r.status);
            return r.text();
        })
        .then(function (text) {
            if (window.NuGenStructureViewer && window.NuGenStructureViewer.load) {
                window.NuGenStructureViewer.load(text, 'Crambin (1CRN) · 46 residues');
            }
            const blob = new Blob([text], { type: 'text/plain' });
            const file = new File([blob], 'sample_protein.pdb', { type: 'text/plain' });
            currentStructureFile = file;
            const fd = new FormData();
            fd.append('file', file);
            return fetch('/api/structure/parse', { method: 'POST', body: fd });
        })
        .then(function (r) { return r.json(); })
        .then(function (data) {
            if (data && data.success) {
                displayStructureInfo(data.structure_info);
                window.structureInfo = data.structure_info;
                const btn = document.getElementById('advancedAnalysisBtn');
                if (btn) btn.disabled = false;
            }
        })
        .catch(function (e) {
            if (window.showAlert) window.showAlert('Failed to load sample: ' + e.message, 'warning');
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('superimposeBtn', '<i class="fas fa-project-diagram me-2"></i>Align Structures');
    });
});

function displaySuperimposeResults(data) {
    const quality = data.rmsd < 2 ? 'Excellent' : data.rmsd < 5 ? 'Good' : 'Fair';
    const tiles = [
        { label: 'RMSD',          value: data.rmsd.toFixed(3) + ' Å' },
        { label: 'Atoms Aligned', value: data.atoms_aligned },
        { label: 'Quality',       value: quality },
    ];
    const tilesHtml = '<div class="rc-stats">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('superimposeResults', {
            title: 'Structure Superposition',
            meta: `RMSD ${data.rmsd.toFixed(3)} Å · ${data.atoms_aligned} atoms`,
            summary: tilesHtml,
            raw: JSON.stringify(data, null, 2),
            workspaceItem: { type: 'superimpose', name: `Superimpose RMSD ${data.rmsd.toFixed(3)}`, data: data },
            downloads: [
                { label: 'JSON', filename: 'superimpose.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('superimposeResults').innerHTML = tilesHtml;
    }
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('geometryBtn', '<i class="fas fa-ruler me-2"></i>Calculate Geometry');
    });
});

function displayGeometry(geometry) {
    const tiles = [
        { label: 'Distances', value: geometry.distances.length },
        { label: 'Angles',    value: geometry.angles.length },
        { label: 'Dihedrals', value: geometry.dihedrals.length },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let accordion = '<div class="accordion" id="geometryAccordion">';
    accordion += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button" data-bs-toggle="collapse" data-bs-target="#distances">Distances (first 10)</button></h2>';
    accordion += '<div id="distances" class="accordion-collapse collapse show" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    accordion += '<table class="table table-sm"><thead><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead><tbody>';
    geometry.distances.slice(0, 10).forEach(d => {
        accordion += `<tr><td>${d.residue1}</td><td>${d.residue2}</td><td>${d.distance}</td></tr>`;
    });
    accordion += '</tbody></table></div></div></div>';

    accordion += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#angles">Angles (first 10)</button></h2>';
    accordion += '<div id="angles" class="accordion-collapse collapse" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    accordion += '<table class="table table-sm"><thead><tr><th>Residues</th><th>Angle (°)</th></tr></thead><tbody>';
    geometry.angles.slice(0, 10).forEach(a => {
        accordion += `<tr><td>${a.residues}</td><td>${a.angle_degrees}</td></tr>`;
    });
    accordion += '</tbody></table></div></div></div>';

    accordion += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#dihedrals">Dihedrals</button></h2>';
    accordion += '<div id="dihedrals" class="accordion-collapse collapse" data-bs-parent="#geometryAccordion"><div class="accordion-body">';
    accordion += '<table class="table table-sm"><thead><tr><th>Residues</th><th>Dihedral (°)</th></tr></thead><tbody>';
    geometry.dihedrals.forEach(d => {
        accordion += `<tr><td>${d.residues}</td><td>${d.dihedral_degrees}</td></tr>`;
    });
    accordion += '</tbody></table></div></div></div></div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('geometryResults', {
            title: 'Structure Geometry',
            meta: `${geometry.distances.length} distances · ${geometry.angles.length} angles · ${geometry.dihedrals.length} dihedrals`,
            summary: tilesHtml + accordion,
            raw: JSON.stringify(geometry, null, 2),
            workspaceItem: { type: 'structure-geometry', name: `Geometry (${geometry.distances.length} distances)`, data: geometry },
            downloads: [
                { label: 'JSON', filename: 'structure-geometry.json', text: JSON.stringify(geometry, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('geometryResults').innerHTML = tilesHtml + accordion;
    }
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('qualityBtn', '<i class="fas fa-check-circle me-2"></i>Analyze Quality');
    });
});

function displayQuality(quality) {
    const tiles = [
        { label: 'Mean B-factor',  value: quality.bfactor_stats.mean ?? '—' },
        { label: 'Median B-factor', value: quality.bfactor_stats.median ?? '—' },
        { label: 'Mean Occupancy', value: quality.occupancy_stats.mean ?? '—' },
        { label: 'Incomplete',     value: quality.residue_completeness.length },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let details = '<div class="row g-3">';
    if (quality.bfactor_stats.mean != null) {
        details += '<div class="col-md-6"><h6 class="small">B-factor Statistics</h6><table class="table table-sm">';
        details += `<tr><td>Mean</td><td>${quality.bfactor_stats.mean}</td></tr>`;
        details += `<tr><td>Median</td><td>${quality.bfactor_stats.median}</td></tr>`;
        details += `<tr><td>Range</td><td>${quality.bfactor_stats.min} - ${quality.bfactor_stats.max}</td></tr>`;
        details += '</table></div>';
    }
    if (quality.occupancy_stats.mean != null) {
        details += '<div class="col-md-6"><h6 class="small">Occupancy Statistics</h6><table class="table table-sm">';
        details += `<tr><td>Mean</td><td>${quality.occupancy_stats.mean}</td></tr>`;
        details += `<tr><td>Range</td><td>${quality.occupancy_stats.min} - ${quality.occupancy_stats.max}</td></tr>`;
        details += '</table></div>';
    }
    details += '</div>';
    if (quality.residue_completeness.length > 0) {
        details += '<div class="alert alert-warning small mt-3 mb-0"><strong>Warning:</strong> ' +
            quality.residue_completeness.length + ' residues have missing backbone atoms</div>';
    }

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('qualityResults', {
            title: 'Structure Quality',
            meta: `B-factor mean ${quality.bfactor_stats.mean ?? '—'} · ${quality.residue_completeness.length} incomplete`,
            summary: tilesHtml,
            details: details,
            raw: JSON.stringify(quality, null, 2),
            workspaceItem: { type: 'structure-quality', name: `Quality (${quality.residue_completeness.length} incomplete)`, data: quality },
            downloads: [
                { label: 'JSON', filename: 'structure-quality.json', text: JSON.stringify(quality, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('qualityResults').innerHTML = tilesHtml + details;
    }
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('contactsBtn', '<i class="fas fa-network-wired me-2"></i>Calculate Contacts');
    });
});

function displayContacts(data) {
    const tiles = [
        { label: 'Contacts', value: data.contact_count },
        { label: 'Cutoff',   value: data.cutoff + ' Å' },
        { label: 'Shown',    value: data.contacts.length },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    data.contacts.forEach(c => { rows += `<tr><td>${c.residue1}</td><td>${c.residue2}</td><td>${c.distance}</td></tr>`; });
    const table = `<div class="table-responsive"><table class="table table-sm table-hover">
        <thead class="table-light"><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead>
        <tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('contactsResults', {
            title: 'Structure Contacts',
            meta: `${data.contact_count} contacts within ${data.cutoff} Å`,
            summary: tilesHtml + table,
            raw: JSON.stringify(data.contacts, null, 2),
            workspaceItem: { type: 'structure-contacts', name: `Contacts (${data.contact_count})`, data: data.contacts },
            downloads: [
                { label: 'JSON', filename: 'structure-contacts.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('contactsResults').innerHTML = tilesHtml + table;
    }

    const chain = (document.getElementById('contactChainId').value || 'A').toUpperCase();
    const resiMap = {};
    data.contacts.forEach(c => {
        [c.residue1, c.residue2].forEach(r => {
            const num = _structureResidueNum(r);
            if (num !== null) resiMap[`${chain}:${num}`] = '#f59e0b';
        });
    });
    applyResidueColors(resiMap);
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('interactionsBtn', '<i class="fas fa-link me-2"></i>Detect Interactions');
    });
});

function displayInteractions(data) {
    const tiles = [
        { label: 'H-Bonds',      value: data.hbond_count },
        { label: 'Salt Bridges', value: data.salt_bridge_count },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let accordion = '<div class="accordion" id="interactionsAccordion">';
    accordion += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button" data-bs-toggle="collapse" data-bs-target="#hbonds">';
    accordion += `Hydrogen Bonds (${data.hbond_count} found, showing first 50)</button></h2>`;
    accordion += '<div id="hbonds" class="accordion-collapse collapse show" data-bs-parent="#interactionsAccordion"><div class="accordion-body">';
    accordion += '<table class="table table-sm"><thead><tr><th>Donor</th><th>Acceptor</th><th>Distance (Å)</th></tr></thead><tbody>';
    data.hbonds.forEach(hb => { accordion += `<tr><td>${hb.donor}</td><td>${hb.acceptor}</td><td>${hb.distance}</td></tr>`; });
    accordion += '</tbody></table></div></div></div>';

    accordion += '<div class="accordion-item"><h2 class="accordion-header"><button class="accordion-button collapsed" data-bs-toggle="collapse" data-bs-target="#saltbridges">';
    accordion += `Salt Bridges (${data.salt_bridge_count} found)</button></h2>`;
    accordion += '<div id="saltbridges" class="accordion-collapse collapse" data-bs-parent="#interactionsAccordion"><div class="accordion-body">';
    accordion += '<table class="table table-sm"><thead><tr><th>Residue 1</th><th>Residue 2</th><th>Distance (Å)</th></tr></thead><tbody>';
    data.salt_bridges.forEach(sb => { accordion += `<tr><td>${sb.residue1}</td><td>${sb.residue2}</td><td>${sb.distance}</td></tr>`; });
    accordion += '</tbody></table></div></div></div></div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('interactionsResults', {
            title: 'Structure Interactions',
            meta: `${data.hbond_count} H-bonds · ${data.salt_bridge_count} salt bridges`,
            summary: tilesHtml + accordion,
            raw: JSON.stringify(data, null, 2),
            workspaceItem: { type: 'structure-interactions', name: `Interactions (${data.hbond_count + data.salt_bridge_count})`, data: data },
            downloads: [
                { label: 'JSON', filename: 'structure-interactions.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('interactionsResults').innerHTML = tilesHtml + accordion;
    }
}

// Helper Functions: showLoading, hideLoading, showAlert are provided globally by static/js/utils.js

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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('dsspBtn', '<i class="fas fa-dna me-2"></i>Analyze Secondary Structure');
        showAlert(friendlyError(error, 'pdb'), 'danger');
    });
});

function displayDsspResults(result) {
    if (!result.dssp_available) {
        document.getElementById('dsspResults').innerHTML =
            '<div class="alert alert-warning mb-0"><strong>DSSP not available or failed.</strong> ' + (result.dssp_error || '') + '</div>';
        return;
    }
    const ss = result.secondary_structure;
    const total = Object.values(ss.counts).reduce((a, b) => a + b, 0);
    const topKeys = Object.entries(ss.counts).sort((a, b) => b[1] - a[1]).slice(0, 4);
    const tiles = [
        { label: 'Residues', value: total },
    ];
    topKeys.forEach(([key, count]) => {
        tiles.push({ label: (ss.mapping[key] || key), value: count });
    });
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let distribution = '<h6 class="mb-2"><i class="fas fa-chart-pie"></i> Secondary Structure Distribution</h6><div class="row g-2 mb-3">';
    for (const [key, count] of Object.entries(ss.counts)) {
        const name = ss.mapping[key] || key;
        distribution += `<div class="col-md-4"><div class="border rounded p-2"><strong>${name} (${key}):</strong> ${count}</div></div>`;
    }
    distribution += '</div>';

    let table = '';
    if (result.residue_details.length > 0) {
        table += '<h6 class="mb-2"><i class="fas fa-list"></i> Residue Details (first 50)</h6>';
        table += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>SS</th><th>Type</th><th>Accessibility</th><th>Phi</th><th>Psi</th></tr></thead><tbody>';
        result.residue_details.forEach(r => {
            table += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.ss}</td><td class="small">${r.ss_name}</td><td>${r.accessibility}</td><td>${r.phi !== null ? r.phi + '°' : '-'}</td><td>${r.psi !== null ? r.psi + '°' : '-'}</td></tr>`;
        });
        table += '</tbody></table></div>';
    }

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('dsspResults', {
            title: 'DSSP Secondary Structure',
            meta: `${total} residues classified`,
            summary: tilesHtml + distribution,
            details: table,
            raw: JSON.stringify(result, null, 2),
            workspaceItem: { type: 'dssp', name: `DSSP (${total} residues)`, data: result },
            downloads: [
                { label: 'JSON', filename: 'dssp.json', text: JSON.stringify(result, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('dsspResults').innerHTML = tilesHtml + distribution + table;
    }

    if (result.residue_details) {
        const resiMap = {};
        result.residue_details.forEach(r => {
            const num = _structureResidueNum(r.residue);
            if (num !== null) resiMap[`${r.chain}:${num}`] = _structureSsColor(r.ss);
        });
        applyResidueColors(resiMap);
    }
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('ramaBtn', '<i class="fas fa-chart-scatter me-2"></i>Calculate Phi/Psi Angles');
        showAlert(friendlyError(error, 'pdb'), 'danger');
    });
});

function displayRamaResults(data) {
    const tiles = [
        { label: 'Residues', value: data.total_residues },
        { label: 'Favored',  value: data.classification.favored },
        { label: 'Allowed',  value: data.classification.allowed },
        { label: 'Outliers', value: data.classification.outliers },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let summary = tilesHtml +
        '<div class="row g-2 mb-3">' +
        `<div class="col-md-4"><div class="border rounded p-2 bg-success text-white text-center"><strong>Favored:</strong> ${data.classification.favored}</div></div>` +
        `<div class="col-md-4"><div class="border rounded p-2 bg-warning text-center"><strong>Allowed:</strong> ${data.classification.allowed}</div></div>` +
        `<div class="col-md-4"><div class="border rounded p-2 bg-danger text-white text-center"><strong>Outliers:</strong> ${data.classification.outliers}</div></div>` +
        '</div>';

    let table = '<h6 class="mb-2"><i class="fas fa-table"></i> Phi/Psi Angles (first 100)</h6>';
    table += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Residue</th><th>Phi (°)</th><th>Psi (°)</th></tr></thead><tbody>';
    data.phi_psi_data.forEach(r => { table += `<tr><td>${r.residue}</td><td>${r.phi}</td><td>${r.psi}</td></tr>`; });
    table += '</tbody></table></div>';

    const tsv = ['residue\tphi\tpsi'];
    data.phi_psi_data.forEach(r => tsv.push([r.residue, r.phi, r.psi].join('\t')));

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('ramaResults', {
            title: 'Ramachandran',
            meta: `${data.total_residues} residues · ${data.classification.outliers} outliers`,
            summary: summary,
            details: table,
            raw: JSON.stringify(data, null, 2),
            workspaceItem: { type: 'ramachandran', name: `Ramachandran (${data.total_residues} residues)`, data: data },
            downloads: [
                { label: 'TSV',  filename: 'ramachandran.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'ramachandran.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('ramaResults').innerHTML = summary + table;
    }

    // Highlight Ramachandran outliers (phi/psi outside favored/allowed) on the viewer.
    const chain = (document.getElementById('ramaChainId').value || 'A').toUpperCase();
    const outliers = {};
    data.phi_psi_data.forEach(r => {
        const phi = r.phi, psi = r.psi;
        const favored =
            (phi >= -180 && phi <= -30 && psi >= -90  && psi <= 30) ||
            (phi >= -180 && phi <= -30 && psi >= 60   && psi <= 180);
        const allowed =
            (phi >= -180 && phi <= 0   && psi >= -180 && psi <= 180);
        if (!favored && !allowed) {
            const num = _structureResidueNum(r.residue);
            if (num !== null) outliers[`${chain}:${num}`] = '#ef4444';
        }
    });
    applyResidueColors(outliers);
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('sasaBtn', '<i class="fas fa-water me-2"></i>Calculate SASA');
        showAlert(friendlyError(error, 'pdb'), 'danger');
    });
});

function displaySasaResults(data) {
    let tiles = [];
    let summary = '';
    if (data.total_sasa !== undefined) {
        tiles.push({ label: 'Total SASA', value: data.total_sasa + ' Ų' });
        tiles.push({ label: 'Chains', value: Object.keys(data.chain_sasa || {}).length });
        tiles.push({ label: 'Residues', value: (data.residue_sasa || []).length });
        summary += '<h6 class="mb-2">Chain SASA:</h6><div class="row g-2 mb-3">';
        for (const [chain, sasa] of Object.entries(data.chain_sasa || {})) {
            summary += `<div class="col-md-3"><div class="border rounded p-2"><strong>Chain ${chain}:</strong> ${sasa} Ų</div></div>`;
        }
        summary += '</div>';
        if (data.residue_sasa && data.residue_sasa.length > 0) {
            summary += '<h6 class="mb-2">Residue SASA (first 50):</h6>';
            summary += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>SASA (Ų)</th></tr></thead><tbody>';
            data.residue_sasa.forEach(r => { summary += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.sasa}</td></tr>`; });
            summary += '</tbody></table></div>';
        }
    } else if (data.method) {
        tiles.push({ label: 'Method', value: data.method });
        tiles.push({ label: 'Chains', value: Object.keys(data.chain_accessibility || {}).length });
        tiles.push({ label: 'Residues', value: (data.residue_accessibility || []).length });
        summary += '<h6 class="mb-2">Chain Accessibility:</h6><div class="row g-2 mb-3">';
        for (const [chain, acc] of Object.entries(data.chain_accessibility || {})) {
            summary += `<div class="col-md-3"><div class="border rounded p-2"><strong>Chain ${chain}:</strong> ${acc}</div></div>`;
        }
        summary += '</div>';
        if (data.residue_accessibility) {
            summary += '<h6 class="mb-2">Relative Accessibility (first 50):</h6>';
            summary += '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr><th>Chain</th><th>Residue</th><th>Accessibility</th></tr></thead><tbody>';
            data.residue_accessibility.forEach(r => { summary += `<tr><td>${r.chain}</td><td>${r.residue}</td><td>${r.relative_accessibility}</td></tr>`; });
            summary += '</tbody></table></div>';
        }
    }
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('sasaResults', {
            title: 'Structure SASA',
            meta: data.total_sasa !== undefined ? `${data.total_sasa} Ų total` : (data.method || 'SASA'),
            summary: tilesHtml + summary,
            raw: JSON.stringify(data, null, 2),
            workspaceItem: { type: 'sasa', name: `SASA analysis`, data: data },
            downloads: [
                { label: 'JSON', filename: 'sasa.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('sasaResults').innerHTML = tilesHtml + summary;
    }

    // Color residues on the 3D viewer by relative solvent accessibility.
    const resiMap = {};
    const rows = (data.residue_sasa && data.residue_sasa.length)
        ? data.residue_sasa.map(r => ({ chain: r.chain, residue: r.residue, v: r.sasa }))
        : (data.residue_accessibility || []).map(r => ({ chain: r.chain, residue: r.residue, v: r.relative_accessibility }));
    if (rows.length) {
        const values = rows.map(r => r.v).filter(v => typeof v === 'number');
        const min = Math.min.apply(null, values);
        const max = Math.max.apply(null, values);
        rows.forEach(r => {
            const num = _structureResidueNum(r.residue);
            if (num === null || typeof r.v !== 'number') return;
            const t = max > min ? (r.v - min) / (max - min) : 0.5;
            resiMap[`${r.chain}:${num}`] = _structureGradient(t);
        });
        applyResidueColors(resiMap);
    }
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
            showAlert(friendlyError(data.error, 'pdb'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('selectionBtn', '<i class="fas fa-crop me-2"></i>Extract Selection');
        showAlert(friendlyError(error, 'pdb'), 'danger');
    });
});

function displaySelectionResults(extraction) {
    const tiles = [
        { label: 'Type',    value: extraction.selection_type },
        { label: 'Value',   value: extraction.selection_value },
        { label: 'Matched', value: extraction.extracted_count },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let table = '';
    if (extraction.details.length > 0) {
        const keys = Object.keys(extraction.details[0]);
        table = '<div class="table-responsive"><table class="table table-sm table-bordered"><thead><tr>';
        keys.forEach(k => { table += `<th>${k}</th>`; });
        table += '</tr></thead><tbody>';
        extraction.details.forEach(item => {
            table += '<tr>';
            keys.forEach(k => { table += `<td>${item[k]}</td>`; });
            table += '</tr>';
        });
        table += '</tbody></table></div>';
    }

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('selectionResults', {
            title: 'Structure Selection',
            meta: `${extraction.extracted_count} ${extraction.selection_type}`,
            summary: tilesHtml + table,
            raw: JSON.stringify(extraction, null, 2),
            workspaceItem: { type: 'structure-selection', name: `Selection ${extraction.selection_type}=${extraction.selection_value}`, data: extraction },
            downloads: [
                { label: 'JSON', filename: 'selection.json', text: JSON.stringify(extraction, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('selectionResults').innerHTML = tilesHtml + table;
    }
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
