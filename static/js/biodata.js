/**
 * BioData JavaScript - Handles CodonTable, IUPACData, and PDBData operations
 */

// Codon Table Form Handler
document.getElementById('codonForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const sequence = document.getElementById('codonSeq').value.trim();
    const tableId = document.getElementById('codonTableId').value;

    if (!sequence) {
        showAlert('Please enter a DNA sequence', 'warning');
        return;
    }

    showLoading('codonBtn');

    fetch('/api/biodata/translate', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({sequence: sequence, table_id: tableId})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('codonBtn', '<i class="fas fa-exchange-alt me-2"></i>Translate');
        if (data.success) {
            displayCodonResults(data);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('codonBtn', '<i class="fas fa-exchange-alt me-2"></i>Translate');
        showAlert('Network error: ' + error.message, 'danger');
    });
});

function displayCodonResults(data) {
    const resultsDiv = document.getElementById('codonResults');
    let html = `
        <p class="small mb-1"><strong>Table:</strong> ${data.table_name} (ID: ${data.table_id})</p>
        <p class="small mb-1"><strong>Original:</strong> ${data.length_original} bp</p>
        <p class="small mb-1"><strong>Translated:</strong> ${data.length_translated} AA</p>
        <div class="mb-2"><label class="small"><strong>Result:</strong></label>
        <textarea class="form-control form-control-sm sequence-display" rows="3" readonly>${data.translated_sequence}</textarea></div>
        <p class="small mb-1"><strong>Start:</strong> ${data.start_codons.map(c => '<code class="badge bg-success">'+c+'</code>').join(' ')}</p>
        <p class="small mb-0"><strong>Stop:</strong> ${data.stop_codons.map(c => '<code class="badge bg-danger">'+c+'</code>').join(' ')}</p>
    `;
    resultsDiv.innerHTML = html;
}

function loadCodonTables() {
    fetch('/api/biodata/codon_tables')
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            displayCodonTables(data.tables);
        }
    });
}

function displayCodonTables(tables) {
    const div = document.getElementById('codonTablesDiv');
    let html = '<div class="table-responsive"><table class="table table-sm table-hover mb-0"><thead><tr><th>ID</th><th>Name</th><th>Start</th><th>Stop</th></tr></thead><tbody>';
    Object.values(tables).forEach(table => {
        html += `<tr><td><span class="badge bg-primary">${table.id}</span></td><td class="small">${table.name}</td>
                 <td>${table.start_codons.map(c => '<code class="badge bg-success">'+c+'</code>').join(' ')}</td>
                 <td>${table.stop_codons.map(c => '<code class="badge bg-danger">'+c+'</code>').join(' ')}</td></tr>`;
    });
    html += '</tbody></table></div>';
    div.innerHTML = html;
}

function loadCodonExample() {
    document.getElementById('codonSeq').value = 'ATGAAATTTAAAGGTCTCGACACCCTGAAGAAAGTTTATGGTGCTATTGGTGGCGGTATTGGTGCTATGGGTATGATTCTGAAAAAACTTGGT';
    document.getElementById('codonTableId').value = '1';
}

// IUPAC Lookup Form Handler
document.getElementById('iupacLookupForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const code = document.getElementById('iupacCode').value.trim();
    const codeType = document.getElementById('iupacCodeType').value;

    if (!code) {
        showAlert('Please enter an IUPAC code', 'warning');
        return;
    }

    showLoading('iupacBtn');

    fetch('/api/biodata/iupac_lookup', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({code: code, type: codeType})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('iupacBtn', '<i class="fas fa-search me-2"></i>Lookup');
        if (data.success) {
            const resultsDiv = document.getElementById('iupacResults');
            resultsDiv.innerHTML = `
                <p class="small mb-1"><strong>Code:</strong> <code class="badge bg-primary">${data.code}</code></p>
                <p class="small mb-1"><strong>Bases:</strong> ${data.bases}</p>
                <p class="small mb-0"><strong>Complement:</strong> <code class="badge bg-info">${data.complement}</code></p>
            `;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('iupacBtn', '<i class="fas fa-search me-2"></i>Lookup');
        showAlert('Network error: ' + error.message, 'danger');
    });
});

function showAllIUPACCodes() {
    const codeType = document.getElementById('iupacCodeType').value;
    fetch('/api/biodata/iupac_codes?type=' + codeType)
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            const resultsDiv = document.getElementById('iupacResults');
            let html = '<div class="table-responsive"><table class="table table-sm mb-0"><thead><tr><th>Code</th><th>Bases</th><th>Complement</th></tr></thead><tbody>';
            for (let [code, bases] of Object.entries(data.codes.values)) {
                const complement = data.codes.complement[code] || 'N';
                html += `<tr><td><code class="badge bg-primary">${code}</code></td><td class="small">${bases}</td><td><code class="badge bg-info">${complement}</code></td></tr>`;
            }
            html += '</tbody></table></div>';
            resultsDiv.innerHTML = html;
        }
    });
}

// Protein Converter Form Handler
document.getElementById('proteinConvertForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const input = document.getElementById('proteinInput').value.trim();
    const convType = document.getElementById('proteinConvType').value;

    if (!input) {
        showAlert('Please enter protein sequence', 'warning');
        return;
    }

    showLoading('proteinBtn');

    fetch('/api/biodata/convert_protein', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({input: input, conversion_type: convType})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('proteinBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert');
        if (data.success) {
            const resultsDiv = document.getElementById('iupacResults');
            resultsDiv.innerHTML = `
                <p class="small mb-1"><strong>Result:</strong></p>
                <textarea class="form-control form-control-sm sequence-display" rows="3" readonly>${data.result}</textarea>
                <p class="small mt-2 mb-0"><strong>Count:</strong> ${data.count} residues</p>
            `;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('proteinBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert');
        showAlert('Network error: ' + error.message, 'danger');
    });
});

function loadProteinExample() {
    document.getElementById('proteinInput').value = 'ACDEFGHIKLMNPQRSTVWY';
    document.getElementById('proteinConvType').value = '1to3';
}

// Molecular Weight Form Handler
document.getElementById('weightForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const sequence = document.getElementById('weightSeq').value.trim();
    const seqType = document.getElementById('weightType').value;
    const weightType = document.getElementById('weightCalcType').value;

    if (!sequence) {
        showAlert('Please enter a sequence', 'warning');
        return;
    }

    showLoading('weightBtn');

    fetch('/api/biodata/molecular_weight', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({sequence: sequence, seq_type: seqType, weight_type: weightType})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('weightBtn', '<i class="fas fa-calculator me-2"></i>Calculate');
        if (data.success) {
            const resultsDiv = document.getElementById('iupacResults');
            let html = `
                <p class="small mb-1"><strong>Molecular Weight:</strong> ${data.weight} Da</p>
                <p class="small mb-1"><strong>Length:</strong> ${data.length}</p>
                <p class="small mb-1"><strong>Type:</strong> ${data.weight_type}</p>
                <p class="small mb-1"><strong>Composition:</strong></p>
                <div class="d-flex flex-wrap gap-1">
            `;
            for (let [unit, count] of Object.entries(data.composition)) {
                html += `<span class="badge bg-secondary">${unit}: ${count}</span>`;
            }
            html += '</div>';
            resultsDiv.innerHTML = html;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('weightBtn', '<i class="fas fa-calculator me-2"></i>Calculate');
        showAlert('Network error: ' + error.message, 'danger');
    });
});

function loadWeightExample() {
    document.getElementById('weightSeq').value = 'ACDEFGHIKLMNPQRSTVWY';
    document.getElementById('weightType').value = 'protein';
    document.getElementById('weightCalcType').value = 'average';
}

// PDB Conversion Functions
function loadPDBConversion(convType) {
    fetch('/api/biodata/pdb_conversions?type=' + convType)
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            const div = document.getElementById('pdbConversionTable');
            let html = '<div class="table-responsive"><table class="table table-sm table-hover mb-0"><thead><tr><th>From</th><th>To</th></tr></thead><tbody>';
            for (let [from, to] of Object.entries(data.conversions)) {
                html += `<tr><td><code class="badge bg-primary">${from}</code></td><td><code class="badge bg-success">${to}</code></td></tr>`;
            }
            html += '</tbody></table></div>';
            div.innerHTML = html;
        }
    });
}

function loadAtomWeights() {
    fetch('/api/biodata/atom_weights')
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            const div = document.getElementById('pdbConversionTable');
            let html = '<div class="table-responsive" style="max-height:400px;overflow-y:auto;"><table class="table table-sm table-hover mb-0"><thead class="sticky-top bg-white"><tr><th>Element</th><th>Symbol</th><th>Weight (Da)</th></tr></thead><tbody>';
            const elements = {H:'Hydrogen',C:'Carbon',N:'Nitrogen',O:'Oxygen',P:'Phosphorus',S:'Sulfur',Na:'Sodium',Mg:'Magnesium',Cl:'Chlorine',K:'Potassium',Ca:'Calcium',Fe:'Iron',Zn:'Zinc',Se:'Selenium'};
            for (let [symbol, name] of Object.entries(elements)) {
                if (data.weights[symbol]) {
                    html += `<tr><td class="small">${name}</td><td><code class="badge bg-primary">${symbol}</code></td><td>${data.weights[symbol]}</td></tr>`;
                }
            }
            html += '</tbody></table></div>';
            div.innerHTML = html;
        }
    });
}

function showPDBInfo() {
    const div = document.getElementById('pdbConversionTable');
    div.innerHTML = `
        <div class="alert alert-info mb-0">
            <h6 class="small"><i class="fas fa-info-circle"></i> PDB Residue Naming</h6>
            <p class="small mb-2">The Protein Data Bank (PDB) uses 3-letter codes for amino acid and nucleic acid residues in structure files.</p>
            <p class="small mb-1"><strong>Common Conversions:</strong></p>
            <ul class="small mb-0">
                <li>Protein: ALA (A), CYS (C), ASP (D), GLU (E), etc.</li>
                <li>DNA: DA (A), DC (C), DG (G), DT (T)</li>
                <li>RNA: A (A), C (C), G (G), U (U)</li>
            </ul>
        </div>
    `;
}
