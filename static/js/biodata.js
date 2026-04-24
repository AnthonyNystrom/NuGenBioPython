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
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('codonBtn', '<i class="fas fa-exchange-alt me-2"></i>Translate');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayCodonResults(data) {
    const tiles = [
        { label: 'Table',      value: `#${data.table_id}` },
        { label: 'DNA Length', value: data.length_original + ' bp' },
        { label: 'Protein',    value: data.length_translated + ' AA' },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const summary = tilesHtml + `
        <p class="small mb-2"><strong>Table:</strong> ${data.table_name}</p>
        <div class="mb-2"><label class="small"><strong>Translated:</strong></label>
        <textarea class="form-control form-control-sm sequence-display" rows="3" readonly>${data.translated_sequence}</textarea></div>
        <p class="small mb-1"><strong>Start codons:</strong> ${data.start_codons.map(c => '<code class="badge bg-success">'+c+'</code>').join(' ')}</p>
        <p class="small mb-0"><strong>Stop codons:</strong> ${data.stop_codons.map(c => '<code class="badge bg-danger">'+c+'</code>').join(' ')}</p>
    `;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('codonResults', {
            title: 'Codon Translation',
            meta: `${data.table_name} · ${data.length_original} bp → ${data.length_translated} AA`,
            summary: summary,
            raw: data.translated_sequence,
            copyText: data.translated_sequence,
            workspaceItem: { type: 'protein', name: `Translation (${data.length_translated} AA)`, data: data.translated_sequence, meta: { kind: 'protein' } },
            downloads: [
                { label: 'FASTA', filename: 'translation.fasta', text: `>translation_table_${data.table_id}\n${data.translated_sequence}`, mime: 'text/plain' },
                { label: 'JSON',  filename: 'translation.json',  text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('codonResults').innerHTML = summary;
    }
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
            displayIUPACLookup(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('iupacBtn', '<i class="fas fa-search me-2"></i>Lookup');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayIUPACLookup(data) {
    const tiles = [
        { label: 'Code',       value: data.code },
        { label: 'Bases',      value: data.bases },
        { label: 'Complement', value: data.complement },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('iupacResults', {
            title: 'IUPAC Lookup',
            meta: `Code ${data.code} → ${data.bases}`,
            summary: tilesHtml,
            raw: JSON.stringify(data, null, 2),
            copyText: data.bases,
            workspaceItem: { type: 'iupac', name: `IUPAC ${data.code}`, data: data },
            downloads: [
                { label: 'JSON', filename: 'iupac.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('iupacResults').innerHTML = tilesHtml;
    }
}

function showAllIUPACCodes() {
    const codeType = document.getElementById('iupacCodeType').value;
    fetch('/api/biodata/iupac_codes?type=' + codeType)
    .then(response => response.json())
    .then(data => {
        if (data.success) displayIUPACTable(data.codes, codeType);
    });
}

function displayIUPACTable(codes, codeType) {
    const entries = Object.entries(codes.values);
    const tiles = [{ label: 'Codes', value: entries.length }, { label: 'Type', value: codeType }];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    const tsv = ['code\tbases\tcomplement'];
    entries.forEach(([code, bases]) => {
        const complement = codes.complement[code] || 'N';
        rows += `<tr><td><code class="badge bg-primary">${code}</code></td><td class="small">${bases}</td><td><code class="badge bg-info">${complement}</code></td></tr>`;
        tsv.push([code, bases, complement].join('\t'));
    });
    const table = `<div class="table-responsive"><table class="table table-sm mb-0">
        <thead><tr><th>Code</th><th>Bases</th><th>Complement</th></tr></thead><tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('iupacResults', {
            title: 'IUPAC Codes',
            meta: `${codeType} · ${entries.length} codes`,
            summary: tilesHtml + table,
            raw: tsv.join('\n'),
            downloads: [
                { label: 'TSV',  filename: `iupac-${codeType}.tsv`,  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: `iupac-${codeType}.json`, text: JSON.stringify(codes, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('iupacResults').innerHTML = tilesHtml + table;
    }
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
            displayProteinConvert(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('proteinBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayProteinConvert(data) {
    const tiles = [
        { label: 'Residues', value: data.count },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    const summary = tilesHtml +
        `<p class="small mb-1"><strong>Converted:</strong></p>` +
        `<textarea class="form-control form-control-sm sequence-display" rows="3" readonly>${data.result}</textarea>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('iupacResults', {
            title: 'Protein Conversion',
            meta: `${data.count} residues`,
            summary: summary,
            raw: data.result,
            copyText: data.result,
            workspaceItem: { type: 'protein', name: `Converted protein (${data.count} residues)`, data: data.result, meta: { kind: 'protein' } },
            downloads: [
                { label: 'Text', filename: 'protein-convert.txt', text: data.result, mime: 'text/plain' },
            ],
        });
    } else {
        document.getElementById('iupacResults').innerHTML = summary;
    }
}

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
            displayMolecularWeight(data);
        } else {
            showAlert(friendlyError(data.error, 'server'), 'danger');
        }
    })
    .catch(error => {
        hideLoading('weightBtn', '<i class="fas fa-calculator me-2"></i>Calculate');
        showAlert(friendlyError(error, 'server'), 'danger');
    });
});

function displayMolecularWeight(data) {
    const tiles = [
        { label: 'Weight', value: data.weight + ' Da' },
        { label: 'Length', value: data.length },
        { label: 'Type',   value: data.weight_type },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let composition = '<p class="small mb-1"><strong>Composition:</strong></p><div class="d-flex flex-wrap gap-1">';
    for (let [unit, count] of Object.entries(data.composition)) {
        composition += `<span class="badge bg-secondary">${unit}: ${count}</span>`;
    }
    composition += '</div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('iupacResults', {
            title: 'Molecular Weight',
            meta: `${data.weight} Da · ${data.length} ${data.weight_type}`,
            summary: tilesHtml + composition,
            raw: JSON.stringify(data, null, 2),
            workspaceItem: { type: 'molecular-weight', name: `${data.weight} Da (${data.length} ${data.weight_type})`, data: data },
            downloads: [
                { label: 'JSON', filename: 'molecular-weight.json', text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('iupacResults').innerHTML = tilesHtml + composition;
    }
}

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
        if (data.success) displayPDBConversionTable(data.conversions, convType);
    });
}

function displayPDBConversionTable(conversions, convType) {
    const entries = Object.entries(conversions);
    const tiles = [{ label: 'Mappings', value: entries.length }, { label: 'Type', value: convType }];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    const tsv = ['from\tto'];
    entries.forEach(([from, to]) => {
        rows += `<tr><td><code class="badge bg-primary">${from}</code></td><td><code class="badge bg-success">${to}</code></td></tr>`;
        tsv.push([from, to].join('\t'));
    });
    const table = `<div class="table-responsive"><table class="table table-sm table-hover mb-0">
        <thead><tr><th>From</th><th>To</th></tr></thead><tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('pdbConversionTable', {
            title: 'PDB Conversion',
            meta: `${convType} · ${entries.length} mappings`,
            summary: tilesHtml + table,
            raw: tsv.join('\n'),
            downloads: [
                { label: 'TSV',  filename: `pdb-${convType}.tsv`,  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: `pdb-${convType}.json`, text: JSON.stringify(conversions, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('pdbConversionTable').innerHTML = tilesHtml + table;
    }
}

function loadAtomWeights() {
    fetch('/api/biodata/atom_weights')
    .then(response => response.json())
    .then(data => {
        if (data.success) displayAtomWeightsTable(data.weights);
    });
}

function displayAtomWeightsTable(weights) {
    const elements = {H:'Hydrogen',C:'Carbon',N:'Nitrogen',O:'Oxygen',P:'Phosphorus',S:'Sulfur',Na:'Sodium',Mg:'Magnesium',Cl:'Chlorine',K:'Potassium',Ca:'Calcium',Fe:'Iron',Zn:'Zinc',Se:'Selenium'};
    const present = Object.entries(elements).filter(([sym]) => weights[sym]);
    const tiles = [{ label: 'Elements', value: present.length }];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    const tsv = ['element\tsymbol\tweight_da'];
    present.forEach(([symbol, name]) => {
        rows += `<tr><td class="small">${name}</td><td><code class="badge bg-primary">${symbol}</code></td><td>${weights[symbol]}</td></tr>`;
        tsv.push([name, symbol, weights[symbol]].join('\t'));
    });
    const table = `<div class="table-responsive" style="max-height:400px;overflow-y:auto;"><table class="table table-sm table-hover mb-0">
        <thead class="sticky-top bg-white"><tr><th>Element</th><th>Symbol</th><th>Weight (Da)</th></tr></thead><tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('pdbConversionTable', {
            title: 'Atomic Weights',
            meta: `${present.length} elements`,
            summary: tilesHtml + table,
            raw: tsv.join('\n'),
            downloads: [
                { label: 'TSV',  filename: 'atom-weights.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'atom-weights.json', text: JSON.stringify(weights, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('pdbConversionTable').innerHTML = tilesHtml + table;
    }
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
