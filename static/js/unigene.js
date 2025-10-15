/**
 * UniGene Analysis JavaScript
 * Handles Parse Multiple and Read Single tabs
 */

// Parse Multiple Tab
document.getElementById('parseForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const fileInput = document.getElementById('parseFile');
    const maxRecords = document.getElementById('parseMaxRecords').value;

    if (!fileInput.files[0]) {
        showAlert('Please select a UniGene file', 'warning');
        return;
    }

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('max_records', maxRecords);

    showLoading('parseBtn');

    fetch('/api/unigene/parse', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');

        if (data.success) {
            displayParseResults(data.records, data.count);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');
        showAlert('Error parsing file: ' + error.message, 'danger');
    });
});

// Read Single Tab
document.getElementById('readForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const fileInput = document.getElementById('readFile');

    if (!fileInput.files[0]) {
        showAlert('Please select a UniGene file', 'warning');
        return;
    }

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);

    showLoading('readBtn');

    fetch('/api/unigene/read', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Record');

        if (data.success) {
            displayReadResults(data.record);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Record');
        showAlert('Error reading file: ' + error.message, 'danger');
    });
});

function displayParseResults(records, count) {
    const resultsDiv = document.getElementById('parseResults');
    const countBadge = document.getElementById('parseCount');

    countBadge.textContent = `${count} clusters`;
    countBadge.style.display = 'inline-block';

    if (records.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted small text-center mb-0">No records found</p>';
        return;
    }

    let html = '<div class="accordion" id="parseAccordion">';

    records.forEach((record, index) => {
        const accordionId = `parseRecord${index}`;
        const isFirst = index === 0;

        html += `
            <div class="accordion-item border mb-1">
                <h2 class="accordion-header" id="heading${index}">
                    <button class="accordion-button ${isFirst ? '' : 'collapsed'} p-2 small" type="button" data-bs-toggle="collapse" data-bs-target="#${accordionId}" aria-expanded="${isFirst}" aria-controls="${accordionId}">
                        <strong>${record.cluster_id}</strong>
                        <span class="badge bg-info ms-2">${record.gene_symbol}</span>
                        <span class="ms-2 text-muted">${record.title}</span>
                    </button>
                </h2>
                <div id="${accordionId}" class="accordion-collapse collapse ${isFirst ? 'show' : ''}" aria-labelledby="heading${index}" data-bs-parent="#parseAccordion">
                    <div class="accordion-body p-2">
                        <div class="row g-2">
                            <div class="col-md-6">
                                <p class="small mb-1"><strong>Cluster ID:</strong> <code>${record.cluster_id}</code></p>
                                <p class="small mb-1"><strong>Gene:</strong> ${record.gene_symbol}</p>
                                <p class="small mb-1"><strong>Species:</strong> ${record.species}</p>
                                <p class="small mb-1"><strong>Chromosome:</strong> ${record.chromosome}</p>
                            </div>
                            <div class="col-md-6">
                                <p class="small mb-1"><strong>Title:</strong> ${record.title}</p>
                                <p class="small mb-1"><strong>Sequences:</strong> <span class="badge bg-secondary">${record.sequence_count}</span></p>
                                <p class="small mb-1"><strong>Tissues:</strong> <span class="badge bg-info">${record.tissue_count}</span></p>
                            </div>
                        </div>

                        ${record.tissues && record.tissues.length > 0 ? `
                            <hr class="my-2">
                            <p class="small mb-1"><strong>Expression Tissues (${record.tissues.length}):</strong></p>
                            <div class="d-flex flex-wrap gap-1">
                                ${record.tissues.slice(0, 10).map(tissue => `<span class="badge bg-secondary">${tissue}</span>`).join('')}
                                ${record.tissues.length > 10 ? `<span class="badge bg-secondary">+${record.tissues.length - 10} more</span>` : ''}
                            </div>
                        ` : ''}
                    </div>
                </div>
            </div>
        `;
    });

    html += '</div>';
    resultsDiv.innerHTML = html;
}

function displayReadResults(record) {
    const resultsDiv = document.getElementById('readResults');

    let html = `
        <div class="card mb-2">
            <div class="card-header p-2 bg-info text-white">
                <h6 class="mb-0"><i class="fas fa-sitemap"></i> ${record.cluster_id} - ${record.gene_symbol}</h6>
            </div>
            <div class="card-body p-2">
                <div class="row g-2 mb-2">
                    <div class="col-md-6">
                        <p class="small mb-1"><strong>Cluster ID:</strong> <code>${record.cluster_id}</code></p>
                        <p class="small mb-1"><strong>Gene Symbol:</strong> ${record.gene_symbol}</p>
                        <p class="small mb-1"><strong>Species:</strong> ${record.species}</p>
                        <p class="small mb-1"><strong>Title:</strong> ${record.title}</p>
                    </div>
                    <div class="col-md-6">
                        <p class="small mb-1"><strong>Chromosome:</strong> ${record.chromosome}</p>
                        <p class="small mb-1"><strong>Cytoband:</strong> ${record.cytoband}</p>
                        <p class="small mb-1"><strong>Gene ID:</strong> ${record.gene_id}</p>
                        <p class="small mb-1"><strong>LocusLink:</strong> ${record.locuslink}</p>
                    </div>
                </div>

                ${record.homol !== 'N/A' ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Homology:</strong> ${record.homol}</p>
                ` : ''}

                ${record.restr_expr !== 'N/A' ? `
                    <p class="small mb-1"><strong>Restricted Expression:</strong> ${record.restr_expr}</p>
                ` : ''}

                ${record.gnm_terminus !== 'N/A' ? `
                    <p class="small mb-1"><strong>Genome Terminus:</strong> ${record.gnm_terminus}</p>
                ` : ''}

                ${record.txmap !== 'N/A' ? `
                    <p class="small mb-1"><strong>TxMap:</strong> ${record.txmap}</p>
                ` : ''}

                <hr class="my-2">
                <p class="small mb-1"><strong>Statistics:</strong></p>
                <div class="row g-2">
                    <div class="col-md-4">
                        <div class="text-center p-2 bg-light rounded">
                            <h5 class="text-info mb-0">${record.sequence_count}</h5>
                            <small>Sequences</small>
                        </div>
                    </div>
                    <div class="col-md-4">
                        <div class="text-center p-2 bg-light rounded">
                            <h5 class="text-info mb-0">${record.tissue_count}</h5>
                            <small>Tissues</small>
                        </div>
                    </div>
                    <div class="col-md-4">
                        <div class="text-center p-2 bg-light rounded">
                            <h5 class="text-info mb-0">${record.protein_similarities.length}</h5>
                            <small>Protein Sim.</small>
                        </div>
                    </div>
                </div>

                ${record.tissue_details && record.tissue_details.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Expression Tissues (${record.tissue_details.length}):</strong></p>
                    <div class="table-responsive" style="max-height: 200px; overflow-y: auto;">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light sticky-top">
                                <tr>
                                    <th>Tissue</th>
                                    <th>Frequency</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.tissue_details.map(tissue => `
                                    <tr>
                                        <td>${tissue.tissue}</td>
                                        <td><span class="badge bg-info">${tissue.frequency}</span></td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                ` : ''}

                ${record.sequences && record.sequences.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Sequences (${record.sequences.length}):</strong></p>
                    <div class="table-responsive" style="max-height: 300px; overflow-y: auto;">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light sticky-top">
                                <tr>
                                    <th>Accession</th>
                                    <th>Type</th>
                                    <th>Clone</th>
                                    <th>End</th>
                                    <th>Library ID</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.sequences.slice(0, 50).map(seq => `
                                    <tr>
                                        <td><code>${seq.acc}</code></td>
                                        <td><span class="badge bg-secondary">${seq.seqtype}</span></td>
                                        <td>${seq.clone}</td>
                                        <td>${seq.end}</td>
                                        <td>${seq.lid}</td>
                                    </tr>
                                `).join('')}
                                ${record.sequences.length > 50 ? `
                                    <tr>
                                        <td colspan="5" class="text-center text-muted">
                                            <small>+${record.sequences.length - 50} more sequences</small>
                                        </td>
                                    </tr>
                                ` : ''}
                            </tbody>
                        </table>
                    </div>
                ` : ''}

                ${record.protein_similarities && record.protein_similarities.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Protein Similarities (${record.protein_similarities.length}):</strong></p>
                    <div class="table-responsive" style="max-height: 300px; overflow-y: auto;">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light sticky-top">
                                <tr>
                                    <th>Organism</th>
                                    <th>Protein ID</th>
                                    <th>Identity %</th>
                                    <th>Alignment Length</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.protein_similarities.map(prot => `
                                    <tr>
                                        <td>${prot.organism}</td>
                                        <td><code>${prot.protein_id}</code></td>
                                        <td><span class="badge bg-success">${prot.percent}%</span></td>
                                        <td>${prot.alignment_length}</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                ` : ''}

                ${record.sts && record.sts.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>STS Markers (${record.sts.length}):</strong></p>
                    <div class="table-responsive">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light">
                                <tr>
                                    <th>Accession</th>
                                    <th>UniSTS</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.sts.map(sts => `
                                    <tr>
                                        <td><code>${sts.acc}</code></td>
                                        <td>${sts.unists}</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                ` : ''}
            </div>
        </div>
    `;

    resultsDiv.innerHTML = html;
}

function loadParseExample() {
    const exampleData = `ID          Hs.1
TITLE       actin beta
GENE        ACTB
CYTOBAND    7p22.1
LOCUSLINK   60
CHROMOSOME  7
GENE_ID     60
EXPRESS     brain| 100
EXPRESS     heart| 85
EXPRESS     liver| 95
SCOUNT      2
SEQUENCE    ACC=NM_001101.5; NID=g195972895; CLONE=IMAGE:123456; END=5'; LID=1234; SEQTYPE=mRNA
SEQUENCE    ACC=BC000123.1; NID=g23456789; CLONE=IMAGE:234567; END=3'; LID=2345; SEQTYPE=EST
PROTSIM     ORG=Mus musculus; PROTGI=123456; PCT=98; ALN=375
STS         ACC=G12345; UNISTS=12345
//
ID          Hs.2
TITLE       glyceraldehyde-3-phosphate dehydrogenase
GENE        GAPDH
CYTOBAND    12p13.31
LOCUSLINK   2597
CHROMOSOME  12
GENE_ID     2597
EXPRESS     brain| 90
EXPRESS     liver| 95
EXPRESS     muscle| 88
SCOUNT      2
SEQUENCE    ACC=NM_002046.7; NID=g123456789; CLONE=IMAGE:345678; END=5'; LID=3456; SEQTYPE=mRNA
SEQUENCE    ACC=AI234567.1; NID=g34567890; CLONE=IMAGE:456789; END=3'; LID=4567; SEQTYPE=EST
PROTSIM     ORG=Rattus norvegicus; PROTGI=234567; PCT=97; ALN=335
STS         ACC=G23456; UNISTS=23456
//
ID          Hs.3
TITLE       tumor protein p53
GENE        TP53
CYTOBAND    17p13.1
LOCUSLINK   7157
CHROMOSOME  17
GENE_ID     7157
EXPRESS     brain| 75
EXPRESS     liver| 80
EXPRESS     lung| 70
SCOUNT      1
SEQUENCE    ACC=NM_000546.6; NID=g456789012; CLONE=IMAGE:567890; END=5'; LID=5678; SEQTYPE=mRNA
PROTSIM     ORG=Mus musculus; PROTGI=345678; PCT=85; ALN=393
STS         ACC=G34567; UNISTS=34567
//`;

    const blob = new Blob([exampleData], { type: 'text/plain' });
    const file = new File([blob], 'example_unigene.data', { type: 'text/plain' });

    try {
        const dataTransfer = new DataTransfer();
        dataTransfer.items.add(file);
        document.getElementById('parseFile').files = dataTransfer.files;
    } catch (e) {
        // Fallback: just parse directly
    }

    const formData = new FormData();
    formData.append('file', file);
    formData.append('max_records', document.getElementById('parseMaxRecords').value);

    showLoading('parseBtn');

    fetch('/api/unigene/parse', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');

        if (data.success) {
            displayParseResults(data.records, data.count);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');
        showAlert('Error loading example: ' + error.message, 'danger');
    });
}

function loadReadExample() {
    const exampleData = `ID          Hs.1
TITLE       actin beta
GENE        ACTB
CYTOBAND    7p22.1
LOCUSLINK   60
CHROMOSOME  7
GENE_ID     60
EXPRESS     brain| 100
EXPRESS     heart| 85
EXPRESS     liver| 95
EXPRESS     muscle| 120
EXPRESS     kidney| 78
SCOUNT      2
SEQUENCE    ACC=NM_001101.5; NID=g195972895; CLONE=IMAGE:123456; END=5'; LID=1234; SEQTYPE=mRNA
SEQUENCE    ACC=BC123456.1; NID=g23456789; CLONE=IMAGE:234567; END=3'; LID=2345; SEQTYPE=EST
PROTSIM     ORG=Mus musculus; PROTGI=123456; PCT=98; ALN=375
PROTSIM     ORG=Rattus norvegicus; PROTGI=234567; PCT=97; ALN=375
STS         ACC=G12345; UNISTS=12345
//`;

    const blob = new Blob([exampleData], { type: 'text/plain' });
    const file = new File([blob], 'example_unigene_single.data', { type: 'text/plain' });

    try {
        const dataTransfer = new DataTransfer();
        dataTransfer.items.add(file);
        document.getElementById('readFile').files = dataTransfer.files;
    } catch (e) {
        // Fallback: just parse directly
    }

    const formData = new FormData();
    formData.append('file', file);

    showLoading('readBtn');

    fetch('/api/unigene/read', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Record');

        if (data.success) {
            displayReadResults(data.record);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Record');
        showAlert('Error loading example: ' + error.message, 'danger');
    });
}
