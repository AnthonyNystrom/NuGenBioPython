/**
 * Advanced Tools JavaScript
 * Handles HMM, Literature Mining, Nexus, SCOP, and Codon Alignment
 */

// ============================================================================
// HMM FUNCTIONS
// ============================================================================

/**
 * Build HMM model
 */
function buildHMM() {
    const hmmType = document.getElementById('hmmType').value;
    const states = parseInt(document.getElementById('states').value);
    const sequence = document.getElementById('sequence').value.trim();

    if (!sequence) {
        showAlert('Please enter a training sequence', 'warning');
        return;
    }

    showLoading('buildHmmBtn');

    fetch('/api/advanced/hmm/build', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            type: hmmType,
            states: states,
            sequence: sequence
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('buildHmmBtn', '<i class="fas fa-cogs"></i> Build HMM');

        if (data.success) {
            displayHmmResults(data.model);
            window.currentHMM = data.model;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('buildHmmBtn', '<i class="fas fa-cogs"></i> Build HMM');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Train HMM using Baum-Welch
 */
function trainHMM() {
    if (!window.currentHMM) {
        showAlert('Please build an HMM first', 'warning');
        return;
    }

    const sequence = document.getElementById('sequence').value.trim();
    const iterations = parseInt(document.getElementById('trainIterations').value) || 10;

    showLoading('trainHmmBtn');

    fetch('/api/advanced/hmm/train', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: sequence,
            iterations: iterations,
            model_id: window.currentHMM.model_id
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('trainHmmBtn', '<i class="fas fa-graduation-cap"></i> Train with Baum-Welch');

        if (data.success) {
            displayTrainingResults(data.training);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('trainHmmBtn', '<i class="fas fa-graduation-cap"></i> Train with Baum-Welch');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Decode sequence using Viterbi
 */
function decodeSequence() {
    if (!window.currentHMM) {
        showAlert('Please build an HMM first', 'warning');
        return;
    }

    const testSequence = document.getElementById('testSequence').value.trim();

    if (!testSequence) {
        showAlert('Please enter a test sequence', 'warning');
        return;
    }

    showLoading('decodeBtn');

    fetch('/api/advanced/hmm/decode', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequence: testSequence,
            model_id: window.currentHMM.model_id
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('decodeBtn', '<i class="fas fa-search"></i> Decode with Viterbi');

        if (data.success) {
            displayDecodingResults(data.decoding);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('decodeBtn', '<i class="fas fa-search"></i> Decode with Viterbi');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Display HMM build results
 */
function displayHmmResults(model) {
    const vizDiv = document.getElementById('hmmViz');
    const analysisDiv = document.getElementById('hmmAnalysis');

    // Visualization
    let vizHtml = `
        <div class="text-center mb-3">
            <h6><i class="fas fa-project-diagram"></i> HMM State Diagram</h6>
            <div class="d-flex justify-content-center align-items-center flex-wrap gap-3 p-3">
    `;

    for (let i = 0; i < model.states; i++) {
        vizHtml += `
            <div class="text-center">
                <div class="rounded-circle bg-primary text-white d-flex align-items-center justify-content-center"
                     style="width: 70px; height: 70px; font-weight: bold; font-size: 1.1rem;">
                    S${i + 1}
                </div>
                <small class="text-muted mt-1 d-block">State ${i + 1}</small>
            </div>
        `;
        if (i < model.states - 1) {
            vizHtml += '<i class="fas fa-arrow-right text-primary fa-2x"></i>';
        }
    }

    vizHtml += `
            </div>
            <p class="text-muted small">Model ID: ${model.model_id}</p>
        </div>
    `;
    vizDiv.innerHTML = vizHtml;

    // Analysis
    let analysisHtml = `
        <div class="row g-3">
            <div class="col-md-6">
                <h6><i class="fas fa-info-circle"></i> Model Parameters</h6>
                <table class="table table-sm table-bordered">
                    <tbody>
                        <tr><th>Type</th><td>${model.type}</td></tr>
                        <tr><th>States</th><td>${model.states}</td></tr>
                        <tr><th>Sequence Length</th><td>${model.sequence_length}</td></tr>
                        <tr><th>Alphabet</th><td>${model.alphabet}</td></tr>
                        <tr><th>Emissions</th><td>${model.emissions_possible}</td></tr>
                        <tr><th>Transitions</th><td>${model.transitions_possible}</td></tr>
                    </tbody>
                </table>
            </div>
            <div class="col-md-6">
                <h6><i class="fas fa-chart-bar"></i> Model Statistics</h6>
                <table class="table table-sm table-bordered">
                    <tbody>
                        <tr><th>Build Time</th><td>${model.training_time}ms</td></tr>
                        <tr><th>Model Created</th><td><span class="badge bg-success">Yes</span></td></tr>
                        <tr><th>Status</th><td><span class="badge bg-info">Ready for Training</span></td></tr>
                    </tbody>
                </table>
            </div>
        </div>
    `;

    analysisDiv.innerHTML = analysisHtml;

    // Show training and decoding sections
    document.getElementById('hmmTraining').style.display = 'block';
    document.getElementById('hmmDecoding').style.display = 'block';
}

/**
 * Display training results
 */
function displayTrainingResults(training) {
    const resultsDiv = document.getElementById('trainingResults');

    let html = `
        <div class="alert alert-success">
            <h6><i class="fas fa-graduation-cap"></i> Training Complete</h6>
            <hr>
            <div class="row">
                <div class="col-md-6">
                    <p class="mb-1"><strong>Algorithm:</strong> Baum-Welch</p>
                    <p class="mb-1"><strong>Iterations:</strong> ${training.iterations}</p>
                    <p class="mb-1"><strong>Converged:</strong> ${training.converged ? 'Yes' : 'No'}</p>
                    <p class="mb-1"><strong>Training Time:</strong> ${training.training_time}ms</p>
                </div>
                <div class="col-md-6">
                    <p class="mb-1"><strong>Final Log-Likelihood:</strong> ${training.final_log_likelihood.toFixed(4)}</p>
                    <p class="mb-1"><strong>Initial Log-Likelihood:</strong> ${training.initial_log_likelihood.toFixed(4)}</p>
                    <p class="mb-1"><strong>Improvement:</strong> ${training.improvement.toFixed(4)}</p>
                </div>
            </div>
        </div>
    `;

    if (training.iteration_history && training.iteration_history.length > 0) {
        html += `
            <div class="mt-3">
                <h6>Convergence History:</h6>
                <div class="table-responsive">
                    <table class="table table-sm table-bordered">
                        <thead class="table-light">
                            <tr>
                                <th>Iteration</th>
                                <th>Log-Likelihood</th>
                                <th>Improvement</th>
                            </tr>
                        </thead>
                        <tbody>
        `;

        training.iteration_history.forEach((item, idx) => {
            html += `
                <tr>
                    <td>${idx + 1}</td>
                    <td>${item.log_likelihood.toFixed(4)}</td>
                    <td>${item.improvement ? item.improvement.toFixed(4) : '-'}</td>
                </tr>
            `;
        });

        html += `
                        </tbody>
                    </table>
                </div>
            </div>
        `;
    }

    resultsDiv.innerHTML = html;
}

/**
 * Display Viterbi decoding results
 */
function displayDecodingResults(decoding) {
    const resultsDiv = document.getElementById('decodingResults');

    let html = `
        <div class="alert alert-success">
            <h6><i class="fas fa-search"></i> Viterbi Decoding Complete</h6>
            <hr>
            <p class="mb-2"><strong>Test Sequence Length:</strong> ${decoding.sequence_length}</p>
            <p class="mb-2"><strong>Log Probability:</strong> ${decoding.log_probability.toFixed(4)}</p>
            <p class="mb-3"><strong>Decoding Time:</strong> ${decoding.decoding_time}ms</p>

            <div class="mt-3">
                <h6>Predicted State Path:</h6>
                <div class="bg-light p-3 rounded">
                    <code style="word-break: break-all;">${decoding.state_path}</code>
                </div>
            </div>
        </div>
    `;

    if (decoding.state_sequence && decoding.state_sequence.length > 0) {
        html += `
            <div class="mt-3">
                <h6>Detailed State Alignment:</h6>
                <div class="table-responsive">
                    <table class="table table-sm table-bordered font-monospace small">
                        <thead class="table-light">
                            <tr>
                                <th>Position</th>
                                <th>Symbol</th>
                                <th>State</th>
                            </tr>
                        </thead>
                        <tbody>
        `;

        const limit = Math.min(decoding.state_sequence.length, 50);
        for (let i = 0; i < limit; i++) {
            const item = decoding.state_sequence[i];
            html += `
                <tr>
                    <td>${i + 1}</td>
                    <td class="text-center">${item.symbol}</td>
                    <td>${item.state}</td>
                </tr>
            `;
        }

        if (decoding.state_sequence.length > 50) {
            html += `<tr><td colspan="3" class="text-center text-muted">... (showing first 50 of ${decoding.state_sequence.length})</td></tr>`;
        }

        html += `
                        </tbody>
                    </table>
                </div>
            </div>
        `;
    }

    resultsDiv.innerHTML = html;
}

/**
 * Load HMM example
 */
function loadHmmExample() {
    const hmmType = document.getElementById('hmmType').value;

    switch(hmmType) {
        case 'sequence':
            document.getElementById('states').value = '3';
            document.getElementById('sequence').value = 'ATCGATCGATCGATCGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGG';
            document.getElementById('testSequence').value = 'ATCGATCGAAATTTCCCGGG';
            break;

        case 'emission':
            document.getElementById('states').value = '2';
            document.getElementById('sequence').value = 'HHHCCCHHCCHHHCCCHHHHCCCCHHHHCCCCHHHHHCCCCC';
            document.getElementById('testSequence').value = 'HHHCCCHHCC';
            break;

        case 'transition':
            document.getElementById('states').value = '2';
            document.getElementById('sequence').value = 'UUUUDDDUUUDDUUUUUDDDDDDUUUUUDDDDUUUUUDDDD';
            document.getElementById('testSequence').value = 'UUUDDD';
            break;

        case 'profile':
            document.getElementById('states').value = '5';
            document.getElementById('sequence').value = 'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEF';
            document.getElementById('testSequence').value = 'ACDEFGHIKLMNPQRSTV';
            break;
    }

    document.getElementById('trainIterations').value = '10';
}

/**
 * Update HMM description based on type
 */
function updateHmmDescription() {
    const hmmType = document.getElementById('hmmType').value;
    const descriptionSpan = document.getElementById('hmmDescription');
    const sequenceTextarea = document.getElementById('sequence');

    const descriptions = {
        'sequence': 'DNA/RNA sequence with repeating patterns for pattern recognition',
        'emission': 'Observation sequence for emission probability modeling (e.g., H=Hot, C=Cold)',
        'transition': 'State transition sequence for market/system modeling (e.g., U=Up, D=Down)',
        'profile': 'Protein sequence for profile HMM to identify conserved domains'
    };

    const placeholders = {
        'sequence': 'Enter DNA/RNA sequence (ATCG)...',
        'emission': 'Enter observation sequence (e.g., HHHCCCHHH for weather states)...',
        'transition': 'Enter state sequence (e.g., UUUDDD for market movements)...',
        'profile': 'Enter protein sequence (single letter amino acids)...'
    };

    descriptionSpan.textContent = descriptions[hmmType] || descriptions['sequence'];
    sequenceTextarea.placeholder = placeholders[hmmType] || placeholders['sequence'];
}

// ============================================================================
// LITERATURE MINING FUNCTIONS
// ============================================================================

/**
 * Search literature
 */
function searchLiterature() {
    const query = document.getElementById('pubmedQuery').value.trim();
    const maxResults = parseInt(document.getElementById('maxResults').value);
    const sortBy = document.getElementById('sortBy').value;

    if (!query) {
        showAlert('Please enter a search query', 'warning');
        return;
    }

    showLoading('searchLiteratureBtn');

    fetch('/api/advanced/literature/search', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            query: query,
            max_results: maxResults,
            sort_by: sortBy
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('searchLiteratureBtn', '<i class="fas fa-search"></i> Search Literature');

        if (data.success) {
            displayLiteratureResults(data.results);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('searchLiteratureBtn', '<i class="fas fa-search"></i> Search Literature');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Display literature search results
 */
function displayLiteratureResults(results) {
    const resultsDiv = document.getElementById('literatureResults');

    if (!results || results.length === 0) {
        resultsDiv.innerHTML = '<div class="alert alert-info">No results found</div>';
        return;
    }

    let html = `
        <div class="mb-3">
            <span class="badge bg-primary">${results.length} articles found</span>
        </div>
    `;

    results.forEach((article, idx) => {
        const collapseId = `article${idx}`;
        const authorsStr = Array.isArray(article.authors) ? article.authors.slice(0, 3).join(', ') : article.authors;
        const shortAuthors = authorsStr.length > 100 ? authorsStr.substring(0, 100) + '...' : authorsStr;

        html += `
            <div class="card mb-2">
                <div class="card-body p-3">
                    <div class="d-flex justify-content-between align-items-start mb-2">
                        <h6 class="mb-0 flex-grow-1">
                            <a href="#" class="text-decoration-none text-dark"
                               data-bs-toggle="collapse"
                               data-bs-target="#${collapseId}">
                                ${article.title}
                            </a>
                        </h6>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/${article.pmid}/"
                           target="_blank"
                           class="btn btn-sm btn-outline-primary ms-2 flex-shrink-0">
                            <i class="fas fa-external-link-alt"></i>
                        </a>
                    </div>

                    <div class="text-muted small mb-2">
                        <strong>Authors:</strong> ${shortAuthors}
                    </div>

                    <div class="text-muted small mb-2">
                        ${article.journal} ${article.year ? `(${article.year})` : ''} • PMID: ${article.pmid}
                    </div>

                    <div class="collapse" id="${collapseId}">
                        <hr class="my-2">
                        <div class="small">
                            <strong>Abstract:</strong>
                            <p class="mb-2 mt-1">${article.abstract || 'No abstract available'}</p>
                        </div>
                        <div class="d-flex gap-2 flex-wrap">
                            <a href="https://pubmed.ncbi.nlm.nih.gov/${article.pmid}/"
                               target="_blank"
                               class="btn btn-sm btn-outline-primary">
                                <i class="fas fa-book me-1"></i> PubMed
                            </a>
                            <a href="https://scholar.google.com/scholar?q=${encodeURIComponent(article.title)}"
                               target="_blank"
                               class="btn btn-sm btn-outline-secondary">
                                <i class="fas fa-graduation-cap me-1"></i> Scholar
                            </a>
                        </div>
                    </div>

                    <div class="text-center mt-2">
                        <a href="#" class="text-decoration-none small"
                           data-bs-toggle="collapse"
                           data-bs-target="#${collapseId}">
                            <i class="fas fa-chevron-down"></i> Toggle details
                        </a>
                    </div>
                </div>
            </div>
        `;
    });

    resultsDiv.innerHTML = html;
}

/**
 * Load literature example
 */
function loadLiteratureExample() {
    document.getElementById('pubmedQuery').value = 'CRISPR gene editing';
    document.getElementById('maxResults').value = '10';
    document.getElementById('sortBy').value = 'relevance';
}

// ============================================================================
// NEXUS FORMAT FUNCTIONS
// ============================================================================

/**
 * Parse Nexus file/data
 */
function parseNexus() {
    const nexusFile = document.getElementById('nexusFile').files[0];
    const nexusData = document.getElementById('nexusData').value.trim();

    if (!nexusFile && !nexusData) {
        showAlert('Please upload a file or paste Nexus data', 'warning');
        return;
    }

    showLoading('parseNexusBtn');

    const formData = new FormData();
    if (nexusFile) {
        formData.append('file', nexusFile);
    } else {
        formData.append('data', nexusData);
    }

    fetch('/api/advanced/nexus/parse', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('parseNexusBtn', '<i class="fas fa-play"></i> Parse Nexus');

        if (data.success) {
            displayNexusResults(data.parsed);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseNexusBtn', '<i class="fas fa-play"></i> Parse Nexus');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Display Nexus parsing results
 */
function displayNexusResults(parsed) {
    const resultsDiv = document.getElementById('nexusResults');

    if (!parsed) {
        resultsDiv.innerHTML = '<div class="alert alert-warning">No data parsed</div>';
        return;
    }

    let html = `
        <div class="alert alert-success">
            <h6><i class="fas fa-check-circle"></i> Nexus File Parsed Successfully</h6>
            <hr>
            <div class="row">
                <div class="col-md-6">
                    <p class="mb-1"><strong>Number of Taxa:</strong> ${parsed.ntax}</p>
                    <p class="mb-1"><strong>Number of Characters:</strong> ${parsed.nchar}</p>
                    <p class="mb-1"><strong>Data Type:</strong> ${parsed.datatype}</p>
                </div>
                <div class="col-md-6">
                    <p class="mb-1"><strong>Has Trees:</strong> ${parsed.has_trees ? 'Yes' : 'No'}</p>
                    <p class="mb-1"><strong>Has Alignments:</strong> ${parsed.has_matrix ? 'Yes' : 'No'}</p>
                </div>
            </div>
        </div>
    `;

    if (parsed.taxlabels && parsed.taxlabels.length > 0) {
        html += `
            <div class="mt-3">
                <h6>Taxa:</h6>
                <div class="d-flex flex-wrap gap-2">
                    ${parsed.taxlabels.map(tax => `<span class="badge bg-secondary">${tax}</span>`).join('')}
                </div>
            </div>
        `;
    }

    if (parsed.matrix && Object.keys(parsed.matrix).length > 0) {
        html += `
            <div class="mt-3">
                <h6>Alignment Preview:</h6>
                <div class="table-responsive">
                    <table class="table table-sm table-bordered font-monospace small">
                        <tbody>
        `;

        for (const [taxon, sequence] of Object.entries(parsed.matrix)) {
            const shortSeq = sequence.length > 50 ? sequence.substring(0, 50) + '...' : sequence;
            html += `
                <tr>
                    <td><strong>${taxon}</strong></td>
                    <td>${shortSeq}</td>
                </tr>
            `;
        }

        html += `
                        </tbody>
                    </table>
                </div>
            </div>
        `;
    }

    resultsDiv.innerHTML = html;
}

/**
 * Load Nexus example
 */
function loadNexusExample() {
    const exampleNexus = `#NEXUS
BEGIN DATA;
DIMENSIONS NTAX=4 NCHAR=15;
FORMAT DATATYPE=DNA MISSING=? GAP=-;
MATRIX
Species1   ATGCTAGCTAGCTAG
Species2   ATGCTAGGTAGCTAG
Species3   ATGCTAGCAAGCTAG
Species4   ATGCTAGCTAGGTAG
;
END;`;

    document.getElementById('nexusData').value = exampleNexus;
}

// ============================================================================
// SCOP CLASSIFICATION FUNCTIONS
// ============================================================================

/**
 * Lookup SCOP classification
 */
function lookupSCOP() {
    const scopId = document.getElementById('scopId').value.trim();
    const scopLevel = document.getElementById('scopLevel').value;

    if (!scopId) {
        showAlert('Please enter a SCOP ID', 'warning');
        return;
    }

    showLoading('lookupScopBtn');

    fetch('/api/advanced/scop/lookup', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            scop_id: scopId,
            level: scopLevel
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('lookupScopBtn', '<i class="fas fa-search"></i> Lookup SCOP');

        if (data.success) {
            displayScopResults(data.classification);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('lookupScopBtn', '<i class="fas fa-search"></i> Lookup SCOP');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Display SCOP classification results
 */
function displayScopResults(classification) {
    const resultsDiv = document.getElementById('scopResults');

    if (!classification) {
        resultsDiv.innerHTML = '<div class="alert alert-warning">No classification data</div>';
        return;
    }

    let html = `
        <div class="alert alert-success">
            <h6><i class="fas fa-sitemap"></i> SCOP Classification Hierarchy</h6>
            <hr>
            <div class="row">
                <div class="col-md-6">
                    <p class="mb-1"><strong>SCOP ID:</strong> <code>${classification.scop_id}</code></p>
                    <p class="mb-1"><strong>Class:</strong> ${classification.class_name || 'N/A'}</p>
                    <p class="mb-1"><strong>Fold:</strong> ${classification.fold_name || 'N/A'}</p>
                    <p class="mb-1"><strong>Superfamily:</strong> ${classification.superfamily_name || 'N/A'}</p>
                </div>
                <div class="col-md-6">
                    <p class="mb-1"><strong>Family:</strong> ${classification.family_name || 'N/A'}</p>
                    <p class="mb-1"><strong>Protein:</strong> ${classification.protein_name || 'N/A'}</p>
                    <p class="mb-1"><strong>Species:</strong> ${classification.species_name || 'N/A'}</p>
                </div>
            </div>
        </div>
    `;

    if (classification.description) {
        html += `
            <div class="mt-3">
                <h6>Description:</h6>
                <p class="text-muted">${classification.description}</p>
            </div>
        `;
    }

    resultsDiv.innerHTML = html;
}

/**
 * Load SCOP example
 */
function loadScopExample() {
    document.getElementById('scopId').value = 'd1dlwa_';
    document.getElementById('scopLevel').value = 'domain';
}

// ============================================================================
// CODON ALIGNMENT FUNCTIONS
// ============================================================================

/**
 * Align codons
 */
function alignCodons() {
    const sequences = document.getElementById('codonSequences').value.trim();
    const geneticCode = document.getElementById('geneticCode').value;
    const alignMethod = document.getElementById('alignMethod').value;

    if (!sequences) {
        showAlert('Please enter sequences', 'warning');
        return;
    }

    showLoading('alignCodonBtn');

    fetch('/api/advanced/codon/align', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sequences: sequences,
            genetic_code: geneticCode,
            method: alignMethod
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('alignCodonBtn', '<i class="fas fa-play"></i> Align Codons');

        if (data.success) {
            displayCodonResults(data.alignment);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('alignCodonBtn', '<i class="fas fa-play"></i> Align Codons');
        showAlert('Network error: ' + error.message, 'danger');
    });
}

/**
 * Display codon alignment results
 */
function displayCodonResults(alignment) {
    const resultsDiv = document.getElementById('codonResults');

    if (!alignment) {
        resultsDiv.innerHTML = '<div class="alert alert-warning">No alignment data</div>';
        return;
    }

    let html = `
        <div class="alert alert-success">
            <h6><i class="fas fa-align-left"></i> Codon Alignment Results</h6>
            <hr>
            <div class="row">
                <div class="col-md-6">
                    <p class="mb-1"><strong>Sequences Aligned:</strong> ${alignment.aligned_count}</p>
                    <p class="mb-1"><strong>Alignment Length:</strong> ${alignment.alignment_length} codons (${alignment.alignment_length * 3} bp)</p>
                </div>
                <div class="col-md-6">
                    <p class="mb-1"><strong>Method:</strong> ${alignment.method}</p>
                    <p class="mb-1"><strong>Genetic Code:</strong> Table ${alignment.genetic_code}</p>
                </div>
            </div>
        </div>
    `;

    if (alignment.alignment && alignment.alignment.length > 0) {
        html += `
            <div class="mt-3">
                <h6>Alignment Preview:</h6>
                <div class="table-responsive">
                    <table class="table table-sm table-bordered font-monospace small">
                        <tbody>
        `;

        alignment.alignment.forEach(seq => {
            const shortSeq = seq.sequence.length > 60 ? seq.sequence.substring(0, 60) + '...' : seq.sequence;
            html += `
                <tr>
                    <td width="150"><strong>${seq.id}</strong></td>
                    <td>${shortSeq}</td>
                </tr>
            `;
        });

        html += `
                        </tbody>
                    </table>
                </div>
            </div>
        `;
    }

    if (alignment.note) {
        html += `<p class="text-muted small mt-2"><em>${alignment.note}</em></p>`;
    }

    resultsDiv.innerHTML = html;
}

/**
 * Load codon alignment example
 */
function loadCodonExample() {
    const exampleSequences = `>Sequence1
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>Sequence2
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>Sequence3
ATGGCGATTGTAATGGGCCGTTGAAAGGGTGCCCGATAG`;

    document.getElementById('codonSequences').value = exampleSequences;
    document.getElementById('geneticCode').value = '1';
    document.getElementById('alignMethod').value = 'default';
}

// ============================================================================
// Initialize on page load
// ============================================================================

document.addEventListener('DOMContentLoaded', function() {
    // Initialize HMM description
    if (document.getElementById('hmmType')) {
        updateHmmDescription();
    }
});
