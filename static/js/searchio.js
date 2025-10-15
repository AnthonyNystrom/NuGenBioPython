// SearchIO Parser - JavaScript for all 7 tabs

// Parse Tab
document.getElementById('parseForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('parseFile');
    const format = document.getElementById('parseFormat').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('format', format);

    showLoading('parseBtn');

    fetch('/api/searchio/parse', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
        if (data.success) {
            displayParseResults(data.results, data.count);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
        showAlert('Error parsing file', 'danger');
    });
});

function displayParseResults(results, count) {
    const resultsDiv = document.getElementById('parseResults');
    const countBadge = document.getElementById('parseCount');

    countBadge.textContent = count + ' queries';
    countBadge.style.display = 'inline-block';

    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted small">No results found</p>';
        return;
    }

    let html = '';
    results.forEach(query => {
        html += `
            <div class="card mb-2">
                <div class="card-header py-1 small">
                    <strong>${query.id}</strong> - ${query.description} (${query.seq_len} bp/aa)
                </div>
                <div class="card-body p-2">
                    <table class="table table-sm small mb-0">
                        <thead><tr><th>Hit ID</th><th>Description</th><th>E-value</th><th>Bit Score</th></tr></thead>
                        <tbody>`;
        query.hits.slice(0, 10).forEach(hit => {
            html += `<tr>
                <td><code class="small">${hit.id}</code></td>
                <td class="small">${hit.description.substring(0, 60)}</td>
                <td class="small">${hit.evalue}</td>
                <td class="small">${hit.bitscore}</td>
            </tr>`;
        });
        html += '</tbody></table></div></div>';
    });
    resultsDiv.innerHTML = html;
}

// Read Tab
document.getElementById('readForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('readFile');
    const format = document.getElementById('readFormat').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('format', format);

    showLoading('readBtn');

    fetch('/api/searchio/read', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Result');
        if (data.success) {
            displayReadResults(data.result);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Result');
        showAlert('Error reading file', 'danger');
    });
});

function displayReadResults(result) {
    const resultsDiv = document.getElementById('readResults');

    let html = `
        <div class="card mb-2">
            <div class="card-header py-1 small bg-primary text-white">
                <strong>${result.id}</strong> - ${result.description}
            </div>
            <div class="card-body p-2 small">
                <p class="mb-1"><strong>Sequence Length:</strong> ${result.seq_len} bp/aa</p>
                <p class="mb-2"><strong>Number of Hits:</strong> ${result.num_hits}</p>
                <strong>Top Hits with HSPs:</strong>
            </div>
        </div>`;

    result.hits.slice(0, 10).forEach(hit => {
        html += `
            <div class="card mb-2">
                <div class="card-header py-1 small">
                    <strong>${hit.id}</strong> - ${hit.description.substring(0, 80)}
                </div>
                <div class="card-body p-2">
                    <p class="small mb-1"><strong>Hit Length:</strong> ${hit.length} | <strong>HSPs:</strong> ${hit.num_hsps}</p>
                    <table class="table table-sm small mb-0">
                        <thead><tr><th>E-value</th><th>Bit Score</th><th>Query Range</th><th>Hit Range</th></tr></thead>
                        <tbody>`;
        hit.hsps.forEach(hsp => {
            html += `<tr>
                <td>${hsp.evalue}</td>
                <td>${hsp.bitscore}</td>
                <td>${hsp.query_start}-${hsp.query_end}</td>
                <td>${hsp.hit_start}-${hsp.hit_end}</td>
            </tr>`;
        });
        html += '</tbody></table></div></div>';
    });

    resultsDiv.innerHTML = html;
}

// Index Tab
document.getElementById('indexForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('indexFile');
    const format = document.getElementById('indexFormat').value;
    const limit = document.getElementById('indexLimit').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('format', format);
    formData.append('limit', limit);

    showLoading('indexBtn');

    fetch('/api/searchio/index', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('indexBtn', '<i class="fas fa-database me-2"></i>Create Index');
        if (data.success) {
            displayIndexResults(data.results);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('indexBtn', '<i class="fas fa-database me-2"></i>Create Index');
        showAlert('Error indexing file', 'danger');
    });
});

function displayIndexResults(results) {
    const resultsDiv = document.getElementById('indexResults');
    const countBadge = document.getElementById('indexCount');

    countBadge.textContent = results.total_queries + ' total queries';
    countBadge.style.display = 'inline-block';

    let html = `<p class="small mb-2"><strong>Total Queries in File:</strong> ${results.total_queries}</p>`;
    html += '<table class="table table-sm small"><thead><tr><th>Query ID</th><th>Description</th><th>Length</th><th>Hits</th></tr></thead><tbody>';

    results.queries.forEach(query => {
        html += `<tr>
            <td><code class="small">${query.id}</code></td>
            <td class="small">${query.description.substring(0, 60)}</td>
            <td class="small">${query.seq_len}</td>
            <td><span class="badge bg-secondary">${query.num_hits}</span></td>
        </tr>`;
    });

    html += '</tbody></table>';
    resultsDiv.innerHTML = html;
}

// Convert Tab
document.getElementById('convertForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('convertFile');
    const inputFormat = document.getElementById('convertInputFormat').value;
    const outputFormat = document.getElementById('convertOutputFormat').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('input_format', inputFormat);
    formData.append('output_format', outputFormat);

    showLoading('convertBtn');

    fetch('/api/searchio/convert', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert Format');
        if (data.success) {
            displayConvertResults(data.result);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert Format');
        showAlert('Error converting file', 'danger');
    });
});

function displayConvertResults(result) {
    const resultsDiv = document.getElementById('convertResults');

    let html = `
        <div class="card mb-2">
            <div class="card-body p-2 small">
                <p class="mb-1"><strong>Queries Converted:</strong> ${result.count}</p>
                <p class="mb-1"><strong>Output Size:</strong> ${result.content_size} bytes</p>
                <strong>Preview:</strong>
                <pre class="bg-light p-2 small mt-1" style="max-height:300px;overflow-y:auto;">${result.content_preview}</pre>
            </div>
        </div>`;

    resultsDiv.innerHTML = html;
}

// Filter Tab
document.getElementById('filterForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('filterFile');
    const format = document.getElementById('filterFormat').value;
    const evalue = document.getElementById('filterEvalue').value;
    const bitscore = document.getElementById('filterBitscore').value;
    const identity = document.getElementById('filterIdentity').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('format', format);
    if (evalue) formData.append('evalue_threshold', evalue);
    if (bitscore) formData.append('bitscore_threshold', bitscore);
    if (identity) formData.append('min_identity', identity);

    showLoading('filterBtn');

    fetch('/api/searchio/filter', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('filterBtn', '<i class="fas fa-filter me-2"></i>Apply Filters');
        if (data.success) {
            displayFilterResults(data.results, data.count);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('filterBtn', '<i class="fas fa-filter me-2"></i>Apply Filters');
        showAlert('Error filtering file', 'danger');
    });
});

function displayFilterResults(results, count) {
    const resultsDiv = document.getElementById('filterResults');
    const countBadge = document.getElementById('filterCount');

    countBadge.textContent = count + ' queries';
    countBadge.style.display = 'inline-block';

    let html = '';
    results.forEach(query => {
        html += `
            <div class="card mb-2">
                <div class="card-header py-1 small">
                    <strong>${query.id}</strong> - ${query.description}
                    <span class="badge bg-success ms-2">${query.filtered_count}/${query.original_count} hits passed</span>
                </div>
                <div class="card-body p-2">
                    <table class="table table-sm small mb-0">
                        <thead><tr><th>Hit ID</th><th>Description</th><th>E-value</th><th>Bit Score</th><th>Identity</th></tr></thead>
                        <tbody>`;
        query.hits.slice(0, 10).forEach(hit => {
            html += `<tr>
                <td><code class="small">${hit.id}</code></td>
                <td class="small">${hit.description.substring(0, 50)}</td>
                <td class="small">${hit.evalue}</td>
                <td class="small">${hit.bitscore}</td>
                <td class="small">${hit.identity}</td>
            </tr>`;
        });
        html += '</tbody></table></div></div>';
    });

    resultsDiv.innerHTML = html;
}

// Write Tab
document.getElementById('writeForm')?.addEventListener('submit', function(e) {
    e.preventDefault();
    const fileInput = document.getElementById('writeFile');
    const inputFormat = document.getElementById('writeInputFormat').value;
    const outputFormat = document.getElementById('writeOutputFormat').value;
    const maxQueries = document.getElementById('writeMaxQueries').value;

    if (!fileInput.files[0]) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('input_format', inputFormat);
    formData.append('output_format', outputFormat);
    formData.append('max_queries', maxQueries);

    showLoading('writeBtn');

    fetch('/api/searchio/write', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('writeBtn', '<i class="fas fa-file-export me-2"></i>Write Output');
        if (data.success) {
            displayWriteResults(data.result);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('writeBtn', '<i class="fas fa-file-export me-2"></i>Write Output');
        showAlert('Error writing file', 'danger');
    });
});

function displayWriteResults(result) {
    const resultsDiv = document.getElementById('writeResults');

    let html = `
        <div class="card mb-2">
            <div class="card-body p-2 small">
                <p class="mb-1"><strong>Queries Written:</strong> ${result.count}</p>
                <p class="mb-1"><strong>Output Size:</strong> ${result.content_size} bytes</p>
                <strong>Content Preview:</strong>
                <pre class="bg-light p-2 small mt-1" style="max-height:400px;overflow-y:auto;">${result.content}</pre>
            </div>
        </div>`;

    resultsDiv.innerHTML = html;
}

// Formats Tab
function loadFormats() {
    showLoading('formatsBtn');

    fetch('/api/searchio/formats', {
        method: 'GET'
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('formatsBtn', '<i class="fas fa-list-alt me-2"></i>Load Supported Formats');
        if (data.success) {
            displayFormats(data.formats);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('formatsBtn', '<i class="fas fa-list-alt me-2"></i>Load Supported Formats');
        showAlert('Error loading formats', 'danger');
    });
}

function displayFormats(formats) {
    const resultsDiv = document.getElementById('formatsResults');

    let html = '';
    for (const [category, formatList] of Object.entries(formats)) {
        html += `<div class="card mb-2">
            <div class="card-header py-1 small bg-primary text-white"><strong>${category}</strong></div>
            <div class="card-body p-2">
                <table class="table table-sm small mb-0">
                    <thead><tr><th>Format ID</th><th>Name</th><th>Description</th><th>Can Write</th></tr></thead>
                    <tbody>`;
        formatList.forEach(fmt => {
            const canWrite = fmt.can_write ? '<span class="badge bg-success">Yes</span>' : '<span class="badge bg-secondary">No</span>';
            html += `<tr>
                <td><code class="small">${fmt.id}</code></td>
                <td class="small">${fmt.name}</td>
                <td class="small">${fmt.description}</td>
                <td>${canWrite}</td>
            </tr>`;
        });
        html += '</tbody></table></div></div>';
    }

    resultsDiv.innerHTML = html;
}

// Example data functions
function loadParseExample() {
    const selectedFormat = document.getElementById('parseFormat').value;

    const examples = {
        'blast-xml': {
            content: `<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.10.0+</BlastOutput_version>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test protein</BlastOutput_query-def>
  <BlastOutput_query-len>50</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test protein</Iteration_query-def>
      <Iteration_query-len>50</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|123456|ref|NP_000001.1|</Hit_id>
          <Hit_def>example protein [Homo sapiens]</Hit_def>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>85.5</Hsp_bit-score>
              <Hsp_evalue>2e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>59</Hsp_hit-to>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>`,
            filename: 'example_blast.xml',
            type: 'text/xml'
        },
        'blast-tab': {
            content: `Query_1\tgi|123456|ref|NP_000001.1|\t90.00\t50\t5\t0\t1\t50\t10\t59\t2e-25\t85.5
Query_1\tgi|789012|ref|NP_000002.1|\t75.00\t48\t12\t0\t3\t50\t15\t62\t5e-18\t65.2`,
            filename: 'example_blast.tab',
            type: 'text/plain'
        },
        'hmmer3-tab': {
            content: `#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                 -          test_profile         -            2.1e-25   85.5   0.0   2.5e-25   85.2   0.0   1.0   1   0   0   1   1   1   1 example sequence 1
seq2                 -          test_profile         -            5.3e-18   65.2   0.1   6.1e-18   64.9   0.1   1.0   1   0   0   1   1   1   1 example sequence 2`,
            filename: 'example_hmmer.tab',
            type: 'text/plain'
        },
        'hmmer3-text': {
            content: `# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       test_profile  [M=100]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence     Description
    ------- ------ -----    ------- ------ -----   ---- --  --------     -----------
    2.1e-25   85.5   0.0    2.5e-25   85.2   0.0    1.0  1  seq1         example sequence 1
    5.3e-18   65.2   0.1    6.1e-18   64.9   0.1    1.0  1  seq2         example sequence 2`,
            filename: 'example_hmmer.txt',
            type: 'text/plain'
        }
    };

    const example = examples[selectedFormat] || examples['blast-xml'];
    const blob = new Blob([example.content], { type: example.type });
    const file = new File([blob], example.filename, { type: example.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('format', selectedFormat);

    showLoading('parseBtn');
    fetch('/api/searchio/parse', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
            if (data.success) displayParseResults(data.results, data.count);
        })
        .catch(error => {
            hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
        });
}

function loadReadExample() {
    // Read is similar to Parse, just use the Parse example logic
    const selectedFormat = document.getElementById('readFormat').value;
    const parseExampleData = getExampleData(selectedFormat);

    const blob = new Blob([parseExampleData.content], { type: parseExampleData.type });
    const file = new File([blob], parseExampleData.filename, { type: parseExampleData.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('format', selectedFormat);

    showLoading('readBtn');
    fetch('/api/searchio/read', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Result');
            if (data.success) displayReadResults(data.result);
        })
        .catch(error => {
            hideLoading('readBtn', '<i class="fas fa-book-open me-2"></i>Read Result');
        });
}

// Helper function to get example data for a format
function getExampleData(format) {
    const examples = {
        'blast-xml': {
            content: `<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.10.0+</BlastOutput_version>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test protein</BlastOutput_query-def>
  <BlastOutput_query-len>50</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test protein</Iteration_query-def>
      <Iteration_query-len>50</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|123456|ref|NP_000001.1|</Hit_id>
          <Hit_def>example protein [Homo sapiens]</Hit_def>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>85.5</Hsp_bit-score>
              <Hsp_evalue>2e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>59</Hsp_hit-to>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>`,
            filename: 'example.xml',
            type: 'text/xml'
        },
        'blast-tab': {
            content: `Query_1\tgi|123456|ref|NP_000001.1|\t90.00\t50\t5\t0\t1\t50\t10\t59\t2e-25\t85.5`,
            filename: 'example.tab',
            type: 'text/plain'
        },
        'hmmer3-tab': {
            content: `seq1                 -          test_profile         -            2.1e-25   85.5   0.0   2.5e-25   85.2   0.0   1.0   1   0   0   1   1   1   1 example sequence 1`,
            filename: 'example.hmmer',
            type: 'text/plain'
        }
    };
    return examples[format] || examples['blast-xml'];
}

function loadIndexExample() {
    const selectedFormat = document.getElementById('indexFormat').value;
    const exampleData = getExampleData(selectedFormat);

    const blob = new Blob([exampleData.content], { type: exampleData.type });
    const file = new File([blob], exampleData.filename, { type: exampleData.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('format', selectedFormat);
    formData.append('limit', '50');

    showLoading('indexBtn');
    fetch('/api/searchio/index', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('indexBtn', '<i class="fas fa-database me-2"></i>Create Index');
            if (data.success) displayIndexResults(data.results);
        })
        .catch(error => {
            hideLoading('indexBtn', '<i class="fas fa-database me-2"></i>Create Index');
        });
}

function loadConvertExample() {
    const inputFormat = document.getElementById('convertInputFormat').value;
    const exampleData = getExampleData(inputFormat);

    const blob = new Blob([exampleData.content], { type: exampleData.type });
    const file = new File([blob], exampleData.filename, { type: exampleData.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('input_format', inputFormat);
    formData.append('output_format', inputFormat); // Same format for safety

    showLoading('convertBtn');
    fetch('/api/searchio/convert', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert Format');
            if (data.success) {
                displayConvertResults(data.result);
            } else {
                document.getElementById('convertResults').innerHTML =
                    `<div class="alert alert-info small mb-0">Note: ${data.error}</div>`;
            }
        })
        .catch(error => {
            hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert Format');
        });
}

function loadFilterExample() {
    const selectedFormat = document.getElementById('filterFormat').value;
    const exampleData = getExampleData(selectedFormat);

    const blob = new Blob([exampleData.content], { type: exampleData.type });
    const file = new File([blob], exampleData.filename, { type: exampleData.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('format', selectedFormat);
    formData.append('evalue_threshold', '1e-20');

    showLoading('filterBtn');
    fetch('/api/searchio/filter', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('filterBtn', '<i class="fas fa-filter me-2"></i>Apply Filters');
            if (data.success) displayFilterResults(data.results, data.count);
        })
        .catch(error => {
            hideLoading('filterBtn', '<i class="fas fa-filter me-2"></i>Apply Filters');
        });
}

function loadWriteExample() {
    const inputFormat = document.getElementById('writeInputFormat').value;

    // Write requires complete BLAST XML - use hardcoded complete version for blast-xml
    let exampleData;
    if (inputFormat === 'blast-xml') {
        exampleData = {
            content: `<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.10.0+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul et al., Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test protein</BlastOutput_query-def>
  <BlastOutput_query-len>50</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test protein</Iteration_query-def>
      <Iteration_query-len>50</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|123456|ref|NP_000001.1|</Hit_id>
          <Hit_def>example protein [Homo sapiens]</Hit_def>
          <Hit_accession>NP_000001</Hit_accession>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>85.5</Hsp_bit-score>
              <Hsp_score>210</Hsp_score>
              <Hsp_evalue>2e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>59</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>45</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
              <Hsp_qseq>MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGET</Hsp_qseq>
              <Hsp_hseq>MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGET</Hsp_hseq>
              <Hsp_midline>MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGET</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>100000000</Statistics_db-num>
          <Statistics_db-len>50000000000</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>`,
            filename: 'example_write.xml',
            type: 'text/xml'
        };
    } else {
        exampleData = getExampleData(inputFormat);
    }

    const blob = new Blob([exampleData.content], { type: exampleData.type });
    const file = new File([blob], exampleData.filename, { type: exampleData.type });
    const formData = new FormData();
    formData.append('file', file);
    formData.append('input_format', inputFormat);
    formData.append('output_format', inputFormat); // Same format
    formData.append('max_queries', '5');

    showLoading('writeBtn');
    fetch('/api/searchio/write', { method: 'POST', body: formData })
        .then(response => response.json())
        .then(data => {
            hideLoading('writeBtn', '<i class="fas fa-file-export me-2"></i>Write Output');
            if (data.success) {
                displayWriteResults(data.result);
            } else {
                document.getElementById('writeResults').innerHTML =
                    `<div class="alert alert-danger small mb-0">Error: ${data.error}</div>`;
            }
        })
        .catch(error => {
            hideLoading('writeBtn', '<i class="fas fa-file-export me-2"></i>Write Output');
            document.getElementById('writeResults').innerHTML =
                `<div class="alert alert-danger small mb-0">Error: ${error.message}</div>`;
        });
}
