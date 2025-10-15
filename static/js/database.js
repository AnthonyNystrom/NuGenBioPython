// Database & Search JavaScript

// Search Tab
document.getElementById('entrezForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const database = document.getElementById('database').value;
    const term = document.getElementById('searchTerm').value.trim();
    const email = document.getElementById('email').value.trim();
    const retmax = document.getElementById('retmax').value;
    const dateFrom = document.getElementById('dateFrom').value;
    const dateTo = document.getElementById('dateTo').value;

    if (!term || !email) {
        return;
    }

    // Add date range to search term if provided
    let finalTerm = term;
    if (dateFrom || dateTo) {
        const fromYear = dateFrom ? new Date(dateFrom).getFullYear() : '1900';
        const toYear = dateTo ? new Date(dateTo).getFullYear() : new Date().getFullYear();
        finalTerm += ` AND ${fromYear}:${toYear}[pdat]`;
    }

    const searchData = {
        database: database,
        term: finalTerm,
        email: email,
        retmax: retmax,
        retstart: 0
    };

    showLoading('searchBtn');

    fetch('/api/database/search', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(searchData)
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('searchBtn', '<i class="fas fa-search me-2"></i>Search Database');
        if (data.success) {
            displaySearchResults(data.results, data.count, database);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('searchBtn', '<i class="fas fa-search me-2"></i>Search Database');
    });
});

function displaySearchResults(results, count, database) {
    const resultsDiv = document.getElementById('searchResults');
    const countBadge = document.getElementById('resultCount');

    countBadge.textContent = `${count} total results (showing ${results.length})`;
    countBadge.style.display = 'inline-block';

    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No results found.</p>';
        return;
    }

    let html = '<div class="table-responsive"><table class="table table-hover table-sm"><thead><tr><th>ID</th><th>Title</th>';
    if (database === 'pubmed') {
        html += '<th>Authors</th><th>Journal</th><th>Date</th>';
    }
    html += '<th>Actions</th></tr></thead><tbody>';

    results.forEach(result => {
        html += `<tr>
            <td><code>${result.id}</code></td>
            <td>${result.title || 'N/A'}</td>`;
        if (database === 'pubmed') {
            const authors = Array.isArray(result.authors) ? result.authors.slice(0, 2).join(', ') : '';
            html += `<td><small>${authors}</small></td>
                    <td><small>${result.journal || ''}</small></td>
                    <td><small>${result.date || ''}</small></td>`;
        }
        html += `<td>
                <button class="btn btn-sm gradient-btn-info" onclick="viewFullRecord('${result.id}', '${database}')">
                    <i class="fas fa-eye"></i> View
                </button>
            </td>`;
        html += '</tr>';
    });

    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function loadExampleSearch() {
    document.getElementById('database').value = 'pubmed';
    document.getElementById('searchTerm').value = 'CRISPR gene editing[title]';
    document.getElementById('retmax').value = '20';
}

// Global Query Tab
document.getElementById('globalQueryForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const term = document.getElementById('globalTerm').value.trim();
    const email = document.getElementById('globalEmail').value.trim();

    if (!term || !email) return;

    showLoading('globalQueryBtn');

    fetch('/api/database/global', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({term: term, email: email})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('globalQueryBtn', '<i class="fas fa-search me-2"></i>Search All Databases');
        if (data.success) {
            displayGlobalResults(data.results);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('globalQueryBtn', '<i class="fas fa-search me-2"></i>Search All Databases');
    });
});

function displayGlobalResults(results) {
    const resultsDiv = document.getElementById('globalResults');

    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No results found in any database.</p>';
        return;
    }

    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>Database</th><th>Count</th></tr></thead><tbody>';

    results.forEach(result => {
        html += `<tr>
            <td><code>${result.database}</code></td>
            <td><span class="badge bg-primary">${result.count.toLocaleString()}</span></td>
        </tr>`;
    });

    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

// Fetch Records Tab
document.getElementById('fetchForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const database = document.getElementById('fetchDatabase').value;
    const idsText = document.getElementById('fetchIds').value.trim();
    const format = document.getElementById('fetchFormat').value;
    const email = document.getElementById('fetchEmail').value.trim();

    if (!idsText || !email) return;

    const ids = idsText.split(',').map(id => id.trim()).filter(id => id);

    showLoading('fetchBtn');

    fetch('/api/database/fetch', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            ids: ids,
            database: database,
            rettype: format,
            email: email
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('fetchBtn', '<i class="fas fa-download me-2"></i>Fetch Records');
        if (data.success) {
            displayFetchedRecords(data.data, data.format);
            window.fetchedData = {data: data.data, format: data.format};
            document.getElementById('downloadFetchedBtn').style.display = 'inline-block';
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('fetchBtn', '<i class="fas fa-download me-2"></i>Fetch Records');
    });
});

function displayFetchedRecords(data, format) {
    const resultsDiv = document.getElementById('fetchResults');
    resultsDiv.innerHTML = `<pre class="mb-0" style="max-height: 500px; overflow-y: auto; font-size: 12px;">${data}</pre>`;
}

function downloadFetchedRecords() {
    if (!window.fetchedData) return;

    const blob = new Blob([window.fetchedData.data], {type: 'text/plain'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `records.${window.fetchedData.format}`;
    a.click();
    URL.revokeObjectURL(url);
}

// Related Records Tab
document.getElementById('linkForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const id = document.getElementById('linkId').value.trim();
    const fromDb = document.getElementById('linkFromDb').value;
    const toDb = document.getElementById('linkToDb').value;
    const email = document.getElementById('linkEmail').value.trim();

    if (!id || !email) return;

    showLoading('linkBtn');

    fetch('/api/database/link', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            id: id,
            from_db: fromDb,
            to_db: toDb,
            email: email
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('linkBtn', '<i class="fas fa-link me-2"></i>Find Related Records');
        if (data.success) {
            displayLinkedRecords(data.results, data.count);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('linkBtn', '<i class="fas fa-link me-2"></i>Find Related Records');
    });
});

function displayLinkedRecords(links, count) {
    const resultsDiv = document.getElementById('linkResults');
    const countBadge = document.getElementById('linkCount');

    countBadge.textContent = `${count} linked records`;
    countBadge.style.display = 'inline-block';

    if (links.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No linked records found.</p>';
        return;
    }

    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>ID</th><th>Title</th></tr></thead><tbody>';

    links.forEach(link => {
        html += `<tr>
            <td><code>${link.id}</code></td>
            <td>${link.title || link.name || 'N/A'}</td>
        </tr>`;
    });

    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

// Database Info Tab
document.getElementById('infoForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const database = document.getElementById('infoDatabase').value;
    const email = document.getElementById('infoEmail').value.trim();

    if (!email) return;

    showLoading('infoBtn');

    fetch('/api/database/info', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            database: database,
            email: email
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
        if (data.success) {
            displayDatabaseInfo(data.database || data.databases, database);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
    });
});

function displayDatabaseInfo(info, database) {
    const resultsDiv = document.getElementById('infoResults');

    let html = '';

    if (!database) {
        // Show list of all databases
        if (Array.isArray(info)) {
            html = '<h6 class="mb-3">Available NCBI Databases</h6>';
            html += '<div class="row g-2">';
            info.forEach(db => {
                html += `<div class="col-md-3"><span class="badge bg-secondary">${db.name}</span></div>`;
            });
            html += '</div>';
        }
    } else {
        // Show specific database info
        if (info) {
            html += `<div class="row g-3">`;
            html += `<div class="col-md-6">
                <h6>Database Details</h6>
                <table class="table table-sm">
                    <tr><th>Name</th><td>${info.name || 'N/A'}</td></tr>
                    <tr><th>Description</th><td>${info.description || 'N/A'}</td></tr>
                    <tr><th>Menu Name</th><td>${info.menu_name || 'N/A'}</td></tr>
                    <tr><th>Record Count</th><td>${info.count ? parseInt(info.count).toLocaleString() : 'N/A'}</td></tr>
                    <tr><th>Last Update</th><td>${info.last_update || 'N/A'}</td></tr>
                </table>
            </div>`;

            if (info.fields && info.fields.length > 0) {
                html += `<div class="col-md-6">
                    <h6>Searchable Fields (first 10)</h6>
                    <ul class="small">`;
                info.fields.forEach(field => {
                    html += `<li><code>${field}</code></li>`;
                });
                html += `</ul></div>`;
            }
            html += `</div>`;
        }
    }

    resultsDiv.innerHTML = html;
}

// View Full Record Modal
let currentRecordData = null;

function viewFullRecord(recordId, database) {
    const email = document.getElementById('email').value || 'test@example.com';

    // Show modal with loading spinner
    const modal = new bootstrap.Modal(document.getElementById('recordModal'));
    document.getElementById('recordModalLabel').textContent = `Record ${recordId} - ${database}`;
    document.getElementById('recordModalContent').innerHTML = `
        <div class="text-center py-5">
            <div class="spinner-border text-primary" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
            <p class="mt-3">Fetching full record...</p>
        </div>`;
    modal.show();

    // Determine format based on database
    const formatMap = {
        'nucleotide': 'gb',
        'protein': 'gp',
        'gene': 'xml',
        'pubmed': 'xml',
        'pmc': 'xml'
    };
    const format = formatMap[database] || 'text';

    // Fetch full record
    fetch('/api/database/fetch', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            ids: [recordId],
            database: database,
            rettype: format,
            email: email
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            currentRecordData = {
                id: recordId,
                database: database,
                format: data.format,
                data: data.data
            };
            displayRecordInModal(data.data, database, recordId);
        } else {
            document.getElementById('recordModalContent').innerHTML = `
                <div class="alert alert-danger">
                    <i class="fas fa-exclamation-triangle"></i> Error: ${data.error}
                </div>`;
        }
    })
    .catch(error => {
        document.getElementById('recordModalContent').innerHTML = `
            <div class="alert alert-danger">
                <i class="fas fa-exclamation-triangle"></i> Error fetching record: ${error}
            </div>`;
    });
}

function displayRecordInModal(recordData, database, recordId) {
    let html = '';

    // Create tabbed interface for different views
    html += `
        <ul class="nav nav-tabs mb-3" role="tablist">
            <li class="nav-item">
                <button class="nav-link active" data-bs-toggle="tab" data-bs-target="#raw-view">Raw Data</button>
            </li>
            <li class="nav-item">
                <button class="nav-link" data-bs-toggle="tab" data-bs-target="#formatted-view">Formatted</button>
            </li>
        </ul>
        <div class="tab-content">
            <div class="tab-pane fade show active" id="raw-view">
                <pre class="border p-3" style="max-height: 500px; overflow-y: auto; font-size: 11px; background: #f8f9fa;">${escapeHtml(recordData)}</pre>
            </div>
            <div class="tab-pane fade" id="formatted-view">
                ${formatRecordData(recordData, database)}
            </div>
        </div>`;

    document.getElementById('recordModalContent').innerHTML = html;
}

function formatRecordData(data, database) {
    // Extract key information from the record
    const lines = data.split('\n');
    let html = '<div class="p-3">';

    if (database === 'pubmed') {
        html += '<h6>PubMed Record</h6>';
        html += '<div class="small">';
        lines.slice(0, 30).forEach(line => {
            if (line.trim()) {
                html += `<p class="mb-1">${escapeHtml(line)}</p>`;
            }
        });
        html += '</div>';
    } else if (database === 'nucleotide' || database === 'protein') {
        html += `<h6>${database === 'nucleotide' ? 'Nucleotide' : 'Protein'} Record</h6>`;
        html += '<div class="small">';
        lines.slice(0, 50).forEach(line => {
            if (line.trim()) {
                html += `<p class="mb-1 font-monospace">${escapeHtml(line)}</p>`;
            }
        });
        html += '</div>';
    } else {
        html += '<div class="small"><pre>' + escapeHtml(data) + '</pre></div>';
    }

    html += '</div>';
    return html;
}

function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// Download Record
document.getElementById('downloadRecordBtn').addEventListener('click', function() {
    if (!currentRecordData) return;

    const { id, database, format, data } = currentRecordData;
    const blob = new Blob([data], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${database}_${id}.${format}`;
    a.click();
    URL.revokeObjectURL(url);
});
