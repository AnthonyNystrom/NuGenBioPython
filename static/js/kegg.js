// KEGG Database JavaScript - All 5 operations

// Utility Functions
function showLoading(btnId) {
    const btn = document.getElementById(btnId);
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner-border spinner-border-sm me-2"></span>Loading...';
}

function hideLoading(btnId, originalHtml) {
    const btn = document.getElementById(btnId);
    btn.disabled = false;
    btn.innerHTML = originalHtml;
}

function showAlert(message, type = 'danger') {
    const alertHtml = `<div class="alert alert-${type} alert-dismissible fade show" role="alert">
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    </div>`;
    const container = document.querySelector('.card-body.p-4');
    container.insertAdjacentHTML('afterbegin', alertHtml);
    setTimeout(() => {
        const alert = container.querySelector('.alert');
        if (alert) alert.remove();
    }, 5000);
}

// Tab 1: Search (kegg_find)
document.getElementById('searchForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const database = document.getElementById('searchDatabase').value;
    const query = document.getElementById('searchQuery').value.trim();
    const organism = document.getElementById('searchOrganism').value.trim();
    if (!query) return;
    showLoading('searchBtn');
    fetch('/api/kegg/search', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({database, query, organism: organism || null})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('searchBtn', '<i class="fas fa-search me-2"></i>Search KEGG');
        if (data.success) {
            displaySearchResults(data.results, data.count);
            document.getElementById('searchCount').textContent = data.count + ' results';
            document.getElementById('searchCount').style.display = 'inline-block';
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('searchBtn', '<i class="fas fa-search me-2"></i>Search KEGG');
        showAlert('Error: ' + error);
    });
});

function displaySearchResults(results, count) {
    const resultsDiv = document.getElementById('searchResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No results found.</p>';
        return;
    }
    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>ID</th><th>Definition</th><th>Actions</th></tr></thead><tbody>';
    results.forEach(result => {
        html += `<tr><td><code>${result.id}</code></td><td>${result.definition}</td><td>`;
        html += `<button class="btn btn-sm gradient-btn-info" onclick="viewEntry('${result.id}')">`;
        html += '<i class="fas fa-eye"></i> View</button></td></tr>';
    });
    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function loadSearchExample() {
    document.getElementById('searchDatabase').value = 'pathway';
    document.getElementById('searchQuery').value = 'cancer';
    document.getElementById('searchOrganism').value = '';
}

// Tab 2: List (kegg_list)
document.getElementById('listForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const database = document.getElementById('listDatabase').value;
    const organism = document.getElementById('listOrganism').value.trim();
    const limit = parseInt(document.getElementById('listLimit').value);
    showLoading('listBtn');
    fetch('/api/kegg/list', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({database, organism: organism || null, limit})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('listBtn', '<i class="fas fa-list me-2"></i>List Entries');
        if (data.success) {
            displayListResults(data.results, data.total, data.displayed);
            document.getElementById('listCount').textContent = data.displayed + ' of ' + data.total + ' entries';
            document.getElementById('listCount').style.display = 'inline-block';
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('listBtn', '<i class="fas fa-list me-2"></i>List Entries');
        showAlert('Error: ' + error);
    });
});

function displayListResults(results, total, displayed) {
    const resultsDiv = document.getElementById('listResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No entries found.</p>';
        return;
    }
    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>ID</th><th>Definition</th><th>Actions</th></tr></thead><tbody>';
    results.forEach(result => {
        html += `<tr><td><code>${result.id}</code></td><td>${result.definition || 'N/A'}</td><td>`;
        html += `<button class="btn btn-sm gradient-btn-info" onclick="viewEntry('${result.id}')">`;
        html += '<i class="fas fa-eye"></i> View</button></td></tr>';
    });
    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function loadListExample() {
    document.getElementById('listDatabase').value = 'pathway';
    document.getElementById('listOrganism').value = '';
    document.getElementById('listLimit').value = '50';
}

// Tab 3: Link (kegg_link)
document.getElementById('linkForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const targetDb = document.getElementById('linkTargetDb').value;
    const sourceDb = document.getElementById('linkSourceDb').value;
    const sourceId = document.getElementById('linkSourceId').value.trim();
    showLoading('linkBtn');
    fetch('/api/kegg/link', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({target_db: targetDb, source_db: sourceDb, source_id: sourceId || null})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('linkBtn', '<i class="fas fa-link me-2"></i>Find Links');
        if (data.success) {
            displayLinkResults(data.results, data.count);
            document.getElementById('linkCount').textContent = data.count + ' links';
            document.getElementById('linkCount').style.display = 'inline-block';
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('linkBtn', '<i class="fas fa-link me-2"></i>Find Links');
        showAlert('Error: ' + error);
    });
});

function displayLinkResults(results, count) {
    const resultsDiv = document.getElementById('linkResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No links found.</p>';
        return;
    }
    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>Source</th><th>Target</th><th>Actions</th></tr></thead><tbody>';
    results.forEach(result => {
        html += `<tr><td><code>${result.source}</code></td><td><code>${result.target}</code></td><td>`;
        html += `<button class="btn btn-sm btn-outline-primary btn-sm" onclick="viewEntry('${result.target}')">`;
        html += '<i class="fas fa-eye"></i></button></td></tr>';
    });
    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function loadLinkExample() {
    document.getElementById('linkTargetDb').value = 'pathway';
    document.getElementById('linkSourceDb').value = 'genes';
    document.getElementById('linkSourceId').value = 'hsa:5594';
}

// Tab 4: Convert (kegg_conv)
document.getElementById('convertForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const targetDb = document.getElementById('convertTargetDb').value;
    const sourceDb = document.getElementById('convertSourceDb').value;
    const ids = document.getElementById('convertIds').value.trim();
    showLoading('convertBtn');
    fetch('/api/kegg/convert', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({target_db: targetDb, source_db: sourceDb, ids: ids || null})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert IDs');
        if (data.success) {
            displayConvertResults(data.results, data.count);
            document.getElementById('convertCount').textContent = data.count + ' conversions';
            document.getElementById('convertCount').style.display = 'inline-block';
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert IDs');
        showAlert('Error: ' + error);
    });
});

function displayConvertResults(results, count) {
    const resultsDiv = document.getElementById('convertResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No conversions found.</p>';
        return;
    }
    let html = '<div class="table-responsive"><table class="table table-hover table-sm">';
    html += '<thead><tr><th>Source ID</th><th>Target ID</th></tr></thead><tbody>';
    results.forEach(result => {
        html += `<tr><td><code>${result.source_id}</code></td><td><code>${result.target_id}</code></td></tr>`;
    });
    html += '</tbody></table></div>';
    resultsDiv.innerHTML = html;
}

function loadConvertExample() {
    document.getElementById('convertTargetDb').value = 'ncbi-geneid';
    document.getElementById('convertSourceDb').value = 'hsa';
    document.getElementById('convertIds').value = 'hsa:5594';
}

// Tab 5: Info (kegg_info)
document.getElementById('infoForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const database = document.getElementById('infoDatabase').value;
    showLoading('infoBtn');
    fetch('/api/kegg/info', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({database: database || null})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
        if (data.success) {
            displayInfoResults(data.info, data.database);
        } else {
            showAlert('Error: ' + data.error);
        }
    })
    .catch(error => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
        showAlert('Error: ' + error);
    });
});

function displayInfoResults(info, database) {
    const resultsDiv = document.getElementById('infoResults');
    let html = '<div class="p-3">';
    html += `<h6>Database: ${database}</h6>`;
    html += '<pre class="border p-3" style="max-height: 400px; overflow-y: auto; font-size: 11px; background: #f8f9fa;">';
    html += info.raw_data;
    html += '</pre></div>';
    resultsDiv.innerHTML = html;
}

// View Entry Details Modal
function viewEntry(entryId) {
    const modal = new bootstrap.Modal(document.getElementById('entryModal'));
    document.getElementById('entryModalLabel').textContent = 'Entry: ' + entryId;
    document.getElementById('entryModalContent').innerHTML = '<div class="text-center py-5"><div class="spinner-border text-primary"></div><p class="mt-3">Fetching entry details...</p></div>';
    modal.show();
    fetch('/api/kegg/get/' + entryId)
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            displayEntryInModal(data);
        } else {
            document.getElementById('entryModalContent').innerHTML = '<div class="alert alert-danger">Error: ' + data.error + '</div>';
        }
    })
    .catch(error => {
        document.getElementById('entryModalContent').innerHTML = '<div class="alert alert-danger">Error: ' + error + '</div>';
    });
}

function displayEntryInModal(data) {
    let html = '';
    if (data.image_url) {
        html += '<div class="text-center mb-3"><img src="' + data.image_url + '" class="img-fluid border" alt="Pathway diagram" style="max-width: 100%;" onerror="this.style.display=\'none\'"></div>';
    }
    html += '<h6>Entry Details</h6>';
    html += '<pre class="border p-3" style="max-height: 500px; overflow-y: auto; font-size: 11px; background: #f8f9fa;">';
    html += data.raw_data;
    html += '</pre>';
    document.getElementById('entryModalContent').innerHTML = html;
}
