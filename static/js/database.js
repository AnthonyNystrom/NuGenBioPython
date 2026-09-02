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

    if (!term) {
        showAlert('Enter a search term first.', 'warning');
        return;
    }
    if (!email) {
        showAlert('NCBI requires an email address on every Entrez request — this lets them contact you if your usage triggers a rate-limit issue. Enter your email above to continue.', 'warning');
        const emailEl = document.getElementById('email');
        if (emailEl) emailEl.focus();
        return;
    }
    if (!/^\S+@\S+\.\S+$/.test(email)) {
        showAlert('Email looks invalid. Please enter a real email address — NCBI requires it for Entrez requests.', 'warning');
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
            showAlert(friendlyError(data.error, 'ncbi'), 'warning');
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

    const tiles = [
        { label: 'Total',    value: count },
        { label: 'Returned', value: results.length },
        { label: 'Database', value: database },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let html = '<div class="table-responsive"><table class="table table-hover table-sm"><thead><tr><th>ID</th><th>Title</th>';
    if (database === 'pubmed') {
        html += '<th>Authors</th><th>Journal</th><th>Date</th>';
    }
    html += '<th>Actions</th></tr></thead><tbody>';

    results.forEach(result => {
        html += `<tr>
            <td><code>${escapeHtml(result.id)}</code></td>
            <td>${escapeHtml(result.title || 'N/A')}</td>`;
        if (database === 'pubmed') {
            const authors = Array.isArray(result.authors) ? result.authors.slice(0, 2).join(', ') : '';
            html += `<td><small>${escapeHtml(authors)}</small></td>
                    <td><small>${escapeHtml(result.journal || '')}</small></td>
                    <td><small>${escapeHtml(result.date || '')}</small></td>`;
        }
        html += `<td>
                <button class="btn-app-sm btn-app-secondary" data-action="viewFullRecord" data-action-args="[&quot;${escapeHtml(result.id)}&quot;, &quot;${escapeHtml(database)}&quot;]">
                    <i class="fas fa-eye"></i> View
                </button>
            </td>`;
        html += '</tr>';
    });

    html += '</tbody></table></div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('searchResults', {
            title: 'Entrez Search',
            meta: `${database} · ${count} total`,
            summary: tilesHtml + html,
            raw: JSON.stringify(results, null, 2),
            workspaceItem: { type: 'entrez-search', name: `${database} search (${results.length})`, data: results },
            downloads: [
                { label: 'JSON', filename: `${database}-search.json`, text: JSON.stringify(results, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        resultsDiv.innerHTML = tilesHtml + html;
    }
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
            showAlert(friendlyError(data.error, 'ncbi'), 'warning');
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

    const total = results.reduce((s, r) => s + (r.count || 0), 0);
    const tiles = [
        { label: 'Databases', value: results.length },
        { label: 'Total Hits', value: total.toLocaleString() },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    results.forEach(r => { rows += `<tr><td><code>${r.database}</code></td><td><span class="badge bg-primary">${r.count.toLocaleString()}</span></td></tr>`; });
    const table = `<div class="table-responsive"><table class="table table-hover table-sm">
        <thead><tr><th>Database</th><th>Count</th></tr></thead><tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('globalResults', {
            title: 'Global Entrez Search',
            meta: `${results.length} databases · ${total.toLocaleString()} total`,
            summary: tilesHtml + table,
            raw: JSON.stringify(results, null, 2),
            workspaceItem: { type: 'entrez-global', name: `Global search (${total.toLocaleString()})`, data: results },
            downloads: [
                { label: 'JSON', filename: 'entrez-global.json', text: JSON.stringify(results, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        resultsDiv.innerHTML = tilesHtml + table;
    }
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
            showAlert(friendlyError(data.error, 'ncbi'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('fetchBtn', '<i class="fas fa-download me-2"></i>Fetch Records');
    });
});

function displayFetchedRecords(data, format) {
    const resultsDiv = document.getElementById('fetchResults');
    const chars = (data || '').length;

    const tiles = [
        { label: 'Format', value: format },
        { label: 'Size',   value: chars.toLocaleString() + ' chars' },
    ];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('fetchResults', {
            title: 'Fetched Records',
            meta: `${format} · ${chars.toLocaleString()} chars`,
            summary: tilesHtml + `<p class="small mb-0 text-muted">See the Raw tab for the fetched content.</p>`,
            raw: data,
            copyText: data,
            workspaceItem: { type: 'entrez-fetch', name: `Fetched ${format}`, data: data },
            downloads: [
                { label: format.toUpperCase(), filename: `records.${format}`, text: data, mime: 'text/plain' },
            ],
        });
    } else {
        resultsDiv.innerHTML = `<pre class="mb-0" style="max-height: 500px; overflow-y: auto; font-size: 12px;">${escapeHtml(data)}</pre>`;
    }
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
            showAlert(friendlyError(data.error, 'ncbi'), 'warning');
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

    const tiles = [{ label: 'Linked Records', value: count }];
    const tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
        `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
    ).join('') + '</div>';

    let rows = '';
    links.forEach(link => { rows += `<tr><td><code>${link.id}</code></td><td>${link.title || link.name || 'N/A'}</td></tr>`; });
    const table = `<div class="table-responsive"><table class="table table-hover table-sm">
        <thead><tr><th>ID</th><th>Title</th></tr></thead><tbody>${rows}</tbody></table></div>`;

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('linkResults', {
            title: 'Entrez Linked Records',
            meta: `${count} linked record${count === 1 ? '' : 's'}`,
            summary: tilesHtml + table,
            raw: JSON.stringify(links, null, 2),
            workspaceItem: { type: 'entrez-link', name: `Entrez links (${count})`, data: links },
            downloads: [
                { label: 'JSON', filename: 'entrez-link.json', text: JSON.stringify(links, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        resultsDiv.innerHTML = tilesHtml + table;
    }
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
            showAlert(friendlyError(data.error, 'ncbi'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
    });
});

function displayDatabaseInfo(info, database) {
    let title, meta, summary, raw, tilesHtml = '';

    if (!database && Array.isArray(info)) {
        const tiles = [{ label: 'Databases', value: info.length }];
        tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
            `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
        ).join('') + '</div>';
        summary = tilesHtml + '<h6 class="mb-3">Available NCBI Databases</h6><div class="row g-2">' +
            info.map(db => `<div class="col-md-3"><span class="badge bg-secondary">${db.name}</span></div>`).join('') +
            '</div>';
        title = 'NCBI Databases';
        meta = `${info.length} databases`;
        raw = JSON.stringify(info, null, 2);
    } else if (info) {
        const tiles = [
            { label: 'Records',    value: info.count ? parseInt(info.count).toLocaleString() : 'N/A' },
            { label: 'Fields',     value: (info.fields || []).length },
            { label: 'Updated',    value: info.last_update || 'N/A' },
        ];
        tilesHtml = '<div class="rc-stats mb-3">' + tiles.map(t =>
            `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`
        ).join('') + '</div>';
        let body = `<div class="row g-3"><div class="col-md-6"><h6>Database Details</h6><table class="table table-sm">
                <tr><th>Name</th><td>${info.name || 'N/A'}</td></tr>
                <tr><th>Description</th><td>${info.description || 'N/A'}</td></tr>
                <tr><th>Menu Name</th><td>${info.menu_name || 'N/A'}</td></tr>
                <tr><th>Record Count</th><td>${info.count ? parseInt(info.count).toLocaleString() : 'N/A'}</td></tr>
                <tr><th>Last Update</th><td>${info.last_update || 'N/A'}</td></tr>
            </table></div>`;
        if (info.fields && info.fields.length > 0) {
            body += `<div class="col-md-6"><h6>Searchable Fields (first 10)</h6><ul class="small">` +
                info.fields.map(f => `<li><code>${f}</code></li>`).join('') + `</ul></div>`;
        }
        body += '</div>';
        summary = tilesHtml + body;
        title = 'Database Info';
        meta = `${database} · ${info.count ? parseInt(info.count).toLocaleString() : 'N/A'} records`;
        raw = JSON.stringify(info, null, 2);
    } else {
        document.getElementById('infoResults').innerHTML = '<p class="text-muted">No info available.</p>';
        return;
    }

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('infoResults', {
            title: title,
            meta: meta,
            summary: summary,
            raw: raw,
            workspaceItem: { type: 'entrez-info', name: meta, data: info },
            downloads: [
                { label: 'JSON', filename: 'entrez-info.json', text: raw, mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('infoResults').innerHTML = summary;
    }
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
    const raw = `<pre class="code-block" style="max-height: 500px; overflow-y: auto; font-size: 11px; margin: 0;">${escapeHtml(recordData)}</pre>`;
    ResultsPanel.render('recordModalContent', [
        { id: 'raw', title: 'Raw Data', content: raw, active: true },
        { id: 'fmt', title: 'Formatted', content: formatRecordData(recordData, database) },
    ], { prefix: 'record-view' });
}

function formatRecordData(data, database) {
    const F = window.NuGenFormatters || {};
    const looksXml = typeof data === 'string' && data.trim().startsWith('<');
    const looksGenBank = typeof data === 'string' && /^LOCUS\b/m.test(data);
    const looksFasta = typeof data === 'string' && data.trim().startsWith('>');

    if ((database === 'pubmed' || database === 'pmc') && looksXml && F.formatPubMedXML) {
        return '<div class="p-2">' + F.formatPubMedXML(data) + '</div>';
    }
    if ((database === 'nucleotide' || database === 'protein' || database === 'gene') && looksGenBank && F.formatGenBankText) {
        return '<div class="p-2">' + F.formatGenBankText(data) + '</div>';
    }
    if (looksFasta && F.formatFasta) {
        return '<div class="p-2">' + F.formatFasta(data) + '</div>';
    }
    if (looksXml && F.formatXMLPretty) {
        return '<div class="p-2">' + F.formatXMLPretty(data) + '</div>';
    }
    return '<div class="p-3 small font-monospace" style="white-space:pre-wrap;">' + escapeHtml(data) + '</div>';
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
