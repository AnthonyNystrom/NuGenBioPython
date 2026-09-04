// KEGG Database JavaScript - All 5 operations
// showLoading, hideLoading, showAlert are provided globally by static/js/utils.js

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
            showAlert(friendlyError(data.error, 'kegg'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('searchBtn', '<i class="fas fa-search me-2"></i>Search KEGG');
        showAlert(friendlyError(error, 'kegg'), 'warning');
    });
});

function displaySearchResults(results, count) {
    const resultsDiv = document.getElementById('searchResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No results found.</p>';
        return;
    }
    _keggRenderTableCard('searchResults', 'KEGG Search', results, count,
        ['ID', 'Definition', 'Actions'],
        r => [`<code>${r.id}</code>`, r.definition || '', `<button class="btn-app-sm btn-app-secondary" data-action="viewEntry" data-action-args="[&quot;${r.id}&quot;]"><i class="fas fa-eye"></i> View</button>`]);
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
            showAlert(friendlyError(data.error, 'kegg'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('listBtn', '<i class="fas fa-list me-2"></i>List Entries');
        showAlert(friendlyError(error, 'kegg'), 'warning');
    });
});

function displayListResults(results, total, displayed) {
    const resultsDiv = document.getElementById('listResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No entries found.</p>';
        return;
    }
    _keggRenderTableCard('listResults', 'KEGG List', results, displayed,
        ['ID', 'Definition', 'Actions'],
        r => [`<code>${r.id}</code>`, r.definition || 'N/A', `<button class="btn-app-sm btn-app-secondary" data-action="viewEntry" data-action-args="[&quot;${r.id}&quot;]"><i class="fas fa-eye"></i> View</button>`],
        `${displayed} of ${total}`);
}

function _keggRenderTableCard(containerId, title, results, count, headers, rowFn, metaOverride) {
    const tiles = [{ label: 'Entries', value: count }];
    const tilesHtml = '<div class="rc-stats mb-3">' +
        tiles.map(t => `<div class="rc-stat"><div class="rc-stat-value">${t.value}</div><div class="rc-stat-label">${t.label}</div></div>`).join('') +
        '</div>';
    let rows = '';
    results.forEach(r => {
        rows += '<tr>' + rowFn(r).map(c => `<td>${c}</td>`).join('') + '</tr>';
    });
    const table = `<div class="table-responsive"><table class="table table-hover table-sm">
        <thead><tr>${headers.map(h => `<th>${h}</th>`).join('')}</tr></thead><tbody>${rows}</tbody></table></div>`;

    const tsv = [headers.join('\t')];
    results.forEach(r => {
        tsv.push(rowFn(r).map(c => String(c).replace(/<[^>]*>/g, '').replace(/\t/g, ' ').trim()).join('\t'));
    });

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount(containerId, {
            title: title,
            meta: metaOverride || `${count} ${count === 1 ? 'entry' : 'entries'}`,
            summary: tilesHtml + table,
            raw: JSON.stringify(results, null, 2),
            workspaceItem: { type: 'kegg', name: `${title} (${count})`, data: results },
            downloads: [
                { label: 'TSV',  filename: 'kegg.tsv',  text: tsv.join('\n'), mime: 'text/tab-separated-values' },
                { label: 'JSON', filename: 'kegg.json', text: JSON.stringify(results, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById(containerId).innerHTML = tilesHtml + table;
    }
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
            showAlert(friendlyError(data.error, 'kegg'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('linkBtn', '<i class="fas fa-link me-2"></i>Find Links');
        showAlert(friendlyError(error, 'kegg'), 'warning');
    });
});

function displayLinkResults(results, count) {
    const resultsDiv = document.getElementById('linkResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No links found.</p>';
        return;
    }
    _keggRenderTableCard('linkResults', 'KEGG Links', results, count,
        ['Source', 'Target', 'Actions'],
        r => [`<code>${r.source}</code>`, `<code>${r.target}</code>`,
              `<button class="btn-app-sm btn-app-secondary" data-action="viewEntry" data-action-args="[&quot;${r.target}&quot;]"><i class="fas fa-eye"></i></button>`]);
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
            showAlert(friendlyError(data.error, 'kegg'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('convertBtn', '<i class="fas fa-exchange-alt me-2"></i>Convert IDs');
        showAlert(friendlyError(error, 'kegg'), 'warning');
    });
});

function displayConvertResults(results, count) {
    const resultsDiv = document.getElementById('convertResults');
    if (results.length === 0) {
        resultsDiv.innerHTML = '<p class="text-muted">No conversions found.</p>';
        return;
    }
    _keggRenderTableCard('convertResults', 'KEGG Convert', results, count,
        ['Source ID', 'Target ID'],
        r => [`<code>${r.source_id}</code>`, `<code>${r.target_id}</code>`]);
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
            showAlert(friendlyError(data.error, 'kegg'), 'warning');
        }
    })
    .catch(error => {
        hideLoading('infoBtn', '<i class="fas fa-info-circle me-2"></i>Get Database Info');
        showAlert(friendlyError(error, 'kegg'), 'warning');
    });
});

function displayInfoResults(info, database) {
    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('infoResults', {
            title: 'KEGG Database Info',
            meta: database,
            summary: `<p class="small text-muted mb-0">See the Raw tab for the full info record.</p>`,
            raw: info.raw_data,
            copyText: info.raw_data,
            workspaceItem: { type: 'kegg-info', name: `KEGG info ${database}`, data: info },
            downloads: [
                { label: 'Text', filename: `kegg-info-${database}.txt`, text: info.raw_data, mime: 'text/plain' },
            ],
        });
    } else {
        document.getElementById('infoResults').innerHTML =
            `<div class="p-3"><h6>Database: ${escapeHtml(database)}</h6><pre class="border p-3 u-s1cc113d5">${escapeHtml(info.raw_data)}</pre></div>`;
    }
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
            document.getElementById('entryModalContent').innerHTML =
                '<div class="alert alert-danger">' + escapeHtml(friendlyError(data.error, 'kegg')) + '</div>';
        }
    })
    .catch(error => {
        document.getElementById('entryModalContent').innerHTML =
            '<div class="alert alert-danger">' + escapeHtml(friendlyError(error, 'kegg')) + '</div>';
    });
}

function displayEntryInModal(data) {
    const entryId = (document.getElementById('entryModalLabel').textContent || '').replace(/^Entry:\s*/, '');
    let image = '';
    if (data.image_url) {
        image = `<div class="text-center mb-3"><img src="${data.image_url}" class="img-fluid border u-maxw100" alt="Pathway diagram" data-hide-on-error="1"></div>`;
    }

    if (typeof ResultsCard !== 'undefined') {
        ResultsCard.mount('entryModalContent', {
            title: 'KEGG Entry',
            meta: entryId,
            summary: image + '<p class="small mb-0 text-muted">See the Raw tab for the full KEGG record.</p>',
            raw: data.raw_data,
            copyText: data.raw_data,
            workspaceItem: { type: 'kegg-entry', name: `KEGG ${entryId}`, data: data },
            downloads: [
                { label: 'Text', filename: `${entryId.replace(/[^a-z0-9]/gi,'_')}.txt`, text: data.raw_data, mime: 'text/plain' },
                { label: 'JSON', filename: `${entryId.replace(/[^a-z0-9]/gi,'_')}.json`, text: JSON.stringify(data, null, 2), mime: 'application/json' },
            ],
        });
    } else {
        document.getElementById('entryModalContent').innerHTML = image +
            '<pre class="border p-3 u-s8dc56dcd">' + escapeHtml(data.raw_data) + '</pre>';
    }
}
