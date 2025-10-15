/**
 * Pathway Analysis - Complete JavaScript
 * Handles Bio.Pathway: Reaction, System, Network, MultiGraph
 */

// Global state
let pathwayReactions = [];
let currentSystem = null;
let currentNetwork = null;

// ============================================================================
// REACTION BUILDER TAB
// ============================================================================

function addReactant() {
    const container = document.getElementById('reactantsContainer');
    const count = container.children.length;

    const div = document.createElement('div');
    div.className = 'row g-2 mb-2 reactant-item';
    div.innerHTML = `
        <div class="col-md-6">
            <input type="text" class="form-control form-control-sm reactant-species"
                   placeholder="Species name (e.g., Glucose)">
        </div>
        <div class="col-md-3">
            <input type="number" class="form-control form-control-sm reactant-coeff"
                   value="1" min="1" max="10">
        </div>
        <div class="col-md-3">
            <button type="button" class="btn btn-sm btn-outline-danger w-100"
                    onclick="this.parentElement.parentElement.remove()">
                <i class="fas fa-times"></i>
            </button>
        </div>
    `;
    container.appendChild(div);
}

function addProduct() {
    const container = document.getElementById('productsContainer');
    const count = container.children.length;

    const div = document.createElement('div');
    div.className = 'row g-2 mb-2 product-item';
    div.innerHTML = `
        <div class="col-md-6">
            <input type="text" class="form-control form-control-sm product-species"
                   placeholder="Species name (e.g., Pyruvate)">
        </div>
        <div class="col-md-3">
            <input type="number" class="form-control form-control-sm product-coeff"
                   value="1" min="1" max="10">
        </div>
        <div class="col-md-3">
            <button type="button" class="btn btn-sm btn-outline-danger w-100"
                    onclick="this.parentElement.parentElement.remove()">
                <i class="fas fa-times"></i>
            </button>
        </div>
    `;
    container.appendChild(div);
}

function addCatalyst() {
    const container = document.getElementById('catalystsContainer');

    const div = document.createElement('div');
    div.className = 'row g-2 mb-2 catalyst-item';
    div.innerHTML = `
        <div class="col-md-9">
            <input type="text" class="form-control form-control-sm catalyst-name"
                   placeholder="Enzyme/Catalyst name">
        </div>
        <div class="col-md-3">
            <button type="button" class="btn btn-sm btn-outline-danger w-100"
                    onclick="this.parentElement.parentElement.remove()">
                <i class="fas fa-times"></i>
            </button>
        </div>
    `;
    container.appendChild(div);
}

document.addEventListener('DOMContentLoaded', function() {
    const reactionForm = document.getElementById('reactionBuilderForm');
    if (reactionForm) {
        reactionForm.addEventListener('submit', function(e) {
            e.preventDefault();

            const reactionName = document.getElementById('reactionName').value.trim();
            const reversible = document.getElementById('reversibleCheck').checked;

            // Collect reactants
            const reactants = {};
            document.querySelectorAll('.reactant-item').forEach(item => {
                const species = item.querySelector('.reactant-species').value.trim();
                const coeff = parseInt(item.querySelector('.reactant-coeff').value);
                if (species) {
                    reactants[species] = -coeff; // Negative for reactants in Bio.Pathway
                }
            });

            // Collect products
            const products = {};
            document.querySelectorAll('.product-item').forEach(item => {
                const species = item.querySelector('.product-species').value.trim();
                const coeff = parseInt(item.querySelector('.product-coeff').value);
                if (species) {
                    products[species] = coeff; // Positive for products
                }
            });

            // Collect catalysts
            const catalysts = [];
            document.querySelectorAll('.catalyst-item').forEach(item => {
                const name = item.querySelector('.catalyst-name').value.trim();
                if (name) {
                    catalysts.push(name);
                }
            });

            if (Object.keys(reactants).length === 0 || Object.keys(products).length === 0) {
                showAlert('Please add at least one reactant and one product', 'warning');
                return;
            }

            // Combine reactants and products for Bio.Pathway format
            const allSpecies = {...reactants, ...products};

            const reaction = {
                name: reactionName,
                reactants: reactants,
                products: products,
                species: allSpecies,
                catalysts: catalysts,
                reversible: reversible
            };

            pathwayReactions.push(reaction);
            updateReactionsList();

            // Clear form
            reactionForm.reset();
            document.getElementById('reactantsContainer').innerHTML = '';
            document.getElementById('productsContainer').innerHTML = '';
            document.getElementById('catalystsContainer').innerHTML = '';

            // Add initial empty fields
            addReactant();
            addProduct();
        });
    }
});

function updateReactionsList() {
    const container = document.getElementById('reactionsList');
    const badge = document.getElementById('reactionsCount');

    badge.textContent = `${pathwayReactions.length} reaction${pathwayReactions.length !== 1 ? 's' : ''}`;

    if (pathwayReactions.length === 0) {
        container.innerHTML = '<p class="text-muted text-center small">No reactions added yet</p>';
        return;
    }

    let html = '';
    pathwayReactions.forEach((rxn, idx) => {
        const reactantsStr = Object.entries(rxn.reactants).map(([s, c]) => `${Math.abs(c)} ${s}`).join(' + ');
        const productsStr = Object.entries(rxn.products).map(([s, c]) => `${c} ${s}`).join(' + ');

        html += `
            <div class="card mb-2">
                <div class="card-body p-2">
                    <div class="d-flex justify-content-between align-items-start mb-1">
                        <h6 class="mb-0 small"><i class="fas fa-flask"></i> ${rxn.name}</h6>
                        <button class="btn btn-sm btn-outline-danger" onclick="removeReaction(${idx})">
                            <i class="fas fa-trash"></i>
                        </button>
                    </div>
                    <div class="small">
                        <span class="text-primary">${reactantsStr}</span>
                        <i class="fas ${rxn.reversible ? 'fa-exchange-alt text-success' : 'fa-arrow-right text-warning'} mx-1"></i>
                        <span class="text-info">${productsStr}</span>
                    </div>
                    ${rxn.catalysts.length > 0 ? `
                        <div class="mt-1">
                            <small class="badge bg-success">${rxn.catalysts.join(', ')}</small>
                        </div>
                    ` : ''}
                </div>
            </div>
        `;
    });

    container.innerHTML = html;
}

function removeReaction(idx) {
    pathwayReactions.splice(idx, 1);
    updateReactionsList();
}

function clearReactions() {
    pathwayReactions = [];
    updateReactionsList();
}

function loadReactionExample() {
    clearReactions();

    pathwayReactions = [
        {
            name: 'Glucose Phosphorylation',
            reactants: {'Glucose': -1, 'ATP': -1},
            products: {'Glucose-6-phosphate': 1, 'ADP': 1},
            species: {'Glucose': -1, 'ATP': -1, 'Glucose-6-phosphate': 1, 'ADP': 1},
            catalysts: ['Hexokinase'],
            reversible: false
        },
        {
            name: 'G6P Isomerization',
            reactants: {'Glucose-6-phosphate': -1},
            products: {'Fructose-6-phosphate': 1},
            species: {'Glucose-6-phosphate': -1, 'Fructose-6-phosphate': 1},
            catalysts: ['Phosphoglucose isomerase'],
            reversible: true
        },
        {
            name: 'F6P Phosphorylation',
            reactants: {'Fructose-6-phosphate': -1, 'ATP': -1},
            products: {'Fructose-1,6-bisphosphate': 1, 'ADP': 1},
            species: {'Fructose-6-phosphate': -1, 'ATP': -1, 'Fructose-1,6-bisphosphate': 1, 'ADP': 1},
            catalysts: ['Phosphofructokinase'],
            reversible: false
        }
    ];

    updateReactionsList();
}

// ============================================================================
// SYSTEM ANALYSIS TAB
// ============================================================================

function analyzeSystem() {
    if (pathwayReactions.length === 0) {
        showAlert('Please add reactions first', 'warning');
        return;
    }

    showLoading('analyzeSystemBtn');

    fetch('/api/pathway/analyze_system', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({reactions: pathwayReactions})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('analyzeSystemBtn', '<i class="fas fa-calculator"></i> Analyze System');

        if (data.success) {
            displaySystemAnalysis(data.analysis);
            currentSystem = data.analysis;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('analyzeSystemBtn', '<i class="fas fa-calculator"></i> Analyze System');
        showAlert('Request failed', 'danger');
    });
}

function displaySystemAnalysis(analysis) {
    const container = document.getElementById('systemAnalysisResults');

    let html = `
        <div class="row g-3 mb-3">
            <div class="col-md-3">
                <div class="card bg-primary text-white">
                    <div class="card-body p-3 text-center">
                        <h3 class="mb-0">${analysis.reaction_count}</h3>
                        <small>Reactions</small>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-info text-white">
                    <div class="card-body p-3 text-center">
                        <h3 class="mb-0">${analysis.species_count}</h3>
                        <small>Species</small>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-success text-white">
                    <div class="card-body p-3 text-center">
                        <h3 class="mb-0">${analysis.reversible_count}</h3>
                        <small>Reversible</small>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card bg-warning text-dark">
                    <div class="card-body p-3 text-center">
                        <h3 class="mb-0">${analysis.irreversible_count}</h3>
                        <small>Irreversible</small>
                    </div>
                </div>
            </div>
        </div>

        <div class="row">
            <div class="col-md-6">
                <h6 class="mb-2"><i class="fas fa-list"></i> All Species</h6>
                <div class="p-3 bg-light rounded">
                    ${analysis.species.map(s => `<span class="badge bg-secondary me-1 mb-1">${s}</span>`).join('')}
                </div>
            </div>
            <div class="col-md-6">
                <h6 class="mb-2"><i class="fas fa-info-circle"></i> System Properties</h6>
                <ul class="list-unstyled small">
                    <li><strong>Reversible reactions:</strong> ${analysis.reversible_reactions.join(', ') || 'None'}</li>
                    <li><strong>Irreversible reactions:</strong> ${analysis.irreversible_reactions.join(', ') || 'None'}</li>
                </ul>
            </div>
        </div>
    `;

    container.innerHTML = html;
}

// ============================================================================
// NETWORK ANALYSIS TAB
// ============================================================================

function analyzeNetwork() {
    if (pathwayReactions.length === 0) {
        showAlert('Please add reactions first', 'warning');
        return;
    }

    showLoading('analyzeNetworkBtn');

    fetch('/api/pathway/analyze_network', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({reactions: pathwayReactions})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('analyzeNetworkBtn', '<i class="fas fa-project-diagram"></i> Analyze Network');

        if (data.success) {
            displayNetworkAnalysis(data.analysis);
            currentNetwork = data.analysis;
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('analyzeNetworkBtn', '<i class="fas fa-project-diagram"></i> Analyze Network');
        showAlert('Request failed', 'danger');
    });
}

function displayNetworkAnalysis(analysis) {
    const container = document.getElementById('networkAnalysisResults');

    let html = `
        <div class="row g-3 mb-3">
            <div class="col-md-4">
                <div class="card border-success">
                    <div class="card-header bg-success text-white py-2">
                        <h6 class="mb-0 small"><i class="fas fa-sign-in-alt"></i> Sources</h6>
                    </div>
                    <div class="card-body p-2">
                        ${analysis.sources.length > 0 ?
                            analysis.sources.map(s => `<span class="badge bg-success me-1 mb-1">${s}</span>`).join('') :
                            '<p class="text-muted small mb-0">No sources found</p>'}
                    </div>
                </div>
            </div>
            <div class="col-md-4">
                <div class="card border-danger">
                    <div class="card-header bg-danger text-white py-2">
                        <h6 class="mb-0 small"><i class="fas fa-sign-out-alt"></i> Sinks</h6>
                    </div>
                    <div class="card-body p-2">
                        ${analysis.sinks.length > 0 ?
                            analysis.sinks.map(s => `<span class="badge bg-danger me-1 mb-1">${s}</span>`).join('') :
                            '<p class="text-muted small mb-0">No sinks found</p>'}
                    </div>
                </div>
            </div>
            <div class="col-md-4">
                <div class="card border-info">
                    <div class="card-header bg-info text-white py-2">
                        <h6 class="mb-0 small"><i class="fas fa-exchange-alt"></i> Intermediates</h6>
                    </div>
                    <div class="card-body p-2">
                        ${analysis.intermediates.length > 0 ?
                            analysis.intermediates.map(s => `<span class="badge bg-info me-1 mb-1">${s}</span>`).join('') :
                            '<p class="text-muted small mb-0">No intermediates</p>'}
                    </div>
                </div>
            </div>
        </div>

        <div class="card mb-3">
            <div class="card-header py-2">
                <h6 class="mb-0"><i class="fas fa-network-wired"></i> Interactions</h6>
            </div>
            <div class="card-body p-2">
                <div class="table-responsive">
                    <table class="table table-sm table-hover">
                        <thead>
                            <tr>
                                <th>Source</th>
                                <th>→</th>
                                <th>Sink</th>
                                <th>Reaction</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${analysis.interactions.map(i => `
                                <tr>
                                    <td><span class="badge bg-primary">${i.source}</span></td>
                                    <td><i class="fas fa-arrow-right text-muted"></i></td>
                                    <td><span class="badge bg-info">${i.sink}</span></td>
                                    <td><small>${i.reaction}</small></td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        ${analysis.species_connections ? `
            <div class="card">
                <div class="card-header py-2">
                    <h6 class="mb-0"><i class="fas fa-sitemap"></i> Species Connections</h6>
                </div>
                <div class="card-body p-2">
                    <div class="row">
                        ${Object.entries(analysis.species_connections).map(([species, conn]) => `
                            <div class="col-md-6 mb-2">
                                <strong>${species}</strong>
                                <ul class="small mb-0">
                                    ${conn.upstream.length > 0 ? `<li>Upstream: ${conn.upstream.join(', ')}</li>` : ''}
                                    ${conn.downstream.length > 0 ? `<li>Downstream: ${conn.downstream.join(', ')}</li>` : ''}
                                </ul>
                            </div>
                        `).join('')}
                    </div>
                </div>
            </div>
        ` : ''}
    `;

    container.innerHTML = html;
}

// ============================================================================
// VISUALIZATION TAB
// ============================================================================

function generateVisualization() {
    if (pathwayReactions.length === 0) {
        showAlert('Please add reactions first', 'warning');
        return;
    }

    showLoading('generateVizBtn');

    fetch('/api/pathway/visualize', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({reactions: pathwayReactions})
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('generateVizBtn', '<i class="fas fa-eye"></i> Generate Visualization');

        if (data.success) {
            displayVisualization(data.visualization);
        } else {
            showAlert('Error: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('generateVizBtn', '<i class="fas fa-eye"></i> Generate Visualization');
        showAlert('Request failed', 'danger');
    });
}

function displayVisualization(viz) {
    const container = document.getElementById('pathwayVisualization');

    if (viz.graph_image) {
        container.innerHTML = `
            <div class="text-center">
                <img src="${viz.graph_image}" class="img-fluid" alt="Pathway Graph"
                     style="max-width: 100%; border: 1px solid #ddd; border-radius: 8px;">
            </div>
            <div class="mt-3 p-3 bg-light rounded">
                <h6><i class="fas fa-info-circle"></i> Graph Information</h6>
                <ul class="small mb-0">
                    <li><strong>Nodes:</strong> ${viz.node_count}</li>
                    <li><strong>Edges:</strong> ${viz.edge_count}</li>
                    <li><strong>Graph Type:</strong> ${viz.graph_type}</li>
                </ul>
            </div>
        `;
    } else {
        container.innerHTML = '<div class="alert alert-warning">Visualization not available</div>';
    }
}

function exportPathway(format) {
    if (pathwayReactions.length === 0) {
        showAlert('Please add reactions first', 'warning');
        return;
    }

    fetch('/api/pathway/export', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            reactions: pathwayReactions,
            format: format
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            // Trigger download
            const blob = new Blob([data.content], {type: 'text/plain'});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `pathway.${format}`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);

            showAlert(`Pathway exported as ${format.toUpperCase()}`, 'success');
        } else {
            showAlert('Export failed: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        showAlert('Export request failed', 'danger');
    });
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

function showLoading(btnId) {
    const btn = document.getElementById(btnId);
    if (btn) {
        btn.disabled = true;
        btn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading...';
    }
}

function hideLoading(btnId, originalHtml) {
    const btn = document.getElementById(btnId);
    if (btn) {
        btn.disabled = false;
        btn.innerHTML = originalHtml;
    }
}

function showAlert(message, type) {
    let alertDiv = document.getElementById('pathwayAlert');
    if (!alertDiv) {
        alertDiv = document.createElement('div');
        alertDiv.id = 'pathwayAlert';
        alertDiv.style.position = 'fixed';
        alertDiv.style.top = '20px';
        alertDiv.style.right = '20px';
        alertDiv.style.zIndex = '9999';
        document.body.appendChild(alertDiv);
    }

    alertDiv.innerHTML = `
        <div class="alert alert-${type} alert-dismissible fade show" role="alert">
            ${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        </div>
    `;

    setTimeout(() => {
        alertDiv.innerHTML = '';
    }, 5000);
}

// Initialize
window.addEventListener('load', function() {
    // Add initial fields
    if (document.getElementById('reactantsContainer')) {
        addReactant();
        addProduct();
    }
});
