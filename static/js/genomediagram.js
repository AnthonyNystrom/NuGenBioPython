/**
 * Genome Diagrams - Complete JavaScript
 * Handles all tabs and functionality
 */

// Global state
let featureCount = 0;
let trackCount = 0;
let graphCount = 0;
let linkCount = 0;
let currentDiagramData = null;

// ============================================================================
// FILE UPLOAD HANDLERS
// ============================================================================

function handleBasicFileUpload(event) {
    handleFileUpload(event, 'basic');
}

function handleMultiTrackFileUpload(event) {
    handleFileUpload(event, 'multitrack');
}

function handleAdvancedFileUpload(event) {
    handleFileUpload(event, 'advanced');
}

function handleFileUpload(event, context) {
    const file = event.target.files[0];
    if (!file) return;

    const formData = new FormData();
    formData.append('file', file);

    // Show loading
    const infoId = `${context}FileInfo`;
    const infoDiv = document.getElementById(infoId);
    if (infoDiv) {
        infoDiv.style.display = 'block';
        infoDiv.innerHTML = '<small class="text-muted"><i class="fas fa-spinner fa-spin"></i> Parsing file...</small>';
    }

    fetch('/api/genomediagram/upload_file', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if (data.success) {
            // Populate with parsed features
            if (context === 'basic') {
                loadFeaturesFromFile(data);
            } else if (context === 'multitrack') {
                loadMultiTrackFeaturesFromFile(data);
            } else if (context === 'advanced') {
                loadAdvancedFeaturesFromFile(data);
            }

            // Show file info
            if (infoDiv) {
                let infoHtml = `<div class="alert alert-success py-1 px-2 small mb-0">`;
                infoHtml += `<i class="fas fa-check-circle"></i> <strong>${data.filename}</strong><br>`;
                infoHtml += `Format: ${data.info.format} | Length: ${data.info.length.toLocaleString()} bp`;
                if (data.info.features) {
                    infoHtml += ` | ${data.info.features} features`;
                }
                infoHtml += `</div>`;
                infoDiv.innerHTML = infoHtml;
            }
        } else {
            if (infoDiv) {
                infoDiv.innerHTML = `<div class="alert alert-danger py-1 px-2 small mb-0">
                    <i class="fas fa-exclamation-circle"></i> ${data.error}
                </div>`;
            }
        }
    })
    .catch(error => {
        if (infoDiv) {
            infoDiv.innerHTML = `<div class="alert alert-danger py-1 px-2 small mb-0">
                <i class="fas fa-exclamation-circle"></i> Upload failed
            </div>`;
        }
    });
}

function loadFeaturesFromFile(data) {
    const container = document.getElementById('featuresContainer');
    container.innerHTML = '';
    featureCount = 0;

    document.getElementById('genomeLength').value = data.info.length;

    data.features.forEach(feature => {
        addFeature();
        const lastFeature = document.querySelector('.feature-item:last-child');
        lastFeature.querySelector('.feature-name').value = feature.name;
        lastFeature.querySelector('.feature-start').value = feature.start;
        lastFeature.querySelector('.feature-end').value = feature.end;
        lastFeature.querySelector('.feature-color').value = feature.color || 'blue';
    });

    updateFeatureSummary();
}

function loadMultiTrackFeaturesFromFile(data) {
    const container = document.getElementById('tracksContainer');
    container.innerHTML = '';
    trackCount = 0;

    document.getElementById('multiTrackGenomeLength').value = data.info.length;

    // Group features by type for different tracks
    const featuresByType = {};
    data.features.forEach(feature => {
        const type = feature.type || 'gene';
        if (!featuresByType[type]) {
            featuresByType[type] = [];
        }
        featuresByType[type].push(feature);
    });

    // Create tracks
    Object.keys(featuresByType).forEach(type => {
        addTrack();
        const lastTrack = document.querySelector('.track-item:last-child');
        lastTrack.querySelector('.track-name').value = type.charAt(0).toUpperCase() + type.slice(1);

        // Add features to track
        const featuresContainer = lastTrack.querySelector('.track-features');
        featuresByType[type].forEach(feature => {
            const featureDiv = document.createElement('div');
            featureDiv.className = 'row g-2 mb-2 track-feature-item p-2 border rounded bg-light';
            featureDiv.innerHTML = `
                <div class="col-3">
                    <input type="text" class="form-control form-control-sm track-feature-name"
                           placeholder="Name" value="${feature.name}">
                </div>
                <div class="col-2">
                    <input type="number" class="form-control form-control-sm track-feature-start"
                           min="1" value="${feature.start}">
                </div>
                <div class="col-2">
                    <input type="number" class="form-control form-control-sm track-feature-end"
                           min="1" value="${feature.end}">
                </div>
                <div class="col-2">
                    <select class="form-select form-select-sm track-feature-color">
                        <option value="${feature.color || 'blue'}" selected>${feature.color || 'blue'}</option>
                        <option value="blue">Blue</option>
                        <option value="red">Red</option>
                        <option value="green">Green</option>
                        <option value="orange">Orange</option>
                        <option value="purple">Purple</option>
                    </select>
                </div>
                <div class="col-3">
                    <button type="button" class="btn btn-sm btn-outline-danger w-100"
                            onclick="this.parentElement.parentElement.remove()">
                        <i class="fas fa-times"></i>
                    </button>
                </div>
            `;
            featuresContainer.appendChild(featureDiv);
        });
    });
}

function loadAdvancedFeaturesFromFile(data) {
    const container = document.getElementById('advancedFeaturesContainer');
    container.innerHTML = '';
    featureCount = 0;

    document.getElementById('advancedGenomeLength').value = data.info.length;

    data.features.forEach(feature => {
        addAdvancedFeature();
        const lastFeature = document.querySelector('.advanced-feature-item:last-child');
        lastFeature.querySelector('.advanced-feature-name').value = feature.name;
        lastFeature.querySelector('.advanced-feature-start').value = feature.start;
        lastFeature.querySelector('.advanced-feature-end').value = feature.end;
        lastFeature.querySelector('.advanced-feature-color').value = feature.color || 'blue';
        lastFeature.querySelector('.advanced-feature-strand').value = feature.strand || '1';
    });
}

// ============================================================================
// BASIC DIAGRAM TAB
// ============================================================================

function addFeature() {
    featureCount++;
    const container = document.getElementById('featuresContainer');

    const featureDiv = document.createElement('div');
    featureDiv.className = 'row g-2 mb-2 p-2 border rounded bg-light feature-item';
    featureDiv.id = `feature-${featureCount}`;

    featureDiv.innerHTML = `
        <div class="col-md-3">
            <input type="text" class="form-control form-control-sm feature-name"
                   placeholder="Feature ${featureCount}" value="Feature ${featureCount}">
        </div>
        <div class="col-md-2">
            <input type="number" class="form-control form-control-sm feature-start"
                   min="1" value="${1000 * featureCount}">
        </div>
        <div class="col-md-2">
            <input type="number" class="form-control form-control-sm feature-end"
                   min="1" value="${1000 * featureCount + 500}">
        </div>
        <div class="col-md-3">
            <select class="form-select form-select-sm feature-color">
                <option value="blue">Blue</option>
                <option value="red">Red</option>
                <option value="green">Green</option>
                <option value="orange">Orange</option>
                <option value="purple">Purple</option>
                <option value="brown">Brown</option>
                <option value="pink">Pink</option>
                <option value="cyan">Cyan</option>
            </select>
        </div>
        <div class="col-md-2">
            <button type="button" class="btn btn-outline-danger btn-sm w-100" onclick="removeFeature(${featureCount})">
                <i class="fas fa-trash"></i>
            </button>
        </div>
    `;

    container.appendChild(featureDiv);
    updateFeatureSummary();
}

function removeFeature(id) {
    const featureDiv = document.getElementById(`feature-${id}`);
    if (featureDiv) {
        featureDiv.remove();
        updateFeatureSummary();
    }
}

function clearFeatures() {
    const container = document.getElementById('featuresContainer');
    container.innerHTML = '';
    featureCount = 0;
    updateFeatureSummary();
}

function updateFeatureSummary() {
    const summaryDiv = document.getElementById('featureSummary');
    const features = document.querySelectorAll('.feature-item');

    if (features.length === 0) {
        summaryDiv.innerHTML = `
            <p class="text-muted text-center small mb-0">
                <i class="fas fa-plus fa-2x mb-2"></i><br>
                Add features to see summary
            </p>
        `;
        return;
    }

    const colorCounts = {};
    features.forEach(feature => {
        const color = feature.querySelector('.feature-color').value;
        colorCounts[color] = (colorCounts[color] || 0) + 1;
    });

    let html = `
        <div class="text-center mb-2">
            <h5 class="mb-0">${features.length}</h5>
            <small class="text-muted">Total Features</small>
        </div>
        <hr class="my-2">
        <p class="small mb-1"><strong>Color Distribution:</strong></p>
    `;

    Object.entries(colorCounts).forEach(([color, count]) => {
        html += `
            <div class="d-flex justify-content-between align-items-center mb-1">
                <span class="small"><i class="fas fa-square" style="color: ${color}"></i> ${color}</span>
                <span class="badge bg-secondary">${count}</span>
            </div>
        `;
    });

    summaryDiv.innerHTML = html;
}

function loadBasicExample() {
    clearFeatures();

    document.getElementById('genomeLength').value = '25000';
    document.getElementById('diagramTitle').value = 'Example Bacterial Operon';

    const exampleFeatures = [
        { name: 'Promoter', start: 1000, end: 1200, color: 'green' },
        { name: 'Gene A', start: 1500, end: 3000, color: 'blue' },
        { name: 'Gene B', start: 3200, end: 4800, color: 'blue' },
        { name: 'Gene C', start: 5000, end: 6500, color: 'blue' },
        { name: 'Terminator', start: 6800, end: 7000, color: 'red' },
        { name: 'Repeat Region', start: 10000, end: 12000, color: 'orange' },
        { name: 'Regulatory Site', start: 15000, end: 15300, color: 'purple' },
        { name: 'Transposon', start: 18000, end: 20000, color: 'brown' }
    ];

    exampleFeatures.forEach(feature => {
        addFeature();
        const lastFeature = document.querySelector('.feature-item:last-child');
        lastFeature.querySelector('.feature-name').value = feature.name;
        lastFeature.querySelector('.feature-start').value = feature.start;
        lastFeature.querySelector('.feature-end').value = feature.end;
        lastFeature.querySelector('.feature-color').value = feature.color;
    });

    updateFeatureSummary();
}

// ============================================================================
// MULTI-TRACK TAB
// ============================================================================

function addTrack() {
    trackCount++;
    const container = document.getElementById('tracksContainer');

    const trackDiv = document.createElement('div');
    trackDiv.className = 'card mb-3 track-item';
    trackDiv.id = `track-${trackCount}`;

    trackDiv.innerHTML = `
        <div class="card-header bg-light py-2">
            <div class="row g-2 align-items-center">
                <div class="col-md-3">
                    <label class="form-label small mb-0">Track Name</label>
                    <input type="text" class="form-control form-control-sm track-name"
                           value="Track ${trackCount}">
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Height</label>
                    <input type="number" class="form-control form-control-sm track-height"
                           min="0.1" max="2" step="0.1" value="0.5">
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Greytrack</label>
                    <select class="form-select form-select-sm track-greytrack">
                        <option value="true">Yes</option>
                        <option value="false">No</option>
                    </select>
                </div>
                <div class="col-md-3">
                    <label class="form-label small mb-0">Track Number</label>
                    <input type="number" class="form-control form-control-sm track-number"
                           min="1" value="${trackCount}">
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">&nbsp;</label>
                    <button type="button" class="btn btn-outline-danger btn-sm w-100"
                            onclick="removeTrack(${trackCount})">
                        <i class="fas fa-trash"></i> Remove
                    </button>
                </div>
            </div>
        </div>
        <div class="card-body p-2">
            <div class="track-features" id="track-${trackCount}-features">
                <!-- Features for this track -->
            </div>
            <button type="button" class="btn btn-sm btn-outline-primary mt-2"
                    onclick="addTrackFeature(${trackCount})">
                <i class="fas fa-plus"></i> Add Feature to Track
            </button>
        </div>
    `;

    container.appendChild(trackDiv);
}

function removeTrack(id) {
    const trackDiv = document.getElementById(`track-${id}`);
    if (trackDiv) {
        trackDiv.remove();
    }
}

function addTrackFeature(trackId) {
    const container = document.getElementById(`track-${trackId}-features`);

    const featureDiv = document.createElement('div');
    featureDiv.className = 'row g-2 mb-2 track-feature-item p-2 border rounded bg-light';

    featureDiv.innerHTML = `
        <div class="col-md-3">
            <input type="text" class="form-control form-control-sm track-feature-name"
                   placeholder="Feature name">
        </div>
        <div class="col-md-2">
            <input type="number" class="form-control form-control-sm track-feature-start"
                   min="1" placeholder="Start">
        </div>
        <div class="col-md-2">
            <input type="number" class="form-control form-control-sm track-feature-end"
                   min="1" placeholder="End">
        </div>
        <div class="col-md-2">
            <select class="form-select form-select-sm track-feature-color">
                <option value="blue">Blue</option>
                <option value="red">Red</option>
                <option value="green">Green</option>
                <option value="orange">Orange</option>
                <option value="purple">Purple</option>
                <option value="brown">Brown</option>
            </select>
        </div>
        <div class="col-md-1">
            <select class="form-select form-select-sm track-feature-strand">
                <option value="1">+</option>
                <option value="-1">-</option>
                <option value="0">None</option>
            </select>
        </div>
        <div class="col-md-2">
            <button type="button" class="btn btn-outline-danger btn-sm w-100"
                    onclick="this.parentElement.parentElement.remove()">
                <i class="fas fa-times"></i>
            </button>
        </div>
    `;

    container.appendChild(featureDiv);
}

function clearTracks() {
    const container = document.getElementById('tracksContainer');
    container.innerHTML = '';
    trackCount = 0;
}

function loadMultiTrackExample() {
    clearTracks();

    document.getElementById('multiTrackGenomeLength').value = '30000';
    document.getElementById('multiTrackTitle').value = 'Multi-Track Genome View';

    // Track 1: Genes
    addTrack();
    let track = document.querySelector('.track-item:last-child');
    track.querySelector('.track-name').value = 'Genes';
    track.querySelector('.track-height').value = '0.5';

    const genesFeatures = [
        { name: 'geneA', start: 1000, end: 2500, color: 'blue', strand: 1 },
        { name: 'geneB', start: 3000, end: 4500, color: 'blue', strand: 1 },
        { name: 'geneC', start: 5000, end: 6000, color: 'blue', strand: -1 }
    ];

    genesFeatures.forEach(f => {
        addTrackFeature(trackCount);
        const lastFeature = track.querySelector('.track-feature-item:last-child');
        lastFeature.querySelector('.track-feature-name').value = f.name;
        lastFeature.querySelector('.track-feature-start').value = f.start;
        lastFeature.querySelector('.track-feature-end').value = f.end;
        lastFeature.querySelector('.track-feature-color').value = f.color;
        lastFeature.querySelector('.track-feature-strand').value = f.strand;
    });

    // Track 2: Regulatory
    addTrack();
    track = document.querySelector('.track-item:last-child');
    track.querySelector('.track-name').value = 'Regulatory Elements';
    track.querySelector('.track-height').value = '0.3';
    track.querySelector('.track-number').value = '2';

    const regFeatures = [
        { name: 'Promoter1', start: 800, end: 1000, color: 'green', strand: 1 },
        { name: 'Promoter2', start: 2800, end: 3000, color: 'green', strand: 1 },
        { name: 'Enhancer', start: 8000, end: 8500, color: 'orange', strand: 0 }
    ];

    regFeatures.forEach(f => {
        addTrackFeature(trackCount);
        const lastFeature = track.querySelector('.track-feature-item:last-child');
        lastFeature.querySelector('.track-feature-name').value = f.name;
        lastFeature.querySelector('.track-feature-start').value = f.start;
        lastFeature.querySelector('.track-feature-end').value = f.end;
        lastFeature.querySelector('.track-feature-color').value = f.color;
        lastFeature.querySelector('.track-feature-strand').value = f.strand;
    });
}

// ============================================================================
// DATA TRACKS TAB
// ============================================================================

function addDataGraph() {
    graphCount++;
    const container = document.getElementById('dataGraphsContainer');

    const graphDiv = document.createElement('div');
    graphDiv.className = 'card mb-3 data-graph-item';
    graphDiv.id = `graph-${graphCount}`;

    graphDiv.innerHTML = `
        <div class="card-header bg-light py-2">
            <div class="row g-2 align-items-center">
                <div class="col-md-3">
                    <label class="form-label small mb-0">Graph Type</label>
                    <select class="form-select form-select-sm graph-type">
                        <option value="gc_content">GC Content</option>
                        <option value="gc_skew">GC Skew</option>
                        <option value="custom">Custom Data</option>
                    </select>
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Style</label>
                    <select class="form-select form-select-sm graph-style">
                        <option value="line">Line</option>
                        <option value="bar">Bar</option>
                        <option value="heat">Heat</option>
                    </select>
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Window</label>
                    <input type="number" class="form-control form-control-sm graph-window"
                           value="1000" min="10">
                </div>
                <div class="col-md-3">
                    <label class="form-label small mb-0">Color</label>
                    <select class="form-select form-select-sm graph-color">
                        <option value="blue">Blue</option>
                        <option value="red">Red</option>
                        <option value="green">Green</option>
                        <option value="purple">Purple</option>
                        <option value="black">Black</option>
                    </select>
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">&nbsp;</label>
                    <button type="button" class="btn btn-outline-danger btn-sm w-100"
                            onclick="removeDataGraph(${graphCount})">
                        <i class="fas fa-trash"></i>
                    </button>
                </div>
            </div>
        </div>
        <div class="card-body p-2">
            <div class="custom-data-section" style="display: none;">
                <label class="form-label small">Custom Data (comma-separated values)</label>
                <textarea class="form-control form-control-sm graph-custom-data" rows="2"
                          placeholder="0.5,0.6,0.4,0.7,0.5..."></textarea>
            </div>
        </div>
    `;

    container.appendChild(graphDiv);

    // Toggle custom data section
    graphDiv.querySelector('.graph-type').addEventListener('change', function(e) {
        const customSection = graphDiv.querySelector('.custom-data-section');
        customSection.style.display = e.target.value === 'custom' ? 'block' : 'none';
    });
}

function removeDataGraph(id) {
    const graphDiv = document.getElementById(`graph-${id}`);
    if (graphDiv) {
        graphDiv.remove();
    }
}

function clearDataGraphs() {
    const container = document.getElementById('dataGraphsContainer');
    container.innerHTML = '';
    graphCount = 0;
}

function loadDataTrackExample() {
    clearDataGraphs();

    document.getElementById('dataGenomeLength').value = '20000';
    document.getElementById('dataTitle').value = 'Genome Data Visualization';
    document.getElementById('dataSequence').value = 'ATGCATGCATGCATGCGCGCGCGCATATATATATGCGCGCGC';

    // Add GC content graph
    addDataGraph();
    let graph = document.querySelector('.data-graph-item:last-child');
    graph.querySelector('.graph-type').value = 'gc_content';
    graph.querySelector('.graph-style').value = 'line';
    graph.querySelector('.graph-window').value = '500';

    // Add GC skew graph
    addDataGraph();
    graph = document.querySelector('.data-graph-item:last-child');
    graph.querySelector('.graph-type').value = 'gc_skew';
    graph.querySelector('.graph-style').value = 'bar';
    graph.querySelector('.graph-color').value = 'purple';
}

// ============================================================================
// ADVANCED FEATURES TAB
// ============================================================================

function addAdvancedFeature() {
    featureCount++;
    const container = document.getElementById('advancedFeaturesContainer');

    const featureDiv = document.createElement('div');
    featureDiv.className = 'card mb-2 advanced-feature-item';
    featureDiv.id = `advanced-feature-${featureCount}`;

    featureDiv.innerHTML = `
        <div class="card-body p-2">
            <div class="row g-2">
                <div class="col-md-2">
                    <label class="form-label small mb-0">Name</label>
                    <input type="text" class="form-control form-control-sm advanced-feature-name"
                           value="Feature ${featureCount}">
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">Start</label>
                    <input type="number" class="form-control form-control-sm advanced-feature-start"
                           min="1" value="${1000 * featureCount}">
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">End</label>
                    <input type="number" class="form-control form-control-sm advanced-feature-end"
                           min="1" value="${1000 * featureCount + 500}">
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">Strand</label>
                    <select class="form-select form-select-sm advanced-feature-strand">
                        <option value="1">+ (forward)</option>
                        <option value="-1">- (reverse)</option>
                        <option value="0">None</option>
                    </select>
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Sigil</label>
                    <select class="form-select form-select-sm advanced-feature-sigil">
                        <option value="BOX">BOX (rectangle)</option>
                        <option value="ARROW">ARROW (directional)</option>
                        <option value="OCTO">OCTO (octagon)</option>
                        <option value="JAGGY">JAGGY (jagged)</option>
                        <option value="BIGARROW">BIGARROW (large)</option>
                    </select>
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">Color</label>
                    <select class="form-select form-select-sm advanced-feature-color">
                        <option value="blue">Blue</option>
                        <option value="red">Red</option>
                        <option value="green">Green</option>
                        <option value="orange">Orange</option>
                        <option value="purple">Purple</option>
                    </select>
                </div>
                <div class="col-md-2">
                    <label class="form-label small mb-0">Label Position</label>
                    <select class="form-select form-select-sm advanced-feature-label-pos">
                        <option value="start">Start</option>
                        <option value="middle">Middle</option>
                        <option value="end">End</option>
                    </select>
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">Label</label>
                    <select class="form-select form-select-sm advanced-feature-label">
                        <option value="true">Yes</option>
                        <option value="false">No</option>
                    </select>
                </div>
                <div class="col-md-1">
                    <button type="button" class="btn btn-outline-danger btn-sm w-100 mt-4"
                            onclick="removeAdvancedFeature(${featureCount})">
                        <i class="fas fa-trash"></i>
                    </button>
                </div>
            </div>
        </div>
    `;

    container.appendChild(featureDiv);
}

function removeAdvancedFeature(id) {
    const featureDiv = document.getElementById(`advanced-feature-${id}`);
    if (featureDiv) {
        featureDiv.remove();
    }
}

function clearAdvancedFeatures() {
    const container = document.getElementById('advancedFeaturesContainer');
    container.innerHTML = '';
    featureCount = 0;
}

function addCrossLink() {
    linkCount++;
    const container = document.getElementById('crossLinksContainer');

    const linkDiv = document.createElement('div');
    linkDiv.className = 'card mb-2 cross-link-item';
    linkDiv.id = `link-${linkCount}`;

    linkDiv.innerHTML = `
        <div class="card-body p-2">
            <div class="row g-2">
                <div class="col-md-5">
                    <label class="form-label small mb-0">Track 1 (Start - End)</label>
                    <div class="input-group input-group-sm">
                        <input type="number" class="form-control form-control-sm link-track1"
                               value="1" min="1">
                        <input type="number" class="form-control form-control-sm link-start1"
                               placeholder="Start" min="1">
                        <input type="number" class="form-control form-control-sm link-end1"
                               placeholder="End" min="1">
                    </div>
                </div>
                <div class="col-md-5">
                    <label class="form-label small mb-0">Track 2 (Start - End)</label>
                    <div class="input-group input-group-sm">
                        <input type="number" class="form-control form-control-sm link-track2"
                               value="2" min="1">
                        <input type="number" class="form-control form-control-sm link-start2"
                               placeholder="Start" min="1">
                        <input type="number" class="form-control form-control-sm link-end2"
                               placeholder="End" min="1">
                    </div>
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">Color</label>
                    <select class="form-select form-select-sm link-color">
                        <option value="lightblue">Light Blue</option>
                        <option value="lightgreen">Light Green</option>
                        <option value="lightcoral">Light Coral</option>
                        <option value="lightgray">Light Gray</option>
                    </select>
                </div>
                <div class="col-md-1">
                    <label class="form-label small mb-0">&nbsp;</label>
                    <button type="button" class="btn btn-outline-danger btn-sm w-100"
                            onclick="removeCrossLink(${linkCount})">
                        <i class="fas fa-trash"></i>
                    </button>
                </div>
            </div>
        </div>
    `;

    container.appendChild(linkDiv);
}

function removeCrossLink(id) {
    const linkDiv = document.getElementById(`link-${id}`);
    if (linkDiv) {
        linkDiv.remove();
    }
}

function clearCrossLinks() {
    const container = document.getElementById('crossLinksContainer');
    container.innerHTML = '';
    linkCount = 0;
}

function loadAdvancedExample() {
    clearAdvancedFeatures();
    clearCrossLinks();

    document.getElementById('advancedGenomeLength').value = '15000';
    document.getElementById('advancedTitle').value = 'Advanced Feature Styling';

    // Add features with different sigils
    const features = [
        { name: 'Gene1', start: 1000, end: 2500, strand: 1, sigil: 'ARROW', color: 'blue', label: true, pos: 'middle' },
        { name: 'Gene2', start: 3000, end: 4000, strand: -1, sigil: 'ARROW', color: 'red', label: true, pos: 'middle' },
        { name: 'Promoter', start: 800, end: 1000, strand: 1, sigil: 'BOX', color: 'green', label: true, pos: 'start' },
        { name: 'Repeat', start: 5000, end: 7000, strand: 0, sigil: 'JAGGY', color: 'orange', label: true, pos: 'middle' }
    ];

    features.forEach(f => {
        addAdvancedFeature();
        const lastFeature = document.querySelector('.advanced-feature-item:last-child');
        lastFeature.querySelector('.advanced-feature-name').value = f.name;
        lastFeature.querySelector('.advanced-feature-start').value = f.start;
        lastFeature.querySelector('.advanced-feature-end').value = f.end;
        lastFeature.querySelector('.advanced-feature-strand').value = f.strand;
        lastFeature.querySelector('.advanced-feature-sigil').value = f.sigil;
        lastFeature.querySelector('.advanced-feature-color').value = f.color;
        lastFeature.querySelector('.advanced-feature-label').value = f.label.toString();
        lastFeature.querySelector('.advanced-feature-label-pos').value = f.pos;
    });
}

// ============================================================================
// FORM SUBMISSIONS
// ============================================================================

document.addEventListener('DOMContentLoaded', function() {

    // Basic Diagram Form
    const basicForm = document.getElementById('basicDiagramForm');
    if (basicForm) {
        basicForm.addEventListener('submit', function(e) {
            e.preventDefault();

            const genomeLength = parseInt(document.getElementById('genomeLength').value);
            const features = [];

            document.querySelectorAll('.feature-item').forEach(featureDiv => {
                const name = featureDiv.querySelector('.feature-name').value;
                const start = parseInt(featureDiv.querySelector('.feature-start').value);
                const end = parseInt(featureDiv.querySelector('.feature-end').value);
                const color = featureDiv.querySelector('.feature-color').value;

                if (name && start && end && start < end && start >= 1 && end <= genomeLength) {
                    features.push({ name, start, end, color });
                }
            });

            if (features.length === 0) {
                showAlert('Please add at least one valid feature', 'warning');
                return;
            }

            showLoading('createBasicDiagramBtn');

            const diagramType = document.getElementById('diagramType').value;
            const pageSize = document.getElementById('pageSize').value;
            const diagramTitle = document.getElementById('diagramTitle').value;

            fetch('/api/genomediagram/create', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    genome_length: genomeLength,
                    features: features,
                    diagram_type: diagramType,
                    page_size: pageSize,
                    title: diagramTitle
                })
            })
            .then(response => response.json())
            .then(data => {
                hideLoading('createBasicDiagramBtn', '<i class="fas fa-chart-area"></i> Create Diagram');

                if (data.success) {
                    displayDiagram(data.diagram, 'basicDiagramDisplay');
                    currentDiagramData = data.diagram;
                } else {
                    showAlert('Error: ' + data.error, 'danger');
                }
            })
            .catch(error => {
                hideLoading('createBasicDiagramBtn', '<i class="fas fa-chart-area"></i> Create Diagram');
                showAlert('Request failed', 'danger');
            });
        });
    }

    // Multi-Track Form
    const multiTrackForm = document.getElementById('multiTrackForm');
    if (multiTrackForm) {
        multiTrackForm.addEventListener('submit', function(e) {
            e.preventDefault();

            const genomeLength = parseInt(document.getElementById('multiTrackGenomeLength').value);
            const tracks = [];

            document.querySelectorAll('.track-item').forEach(trackDiv => {
                const trackName = trackDiv.querySelector('.track-name').value;
                const trackHeight = parseFloat(trackDiv.querySelector('.track-height').value);
                const trackGreytrack = trackDiv.querySelector('.track-greytrack').value === 'true';
                const trackNumber = parseInt(trackDiv.querySelector('.track-number').value);

                const features = [];
                trackDiv.querySelectorAll('.track-feature-item').forEach(featureDiv => {
                    const name = featureDiv.querySelector('.track-feature-name').value;
                    const start = parseInt(featureDiv.querySelector('.track-feature-start').value);
                    const end = parseInt(featureDiv.querySelector('.track-feature-end').value);
                    const color = featureDiv.querySelector('.track-feature-color').value;
                    const strand = parseInt(featureDiv.querySelector('.track-feature-strand').value);

                    if (name && start && end && start < end) {
                        features.push({ name, start, end, color, strand });
                    }
                });

                if (features.length > 0) {
                    tracks.push({
                        name: trackName,
                        height: trackHeight,
                        greytrack: trackGreytrack,
                        track_number: trackNumber,
                        features: features
                    });
                }
            });

            if (tracks.length === 0) {
                showAlert('Please add at least one track with features', 'warning');
                return;
            }

            showLoading('createMultiTrackBtn');

            const diagramType = document.getElementById('multiTrackDiagramType').value;
            const pageSize = document.getElementById('multiTrackPageSize').value;
            const title = document.getElementById('multiTrackTitle').value;

            fetch('/api/genomediagram/create_multitrack', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    genome_length: genomeLength,
                    tracks: tracks,
                    diagram_type: diagramType,
                    page_size: pageSize,
                    title: title
                })
            })
            .then(response => response.json())
            .then(data => {
                hideLoading('createMultiTrackBtn', '<i class="fas fa-layer-group"></i> Create Multi-Track Diagram');

                if (data.success) {
                    displayDiagram(data.diagram, 'multiTrackDisplay');
                    currentDiagramData = data.diagram;
                } else {
                    showAlert('Error: ' + data.error, 'danger');
                }
            })
            .catch(error => {
                hideLoading('createMultiTrackBtn', '<i class="fas fa-layer-group"></i> Create Multi-Track Diagram');
                showAlert('Request failed', 'danger');
            });
        });
    }

    // Data Tracks Form
    const dataForm = document.getElementById('dataTracksForm');
    if (dataForm) {
        dataForm.addEventListener('submit', function(e) {
            e.preventDefault();

            const genomeLength = parseInt(document.getElementById('dataGenomeLength').value);
            const sequence = document.getElementById('dataSequence').value.toUpperCase().replace(/\s/g, '');
            const graphs = [];

            document.querySelectorAll('.data-graph-item').forEach(graphDiv => {
                const graphType = graphDiv.querySelector('.graph-type').value;
                const graphStyle = graphDiv.querySelector('.graph-style').value;
                const graphWindow = parseInt(graphDiv.querySelector('.graph-window').value);
                const graphColor = graphDiv.querySelector('.graph-color').value;

                let customData = null;
                if (graphType === 'custom') {
                    const dataText = graphDiv.querySelector('.graph-custom-data').value;
                    if (dataText) {
                        customData = dataText.split(',').map(v => parseFloat(v.trim())).filter(v => !isNaN(v));
                    }
                }

                graphs.push({
                    type: graphType,
                    style: graphStyle,
                    window: graphWindow,
                    color: graphColor,
                    custom_data: customData
                });
            });

            if (graphs.length === 0) {
                showAlert('Please add at least one data graph', 'warning');
                return;
            }

            if (!sequence && graphs.some(g => g.type !== 'custom')) {
                showAlert('Please provide a DNA sequence for GC content/skew calculation', 'warning');
                return;
            }

            showLoading('createDataTracksBtn');

            const title = document.getElementById('dataTitle').value;
            const diagramType = document.getElementById('dataDiagramType').value;

            fetch('/api/genomediagram/create_data_tracks', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    genome_length: genomeLength,
                    sequence: sequence,
                    graphs: graphs,
                    title: title,
                    diagram_type: diagramType
                })
            })
            .then(response => response.json())
            .then(data => {
                hideLoading('createDataTracksBtn', '<i class="fas fa-chart-line"></i> Create Data Tracks');

                if (data.success) {
                    displayDiagram(data.diagram, 'dataTracksDisplay');
                    currentDiagramData = data.diagram;
                } else {
                    showAlert('Error: ' + data.error, 'danger');
                }
            })
            .catch(error => {
                hideLoading('createDataTracksBtn', '<i class="fas fa-chart-line"></i> Create Data Tracks');
                showAlert('Request failed', 'danger');
            });
        });
    }

    // Advanced Features Form
    const advancedForm = document.getElementById('advancedFeaturesForm');
    if (advancedForm) {
        advancedForm.addEventListener('submit', function(e) {
            e.preventDefault();

            const genomeLength = parseInt(document.getElementById('advancedGenomeLength').value);
            const features = [];

            document.querySelectorAll('.advanced-feature-item').forEach(featureDiv => {
                const name = featureDiv.querySelector('.advanced-feature-name').value;
                const start = parseInt(featureDiv.querySelector('.advanced-feature-start').value);
                const end = parseInt(featureDiv.querySelector('.advanced-feature-end').value);
                const strand = parseInt(featureDiv.querySelector('.advanced-feature-strand').value);
                const sigil = featureDiv.querySelector('.advanced-feature-sigil').value;
                const color = featureDiv.querySelector('.advanced-feature-color').value;
                const label = featureDiv.querySelector('.advanced-feature-label').value === 'true';
                const labelPos = featureDiv.querySelector('.advanced-feature-label-pos').value;

                if (name && start && end && start < end) {
                    features.push({
                        name, start, end, strand, sigil, color, label, label_position: labelPos
                    });
                }
            });

            const crossLinks = [];
            document.querySelectorAll('.cross-link-item').forEach(linkDiv => {
                const track1 = parseInt(linkDiv.querySelector('.link-track1').value);
                const start1 = parseInt(linkDiv.querySelector('.link-start1').value);
                const end1 = parseInt(linkDiv.querySelector('.link-end1').value);
                const track2 = parseInt(linkDiv.querySelector('.link-track2').value);
                const start2 = parseInt(linkDiv.querySelector('.link-start2').value);
                const end2 = parseInt(linkDiv.querySelector('.link-end2').value);
                const color = linkDiv.querySelector('.link-color').value;

                if (start1 && end1 && start2 && end2) {
                    crossLinks.push({
                        track1, start1, end1, track2, start2, end2, color
                    });
                }
            });

            if (features.length === 0) {
                showAlert('Please add at least one feature', 'warning');
                return;
            }

            showLoading('createAdvancedBtn');

            const title = document.getElementById('advancedTitle').value;
            const diagramType = document.getElementById('advancedDiagramType').value;

            fetch('/api/genomediagram/create_advanced', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    genome_length: genomeLength,
                    features: features,
                    cross_links: crossLinks,
                    title: title,
                    diagram_type: diagramType
                })
            })
            .then(response => response.json())
            .then(data => {
                hideLoading('createAdvancedBtn', '<i class="fas fa-magic"></i> Create Advanced Diagram');

                if (data.success) {
                    displayDiagram(data.diagram, 'advancedDisplay');
                    currentDiagramData = data.diagram;
                } else {
                    showAlert('Error: ' + data.error, 'danger');
                }
            })
            .catch(error => {
                hideLoading('createAdvancedBtn', '<i class="fas fa-magic"></i> Create Advanced Diagram');
                showAlert('Request failed', 'danger');
            });
        });
    }
});

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

function displayDiagram(diagramData, containerId) {
    const displayDiv = document.getElementById(containerId);

    displayDiv.innerHTML = `
        <div>
            <img src="${diagramData}" class="img-fluid" alt="Genome Diagram"
                 style="max-width: 100%; border: 1px solid #ddd; border-radius: 8px;">
        </div>
    `;
}

function exportDiagram(format) {
    if (!currentDiagramData) {
        showAlert('Please create a diagram first', 'warning');
        return;
    }

    showLoading('exportBtn');

    fetch('/api/genomediagram/export', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
            diagram_data: currentDiagramData,
            format: format
        })
    })
    .then(response => response.json())
    .then(data => {
        hideLoading('exportBtn', '<i class="fas fa-download"></i> Export Diagram');

        if (data.success) {
            // Trigger download
            const link = document.createElement('a');
            link.href = data.file_data;
            link.download = `genome_diagram.${format}`;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);

            showAlert(`Diagram exported as ${format.toUpperCase()}`, 'success');
        } else {
            showAlert('Export failed: ' + data.error, 'danger');
        }
    })
    .catch(error => {
        hideLoading('exportBtn', '<i class="fas fa-download"></i> Export Diagram');
        showAlert('Export request failed', 'danger');
    });
}

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
    // Create alert if doesn't exist
    let alertDiv = document.getElementById('genomeDiagramAlert');
    if (!alertDiv) {
        alertDiv = document.createElement('div');
        alertDiv.id = 'genomeDiagramAlert';
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

// Initialize on page load
window.addEventListener('load', function() {
    // Add one feature to basic tab
    if (document.getElementById('featuresContainer')) {
        addFeature();
    }
});
