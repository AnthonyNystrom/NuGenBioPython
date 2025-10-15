/**
 * SwissProt/UniProt Parser JavaScript
 * Handles Parse Multiple and Read Single tabs
 */

// Parse Multiple Tab
document.getElementById('parseForm').addEventListener('submit', function(e) {
    e.preventDefault();

    const fileInput = document.getElementById('parseFile');
    const maxRecords = document.getElementById('parseMaxRecords').value;

    if (!fileInput.files[0]) {
        showAlert('Please select a SwissProt file', 'warning');
        return;
    }

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('max_records', maxRecords);

    showLoading('parseBtn');

    fetch('/api/swissprot/parse', {
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
        showAlert('Please select a SwissProt file', 'warning');
        return;
    }

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);

    showLoading('readBtn');

    fetch('/api/swissprot/read', {
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

    countBadge.textContent = `${count} records`;
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
                        <strong>${record.entry_name}</strong>
                        <span class="badge bg-warning ms-2">${record.sequence_length} AA</span>
                        <span class="ms-2 text-muted">${record.description}</span>
                    </button>
                </h2>
                <div id="${accordionId}" class="accordion-collapse collapse ${isFirst ? 'show' : ''}" aria-labelledby="heading${index}" data-bs-parent="#parseAccordion">
                    <div class="accordion-body p-2">
                        <div class="row g-2">
                            <div class="col-md-6">
                                <p class="small mb-1"><strong>Accessions:</strong> ${record.accessions.join(', ')}</p>
                                <p class="small mb-1"><strong>Gene:</strong> ${record.gene_name || 'N/A'}</p>
                                <p class="small mb-1"><strong>Organism:</strong> ${record.organism}</p>
                                <p class="small mb-1"><strong>Length:</strong> ${record.sequence_length} AA</p>
                            </div>
                            <div class="col-md-6">
                                <p class="small mb-1"><strong>Keywords:</strong></p>
                                <div class="d-flex flex-wrap gap-1">
                                    ${record.keywords.slice(0, 5).map(kw => `<span class="badge bg-secondary">${kw}</span>`).join('')}
                                    ${record.keywords.length > 5 ? `<span class="badge bg-secondary">+${record.keywords.length - 5} more</span>` : ''}
                                </div>
                            </div>
                        </div>

                        ${record.features && record.features.length > 0 ? `
                            <hr class="my-2">
                            <p class="small mb-1"><strong>Features (${record.features.length}):</strong></p>
                            <div class="table-responsive">
                                <table class="table table-sm table-bordered small mb-0">
                                    <thead>
                                        <tr>
                                            <th>Type</th>
                                            <th>Location</th>
                                            <th>Description</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        ${record.features.slice(0, 5).map(feature => `
                                            <tr>
                                                <td><span class="badge bg-info">${feature.type}</span></td>
                                                <td>${feature.location}</td>
                                                <td>${feature.description || 'N/A'}</td>
                                            </tr>
                                        `).join('')}
                                        ${record.features.length > 5 ? `
                                            <tr>
                                                <td colspan="3" class="text-center text-muted">
                                                    <small>+${record.features.length - 5} more features</small>
                                                </td>
                                            </tr>
                                        ` : ''}
                                    </tbody>
                                </table>
                            </div>
                        ` : ''}

                        ${record.cross_references && record.cross_references.length > 0 ? `
                            <hr class="my-2">
                            <p class="small mb-1"><strong>Cross-References (${record.cross_references.length}):</strong></p>
                            <div class="d-flex flex-wrap gap-1">
                                ${record.cross_references.slice(0, 10).map(ref => `
                                    <span class="badge bg-secondary">${ref.database}: ${ref.id}</span>
                                `).join('')}
                                ${record.cross_references.length > 10 ? `
                                    <span class="badge bg-secondary">+${record.cross_references.length - 10} more</span>
                                ` : ''}
                            </div>
                        ` : ''}

                        ${record.comments && record.comments.length > 0 ? `
                            <hr class="my-2">
                            <p class="small mb-1"><strong>Comments (${record.comments.length}):</strong></p>
                            ${record.comments.slice(0, 3).map(comment => `
                                <div class="alert alert-info p-2 mb-1 small">
                                    ${comment}
                                </div>
                            `).join('')}
                            ${record.comments.length > 3 ? `
                                <p class="text-muted small mb-0">+${record.comments.length - 3} more comments</p>
                            ` : ''}
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
            <div class="card-header p-2 bg-warning text-white">
                <h6 class="mb-0"><i class="fas fa-protein"></i> ${record.entry_name}</h6>
            </div>
            <div class="card-body p-2">
                <div class="row g-2">
                    <div class="col-md-6">
                        <p class="small mb-1"><strong>Accessions:</strong> ${record.accessions.join(', ')}</p>
                        <p class="small mb-1"><strong>Data Class:</strong> ${record.data_class || 'N/A'}</p>
                        <p class="small mb-1"><strong>Molecule Type:</strong> ${record.molecule_type || 'N/A'}</p>
                        <p class="small mb-1"><strong>Gene:</strong> ${record.gene_name || 'N/A'}</p>
                    </div>
                    <div class="col-md-6">
                        <p class="small mb-1"><strong>Organism:</strong> ${record.organism}</p>
                        <p class="small mb-1"><strong>Length:</strong> ${record.sequence_length} AA</p>
                        <p class="small mb-1"><strong>Description:</strong> ${record.description}</p>
                    </div>
                </div>

                <hr class="my-2">
                <p class="small mb-1"><strong>Keywords (${record.keywords.length}):</strong></p>
                <div class="d-flex flex-wrap gap-1 mb-2">
                    ${record.keywords.map(kw => `<span class="badge bg-secondary">${kw}</span>`).join('')}
                </div>

                <hr class="my-2">
                <p class="small mb-1"><strong>Sequence:</strong></p>
                <div class="bg-light p-2 rounded" style="font-family: monospace; font-size: 12px; word-break: break-all; max-height: 200px; overflow-y: auto;">
                    ${record.sequence}
                </div>

                ${record.features && record.features.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Features (${record.features.length}):</strong></p>
                    <div class="table-responsive" style="max-height: 300px; overflow-y: auto;">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light sticky-top">
                                <tr>
                                    <th>Type</th>
                                    <th>Location</th>
                                    <th>Description</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.features.map(feature => `
                                    <tr>
                                        <td><span class="badge bg-info">${feature.type}</span></td>
                                        <td>${feature.location}</td>
                                        <td>${feature.description || 'N/A'}</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                ` : ''}

                ${record.cross_references && record.cross_references.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Cross-References (${record.cross_references.length}):</strong></p>
                    <div class="table-responsive" style="max-height: 300px; overflow-y: auto;">
                        <table class="table table-sm table-bordered small mb-0">
                            <thead class="table-light sticky-top">
                                <tr>
                                    <th>Database</th>
                                    <th>ID</th>
                                    <th>Info</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${record.cross_references.map(ref => `
                                    <tr>
                                        <td><span class="badge bg-secondary">${ref.database}</span></td>
                                        <td>${ref.id}</td>
                                        <td>${ref.info || 'N/A'}</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                ` : ''}

                ${record.comments && record.comments.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>Comments (${record.comments.length}):</strong></p>
                    <div style="max-height: 300px; overflow-y: auto;">
                        ${record.comments.map(comment => `
                            <div class="alert alert-info p-2 mb-1 small">
                                <strong>${comment.type}:</strong> ${comment.text}
                            </div>
                        `).join('')}
                    </div>
                ` : ''}

                ${record.references && record.references.length > 0 ? `
                    <hr class="my-2">
                    <p class="small mb-1"><strong>References (${record.references.length}):</strong></p>
                    <div style="max-height: 300px; overflow-y: auto;">
                        ${record.references.map((ref, idx) => `
                            <div class="alert alert-secondary p-2 mb-1 small">
                                <strong>[${idx + 1}]</strong> ${ref.authors || 'N/A'}<br>
                                <em>${ref.title || 'N/A'}</em><br>
                                ${ref.journal || 'N/A'}
                            </div>
                        `).join('')}
                    </div>
                ` : ''}
            </div>
        </div>
    `;

    resultsDiv.innerHTML = html;
}

function loadParseExample() {
    const swissExample = `ID   CYB_BOVIN               Reviewed;         378 AA.
AC   P00157;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DT   21-JUL-1986, sequence version 1.
DT   11-DEC-2019, entry version 151.
DE   RecName: Full=Cytochrome b;
GN   Name=CYTB; Synonyms=COB;
OS   Bos taurus (Bovine).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
OC   Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia;
OC   Pecora; Bovidae; Bovinae; Bos.
KW   Complete proteome; Electron transport; Heme; Iron; Membrane;
KW   Metal-binding; Mitochondrion; Mitochondrion inner membrane;
KW   Respiratory chain; Transmembrane; Transmembrane helix; Transport;
KW   Ubiquinone.
DR   EMBL; V00159; CAA23506.1; -; Genomic_DNA.
DR   PIR; A00729; CYBO.
DR   GO; GO:0005739; C:mitochondrion; IEA:UniProtKB-SubCell.
CC   -!- FUNCTION: Component of the ubiquinol-cytochrome c reductase
CC       complex (complex III or cytochrome b-c1 complex).
SQ   SEQUENCE   378 AA;  42851 MW;  61FE8C26B13BAC23 CRC64;
     MTNIRKSHPL MKLVDLPGVN PSLMQNIVIH SDQTTMWLIW NKKVALYFFP MLILTGMLVF
     LQTNGKLNRD TVGSNIIVPS NIKVSSIALW LWGGFSVDKA TLNRFFAFHF ILPFTMVALA
     GIHLTFLHES GSNNPLGFTS DSDKIPFHPY YTVKDLLGIL ILILLLLLLA LLSPDMLGDP
     DNYMPADPLN TPLHIKPEWF FLFAYAILRS VPNKLGGVLA LFLSIVILLA MPFLHTSKHR
     GMMFRPLSQA LFWILVADLL VLTWIGSQPV EYPYTIIGQM ASILYFSILA FLPIAAGIMM
     VLLSGVASAL LKDLLLFGGM LKVGYSSQLE LEYTWMLQAK ADIWAVGIQE QEGIQQKLYM
     LGIMVGFAFF LVLMPLRAG
//
ID   INS_HUMAN               Reviewed;         110 AA.
AC   P01308;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DT   21-JUL-1986, sequence version 1.
DT   11-DEC-2019, entry version 221.
DE   RecName: Full=Insulin;
GN   Name=INS;
OS   Homo sapiens (Human).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
OC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
OC   Catarrhini; Hominidae; Homo.
KW   Amyloidosis; Diabetes mellitus; Disease mutation; Hormone;
KW   Pharmaceutical; Secreted; Signal.
DR   EMBL; V00565; CAA23839.1; -; Genomic_DNA.
DR   PIR; A00022; INHM.
DR   GO; GO:0005576; C:extracellular region; IEA:UniProtKB-SubCell.
CC   -!- FUNCTION: Insulin decreases blood glucose concentration by
CC       regulating carbohydrate metabolism in muscle, fat and liver.
SQ   SEQUENCE   110 AA;  12171 MW;  C514A3AA7AFF7783 CRC64;
     MALWMRLLPL LALLALWGPD PAAAFVNQHL CGSHLVEALY LVCGERGFFY TPKTRREAED
     LQVGQVELGG GPGAGSLQPL ALEGSLQKRG IVEQCCTSIC SLYQLENYCN
//
ID   HBA_HUMAN               Reviewed;         142 AA.
AC   P69905;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DT   23-JAN-2007, sequence version 2.
DT   11-DEC-2019, entry version 160.
DE   RecName: Full=Hemoglobin subunit alpha;
GN   Name=HBA1; Synonyms=HBA;
OS   Homo sapiens (Human).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
OC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
OC   Catarrhini; Hominidae; Homo.
KW   Heme; Iron; Metal-binding; Oxygen transport; Transport.
DR   EMBL; V00488; CAA23751.1; -; Genomic_DNA.
DR   PIR; A01543; HAHU.
DR   GO; GO:0005840; C:ribosome; IEA:UniProtKB-SubCell.
CC   -!- FUNCTION: Involved in oxygen transport from the lung to the
CC       various peripheral tissues.
SQ   SEQUENCE   142 AA;  15258 MW;  88E7BFD45F7D5B54 CRC64;
     MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS HGSAQVKGHG
     KKVADALTNA VAHVDDMPNA LSALSDLHAH KLRVDPVNFK LLSHCLLVTL AAHLPAEFTP
     AVHASLDKFL ASVSTVLTSK YR
//`;

    const blob = new Blob([swissExample], { type: 'text/plain' });
    const file = new File([blob], 'example_cytb_bovin.dat', { type: 'text/plain' });

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

    fetch('/api/swissprot/parse', {
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
    const swissExample = `ID   CYB_BOVIN               Reviewed;         378 AA.
AC   P00157;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DT   21-JUL-1986, sequence version 1.
DT   11-DEC-2019, entry version 151.
DE   RecName: Full=Cytochrome b;
GN   Name=CYTB; Synonyms=COB;
OS   Bos taurus (Bovine).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
OC   Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia;
OC   Pecora; Bovidae; Bovinae; Bos.
KW   Complete proteome; Electron transport; Heme; Iron; Membrane;
KW   Metal-binding; Mitochondrion; Mitochondrion inner membrane;
KW   Respiratory chain; Transmembrane; Transmembrane helix; Transport;
KW   Ubiquinone.
DR   EMBL; V00159; CAA23506.1; -; Genomic_DNA.
DR   PIR; A00729; CYBO.
DR   GO; GO:0005739; C:mitochondrion; IEA:UniProtKB-SubCell.
CC   -!- FUNCTION: Component of the ubiquinol-cytochrome c reductase
CC       complex (complex III or cytochrome b-c1 complex), which is a
CC       respiratory chain that generates an electrochemical potential
CC       coupled to ATP synthesis.
CC   -!- SUBCELLULAR LOCATION: Mitochondrion inner membrane; Multi-pass
CC       membrane protein.
RN   [1]
RP   NUCLEOTIDE SEQUENCE [GENOMIC DNA].
RA   Anderson S., de Bruijn M.H.L., Coulson A.R., Eperon I.C., Sanger F.,
RA   Young I.G.;
RT   "Complete sequence of bovine mitochondrial DNA. Conserved features of
RT   the mammalian mitochondrial genome.";
RL   J. Mol. Biol. 156:683-717(1982).
SQ   SEQUENCE   378 AA;  42851 MW;  61FE8C26B13BAC23 CRC64;
     MTNIRKSHPL MKLVDLPGVN PSLMQNIVIH SDQTTMWLIW NKKVALYFFP MLILTGMLVF
     LQTNGKLNRD TVGSNIIVPS NIKVSSIALW LWGGFSVDKA TLNRFFAFHF ILPFTMVALA
     GIHLTFLHES GSNNPLGFTS DSDKIPFHPY YTVKDLLGIL ILILLLLLLA LLSPDMLGDP
     DNYMPADPLN TPLHIKPEWF FLFAYAILRS VPNKLGGVLA LFLSIVILLA MPFLHTSKHR
     GMMFRPLSQA LFWILVADLL VLTWIGSQPV EYPYTIIGQM ASILYFSILA FLPIAAGIMM
     VLLSGVASAL LKDLLLFGGM LKVGYSSQLE LEYTWMLQAK ADIWAVGIQE QEGIQQKLYM
     LGIMVGFAFF LVLMPLRAG
//`;

    const blob = new Blob([swissExample], { type: 'text/plain' });
    const file = new File([blob], 'example_cytb_bovin.dat', { type: 'text/plain' });

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

    fetch('/api/swissprot/read', {
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
