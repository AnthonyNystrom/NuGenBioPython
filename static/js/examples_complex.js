/* Complex example data for every major tool section.
 *
 * Each function is paired with an existing simple `loadXxxExample` fn in
 * the tool's own JS. These complex examples use realistic biological data
 * (human gene snippets, multi-sequence alignments, richer matrices, etc.)
 * so users can see the tool perform on substantive input.
 *
 * Function naming: `loadXxxExampleComplex` — paired with `loadXxxExample`
 * (simple) or the tool's existing loader.
 */
(function () {
    'use strict';

    // -------- Sequence analysis --------
    // Human TP53 exon 4 (~500bp) — varied GC, real codons, has ORFs.
    window.loadExampleComplex = function () {
        var seqType = document.getElementById('sequenceType');
        var seqInput = document.getElementById('sequenceInput');
        if (seqType) seqType.value = 'dna';
        if (seqInput) {
            seqInput.value =
                'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA' +
                'GACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATG' +
                'GATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCA' +
                'GATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCT' +
                'ACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAG' +
                'AAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAG' +
                'TCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACC' +
                'TGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATG';
        }
    };

    // -------- Pairwise alignment --------
    // Two divergent mammalian HBB fragments (human vs mouse beta-globin) —
    // shows gaps, mismatches, a realistic protein-coding comparison.
    window.loadExampleAlignmentComplex = function () {
        var s1 = document.getElementById('sequence1');
        var s2 = document.getElementById('sequence2');
        if (s1) s1.value =
            'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK' +
            'VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG' +
            'KEFTPPVQAAYQKVVAGVANALAHKYH';
        if (s2) s2.value =
            'MVHLTDAEKAAVSGLWGKVNADEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPK' +
            'VKAHGKKVITAFNEGLKNLDNLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLG' +
            'KEFTPCAQAAFQKVVAGVASALAHKYH';
        var mode = document.getElementById('alignmentMode');
        if (mode) mode.value = 'global';
        var matrix = document.getElementById('substitutionMatrix');
        if (matrix) { matrix.value = 'BLOSUM62'; matrix.dispatchEvent(new Event('change')); }
    };

    // -------- BLAST --------
    // Full human beta-globin mRNA (~600bp). Complex = actual searchable gene,
    // expected to return high-identity hits against nt.
    window.loadBlastExampleComplex = function () {
        var seq = document.getElementById('querySequence');
        if (seq) seq.value =
            'ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGA' +
            'GGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGC' +
            'AGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATG' +
            'CTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGC' +
            'TCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT' +
            'CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA' +
            'CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCA' +
            'CTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACT';
        var prog = document.getElementById('blastProgram');
        var db = document.getElementById('database');
        var max = document.getElementById('maxHits');
        if (prog) prog.value = 'blastn';
        if (db) db.value = 'nt';
        if (max) max.value = '25';
        var file = document.getElementById('sequenceFile');
        if (file) file.value = '';
        if (window.showAlert) window.showAlert('Loaded: human β-globin mRNA (full coding region).', 'success');
    };

    // -------- Motif finder --------
    // TATA-box variants. Each sequence contains TATAAA-like consensus at
    // different positions — the motif build + search should find them.
    window.loadMotifExampleComplex = function () {
        var mot = document.getElementById('motifSequences');
        if (mot) mot.value = [
            'TATAAA',
            'TATATAA',
            'TATAAAA',
            'TATAAATA',
            'CTATAAAA',
            'TATATATA',
            'GTATAAAT',
            'ATATAAAG',
            'TATAAAAG',
            'CTATAAAT',
        ].join('\n');
        var search = document.getElementById('searchSequence');
        if (search) search.value =
            // Synthetic promoter region with embedded TATA-like boxes
            'CCGGCGCGTGACCAGTTCCGCTATAAAAGCGCGCCGCGACGTGGATTCTATAAAGCCGCG' +
            'TGAGCTATAAATCCCGCGAGTACGCGTGCCTATATAAAGCCCGCGAGTGAGCGTGCTATA' +
            'TATCCCGCGTATGAGCGTGAGTCGTATAAAGCCCGCGCCCCGAGTGCCTATATATGAGCG';
    };

    // -------- Structure --------
    // Complex: fetch a real mid-sized PDB direct from RCSB. Porcine insulin
    // (4INS) is ~800 atoms across two chains with disulfide bridges and
    // clear secondary structure, so the cartoon view is immediately
    // informative. Uses the same client-side fetch path as the manual
    // "Fetch" button — no backend endpoint required.
    window.loadExampleComplexStructure = function () {
        var pdbId = '4INS';
        var labelText = 'Porcine insulin (4INS)';
        if (window.showLoading) window.showLoading('parseStructureBtn');

        fetch('https://files.rcsb.org/download/' + pdbId + '.pdb')
            .then(function (r) {
                if (!r.ok) throw new Error('PDB ' + pdbId + ' not found (HTTP ' + r.status + ')');
                return r.text();
            })
            .then(function (pdbText) {
                if (window.NuGenStructureViewer && window.NuGenStructureViewer.load) {
                    window.NuGenStructureViewer.load(pdbText, labelText);
                }
                var blob = new Blob([pdbText], { type: 'text/plain' });
                var file = new File([blob], pdbId + '.pdb', { type: 'text/plain' });
                window.currentStructureFile = file;
                var fd = new FormData();
                fd.append('file', file);
                return fetch('/api/structure/parse', { method: 'POST', body: fd });
            })
            .then(function (r) { return r.json(); })
            .then(function (data) {
                if (window.hideLoading) window.hideLoading('parseStructureBtn', '<i class="fas fa-cube me-2"></i>Parse');
                if (data && data.success && typeof displayStructureInfo === 'function') {
                    displayStructureInfo(data.structure_info);
                    window.structureInfo = data.structure_info;
                    var btn = document.getElementById('advancedAnalysisBtn');
                    if (btn) btn.disabled = false;
                }
            })
            .catch(function (err) {
                if (window.hideLoading) window.hideLoading('parseStructureBtn', '<i class="fas fa-cube me-2"></i>Parse');
                if (window.showAlert) window.showAlert(window.friendlyError ? window.friendlyError(err, 'pdb') : ('Failed: ' + err.message), 'danger');
            });
    };

    // -------- GenomeDiagram --------
    // Helper: assign a color to a <select>, falling back to 'blue' when the
    // requested color isn't in the select's option list (otherwise the browser
    // silently sets value to '' and the backend raises "Invalid color value ''").
    function setColorSafe(select, color) {
        if (!select) return;
        var opts = Array.from(select.options).map(function (o) { return o.value; });
        select.value = opts.indexOf(color) >= 0 ? color : 'blue';
    }

    // Complex basic: real pUC19 plasmid (2686 bp) with 11 features at their
    // actual coordinates. Covers most of the plasmid so the diagram is dense
    // and visually informative rather than leaving large empty stretches.
    window.loadBasicExampleComplexGenome = function () {
        if (typeof clearFeatures === 'function') clearFeatures();
        var len = document.getElementById('genomeLength');
        var title = document.getElementById('diagramTitle');
        if (len) len.value = '2686';
        if (title) title.value = 'pUC19 plasmid (real coordinates)';
        var features = [
            { name: 'AmpR signal pep', start: 1,    end: 80,   color: 'pink' },
            { name: 'AmpR (β-lactamase)', start: 80,   end: 940,  color: 'red' },
            { name: 'AmpR promoter', start: 941,  end: 1045, color: 'green' },
            { name: 'ori',            start: 1151, end: 1739, color: 'brown' },
            { name: 'bom',            start: 1887, end: 2029, color: 'gray' },
            { name: 'M13 Fwd primer', start: 2183, end: 2200, color: 'purple' },
            { name: 'lac promoter',   start: 2209, end: 2316, color: 'green' },
            { name: 'lac operator',   start: 2317, end: 2333, color: 'orange' },
            { name: 'MCS',            start: 2346, end: 2402, color: 'cyan' },
            { name: 'lacZ alpha',     start: 2358, end: 2613, color: 'blue' },
            { name: 'M13 Rev primer', start: 2627, end: 2644, color: 'purple' }
        ];
        features.forEach(function (f) {
            if (typeof addFeature === 'function') addFeature();
            var row = document.querySelector('.feature-item:last-child');
            if (!row) return;
            row.querySelector('.feature-name').value = f.name;
            row.querySelector('.feature-start').value = f.start;
            row.querySelector('.feature-end').value = f.end;
            setColorSafe(row.querySelector('.feature-color'), f.color);
        });
        // Prefer circular view for a real plasmid — much more accurate.
        var typeSel = document.getElementById('diagramType');
        if (typeSel) typeSel.value = 'circular';
        if (typeof updateFeatureSummary === 'function') updateFeatureSummary();
    };

    // Complex advanced: E. coli lac + trp + rRNA regions (25 kb) using both
    // strands and varied sigils so the rendered diagram is rich and readable.
    window.loadAdvancedExampleComplexGenome = function () {
        if (typeof clearAdvancedFeatures === 'function') clearAdvancedFeatures();
        var len = document.getElementById('advancedGenomeLength');
        var title = document.getElementById('advancedTitle');
        if (len) len.value = '25000';
        if (title) title.value = 'E. coli operon cluster (complex)';
        // Genes use ARROW sigils + proper strand; regulatory elements use BOX.
        var features = [
            { name: 'lacI',            start: 1000,  end: 2080,  color: 'darkblue', strand: -1, sigil: 'ARROW' },
            { name: 'CAP site',        start: 2100,  end: 2140,  color: 'orange',   strand: 0,  sigil: 'BOX' },
            { name: 'lac promoter',    start: 2150,  end: 2250,  color: 'green',    strand: 1,  sigil: 'BOX' },
            { name: 'lacZ',            start: 2270,  end: 5300,  color: 'blue',     strand: 1,  sigil: 'ARROW' },
            { name: 'lacY',            start: 5350,  end: 6600,  color: 'blue',     strand: 1,  sigil: 'ARROW' },
            { name: 'lacA',            start: 6650,  end: 7260,  color: 'blue',     strand: 1,  sigil: 'ARROW' },
            { name: 'lac terminator',  start: 7280,  end: 7400,  color: 'red',      strand: 0,  sigil: 'BOX' },
            { name: 'trp leader',      start: 9960,  end: 10000, color: 'teal',     strand: 1,  sigil: 'BOX' },
            { name: 'trp operator',    start: 10000, end: 10050, color: 'orange',   strand: 0,  sigil: 'BOX' },
            { name: 'trpE',            start: 10100, end: 11500, color: 'purple',   strand: 1,  sigil: 'ARROW' },
            { name: 'trpD',            start: 11550, end: 13090, color: 'purple',   strand: 1,  sigil: 'ARROW' },
            { name: 'trpC',            start: 13150, end: 14700, color: 'purple',   strand: 1,  sigil: 'ARROW' },
            { name: 'trpB',            start: 14750, end: 15950, color: 'purple',   strand: 1,  sigil: 'ARROW' },
            { name: 'trpA',            start: 16000, end: 16800, color: 'purple',   strand: 1,  sigil: 'ARROW' },
            { name: 'rrsB (16S rRNA)', start: 18000, end: 19540, color: 'brown',    strand: 1,  sigil: 'ARROW' },
            { name: 'rrfB (5S rRNA)',  start: 19600, end: 19720, color: 'brown',    strand: 1,  sigil: 'ARROW' },
            { name: 'IS5 transposon',  start: 21000, end: 22200, color: 'gray',     strand: -1, sigil: 'ARROW' },
            { name: 'tRNA-Ala',        start: 22900, end: 22977, color: 'pink',     strand: -1, sigil: 'ARROW' },
            { name: 'tRNA-Ile',        start: 23100, end: 23177, color: 'pink',     strand: -1, sigil: 'ARROW' },
            { name: 'tRNA-Lys',        start: 23300, end: 23377, color: 'pink',     strand: -1, sigil: 'ARROW' },
            { name: 'tRNA-Val',        start: 23500, end: 23577, color: 'pink',     strand: -1, sigil: 'ARROW' },
            { name: 'aspS gene',       start: 24000, end: 24900, color: 'cyan',     strand: -1, sigil: 'ARROW' }
        ];
        features.forEach(function (f) {
            if (typeof addAdvancedFeature === 'function') addAdvancedFeature();
            var row = document.querySelector('.advanced-feature-item:last-child');
            if (!row) return;
            var setVal = function (sel, v) { var el = row.querySelector(sel); if (el) el.value = v; };
            setVal('.advanced-feature-name',   f.name);
            setVal('.advanced-feature-start',  f.start);
            setVal('.advanced-feature-end',    f.end);
            setVal('.advanced-feature-strand', f.strand);
            setVal('.advanced-feature-sigil',  f.sigil);
            setColorSafe(row.querySelector('.advanced-feature-color'), f.color);
        });
    };

    // -------- Restriction (Basic) --------
    // pUC19-like multi-site sequence (~600bp) with many common enzyme sites.
    window.loadBasicExampleComplex = function () {
        var seq = document.getElementById('basicSequence');
        if (seq) seq.value =
            'GAATTCAAGCTTATCGATCGAATTCCTGCAGGGATCCAAGCTTTCTAGATGCATGCCTGCAGGAATTC' +
            'GGTACCGCGGCCGCTCTAGAACTAGTGGATCCCCCGGGCTGCAGGAATTCGATATCAAGCTTATCGATA' +
            'CCGTCGACCTCGAGGGGGGGCCCGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTGCGCGCTTG' +
            'GCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGA' +
            'GCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGC' +
            'TCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGG' +
            'AGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGG' +
            'CTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCA';
        // Pre-select common enzymes (trigger click on "Common" button)
        setTimeout(function () {
            var btn = document.querySelector('[data-action="selectCommonBasicEnzymes"]');
            if (btn) btn.click();
        }, 150);
    };

    // -------- Restriction advanced --------
    window.loadAdvancedExampleComplex = function () {
        var seq = document.getElementById('advancedSequence');
        if (seq) seq.value =
            // Long complex sequence with many potential cuts
            'GAATTCGCGGCCGCTCTAGAACTAGTGGATCCCCCGGGCTGCAGGAATTCGATATCAAGCTTATCGATA' +
            'CCGTCGACCTCGAGGGGGGGCCCGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTCTGAGAGTG' +
            'CACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCA' +
            'TTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAA' +
            'GGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACG';
        var filt = document.getElementById('filterType');
        if (filt) filt.value = 'unique';
        var allEnz = document.getElementById('useAllEnzymes');
        if (allEnz) allEnz.checked = true;
        setTimeout(function () {
            var f = document.getElementById('advancedAnalysisForm');
            if (f) f.dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        }, 400);
    };

    // -------- Restriction map --------
    window.loadMapExampleComplex = function () {
        var seq = document.getElementById('mapSequence');
        var enz = document.getElementById('mapEnzymes');
        if (seq) seq.value =
            'GAATTCAAGCTTATCGATCGAATTCCTGCAGGGATCCAAGCTTTCTAGATGCATGCCTGCAGGAATTC' +
            'GGTACCGCGGCCGCTCTAGAACTAGTGGATCCCCCGGGCTGCAGGAATTCGATATCAAGCTTATCGATA' +
            'CCGTCGACCTCGAGGGGGGGCCCGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTGCGCGCTT';
        if (enz) enz.value = 'EcoRI, BamHI, HindIII, PstI, XhoI, SalI, NotI, XbaI, SmaI, KpnI';
        setTimeout(function () {
            var f = document.getElementById('mapForm');
            if (f) f.dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        }, 400);
    };

    // -------- Restriction compatible --------
    window.loadCompatibleExampleComplex = function () {
        var enz = document.getElementById('compatibleEnzymes');
        if (enz) enz.value = 'BamHI, BglII, XbaI, SpeI, NheI, AvrII, EcoRI, MfeI, NcoI, PciI, ApaI, PspOMI';
        setTimeout(function () {
            var f = document.getElementById('compatibleEndsForm');
            if (f) f.dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        }, 300);
    };

    // -------- Features ORF --------
    // Real beta-actin-like coding region with multiple frames/ORFs
    window.loadORFExampleComplex = function () {
        var seq = document.getElementById('orfSequence');
        if (seq) seq.value =
            'ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGC' +
            'TTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCAC' +
            'CAGGGCGTGATGGTGGGCATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCCCAGAGC' +
            'AAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGAC' +
            'GACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAG' +
            'CACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACC' +
            'CAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTA' +
            'TCCCTGTACGCCTCTGGCCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACC' +
            'CACACTGTGCCCATCTACGAGGGGTATGCCCTCCCCCATGCCATCCTGCGTCTGGACCTG';
        var minLen = document.getElementById('minLength');
        if (minLen) minLen.value = '60';
        var strand = document.getElementById('strand');
        if (strand) strand.value = 'both';
    };

    // -------- Features extract --------
    window.loadExtractExampleComplex = function () {
        var seq = document.getElementById('extractSequence');
        if (seq) seq.value =
            'GAATTCATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAG' +
            'GCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCC' +
            'AGGCACCAGGGCGTGATGGTGGGCATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCC' +
            'CAGAGCAAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGTAA';
        var start = document.getElementById('extractStart');
        var end = document.getElementById('extractEnd');
        var strand = document.getElementById('extractStrand');
        var trans = document.getElementById('translateCheck');
        if (start) start.value = '7';
        if (end) end.value = '205';
        if (strand) strand.value = '1';
        if (trans) trans.checked = true;
    };

    // -------- Features compound --------
    window.loadCompoundExampleComplex = function () {
        var seq = document.getElementById('compoundSequence');
        if (seq) seq.value =
            'ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGC' +
            'TTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCAC' +
            'CAGGGCGTGATGGTGGGCATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCCCAGAGC' +
            'AAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGAC' +
            'GACATGGAGAAAATCTGGCACCACACCTTCTAA';
        var loc = document.getElementById('compoundLocations');
        if (loc) loc.value = '1,60,1\n65,130,1\n135,200,-1\n210,260,1';
    };

    // -------- Phylogeny parse --------
    // Newick with bootstrap values, 15 mammal taxa
    window.loadExampleTreeComplex = function () {
        var sel = document.getElementById('treeFormat');
        if (sel) sel.value = 'newick';
        var inp = document.getElementById('treeString');
        if (inp) inp.value =
            '(((((((human:0.006550,chimp:0.006655)100:0.033537,(gorilla:0.006469,orangutan:0.019601)88:0.001015)92:0.004203,' +
            '(gibbon:0.024680,siamang:0.028370)83:0.005402)76:0.007212,(macaque:0.037470,baboon:0.037470)95:0.008030)61:0.012870,' +
            '((squirrel_monkey:0.042300,marmoset:0.042300)90:0.024270,tamarin:0.066570)78:0.019100)55:0.026430,' +
            '(lemur:0.089800,loris:0.089800)85:0.025260)99:0.048440,(mouse:0.145500,rat:0.145500)100:0.019720):0;';
    };

    // -------- Phylogeny build --------
    // Larger alignment — 6 sequences, 180bp
    window.loadExampleAlignmentComplexPhylo = function () {
        var af = document.getElementById('alignmentFormat');
        if (af) af.value = 'fasta';
        var inp = document.getElementById('alignmentString');
        if (inp) inp.value = [
            '>Human',
            'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
            '>Chimp',
            'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
            '>Gorilla',
            'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLNHLDNLKGTFSQLSELHCDKLHVDPENFKLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
            '>Cow',
            'MLTAEEKAAVTAFWGKVKVDEVGGEALGRLLVVYPWTQRFFESFGDLSTADAVMNNPKVKAHGKKVLDSFSNGMKHLDDLKGTFAALSELHCDKLHVDPENFRLLGNVLVVVLARNFGKEFTPVLQADFQKVVAGVANALAHRYH',
            '>Mouse',
            'MVHLTDAEKSAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRYFDKFGNLSSASAIMGNPKVKAHGAKVITAFNDGLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH',
            '>Zebrafish',
            'MVEWTDAERTAILGLWGKLNIDEIGPQALARCLIVYPWTQRYFATFGNLSSPAAIMGNPKVAAHGRTVMGGLERAIKNMDNIKNTYAALSVMHSEKLHVDPDNFRLLADCITVCAAMKFGPSGFNADVQEAWQKFLAVVVSALCRQYH',
        ].join('\n');
        var method = document.getElementById('buildMethod');
        var model = document.getElementById('distanceModel');
        if (method) method.value = 'upgma';
        if (model) model.value = 'identity';
    };

    // -------- Clustering --------
    // 10 samples × 5 features with clear cluster structure
    window.loadClusteringExampleComplex = function () {
        var mat = document.getElementById('dataMatrix');
        var names = document.getElementById('sampleNames');
        if (mat) mat.value = [
            '12.3, 8.2, 15.7, 3.1, 9.8',
            '11.8, 8.5, 16.2, 3.3, 10.1',
            '13.1, 7.9, 15.3, 2.8, 9.4',
            '12.7, 8.1, 15.9, 3.0, 9.7',
            '2.4, 22.5, 4.2, 18.3, 1.2',
            '2.1, 23.1, 4.5, 18.8, 1.5',
            '2.7, 22.8, 3.9, 17.9, 1.3',
            '8.1, 5.2, 9.3, 10.4, 14.7',
            '7.9, 4.9, 8.8, 10.1, 15.2',
            '8.4, 5.4, 9.5, 10.7, 14.4',
        ].join('\n');
        if (names) names.value = 'BRCA1,BRCA2,BRCA3,BRCA4,HER2_1,HER2_2,HER2_3,Ctrl_1,Ctrl_2,Ctrl_3';
        var method = document.getElementById('clusterMethod');
        var k = document.getElementById('nClusters');
        if (method) { method.value = 'kmeans'; method.dispatchEvent(new Event('change')); }
        if (k) k.value = '3';
    };

    // -------- Popgen --------
    // Complex: 4 populations × 10 loci (backend generates via /popgen/load_example_complex).
    window.loadPopgenExampleComplex = function () {
        if (window.showLoading) window.showLoading('parsePopgenBtn');
        fetch('/api/popgen/load_example_complex', { method: 'POST' })
            .then(function (r) { return r.json(); })
            .then(function (data) {
                if (window.hideLoading) window.hideLoading('parsePopgenBtn', '<i class="fas fa-dna"></i> Parse');
                if (data && data.success && typeof displayPopgenResults === 'function') {
                    displayPopgenResults(data.results);
                    if (window.showAlert) window.showAlert('Loaded: 4 populations × 10 loci example.', 'success');
                } else if (window.showAlert) {
                    window.showAlert(window.friendlyError ? window.friendlyError(data && data.error, 'server') : 'Failed to load example', 'warning');
                }
            })
            .catch(function (err) {
                if (window.hideLoading) window.hideLoading('parsePopgenBtn', '<i class="fas fa-dna"></i> Parse');
                if (window.showAlert) window.showAlert(window.friendlyError ? window.friendlyError(err, 'server') : 'Failed to load example', 'danger');
            });
    };

    // -------- SwissProt / UniGene --------
    // Complex: multi-record example (ALB + HBB). Uses the same Blob + FormData
    // flow as the simple loadParseExample — bypasses the file input entirely.
    window.loadParseExampleComplex = function () {
        var parseFile = document.getElementById('parseFile');
        if (!parseFile) return;  // page doesn't host this form
        var text =
`ID   ALBU_HUMAN              Reviewed;         609 AA.
AC   P02768; A6NF33; B3GQS0; B7WNR0;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DT   11-DEC-2019, entry version 212.
DE   RecName: Full=Albumin;
DE   AltName: Full=Serum albumin;
GN   Name=ALB;
OS   Homo sapiens (Human).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia.
KW   3D-structure; Acetylation; Alternative splicing;
KW   Direct protein sequencing; Disulfide bond; Glycation; Glycoprotein;
KW   Metal-binding; Phosphoprotein; Plasma; Reference proteome;
KW   Secreted; Signal; Transport.
DR   EMBL; V00494; CAA23758.1; -; mRNA.
DR   PIR; A93762; ABHU.
DR   GO; GO:0005615; C:extracellular space; IEA:UniProtKB-SubCell.
CC   -!- FUNCTION: Serum albumin, the main protein of plasma, has a good
CC       binding capacity for water, Ca(2+), Na(+), K(+), fatty acids,
CC       hormones, bilirubin and drugs.
FT   CHAIN           25..609
FT                   /note="Albumin"
FT   DISULFID        77..86
FT   DISULFID        99..115
SQ   SEQUENCE   609 AA;  69367 MW;  B6F6C2F3C8A8A8AB CRC64;
     MKWVTFISLL FLFSSAYSRG VFRRDAHKSE VAHRFKDLGE ENFKALVLIA FAQYLQQCPF
     EDHVKLVNEV TEFAKTCVAD ESAENCDKSL HTLFGDKLCT VATLRETYGE MADCCAKQEP
//
ID   HBB_HUMAN               Reviewed;         147 AA.
AC   P68871; A4GX73; B2ZUE0; P02023; Q13852;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=Hemoglobin subunit beta;
DE   AltName: Full=Beta-globin;
GN   Name=HBB;
OS   Homo sapiens (Human).
OC   Eukaryota; Metazoa; Chordata; Mammalia.
KW   3D-structure; Acetylation; Antimicrobial; Blood group antigen;
KW   Defense response; Disease variant; Glycation; Heme; Iron;
KW   Metal-binding; Oxygen transport; Phosphoprotein; Transport.
DR   EMBL; U01317; AAB59424.1; -; Genomic_DNA.
DR   GO; GO:0005833; C:hemoglobin complex; IDA:UniProtKB.
CC   -!- FUNCTION: Involved in oxygen transport from the lung to the
CC       various peripheral tissues.
FT   CHAIN           2..147
FT                   /note="Hemoglobin subunit beta"
FT   METAL           93
FT                   /note="Iron (heme distal ligand)"
SQ   SEQUENCE   147 AA;  15998 MW;  A31F6F1D08A6CE4A CRC64;
     MVHLTPEEKS AVTALWGKVN VDEVGGEALG RLLVVYPWTQ RFFESFGDLS TPDAVMGNPK
     VKAHGKKVLG AFSDGLAHLD NLKGTFATLS ELHCDKLHVD PENFRLLGNV LVCVLAHHFG
     KEFTPPVQAA YQKVVAGVAN ALAHKYH
//`;
        var blob = new Blob([text], { type: 'text/plain' });
        var file = new File([blob], 'example_multi.dat', { type: 'text/plain' });
        try {
            var dt = new DataTransfer();
            dt.items.add(file);
            parseFile.files = dt.files;
        } catch (e) { /* fall through to direct post */ }

        var fd = new FormData();
        fd.append('file', file);
        var maxEl = document.getElementById('parseMaxRecords');
        if (maxEl) fd.append('max_records', maxEl.value);

        if (window.showLoading) window.showLoading('parseBtn');
        fetch('/api/swissprot/parse', { method: 'POST', body: fd })
            .then(function (r) { return r.json(); })
            .then(function (data) {
                if (window.hideLoading) window.hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');
                if (data && data.success && typeof displayParseResults === 'function') {
                    displayParseResults(data.records, data.count);
                } else if (window.showAlert) {
                    window.showAlert(window.friendlyError ? window.friendlyError(data && data.error, 'uniprot') : 'Failed', 'warning');
                }
            })
            .catch(function (err) {
                if (window.hideLoading) window.hideLoading('parseBtn', '<i class="fas fa-list me-2"></i>Parse Records');
                if (window.showAlert) window.showAlert(window.friendlyError ? window.friendlyError(err, 'uniprot') : 'Failed', 'danger');
            });
    };

    // -------- SearchIO --------
    // Complex: multi-query, multi-hit BLAST XML — parses immediately.
    window.loadParseExampleComplexSearchIO = function () {
        var sel = document.getElementById('parseFormat');
        if (sel) sel.value = 'blast-xml';
        var xml = `<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.14.0+</BlastOutput_version>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>sp|P68871|HBB_HUMAN hemoglobin subunit beta</BlastOutput_query-def>
  <BlastOutput_query-len>147</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>HBB_HUMAN</Iteration_query-def>
      <Iteration_query-len>147</Iteration_query-len>
      <Iteration_hits>
        <Hit><Hit_num>1</Hit_num><Hit_id>sp|P02088|HBB_MOUSE</Hit_id><Hit_def>Hemoglobin beta chain (Mouse)</Hit_def><Hit_len>147</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>271.2</Hsp_bit-score><Hsp_evalue>1e-92</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>147</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>147</Hsp_hit-to><Hsp_identity>116</Hsp_identity><Hsp_positive>128</Hsp_positive><Hsp_align-len>147</Hsp_align-len></Hsp></Hit_hsps>
        </Hit>
        <Hit><Hit_num>2</Hit_num><Hit_id>sp|P02091|HBB_RAT</Hit_id><Hit_def>Hemoglobin beta-1 chain (Rat)</Hit_def><Hit_len>147</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>268.8</Hsp_bit-score><Hsp_evalue>4e-91</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>147</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>147</Hsp_hit-to><Hsp_identity>115</Hsp_identity><Hsp_positive>127</Hsp_positive><Hsp_align-len>147</Hsp_align-len></Hsp></Hit_hsps>
        </Hit>
        <Hit><Hit_num>3</Hit_num><Hit_id>sp|P02113|HBB_CHICK</Hit_id><Hit_def>Hemoglobin beta chain (Chicken)</Hit_def><Hit_len>146</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>215.7</Hsp_bit-score><Hsp_evalue>6e-70</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>146</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>146</Hsp_hit-to><Hsp_identity>98</Hsp_identity><Hsp_positive>115</Hsp_positive><Hsp_align-len>146</Hsp_align-len></Hsp></Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>Query_2</Iteration_query-ID>
      <Iteration_query-def>sp|P68871|HBA_HUMAN hemoglobin subunit alpha</Iteration_query-def>
      <Iteration_query-len>142</Iteration_query-len>
      <Iteration_hits>
        <Hit><Hit_num>1</Hit_num><Hit_id>sp|P01942|HBA_MOUSE</Hit_id><Hit_def>Hemoglobin alpha chain (Mouse)</Hit_def><Hit_len>142</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>268.1</Hsp_bit-score><Hsp_evalue>9e-91</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>142</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>142</Hsp_hit-to><Hsp_identity>119</Hsp_identity><Hsp_positive>129</Hsp_positive><Hsp_align-len>142</Hsp_align-len></Hsp></Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>`;
        var blob = new Blob([xml], { type: 'text/xml' });
        var file = new File([blob], 'example_multiquery.xml', { type: 'text/xml' });
        var fd = new FormData();
        fd.append('file', file);
        fd.append('format', 'blast-xml');

        if (window.showLoading) window.showLoading('parseBtn');
        fetch('/api/searchio/parse', { method: 'POST', body: fd })
            .then(function (r) { return r.json(); })
            .then(function (data) {
                if (window.hideLoading) window.hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
                if (data && data.success && typeof displayParseResults === 'function') {
                    displayParseResults(data.results, data.count);
                } else if (window.showAlert) {
                    window.showAlert(window.friendlyError ? window.friendlyError(data && data.error, 'server') : 'Failed', 'warning');
                }
            })
            .catch(function (err) {
                if (window.hideLoading) window.hideLoading('parseBtn', '<i class="fas fa-file-import me-2"></i>Parse Results');
                if (window.showAlert) window.showAlert(window.friendlyError ? window.friendlyError(err, 'server') : 'Failed', 'danger');
            });
    };

    // -------- Biodata codon --------
    // Full beta-actin ORF — table 1, lots of start/stop possibilities
    window.loadCodonExampleComplex = function () {
        var seq = document.getElementById('codonSeq');
        var tab = document.getElementById('codonTableId');
        if (seq) seq.value =
            'ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGC' +
            'TTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCAC' +
            'CAGGGCGTGATGGTGGGCATGGGTCAGAAGGACTCCTATGTGGGTGACGAGGCCCAGAGC' +
            'AAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGAC' +
            'GACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAG' +
            'CACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACC' +
            'CAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTA' +
            'TCCCTGTACGCCTCTGGCCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACC' +
            'CACACTGTGCCCATCTACGAGGGGTATGCCCTCCCCCATGCCATCCTGCGTCTGGACCTG' +
            'GCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTC' +
            'ACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCC' +
            'CTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTAC' +
            'GAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCA' +
            'CTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAAC' +
            'TCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCT' +
            'GGCGGCACCACCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGATCACTGCCCTG' +
            'GCACCCAGCACAATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCCGTGTGG' +
            'ATCGGCGGCTCCATCCTGGCCTCGCTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAG' +
            'GAGTATGACGAGTCCGGCCCCTCCATCGTCCACCGCAAATGCTTCTAG';
        if (tab) tab.value = '1';
    };

    // -------- Biodata protein convert --------
    window.loadProteinExampleComplex = function () {
        var inp = document.getElementById('proteinInput');
        var typ = document.getElementById('proteinConvType');
        if (inp) inp.value = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH';
        if (typ) typ.value = '1to3';
    };

    // -------- Biodata molecular weight --------
    window.loadWeightExampleComplex = function () {
        var s = document.getElementById('weightSeq');
        var t = document.getElementById('weightType');
        var w = document.getElementById('weightCalcType');
        if (s) s.value = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH';
        if (t) t.value = 'protein';
        if (w) w.value = 'average';
    };

    // -------- HMM --------
    window.loadHmmExampleComplex = function () {
        var states = document.getElementById('states');
        var seq = document.getElementById('sequence');
        var type = document.getElementById('hmmType');
        if (states) states.value = '4';
        if (type) type.value = 'sequence';
        if (seq) seq.value = 'ACGTACGTACGTACGTACGTACGTACGTACGT' +
            'AAACCCGGGTTTAAACCCGGGTTTAAACCCGG' +
            'GTATATATACGTACGTACGTACGTATATATAT' +
            'ACGTACGTACGTAAAACCCCGGGGTTTTACGT';
    };

    // -------- Entrez (database) --------
    window.loadExampleSearchComplex = function () {
        var db = document.getElementById('database');
        var term = document.getElementById('searchTerm');
        var max = document.getElementById('retmax');
        if (db) db.value = 'pubmed';
        if (term) term.value = '(CRISPR[Title]) AND (gene therapy[MeSH] OR "cas9"[All Fields]) AND 2023:2024[dp]';
        if (max) max.value = '50';
    };

    // -------- KEGG --------
    window.loadSearchExampleComplex = function () {
        var db = document.getElementById('searchDatabase');
        var q = document.getElementById('searchQuery');
        if (db) db.value = 'pathway';
        if (q) q.value = 'metabolism';
    };

    window.loadListExampleComplex = function () {
        var db = document.getElementById('listDatabase');
        var org = document.getElementById('listOrganism');
        var lim = document.getElementById('listLimit');
        if (db) db.value = 'pathway';
        if (org) org.value = 'hsa';
        if (lim) lim.value = '100';
    };

    // -------- Pathway --------
    // Complex: full 8-reaction glycolysis pathway. Replaces the reactions list
    // (simple loader clears and sets 3 reactions; this replaces with 8).
    window.loadReactionExampleComplex = function () {
        if (typeof window.clearReactions === 'function') window.clearReactions();
        if (!Array.isArray(window.pathwayReactions)) return;
        var rxns = [
            { name: 'Hexokinase', reactants: {'Glucose': -1, 'ATP': -1}, products: {'Glucose-6-P': 1, 'ADP': 1}, catalysts: ['Hexokinase'], reversible: false },
            { name: 'Phosphoglucose isomerase', reactants: {'Glucose-6-P': -1}, products: {'Fructose-6-P': 1}, catalysts: ['PGI'], reversible: true },
            { name: 'Phosphofructokinase', reactants: {'Fructose-6-P': -1, 'ATP': -1}, products: {'Fructose-1,6-BP': 1, 'ADP': 1}, catalysts: ['PFK-1'], reversible: false },
            { name: 'Aldolase', reactants: {'Fructose-1,6-BP': -1}, products: {'G3P': 1, 'DHAP': 1}, catalysts: ['Aldolase'], reversible: true },
            { name: 'Triose-P isomerase', reactants: {'DHAP': -1}, products: {'G3P': 1}, catalysts: ['TPI'], reversible: true },
            { name: 'GAPDH', reactants: {'G3P': -1, 'NAD+': -1, 'Pi': -1}, products: {'1,3-BPG': 1, 'NADH': 1}, catalysts: ['GAPDH'], reversible: true },
            { name: 'Pyruvate kinase', reactants: {'PEP': -1, 'ADP': -1}, products: {'Pyruvate': 1, 'ATP': 1}, catalysts: ['Pyruvate kinase'], reversible: false },
            { name: 'Lactate dehydrogenase', reactants: {'Pyruvate': -1, 'NADH': -1}, products: {'Lactate': 1, 'NAD+': 1}, catalysts: ['LDH'], reversible: true }
        ];
        // Copy stoichiometry into the species field to match simple loader shape
        rxns.forEach(function (r) {
            r.species = Object.assign({}, r.reactants, r.products);
        });
        window.pathwayReactions.length = 0;
        rxns.forEach(function (r) { window.pathwayReactions.push(r); });
        if (typeof window.updateReactionsList === 'function') window.updateReactionsList();
    };

    // -------- Literature (PubMed) --------
    window.loadLiteratureExampleComplex = function () {
        var q = document.getElementById('pubmedQuery');
        var m = document.getElementById('maxResults');
        if (q) q.value = '(BRCA1[Title]) AND (breast cancer[MeSH]) AND (2020:2024[dp])';
        if (m) m.value = '20';
    };

    // -------- Nexus --------
    window.loadNexusExampleComplex = function () {
        var d = document.getElementById('nexusData');
        if (d) d.value =
`#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=8;
    TAXLABELS Human Chimp Gorilla Orangutan Gibbon Macaque Baboon Marmoset;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=30;
    FORMAT DATATYPE=DNA MISSING=? GAP=-;
    MATRIX
        Human     ATCGATCGATCGATCGATCGATCGATCGAT
        Chimp     ATCGATCGATGGATCGATCGATCGATCGAT
        Gorilla   ATCGATCGATCGATCGATCGATCGTTCGAT
        Orangutan ATCGATCGCTCGATCGATCGATCGATCGAT
        Gibbon    ATCGATCGATCGATCGATCGATCCATCGAT
        Macaque   TTCGATCGATCGATCGATCGATCGATCGAT
        Baboon    TTCGATCGATCGATCGATCGATCGATCGAT
        Marmoset  ATCGATCGATAGATCGATCGATCGATCGAT
    ;
END;

BEGIN TREES;
    TREE tree1 = ((Human,Chimp),(Gorilla,Orangutan),((Gibbon,Macaque),(Baboon,Marmoset)));
END;`;
    };

    // -------- SCOP --------
    window.loadScopExampleComplex = function () {
        var id = document.getElementById('scopId');
        var lvl = document.getElementById('scopLevel');
        if (id) id.value = 'd1n0wa_';
        if (lvl) lvl.value = 'superfamily';
    };

    // -------- Codon alignment (advanced HMM page) --------
    window.loadCodonExampleComplexAdv = function () {
        var seqs = document.getElementById('codonSequences');
        if (seqs) seqs.value = [
            '>Human',
            'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC',
            '>Chimp',
            'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC',
            '>Gorilla',
            'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTCACTGCCCTGTGGGGCAAGGTGAAC',
            '>Mouse',
            'ATGGTGCACCTGACTGATGCTGAGAAGTCTGCTGTCTCTTGCCTGTGGGGAAAGGTGAAC',
        ].join('\n');
        var gc = document.getElementById('geneticCode');
        var meth = document.getElementById('alignMethod');
        if (gc) gc.value = '1';
        if (meth) meth.value = 'simple';
    };

})();
