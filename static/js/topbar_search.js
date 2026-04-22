// Topbar tool search — fuzzy-match across the known tool catalog and
// navigate on click / Enter. Empty/escape closes the dropdown.
(function () {
    'use strict';

    const TOOLS = [
        { label: 'Dashboard',           href: '/',              icon: 'fa-home',            keywords: 'home dashboard overview' },
        { label: 'Sequences hub',       href: '/sequences',     icon: 'fa-dna',             keywords: 'sequence dna rna protein' },
        { label: 'Sequence analysis',   href: '/sequence',      icon: 'fa-dna',             keywords: 'sequence analyze gc content composition' },
        { label: 'SeqIO',               href: '/seqio',         icon: 'fa-file-code',       keywords: 'seqio fasta genbank format parse' },
        { label: 'Features',            href: '/features',      icon: 'fa-puzzle-piece',    keywords: 'orf features genbank annotate extract compound' },
        { label: 'Compare hub',         href: '/compare',       icon: 'fa-align-center',    keywords: 'compare alignment blast searchio' },
        { label: 'Pairwise alignment',  href: '/alignment',     icon: 'fa-align-center',    keywords: 'align pairwise needleman smith waterman' },
        { label: 'BLAST',               href: '/blast',         icon: 'fa-rocket',          keywords: 'blast ncbi search homology' },
        { label: 'SearchIO',            href: '/searchio',      icon: 'fa-search',          keywords: 'searchio parser blast results' },
        { label: 'Patterns hub',        href: '/patterns',      icon: 'fa-shapes',          keywords: 'pattern motif restriction' },
        { label: 'Motif finder',        href: '/motifs',        icon: 'fa-fingerprint',     keywords: 'motif jaspar meme transfac search' },
        { label: 'Restriction',         href: '/restriction',   icon: 'fa-cut',             keywords: 'restriction enzyme digest map cut' },
        { label: 'Structure',           href: '/structure',     icon: 'fa-cube',            keywords: 'structure pdb 3dmol geometry quality dssp ramachandran sasa' },
        { label: 'Phylogeny hub',       href: '/phylogeny',     icon: 'fa-project-diagram', keywords: 'phylogeny tree clustering popgen' },
        { label: 'Phylo trees',         href: '/phylo',         icon: 'fa-tree',            keywords: 'phylo tree newick nexus phyloxml' },
        { label: 'Clustering',          href: '/clustering',    icon: 'fa-circle-nodes',    keywords: 'cluster kmeans hierarchical dbscan pca som' },
        { label: 'Population genetics', href: '/popgen',        icon: 'fa-users',           keywords: 'popgen population hardy weinberg fst allele' },
        { label: 'External data',       href: '/data',          icon: 'fa-database',        keywords: 'database external entrez kegg pathway' },
        { label: 'Entrez / NCBI',       href: '/database',      icon: 'fa-server',          keywords: 'entrez ncbi pubmed nucleotide protein search fetch' },
        { label: 'KEGG',                href: '/kegg',          icon: 'fa-project-diagram', keywords: 'kegg pathway gene reaction enzyme' },
        { label: 'Pathway',             href: '/pathway',       icon: 'fa-network-wired',   keywords: 'pathway reaction sbml metabolic network system' },
        { label: 'SwissProt',           href: '/swissprot',     icon: 'fa-flask',           keywords: 'swissprot uniprot protein record annotation' },
        { label: 'UniGene',             href: '/unigene',       icon: 'fa-sitemap',         keywords: 'unigene gene cluster expression' },
        { label: 'Biodata tables',      href: '/biodata',       icon: 'fa-table',           keywords: 'biodata codon iupac protein molecular weight' },
        { label: 'HMM',                 href: '/hmm',           icon: 'fa-circle-nodes',    keywords: 'hmm hidden markov model literature nexus scop codon' },
        { label: 'Genome diagram',      href: '/genomediagram', icon: 'fa-image',           keywords: 'genome diagram circular linear plot' },
    ];

    const input = document.getElementById('topbarSearchInput');
    const list  = document.getElementById('topbarSearchResults');
    if (!input || !list) return;

    let activeIndex = -1;
    let lastResults = [];

    function score(tool, q) {
        const hay = (tool.label + ' ' + tool.keywords).toLowerCase();
        const labelHay = tool.label.toLowerCase();
        if (labelHay === q) return 100;
        if (labelHay.startsWith(q)) return 80;
        if (hay.includes(' ' + q) || hay.startsWith(q)) return 60;
        if (hay.includes(q)) return 40;
        return 0;
    }

    function render(q) {
        if (!q) { hide(); return; }
        const ranked = TOOLS
            .map(t => ({ tool: t, s: score(t, q) }))
            .filter(r => r.s > 0)
            .sort((a, b) => b.s - a.s)
            .slice(0, 8);
        lastResults = ranked.map(r => r.tool);

        if (lastResults.length === 0) {
            list.innerHTML = '<li class="ts-empty">No matching tools.</li>';
        } else {
            list.innerHTML = lastResults.map((t, i) =>
                `<li role="option" data-href="${t.href}" data-idx="${i}"${i === 0 ? ' aria-selected="true"' : ''}>
                    <i class="fas ${t.icon}" aria-hidden="true"></i>
                    <span>${t.label}</span>
                    <span class="ts-hint">${t.href}</span>
                </li>`
            ).join('');
            activeIndex = 0;
        }
        list.hidden = false;
    }

    function hide() {
        list.hidden = true;
        list.innerHTML = '';
        activeIndex = -1;
        lastResults = [];
    }

    function navigate(i) {
        const t = lastResults[i];
        if (t) window.location.href = t.href;
    }

    function setActive(i) {
        if (!lastResults.length) return;
        activeIndex = (i + lastResults.length) % lastResults.length;
        list.querySelectorAll('li').forEach(function (li, idx) {
            li.setAttribute('aria-selected', idx === activeIndex ? 'true' : 'false');
        });
    }

    input.addEventListener('input', function () { render(input.value.trim().toLowerCase()); });
    input.addEventListener('focus',  function () { if (input.value.trim()) render(input.value.trim().toLowerCase()); });

    input.addEventListener('keydown', function (e) {
        if (list.hidden) return;
        if (e.key === 'ArrowDown') { e.preventDefault(); setActive(activeIndex + 1); }
        else if (e.key === 'ArrowUp') { e.preventDefault(); setActive(activeIndex - 1); }
        else if (e.key === 'Enter')   { e.preventDefault(); navigate(activeIndex); }
        else if (e.key === 'Escape')  { e.preventDefault(); hide(); input.blur(); }
    });

    list.addEventListener('mousedown', function (e) {
        const li = e.target.closest('li[data-href]');
        if (li) { e.preventDefault(); navigate(parseInt(li.dataset.idx, 10)); }
    });

    document.addEventListener('click', function (e) {
        if (e.target !== input && !list.contains(e.target)) hide();
    });

    // Cmd/Ctrl-K to focus
    document.addEventListener('keydown', function (e) {
        if ((e.metaKey || e.ctrlKey) && e.key === 'k') {
            e.preventDefault();
            input.focus();
            input.select();
        }
    });
})();
