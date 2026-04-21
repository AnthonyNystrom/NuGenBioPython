"""
Inventory of URLs, example buttons, static sample files, and API endpoints
that must keep working through the UI redesign.

This module is imported directly by every test file — the inventory IS the
test oracle, not a human-maintained markdown doc.

Design rules:
- Any page URL that renders HTML today belongs in PAGE_URLS.
- Any file under static/ that a "Load Example" button (or a backend route)
  depends on belongs in STATIC_SAMPLES.
- Any endpoint with a plausible "Load Example + Run" flow belongs in
  API_EXAMPLES with a concrete payload. Endpoints that hit external services
  (NCBI, KEGG, BLAST remote) are listed in EXTERNAL_API_ENDPOINTS and NOT
  exercised in the automated suite — Playwright coverage in phase 0.6
  exercises them with mocks or against live services on a separate schedule.
- LOAD_EXAMPLE_BUTTONS documents every UI example trigger for Playwright
  automation in phase 0.6. It does not drive the API tests directly.
"""

# ---------------------------------------------------------------------------
# Page URLs
# ---------------------------------------------------------------------------
# Every GET URL that renders an HTML page today. Phase 3 hub consolidation
# must keep every one of these returning 200.

PAGE_URLS = [
    "/",
    "/sequence",
    "/features",
    "/seqio",
    "/alignment",
    "/phylo",
    "/structure",
    "/database",
    "/motifs",
    "/restriction",
    "/clustering",
    "/blast",
    "/kegg",
    "/genomediagram",
    "/popgen",
    "/pathway",
    "/unigene",
    "/hmm",
    "/searchio",
    "/swissprot",
    "/biodata",
    # Phase-3 hub URLs — new canonical destinations. The old URLs above
    # resolve to hub views with the correct tab preselected; every old URL
    # must continue to render 200 for backward compatibility.
    "/patterns",
]

# Hub → active_tab mapping. Used by tests to verify that old URLs still
# preselect their correct tab inside the new hub view.
HUB_TAB_PRESELECTION = [
    # (url, hub_name, expected_active_tab)
    ("/patterns",    "patterns", "motifs"),
    ("/motifs",      "patterns", "motifs"),
    ("/restriction", "patterns", "restriction"),
]

# ---------------------------------------------------------------------------
# Static sample / example files
# ---------------------------------------------------------------------------
# Files under static/ that the app or its Load Example buttons depend on.
# Paths are relative to the repo root.

STATIC_SAMPLES = [
    "static/sample_sequences.fasta",
    "static/sample_sequences.tab",
    "static/sample_sequence.gb",
    "static/sample_sequence.embl",
    "static/sample_sequence.swiss",
    "static/sample_alignment.stockholm",
    "static/sample_alignment.clustal",
    "static/sample_alignment.nexus",
    "static/sample_alignment.phylip",
    "static/sample_protein.pdb",
    "static/sample_assembly.ace",
    "static/example_alignment.fasta",
    "static/example_popgen.gen",
]

# Files that Flask auto-generates on first boot (config.py). We assert the
# server creates these, not that they pre-exist in the repo.
STATIC_AUTO_GENERATED = [
    "uploads/example_popgen.gen",
]


# ---------------------------------------------------------------------------
# API endpoints with self-contained example payloads
# ---------------------------------------------------------------------------
# Each entry can be POSTed to the given endpoint with `payload` and must
# return JSON with `success: True` and every key in `expect_keys`.
#
# Keep payloads MINIMAL but representative. The point is not to exhaustively
# test every code path; it's to detect regressions in the Load Example flow.

# Small reusable fixtures
_DNA_SHORT = "ATGGCAGCATTCCAGCATGCAATAACACCGACCGCAATGGCAAATAAAG"
_DNA_MED = (
    "ATGCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCGTAGCTAGCTAGC"
    "ATGCGATCGATCGTAGCATGCATGCATGCGTACGTAGCTAGCTACGTAGC"
)
_PROTEIN_SHORT = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLR"
_FASTA_TWO = ">s1\nATGCATGCATGC\n>s2\nGGGGCCCCGGGG\n"
_FASTQ_ONE = "@r1\nACGTACGT\n+\nIIIIIIII\n"

API_EXAMPLES = [
    # --- Sequence analysis ----------------------------------------------------
    {
        "endpoint": "/api/sequence/analyze",
        "payload": {"sequence": _DNA_SHORT, "type": "dna"},
        "expect_keys": ["analysis"],
        "notes": "mirror of sequence.html Load Example (DNA)",
    },
    {
        "endpoint": "/api/sequence/gc_analysis",
        "payload": {"sequence": _DNA_MED, "window_size": 10},
        "expect_keys": ["gc_content"],
    },
    {
        "endpoint": "/api/sequence/find_orfs",
        "payload": {"sequence": _DNA_MED, "min_length": 30},
        "expect_keys": ["orfs"],
    },
    {
        "endpoint": "/api/sequence/melting_temp",
        "payload": {"sequence": _DNA_SHORT},
        "expect_keys": ["tm_wallace", "tm_gc"],
    },
    {
        "endpoint": "/api/sequence/codon_usage_analysis",
        "payload": {"sequence": _DNA_MED},
        "expect_keys": [],  # route wraps a variable shape; we only assert success
    },
    # --- SeqIO (JSON-only variants; file-upload variants handled in test_api_uploads) --
    {
        "endpoint": "/api/seqio/convert",
        "payload": {
            "sequences": [{"id": "a", "sequence": "ATGCATGC"}],
            "output_format": "fasta",
        },
        "expect_keys": ["converted", "count"],
    },
    {
        "endpoint": "/api/seqio/statistics",
        "payload": {
            "sequences": [
                {"id": "a", "sequence": "ATGCATGC"},
                {"id": "b", "sequence": "GGGGCCCC"},
            ]
        },
        "expect_keys": ["statistics"],
    },
    {
        "endpoint": "/api/seqio/filter",
        "payload": {
            "sequences": [
                {"id": "short", "sequence": "ATG"},
                {"id": "long", "sequence": "ATGCATGCATGCATGC"},
            ],
            "filter_type": "length",
            "min_length": 5,
        },
        "expect_keys": ["filtered", "filtered_count"],
    },
    {
        "endpoint": "/api/seqio/sort",
        "payload": {
            "sequences": [{"id": "b", "sequence": "G"}, {"id": "a", "sequence": "A"}],
            "sort_by": "id",
        },
        "expect_keys": ["sorted"],
    },
    {
        "endpoint": "/api/seqio/slice",
        "payload": {
            "sequences": [{"id": "a", "sequence": "ATGCATGC"}],
            "start": 2,
            "end": 5,
        },
        "expect_keys": ["sliced"],
    },
    # --- Alignment ----------------------------------------------------------
    {
        "endpoint": "/api/alignment/pairwise",
        "payload": {
            "sequence1": "ACGTACGTACGT",
            "sequence2": "ACGTACCGTACGT",
            "match_score": 2,
            "mismatch_score": -1,
            "gap_open": -2,
            "gap_extend": -0.5,
            "mode": "global",
        },
        "expect_keys": ["alignment", "score", "statistics"],
        "notes": "mirror of alignment.html Load Pairwise Example",
    },
    {
        "endpoint": "/api/alignment/multiple",
        "payload": {
            "sequences": [
                "ACGTACGTACGT",
                "ACGTACCGTACGT",
                "ACGTACGTCCGT",
            ]
        },
        "expect_keys": ["alignments"],
    },
    {
        "endpoint": "/api/alignment/matrices",
        "method": "GET",
        "expect_keys": ["matrices"],
    },
    # --- Restriction --------------------------------------------------------
    {
        "endpoint": "/api/restriction/analyze",
        "payload": {
            "sequence": "GAATTCGGGGGGAAGCTT",
            "enzymes": ["EcoRI", "HindIII"],
        },
        "expect_keys": ["analysis"],
    },
    {
        "endpoint": "/api/restriction/list_enzymes",
        "method": "GET",
        "expect_keys": ["enzymes"],
    },
    # --- Features -----------------------------------------------------------
    {
        "endpoint": "/api/features/orf_find",
        "payload": {"sequence": _DNA_MED, "min_length": 30, "strand": "both"},
        "expect_keys": ["orfs"],
    },
    {
        "endpoint": "/api/features/create",
        "payload": {
            "sequence": _DNA_SHORT,
            "feature_type": "gene",
            "start": 0,
            "end": 30,
            "strand": 1,
        },
        "expect_keys": ["feature"],
    },
    {
        "endpoint": "/api/features/extract",
        "payload": {
            "sequence": _DNA_SHORT,
            "start": 0,
            "end": 15,
            "strand": 1,
        },
        "expect_keys": [],
    },
    # --- Motifs -------------------------------------------------------------
    {
        "endpoint": "/api/motifs/create",
        "payload": {
            "sequences": [
                "ACGT",
                "AGGT",
                "ACCT",
                "ATGT",
            ]
        },
        "expect_keys": ["motif"],
    },
    # --- Clustering ---------------------------------------------------------
    {
        "endpoint": "/api/clustering/analyze",
        "payload": {
            "matrix": [[1, 2], [1, 3], [10, 10], [11, 11]],
            "method": "kmeans",
            "n_clusters": 2,
        },
        "expect_keys": ["results"],
    },
    # --- PopGen (pure-python path, no external binary) ----------------------
    {
        "endpoint": "/api/popgen/load_example",
        "payload": {},
        "expect_keys": ["results"],
        "notes": "backend-served example used by the PopGen page",
    },
    # --- Advanced: protparam ------------------------------------------------
    {
        "endpoint": "/api/sequence/protparam",
        "payload": {"sequence": _PROTEIN_SHORT},
        "expect_keys": ["analysis"],
    },
]


# ---------------------------------------------------------------------------
# File-upload API endpoints (multipart)
# ---------------------------------------------------------------------------
# Listed separately because they use a different calling convention.
# Each entry specifies which fixture file (from STATIC_SAMPLES or inline)
# to upload.

API_UPLOAD_EXAMPLES = [
    {
        "endpoint": "/api/seqio/parse",
        "fixture_bytes": _FASTA_TWO.encode(),
        "fixture_name": "example.fasta",
        "form": {"format": "fasta"},
        "expect_keys": ["sequences"],
    },
    {
        "endpoint": "/api/seqio/parse_fastq",
        "fixture_bytes": _FASTQ_ONE.encode(),
        "fixture_name": "example.fastq",
        "form": {"format": "fastq"},
        "expect_keys": ["sequences"],
    },
    {
        "endpoint": "/api/structure/parse",
        "fixture_file": "static/sample_protein.pdb",
        "fixture_name": "sample_protein.pdb",
        "form": {},
        "expect_keys": ["structure_info"],
    },
    {
        "endpoint": "/api/alignment/parse_file",
        "fixture_file": "static/sample_alignment.stockholm",
        "fixture_name": "sample_alignment.stockholm",
        "form": {"format": "stockholm"},
        "expect_keys": ["sequences"],
    },
    {
        "endpoint": "/api/features/parse_genbank",
        "fixture_file": "static/sample_sequence.gb",
        "fixture_name": "sample_sequence.gb",
        "form": {"format": "genbank"},
        "expect_keys": ["features"],
    },
]


# ---------------------------------------------------------------------------
# External API endpoints (skipped in CI; Playwright can exercise with mocks)
# ---------------------------------------------------------------------------

EXTERNAL_API_ENDPOINTS = [
    "/api/blast/search",            # NCBI BLAST
    "/api/database/search",         # NCBI Entrez
    "/api/database/fetch",          # NCBI Entrez
    "/api/database/get_db_info",    # NCBI Entrez
    "/api/kegg/find",               # KEGG REST
    "/api/kegg/get",                # KEGG REST
    "/api/kegg/link",               # KEGG REST
    "/api/kegg/convert",            # KEGG REST
    "/api/pathway/reaction",        # Reactome / KEGG
    "/api/swissprot/fetch",         # UniProt
    "/api/unigene/fetch",           # NCBI UniGene
    "/api/biodata/pubmed",          # NCBI
]


# ---------------------------------------------------------------------------
# Load Example button metadata (for Playwright automation in phase 0.6)
# ---------------------------------------------------------------------------
# NOT consumed by pytest. Documents what each Load Example button does so
# the Playwright harness can click, run, and assert.
#
# Keyed by onclick function name (stable hook) rather than button text (which
# is often just "Example" with an icon and unreliable to match on).
# Extracted from templates via `grep -rho 'onclick="load[A-Za-z]*[Ee]xample'`.

LOAD_EXAMPLE_BUTTONS = [
    # /alignment
    {"page": "/alignment", "onclick_fn": "loadExampleAlignment"},
    {"page": "/alignment", "onclick_fn": "loadExampleMSA"},
    {"page": "/alignment", "onclick_fn": "loadConservationExample"},
    {"page": "/alignment", "onclick_fn": "loadFileIOExample"},
    {"page": "/alignment", "onclick_fn": "loadAllAlignmentsExample"},
    {"page": "/alignment", "onclick_fn": "loadCoordinatesExample"},
    {"page": "/alignment", "onclick_fn": "loadDetailedStatsExample"},
    {"page": "/alignment", "onclick_fn": "loadIdentityMatrixExample"},
    {"page": "/alignment", "onclick_fn": "loadCodonAwareExample"},
    # /biodata
    {"page": "/biodata", "onclick_fn": "loadCodonExample"},
    {"page": "/biodata", "onclick_fn": "loadProteinExample"},
    {"page": "/biodata", "onclick_fn": "loadWeightExample"},
    # /blast
    {"page": "/blast", "onclick_fn": "loadBlastExample"},
    # /clustering
    {"page": "/clustering", "onclick_fn": "loadClusteringExample"},
    # /database
    {"page": "/database", "onclick_fn": "loadExampleSearch"},
    # /features
    {"page": "/features", "onclick_fn": "loadORFExample"},
    {"page": "/features", "onclick_fn": "loadCreateExample"},
    {"page": "/features", "onclick_fn": "loadExtractExample"},
    {"page": "/features", "onclick_fn": "loadCompoundExample"},
    {"page": "/features", "onclick_fn": "loadAnnotateExample"},
    # /genomediagram
    {"page": "/genomediagram", "onclick_fn": "loadBasicExample"},
    {"page": "/genomediagram", "onclick_fn": "loadMultiTrackExample"},
    {"page": "/genomediagram", "onclick_fn": "loadDataTrackExample"},
    {"page": "/genomediagram", "onclick_fn": "loadAdvancedExample"},
    # /hmm (Advanced)
    {"page": "/hmm", "onclick_fn": "loadHmmExample"},
    {"page": "/hmm", "onclick_fn": "loadLiteratureExample"},
    {"page": "/hmm", "onclick_fn": "loadNexusExample"},
    {"page": "/hmm", "onclick_fn": "loadScopExample"},
    {"page": "/hmm", "onclick_fn": "loadCodonExample"},
    # /kegg
    {"page": "/kegg", "onclick_fn": "loadSearchExample"},
    {"page": "/kegg", "onclick_fn": "loadListExample"},
    {"page": "/kegg", "onclick_fn": "loadLinkExample"},
    {"page": "/kegg", "onclick_fn": "loadConvertExample"},
    # /motifs
    {"page": "/motifs", "onclick_fn": "loadMotifExample"},
    # /pathway
    {"page": "/pathway", "onclick_fn": "loadReactionExample"},
    # /phylo
    {"page": "/phylo", "onclick_fn": "loadExampleTree"},
    # /popgen
    {"page": "/popgen", "onclick_fn": "loadPopgenExample"},
    # /restriction
    {"page": "/restriction", "onclick_fn": "loadBasicExample"},
    {"page": "/restriction", "onclick_fn": "loadAdvancedExample"},
    {"page": "/restriction", "onclick_fn": "loadMapExample"},
    {"page": "/restriction", "onclick_fn": "loadCompatibleExample"},
    # /searchio
    {"page": "/searchio", "onclick_fn": "loadParseExample"},
    {"page": "/searchio", "onclick_fn": "loadReadExample"},
    {"page": "/searchio", "onclick_fn": "loadIndexExample"},
    {"page": "/searchio", "onclick_fn": "loadConvertExample"},
    {"page": "/searchio", "onclick_fn": "loadFilterExample"},
    {"page": "/searchio", "onclick_fn": "loadWriteExample"},
    # /sequence
    {"page": "/sequence", "onclick_fn": "loadExample"},
    # /structure
    {"page": "/structure", "onclick_fn": "loadExample"},
    {"page": "/structure", "onclick_fn": "loadSuperimposeExample"},
    {"page": "/structure", "onclick_fn": "loadGeometryExample"},
    {"page": "/structure", "onclick_fn": "loadQualityExample"},
    {"page": "/structure", "onclick_fn": "loadContactsExample"},
    {"page": "/structure", "onclick_fn": "loadInteractionsExample"},
    {"page": "/structure", "onclick_fn": "loadDsspExample"},
    {"page": "/structure", "onclick_fn": "loadRamaExample"},
    {"page": "/structure", "onclick_fn": "loadSasaExample"},
    {"page": "/structure", "onclick_fn": "loadSelectionExample"},
    # /swissprot
    {"page": "/swissprot", "onclick_fn": "loadParseExample"},
    {"page": "/swissprot", "onclick_fn": "loadReadExample"},
    # /unigene
    {"page": "/unigene", "onclick_fn": "loadParseExample"},
    {"page": "/unigene", "onclick_fn": "loadReadExample"},
]
