// Headless test runner for static/js/formatters.js.
// Invoked by tests/test_formatters.py, which ensures node + jsdom exist.
//
// To run manually:
//   npm install --prefix tests/formatters jsdom
//   node tests/formatters/run.js
'use strict';

const path = require('path');
const fs = require('fs');
const { JSDOM } = require(path.join(__dirname, 'node_modules', 'jsdom'));

// url must be set so jsdom gives a non-opaque origin (required by
// sessionStorage / localStorage access).
const dom = new JSDOM('<!doctype html><html><body></body></html>', {
    url: 'http://localhost/',
});
global.window = dom.window;
global.document = dom.window.document;
global.DOMParser = dom.window.DOMParser;
// workspace.js accesses sessionStorage as a bare global (mirrors browser
// semantics where `window.sessionStorage` and `sessionStorage` are the same).
global.sessionStorage = dom.window.sessionStorage;
global.Event = dom.window.Event;

const REPO_ROOT = path.resolve(__dirname, '..', '..');

// Stub sessionStorage (jsdom provides it, but ensure isolation across runs)
global.window.sessionStorage.clear();

const utilsSrc = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'utils.js'), 'utf8');
eval(utilsSrc.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const RP = global.window.ResultsPanel;

const src = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'formatters.js'), 'utf8');
eval(src.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const F = global.window.NuGenFormatters;

const wsSrc = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'workspace.js'), 'utf8');
eval(wsSrc.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const W = global.window.Workspace;

const rcSrc = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'results_card.js'), 'utf8');
eval(rcSrc.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const RC = global.window.ResultsCard;

let pass = 0, fail = 0;
function t(ok, msg) {
    if (ok) { pass++; console.log('PASS', msg); }
    else { fail++; console.log('FAIL', msg); }
}
function run(name, html, ...needles) {
    const ok = typeof html === 'string' && html.length > 10 &&
        needles.every(function (n) { return html.includes(n); });
    t(ok, name + (ok ? '' : ' (missing: ' + needles.filter(function (n) { return !html.includes(n); }).join(', ') + ')'));
}

// ---------- PubMed ----------
const richPubMed =
'<?xml version="1.0"?>\n' +
'<PubmedArticleSet><PubmedArticle>' +
'  <MedlineCitation>' +
'    <PMID Version="1">35313109</PMID>' +
'    <Article>' +
'      <Journal>' +
'        <JournalIssue>' +
'          <Volume>16</Volume><Issue>4</Issue>' +
'          <PubDate><Year>2022</Year><Month>Apr</Month></PubDate>' +
'        </JournalIssue>' +
'        <Title>Expert review of respiratory medicine</Title>' +
'        <ISOAbbreviation>Expert Rev Respir Med</ISOAbbreviation>' +
'      </Journal>' +
'      <ArticleTitle>CRISPR gene editing - what are the possibilities for respiratory medicine?</ArticleTitle>' +
'      <Pagination><MedlinePgn>371-374</MedlinePgn></Pagination>' +
'      <Abstract>' +
'        <AbstractText Label="INTRODUCTION">CRISPR-Cas9 is revolutionizing gene therapy.</AbstractText>' +
'        <AbstractText Label="AREAS COVERED">This editorial reviews the possibilities for respiratory diseases.</AbstractText>' +
'      </Abstract>' +
'      <AuthorList>' +
'        <Author>' +
'          <LastName>Harrison</LastName><Initials>PT</Initials>' +
'          <AffiliationInfo><Affiliation>University College Cork, Cork, Ireland.</Affiliation></AffiliationInfo>' +
'          <Identifier Source="ORCID">0000-0003-0000-0000</Identifier>' +
'        </Author>' +
'      </AuthorList>' +
'      <Language>eng</Language>' +
'      <PublicationTypeList>' +
'        <PublicationType>Editorial</PublicationType>' +
'        <PublicationType>Review</PublicationType>' +
'      </PublicationTypeList>' +
'      <ChemicalList><Chemical><NameOfSubstance>CRISPR-Cas9 nuclease</NameOfSubstance></Chemical></ChemicalList>' +
'      <MeshHeadingList>' +
'        <MeshHeading><DescriptorName MajorTopicYN="N">CRISPR-Cas Systems</DescriptorName></MeshHeading>' +
'        <MeshHeading><DescriptorName MajorTopicYN="Y">Gene Editing</DescriptorName></MeshHeading>' +
'        <MeshHeading><DescriptorName MajorTopicYN="Y">Genetic Therapy</DescriptorName><QualifierName>methods</QualifierName></MeshHeading>' +
'      </MeshHeadingList>' +
'      <KeywordList><Keyword>CRISPR</Keyword><Keyword>cystic fibrosis</Keyword></KeywordList>' +
'    </Article>' +
'    <MedlineJournalInfo><Country>England</Country><MedlineTA>Expert Rev Respir Med</MedlineTA></MedlineJournalInfo>' +
'  </MedlineCitation>' +
'  <PubmedData>' +
'    <History><PubMedPubDate PubStatus="pubmed"><Year>2022</Year><Month>3</Month><Day>22</Day></PubMedPubDate></History>' +
'    <ArticleIdList>' +
'      <ArticleId IdType="pubmed">35313109</ArticleId>' +
'      <ArticleId IdType="doi">10.1080/17476348.2022.2056021</ArticleId>' +
'    </ArticleIdList>' +
'    <ReferenceList>' +
'      <Reference><Citation>Smith A, et al.</Citation><ArticleIdList><ArticleId IdType="pubmed">12345678</ArticleId></ArticleIdList></Reference>' +
'    </ReferenceList>' +
'  </PubmedData>' +
'</PubmedArticle></PubmedArticleSet>';

const pubHtml = F.formatPubMedXML(richPubMed);
run('PubMed title', pubHtml, 'CRISPR gene editing');
run('PubMed journal+date', pubHtml, 'Expert review of respiratory medicine', '2022');
run('PubMed PMID', pubHtml, '35313109');
run('PubMed DOI scoped (not reference DOI)', pubHtml, '10.1080/17476348.2022.2056021');
run('PubMed author', pubHtml, 'Harrison PT');
run('PubMed ORCID badge', pubHtml, 'ORCID');
run('PubMed affiliation', pubHtml, 'University College Cork');
run('PubMed abstract labelled', pubHtml, 'INTRODUCTION');
run('PubMed pub type', pubHtml, 'Editorial');
run('PubMed keyword', pubHtml, 'cystic fibrosis');
run('PubMed MeSH major marker', pubHtml, '*');
run('PubMed MeSH qualifier', pubHtml, 'methods');
run('PubMed chemical', pubHtml, 'CRISPR-Cas9 nuclease');
run('PubMed language', pubHtml, 'eng');
run('PubMed history', pubHtml, 'pubmed');
run('PubMed reference link', pubHtml, 'pubmed.ncbi.nlm.nih.gov/12345678');

// ---------- GenBank ----------
const gb =
'LOCUS       TEST                 100 bp    DNA     linear   UNK 01-JAN-2025\n' +
'DEFINITION  Test seq.\n' +
'ACCESSION   T001\n' +
'VERSION     T001.1\n' +
'FEATURES             Location/Qualifiers\n' +
'     gene            1..100\n' +
'                     /gene="G"\n' +
'                     /note="full qualifier not truncated"\n' +
'     CDS             1..100\n' +
'                     /protein_id="P1"\n' +
'ORIGIN\n' +
'        1 acgtacgtac gtacgtacgt\n' +
'//';
const gbHtml = F.formatGenBankText(gb);
run('GenBank definition', gbHtml, 'Test seq');
run('GenBank qualifier full', gbHtml, 'full qualifier not truncated');
run('GenBank features table', gbHtml, 'gene', 'CDS');
run('GenBank sequence', gbHtml, 'ACGTACGTAC');

// ---------- FASTA / Newick / XML / Pairwise ----------
run('FASTA', F.formatFasta('>s1\nATGC\n>s2\nGGGG\n'), 's1', 's2');
run('Newick', F.formatNewickOutline('(A,B);'), 'A', 'B');
run('XML pretty', F.formatXMLPretty('<a><b/></a>'), 'a');
run('Pairwise', F.formatPairwiseAlignment('target 0 ACGT 4\n       ||||\nquery  0 ACGT 4'), 'Legend');

// ---------- Graceful failure ----------
run('Empty handled', F.formatPubMedXML(''), 'Empty');
run('Malformed handled', F.formatGenBankText('not genbank'), 'Not a GenBank');

// ---------- ResultsPanel helper ----------
const rpHtml = RP.tabs([
    { id: 'a', title: 'First',  content: '<div>ONE</div>', active: true },
    { id: 'b', title: 'Second', content: '<div>TWO</div>' },
], { prefix: 'rp-test' });
run('ResultsPanel tab titles',   rpHtml, 'First', 'Second');
run('ResultsPanel tab contents', rpHtml, 'ONE', 'TWO');
run('ResultsPanel active class', rpHtml, 'class="nav-link active"');
run('ResultsPanel prefix',       rpHtml, 'id="rp-test-a"', 'data-bs-target="#rp-test-a"');
// XSS-safe title
const rpXss = RP.tabs([{ id: 'x', title: '<script>x</script>', content: '' }]);
t(!rpXss.includes('<script>x</script>'), 'ResultsPanel title XSS-safe');
t(rpXss.includes('&lt;script&gt;'), 'ResultsPanel title HTML-encoded');

// ---------- Workspace helper ----------
W.clear();
t(W.count() === 0, 'Workspace empty after clear()');

// Add / list / get
const entry = W.add({ type: 'sequence', name: 'seq1', data: 'ATGCATGC' });
t(!!entry.id, 'Workspace.add returns entry with id');
t(W.count() === 1, 'Workspace count=1 after add');
t(W.count('sequence') === 1, 'Workspace.count("sequence") filters by type');
t(W.list('sequence')[0].name === 'seq1', 'Workspace.list returns the right item');
t(W.get(entry.id).data === 'ATGCATGC', 'Workspace.get by id');

// Dedup on identical (type + data)
const dup = W.add({ type: 'sequence', name: 'seq1-again', data: 'ATGCATGC' });
t(dup.id === entry.id, 'Workspace.add dedupes on same (type, data)');
t(W.count() === 1, 'Workspace still has 1 item after dup add');

// Different data → distinct entry
W.add({ type: 'sequence', name: 'seq2', data: 'GGGGCCCC' });
t(W.count('sequence') === 2, 'Workspace stores distinct-data items');

// Type filtering
W.add({ type: 'tree', name: 't1', data: '((A,B),(C,D));' });
t(W.count('tree') === 1, 'Workspace stores items of different types');
t(W.list('sequence').length === 2, 'type filter unaffected by other types');

// Remove
W.remove(entry.id);
t(W.count() === 2, 'Workspace.remove reduces count');
t(W.get(entry.id) === null, 'removed item no longer retrievable');

// Persist to sessionStorage
const raw = global.window.sessionStorage.getItem('nugen_workspace_v1');
t(!!raw, 'Workspace persisted to sessionStorage');
const parsed = JSON.parse(raw);
t(parsed.items.length === W.count(), 'persisted count matches in-memory count');

// onChange listener
let events = 0;
const off = W.onChange(function () { events += 1; });
W.add({ type: 'sequence', name: 'seq3', data: 'AAAA' });
t(events >= 1, 'onChange listener fired on add');
off();

// Cap at 50 items
W.clear();
for (let i = 0; i < 60; i++) W.add({ type: 'sequence', name: 's' + i, data: 'X' + i });
t(W.count() === 50, 'Workspace caps at 50 items');

// Clear
W.clear();
t(W.count() === 0, 'Workspace.clear() empties all items');

// ---------- ResultsCard ----------
global.document.body.innerHTML = '<div id="rc-test"></div>';
RC.mount('rc-test', {
    title: 'Demo',
    meta: '<span class="meta">42 items</span>',
    summary: '<div>SUMMARY_CONTENT</div>',
    details: '<div>DETAILS_CONTENT</div>',
    raw: 'RAW_TEXT_BODY',
    downloads: [
        { label: 'A.txt', filename: 'a.txt', text: 'hello' },
        { label: 'B.json', filename: 'b.json', text: '{}', mime: 'application/json' },
    ],
    workspaceItem: { type: 'sequence', name: 'from rc', data: 'ATGC' },
});
const rcHtml = global.document.getElementById('rc-test').innerHTML;
run('ResultsCard renders title',    rcHtml, 'Demo');
run('ResultsCard renders meta',     rcHtml, '42 items');
run('ResultsCard Summary tab',      rcHtml, 'SUMMARY_CONTENT');
run('ResultsCard Details tab',      rcHtml, 'DETAILS_CONTENT');
run('ResultsCard Raw tab (escaped)',rcHtml, 'RAW_TEXT_BODY');
run('ResultsCard Save button',      rcHtml, 'data-rc-save');
run('ResultsCard Copy button',      rcHtml, 'data-rc-copy');
run('ResultsCard Download menu',    rcHtml, 'data-rc-download="0"', 'data-rc-download="1"');

// Save wires to Workspace.add
W.clear();
global.document.querySelector('[data-rc-save]').click();
t(W.count() === 1, 'ResultsCard Save wires to Workspace.add');
t(W.list()[0].name === 'from rc', 'Workspace entry carries the configured name');

// Minimal spec — only title + summary
global.document.body.innerHTML = '<div id="rc-min"></div>';
RC.mount('rc-min', { title: 'Only summary', summary: '<p>just this</p>' });
const rcMin = global.document.getElementById('rc-min').innerHTML;
t(rcMin.includes('just this'), 'ResultsCard works with only Summary');
t(!rcMin.includes('data-rc-save'), 'ResultsCard omits Save when no workspaceItem');
t(!rcMin.includes('data-rc-download'), 'ResultsCard omits Download when no downloads');

// Polish-1: showLoading / showError APIs
global.document.body.innerHTML = '<div id="rc-load"></div>';
RC.showLoading('rc-load', { title: 'Running BLAST…' });
const rcLoad = global.document.getElementById('rc-load').innerHTML;
t(rcLoad.includes('results-card-loading'), 'showLoading marks card with loading class');
t(rcLoad.includes('aria-busy="true"'), 'showLoading sets aria-busy');
t(rcLoad.includes('Running BLAST'), 'showLoading shows custom title');
t(rcLoad.includes('rc-skeleton-row'), 'showLoading renders skeleton rows');

global.document.body.innerHTML = '<div id="rc-err"></div>';
RC.showError('rc-err', { title: 'BLAST failed', message: 'NCBI returned 503' });
const rcErr = global.document.getElementById('rc-err').innerHTML;
t(rcErr.includes('results-card-error'), 'showError marks card with error class');
t(rcErr.includes('role="alert"'), 'showError sets role=alert');
t(rcErr.includes('BLAST failed'), 'showError shows the title');
t(rcErr.includes('NCBI returned 503'), 'showError shows the message');

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
