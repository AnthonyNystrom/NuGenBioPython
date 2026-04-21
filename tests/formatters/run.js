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

const dom = new JSDOM('<!doctype html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.DOMParser = dom.window.DOMParser;

const REPO_ROOT = path.resolve(__dirname, '..', '..');
const utilsSrc = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'utils.js'), 'utf8');
eval(utilsSrc.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const RP = global.window.ResultsPanel;

const src = fs.readFileSync(path.join(REPO_ROOT, 'static', 'js', 'formatters.js'), 'utf8');
eval(src.replace(/\}\)\(window\);\s*$/, '})(global.window);'));
const F = global.window.NuGenFormatters;

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

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
