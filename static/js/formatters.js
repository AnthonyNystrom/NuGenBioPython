// Readable formatters for raw bio-text outputs.
// Exposes window.NuGenFormatters; each formatter returns an HTML string.
(function (global) {
    'use strict';

    function esc(s) {
        return String(s == null ? '' : s)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }

    function empty(msg) {
        return '<div class="alert alert-secondary mb-0">' + esc(msg || 'Nothing to format.') + '</div>';
    }

    // --- PubMed XML ------------------------------------------------------
    function formatPubMedXML(xml) {
        if (!xml) return empty('Empty record.');
        let doc;
        try {
            doc = new DOMParser().parseFromString(xml, 'text/xml');
            if (doc.querySelector('parsererror')) throw new Error('parse error');
        } catch (e) {
            return empty('Could not parse as PubMed XML.');
        }

        const articles = doc.querySelectorAll('PubmedArticle');
        if (articles.length === 0) return empty('No PubmedArticle found.');

        let html = '';
        articles.forEach(function (art) {
            const pmid = text(art, 'PMID');
            const title = text(art, 'ArticleTitle');
            const vernacularTitle = text(art, 'VernacularTitle');

            // Abstract (all labelled sections)
            const abstractParts = art.querySelectorAll('Abstract AbstractText');
            const abstract = Array.from(abstractParts).map(function (n) {
                const label = n.getAttribute('Label');
                return '<div class="mb-1">' + (label ? '<strong>' + esc(label) + ':</strong> ' : '') + esc(n.textContent) + '</div>';
            }).join('');
            const otherAbstract = Array.from(art.querySelectorAll('OtherAbstract AbstractText')).map(function (n) {
                return '<div class="mb-1">' + esc(n.textContent) + '</div>';
            }).join('');

            // Journal
            const journal = text(art, 'Journal Title');
            const isoAbbrev = text(art, 'Journal ISOAbbreviation');
            const year = text(art, 'JournalIssue PubDate Year') || text(art, 'JournalIssue PubDate MedlineDate').slice(0, 4);
            const month = text(art, 'JournalIssue PubDate Month');
            const day = text(art, 'JournalIssue PubDate Day');
            const volume = text(art, 'JournalIssue Volume');
            const issue = text(art, 'JournalIssue Issue');
            const pages = text(art, 'Pagination MedlinePgn');

            // Article IDs — scope to the article's direct PubmedData > ArticleIdList so
            // we don't pick up per-Reference IDs (which also use <ArticleIdList>).
            const articleIdList = art.querySelector(':scope > PubmedData > ArticleIdList');
            const idNodes = articleIdList
                ? Array.from(articleIdList.querySelectorAll(':scope > ArticleId'))
                : [];
            const ids = {};
            idNodes.forEach(function (n) { ids[n.getAttribute('IdType')] = n.textContent; });
            const doiValue = ids.doi || '';
            const pmcValue = ids.pmc || '';

            // Language, country
            const language = Array.from(art.querySelectorAll('Language')).map(n => n.textContent).join(', ');
            const country = text(art, 'MedlineJournalInfo Country');
            const medlineTA = text(art, 'MedlineJournalInfo MedlineTA');

            // Publication types
            const pubTypes = Array.from(art.querySelectorAll('PublicationTypeList PublicationType')).map(function (n) {
                return esc(n.textContent);
            });

            // Authors with affiliations
            const authors = Array.from(art.querySelectorAll('AuthorList Author')).map(function (a) {
                const last = text(a, 'LastName');
                const initials = text(a, 'Initials');
                const collective = text(a, 'CollectiveName');
                const affiliations = Array.from(a.querySelectorAll('AffiliationInfo Affiliation'))
                    .map(n => esc(n.textContent)).filter(Boolean);
                const orcid = Array.from(a.querySelectorAll('Identifier'))
                    .filter(n => n.getAttribute('Source') === 'ORCID')
                    .map(n => esc(n.textContent))[0];
                const name = collective ? esc(collective) : esc(last + (initials ? ' ' + initials : ''));
                return { name: name, affiliations: affiliations, orcid: orcid };
            }).filter(a => a.name);

            // Keywords
            const keywords = Array.from(art.querySelectorAll('KeywordList Keyword')).map(function (n) {
                return esc(n.textContent);
            }).filter(Boolean);

            // MeSH headings
            const meshes = Array.from(art.querySelectorAll('MeshHeadingList MeshHeading')).map(function (mh) {
                const desc = text(mh, 'DescriptorName');
                const majorTopic = mh.querySelector('DescriptorName') && mh.querySelector('DescriptorName').getAttribute('MajorTopicYN') === 'Y';
                const quals = Array.from(mh.querySelectorAll('QualifierName')).map(function (q) {
                    return esc(q.textContent) + (q.getAttribute('MajorTopicYN') === 'Y' ? '*' : '');
                });
                return { name: esc(desc), major: majorTopic, quals: quals };
            });

            // Chemicals
            const chemicals = Array.from(art.querySelectorAll('ChemicalList Chemical NameOfSubstance')).map(function (n) {
                return esc(n.textContent);
            });

            // Grants
            const grants = Array.from(art.querySelectorAll('GrantList Grant')).map(function (g) {
                const gid = text(g, 'GrantID');
                const agency = text(g, 'Agency');
                const gcountry = text(g, 'Country');
                return [gid, agency, gcountry].filter(Boolean).map(esc).join(' — ');
            }).filter(Boolean);

            // Publication dates history
            const history = Array.from(art.querySelectorAll('History PubMedPubDate')).map(function (d) {
                const status = d.getAttribute('PubStatus');
                const y = text(d, 'Year'), m = text(d, 'Month'), dd = text(d, 'Day');
                return esc(status) + ': ' + [y, m, dd].filter(Boolean).map(esc).join('-');
            });

            // References
            const references = Array.from(art.querySelectorAll('ReferenceList Reference')).map(function (r) {
                const cite = text(r, 'Citation');
                const rIds = Array.from(r.querySelectorAll('ArticleId')).map(function (id) {
                    const t = id.getAttribute('IdType');
                    const v = id.textContent;
                    if (t === 'pubmed') return '<a href="https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(v) + '/" target="_blank" rel="noopener">PMID ' + esc(v) + '</a>';
                    if (t === 'pmc') return '<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/' + encodeURIComponent(v) + '/" target="_blank" rel="noopener">' + esc(v) + '</a>';
                    if (t === 'doi') return '<a href="https://doi.org/' + encodeURIComponent(v) + '" target="_blank" rel="noopener">DOI</a>';
                    return esc(t) + ':' + esc(v);
                }).join(' ');
                return esc(cite) + (rIds ? ' <span class="small text-muted">[' + rIds + ']</span>' : '');
            });

            // ---- Render ----
            html += '<div class="card mb-3"><div class="card-body">';
            if (title) html += '<h5 class="card-title mb-2">' + esc(title) + '</h5>';
            if (vernacularTitle && vernacularTitle !== title) {
                html += '<p class="text-muted fst-italic mb-2">' + esc(vernacularTitle) + '</p>';
            }

            const meta = [];
            if (journal || isoAbbrev) {
                let j = '<em>' + esc(journal || isoAbbrev) + '</em>';
                const dateParts = [year, month, day].filter(Boolean).join(' ');
                if (dateParts) j += '. ' + esc(dateParts);
                if (volume) j += ';' + esc(volume);
                if (issue) j += '(' + esc(issue) + ')';
                if (pages) j += ':' + esc(pages);
                meta.push(j);
            }
            if (pmid) meta.push('PMID: <code>' + esc(pmid) + '</code>');
            if (pmcValue) meta.push('PMC: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/' + encodeURIComponent(pmcValue) + '/" target="_blank" rel="noopener">' + esc(pmcValue) + '</a>');
            if (doiValue) meta.push('DOI: <a href="https://doi.org/' + encodeURIComponent(doiValue) + '" target="_blank" rel="noopener">' + esc(doiValue) + '</a>');
            if (meta.length) html += '<p class="text-muted small mb-2">' + meta.join(' &middot; ') + '</p>';

            if (pubTypes.length) {
                html += '<div class="mb-2">' +
                    pubTypes.map(t => '<span class="badge bg-info-subtle text-info-emphasis border border-info-subtle me-1">' + t + '</span>').join('') +
                    '</div>';
            }

            if (authors.length) {
                html += '<p class="mb-2"><strong>Authors:</strong> ' +
                    authors.map(a => a.name + (a.orcid ? ' <span class="badge bg-success-subtle text-success-emphasis border border-success-subtle small">ORCID</span>' : '')).join(', ') +
                    '</p>';
                const withAff = authors.filter(a => a.affiliations.length);
                if (withAff.length) {
                    html += '<details class="mb-2"><summary class="small text-muted">Affiliations (' + withAff.length + ')</summary>' +
                        '<ol class="small mt-1">' +
                        withAff.map(a => '<li><strong>' + a.name + '</strong>: ' + a.affiliations.join('; ') + '</li>').join('') +
                        '</ol></details>';
                }
            }

            if (abstract) {
                html += '<div class="mb-3"><strong>Abstract</strong><div class="small mt-1" style="line-height:1.5">' + abstract + '</div></div>';
            }
            if (otherAbstract) {
                html += '<details class="mb-3"><summary><strong>Other abstract</strong></summary><div class="small mt-1" style="line-height:1.5">' + otherAbstract + '</div></details>';
            }

            if (keywords.length) {
                html += '<div class="mb-2"><strong>Keywords:</strong> ' +
                    keywords.map(k => '<span class="badge bg-light text-dark border me-1 mb-1">' + k + '</span>').join('') +
                    '</div>';
            }

            if (meshes.length) {
                html += '<div class="mb-2"><strong>MeSH headings</strong> <span class="text-muted small">(' + meshes.length + ')</span><div class="mt-1">' +
                    meshes.map(function (m) {
                        const qs = m.quals.length ? ' <span class="text-muted small">/ ' + m.quals.join(' / ') + '</span>' : '';
                        const cls = m.major ? 'bg-primary-subtle text-primary-emphasis border border-primary-subtle' : 'bg-light text-dark border';
                        return '<span class="badge ' + cls + ' me-1 mb-1">' + m.name + (m.major ? '*' : '') + qs + '</span>';
                    }).join('') +
                    '</div></div>';
            }

            if (chemicals.length) {
                html += '<div class="mb-2"><strong>Chemicals:</strong> ' +
                    chemicals.map(c => '<span class="badge bg-light text-dark border me-1 mb-1">' + c + '</span>').join('') +
                    '</div>';
            }

            if (grants.length) {
                html += '<details class="mb-2"><summary><strong>Grants</strong> <span class="text-muted small">(' + grants.length + ')</span></summary>' +
                    '<ul class="small mt-1 mb-0">' + grants.map(g => '<li>' + g + '</li>').join('') + '</ul></details>';
            }

            const sideMeta = [];
            if (language) sideMeta.push('<strong>Language:</strong> ' + esc(language));
            if (country) sideMeta.push('<strong>Country:</strong> ' + esc(country));
            if (medlineTA) sideMeta.push('<strong>Medline TA:</strong> <code>' + esc(medlineTA) + '</code>');
            if (sideMeta.length) html += '<p class="small text-muted mb-2">' + sideMeta.join(' &middot; ') + '</p>';

            if (history.length) {
                html += '<details class="mb-2"><summary class="small text-muted">Publication history</summary>' +
                    '<ul class="small mt-1 mb-0">' + history.map(h => '<li>' + h + '</li>').join('') + '</ul></details>';
            }

            if (references.length) {
                html += '<details class="mb-0"><summary><strong>References</strong> <span class="text-muted small">(' + references.length + ')</span></summary>' +
                    '<ol class="small mt-2 mb-0" style="line-height:1.5">' +
                    references.map(r => '<li class="mb-1">' + r + '</li>').join('') +
                    '</ol></details>';
            }

            html += '</div></div>';
        });
        return html;
    }

    function text(ctx, selector) {
        // simple space-separated selector chain: descendants
        const parts = selector.split(/\s+/);
        let cur = ctx;
        for (let i = 0; i < parts.length; i++) {
            if (!cur) return '';
            cur = cur.querySelector(parts[i]);
        }
        return cur && cur.textContent ? cur.textContent.trim() : '';
    }

    // --- GenBank / GenPept text ------------------------------------------
    function formatGenBankText(gb) {
        if (!gb || !/^LOCUS\b/m.test(gb)) return empty('Not a GenBank/GenPept record.');
        const meta = {};
        const metaFields = ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE', '  ORGANISM'];
        metaFields.forEach(function (field) {
            const trimmed = field.trim();
            const re = new RegExp('^' + field.replace(/\s/g, '\\s') + '\\s+([\\s\\S]*?)(?=\\n[A-Z ]{2,}\\s{2,}|\\nFEATURES|\\nORIGIN|\\n//)', 'm');
            const m = gb.match(re);
            if (m) meta[trimmed] = m[1].replace(/\n\s+/g, ' ').trim();
        });

        // FEATURES block
        const featuresMatch = gb.match(/^FEATURES[^\n]*\n([\s\S]*?)(?=^ORIGIN|^\/\/|^CONTIG)/m);
        const features = [];
        if (featuresMatch) {
            const block = featuresMatch[1];
            const lines = block.split('\n');
            let current = null;
            lines.forEach(function (line) {
                if (/^ {5}\S/.test(line)) {
                    if (current) features.push(current);
                    const parts = line.trim().split(/\s+/);
                    current = { type: parts[0], location: parts.slice(1).join(' '), qualifiers: {} };
                } else if (current && /^ {21}/.test(line)) {
                    const q = line.trim();
                    const m = q.match(/^\/(\w+)=?"?([^"]*)"?$/);
                    if (m) {
                        const k = m[1];
                        const v = m[2] || '';
                        current.qualifiers[k] = current.qualifiers[k]
                            ? current.qualifiers[k] + ' ' + v
                            : v;
                    } else if (current._last) {
                        current.qualifiers[current._last] += ' ' + q.replace(/"$/, '');
                    }
                }
            });
            if (current) features.push(current);
        }

        // ORIGIN sequence
        const originMatch = gb.match(/^ORIGIN[^\n]*\n([\s\S]*?)(?=^\/\/)/m);
        const sequence = originMatch
            ? originMatch[1].replace(/[\s0-9]/g, '').toUpperCase()
            : '';

        let html = '<div class="card mb-3"><div class="card-body">';
        if (meta.DEFINITION) html += '<h5 class="card-title mb-2">' + esc(meta.DEFINITION) + '</h5>';
        const metaRow = [];
        if (meta.LOCUS) {
            const locusParts = meta.LOCUS.split(/\s+/);
            metaRow.push('<strong>Locus:</strong> <code>' + esc(locusParts[0]) + '</code>');
            if (locusParts[1]) metaRow.push('<strong>Length:</strong> ' + esc(locusParts[1]) + ' ' + esc(locusParts[2] || ''));
        }
        if (meta.ACCESSION) metaRow.push('<strong>Accession:</strong> <code>' + esc(meta.ACCESSION.split(/\s+/)[0]) + '</code>');
        if (meta.VERSION) metaRow.push('<strong>Version:</strong> <code>' + esc(meta.VERSION.split(/\s+/)[0]) + '</code>');
        if (meta.ORGANISM) metaRow.push('<strong>Organism:</strong> <em>' + esc(meta.ORGANISM.split('\n')[0]) + '</em>');
        if (meta.SOURCE) metaRow.push('<strong>Source:</strong> ' + esc(meta.SOURCE));
        if (metaRow.length) html += '<p class="small mb-3">' + metaRow.join(' &middot; ') + '</p>';
        if (meta.KEYWORDS && meta.KEYWORDS !== '.') html += '<p class="small text-muted mb-2"><strong>Keywords:</strong> ' + esc(meta.KEYWORDS) + '</p>';

        if (features.length) {
            html += '<details class="mb-2" open><summary><strong>Features</strong> <span class="text-muted small">(' + features.length + ')</span></summary>';
            html += '<div class="table-responsive mt-2"><table class="table table-sm table-striped mb-0"><thead><tr><th style="width:120px">Type</th><th style="width:200px">Location</th><th>Qualifiers</th></tr></thead><tbody>';
            features.forEach(function (f) {
                const q = Object.keys(f.qualifiers).map(function (k) {
                    const v = String(f.qualifiers[k]).replace(/^"|"$/g, '');
                    return '<span class="text-muted">/' + esc(k) + '</span>="' + esc(v) + '"';
                }).join('<br>');
                html += '<tr><td><code>' + esc(f.type) + '</code></td><td><code class="small">' + esc(f.location) + '</code></td><td class="small">' + q + '</td></tr>';
            });
            html += '</tbody></table></div></details>';
        }

        if (sequence) {
            const chunks = sequence.match(/.{1,60}/g) || [];
            html += '<details><summary><strong>Sequence</strong> <span class="text-muted small">(' + sequence.length + ' bp/aa)</span></summary>';
            html += '<pre class="mt-2 p-2 border rounded small" style="max-height:400px; overflow:auto; font-size:11px; background:#f8fafc; line-height:1.4;">' +
                esc(chunks.join('\n')) + '</pre></details>';
        }

        html += '</div></div>';
        return html;
    }

    // --- Generic XML pretty-print ----------------------------------------
    function formatXMLPretty(xml) {
        if (!xml) return empty('Empty XML.');
        const PADDING = '  ';
        const reg = /(>)(<)(\/*)/g;
        const xmlLine = String(xml).replace(reg, '$1\n$2$3').replace(/^\s*\n/gm, '');
        let pad = 0;
        const lines = xmlLine.split('\n');
        const out = [];
        lines.forEach(function (node) {
            let indent = 0;
            if (node.match(/^<\/\w/)) { pad -= 1; }
            else if (node.match(/^<\w[^>]*[^\/]>.*$/) && !node.match(/^<\w[^>]*\/\s*>/)) { indent = 1; }
            out.push(PADDING.repeat(Math.max(0, pad)) + node);
            pad += indent;
        });
        // Color tags
        const colored = out.map(function (line) {
            return esc(line)
                .replace(/&lt;(\/?\w[\w:-]*)/g, '&lt;<span style="color:#0284c7">$1</span>')
                .replace(/(\w+)=&quot;([^&]*)&quot;/g, '<span style="color:#b45309">$1</span>=<span style="color:#15803d">&quot;$2&quot;</span>');
        }).join('\n');
        return '<pre class="p-3 border rounded small" style="max-height:500px; overflow:auto; font-size:11px; background:#f8fafc; line-height:1.4">' + colored + '</pre>';
    }

    // --- FASTA / sequence-text records ----------------------------------
    function formatFasta(text) {
        if (!text) return empty('Empty FASTA.');
        const records = [];
        const lines = String(text).split('\n');
        let cur = null;
        lines.forEach(function (line) {
            if (line.startsWith('>')) {
                if (cur) records.push(cur);
                const header = line.slice(1).trim();
                const sp = header.indexOf(' ');
                cur = {
                    id: sp === -1 ? header : header.slice(0, sp),
                    description: sp === -1 ? '' : header.slice(sp + 1),
                    seq: ''
                };
            } else if (cur) {
                cur.seq += line.trim();
            }
        });
        if (cur) records.push(cur);
        if (!records.length) return empty('No FASTA records found.');

        let html = '<div class="table-responsive"><table class="table table-sm table-striped mb-0"><thead><tr>' +
            '<th>ID</th><th>Description</th><th style="width:80px">Length</th><th>Preview</th></tr></thead><tbody>';
        records.forEach(function (r) {
            html += '<tr>' +
                '<td><code>' + esc(r.id) + '</code></td>' +
                '<td class="small">' + esc(r.description) + '</td>' +
                '<td><span class="badge bg-light text-dark border">' + r.seq.length + '</span></td>' +
                '<td class="font-monospace small text-muted">' + esc(r.seq.slice(0, 60)) + (r.seq.length > 60 ? '…' : '') + '</td>' +
                '</tr>';
        });
        html += '</tbody></table></div>';
        html += '<p class="text-muted small mt-2 mb-0">' + records.length + ' record' + (records.length === 1 ? '' : 's') + '</p>';
        return html;
    }

    // --- Newick outline --------------------------------------------------
    function formatNewickOutline(nwk) {
        if (!nwk) return empty('Empty tree.');
        const s = String(nwk).trim().replace(/;$/, '');
        let depth = 0;
        let out = '';
        let token = '';
        function flushToken() {
            if (token.trim()) {
                const parts = token.split(':');
                const name = parts[0].trim() || '<em class="text-muted">(internal)</em>';
                const bl = parts[1] ? ' <span class="text-muted small">:' + esc(parts[1]) + '</span>' : '';
                out += '  '.repeat(depth) + '• ' + esc(parts[0].trim() || '') + bl + '\n';
            }
            token = '';
        }
        for (let i = 0; i < s.length; i++) {
            const ch = s[i];
            if (ch === '(') {
                flushToken();
                out += '  '.repeat(depth) + '┐\n';
                depth++;
            } else if (ch === ')') {
                flushToken();
                depth = Math.max(0, depth - 1);
            } else if (ch === ',') {
                flushToken();
            } else {
                token += ch;
            }
        }
        flushToken();
        return '<pre class="p-3 border rounded small" style="max-height:500px; overflow:auto; background:#f8fafc;">' + out + '</pre>';
    }

    // --- Pairwise alignment (BioPython str format) -----------------------
    function formatPairwiseAlignment(txt) {
        if (!txt) return empty('Empty alignment.');
        const lines = String(txt).split('\n').filter(function (l) { return l.length > 0; });
        // BioPython emits repeating 3-line blocks: target, middle (match indicators), query
        // Each line looks like: "target  0 ACGT-CG 7"
        // We parse positional sequence chunks, color them, then reassemble.
        const target = [];
        const middle = [];
        const query = [];
        for (let i = 0; i + 2 < lines.length; i += 3) {
            target.push(extractSeq(lines[i]));
            middle.push(extractMiddle(lines[i + 1]));
            query.push(extractSeq(lines[i + 2]));
        }
        function extractSeq(line) {
            const parts = line.split(/\s+/);
            if (parts.length >= 3) return parts[parts.length - 2] || '';
            return line;
        }
        function extractMiddle(line) {
            // middle line has no leading label; just the indicator string between the number columns
            return line.replace(/^\s+/, '');
        }

        const joinedT = target.join('');
        const joinedQ = query.join('');
        const width = 60;
        let html = '<div class="mb-2 small text-muted">Legend: ' +
            '<span style="background:#d1fae5;padding:0 4px;border-radius:2px">match</span> ' +
            '<span style="background:#fee2e2;padding:0 4px;border-radius:2px">mismatch</span> ' +
            '<span style="background:#e5e7eb;padding:0 4px;border-radius:2px">gap</span></div>';
        html += '<pre class="p-3 border rounded small font-monospace" style="max-height:500px; overflow:auto; background:#ffffff; line-height:1.4">';
        for (let pos = 0; pos < joinedT.length; pos += width) {
            const tChunk = joinedT.slice(pos, pos + width);
            const qChunk = joinedQ.slice(pos, pos + width);
            let tRow = '';
            let qRow = '';
            for (let i = 0; i < tChunk.length; i++) {
                const tc = tChunk[i];
                const qc = qChunk[i] || ' ';
                let bg = '';
                if (tc === '-' || qc === '-') bg = '#e5e7eb';
                else if (tc === qc) bg = '#d1fae5';
                else bg = '#fee2e2';
                tRow += '<span style="background:' + bg + '">' + esc(tc) + '</span>';
                qRow += '<span style="background:' + bg + '">' + esc(qc) + '</span>';
            }
            const label = String(pos + 1).padStart(6, ' ');
            html += esc(label) + ' ' + tRow + '\n' + '       ' + qRow + '\n\n';
        }
        html += '</pre>';
        return html;
    }

    global.NuGenFormatters = {
        escapeHtml: esc,
        formatPubMedXML: formatPubMedXML,
        formatGenBankText: formatGenBankText,
        formatXMLPretty: formatXMLPretty,
        formatFasta: formatFasta,
        formatNewickOutline: formatNewickOutline,
        formatPairwiseAlignment: formatPairwiseAlignment
    };
})(window);
