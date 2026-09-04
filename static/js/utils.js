// Shared UI helpers for NuGenBioPython frontend.
// Loaded from base.html before all other page scripts, so these globals are
// available on every page.
(function (global) {
    'use strict';

    // D13: buttonId → results container id convention. Every tool that
    // follows the pattern `<button id="xxxBtn">` with results rendered
    // into `#xxxResults` (or a mapped container) gets a loading skeleton
    // automatically. Tools can opt out by not registering their button.
    var LOADING_TARGETS = {
        alignBtn: 'alignmentResults',
        analyzeBtn: 'resultsArea',
        blastBtn: 'blastResults',
        basicAnalyzeBtn: 'basicResults',
        advancedAnalyzeBtn: 'advancedResults',
        generateMapBtn: 'mapVisualization',
        clusterBtn: 'clusterAssignments',
        parsePopgenBtn: 'popgenResults',
        buildTreeBtn: 'treeInfo',
        parseTreeBtn: 'treeInfo',
        searchBtn: 'searchResults',
        listBtn: 'listResults',
        linkBtn: 'linkResults',
        convertBtn: 'convertResults',
        infoBtn: 'infoResults',
        parseBtn: 'parseResults',
        readBtn: 'readResults',
        indexBtn: 'indexResults',
        filterBtn: 'filterResults',
        writeBtn: 'writeResults',
        codonBtn: 'codonResults',
        orfBtn: 'orfResults',
        extractBtn: 'extractResults',
        compoundBtn: 'compoundResults',
        annotateSubmitBtn: 'annotateResults',
        createMotifBtn: 'createLogoView',
        searchMotifBtn: 'searchResults',
        superimposeBtn: 'superimposeResults',
        geometryBtn: 'geometryResults',
        qualityBtn: 'qualityResults',
        contactsBtn: 'contactsResults',
        interactionsBtn: 'interactionsResults',
        dsspBtn: 'dsspResults',
        ramaBtn: 'ramaResults',
        sasaBtn: 'sasaResults',
        selectionBtn: 'selectionResults',
        analyzeSystemBtn: 'systemAnalysisResults',
        analyzeNetworkBtn: 'networkAnalysisResults',
        generateVizBtn: 'pathwayVisualization',
        buildHmmBtn: 'hmmAnalysis',
        searchLiteratureBtn: 'literatureResults',
        parseNexusBtn: 'nexusResults',
        lookupScopBtn: 'scopResults',
        alignCodonBtn: 'codonResults',
        // Expanded coverage (B3): every fetching button maps to its result
        // container so page-load skeletons + in-flight auto-timeouts fire.
        advProtParamBtn: 'advProtResults',
        batchBtn: 'batchResults',
        calculateStatsBtn: 'statsResults',
        checksumBtn: 'checksumResults',
        consensusBtn: 'consensusResults',
        createAdvancedBtn: 'advancedDiagramDisplay',
        createBasicDiagramBtn: 'basicDiagramDisplay',
        createBtn: 'createResults',
        createDataTracksBtn: 'dataTracksDisplay',
        createMultiTrackBtn: 'multiTrackDisplay',
        decodeBtn: 'decodingResults',
        fastqBtn: 'fastqResults',
        fetchBtn: 'fetchResults',
        formatsBtn: 'formatsResults',
        globalQueryBtn: 'globalResults',
        iupacBtn: 'iupacResults',
        meltingTempBtn: 'meltingTempResults',
        molWeightBtn: 'molWeightResults',
        msaAlignBtn: 'msaResults',
        parseStructureBtn: 'structureInfo',
        protParamBtn: 'protParamResults',
        proteinBtn: 'proteinResults',
        refDataBtn: 'refDataResults',
        sixFrameBtn: 'sixFrameResults',
        sliceBtn: 'sliceResults',
        sortBtn: 'sortResults',
        trainHmmBtn: 'trainingResults',
        trimBtn: 'trimResults',
        weightBtn: 'weightResults',
    };

    // How long to wait before a loading state auto-degrades to an error.
    // Covers external-service hangs (KEGG / Entrez / BLAST that stall).
    var LOADING_TIMEOUT_MS = 60000;
    // Some buttons kick off requests that legitimately take longer.
    // BLAST queries against NCBI commonly take 30-90s.
    var PER_BTN_TIMEOUT_MS = {
        blastBtn: 180000,          // 3 min — BLAST against nt/nr is slow
        generateVizBtn: 90000,     // pathway graphviz render
    };
    var _loadingTimers = {};

    function showLoading(btnId, label) {
        var btn = document.getElementById(btnId);
        if (!btn) return;
        if (!btn.dataset.originalHtml) {
            btn.dataset.originalHtml = btn.innerHTML;
        }
        btn.disabled = true;
        var text = label || 'Processing...';
        btn.innerHTML = '<span class="spinner-border spinner-border-sm me-2"></span>' + escapeHtml(text);

        var targetId = LOADING_TARGETS[btnId];
        if (targetId && global.ResultsCard && global.ResultsCard.showLoading) {
            global.ResultsCard.showLoading(targetId, { title: text });
        }

        // Auto-degrade to error if the request never resolves.
        var timeout = PER_BTN_TIMEOUT_MS[btnId] || LOADING_TIMEOUT_MS;
        if (_loadingTimers[btnId]) clearTimeout(_loadingTimers[btnId]);
        _loadingTimers[btnId] = setTimeout(function () {
            if (!btn.disabled) return;  // already completed
            hideLoading(btnId);
            var secs = Math.round(timeout / 1000);
            var msg = 'The request took longer than ' + secs +
                      's — likely a slow external service. Try again in a moment.';
            if (targetId && global.ResultsCard && global.ResultsCard.showError) {
                global.ResultsCard.showError(targetId, { title: 'Request timed out', message: msg });
            }
            showAlert(msg, 'warning');
        }, timeout);
    }

    function hideLoading(btnId, originalHtml) {
        var btn = document.getElementById(btnId);
        if (!btn) return;
        btn.disabled = false;
        if (_loadingTimers[btnId]) {
            clearTimeout(_loadingTimers[btnId]);
            delete _loadingTimers[btnId];
        }
        var restore = originalHtml || btn.dataset.originalHtml;
        if (restore !== undefined) {
            btn.innerHTML = restore;
        }
        var targetId = LOADING_TARGETS[btnId];
        if (targetId) {
            var container = document.getElementById(targetId);
            if (container && container.querySelector('.results-card-loading')) {
                container.innerHTML = '';
            }
        }
    }

    function showAlert(message, type) {
        type = type || 'info';

        var alertDiv = document.createElement('div');
        alertDiv.className = 'alert alert-' + type + ' alert-dismissible fade show';
        alertDiv.setAttribute('role', 'alert');

        var msgNode = document.createElement('span');
        msgNode.textContent = message == null ? '' : String(message);
        alertDiv.appendChild(msgNode);

        var closeBtn = document.createElement('button');
        closeBtn.type = 'button';
        closeBtn.className = 'btn-close';
        closeBtn.setAttribute('data-bs-dismiss', 'alert');
        closeBtn.setAttribute('aria-label', 'Close');
        alertDiv.appendChild(closeBtn);

        // D4: dedicated toast target, no longer depends on Bootstrap card
        // chrome. Falls through to <main> and body.
        var container =
            document.getElementById('app-alerts') ||
            document.getElementById('main-content') ||
            document.querySelector('.main-content') ||
            document.body;
        container.insertBefore(alertDiv, container.firstChild);

        setTimeout(function () {
            if (alertDiv.parentNode) alertDiv.remove();
        }, 5000);
    }

    // fetch wrapper that checks response.ok and surfaces server error payload.
    // Returns parsed JSON, or throws Error with a useful message.
    function fetchJSON(url, opts) {
        return fetch(url, opts).then(function (res) {
            var ct = res.headers.get('content-type') || '';
            var isJson = ct.indexOf('application/json') !== -1;
            var bodyPromise = isJson ? res.json().catch(function () { return null; })
                                     : res.text().then(function (t) { return { _text: t }; });
            return bodyPromise.then(function (body) {
                if (!res.ok) {
                    var msg = 'HTTP ' + res.status;
                    if (body && typeof body === 'object') {
                        if (body.error) msg = body.error;
                        else if (body._text) msg = body._text.slice(0, 200);
                    }
                    var err = new Error(msg);
                    err.status = res.status;
                    err.body = body;
                    throw err;
                }
                return body;
            });
        });
    }

    function escapeHtml(s) {
        return String(s)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }

    // ---- ResultsPanel --------------------------------------------------
    // Produce consistent tabbed result scaffolding so the 3+ pages that
    // already hand-rolled "Formatted / Raw" tabs (database modal,
    // alignment result, seqio convert modal, phylo Newick export) can share
    // one implementation. Each tab's `content` is trusted HTML — callers
    // pass already-sanitized markup from formatters.js or escapeHtml wrap.
    //
    // Usage:
    //   ResultsPanel.tabs([
    //     { id: "fmt", title: "Formatted", content: "<div>…</div>", active: true },
    //     { id: "raw", title: "Raw",       content: "<pre>…</pre>" },
    //   ])
    //   → returns an HTML string ready to drop into innerHTML.
    //
    // Every tab id is namespaced with a unique prefix so multiple panels
    // on the same page don't collide.
    let _panelSeq = 0;
    function buildResultsPanelTabs(tabs, opts) {
        opts = opts || {};
        _panelSeq += 1;
        const prefix = opts.prefix || ('rp' + _panelSeq);

        // Normalize: ensure exactly one tab is active (first by default)
        let hasActive = false;
        const normalized = tabs.map(function (t, i) {
            const active = !!t.active;
            if (active) hasActive = true;
            return {
                id: prefix + '-' + (t.id || ('tab' + i)),
                title: t.title,
                content: t.content == null ? '' : t.content,
                active: active,
            };
        });
        if (!hasActive && normalized.length) normalized[0].active = true;

        let nav = '<ul class="nav nav-tabs mb-2" role="tablist">';
        let panes = '<div class="tab-content">';
        normalized.forEach(function (t) {
            const actCls = t.active ? ' active' : '';
            const showCls = t.active ? ' show active' : '';
            nav += '<li class="nav-item" role="presentation">' +
                '<button class="nav-link' + actCls + '" type="button" role="tab"' +
                ' data-bs-toggle="tab" data-bs-target="#' + t.id + '">' +
                escapeHtml(t.title) + '</button></li>';
            panes += '<div class="tab-pane fade' + showCls + '" id="' + t.id +
                '" role="tabpanel">' + t.content + '</div>';
        });
        nav += '</ul>';
        panes += '</div>';
        return nav + panes;
    }

    // Helper: render directly into a container element by id
    function renderResultsPanel(containerId, tabs, opts) {
        const el = typeof containerId === 'string'
            ? document.getElementById(containerId)
            : containerId;
        if (!el) return;
        el.innerHTML = buildResultsPanelTabs(tabs, opts);
    }

    const ResultsPanel = {
        tabs: buildResultsPanelTabs,
        render: renderResultsPanel,
    };

    // ---- Global keyboard shortcuts -------------------------------------
    // Phase-5 polish: cheap UX wins that don't collide with typing.
    //
    //   /        → focus the first visible textarea/input (scroll into view)
    //   g then h → go to Dashboard (vim-style chord)
    //   Shift-? → open a small help overlay listing the bindings
    //   Esc     → close the overlay (also handled elsewhere for sidebar/panel)
    //
    // All shortcuts are skipped while an editable field has focus, so they
    // never steal keystrokes from the user's actual typing.
    function isEditable(el) {
        if (!el) return false;
        const tag = el.tagName;
        return tag === 'TEXTAREA' || tag === 'INPUT' || tag === 'SELECT' ||
               el.isContentEditable;
    }
    function focusFirstInput() {
        const candidates = Array.from(document.querySelectorAll(
            'textarea:not([readonly]), input[type="text"], input[type="search"], input:not([type])'
        ));
        const visible = candidates.find(function (el) {
            return el.offsetParent !== null && !el.disabled;
        });
        if (visible) {
            visible.focus({ preventScroll: false });
            visible.scrollIntoView({ block: 'center', behavior: 'smooth' });
        }
    }
    function toggleShortcutsHelp() {
        let overlay = document.getElementById('ng-shortcuts-help');
        if (overlay) {
            overlay.remove();
            // Return focus to the ? button in the topbar
            var shortcutsTrigger = document.querySelector('[data-action="toggleShortcutsHelp"]');
            if (shortcutsTrigger) shortcutsTrigger.focus();
            return;
        }
        overlay = document.createElement('div');
        overlay.id = 'ng-shortcuts-help';
        overlay.className = 'ng-shortcuts-overlay';
        overlay.innerHTML =
            '<div class="ng-shortcuts-card" role="dialog" aria-label="Keyboard shortcuts">' +
                '<div class="ng-shortcuts-title">Keyboard shortcuts</div>' +
                '<dl>' +
                    '<dt><kbd>⌘K</kbd> / <kbd>Ctrl+K</kbd></dt><dd>Focus the topbar search</dd>' +
                    '<dt><kbd>/</kbd></dt><dd>Focus the first input on the page</dd>' +
                    '<dt><kbd>⌘Enter</kbd></dt><dd>Submit the active form</dd>' +
                    '<dt><kbd>⌘S</kbd></dt><dd>Save the most-recent result to Workspace</dd>' +
                    '<dt><kbd>g</kbd> <kbd>h</kbd></dt><dd>Go to Dashboard</dd>' +
                    '<dt><kbd>g</kbd> <kbd>s</kbd></dt><dd>Go to Sequences</dd>' +
                    '<dt><kbd>g</kbd> <kbd>b</kbd></dt><dd>Go to BLAST</dd>' +
                    '<dt><kbd>g</kbd> <kbd>a</kbd></dt><dd>Go to Alignment</dd>' +
                    '<dt><kbd>g</kbd> <kbd>t</kbd></dt><dd>Go to Structure</dd>' +
                    '<dt><kbd>g</kbd> <kbd>w</kbd></dt><dd>Toggle Workspace drawer</dd>' +
                    '<dt><kbd>?</kbd></dt><dd>Show / hide this help</dd>' +
                    '<dt><kbd>Esc</kbd></dt><dd>Close panels and overlays</dd>' +
                '</dl>' +
                '<button type="button" class="btn-app-sm btn-app-secondary">Close</button>' +
            '</div>';
        overlay.addEventListener('click', function (e) {
            if (e.target === overlay) overlay.remove();
        });
        overlay.querySelector('button').addEventListener('click', function () { overlay.remove(); });
        document.body.appendChild(overlay);

        // D25: focus-trap within overlay — Tab stays inside, Esc closes.
        var firstFocusable = overlay.querySelector('button');
        if (firstFocusable) firstFocusable.focus();
        overlay.addEventListener('keydown', function (e) {
            if (e.key !== 'Tab') return;
            var focusables = overlay.querySelectorAll('button, [href], [tabindex]:not([tabindex="-1"])');
            if (!focusables.length) return;
            var first = focusables[0];
            var last = focusables[focusables.length - 1];
            if (e.shiftKey && document.activeElement === first) {
                e.preventDefault();
                last.focus();
            } else if (!e.shiftKey && document.activeElement === last) {
                e.preventDefault();
                first.focus();
            }
        });
    }

    // ---- data-action event delegation ---------------------------------
    // Phase-5 polish: lets us remove inline onclick= handlers page-wide so
    // we can tighten CSP later. A `<button data-action="foo">` is
    // equivalent to `<button data-action="foo">` — the global `foo` function
    // still lives on window. If it takes no arguments, use data-action
    // alone; for arguments pass data-action-args (JSON).
    //
    //     <button data-action="loadExample">Example</button>
    //     <button data-action="showDetails" data-action-args='[123,"x"]'>…</button>
    document.addEventListener('click', function (e) {
        const btn = e.target.closest('[data-action]');
        if (!btn) return;
        const name = btn.dataset.action;
        const fn = window[name];
        if (typeof fn !== 'function') return;
        let args = [];
        if (btn.dataset.actionArgs) {
            try { args = JSON.parse(btn.dataset.actionArgs); } catch (_) {}
            if (!Array.isArray(args)) args = [args];
        }
        fn.apply(btn, args);
    });

    // ---- non-click events ----------------------------------------------
    // onchange / onsubmit / oninput are inline event handlers too, so the
    // nonce-based CSP blocks them exactly like onclick did. These mirror the
    // click dispatcher:
    //
    //     <select data-action-change="loadEnzymeBrowser">
    //     <input  data-action-change="loadSequenceFile" data-action-args="[1]">
    //     <form   data-action-submit="buildHMM">
    //
    // The DOM event is appended as the final argument, so a handler written as
    // `fn(event)` keeps working and one written as `fn(1)` simply ignores it.
    // Submit handlers get preventDefault() called for them, which is what
    // every migrated onsubmit did by hand.
    function delegate(eventName, attribute, preventDefault) {
        document.addEventListener(eventName, function (e) {
            var el = e.target && e.target.closest
                ? e.target.closest('[' + attribute + ']')
                : null;
            if (!el) return;
            var name = el.getAttribute(attribute);
            var fn = global[name];
            if (typeof fn !== 'function') return;
            if (preventDefault) e.preventDefault();
            var args = [];
            if (el.dataset.actionArgs) {
                try { args = JSON.parse(el.dataset.actionArgs); } catch (_) {}
                if (!Array.isArray(args)) args = [args];
            }
            fn.apply(el, args.concat([e]));
        }, eventName === 'submit');
    }

    delegate('change', 'data-action-change', false);
    delegate('input', 'data-action-input', false);
    delegate('submit', 'data-action-submit', true);

    // Images that should quietly disappear when they fail to load — the old
    // inline onerror="this.style.display='none'". Error events on <img> do
    // not bubble, so this listens in the capture phase.
    document.addEventListener('error', function (e) {
        var el = e.target;
        if (el && el.tagName === 'IMG' && el.hasAttribute('data-hide-on-error')) {
            el.style.display = 'none';
        }
    }, true);



    // D15: real keyboard shortcuts catalog. Cmd/Ctrl-Enter submits,
    // Cmd/Ctrl-S saves latest result, g-chords navigate.
    let chord = null;
    let chordTimer = null;
    var CHORD_NAV = {
        h: '/',
        s: '/sequence',
        a: '/alignment',
        b: '/blast',
        m: '/motifs',
        r: '/restriction',
        t: '/structure',
        p: '/phylogeny',
        d: '/data',
        k: '/kegg',
        w: null,  // special: toggle workspace drawer
    };

    document.addEventListener('keydown', function (e) {
        if (e.key === 'Escape') {
            const overlay = document.getElementById('ng-shortcuts-help');
            if (overlay) overlay.remove();
            return;
        }

        // Cmd/Ctrl-Enter: submit active form (works from inside textarea/input)
        if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
            var focused = document.activeElement;
            var form = focused && focused.form;
            if (!form) {
                form = document.querySelector('.tool-input-pane form, .tool-page form, form[id]');
            }
            if (form) {
                e.preventDefault();
                var submit = form.querySelector('[type="submit"]');
                if (submit) submit.click();
                else form.dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
            }
            return;
        }

        // Cmd/Ctrl-S: save latest result to Workspace (first Save button visible)
        if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 's') {
            var saveBtn = document.querySelector('.results-card [data-rc-save]');
            if (saveBtn) {
                e.preventDefault();
                saveBtn.click();
            }
            return;
        }

        if (isEditable(document.activeElement)) return;
        if (e.ctrlKey || e.metaKey || e.altKey) return;

        if (e.key === '/') {
            e.preventDefault();
            focusFirstInput();
            return;
        }
        if (e.key === '?') {
            e.preventDefault();
            toggleShortcutsHelp();
            return;
        }
        // Chord handling: g then [target]
        if (chord === 'g') {
            clearTimeout(chordTimer);
            chord = null;
            var dest = CHORD_NAV[e.key.toLowerCase()];
            if (e.key.toLowerCase() === 'w') {
                var wsBtn = document.getElementById('topbarWorkspaceBtn');
                if (wsBtn) { e.preventDefault(); wsBtn.click(); }
                return;
            }
            if (dest) {
                e.preventDefault();
                window.location.href = dest;
                return;
            }
        } else if (e.key === 'g') {
            chord = 'g';
            chordTimer = setTimeout(function () { chord = null; }, 1200);
        }
    });

    // Per-URL timeout overrides — endpoints known to run long against external
    // services. Must be >= the matching PER_BTN_TIMEOUT_MS entry so the button
    // skeleton and the fetch abort are aligned.
    var URL_TIMEOUT_MS = [
        { match: '/api/blast/',            timeout: 180000 }, // NCBI BLAST, 3 min
        { match: '/api/pathway/visualize', timeout: 90000  }, // graphviz render
    ];
    function resolveTimeoutForUrl(resource) {
        var url = typeof resource === 'string' ? resource : (resource && resource.url) || '';
        for (var i = 0; i < URL_TIMEOUT_MS.length; i++) {
            if (url.indexOf(URL_TIMEOUT_MS[i].match) >= 0) return URL_TIMEOUT_MS[i].timeout;
        }
        return LOADING_TIMEOUT_MS;
    }

    // Global fetch timeout. Wraps window.fetch so any request that stalls
    // against an external service (NCBI, KEGG, UniProt) aborts and surfaces a
    // real error instead of hanging the UI forever. Endpoints matched by
    // URL_TIMEOUT_MS get a longer ceiling (BLAST: 180s, pathway viz: 90s).
    // Guarded so test environments without fetch (jsdom) don't crash.
    if (typeof global.fetch === 'function' && typeof AbortController === 'function' && !global.__fetchPatched) {
        global.__fetchPatched = true;
        var origFetch = global.fetch.bind(global);
        global.fetch = function (resource, opts) {
            opts = opts || {};
            if (opts.signal) return origFetch(resource, opts);
            var controller = new AbortController();
            var timeoutMs = resolveTimeoutForUrl(resource);
            var to = setTimeout(function () { controller.abort(); }, timeoutMs);
            opts.signal = controller.signal;
            return origFetch(resource, opts).finally(function () { clearTimeout(to); });
        };
    }

    // B5: keep `aria-selected` in sync with Bootstrap tab state across the app,
    // so every `role="tab"` button has the correct state without per-template
    // markup. Covers both initial render and click/keyboard-driven changes.
    function syncTabAriaSelected(root) {
        var scope = root || document;
        scope.querySelectorAll('[role="tab"]').forEach(function (tab) {
            var isActive = tab.classList.contains('active');
            tab.setAttribute('aria-selected', isActive ? 'true' : 'false');
        });
    }
    document.addEventListener('DOMContentLoaded', function () { syncTabAriaSelected(); });
    // Bootstrap fires shown.bs.tab on activation; we listen at document level.
    document.addEventListener('shown.bs.tab', function () { syncTabAriaSelected(); });
    // Fallback for clicks on tabs that don't use Bootstrap's API.
    document.addEventListener('click', function (e) {
        var tab = e.target.closest('[role="tab"]');
        if (!tab) return;
        // Let Bootstrap flip .active first, then resync on the next tick.
        setTimeout(syncTabAriaSelected, 0);
    });

    // P4: map raw errors from external services to user-friendly messages.
    // Takes an error object/string + a service name; returns a human-readable
    // sentence suitable for ResultsCard.showError() / showAlert().
    function friendlyError(error, service) {
        var svc = service || 'the service';
        var svcName = {
            kegg:    'KEGG',
            ncbi:    'NCBI',
            blast:   'NCBI BLAST',
            uniprot: 'UniProt',
            unigene: 'NCBI UniGene',
            pubmed:  'PubMed',
            pdb:     'RCSB PDB',
            server:  'The server',
        }[service] || svc;
        var raw = error && error.message ? error.message : String(error || '');
        var lowered = raw.toLowerCase();
        if (lowered.indexOf('abort') >= 0 || lowered.indexOf('timeout') >= 0) {
            return svcName + ' is slow to respond right now. Try again in a moment.';
        }
        if (lowered.indexOf('failed to fetch') >= 0 || lowered.indexOf('networkerror') >= 0 || lowered.indexOf('network request failed') >= 0) {
            return "Can't reach " + svcName + '. Check your internet connection or try again later.';
        }
        if (raw.match(/\b429\b/) || lowered.indexOf('rate') >= 0) {
            return svcName + ' is rate-limiting your requests. Wait a minute and retry.';
        }
        if (raw.match(/\b5\d\d\b/)) {
            return svcName + ' returned a server error. The service may be temporarily down.';
        }
        if (raw.match(/\b404\b/)) {
            return svcName + " returned 404 — the record doesn't exist or the ID is wrong.";
        }
        if (raw.match(/\b401\b/) || raw.match(/\b403\b/)) {
            return svcName + ' rejected the request. Some endpoints require an email parameter.';
        }
        // Unknown — trim to readable length
        return raw.slice(0, 200) || (svcName + ' returned an unexpected error.');
    }

    // ---- Theme-aware server rendering ----------------------------------
    // Server-rendered figures (genome diagrams, phylo trees, restriction
    // maps, dendrograms, pathway graphs, sequence logos) are delivered as
    // opaque `<img src="data:image/svg+xml;base64,...">` URLs, so CSS cannot
    // restyle them for dark mode — the renderer itself has to know the theme.
    //
    // Rather than touch every fetch() call site (there are ~158 across page
    // scripts and inline template blocks), tag same-origin /api/ requests
    // with the active theme here. Server side, utils.plot_helpers.
    // resolve_theme() reads this header and recolors the figure chrome.
    function currentTheme() {
        return document.documentElement.getAttribute('data-theme') || 'light';
    }

    function installThemeHeader() {
        if (!global.fetch || global.__nugenThemeFetchPatched) return;
        global.__nugenThemeFetchPatched = true;
        var nativeFetch = global.fetch.bind(global);
        global.fetch = function (input, init) {
            try {
                var url = (typeof input === 'string')
                    ? input
                    : (input && input.url) || '';
                // Same-origin API calls only — never leak the header to the
                // external biology services the app also talks to.
                var isAbsolute = /^[a-z][a-z0-9+.-]*:\/\//i.test(url);
                if (!isAbsolute && url.indexOf('/api/') !== -1) {
                    init = init ? Object.assign({}, init) : {};
                    var headers = new Headers(
                        init.headers ||
                        (typeof input === 'object' && input && input.headers) ||
                        {}
                    );
                    headers.set('X-Theme', currentTheme());
                    init.headers = headers;
                }
            } catch (e) {
                // Theming must never be able to block a request.
            }
            return nativeFetch(input, init);
        };
    }

    installThemeHeader();

    // ---- Figure export -------------------------------------------------
    // Server-rendered diagrams arrive as `data:` URLs. Before this, only the
    // phylo page offered a download at all, and it saved SVG payloads under a
    // .png filename — the file opened as garbage in image viewers. This gives
    // every diagram one control that gets the format right: keep the vector
    // SVG (ideal for figures and further editing), or rasterize to a
    // high-resolution PNG for slides and manuscripts.

    function figureMimeType(dataUrl) {
        var m = /^data:([^;,]+)/.exec(String(dataUrl || ''));
        return m ? m[1].toLowerCase() : '';
    }

    function figureExtension(dataUrl) {
        var mime = figureMimeType(dataUrl);
        if (mime.indexOf('svg') !== -1) return 'svg';
        if (mime.indexOf('png') !== -1) return 'png';
        if (mime.indexOf('jpeg') !== -1 || mime.indexOf('jpg') !== -1) return 'jpg';
        return 'img';
    }

    function figureFilename(basename, dataUrl) {
        var safe = String(basename || 'figure')
            .replace(/[^A-Za-z0-9._-]+/g, '_')
            .replace(/^_+|_+$/g, '') || 'figure';
        return safe + '.' + figureExtension(dataUrl);
    }

    function triggerDownload(href, filename) {
        var link = document.createElement('a');
        link.href = href;
        link.download = filename;
        link.style.display = 'none';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }

    // Save the figure exactly as the server produced it.
    function downloadFigure(dataUrl, basename) {
        if (!dataUrl) return;
        triggerDownload(dataUrl, figureFilename(basename, dataUrl));
    }

    // Rasterize a figure to PNG at `scale`x its intrinsic size. Vector SVG
    // has no intrinsic pixel size in some browsers, so fall back to a
    // sensible publication-width canvas when the image reports 0.
    function downloadFigureAsPng(dataUrl, basename, scale) {
        if (!dataUrl) return;
        if (figureExtension(dataUrl) === 'png') {
            downloadFigure(dataUrl, basename);
            return;
        }
        scale = scale || 2;
        var img = new Image();
        img.onload = function () {
            var w = (img.naturalWidth || img.width || 1200) * scale;
            var h = (img.naturalHeight || img.height || 800) * scale;
            var canvas = document.createElement('canvas');
            canvas.width = Math.max(1, Math.round(w));
            canvas.height = Math.max(1, Math.round(h));
            var ctx = canvas.getContext('2d');
            // Diagrams are drawn on an opaque ground server-side, but a
            // transparent PNG on a dark slide is unreadable — paint the
            // current surface color behind it.
            ctx.fillStyle = currentTheme() === 'dark' ? '#141c2e' : '#ffffff';
            ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
            try {
                triggerDownload(canvas.toDataURL('image/png'),
                                figureFilename(basename, 'data:image/png'));
            } catch (e) {
                // Canvas export blocked — fall back to the original asset.
                downloadFigure(dataUrl, basename);
            }
        };
        img.onerror = function () { downloadFigure(dataUrl, basename); };
        img.src = dataUrl;
    }

    // Markup for a pair of export buttons. Uses data-action so the handlers
    // stay out of inline onclick attributes (see the CSP note in app.py).
    function figureExportControls(dataUrl, basename) {
        if (!dataUrl) return '';
        var args = JSON.stringify([dataUrl, basename || 'figure']);
        var esc = escapeHtml(args);
        return '' +
            '<div class="figure-export-controls d-flex gap-2 mt-2 justify-content-end">' +
              '<button type="button" class="btn-app-ghost btn-app-sm" ' +
                'data-action="downloadFigure" data-action-args="' + esc + '" ' +
                'title="Download the original vector figure">' +
                '<i class="fas fa-download"></i> SVG' +
              '</button>' +
              '<button type="button" class="btn-app-ghost btn-app-sm" ' +
                'data-action="downloadFigureAsPng" data-action-args="' + esc + '" ' +
                'title="Download a high-resolution raster copy">' +
                '<i class="fas fa-image"></i> PNG' +
              '</button>' +
            '</div>';
    }

    var FigureExport = {
        mimeType: figureMimeType,
        extension: figureExtension,
        filename: figureFilename,
        download: downloadFigure,
        downloadPng: downloadFigureAsPng,
        controls: figureExportControls,
        inline: inlineFigure,
        inlineControls: inlineFigureControls,
        svgToDataUrl: svgElementToDataUrl,
    };

    // ---- Generic DOM actions for the data-action dispatcher -------------
    // A handful of inline handlers did DOM work directly rather than calling a
    // named function ("this.parentElement.parentElement.remove()"). Inline
    // event handlers cannot be covered by a CSP nonce, so they have to become
    // named actions before 'unsafe-inline' can be dropped from script-src.
    // The dispatcher calls fn.apply(btn, args), so `this` is the button.

    // Removes the row/group a button sits in — the old
    // `this.parentElement.parentElement.remove()`.
    function disableConsensusAndTrimButtons() {
        ['consensusBtn', 'trimBtn'].forEach(function (id) {
            var el = document.getElementById(id);
            if (el) el.disabled = true;
        });
    }

    function removeParentRow() {
        var el = this && this.parentElement && this.parentElement.parentElement;
        if (el && el.remove) el.remove();
    }

    function removeElementById(id) {
        var el = document.getElementById(id);
        if (el) el.remove();
    }

    // ---- data-css: dynamic styles under a strict CSP ---------------------
    // A style="..." attribute in markup is refused when style-src drops
    // 'unsafe-inline' (verified in Chromium: markup style= and
    // setAttribute('style', ...) are blocked, while el.style.prop = x and
    // el.style.cssText are exempt — the CSSOM is not covered by style-src).
    //
    // Most inline styles became utility classes, but a handful carry a value
    // only known at render time (a cluster's colour, a bar's width). Those
    // are emitted as data-css="..." and applied here through the CSSOM.
    // A MutationObserver catches them however they are inserted, so callers
    // building HTML strings need no extra step.
    function applyDataCss(root) {
        if (!root || !root.querySelectorAll) return;
        if (root.hasAttribute && root.hasAttribute('data-css')) {
            try { root.style.cssText = root.getAttribute('data-css'); } catch (e) {}
        }
        root.querySelectorAll('[data-css]').forEach(function (el) {
            try { el.style.cssText = el.getAttribute('data-css'); } catch (e) {}
        });
    }

    function watchDataCss() {
        applyDataCss(document.body || document.documentElement);
        if (!global.MutationObserver) return;
        new MutationObserver(function (records) {
            records.forEach(function (rec) {
                rec.addedNodes.forEach(function (node) {
                    if (node.nodeType === 1) applyDataCss(node);
                });
            });
        }).observe(document.documentElement, { childList: true, subtree: true });
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', watchDataCss);
    } else {
        watchDataCss();
    }

    // ---- Inline-SVG figures ---------------------------------------------
    // Figures are injected into the page rather than shipped as opaque
    // <img src="data:...">, so labels are selectable text and CSS can reach
    // inside them. Export therefore has to serialise the live element rather
    // than hand back a URL it was given.

    function findFigureSvg(from) {
        if (!from) return null;
        const scope = from.closest('.figure-container') || from.parentElement
            || document;
        return scope.querySelector('svg.figure-svg')
            || document.querySelector('svg.figure-svg');
    }

    function svgElementToDataUrl(svg) {
        const clone = svg.cloneNode(true);
        // A standalone file needs the namespace declared; inline it is implied.
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
        const markup = new XMLSerializer().serializeToString(clone);
        return 'data:image/svg+xml;base64,' +
            btoa(unescape(encodeURIComponent(markup)));
    }

    // data-action handlers, bound to the button that was clicked.
    function downloadInlineSvg(basename) {
        const svg = findFigureSvg(this);
        if (svg) downloadFigure(svgElementToDataUrl(svg), basename || 'figure');
    }

    function downloadInlineSvgAsPng(basename) {
        const svg = findFigureSvg(this);
        if (svg) downloadFigureAsPng(svgElementToDataUrl(svg),
                                     basename || 'figure');
    }

    // Export controls for an inline figure. Unlike FigureExport.controls(),
    // which is handed a data URL, these resolve the figure at click time.
    function inlineFigureControls(basename) {
        const name = escapeHtml(basename || 'figure');
        const args = escapeHtml(JSON.stringify([basename || 'figure']));
        return '' +
            '<div class="figure-export-controls d-flex gap-2 mt-2 justify-content-end">' +
              '<button type="button" class="btn-app-ghost btn-app-sm" ' +
                'data-action="downloadInlineSvg" data-action-args="' + args + '" ' +
                'title="Download the vector figure">' +
                '<i class="fas fa-download"></i> SVG</button>' +
              '<button type="button" class="btn-app-ghost btn-app-sm" ' +
                'data-action="downloadInlineSvgAsPng" data-action-args="' + args + '" ' +
                'title="Download a high-resolution raster copy">' +
                '<i class="fas fa-image"></i> PNG</button>' +
            '</div>';
    }

    // Wrap a server-rendered figure in a container with export controls.
    //
    // Accepts either inline SVG markup (the normal case) or a data: URL.
    // The fallback matters: the sequence-logo renderer drops to a PNG when
    // logomaker is unavailable, and a PNG cannot be inlined as SVG.
    function inlineFigure(figure, basename, altText) {
        if (!figure) return '';
        const isSvgMarkup = String(figure).trim().slice(0, 4) === '<svg';
        if (isSvgMarkup) {
            return '<div class="figure-container">' + figure +
                   inlineFigureControls(basename) + '</div>';
        }
        const alt = escapeHtml(altText || basename || 'Figure');
        return '<div class="figure-container">' +
               '<img src="' + figure + '" class="img-fluid u-maxw100" alt="' +
               alt + '">' + figureExportControls(figure, basename) + '</div>';
    }

    global.showLoading = showLoading;
    global.hideLoading = hideLoading;
    global.showAlert = showAlert;
    global.fetchJSON = fetchJSON;
    global.ResultsPanel = ResultsPanel;
    global.friendlyError = friendlyError;
    global.escapeHtml = escapeHtml;
    global.currentTheme = currentTheme;
    global.removeParentRow = removeParentRow;
    global.removeElementById = removeElementById;
    global.applyDataCss = applyDataCss;
    global.downloadInlineSvg = downloadInlineSvg;
    global.downloadInlineSvgAsPng = downloadInlineSvgAsPng;
    global.disableConsensusAndTrimButtons = disableConsensusAndTrimButtons;
    global.FigureExport = FigureExport;
    global.downloadFigure = downloadFigure;
    global.downloadFigureAsPng = downloadFigureAsPng;
    global.toggleShortcutsHelp = toggleShortcutsHelp;
    global.NuGenUtils = {
        showLoading: showLoading,
        hideLoading: hideLoading,
        showAlert: showAlert,
        fetchJSON: fetchJSON,
        escapeHtml: escapeHtml,
        ResultsPanel: ResultsPanel,
        FigureExport: FigureExport,
    };
})(window);
