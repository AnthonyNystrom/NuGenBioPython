// Shared UI helpers for NuGenBioPython frontend.
// Loaded from base.html before all other page scripts, so these globals are
// available on every page.
(function (global) {
    'use strict';

    function showLoading(btnId, label) {
        var btn = document.getElementById(btnId);
        if (!btn) return;
        if (!btn.dataset.originalHtml) {
            btn.dataset.originalHtml = btn.innerHTML;
        }
        btn.disabled = true;
        var text = label || 'Processing...';
        btn.innerHTML = '<span class="spinner-border spinner-border-sm me-2"></span>' + escapeHtml(text);
    }

    function hideLoading(btnId, originalHtml) {
        var btn = document.getElementById(btnId);
        if (!btn) return;
        btn.disabled = false;
        var restore = originalHtml || btn.dataset.originalHtml;
        if (restore !== undefined) {
            btn.innerHTML = restore;
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

        var container =
            document.querySelector('.card-body.p-4') ||
            document.querySelector('.card-body') ||
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
        if (overlay) { overlay.remove(); return; }
        overlay = document.createElement('div');
        overlay.id = 'ng-shortcuts-help';
        overlay.className = 'ng-shortcuts-overlay';
        overlay.innerHTML =
            '<div class="ng-shortcuts-card" role="dialog" aria-label="Keyboard shortcuts">' +
                '<div class="ng-shortcuts-title">Keyboard shortcuts</div>' +
                '<dl>' +
                    '<dt><kbd>/</kbd></dt><dd>Focus the first input on the page</dd>' +
                    '<dt><kbd>g</kbd> <kbd>h</kbd></dt><dd>Go to Dashboard</dd>' +
                    '<dt><kbd>?</kbd></dt><dd>Show / hide this help</dd>' +
                    '<dt><kbd>Esc</kbd></dt><dd>Close panels and overlays</dd>' +
                '</dl>' +
                '<button type="button" class="btn btn-sm btn-secondary">Close</button>' +
            '</div>';
        overlay.addEventListener('click', function (e) {
            if (e.target === overlay) overlay.remove();
        });
        overlay.querySelector('button').addEventListener('click', function () { overlay.remove(); });
        document.body.appendChild(overlay);
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

    let chord = null;
    let chordTimer = null;
    document.addEventListener('keydown', function (e) {
        // Always allow Esc to close the help overlay
        if (e.key === 'Escape') {
            const overlay = document.getElementById('ng-shortcuts-help');
            if (overlay) overlay.remove();
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
        // Chord handling: g then h → /
        if (chord === 'g') {
            clearTimeout(chordTimer);
            chord = null;
            if (e.key === 'h') {
                e.preventDefault();
                window.location.href = '/';
                return;
            }
        } else if (e.key === 'g') {
            chord = 'g';
            chordTimer = setTimeout(function () { chord = null; }, 1200);
        }
    });

    global.showLoading = showLoading;
    global.hideLoading = hideLoading;
    global.showAlert = showAlert;
    global.fetchJSON = fetchJSON;
    global.ResultsPanel = ResultsPanel;
    global.NuGenUtils = {
        showLoading: showLoading,
        hideLoading: hideLoading,
        showAlert: showAlert,
        fetchJSON: fetchJSON,
        escapeHtml: escapeHtml,
        ResultsPanel: ResultsPanel,
    };
})(window);
