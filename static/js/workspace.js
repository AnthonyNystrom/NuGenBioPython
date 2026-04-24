// Workspace — shared sessionStorage clipboard for cross-tool workflow.
//
// Adds a collapsible right-rail panel listing items the user has saved
// (sequences, structures, trees, etc.). Any tool can push items in via
// Workspace.add(...) or pull them out via Workspace.list(...).
//
// Phase-4 contract: additive. Never removes or mutates existing DOM;
// defaults collapsed so pages look identical on first load for users who
// don't opt in. Storage is sessionStorage (per-tab, not persisted).
(function (global) {
    'use strict';

    const STORAGE_KEY = 'nugen_workspace_v1';
    const listeners = [];
    let state = null;   // lazily loaded

    // ---- Storage ------------------------------------------------------
    function load() {
        if (state !== null) return state;
        try {
            const raw = sessionStorage.getItem(STORAGE_KEY);
            state = raw ? JSON.parse(raw) : { items: [] };
        } catch (e) {
            state = { items: [] };
        }
        if (!Array.isArray(state.items)) state.items = [];
        return state;
    }
    function save() {
        try {
            sessionStorage.setItem(STORAGE_KEY, JSON.stringify(state));
        } catch (e) {
            // Quota exceeded; silently drop.
        }
        listeners.forEach(function (fn) {
            try { fn(state); } catch (_) {}
        });
    }

    function makeId() {
        return 'ws-' + Date.now().toString(36) + '-' + Math.random().toString(36).slice(2, 8);
    }

    // ---- Public API ---------------------------------------------------
    function add(item) {
        if (!item || !item.type) throw new Error('Workspace.add: item.type required');
        const s = load();
        const entry = {
            id: item.id || makeId(),
            type: item.type,
            name: String(item.name || '(untitled)'),
            data: item.data == null ? '' : item.data,
            meta: item.meta || {},
            createdAt: Date.now(),
        };
        // Deduplicate by (type + data) — if identical data already exists, don't add
        const dup = s.items.find(function (it) {
            return it.type === entry.type && it.data === entry.data;
        });
        if (dup) return dup;
        s.items.unshift(entry);
        // Cap at 50 items to stay under sessionStorage quota comfortably
        if (s.items.length > 50) s.items.length = 50;
        save();
        return entry;
    }

    function list(type) {
        const s = load();
        if (!type) return s.items.slice();
        return s.items.filter(function (i) { return i.type === type; });
    }

    function get(id) {
        return load().items.find(function (i) { return i.id === id; }) || null;
    }

    function remove(id) {
        const s = load();
        const before = s.items.length;
        s.items = s.items.filter(function (i) { return i.id !== id; });
        if (s.items.length !== before) save();
    }

    function clearAll() {
        state = { items: [] };
        save();
    }

    function count(type) {
        return list(type).length;
    }

    function onChange(fn) {
        listeners.push(fn);
        return function off() {
            const i = listeners.indexOf(fn);
            if (i >= 0) listeners.splice(i, 1);
        };
    }

    // ---- Convenience for callers --------------------------------------
    // Populate a <textarea> or <input> with the item's data. Covers the
    // "Load from workspace" → "fill this form field" flow.
    function fillInto(targetEl, item) {
        if (!targetEl || !item) return;
        if ('value' in targetEl) {
            targetEl.value = item.data;
            targetEl.dispatchEvent(new Event('input', { bubbles: true }));
            targetEl.dispatchEvent(new Event('change', { bubbles: true }));
        } else {
            targetEl.textContent = item.data;
        }
    }

    // ---- UI: right-rail panel -----------------------------------------
    const TYPE_LABELS = {
        sequence:  { icon: 'fa-dna',            title: 'Sequences'  },
        alignment: { icon: 'fa-align-center',   title: 'Alignments' },
        structure: { icon: 'fa-cube',           title: 'Structures' },
        tree:      { icon: 'fa-project-diagram',title: 'Trees'      },
        motif:     { icon: 'fa-puzzle-piece',   title: 'Motifs'     },
    };

    function setPanelOpen(panel, open) {
        panel.classList.toggle('ws-collapsed', !open);
        const topBtn = document.getElementById('topbarWorkspaceBtn');
        if (topBtn) {
            topBtn.setAttribute('aria-expanded', open ? 'true' : 'false');
            topBtn.classList.toggle('is-active', open);
        }
    }

    function ensurePanel() {
        let panel = document.getElementById('ws-panel');
        if (panel) return panel;

        // Panel is a right-side drawer; toggle lives on the topbar
        // (#topbarWorkspaceBtn). The legacy floating .ws-toggle hook is
        // preserved for back-compat but hidden via shell.css.
        panel = document.createElement('aside');
        panel.id = 'ws-panel';
        panel.className = 'ws-panel ws-collapsed';
        panel.setAttribute('aria-label', 'Workspace');
        panel.innerHTML =
            '<div class="ws-body">' +
                '<div class="ws-header d-flex align-items-center justify-content-between">' +
                    '<strong><i class="fas fa-layer-group me-2 text-primary"></i>Workspace</strong>' +
                    '<div>' +
                        '<button type="button" class="btn-app-secondary ws-clear btn-link btn-app-sm p-0 me-3" title="Clear all">Clear</button>' +
                        '<button type="button" class="btn-app-secondary ws-close btn-link btn-app-sm p-0" title="Close" aria-label="Close">' +
                            '<i class="fas fa-times"></i>' +
                        '</button>' +
                    '</div>' +
                '</div>' +
                '<div class="ws-hint small text-muted mb-2">' +
                    'Save sequences, structures, or trees here and reuse them in other tools.' +
                '</div>' +
                '<div class="ws-list"></div>' +
            '</div>';
        document.body.appendChild(panel);

        // The in-topbar button is the primary opener; wire it once on render.
        const topBtn = document.getElementById('topbarWorkspaceBtn');
        if (topBtn) {
            topBtn.addEventListener('click', function () {
                setPanelOpen(panel, panel.classList.contains('ws-collapsed'));
            });
        }
        panel.querySelector('.ws-close').addEventListener('click', function () {
            setPanelOpen(panel, false);
        });
        // Esc closes the drawer
        document.addEventListener('keydown', function (e) {
            if (e.key === 'Escape' && !panel.classList.contains('ws-collapsed')) {
                setPanelOpen(panel, false);
            }
        });
        panel.querySelector('.ws-clear').addEventListener('click', function () {
            if (confirm('Clear all workspace items?')) clearAll();
        });
        return panel;
    }

    function renderPanel() {
        const panel = ensurePanel();
        const s = load();
        const listEl = panel.querySelector('.ws-list');

        // Reflect item count in the topbar badge (R1 primary location) and
        // keep the legacy in-panel badge (#ws-badge) in sync if present.
        const topBadge = document.getElementById('topbarWorkspaceBadge');
        const legacyBadge = panel.querySelector('#ws-badge');
        const n = s.items.length;
        [topBadge, legacyBadge].forEach(function (b) {
            if (!b) return;
            if (n) { b.hidden = false; b.textContent = String(n); }
            else   { b.hidden = true; }
        });

        if (!s.items.length) {
            listEl.innerHTML =
                '<div class="ws-empty">' +
                    '<i class="fas fa-layer-group ws-empty-icon" aria-hidden="true"></i>' +
                    '<h6 class="ws-empty-title">Workspace is empty</h6>' +
                    '<p class="ws-empty-text small">' +
                        'Run any analysis (BLAST, alignment, restriction, structure…) and click ' +
                        '<strong><i class="fas fa-layer-group"></i> Save</strong> on the result card. ' +
                        'Saved items live here across pages so you can reuse them.' +
                    '</p>' +
                    '<a class="btn-app-sm btn-app-secondary" href="/sequence">' +
                        '<i class="fas fa-play me-1"></i> Try sequence analysis' +
                    '</a>' +
                '</div>';
            return;
        }

        // Group by type
        const byType = {};
        s.items.forEach(function (it) {
            (byType[it.type] = byType[it.type] || []).push(it);
        });

        let html = '';
        Object.keys(byType).forEach(function (type) {
            const spec = TYPE_LABELS[type] || { icon: 'fa-tag', title: type };
            html += '<div class="ws-group mb-2">';
            html += '<div class="ws-group-title small fw-semibold mb-1">' +
                '<i class="fas ' + spec.icon + ' me-1"></i>' + spec.title +
                ' <span class="text-muted">(' + byType[type].length + ')</span></div>';
            byType[type].forEach(function (it) {
                const preview = String(it.data || '').slice(0, 40).replace(/\s+/g, ' ');
                html += '<div class="ws-item" data-ws-id="' + it.id + '">' +
                    '<div class="ws-item-name text-truncate" title="' + escapeAttr(it.name) + '">' +
                        escapeHtml(it.name) +
                    '</div>' +
                    '<div class="ws-item-preview small text-muted text-truncate">' +
                        escapeHtml(preview) +
                    '</div>' +
                    '<div class="ws-item-actions">' +
                        '<button type="button" class="btn-app-secondary btn-app-sm btn-link p-0 me-2" data-ws-action="copy">Copy</button>' +
                        '<button type="button" class="btn-app-secondary btn-app-sm btn-link p-0 me-2" data-ws-action="use">Use</button>' +
                        '<button type="button" class="btn-app-secondary btn-app-sm btn-link p-0 text-danger" data-ws-action="remove">Remove</button>' +
                    '</div>' +
                '</div>';
            });
            html += '</div>';
        });
        listEl.innerHTML = html;

        // Wire per-item buttons
        listEl.querySelectorAll('.ws-item').forEach(function (row) {
            const id = row.dataset.wsId;
            row.querySelector('[data-ws-action="remove"]').addEventListener('click', function () {
                remove(id);
            });
            row.querySelector('[data-ws-action="copy"]').addEventListener('click', function () {
                const item = get(id);
                if (item && navigator.clipboard) {
                    navigator.clipboard.writeText(String(item.data || ''));
                }
            });
            row.querySelector('[data-ws-action="use"]').addEventListener('click', function () {
                const item = get(id);
                if (!item) return;
                // Fill the currently-focused textarea/input, or the first
                // sequence textarea on the page as a fallback.
                const active = document.activeElement;
                let target = null;
                if (active && (active.tagName === 'TEXTAREA' || active.tagName === 'INPUT')) {
                    target = active;
                } else {
                    target = document.querySelector('textarea.sequence-display') ||
                             document.querySelector('textarea');
                }
                if (target) {
                    fillInto(target, item);
                    if (global.showAlert) showAlert('Loaded "' + item.name + '" into form', 'success');
                } else {
                    if (global.showAlert) showAlert('No input field found on this page', 'warning');
                }
            });
        });
    }

    function escapeHtml(s) {
        return String(s == null ? '' : s)
            .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
    }
    function escapeAttr(s) { return escapeHtml(s); }

    // Re-render whenever the state changes
    onChange(renderPanel);

    // Initialize after the DOM is ready (so we can append to body)
    function init() {
        ensurePanel();
        renderPanel();
    }
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }

    global.Workspace = { add: add, list: list, get: get, remove: remove, clear: clearAll, count: count, onChange: onChange, fillInto: fillInto, renderPanel: renderPanel };
})(window);
