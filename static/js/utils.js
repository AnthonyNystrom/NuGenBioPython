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

    global.showLoading = showLoading;
    global.hideLoading = hideLoading;
    global.showAlert = showAlert;
    global.fetchJSON = fetchJSON;
    global.NuGenUtils = { showLoading: showLoading, hideLoading: hideLoading, showAlert: showAlert, fetchJSON: fetchJSON, escapeHtml: escapeHtml };
})(window);
