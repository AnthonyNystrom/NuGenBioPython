"""
Capture visual baseline screenshots for every page + every Load Example flow.

Baselines live under `tests/visual/baseline/`. Re-run this script whenever
you accept a visual change; commit the new PNGs alongside the code.

One-time setup
--------------
    pip install playwright pytest-playwright pixelmatch-python pillow
    playwright install chromium

Run
---
    # In one terminal: start the Flask app
    FLASK_DEBUG=false SECRET_KEY=test /opt/anaconda3/envs/nugenbio/bin/python app.py
    # In another:
    /opt/anaconda3/envs/nugenbio/bin/python tests/visual/take_baseline.py

Workflow
--------
1. Before any Phase 1+ UI change: rerun with the current design, commit PNGs.
2. Make the UI change.
3. `pytest tests/visual/test_visual_regression.py` — see diffs.
4. Human reviews the diff images in `tests/visual/diff/`.
5. If diffs are intended, rerun this script to replace baseline; commit
   both the code and the updated baselines in the same PR.

This is the gate that protects Load Example flows across phases 1–4.
"""
import os
import sys
import time

BASE_URL = os.environ.get("APP_URL", "http://127.0.0.1:9000")
VIEWPORT = {"width": 1440, "height": 900}
OUTDIR = os.path.join(os.path.dirname(__file__), "baseline")


def main():
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sys.exit(
            "playwright not installed. See top-of-file instructions:\n"
            "  pip install playwright pytest-playwright\n"
            "  playwright install chromium"
        )

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
    from tests.fixtures.inventory import PAGE_URLS, LOAD_EXAMPLE_BUTTONS

    os.makedirs(OUTDIR, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport=VIEWPORT)
        page = ctx.new_page()

        # 1 — Every page URL
        for url in PAGE_URLS:
            slug = url.strip("/").replace("/", "_") or "index"
            dest = os.path.join(OUTDIR, f"page__{slug}.png")
            print(f"-> {url}")
            page.goto(BASE_URL + url, wait_until="networkidle")
            # Wait for webfonts (Inter) to finish loading before snapshot.
            # Without this, first-visit captures can be 5-15 px taller than
            # warm captures because fallback-font metrics differ subtly.
            try:
                page.evaluate("document.fonts.ready")
                page.wait_for_function(
                    "document.fonts && document.fonts.status === 'loaded'",
                    timeout=5000,
                )
            except Exception:
                pass
            time.sleep(0.5)
            page.screenshot(path=dest, full_page=True)

        # 2 — Every Load Example button result
        # Buttons are located by their onclick function name (stable hook)
        # rather than by visible text (often just "Example" with an icon).
        for entry in LOAD_EXAMPLE_BUTTONS:
            slug = entry["page"].strip("/").replace("/", "_") or "index"
            fn = entry["onclick_fn"]
            dest = os.path.join(OUTDIR, f"example__{slug}__{fn}.png")
            print(f"-> {entry['page']} | {fn}")
            try:
                page.goto(BASE_URL + entry["page"], wait_until="networkidle")
                # Try to activate the containing tab pane first (many example
                # buttons live inside inactive tabs). Finds the enclosing
                # .tab-pane and clicks its matching .nav-link.
                selector = f'[onclick^="{fn}("]'
                page.evaluate(
                    """(sel) => {
                        const btn = document.querySelector(sel);
                        if (!btn) return;
                        const pane = btn.closest('.tab-pane');
                        if (pane && pane.id) {
                            const tab = document.querySelector(
                                `[data-bs-target="#${pane.id}"], [href="#${pane.id}"]`
                            );
                            if (tab) tab.click();
                        }
                    }""",
                    selector,
                )
                time.sleep(0.3)
                # Click the button; force=True because the load* handlers only
                # manipulate form values and don't require visibility.
                page.locator(selector).first.click(timeout=3000, force=True)
                time.sleep(1.5)
                page.screenshot(path=dest, full_page=True)
            except Exception as e:
                print(f"   skipped: {str(e).splitlines()[0][:100]}")

        browser.close()
    print(f"\nBaselines written to {OUTDIR}")


if __name__ == "__main__":
    main()
