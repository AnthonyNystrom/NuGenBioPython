"""Capture the README screenshots.

Unlike `take_baseline.py` (which snaps every page for regression testing),
this script produces a curated set of polished screenshots for the README.
Each shot loads an example, configures the form if needed, submits, and
then scrolls so the result is in frame.

Run
---
    # In one terminal:
    FLASK_DEBUG=false SECRET_KEY=test /opt/anaconda3/envs/nugenbio/bin/python app.py
    # In another:
    /opt/anaconda3/envs/nugenbio/bin/python tests/visual/take_readme_shots.py
"""
import os
import sys
import time

BASE_URL = os.environ.get("APP_URL", "http://127.0.0.1:9000")
VIEWPORT = {"width": 1440, "height": 900}
OUTDIR = os.path.join(os.path.dirname(__file__), "..", "..",
                       "docs", "screenshots")
OUTDIR = os.path.abspath(OUTDIR)


# Each shot has these knobs:
#   url                - page to open
#   setup              - JS fn name loaded via data-action or window[fn]()
#   prep_js            - raw JS to run before submitting (e.g. pick a radio)
#   open_tab           - data-bs-target selector to activate a tab first
#   submit             - selector to click (form submit) after setup
#   wait               - seconds after submit before screenshot
#   scroll_to          - element selector to scroll into view before shot
#   clip               - selector for element.screenshot() (overrides full page)
#   full               - bool, full-page screenshot (default False, viewport)
SHOTS = [
    {
        "name": "01_dashboard.png",
        "url": "/",
        "wait": 0.6,
        "full": True,
    },
    {
        "name": "02_structure_viewer.png",
        "url": "/structure",
        "setup": "loadExample",
        "wait": 3.0,   # 3Dmol cartoon render
        "scroll_to": "#view3dPanel",
    },
    {
        "name": "03_genome_diagram.png",
        "url": "/genomediagram",
        "setup": "loadBasicExampleComplexGenome",
        "submit": "#createBasicDiagramBtn",
        "wait": 3.0,
        "scroll_to": "#basicDiagramDisplay",
    },
    {
        "name": "04_sequence_logo.png",
        "url": "/motifs",
        "setup": "loadMotifExample",
        "submit": "#createMotifBtn",
        "wait": 2.5,
        "scroll_to": "#createLogoView",
    },
    {
        "name": "05_phylo_tree.png",
        "url": "/phylo",
        "setup": "loadExampleTreeComplex",
        # "Parse Tree" only handles file uploads. The textarea goes through
        # the separate "Parse String" button.
        "submit": "[data-action=\"parseTreeString\"]",
        "wait": 3.0,
        "scroll_to": "#treeVisualization",
    },
    {
        "name": "06_pathway_network.png",
        "url": "/pathway",
        "setup": "loadReactionExampleComplex",
        "open_tab": "#visualization",
        "submit": "#generateVizBtn",
        "wait": 3.5,
        "scroll_to": "#pathwayVisualization",
    },
    {
        "name": "07_restriction_map.png",
        "url": "/restriction",
        "setup": "loadMapExampleComplex",
        "submit": "#generateMapBtn",
        "wait": 2.5,
        "scroll_to": "#mapVisualization",
    },
    {
        "name": "08_clustering_dendrogram.png",
        "url": "/clustering",
        "setup": "loadClusteringExampleComplex",
        "prep_js": """
            // Switch to Hierarchical so the submit produces a dendrogram.
            const sel = document.getElementById('clusterMethod');
            if (sel) { sel.value = 'hierarchical'; sel.dispatchEvent(new Event('change')); }
        """,
        "submit": "#clusterBtn",
        "wait": 3.0,
        "scroll_to": "#clusterVisualization",
    },
]


def activate_tab_of(page, selector):
    """If the button lives inside a hidden .tab-pane, click the nav tab
    that owns it so the button becomes visible + clickable."""
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


def click_setup(page, fn):
    selector = f'[data-action="{fn}"]'
    if page.locator(selector).count():
        activate_tab_of(page, selector)
        time.sleep(0.25)
        page.locator(selector).first.click(force=True, timeout=4000)
    else:
        # Fallback: most load* handlers are exposed on window.
        page.evaluate(f"typeof {fn} === 'function' && {fn}()")


def main():
    from playwright.sync_api import sync_playwright

    os.makedirs(OUTDIR, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport=VIEWPORT, device_scale_factor=2)
        page = ctx.new_page()

        for shot in SHOTS:
            print(f"-> {shot['name']}  ({shot['url']})")
            page.goto(BASE_URL + shot["url"], wait_until="networkidle")
            # Wait for fonts to settle so first-vs-warm captures match.
            try:
                page.wait_for_function(
                    "document.fonts && document.fonts.status === 'loaded'",
                    timeout=4000,
                )
            except Exception:
                pass
            time.sleep(0.4)

            if shot.get("setup"):
                try:
                    click_setup(page, shot["setup"])
                    time.sleep(0.4)
                except Exception as e:
                    print(f"   (setup {shot['setup']} failed: {e})")

            if shot.get("prep_js"):
                try:
                    page.evaluate(shot["prep_js"])
                    time.sleep(0.2)
                except Exception as e:
                    print(f"   (prep_js failed: {e})")

            if shot.get("open_tab"):
                # e.g. "#visualization" — click its nav-link
                tab_sel = shot["open_tab"]
                try:
                    page.evaluate(
                        """(id) => {
                            const tab = document.querySelector(
                                `[data-bs-target="${id}"], [href="${id}"]`);
                            if (tab) tab.click();
                        }""", tab_sel)
                    time.sleep(0.4)
                except Exception as e:
                    print(f"   (open_tab {tab_sel} failed: {e})")

            if shot.get("submit"):
                sel = shot["submit"]
                try:
                    activate_tab_of(page, sel)
                    time.sleep(0.2)
                    page.locator(sel).first.click(force=True, timeout=4000)
                except Exception as e:
                    print(f"   (submit {sel} failed: {e})")

            time.sleep(shot["wait"])

            if shot.get("scroll_to"):
                try:
                    page.evaluate(
                        """(sel) => {
                            const el = document.querySelector(sel);
                            if (el) el.scrollIntoView(
                                { block: 'start', behavior: 'instant' });
                        }""",
                        shot["scroll_to"],
                    )
                    # Nudge up so the element isn't flush against the header.
                    page.evaluate("window.scrollBy(0, -80)")
                    time.sleep(0.3)
                except Exception as e:
                    print(f"   (scroll_to failed: {e})")

            dest = os.path.join(OUTDIR, shot["name"])
            page.screenshot(path=dest, full_page=bool(shot.get("full")))
            print(f"   saved {dest}")

        browser.close()

    print(f"\nAll screenshots in: {OUTDIR}")


if __name__ == "__main__":
    main()
