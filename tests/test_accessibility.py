"""Form controls must have accessible names.

64 inputs, selects and textareas relied on a placeholder or nothing at all.
A placeholder is not an accessible name: it disappears on focus and screen
readers treat it inconsistently, so a user navigating by keyboard hears
"edit text, blank" with no indication of what the field is for. Most were in
dynamically-built rows (genome diagram tracks, pathway reactants, feature
editors), which is exactly where a user is doing repetitive data entry.
"""
import pathlib
import re

import pytest

CONTROL = re.compile(r'<(input|select|textarea)\b[^>]*>', re.S)


def _source_files(repo_root):
    yield from pathlib.Path(repo_root, "templates").rglob("*.html")
    yield from pathlib.Path(repo_root, "static", "js").glob("*.js")


def _is_comment_line(text, offset):
    line_start = text.rfind("\n", 0, offset) + 1
    return text[line_start:offset].lstrip().startswith(("//", "*", "#"))


def test_every_form_control_has_an_accessible_name(repo_root):
    offenders = []
    for path in _source_files(repo_root):
        text = path.read_text()
        for match in CONTROL.finditer(text):
            tag = match.group(0)
            if 'type="hidden"' in tag:
                continue
            if "aria-label" in tag or re.search(r"\bid=", tag):
                continue
            if _is_comment_line(text, match.start()):
                continue          # prose describing markup, not markup
            line = text[:match.start()].count("\n") + 1
            offenders.append(f"{path.name}:{line}  {' '.join(tag.split())[:60]}")
    assert not offenders, (
        "these controls have no accessible name; add aria-label, or an id "
        f"with a matching <label for>: {offenders}")


def test_no_label_is_derived_from_sample_data(repo_root):
    """A placeholder showing example values is not a name.

    The first automated pass produced aria-label="0.5,0.6,0.4,0.7,0.5..."
    from a placeholder that was demo data.
    """
    offenders = []
    for path in _source_files(repo_root):
        for label in re.findall(r'aria-label="([^"]*)"', path.read_text()):
            stripped = label.replace(",", "").replace(".", "").replace(" ", "")
            if stripped.isdigit() or label.count(",") >= 3:
                offenders.append(f"{path.name}: {label!r}")
    assert not offenders, f"aria-label looks like sample data: {offenders}"


def test_labels_are_not_empty(repo_root):
    offenders = []
    for path in _source_files(repo_root):
        for label in re.findall(r'aria-label="([^"]*)"', path.read_text()):
            if not label.strip():
                offenders.append(path.name)
    assert not offenders, f"empty aria-label in: {offenders}"


@pytest.mark.parametrize("url", ["/genomediagram", "/features", "/pathway"])
def test_pages_with_dynamic_rows_still_render(client, url):
    assert client.get(url).status_code == 200
