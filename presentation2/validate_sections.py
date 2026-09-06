"""presentation2/validate_sections.py -- validates every
presentation2/sections/*.json against the schema in presentation2/SCHEMA.md.

Checks: required keys (section-level and per-slide), layout is one of the
nine allowed values, at most 6 bullets each at most ~18 words, at most 3
equations per slide (each with latex + caption), figure.path/figure.script
exist on disk (the figure.path check can be skipped with check_figures=False
for a fast pre-figure-generation pass -- see build.py), and every
`repo_numbers` entry's value appears (within 1e-6) as a numeric token in the
text of the file it names.

CLI: python presentation2/validate_sections.py [--include-sample]
Exit code 0 iff every section file is valid. Prints a per-section slide
count + minutes line and a total.

Also importable: `validate_all(...)` returns (ok, sections) so build.py can
reuse the parsed, validated section data without re-reading the files.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECTIONS_DIR = Path(__file__).resolve().parent / "sections"

VALID_LAYOUTS = {
    "title", "bullets", "bullets+figure", "figure", "equation",
    "two-column", "table", "section-divider", "quote",
}
MAX_BULLETS = 6
MAX_BULLET_WORDS = 18
MAX_EQUATIONS = 3
REPO_NUMBER_TOL = 1e-6

_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _word_count(s: str) -> int:
    return len(s.split())


def _number_in_text(text: str, value: float, tol: float) -> bool:
    for tok in _NUMBER_RE.findall(text):
        try:
            num = float(tok)
        except ValueError:
            continue
        if abs(num - value) <= tol:
            return True
    return False


class _Validator:
    def __init__(self, check_figures: bool):
        self.check_figures = check_figures
        self.errors: list[str] = []

    def err(self, where: str, msg: str) -> None:
        self.errors.append(f"{where}: {msg}")

    def check_repo_number(self, where: str, key: str, entry) -> None:
        if not isinstance(entry, dict) or "value" not in entry or "file" not in entry:
            self.err(where, f"repo_numbers.{key} needs 'value' and 'file'")
            return
        value = entry["value"]
        file_rel = entry["file"]
        target = ROOT / file_rel
        if not target.exists():
            self.err(where, f"repo_numbers.{key} names missing file '{file_rel}'")
            return
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            self.err(where, f"repo_numbers.{key}.value must be numeric")
            return
        try:
            text = target.read_text(encoding="utf-8", errors="strict")
        except (UnicodeDecodeError, OSError):
            return  # not parseable text -- schema exempts this case
        if not _number_in_text(text, float(value), REPO_NUMBER_TOL):
            self.err(
                where,
                f"repo_numbers.{key} value {value} not found in '{file_rel}' "
                f"(tol {REPO_NUMBER_TOL})",
            )

    def validate_slide(self, slide: dict, section_id: str) -> None:
        sid = slide.get("id", "<no id>")
        where = f"{section_id}/{sid}"

        for req in ("id", "layout", "title", "notes", "sources"):
            if req not in slide:
                self.err(where, f"missing required key '{req}'")

        layout = slide.get("layout")
        if layout is not None and layout not in VALID_LAYOUTS:
            self.err(where, f"invalid layout '{layout}' (allowed: {sorted(VALID_LAYOUTS)})")

        sources = slide.get("sources")
        if sources is not None and (not isinstance(sources, list) or len(sources) < 1):
            self.err(where, "sources must be a non-empty list")

        bullets = slide.get("bullets")
        if bullets is not None:
            if len(bullets) > MAX_BULLETS:
                self.err(where, f"{len(bullets)} bullets > max {MAX_BULLETS}")
            for b in bullets:
                wc = _word_count(b)
                if wc > MAX_BULLET_WORDS:
                    self.err(where, f"bullet has {wc} words > ~{MAX_BULLET_WORDS}: '{b[:50]}'")

        equations = slide.get("equations")
        if equations is not None:
            if len(equations) > MAX_EQUATIONS:
                self.err(where, f"{len(equations)} equations > max {MAX_EQUATIONS}")
            for eq in equations:
                if "latex" not in eq or "caption" not in eq:
                    self.err(where, "equation missing 'latex' or 'caption'")

        figure = slide.get("figure")
        if figure is not None:
            for req in ("path", "caption", "script"):
                if req not in figure:
                    self.err(where, f"figure missing '{req}'")
            if "script" in figure and not (ROOT / figure["script"]).exists():
                self.err(where, f"figure script does not exist: {figure['script']}")
            if self.check_figures and "path" in figure and not (ROOT / figure["path"]).exists():
                self.err(where, f"figure path does not exist: {figure['path']} (run its script first)")

        columns = slide.get("columns")
        if columns is not None:
            for col in columns:
                if "heading" not in col or "bullets" not in col:
                    self.err(where, "column missing 'heading' or 'bullets'")

        table = slide.get("table")
        if table is not None:
            if "header" not in table or "rows" not in table:
                self.err(where, "table missing 'header' or 'rows'")

        repo_numbers = slide.get("repo_numbers")
        if repo_numbers is not None:
            for key, entry in repo_numbers.items():
                self.check_repo_number(where, key, entry)


def validate_all(section_dir=None, include_sample: bool = False, check_figures: bool = True):
    """Validate every presentation2/sections/NN_*.json file. Returns
    (ok, sections): `sections` is the list of parsed section dicts that
    passed structural parsing (section "00" excluded unless include_sample).
    """
    section_dir = Path(section_dir) if section_dir else SECTIONS_DIR
    validator = _Validator(check_figures=check_figures)

    paths = sorted(section_dir.glob("[0-9][0-9]_*.json"))
    sections = []
    total_minutes = 0
    total_slides = 0

    for path in paths:
        if path.stem[:2] == "00" and not include_sample:
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            validator.err(path.name, f"invalid JSON: {e}")
            continue

        for req in ("section", "title", "minutes", "slides"):
            if req not in data:
                validator.err(path.name, f"missing required key '{req}'")
        if "slides" not in data:
            continue

        seen_ids: set[str] = set()
        for slide in data["slides"]:
            validator.validate_slide(slide, data.get("section", path.stem))
            sid = slide.get("id")
            if sid in seen_ids:
                validator.err(path.name, f"duplicate slide id '{sid}'")
            seen_ids.add(sid)

        n = len(data["slides"])
        minutes = data.get("minutes", 0)
        print(f"  {path.name}: {n} slides, {minutes} min")
        total_slides += n
        total_minutes += minutes if isinstance(minutes, (int, float)) else 0
        sections.append(data)

    print(f"TOTAL: {len(sections)} sections, {total_slides} slides, {total_minutes} minutes")

    ok = len(validator.errors) == 0
    if not ok:
        print("\nERRORS:")
        for e in validator.errors:
            print(f"  {e}")
    return ok, sections


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-sample", action="store_true",
                         help="also validate presentation2/sections/00_*.json")
    parser.add_argument("--skip-figure-check", action="store_true",
                         help="skip the figure.path-exists check (use before figure scripts run)")
    parser.add_argument("--sections-dir", default=None,
                         help="directory of NN_*.json section files (default: presentation2/sections)")
    args = parser.parse_args()
    ok, _ = validate_all(section_dir=args.sections_dir, include_sample=args.include_sample,
                          check_figures=not args.skip_figure_check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
