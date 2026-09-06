"""presentation2/validate_sections.py -- validates every
presentation2/sections/*.json against the schema in presentation2/SCHEMA.md.

Checks: required keys (section-level and per-slide), layout is one of the
eleven allowed values, at most 6 bullets each at most ~18 words, at most 3
equations per slide (each with latex + caption), speaker notes word count
(115-255 words, hard bounds), figure.path/figure.script exist on disk and
figure.path is exactly 1600x900 or 1200x1200 px, read straight from the PNG
header (both checks skip with check_figures=False for a fast
pre-figure-generation pass -- see build.py), and every `repo_numbers` entry's
value against its source, in one of two forms:
  (a) literal: {"value": ..., "file": ...} -- the value must appear (within
      `tolerance`, default 1e-6, absolute) as a numeric token in the text of
      the named file.
  (b) computed: {"value": ..., "how": "<python expression>"} -- the
      expression is evaluated in a subprocess (`python -c ...`) with the
      repository root on sys.path and cwd, given a 30s timeout, and the
      result is compared to `value` within `tolerance` (default 1e-6,
      relative). A `"file"` alongside `how` is documentation only (not
      checked). Mismatch, timeout, or an exception in the expression is
      reported naming the slide, key, and expected/computed values.

CLI: python presentation2/validate_sections.py [--include-sample]
     python presentation2/validate_sections.py --section <path/to.json>
`--section` validates exactly one section JSON file (ignoring
--sections-dir/--include-sample), for spot-checking a single file while
others are still in progress.
Exit code 0 iff every section file is valid. Prints a per-section slide
count + minutes line and a total.

Also importable: `validate_all(...)` returns (ok, sections) so build.py can
reuse the parsed, validated section data without re-reading the files.
"""
from __future__ import annotations

import argparse
import json
import re
import struct
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECTIONS_DIR = Path(__file__).resolve().parent / "sections"

VALID_LAYOUTS = {
    "title", "bullets", "bullets+figure", "figure", "equation",
    "equation+figure", "two-column", "table", "bullets+table",
    "section-divider", "quote",
}
MAX_BULLETS = 6
MAX_BULLET_WORDS = 18
MAX_EQUATIONS = 3
MIN_NOTES_WORDS = 115
MAX_NOTES_WORDS = 255
ALLOWED_FIGURE_SIZES = {(1600, 900), (1200, 1200)}
REPO_NUMBER_TOL = 1e-6
HOW_TIMEOUT_S = 30

_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


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


def _png_dimensions(path: Path):
    """Read (width, height) straight from a PNG's IHDR chunk (bytes 16:24 of
    the file), no image library needed. Returns None if `path` is not a
    well-formed PNG."""
    try:
        with path.open("rb") as f:
            header = f.read(24)
    except OSError:
        return None
    if len(header) < 24 or header[:8] != _PNG_SIGNATURE:
        return None
    width, height = struct.unpack(">II", header[16:24])
    return width, height


class _Validator:
    def __init__(self, check_figures: bool):
        self.check_figures = check_figures
        self.errors: list[str] = []

    def err(self, where: str, msg: str) -> None:
        self.errors.append(f"{where}: {msg}")

    def check_repo_number(self, where: str, key: str, entry) -> None:
        if not isinstance(entry, dict) or "value" not in entry:
            self.err(where, f"repo_numbers.{key} needs 'value'")
            return
        value = entry["value"]
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            self.err(where, f"repo_numbers.{key}.value must be numeric")
            return
        tolerance = entry.get("tolerance", REPO_NUMBER_TOL)
        if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool):
            self.err(where, f"repo_numbers.{key}.tolerance must be numeric")
            return

        if "how" in entry:
            self.check_computed_repo_number(where, key, entry, float(value), float(tolerance))
            return

        if "file" not in entry:
            self.err(where, f"repo_numbers.{key} needs 'file' (or 'how')")
            return
        file_rel = entry["file"]
        target = ROOT / file_rel
        if not target.exists():
            self.err(where, f"repo_numbers.{key} names missing file '{file_rel}'")
            return
        try:
            text = target.read_text(encoding="utf-8", errors="strict")
        except (UnicodeDecodeError, OSError):
            return  # not parseable text -- schema exempts this case
        if not _number_in_text(text, float(value), float(tolerance)):
            self.err(
                where,
                f"repo_numbers.{key} value {value} not found in '{file_rel}' "
                f"(tol {tolerance})",
            )

    def check_computed_repo_number(self, where: str, key: str, entry, expected: float,
                                    tolerance: float) -> None:
        how = entry["how"]
        if not isinstance(how, str) or not how.strip():
            self.err(where, f"repo_numbers.{key}.how must be a non-empty python expression")
            return
        code = f"import sys; sys.path.insert(0, '.'); print(repr(float({how})))"
        try:
            result = subprocess.run(
                [sys.executable, "-c", code],
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=HOW_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            self.err(
                where,
                f"repo_numbers.{key}.how timed out after {HOW_TIMEOUT_S}s (how={how!r})",
            )
            return

        if result.returncode != 0:
            stderr_lines = result.stderr.strip().splitlines()
            last_line = stderr_lines[-1] if stderr_lines else "unknown error"
            self.err(
                where,
                f"repo_numbers.{key}.how raised an exception evaluating {how!r}: {last_line}",
            )
            return

        try:
            computed = float(result.stdout.strip())
        except ValueError as e:
            self.err(
                where,
                f"repo_numbers.{key}.how produced unparseable output {result.stdout!r} "
                f"(how={how!r}): {e}",
            )
            return

        denom = abs(expected) if expected != 0 else 1.0
        rel_err = abs(computed - expected) / denom
        if rel_err > tolerance:
            self.err(
                where,
                f"repo_numbers.{key}.how mismatch: expected {expected}, computed {computed} "
                f"(rel err {rel_err:.3g} > tol {tolerance}, how={how!r})",
            )

    def validate_slide(self, slide: dict, section_id: str) -> None:
        sid = slide.get("id", "<no id>")
        where = f"{section_id}/{sid}"

        for req in ("id", "layout", "title", "notes", "sources"):
            if req not in slide:
                self.err(where, f"missing required key '{req}'")

        notes = slide.get("notes")
        if isinstance(notes, str):
            wc = _word_count(notes)
            if wc < MIN_NOTES_WORDS or wc > MAX_NOTES_WORDS:
                self.err(
                    where,
                    f"notes has {wc} words, must be {MIN_NOTES_WORDS}-{MAX_NOTES_WORDS}",
                )

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
            if self.check_figures and "path" in figure:
                target = ROOT / figure["path"]
                if not target.exists():
                    self.err(where, f"figure path does not exist: {figure['path']} (run its script first)")
                else:
                    dims = _png_dimensions(target)
                    if dims not in ALLOWED_FIGURE_SIZES:
                        got = f"{dims[0]}x{dims[1]}" if dims else "unreadable/not a PNG"
                        self.err(
                            where,
                            f"figure {figure['path']} is {got} px, must be exactly "
                            f"1600x900 or 1200x1200",
                        )

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


def validate_all(section_dir=None, include_sample: bool = False, check_figures: bool = True,
                  only_file=None):
    """Validate every presentation2/sections/NN_*.json file. Returns
    (ok, sections): `sections` is the list of parsed section dicts that
    passed structural parsing (section "00" excluded unless include_sample).

    If `only_file` is given, validate exactly that one file instead of
    globbing `section_dir`, and never skip it for being a "00" section.
    """
    validator = _Validator(check_figures=check_figures)

    if only_file is not None:
        paths = [Path(only_file)]
    else:
        section_dir = Path(section_dir) if section_dir else SECTIONS_DIR
        paths = sorted(section_dir.glob("[0-9][0-9]_*.json"))
    sections = []
    total_minutes = 0
    total_slides = 0

    for path in paths:
        if only_file is None and path.stem[:2] == "00" and not include_sample:
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
    parser.add_argument("--section", default=None,
                         help="validate exactly one section JSON file "
                              "(ignores --sections-dir/--include-sample)")
    args = parser.parse_args()
    if args.section:
        ok, _ = validate_all(only_file=args.section, check_figures=not args.skip_figure_check)
    else:
        ok, _ = validate_all(section_dir=args.sections_dir, include_sample=args.include_sample,
                              check_figures=not args.skip_figure_check)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
