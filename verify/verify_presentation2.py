"""Acceptance checks for the presentation2 renderer pipeline (build.py,
validate_sections.py, mathimg.py, render_pptx.js, render_html.py).

This runs the build against an isolated temp copy of the schema-conformance
sample section only (presentation2/sections/00_sample.json), never against
the live presentation2/sections/ directory -- that directory is being
written concurrently by other section-content workers, and this check must
be deterministic regardless of the state of their in-progress files.

What this confirms:
  * `python presentation2/build.py --include-sample --sections-dir <isolated>`
    exits 0 (validate -> figures -> equations -> pptx -> html -> office
    validate, per presentation2/SCHEMA.md's renderer contract);
  * the produced pptx opens with python-pptx, its slide count equals the
    sample section's slide count, every slide's rendered text contains its
    JSON `title`, every slide's notes slide matches its JSON `notes`
    verbatim, and every slide carrying an `equations` array -- on any
    layout, not only `equation` -- contains at least as many picture shapes
    as it has `equations` (exercises the new `equation+figure` and
    `bullets+table` layouts, and a `bullets+figure` slide with an equation);
  * every `equation`/`equation+figure` slide that also carries `bullets`
    renders each bullet's text in both the pptx slide text and the HTML
    page (SCHEMA.md: bullets render on every layout that carries them,
    including those two);
  * out/presentation2/index.html has the same number of `<section` elements
    as slides, embeds no external `<script src=`, and stays under 16 MB;
  * validate_sections.py rejects speaker notes shorter than the 115-255 word
    hard bound and a figure PNG whose pixel size is not exactly 1600x900 or
    1200x1200 (negative tests using temp copies).

Sandbox note: this creates a tempfile.TemporaryDirectory(); inside a
read-only or workspace-write sandbox that can raise PermissionError. If that
is the ONLY failure, report `TESTS: pass (verify_presentation2 skipped:
sandbox temp dir)` per CLAUDE.md and let the orchestrator re-run outside
the sandbox.

CLI: python verify/verify_presentation2.py
Exit code 0 iff every check passes; prints 'N/N presentation2 checks passed'.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE

ROOT = Path(__file__).resolve().parents[1]
PRESENTATION2 = ROOT / "presentation2"
SAMPLE_JSON = PRESENTATION2 / "sections" / "00_sample.json"
OUT_DIR = ROOT / "out" / "presentation2"
PPTX_PATH = OUT_DIR / "qd_physics_2h.pptx"
HTML_PATH = OUT_DIR / "index.html"

MAX_HTML_BYTES = 16 * 1024 * 1024

checks: list[bool] = []


def ck(ok: bool, name: str) -> None:
    checks.append(bool(ok))
    print(("PASS" if ok else "FAIL"), name)


def check_repo_numbers_how(sample: dict) -> None:
    """repo_numbers form (b): validate_sections.py --section must accept the
    sample's computed ('how') entries as-is, and must reject a copy where a
    'how' entry's expected value has been tampered with."""
    result_ok = subprocess.run(
        [sys.executable, str(PRESENTATION2 / "validate_sections.py"),
         "--section", str(SAMPLE_JSON)],
        cwd=str(ROOT), capture_output=True, text=True,
    )
    ck(result_ok.returncode == 0,
       "validate_sections.py --section 00_sample.json (correct 'how' values) exits 0")

    mutated = json.loads(json.dumps(sample))  # deep copy
    tampered = False
    for slide in mutated["slides"]:
        for entry in slide.get("repo_numbers", {}).values():
            if "how" in entry:
                entry["value"] = entry["value"] + 1000.0  # deliberately wrong
                tampered = True
    ck(tampered, "sample fixture has at least one 'how' entry to tamper with")

    with tempfile.TemporaryDirectory() as tmp:
        bad_path = Path(tmp) / SAMPLE_JSON.name
        bad_path.write_text(json.dumps(mutated), encoding="utf-8")
        result_bad = subprocess.run(
            [sys.executable, str(PRESENTATION2 / "validate_sections.py"),
             "--section", str(bad_path)],
            cwd=str(ROOT), capture_output=True, text=True,
        )
    ck(result_bad.returncode != 0,
       "validate_sections.py --section rejects a tampered 'how' expected value")


def check_notes_length_validation(sample: dict) -> None:
    """validate_sections.py must reject a slide whose speaker notes are far
    shorter than the schema's 115-255 word hard bound."""
    mutated = json.loads(json.dumps(sample))
    mutated["slides"][0]["notes"] = "Too short to pass the word-count check."

    with tempfile.TemporaryDirectory() as tmp:
        bad_path = Path(tmp) / SAMPLE_JSON.name
        bad_path.write_text(json.dumps(mutated), encoding="utf-8")
        result = subprocess.run(
            [sys.executable, str(PRESENTATION2 / "validate_sections.py"),
             "--section", str(bad_path)],
            cwd=str(ROOT), capture_output=True, text=True,
        )
    ck(result.returncode != 0,
       "validate_sections.py --section rejects notes shorter than 115 words")


_FAKE_PNG_PATH = PRESENTATION2 / "figures" / "out" / "_verify_tmp_badsize.png"


def _write_fake_png(width: int, height: int, path: Path) -> None:
    """A PNG whose first 24 bytes (signature + IHDR width/height) are well
    formed -- all validate_sections.py's `_png_dimensions` reads -- but whose
    body is otherwise not a real, decodable image; good enough for a check
    that never asks a real image library to open the file."""
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr_len = (13).to_bytes(4, "big")
    ihdr_body = (
        width.to_bytes(4, "big") + height.to_bytes(4, "big") + bytes([8, 6, 0, 0, 0])
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(sig + ihdr_len + b"IHDR" + ihdr_body + b"\x00\x00\x00\x00")


def check_figure_size_validation(sample: dict) -> None:
    """validate_sections.py must reject a figure PNG whose pixel size is not
    exactly 1600x900 or 1200x1200."""
    mutated = json.loads(json.dumps(sample))
    mutated_slide = next(s for s in mutated["slides"] if s.get("figure"))
    mutated_slide["figure"]["path"] = "presentation2/figures/out/_verify_tmp_badsize.png"

    _write_fake_png(800, 600, _FAKE_PNG_PATH)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            bad_path = Path(tmp) / SAMPLE_JSON.name
            bad_path.write_text(json.dumps(mutated), encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(PRESENTATION2 / "validate_sections.py"),
                 "--section", str(bad_path)],
                cwd=str(ROOT), capture_output=True, text=True,
            )
    finally:
        if _FAKE_PNG_PATH.exists():
            _FAKE_PNG_PATH.unlink()
    ck(result.returncode != 0,
       "validate_sections.py --section rejects a figure that is not 1600x900 or 1200x1200 px")


def main() -> int:
    sample = json.loads(SAMPLE_JSON.read_text(encoding="utf-8"))
    slides_by_id = {s["id"]: s for s in sample["slides"]}
    expected_n_slides = len(sample["slides"])

    check_repo_numbers_how(sample)
    check_notes_length_validation(sample)
    check_figure_size_validation(sample)

    with tempfile.TemporaryDirectory() as tmp:
        sections_dir = Path(tmp)
        shutil.copy(SAMPLE_JSON, sections_dir / SAMPLE_JSON.name)

        result = subprocess.run(
            [sys.executable, str(PRESENTATION2 / "build.py"),
             "--include-sample", "--sections-dir", str(sections_dir)],
            cwd=str(ROOT),
        )
        ck(result.returncode == 0, "presentation2/build.py --include-sample exits 0")
        if result.returncode != 0:
            print(f"{sum(checks)}/{len(checks)} presentation2 checks passed")
            return 1

    prs = Presentation(str(PPTX_PATH))
    slides = list(prs.slides)
    ck(len(slides) == expected_n_slides,
       f"pptx slide count == sample slide count ({len(slides)} vs {expected_n_slides})")

    html_text = HTML_PATH.read_text(encoding="utf-8")

    for i, pptx_slide in enumerate(slides):
        json_slide = sample["slides"][i] if i < len(sample["slides"]) else None
        label = json_slide["id"] if json_slide else f"index {i}"

        slide_text = "\n".join(
            shape.text_frame.text for shape in pptx_slide.shapes if shape.has_text_frame
        )
        title_ok = bool(json_slide) and json_slide["title"] in slide_text
        ck(title_ok, f"slide {label}: rendered text contains its JSON title")

        notes_text = ""
        if pptx_slide.has_notes_slide:
            notes_text = pptx_slide.notes_slide.notes_text_frame.text
        notes_ok = bool(json_slide) and notes_text == json_slide["notes"]
        ck(notes_ok, f"slide {label}: notes slide matches its JSON notes verbatim")

        if json_slide and json_slide.get("equations"):
            # Every layout that carries an `equations` array must render
            # them as pictures, not only the `equation` layout (SCHEMA.md).
            n_pictures = sum(
                1 for shape in pptx_slide.shapes if shape.shape_type == MSO_SHAPE_TYPE.PICTURE
            )
            n_equations = len(json_slide["equations"])
            ck(n_pictures >= n_equations,
               f"slide {label} ({json_slide['layout']}): "
               f"{n_pictures} pictures >= {n_equations} equations")

        if json_slide and json_slide.get("layout") in ("equation", "equation+figure") \
                and json_slide.get("bullets"):
            # SCHEMA.md: bullets render on every layout that carries them,
            # including "equation" and "equation+figure" -- assert each
            # bullet's text actually shows up in both renderers' output.
            for bullet in json_slide["bullets"]:
                ck(bullet in slide_text,
                   f"slide {label} ({json_slide['layout']}): pptx text contains bullet "
                   f"'{bullet[:40]}'")
                ck(bullet in html_text,
                   f"slide {label} ({json_slide['layout']}): html contains bullet "
                   f"'{bullet[:40]}'")

    n_sections = len(re.findall(r"<section", html_text))
    ck(n_sections == expected_n_slides,
       f"index.html has the same number of <section> elements as slides "
       f"({n_sections} vs {expected_n_slides})")
    ck("<script src=" not in html_text, "index.html has no external <script src=...>")
    ck(not re.search(r"""(?:src|href)\s*=\s*["']https?://""", html_text),
       "index.html has no external http(s) src/href reference")
    html_bytes = HTML_PATH.stat().st_size
    ck(html_bytes < MAX_HTML_BYTES, f"index.html is under 16 MB (got {html_bytes} bytes)")

    print(f"{sum(checks)}/{len(checks)} presentation2 checks passed")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PermissionError as e:
        print(f"SKIPPED (sandbox temp dir): {e}")
        print("See CLAUDE.md 'Sandbox note for workers' -- re-run outside the sandbox.")
        sys.exit(0)
