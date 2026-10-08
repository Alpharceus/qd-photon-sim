"""Acceptance checks for the presentation2 renderer pipeline (build.py,
validate_sections.py, mathimg.py, render_pptx.js, render_html.py,
figures/deckstyle.py).

This runs the build against an isolated copy of the schema-conformance
sample section only (presentation2/sections/00_sample.json), never against
the live presentation2/sections/ directory -- that directory is being
written concurrently by other section-content workers, and this check must
be deterministic regardless of the state of their in-progress files. It
builds into `--out-dir` (default out/presentation2_selftest) and refuses to
run against the deliverable directory out/presentation2.

What this confirms:
  * `build.py` exits 0 for BOTH variants of the sample deck -- with reveal
    (into <out-dir>) and with --no-reveal (into <out-dir>/no_reveal); each
    build runs validate -> figures -> equations -> pptx -> html -> the pptx
    skill's office validate.py, which is also re-run here as its own check;
  * the produced pptx opens with python-pptx, its slide count equals the
    sample section's slide count, every slide's rendered text contains its
    JSON `title`, every slide's notes slide matches its JSON `notes`
    verbatim, and every slide carrying an `equations` array -- on any
    layout, not only `equation` -- contains at least as many picture shapes
    as it has `equations`;
  * every `equation`/`equation+figure` slide that also carries `bullets`
    renders each bullet's text in both the pptx slide text and the HTML
    page, and every reveal slide renders every bullet in the pptx;
  * reveal (SCHEMA.md "Reveal and build-up frames"): validate_sections.py
    rejects each broken reveal/frames rule on a mutated copy (length
    mismatch, non-contiguous steps, step > 8, reveal on a title slide,
    [click] count mismatch, [click] on a slide without reveal, frames count,
    figure.path != frames[-1], non-increasing frame steps, unknown key, bad
    reveal_effect, wrong-size frame PNG) with the specific error;
  * the injected pptx timing XML: every slide part is well-formed XML with
    unique shape ids; a reveal slide has exactly one <p:timing> as the last
    child of <p:sld>, one click group per reveal step (== max step), click
    group g animates exactly the shapes named rv-g-*, every spTgt/bldP spid
    resolves to a shape on that slide (bldP only to text shapes), cTn ids
    are unique, and the presetID matches reveal_effect; non-reveal slides
    have no <p:timing>;
  * the --no-reveal pptx contains no p:timing anywhere and no rv- shapes,
    and its HTML has no data-reveal attribute;
  * index.html: same number of `<section` elements as slides; each reveal
    slide carries exactly the expected data-reveal attributes (steps 1..max);
    a frame stack holds every frame with the complete graph last; exactly
    one <script>, inline, and it parses (`node --check`); no external
    resource of any kind (script src, link, @import, http(s) or
    protocol-relative src/href/url(), non-data img src); and with every
    script REMOVED the page still shows all content: no CSS rule hides
    anything outside the body.presenting scope (except non-final frames and
    the JS-created notes overlay), no element is hidden inline, and every
    bullet / column / table cell / caption text is present; under 16 MB;
  * every equation image the sample uses (equations and math-only bullets)
    was written by the current mathimg cache version and its inked content
    lies strictly inside the canvas with >= EQ_MIN_MARGIN_PX white pixels on
    every side (a cropped equation touches an edge);
  * HTML dark theme: every slide's --accent-dark reaches WCAG AA contrast
    (>= 4.5:1, WCAG 2.x relative-luminance formula, computed here
    independently) against the dark --bg, and the dark-theme CSS (both the
    prefers-color-scheme block and [data-theme="dark"]) colours the
    accent-coloured titles (h2, .big-title, .quote) with var(--accent-dark);
  * validate_sections.py rejects speaker notes shorter than the 115-255 word
    hard bound, a figure PNG whose pixel size is not exactly 1600x900 or
    1200x1200, and a tampered `how` value (negative tests using temp copies);
    it counts only the AUTHORED notes: a long appended "Sources:" block
    neither rescues too-short notes nor pushes valid notes over the bound;
  * layout (fix round 2; the build fixture is the sample with a subtitle
    added to every titled slide, so every layout is exercised with one):
    on every pptx slide no two content shapes overlap (build-up frames at
    the identical box excepted), every content shape ends above the one-line
    source band, the source band is ONE line "Source: <file> @ <hash>; ..."
    and no slide text contains "(undefined)"; a slide without repo_numbers
    gets the band exactly when its `sources` name a repository file (the
    fixture adds one to an equation slide); every subtitle is rendered in
    the pptx and the HTML; table rows use >= 12 pt; the pptx notes are the
    JSON notes verbatim plus, for a slide with repo_numbers, a "Sources:"
    block with one line per entry (the HTML notes carry the same block);
  * inline math: a $...$ span in running text is rendered as Unicode text
    (mathimg.inline_runs, subscripts as sub/superscript runs), never as raw
    LaTeX, in both the pptx and the HTML.

Sandbox note: this creates temporary fixture files; inside a
read-only or workspace-write sandbox that can raise PermissionError. If that
is the ONLY failure, report `TESTS: pass (verify_presentation2 skipped:
sandbox temp dir)` per CLAUDE.md and let the orchestrator re-run outside
the sandbox.

CLI: python verify/verify_presentation2.py [--out-dir out/presentation2_selftest]
Exit code 0 iff every check passes; prints 'N/N presentation2 checks passed'.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from lxml import etree
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE

ROOT = Path(__file__).resolve().parents[1]
PRESENTATION2 = ROOT / "presentation2"
SAMPLE_JSON = PRESENTATION2 / "sections" / "00_sample.json"
DELIVERABLE_DIR = ROOT / "out" / "presentation2"
DEFAULT_OUT_DIR = ROOT / "out" / "presentation2_selftest"
OFFICE_VALIDATE = Path.home() / ".claude" / "skills" / "pptx" / "scripts" / "office" / "validate.py"

MAX_HTML_BYTES = 16 * 1024 * 1024
NS = {
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
}
P = "{%s}" % NS["p"]
PRESET_ID = {"fade": "10", "appear": "1"}

sys.path.insert(0, str(PRESENTATION2))
import build  # noqa: E402
import mathimg  # noqa: E402
import validate_sections  # noqa: E402

EMU = 914400
SOURCE_BAND_TOP_IN = 6.74   # render_pptx.js SOURCE_Y: content must end above it
FOOTER_TOP_IN = 7.07        # render_pptx.js FOOTER_Y

EQ_MIN_MARGIN_PX = 4
WCAG_AA = 4.5

checks: list[bool] = []
try:  # check names carry Unicode (inline-math cases); never die on a cp1252 console
    sys.stdout.reconfigure(errors="backslashreplace")
except (AttributeError, ValueError):
    pass


def ck(ok: bool, name: str) -> None:
    checks.append(bool(ok))
    print(("PASS" if ok else "FAIL"), name)


# ---------------------------------------------------------------------------
# validator negative tests
# ---------------------------------------------------------------------------

def run_validator(section: dict, scratch: Path, name: str = "00_sample.json"):
    scratch.mkdir(parents=True, exist_ok=True)
    path = scratch / name
    path.write_text(json.dumps(section), encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(PRESENTATION2 / "validate_sections.py"), "--section", str(path)],
        cwd=str(ROOT), capture_output=True, text=True,
    )


def check_repo_numbers_how(sample: dict, scratch: Path) -> None:
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
    if result_ok.returncode != 0:
        print(result_ok.stdout[-2000:])

    mutated = json.loads(json.dumps(sample))  # deep copy
    tampered = False
    for slide in mutated["slides"]:
        for entry in slide.get("repo_numbers", {}).values():
            if "how" in entry:
                entry["value"] = entry["value"] + 1000.0  # deliberately wrong
                tampered = True
    ck(tampered, "sample fixture has at least one 'how' entry to tamper with")
    result_bad = run_validator(mutated, scratch)
    ck(result_bad.returncode != 0,
       "validate_sections.py --section rejects a tampered 'how' expected value")


def check_notes_length_validation(sample: dict, scratch: Path) -> None:
    """validate_sections.py must reject a slide whose speaker notes are far
    shorter than the schema's 115-255 word hard bound."""
    mutated = json.loads(json.dumps(sample))
    mutated["slides"][0]["notes"] = "Too short to pass the word-count check."
    result = run_validator(mutated, scratch)
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


def check_figure_size_validation(sample: dict, scratch: Path) -> None:
    """validate_sections.py must reject a figure PNG whose pixel size is not
    exactly 1600x900 or 1200x1200."""
    mutated = json.loads(json.dumps(sample))
    mutated_slide = next(s for s in mutated["slides"] if s.get("figure") and not s.get("reveal"))
    mutated_slide["figure"]["path"] = "presentation2/figures/out/_verify_tmp_badsize.png"
    _write_fake_png(800, 600, _FAKE_PNG_PATH)
    try:
        result = run_validator(mutated, scratch)
    finally:
        if _FAKE_PNG_PATH.exists():
            _FAKE_PNG_PATH.unlink()
    ck(result.returncode != 0,
       "validate_sections.py --section rejects a figure that is not 1600x900 or 1200x1200 px")


def _slide(section: dict, sid: str) -> dict:
    return next(s for s in section["slides"] if s["id"] == sid)


def check_reveal_validation(sample: dict, scratch: Path) -> None:
    """Each SCHEMA.md reveal / frames rule, broken one at a time on a copy of
    the sample (repo_numbers slides dropped for speed; they are covered
    above). The validator must exit non-zero AND name the broken rule."""
    base = json.loads(json.dumps(sample))
    base["slides"] = [s for s in base["slides"] if "repo_numbers" not in s]
    reveal_fig = next(s["id"] for s in base["slides"]
                      if s.get("reveal", {}).get("figure_frames") is not None)
    reveal_cols = next(s["id"] for s in base["slides"] if s.get("reveal", {}).get("columns"))
    plain = next(s["id"] for s in base["slides"]
                 if "reveal" not in s and s["layout"] == "bullets")
    title = next(s["id"] for s in base["slides"] if s["layout"] == "title")

    ok_result = run_validator(base, scratch)
    ck(ok_result.returncode == 0, "reveal fixture copy (unmutated) validates cleanly")
    if ok_result.returncode != 0:
        print(ok_result.stdout[-2000:])

    def case(name: str, expect: str, mutate) -> None:
        m = json.loads(json.dumps(base))
        mutate(m)
        r = run_validator(m, scratch)
        ok = r.returncode != 0 and expect in r.stdout
        ck(ok, f"reveal validation rejects {name} (expects '{expect}')")
        if not ok:
            print(r.stdout[-1500:])

    case("an array whose length differs from its target", "reveal.bullets has 2 entries",
         lambda m: _slide(m, reveal_fig)["reveal"].__setitem__("bullets", [0, 1]))
    case("non-contiguous steps", "contiguous from 1",
         lambda m: _slide(m, reveal_cols)["reveal"].__setitem__("columns", [1, 3]))
    case("a step above 8", "integers 0-8",
         lambda m: _slide(m, reveal_cols)["reveal"].__setitem__("columns", [1, 9]))
    case("reveal on a title slide", "not allowed on the 'title' layout",
         lambda m: _slide(m, title).__setitem__("reveal", {"bullets": [1]}))

    def drop_click(m):
        s = _slide(m, reveal_fig)
        s["notes"] = s["notes"].replace("[click]", "", 1)
    case("a [click] count that differs from max(step)", "markers, expected", drop_click)

    def click_on_plain(m):
        s = _slide(m, plain)
        s["notes"] = s["notes"] + " [click]"
    case("[click] markers on a slide without reveal", "has no reveal block", click_on_plain)

    def one_frame(m):
        f = _slide(m, reveal_fig)["figure"]
        f["frames"] = f["frames"][-1:]
        _slide(m, reveal_fig)["reveal"]["figure_frames"] = [0]
    case("a single figure frame", "must be 2-4", one_frame)

    def path_not_last(m):
        f = _slide(m, reveal_fig)["figure"]
        f["path"] = f["frames"][0]
    case("figure.path != frames[-1]", "figure.path must equal figure.frames[-1]", path_not_last)
    case("non-increasing frame steps", "strictly increasing",
         lambda m: _slide(m, reveal_fig)["reveal"].__setitem__("figure_frames", [0, 3, 2]))
    case("an unknown reveal key", "unknown reveal key",
         lambda m: _slide(m, reveal_cols)["reveal"].__setitem__("rows", [1]))
    case("an invalid reveal_effect", "reveal_effect 'wipe'",
         lambda m: _slide(m, reveal_cols).__setitem__("reveal_effect", "wipe"))

    def bad_frame(m):
        _slide(m, reveal_fig)["figure"]["frames"][0] = "presentation2/figures/out/_verify_tmp_badsize.png"
    _write_fake_png(800, 600, _FAKE_PNG_PATH)
    try:
        case("a figure frame that is not 1600x900 or 1200x1200", "figure frame 0", bad_frame)
    finally:
        if _FAKE_PNG_PATH.exists():
            _FAKE_PNG_PATH.unlink()


# ---------------------------------------------------------------------------
# builds
# ---------------------------------------------------------------------------

def run_build(sections_dir: Path, out_dir: Path, no_reveal: bool) -> int:
    cmd = [sys.executable, str(PRESENTATION2 / "build.py"), "--include-sample",
           "--sections-dir", str(sections_dir), "--out-dir", str(out_dir)]
    if no_reveal:
        cmd.append("--no-reveal")
    r = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
    return r.returncode


def office_validate(pptx: Path, label: str) -> None:
    if not OFFICE_VALIDATE.exists():
        ck(False, f"{label}: office validate.py present at {OFFICE_VALIDATE}")
        return
    r = subprocess.run([sys.executable, str(OFFICE_VALIDATE), str(pptx)],
                       cwd=str(ROOT), capture_output=True, text=True)
    ck(r.returncode == 0, f"{label}: pptx skill office/validate.py passes")
    if r.returncode != 0:
        print(r.stdout[-2000:])


# ---------------------------------------------------------------------------
# pptx checks
# ---------------------------------------------------------------------------

def slide_parts(pptx: Path) -> dict:
    with zipfile.ZipFile(pptx) as z:
        return {n: z.read(n) for n in z.namelist()}


def check_pptx_common(prs, sample: dict, html_text: str, label: str) -> None:
    slides = list(prs.slides)
    n = len(sample["slides"])
    ck(len(slides) == n, f"{label}: pptx slide count == sample slide count ({len(slides)} vs {n})")
    for i, pptx_slide in enumerate(slides):
        json_slide = sample["slides"][i] if i < n else None
        sid = json_slide["id"] if json_slide else f"index {i}"
        slide_text = "\n".join(
            shape.text_frame.text for shape in pptx_slide.shapes if shape.has_text_frame
        )
        ck(bool(json_slide) and json_slide["title"] in slide_text,
           f"{label} slide {sid}: rendered text contains its JSON title")
        notes_text = pptx_slide.notes_slide.notes_text_frame.text if pptx_slide.has_notes_slide else ""
        if not json_slide:
            ck(False, f"{label} slide {sid}: has a JSON slide")
            continue
        rn = json_slide.get("repo_numbers") or {}
        if rn:
            head = json_slide["notes"] + "\n\n" + build.NOTES_SOURCES_HEADING + "\n"
            tail = notes_text[len(head):].split("\n") if notes_text.startswith(head) else []
            ok = len(tail) == len(rn) and all(line.startswith(f"{k} = ") for line, k in zip(tail, rn))
            ck(ok, f"{label} slide {sid}: notes = JSON notes verbatim + 'Sources:' + one line per repo_numbers entry")
        else:
            ck(notes_text == json_slide["notes"],
               f"{label} slide {sid}: notes slide matches its JSON notes verbatim (no repo_numbers)")
        ck(validate_sections.authored_notes(notes_text) .strip() == json_slide["notes"].strip(),
           f"{label} slide {sid}: validate_sections.authored_notes() recovers exactly the authored notes")
        check_slide_layout(pptx_slide, json_slide, slide_text, label)
        if json_slide.get("equations"):
            n_pictures = sum(1 for s in pptx_slide.shapes if s.shape_type == MSO_SHAPE_TYPE.PICTURE)
            n_eq = len(json_slide["equations"])
            ck(n_pictures >= n_eq,
               f"{label} slide {sid} ({json_slide['layout']}): {n_pictures} pictures >= {n_eq} equations")
        if json_slide.get("layout") in ("equation", "equation+figure") and json_slide.get("bullets"):
            for bullet in json_slide["bullets"]:
                ck(mathimg.inline_plain(bullet) in slide_text,
                   f"{label} slide {sid} ({json_slide['layout']}): pptx text contains bullet '{bullet[:40]}'")
                ck(bullet in html_text,
                   f"{label} slide {sid} ({json_slide['layout']}): html contains bullet '{bullet[:40]}'")
        if json_slide.get("reveal"):
            text_bullets = [b for b in json_slide.get("bullets", []) if not re.fullmatch(r"\$[^$]+\$", b.strip())]
            missing = [b for b in text_bullets if mathimg.inline_plain(b) not in slide_text]
            ck(not missing, f"{label} slide {sid} (reveal): every text bullet is in the pptx text")


def _box(shape):
    return (shape.left / EMU, shape.top / EMU, (shape.left + shape.width) / EMU,
            (shape.top + shape.height) / EMU)


def _overlap(a, b, tol=0.01) -> bool:
    return (min(a[2], b[2]) - max(a[0], b[0]) > tol) and (min(a[3], b[3]) - max(a[1], b[1]) > tol)


def check_slide_layout(pptx_slide, js: dict, slide_text: str, label: str) -> None:
    """Fix round 2 layout rules on one rendered slide (geometry of the shape
    boxes the renderer wrote; PowerPoint draws text inside them)."""
    sid = js["id"]
    shapes = list(pptx_slide.shapes)
    footer, source, content = [], [], []
    for sh in shapes:
        top = sh.top / EMU
        txt = sh.text_frame.text if sh.has_text_frame else ""
        if top >= FOOTER_TOP_IN - 0.01:
            footer.append(sh)
        elif top >= SOURCE_BAND_TOP_IN - 0.01 and txt.startswith("Source: "):
            source.append(sh)
        else:
            content.append(sh)
    ck("(undefined)" not in slide_text,
       f"{label} slide {sid}: no '(undefined)' in the slide text")
    if js.get("repo_numbers"):
        line = source[0].text_frame.text if len(source) == 1 else ""
        files = {build.entry_source(e) for e in js["repo_numbers"].values()}
        ck(len(source) == 1 and "\n" not in line and "\x0b" not in line
           and (source[0].height / EMU) <= 0.3 and any(f in line for f in files),
           f"{label} slide {sid}: one-line source band names a repo_numbers file ({line[:70]!r})")
    elif build.source_files(js):
        # band-only provenance: no repo_numbers, but `sources` names repo files
        line = source[0].text_frame.text if len(source) == 1 else ""
        ck(len(source) == 1 and "\n" not in line and "\x0b" not in line
           and (source[0].height / EMU) <= 0.3 and any(f in line for f in build.source_files(js)),
           f"{label} slide {sid}: one-line source band names a file from its sources ({line[:70]!r})")
    else:
        ck(not source, f"{label} slide {sid}: no source band without repo_numbers or file sources")
    boxes = [(_box(sh), sh) for sh in content]
    low = [sh.name for b, sh in boxes if b[3] > SOURCE_BAND_TOP_IN + 0.01]
    ck(not low, f"{label} slide {sid}: every content shape ends above the source band ({low[:3]})")
    clashes = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (a, sa), (b, sb) = boxes[i], boxes[j]
            same_box = all(abs(x - y) < 0.005 for x, y in zip(a, b))
            both_pics = sa.shape_type == MSO_SHAPE_TYPE.PICTURE and sb.shape_type == MSO_SHAPE_TYPE.PICTURE
            if same_box and both_pics:
                continue  # build-up frames stack at the identical box by design
            if _overlap(a, b):
                clashes.append((sa.name, sb.name))
    ck(not clashes, f"{label} slide {sid}: no two content shapes overlap ({clashes[:2]})")
    if js.get("subtitle"):
        ck(mathimg.inline_plain(js["subtitle"]) in slide_text,
           f"{label} slide {sid} ({js['layout']}): subtitle rendered")
    for sh in content:
        if sh.has_table:
            sizes = [r.font.size.pt for row in sh.table.rows for c in row.cells
                     for para in c.text_frame.paragraphs for r in para.runs if r.font.size is not None]
            ck(bool(sizes) and min(sizes) >= 12,
               f"{label} slide {sid}: table text >= 12 pt (min {min(sizes) if sizes else None})")
    raw_math = re.findall(r"\$[^$]*\\[A-Za-z]+[^$]*\$", slide_text)
    ck(not raw_math, f"{label} slide {sid}: no raw inline LaTeX in the slide text ({raw_math[:1]})")


def check_notes_sources_validation(sample: dict, scratch: Path) -> None:
    """The appended "Sources:" list never counts toward the notes bound."""
    filler = "\n\n" + build.NOTES_SOURCES_HEADING + "\n" + "\n".join(
        f"k{i} = {i} (out/rt_edge/verdict.md @ 0000000)" for i in range(120))
    mutated = json.loads(json.dumps(sample))
    mutated["slides"][1]["notes"] = mutated["slides"][1]["notes"] + filler
    r = run_validator(mutated, scratch)
    ck(r.returncode == 0,
       "validate_sections.py accepts valid notes followed by a 480-word appended 'Sources:' block")
    if r.returncode != 0:
        print(r.stdout[-1500:])
    mutated = json.loads(json.dumps(sample))
    mutated["slides"][1]["notes"] = "Too short to pass the word-count check." + filler
    r = run_validator(mutated, scratch)
    ck(r.returncode != 0,
       "validate_sections.py rejects too-short authored notes even with a long 'Sources:' block")


def check_inline_math() -> None:
    """mathimg.inline_runs: known conversions (the expectation is written
    here, not produced by the code)."""
    cases = [
        (r"critical angle $\theta_c = \arcsin(1/n) = 16.6^\circ$.",
         [("critical angle θ", ""), ("c", "sub"), (" = arcsin(1/n) = 16.6°.", "")]),
        (r"$n_{eff} = 3.25$", [("n", ""), ("eff", "sub"), (" = 3.25", "")]),
        (r"$0.587\ \mu m^2$", [("0.587 μm", ""), ("2", "sup")]),
        (r"$e^{-\alpha L}$", [("e", ""), ("−αL", "sup")]),
        ("plain text, no math", [("plain text, no math", "")]),
        # braced scripts OUTSIDE $...$ are markup too (final polish round):
        # the 02-11 linewidth law must never show a raw "^{...}".
        ("b/(e^{E_LO/kT} − 1)", [("b/(e", ""), ("E_LO/kT", "sup"), (" − 1)", "")]),
        ("Loading 1 − e^{-µ}, J_{0,SRH}", [("Loading 1 − e", ""), ("−µ", "sup"), (", J", ""), ("0,SRH", "sub")]),
        ("cm^-2 and E_X stay literal", [("cm^-2 and E_X stay literal", "")]),
    ]
    for text, want in cases:
        got = mathimg.inline_runs(text)
        ck(got == want, f"inline math {text!r} -> {want} (got {got})")
    s = "b/(e^{E_LO/kT} − 1)"
    ck(s in mathimg.inline_map([s]) and "^{" not in mathimg.inline_plain(s),
       "inline_map carries a braced-script-only string (render_pptx.js converts it; no raw '^{')")
    import render_html  # noqa: E402  (same helper the HTML build uses)
    ck(render_html.inline_text(s) == "b/(e<sup>E_LO/kT</sup> − 1)",
       f"render_html.inline_text renders the braced script as <sup> (got {render_html.inline_text(s)!r})")


def check_band_only_provenance() -> None:
    """A slide without repo_numbers whose `sources` name repository files gets
    the one-line source band (no notes 'Sources:' block); literature-only
    sources give no band."""
    lit = {"id": "x-1", "sources": ["Sze & Ng, Physics of Semiconductor Devices 3rd ed."], "notes": "n"}
    ck(build.slide_provenance(lit, {}) is None, "literature-only sources: no source band")
    band = {"id": "x-2", "notes": "n",
            "sources": ["Sze & Ng, ch. 2", "fsim_core/transport.py: Diode.iv", "not/a/file.py: x"]}
    prov = build.slide_provenance(band, {"fsim_core/transport.py": "abc1234"})
    ck(bool(prov) and prov["line"] == "Source: fsim_core/transport.py @ abc1234" and prov["full"] == [],
       f"file sources without repo_numbers: band 'Source: fsim_core/transport.py @ <hash>' only ({prov})")
    ck(build.notes_with_sources("n", prov) == "n", "band-only provenance adds no 'Sources:' block to the notes")


def check_timing(parts: dict, sample: dict) -> None:
    """Well-formedness and wiring of the injected <p:timing> blocks."""
    for i, js in enumerate(sample["slides"]):
        sid = js["id"]
        raw = parts.get(f"ppt/slides/slide{i + 1}.xml")
        if raw is None:
            ck(False, f"timing {sid}: slide part slide{i + 1}.xml exists")
            continue
        try:
            root = etree.fromstring(raw)
        except etree.XMLSyntaxError as e:
            ck(False, f"timing {sid}: slide XML is well-formed ({e})")
            continue
        shapes = {}  # id -> (element local name, cNvPr name)
        ids = []
        for c in root.iter(P + "cNvPr"):
            ids.append(c.get("id"))
            shapes[c.get("id")] = (etree.QName(c.getparent().getparent()).localname, c.get("name", ""))
        ids_ok = len(ids) == len(set(ids))
        timings = root.findall(P + "timing")
        rv_named = {sid_: nm for sid_, (_, nm) in shapes.items() if nm.startswith("rv-")}

        if not js.get("reveal"):
            ck(ids_ok and not timings and not rv_named,
               f"timing {sid}: no reveal -> no <p:timing>, no rv- shapes, unique shape ids")
            continue

        max_step = validate_sections.max_reveal_step(js)
        ck(ids_ok, f"timing {sid}: shape ids unique ({len(ids)} ids)")
        ck(len(timings) == 1 and root[-1].tag == P + "timing"
           and root[-2].tag == P + "clrMapOvr",
           f"timing {sid}: exactly one <p:timing>, last child of <p:sld>, after <p:clrMapOvr>")
        if len(timings) != 1:
            continue
        t = timings[0]
        ctn_ids = [c.get("id") for c in t.iter(P + "cTn")]
        ck(len(ctn_ids) == len(set(ctn_ids)), f"timing {sid}: {len(ctn_ids)} cTn ids unique")
        main = [c for c in t.iter(P + "cTn") if c.get("nodeType") == "mainSeq"]
        clicks = main[0].find(P + "childTnLst").findall(P + "par") if main else []
        ck(len(clicks) == max_step,
           f"timing {sid}: {len(clicks)} click groups == max reveal step {max_step}")
        spids_all = [e.get("spid") for e in t.iter(P + "spTgt")]
        ck(bool(spids_all) and all(s in shapes for s in spids_all),
           f"timing {sid}: all {len(spids_all)} spTgt spids resolve to shapes on the slide")
        group_ok = True
        for g, par in enumerate(clicks, start=1):
            targeted = {e.get("spid") for e in par.iter(P + "spTgt")}
            expected = {s for s, nm in rv_named.items() if nm.startswith(f"rv-{g}-")}
            group_ok &= targeted == expected and bool(expected)
            effects = [c for c in par.iter(P + "cTn") if c.get("presetClass") == "entr"]
            group_ok &= bool(effects) and effects[0].get("nodeType") == "clickEffect"
            group_ok &= all(e.get("nodeType") == "withEffect" for e in effects[1:])
        ck(group_ok and set().union(*[{e.get("spid") for e in p_.iter(P + "spTgt")} for p_ in clicks]) == set(rv_named),
           f"timing {sid}: click group g animates exactly the rv-g-* shapes, first is clickEffect")
        effect = js.get("reveal_effect", "fade")
        presets = {c.get("presetID") for c in t.iter(P + "cTn") if c.get("presetClass") == "entr"}
        ck(presets == {PRESET_ID[effect]},
           f"timing {sid}: presetID {sorted(presets)} matches reveal_effect '{effect}'")
        bld = [e.get("spid") for e in t.iter(P + "bldP")]
        ck(all(s in shapes and shapes[s][0] == "sp" for s in bld),
           f"timing {sid}: {len(bld)} bldP spids resolve to text shapes")


def check_no_reveal_pptx(parts: dict) -> None:
    xml_parts = {n: b for n, b in parts.items() if n.endswith(".xml")}
    ck(not any(b"p:timing" in b or b"<timing" in b for b in xml_parts.values()),
       "no-reveal pptx: no p:timing in any package part")
    ck(not any(b'name="rv-' in b for b in xml_parts.values()),
       "no-reveal pptx: no reveal-named (rv-) shapes")
    dup = False
    for n, b in xml_parts.items():
        if n.startswith("ppt/slides/slide"):
            ids = re.findall(rb'<p:cNvPr id="(\d+)"', b)
            dup |= len(ids) != len(set(ids))
    ck(not dup, "no-reveal pptx: shape ids unique on every slide")


# ---------------------------------------------------------------------------
# html checks
# ---------------------------------------------------------------------------

def section_html(html_text: str, sid: str) -> str:
    m = re.search(r'<section class="slide[^"]*" id="%s".*?</section>' % re.escape(sid), html_text, re.S)
    return m.group(0) if m else ""


def expected_reveal_attrs(slide: dict) -> list:
    """The data-reveal values render_html.py must write for `slide`, sorted."""
    r = slide.get("reveal") or {}
    vals = []
    for key in ("bullets", "equations", "columns"):
        vals += [v for v in r.get(key, []) if v > 0]
    rows = r.get("table_rows", [])
    vals += [v for v in rows if v > 0]
    if rows and rows[0] > 0 and (slide.get("table") or {}).get("header"):
        vals.append(rows[0])  # the header row enters with the first group
    frames = r.get("figure_frames")
    if frames is not None and (slide.get("figure") or {}).get("frames"):
        vals += [v for v in frames if v > 0]
        if frames[0] > 0 and slide["figure"].get("caption"):
            vals.append(frames[0])  # caption enters with the first frame
    elif r.get("figure", 0) > 0 and slide.get("figure"):
        vals.append(r["figure"])
    return sorted(vals)


def css_rules(css: str):
    """(selector, declarations) for every innermost rule, @media flattened."""
    rules, stack, buf = [], [], ""
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    for ch in css:
        if ch == "{":
            stack.append(buf.strip())
            buf = ""
        elif ch == "}":
            if stack:
                sel = stack.pop()
                if buf.strip():
                    rules.append((sel, buf.strip()))
            buf = ""
        else:
            buf += ch
    return rules


HIDING = re.compile(r"display\s*:\s*none|visibility\s*:\s*hidden|opacity\s*:\s*0(?![.\d])")
ALLOWED_HIDING_SELECTORS = {
    ".frames > img:not(:last-child)",  # scroll mode shows the complete (last) frame
    ".rv-notes",                       # the notes overlay the script creates
    ".present-btn, .rv-notes",         # print: hide the script-created controls
}


def check_html(html_text: str, sample: dict, label: str, reveal: bool) -> None:
    n = len(sample["slides"])
    n_sections = len(re.findall(r"<section", html_text))
    ck(n_sections == n, f"{label}: index.html has one <section> per slide ({n_sections} vs {n})")

    scripts = re.findall(r"<script\b([^>]*)>(.*?)</script>", html_text, re.S | re.I)
    ck(len(scripts) == 1 and "src" not in scripts[0][0].lower(),
       f"{label}: exactly one <script>, inline (present mode), no src")
    if scripts:
        js_path = HTML_SCRATCH / f"present_{label.replace(' ', '_')}.js"
        js_path.parent.mkdir(parents=True, exist_ok=True)
        js_path.write_text(scripts[0][1], encoding="utf-8")
        r = subprocess.run(["node", "--check", str(js_path)], capture_output=True, text=True)
        ck(r.returncode == 0, f"{label}: present-mode script parses (node --check)")

    ext = []
    ext += re.findall(r"<script[^>]+src\s*=", html_text, re.I)
    ext += re.findall(r"<link\b", html_text, re.I)
    ext += re.findall(r"@import", html_text, re.I)
    ext += re.findall(r"""(?:src|href)\s*=\s*["']\s*(?:https?:)?//""", html_text, re.I)
    ext += re.findall(r"""url\(\s*["']?\s*(?:https?:)?//""", html_text, re.I)
    ext += [s for s in re.findall(r"""<img[^>]*\ssrc\s*=\s*["']([^"']{0,40})""", html_text, re.I)
            if not s.startswith("data:")]
    ext += re.findall(r"<(?:iframe|object|embed)\b", html_text, re.I)
    ck(not ext, f"{label}: index.html references no external resource ({ext[:3]})")

    # ---- scripts removed: all content still shown ----
    nojs = re.sub(r"<script\b.*?</script>", "", html_text, flags=re.S | re.I)
    style = "".join(re.findall(r"<style>(.*?)</style>", nojs, re.S))
    bad = [sel for sel, decl in css_rules(style)
           if HIDING.search(decl) and "presenting" not in sel and sel not in ALLOWED_HIDING_SELECTORS]
    ck(not bad, f"{label}: no-JS: every CSS rule that hides content is scoped to body.presenting ({bad})")
    inline_hidden = re.findall(r"""<[^>]+\s(?:hidden(?:\s|>|=)|style\s*=\s*["'][^"']*(?:display\s*:\s*none|visibility\s*:\s*hidden))""", nojs, re.I)
    state = re.findall(r"""class\s*=\s*["'][^"']*\b(?:rv-hidden|presenting|current)\b""", nojs)
    ck(not inline_hidden and not state,
       f"{label}: no-JS: no element hidden inline or pre-set to a present-mode state")
    inline_free = re.sub(r'</?(?:sub|sup)>|<span class="mathinline">|</span>', "", nojs)
    body_text = re.sub(r"<[^>]+>", " ", inline_free)
    import html as _html
    body_text = " ".join(_html.unescape(body_text).split())
    missing = []
    for s in sample["slides"]:
        items = list(s.get("bullets", []))
        for col in s.get("columns", []):
            items += [col["heading"], *col["bullets"]]
        table = s.get("table") or {}
        items += list(table.get("header", [])) + [c for row in table.get("rows", []) for c in row]
        for eq in s.get("equations", []):
            items.append(eq["caption"])
        if s.get("figure", {}).get("caption"):
            items.append(s["figure"]["caption"])
        if s.get("subtitle"):
            items.append(s["subtitle"])
        missing += [x for x in items if " ".join(mathimg.inline_plain(str(x)).split()) not in body_text]
    ck(not missing, f"{label}: no-JS: every bullet/column/table/caption/subtitle text is in the page ({missing[:3]})")
    ck("$\\" not in body_text and "(undefined)" not in body_text and "(None)" not in body_text,
       f"{label}: no raw inline LaTeX and no '(undefined)' in the page text")
    for s in sample["slides"]:
        sec = section_html(html_text, s["id"])
        if s.get("subtitle") and s["layout"] not in ("quote",):
            ck('<p class="subtitle">' in sec, f"{label} slide {s['id']} ({s['layout']}): subtitle rendered")
        n_prov = len(re.findall(r'<span class="provenance"', sec))
        n_src = len(re.findall(r'<div class="notes-sources">', sec))
        want = 1 if s.get("repo_numbers") else 0
        want_prov = 1 if (s.get("repo_numbers") or build.source_files(s)) else 0
        ck(n_prov == want_prov and n_src == want,
           f"{label} slide {s['id']}: {want_prov} one-line provenance + {want} notes 'Sources:' block ({n_prov}/{n_src})")

    html_bytes = len(html_text.encode("utf-8"))
    ck(html_bytes < MAX_HTML_BYTES, f"{label}: index.html is under 16 MB (got {html_bytes} bytes)")

    if not reveal:
        ck('data-reveal="' not in html_text, f"{label}: no data-reveal attribute (reveal disabled)")
        return
    for s in sample["slides"]:
        sec = section_html(html_text, s["id"])
        got = sorted(int(v) for v in re.findall(r'data-reveal="(\d+)"', sec))
        want = expected_reveal_attrs(s)
        ck(got == want, f"{label} slide {s['id']}: data-reveal steps {got} == expected {want}")
        if s.get("reveal"):
            mx = validate_sections.max_reveal_step(s)
            ck(sorted(set(got)) == list(range(1, mx + 1)),
               f"{label} slide {s['id']}: data-reveal covers steps 1..{mx}")
        frames = (s.get("figure") or {}).get("frames")
        if s.get("reveal", {}).get("figure_frames") is not None and frames:
            stack = re.search(r'<div class="frames">(.*?)</div>', sec, re.S)
            imgs = re.findall(r"<img\b[^>]*>", stack.group(1)) if stack else []
            order = [int(x) for x in re.findall(r'data-frame="(\d+)"', stack.group(1))] if stack else []
            ck(len(imgs) == len(frames) and order == list(range(len(frames))),
               f"{label} slide {s['id']}: frame stack holds {len(frames)} frames, complete graph last")


# ---------------------------------------------------------------------------
# equation images and dark-theme contrast
# ---------------------------------------------------------------------------

def check_equation_images(sample: dict) -> None:
    """No equation PNG is cropped: the non-white content bbox must sit
    strictly inside the image with a white margin on all four sides."""
    from PIL import Image, ImageChops
    latexes = []
    for s in sample["slides"]:
        latexes += [eq["latex"] for eq in s.get("equations", [])]
        latexes += [b.strip()[1:-1] for b in s.get("bullets", []) if re.fullmatch(r"\$[^$]+\$", b.strip())]
    for latex in latexes:
        path = mathimg.cache_path(latex)
        if not path.exists():
            ck(False, f"equation image exists for {latex!r}")
            continue
        ck(mathimg.cached_version(path) == mathimg.CACHE_VERSION,
           f"equation {path.name}: written by mathimg cache version {mathimg.CACHE_VERSION}")
        with Image.open(path) as im:
            gray = im.convert("L")
            w, h = gray.size
            ink = gray.point(lambda v: 255 if v < 250 else 0)
            bbox = ink.getbbox()
        if bbox is None:
            ck(False, f"equation {path.name}: has visible content")
            continue
        left, top, right, bottom = bbox[0], bbox[1], w - bbox[2], h - bbox[3]
        m = min(left, top, right, bottom)
        ck(m >= EQ_MIN_MARGIN_PX,
           f"equation {path.name} ({latex[:30]!r}): content inside {w}x{h} canvas, "
           f"margins L{left} T{top} R{right} B{bottom} >= {EQ_MIN_MARGIN_PX} px")


def _lum(hex6: str) -> float:
    def lin(c):
        c /= 255.0
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (int(hex6[i:i + 2], 16) for i in (0, 2, 4))
    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def contrast(a: str, b: str) -> float:
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def check_dark_contrast(html_text: str, label: str) -> None:
    style = "".join(re.findall(r"<style>(.*?)</style>", html_text, re.S))
    m = re.search(r':root\[data-theme="dark"\]\s*\{[^}]*--bg:\s*#([0-9A-Fa-f]{6})', style)
    bg = m.group(1) if m else None
    ck(bg is not None, f"{label}: dark-theme --bg found ({bg})")
    if bg is None:
        return
    pairs = re.findall(r'--accent:#([0-9A-Fa-f]{6}); --accent-dark:#([0-9A-Fa-f]{6});', html_text)
    n_sections = len(re.findall(r'<section class="slide', html_text))
    ck(len(pairs) == n_sections, f"{label}: every slide sets --accent-dark ({len(pairs)}/{n_sections})")
    worst = min((contrast(d, bg), a, d) for a, d in pairs) if pairs else (0, "", "")
    ck(worst[0] >= WCAG_AA,
       f"{label}: dark-theme title colour contrast >= {WCAG_AA}:1 on #{bg} "
       f"(worst {worst[0]:.2f} for accent #{worst[1]} -> #{worst[2]})")
    rules = css_rules(style)
    for scope in (':root:not([data-theme="light"])', ':root[data-theme="dark"]'):
        covered = set()
        for sel, decl in rules:
            if "var(--accent-dark)" in decl and "color" in decl:
                for part in sel.split(","):
                    part = part.strip()
                    if part.startswith(scope):
                        covered |= {t for t in ("h2", ".big-title", ".quote") if part.endswith(t)}
        ck(covered == {"h2", ".big-title", ".quote"},
           f"{label}: dark CSS ({scope}) colours h2/.big-title/.quote with var(--accent-dark) ({sorted(covered)})")


HTML_SCRATCH = DEFAULT_OUT_DIR / "_scratch"


def main() -> int:
    global HTML_SCRATCH
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR),
                    help="self-test build directory (default out/presentation2_selftest)")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir = (ROOT / out_dir).resolve() if not out_dir.is_absolute() else out_dir.resolve()
    if out_dir == DELIVERABLE_DIR.resolve() or DELIVERABLE_DIR.resolve() in out_dir.parents:
        print(f"refusing to run against the deliverable directory {DELIVERABLE_DIR}; "
              "use --out-dir out/presentation2_selftest")
        return 2
    no_reveal_dir = out_dir / "no_reveal"
    scratch = out_dir / "_scratch"
    sections_dir = out_dir / "_sections"
    HTML_SCRATCH = scratch

    sample = json.loads(SAMPLE_JSON.read_text(encoding="utf-8"))

    check_repo_numbers_how(sample, scratch)
    check_notes_length_validation(sample, scratch)
    check_figure_size_validation(sample, scratch)
    check_reveal_validation(sample, scratch)
    check_notes_sources_validation(sample, scratch)
    check_inline_math()
    check_band_only_provenance()

    # Build fixture: the sample with a subtitle on every titled slide that has
    # none, so the subtitle band is exercised on every layout (fix round 2).
    # Every later check reads this fixture, not the bare sample.
    sample = json.loads(json.dumps(sample))
    for s in sample["slides"]:
        if not s.get("subtitle") and s["layout"] not in ("title", "section-divider", "quote"):
            s["subtitle"] = f"Fixture subtitle line for the {s['layout']} layout"
    # band-only provenance end to end: the first equation slide without
    # repo_numbers names a repository file in its sources.
    band_slide = next(s for s in sample["slides"] if s["layout"] == "equation" and not s.get("repo_numbers"))
    band_slide["sources"] = list(band_slide["sources"]) + ["fsim_core/transport.py (fixture: band-only provenance)"]
    if sections_dir.exists():
        shutil.rmtree(sections_dir)
    sections_dir.mkdir(parents=True)
    (sections_dir / SAMPLE_JSON.name).write_text(json.dumps(sample, ensure_ascii=False, indent=1), encoding="utf-8")

    rc_reveal = run_build(sections_dir, out_dir, no_reveal=False)
    ck(rc_reveal == 0, "build.py --include-sample (reveal) exits 0 (incl. its office-validate step)")
    rc_plain = run_build(sections_dir, no_reveal_dir, no_reveal=True)
    ck(rc_plain == 0, "build.py --include-sample --no-reveal exits 0 (incl. its office-validate step)")
    if rc_reveal != 0 or rc_plain != 0:
        print(f"{sum(checks)}/{len(checks)} presentation2 checks passed")
        return 1

    pptx_reveal = out_dir / "qd_physics_2h.pptx"
    pptx_plain = no_reveal_dir / "qd_physics_2h.pptx"
    office_validate(pptx_reveal, "reveal")
    office_validate(pptx_plain, "no-reveal")

    html_reveal = (out_dir / "index.html").read_text(encoding="utf-8")
    html_plain = (no_reveal_dir / "index.html").read_text(encoding="utf-8")

    check_pptx_common(Presentation(str(pptx_reveal)), sample, html_reveal, "reveal")
    check_pptx_common(Presentation(str(pptx_plain)), sample, html_plain, "no-reveal")
    check_timing(slide_parts(pptx_reveal), sample)
    check_no_reveal_pptx(slide_parts(pptx_plain))
    check_html(html_reveal, sample, "reveal html", reveal=True)
    check_html(html_plain, sample, "no-reveal html", reveal=False)
    check_equation_images(sample)
    check_dark_contrast(html_reveal, "reveal html")
    check_dark_contrast(html_plain, "no-reveal html")

    print(f"{sum(checks)}/{len(checks)} presentation2 checks passed")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PermissionError as e:
        print(f"SKIPPED (sandbox temp dir): {e}")
        print("See CLAUDE.md 'Sandbox note for workers' -- re-run outside the sandbox.")
        sys.exit(0)
