"""Acceptance checks for the qd-photon-sim / rt-edge-emitter presentation
built by scripts/make_presentation.py: out/presentation/qd_edge_sps.pptx and
out/presentation/index.html.

This is a structural/consistency check, not a physics check (the physics
numbers on the slides are read from out/rt_edge/verdict.md, out/rt_edge/
sweep.csv, and fsim_core at build time by make_presentation.py itself, and
the rest of the verify/ suite already covers the physics those files
report). What this file confirms:
  * the pptx opens and has a plausible slide count, a title on every slide,
    enough slides carrying a picture, and enough slides with speaker notes;
  * every embedded picture's bytes match a real file under
    out/presentation/figures/ or out/rt_edge/ (content-hash match, not path
    metadata, since python-pptx does not retain the original source path of
    an embedded image);
  * out/presentation/index.html mirrors the same slide count, embeds no
    external script or network reference, and stays under the 16 MB
    artifact limit;
  * the deck's results slide states the same VERDICT line as
    out/rt_edge/verdict.md, verbatim (substring match), and so does the
    HTML mirror.

CLI: python verify/verify_presentation.py
Exit code 0 iff every check passes; prints 'N/N presentation checks passed'.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "out" / "presentation"
PPTX_PATH = OUT_DIR / "qd_edge_sps.pptx"
HTML_PATH = OUT_DIR / "index.html"
FIG_DIR = OUT_DIR / "figures"
RT_DIR = ROOT / "out" / "rt_edge"
VERDICT_MD = RT_DIR / "verdict.md"

MAX_HTML_BYTES = 16 * 1024 * 1024

checks = []


def ck(ok: bool, name: str):
    checks.append(bool(ok))
    if not ok:
        print("FAIL", name)


def sha1_of(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def main() -> int:
    verdict_text = VERDICT_MD.read_text(encoding="utf-8")
    verdict_line = re.search(r"^VERDICT:.*$", verdict_text, re.MULTILINE).group(0).strip()

    allowed_hashes = set()
    for candidate_dir in (FIG_DIR, RT_DIR):
        for f in candidate_dir.glob("*.png"):
            allowed_hashes.add(sha1_of(f))

    prs = Presentation(str(PPTX_PATH))
    slides = list(prs.slides)
    n_slides = len(slides)
    ck(14 <= n_slides <= 20, f"slide count in [14, 20] (got {n_slides})")

    n_with_picture = 0
    n_with_notes = 0
    verdict_seen_pptx = False
    linewidth_seen_results_pptx = False
    anchor_seen_results_pptx = False

    for i, slide in enumerate(slides, start=1):
        title_shape = slide.shapes.title
        title_text = title_shape.text.strip() if title_shape is not None else ""
        ck(bool(title_text), f"slide {i} has a non-empty title")

        slide_text_parts = [title_text]
        has_picture = False
        for shape in slide.shapes:
            if shape.has_text_frame:
                slide_text_parts.append(shape.text_frame.text)
            if shape.shape_type == MSO_SHAPE_TYPE.PICTURE:
                has_picture = True
                blob_hash = sha1_of_bytes(shape.image.blob)
                ck(blob_hash in allowed_hashes,
                   f"slide {i} picture matches a real file under "
                   "out/presentation/figures or out/rt_edge")
        if has_picture:
            n_with_picture += 1

        if slide.has_notes_slide and slide.notes_slide.notes_text_frame.text.strip():
            n_with_notes += 1

        if verdict_line in "\n".join(slide_text_parts):
            verdict_seen_pptx = True
        if title_text.startswith("Results:") and "linewidth" in "\n".join(slide_text_parts).lower():
            linewidth_seen_results_pptx = True
        if title_text.startswith("Results:") and "6.5 meV" in "\n".join(slide_text_parts):
            anchor_seen_results_pptx = True

    ck(n_with_picture >= 8, f"at least 8 slides contain a picture (got {n_with_picture})")
    ck(n_with_notes >= 10, f"at least 10 slides have non-empty speaker notes (got {n_with_notes})")
    ck(verdict_seen_pptx, "the pptx's results slide states the verdict.md VERDICT line verbatim")
    ck(linewidth_seen_results_pptx, "the pptx's results slide includes the phrase 'linewidth'")
    ck(anchor_seen_results_pptx, "the pptx's results slide mentions the 6.5 meV anchor result")

    html = HTML_PATH.read_text(encoding="utf-8")
    n_sections = len(re.findall(r"<section", html))
    ck(n_sections == n_slides,
       f"index.html has the same number of <section> elements as slides "
       f"({n_sections} vs {n_slides})")
    ck("<script src=" not in html, "index.html has no external <script src=...>")
    ck(not re.search(r"""(?:src|href)\s*=\s*["']https?://""", html),
       "index.html has no external http(s) src/href reference")
    html_bytes = HTML_PATH.stat().st_size
    ck(html_bytes < MAX_HTML_BYTES, f"index.html is under 16 MB (got {html_bytes} bytes)")
    ck(verdict_line in html, "index.html states the verdict.md VERDICT line verbatim")
    ck("linewidth" in html.lower(), "index.html includes the phrase 'linewidth'")
    ck("6.5 meV" in html, "index.html mentions the 6.5 meV anchor result")

    print(f"{sum(checks)}/{len(checks)} presentation checks passed")
    return 0 if all(checks) else 1


def sha1_of_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()


if __name__ == "__main__":
    sys.exit(main())
