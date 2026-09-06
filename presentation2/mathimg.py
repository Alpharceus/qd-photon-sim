"""presentation2/mathimg.py -- render a LaTeX string to a PNG using matplotlib
mathtext (usetex=False; no TeX install required). Used by build.py to turn
every slide's `equations[].latex` into an image that render_pptx.js and
render_html.py both place, and directly by render_html.py (same process).

Cached by content hash under presentation2/figures/out/eq/eq_<hash>.png, so a
rebuild only re-renders equations that changed. render_pptx.js (Node) needs
to find the same cached file without re-running this module, so it mirrors
the hash below (sha256 of the raw latex string, hex, first 16 chars) --
keep the two in sync if this changes.

CLI: python presentation2/mathimg.py "<latex>"   -> prints the output path.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CACHE_DIR = Path(__file__).resolve().parent / "figures" / "out" / "eq"
FONT_SIZE = 40
DPI = 300
PAD_PX = 16


def _cache_key(latex: str) -> str:
    return hashlib.sha256(latex.encode("utf-8")).hexdigest()[:16]


def cache_path(latex: str) -> Path:
    return CACHE_DIR / f"eq_{_cache_key(latex)}.png"


def render_equation(latex: str, force: bool = False) -> Path:
    """Render `latex` (no surrounding $ needed) to a white-background PNG
    at 300 dpi, tightly cropped to the glyphs plus a small pad. Returns the
    cached path; skips rendering if it already exists (unless force=True)."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = cache_path(latex)
    if out_path.exists() and not force:
        return out_path

    text = f"${latex}$"

    # Measure the rendered text first (a throwaway figure), then build a
    # figure sized exactly to it -- avoids acres of white margin or clipping.
    probe = plt.figure()
    t = probe.text(0, 0, text, fontsize=FONT_SIZE)
    probe.canvas.draw()
    bbox = t.get_window_extent(renderer=probe.canvas.get_renderer())
    plt.close(probe)

    w_in = (bbox.width + 2 * PAD_PX) / DPI
    h_in = (bbox.height + 2 * PAD_PX) / DPI

    fig = plt.figure(figsize=(w_in, h_in), dpi=DPI)
    fig.patch.set_facecolor("white")
    fig.text(0.5, 0.5, text, fontsize=FONT_SIZE, color="black", ha="center", va="center")
    fig.savefig(out_path, dpi=DPI, facecolor="white")
    plt.close(fig)
    return out_path


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: mathimg.py <latex>", file=sys.stderr)
        sys.exit(1)
    print(render_equation(sys.argv[1]))


if __name__ == "__main__":
    main()
