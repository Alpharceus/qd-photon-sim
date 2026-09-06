"""presentation2/render_html.py -- self-contained HTML renderer per the
contract in presentation2/SCHEMA.md. Reads every presentation2/sections/
NN_*.json (section "00" skipped unless --include-sample) and writes
out/presentation2/index.html: one <section> per slide, images embedded as
base64 data URIs, inline CSS only, no external scripts or network
references, a sticky table of contents, and a prefers-color-scheme dark
mode. MathJax-free: equations are images (from mathimg.py's cache); inline
$...$ text is shown as a monospace span, never rendered as math.

Never special-cases section content -- every layout branch below works from
the generic slide object the schema defines.

CLI: python presentation2/render_html.py [--include-sample] [--sections-dir DIR]
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EQ_DIR = HERE / "figures" / "out" / "eq"
OUT_DIR = ROOT / "out" / "presentation2"
OUT_HTML = OUT_DIR / "index.html"

ACCENT_PALETTE = [
    "1E2761", "065A82", "2C5F2D", "6D2E46",
    "028090", "B85042", "36454F", "990011",
]

INLINE_MATH_RE = re.compile(r"\$[^$]+\$")


def accent_for(section_id: str) -> str:
    try:
        n = int(section_id)
    except (TypeError, ValueError):
        n = 0
    return ACCENT_PALETTE[n % len(ACCENT_PALETTE)]


def equation_image_path(latex: str) -> Path:
    h = hashlib.sha256(latex.encode("utf-8")).hexdigest()[:16]
    return EQ_DIR / f"eq_{h}.png"


def data_uri(path: Path) -> str:
    data = path.read_bytes()
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def esc(text) -> str:
    return html.escape(str(text), quote=True)


def inline_text(text: str) -> str:
    """Escape `text` for HTML, then re-wrap any $...$ span in a monospace
    <span> -- text only, never rendered as math (MathJax-free contract)."""
    escaped = esc(text)

    def repl(m: "re.Match[str]") -> str:
        return f'<span class="mathinline">{m.group(0)}</span>'

    # esc() already escaped `$`? No -- html.escape does not touch `$`, so the
    # regex still matches the escaped string correctly.
    return INLINE_MATH_RE.sub(repl, escaped)


def load_sections(sections_dir: Path, include_sample: bool):
    sections = []
    for path in sorted(sections_dir.glob("[0-9][0-9]_*.json")):
        if path.stem[:2] == "00" and not include_sample:
            continue
        sections.append(json.loads(path.read_text(encoding="utf-8")))
    return sections


def render_bullets(bullets) -> str:
    items = "".join(f"<li>{inline_text(b)}</li>" for b in bullets)
    return f'<ul class="bullets">{items}</ul>'


def render_figure(figure) -> str:
    img_path = ROOT / figure["path"]
    if not img_path.exists():
        return ""
    src = data_uri(img_path)
    caption = f'<figcaption>{esc(figure.get("caption", ""))}</figcaption>' if figure.get("caption") else ""
    return f'<figure class="fig"><img src="{src}" alt="{esc(figure.get("caption", ""))}">{caption}</figure>'


def render_equation_block(eq) -> str:
    img_path = equation_image_path(eq["latex"])
    body = ""
    if img_path.exists():
        src = data_uri(img_path)
        body = f'<img src="{src}" alt="{esc(eq["latex"])}">'
    caption = f'<figcaption>{esc(eq.get("caption", ""))}</figcaption>' if eq.get("caption") else ""
    return f'<figure class="eq">{body}{caption}</figure>'


def render_table(table) -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in table.get("header", []))
    rows = "".join(
        "<tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>"
        for row in table.get("rows", [])
    )
    return (
        '<div class="table-wrap"><table><thead><tr>'
        f"{head}</tr></thead><tbody>{rows}</tbody></table></div>"
    )


def render_provenance(repo_numbers) -> str:
    if not repo_numbers:
        return ""
    parts = [f"{k} = {v.get('value')} ({esc(v.get('file'))})" for k, v in repo_numbers.items()]
    return f'<p class="provenance">Source: {" | ".join(parts)}</p>'


def render_sources(sources) -> str:
    if not sources:
        return ""
    items = "".join(f"<li>{esc(s)}</li>" for s in sources)
    return f'<ul class="sources">{items}</ul>'


def render_equations(equations) -> str:
    return "".join(render_equation_block(eq) for eq in equations)


def render_equations_extra(equations) -> str:
    """The generic "equations render on every layout" block appended below
    a layout's own content when that layout does not already make the
    equations array its primary content (equation, equation+figure)."""
    if not equations:
        return ""
    return f'<div class="equations equations-extra">{render_equations(equations)}</div>'


def render_slide(slide, section, accent) -> str:
    layout = slide["layout"]
    title = esc(slide.get("title", ""))
    subtitle = f'<p class="subtitle">{esc(slide["subtitle"])}</p>' if slide.get("subtitle") else ""
    body = ""

    if layout == "title":
        body = f'<h1 class="big-title">{title}</h1>{subtitle}'
    elif layout == "section-divider":
        body = f'<h1 class="big-title">{title}</h1>{subtitle}'
    elif layout == "quote":
        body = f'<blockquote class="quote">"{title}"</blockquote>{subtitle}'
    elif layout == "bullets":
        body = f"<h2>{title}</h2>{subtitle}"
        if slide.get("bullets"):
            body += render_bullets(slide["bullets"])
    elif layout == "bullets+figure":
        body = f'<h2>{title}</h2><div class="split">'
        body += f'<div class="split-text">{render_bullets(slide.get("bullets", []))}</div>'
        if slide.get("figure"):
            body += f'<div class="split-fig">{render_figure(slide["figure"])}</div>'
        body += "</div>"
    elif layout == "figure":
        body = f"<h2>{title}</h2>"
        if slide.get("figure"):
            body += render_figure(slide["figure"])
    elif layout == "equation":
        body = f'<h2>{title}</h2><div class="equations">{render_equations(slide.get("equations", []))}</div>'
    elif layout == "equation+figure":
        body = f'<h2>{title}</h2><div class="split">'
        body += f'<div class="split-text equations">{render_equations(slide.get("equations", []))}</div>'
        if slide.get("figure"):
            body += f'<div class="split-fig">{render_figure(slide["figure"])}</div>'
        body += "</div>"
    elif layout == "two-column":
        body = f'<h2>{title}</h2><div class="columns">'
        for col in slide.get("columns", []):
            body += (
                '<div class="column">'
                f'<h3>{esc(col.get("heading", ""))}</h3>'
                f'{render_bullets(col.get("bullets", []))}'
                "</div>"
            )
        body += "</div>"
    elif layout == "table":
        body = f'<h2>{title}</h2>{render_table(slide.get("table", {}))}'
    elif layout == "bullets+table":
        body = f"<h2>{title}</h2>{subtitle}"
        if slide.get("bullets"):
            body += render_bullets(slide["bullets"])
        body += render_table(slide.get("table", {}))
    else:
        body = f"<h2>{title}</h2>"

    # Equations are rendered on EVERY layout that carries an equations
    # array (SCHEMA.md), not only "equation" and "equation+figure", which
    # already made the array their primary content above.
    if layout not in ("equation", "equation+figure"):
        body += render_equations_extra(slide.get("equations"))

    dark = layout == "section-divider"
    slide_style = f'style="--accent:#{accent};"' + (' data-dark="1"' if dark else "")

    notes = esc(slide.get("notes", ""))
    sources = render_sources(slide.get("sources"))
    provenance = render_provenance(slide.get("repo_numbers"))

    return f"""
<section class="slide layout-{esc(layout)}" id="{esc(slide['id'])}" {slide_style}>
  <div class="slide-body">{body}</div>
  <footer class="slide-footer">
    <span class="section-name">{esc(section.get('title', ''))}</span>
    {provenance}
    <span class="slide-id">{esc(slide['id'])}</span>
  </footer>
  {sources}
  <details class="notes"><summary>Speaker notes</summary><p>{notes}</p></details>
</section>
""".strip()


def render_toc(sections) -> str:
    items = []
    for section in sections:
        accent = accent_for(section.get("section", "0"))
        slide_links = "".join(
            f'<a href="#{esc(s["id"])}">{esc(s["id"])}: {esc(s.get("title", ""))}</a>'
            for s in section.get("slides", [])
        )
        items.append(
            f'<div class="toc-section" style="--accent:#{accent};">'
            f'<div class="toc-section-title">{esc(section.get("title", ""))}</div>'
            f'<nav class="toc-slides">{slide_links}</nav></div>'
        )
    return f'<nav class="toc">{"".join(items)}</nav>'


CSS = """
:root {
  --bg: #FFFFFF; --ink: #212121; --muted: #6B6B6B; --panel: #F2F2F2; --border: #DDDDDD;
}
@media (prefers-color-scheme: dark) {
  :root { --bg: #16181D; --ink: #EDEDED; --muted: #A0A4AE; --panel: #21242B; --border: #34384040; }
}
* { box-sizing: border-box; }
html, body { margin: 0; padding: 0; background: var(--bg); color: var(--ink);
  font-family: Calibri, Arial, sans-serif; overflow-x: hidden; }
.layout { display: flex; align-items: flex-start; }
.toc { position: sticky; top: 0; align-self: flex-start; width: 260px; max-height: 100vh;
  overflow-y: auto; padding: 1rem; border-right: 1px solid var(--border); flex: 0 0 260px; }
.toc-section { margin-bottom: 1rem; border-left: 4px solid var(--accent); padding-left: 0.5rem; }
.toc-section-title { font-weight: bold; margin-bottom: 0.25rem; }
.toc-slides { display: flex; flex-direction: column; }
.toc-slides a { color: var(--muted); text-decoration: none; font-size: 0.8rem; padding: 0.1rem 0; }
.toc-slides a:hover { color: var(--ink); }
.deck { flex: 1 1 auto; min-width: 0; padding: 1rem 2rem; }
.slide { max-width: 960px; margin: 0 auto 3rem auto; padding: 1.5rem 2rem; border-radius: 8px;
  border: 1px solid var(--border); background: var(--bg); }
.slide[data-dark="1"] { background: var(--accent); color: #FFFFFF; }
.slide[data-dark="1"] .slide-footer, .slide[data-dark="1"] .subtitle { color: #DDDDDD; }
.slide h1, .slide h2, .slide h3 { font-family: Cambria, Georgia, serif; }
.slide h2 { color: var(--accent); }
.slide[data-dark="1"] h1 { color: #FFFFFF; }
.big-title { font-size: 2.2rem; text-align: center; color: var(--accent); }
.slide[data-dark="1"] .big-title { color: #FFFFFF; }
.subtitle { text-align: center; font-style: italic; color: var(--muted); }
.quote { font-size: 1.6rem; font-style: italic; font-weight: bold; color: var(--accent);
  text-align: center; border: none; margin: 1rem 0; }
.bullets { line-height: 1.6; }
.bullets li { margin-bottom: 0.4rem; }
.split { display: flex; gap: 1.5rem; }
.split-text, .split-fig { flex: 1 1 0; min-width: 0; }
.fig, .eq { text-align: center; margin: 0.5rem 0; }
.fig img, .eq img { max-width: 100%; height: auto; }
.fig figcaption, .eq figcaption { font-size: 0.85rem; font-style: italic; color: var(--muted); margin-top: 0.3rem; }
.equations { display: flex; flex-direction: column; gap: 1rem; }
.equations-extra { margin-top: 1rem; padding-top: 1rem; border-top: 1px solid var(--border); }
.columns { display: flex; gap: 1.5rem; }
.column { flex: 1 1 0; min-width: 0; }
.table-wrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid var(--border); padding: 0.4rem 0.6rem; text-align: left; font-size: 0.9rem; }
th { background: var(--accent); color: #FFFFFF; }
.mathinline { font-family: "Courier New", Consolas, monospace; background: var(--panel); padding: 0 0.2rem; }
.slide-footer { display: flex; justify-content: space-between; gap: 1rem; font-size: 0.75rem;
  color: var(--muted); border-top: 1px solid var(--border); margin-top: 1rem; padding-top: 0.4rem; }
.provenance { font-style: italic; }
.sources { font-size: 0.75rem; color: var(--muted); margin: 0.5rem 0 0 0; padding-left: 1.2rem; }
.notes { font-size: 0.8rem; color: var(--muted); margin-top: 0.5rem; }
.notes summary { cursor: pointer; }
@media (max-width: 800px) {
  .layout { flex-direction: column; }
  .toc { position: static; width: auto; max-height: none; border-right: none; border-bottom: 1px solid var(--border); }
  .split, .columns { flex-direction: column; }
}
"""


def build_html(sections) -> str:
    toc = render_toc(sections)
    slides_html = []
    for section in sections:
        accent = accent_for(section.get("section", "0"))
        for slide in section.get("slides", []):
            slides_html.append(render_slide(slide, section, accent))
    body = "\n".join(slides_html)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>How quantum dots emit single photons</title>
<style>{CSS}</style>
</head>
<body>
<div class="layout">
{toc}
<main class="deck">
{body}
</main>
</div>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-sample", action="store_true")
    parser.add_argument("--sections-dir", default=None)
    args = parser.parse_args()

    sections_dir = Path(args.sections_dir) if args.sections_dir else HERE / "sections"
    sections = load_sections(sections_dir, args.include_sample)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    html_text = build_html(sections)
    OUT_HTML.write_text(html_text, encoding="utf-8")
    size_mb = OUT_HTML.stat().st_size / (1024 * 1024)
    n_slides = sum(len(s.get("slides", [])) for s in sections)
    print(f"wrote {OUT_HTML} ({n_slides} slides, {size_mb:.2f} MB)")
    if size_mb > 16:
        print("ERROR: index.html exceeds the 16 MB artifact limit", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
