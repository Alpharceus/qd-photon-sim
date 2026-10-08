"""Build a self-contained HTML graphs page for the nitride nanowire sweep.

Reads out/nitride_nanowire/full/ (results.md, manifest.json, every *.png) plus
the Q1-Q7 sources and the strain-anchor review note, and writes ONE file,
<out-dir>/graphs.html, with every PNG embedded as a base64 data URI. The page
makes no external requests (no http(s) src/href, no web fonts, no scripts).

Idempotent: the output depends only on the input files and the git HEAD; no
timestamps are written. Nothing other than graphs.html is written.

Usage:
    python scripts/build_nanowire_graphs_page.py [--out-dir out/nitride_nanowire/full]

This script reports numbers; it computes none except the plain difference
predicted - measured (nm) in the strain-anchor table, a [DR] reduction of the
quoted values.
"""
from __future__ import annotations

import argparse
import base64
import html
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUT = Path("out") / "nitride_nanowire" / "full"
BRIEF = REPO / ".workers" / "briefs" / "simulator-physics-audit.md"
COHERENCE = REPO / ".workers" / "review" / "nitride-nanowire-physics-opus-coherence-findings.md"
STRAIN_REVIEW = REPO / ".workers" / "review" / "audit-nitride-strain-mass-anchors.md"
MEMORY_DIR = Path.home() / ".claude" / "projects" / "C--Users-rama4-Downloads-quantum-dot" / "memory"
MEMORY_FILES = [MEMORY_DIR / "nanowire-next-step.md", MEMORY_DIR / "qd-photon-sim-handoff.md"]

SHARED_QUAL = " Qualification (applies to every figure this run):"

VERDICT_COLS = [
    "family", "regime", "strain_bound", "rep_rate_hz", "idealized_status",
    "eligible", "paired_optical_pass", "quality_pass", "hardware_qualified",
    "rti_qualified", "coverage", "invalid",
]
BEST_COLS = [
    "family", "value", "row_id", "core_radius_nm", "height_nm", "x_in", "T_hs",
    "rep_rate_hz", "strain_bound", "regime", "delivered_flux",
    "photons_per_cycle_delivered", "quality_pass", "hardware_qualified", "rti_qualified",
]

E = html.escape


def read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return None


def kv_tokens(line: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in re.finditer(r"([A-Za-z_][A-Za-z0-9_]*)=(\S+)", line)}


def git_head() -> str | None:
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                           capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    out = r.stdout.strip()
    return out if r.returncode == 0 and out else None


def md_sections(md: str) -> dict[str, str]:
    """Map '## heading' -> body text (up to the next '## ' heading)."""
    sections: dict[str, str] = {}
    cur = None
    buf: list[str] = []
    for line in md.splitlines():
        if line.startswith("## "):
            if cur is not None:
                sections[cur] = "\n".join(buf).strip()
            cur, buf = line[3:].strip(), []
        elif cur is not None:
            buf.append(line)
    if cur is not None:
        sections[cur] = "\n".join(buf).strip()
    return sections


def table(headers: list[str], rows: list[list[str]], cls: str = "") -> str:
    h = "".join(f"<th>{E(c)}</th>" for c in headers)
    body = "".join("<tr>" + "".join(f"<td>{E(c)}</td>" for c in r) + "</tr>" for r in rows)
    return (f'<div class="tw"><table class="{cls}"><thead><tr>{h}</tr></thead>'
            f"<tbody>{body}</tbody></table></div>")


def parse_md_table(lines: list[str]) -> tuple[list[str], list[list[str]]]:
    rows = [[c.strip() for c in ln.strip().strip("|").split("|")] for ln in lines]
    rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
    return (rows[0], rows[1:]) if rows else ([], [])


# ---------------------------------------------------------------- sections

def section_headline(md: str) -> str:
    verdicts = [ln[len("VERDICT: "):] for ln in md.splitlines() if ln.startswith("VERDICT: ")]
    bests = [ln for ln in md.splitlines() if ln.startswith("BEST_PASSING_FLUX")]
    vrows = [[kv_tokens(v).get(c, "") for c in VERDICT_COLS] for v in verdicts]
    brows = []
    for b in bests:
        kv = kv_tokens(b)
        brows.append([b.split(" ", 1)[0]] + [kv.get(c, "") for c in BEST_COLS])
    out = ["<section id='headline'><h2>Headline verdicts</h2>",
           f"<p class='note'>{len(verdicts)} VERDICT lines and {len(bests)} BEST rows, "
           "parsed verbatim from results.md. Full raw lines are in the collapsible blocks.</p>",
           table(VERDICT_COLS, vrows, "dense"),
           "<details><summary>Raw VERDICT lines</summary><pre>"
           + E("\n".join("VERDICT: " + v for v in verdicts)) + "</pre></details>",
           "<h3>Best passing flux</h3>",
           table(["line"] + BEST_COLS, brows, "dense"),
           "<details><summary>Raw BEST lines</summary><pre>" + E("\n".join(bests))
           + "</pre></details></section>"]
    return "\n".join(out)


def section_figures(out_dir: Path, manifest: dict) -> tuple[str, int]:
    mapping = manifest.get("plot_row_mapping", {}) if isinstance(manifest, dict) else {}
    pngs = sorted(out_dir.glob("*.png"), key=lambda p: p.name)
    shared = None
    parts = ["<section id='figures'><h2>Figures</h2>"]
    for p in pngs:
        cap = ""
        entry = mapping.get(p.name)
        if isinstance(entry, dict):
            cap = str(entry.get("__caption__", ""))
        if SHARED_QUAL in cap:
            cap, rest = cap.split(SHARED_QUAL, 1)
            shared = rest.strip()
        data = base64.b64encode(p.read_bytes()).decode("ascii")
        parts.append(
            f"<figure><img alt='{E(p.stem)}' src='data:image/png;base64,{data}'>"
            f"<figcaption><b>{E(p.name)}</b>. {E(cap.strip()) or 'No caption in manifest.json.'}"
            "</figcaption></figure>")
    if shared:
        parts.insert(1, f"<p class='note'><b>Applies to every figure:</b> {E(shared)}</p>")
    parts.append("</section>")
    return "\n".join(parts), len(pngs)


def section_questions() -> str:
    brief = read_text(BRIEF) or ""
    coh = read_text(COHERENCE) or ""
    answers: dict[str, tuple[str, str]] = {}
    for m in re.finditer(r"^- Nanowire (Q[1-7])(?: \(([^)]*)\))?: (.+)$", brief, re.M):
        answers[m.group(1)] = (m.group(2) or "", m.group(3).strip())
    questions: dict[str, tuple[str, str]] = {}
    for m in re.finditer(r"^(Q[1-7]) (.*) \[(.*)\]\.?\s*$", coh, re.M):
        questions[m.group(1)] = (m.group(2).strip(), m.group(3).strip())
    mem_notes = []
    for mf in MEMORY_FILES:
        txt = read_text(mf)
        if txt is None:
            mem_notes.append(f"{mf.name}: not readable")
        elif re.search(r"answers? to (the open questions )?Q1-Q7", txt, re.I):
            mem_notes.append(f"{mf.name}: lists Q1-Q7 as awaiting user answers at the round-3 close (2026-09-22)")
        else:
            mem_notes.append(f"{mf.name}: read, no Q1-Q7 status found")
    rows = []
    for i in range(1, 8):
        q = f"Q{i}"
        qtext, default = questions.get(q, ("question text not found", ""))
        if q in answers:
            topic, ans = answers[q]
            status = "answered (user)"
            detail = (f"{topic}: " if topic else "") + ans
        else:
            status = "open"
            detail = (f"orchestrator default applied, not a user answer: {default}" if default
                      else "no user answer found")
        rows.append([q, qtext, status, detail])
    src = [f"user decisions: {BRIEF.relative_to(REPO).as_posix()}"
           + ("" if brief else " (not readable)"),
           f"question text and defaults: {COHERENCE.relative_to(REPO).as_posix()}"
           + ("" if coh else " (not readable)")] + mem_notes
    return ("<section id='questions'><h2>Open questions Q1-Q7</h2>"
            "<p class='note'>A status is 'answered' only when the audit brief's 'User decisions "
            "already taken' records a user answer; every other question is 'open'.</p>"
            + table(["Q", "question", "status", "answer or default"], rows)
            + "<p class='src'>Sources: " + E("; ".join(src)) + "</p></section>")


def section_strain(md: str) -> str:
    rev = read_text(STRAIN_REVIEW)
    rows_out = []
    heads: list[str] = []
    if rev:
        lines = rev.splitlines()
        try:
            start = next(i for i, ln in enumerate(lines) if ln.startswith("### Held-out nanowire anchors"))
        except StopIteration:
            start = None
        if start is not None:
            tbl = []
            for ln in lines[start + 1:]:
                if ln.startswith("|"):
                    tbl.append(ln)
                elif tbl:
                    break
            heads, rows_out = parse_md_table(tbl)
    cur = {}
    m = re.search(r"2013 relaxed predicted lambda_nm=([0-9.]+)", md)
    if m:
        cur[("2013", "relaxed")] = float(m.group(1))
    m = re.search(r"2014 relaxed predicted lambda_nm=([0-9.]+) / unrelaxed predicted lambda_nm=([0-9.]+)", md)
    if m:
        cur[("2014", "relaxed")] = float(m.group(1))
        cur[("2014", "unrelaxed")] = float(m.group(2))
    headers = ["anchor", "endpoint", "before Phase B E_X (eV) / nm", "after Phase B E_X (eV) / nm",
               "current results.md nm", "current - measured (nm)", "measured"]
    rows = []
    for r in rows_out:
        if len(r) < 5:
            continue
        anchor, endpoint, before, after, measured = r[:5]
        year = "2013" if "2013" in anchor else "2014" if "2014" in anchor else ""
        ep = "unrelaxed" if "unrelaxed" in endpoint else "relaxed"
        c = cur.get((year, ep))
        mm = re.search(r"([0-9]+(?:\.[0-9]+)?)", measured)
        delta = f"{c - float(mm.group(1)):+.1f}" if (c is not None and mm) else "n/a"
        rows.append([anchor, endpoint, before, after,
                     f"{c:.2f}" if c is not None else "not reported in results.md", delta, measured])
    if not rows:
        body = "<p>Strain-anchor review table not found; nothing to show.</p>"
    else:
        body = table(headers, rows)
    return ("<section id='strain'><h2>Strain anchors: Deshpande 2013 replay and 2014 comparison</h2>"
            "<p class='note'>Before and after Phase B (Yan 2014 anisotropic deformation potentials and "
            "Rinke 2008 electron-mass axis swap) are the E_X values from the review note; the current column "
            "is parsed from results.md's 'Opposite-endpoint anchor match' section. The two anchors are matched "
            "by opposite strain endpoints and are never averaged. 'current - measured' is a plain difference [DR].</p>"
            + body
            + f"<p class='src'>Sources: {E(STRAIN_REVIEW.relative_to(REPO).as_posix())}; results.md</p></section>")


CSS = """
:root{--bg:#fbfaf7;--fg:#1d1d1b;--mut:#5d5c57;--line:#dcd9d0;--card:#ffffff;--acc:#2b5d8a;--th:#f0eee8}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15161a;--fg:#e6e4de;--mut:#a09e97;--line:#34353b;--card:#1d1e23;--acc:#8fb8e0;--th:#24252b}}
:root[data-theme="dark"]{--bg:#15161a;--fg:#e6e4de;--mut:#a09e97;--line:#34353b;--card:#1d1e23;--acc:#8fb8e0;--th:#24252b}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
.eyebrow{font:12px/1.4 ui-monospace,Consolas,monospace;color:var(--mut);letter-spacing:.04em;text-transform:uppercase}
h1{font-size:1.7rem;margin:.3em 0 .4em}
h2{font-size:1.25rem;margin:2.2em 0 .6em;padding-top:.6em;border-top:1px solid var(--line)}
h3{font-size:1.05rem;margin:1.4em 0 .5em}
.note,.src{color:var(--mut);font-size:.9rem}
.tw{overflow-x:auto;border:1px solid var(--line);border-radius:6px;background:var(--card);margin:.6em 0}
table{border-collapse:collapse;width:100%;font-size:.85rem}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:var(--th);position:sticky;top:0}
table.dense td,table.dense th{white-space:nowrap;font-family:ui-monospace,Consolas,monospace;font-size:.78rem}
figure{margin:1.4em 0;background:var(--card);border:1px solid var(--line);border-radius:6px;padding:10px}
figure img{display:block;max-width:100%;height:auto;margin:0 auto;background:#fff;border-radius:4px}
figcaption{font-size:.88rem;color:var(--mut);margin-top:8px}
pre{white-space:pre-wrap;word-break:break-word;font-size:.75rem;background:var(--card);border:1px solid var(--line);padding:8px;border-radius:6px}
details summary{cursor:pointer;color:var(--acc);font-size:.9rem}
nav a{color:var(--acc);margin-right:14px;font-size:.9rem}
"""


def build(out_dir: Path) -> Path:
    md = read_text(out_dir / "results.md")
    if md is None:
        raise SystemExit(f"results.md not found in {out_dir}")
    try:
        manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        manifest = {}
    head = git_head()
    gen = manifest.get("report_only_generation_commit") if isinstance(manifest, dict) else None
    dirty = manifest.get("report_only_generation_worktree_dirty") if isinstance(manifest, dict) else None
    commit = head or (gen[:7] if gen else "unknown")
    eyebrow = f"qd-photon-sim / nitride nanowire / commit {commit}"
    if gen:
        eyebrow += f" / results generated at {gen[:7]}" + (" (dirty worktree)" if dirty else "")
    intro = md_sections("## _intro\n" + md.split("\n## ", 1)[0]).get("_intro", "")
    intro = "\n".join(ln for ln in intro.splitlines() if not ln.startswith("# ")).strip()
    fig_html, n_png = section_figures(out_dir, manifest)
    page = "\n".join([
        "<!DOCTYPE html>",
        "<html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Nanowire Sweep Graphs</title>",
        f"<style>{CSS}</style></head><body><main>",
        f"<div class='eyebrow'>{E(eyebrow)}</div>",
        "<h1>Nitride nanowire sweep: graphs and verdicts</h1>",
        f"<p class='note'>{E(intro)}</p>",
        "<nav><a href='#headline'>Verdicts</a><a href='#figures'>Figures</a>"
        "<a href='#questions'>Q1-Q7</a><a href='#strain'>Strain anchors</a></nav>",
        section_headline(md),
        fig_html,
        section_questions(),
        section_strain(md),
        f"<p class='src'>Built by scripts/build_nanowire_graphs_page.py from {E(out_dir.name)}/ "
        f"({n_png} PNGs embedded).</p>",
        "</main></body></html>",
        "",
    ])
    target = out_dir / "graphs.html"
    data = page.encode("utf-8")
    if not target.exists() or target.read_bytes() != data:
        target.write_bytes(data)
    return target


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    args = ap.parse_args(argv)
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO / out_dir
    target = build(out_dir)
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
