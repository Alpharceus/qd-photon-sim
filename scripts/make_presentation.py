"""Builds the non-technical qd-photon-sim / rt-edge-emitter presentation from
real repository outputs (docs/rt_edge_contract.md's acceptance sweep,
out/rt_edge/*) plus a handful of explanatory diagrams this script draws itself
with matplotlib, using numbers read straight out of fsim_core at build time.

WHY THIS FILE IS SHAPED THIS WAY. Every number that ends up on a slide has to
trace to something the repository actually produced: the acceptance sweep
(out/rt_edge/verdict.md, out/rt_edge/sweep.csv), the evidence ledger
(out/rt_edge/evidence.json / verify/data/rt_edge_anchors.yaml), or a live call
into fsim_core (dot_levels, linewidth, transport, waveguide, materials). None
of the physics numbers below are typed in by hand; they are parsed out of the
markdown/CSV/JSON artifacts or computed fresh from the modules. The only
hand-written text is narrative framing (why single photons matter, why 300 K
matters) and the citation strings, which are themselves copied verbatim from
the repository's own evidence ledger, not invented.

OUTPUTS.
  out/presentation/qd_edge_sps.pptx  -- 16:9 python-pptx deck, one idea per
      slide, large type, short bullets, full explanation + citations in the
      speaker notes.
  out/presentation/index.html        -- the same slide sequence as scrollable
      sections, same figures embedded as base64 data URIs, no external
      scripts, light theme with a prefers-color-scheme dark variant.
  out/presentation/figures/*.png     -- the 8 matplotlib figures this script
      draws (band diagram, p-i-n cartoon, exciton/biexciton overlap, ridge
      waveguide mode, g2(0) dip, headline-model g2 vs linewidth, collected
      flux vs heat-sink temperature, finite-pulse vs static (non-headline)
      g2); the only other figure used is the existing repository file
      out/rt_edge/gui-smoke.png.

Short-deck refresh (2026-09-24): 19 slides carrying the current VERDICT line,
flux clearing the floor at 230 K,
pulsed g2 ~0.98 from re-excitation within the pump pulse, no Purcell
recovery from a planar cavity (+1.3 %), the post-H4 beta, the "static
(non-headline)" label on every card-level (drive.finite_pulse: false)
number, and the commit hash on the verdict slides. Headline numbers are
taken only from the verdict's own headline-model rows of sweep.csv.

CLI: python scripts/make_presentation.py [--gui | --no-gui]
  Default (and --no-gui): use the committed out/rt_edge/gui-smoke.png as-is,
  so the script writes nothing outside out/presentation/. --gui first
  regenerates that screenshot (fsim_gui/designer.py --frames 10 --screenshot
  ...). Outputs are otherwise rebuilt deterministically: no timestamps or
  randomness are introduced by this script itself (the deck does quote the
  read-only git commit hashes of verdict.md and HEAD).
"""
from __future__ import annotations

import argparse
import base64
import csv
import html as html_lib
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from pptx import Presentation
from pptx.util import Emu, Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import dot_levels as DL
from fsim_core import linewidth as LW
from fsim_core import materials as MAT
from fsim_core import transport as TR
from fsim_core import waveguide as WG

RT_DIR = ROOT / "out" / "rt_edge"
VERDICT_MD = RT_DIR / "verdict.md"
SWEEP_CSV = RT_DIR / "sweep.csv"
EVIDENCE_JSON = RT_DIR / "evidence.json"
GUI_PNG = RT_DIR / "gui-smoke.png"
ENVELOPE_PNG = RT_DIR / "envelope.png"
ANCHORS_YAML = ROOT / "verify" / "data" / "rt_edge_anchors.yaml"
GAASP_CARD = ROOT / "cards" / "edge-inp-gaasp-design.yaml"

OUT_DIR = ROOT / "out" / "presentation"
FIG_DIR = OUT_DIR / "figures"
PPTX_PATH = OUT_DIR / "qd_edge_sps.pptx"
HTML_PATH = OUT_DIR / "index.html"

SLIDE_W_IN = 13.333
SLIDE_H_IN = 7.5
FIG_DPI = 150


# ===================================================================
# Reading real repository outputs
# ===================================================================

def split_markdown_sections(md: str) -> dict:
    """Split on '## ' headers; '_top' holds everything above the first one
    (the title line and the fenced VERDICT/CARD blocks)."""
    parts = re.split(r"\n## ", md)
    sections = {"_top": parts[0]}
    for part in parts[1:]:
        header, _, body = part.partition("\n")
        sections[header.strip()] = body
    return sections


def parse_card_line(line: str) -> dict:
    """'CARD: edge-inp-gaasp-design role=primary g2_pulsed_min=0.4117 ...'
    -> {'card_id': 'edge-inp-gaasp-design', 'role': 'primary',
        'g2_pulsed_min': 0.4117, ...} (numeric fields parsed as float)."""
    tokens = line.split()
    kv = {"card_id": tokens[1]}
    for tok in tokens[2:]:
        key, _, value = tok.partition("=")
        try:
            kv[key] = float(value)
        except ValueError:
            kv[key] = value
    return kv


def facet_model_note(text: str) -> str:
    """Pull the 'Independent check on eta_facet' forward/back-solved line out
    of a verdict.md's "Best diagnostic-g2 row and brightness decomposition"
    section. pkg5-fix2, item 6: that paragraph no longer contains the
    literal substring "facet_factor" (it now reads e.g. "forward = 0.810767
    via ray-series-midpoint, back-solved = ..."), so match on "forward ="
    and "via " instead."""
    sections = split_markdown_sections(text)
    facet_section = sections.get("Best diagnostic-g2 row and brightness decomposition", "")
    return next((line.strip() for line in facet_section.splitlines()
                if "forward =" in line and "via " in line), "")


def parse_verdict_md(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    sections = split_markdown_sections(text)
    verdict_line = re.search(r"^VERDICT:.*$", text, re.MULTILINE).group(0).strip()
    card_lines = [m.strip() for m in re.findall(r"^CARD:.*$", text, re.MULTILINE)]
    cards = {c["card_id"]: c for c in (parse_card_line(l) for l in card_lines)}
    ceiling = sections["Literature ceiling"].strip()
    evidence_body = sections["Evidence gate"]
    missing_claims = re.findall(r"^\s*-\s*`([^`]+)`:\s*(.+)$", evidence_body, re.MULTILINE)
    fail_body = sections["Fail reasons"]
    fail_reasons = re.findall(r"^-\s*`([^`]+)`", fail_body, re.MULTILINE)
    generated_m = re.search(r"^Generated ([^;]+);", text, re.MULTILINE)
    generated = generated_m.group(1).strip() if generated_m else ""
    metrics = {}
    for token in verdict_line.removeprefix("VERDICT: ").split():
        key, sep, value = token.partition("=")
        if sep:
            metrics[key] = value
    facet_note = facet_model_note(text)
    anchor_section = sections.get(
        "Verified 6.5 meV anchor (Chatzarakis et al., Phys. Rev. Applied 20, 034011, 2023)", "")
    anchor_lines = [line.strip() for line in anchor_section.splitlines()
                    if "at gamma300 = 6.5 meV" in line]
    # P10 short-deck refresh: the per-temperature table, the flux floor, the
    # model-sensitivity table and the convention-matched Reischle comparison,
    # all parsed from the same committed verdict.md (never typed in here).
    per_T = {}
    for m in re.finditer(r"^\| (\d+) \| (\d+)/(\d+) \| (\d+)/\d+ \| [^|]+ \| ([^|]+) \| ([^|]+) \|",
                         sections.get("Per-temperature acceptance", ""), re.MULTILINE):
        per_T[int(m.group(1))] = dict(eligible=int(m.group(2)), total=int(m.group(3)),
                                      passes=int(m.group(4)), g2_min=m.group(5).strip(),
                                      flux_max=m.group(6).strip())
    floor_m = re.search(r"eligibility floor \[A\]: ([0-9.eE+]+) photons/s",
                        sections.get("Coverage", ""))
    model_sens = {}
    for m in re.finditer(r"^\| (finite_pulse:\w+,tau_cap_density:\w+)[^|]* \| (\d+)/(\d+) \| (\d+) \| "
                         r"([^|]+) \| ([^|]+) \|", sections.get("Model sensitivity", ""),
                         re.MULTILINE):
        model_sens[m.group(1)] = dict(eligible=int(m.group(2)), total=int(m.group(3)),
                                      passes=int(m.group(4)), g2_min=m.group(5).strip(),
                                      flux_max=m.group(6).strip())

    def _grab(pattern):
        m = re.search(pattern, ceiling)
        return m.group(1) if m else "n/a"

    reischle = dict(raw=_grab(r"g2\(0\) = ([0-9.]+) raw"),
                    corrected=_grab(r"([0-9.]+) after background correction"),
                    deconv=_grab(r"g2_b\(0\) = ([0-9.]+) \[V\]"),
                    deconv_err=_grab(r"g2_b\(0\) = [0-9.]+ \[V\] \+/- ([0-9.]+)"),
                    factor=_grab(r"about ([0-9.]+)x better"))
    return dict(verdict_line=verdict_line, card_lines=card_lines, cards=cards, ceiling=ceiling,
                missing_claims=missing_claims, fail_reasons=fail_reasons,
                generated=generated, metrics=metrics,
                gamma300_threshold=metrics.get("gamma300_threshold", ""),
                flux_margin=metrics.get("flux_margin", ""), facet_model_note=facet_note,
                gamma300_anchor_lines=anchor_lines, per_T=per_T,
                flux_floor=float(floor_m.group(1)) if floor_m else float("nan"),
                model_sensitivity=model_sens, reischle=reischle)


def headline_model_rows(rows: list[dict], verdict: dict) -> list[dict]:
    """sweep.csv rows of the verdict's own headline model (the VERDICT
    line's `model=finite_pulse:<b>,tau_cap_density:<b>` token), so every
    headline number on a slide comes from the same model the verdict gates
    on. The other three model combinations are non-headline diagnostics."""
    token = verdict["metrics"].get("model", "finite_pulse:true,tau_cap_density:false")
    fp = "finite_pulse:true" in token
    tcd = "tau_cap_density:true" in token
    return [r for r in rows
            if (r.get("model_finite_pulse") == "True") == fp
            and (r.get("model_tau_cap_density") == "True") == tcd]


def static_model_rows(rows: list[dict]) -> list[dict]:
    """The card-level static-loading model (drive.finite_pulse: false,
    tau_cap density law off): what both edge cards ship with, so every number
    from these rows is labelled "static (non-headline)"."""
    return [r for r in rows if r.get("model_finite_pulse") == "False"
            and r.get("model_tau_cap_density") == "False"]


STATIC_LABEL = "static (non-headline)"

SWEEP_KEY_COLS = ("card_id", "delta_xx_meV", "gamma300_meV", "irf_ps", "T_hs_K",
                  "emission_NA", "emission_R_back", "emission_L_um")


def best_row_comparison(rows: list[dict], verdict: dict) -> dict:
    """The headline model's lowest-g2 ELIGIBLE row and the static-loading
    row at the identical configuration -- the difference between the two is
    the re-excitation effect fsim_core.pulse_counting models (a dot re-fills
    and re-emits within one finite pump pulse)."""
    hl = [r for r in headline_model_rows(rows, verdict) if r["eligible_row"] == "True"]
    best = min(hl, key=lambda r: float(r["g2_pulsed"]))
    key = tuple(best[c] for c in SWEEP_KEY_COLS)
    static = next(r for r in static_model_rows(rows)
                  if tuple(r[c] for c in SWEEP_KEY_COLS) == key)
    return dict(best=best, static=static, eligible_rows=hl,
                eligible_cards=sorted(set(r["card_id"] for r in hl)),
                eligible_T=sorted(set(float(r["T_hs_K"]) for r in hl)))


def git_hashes() -> dict:
    """Read-only git lookups: the commit that last wrote out/rt_edge/
    verdict.md (the commit the verdict was computed and committed at) and
    the current HEAD. Falls back to 'unknown' outside a git checkout."""
    def run(args):
        try:
            out = subprocess.run(["git", *args], cwd=str(ROOT), capture_output=True,
                                 text=True, check=True).stdout.strip()
            return out or "unknown"
        except (OSError, subprocess.CalledProcessError):
            return "unknown"
    return dict(verdict=run(["log", "-1", "--format=%h", "--", "out/rt_edge/verdict.md"]),
                head=run(["rev-parse", "--short", "HEAD"]))


def anchor_values(path: Path) -> dict:
    ledger = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {a["id"]: float(a["value"]) for a in ledger["anchors"] if a.get("value") is not None}


def planar_cavity_numbers() -> dict:
    """Audit H5, 'no Purcell recovery': the angle-integrated TOTAL emission-
    rate multiplier of an in-plane dipole in the 300 K-tracked planar DBR
    lambda-cavity, evaluated live with fsim_core.dbr exactly as
    scripts/run_rt_campaign.py's rt_cavity() does. That script runs its whole
    campaign at import time, so its design constants (N_H, N_L, N_C, the
    emitter_energy E0, the pair counts) are read out of its source text with
    ast/regex instead of being copied here. The 1-D on-axis LDOS and the
    120 K total rate come from the committed out/tier_geometry/
    geometry_sweep.csv."""
    import ast
    from fsim_core.cavity import emitter_energy
    from fsim_core.dbr import cavity_stack, planar_total_rate

    src = (ROOT / "scripts" / "run_rt_campaign.py").read_text(encoding="utf-8")
    consts = {}
    for node in ast.parse(src).body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Tuple)):
            names = [t.id for t in node.targets[0].elts]
            if names == ["N_H", "N_L", "N_C"]:
                consts = dict(zip(names, ast.literal_eval(node.value)))
    e0 = float(re.search(r"emitter_energy\(300\.0,\s*([0-9.]+)\)", src).group(1))
    pairs = int(re.search(r"^cav_rt = rt_cavity\((\d+)\)", src, re.MULTILINE).group(1))
    extra = int(re.search(r"cavity_stack\(N_H, N_L, N_C, pairs, pairs \+ (\d+), lam0\)", src).group(1))
    lam0 = 1239.841984 / emitter_energy(300.0, e0)   # hc in eV nm, as run_rt_campaign.py
    top, sp, bot = cavity_stack(consts["N_H"], consts["N_L"], consts["N_C"], pairs,
                                pairs + extra, lam0)
    f_total_300 = float(planar_total_rate(top, bot, consts["N_C"], sp.d_nm, 0.5 * sp.d_nm,
                                          lam0, n_out=consts["N_H"]))
    geo = load_sweep_rows(ROOT / "out" / "tier_geometry" / "geometry_sweep.csv")
    cav = [r for r in geo if r["cavity_top_pairs"] != "none"]
    return dict(pairs=pairs, lam0=lam0, f_total_300=f_total_300,
                pct_300=100.0 * (f_total_300 - 1.0),
                onaxis_max=max(float(r["F_planar_1d_onaxis"]) for r in cav),
                onaxis_min=min(float(r["F_planar_1d_onaxis"]) for r in cav),
                f_total_120_max=max(float(r["F_total_planar"]) for r in cav
                                    if float(r["T_hs_K"]) == 120.0))


def load_sweep_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def first_row_for_card(rows: list[dict], card_id: str) -> dict:
    for row in rows:
        if row["card_id"] == card_id:
            return row
    raise KeyError(card_id)


def load_anchor_sources(path: Path) -> dict:
    """claim -> ordered list of 'Author et al., Journal Vol, Page (Year)'
    citation strings, straight out of the evidence ledger (skips the
    unretrieved-independent-source placeholders)."""
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    by_claim: dict[str, list[str]] = {}
    for anchor in doc["anchors"]:
        src = (anchor.get("source") or "").strip()
        if not src or "not retrieved" in src.lower():
            continue
        by_claim.setdefault(anchor["claim"], [])
        if src not in by_claim[anchor["claim"]]:
            by_claim[anchor["claim"]].append(src)
    return by_claim


# ===================================================================
# Live numbers from fsim_core (never typed in by hand)
# ===================================================================

def gaasp_band_levels() -> DL.DotLevels:
    """The design-card composition, dot_levels.class_presets()'s
    'InP/GaAsP0.4/AlGaAs0.4 on GaAs' (GaAs0.60P0.40 well, matches
    cards/edge-inp-gaasp-design.yaml's ret.system)."""
    presets = DL.class_presets()
    return DL.levels(presets["InP/GaAsP0.4/AlGaAs0.4 on GaAs"])


def material_switch_numbers() -> dict:
    InP = MAT.binary("InP")
    gaasp = MAT.GaAsP(0.4)
    mismatch = MAT.mismatch(gaasp, InP)
    eg_inp = MAT.bandgap(InP, 300.0)
    eg_gaasp = MAT.bandgap(gaasp, 300.0)
    return dict(mismatch_pct=100.0 * mismatch, eg_inp_eV=eg_inp, eg_gaasp_eV=eg_gaasp)


def diode_numbers(sweep_rows: list[dict]) -> dict:
    """sweep_rows should be the headline-model rows (headline_model_rows),
    so the quoted V_j is not a static (non-headline) row's."""
    diode = TR.hkust_preset()
    row = first_row_for_card(sweep_rows, "edge-inp-gaasp-design")
    return dict(diode=diode, V_j=float(row["V_j_pulsed_V"]), I_uA=float(row["I_uA"]),
                d_i_nm=diode.d_i_nm, d_active_nm=diode.d_active_nm)


def waveguide_numbers(sweep_rows: list[dict], card_path: Path) -> dict:
    """The actual edge-emission result for the primary design card, built the
    same way fsim_core.device.evaluate() builds it internally (device.py's
    private _resolve_edge / _retention_system helpers, imported read-only
    here rather than duplicated): the shared ret.system stack (InP dot /
    GaAs0.6P0.4 well / (Al0.50Ga0.50)0.51In0.49P barrier) at the card's own
    668 nm emission wavelength. betas/etas come straight from sweep.csv's
    own edge_beta/edge_eta_total columns (both design cards)."""
    from fsim_core.device import DeviceDesign, _resolve_edge, _retention_system

    design = DeviceDesign.load(str(card_path))
    system = _retention_system(design.ret)
    edge, lambda_nm = _resolve_edge(design.ret, design.emission, 300.0)
    g = system.geometry

    def n_at(mat):
        return MAT.refractive_index(mat.label, lambda_nm)

    layers = [
        WG.Layer("barrier_lower", n_at(system.barrier), design.emission.cladding_nm),
        WG.Layer("matrix_lower", n_at(system.matrix), design.emission.core_half_nm),
        WG.Layer("dot", n_at(system.dot), max(g.height_nm, 0.1), True),
        WG.Layer("matrix_upper", n_at(system.matrix), design.emission.core_half_nm),
        WG.Layer("barrier_upper", n_at(system.barrier), design.emission.cladding_nm),
    ]
    mode = WG.slab_modes(layers, lambda_nm, "TE", 1)[0]
    betas = sorted(set(round(100.0 * float(r["edge_beta"]), 2) for r in sweep_rows))
    return dict(stack=layers, mode=mode, edge=edge, lambda_nm=lambda_nm, beta_pct_range=betas)


def linewidth_numbers() -> dict:
    p = LW.LinewidthParams(gamma300=16.0)
    return dict(params=p, gamma_4K=LW.gamma_anchor(4.0, p), gamma_300K=LW.gamma_anchor(300.0, p))


def evidence_numbers(path: Path) -> dict:
    import json
    doc = json.loads(path.read_text(encoding="utf-8"))
    return dict(passed=doc["checks_passed"], total=doc["checks_total"],
                hallucination_ok=doc.get("hallucination_tests_passed", doc["all_checks_passed"]))


# ===================================================================
# Matplotlib diagrams (the 3-5 explanatory figures this script draws itself)
# ===================================================================

LIGHT = dict(fig="#ffffff", ax="#1b1f24", grid="#d8dee4",
             cb="#1f6feb", vb="#c2410c", accent="#8250df", muted="#57606a")


def _style_ax(ax):
    ax.set_facecolor(LIGHT["fig"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors=LIGHT["ax"], labelsize=9)
    ax.xaxis.label.set_color(LIGHT["ax"])
    ax.yaxis.label.set_color(LIGHT["ax"])
    ax.title.set_color(LIGHT["ax"])


def _trap_step(xs, depth, step):
    """A matrix(0) | dot(-depth) | matrix(0) step profile, then a rise of
    +step in the outer barrier regions -- shared shape for the electron and
    hole trap panels below (both dot_levels.DotLevels quantities V_e/V_h,
    dE_e_matrix/dE_h_matrix are already stored as positive 'depth', so both
    carriers can be drawn with the same sign convention)."""
    y = [step, step, 0.0, 0.0, -depth, -depth, 0.0, 0.0, step, step]
    x = [xs[0], xs[1], xs[1], xs[2], xs[2], xs[3], xs[3], xs[4], xs[4], xs[5]]
    return x, y


def draw_band_diagram(path: Path, lv: DL.DotLevels) -> None:
    """Band diagram of the InP dot in the GaAs0.6P0.4 well: two mirrored trap
    panels (electron confinement, hole confinement), both from a live call to
    fsim_core.dot_levels.levels(). Plotted as 'depth into the well' rather
    than on one shared, true-scale conduction-to-valence-band axis, because
    the ~1.9 eV band gap would otherwise squeeze the ~30-100 meV confinement
    detail into an unreadable sliver; dot_levels.DotLevels already reports
    V_e/V_h/E_e/E_h/dE_e_matrix/dE_h_matrix in exactly this positive-depth
    convention, so no sign or scale is invented for the plot."""
    xs = [-45, -15, -2, 2, 15, 45]
    e_x, e_y = _trap_step(xs, lv.V_e, lv.V_e_mb)
    h_x, h_y = _trap_step(xs, lv.V_h, lv.V_h_mb)
    e_level = -(lv.V_e - lv.E_e)     # = -dE_e_matrix: confined electron level
    h_level = -(lv.V_h - lv.E_h)     # = -dE_h_matrix: confined hole level

    fig, (ax_e, ax_h) = plt.subplots(2, 1, figsize=(8.4, 6.4), dpi=FIG_DPI, sharex=True)
    for ax, x, y, level, depth, carrier, color, esc in (
            (ax_e, e_x, e_y, e_level, lv.V_e, "e-", LIGHT["cb"], lv.dE_e_matrix),
            (ax_h, h_x, h_y, h_level, lv.V_h, "h+", LIGHT["vb"], lv.dE_h_matrix)):
        ax.plot(x, y, color=color, lw=2.5)
        ax.fill_between(x, y, min(y) - 20, color=color, alpha=0.10)
        ax.hlines(level, -2, 2, color=color, lw=1.8, linestyle="--")
        ax.annotate(f"{carrier} confined level\nescape energy {esc:.0f} meV",
                    xy=(2.3, level), xytext=(6, level), fontsize=8.5, color=color, va="center")
        ax.set_ylabel(f"{'electron' if carrier=='e-' else 'hole'} depth (meV)")
        ax.set_ylim(-depth - 25, lv.V_e_mb + 30 if carrier == "e-" else lv.V_h_mb + 30)
        _style_ax(ax)

    ax_e.text(-30, 15, "barrier\n(Al,Ga)InP", ha="center", fontsize=8, color=LIGHT["muted"])
    ax_e.text(-8, 15, "GaAs0.6P0.4\nwell", ha="center", fontsize=8, color=LIGHT["muted"])
    ax_e.annotate("InP dot", xy=(0, -lv.V_e), xytext=(9, -lv.V_e - 15), fontsize=8,
                 color="#0b2545", fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color="#0b2545", lw=0.8))
    ax_e.text(30, 15, "barrier\n(Al,Ga)InP", ha="center", fontsize=8, color=LIGHT["muted"])
    ax_e.set_title(f"Electron trap: {lv.V_e:.0f} meV deep, confined level {lv.E_e:.0f} meV up "
                   "(conduction band)")
    ax_h.set_title(f"Hole trap: {lv.V_h:.0f} meV deep, confined level {lv.E_h:.0f} meV up "
                   "(valence band, hole-energy convention)")
    ax_h.set_xlabel("growth position (nm, illustrative)")
    fig.suptitle(f"An InP dot in a GaAs0.6P0.4 well confines one electron and one hole\n"
                 f"recombination emits {lv.lambda_nm:.0f} nm  (E_X = {lv.E_X_eV:.3f} eV)",
                 fontsize=11.5, color=LIGHT["ax"])
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def draw_pin_diode(path: Path, diode: TR.Diode, V_j: float, I_uA: float) -> None:
    """p-i-n diode cartoon: real layer names and thicknesses from
    fsim_core.transport.hkust_preset()."""
    d_active = diode.d_active_nm
    d_barrier = 0.5 * (diode.d_i_nm - d_active)
    p_w, n_w = 25.0, 25.0
    edges = [0.0, p_w, p_w + d_barrier, p_w + d_barrier + d_active,
             p_w + 2 * d_barrier + d_active, p_w + 2 * d_barrier + d_active + n_w]
    labels = [f"p-cladding\n{diode.p_cladding.label}", f"barrier\n{diode.barrier.label}",
              f"active\n{diode.active.label}\n+ InP dot",
              f"barrier\n{diode.barrier.label}", f"n-cladding\n{diode.n_cladding.label}"]
    colors = ["#c2410c", "#f0b429", "#1f6feb", "#f0b429", "#57606a"]

    fig, ax = plt.subplots(figsize=(8.4, 4.4), dpi=FIG_DPI)
    for i, label in enumerate(labels):
        left, right = edges[i], edges[i + 1]
        ax.axvspan(left, right, color=colors[i], alpha=0.28 if i in (1, 3) else 0.4)
        ax.text(0.5 * (left + right), 1.15, label, ha="center", va="bottom", fontsize=8,
                color=LIGHT["ax"])
    dot_x = 0.5 * (edges[2] + edges[3])
    ax.scatter([dot_x], [0.5], s=140, color="#0b2545", zorder=5)
    ax.text(dot_x, 0.36, "InP dot", ha="center", fontsize=8, color=LIGHT["ax"])

    ax.annotate("", xy=(edges[2] + 0.15 * d_active, 0.75), xytext=(edges[0] + 2, 0.75),
                arrowprops=dict(arrowstyle="-|>", color=LIGHT["vb"], lw=2))
    ax.text(0.5 * (edges[0] + edges[2]), 0.85, "holes +", color=LIGHT["vb"], fontsize=9, ha="center")
    ax.annotate("", xy=(edges[3] - 0.15 * d_active, 0.25), xytext=(edges[5] - 2, 0.25),
                arrowprops=dict(arrowstyle="-|>", color=LIGHT["cb"], lw=2))
    ax.text(0.5 * (edges[3] + edges[5]), 0.13, "electrons -", color=LIGHT["cb"], fontsize=9, ha="center")

    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, 1.3)
    ax.set_yticks([])
    ax.set_xlabel(f"growth position (nm); intrinsic region {diode.d_i_nm:.0f} nm, "
                  f"active layer {d_active:.0f} nm")
    ax.set_title(f"Forward bias pushes carriers into the dot  "
                 f"(V_j ≈ {V_j:.2f} V at I = {I_uA*1000:.1f} nA, 300 K)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def draw_spectral_overlap(path: Path, p: LW.LinewidthParams, delta_xx_lo: float,
                          delta_xx_hi: float) -> None:
    """Exciton (X) and biexciton (XX) lines at 4 K vs 300 K, widths from
    fsim_core.linewidth.gamma_anchor, splitting from the design cards'
    dot.delta_xx range (4-7 meV, docs/rt_edge_contract.md)."""
    def lorentzian(e, e0, gamma):
        hw = 0.5 * gamma
        return (hw ** 2) / ((e - e0) ** 2 + hw ** 2)

    e = np.linspace(-15, 30, 2000)
    g4, g300 = LW.gamma_anchor(4.0, p), LW.gamma_anchor(300.0, p)
    delta_mid = 0.5 * (delta_xx_lo + delta_xx_hi)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), dpi=FIG_DPI, sharey=True)
    for ax, T, gamma in zip(axes, [4, 300], [g4, g300]):
        x_line = lorentzian(e, 0.0, gamma)
        xx_line = lorentzian(e, delta_mid, gamma)
        ax.plot(e, x_line, color=LIGHT["cb"], lw=1.8, label="exciton (X)")
        ax.plot(e, xx_line, color=LIGHT["vb"], lw=1.8, label="biexciton (XX)")
        ax.fill_between(e, np.minimum(x_line, xx_line), color=LIGHT["accent"], alpha=0.35)
        ax.set_title(f"T = {T} K   Γ ≈ {gamma:.2f} meV")
        ax.set_xlabel("energy from X line (meV)")
        _style_ax(ax)
    axes[0].set_ylabel("intensity (normalized)")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(f"X-XX splitting {delta_xx_lo:.0f}-{delta_xx_hi:.0f} meV vs 300 K linewidth "
                 f"class range {LW.GAMMA300_CLASS_RANGE[0]:.0f}-{LW.GAMMA300_CLASS_RANGE[1]:.0f} meV",
                 fontsize=10, color=LIGHT["ax"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def draw_waveguide_mode(path: Path, stack: list[WG.Layer], mode: WG.SlabMode,
                        edge: WG.EdgeResult) -> None:
    """Ridge waveguide cross-section with the real vertical mode profile from
    fsim_core.waveguide.slab_modes on the primary design card's own edge-
    emission stack (barrier / matrix / dot / matrix / barrier, the same 5
    layers fsim_core.device._resolve_edge builds for cards/edge-inp-gaasp-
    design.yaml)."""
    edges = [0.0]
    for layer in stack:
        edges.append(edges[-1] + layer.thickness_nm)
    center = edges[-1] / 2.0
    window = 260.0
    z = mode.z_nm - center
    intensity = mode.field ** 2
    intensity = intensity / intensity.max()

    fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=FIG_DPI)
    palette = {"barrier_lower": "#8250df", "matrix_lower": "#1f6feb", "dot": "#0b2545",
               "matrix_upper": "#1f6feb", "barrier_upper": "#8250df"}
    for i, layer in enumerate(stack):
        left, right = edges[i] - center, edges[i + 1] - center
        if right < -window or left > window:
            continue
        ax.axvspan(max(left, -window), min(right, window),
                  color=palette.get(layer.name, "#c8d1da"), alpha=0.30)
    ax.plot(z, intensity, color="#d1242f", lw=2.0, label="guided mode |E|² (normalized)")
    ax.set_xlim(-window, window)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("vertical position (nm, dot at 0; claddings continue off-frame at ±1000 nm)")
    ax.set_ylabel("mode intensity (normalized)")
    ax.set_title(f"n_eff = {edge.n_eff:.3f},  β = {100*edge.beta:.2f}%,  "
                 f"facet transmission = {100*edge.T_facet:.0f}%")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def draw_g2_dip(path: Path, g2_best: float, g2_ceiling: float) -> None:
    """Schematic antibunching dip: the depth (g2 at zero delay) is the real
    pooled pulsed minimum from out/rt_edge/verdict.md (headline finite-pulse
    model); the best published 300 K single-dot value (optically pumped,
    Laferriere et al. 2023, read from the evidence ledger) is drawn for
    comparison. The delay axis is illustrative (the sweep reports the scalar
    g2(0), not a measured g2(tau) trace)."""
    tau = np.linspace(-4, 4, 800)

    def dip(g0, width):
        return 1.0 - (1.0 - g0) * np.exp(-np.abs(tau) / width)

    fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=FIG_DPI)
    ax.plot(tau, dip(g2_best, 0.7), color=LIGHT["cb"], lw=2.2,
            label=f"this design's best eligible corner (finite-pulse headline), g2(0) = {g2_best:.4f}")
    ax.plot(tau, dip(g2_ceiling, 1.0), color=LIGHT["muted"], lw=1.8, linestyle="--",
            label=f"best published 300 K dot (optical pump), g2(0) = {g2_ceiling:.2f}")
    ax.axhline(0.5, color=LIGHT["vb"], lw=1.2, linestyle=":", label="single-photon threshold (0.5)")
    ax.axhline(1.0, color=LIGHT["ax"], lw=1.0, linestyle=":", alpha=0.6)
    ax.text(3.6, 1.02, "ordinary lamp", fontsize=8, color=LIGHT["ax"], ha="right")
    ax.text(0, -0.08, "0 = perfect single photon", fontsize=8, color=LIGHT["ax"], ha="center")
    ax.set_xlabel("delay between photon arrivals (arbitrary units)")
    ax.set_ylabel("g2(delay)")
    ax.set_ylim(-0.15, 1.2)
    ax.set_title("The antibunching dip: how often two photons arrive together")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def draw_linewidth_sweep(path: Path, rows: list[dict]) -> dict:
    """Plot the CSV's pulsed g2 landscape at the largest X-XX splitting, one
    line per heat-sink temperature. `rows` must be the headline-model rows
    only (headline_model_rows): mixing in the three non-headline model
    combinations would put their lower g2 on a headline figure. The
    diagnostic column diag_g2_pulsed is used so that rows below the flux
    floor still show their physics value.

    The two reference linewidths are read from the shared evidence ledger,
    rather than copied into this presentation source.
    """
    anchor = anchor_values(ANCHORS_YAML)
    best_split = max(float(r["delta_xx_meV"]) for r in rows)
    selected = [r for r in rows if float(r["delta_xx_meV"]) == best_split]
    by_T: dict[float, dict[float, list[float]]] = {}
    for row in selected:
        g = float(row["diag_g2_pulsed"])
        if not np.isfinite(g):
            continue
        by_T.setdefault(float(row["T_hs_K"]), {}).setdefault(
            float(row["gamma300_meV"]), []).append(g)
    low = anchor["chatzarakis23-gamma300-class"]
    mid = anchor["matsuda01-gamma300-class"]
    fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=FIG_DPI)
    colors = [LIGHT["cb"], LIGHT["accent"], "#1a7f37", LIGHT["muted"]]
    all_g2 = []
    for color, T in zip(colors, sorted(by_T)):
        gam = sorted(by_T[T])
        g2 = [min(by_T[T][x]) for x in gam]
        all_g2 += g2
        ax.plot(gam, g2, "o-", color=color, lw=2.0, label=f"T_hs = {T:g} K")
    ax.axhline(0.5, color=LIGHT["vb"], lw=1.5, linestyle="--",
               label="pulsed g2 threshold = 0.5")
    for value, label, color in ((low, f"verified anchor {low:g} meV", LIGHT["vb"]),
                                (mid, f"verified anchor {mid:g} meV", LIGHT["accent"])):
        ax.axvline(value, color=color, lw=1.2, linestyle=":", label=label)
    ax.set_xlabel("300 K exciton linewidth gamma300 (meV)")
    ax.set_ylabel("pulsed intrinsic g2(0), best lever combination")
    ax.set_title(f"Headline finite-pulse model, X-XX splitting {best_split:g} meV: "
                 "g2 stays near 1 at every linewidth")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=8, frameon=False, loc="lower right", ncol=2)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)
    gammas = sorted({float(r["gamma300_meV"]) for r in rows})
    return {"split": best_split, "low": low, "mid": mid, "gamma_lo": gammas[0],
            "gamma_hi": gammas[-1], "g2_lo": min(all_g2), "g2_hi": max(all_g2)}


def draw_flux_vs_T(path: Path, rows: list[dict], floor: float) -> dict:
    """Maximum collected pulsed flux per heat-sink temperature and card,
    headline-model rows of out/rt_edge/sweep.csv, against the collected-flux
    eligibility floor parsed from verdict.md."""
    by = {}
    for r in rows:
        by.setdefault(r["card_id"], {}).setdefault(float(r["T_hs_K"]), []).append(
            float(r["collected_flux_pulsed_s"]))
    fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=FIG_DPI)
    for color, card in zip((LIGHT["cb"], LIGHT["vb"]), sorted(by)):
        T = sorted(by[card])
        f = [max(by[card][t]) for t in T]
        ax.plot(T, f, "o-", color=color, lw=2.2, label=card)
        for t, v in zip(T, f):
            ax.annotate(f"{v:.0f}", xy=(t, v), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=8, color=color)
    ax.axhline(floor, color=LIGHT["ax"], lw=1.4, linestyle="--",
               label=f"eligibility floor {floor:g} photons/s")
    ax.set_yscale("log")
    ax.set_xlabel("heat-sink temperature T_hs (K)")
    ax.set_ylabel("max collected pulsed flux (photons/s)")
    ax.set_title("Headline finite-pulse model: the floor is cleared only at the coldest T_hs")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)
    return by


def draw_reexcitation(path: Path, headline: list[dict], static: list[dict], card: str) -> None:
    """Lowest pulsed intrinsic g2(0) per heat-sink temperature on one card,
    headline finite-pulse model (re-excitation within the pump pulse kept,
    fsim_core.pulse_counting) versus the card-level static-loading model
    (labelled static (non-headline)). Both from out/rt_edge/sweep.csv's
    diag_g2_pulsed column (valid rows, including rows below the flux floor)."""
    def curve(rows):
        by = {}
        for r in rows:
            g = float(r["diag_g2_pulsed"])
            if r["card_id"] == card and np.isfinite(g):
                by.setdefault(float(r["T_hs_K"]), []).append(g)
        T = sorted(by)
        return T, [min(by[t]) for t in T]

    fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=FIG_DPI)
    for rows, color, style, label, dy in (
            (headline, LIGHT["cb"], "o-", "finite-pulse model (headline): re-excitation kept", 7),
            (static, LIGHT["muted"], "s--", f"static loading, {STATIC_LABEL}", -14)):
        T, g = curve(rows)
        ax.plot(T, g, style, color=color, lw=2.2, label=label)
        for t, v in zip(T, g):
            ax.annotate(f"{v:.3f}", xy=(t, v), xytext=(0, dy), textcoords="offset points",
                        ha="center", fontsize=8, color=color)
    ax.axhline(0.5, color=LIGHT["vb"], lw=1.5, linestyle=":", label="single-photon threshold 0.5")
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("heat-sink temperature T_hs (K)")
    ax.set_ylabel("lowest pulsed intrinsic g2(0)")
    ax.set_title(f"{card}: re-excitation during the pulse pushes g2 toward 1")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(path, facecolor=LIGHT["fig"])
    plt.close(fig)


def verdict_brightness_factors() -> dict:
    """Read the verdict's brightness-decomposition table for result slides.
    pkg5b: when no row is eligible (n_eligible == 0), run_rt_edge.py's
    write_markdown prints a "Why no row is eligible" heading in place of
    "Best diagnostic-g2 row and brightness decomposition" ahead of the SAME
    factor table -- try both headings rather than crashing on whichever one
    the current run did not print."""
    text = VERDICT_MD.read_text(encoding="utf-8")
    sections = split_markdown_sections(text)
    section = sections.get("Best diagnostic-g2 row and brightness decomposition")
    if section is None:
        section = sections.get("Why no row is eligible", "")
    return {name.strip(): value.strip() for name, value in
            re.findall(r"^\| ([^|]+) \| ([^|]+) \|$", section, re.MULTILINE)}


# ===================================================================
# Slide content
# ===================================================================

class Slide:
    __slots__ = ("title", "bullets", "notes", "image", "verbatim")

    def __init__(self, title, bullets, notes, image=None, verbatim=None):
        self.title = title
        self.bullets = bullets
        self.notes = notes
        self.image = image          # Path or None
        self.verbatim = verbatim    # optional monospace data line(s), list[str]


def build_slides(verdict: dict, sweep_rows: list[dict], lv: DL.DotLevels,
                 mat_sw: dict, diode_info: dict, wg_info: dict, lw_info: dict,
                 evidence: dict, anchors: dict, figures: dict, linewidth_sweep: dict,
                 brightness: dict, cmp: dict, planar: dict, hashes: dict,
                 anchor_vals: dict) -> list[Slide]:
    """The 19-slide short deck. Every number is parsed from out/rt_edge/* and
    out/tier_geometry/* or computed live from fsim_core; nothing physical is
    typed in here."""
    def cite(claim, idx=0, default=""):
        srcs = anchors.get(claim, [])
        return srcs[idx] if idx < len(srcs) else default

    m = verdict["metrics"]
    beta_lo, beta_hi = min(wg_info["beta_pct_range"]), max(wg_info["beta_pct_range"])
    beta_range = (f"{beta_lo:.2f}%" if f"{beta_lo:.2f}" == f"{beta_hi:.2f}"
                  else f"{beta_lo:.2f}-{beta_hi:.2f}%")
    flux_margin = verdict["flux_margin"]
    floor = verdict["flux_floor"]
    anchor_text = " ".join(verdict["gamma300_anchor_lines"])
    loading_key = next(key for key in brightness if key.startswith("loading = 1 - e^-mu"))
    mean_counts_key = next((key for key in brightness
                            if key.startswith("finite_pulse_mean_counts")), None)
    eta_total_key = next(key for key in brightness if key.startswith("eta_total"))
    beta_key = next(key for key in brightness if key.startswith("beta"))
    best, static = cmp["best"], cmp["static"]
    static_g2 = float(static["diag_g2_pulsed"])
    static_flux = float(static["collected_flux_pulsed_s"])
    pulse_ps = 1000.0 * float(best["pulse_width_ns"])
    rep_mhz = float(best["rep_rate_hz"]) / 1e6
    elig_T = ", ".join(f"{t:g}" for t in cmp["eligible_T"])
    elig_cards = ", ".join(cmp["eligible_cards"])
    per_T = verdict["per_T"]
    flux_by_T = ", ".join(f"{per_T[t]['flux_max']} at {t} K" for t in sorted(per_T))
    n_missing = len(verdict["missing_claims"])
    fail_reasons = ", ".join(verdict["fail_reasons"])
    static_sens = verdict["model_sensitivity"].get("finite_pulse:false,tau_cap_density:false", {})
    reis = verdict["reischle"]
    chatz_g2 = anchor_vals["chatzarakis23-g2-temperature"]
    laf_g2 = anchor_vals["laferriere23-g2-temperature"]
    pct = planar["pct_300"]
    commit_line = (f"verdict.md computed and committed at {hashes['verdict']} "
                   f"(deck built at HEAD {hashes['head']})")

    slides = []

    slides.append(Slide(
        "qd-photon-sim",
        ["Branch: rt-edge-emitter",
         "A physics simulator that designs and checks a quantum-dot single-photon source",
         "This deck reports what the simulator found, including where the design falls short"],
        "Project: qd-photon-sim, branch rt-edge-emitter. This is a computational design study "
        "produced by the simulator's own pipeline and acceptance sweep; no author or institution "
        "is claimed beyond the repository itself. The question it answers is narrow and practical: "
        "can an electrically driven, edge-emitting indium-phosphide quantum dot emit single photons "
        "without a cryostat, at a heat-sink temperature of 230 K or warmer? The short answer, "
        "which the last third of the deck builds up step by step, is no: the corrected model "
        "collects enough light at 230 K but the photons are not single. Every number quoted "
        "later in this deck traces to a real file the repository generated "
        f"({VERDICT_MD.name}, {SWEEP_CSV.name}, the geometry sweep) or to a live call into the "
        "simulator's physics modules made while building this deck, so the slides refresh "
        "automatically whenever the sweep is re-run."))

    slides.append(Slide(
        "What is a quantum dot?",
        ["A speck of semiconductor a few nanometres across, grown inside a larger crystal",
         "Electrons and holes get trapped inside it, like an atom -- an ‘artificial atom’",
         "The trapped charges recombine and emit exactly one photon at a time"],
        "Self-assembled quantum dots form during epitaxial crystal growth (Stranski-Krastanov "
        "growth: a thin strained layer of one semiconductor spontaneously balls up into "
        "nanometre-scale islands on top of another). Because the dot material has a smaller "
        "electronic band gap than the crystal around it, an electron and a hole that wander into "
        "the dot get trapped there, just as an electron is trapped around a nucleus in a real "
        "atom -- hence 'artificial atom'. When the trapped electron and hole recombine they give "
        "up their energy as one photon. Because the dot holds only a few discrete states, it "
        "cannot emit a second identical photon until it has been refilled, and that refilling "
        "delay is what makes the output arrive one photon at a time. The rest of this deck asks "
        "whether that delay survives at room temperature, where heat shakes carriers out of the "
        "dot and broadens its emission lines."))

    slides.append(Slide(
        "Why one photon at a time matters",
        ["Quantum key distribution needs single photons -- extra photons can be intercepted and copied",
         "One photon per pulse makes the resulting encryption key provably unforgeable",
         "g2(0) measures this directly: 0 is perfect, 1 is an ordinary lamp"],
        "In quantum key distribution, security depends on each pulse carrying at most one photon: "
        "if a pulse ever contains two photons, an eavesdropper can silently siphon off the extra "
        "one without being detected, breaking the security proof. The standard way to measure how "
        "close a source comes to 'exactly one photon' is the second-order correlation function "
        "g2(0): a perfect single-photon source gives g2(0) = 0, while an ordinary lamp (many "
        "photons, uncorrelated) gives g2(0) = 1. This design study's target is g2(0) < 0.5, the "
        "usual working definition of 'single-photon' behaviour. Keep that scale in mind: every "
        "later slide that reports a g2 value is really saying how far the device sits between "
        "a perfect single-photon emitter and a light bulb, and the headline value this deck "
        "arrives at sits very close to the light-bulb end."))

    slides.append(Slide(
        "Why room temperature and electrical driving matter",
        ["Every working single-photon quantum dot today needs a cryostat and a laser",
         "A source that runs at 300 K on a simple diode current plugs in like an LED",
         "That is the difference between a lab demonstration and a deployable device"],
        "Essentially all published single-photon quantum dots operate in a liquid-helium cryostat "
        "(4-80 K) and are excited optically by a separate pump laser tuned onto the dot. A source "
        "that instead runs at room temperature (300 K) and is turned on and off with an ordinary "
        "p-i-n diode current -- the same physics as an LED -- removes both the cryostat and the "
        "laser, which is what would actually make a single-photon source practical to deploy "
        "outside a physics lab. The acceptance rule used here is already relaxed from that ideal: "
        "it accepts any heat-sink temperature from 230 K upward, which a small thermoelectric "
        "cooler can reach, so a pass at 230 K would still count as a deployable, cryostat-free "
        "device. Even with that relaxation, the design fails, as the verdict slide shows."))

    slides.append(Slide(
        "Why edge (in-plane) emission matters",
        ["Most quantum-dot sources shine light straight up, out of the top of the chip",
         "Edge emission sends photons sideways, directly into a waveguide or an optical fibre",
         "That is how the light gets onto a photonic chip or into a communication line"],
        "The vertical single-photon sources in the literature (micropillars, DBR cavities, "
        "resonant-cavity LEDs) collect light out of the top surface, which suits a free-space "
        "detector but not a photonic circuit. An edge-emitting, in-plane design instead guides "
        "the light sideways through a ridge waveguide etched into the chip, so it can be coupled "
        "directly into an on-chip waveguide network or into an optical fibre butted against the "
        "cleaved facet -- the geometry a real photonic-integrated single-photon source would need. "
        "The price of this geometry is collection efficiency: only the few percent of the dot's "
        "light that couples into the guided mode ever reaches the facet, so brightness, not only "
        "purity, becomes a gate the design has to clear. The brightness slide later shows "
        "exactly how much of the light survives that path."))

    slides.append(Slide(
        "What this software does",
        ["A physics calculator, not a lab measurement: predicts g2(0), brightness and the operating window of a design",
         "Built from published material data, not fitted to make an answer come out right",
         "It also checks its own literature citations for errors"],
        "qd-photon-sim is a simulator: given a 'design card' describing a dot, its host materials, "
        "the drive electronics and the output waveguide, it predicts the single-photon purity "
        "g2(0), the collected brightness, and the temperature/current window over which the design "
        "would work -- all built up from published material constants (lattice constants, band "
        "gaps, effective masses, linewidths) rather than tuned to produce a desired answer. A "
        "companion module (verify/verify_rt_edge_papers.py) independently re-checks every cited "
        "literature number for unit errors, out-of-range values, and duplicated or unsourced "
        f"citations; the current evidence ledger passes {evidence['passed']}/{evidence['total']} "
        "such checks. The simulator is therefore only as good as its inputs, and the deck flags "
        "every input that is an estimate or an assumption rather than a measurement, so the "
        "audience can see which conclusions rest on data."))

    slides.append(Slide(
        "How it works: the pipeline",
        ["Materials → dot levels → diode injection → heating → line overlap → filtering → waveguide → g2",
         "Every step keeps a literature tag: verified, derived, class estimate, or assumption",
         "Every conclusion is a range swept over the uncertain inputs, not one confident number"],
        "The pipeline: (1) materials.py supplies band gaps, lattice constants and band offsets from "
        "Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001); (2) dot_levels.py solves the confined "
        "electron/hole levels of the dot-in-a-well-in-a-barrier stack; (3) transport.py "
        "models the p-i-n diode and carrier injection into the dot; (4) thermal.py raises the "
        "junction temperature under drive; (5) linewidth.py gives the exciton/biexciton "
        "linewidths at that temperature; (6) the spectral overlap sets how much biexciton light "
        "leaks through a filter tuned to the exciton; (7) waveguide.py computes how much light the "
        "edge waveguide collects; (8) the device evaluator combines all of this into g2(0) and "
        "brightness, and for pulsed drive pulse_counting.py follows the dot through the whole "
        "finite pump pulse. Every number carries a provenance tag ([V] verified, [DR] derived, "
        "[E] class estimate, [A] assumption), and the acceptance sweep varies the uncertain "
        "inputs over their literature ranges rather than picking one value."))

    slides.append(Slide(
        "Anatomy of the artificial atom",
        [f"An InP dot inside a GaAs0.6P0.4 well traps one electron {lv.E_e:.0f} meV deep and one hole {lv.E_h:.0f} meV deep",
         f"Recombination of that trapped pair is the single photon, here at {lv.lambda_nm:.0f} nm"],
        "This band diagram is drawn from a live call to fsim_core.dot_levels.levels() on the "
        "simulator's own 'InP/GaAsP0.4/AlGaAs0.4 on GaAs' design-card class: a 4 nm InP dot in a "
        f"GaAs0.6P0.4 well gives an electron confinement V_e = {lv.V_e:.0f} meV and a hole "
        f"confinement V_h = {lv.V_h:.0f} meV (band offsets from materials.offsets, Vurgaftman et "
        f"al. 2001 model-solid theory); the bound levels sit E_e = {lv.E_e:.0f} meV and "
        f"E_h = {lv.E_h:.0f} meV inside those wells, giving an exciton transition energy "
        f"E_X = {lv.E_X_eV:.3f} eV ({lv.lambda_nm:.0f} nm). This is dot_levels.py's separable-disk "
        "model, honestly documented as accurate to tens of meV for this geometry class, not a "
        "full k.p calculation. The escape energies printed on the figure matter later: the "
        "shallower they are, the faster heat empties the dot and the faster it has to be "
        "refilled, which is what drives the re-excitation problem on the verdict slides.",
        image=figures["band"]))

    slides.append(Slide(
        "The material swap: fixing an inverted design",
        ["The literal 'InP cladding, GaAsP active' is backwards: InP has the smaller band gap",
         f"GaAsP is {mat_sw['mismatch_pct']:.1f}% lattice-mismatched to InP -- it cannot grow as a cladding",
         "Instead: InP dots inside a GaAsP well, clad by (Al,Ga)InP -- what real labs grow"],
        "The starting brief asked for 'InP cladding, GaAsP active region'. Checked against the "
        f"materials database this is physically inverted: InP's 300 K band gap is "
        f"{mat_sw['eg_inp_eV']:.2f} eV versus GaAs0.6P0.4's {mat_sw['eg_gaasp_eV']:.2f} eV "
        "(fsim_core.materials.bandgap), so InP is the narrower-gap material and cannot act as a "
        "confining cladding around a wider-gap GaAsP active layer. Separately, GaAsP is "
        f"{mat_sw['mismatch_pct']:.1f}% lattice-mismatched to InP (fsim_core.materials.mismatch, "
        "Vurgaftman et al. 2001), too large to grow a coherent strained layer of useful "
        "thickness. The platform that does exist in the literature is InP self-assembled quantum "
        "dots grown inside a GaAs1-xPx quantum well, clad by (Al,Ga)InP or AlGaAs, emitting in "
        "the 660-755 nm range (Gu et al., Optics Express 33, 23732 (2025), HKUST class; Reischle "
        "et al., Optics Express 16, 12771 (2008), Stuttgart InP/GaInP class). This study uses "
        "that platform as both its primary (GaAsP-well) and fallback (GaInP-well) design cards.",
        image=figures["band"]))

    slides.append(Slide(
        "Turning on the light: the p-i-n diode",
        ["A p-i-n diode: p-doped and n-doped layers around an intrinsic core",
         "Electrons and holes are pushed in from opposite sides and meet at the dot",
         f"At this design's operating current the junction sits near {diode_info['V_j']:.2f} V forward bias"],
        "Electrical driving replaces the pump laser used by almost every other single-photon "
        "quantum-dot demonstration. fsim_core.transport.hkust_preset() builds the actual p-i-n "
        f"stack used here: {diode_info['diode'].p_cladding.label} p- and n-claddings around a "
        f"{diode_info['diode'].d_i_nm:.0f} nm intrinsic region, with the {diode_info['diode'].active.label} "
        f"active layer ({diode_info['diode'].d_active_nm:.0f} nm) at its centre, following Gu et "
        "al.'s (Optics Express 33, 23732, 2025) own electrically pumped device. Forward bias tilts "
        "the bands so electrons flow in from the n-side and holes from the p-side and both are "
        "captured by the dot in the middle. At the single-dot-scale drive current this design uses "
        f"(I = {diode_info['I_uA']*1000:.1f} nA), fsim_core.transport's own I-V solver puts the "
        f"junction voltage at V_j ≈ {diode_info['V_j']:.2f} V (first headline-model sweep row of "
        "the primary card). In pulsed operation this current is switched on for "
        f"{pulse_ps:.0f} ps at {rep_mhz:.0f} MHz; how fast the dot refills during that window is "
        "the crux of the verdict.", image=figures["pin"]))

    slides.append(Slide(
        "Two lines, one at a time: exciton and biexciton",
        ["One trapped pair (exciton, X) and two trapped pairs (biexciton, XX) emit at slightly different colours",
         "At 4 K the two lines are razor-sharp and clearly separate",
         f"At 300 K thermal broadening ({LW.GAMMA300_CLASS_RANGE[0]:.0f}-{LW.GAMMA300_CLASS_RANGE[1]:.0f} meV) swallows the 3-7 meV gap between them"],
        "Filtering out only the exciton (single-pair) line and rejecting the biexciton (two-pair) "
        "line is how a quantum dot is made to emit one photon at a time: if the two lines are "
        "spectrally distinct, a narrow filter passes X and blocks XX. fsim_core.linewidth.py's "
        "literature-anchored model (Matsuda et al., Phys. Rev. B 63, 121304(R) (2001); Laferriere "
        "et al., Nano Letters 23, 962 (2023); Chatzarakis et al., Phys. Rev. Applied 20, 034011 "
        f"(2023)) gives Γ(4 K) ≈ {lw_info['gamma_4K']:.2f} meV, far narrower than the "
        "X-XX splitting; but Γ(300 K) is anchored to a class range (no InP/GaAsP single-dot "
        "300 K linewidth has been published), comparable to or larger than the X-XX splittings "
        "seen across the InP-dot literature (Reischle et al. 2008; Bommer et al., J. Appl. Phys. "
        "110, 063108 (2011)). At 300 K the two lines overlap for essentially every published "
        "dot. This slide explains why purity is hard at room temperature; the results slides "
        "show that, in the corrected model, it is not even the dominant problem.",
        image=figures["overlap"]))

    slides.append(Slide(
        "Getting light out: the edge waveguide",
        ["A ridge waveguide guides light sideways, to the chip edge, not straight up",
         f"After the audit's beta fix (H4), the mode captures about {beta_range} of the dot's light",
         f"An uncoated facet lets about {100*wg_info['edge'].T_facet:.0f}% of that light escape the chip"],
        "fsim_core.waveguide.py solves the guided optical mode of the ridge stack (an "
        "(Al,Ga)InP core around the GaAsP well, on a GaAs substrate) with a scalar effective-index "
        "model (Lecamp, Lalanne & Hugonin, Phys. Rev. Lett. 99, 023902 (2007)). For this design's "
        f"stack at {wg_info['lambda_nm']:.0f} nm the solved mode has n_eff = {wg_info['edge'].n_eff:.2f} "
        f"and a beta factor (fraction of the dot's spontaneous emission that couples into the "
        f"guided mode) of {100*wg_info['edge'].beta:.2f}%, ranging {beta_range} "
        "across the acceptance-sweep grid. Audit finding H4 corrected the mode area used in the "
        "guided Purcell factor to the energy-normalised area at the dot, which roughly doubled "
        "beta; the verdict's best row now reports beta = "
        f"{brightness[beta_key]}. An uncoated facet then transmits "
        f"{100*wg_info['edge'].T_facet:.0f}% of the guided light out of the chip (Fresnel "
        "reflection). Per Lemma 1 of docs/rt_edge_contract.md, none of this changes the "
        "intrinsic multiphoton probability; it only sets how much emitted light is collected.",
        image=figures["wg"]))

    slides.append(Slide(
        "g2(0): the number that matters",
        ["g2(0): how often two photons arrive together. 0 = perfect single photon, 1 = ordinary lamp",
         "0.5 is the usual single-photon threshold",
         f"Headline finite-pulse model: the best eligible corner reaches only {m['g2_min']} -- close to a lamp"],
        "g2(0) is the only equation-like quantity this deck relies on, and even it is explained "
        "rather than derived: it is the probability of detecting two photons at (almost) the same "
        "instant, normalised so that a classical, many-photon source gives 1 and a perfect "
        "single-photon emitter gives 0. The dashed curve is the best published 300 K single-dot "
        f"value in the evidence ledger, g2(0) = {laf_g2:g} (Laferriere et al. 2023, optically "
        "pumped InAsP/InP nanowire dot). The solid curve is this design's best eligible corner "
        f"under the headline finite-pulse model: g2(0) = {m['g2_min']} (pooled median over the "
        f"eligible rows {m['g2_median_eligible']}). That is worse than the published optical "
        "result and nowhere near the 0.5 threshold. The delay axis of the figure is "
        "illustrative: the sweep reports the scalar g2(0), not a measured g2(tau) trace, so only "
        "the depth of each dip carries data. The next slides show where that value comes from.",
        image=figures["g2"]))

    slides.append(Slide(
        "Results: the room-temperature acceptance sweep",
        [f"Flux clears the {floor:g} photons/s floor only at T_hs = {elig_T} K: max {m['flux_max']} photons/s, "
         f"{m['eligible']} rows eligible, all on {elig_cards}",
         f"Pulsed g2(0) at those rows: min {m['g2_min']}, median {m['g2_median_eligible']} -- far above 0.5",
         f"Across the swept 300 K linewidth ({linewidth_sweep['gamma_lo']:g}-{linewidth_sweep['gamma_hi']:g} meV) "
         f"g2 stays at {linewidth_sweep['g2_lo']:.3f}-{linewidth_sweep['g2_hi']:.3f}",
         "The 6.5 meV linewidth anchor fails on both cards"],
        "The VERDICT line printed at the bottom of this slide is copied verbatim from "
        f"out/rt_edge/verdict.md ({commit_line}). The sweep was generated by "
        f"scripts/run_rt_edge.py ({verdict['generated']}) over the full contract grid "
        "(docs/rt_edge_contract.md): X-XX splitting, 300 K linewidth, detector response, "
        "heat-sink temperature and the collection levers, for both the primary (GaAsP-well) and "
        "fallback (GaInP-well) cards, with the corrected finite-pulse model as the headline. "
        f"Only {m['eligible']} rows clear the collected-flux floor, all at the coldest heat sink; "
        f"there the lowest pulsed g2 is {m['g2_min']}. The figure plots, from sweep.csv at the "
        f"widest {linewidth_sweep['split']:g} meV splitting, the lowest g2 at each linewidth and "
        "temperature: the curves are flat and near 1, so a narrower line would not rescue the "
        f"design. At the 6.5 meV anchor: {anchor_text} The verified class anchors "
        f"{linewidth_sweep['low']:g} and {linewidth_sweep['mid']:g} meV come from the evidence ledger.",
        image=figures["linewidth"], verbatim=[verdict["verdict_line"], commit_line]))

    slides.append(Slide(
        "Brightness: the flux floor is cleared only at 230 K",
        [f"Flux margin (flux_max / floor) is {flux_margin}: above 1, so the floor is cleared",
         f"Flux falls as the heat sink warms: {flux_by_T} photons/s",
         f"Post-H4 waveguide coupling at the best row: beta = {brightness[beta_key]}",
         f"No Purcell recovery: a {planar['pairs']}-pair planar DBR cavity adds only +{pct:.1f}% to the total emission rate at 300 K"],
        "The figure plots the maximum collected pulsed flux per heat-sink temperature and card "
        "from the headline-model rows of out/rt_edge/sweep.csv, against the floor parsed from "
        "verdict.md. The best row's flux chain, read from verdict.md: "
        f"finite_pulse_mean_counts = {brightness.get(mean_counts_key, 'n/a')}, "
        f"eta_total = {brightness[eta_total_key]}, repetition rate = {brightness['rep rate (Hz)']} "
        f"Hz, giving {brightness['reported collected_flux_pulsed_s']}. The dominant limiter is "
        f"confinement retention S = {brightness['S (confinement retention)']}: heat empties the "
        f"dot. (Static factors loading = {brightness[loading_key]} and t_X = "
        f"{brightness['t_X (spectral transmission)']} are reference only.) A cavity cannot fix "
        f"this: its 1-D on-axis LDOS reaches {planar['onaxis_min']:.0f}-{planar['onaxis_max']:.0f}x "
        "(out/tier_geometry/geometry_sweep.csv), but the angle-integrated total rate, computed "
        f"live with fsim_core.dbr.planar_total_rate at {planar['lam0']:.0f} nm, is "
        f"{planar['f_total_300']:.4f} at 300 K (+{pct:.1f}%; the 120 K geometry sweep gives at "
        f"most {planar['f_total_120_max']:.4f}). A planar cavity redistributes emission, it does "
        "not speed it up.", image=figures["flux"]))

    slides.append(Slide(
        "Validation: 80 K and 230 K comparisons",
        [f"Reischle 2008 (80 K, electrical): g2(0) = {reis['raw']} raw, {reis['deconv']} +/- {reis['deconv_err']} "
         f"IRF-deconvolved, {reis['corrected']} background-corrected",
         f"Convention-matched, that 80 K device is about {reis['factor']}x better than this design's best 230 K corner",
         f"Chatzarakis 2023 measured g2(0) = {chatz_g2:g} at 230 K (optical pump); the headline model gives {m['g2_min']} here",
         f"Evidence ledger: {evidence['passed']}/{evidence['total']} automated checks pass; "
         f"{n_missing} claims still lack a second source"],
        "The comparisons on this slide are read from verdict.md's literature-ceiling section and "
        "from the shared evidence ledger (verify/data/rt_edge_anchors.yaml). The convention-"
        "matched comparison uses Reischle's IRF-deconvolved but background-included value, "
        "because the sweep's pulsed g2 is likewise background-included and never IRF-convolved. "
        "The 230 K Chatzarakis value is an optically pumped InAs/GaAs measurement, so it is "
        "context, not a like-for-like test. verify/verify_rt_edge_papers.py is a self-contained "
        "hallucination check: it re-derives every anchored claim from the ledger, confirms each "
        "has a real DOI or figure locator, checks units and magnitude, and distinguishes raw "
        "from background-corrected numbers. It currently passes "
        f"{evidence['passed']}/{evidence['total']} checks. {n_missing} claims still lack the "
        "required second, independent verified source and are reported as incomplete rather "
        "than silently accepted: " +
        "; ".join(f"{claim} ({reason})" for claim, reason in verdict["missing_claims"]) + ". " +
        f"For context, the 80 K claim's only source is {cite('reischle2008_g2_80K')}."))

    slides.append(Slide(
        "The design GUI",
        ["Engineers explore designs interactively -- sliders for the dot, the drive, the waveguide",
         "The same evaluate() pipeline runs live behind the GUI as behind this sweep",
         f"GUI readouts are card-level numbers: {STATIC_LABEL}, not the finite-pulse headline"],
        "fsim_gui/designer.py is a live GUI (Dear PyGui) over the identical fsim_core.device "
        "evaluator used by the acceptance sweep and by verify/verify_rt_edge_cards.py: moving a "
        "slider re-runs the full materials-to-g2 pipeline and redraws the results panel and g2(0) "
        "curve immediately, so a designer can explore the same trade-offs this deck reports. The "
        "screenshot on this slide (out/rt_edge/gui-smoke.png) is captured from a running instance "
        "of the tool via fsim_gui/designer.py --frames 10 --screenshot, not a mockup. One caution: "
        "the GUI evaluates a design card as shipped, and both edge cards ship with the static "
        "per-pulse loading model (drive.finite_pulse: false). Every g2 the GUI shows for them is "
        f"therefore a card-level value, labelled {STATIC_LABEL}; the headline finite-pulse model "
        "that the verdict gates on gives g2 near 1 at the same points. Use the GUI to explore "
        "trends, and the verdict for the answer.", image=GUI_PNG))

    slides.append(Slide(
        "The honest verdict",
        [f"1. Flux clears the floor at {elig_T} K: {m['flux_max']} > {floor:g} photons/s",
         f"2. But pulsed g2(0) = {m['g2_min']}: the dot is re-excited during the {pulse_ps:.0f} ps pulse",
         f"3. Same corner, static loading: g2 = {static_g2:.4f}, {STATIC_LABEL}",
         f"4. No Purcell recovery: +{pct:.1f}% total rate from a planar cavity",
         f"5. VERDICT: FAIL ({fail_reasons}) at commit {hashes['verdict']}"],
        "This is the verdict chain. First, the corrected collection model (audit H4) lifts the "
        f"collected flux above the floor, but only at {elig_T} K and only on {elig_cards}. "
        f"Second, at exactly those rows the pulsed g2 is {m['g2_min']}. Third, the cause: at "
        "these temperatures carriers escape the dot quickly, so during a "
        f"{pulse_ps:.0f} ps current pulse the dot empties and refills and can emit more than "
        "once per pulse (fsim_core.pulse_counting, which keeps re-excitation inside the pulse). "
        "The figure compares the two models from sweep.csv: at the best corner the legacy "
        f"static-loading model gives {static_g2:.4f} with {static_flux:.0f} photons/s, below the "
        f"floor (a card-level, {STATIC_LABEL} number; that model has "
        f"{static_sens.get('eligible', 'n/a')}/{static_sens.get('total', 'n/a')} eligible rows). "
        f"Fourth, a planar cavity adds only +{pct:.1f}% to the emission rate, so it cannot outrun "
        f"the escape. Fifth, the verdict fails on {fail_reasons}. {commit_line}. Closing the "
        "evidence gaps alone would not flip this sweep to PASS.",
        image=figures["reexc"], verbatim=[verdict["verdict_line"], commit_line]))

    slides.append(Slide(
        "Limitations and what would have to change",
        [f"Re-excitation within the {pulse_ps:.0f} ps pulse sets g2: slower escape or a shorter pulse is needed",
         "Filter window follows the linewidth, fixing t_X=0.5 by convention",
         "Residual background is transferred from 80 K data to this 300 K platform",
         "Pulsed drive is required: the CW dip is narrower than detector response"],
        "What would have to be true for a future revision to pass: (1) much slower thermal "
        "escape, from a deeper confinement barrier, so that the dot is not emptied and refilled "
        f"within one {pulse_ps:.0f} ps pulse, or a drive pulse short compared with the escape "
        "time; (2) a measured, not class-proxy, 300 K single-dot linewidth for an InP/GaAsP or "
        "InP/GaInP dot; (3) a measured X-XX splitting, since the sweep only brackets it; and "
        f"(4) an independent second source for each of the {n_missing} evidence claims still "
        "open in out/rt_edge/evidence.json. The modelling conventions also matter: the filter "
        "window tracks the linewidth, so t_X is 0.5 by construction and flux cannot depend on "
        "gamma300; the residual background is borrowed from an 80 K device; and because the CW "
        "antibunching dip is narrower than any detector's timing response, only pulsed, gated "
        "operation can demonstrate antibunching. None of these are engineering details the "
        "simulator can resolve on its own -- they are measurements and material developments "
        "that have not yet been published."))

    return slides


# ===================================================================
# python-pptx rendering
# ===================================================================

def fit_picture(slide, image_path: Path, left: Emu, top: Emu, max_w: Emu, max_h: Emu):
    pic = slide.shapes.add_picture(str(image_path), left, top, width=max_w)
    if pic.height > max_h:
        pic._element.getparent().remove(pic._element)
        pic = slide.shapes.add_picture(str(image_path), left, top, height=max_h)
    pic.left = int(left + (max_w - pic.width) / 2)
    pic.top = int(top + (max_h - pic.height) / 2)
    return pic


def set_title(slide, text: str):
    title = slide.shapes.title
    title.left, title.top = Inches(0.5), Inches(0.25)
    title.width, title.height = Inches(SLIDE_W_IN - 1.0), Inches(1.0)
    title.text_frame.text = text
    run = title.text_frame.paragraphs[0].runs[0]
    run.font.size = Pt(36)
    run.font.bold = True
    run.font.color.rgb = RGBColor(0x1B, 0x1F, 0x24)
    return title


def add_bullets(slide, bullets, left, top, width, height, font_pt=20):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    for i, text in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = f"• {text}"
        p.font.size = Pt(font_pt)
        p.font.color.rgb = RGBColor(0x1B, 0x1F, 0x24)
        p.space_after = Pt(10)
    return box


def build_pptx(slides: list[Slide], path: Path):
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W_IN)
    prs.slide_height = Inches(SLIDE_H_IN)
    layout = prs.slide_layouts[5]   # "Title Only"

    for s in slides:
        slide = prs.slides.add_slide(layout)
        set_title(slide, s.title)

        content_top = Inches(1.35)
        # P10: a verbatim data block (the full VERDICT line, ~700 characters,
        # plus the commit line) gets its own full-width strip at the bottom
        # of the slide in small monospace, so it never overflows a column.
        mono_h = Inches(1.25) if s.verbatim else 0
        content_h = Inches(SLIDE_H_IN - 1.35 - 0.3) - mono_h
        if s.image is not None:
            img_w = Inches(6.6)
            fit_picture(slide, s.image, Inches(SLIDE_W_IN - 0.5 - 6.6), content_top,
                       img_w, content_h)
            bullets_w = Inches(5.7)
        else:
            bullets_w = Inches(SLIDE_W_IN - 1.6)
        bullets_left = Inches(0.6) if s.image is not None else Inches(0.8)
        n_chars = sum(len(b) for b in s.bullets)
        if s.image is not None:
            font_pt = 18 if n_chars < 300 else (16 if n_chars < 380 else 14)
        else:
            font_pt = 22 if n_chars < 420 else 18
        add_bullets(slide, s.bullets, bullets_left, content_top, bullets_w, content_h,
                   font_pt=font_pt)
        if s.verbatim:
            box = slide.shapes.add_textbox(Inches(0.5), content_top + content_h + Inches(0.1),
                                           Inches(SLIDE_W_IN - 1.0), mono_h - Inches(0.1))
            tf = box.text_frame
            tf.word_wrap = True
            for j, text in enumerate(s.verbatim):
                p = tf.paragraphs[0] if j == 0 else tf.add_paragraph()
                p.text = text
                p.font.name = "Courier New"
                p.font.size = Pt(8)
                p.font.color.rgb = RGBColor(0x57, 0x60, 0x6A)

        notes = slide.notes_slide.notes_text_frame
        notes.text = s.notes

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    prs.save(str(path))


# ===================================================================
# Self-contained HTML rendering
# ===================================================================

def data_uri(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


HTML_HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>qd-photon-sim: rt-edge-emitter</title>
<style>
:root {
  --bg: #ffffff; --fg: #1b1f24; --muted: #57606a; --card: #f6f8fa;
  --accent: #1f6feb; --border: #d8dee4; --mono-bg: #eef1f4;
}
@media (prefers-color-scheme: dark) {
  :root { --bg: #0d1117; --fg: #e6edf3; --muted: #9198a1; --card: #161b22;
          --accent: #58a6ff; --border: #30363d; --mono-bg: #1c2128; }
}
* { box-sizing: border-box; }
html, body { margin: 0; padding: 0; background: var(--bg); color: var(--fg);
             font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif;
             overflow-x: hidden; }
section { max-width: 980px; margin: 0 auto; padding: 3rem 1.5rem;
          border-bottom: 1px solid var(--border); }
section:last-child { border-bottom: none; }
h1 { font-size: 2.2rem; margin: 0 0 1rem; }
h2 { font-size: 1.8rem; margin: 0 0 1.2rem; color: var(--fg); }
.subtitle { color: var(--muted); font-size: 1.1rem; }
ul { font-size: 1.15rem; line-height: 1.6; padding-left: 1.3rem; }
li { margin-bottom: 0.5rem; }
img { display: block; max-width: 100%; height: auto; margin: 1.2rem auto;
      border-radius: 6px; border: 1px solid var(--border); }
.notes { margin-top: 1.5rem; padding: 1rem 1.2rem; background: var(--card);
         border-radius: 8px; font-size: 0.95rem; color: var(--muted); line-height: 1.55; }
.notes summary { cursor: pointer; color: var(--accent); font-weight: 600; }
.mono { font-family: Consolas, Menlo, monospace; background: var(--mono-bg);
        padding: 0.6rem 0.9rem; border-radius: 6px; white-space: pre-wrap;
        overflow-wrap: anywhere; font-size: 0.85rem; margin: 1rem 0; }
.slideno { color: var(--muted); font-size: 0.85rem; margin-bottom: 0.4rem; }
</style>
</head>
<body>
"""

HTML_TAIL = """
</body>
</html>
"""


def build_html(slides: list[Slide], path: Path):
    parts = [HTML_HEAD]
    total = len(slides)
    for i, s in enumerate(slides, start=1):
        parts.append(f'<section id="slide-{i}">')
        parts.append(f'<div class="slideno">Slide {i} / {total}</div>')
        parts.append(f"<h2>{html_lib.escape(s.title)}</h2>")
        if s.image is not None:
            parts.append(f'<img src="{data_uri(s.image)}" alt="{html_lib.escape(s.title)}">')
        parts.append("<ul>")
        for b in s.bullets:
            parts.append(f"<li>{html_lib.escape(b)}</li>")
        parts.append("</ul>")
        if s.verbatim:
            for v in s.verbatim:
                parts.append(f'<div class="mono">{html_lib.escape(v)}</div>')
        parts.append("<details class=\"notes\"><summary>Speaker notes</summary><div>")
        parts.append(html_lib.escape(s.notes).replace("\n", "<br>"))
        parts.append("</div></details>")
        parts.append("</section>")
    parts.append(HTML_TAIL)
    path.write_text("".join(parts), encoding="utf-8")


# ===================================================================
# Main
# ===================================================================

def regenerate_gui_screenshot():
    GUI_PNG.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(ROOT / "fsim_gui" / "designer.py"), "--frames", "10",
           "--screenshot", str(GUI_PNG)]
    subprocess.run(cmd, cwd=str(ROOT), check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # P10: regenerating the GUI screenshot writes out/rt_edge/gui-smoke.png,
    # which is outside this script's own output directory; it is now opt-in
    # (--gui). --no-gui is still accepted (the default behaviour) so older
    # command lines keep working.
    parser.add_argument("--gui", action="store_true",
                       help="regenerate out/rt_edge/gui-smoke.png first (opt-in)")
    parser.add_argument("--no-gui", action="store_true",
                       help="use the committed out/rt_edge/gui-smoke.png as-is (default)")
    args = parser.parse_args()

    if args.gui and not args.no_gui:
        regenerate_gui_screenshot()

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    verdict = parse_verdict_md(VERDICT_MD)
    sweep_rows = load_sweep_rows(SWEEP_CSV)
    headline_rows = headline_model_rows(sweep_rows, verdict)
    anchors = load_anchor_sources(ANCHORS_YAML)
    anchor_vals = anchor_values(ANCHORS_YAML)
    evidence = evidence_numbers(EVIDENCE_JSON)

    lv = gaasp_band_levels()
    mat_sw = material_switch_numbers()
    diode_info = diode_numbers(headline_rows)
    wg_info = waveguide_numbers(headline_rows, GAASP_CARD)
    lw_info = linewidth_numbers()
    cmp = best_row_comparison(sweep_rows, verdict)
    planar = planar_cavity_numbers()
    hashes = git_hashes()

    figures = {
        "band": FIG_DIR / "band_diagram.png",
        "pin": FIG_DIR / "pin_diode.png",
        "overlap": FIG_DIR / "spectral_overlap.png",
        "wg": FIG_DIR / "waveguide_mode.png",
        "g2": FIG_DIR / "g2_dip.png",
        "linewidth": FIG_DIR / "linewidth_sweep.png",
        "flux": FIG_DIR / "flux_vs_T.png",
        "reexc": FIG_DIR / "reexcitation_g2.png",
    }
    draw_band_diagram(figures["band"], lv)
    draw_pin_diode(figures["pin"], diode_info["diode"], diode_info["V_j"], diode_info["I_uA"])
    draw_spectral_overlap(figures["overlap"], lw_info["params"], 4.0, 7.0)
    draw_waveguide_mode(figures["wg"], wg_info["stack"], wg_info["mode"], wg_info["edge"])
    draw_g2_dip(figures["g2"], float(verdict["metrics"]["g2_min"]),
                anchor_vals["laferriere23-g2-temperature"])
    linewidth_sweep = draw_linewidth_sweep(figures["linewidth"], headline_rows)
    draw_flux_vs_T(figures["flux"], headline_rows, verdict["flux_floor"])
    draw_reexcitation(figures["reexc"], headline_rows, static_model_rows(sweep_rows),
                      cmp["best"]["card_id"])
    brightness = verdict_brightness_factors()

    slides = build_slides(verdict, headline_rows, lv, mat_sw, diode_info, wg_info, lw_info,
                          evidence, anchors, figures, linewidth_sweep, brightness,
                          cmp, planar, hashes, anchor_vals)

    # Speaker-notes length rule of the deck (120-250 words per slide).
    for i, s in enumerate(slides, start=1):
        n_words = len(s.notes.split())
        if not 120 <= n_words <= 250:
            print(f"WARNING slide {i} notes have {n_words} words (rule: 120-250)")

    build_pptx(slides, PPTX_PATH)
    build_html(slides, HTML_PATH)
    print(f"slides: {len(slides)}")
    print(f"wrote {PPTX_PATH}")
    print(f"wrote {HTML_PATH}")


if __name__ == "__main__":
    main()
