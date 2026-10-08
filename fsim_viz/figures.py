"""fsim-viz: figure factory. Consumes fsim-core results objects only (three-layer
rule) and writes, for every figure, the underlying CSV alongside PDF/SVG/PNG --
the group re-renders in Origin, so data export is first-class.
"""
from __future__ import annotations

import csv
import json
import textwrap
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects  # noqa: E402  (label halos)
import matplotlib.pyplot as plt
import numpy as np

import fsim_theme
from fsim_core.card import Tag
from fsim_core.fitting import FitResult, V_A_TOL
from fsim_core.integrator import g2_of_T, solve_Tc

# Provenance is metadata, not hue (studio-02 / 02-charts section 3): tags are
# drawn as line-form chips ([V] solid, [DR] solid+inset, [E] outline, [A]
# dashed outline). The name is kept for importers; it now maps onto the
# ordinal tag ramp (V darkest -> A lightest) instead of green/amber/red.
_ORD = fsim_theme.load_tokens()["color"]["light"]["ordinal_tag"]
TAG_COLORS = {Tag.V: _ORD[0], Tag.DR: _ORD[1], Tag.E: _ORD[2], Tag.A: _ORD[3]}


# ------------------------------------------------------------------ theme helpers

def _rgba_tuple(css: str) -> tuple:
    """'rgba(82,81,78,0.10)' -> matplotlib RGBA tuple."""
    r, g, b, a = (float(v) for v in css[css.index("(") + 1:css.index(")")].split(","))
    return (r / 255.0, g / 255.0, b / 255.0, a)


def _palette(mode: str = "light") -> SimpleNamespace:
    tok = fsim_theme.load_tokens()
    c = tok["color"][mode]
    m = tok["mark"]
    return SimpleNamespace(
        s=c["series"], ink1=c["ink-1"], ink2=c["ink-2"], muted=c["muted"], ref=c["ref"],
        ref_wash=_rgba_tuple(c["ref-wash"]), surface=c["surface"], grid=c["grid"],
        status=c["status"], div=c["div"], ordinal=c["ordinal_tag"],
        band=m["band_opacity"], band_overlap=m["overlap_band_opacity"],
        edge_w=m["edge_width"], edge_a=m["edge_opacity"], ref_w=0.8)


@contextmanager
def _themed(mode: str = "light"):
    """Apply the fsim matplotlib style for one figure function only: the
    rcParams are restored on exit, so importing/using fsim_viz never leaks
    style into a caller's own plots."""
    with plt.rc_context():
        fsim_theme.apply_matplotlib(mode)
        yield _palette(mode)


def _ref_h(ax, y, label, P, *, va="bottom"):
    """Threshold reference: 1px `ref` ink with a direct right-edge label."""
    ax.axhline(y, color=P.ref, lw=P.ref_w, zorder=1.5)
    if label:
        ax.annotate(label, xy=(1.0, y), xycoords=("axes fraction", "data"),
                    xytext=(-3, 2 if va == "bottom" else -2), textcoords="offset points",
                    ha="right", va=va, fontsize=8, color=P.ink2)


def _ref_v(ax, x, label, P):
    """Vertical reference (T_c, T_target, ...): `ref` ink rule, label at top."""
    ax.axvline(x, color=P.ref, lw=P.ref_w, zorder=1.5)
    if label:
        ax.annotate(label, xy=(x, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                    fontsize=8, color=P.ink1)


def _band(ax, x, lo, hi, color, P, *, overlap=False, label=None, edge_ls=("-", "-"),
          edge_labels=(None, None)):
    """Envelope rule: fill in the series hue at 12% (8% when bands overlap)
    with lo/hi edge lines at 1.25 pt / 60% -- never an unlabeled blend."""
    ax.fill_between(x, lo, hi, color=color, alpha=P.band_overlap if overlap else P.band,
                    lw=0, label=label)
    ax.plot(x, lo, color=color, lw=P.edge_w, alpha=P.edge_a, ls=edge_ls[0], label=edge_labels[0])
    ax.plot(x, hi, color=color, lw=P.edge_w, alpha=P.edge_a, ls=edge_ls[1], label=edge_labels[1])


def _tc_interval(ax, lo, hi, P, *, unit="K"):
    """T_c as an interval: ref-wash span labelled 'T_c lo-hi K', or the
    'not crossed in range' label when either end is NaN."""
    if lo == lo and hi == hi:
        ax.axvspan(lo, hi, color=P.ref_wash, lw=0, zorder=0.5)
        ax.annotate(f"T$_c$ {lo:.0f}–{hi:.0f} {unit}", xy=(0.5 * (lo + hi), 1.0),
                    xycoords=("data", "axes fraction"), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom", fontsize=8,
                    color=P.ink1)
    else:
        ax.annotate("T$_c$ not crossed in range", xy=(0.99, 0.97), xycoords="axes fraction",
                    ha="right", va="top", fontsize=8, color=P.ink2)


_CHIP_STYLE = {  # line-form grammar: face, edge, linestyle, text colour key
    "V": ("ink1", "ink1", "-", "surface"),
    "DR": ("ink2", "ink1", "-", "surface"),
    "E": (None, "ink2", "-", "ink1"),
    "A": (None, "ink2", "--", "ink1"),
}


def _tag_of(tag) -> str:
    return tag.name if isinstance(tag, Tag) else str(tag).strip("[]").upper()


def _chip(fig, x, y, text, tag, P, *, ha="left"):
    fc, ec, ls, tc = _CHIP_STYLE[_tag_of(tag)]
    fig.text(x, y, text, ha=ha, va="bottom", fontsize=8.5, family="monospace",
             color=getattr(P, tc),
             bbox=dict(boxstyle="square,pad=0.35", fc=getattr(P, fc) if fc else "none",
                       ec=getattr(P, ec), ls=ls, lw=0.9))


def _status_chip(fig, x, y, ok: bool, word: str, P, *, ha="right"):
    """Verdicts use the status pair with an icon AND a word, never hue alone."""
    col = P.status["pass"] if ok else P.status["fail"]
    fig.text(x, y, ("✓ " if ok else "✗ ") + word, ha=ha, va="bottom", fontsize=8.5,
             color=P.ink1, fontweight=600,
             bbox=dict(boxstyle="square,pad=0.35", fc=(*matplotlib.colors.to_rgb(col), 0.14),
                       ec=col, lw=1.0))


def _status_ax(ax, ok: bool, word: str, P):
    """Axes-level verdict chip (icon + word), top-right above the plot."""
    col = P.status["pass"] if ok else P.status["fail"]
    ax.text(1.0, 1.02, ("✓ " if ok else "✗ ") + word, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.5, fontweight=600, color=P.ink1,
            bbox=dict(boxstyle="square,pad=0.3", fc=(*matplotlib.colors.to_rgb(col), 0.14),
                      ec=col, lw=1.0))


def _title(fig, text, tag, P, *, y=1.0, width_chars=None, status=None):
    """Figure title in ink (wrapped to the figure width) with the tag-chain
    chip on the line above it, left; an optional (ok, word) verdict chip on
    the right of the same line."""
    if width_chars is None:
        width_chars = int(fig.get_figwidth() * 11.5)
    lines = textwrap.wrap(text, width_chars) or [""]
    fig.suptitle("\n".join(lines), fontsize=10, color=P.ink1, y=y, va="bottom")
    chip_y = y + (0.175 * len(lines) + 0.06) / fig.get_figheight()  # inches -> fig frac
    _chip(fig, 0.01, chip_y, f"tag chain [{_tag_of(tag)}]", tag, P)
    if status is not None:
        _status_chip(fig, 0.99, chip_y, status[0], status[1], P)


def _data_marks(P):
    """Measured ([V]) points: ink markers with a surface ring + error bars."""
    return dict(color=P.ink1, mfc=P.ink1, mec=P.surface, mew=1.0, ecolor=P.ink2,
                elinewidth=0.9, capsize=3, zorder=5)


def _write_csv(path: Path, header: list[str], rows) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        wtr = csv.writer(f)
        wtr.writerow(header)
        wtr.writerows(rows)


def phase0_bundle(fit: FitResult, card_path: Path, outdir: Path) -> dict:
    """One-click report bundle: figure (pdf/svg/png) + CSV per panel + fitted
    params + the exact card, so the plot is regenerable from the bundle."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    Ts_d = np.array([r["T"] for r in fit.data])
    g2_d = np.array([r["g2"] for r in fit.data])
    err_d = np.array([r["err"] for r in fit.data])
    upper = np.array([r.get("bound") == "upper" for r in fit.data])
    ws = np.array([r["w"] for r in fit.data])
    dxs = np.array([r["dx"] for r in fit.data])

    # model curve with the published windows interpolated in T
    p = dict(fit.params)
    p["w"] = lambda T: np.interp(T, Ts_d, ws)
    p["dx"] = lambda T: np.interp(T, Ts_d, dxs)
    Ts = np.linspace(min(Ts_d) - 18, max(Ts_d) + 60, 300)
    pts = [g2_of_T(T, p) for T in Ts]
    g2_m = np.array([q.g2 for q in pts])
    eps_m = np.array([q.eps for q in pts])
    rho_m = np.array([q.rho for q in pts])
    gam_m = np.array([q.gamma for q in pts])
    Tc = solve_Tc(p)

    # ---- CSV exports (one per plotted series group)
    _write_csv(outdir / "fit_data_points.csv",
               ["T_K", "g2", "err", "bound", "w_meV", "dx_meV", "model_g2", "residual"],
               ((r["T"], r["g2"], r["err"], r.get("bound", "value"), r["w"], r["dx"], m, res)
                for r, m, res in zip(fit.data, fit.model_g2, fit.residuals)))
    _write_csv(outdir / "fit_model_curve.csv",
               ["T_K", "g2_model", "eps", "rho", "gamma_meV",
                "w_interp_meV_extrap_clamped", "dx_interp_meV_extrap_clamped"],
               zip(Ts, g2_m, eps_m, rho_m, gam_m,
                   np.interp(Ts, Ts_d, ws), np.interp(Ts, Ts_d, dxs)))
    with open(outdir / "fit_params.json", "w", encoding="utf-8") as f:
        json.dump({"params": fit.params, "Tc_K": None if np.isnan(Tc) else Tc,
                   "tag_chain": fit.tag.label, "passed_pm0.03": fit.passed,
                   "cost": fit.cost, "notes": fit.notes}, f, indent=2)
    (outdir / "card_snapshot.yaml").write_text(
        Path(card_path).read_text(encoding="utf-8"), encoding="utf-8")

    # ---- figure
    with _themed("light") as P:
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

        ax = axes[0]
        val = ~upper
        ax.errorbar(Ts_d[val], g2_d[val], yerr=err_d[val], fmt="o", ms=6,
                    label="data (Fig. 5) [V]", **_data_marks(P))
        if upper.any():
            ax.errorbar(Ts_d[upper], g2_d[upper], yerr=err_d[upper], fmt="v", ms=7,
                        uplims=True, label="upper bound (≤)", **_data_marks(P))
        ax.plot(Ts, g2_m, color=P.s[0], lw=2, label="F-series model")
        _ref_h(ax, 0.5, "$g^{(2)}(0)$ = 0.5 ceiling", P)
        if np.isfinite(Tc):
            _ref_v(ax, Tc, f"T$_c$ = {Tc:.0f} K", P)
        ax.set_xlabel("T (K)")
        ax.set_ylabel("$g^{(2)}(0)$")
        ax.set_ylim(0, max(1.0, float(np.nanmax(g2_m)) * 1.02))
        ax.legend(fontsize=8, loc="upper left")
        ax.set_title("V-a fit: $g^{(2)}(T)$, published windows", fontsize=10, pad=14)

        ax = axes[1]
        ax.plot(Ts, eps_m, color=P.s[1], lw=2, label=r"$\varepsilon(T)=t_{XX}/t_X$")
        ax.plot(Ts, rho_m**2, color=P.s[2], lw=2, label=r"$\rho(T)^2$")
        ax.plot(Ts, rho_m**2 * (1 - eps_m), color=P.ink2, lw=1.4, ls="--",
                label=r"$\rho^2(1-\varepsilon)$")
        _ref_h(ax, 0.5, "1/2 master ceiling", P)
        ax.set_xlabel("T (K)")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8)
        ax.set_title("decomposition (master ceiling at 1/2)", fontsize=10, pad=14)

        ax = axes[2]
        ax.plot(Ts, gam_m, color=P.s[0], lw=2, label="model Γ(T)")
        gA, gAv, gAt, T_anchor = fit.gamma_anchor
        ax.errorbar([T_anchor], [gAv], yerr=[gAt], fmt="s", ms=6,
                    label="published anchor [V]", **_data_marks(P))
        ax.plot([T_anchor], [gA], "o", mfc="none", mec=P.s[0], ms=10, mew=1.8,
                label="model @ anchor", zorder=6)
        ax.set_xlabel("T (K)")
        ax.set_ylabel(r"$\Gamma$ (meV)")
        ax.legend(fontsize=8)
        ax.set_title(r"$\Gamma(T)$ vs anchor", fontsize=10, pad=14)

        status = "PASS" if fit.passed else "FAIL"
        fig.tight_layout()
        _title(fig, f"{Path(card_path).stem}: V-a fit [{status} at ±{V_A_TOL}]   "
                    f"tag chain {fit.tag.label}" + ("  — PLACEHOLDER DATA" if fit.notes else ""),
               fit.tag, P, y=1.0, status=(fit.passed, f"{status} at ±{V_A_TOL}"))
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"phase0_fit.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)

    return {"outdir": str(outdir), "Tc": Tc}


def phase1_bundle(drive: dict, thermal: dict, outdir: Path) -> dict:
    """Phase-1 report bundle: F1b drive-factor curve + Delta T_J envelope maps
    (per-figure CSVs, per the plan). Inputs are plain result dicts computed by
    fsim-core callers (three-layer rule: no physics here)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mus = drive["mu"]
    _write_csv(outdir / "drive_factor.csv",
               ["mu"] + [f"f_eps={e}" for e in drive["curves"]],
               zip(mus, *drive["curves"].values()))
    d_um = thermal["diam_um"]
    for tpl, per_I in thermal["templates"].items():
        _write_csv(outdir / f"tj_map_{tpl.lower().replace('/', '_')}.csv",
                   ["diameter_um"] + [f"dTj_K_{I}uA_lo,dTj_K_{I}uA_hi".split(",")[i]
                                      for I in thermal["currents_uA"] for i in (0, 1)],
                   zip(d_um, *[col for I in thermal["currents_uA"]
                               for col in (per_I[I][0], per_I[I][1])]))

    with _themed("light") as P:
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

        ax = axes[0]
        for i, (e, f) in enumerate(drive["curves"].items()):
            ax.plot(mus, f, lw=2, color=P.s[i % len(P.s)], label=rf"$\varepsilon$ = {e}")
        ax.axvline(drive["mu_op"], color=P.ref, lw=P.ref_w)
        ax.annotate(rf"$\mu_{{op}}$ = {drive['mu_op']}", (drive["mu_op"], 1.85),
                    xytext=(3, 0), textcoords="offset points",
                    fontsize=8, ha="left", color=P.ink1)
        _ref_h(ax, 2.0, "2", P, va="top")
        ax.set_xscale("log")
        ax.set_xlabel(r"mean loading $\mu$")
        ax.set_ylabel(r"$g^{(2)}(\mu)\,/\,\varepsilon$")
        ax.set_ylim(0.95, 2.1)
        ax.legend(fontsize=8)
        ax.set_title("F1b: finite-$\\mu$ drive penalty (cap-2)", fontsize=10)

        for ax, tpl in zip(axes[1:], thermal["templates"]):
            per_I = thermal["templates"][tpl]
            for k, I in enumerate(thermal["currents_uA"]):
                color = P.s[k % len(P.s)]
                lo, hi = per_I[I]
                lo = np.where(np.isinf(lo), np.nan, lo)  # runaway region: gap, not a line
                hi = np.where(np.isinf(hi), np.nan, hi)
                # several current bands share one axis -> overlap rule (8% + edges)
                _band(ax, d_um, lo, hi, color, P, overlap=True,
                      label=f"{I} µA" + (" (runaway ←)" if np.isnan(hi).any() else ""))
            ax.set_xlabel("mesa diameter (µm)")
            ax.set_ylabel(r"$\Delta T_J$ (K)")
            ax.set_yscale("log")
            ax.legend(fontsize=8, title="CW drive (band: epi-k range)", title_fontsize=8)
            ax.set_title(f"{tpl}: $T_J-T_{{hs}}$ @ {thermal['T_hs']:.0f} K "
                         f"(band: epi-k range)", fontsize=10)

        fig.tight_layout()
        _title(fig, "Phase 1 — Module D drive penalty and Module A junction-heating envelopes "
                    "(tag chain [A]: requirement envelopes, not predictions)", "A", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"phase1_drive_thermal.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)
    return {"outdir": str(outdir)}


def vb_figure(vb: dict, outdir: Path) -> None:
    """V-b figure: measured Laferriere g2(T) points vs the reachable envelope of
    the F1 spectral model over the swept (unpublished) parameters."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    Ts = [77.0, 220.0, 300.0]
    lo = [min(vb["reach"][T][d][0] for d in vb["reach"][T]) for T in Ts]
    hi = [max(vb["reach"][T][d][1] for d in vb["reach"][T]) for T in Ts]
    _write_csv(outdir / "vb_envelope.csv",
               ["T_K", "g2_measured", "err", "reachable_lo", "reachable_hi"],
               ((T, r["g2"], r["err"], l, h)
                for T, r, l, h in zip(Ts, vb["data"], lo, hi)))

    with _themed("light") as P:
        fig, ax = plt.subplots(figsize=(5.6, 3.9))
        _band(ax, Ts, lo, hi, P.s[0], P, label="F1 reachable envelope (swept [E] inputs)")
        ax.errorbar(Ts, [r["g2"] for r in vb["data"]], yerr=[r["err"] for r in vb["data"]],
                    fmt="o", ms=7, label="measured (Fig. 5) [V]", **_data_marks(P))
        _ref_h(ax, 0.5, "Theorem-0 saturated bound (ε→1, μ→∞)", P)
        ax.set_xlabel("T (K)")
        ax.set_ylabel("$g^{(2)}(0)$")
        ax.set_ylim(0, 1.0)
        ax.legend(fontsize=8, loc="upper left")
        ok = bool(vb["covered_delta"])
        ax.set_title("V-b (Laferrière): ε→1 limit under published windows", fontsize=10,
                     pad=22)
        _status_ax(ax, ok, "CONSISTENT" if ok else "NOT COVERED", P)
        fig.tight_layout()
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"vb_laferriere.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def cavity_design_figure(cv: dict, outdir: Path) -> None:
    """Phase-2 cavity deliverable: tracking detuning, cavity-filter eps map,
    and F_eff/F_P vs kappa."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _write_csv(outdir / "cavity_tracking_detuning.csv",
               ["T_K", "detuning_meV"], zip(cv["Ts"], cv["det_meV"]))
    _write_csv(outdir / "cavity_eps_map.csv",
               ["kappa_meV"] + [f"eps_delta_{d:.1f}" for d in cv["deltas"]],
               ((k, *row) for k, row in zip(cv["kappas"], cv["eps_band"])))

    with _themed("light") as P:
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
        ax = axes[0]
        ax.plot(cv["Ts"], cv["det_meV"], color=P.s[0], lw=2)
        _ref_h(ax, 0, "zero detuning", P)
        _ref_v(ax, cv["T_target"], f"$T_{{target}}$ = {cv['T_target']:.0f} K", P)
        ax.set_xlabel("T (K)")
        ax.set_ylabel("X–mode detuning (meV)")
        ax.set_title("mode-tracking rule (F6 iv): mode red of\ncryogenic line, zero at "
                     "$T_{target}$", fontsize=9, pad=14)

        ax = axes[1]
        eps_lo, eps_hi = cv["eps_band"].min(axis=1), cv["eps_band"].max(axis=1)
        _band(ax, cv["kappas"], eps_lo, eps_hi, P.s[0], P, edge_ls=("-", "--"),
              edge_labels=(rf"$\Delta$ = {cv['deltas'][-1]:.1f} meV",
                           rf"$\Delta$ = {cv['deltas'][0]:.1f} meV"))
        ax.set_xlabel(r"cavity $\kappa$ (meV)")
        ax.set_ylabel(r"$\varepsilon$ at $T_{target}$ (cavity filter)")
        ax.legend(fontsize=8)
        ax.set_title(rf"cavity-only $\varepsilon$; $\Gamma(T_t)$ = {cv['gam_t']:.1f} meV "
                     "[A proxy]", fontsize=9, pad=14)

        ax = axes[2]
        ax.plot(cv["kappas"], cv["Feff"], color=P.s[0], lw=2)
        ax.set_xlabel(r"cavity $\kappa$ (meV)")
        ax.set_ylabel(r"$F_{eff}/F_P = \kappa/(\kappa+\Gamma)$")
        ax.set_title("spectral-overlap Purcell penalty (F6 iii)", fontsize=9, pad=14)

        fig.tight_layout()
        _title(fig, "Phase 2 — Module B design rules (tag chain [A]: requirement envelopes; "
                    "κ, F_P, G await MEEP/COMSOL)", "A", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"phase2_cavity_design.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def phase3_packet_figure(m: dict, s: dict, a: dict, r: dict, outdir: Path) -> None:
    """Phase-3 packet: (ii) (Delta,rho) T_c map, (i) staged envelope,
    (iii) F5 aperture rules, (iv) measurement ranking. CSV per panel."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    _write_csv(outdir / "map_tc_delta_rho.csv",
               ["rho\\delta_meV"] + [f"{d:.2f}" for d in m["deltas"]],
               ((f"{rho:.3f}", *row) for rho, row in zip(m["rhos"], m["Tc"])))
    g_lo, g_hi = m["g_band"]
    _write_csv(outdir / "map_rho_required_300K.csv",
               ["delta_meV", f"rho_req_gamma{g_lo}", f"rho_req_gamma{g_hi}"],
               zip(m["deltas"], m["rho_req"][g_lo], m["rho_req"][g_hi]))
    for T_hs, e in s["envelopes"].items():
        _write_csv(outdir / f"staged_envelope_{T_hs:.0f}K.csv",
                   ["delta_meV", "g2_lo", "g2_hi"],
                   ((d, lo, hi) for d, (lo, hi) in zip(s["deltas"], e["band"])))
    _write_csv(outdir / "aperture_g2_penalty.csv",
               ["density_cm2\\diam_um"] + [f"{d:.2f}" for d in a["diams"]],
               ((f"{n:.2e}", *row) for n, row in zip(a["dens"], a["g2pen"])))
    _write_csv(outdir / "measurement_ranking.csv",
               ["input", "envelope_narrowing_K"],
               ((k, v) for k, v in r["ranked"]))

    with _themed("light") as P:
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

        ax = axes[0, 0]
        # sequential blue ramp in 50 K bands (no rainbow), ink contours with a halo
        cf = ax.contourf(m["deltas"], m["rhos"], m["Tc"], levels=np.arange(50, 450, 50),
                         cmap="fsim_seq")
        cs = ax.contour(m["deltas"], m["rhos"], m["Tc"], levels=[77, 120, 200],
                        colors=P.ink1, linewidths=0.9)
        lbl = ax.clabel(cs, fmt="%.0f K", fontsize=7, colors=P.ink1)
        cs300 = ax.contour(m["deltas"], m["rhos"], m["Tc"], levels=[300],
                           colors=P.ink1, linewidths=2.0)
        lbl += ax.clabel(cs300, fmt={300: "300 K target"}, fontsize=7.5, colors=P.ink1)
        for t in lbl:
            t.set_path_effects([matplotlib.patheffects.withStroke(linewidth=2.5,
                                                                   foreground=P.surface)])
        _band(ax, m["deltas"], m["rho_req"][g_lo], m["rho_req"][g_hi], P.s[1], P,
              edge_ls=("-", "--"),
              label=r"$\rho$ required for 300 K [$\Gamma$(300) 6–7 meV]")
        ax.axvline(5.0, color=P.ref, lw=P.ref_w)
        ax.annotate(">50% of (211)B dots", (5.0, 0.615), xytext=(3, 0),
                    textcoords="offset points", rotation=90, fontsize=7, color=P.ink1)
        for i, p in enumerate(m["points"]):
            ax.plot(p["delta"], p["rho"], "o", mfc=P.surface, mec=P.ink1, mew=1.2, ms=6,
                    zorder=6)
            # alternate leader-line offsets so neighbouring labels never collide
            dy = (14 if i % 2 == 0 else -16) + 4 * (i // 2 % 2)
            ax.annotate(p["name"], (p["delta"], p["rho"]), xytext=(10, dy),
                        textcoords="offset points", fontsize=6.5, color=P.ink1,
                        arrowprops=dict(arrowstyle="-", color=P.ink2, lw=0.6,
                                        shrinkA=0, shrinkB=3),
                        path_effects=[matplotlib.patheffects.withStroke(
                            linewidth=2.5, foreground=P.surface)], zorder=7)
        plt.colorbar(cf, ax=ax, label="$T_c$ (K)")
        ax.set_xlabel(r"$\Delta_{XX}$ (meV)")
        ax.set_ylabel(r"signal purity $\rho$")
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(False)
        ax.set_title("(ii) master-ceiling map: $T_c(\\Delta,\\rho)$, narrow-filter bound",
                     fontsize=9)

        ax = axes[0, 1]
        for k, (T_hs, e) in enumerate(s["envelopes"].items()):
            band = e["band"]
            _band(ax, s["deltas"], band[:, 0], band[:, 1], P.s[k % len(P.s)], P,
                  overlap=len(s["envelopes"]) > 1,
                  label=f"$T_{{hs}}$ = {T_hs:.0f} K (band; upper edge = worst case)")
        _ref_h(ax, 0.5, "$g^{(2)}(0)$ = 0.5 ceiling", P)
        _ref_h(ax, 0.1, "$g^{(2)}(0)$ = 0.1", P)
        ax.set_xlabel(r"$\Delta_{XX}$ (meV)  [A: unmeasured on InP/GaAsP]")
        ax.set_ylabel("$g^{(2)}(0)$ envelope")
        ax.set_yscale("log")
        ax.legend(fontsize=7)
        ax.set_title("(i) staged device envelope, electrical, w=$\\Gamma$ convention",
                     fontsize=9)

        ax = axes[1, 0]
        cf = ax.contourf(a["diams"], a["dens"], a["g2pen"],
                         levels=[0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0], cmap="fsim_seq")
        cs = ax.contour(a["diams"], a["dens"], a["Nw"], levels=[1.0], colors=P.ink1,
                        linewidths=1.5)
        for t in ax.clabel(cs, fmt="$N_w$=%.0f", fontsize=7, colors=P.ink1):
            t.set_path_effects([matplotlib.patheffects.withStroke(linewidth=2.5,
                                                                   foreground=P.surface)])
        ax.set_yscale("log")
        ax.grid(False)
        plt.colorbar(cf, ax=ax, label="F5 aperture $g^{(2)}$ penalty")
        ax.set_xlabel("aperture diameter (µm)")
        ax.set_ylabel("QD density (cm$^{-2}$)")
        ax.set_title(f"(iii) F5 aperture/density rules (w = {a['w']:.1f} meV)", fontsize=9)

        ax = axes[1, 1]
        items = r["ranked"][::-1]
        vals = [v for _, v in items]
        ax.barh(["\n".join(textwrap.wrap(str(r["labels"][k]), 28)) for k, _ in items], vals,
                color=P.s[0], height=0.5)
        vmax = max(vals) if vals else 1.0
        for yi, v in enumerate(vals):
            ax.annotate(f"{v:.0f} K", (v, yi), xytext=(3, 0), textcoords="offset points",
                        va="center", ha="left", fontsize=7.5, color=P.ink2,
                        family="monospace")
        ax.set_xlim(0, vmax * 1.15 if vmax > 0 else 1.0)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel(f"$T_c$ envelope narrowing (K) out of {r['width_full']:.0f} K total")
        ax.set_title("(iv) measurement-priority ranking (in-house queue)", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)

        fig.tight_layout()
        _title(fig, "Phase 3 — requirement envelopes and design rules "
                    "(tag chain [A]: every panel inherits unmeasured inputs; "
                    "see the packet note for the chains)", "A", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"phase3_packet.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def spec_figure(staged_rows: list, k300_rows: list, be_curves: dict, outdir: Path) -> None:
    """Spec-mode deliverable (fsim_core.spec): 2 panels summarizing the
    design-target spec sheets.
      left:  kappa_max vs Delta_XX, one line per (T_op, target) staged group,
             plus the 300 K route's kappa_max vs Delta_XX per Gamma(300) anchor.
      right: G_required vs injection background b_e for a few representative
             design points (staged + 300 K), np.inf entries (rho_req >= 1,
             no finite gain closes the gap) simply stop the curve.
    Consumes fsim_core.spec results only (three-layer rule); every array here
    was computed in fsim_core/scripts, not derived in this module.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    _write_csv(outdir / "spec_kappa_vs_delta.csv",
               ["group", "delta_meV", "kappa_max_meV"],
               ((f"staged T_op={r['T_op']:.0f}K t={r['target_g2']:.2f}", r["delta_xx"], r["kappa_max"])
                for r in staged_rows))
    _write_csv(outdir / "spec_G_req_vs_be.csv",
               ["label", "b_e", "G_required"],
               ((label, be, g) for label, (bes, gs) in be_curves.items() for be, g in zip(bes, gs)))

    with _themed("light") as P:
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

        ax = axes[0]
        groups = {}
        for r in staged_rows:
            groups.setdefault((r["T_op"], r["target_g2"]), []).append(r)
        for (i, ((T_op, target), rs)) in enumerate(sorted(groups.items())):
            rs = sorted(rs, key=lambda r: r["delta_xx"])
            deltas = [r["delta_xx"] for r in rs]
            kmax = [r["kappa_max"] for r in rs]
            ax.plot(deltas, kmax, "o-", color=P.s[i % len(P.s)], lw=1.8, ms=4,
                    label=f"staged $T_{{op}}$={T_op:.0f} K, target={target:.2f}")
        k_groups = {}
        for k in k300_rows:
            k_groups.setdefault(k["gamma300"], []).append(k)
        for i, (g300, ks) in enumerate(sorted(k_groups.items())):
            ks = sorted(ks, key=lambda k: k["delta"])
            deltas = [k["delta"] for k in ks]
            kmax = [k["kappa_max"] for k in ks]
            ax.plot(deltas, kmax, "s--", color=P.ink2 if i else P.ink1, lw=1.4, ms=4,
                    label=rf"300 K, $\Gamma$(300)={g300:.1f} meV")
        ax.set_xlabel(r"$\Delta_{XX}$ (meV)")
        ax.set_ylabel(r"$\kappa_{max}$ (meV)  [cavity linewidth ceiling]")
        ax.legend(fontsize=7)
        ax.set_title("kappa ceiling vs splitting (F6 i inverted)", fontsize=9)

        ax = axes[1]
        for i, (label, (bes, gs)) in enumerate(sorted(be_curves.items())):
            finite = np.isfinite(gs)
            ax.plot(np.asarray(bes)[finite], gs[finite], lw=1.8,
                    color=P.s[i % len(P.s)], label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"injection background $b_e$ (signal units)")
        ax.set_ylabel(r"$G_{required}$ (collection gain)")
        ax.legend(fontsize=7)
        ax.set_title("cavity gain requirement vs background (F6 ii inverted)", fontsize=9)

        fig.tight_layout()
        _title(fig, "Spec mode — design-target requirement sheets (tag chain [A]: "
                    "necessary conditions from the validated model, not existence proofs "
                    "-- kappa/V~ achievability is MEEP's question, R_th/mesa is COMSOL's)",
               "A", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"spec_figure.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def zhao_fit_figure(rows: list, curve: dict, I_th: float, outdir: Path) -> None:
    """Module-D (fsim_core.sde) Tier-2 validation figure for cards/zhao.yaml:
    left panel is the headline sim-vs-published g2(0) comparison at I=4 I_th
    (normal/quiet), right panel is the g2 vs I/I_th curve for both pump
    modes with the I=4 I_th published points overlaid. CSVs are written by
    the caller (scripts/run_zhao_fit.py), matching the run_spec.py/
    spec_figure convention -- this function only plots. [E]-tag caption note:
    the simulation is fsim_core.sde's literature-class Tier-2 fit (see that
    module's and run_zhao_fit.py's docstrings for the tuned-parameter table
    and the honestly-reported quiet-pump near-miss -- do not read a visually
    small gap on this plot as a precision match; see the printed sigma
    table)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with _themed("light") as P:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        # model (sim) = filled slot colour; measured (published) = hollow ink marker
        colors = {"normal": P.s[0], "quiet": P.s[1]}
        pub_marker = {"normal": "o", "quiet": "s"}
        g2_max = 1.0

        ax = axes[0]
        for i, r in enumerate(rows):
            x = i
            ax.errorbar([x - 0.08], [r["sim_g2"]], yerr=[r["sim_se"]], fmt="o", ms=8,
                        color=colors[r["pump"]], capsize=4,
                        label=f"{r['pump']} (sim)" if i < 2 else None)
            ax.errorbar([x + 0.08], [r["published_g2"]], yerr=[r["published_err"]],
                        fmt=pub_marker.get(r["pump"], "s"), ms=8, color=P.ink1, mfc=P.surface,
                        mec=P.ink1, ecolor=P.ink2, capsize=4,
                        label=f"{r['pump']} (published) [V]" if i < 2 else None)
            g2_max = max(g2_max, r["sim_g2"] + r["sim_se"], r["published_g2"] + r["published_err"])
        _ref_h(ax, 1.0, "$g^{(2)}(0)$ = 1", P)
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([r["pump"] for r in rows])
        ax.set_ylabel("$g^{(2)}(0)$")
        ax.set_ylim(0, g2_max * 1.05)
        ax.set_title("I = 4 $I_{th}$: sim (filled) vs Zhao et al. (open)", fontsize=9.5)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, fontsize=7, loc="lower right")

        ax = axes[1]
        for label, mk in (("normal", "o-"), ("quiet", "s-")):
            c = curve[label]
            g2 = np.asarray(c["g2"])
            se = np.asarray(c["se"])
            ax.plot(c["I"], g2, mk, color=colors[label], lw=1.8, ms=5, label=f"{label} (sim)")
            ax.fill_between(c["I"], g2 - se, g2 + se, color=colors[label], alpha=0.10, lw=0,
                            label=f"{label} ±SE")
            g2_max = max(g2_max, float(np.nanmax(g2 + se)))
        for r in rows:
            ax.errorbar([4.0], [r["published_g2"]], yerr=[r["published_err"]],
                        fmt=pub_marker.get(r["pump"], "s"), ms=7, color=P.ink1, mfc=P.surface,
                        mec=P.ink1, ecolor=P.ink2, capsize=4, zorder=6)
        _ref_h(ax, 1.0, "$g^{(2)}(0)$ = 1", P)
        ax.set_xlabel(r"$I / I_{th}$")
        ax.set_ylabel("$g^{(2)}(0)$")
        ax.set_ylim(0, g2_max * 1.05)
        ax.legend(fontsize=7.5, loc="lower right")
        ax.set_title(f"g$^{{(2)}}$(0) vs pump strength ($I_{{th}}$ raw axis = {I_th:.3f})",
                     fontsize=9.5)

        fig.tight_layout()
        _title(fig, "cards/zhao.yaml Tier-2 fit -- fsim_core.sde stochastic drive engine "
                    "(tag chain [E]: literature-class tuned rate-equation model; "
                    "see run_zhao_fit.py for the honest verdict)", "E", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"zhao_fit_figure.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def vc_reischle_figure(vc: dict, outdir: Path) -> None:
    """V-c consistency figure: digitized rho(w) envelope vs the rho required by
    the measured g2 under the trion (eps=0) F-series prediction g2 = 1-rho^2."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ws, lo, hi = vc["w_meV"], vc["rho_lo"], vc["rho_hi"]
    _write_csv(outdir / "vc_rho_window.csv",
               ["w_meV", "rho_lo", "rho_hi", "rho_required", "rho_required_err"],
               ((w, l, h, vc["rho_req"], vc["rho_req_err"]) for w, l, h in zip(ws, lo, hi)))

    with _themed("light") as P:
        fig, ax = plt.subplots(figsize=(5.4, 3.8))
        _band(ax, ws, lo, hi, P.s[0], P,
              label=r"digitized $\rho(w)$ (zero-level conventions)")
        ax.axhspan(vc["rho_req"] - vc["rho_req_err"], vc["rho_req"] + vc["rho_req_err"],
                   color=P.ref_wash, lw=0)
        ax.axhline(vc["rho_req"], color=P.ref, lw=1.2,
                   label=rf"$\rho$ required by $g^{{(2)}}$ = {vc['g2']}$\pm${vc['g2_err']}")
        ax.set_xlabel("detection window $w$ (meV)")
        ax.set_ylabel(r"signal fraction $\rho$")
        ax.set_ylim(0.7, 1.0)
        ax.legend(fontsize=8, loc="lower left")
        verdict = "CONSISTENT" if vc["consistent"] else "INCONSISTENT"
        ax.set_title("V-c (Reischle, 100 MHz, ~40 K): trion, $\\varepsilon\\approx 0$",
                     fontsize=10, pad=22)
        _status_ax(ax, bool(vc["consistent"]), verdict, P)
        fig.tight_layout()
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"vc_reischle.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)
