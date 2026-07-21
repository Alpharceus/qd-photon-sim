"""fsim-viz: figure factory. Consumes fsim-core results objects only (three-layer
rule) and writes, for every figure, the underlying CSV alongside PDF/SVG/PNG --
the group re-renders in Origin, so data export is first-class.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fsim_core.card import Tag
from fsim_core.fitting import FitResult, V_A_TOL
from fsim_core.integrator import g2_of_T, solve_Tc

TAG_COLORS = {Tag.V: "#1a9641", Tag.DR: "#e3a21a", Tag.E: "#e3a21a", Tag.A: "#d7191c"}


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
    tagc = TAG_COLORS[fit.tag]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

    ax = axes[0]
    val = ~upper
    ax.errorbar(Ts_d[val], g2_d[val], yerr=err_d[val], fmt="o", color="#333333",
                capsize=3, label="data (Fig. 5)", zorder=5)
    if upper.any():
        ax.errorbar(Ts_d[upper], g2_d[upper], yerr=err_d[upper], fmt="v",
                    color="#333333", capsize=3, uplims=True, zorder=5,
                    label="upper bound")
    ax.plot(Ts, g2_m, color="#2166ac", lw=2, label="F-series model")
    ax.axhline(0.5, color="#888888", ls=":", lw=1)
    if np.isfinite(Tc):
        ax.axvline(Tc, color="#d7191c", ls="--", lw=1.2)
        ax.annotate(f"$T_c$ = {Tc:.0f} K", (Tc, 0.52), color="#d7191c", fontsize=9,
                    ha="right", rotation=90)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("$g^{(2)}(0)$")
    ax.set_ylim(0, max(0.7, g2_m.max() * 1.05))
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("V-a fit: $g^{(2)}(T)$, published windows", fontsize=10)

    ax = axes[1]
    ax.plot(Ts, eps_m, color="#5e3c99", lw=2, label=r"$\varepsilon(T)=t_{XX}/t_X$")
    ax.plot(Ts, rho_m**2, color="#e66101", lw=2, label=r"$\rho(T)^2$")
    ax.plot(Ts, rho_m**2 * (1 - eps_m), color="#333333", lw=1.4, ls="--",
            label=r"$\rho^2(1-\varepsilon)$")
    ax.axhline(0.5, color="#888888", ls=":", lw=1)
    ax.set_xlabel("T (K)")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("decomposition (master ceiling at 1/2)", fontsize=10)

    ax = axes[2]
    ax.plot(Ts, gam_m, color="#2166ac", lw=2)
    gA, gAv, gAt, T_anchor = fit.gamma_anchor
    ax.errorbar([T_anchor], [gAv], yerr=[gAt], fmt="s", color="#d7191c",
                capsize=4, label="published anchor")
    ax.plot([T_anchor], [gA], "x", color="#2166ac", ms=9, mew=2, label="model @ anchor")
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$\Gamma$ (meV)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(r"$\Gamma(T)$ vs anchor", fontsize=10)

    status = "PASS" if fit.passed else "FAIL"
    fig.suptitle(
        f"{Path(card_path).stem}: V-a fit [{status} at ±{V_A_TOL}]   "
        f"tag chain {fit.tag.label}" + ("  — PLACEHOLDER DATA" if fit.notes else ""),
        fontsize=11, color=tagc, y=1.02,
    )
    fig.tight_layout()
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

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

    ax = axes[0]
    for e, f in drive["curves"].items():
        ax.plot(mus, f, lw=2, label=rf"$\varepsilon$ = {e}")
    ax.axvline(drive["mu_op"], color="#888888", ls="--", lw=1)
    ax.annotate(rf"$\mu_{{op}}$ = {drive['mu_op']}", (drive["mu_op"], 1.85),
                fontsize=8, ha="left", color="#555555")
    ax.axhline(2.0, color="#888888", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel(r"mean loading $\mu$")
    ax.set_ylabel(r"$g^{(2)}(\mu)\,/\,\varepsilon$")
    ax.set_ylim(0.95, 2.1)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("F1b: finite-$\\mu$ drive penalty (cap-2)", fontsize=10)

    for ax, tpl in zip(axes[1:], thermal["templates"]):
        per_I = thermal["templates"][tpl]
        for I, color in zip(thermal["currents_uA"], ("#2166ac", "#e66101", "#d7191c")):
            lo, hi = per_I[I]
            lo = np.where(np.isinf(lo), np.nan, lo)  # runaway region: gap, not a line
            hi = np.where(np.isinf(hi), np.nan, hi)
            ax.fill_between(d_um, lo, hi, alpha=0.25, color=color, lw=0)
            ax.plot(d_um, hi, color=color, lw=1.5,
                    label=f"{I} µA" + (" (runaway ←)" if np.isnan(hi).any() else ""))
        ax.set_xlabel("mesa diameter (µm)")
        ax.set_ylabel(r"$\Delta T_J$ (K)")
        ax.set_yscale("log")
        ax.legend(frameon=False, fontsize=8, title="CW drive", title_fontsize=8)
        ax.set_title(f"{tpl}: $T_J-T_{{hs}}$ @ {thermal['T_hs']:.0f} K "
                     f"(band: epi-k range)", fontsize=10)

    fig.suptitle("Phase 1 — Module D drive penalty and Module A junction-heating envelopes "
                 "(tag chain [A]: requirement envelopes, not predictions)",
                 fontsize=10, color="#d7191c", y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"phase1_drive_thermal.{ext}", bbox_inches="tight", dpi=200)
    plt.close(fig)
    return {"outdir": str(outdir)}
