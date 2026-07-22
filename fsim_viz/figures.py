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

    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.fill_between(Ts, lo, hi, alpha=0.25, color="#2166ac", lw=0,
                    label="F1 reachable envelope (swept [E] inputs)")
    ax.errorbar(Ts, [r["g2"] for r in vb["data"]], yerr=[r["err"] for r in vb["data"]],
                fmt="o", color="#333333", capsize=4, ms=8, label="measured (Fig. 5)", zorder=5)
    ax.axhline(0.5, color="#d7191c", ls=":", lw=1.5)
    ax.annotate("Theorem-0 saturated bound (ε→1, μ→∞)", (90, 0.51), fontsize=8,
                color="#d7191c")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("$g^{(2)}(0)$")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ok = bool(vb["covered_delta"])
    ax.set_title("V-b (Laferrière): ε→1 limit under published windows — "
                 + ("CONSISTENT" if ok else "NOT COVERED"),
                 fontsize=10, color="#1a9641" if ok else "#d7191c")
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

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axes[0]
    ax.plot(cv["Ts"], cv["det_meV"], color="#2166ac", lw=2)
    ax.axhline(0, color="#888888", ls=":", lw=1)
    ax.axvline(cv["T_target"], color="#d7191c", ls="--", lw=1.2)
    ax.annotate(f"$T_{{target}}$ = {cv['T_target']:.0f} K",
                (cv["T_target"], cv["det_meV"].max() * 0.7),
                color="#d7191c", fontsize=9, rotation=90, ha="right")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("X–mode detuning (meV)")
    ax.set_title("mode-tracking rule (F6 iv): mode red of\ncryogenic line, zero at $T_{target}$",
                 fontsize=9)

    ax = axes[1]
    ax.fill_between(cv["kappas"], cv["eps_band"].min(axis=1), cv["eps_band"].max(axis=1),
                    alpha=0.3, color="#5e3c99", lw=0)
    ax.plot(cv["kappas"], cv["eps_band"].min(axis=1), color="#5e3c99", lw=1.5,
            label=rf"$\Delta$ = {cv['deltas'][-1]:.1f} meV")
    ax.plot(cv["kappas"], cv["eps_band"].max(axis=1), color="#5e3c99", lw=1.5, ls="--",
            label=rf"$\Delta$ = {cv['deltas'][0]:.1f} meV")
    ax.set_xlabel(r"cavity $\kappa$ (meV)")
    ax.set_ylabel(r"$\varepsilon$ at $T_{target}$ (cavity filter)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(rf"cavity-only $\varepsilon$; $\Gamma(T_t)$ = {cv['gam_t']:.1f} meV [A proxy]",
                 fontsize=9)

    ax = axes[2]
    ax.plot(cv["kappas"], cv["Feff"], color="#e66101", lw=2)
    ax.set_xlabel(r"cavity $\kappa$ (meV)")
    ax.set_ylabel(r"$F_{eff}/F_P = \kappa/(\kappa+\Gamma)$")
    ax.set_title("spectral-overlap Purcell penalty (F6 iii)", fontsize=9)

    fig.suptitle("Phase 2 — Module B design rules (tag chain [A]: requirement envelopes; "
                 "κ, F_P, G await MEEP/COMSOL)", fontsize=10, color="#d7191c", y=1.03)
    fig.tight_layout()
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

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

    ax = axes[0, 0]
    cf = ax.contourf(m["deltas"], m["rhos"], m["Tc"], levels=np.arange(50, 425, 25),
                     cmap="viridis")
    cs = ax.contour(m["deltas"], m["rhos"], m["Tc"], levels=[77, 120, 200, 300],
                    colors="white", linewidths=1.2)
    ax.clabel(cs, fmt="%.0f K", fontsize=7)
    ax.plot(m["deltas"], m["rho_req"][g_lo], color="#d7191c", lw=2)
    ax.plot(m["deltas"], m["rho_req"][g_hi], color="#d7191c", lw=2, ls="--")
    ax.fill_between(m["deltas"], m["rho_req"][g_lo], m["rho_req"][g_hi],
                    color="#d7191c", alpha=0.25, lw=0,
                    label=r"$\rho$ required for 300 K [$\Gamma$(300) 6–7 meV]")
    ax.axvline(5.0, color="#e3a21a", ls=":", lw=1.5)
    ax.annotate(">50% of (211)B dots", (5.05, 0.615), rotation=90, fontsize=7,
                color="#e3a21a")
    for p in m["points"]:
        ax.plot(p["delta"], p["rho"], "o", color="#ffffff", mec="#333333", ms=7)
        ax.annotate(p["name"], (p["delta"] + 0.15, p["rho"] - 0.012), fontsize=6.5)
    plt.colorbar(cf, ax=ax, label="$T_c$ (K)")
    ax.set_xlabel(r"$\Delta_{XX}$ (meV)")
    ax.set_ylabel(r"signal purity $\rho$")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.set_title("(ii) master-ceiling map: $T_c(\\Delta,\\rho)$, narrow-filter bound",
                 fontsize=9)

    ax = axes[0, 1]
    colors = {77.0: "#2166ac", 120.0: "#e66101"}
    for T_hs, e in s["envelopes"].items():
        band = e["band"]
        ax.fill_between(s["deltas"], band[:, 0], band[:, 1], alpha=0.3,
                        color=colors[T_hs], lw=0)
        ax.plot(s["deltas"], band[:, 1], color=colors[T_hs], lw=1.8,
                label=f"$T_{{hs}}$ = {T_hs:.0f} K (worst case)")
    ax.axhline(0.5, color="#888888", ls=":", lw=1)
    ax.axhline(0.1, color="#888888", ls=":", lw=1)
    ax.set_xlabel(r"$\Delta_{XX}$ (meV)  [A: unmeasured on InP/GaAsP]")
    ax.set_ylabel("$g^{(2)}(0)$ envelope")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("(i) staged device envelope, electrical, w=$\\Gamma$ convention",
                 fontsize=9)

    ax = axes[1, 0]
    cf = ax.contourf(a["diams"], a["dens"], a["g2pen"],
                     levels=[0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0], cmap="magma_r")
    cs = ax.contour(a["diams"], a["dens"], a["Nw"], levels=[1.0], colors="#2166ac",
                    linewidths=1.5)
    ax.clabel(cs, fmt="$N_w$=%.0f", fontsize=7)
    ax.set_yscale("log")
    plt.colorbar(cf, ax=ax, label="F5 aperture $g^{(2)}$ penalty")
    ax.set_xlabel("aperture diameter (µm)")
    ax.set_ylabel("QD density (cm$^{-2}$)")
    ax.set_title(f"(iii) F5 aperture/density rules (w = {a['w']:.1f} meV)", fontsize=9)

    ax = axes[1, 1]
    items = r["ranked"][::-1]
    ax.barh([r["labels"][k] for k, _ in items], [v for _, v in items],
            color="#2166ac")
    ax.set_xlabel(f"$T_c$ envelope narrowing (K) out of {r['width_full']:.0f} K total")
    ax.set_title("(iv) measurement-priority ranking (in-house queue)", fontsize=9)
    ax.tick_params(axis="y", labelsize=7)

    fig.suptitle("Phase 3 — requirement envelopes and design rules "
                 "(tag chain [A]: every panel inherits unmeasured inputs; "
                 "see the packet note for the chains)",
                 fontsize=10, color="#d7191c", y=1.0)
    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"phase3_packet.{ext}", bbox_inches="tight", dpi=200)
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

    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.fill_between(ws, lo, hi, alpha=0.3, color="#2166ac", lw=0,
                    label=r"digitized $\rho(w)$ (zero-level conventions)")
    ax.plot(ws, lo, color="#2166ac", lw=1.2)
    ax.plot(ws, hi, color="#2166ac", lw=1.2)
    ax.axhspan(vc["rho_req"] - vc["rho_req_err"], vc["rho_req"] + vc["rho_req_err"],
               color="#d7191c", alpha=0.35, lw=0)
    ax.axhline(vc["rho_req"], color="#d7191c", lw=1.5,
               label=rf"$\rho$ required by $g^{{(2)}}$ = {vc['g2']}$\pm${vc['g2_err']}")
    ax.set_xlabel("detection window $w$ (meV)")
    ax.set_ylabel(r"signal fraction $\rho$")
    ax.set_ylim(0.7, 1.0)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    verdict = "CONSISTENT" if vc["consistent"] else "INCONSISTENT"
    ax.set_title(f"V-c (Reischle, 100 MHz, ~40 K): trion, $\\varepsilon\\approx 0$ — "
                 f"{verdict}", fontsize=10,
                 color="#1a9641" if vc["consistent"] else "#d7191c")
    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"vc_reischle.{ext}", bbox_inches="tight", dpi=200)
    plt.close(fig)
