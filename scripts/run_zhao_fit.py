"""Module D validation driver (update-plan Task 4): fits fsim_core.sde's
stochastic three-level QD-laser drive engine to cards/zhao.yaml's published
g2(0) values (Zhao et al., PRResearch 6, L032021 (2024) / arXiv:2309.09703
[V, read]) for normal (shot-noise-limited) vs quiet (sub-Poissonian,
high-impedance-drive) electrical pumping at I = 4 x I_th, and writes the
out/zhao/ report bundle (CSV + figure, house convention -- see
fsim_viz.figures.zhao_fit_figure).

TUNED-PARAMETER TABLE (Tier-2 fit; every entry [E] unless noted -- see
fsim_core/sde.py's QDLaserParams field comments for the literature band
behind each one):

    N_dots       = 2.0e4     (photon/carrier-number scale; sets the g2-1 magnitude)
    E_a_ES_meV   = 200.0     (RS-ES escape suppression; pushed to the deep end of a widened 40-250 meV band)
    E_a_GS_meV   = 150.0     (ES-GS escape suppression; same widened band)
    eps_gain     = 2.5e-3    (nonlinear gain compression; damps the relaxation
                              oscillation -- REQUIRED for the model to be
                              numerically usable at all; without it S(t) rings
                              with near-undamped amplitude swings of order the
                              mean itself)
    tau_p, g0, beta_sp, tau_cap, tau_relax, tau_*_spon, degeneracies: literature
    defaults (see fsim_core/sde.py), not separately re-tuned -- varying them
    inside their stated bands did not change the qualitative finding below.

HONEST RESULT (do not be misled by a PASS on the normal-pump row): the
"normal" (F_pump=1) pump-noise baseline reproduces cards/zhao.yaml's
1.0224 +- 0.003 point closely (this module's deterministic threshold physics
and F=1 shot-noise floor are well tuned). The "quiet" (F_pump=0.05-0.1) row
does NOT reproduce the published sub-unity g2(0) = 0.9823 +- 0.006 -- this
model's g2 stays measurably above 1 for every F_pump tried, including the
unphysical extreme F_pump=0 (a perfectly silent pump). fsim_core/sde.py's
module docstring ("TUNING NEAR-MISS") gives the full numerical investigation
(12 tuning experiments) and the mechanism: pump-uncorrelated (F=1) shot noise
injected at each of the three RS->ES->GS->S relay hops dilutes the pump's
regulation signal by roughly half per hop, and the direct GS->S stimulated-
emission/loss shot noise (itself F=1, comparable in magnitude to R_pump) mops
up what little survives -- under 1% of Var(S) is pump-attributable by the
time it reaches the photon field. This is reported here without fudging, per
the plan's own instruction: a documented near-miss, not a forced pass.
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card
from fsim_core.sde import QDLaserParams, find_threshold, g2_vs_pump
from fsim_viz.figures import zhao_fit_figure

T_END = 2.0e-9      # s -- ~1000 photon lifetimes (tau_p=2ps default)
DT = 1.0e-13         # s
N_RUNS_MAIN = 60      # at I=4 I_th: SE << the 0.003/0.006 published error bars
N_RUNS_CURVE = 30     # at the extra curve points (1.2, 2 I_th): visual only
F_PUMP_NORMAL = 1.0
F_PUMP_QUIET = 0.08   # [E: 0.05-0.1 residual-suppression band, per the task brief]
CURVE_I_OVER_ITH = (1.2, 2.0, 4.0)


def _write_csv(path: Path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    t_start = time.time()
    card = load_card(ROOT / "cards" / "zhao.yaml")
    pub = {r["pump"]: r for r in card.datasets["g2_vs_pump"].rows}
    print(f"cards/zhao.yaml: {card.meta['source']}")
    print(f"  published: normal g2={pub['normal']['g2']} +- {pub['normal']['err']}   "
          f"quiet g2={pub['quiet']['g2']} +- {pub['quiet']['err']}   "
          f"(both at I/I_th = {pub['normal']['I_over_Ith']})")

    params = QDLaserParams()
    I_th = find_threshold(params)
    print(f"\nfind_threshold(params) = {I_th:.4f} (raw I_over_Ith axis; 1.0 = params.I_ref by calibration)")

    I_op = 4.0 * I_th
    print(f"\nrunning at I = 4 x I_th (I_over_Ith axis value {I_op:.4f}), "
          f"n_runs={N_RUNS_MAIN}, t_end={T_END*1e9:.1f} ns, dt={DT*1e12:.2f} ps ...")

    results = {}
    for label, F_pump in (("normal", F_PUMP_NORMAL), ("quiet", F_PUMP_QUIET)):
        r = g2_vs_pump(params, I_op, n_runs=N_RUNS_MAIN, t_end=T_END, dt=DT,
                        seed=42, F_pump=F_pump)
        results[label] = r

    print(f"\n{'pump':8s} {'F_pump':8s} {'sim g2':>10s} {'sim SE':>9s} {'published':>11s} "
          f"{'pub err':>8s} {'diff':>8s} {'sigma':>7s} {'verdict':>8s}")
    rows = []
    all_pass = True
    for label, F_pump in (("normal", F_PUMP_NORMAL), ("quiet", F_PUMP_QUIET)):
        r = results[label]
        p = pub[label]
        diff = r["g2_mean"] - p["g2"]
        combined_sigma = float(np.hypot(r["g2_se"], p["err"]))
        within_1sigma = abs(diff) <= combined_sigma
        all_pass = all_pass and within_1sigma
        verdict = "PASS" if within_1sigma else "FAIL"
        print(f"{label:8s} {F_pump:8.3f} "
              f"{r['g2_mean']:10.4f} {r['g2_se']:9.4f} {p['g2']:11.4f} {p['err']:8.4f} "
              f"{diff:8.4f} {abs(diff)/combined_sigma:7.2f} {verdict:>8s}")
        rows.append(dict(pump=label, F_pump=(F_PUMP_NORMAL if label == "normal" else F_PUMP_QUIET),
                          sim_g2=r["g2_mean"], sim_se=r["g2_se"], published_g2=p["g2"],
                          published_err=p["err"], diff=diff, combined_sigma=combined_sigma,
                          sigma_ratio=abs(diff) / combined_sigma, verdict=verdict))

    print()
    if all_pass:
        print("ZHAO TIER-2: PASS -- both normal and quiet g2(0) reproduced within combined 1-sigma.")
    else:
        print("ZHAO TIER-2: FAIL -- see per-row verdicts above and the residuals table. This is an "
              "HONEST near-miss, not a fudge: fsim_core/sde.py's module docstring ('TUNING NEAR-MISS') "
              "and this script's own module docstring document the 12-experiment tuning investigation "
              "and the mechanism (each carrier-relay hop's own pump-uncorrelated shot noise dilutes "
              "the pump-statistics signal by more than an order of magnitude before it reaches the "
              "photon field). The normal-pump (Poisson, F_pump=1) row is the honest test of this "
              "model's baseline physics; the quiet-pump row is the one that falls short.")

    # ---- g2 vs I/I_th curve (both pump modes) -- figure + CSV
    print(f"\ng2 vs I/I_th curve at {CURVE_I_OVER_ITH} (n_runs={N_RUNS_CURVE}) ...")
    curve = {"normal": {"I": [], "g2": [], "se": []}, "quiet": {"I": [], "g2": [], "se": []}}
    for I_over_Ith in CURVE_I_OVER_ITH:
        for label, F_pump in (("normal", F_PUMP_NORMAL), ("quiet", F_PUMP_QUIET)):
            if I_over_Ith == 4.0:
                r = results[label]  # reuse the main-comparison run, already computed
            else:
                r = g2_vs_pump(params, I_over_Ith * I_th, n_runs=N_RUNS_CURVE, t_end=T_END, dt=DT,
                                seed=42, F_pump=F_pump)
            curve[label]["I"].append(I_over_Ith)
            curve[label]["g2"].append(r["g2_mean"])
            curve[label]["se"].append(r["g2_se"])
            print(f"  I/I_th={I_over_Ith:4.1f}  {label:7s}: g2={r['g2_mean']:.4f} +- {r['g2_se']:.4f}")

    outdir = ROOT / "out" / "zhao"
    outdir.mkdir(parents=True, exist_ok=True)
    _write_csv(outdir / "zhao_fit_comparison.csv",
               list(rows[0].keys()), (list(r.values()) for r in rows))
    _write_csv(outdir / "zhao_g2_vs_pump_curve.csv",
               ["I_over_Ith", "g2_normal", "se_normal", "g2_quiet", "se_quiet"],
               zip(curve["normal"]["I"], curve["normal"]["g2"], curve["normal"]["se"],
                   curve["quiet"]["g2"], curve["quiet"]["se"]))
    zhao_fit_figure(rows, curve, I_th, outdir)

    print(f"\nzhao bundle -> {outdir}  (zhao_fit_comparison.csv, zhao_g2_vs_pump_curve.csv, "
          f"zhao_fit_figure.[pdf/svg/png])")
    print(f"\ntotal wall time: {time.time()-t_start:.1f} s")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
