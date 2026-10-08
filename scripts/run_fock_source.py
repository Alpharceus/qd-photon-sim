"""Sweep the cryogenic Fock-state source card and write the boson-sampling
figure of merit.

Run:  python scripts/run_fock_source.py          (writes ONLY to out/fock_source/)
      python scripts/run_fock_source.py --quick --out-dir <dir>
          reduced sweep into another directory; used by verify_resonant_source (m4)
          to test determinism without touching out/fock_source/.

Inputs: cards/fock-source-4K-design.yaml (fsim_core.resonant_source.load_fock_card).
Every [A] range of the card is swept (lo / mid / hi or lo / hi, see AXES) and
never averaged; the four efficiencies (and beta0) enter eta monotonically and
are swept as the two joint corners "pess" (all lo) and "opt" (all hi), which
bracket eta exactly. The ZPL-filtered output is swept over off / on at the
two phonon alpha values; the dephasing mode (enhanced / dot, see
resonant_source.source_rates) is swept as a categorical axis.

Phonon treatment (physics review 2026-10-08, finding 1): ZPL output "off" keeps
the sideband photons in eta but multiplies the off-diagonal HOM by Z(T)^2
(Iles-Smith et al., Nat. Photonics 11, 521 (2017)); "on" multiplies eta by Z
and leaves the HOM unpenalised. Both are applied through
rs.hom_phonon_factor / rs.end_to_end_eta.

Outputs in out/fock_source/ (no timestamps: a second run is byte-identical)
  sweep.csv           one row per (physics point, efficiency scenario, ZPL)
  summary.json        verdict numbers and best corners
  results.md          short report with a VERDICT: line
  bs_<best|central>_N<N>.json   per-run export for the Julia MZI-mesh code

JSON schema of bs_*.json  (consumed by the Julia MZI-mesh simulator)
  schema        "fock_source/bs/1"
  N             int, number of photons entering the network
  f_rep_MHz     float, repetition rate
  g2            list[N] of float, per-photon g2(0) (identical: one source)
  eta           float, end-to-end per-photon efficiency
                eta_source * eta_demux * eta_circuit * eta_det
  eta_components dict(source, demux, circuit, det, n_ph, beta, eta_out, Z)
  hom_matrix    N x N nested list, pairwise HOM indistinguishability I from the
                spectral-diffusion model (F3 hom_vs_separation at
                dt_ij = |i-j|/f_rep); symmetric, unit diagonal by definition
  hom_matrix_pulse_corrected  same times the single-pulse factor
                (re-excitation / cascade cap); the Julia code should use this
                one if it wants the conservative value
  pulse_factor  float
  R_N_Hz        float, (f_rep/N) eta^N
  provenance    dict(card entry -> tag A|E) for every input
  config        dict(card entry -> value) of the evaluated corner
  note          str

Decisions (see also the module docstring of fsim_core.resonant_source):
  * VERDICT thresholds are [A]: per-photon g2(0) <= G2_MAX and the CONSERVATIVE
    mean pairwise pulse-corrected I - g2(0) >= I_MIN (the work-order's I >= 0.9
    anchor; the two-photon component lowers the measured visibility by O(g2),
    review finding 4). The uncorrected I columns are kept in the CSV.
  * Best corner per N = the qualifying row with the largest R_N; if no row
    qualifies, the row with the largest conservative I (then the lowest g2) is
    reported and flagged.
  * A fixed-seed Gillespie Monte Carlo cross-checks g2 of the minimum-g2 resonant
    corner (rectangular-pulse equivalent, gamma* = 0).
"""
import argparse
import csv
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import resonant_source as rs  # noqa: E402

CARD = ROOT / "cards" / "fock-source-4K-design.yaml"
OUT = ROOT / "out" / "fock_source"
NS = (4, 8, 20)
G2_MAX = 0.05          # [A] per-photon g2(0) requirement
I_MIN = 0.9            # [A] work-order anchor: I >= 0.9
MC_SEED = 20261008
MC_TRAJ = 1_000_000

# swept physics axes: card entry -> number of points over its range
AXES = {"gamma0_ns": 2, "F_P": 3, "kappa_ueV": 2, "pulse_fwhm_ps": 2, "hom_width_multiple": 3,
        "sigma_sd_ueV": 3, "tau_sd_ns": 2}
TPE_AXES = {"delta_xx_meV": 2}
EFF_KEYS = ("beta0", "eta_out", "eta_demux", "eta_circuit", "eta_det")
COMBOS = (("X+", "resonant_pi"), ("X", "resonant_pi"), ("X", "tpe"))


def f6(x):
    return repr(float(x))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="reduced sweep (verification)")
    ap.add_argument("--out-dir", default=None, help="output directory (default: out/fock_source)")
    args = ap.parse_args(argv)
    out_dir = Path(args.out_dir) if args.out_dir else OUT
    axes_pts = dict(AXES)
    if args.quick:
        axes_pts.update(gamma0_ns=1, F_P=2, kappa_ueV=1, pulse_fwhm_ps=1, hom_width_multiple=2,
                        sigma_sd_ueV=2, tau_sd_ns=2)
    mc_traj = 100_000 if args.quick else MC_TRAJ
    card = rs.load_fock_card(CARD)
    out_dir.mkdir(parents=True, exist_ok=True)
    shapes = list(card["pulse_shape"].points())
    if args.quick:
        shapes = shapes[:1]
        tpe_axes = {"delta_xx_meV": 1}
    else:
        tpe_axes = dict(TPE_AXES)
    alphas = card["phonon_alpha_ps2"].points(2)
    zpl_choices = ([(f"off_alpha{a:g}", False, a) for a in alphas]
                   + [(f"on_alpha{a:g}", True, a) for a in alphas])
    dm_choices = list(card["dephasing_mode"].points())
    eff = {"pess": {k: card[k].points(2)[0] for k in EFF_KEYS},
           "opt": {k: card[k].points(2)[1] for k in EFF_KEYS}}
    off_lo, off_hi = card["remote_offset_ueV"].points(2)
    rows, cfgs = [], []
    fixed = {n: card[n].central() for n in ("T_K", "cav_detune_ueV", "gamma_xx_ratio", "f_rep_MHz")}
    f_rep = fixed["f_rep_MHz"]
    for (trans, exc), shape, dmode in product(COMBOS, shapes, dm_choices):
        axes = dict(axes_pts)
        if exc == "tpe":
            axes.update(tpe_axes)
        names = list(axes)
        grids = [card[n].points(axes[n]) for n in names]
        for vals in product(*grids):
            ov = dict(zip(names, vals))
            ov.update(transition=trans, excitation=exc, pulse_shape=shape, dephasing_mode=dmode)
            # zpl_output=True here: the matrices come back WITHOUT the Z^2 penalty; the
            # per-row factor rs.hom_phonon_factor (Z^2 off, 1 on) is applied below
            src0 = rs.source_from_card(card, zpl_output=True, **ov)
            r = rs.bs_figure_of_merit(src0, max(NS))             # one call per physics point
            Hfull = np.array(r["hom_matrix"])
            Hcfull = np.array(r["hom_matrix_pulse_corrected"])
            mI0, mIc0 = {}, {}
            for N in NS:
                off = ~np.eye(N, dtype=bool)
                mI0[N] = float(Hfull[:N, :N][off].mean())
                mIc0[N] = float(Hcfull[:N, :N][off].mean())
            rates = rs.source_rates(src0)
            sig = rs._rate_from_ueV(src0.sigma_sd_ueV)
            I_rem0 = [rs.hom_remote(rates["Gamma"], rates["gamma_star"], sig, rs._rate_from_ueV(o))
                      for o in (off_lo, off_hi)]
            g2, n_ph = r["g2"][0], r["eta_components"]["n_ph"]
            for ename, ev in eff.items():
                for zname, zon, alpha in zpl_choices:
                    ov2 = dict(ov, **ev)
                    ov2["zpl_output"] = zon
                    ov2["phonon_alpha_ps2"] = alpha
                    src = rs.source_from_card(card, **ov2)
                    eta, _ = rs.end_to_end_eta(src, n_ph)
                    zf = rs.hom_phonon_factor(src)
                    mI = {N: mI0[N] * zf for N in NS}
                    mIc = {N: mIc0[N] * zf for N in NS}
                    I_rem = [v * zf for v in I_rem0]
                    row = {"transition": trans, "excitation": exc, "pulse_shape": shape,
                           "dephasing_mode": dmode}
                    for n in AXES:
                        row[n] = ov[n]
                    row["delta_xx_meV"] = ov.get("delta_xx_meV", "")
                    row.update(eff_scenario=ename, zpl=zname, Gamma_ns=rates["Gamma"],
                               gamma_star_ns=rates["gamma_star"], g2=g2, n_ph=n_ph,
                               pulse_factor=r["pulse_factor"], hom_phonon_factor=zf, eta=eta)
                    for N in NS:
                        row[f"I_mean_N{N}"] = mI[N]
                        row[f"I_mean_corr_N{N}"] = mIc[N]
                        row[f"I_cons_N{N}"] = mIc[N] - g2
                        row[f"R_N{N}_Hz"] = rs.coincidence_rate_Hz(f_rep, N, eta)
                    row["I_remote_offset_lo"], row["I_remote_offset_hi"] = I_rem
                    rows.append(row)
                    cfgs.append(ov2)

    cols = list(rows[0])
    with open(out_dir / "sweep.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for row in rows:
            w.writerow([f6(row[c]) if isinstance(row[c], (float, np.floating)) else row[c] for c in cols])

    summary = {"rows": len(rows), "G2_MAX": G2_MAX, "I_MIN": I_MIN, "thresholds_tag": "A"}
    g2_min_i = int(np.argmin([r["g2"] for r in rows]))
    summary["g2_min"] = rows[g2_min_i]["g2"]
    summary["g2_min_config"] = cfgs[g2_min_i]
    best = {}
    for N in NS:
        Imax_i = int(np.argmax([r[f"I_mean_N{N}"] for r in rows]))
        Icmax_i = int(np.argmax([r[f"I_mean_corr_N{N}"] for r in rows]))
        Icons_i = int(np.argmax([r[f"I_cons_N{N}"] for r in rows]))
        ok = [i for i, r in enumerate(rows) if r["g2"] <= G2_MAX and r[f"I_cons_N{N}"] >= I_MIN]
        if ok:
            bi = max(ok, key=lambda i: rows[i][f"R_N{N}_Hz"])
            flag = "qualifying"
        else:
            bi = max(range(len(rows)), key=lambda i: (rows[i][f"I_cons_N{N}"], -rows[i]["g2"]))
            flag = "no-qualifying-corner"
        best[N] = bi
        summary[f"N{N}"] = {"I_mean_max": rows[Imax_i][f"I_mean_N{N}"],
                            "I_mean_corr_max": rows[Icmax_i][f"I_mean_corr_N{N}"],
                            "I_cons_max": rows[Icons_i][f"I_cons_N{N}"],
                            "qualifying_rows": len(ok), "best_flag": flag,
                            "best_row": rows[bi], "best_config": cfgs[bi],
                            "R_N_Hz_best": rows[bi][f"R_N{N}_Hz"]}

    # per-run JSON exports (best corner and the card's central point)
    central = rs.source_from_card(card)
    for N in NS:
        for tag, src in (("best", rs.source_from_card(card, **cfgs[best[N]])), ("central", central)):
            out = rs.bs_figure_of_merit(src, N)
            out["schema"] = "fock_source/bs/1"
            out["config"] = {k: getattr(src, k) for k in rs._SOURCE_FIELDS}
            out["note"] = ("per-photon g2 and eta of a single time-demultiplexed source; "
                           "hom_matrix = F3 spectral-diffusion model, hom_matrix_pulse_corrected adds the "
                           "single-pulse (re-excitation / cascade) factor")
            with open(out_dir / f"bs_{tag}_N{N}.json", "w", encoding="utf-8") as fh:
                json.dump(out, fh, indent=1, sort_keys=True)
                fh.write("\n")

    # fixed-seed second method: Gillespie MC of the minimum-g2 resonant corner (rect equivalent, gamma* = 0)
    res_rows = [i for i, r in enumerate(rows) if r["excitation"] == "resonant_pi"]
    bi = min(res_rows, key=lambda i: rows[i]["g2"])
    Gb, tb = rows[bi]["Gamma_ns"], cfgs[bi]["pulse_fwhm_ps"] * 1e-3
    tl = rs.TwoLevel(Gb)
    det = rs.pulsed_photon_statistics(tl, "rect", np.pi, tb, f_rep, periodic=False)
    mc = rs.pulsed_counting_mc(tl, np.pi, tb, f_rep, n_traj=mc_traj, seed=MC_SEED)
    z = (mc["g2"] - det["g2"]) / mc["g2_se"]
    summary["mc_check"] = {"Gamma_ns": Gb, "tau_p_ns": tb, "g2_det": det["g2"], "g2_mc": mc["g2"],
                           "g2_mc_se": mc["g2_se"], "z": z, "seed": MC_SEED, "n_traj": mc_traj}
    with open(out_dir / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=float)
        fh.write("\n")

    # report
    npass = [N for N in NS if summary[f"N{N}"]["best_flag"] == "qualifying"]
    verdict = "PASS" if len(npass) == len(NS) else ("PARTIAL" if npass else "FAIL")
    vl = (f"VERDICT: {verdict} g2_min={summary['g2_min']:.4g} "
          + " ".join(f"I_mean_corr_max_N{N}={summary[f'N{N}']['I_mean_corr_max']:.4f}" for N in NS) + " "
          + " ".join(f"I_cons_max_N{N}={summary[f'N{N}']['I_cons_max']:.4f}" for N in NS) + " "
          + " ".join(f"R_N{N}_best_Hz={summary[f'N{N}']['R_N_Hz_best']:.3e}" for N in NS)
          + f" qualifying_N={','.join(map(str, npass)) or 'none'} "
          + " ".join(f"qualifying_rows_N{N}={summary[f'N{N}']['qualifying_rows']}" for N in NS)
          + f" rows={len(rows)}")
    lines = ["# Cryogenic Fock-state source sweep", "",
             f"Card: `cards/fock-source-4K-design.yaml` (all entries [A]/[E]); {len(rows)} rows in `sweep.csv`.",
             f"Thresholds [A]: g2(0) <= {G2_MAX}, conservative mean pairwise corrected I - g2(0) >= {I_MIN}.", "",
             "```", vl, "```", "",
             f"min g2(0) = {summary['g2_min']:.4g} at {json.dumps(summary['g2_min_config'], sort_keys=True, default=float)}", ""]
    for N in NS:
        s = summary[f"N{N}"]
        lines += [f"## N = {N}",
                  f"- max mean pairwise I (spectral diffusion only) = {s['I_mean_max']:.4f}; "
                  f"with the single-pulse factor = {s['I_mean_corr_max']:.4f}; "
                  f"conservative I - g2 = {s['I_cons_max']:.4f}",
                  f"- qualifying rows: {s['qualifying_rows']} / {len(rows)}; best corner is {s['best_flag']}",
                  f"- R_N at the best corner = {s['R_N_Hz_best']:.4e} Hz "
                  f"(eta = {s['best_row']['eta']:.4f}, g2 = {s['best_row']['g2']:.4g}, "
                  f"I_corr = {s['best_row'][f'I_mean_corr_N{N}']:.4f}, "
                  f"I_cons = {s['best_row'][f'I_cons_N{N}']:.4f})",
                  f"- best config: `{json.dumps(s['best_config'], sort_keys=True, default=float)}`", ""]
    m = summary["mc_check"]
    lines += ["## Second method", f"Gillespie MC (seed {m['seed']}, {m['n_traj']} trajectories) of the minimum-g2 resonant "
              f"corner, rectangular-pulse equivalent: g2 = {m['g2_mc']:.5f} +- {m['g2_mc_se']:.5f} vs deterministic "
              f"{m['g2_det']:.5f} (z = {m['z']:+.2f}).", "",
              "## Reading the numbers",
              "- PASS / PARTIAL / FAIL refers to the swept ASSUMED ranges: PASS means at least one corner of the ranges "
              "meets the [A] thresholds, not that a device does; `qualifying_rows_N*` in the VERDICT line counts the corners.",
              "- Ranges are swept, never averaged; the efficiencies are bracketed by the joint corners `pess` and `opt`.",
              "- `R_N = (f_rep/N) eta^N` [DR] is for time demultiplexing of ONE source; it falls super-exponentially with N.",
              "- Unfiltered output (`zpl=off_*`) multiplies the off-diagonal HOM by Z(T)^2 (sideband photons are kept in eta but are "
              "distinguishable); ZPL-filtered output (`zpl=on_*`) multiplies eta by Z instead.",
              "- `dephasing_mode`: enhanced = gamma* scales with the Purcell-enhanced rate (I = 1/m independent of F_P by construction); "
              "dot = gamma* is a bare-dot property. Swept, never averaged.",
              "- The verdict uses the conservative I - g2(0): the two-photon component lowers the measured visibility by O(g2).",
              "- The pulse factor couples re-excitation (resonant) or the cascade cap (TPE) to I by an [A] factorisation.",
              "- The V-a Gamma_0 = 119 ueV is not used; the homogeneous width is a multiple of the transform limit.", ""]
    (out_dir / "results.md").write_text("\n".join(lines), encoding="utf-8")
    print(vl)
    return 0


if __name__ == "__main__":
    sys.exit(main())
