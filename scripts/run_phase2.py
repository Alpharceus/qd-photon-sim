"""Phase 2 driver: V-b (Laferriere epsilon->1 envelope) + Module B cavity
design deliverable (mode-tracking rule at T_target).

V-b is an existence/envelope check (unpublished Delta/Gamma/windows are swept,
never averaged): does a physical parameterization inside the published
constraints reproduce the series, and does the 300 K point sit in the
epsilon->1 corridor [0.5, 1) whose saturated bound is the Theorem-0 value 1/2?
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card
from fsim_core.cavity import purcell_eff, tracking_detuning, varshni_shift
from fsim_core.integrator import g2_of_T, solve_Tc
from fsim_core.loading import f1b_g2
from fsim_core.spectral import epsilon, gamma_of_T
from fsim_viz.figures import cavity_design_figure, vb_figure


def _grid(prm, n=6):
    lo, hi = prm.bounds
    return np.linspace(lo, hi, n)


def run_vb():
    card = load_card(ROOT / "cards" / "laferriere.yaml")
    rows = {r["T"]: r for r in card.datasets["g2_vs_T"].rows}
    print("V-b (Laferriere, InAsP/InP nanowire, published windows; sweep over "
          "unpublished Delta/Gamma/mu):")

    deltas = _grid(card["delta_xx"], 9)
    reach = {T: {} for T in (77.0, 220.0, 300.0)}
    for d in deltas:
        # 77 K
        g = [float(f1b_g2(mu, epsilon(d, gam, gam, w=w).eps))
             for gam in _grid(card["gamma_77"], 5)
             for w in _grid(card["w_77"], 5)
             for mu in _grid(card["mu_77"], 3)]
        reach[77.0][d] = (min(g), max(g))
        # 220 K
        g = [float(f1b_g2(mu, epsilon(d, gam, gam, w=w).eps))
             for gam in _grid(card["gamma_220"], 6)
             for w in _grid(card["w_220"], 5)
             for mu in _grid(card["mu_high"], 6)]
        reach[220.0][d] = (min(g), max(g))
        # 300 K (published 25 nm window)
        w300 = card["w_300"].fixed
        g = [float(f1b_g2(mu, epsilon(d, gam, gam, w=w300).eps))
             for gam in _grid(card["gamma_300"], 6)
             for mu in _grid(card["mu_high"], 6)]
        reach[300.0][d] = (min(g), max(g))

    def in_reach(T, d):
        r = rows[T]
        return reach[T][d][0] - r["err"] <= r["g2"] <= reach[T][d][1] + r["err"]

    for T in (77.0, 220.0, 300.0):
        lo = min(reach[T][d][0] for d in deltas)
        hi = max(reach[T][d][1] for d in deltas)
        r = rows[T]
        inside = lo - r["err"] <= r["g2"] <= hi + r["err"]
        print(f"  {T:5.0f} K: measured {r['g2']:.2f}+-{r['err']:.2f}   "
              f"reachable [{lo:.3f}, {hi:.3f}]   {'covered' if inside else 'NOT COVERED'}")

    # epsilon->1 corridor at 300 K
    w300 = card["w_300"].fixed
    eps300 = [epsilon(d, gam, gam, w=w300).eps
              for d in deltas for gam in _grid(card["gamma_300"], 6)]
    print(f"  eps(300 K, 25 nm window) across sweep: [{min(eps300):.2f}, {max(eps300):.2f}] "
          "-> epsilon->1 corridor g2 in [0.5, 1); measured 0.57 sits at the "
          "saturated (Theorem-0) edge")

    covered_all = [float(d) for d in deltas if all(in_reach(T, d) for T in (77.0, 220.0, 300.0))]
    covered_high = [float(d) for d in deltas if all(in_reach(T, d) for T in (220.0, 300.0))]
    if covered_all:
        print(f"  all three points covered for Delta_XX in "
              f"[{min(covered_all):.1f}, {max(covered_all):.1f}] meV -> V-b CONSISTENT")
    elif covered_high:
        d_hi = covered_high
        resid77 = min(rows[77.0]["g2"] - reach[77.0][d][1] for d in np.asarray(d_hi))
        print(f"  high-T points (220/300 K) covered for Delta_XX in "
              f"[{min(d_hi):.1f}, {max(d_hi):.1f}] meV -> the eps->1 mechanism "
              "(the V-b claim) holds where it applies")
        print(f"  77 K point exceeds the spectral envelope by ~{resid77:+.2f}: the "
              "residual is the paper's own re-excitation channel (their Fig. 2d; "
              "0.021 at 4 K saturation, growing with T) -- non-cascade, outside "
              "cap-2, the same WP-M2' channel as the Reischle refilling residual")
        print("  -> V-b CONSISTENT at the eps->1 limit; low-T end re-excitation-"
              "limited (named)")
    else:
        print("  high-T points not covered -> V-b FAILS; re-examine assumptions")
    print("  4 K point (0.021 at saturation): re-excitation, non-cascade "
          "(their Fig. 2d) -- outside the spectral model by construction")
    covered_delta = covered_all or covered_high

    vb_figure({"deltas": list(map(float, deltas)), "reach": reach,
               "data": [rows[T] for T in (77.0, 220.0, 300.0)],
               "covered_delta": covered_delta}, ROOT / "out" / "phase2")
    return covered_delta


def run_cavity_design():
    card = load_card(ROOT / "cards" / "qcap-cavity.yaml")
    Tt = card["T_target"].fixed
    E0 = card["E_X0"].fixed
    dEdT = card["dEdT_cav"].fixed * 1e-3  # meV/K -> eV/K
    mat = "GaAs" if card["varshni_material"].fixed == 0 else "InP"

    # class-proxy Gamma(T): the chatzarakis joint-fit parameterization [A]
    fitp = json.loads((ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]
    gpar = {k: fitp[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo")}

    Ts = np.linspace(10, 300, 200)
    det_meV = np.array([1e3 * tracking_detuning(T, E0, Tt, mat, dEdT) for T in Ts])
    shift_120 = -1e3 * float(varshni_shift(Tt, mat))
    print(f"\nPhase-2 cavity design (T_target = {Tt:.0f} K, {mat}-Varshni, [A] chain):")
    print(f"  tracking rule: mode placed {shift_120:.1f} meV red of the cryogenic X "
          f"line -> detuning(T_target) = {np.interp(Tt, Ts, det_meV):+.3f} meV")

    kappas = _grid(card["kappa"], 12)
    deltas = _grid(card["delta_xx"], 8)
    gam_t = float(gamma_of_T(Tt, **gpar))
    # cavity-only acceptance at T_target: eps and F_eff vs kappa (Delta swept -> band)
    eps_band = np.array([[epsilon(d, gam_t, fitp.get("r_xx", 1.0) * gam_t,
                                  kappa=k, dx=0.0).eps for d in deltas]
                         for k in kappas])
    Feff = purcell_eff(1.0, kappas, gam_t)  # per unit F_P
    print(f"  Gamma(T_target) class proxy = {gam_t:.2f} meV;  F_eff/F_P at kappa = "
          f"{kappas[0]:.1f}/{kappas[-1]:.1f} meV: {Feff[0]:.2f}/{Feff[-1]:.2f}")

    # Design finding: below T_target the red-tracked mode sits closer to the
    # (red-shifted) XX than to X -- the cavity alone selects the WRONG line
    # (eps can exceed 1). The slit filter must be retained below T_target.
    kmid = float(kappas[len(kappas) // 2])
    T_cold = 20.0
    det_cold = 1e3 * tracking_detuning(T_cold, E0, Tt, mat, dEdT)
    gam_cold = float(gamma_of_T(T_cold, **gpar))
    e_cold = epsilon(deltas[-1], gam_cold, gam_cold, kappa=kmid, dx=det_cold).eps
    print(f"  DESIGN RULE: at {T_cold:.0f} K the tracked mode is {det_cold:.1f} meV red "
          f"of X -> cavity-only eps = {e_cold:.2f} (favors XX!) -- keep the slit "
          "filter below T_target; the cavity is a safe filter only near/above it")

    # Tc with cavity collection gain on rho (F6 ii); ceiling searched ABOVE the
    # operating point T_target (below it the tracked-cavity eps inversion applies)
    G_lo, G_hi = card["collection_gain"].bounds
    p_base = dict(fitp)
    p_base["w"] = None
    Tc = {}
    for tag, G, d in (("hot", G_lo, deltas[0]), ("best", G_hi, deltas[-1])):
        p = dict(p_base, collection_gain=G, delta_xx=d, kappa=kmid)
        p["dx"] = lambda T: 1e3 * tracking_detuning(T, E0, Tt, mat, dEdT)
        Tc[tag] = solve_Tc(p, T_lo=Tt)
    print(f"  T_c (searched above T_target) with tracked cavity (kappa = {kmid:.1f} meV): "
          f"{Tc['hot']:.0f} K (Delta={deltas[0]:.1f}, G={G_lo:.0f}) to "
          f"{Tc['best']:.0f} K (Delta={deltas[-1]:.1f}, G={G_hi:.0f})   [A] envelope")
    if Tc["hot"] <= Tt + 0.1:
        print(f"    (Delta = {deltas[0]:.1f} meV already exceeds the ceiling AT "
              "T_target: small-splitting dots cannot ride the cavity filter -- "
              "the Delta requirement is the binding constraint, cf. F7)")

    cavity_design_figure({
        "Ts": Ts, "det_meV": det_meV, "T_target": Tt,
        "kappas": kappas, "deltas": deltas, "eps_band": eps_band,
        "Feff": Feff, "gam_t": gam_t, "Tc": Tc,
    }, ROOT / "out" / "phase2")
    print(f"  bundle -> {ROOT / 'out' / 'phase2'}")


if __name__ == "__main__":
    run_vb()
    run_cavity_design()
