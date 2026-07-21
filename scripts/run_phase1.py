"""Phase 1 driver: drive realism (Module D) + analytic thermal (Module A).

1. V-a recheck with the F1b operating point (mu_op) folded in.
2. F1b drive-penalty curve (figure + CSV).
3. Delta T_J envelope maps for the staged-device mesa on GaAs and GaAs/Si
   (the section-4 deliverable; epi-k range swept, never averaged).
4. V-c (Reischle): runs only when the card carries real data.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card
from fsim_core.fitting import V_A_TOL, fit_phase0
from fsim_core.loading import f1b_g2
from fsim_core.thermal import Layer, Stack, t_junction
from fsim_viz.figures import phase0_bundle, phase1_bundle


def va_recheck():
    card_path = ROOT / "cards" / "chatzarakis.yaml"
    card = load_card(card_path)
    fit = fit_phase0(card)
    print(f"V-a recheck with F1b (mu_op = {card['mu_op'].fixed}):")
    worst = float(np.max(np.abs(fit.residuals)))
    print(f"  max |resid| = {worst:.3f}  ->  {'PASS' if fit.passed else 'FAIL'} at +-{V_A_TOL}")
    gA, gAv, gAt, T_anchor = fit.gamma_anchor
    print(f"  Gamma({T_anchor:.0f} K) = {gA:.2f} meV (anchor {gAv}+-{gAt});"
          f"  E_a = {fit.params['E_a']:.1f} meV")
    res = phase0_bundle(fit, card_path, ROOT / "out" / "phase0")
    print(f"  Tc = {res['Tc']:.1f} K  (phase0 bundle regenerated with mu_op)")
    return fit


def drive_curves(mu_op):
    mus = np.logspace(-2, 1.5, 200)
    return {"mu": mus,
            "curves": {e: f1b_g2(mus, e) / e for e in (0.01, 0.1, 0.3)},
            "mu_op": mu_op}


def thermal_maps():
    card = load_card(ROOT / "cards" / "qcap-staged.yaml")
    d_lo, d_hi = card["mesa_diameter"].bounds
    d_um = np.linspace(d_lo, d_hi, 60)
    k_lo, k_hi = card["epi_k"].bounds
    t_epi = card["epi_thickness"].fixed * 1e-6
    t_buf = card["buffer_thickness"].fixed * 1e-6
    k_buf = card["buffer_k"].fixed
    V = card["V_drive"].fixed
    duty = card["duty"].fixed
    T_hs = card["T_heatsink"].fixed
    currents = [10.0, 50.0, 200.0]  # uA, spanning the I_drive card range

    a_epi = card["epi_alpha"].fixed

    def stacks(k_epi):
        pillar = [Layer(t_epi, k_epi, a_epi, spread=False)]  # mesa etched through epi
        return {
            "GaAs": Stack(layers=pillar, k_sub300=55.0, sub_alpha=1.25),
            "GaAs-Si": Stack(layers=pillar + [Layer(t_buf, k_buf, 1.25, spread=True)],
                             k_sub300=148.0, sub_alpha=1.35),
        }

    templates = {}
    for tpl in ("GaAs", "GaAs-Si"):
        per_I = {}
        for I in currents:
            P = duty * I * 1e-6 * V
            dT = {}
            for k_epi in (k_lo, k_hi):
                st = stacks(k_epi)[tpl]
                dT[k_epi] = np.array([
                    t_junction(P, 0.5 * d * 1e-6, st, T_hs) - T_hs for d in d_um
                ])
            # k range swept -> envelope band (low k = hot bound)
            per_I[I] = (np.minimum(dT[k_lo], dT[k_hi]), np.maximum(dT[k_lo], dT[k_hi]))
        templates[tpl] = per_I
    return {"diam_um": d_um, "currents_uA": currents, "templates": templates,
            "T_hs": T_hs}, card


def main():
    fit = va_recheck()

    drive = drive_curves(mu_op=0.33)
    thermal, qcard = thermal_maps()

    print("\nDelta T_J envelope (hot bound = low epi-k), T_hs = "
          f"{thermal['T_hs']:.0f} K, tag chain [A]:")
    for tpl, per_I in thermal["templates"].items():
        d = thermal["diam_um"]
        for I in thermal["currents_uA"]:
            hi = per_I[I][1]
            top = "RUNAWAY" if np.isinf(hi[0]) else f"{hi[0]:6.1f} K"
            print(f"  {tpl:8s} {I:5.0f} uA:  dTj({d[0]:.2f} um) up to {top}"
                  f"   dTj({d[-1]:.1f} um) up to {hi[-1]:5.2f} K")

    res = phase1_bundle(drive, thermal, ROOT / "out" / "phase1")
    print(f"\nphase1 bundle -> {res['outdir']}")

    rcard = load_card(ROOT / "cards" / "reischle.yaml")
    ph = rcard.placeholders()
    if ph:
        print(f"\nV-c (reischle): BLOCKED -- placeholder entries: {', '.join(ph)}")
        print("   Drop the Reischle APL 97, 143513 (2010) PDF into the project root "
              "and re-run to close V-c.")
    return fit


if __name__ == "__main__":
    main()
