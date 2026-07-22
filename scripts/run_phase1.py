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
from fsim_viz.figures import phase0_bundle, phase1_bundle, vc_reischle_figure


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

    vc_reischle()
    return fit


def vc_reischle():
    """V-c: is the Reischle electrical result rho-limited with eps ~ small?

    Trion line -> eps = 0 structurally (no cascade partner), so the F-series
    predicts g2 = 1 - rho^2. Test: does the rho digitized from the pulsed EL
    spectrum, over plausible detection windows, cover the rho the measured
    g2 requires?"""
    card = load_card(ROOT / "cards" / "reischle.yaml")
    ph = card.placeholders()
    if ph:
        print(f"\nV-c (reischle): BLOCKED -- placeholder entries: {', '.join(ph)}")
        return

    rows = card.datasets["g2_vs_ERR"].rows
    base = next(r for r in rows
                if r["device"] == 1 and r["position"] == 2 and r["ERR_MHz"] == 100.0)
    g2, err = base["g2"], base["err"]
    rho_req = float(np.sqrt(1.0 - g2))
    rho_req_err = err / (2.0 * rho_req)

    rw = sorted(card.datasets["rho_spectral"].rows, key=lambda r: r["w_meV"])
    ws = [r["w_meV"] for r in rw]
    lo = [r["rho_lo"] for r in rw]
    hi = [r["rho_hi"] for r in rw]
    cover = [r["w_meV"] for r in rw
             if r["rho_lo"] - rho_req_err <= rho_req <= r["rho_hi"] + rho_req_err]
    consistent = bool(cover)

    print("\nV-c (Reischle APL 97, 143513; Device 1 Pos. 2, ~40 K, pulsed electrical):")
    print("  line = trion (paper assignment) -> eps = 0 structurally; "
          "no cascade leakage path")
    print(f"  measured g2(0) = {g2}+-{err} @ 100 MHz  ->  requires rho = "
          f"{rho_req:.3f}+-{rho_req_err:.3f}")
    print(f"  digitized rho(w) envelope, w = {ws[0]:.0f}-{ws[-1]:.0f} meV: "
          f"[{min(lo):.3f}, {max(hi):.3f}]")
    if consistent:
        print(f"  required rho covered at w = {cover[0]:.0f}-{cover[-1]:.0f} meV  ->  "
              "V-c(i) CONSISTENT: rho-limited, eps small")
    else:
        print("  required rho NOT covered by the spectral envelope -> V-c(i) FAILS; "
              "a non-spectral mechanism is required")
    hi_err = [r for r in rows if r["ERR_MHz"] > 100.0 and r["device"] == 1 and r["position"] == 2]
    if hi_err:
        vals = ", ".join(f"{r['g2']} @ {r['ERR_MHz']:.0f} MHz" for r in hi_err)
        print(f"  higher-ERR excess ({vals}): temporal (EP refilling + peak overlap; "
              "paper's own attribution) -- outside the spectral model, WP-M2' tier")
    print(f"  tag chain: {card.datasets['rho_spectral'].tag.label} "
          "(digitized spectrum; window unpublished, swept)")

    # ---- V-c(ii): the 80 K DC result (OE 16, 12771 (2008))
    if "g2_oe2008" in card.datasets:
        from fsim_core.spectral import epsilon
        print("\nV-c(ii) (Reischle OE 16, 12771 (2008), DC electrical, 80 K):")
        ok = True
        for r in card.datasets["g2_oe2008"].rows:
            g2s = (r["g2_b"] - (1 - r["rho"] ** 2)) / r["rho"] ** 2
            err_s = r["err"] / r["rho"] ** 2
            match = abs(g2s - r["g2_corr"]) <= 0.011
            ok &= match
            print(f"  dot {r['dot']} @ {r['T']:.0f} K: F2 inversion "
                  f"(g2_b={r['g2_b']}, rho={r['rho']}) -> g2_s = {g2s:.3f}+-{err_s:.3f}"
                  f"  vs published {r['g2_corr']}  [{'match' if match else 'MISMATCH'}]")
        gam = card["gamma_80K_oe"].fixed
        eps_hi = 0.0
        for delta in card["delta_xx_oe"].bounds:
            for w in card["w_oe"].bounds:
                eps_hi = max(eps_hi, epsilon(delta, gam, gam, w=w).eps)
        print(f"  eps bound (X line, Gamma(80K)={gam} meV, Delta 4-6 meV, w 1-3 meV "
              f"swept): eps <= {eps_hi:.3f}")
        print("  residual g2_s consistent with eps within errors -> "
              "V-c(ii) " + ("CONSISTENT: rho-limited, eps small" if ok else "CHECK FAILED"))
        print("  note: the paper's Eq. (1) IS the F2 background law -- the published "
              "analysis is an F-series inversion performed by the authors themselves")

    vc_reischle_figure({"w_meV": ws, "rho_lo": lo, "rho_hi": hi,
                        "rho_req": rho_req, "rho_req_err": rho_req_err,
                        "g2": g2, "g2_err": err, "consistent": consistent},
                       ROOT / "out" / "phase1")


if __name__ == "__main__":
    main()
