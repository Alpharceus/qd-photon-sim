"""D1 (devlib-r2-pha2): write the opt-in IBM-calibrated V-a fit bundle fsim_core/data/phase0_fit_params_ibm.json.

    OMP_NUM_THREADS=2 python scripts/fit_ibm_bundle.py [--write]

The fit is fitting.fit_phase0(cards/chatzarakis.yaml, lineshape="ibm"): same blocks and priors as the Lorentzian fit,
eps from the IBM line with the ZPL narrowed so FWHM_IBM = Gamma(T) (no double-counted phonon broadening), and the
78 K published INEQUALITY weighted with err/3 [A].  Class (C) calibration; nothing is held out.  The Lorentzian
bundle (phase0_fit_params.json) and every legacy pin are untouched.  Without --write the result is only printed.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card  # noqa: E402
from fsim_core.fitting import UPPER_BOUND_ERR_DIVISOR, V_A_TOL, fit_phase0  # noqa: E402

OUT = ROOT / "fsim_core" / "data" / "phase0_fit_params_ibm.json"


def main():
    card = load_card(ROOT / "cards" / "chatzarakis.yaml")
    fit = fit_phase0(card, lineshape="ibm")
    ibm = fit.aux["ibm"]
    print(f"{'T':>6} {'w':>6} {'dx':>6} {'data':>7} {'model':>7} {'resid':>7}  gz[meV]")
    for r, m, e, gz in zip(fit.data, fit.model_g2, fit.residuals, ibm["gamma_zpl_meV"]):
        print(f"{r['T']:>6.0f} {r['w']:>6.2f} {r['dx']:>6.2f} {r['g2']:>7.3f} {m:>7.4f} {e:>+7.4f}  {gz:.3f}")
    print("params:", {k: float(v) for k, v in fit.params.items()})
    print("aux:", {k: v for k, v in fit.aux.items() if k != "ibm"}, "cost", fit.cost)
    print("passed +-0.03:", fit.passed)
    bundle = {
        "params": {k: float(v) for k, v in fit.params.items()},
        "lineshape": "ibm",
        "tag_chain": fit.tag.label,
        "passed_pm0.03": bool(fit.passed),
        "cost": float(fit.cost),
        "model_g2": {f"{int(r['T'])}": float(m) for r, m in zip(fit.data, fit.model_g2)},
        "residuals": {f"{int(r['T'])}": float(e) for r, e in zip(fit.data, fit.residuals)},
        "gamma_zpl_meV": {f"{int(t)}": g for t, g in zip(ibm["T"], ibm["gamma_zpl_meV"])},
        "weights": {"upper_bound_err_divisor": UPPER_BOUND_ERR_DIVISOR,
                    "note": "[A] the published 78 K inequality g2 < 0.02 is penalised with err/3 (one-sided), "
                            "not +-err (r2-physics-a Q3); all other blocks and priors identical to the "
                            "Lorentzian fit"},
        "aux": {k: v for k, v in fit.aux.items() if k != "ibm"},
        "ibm": {k: ibm[k] for k in ("phonon", "fwhm_root_bracketed", "spline_vs_exact_max_abs_dg2")},
        "V_A_TOL": V_A_TOL,
        "notes": fit.notes + ["Tc_K not computed: integrator.solve_Tc is the Lorentzian eps path (Tc is a "
                              "Lorentzian-bundle quantity)"],
        "source": "scripts/fit_ibm_bundle.py -> fitting.fit_phase0(lineshape='ibm'); class (C) calibration",
    }
    if "--write" in sys.argv:
        OUT.write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
        print("wrote", OUT)


if __name__ == "__main__":
    main()
