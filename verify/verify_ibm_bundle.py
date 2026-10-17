"""Verifier for the opt-in IBM-calibrated V-a fit bundle (devlib-r2-pha2, D1).

Anchors are numbers the code did not produce: the closed Lorentzian FWHM / transmission in the weak-coupling limit, and
the r2-physics-a Q3 probe tables (IBM FWHM at ZPL = fitted Gamma, IBM eps, the b0-only recomposition), which Fable
computed with a separate script.  Class (C) calibration: the bundle is evaluated on the data it was fit to; the
verifier checks wiring, limits and bit-identity of the Lorentzian default, not held-out validity.

Exit 0 iff all checks pass; prints "N/N ibm_bundle checks passed".
"""
import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import fitting as F  # noqa: E402
from fsim_core.card import load_card  # noqa: E402
from fsim_core.device import DeviceDesign, class_proxy_params, evaluate  # noqa: E402
from fsim_core.fit_bundle import fit_params_path, load_fit_params  # noqa: E402
from fsim_core.integrator import g2_from, retention  # noqa: E402
from fsim_core.loading import f1b_g2  # noqa: E402
from fsim_core.qd_gf import PhononParams, ibm_transmission  # noqa: E402
from fsim_core.spectral import epsilon, gamma_of_T  # noqa: E402

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)


card = load_card(ROOT / "cards" / "chatzarakis.yaml")
rows = sorted(card.datasets["g2_vs_T"].rows, key=lambda r: r["T"])
pp = PhononParams(**F.chatz2023_phonon())
LOR = load_fit_params()["params"]
IBM = load_fit_params("ibm")
BP = IBM["params"]


def gam_of(P, T):
    return float(gamma_of_T(T, P["gamma0"], P["a_ac"], P["b_lo"], P["E_lo"]))


# 1. estimator anchor: weak-coupling limit, FWHM of the IBM line = the Lorentzian FWHM (closed form) [DR]
weak = PhononParams(alpha_ps2=1e-9)
errs = [abs(F.ibm_fwhm(weak, 78.0, g) - g) / g for g in (0.2, 1.0, 5.0)]
check("fwhm estimator: weak-coupling limit equals the Lorentzian FWHM", max(errs) < 2e-3, f"max rel err {max(errs):.1e}")

# 2. weak coupling: narrowed ZPL = Gamma and eps = the closed-form spectral.epsilon
r0 = rows[3]
gam = 4.0
e, _, gz = F.ibm_exact_eps(weak, r0["T"], r0["w"], r0["dx"], 5.9, gam, 0.72)
e_cf = epsilon(5.9, gam, 0.72 * gam, w=r0["w"], dx=r0["dx"]).eps
check("weak coupling: gamma_zpl = Gamma and eps = spectral.epsilon closed form",
      abs(gz - gam) < 1e-2 * gam and abs(e - e_cf) < 2e-3 * e_cf, f"gz {gz:.3f} eps {e:.4f} vs {e_cf:.4f}")

# 3-5. r2-physics-a Q3 probe tables at the LORENTZIAN bundle (ZPL = fitted Gamma, card chatzarakis2023 phonons)
gam_T = [gam_of(LOR, r["T"]) for r in rows]
fw_ref = [1.34, 3.45, 4.87, 5.82, 7.73, 8.70]   # [DR] probe, r2-physics-a Q3 (meV, 3 digits)
fw = [F.ibm_fwhm(pp, r["T"], g) for r, g in zip(rows, gam_T)]
d_fw = max(abs(a - b) for a, b in zip(fw, fw_ref))
check("IBM FWHM at ZPL = fitted Gamma matches the Q3 probe table (1.34/3.45/4.87/5.82/7.73/8.70)", d_fw < 0.02,
      f"max |d| {d_fw:.3f} meV")
eps_full = []
for r, g in zip(rows, gam_T):
    tx = ibm_transmission(r["dx"], pp, r["T"], g, w_meV=r["w"])
    txx = ibm_transmission(r["dx"] - LOR["delta_xx"], pp, r["T"], LOR["r_xx"] * g, w_meV=r["w"])
    eps_full.append(txx / tx)
eps_ref = [0.0088, 0.048, 0.090, 0.124, 0.196, 0.232]   # [DR] probe, Q3 (IBM eps with the double-counted ZPL)
d_eps = max(abs(a - b) for a, b in zip(eps_full, eps_ref))
check("IBM eps at ZPL = fitted Gamma matches the Q3 probe (0.0088/0.048/0.090/0.124/0.196/0.232)", d_eps < 1.5e-3,
      f"max |d| {d_eps:.4f}")
g2_ref = [0.029, 0.079, 0.134, 0.177, 0.285, 0.361]     # [DR] probe, Q3: b0 = 0.009, retention/Gamma fixed, full-width ZPL
g2_mine = []
for r, e_ in zip(rows, eps_full):
    S = float(retention(r["T"], LOR["a_esc"], LOR["E_a"], LOR["b_p"], LOR["E_b"]))
    g2_mine.append(float(g2_from(float(f1b_g2(0.33, e_)), S / (S + 0.009 + LOR["beta"] * (1 - S)))))
d_g2 = max(abs(a - b) for a, b in zip(g2_mine, g2_ref))
check("b0-only recomposition (Q3 closed-form arithmetic, b0 = 0.009) matches the six probe numbers", d_g2 < 0.006,
      f"max |d| {d_g2:.4f}")

# 6. physics sign: the narrowed ZPL is narrower than Gamma; FWHM_IBM(narrowed ZPL) = Gamma
gz_b = IBM["gamma_zpl_meV"]
gam_b = [gam_of(BP, r["T"]) for r in rows]
check("narrowed ZPL < fitted Gamma at every T (the full-width ZPL double counts phonon broadening)",
      all(gz_b[str(int(r["T"]))] < g for r, g in zip(rows, gam_b)))
check("FWHM_IBM(narrowed ZPL) = Gamma(T) at the bundle solution (all six T, 5 ueV)",
      all(abs(F.ibm_fwhm(pp, r["T"], gz_b[str(int(r["T"]))]) - g) < 5e-3 for r, g in zip(rows, gam_b)))

# 7. bundle self-consistency: recompute g2 from the stored parameters with exact eps (no spline)
g2_re = []
for r, g in zip(rows, gam_b):
    e, _, _ = F.ibm_exact_eps(pp, r["T"], r["w"], r["dx"], BP["delta_xx"], g, BP["r_xx"])
    S = float(retention(r["T"], BP["a_esc"], BP["E_a"], BP["b_p"], BP["E_b"]))
    g2_re.append(float(g2_from(float(f1b_g2(BP["mu"], e)), S / (S + BP["b0"] + BP["beta"] * (1 - S)))))
stored = [IBM["model_g2"][str(int(r["T"]))] for r in rows]
check("stored bundle model_g2 reproduced from the stored parameters (exact eps)",
      max(abs(a - b) for a, b in zip(g2_re, stored)) < 1e-6)
check("spline-vs-exact error of the fit loop is small (< 1e-3 in g2)", IBM["ibm"]["spline_vs_exact_max_abs_dg2"] < 1e-3,
      "%.1e" % IBM["ibm"]["spline_vs_exact_max_abs_dg2"])
check("78 K published INEQUALITY met within its +-0.01 (<= 0.03); weighted err/3 [A] recorded",
      stored[0] <= 0.03 and IBM["weights"]["upper_bound_err_divisor"] == 3.0, f"g2(78) = {stored[0]:.4f}")

# 8. device recomposition on the bundle == the fit (IBM path wiring, opt-in attribute)
dev = yaml.safe_load((ROOT / "cards" / "chatzarakis2023.yaml").read_text())["device"]


def run_dev(bundle, pt):
    d = DeviceDesign()
    d.drive.V = d.drive.I_uA = d.drive.duty = d.drive.b_e = 0.0
    d.drive.mode = "PL"
    d.cavity.enabled = False
    d.filter.enabled, d.filter.auto_w = True, False
    d.dot.lineshape = "ibm"
    d.dot.phonon = {**dev["geometry"], "alpha_ps2": dev["phonon_alpha_ps2"]}
    if bundle is not None:
        d.dot.fit_bundle = bundle
        if bundle == "ibm":
            d.dot.r_xx = BP["r_xx"]
    d.thermal.T_hs, d.filter.w, d.filter.dx = pt["T"], pt["w"], pt["dx"]
    d.drive.mu, d.dot.delta_xx = dev["mu"], dev["delta_xx"]
    return float(evaluate(d, T_grid=[pt["T"]])["scalars"]["g2_op"])


dev_ibm = [run_dev("ibm", r) for r in rows]
d_dev = max(abs(a - b) for a, b in zip(dev_ibm, stored))
check("device.evaluate on the IBM bundle reproduces the fitted g2 at all six points (< 2e-4)", d_dev < 2e-4,
      f"max |d| {d_dev:.1e}")
leg_none = [run_dev(None, r) for r in rows]
leg_lor = [run_dev("lorentzian", r) for r in rows]
check("fit_bundle 'lorentzian' == attribute absent: bit-identical legacy device path", leg_none == leg_lor)
check("class_proxy_params() default is the Lorentzian bundle; 'ibm' differs",
      class_proxy_params()["b0"] == LOR["b0"] and class_proxy_params("ibm")["b0"] == BP["b0"] and BP["b0"] != LOR["b0"])

# 9. Lorentzian fit untouched: default fit_phase0 reproduces the tracked bundle
fit0 = F.fit_phase0(card)
snap = json.loads(fit_params_path().read_text(encoding="utf-8"))["params"]
check("fit_phase0 default (Lorentzian) reproduces the tracked bundle (1e-6 rel)",
      all(abs(float(fit0.params[k]) - snap[k]) <= 1e-6 * max(abs(snap[k]), 1e-12) + 1e-12 for k in snap if k in fit0.params))
check("bundle lineshape tag and cost recorded", IBM["lineshape"] == "ibm" and IBM["cost"] > 0)

n = len(RESULTS)
print(f"{sum(RESULTS)}/{n} ibm_bundle checks passed")
sys.exit(0 if all(RESULTS) else 1)
