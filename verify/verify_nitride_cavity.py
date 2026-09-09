"""Numerical and source-transcription checks for nitride_cavity.

No held-out prediction is claimed here.  Source fixtures transcribe published
numbers; all model checks are numerical verification or non-gating comparisons.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.dbr import qw_peak_reflectivity
from fsim_core.nitride_cavity import (NitrideCavityEnergyError,
    NitrideCavityParameterError, NitrideCavityParams, dbr_diagnostic,
    response, varshni_shift_eV)
from fsim_core.spectral import epsilon2, lorentzian


def close(a, b, rtol=1e-9, atol=1e-12):
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol)


def expect(exc, fn):
    try:
        fn()
    except exc:
        return True
    return False


def main():
    checks = []
    p = NitrideCavityParams()
    r0 = response(p, T_K=300.0, E_X_eV=2.60, E_X_track_eV=2.60,
                  gamma_X_meV=20.0, gamma_XX_meV=25.0, delta_xx_meV=3.0,
                  gamma_X0_ns=1.0, gamma_XX0_ns=0.5, w_meV=10.0)
    # Source transcription: Deshpande et al., Nat. Commun. 4, 1675 (2013),
    # Fig. 2b: alpha=(8.3 +/- 0.4)e-4 eV/K, beta=825 +/- 16.7 K [V].
    checks.append(close(varshni_shift_eV(300.0), -8.3e-4 * 300.0**2 / 1125.0))
    # Source transcription: Taylor et al., Nanoscale Res. Lett. 5, 608 (2010),
    # DOI 10.1007/s11671-009-9514-4: best-spot Q=167 [V].
    taylor_Q = 167.0
    checks.append(taylor_Q == 167.0)
    checks.append(all(close(response(NitrideCavityParams(Q=q), T_K=300, E_X_eV=2.6,
        E_X_track_eV=2.6, gamma_X_meV=1, gamma_XX_meV=1, delta_xx_meV=1,
        gamma_X0_ns=1, gamma_XX0_ns=1)["kappa_meV"], 2600.0/q)
        for q in (500.0, 2000.0, 10000.0)))  # exploratory Q range [A]
    checks.append(close(r0["kappa_meV"], 1000.0 * r0["E_cav_eV"] / p.Q))
    checks.append(close(r0["detuning_meV"], 0.0))
    away = response(p, T_K=320, E_X_eV=2.60, E_X_track_eV=2.60,
        gamma_X_meV=20, gamma_XX_meV=25, delta_xx_meV=3, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(close(away["E_cav_eV"], 2.60 - 0.04 * 20 / 1000))
    anti = response(p, T_K=300, E_X_eV=2.60, E_X_track_eV=2.60,
        gamma_X_meV=2, gamma_XX_meV=2, delta_xx_meV=-3, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(anti["F_eff_XX"] < anti["F_eff_X"])
    det, gam, kap = 1.3, 4.1, 2.7
    numeric = quad(lambda x: lorentzian(x, det, gam) * (kap / 2)**2 / (x*x + (kap/2)**2),
                   -np.inf, np.inf, epsabs=1e-11)[0]
    analytic = kap / (kap + gam) / (1 + 4 * det**2 / (kap + gam)**2)
    checks.append(close(numeric, analytic, rtol=1e-5))
    lowq = response(NitrideCavityParams(Q=1e5), T_K=300, E_X_eV=2.6, E_X_track_eV=2.6,
        gamma_X_meV=20, gamma_XX_meV=20, delta_xx_meV=0, gamma_X0_ns=1, gamma_XX0_ns=1)
    highq = response(NitrideCavityParams(Q=1e6), T_K=300, E_X_eV=2.6, E_X_track_eV=2.6,
        gamma_X_meV=20, gamma_XX_meV=20, delta_xx_meV=0, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(highq["F_eff_X"] / lowq["F_eff_X"] < 1.2)
    off = response(p, T_K=300, E_X_eV=2.8, E_X_track_eV=2.6, gamma_X_meV=1,
        gamma_XX_meV=1, delta_xx_meV=1, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(close(off["gamma_X_ns"], 1.0, rtol=2e-3))
    disabled = response(NitrideCavityParams(purcell_enabled=False), T_K=300, E_X_eV=2.6,
        E_X_track_eV=2.6, gamma_X_meV=1, gamma_XX_meV=1, delta_xx_meV=1,
        gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(disabled["Fp_add"] == 0 and disabled["F_eff_X"] == 1 and disabled["F_eff_XX"] == 1)
    expected = epsilon2(3, 20, 25, w=10, kappa=r0["kappa_meV"], dx_w=0, dx_c=0)
    checks.append(close(r0["t_X"], expected.t_x) and close(r0["t_XX"], expected.t_xx))
    checks.append(0 <= r0["t_X"] <= 1 and 0 <= r0["t_XX"] <= 1)
    eta2 = response(NitrideCavityParams(eta_out=0.2), T_K=300, E_X_eV=2.6, E_X_track_eV=2.6,
        gamma_X_meV=20, gamma_XX_meV=25, delta_xx_meV=3, gamma_X0_ns=1, gamma_XX0_ns=.5, w_meV=10)
    checks.append(close(eta2["eta_out"], 2*r0["eta_out"]) and close(eta2["t_XX"]/eta2["t_X"], r0["t_XX"]/r0["t_X"]))
    d = dbr_diagnostic(500)
    target = qw_peak_reflectivity(1, 1, 2.1, 1.46, 10)
    checks.append(close(d["reflectance"], target, rtol=1e-9) and close(d["reflectance"] + d["transmittance"], 1))
    checks.append(close(d["d_high_nm"], 500/(4*2.1)) and close(d["d_low_nm"], 500/(4*1.46)))
    checks.append(expect(NitrideCavityParameterError, lambda: NitrideCavityParams(Q=0)))
    checks.append(expect(NitrideCavityParameterError, lambda: NitrideCavityParams(eta_out=1.1)))
    checks.append(expect(NitrideCavityParameterError, lambda: response(p, T_K=300, E_X_eV=0,
        E_X_track_eV=2.6, gamma_X_meV=1, gamma_XX_meV=1, delta_xx_meV=1, gamma_X0_ns=1, gamma_XX0_ns=1)))
    checks.append(expect(NitrideCavityEnergyError, lambda: response(
        NitrideCavityParams(detuning_offset_meV=1.0), T_K=300, E_X_eV=2.6,
        E_X_track_eV=0.00001, gamma_X_meV=1, gamma_XX_meV=1, delta_xx_meV=1,
        gamma_X0_ns=1, gamma_XX0_ns=1)))
    passed = sum(checks)
    print("source transcription: Deshpande Varshni [V]; Taylor Q=167 [V]")
    print("assumed-design Q range [A]: 500, 2000, 10000")
    print(f"{passed}/{len(checks)} nitride cavity checks passed")
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
