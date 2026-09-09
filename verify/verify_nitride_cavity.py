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
from fsim_core.dbr import power_RT, qw_peak_reflectivity, quarter_wave_stack
from fsim_core.nitride_cavity import (NitrideCavityEnergyError,
    NitrideCavityParameterError, NitrideCavityParams, _line_overlap,
    dbr_diagnostic, response, varshni_shift_eV)
from fsim_core.spectral import cavity_transmission, epsilon2, lorentzian


def close(a, b, rtol=1e-9, atol=1e-12):
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol)


def pinned(a, b):
    """Strict pin: literal must reproduce the model output to 1e-12 absolute,
    not a relative tolerance -- catches an unintended physics-result change."""
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-12)


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
    checks.append(all(close(response(NitrideCavityParams(Q=q), T_K=300, E_X_eV=2.6,
        E_X_track_eV=2.6, gamma_X_meV=1, gamma_XX_meV=1, delta_xx_meV=1,
        gamma_X0_ns=1, gamma_XX0_ns=1)["kappa_meV"], 2600.0/q)
        for q in (500.0, 2000.0, 10000.0)))  # exploratory Q range [A]

    # Regression pin: commit 191ff57's response() at the r0 fixture, to
    # 1e-12 absolute.  This hardening round must not move these numbers;
    # a change here means the physics changed, not just its verification.
    checks.append(pinned(r0["Fp_add"], 75.99088773175333))
    checks.append(pinned(r0["F_eff_X"], 5.637941504754898))
    checks.append(pinned(r0["F_eff_XX"], 4.5703781164350845))
    checks.append(pinned(r0["t_X"], 0.0586512835666765))
    checks.append(pinned(r0["t_XX"], 0.04474200146334138))

    checks.append(close(r0["kappa_meV"], 1000.0 * r0["E_cav_eV"] / p.Q))
    checks.append(close(r0["detuning_meV"], 0.0))
    away = response(p, T_K=320, E_X_eV=2.60, E_X_track_eV=2.60,
        gamma_X_meV=20, gamma_XX_meV=25, delta_xx_meV=3, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(close(away["E_cav_eV"], 2.60 - 0.04 * 20 / 1000))
    anti = response(p, T_K=300, E_X_eV=2.60, E_X_track_eV=2.60,
        gamma_X_meV=2, gamma_XX_meV=2, delta_xx_meV=-3, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(anti["F_eff_XX"] < anti["F_eff_X"])

    # Nonzero-detuning XX-binding-sign check.  At zero cavity detuning the
    # broadened overlap is even in line detuning, so a sign bug
    # (detuning_xx = detuning + delta_xx instead of detuning - delta_xx)
    # would pass the anti-binding check above undetected.  Here E_X sits
    # 10 meV above the cavity and delta_xx=+10 meV puts XX exactly on
    # resonance (detuning_xx=0) while X stays 10 meV detuned, so a wrong
    # sign would swap which line is enhanced far more strongly.
    xxsign = response(NitrideCavityParams(), T_K=300.0, E_X_eV=2.71,
        E_X_track_eV=2.70, gamma_X_meV=2.0, gamma_XX_meV=2.0,
        delta_xx_meV=10.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    kappa_xs = 1000.0 * xxsign["E_cav_eV"] / NitrideCavityParams().Q
    fp_xs = ((3.0 / (4.0 * math.pi**2)) * NitrideCavityParams().Q
             / NitrideCavityParams().mode_volume_norm
             * NitrideCavityParams().spatial_overlap)
    fx_expected = 1.0 + fp_xs * cavity_transmission(xxsign["detuning_meV"], 2.0, kappa_xs)
    fxx_expected = 1.0 + fp_xs * cavity_transmission(
        xxsign["detuning_meV"] - 10.0, 2.0, kappa_xs)
    checks.append(close(xxsign["F_eff_X"], fx_expected)
                  and close(xxsign["F_eff_XX"], fxx_expected))
    # Anchored to the reviewer-quoted reference numbers for this fixture.
    checks.append(math.isclose(xxsign["F_eff_XX"], 31.6, rel_tol=2e-2)
                  and math.isclose(xxsign["F_eff_X"], 1.84, rel_tol=2e-2)
                  and xxsign["F_eff_XX"] > xxsign["F_eff_X"])

    # Purcell magnitude: Fp_add against the closed form (3/(4pi^2))*Q/V*overlap
    # [DR, Purcell, Phys. Rev. 69, 681 (1946)] at probe values distinct from
    # every default, so a corrupted exponent (e.g. mode_volume_norm**2) or a
    # dropped factor cannot hide behind a coincidental default-value match.
    q_probe, v_probe, overlap_probe = 733.0, 3.5, 0.62
    fp_closed = (3.0 / (4.0 * math.pi**2)) * q_probe / v_probe * overlap_probe
    r_purcell = response(NitrideCavityParams(Q=q_probe, mode_volume_norm=v_probe,
        spatial_overlap=overlap_probe), T_K=300, E_X_eV=2.6, E_X_track_eV=2.6,
        gamma_X_meV=20, gamma_XX_meV=20, delta_xx_meV=0, gamma_X0_ns=1, gamma_XX0_ns=1)
    checks.append(close(r_purcell["Fp_add"], fp_closed))

    # Broadened-overlap quadrature check: calls the module's own
    # _line_overlap (which itself delegates to spectral.cavity_transmission,
    # not a formula retyped here) and compares to independent numerical
    # convolution of the Lorentzian line with the Lorentzian cavity filter.
    det, gam, kap = 1.3, 4.1, 2.7
    numeric = quad(lambda x: lorentzian(x, det, gam) * (kap / 2)**2 / (x*x + (kap/2)**2),
                   -np.inf, np.inf, epsabs=1e-11)[0]
    checks.append(close(numeric, _line_overlap(det, gam, kap), rtol=1e-5))

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

    # Anchor checks: Taylor et al., Nanoscale Res. Lett. 5, 608 (2010), DOI
    # 10.1007/s11671-009-9514-4, best-spot Q=167 [V]; Dartsch et al.,
    # J. Cryst. Growth (2011), DOI 10.1016/j.jcrysgro.2010.10.091, Q=220
    # [V abstract-only]; Sebald et al., pss(b) 248, 1777 (2011), DOI
    # 10.1002/pssb.201147144, Q=260 [V abstract-only].  None of the three
    # reports observed Purcell enhancement, so each is checked here as a
    # Purcell-disabled reference case: kappa=E/Q must hold at the reported
    # Q, and the disabled channel must leave F_eff_X=1 (no held-out claim;
    # these are cited best-spot/abstract Q values, not a reproduced count).
    for anchor_q, anchor_label in ((167.0, "Taylor"), (220.0, "Dartsch"), (260.0, "Sebald")):
        anchor_r = response(NitrideCavityParams(Q=anchor_q, purcell_enabled=False),
            T_K=300, E_X_eV=2.6, E_X_track_eV=2.6, gamma_X_meV=20, gamma_XX_meV=20,
            delta_xx_meV=0, gamma_X0_ns=1, gamma_XX0_ns=1)
        checks.append(close(anchor_r["kappa_meV"], 1000.0 * anchor_r["E_cav_eV"] / anchor_q))
        checks.append(anchor_r["Fp_add"] == 0.0 and anchor_r["F_eff_X"] == 1.0)

    d = dbr_diagnostic(500)
    target = qw_peak_reflectivity(1, 1, 2.1, 1.46, 10)
    checks.append(close(d["reflectance"], target, rtol=1e-9) and close(d["reflectance"] + d["transmittance"], 1))
    checks.append(close(d["d_high_nm"], 500/(4*2.1)) and close(d["d_low_nm"], 500/(4*1.46)))
    # d["transmittance"] must be the transfer-matrix T from dbr.power_RT on
    # this same quarter-wave stack, not the vacuous 1-reflectance shortcut.
    layers = quarter_wave_stack(2.1, 1.46, 10, 500)
    _, t_direct = power_RT(layers, 500, 1.0, 1.0)
    checks.append(close(d["transmittance"], t_direct, rtol=1e-9))
    # dbr_diagnostic's n_in/n_out now take effect (previously hardcoded to
    # 1); a GaN-like exit index changes R+T=1 stays exact and R actually moves.
    d_gan = dbr_diagnostic(500, n_out=2.4)
    checks.append(close(d_gan["reflectance"] + d_gan["transmittance"], 1))
    checks.append(d_gan["reflectance"] != d["reflectance"])

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
    print("observed cavity Q examples [V], no Purcell reported: Taylor 167, Dartsch 220, Sebald 260")
    print(f"{passed}/{len(checks)} nitride cavity checks passed")
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
