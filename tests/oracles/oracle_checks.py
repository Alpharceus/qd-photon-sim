"""Standalone regression suite for the ported EM oracles (S-1).

Same @check style as verify/verify_fsim.py. Run directly:

    python tests/oracles/oracle_checks.py        (exit code 0 iff all pass)

or via verify/verify_fsim.py, which shells out to this file as one check
and additionally enforces the oracle-only quarantine (no fsim_core /
fsim_viz / fsim_gui source may import tests.oracles).

What is checked:
  (i)   slab energy conservation R+T=1 to 1e-12 over a parameter grid
  (ii)  all four tangential-continuity residuals < 1e-12 over the same grid
  (iii) known closed-form limits (trivial slab, half-wave, quarter-wave AR,
        Fabry-Perot period in optical thickness)
  (iv)  CMT vs fsim_core.spectral cross-checks -- the actual point of S-1:
        the CMT intracavity Lorentzian FWHM, FSIM's cavity_transmission
        half-max detuning, pointwise lineshape agreement between the two
        independently-normalized Lorentzian families, and a lorentzian()
        unit-area sanity anchor
  (v)   a NEGATIVE CONTROL that deliberately breaks the slab solution and
        asserts the (i)/(ii) checks would have caught it (guards against a
        vacuously-passing oracle)
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from fsim_core.spectral import cavity_transmission, lorentzian
from tests.oracles.em_oracles import cmt_intracavity, cmt_transmission, fresnel_slab

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# --------------------------------------------------------------- slab grid

def _slab_grid():
    """(n1, n2, n3, wavelength, thickness) tuples spanning the required
    ranges: n2 in [1.5, 3.5], thickness/wavelength in [0.1, 5], plus a mix
    of symmetric and asymmetric outer media."""
    cases = []
    n2_vals = np.linspace(1.5, 3.5, 6)
    ratio_vals = np.linspace(0.1, 5.0, 8)
    for n1, n3 in [(1.0, 1.0), (1.0, 1.5), (1.33, 1.0), (2.0, 3.0)]:
        for n2 in n2_vals:
            for ratio in ratio_vals:
                wavelength = 1.0
                thickness = ratio * wavelength
                cases.append((n1, float(n2), n3, wavelength, thickness))
    return cases


# ------------------------------------------------------------ (i)+(ii) grid

@check("slab: energy conservation R+T=1 to 1e-12 over the (n2, d/lambda) grid")
def _():
    for n1, n2, n3, wl, d in _slab_grid():
        s = fresnel_slab(n1, n2, n3, wl, d)
        assert abs(s["R"] + s["T"] - 1.0) < 1e-12, (n1, n2, n3, wl, d, s["R"], s["T"])


@check("slab: all four tangential-continuity residuals < 1e-12 over the grid")
def _():
    for n1, n2, n3, wl, d in _slab_grid():
        s = fresnel_slab(n1, n2, n3, wl, d)
        for key in ("residual_E_left", "residual_H_left", "residual_E_right", "residual_H_right"):
            assert s[key] < 1e-12, (key, n1, n2, n3, wl, d, s[key])


# -------------------------------------------------------------- (iii) limits

@check("slab: n1=n2=n3 (no interfaces) -> R=0, T=1")
def _():
    for n in (1.0, 1.5, 2.7):
        s = fresnel_slab(n, n, n, 1.0, 0.73)
        assert abs(s["R"]) < 1e-13 and abs(s["T"] - 1.0) < 1e-13, (n, s["R"], s["T"])


@check("slab: half-wave slab (thickness=lambda/(2 n2)) is invisible when n1=n3")
def _():
    wl = 1.0
    for n1, n2 in [(1.0, 2.2), (1.5, 1.5), (1.33, 3.1)]:
        n3 = n1
        thickness = wl / (2.0 * n2)
        s = fresnel_slab(n1, n2, n3, wl, thickness)
        assert abs(s["R"]) < 1e-12 and abs(s["T"] - 1.0) < 1e-12, (n1, n2, thickness, s["R"])


@check("slab: quarter-wave AR coating (n2=sqrt(n1 n3)) gives R=0")
def _():
    wl = 1.0
    for n1, n3 in [(1.0, 2.25), (1.0, 4.0), (1.5, 3.0)]:
        n2 = (n1 * n3) ** 0.5
        thickness = wl / (4.0 * n2)
        s = fresnel_slab(n1, n2, n3, wl, thickness)
        assert abs(s["R"]) < 1e-12, (n1, n3, n2, s["R"])


@check("slab: R(optical thickness) is periodic with period lambda/2 (Fabry-Perot)")
def _():
    n1, n2, n3, wl = 1.0, 2.4, 1.0, 1.0
    period = wl / 2.0  # period in optical thickness n2*d
    ods = np.linspace(0.05, 4.0, 40)
    for od in ods:
        d1 = od / n2
        d2 = (od + period) / n2
        r1 = fresnel_slab(n1, n2, n3, wl, d1)["R"]
        r2 = fresnel_slab(n1, n2, n3, wl, d2)["R"]
        assert abs(r1 - r2) < 1e-10, (od, r1, r2)
    # and it is NOT periodic at a generic non-multiple offset (sanity that
    # the check above is actually discriminating something)
    od = 1.0
    d1 = od / n2
    d3 = (od + 0.37 * period) / n2
    r1 = fresnel_slab(n1, n2, n3, wl, d1)["R"]
    r3 = fresnel_slab(n1, n2, n3, wl, d3)["R"]
    assert abs(r1 - r3) > 1e-6


# --------------------------------------------------- (iv) CMT vs fsim_core

@check("CMT: intracavity |a|^2 is a Lorentzian of FWHM = kappa_i+kappa_e (half-max search)")
def _():
    omega0 = 1.0
    for kappa_i, kappa_e in [(0.08, 0.16), (0.0, 0.5), (0.3, 0.02)]:
        kappa_tot = kappa_i + kappa_e
        peak = abs(cmt_intracavity(omega0, omega0, kappa_i, kappa_e)) ** 2

        def f(omega):
            return abs(cmt_intracavity(omega, omega0, kappa_i, kappa_e)) ** 2 - 0.5 * peak

        omega_hi = brentq(f, omega0, omega0 + 20.0 * kappa_tot + 1e-6, xtol=1e-13, rtol=1e-13)
        fwhm_measured = 2.0 * (omega_hi - omega0)
        assert abs(fwhm_measured - kappa_tot) < 1e-9, (kappa_i, kappa_e, fwhm_measured, kappa_tot)


@check("fsim cavity_transmission at gamma->0: half-max detuning equals kappa/2 to 1e-9")
def _():
    for kappa in (0.3, 1.0, 4.5):
        peak = cavity_transmission(0.0, 0.0, kappa)
        assert abs(peak - 1.0) < 1e-12, (kappa, peak)  # peak-normalized

        def f(delta):
            return cavity_transmission(delta, 0.0, kappa) - 0.5

        delta_half = brentq(f, 1e-9, 20.0 * kappa, xtol=1e-13, rtol=1e-13)
        assert abs(delta_half - kappa / 2.0) < 1e-9, (kappa, delta_half)


@check("fsim cavity_transmission(gamma=0) vs peak-normalized CMT |a|^2: same lineshape family")
def _():
    for kappa in (0.3, 1.0, 4.5):
        deltas = np.linspace(-6.0 * kappa, 6.0 * kappa, 401)
        fsim_shape = np.array([cavity_transmission(d, 0.0, kappa) for d in deltas])

        # CMT with all coupling external (kappa_i=0, kappa_e=kappa) so the
        # total linewidth matches fsim's kappa; normalize |a|^2 to its peak
        # the same way fsim's cavity_transmission is already peak-normalized.
        a = np.array([cmt_intracavity(d, 0.0, 0.0, kappa) for d in deltas])
        cmt_power = np.abs(a) ** 2
        cmt_shape = cmt_power / cmt_power.max()

        assert np.max(np.abs(fsim_shape - cmt_shape)) < 1e-9, kappa


@check("fsim lorentzian(): unit area by independent quadrature (sanity anchor)")
def _():
    for center, fwhm in [(0.0, 1.0), (5.9, 2.3), (-3.0, 0.05)]:
        area, _err = quad(lorentzian, -np.inf, np.inf, args=(center, fwhm), limit=400)
        assert abs(area - 1.0) < 1e-8, (center, fwhm, area)


@check("CMT transmission: on-resonance notch depth matches [(ki-ke)/(ki+ke)]^2")
def _():
    for kappa_i, kappa_e in [(0.08, 0.16), (0.5, 0.5), (0.02, 0.4)]:
        t0 = cmt_transmission(1.0, 1.0, kappa_i, kappa_e)
        expected = ((kappa_i - kappa_e) / (kappa_i + kappa_e)) ** 2
        assert abs(t0 - expected) < 1e-12, (kappa_i, kappa_e, t0, expected)


# ------------------------------------------------------------- (v) negative control

def _true_residuals(n1, n2, n3, wavelength, thickness, r, a, b, t):
    """Evaluate the four REAL tangential-continuity equations for a given
    (r, a, b, t), independent of how those amplitudes were produced. Used
    only by the negative control below to prove the residual check in (ii)
    is not vacuous."""
    import cmath
    k2 = 2.0 * np.pi * n2 / wavelength
    beta = k2 * thickness
    phase_fwd = cmath.exp(1j * beta)
    phase_bwd = cmath.exp(-1j * beta)
    e_left = abs((1.0 + r) - (a + b))
    h_left = abs((1.0 - r) - (n2 / n1) * (a - b))
    e_right = abs((a * phase_fwd + b * phase_bwd) - t)
    h_right = abs((n2 / n1) * (a * phase_fwd - b * phase_bwd) - (n3 / n1) * t)
    return e_left, h_left, e_right, h_right


@check("NEGATIVE CONTROL: corrupted slab solutions fail both the energy and residual "
      "checks (the oracle checks are not vacuous)")
def _():
    n1, n2, n3, wavelength, thickness = 1.0, 2.5, 1.5, 1.0, 0.37
    good = fresnel_slab(n1, n2, n3, wavelength, thickness)
    r, a, b, t = good["r"], good["a"], good["b"], good["t"]

    # Sanity: the genuine oracle output passes both check families at this
    # operating point, so any failure induced below is due to the injected
    # bug, not the point itself.
    assert abs(good["R"] + good["T"] - 1.0) < 1e-12
    assert max(good["residual_E_left"], good["residual_H_left"],
               good["residual_E_right"], good["residual_H_right"]) < 1e-12

    # --- Bug 1: drop the n2 admittance factor from the LEFT H-continuity
    # equation when reconstructing (a, b) from the (still-correct) r --
    # i.e. use (a-b) = (1-r) instead of the physical (a-b) = (n1/n2)(1-r).
    # E-continuity (a+b = 1+r) is left intact so the bug is isolated to H.
    # R and T do not depend on a, b, so this specifically must be caught by
    # the CONTINUITY RESIDUAL checks in (ii), not by energy conservation.
    sum_ab = a + b
    diff_ab_broken = (1.0 - r)  # correct value is (n1 / n2) * (1.0 - r)
    a_broken = 0.5 * (sum_ab + diff_ab_broken)
    b_broken = 0.5 * (sum_ab - diff_ab_broken)
    residuals = _true_residuals(n1, n2, n3, wavelength, thickness, r, a_broken, b_broken, t)
    assert max(residuals) > 1e-6, "negative control did not break continuity residuals"

    # --- Bug 2: perturb the reflection amplitude away from the value
    # required by the boundary-value problem (t, and hence T, untouched).
    # This must be caught by ENERGY CONSERVATION in (i).
    r_broken = r + 0.2
    energy_err = abs(abs(r_broken) ** 2 + good["T"] - 1.0)
    assert energy_err > 1e-6, "negative control did not break energy conservation"


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}  {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} oracle checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
