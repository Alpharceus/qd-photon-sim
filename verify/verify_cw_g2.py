"""Regression suite for fsim_core.cw_g2 (CW g2 of the filtered cascade,
background law, IRF convolution / deconvolution).

Run: python verify/verify_cw_g2.py   (exit code 0 iff all pass)
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.cw_g2 import (
    convolve_irf,
    cw_report,
    cw_vs_pulsed_note,
    deconvolve_dip,
    dip_convolved,
    dip_convolved_exp,
    escape_rates_from_retention,
    g2_cw,
    g2_cw_zero,
    g2_cw_zero_low_pump,
    g2_two_level,
    g2_with_background,
    generator,
    intrinsic_from_raw,
    irf_width,
    propagate,
    raw_g2_0,
    steady_state,
)

# cw_report's fit_dip transiently probes near-zero decay times while
# curve_fit converges on the stiff 300 K case below (i); the fitted A_dip/
# tau_dip are not what check (i) asserts on.
warnings.filterwarnings("ignore", category=RuntimeWarning, module="fsim_core.cw_g2")
from fsim_core.integrator import g2_from, retention
from fsim_core.loading import f1b_g2

CHECKS = []
RESULTS = {}   # numbers to print in the summary


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


GX, GXX = 1.0, 2.0          # 1/ns: tau_X = 1 ns, tau_XX = 0.5 ns (cascade ratio)


# ------------------------------------------------------------------ preliminaries

@check("steady state: closed form is the null vector of M, normalized")
def _():
    for r, kx, kxx, p in [(0.3, 0.0, 0.0, 1.0), (2.0, 0.4, 0.1, 0.7), (1e-3, 0.0, 0.0, 1.0)]:
        M = generator(r, GX, GXX, kx, kxx, p)
        P = steady_state(r, GX, GXX, kx, kxx, p)
        assert abs(P.sum() - 1.0) < 1e-14
        assert np.max(np.abs(M @ P)) < 1e-13, (r, kx, kxx, p)
        ns = null_space(M)
        assert ns.shape[1] == 1
        v = ns[:, 0] / ns[:, 0].sum()
        assert np.max(np.abs(v - P)) < 1e-12


@check("propagate: eigen path matches scipy expm at 30 taus")
def _():
    M = generator(0.8, GX, GXX, 0.2, 0.05, 1.0)
    taus = np.linspace(0.0, 8.0, 30)
    for v0 in (np.array([1.0, 0, 0]), np.array([0, 1.0, 0])):
        a = propagate(M, v0, taus, "eig")
        b = propagate(M, v0, taus, "expm")
        assert np.max(np.abs(a - b)) < 1e-12


@check("escape rates: S = gamma/(gamma+k) reproduces integrator.retention")
def _():
    for T in (77.0, 150.0, 250.0):
        kX, kXX = escape_rates_from_retention(GX, 1e5, 240.0, 100.0, 35.0, T, gamma_XX_ns=GXX)
        S = float(retention(T, 1e5, 240.0, 100.0, 35.0))
        assert abs(GX / (GX + kX) - S) < 1e-12
        assert abs(GXX / (GXX + kXX) - S) < 1e-12


# ------------------------------------------------------------------ (a) two-level

@check("(a) two-level limit (pump_ratio=0): g2 = 1 - exp(-(r+gamma_X+k_X)|tau|) to 1e-6")
def _():
    tau = np.linspace(-6.0, 6.0, 1201)
    for r, kx in [(0.5, 0.0), (2.0, 0.3), (0.05, 1.0)]:
        g = g2_cw(tau, r, GX, GXX, t_X=0.7, t_XX=0.2, k_X=kx, k_XX=0.1, pump_ratio=0.0)
        ref = g2_two_level(tau, r, GX, kx)
        assert np.max(np.abs(g - ref)) < 1e-6, (r, kx, np.max(np.abs(g - ref)))


# ------------------------------------------------------------- (b) stationarity

@check("(b) stationarity: g2 -> 1 within 1e-6 at tau = 50/min(rate); G -> I_det^2")
def _():
    for r, kx, kxx, txx in [(0.3, 0.0, 0.0, 0.3), (1.5, 0.2, 0.1, 0.05), (0.02, 0.0, 0.0, 0.5)]:
        rates = [r, GX + kx, GXX + kxx]
        tau_far = 50.0 / min(rates)
        g, parts = g2_cw([tau_far], r, GX, GXX, 1.0, txx, kx, kxx, return_parts=True)
        assert abs(g[0] - 1.0) < 1e-6, (r, g[0])
        assert abs(parts["G"][0] / parts["I_det"] ** 2 - 1.0) < 1e-6


# ---------------------------------------------------------------- (c) eps = 0

@check("(c) eps = 0 gives g2_dot(0) = 0 exactly; t_XX > 0 gives > 0, increasing in t_XX")
def _():
    r = 0.7
    g0 = g2_cw([0.0], r, GX, GXX, 1.0, 0.0)[0]
    assert g0 == 0.0, g0
    assert g2_cw_zero(r, GX, GXX, 0.0) == 0.0
    prev = 0.0
    for txx in np.linspace(0.02, 0.5, 13):
        g = g2_cw([0.0], r, GX, GXX, 1.0, txx)[0]
        assert g > prev, (txx, g, prev)
        prev = g
    # depends on the filter only through eps: scale both transmissions
    g1 = g2_cw([0.0, 0.3, 1.0], r, GX, GXX, 1.0, 0.3)
    g2 = g2_cw([0.0, 0.3, 1.0], r, GX, GXX, 0.4, 0.12)
    assert np.max(np.abs(g1 - g2)) < 1e-12


# ----------------------------------------------------------- (d) low-pump limit

@check("(d) low pump r = 1e-3 gamma_X: g2_dot(0) = eps p S_XX/S_X within 2% (and exact form to 1e-10)")
def _():
    eps = 0.3
    for kx, kxx, p in [(0.0, 0.0, 1.0), (0.5, 0.2, 1.0), (0.1, 0.1, 0.6)]:
        r = 1e-3 * GX
        num = g2_cw([0.0], r, GX, GXX, 1.0, eps, kx, kxx, p)[0]
        exact = g2_cw_zero(r, GX, GXX, eps, kx, kxx, p)
        low = g2_cw_zero_low_pump(GX, GXX, eps, kx, kxx, p)
        assert abs(num - exact) < 1e-10 * max(1.0, exact), (num, exact)
        assert abs(num / low - 1.0) < 0.02, (num, low)
        RESULTS.setdefault("d_low_pump", []).append((kx, kxx, p, num, low))


# ----------------------------------------------------- (e) cascade bunching signature

@check("(e) eps=0.3, r=gamma_X: CW g2_dot(0) exceeds pulsed f1b_g2(mu, eps) for all mu in [0.1, 1]")
def _():
    eps, r = 0.3, GX
    cw0 = g2_cw([0.0], r, GX, GXX, 1.0, eps)[0]
    cw0_closed = g2_cw_zero(r, GX, GXX, eps)
    assert abs(cw0 - cw0_closed) < 1e-12
    mus = np.linspace(0.1, 1.0, 10)
    pulsed = np.array([float(f1b_g2(m, eps)) for m in mus])
    RESULTS["e_cw0"] = cw0
    RESULTS["e_pulsed"] = list(zip(mus, pulsed))
    RESULTS["e_pulsed_sat"] = 2 * eps / (1 + eps) ** 2
    assert np.all(cw0 > pulsed), (cw0, pulsed)
    # the excess is a bunching feature of width ~ tau_X: g2 at |tau| ~ 0.2 ns
    # is still above the pulsed value; the curve is a dip + narrow peak
    tau = np.linspace(-5, 5, 2001)
    g = g2_cw(tau, r, GX, GXX, 1.0, eps)
    assert np.all(g[np.abs(tau) < 0.2] > pulsed.max())
    # far from tau=0 the curve relaxes to 1 from below or above -- no NaNs
    assert np.all(np.isfinite(g))


# --------------------------------------------------------------- (f) background

@check("(f) background law g2_meas(0) = 1 - rho^2 (1 - g2_dot(0)) to 1e-12 (and pointwise)")
def _():
    tau = np.linspace(-4, 4, 801)
    g = g2_cw(tau, 0.6, GX, GXX, 1.0, 0.2, 0.1, 0.05)
    for rho in (1.0, 0.94, 0.88, 0.5):
        gm = g2_with_background(g, rho)
        i0 = np.argmin(np.abs(tau))
        assert abs(gm[i0] - (1 - rho**2 * (1 - g[i0]))) < 1e-12
        assert abs(gm[i0] - g2_from(g[i0], rho)) < 1e-12
        assert np.max(np.abs(gm - (1 - rho**2 * (1 - g)))) < 1e-12


# --------------------------------------------------------------------- (g) IRF

@check("(g1) dip_convolved (Gaussian) vs numeric convolution within 1e-4 (A=0.9, tau_d=0.5 ns, fwhm=0.5 ns)")
def _():
    A, td, fwhm = 0.9, 0.5, 500.0
    tau = np.arange(-10.0, 10.0 + 1e-9, 0.001)
    g = 1 - A * np.exp(-np.abs(tau) / td)
    num = raw_g2_0(tau, g, fwhm, "gaussian")
    closed = dip_convolved(A, td, irf_width(fwhm, "gaussian"))
    RESULTS["g1"] = (num, closed)
    assert abs(num - closed) < 1e-4, (num, closed)
    # also the full curve: convolved curve equals the kernel-weighted average everywhere (spot check at tau=0.3)
    conv = convolve_irf(tau, g, fwhm, "gaussian")
    assert abs(np.interp(0.0, tau, conv) - num) < 1e-15


@check("(g1x) dip_convolved_exp (two-sided exponential IRF) vs numeric convolution within 1e-4")
def _():
    A, td, fwhm = 0.85, 1.0, 500.0
    tau = np.arange(-20.0, 20.0 + 1e-9, 0.001)
    g = 1 - A * np.exp(-np.abs(tau) / td)
    num = raw_g2_0(tau, g, fwhm, "exponential")
    closed = dip_convolved_exp(A, td, irf_width(fwhm, "exponential"))
    assert abs(num - closed) < 1e-4, (num, closed)


@check("(g2) deconvolve_dip(dip_convolved(A)) = A within 1e-9 (both IRF shapes)")
def _():
    for A, td, fwhm in [(0.9, 0.5, 500.0), (0.85, 1.0, 500.0), (0.3, 2.0, 120.0)]:
        raw = dip_convolved(A, td, irf_width(fwhm, "gaussian"))
        assert abs(deconvolve_dip(raw, td, fwhm, "gaussian") - A) < 1e-9
        assert abs(intrinsic_from_raw(raw, td, fwhm, "gaussian") - (1 - A)) < 1e-9
        raw_e = dip_convolved_exp(A, td, irf_width(fwhm, "exponential"))
        assert abs(deconvolve_dip(raw_e, td, fwhm, "exponential") - A) < 1e-9


@check("(g3) Reischle anchor: A=0.85 (g2=0.15), tau_d ~1 ns, 500 ps IRF -> raw g2(0) in [0.3, 0.5]")
def _():
    A, td = 0.85, 1.0
    # Reischle 2008 fitted with IRF C exp(-|tau|/0.5 ns): exponential, time constant 0.5 ns
    raw_exp = dip_convolved_exp(A, td, 0.5)
    # the same 0.5 ns quoted as a Gaussian FWHM, and as the exponential's FWHM (2 ln2 * 0.5 ns = 0.69 ns)
    raw_gauss_fwhm500 = dip_convolved(A, td, irf_width(500.0, "gaussian"))
    raw_gauss_fwhm693 = dip_convolved(A, td, irf_width(2 * np.log(2) * 500.0, "gaussian"))
    # implied recovery time from their own numbers (raw 0.41, deconvolved 0.15, exp IRF 0.5 ns)
    td_implied = 0.5 * (1 - 0.41) / (A - (1 - 0.41))
    RESULTS["g3"] = dict(raw_exp_tau05=raw_exp, raw_gauss_fwhm500=raw_gauss_fwhm500,
                         raw_gauss_fwhm693=raw_gauss_fwhm693, td_implied=td_implied)
    assert 0.3 <= raw_exp <= 0.5, raw_exp
    assert 0.3 <= raw_gauss_fwhm693 <= 0.5, raw_gauss_fwhm693
    # the same anchor through the full pipeline: cw_report with a rho-limited dot
    # (eps small, rho chosen so g2_meas(0) = 0.15), exponential IRF of 0.5 ns time constant
    eps, r = 0.03, 0.4
    g0 = g2_cw_zero(r, GX, GXX, eps)
    rho = np.sqrt((1 - 0.15) / (1 - g0))
    rep = cw_report(r, GX, GXX, 0.0, 0.0, 1.0, eps, rho, 2 * np.log(2) * 500.0,
                    tau_max_ns=15.0, irf_shape="exponential")
    RESULTS["g3_pipeline"] = dict(g2_meas0=rep["g2_meas0"], g2_raw0=rep["g2_raw0"],
                                  tau_dip=rep["tau_dip"], A_dip=rep["A_dip"],
                                  back=rep["g2_intrinsic_from_raw"])
    assert abs(rep["g2_meas0"] - 0.15) < 1e-9
    assert 0.3 <= rep["g2_raw0"] <= 0.5, rep["g2_raw0"]
    # deconvolving the raw number with the fitted single-exponential dip time
    # recovers the intrinsic value only up to the MODEL ERROR of the single-
    # exponential assumption (the true dip has two rates, r+gamma_X and the
    # XX branch): the bias is reported below and bounded here at 0.05.
    RESULTS["g3_pipeline"]["bias"] = rep["g2_intrinsic_from_raw"] - 0.15
    assert abs(rep["g2_intrinsic_from_raw"] - 0.15) < 0.05, rep["g2_intrinsic_from_raw"]


# --------------------------------------------------------------- (h) monotone

@check("(h) raw g2(0) increases monotonically with IRF fwhm (closed form and full pipeline)")
def _():
    A, td = 0.85, 1.0
    prev = -1.0
    for fwhm in np.linspace(50.0, 1500.0, 30):
        v = dip_convolved(A, td, irf_width(fwhm, "gaussian"))
        assert v > prev
        prev = v
    prev = -1.0
    for fwhm in (100.0, 250.0, 500.0, 800.0, 1200.0):
        rep = cw_report(0.5, GX, GXX, 0.1, 0.05, 1.0, 0.1, 0.9, fwhm, tau_max_ns=12.0)
        assert rep["g2_raw0"] > prev, (fwhm, rep["g2_raw0"], prev)
        prev = rep["g2_raw0"]


# --------------------------------------------------------- (i) stiff 300 K grid

_STIFF_PARAMS = dict(r_ns=10.0, gX=1.0, gXX=2.0, kX=1e3, kXX=1e3, t_X=1.0, t_XX=0.2, rho=0.95)


@check("(i1) stiff 300K-class escape rates (r_ns=10, k_X=k_XX=1e3/ns): "
       "cw_report completes in < 2 s")
def _():
    p = _STIFF_PARAMS
    t0 = time.perf_counter()
    rep = cw_report(p["r_ns"], p["gX"], p["gXX"], p["kX"], p["kXX"], p["t_X"], p["t_XX"],
                    p["rho"], irf_fwhm_ps=100.0, tau_max_ns=10.0)
    elapsed = time.perf_counter() - t0
    RESULTS["i_elapsed"] = elapsed
    RESULTS["i_g2_dot0"] = rep["g2_dot0"]
    assert elapsed < 2.0, elapsed


@check("(i2) stiff case's g2_cw0 matches the closed-form low-pump limit "
       "eps p S_XX/S_X within 5% (grid-bounding must not degrade propagate())")
def _():
    p = _STIFF_PARAMS
    eps = p["t_XX"] / p["t_X"]
    rep = cw_report(p["r_ns"], p["gX"], p["gXX"], p["kX"], p["kXX"], p["t_X"], p["t_XX"],
                    p["rho"], irf_fwhm_ps=100.0, tau_max_ns=10.0)
    low = g2_cw_zero_low_pump(p["gX"], p["gXX"], eps, p["kX"], p["kXX"], 1.0)
    RESULTS["i_low"] = low
    assert abs(rep["g2_dot0"] / low - 1.0) < 0.05, (rep["g2_dot0"], low)
    # the numeric curve (propagate + the tau grid), not just the closed-form
    # field, must agree too -- propagate() is exact per tau regardless of grid
    # density, so bounding the grid must not have degraded g2_dot itself.
    g_numeric = float(np.interp(0.0, rep["curves"]["tau"], rep["curves"]["g2_dot"]))
    assert abs(g_numeric / low - 1.0) < 0.05, (g_numeric, low)


# ------------------------------------------------------------ report + note sanity

@check("cw_report: fields present, curves symmetric, raw curve -> g2_meas edge, note derivation present")
def _():
    rep = cw_report(0.8, GX, GXX, 0.1, 0.05, 0.8, 0.16, 0.92, 400.0, tau_max_ns=10.0)
    for k in ("P_ss", "I_X", "I_XX", "g2_dot0", "g2_meas0", "g2_raw0", "tau_dip", "curves", "notes"):
        assert k in rep, k
    c = rep["curves"]
    assert np.max(np.abs(c["g2_dot"] - c["g2_dot"][::-1])) < 1e-12
    assert abs(c["g2_raw"][-1] - c["g2_meas"][-1]) < 1e-6
    assert abs(rep["g2_dot0"] - np.interp(0.0, c["tau"], c["g2_dot"])) < 1e-10
    assert abs(rep["eps"] - 0.2) < 1e-15
    note = cw_vs_pulsed_note()
    assert "(1 + a + a b)" in note and "eps p S_XX/S_X" in note


@check("finding 8a: f1b_g2 small-mu drive-factor slope f = 1 + mu(1/3-eps) + O(mu^2)")
def _():
    # f(mu, eps) := f1b_g2(mu, eps)/eps  [docs: cw_g2.py:65-68, 418-429, 454-462]
    # mu=1e-4 residual is 2.4e-10 against tol 1e-9 (margin ~4x, under the 10x
    # threshold that would call for tightening to mu=1e-5 per spec). Tried
    # that: it does NOT help. loading.loading_probs computes
    # P2 = 1 - P0 - P1 by cancelling two near-1 float64 terms, so float
    # noise in P2 (~1e-16 absolute) swamps the true O(mu^2) term once mu is
    # this small. Observed residual at mu=1e-5 is 1.1e-7 -- WORSE, not
    # better -- confirmed against mpmath at 50 dps, whose noise-free
    # residual at mu=1e-5 is -3.3e-13 (i.e. the 1.1e-7 is float64
    # cancellation noise in f1b_g2, not a sign the expansion is wrong).
    # mu=1e-4 / tol=1e-9 is kept: it is the point where the float64
    # implementation still resolves the true O(mu^2) term.
    mu, eps = 1e-4, 0.2
    resid = f1b_g2(mu, eps) / eps - (1 + mu * (1 / 3 - eps))
    assert abs(resid) < 1e-9, resid


@check("finding 8b: f1b_g2 large-mu limit at eps=1: f(50,1.0) -> 2/(1+eps)^2 = 0.5")
def _():
    assert abs(f1b_g2(50.0, 1.0) / 1.0 - 0.5) < 1e-9


@check("finding 8c: the 1/3 < eps < sqrt(2)-1 dip -- f(50,0.35) > 1 but f(1e-3,0.35) < 1")
def _():
    assert f1b_g2(50.0, 0.35) / 0.35 > 1
    assert f1b_g2(1e-3, 0.35) / 0.35 < 1


@check("finding 8d: eps > sqrt(2)-1 = 0.41421356 -- f(50,0.45) < 1 (large-mu limit 2/1.45^2 = 0.951)")
def _():
    assert f1b_g2(50.0, 0.45) / 0.45 < 1


@check("finding 8e: f1b_g2 large-mu limit at eps=0.2: f(1e4,0.2) -> 2/(1+eps)^2")
def _():
    assert abs(f1b_g2(1e4, 0.2) / 0.2 - 2 / 1.2 ** 2) < 1e-6


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
    print(f"\n{n - failed}/{n} checks passed")
    if "e_cw0" in RESULTS:
        print(f"\n(e) CW g2_dot(0) [eps=0.3, r=gamma_X, gamma_XX=2 gamma_X] = {RESULTS['e_cw0']:.4f}")
        print("    pulsed f1b_g2(mu, 0.3): " + ", ".join(f"mu={m:.1f}:{v:.4f}" for m, v in RESULTS["e_pulsed"]))
        print(f"    pulsed saturation 2 eps/(1+eps)^2 = {RESULTS['e_pulsed_sat']:.4f}")
    if "d_low_pump" in RESULTS:
        for kx, kxx, p, num, low in RESULTS["d_low_pump"]:
            print(f"(d) r=1e-3 gamma_X, k_X={kx}, k_XX={kxx}, p={p}: numeric {num:.5f}  limit eps p S_XX/S_X = {low:.5f}")
    if "g1" in RESULTS:
        print(f"(g1) Gaussian dip: numeric {RESULTS['g1'][0]:.6f}  closed {RESULTS['g1'][1]:.6f}")
    if "i_elapsed" in RESULTS:
        print(f"(i) stiff 300K-class case (r_ns=10, k_X=k_XX=1e3/ns): cw_report took "
              f"{RESULTS['i_elapsed']:.3f} s; g2_dot0 = {RESULTS['i_g2_dot0']:.4f} vs "
              f"low-pump limit {RESULTS['i_low']:.4f}")
    if "g3" in RESULTS:
        g = RESULTS["g3"]
        print(f"(g3) Reischle anchor A=0.85, tau_d=1 ns: raw g2(0) = {g['raw_exp_tau05']:.3f} "
              f"[exp IRF, tau_i=0.5 ns, their fit model]; {g['raw_gauss_fwhm693']:.3f} "
              f"[Gaussian, FWHM = 2 ln2 * 0.5 = 0.69 ns]; {g['raw_gauss_fwhm500']:.3f} "
              f"[Gaussian, FWHM 0.5 ns]; their 0.41 -> 0.15 implies tau_d = {g['td_implied']:.2f} ns")
        p = RESULTS["g3_pipeline"]
        print(f"     pipeline: g2_meas(0)={p['g2_meas0']:.3f} -> raw {p['g2_raw0']:.3f}; fitted dip "
              f"A={p['A_dip']:.3f}, tau_dip={p['tau_dip']:.3f} ns; deconvolved back to {p['back']:.3f} "
              f"(single-exponential model bias {p['bias']:+.3f})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
