"""Regression suite for fsim_core.pulse_counting (finite-pulse photon-
counting moment hierarchy, peer-review-triage.md finding 1).

Every check anchors on a value this module did NOT produce for itself:
(a)/(b) are published/hand numbers (the reproduction appendix of
../quantum-dot-peer-review-2026-09-06.md, reproduced exactly by the Opus
triage in .workers/review/peer-review-triage.md); (c) is an independent
scipy.integrate.solve_ivp integration of the same 9-vector; (d) is the
closed-form Lemma 1 (collection-efficiency invariance) stated in
pulse_counting.py's own docstring.

Check (b) uses EXPLICIT, hardcoded rates rather than DeviceDesign.load() +
evaluate(): pr-pkg1-capture-escape changed the confinement retention
prefactor (card density 3e8 cm^-2 instead of dot_levels' 1e10 default), so
k_X/k_XX (and, empirically, mu_resolved/r_dot too -- other fixes have moved
the card's operating point since the review) at the gainp favourable corner
are no longer the PRE-package-1 values the triage's reproduction reported.
mu_resolved, static g2 and k_X are taken verbatim from that reproduction;
eps is recovered from them by inverting loading.f1b_g2(mu, eps) = static_g2
(the closed form is monotonic in eps on (0, 1), so the root is unique) --
this reconstructs the exact (r_ns, gamma_X_ns, eps, k_X, k_XX) the triage's
own script used, without going back through the (now-drifted) card.

Run: python verify/verify_pulse_counting.py   (exit code 0 iff all pass)
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import cw_g2
from fsim_core.loading import f1b_g2
from fsim_core.pulse_counting import _augmented, _periodic_steady_state, pulse_g2

CHECKS = []
RESULTS = {}


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# ------------------------------------------------------- (a) instantaneous-pulse limit

@check("(a) no-escape, tau_on -> 1e-6 ns limit reproduces loading.f1b_g2(1.0, 0.1) "
       "to rtol 1e-4 (the static per-pulse loading distribution is exact only in "
       "this instantaneous-pulse limit)")
def _():
    mu, eps = 1.0, 0.1
    tau_on_ns = 1e-6
    r_ns = mu / tau_on_ns
    # tau_dark long enough (200 X-lifetimes at gamma_X_ns=1.0) that the
    # periodic steady state returns to the ground state, matching F1b's
    # "isolated instantaneous pulse from an empty dot" assumption.
    tau_dark_ns = 200.0
    result = pulse_g2(r_ns, gamma_X_ns=1.0, gamma_XX_ns=2.0, k_X=0.0, k_XX=0.0,
                      t_X=1.0, t_XX=eps, tau_on_ns=tau_on_ns, tau_dark_ns=tau_dark_ns)
    target = f1b_g2(mu, eps)
    RESULTS["a"] = (result["g2"], target)
    assert result["converged"]
    assert abs(result["g2"] - target) / target < 1e-4, (result["g2"], target)


# ------------------------------------------------- (b) gainp favourable corner (triage)

# Reproduction of the finite-pulse appendix, ../quantum-dot-peer-review-2026-
# 09-06.md "Reproduction appendix", verbatim from .workers/review/peer-
# review-triage.md's "Reproduction of the finite-pulse appendix" table
# (gamma300=6.0, delta_xx=8.0, NA=0.8, R_back=0.95, L=250 um; PRE-package-1
# escape rates -- see the module docstring above). gamma_X_ns=1.0,
# gamma_XX_ns=2.0 (cw_g2.escape_rates_from_retention's own gamma_XX_ns=
# 2*gamma_X_ns default, exactly how the appendix's own
# `generator(rate, 1.0, 2.0, kx, kxx)` call used it); k_XX = 2*k_X for the
# same reason (both k_X, k_XX share the ONE Arrhenius escape-attempt factor
# esc(T), scaled by gamma_X_ns=1.0 / gamma_XX_ns=2.0 respectively).
_TAU_ON_NS = 0.1
_TAU_DARK_NS = 12.4
_GAINP_CASES = [
    dict(T_K=230.0, mu=1.0070875689, static_g2=0.0588386352,
        k_X=75.2708, k_XX=2.0 * 75.2708, dynamic_g2=0.7951419634),
    dict(T_K=300.0, mu=0.9999999987, static_g2=0.1612177631,
        k_X=215.5697, k_XX=2.0 * 215.5697, dynamic_g2=0.9255729538),
]


@check("(b) gainp favourable corner: finite-pulse dot g2 matches the "
       "reproduction appendix to +/- 1e-6 at 230 K and 300 K (explicit "
       "PRE-package-1 rates, not the live card)")
def _():
    rows = []
    for c in _GAINP_CASES:
        eps = brentq(lambda e: f1b_g2(c["mu"], e) - c["static_g2"], 1e-9, 1.0 - 1e-9)
        r_ns = c["mu"] / _TAU_ON_NS
        result = pulse_g2(r_ns, gamma_X_ns=1.0, gamma_XX_ns=2.0, k_X=c["k_X"], k_XX=c["k_XX"],
                          t_X=1.0, t_XX=eps, tau_on_ns=_TAU_ON_NS, tau_dark_ns=_TAU_DARK_NS)
        rows.append((c["T_K"], eps, result["g2"], result["mean_counts"], c["dynamic_g2"]))
        assert result["converged"]
        assert abs(result["g2"] - c["dynamic_g2"]) < 1e-6, (c["T_K"], result["g2"], c["dynamic_g2"])
    RESULTS["b"] = rows


# ------------------------------------------------------ (c) independent integrator

@check("(c) an independent scipy.integrate.solve_ivp (RK45, rtol 1e-11, "
       "atol 1e-14) integration of the same 9-vector agrees with expm to "
       "1e-9 in g2 at the 230 K gainp parameter set")
def _():
    c = _GAINP_CASES[0]
    eps = brentq(lambda e: f1b_g2(c["mu"], e) - c["static_g2"], 1e-9, 1.0 - 1e-9)
    r_ns = c["mu"] / _TAU_ON_NS
    gamma_X_ns, gamma_XX_ns = 1.0, 2.0
    M_on = cw_g2.generator(r_ns, gamma_X_ns, gamma_XX_ns, c["k_X"], c["k_XX"])
    M_off = cw_g2.generator(0.0, gamma_X_ns, gamma_XX_ns, c["k_X"], c["k_XX"])
    J = np.zeros((3, 3))
    J[0, 1] = 1.0 * gamma_X_ns
    J[1, 2] = eps * gamma_XX_ns
    A_on, A_off = _augmented(M_on, J), _augmented(M_off, J)
    p_ss, converged = _periodic_steady_state(M_on, M_off, _TAU_ON_NS, _TAU_DARK_NS)
    assert converged
    v0 = np.zeros(9)
    v0[0:3] = p_ss

    def rhs(_t, v, A):
        return A @ v

    v1 = solve_ivp(rhs, [0.0, _TAU_ON_NS], v0, args=(A_on,),
                   method="RK45", rtol=1e-11, atol=1e-14).y[:, -1]
    v2 = solve_ivp(rhs, [0.0, _TAU_DARK_NS], v1, args=(A_off,),
                   method="RK45", rtol=1e-11, atol=1e-14).y[:, -1]
    g2_ivp = float(v2[6:9].sum() / v2[3:6].sum() ** 2)

    v_expm = expm(A_on * _TAU_ON_NS) @ v0
    v_expm = expm(A_off * _TAU_DARK_NS) @ v_expm
    g2_expm = float(v_expm[6:9].sum() / v_expm[3:6].sum() ** 2)

    RESULTS["c"] = (g2_ivp, g2_expm)
    assert abs(g2_ivp - g2_expm) < 1e-9, (g2_ivp, g2_expm)


# ------------------------------------------------------------- (d) Lemma 1

@check("(d) g2 is invariant under scaling t_X and t_XX by a common factor "
       "(collection efficiency cancels, Lemma 1) to 1e-12")
def _():
    r_ns, gamma_X_ns, gamma_XX_ns = 5.4, 1.0, 2.0
    k_X, k_XX = 3.0, 6.0
    t_X, t_XX = 0.5, 0.05
    g2_ref = pulse_g2(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, t_X, t_XX,
                      _TAU_ON_NS, _TAU_DARK_NS)["g2"]
    for factor in (0.37, 2.0, 1e-3, 10.0):
        g2_scaled = pulse_g2(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX,
                             t_X * factor, t_XX * factor, _TAU_ON_NS, _TAU_DARK_NS)["g2"]
        assert abs(g2_scaled - g2_ref) < 1e-12, (factor, g2_scaled, g2_ref)


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
    print(f"\n{n - failed}/{n} pulse_counting checks passed")
    if "a" in RESULTS:
        got, target = RESULTS["a"]
        print(f"(a) instantaneous-pulse limit: pulse_g2 = {got:.8f}  f1b_g2(1.0, 0.1) = {target:.8f}")
    if "b" in RESULTS:
        for T_K, eps, g2, mean_counts, target in RESULTS["b"]:
            print(f"(b) {T_K:.0f} K: eps={eps:.6f}  finite-pulse dot g2 = {g2:.10f} "
                  f"(reference {target:.10f})  mean_counts = {mean_counts:.6e}")
    if "c" in RESULTS:
        g2_ivp, g2_expm = RESULTS["c"]
        print(f"(c) solve_ivp g2 = {g2_ivp:.12f}  expm g2 = {g2_expm:.12f}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
