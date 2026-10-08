"""Regression suite for fsim_core.pulse_counting (finite-pulse photon-
counting moment hierarchy).

Every check anchors on a value this module did NOT produce for itself:
(a)/(b) are published/hand numbers (the reproduction appendix, reproduced
exactly by peer review); (c) is an independent scipy.integrate.solve_ivp
integration of the same 9-vector; (d) is the closed-form Lemma 1
(collection-efficiency invariance) stated in pulse_counting.py's own
docstring.

Audit D2/D7 additions (audit-pulse-g2-label). (f) the adjacent-peak-
normalised g2_adj = <m(m-1)>/<m_n m_(n+1)> must equal the long-delay-
normalised g2 = <m(m-1)>/<m>^2 when the period is long enough that nothing
carries over (identity: m_n, m_(n+1) independent => <m_n m_(n+1)> = <m>^2).
(g) with carry-over, g2_adj must match an independent event-by-event
(Gillespie) Monte Carlo of the same incoherent three-state chain, which
estimates <m(m-1)> and <m_n m_(n+1)> directly from simulated consecutive
periods -- a number the moment propagation did not produce. (h) the
two-level incoherent rectangular-pump limits derived in the audit (C9,
C13): strong fast pump g2 -> (2 G tau + (G tau)^2)/(1 + G tau)^2 (m = 1 + K,
K ~ Poisson(G tau), i.e. ~ 2 G tau), weak pump g2 -> G tau/3 [DR].

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
from fsim_core.drive_mech import reexc_g2
from fsim_core.loading import f1b_g2
from fsim_core.pulse_counting import (_augmented, _periodic_steady_state,
                                      deterministic_cycle_g2, pulse_g2)

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

# Reproduction of the finite-pulse appendix from ../quantum-dot-peer-review-
# 2026-09-06.md "Reproduction appendix"
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


# ------------------------------------------- (e) cross-check against drive_mech.reexc_g2

# pulse_counting.pulse_g2 and drive_mech.reexc_g2 solve the SAME closed
# counting-moment hierarchy (module docstrings cross-reference each other),
# but reexc_g2 has no escape channels (k_X=k_XX=0 always), counts the X
# line only (t_XX=0 here disables XX counting to match), and has no gate
# (a single pulse from mu0=0, decayed out over a generous tail, rather than
# a periodic steady state read out through a gate). To make pulse_g2's
# periodic construction reproduce that single-shot picture, tau_dark_ns is
# set to the SAME generous multiple of the X lifetime as reexc_g2's own
# t_end_factor tail (_REEXC_TAIL_FACTOR lifetimes -- residual leak into the
# next period is exp(-40) ~ 4e-18, far below the 1e-6 tolerance), so the
# periodic steady state is (to float precision) the ground state reexc_g2
# starts from, and pulse_g2's whole-period (gate_ns=None) counting window
# matches reexc_g2's full decay-tail integration window.
_REEXC_CASES = [
    dict(r_ns=0.5, tau_on_ns=0.3, gamma_X_ns=1.0),
    dict(r_ns=3.0, tau_on_ns=1.0, gamma_X_ns=1.0),
    dict(r_ns=1.0, tau_on_ns=2.0, gamma_X_ns=0.5),
]
_REEXC_TAIL_FACTOR = 40.0


@check("(e) cross-check against drive_mech.reexc_g2 (independent counting-moment "
      "hierarchy, no escape/XX channel/gate) at three parameter sets, k=0/t_XX=0, "
      "tolerance 1e-6")
def _():
    rows = []
    for c in _REEXC_CASES:
        gamma_X_ns = c["gamma_X_ns"]
        gamma_XX_ns = 2.0 * gamma_X_ns
        tau_x_ps = 1.0 / gamma_X_ns
        tail_ns = _REEXC_TAIL_FACTOR / gamma_X_ns
        pc = pulse_g2(c["r_ns"], gamma_X_ns, gamma_XX_ns, k_X=0.0, k_XX=0.0,
                      t_X=1.0, t_XX=0.0, tau_on_ns=c["tau_on_ns"], tau_dark_ns=tail_ns)
        rg = reexc_g2(c["r_ns"], c["tau_on_ns"], tau_x_ps, t_end_factor=_REEXC_TAIL_FACTOR)
        rows.append((c, pc["g2"], rg["g2"]))
        assert pc["converged"]
        assert abs(pc["g2"] - rg["g2"]) < 1e-6, (c, pc["g2"], rg["g2"])
    RESULTS["e"] = rows


# ------------------------------------ (f) adjacent-peak g2 == g2 without carry-over

@check("(f) no inter-period carry-over (dark window >> 1/gamma): g2_adj = "
       "<m(m-1)>/<m_n m_(n+1)> equals the long-delay-normalised g2 to 1e-6, "
       "pulse_g2 and deterministic_cycle_g2 (identity <m_n m_(n+1)> = <m>^2)")
def _():
    rows = []
    # (a)'s instantaneous-pulse case (200 X lifetimes of dark time).
    r = pulse_g2(1.0 / 1e-6, 1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 1e-6, 200.0, adjacent=True)
    rows.append(("pulse_g2 (a)-case", r["g2"], r["g2_adj"]))
    # gainp 230 K corner of (b) (k_X = 75/ns: carry-over ~ exp(-76*12.4)).
    c = _GAINP_CASES[0]
    eps = brentq(lambda e: f1b_g2(c["mu"], e) - c["static_g2"], 1e-9, 1.0 - 1e-9)
    r = pulse_g2(c["mu"] / _TAU_ON_NS, 1.0, 2.0, c["k_X"], c["k_XX"], 1.0, eps,
                 _TAU_ON_NS, _TAU_DARK_NS, adjacent=True)
    rows.append(("pulse_g2 gainp 230 K", r["g2"], r["g2_adj"]))
    # (g)'s finite-pulse rates with a 200 ns dark window (nontrivial g2).
    r = pulse_g2(2.0, 1.0, 2.0, 0.3, 0.6, 0.8, 0.3, 1.0, 200.0, adjacent=True)
    rows.append(("pulse_g2 (g)-rates, 200 ns dark", r["g2"], r["g2_adj"]))
    factors = [r["adjacent_peak_factor"]]
    # deterministic cycle, period = 200 X lifetimes.
    r = deterministic_cycle_g2(0.306, 0.612, 0.1, 0.2, 1.0, 0.3, 200.0 / 0.306,
                               adjacent=True)
    # (one pair into an empty dot gives g2 = 0 exactly; the factor is the test)
    rows.append(("deterministic_cycle_g2 long period", r["g2"], r["g2_adj"]))
    factors.append(r["adjacent_peak_factor"])
    RESULTS["f"] = rows
    for fac in factors:
        assert abs(fac - 1.0) < 1e-6, factors
    for name, g2, g2_adj in rows:
        assert abs(g2_adj - g2) < 1e-6, (name, g2, g2_adj)
    # the adjacent option must not perturb the default outputs (bit-identical)
    base = pulse_g2(2.0, 1.0, 2.0, 0.3, 0.6, 0.8, 0.3, 1.0, 0.5)
    ext = pulse_g2(2.0, 1.0, 2.0, 0.3, 0.6, 0.8, 0.3, 1.0, 0.5, adjacent=True)
    assert base["g2"] == ext["g2"] and base["mean_counts"] == ext["mean_counts"]


# ------------------------------- (g) adjacent-peak g2 vs independent Monte Carlo

def _mc_counts(rates_on, rates_off, t_X, t_XX, tau_on, tau_dark, load, n_chains,
               n_periods, burn, rng):
    """Vectorised Gillespie simulation of the incoherent (0, X, XX) chain.

    rates_* = (r, r2, gX, kX, gXX, kXX); a window with tau = 0 is skipped;
    load=True applies the deterministic one-pair load s -> min(s+1, 2) at the
    start of each period. Radiative X (XX) jumps are detected with
    probability t_X (t_XX). Returns counts[n_chains, n_periods] (after burn)."""
    s = np.zeros(n_chains, dtype=np.int64)
    out = np.zeros((n_chains, n_periods), dtype=np.int64)

    def window(rates, T, m):
        r, r2, gX, kX, gXX, kXX = rates
        t = np.zeros(n_chains)
        alive = np.ones(n_chains, dtype=bool)
        while alive.any():
            lam = np.where(s == 0, r, np.where(s == 1, r2 + gX + kX, gXX + kXX))
            with np.errstate(divide="ignore"):
                dt = rng.exponential(1.0, n_chains) / lam
            t = t + np.where(alive, dt, 0.0)
            alive &= (lam > 0) & (t < T)
            u = rng.random(n_chains) * lam
            det = rng.random(n_chains)
            up0 = alive & (s == 0)
            upX = alive & (s == 1) & (u < r2)
            radX = alive & (s == 1) & (u >= r2) & (u < r2 + gX)
            nrX = alive & (s == 1) & (u >= r2 + gX)
            radXX = alive & (s == 2) & (u < gXX)
            nrXX = alive & (s == 2) & (u >= gXX)
            m += (radX & (det < t_X)) + (radXX & (det < t_XX))
            s[up0] = 1
            s[upX] = 2
            s[radX | nrX] = 0
            s[radXX | nrXX] = 1

    for n in range(burn + n_periods):
        m = np.zeros(n_chains, dtype=np.int64)
        if load:
            np.minimum(s + 1, 2, out=s)
        if tau_on > 0:
            window(rates_on, tau_on, m)
        if tau_dark > 0:
            window(rates_off, tau_dark, m)
        if n >= burn:
            out[:, n - burn] = m
    return out


def _mc_g2_adj(counts, n_groups=40):
    """Ratio estimator <m(m-1)>/<m_n m_(n+1)> and <m(m-1)>/<m>^2 with a
    standard error from n_groups independent groups of chains."""
    a = counts.astype(float)
    fact = (a * (a - 1.0)).mean(axis=1)
    adj = (a[:, :-1] * a[:, 1:]).mean(axis=1)
    mean = a.mean(axis=1)
    groups = np.array_split(np.arange(counts.shape[0]), n_groups)
    ga = [fact[g].mean() / adj[g].mean() for g in groups]
    gl = [fact[g].mean() / mean[g].mean() ** 2 for g in groups]
    return (fact.mean() / adj.mean(), np.std(ga, ddof=1) / np.sqrt(n_groups),
            fact.mean() / mean.mean() ** 2, np.std(gl, ddof=1) / np.sqrt(n_groups))


@check("(g) with inter-period carry-over, g2_adj agrees with an independent "
       "Gillespie Monte Carlo of <m(m-1)>/<m_n m_(n+1)> within 4 sigma "
       "(finite pulse: r=2, gX=1, gXX=2, kX=0.3, kXX=0.6, tX=0.8, tXX=0.3, "
       "on 1 ns / dark 0.5 ns; deterministic cycle: gX=0.306, kX=0.1, "
       "tX=1, tXX=1, period 2 ns), and the MC resolves g2_adj != g2")
def _():
    rng = np.random.default_rng(20260923)
    rows = []
    # finite rectangular pulse
    pr = dict(r=2.0, gX=1.0, gXX=2.0, kX=0.3, kXX=0.6, tX=0.8, tXX=0.3, on=1.0, dark=0.5)
    ex = pulse_g2(pr["r"], pr["gX"], pr["gXX"], pr["kX"], pr["kXX"], pr["tX"], pr["tXX"],
                  pr["on"], pr["dark"], adjacent=True)
    counts = _mc_counts((pr["r"], pr["r"], pr["gX"], pr["kX"], pr["gXX"], pr["kXX"]),
                        (0.0, 0.0, pr["gX"], pr["kX"], pr["gXX"], pr["kXX"]),
                        pr["tX"], pr["tXX"], pr["on"], pr["dark"], False,
                        4000, 150, 20, rng)
    rows.append(("pulse_g2",) + (ex["g2_adj"], ex["g2"], ex["adjacent_peak_factor"])
                + _mc_g2_adj(counts))
    # deterministic one-pair cycle
    dc = dict(gX=0.306, gXX=0.612, kX=0.1, kXX=0.2, tX=1.0, tXX=1.0, T=2.0)
    ex = deterministic_cycle_g2(dc["gX"], dc["gXX"], dc["kX"], dc["kXX"], dc["tX"],
                                dc["tXX"], dc["T"], adjacent=True)
    counts = _mc_counts(None, (0.0, 0.0, dc["gX"], dc["kX"], dc["gXX"], dc["kXX"]),
                        dc["tX"], dc["tXX"], 0.0, dc["T"], True, 4000, 150, 20, rng)
    rows.append(("deterministic_cycle_g2",) + (ex["g2_adj"], ex["g2"],
                                               ex["adjacent_peak_factor"])
                + _mc_g2_adj(counts))
    RESULTS["g"] = rows
    for name, g2_adj, g2, _f, mc_adj, se_adj, mc_g2, se_g2 in rows:
        assert abs(mc_adj - g2_adj) < 4.0 * se_adj, (name, g2_adj, mc_adj, se_adj)
        assert abs(mc_g2 - g2) < 4.0 * se_g2, (name, g2, mc_g2, se_g2)
        # carry-over is large enough that the two normalisations are
        # statistically distinguishable (the check is not vacuous)
        assert abs(g2_adj - g2) > 4.0 * se_adj, (name, g2_adj, g2, se_adj)


# ------------------------------------- (h) incoherent two-level rectangular limits

@check("(h) two-level incoherent rectangular pump (pump_ratio=0, t_XX=0, k=0): "
       "strong fast pump g2 -> (2 G tau + (G tau)^2)/(1 + G tau)^2 to rtol 2e-3, "
       "weak short pump g2 -> G tau/3 to rtol 1e-2 (audit C9/C13, D7)")
def _():
    G, tau = 1.0, 0.1
    strong = pulse_g2(1e4, G, 2.0 * G, 0.0, 0.0, 1.0, 0.0, tau, 200.0,
                      pump_ratio=0.0)["g2"]
    lam = G * tau
    strong_ref = (2.0 * lam + lam ** 2) / (1.0 + lam) ** 2
    tau_w = 1e-3
    weak = pulse_g2(1e-3, G, 2.0 * G, 0.0, 0.0, 1.0, 0.0, tau_w, 200.0,
                    pump_ratio=0.0)["g2"]
    weak_ref = G * tau_w / 3.0
    RESULTS["h"] = (strong, strong_ref, weak, weak_ref)
    assert abs(strong - strong_ref) / strong_ref < 2e-3, (strong, strong_ref)
    assert abs(weak - weak_ref) / weak_ref < 1e-2, (weak, weak_ref)


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
    if "e" in RESULTS:
        for c, g2_pc, g2_rg in RESULTS["e"]:
            print(f"(e) r_ns={c['r_ns']:g} tau_on_ns={c['tau_on_ns']:g} "
                  f"gamma_X_ns={c['gamma_X_ns']:g}: pulse_g2 = {g2_pc:.10f}  "
                  f"reexc_g2 = {g2_rg:.10f}")
    if "f" in RESULTS:
        for name, g2, g2_adj in RESULTS["f"]:
            print(f"(f) {name}: g2 = {g2:.10f}  g2_adj = {g2_adj:.10f}")
    if "g" in RESULTS:
        for name, g2_adj, g2, fac, mc_adj, se_adj, mc_g2, se_g2 in RESULTS["g"]:
            print(f"(g) {name}: g2_adj = {g2_adj:.6f} (MC {mc_adj:.6f} +/- {se_adj:.6f}); "
                  f"g2 = {g2:.6f} (MC {mc_g2:.6f} +/- {se_g2:.6f}); "
                  f"adjacent_peak_factor = {fac:.6f}")
    if "h" in RESULTS:
        strong, strong_ref, weak, weak_ref = RESULTS["h"]
        print(f"(h) strong: g2 = {strong:.6f} (ref {strong_ref:.6f}); "
              f"weak: g2 = {weak:.6e} (ref G tau/3 = {weak_ref:.6e})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
