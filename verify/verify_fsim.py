"""FSIM regression suite (Yuki standard: units, limits, second methods).

Every closed form in Module C/E is checked against an independent method:
analytic limits, numeric quadrature, and Monte-Carlo sampling.
Run: python verify/verify_fsim.py   (exit code 0 iff all pass)
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.integrator import g2_from, g2_of_T, retention, rho_of_T, solve_Tc
from fsim_core.spectral import (
    KB,
    cavity_transmission,
    combined_transmission,
    epsilon,
    epsilon_narrow_filter,
    gamma_of_T,
    lorentzian,
    mc_transmission,
    optimal_width,
    tophat_transmission,
)

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


BASE_P = dict(delta_xx=5.9, gamma0=0.5, a_ac=2e-3, b_lo=25.0, E_lo=36.6,
              w=4.0, dx=-1.0, a_esc=1e5, E_a=240.0, b_p=100.0, E_b=35.0,
              b0=0.005, beta=0.5)


# ---------------------------------------------------------------- top-hat closed form

@check("tophat: wide-window limit t->1")
def _():
    assert abs(tophat_transmission(0.0, 1.0, 1e6) - 1.0) < 1e-5


@check("tophat: narrow-window limit t ~ w * L(0)")
def _():
    w = 1e-6
    t = tophat_transmission(0.3, 1.0, w)
    assert abs(t - w * lorentzian(0.0, 0.3, 1.0)) < 1e-9


@check("tophat: detuning symmetry t(+d)=t(-d)")
def _():
    assert abs(tophat_transmission(1.7, 0.9, 2.3) - tophat_transmission(-1.7, 0.9, 2.3)) < 1e-14


@check("tophat: closed form vs numeric quadrature (20 random cases)")
def _():
    rng = np.random.default_rng(1)
    for _ in range(20):
        d, g, w = rng.uniform(0, 10), rng.uniform(0.05, 8), rng.uniform(0.1, 20)
        num, _err = quad(lorentzian, -w / 2, w / 2, args=(d, g), limit=200)
        assert abs(tophat_transmission(d, g, w) - num) < 1e-7, (d, g, w)


@check("tophat: MC second method (1e6 samples, 4 sigma)")
def _():
    for d, g, w in [(0.0, 1.0, 2.0), (5.9, 6.0, 2.0), (3.0, 0.5, 1.0)]:
        t, s = mc_transmission(d, g, w=w)
        assert abs(t - tophat_transmission(d, g, w)) < 4 * s + 1e-4, (d, g, w)


# ----------------------------------------------------------------- cavity closed form

@check("cavity: on-resonance value kappa/(gamma+kappa)")
def _():
    g, k = 1.3, 0.7
    assert abs(cavity_transmission(0.0, g, k) - k / (g + k)) < 1e-14


@check("cavity: closed form vs numeric quadrature (20 random cases)")
def _():
    rng = np.random.default_rng(2)
    for _ in range(20):
        d, g, k = rng.uniform(0, 10), rng.uniform(0.05, 8), rng.uniform(0.05, 8)
        num, _err = quad(
            lambda x: lorentzian(x, d, g) * (k / 2) ** 2 / (x**2 + (k / 2) ** 2),
            -np.inf, np.inf, limit=400,
        )
        assert abs(cavity_transmission(d, g, k) - num) < 1e-6, (d, g, k)


@check("cavity: MC second method")
def _():
    for d, g, k in [(0.0, 1.0, 1.0), (5.0, 2.0, 0.5)]:
        t, s = mc_transmission(d, g, kappa=k)
        assert abs(t - cavity_transmission(d, g, k)) < 4 * s + 1e-4


@check("combined tophat*cavity: numeric vs MC")
def _():
    d, g, w, k = 2.0, 1.0, 3.0, 1.5
    t, s = mc_transmission(d, g, w=w, kappa=k)
    assert abs(t - combined_transmission(d, g, w, k)) < 4 * s + 1e-4


# ---------------------------------------------------------------------------- epsilon

@check("epsilon: F1 identity eps=1 at zero splitting, equal widths (any dx)")
def _():
    assert abs(epsilon(0.0, 1.0, 1.0, w=2.0, dx=-0.8).eps - 1.0) < 1e-14


@check("epsilon: mirror symmetry (dx, delta) == (-dx, -delta)")
def _():
    a = epsilon(5.9, 1.0, 1.5, w=3.0, dx=-1.2)
    b = epsilon(-5.9, 1.0, 1.5, w=3.0, dx=1.2)
    assert abs(a.eps - b.eps) < 1e-12 and abs(a.t_x - b.t_x) < 1e-12


@check("epsilon: blue-shifted window reduces red-XX leakage")
def _():
    e0 = epsilon(5.9, 2.0, 2.0, w=4.0, dx=0.0).eps
    e1 = epsilon(5.9, 2.0, 2.0, w=4.0, dx=-1.5).eps
    assert e1 < e0


@check("epsilon: narrow-filter bound is the w->0 limit (dx=0)")
def _():
    d, gx, gxx = 4.0, 1.0, 1.5
    r = epsilon(d, gx, gxx, w=1e-7)
    assert abs(r.eps - epsilon_narrow_filter(d, gx, gxx)) < 1e-6


@check("epsilon: monotone non-decreasing in w for centered filter")
def _():
    d, g = 5.9, 2.0
    ws = np.linspace(0.05, 30, 200)
    e = [epsilon(d, g, g, w=w).eps for w in ws]
    assert np.all(np.diff(e) > -1e-12)


@check("optimal_width: eps at w_opt equals the target")
def _():
    d, g = 5.9, 2.0
    tgt = 0.05
    w = optimal_width(d, g, g, tgt)
    assert abs(epsilon(d, g, g, w=w).eps - tgt) < 1e-9


# ------------------------------------------------------------------------ Module E

@check("Gamma(T): T->0 recovers Gamma0; monotone increasing")
def _():
    g = gamma_of_T(np.array([1e-3, 50, 150, 300]), 0.5, 1e-3, 20.0, 36.6)
    assert abs(g[0] - (0.5 + 1e-3 * 1e-3)) < 1e-12  # LO term frozen out at T->0
    assert np.all(np.diff(g) > 0)


@check("retention: frozen limit S->1; monotone decreasing in T")
def _():
    S = retention(np.array([1.0, 78, 150, 230, 300]), 1e5, 240.0, 100.0, 35.0)
    assert abs(S[0] - 1.0) < 1e-12 and np.all(np.diff(S) < 0)


@check("rho: frozen limit 1/(1+b0); high-T below it; beta=0 recovers S/(S+b0)")
def _():
    b0 = 0.02
    assert abs(rho_of_T(1.0, 1e5, 240.0, 100.0, 35.0, b0, 0.5) - 1 / (1 + b0)) < 1e-10
    assert rho_of_T(300.0, 1e5, 240.0, 100.0, 35.0, b0, 0.5) < 1 / (1 + b0)
    S = retention(200.0, 1e5, 240.0, 100.0, 35.0)
    assert abs(rho_of_T(200.0, 1e5, 240.0, 100.0, 35.0, b0, 0.0) - S / (S + b0)) < 1e-12


@check("background law identities: rho=1 -> g2=eps; eps=1 -> g2=1; eps=0 -> 1-rho^2")
def _():
    assert abs(g2_from(0.3, 1.0) - 0.3) < 1e-14
    assert abs(g2_from(1.0, 0.7) - 1.0) < 1e-14
    assert abs(g2_from(0.0, 0.9) - (1 - 0.81)) < 1e-14


@check("master ceiling: g2(Tc) = 0.5 at the solved Tc (incl. callable windows)")
def _():
    p = dict(BASE_P)
    p["w"] = lambda T: np.interp(T, [78, 230], [1.4, 7.66])
    p["dx"] = lambda T: np.interp(T, [78, 230], [-0.37, -3.41])
    Tc = solve_Tc(p)
    assert np.isfinite(Tc) and abs(g2_of_T(Tc, p).g2 - 0.5) < 1e-6


@check("fit round-trip: synthetic g2(T) from known params recovered within +-0.03")
def _():
    from fsim_core.card import Card, DataSet, Param, Tag
    from fsim_core.fitting import fit_phase0

    truth = dict(BASE_P)
    rows = []
    for T, w, dx in [(78, 1.40, -0.37), (120, 3.19, -0.70), (150, 4.12, -1.27),
                     (170, 4.74, -1.65), (210, 6.69, -2.82), (230, 7.66, -3.41)]:
        p = {**truth, "w": w, "dx": dx}
        rows.append({"T": float(T), "g2": g2_of_T(T, p).g2, "err": 0.02,
                     "w": w, "dx": dx})

    def P(name, **kw):
        kw.setdefault("unit", "meV"); kw.setdefault("tag", Tag.A); kw.setdefault("source", "synthetic")
        return Param(name=name, **kw)

    g_truth_250 = float(gamma_of_T(250.0, truth["gamma0"], truth["a_ac"],
                                   truth["b_lo"], truth["E_lo"]))
    card = Card(
        name="synthetic", meta={"name": "synthetic"},
        params={
            "delta_xx": P("delta_xx", value=truth["delta_xx"]),
            "E_lo": P("E_lo", value=truth["E_lo"]),
            "E_b": P("E_b", value=truth["E_b"]),
            "gamma0": P("gamma0", range=(0.05, 2.0)),
            "a_ac": P("a_ac", range=(0.0, 0.02), unit="meV/K"),
            "b_lo": P("b_lo", range=(1.0, 80.0)),
            "a_esc": P("a_esc", range=(1e3, 1e10), unit="1"),
            "E_a": P("E_a", range=(200.0, 290.0)),
            "b_p": P("b_p", range=(1.0, 1e4), unit="1"),
            "b0": P("b0", range=(1e-5, 0.1), unit="1"),
            "beta": P("beta", range=(1e-4, 10.0), unit="1"),
            "gamma_anchor": P("gamma_anchor", value=g_truth_250),
            "gamma_anchor_tol": P("gamma_anchor_tol", value=0.5),
            "gamma_anchor_T": P("gamma_anchor_T", value=250.0, unit="K"),
            "E_a_prior": P("E_a_prior", value=240.0),
            "E_a_prior_tol": P("E_a_prior_tol", value=20.0),
            "b_p_prior_log10": P("b_p_prior_log10", value=2.0, unit="log10(1)"),
            "b_p_prior_log10_tol": P("b_p_prior_log10_tol", value=0.7, unit="log10(1)"),
        },
        datasets={"g2_vs_T": DataSet("g2_vs_T", Tag.A, "synthetic",
                                     ["T", "g2", "err", "w", "dx"], rows)},
    )
    fit = fit_phase0(card, n_starts=16)
    assert fit.passed, f"max |resid| = {np.max(np.abs(fit.residuals)):.4f}"


# --------------------------------------------------------------- Module D (loading)

@check("F1b: mu->0 recovers the F1 identity g2=eps")
def _():
    from fsim_core.loading import f1b_g2
    for eps in (0.01, 0.1, 0.3):
        assert abs(f1b_g2(1e-6, eps) / eps - 1.0) < 1e-4, eps


@check("F1b: saturation limit 2 eps/(1+eps)^2; monotone increasing in mu")
def _():
    from fsim_core.loading import f1b_g2
    eps = 0.1
    assert abs(f1b_g2(50.0, eps) - 2 * eps / (1 + eps) ** 2) < 1e-9
    mus = np.linspace(0.01, 10, 300)
    g = f1b_g2(mus, eps)
    assert np.all(np.diff(g) > -1e-14)


@check("F1b: drive factor bounded in (1, 2)")
def _():
    from fsim_core.loading import drive_factor
    for mu in (0.1, 0.33, 1.0, 3.0, 20.0):
        f = drive_factor(mu, 0.05)
        assert 1.0 < f < 2.0, (mu, f)


@check("F1b: MC second method (finite t_x, three operating points)")
def _():
    from fsim_core.loading import f1b_g2, mc_pulsed_g2
    eps = 0.12
    for mu in (0.33, 1.0, 3.0):
        g2, se = mc_pulsed_g2(mu, 0.6, 0.6 * eps, seed=3)
        assert abs(g2 - f1b_g2(mu, eps)) < 4 * se + 2e-3, mu


@check("background law exact for Poissonian background over the cascade (MC)")
def _():
    from fsim_core.loading import f1b_g2, loading_probs, mc_pulsed_g2
    mu, t_x, eps, b = 0.8, 0.5, 0.15, 0.06
    _, P1, P2 = loading_probs(mu)
    s = t_x * (P1 + P2 * (1 + eps))          # mean signal photons per pulse
    rho = s / (s + b)
    expected = 1 - rho**2 * (1 - f1b_g2(mu, eps))
    g2, se = mc_pulsed_g2(mu, t_x, t_x * eps, b=b, seed=4)
    assert abs(g2 - expected) < 4 * se + 2e-3, (g2, expected)


@check("injection channels: compose additively; E_act=0 is T-independent")
def _():
    from fsim_core.loading import BackgroundChannel, b_injection
    leak = BackgroundChannel("leak", A=0.01, m=1.0)
    wl = BackgroundChannel("WL-EL", A=0.5, m=2.0, E_act=100.0)
    assert abs(b_injection([leak], 2.0, 80) - 0.02) < 1e-12
    assert abs(b_injection([leak], 2.0, 300) - b_injection([leak], 2.0, 80)) < 1e-15
    tot = b_injection([leak, wl], 2.0, 200)
    assert tot > b_injection([leak], 2.0, 200)
    assert b_injection([wl], 2.0, 300) > b_injection([wl], 2.0, 100)


# --------------------------------------------------------------- Module A (thermal)

@check("thermal: FD disk solve matches uniform-flux closed form within 8%")
def _():
    from fsim_core.thermal import fd_disk_rth
    a, k = 1e-6, 55.0
    R_dom = 12 * a
    R_fd = fd_disk_rth(a, k, R_dom=R_dom, Z_dom=R_dom, n=240)
    R_fd += 1.0 / (2 * np.pi * k * R_dom)  # far-field resistance cut off by the finite domain
    R_uf = 8.0 / (3 * np.pi**2 * k * a)
    assert abs(R_fd - R_uf) / R_uf < 0.08, (R_fd, R_uf)


@check("thermal: stack limits (bare disk, pillar 1D term, cone < pillar)")
def _():
    from fsim_core.thermal import Layer, Stack, stack_rth
    a, k = 0.5e-6, 55.0
    bare = Stack(layers=[], k_sub300=k, sub_alpha=0.0)
    assert abs(stack_rth(a, bare, 300.0) - 1 / (4 * k * a)) < 1e-9
    t = 1.5e-6
    pil = Stack(layers=[Layer(t, 15.0, 0.0, spread=False)], k_sub300=k, sub_alpha=0.0)
    con = Stack(layers=[Layer(t, 15.0, 0.0, spread=True)], k_sub300=k, sub_alpha=0.0)
    r_pil = stack_rth(a, pil, 300.0) - 1 / (4 * k * a)
    assert abs(r_pil - t / (15.0 * np.pi * a**2)) < 1e-9
    assert stack_rth(a, con, 300.0) < stack_rth(a, pil, 300.0)


@check("thermal: R_th decreasing in mesa radius; bare Si < bare GaAs substrate")
def _():
    from fsim_core.thermal import Layer, Stack, stack_rth
    gaas = Stack(layers=[Layer(1.5e-6, 15.0, spread=False)], k_sub300=55.0)
    radii = np.linspace(0.15e-6, 1.5e-6, 30)
    R = [stack_rth(a, gaas) for a in radii]
    assert np.all(np.diff(R) < 0)
    a = 0.5e-6
    bare_si = Stack(layers=[], k_sub300=148.0, sub_alpha=1.35)
    bare_ga = Stack(layers=[], k_sub300=55.0)
    assert stack_rth(a, bare_si) < stack_rth(a, bare_ga)


@check("thermal: k(T) runaway reported as inf, not extrapolated garbage")
def _():
    from fsim_core.thermal import Layer, Stack, t_junction
    st = Stack(layers=[Layer(1.5e-6, 8.0, 1.25)], k_sub300=55.0)
    Tj = t_junction(400e-6, 0.15e-6, st, 77.0)
    assert np.isinf(Tj)


@check("thermal: t_junction fixed point self-consistent; P=0 -> T_hs")
def _():
    from fsim_core.thermal import Layer, Stack, stack_rth, t_junction
    st = Stack(layers=[Layer(1.5e-6, 15.0)], k_sub300=55.0)
    a, P, T_hs = 0.3e-6, 200e-6, 77.0
    assert t_junction(0.0, a, st, T_hs) == T_hs
    Tj = t_junction(P, a, st, T_hs)
    assert Tj > T_hs
    assert abs(Tj - (T_hs + P * stack_rth(a, st, 0.5 * (T_hs + Tj)))) < 1e-2


@check("electrical separation: g2_electrical == g2_of_T at T_j with I channels")
def _():
    from fsim_core.integrator import g2_electrical
    from fsim_core.loading import BackgroundChannel
    from fsim_core.thermal import Layer, Stack, t_junction
    p = dict(BASE_P)
    p["mu"] = 0.5
    p["channels"] = [BackgroundChannel("leak", A=0.001, m=1.0)]
    st = Stack(layers=[Layer(1.5e-6, 15.0)], k_sub300=55.0)
    T_hs, I, V, a = 77.0, 50e-6, 1.9, 0.4e-6
    pt = g2_electrical(T_hs, I, V, p, a, st)
    Tj = t_junction(I * V, a, st, T_hs)
    ref = g2_of_T(Tj, {**p, "I": I})
    assert abs(pt.T - Tj) < 1e-9 and abs(pt.g2 - ref.g2) < 1e-12


# --------------------------------------------- F-series unification (foundation doc)

@check("Theorem-0 slice: eps=1 gives g2(mu->0)->1 and g2(mu->inf)->0.5")
def _():
    from fsim_core.loading import f1b_g2
    assert abs(f1b_g2(1e-4, 1.0) - 1.0) < 1e-3
    assert abs(f1b_g2(200.0, 1.0) - 0.5) < 1e-6


@check("F1b expansion: g2(mu) ~ eps[1 + mu(1/3 - eps)] at small mu (foundation form)")
def _():
    from fsim_core.loading import f1b_g2
    for eps in (0.05, 0.3979, 0.6):
        mu = 0.02
        exact = float(f1b_g2(mu, eps))
        approx = eps * (1 + mu * (1.0 / 3.0 - eps))
        assert abs(exact - approx) < 2e-4 * max(eps, 0.1), (eps, exact, approx)


@check("F5 aperture lemma vs MC at (1,1), (1,0.3), (1,1,1), (1,0.5,0.1); equal-p 1-1/N")
def _():
    from fsim_core.loading import aperture_g2, mc_aperture_g2
    assert abs(aperture_g2([0.7] * 4) - (1 - 0.25)) < 1e-12
    for ps in ([1.0, 1.0], [1.0, 0.3], [1.0, 1.0, 1.0], [1.0, 0.5, 0.1]):
        # p=1 emitters fire every pulse: closed form still applies
        g2, se = mc_aperture_g2(ps, seed=7)
        assert abs(g2 - aperture_g2(ps)) < 4 * se + 1e-3, ps


# --------------------------------------------------------------- Module B (cavity)

@check("mode-tracking rule: detuning vanishes at T_target, positive (blue) below it")
def _():
    from fsim_core.cavity import tracking_detuning
    Tt = 120.0
    assert abs(tracking_detuning(Tt, 1.88, Tt)) < 1e-12
    assert tracking_detuning(20.0, 1.88, Tt) > 0
    assert tracking_detuning(300.0, 1.88, Tt) < 0


@check("F_eff limits: kappa>>Gamma -> F_P; kappa<<Gamma -> F_P*kappa/Gamma")
def _():
    from fsim_core.cavity import purcell_eff
    assert abs(purcell_eff(10.0, 1e4, 6.5) - 10.0) < 1e-2
    assert abs(purcell_eff(10.0, 1e-3, 6.5) - 10.0 * 1e-3 / 6.5) < 1e-5


@check("cavity filter suppresses XX vs no-filter (F6 i)")
def _():
    e_cav = epsilon(5.9, 2.0, 2.0, kappa=0.5).eps
    e_none = 1.0  # no filter: both lines fully collected
    assert e_cav < 0.1 < e_none


@check("Lemma 1: SiN beta is brightness-only -- never referenced in the g2 path")
def _():
    root = Path(__file__).resolve().parents[1]
    for f in ("integrator.py", "fitting.py", "spectral.py", "loading.py"):
        src = (root / "fsim_core" / f).read_text(encoding="utf-8")
        assert "sin_beta" not in src and "beta_sin" not in src, f
    from fsim_core.cavity import sin_beta_brightness
    assert abs(sin_beta_brightness(0.8, 0.6) - 0.48) < 1e-12


@check("collection gain raises rho but leaves eps untouched (F6 ii)")
def _():
    p = dict(BASE_P)
    a = g2_of_T(200.0, p)
    b = g2_of_T(200.0, {**p, "collection_gain": 10.0})
    assert b.rho > a.rho and abs(b.eps - a.eps) < 1e-15 and b.g2 < a.g2


# -------------------------------------------------------------- device assembly

@check("device: design YAML round-trip preserves every block")
def _():
    import tempfile

    from fsim_core.device import DeviceDesign

    d = DeviceDesign(name="rt-test")
    d.dot.delta_xx = 4.2
    d.drive.I_uA = 33.0
    d.cavity.enabled = True
    d.thermal.layers.append({"name": "buffer", "t_um": 2.0, "k300": 44.0,
                             "alpha": 1.25, "spread": True})
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "d.yaml"
        d.save(p)
        e = DeviceDesign.load(p)
    assert e.dot.delta_xx == 4.2 and e.drive.I_uA == 33.0
    assert e.cavity.enabled and len(e.thermal.layers) == 2


@check("device: evaluate() composes the same chain as the module-level calls")
def _():
    from fsim_core.device import DeviceDesign, class_proxy_params, evaluate
    from fsim_core.loading import BackgroundChannel, b_injection
    from fsim_core.spectral import gamma_of_T
    from fsim_core.thermal import Layer, Stack, t_junction

    d = DeviceDesign()
    d.filter.auto_w = False
    d.filter.w = 3.0
    res = evaluate(d, T_grid=[d.thermal.T_hs])
    s = res["scalars"]

    pr = class_proxy_params()
    st = Stack(layers=[Layer(1.5e-6, 15.0, 0.5)], k_sub300=55.0)
    Tj = t_junction(d.drive.duty * d.drive.I_uA * 1e-6 * d.drive.V, 0.5e-6,
                    st, d.thermal.T_hs)
    assert abs(s["T_j_op"] - Tj) < 1e-6
    gam = float(gamma_of_T(Tj, pr["gamma0"], pr["a_ac"], pr["b_lo"], pr["E_lo"]))
    spec = epsilon(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=3.0)
    assert abs(s["eps_op"] - spec.eps) < 1e-12
    S = float(retention(Tj, pr["a_esc"], pr["E_a"], pr["b_p"], pr["E_b"]))
    chan = [BackgroundChannel("injection", A=d.drive.b_e, m=d.drive.b_e_m,
                              E_act=d.drive.b_e_Eact, I_ref=d.drive.I_ref_uA)]
    B = pr["b0"] + pr["beta"] * (1 - S) + float(b_injection(chan, d.drive.I_uA, Tj))
    assert abs(s["rho_op"] - S / (S + B)) < 1e-12
    from fsim_core.loading import f1b_g2 as _f
    assert abs(s["g2_op"] - (1 - (S / (S + B)) ** 2
                             * (1 - float(_f(d.drive.mu, spec.eps))))) < 1e-12


@check("device: cavity gain raises rho; runaway design reports inf T_j")
def _():
    from fsim_core.device import DeviceDesign, evaluate

    d = DeviceDesign()
    base = evaluate(d, T_grid=[77.0])["scalars"]
    d.cavity.enabled = True
    cav = evaluate(d, T_grid=[77.0])["scalars"]
    assert cav["rho_op"] > base["rho_op"]

    hot = DeviceDesign()
    hot.thermal.mesa_diameter_um = 0.15
    hot.drive.I_uA = 400.0
    hot.thermal.layers[0]["alpha"] = 1.25
    s = evaluate(hot, T_grid=[77.0])["scalars"]
    assert s["runaway"] and np.isinf(s["T_j_op"])


# ------------------------------------------------------------ three-layer rule

@check("three-layer rule: fsim_core and fsim_viz never import GUI frameworks")
def _():
    root = Path(__file__).resolve().parents[1]
    for pkg in ("fsim_core", "fsim_viz"):
        for f in (root / pkg).glob("*.py"):
            src = f.read_text(encoding="utf-8")
            for banned in ("streamlit", "plotly"):
                assert banned not in src, f"{f.name} references {banned}"


@check("units: KB*300K ~ 25.85 meV (meV/K convention holds)")
def _():
    assert abs(KB * 300.0 - 25.852) < 0.01


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
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
