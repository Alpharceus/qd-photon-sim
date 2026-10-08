"""Regression suite for fsim_core.resonant_source (Fock-state source tier).

Every check compares the module against a number it did not produce: closed
forms, identities, an independent second method (Gillespie Monte Carlo), the
legacy lindblad.pulsed_counting, or published values. Groups: (r) F1 resonant
two-level source, (t) F2 two-photon excitation, (s) F3 spectral diffusion, (m)
F4 boson-sampling figure of merit and card.

Run: python verify/verify_resonant_source.py   (exit code 0 iff all pass)
"""
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import lindblad as lb
from fsim_core import resonant_source as rs

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def relerr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


F_REP = 80.0   # MHz, the work-order table
# [V] work-order table (values produced by lindblad.pulsed_counting, printed
# to 4 decimals): (Gamma ns^-1, tau_p ps, g2(0), photons/pulse).
R1_TABLE = [(1.0, 2.0, 0.0005, 1.0002), (1.0, 10.0, 0.0025, 1.0012),
            (1.0, 50.0, 0.0122, 1.0061), (10.0, 10.0, 0.0238, 1.0119),
            (10.0, 50.0, 0.0992, 1.0485)]


# ------------------------------------------------------------------ F1

@check("(r0) pulse_envelope: area = int Omega dt to 1e-12 (rect, Gauss); Gauss intensity FWHM = duration")
def _r0():
    ok = True
    for shape in ("rect", "gauss"):
        for area in (0.3, np.pi, 5.0):
            for dur in (0.002, 0.05, 1.3):
                win = rs.pulse_window(shape, dur)
                pts = [dur] if shape == "rect" else [0.5 * win]
                val, _ = quad(lambda t: float(rs.pulse_envelope(t, shape, area, dur)),
                              0.0, win, points=pts, epsabs=0, epsrel=1e-13, limit=400)
                ok &= abs(val - area) / area < 1e-12
    dur = 0.37
    win = rs.pulse_window("gauss", dur)
    pk = float(rs.pulse_envelope(0.5 * win, "gauss", np.pi, dur)) ** 2
    hi = float(rs.pulse_envelope(0.5 * win + 0.5 * dur, "gauss", np.pi, dur)) ** 2
    lo = float(rs.pulse_envelope(0.5 * win - 0.5 * dur, "gauss", np.pi, dur)) ** 2
    ok &= abs(hi / pk - 0.5) < 1e-12 and abs(lo / pk - 0.5) < 1e-12
    return bool(ok)


@check("(r1) rect pi pulse: pulsed_photon_statistics == lindblad.pulsed_counting (rel 1e-10) and the work-order table (4 dp)")
def _r1():
    worst, ok = 0.0, True
    for G, tp_ps, g2_tab, n_tab in R1_TABLE:
        tau = tp_ps * 1e-3
        mine = rs.pulsed_photon_statistics(rs.TwoLevel(G), "rect", np.pi, tau, F_REP)
        s_on = lb.build_system(levels=2, gamma_X_ns=G, omega_ns=np.pi / tau)
        s_off = lb.build_system(levels=2, gamma_X_ns=G)
        ref = lb.pulsed_counting(s_on.L, s_off.L, lb.counting_jump(s_off), 2, tau,
                                 1.0e3 / F_REP - tau)
        e = max(relerr(mine["g2"], ref["g2"]), relerr(mine["mean_photons"], ref["mean_counts"]))
        worst = max(worst, e)
        ok &= e < 1e-10
        ok &= abs(round(mine["g2"], 4) - g2_tab) < 1e-9 and abs(round(mine["mean_photons"], 4) - n_tab) < 1e-9
    print(f"    worst rel diff vs lindblad.pulsed_counting = {worst:.1e}")
    return bool(ok)


@check("(r2) small-Gamma tau_p slope: g2/(Gamma tau_p) in [0.245, 0.255] at 1e-3, -> 1/4 (DR)")
def _r2():
    ok = True
    for gt in (1e-3, 1e-4):
        r = rs.pulsed_photon_statistics(rs.TwoLevel(1.0), "rect", np.pi, gt, F_REP)
        slope = r["g2"] / gt
        print(f"    Gamma*tau_p = {gt:g}: g2/(Gamma tau_p) = {slope:.6f}")
        ok &= 0.245 <= slope <= 0.255
    ok &= abs(slope - 0.25) < 1e-3
    return bool(ok)


@check("(r3) area theorem at Gamma = 0: P_e(end) = sin^2(theta/2) to 1e-10 (rect, Gauss; pi/2, pi, 3pi/2, 2pi)")
def _r3():
    worst = 0.0
    for shape in ("rect", "gauss"):
        for th in (0.5 * np.pi, np.pi, 1.5 * np.pi, 2.0 * np.pi):
            rho, _ = rs.propagate_pulse(rs.TwoLevel(0.0), shape, th, 0.04)
            worst = max(worst, abs(rho[1, 1].real - np.sin(th / 2) ** 2))
    print(f"    worst |P_e - sin^2(theta/2)| = {worst:.1e}")
    return worst < 1e-10


@check("(r4) HOM limits: I > 0.999 at Gamma tau_p = 1e-3 (gamma* = 0); with gamma* > 0, -> Gamma/(Gamma+gamma*) rel 1e-3")
def _r4():
    I0, info = rs.pulsed_indistinguishability(rs.TwoLevel(1.0), "rect", np.pi, 1e-3, return_info=True)
    print(f"    gamma*=0: I = {I0:.6f} (quadrature converged={info['converged']}, n={info['n_quad']})")
    ok = I0 > 0.999 and info["converged"]
    for gs in (0.5, 2.0, 10.0):
        I = rs.pulsed_indistinguishability(rs.TwoLevel(1.0, gamma_star_ns=gs), "rect", np.pi, 1e-3)
        ref = 1.0 / (1.0 + gs)
        print(f"    gamma*={gs:g}: I = {I:.6f} vs Gamma/(Gamma+gamma*) = {ref:.6f}")
        ok &= relerr(I, ref) < 1e-3
    return bool(ok)


@check("(r5) re-excitation lowers I: strictly decreasing over Gamma tau_p in {1e-3, 1e-2, 0.1, 0.5}, gamma* = 0; quadrature converged")
def _r5():
    vals, conv = [], True
    for gt in (1e-3, 1e-2, 0.1, 0.5):
        I, info = rs.pulsed_indistinguishability(rs.TwoLevel(1.0), "rect", np.pi, gt, return_info=True)
        vals.append(I)
        conv &= info["converged"]
        print(f"    Gamma*tau_p = {gt:g}: I = {I:.6f} (conv_err {info['conv_err']:.1e}, n={info['n_quad']})")
    return bool(conv and all(a > b for a, b in zip(vals, vals[1:])))


@check("(r6) second method: Gillespie Monte Carlo g2(0) at Gamma tau_p = 0.1 within 3 standard errors (seed fixed)")
def _r6():
    tl = rs.TwoLevel(1.0)
    det = rs.pulsed_photon_statistics(tl, "rect", np.pi, 0.1, F_REP, periodic=False)
    mc = rs.pulsed_counting_mc(tl, np.pi, 0.1, F_REP, n_traj=2_000_000, seed=20261008)
    z = (mc["g2"] - det["g2"]) / mc["g2_se"]
    zn = (mc["mean"] - det["mean_photons"]) / mc["mean_se"]
    print(f"    deterministic g2 = {det['g2']:.6f}, MC g2 = {mc['g2']:.6f} +- {mc['g2_se']:.6f} (z = {z:+.2f}); "
          f"photons z = {zn:+.2f}")
    return abs(z) < 3.0 and abs(zn) < 3.0


@check("(r7) three-region HOM assembly: free decay from |e> (area 0, 2 ns window) gives Gamma/(Gamma+gamma*) for any window split (rel 1e-7)")
def _r7():
    rho = np.zeros((2, 2), dtype=complex)
    rho[1, 1] = 1.0
    ok = True
    for gs in (0.0, 0.7):
        for dur in (0.3, 2.0):
            I = rs.pulsed_indistinguishability(rs.TwoLevel(1.0, gamma_star_ns=gs), "rect", 0.0, dur, rho0=rho)
            ok &= relerr(I, 1.0 / (1.0 + gs)) < 1e-7
    return bool(ok)


@check("(r8) scaling: g2 and I depend on the pulse only through Gamma tau_p (Gamma=2,tau) == (Gamma=1,2 tau), rel 1e-9; rect and Gauss")
def _r8():
    ok = True
    for shape in ("rect", "gauss"):
        a = rs.pulsed_photon_statistics(rs.TwoLevel(2.0), shape, np.pi, 0.02, 1.0)
        b = rs.pulsed_photon_statistics(rs.TwoLevel(1.0), shape, np.pi, 0.04, 1.0)
        ok &= relerr(a["g2"], b["g2"]) < 1e-9 and relerr(a["mean_photons"], b["mean_photons"]) < 1e-9
    Ia = rs.pulsed_indistinguishability(rs.TwoLevel(2.0), "rect", np.pi, 0.02)
    Ib = rs.pulsed_indistinguishability(rs.TwoLevel(1.0), "rect", np.pi, 0.04)
    ok &= relerr(Ia, Ib) < 1e-7
    return bool(ok)


@check("(r9) Purcell convention: Gamma = Gamma_0 (1 + F_P kappa/(line+kappa)) on resonance; detuned Lorentzian overlap; F_P = 0 -> Gamma_0")
def _r9():
    ok = True
    g0, F, line, kap = 1.0, 20.0, 0.2, 5.0
    ok &= relerr(rs.purcell_rate(g0, F, 0.0, line, kap), g0 * (1.0 + F * kap / (line + kap))) < 1e-14
    d = 1.7
    L = kap * (line + kap) / 4.0 / (d ** 2 + ((line + kap) / 2.0) ** 2)
    ok &= relerr(rs.purcell_rate(g0, F, d, line, kap), g0 * (1.0 + F * L)) < 1e-14
    ok &= rs.purcell_rate(g0, 0.0, 0.0, line, kap) == g0
    return bool(ok)


# ------------------------------------------------------------------ F2

def _fit_rabi(times, p):
    from scipy.optimize import curve_fit
    f = lambda t, w, a: a * np.sin(w * t / 2.0) ** 2
    w0 = times[np.argmax(p)]
    popt, _ = curve_fit(f, times, p, p0=[np.pi / (times[-1] / 6.0), 1.0], xtol=1e-14, ftol=1e-14)
    return float(popt[0])


# (Omega, delta) in rad/ns: the three work-order pairs plus smaller ratios.
T1_CASES = [(1.0, 20.0), (0.5, 20.0), (1.0, 50.0), (0.5, 50.0), (0.25, 20.0)]


@check("(t1) closed ladder: fitted Omega_eff == Omega^2/(2 delta) to 1e-3 for Omega/delta <= 0.04 and == exact (R-delta)/2 to 1e-6 "
       "for all, incl. the work-order pairs; the 0.05 deviation is the predicted -(Omega/delta)^2/2")
def _t1():
    ok = True
    for om, de in T1_CASES:
        lam = 0.5 * (np.sqrt(de ** 2 + 2 * om ** 2) - de)           # exact, derivation in the module comment
        times = np.linspace(0.0, 6.0 * np.pi / lam, 1500)
        P = rs.tpe_population_series(rs.Cascade(0.0, delta_ns=de), om, times)
        w = _fit_rabi(times, P[:, 2])
        lead = om ** 2 / (2.0 * de)
        e_lead, e_ex = w / lead - 1.0, w / lam - 1.0
        ratio = om / de
        print(f"    Omega={om:g}, delta={de:g}: Omega/delta={ratio:.4f}, fit/lead-1={e_lead:+.3e} "
              f"(predicted {-ratio ** 2 / 2:+.3e}), fit/exact-1={e_ex:+.1e}, peak P_XX={P[:, 2].max():.5f}")
        ok &= abs(e_ex) < 1e-6
        ok &= abs(e_lead - (-ratio ** 2 / 2.0)) < 0.02 * ratio ** 2 / 2.0 + 1e-6   # first correction, 2% of itself
        if ratio <= 0.04 + 1e-12:
            ok &= abs(e_lead) < 1e-3
        ok &= P[:, 2].max() > 1.0 - 5e-4
    return bool(ok)


@check("(t1b) sign convention: populations are unchanged under delta -> -delta and Omega -> -Omega (closed ladder, 1e-12)")
def _t1b():
    times = np.linspace(0.0, 400.0, 41)
    base = rs.tpe_population_series(rs.Cascade(0.0, delta_ns=20.0), 1.0, times)
    a = rs.tpe_population_series(rs.Cascade(0.0, delta_ns=-20.0), 1.0, times)
    b = rs.tpe_population_series(rs.Cascade(0.0, delta_ns=20.0), -1.0, times)
    return float(max(np.max(np.abs(a - base)), np.max(np.abs(b - base)))) < 1e-12


@check("(t2) cascade cap: I_{XX->X}/I_X == Gamma_XX/(Gamma_XX+Gamma_X) to 1e-9 on the Q2 table from a TPE-prepared |XX> "
       "(adiabatic Gaussian, Gamma tau_p -> 0); <= 2/3 for gamma_XX = 2 gamma_X without escape")
def _t2():
    delta, dur = 50.0, 1.0
    th = rs.tpe_pi_area(delta, "gauss", dur)
    # independent adiabatic-phase estimate (module comment): phi = int lambda_-(Omega(t)) dt at the
    # amplitude-refined area. The exact condition is phi = pi up to the 2nd-order adiabatic correction,
    # which must be < 1e-3 relative at delta = 50 and fall ~ 1/delta^2 (factor ~4 per doubling).
    def phase_deficit(de):
        t_a = rs.tpe_pi_area(de, "gauss", dur)
        win = rs.pulse_window("gauss", dur)
        f = lambda t: rs.tpe_exact_rabi(float(rs.pulse_envelope(t, "gauss", t_a, dur)), de)
        phi, _ = quad(f, 0.0, win, points=[0.5 * win], epsabs=0, epsrel=1e-12, limit=400)
        return abs(phi - np.pi)
    d50, d100, d200 = phase_deficit(50.0), phase_deficit(100.0), phase_deficit(200.0)
    print(f"    TPE-pi area (delta=50) = {th:.10f} rad; |adiabatic phase - pi| = {d50:.2e} (delta 50), {d100:.2e} (100), {d200:.2e} (200)")
    ok = d50 / np.pi < 1e-3 and 3.3 < d50 / d100 < 4.2 and 3.3 < d100 / d200 < 4.2
    table = [(1, 2, 0, 0, 0), (1, 2, 0, 0, 5), (1, 2, 10, 20, 5), (1, 3, 0, 0, 2), (1, 2, 5, 0, 1)]
    for gx, gxx, kx, kxx, gs in table:
        c = rs.Cascade(gx, gxx, kx, kxx, deph=((gs, (0, 1, 0)),) if gs else (), delta_ns=delta)
        r = rs.tpe_cascade_hom(c, "gauss", th, dur)
        s = c.system()
        rx = np.zeros((3, 3), dtype=complex)
        rx[1, 1] = 1.0
        IX = lb.indistinguishability(s.L, 3, rx, np.sqrt(gx) * s.ops["c_X"])
        GX, GXX = gx + kx, gxx + kxx
        cap = GXX / (GXX + GX)
        print(f"    (gX,gXX,kX,kXX,g*)=({gx},{gxx},{kx},{kxx},{gs}): I_X = {IX:.6f} (closed form {GX / (GX + gs):.6f}), "
              f"ratio = {r['I'] / IX:.9f} vs {cap:.9f}, P_XX = {r['P'][2]:.12f}")
        ok &= relerr(r["I"] / IX, cap) < 1e-9 and relerr(IX, GX / (GX + gs)) < 1e-9
        ok &= r["P"][2] > 1.0 - 1e-9
    c = rs.Cascade(1.0, 2.0, delta_ns=delta)
    ok &= rs.tpe_cascade_hom(c, "gauss", th, dur)["I"] <= 2.0 / 3.0 + 1e-9
    return bool(ok)


@check("(t3) re-excitation suppression: X-line g2(0) under TPE < resonant pi-pulse g2(0) at Gamma tau_p in {0.05, 0.1} "
       "(same Gamma_X, same Gaussian FWHM; also below the rectangular value)")
def _t3():
    ok = True
    delta = 2000.0                      # rad/ns, ~1.3 meV = half a 2.6 meV biexciton binding energy (E-class, InGaAs-like)
    for gt in (0.05, 0.1):
        th = rs.tpe_pi_area(delta, "gauss", gt)
        tpe = rs.tpe_photon_statistics(rs.Cascade(1.0, 2.0, delta_ns=delta), "gauss", th, gt)
        res_g = rs.pulsed_photon_statistics(rs.TwoLevel(1.0), "gauss", np.pi, gt)
        res_r = rs.pulsed_photon_statistics(rs.TwoLevel(1.0), "rect", np.pi, gt)
        print(f"    Gamma*tau_p = {gt}: g2 TPE = {tpe['g2']:.5f}; resonant pi Gauss = {res_g['g2']:.5f}, rect = {res_r['g2']:.5f}")
        ok &= tpe["g2"] < res_g["g2"] and tpe["g2"] < res_r["g2"]
    return bool(ok)


@check("(t4) TPE propagation convergence: midpoint scheme is 2nd order (error ratio 3-5 per halving) and the production grid is within 1e-4 of the "
       "Richardson limit")
def _t4():
    delta, gt = 500.0, 0.05
    th = rs.tpe_pi_area(delta, "gauss", gt)
    c = rs.Cascade(1.0, 2.0, delta_ns=delta)
    pc = rs._pieces_cascade(c)
    n0 = rs._default_steps("gauss", gt, None, pc.max_dt)
    g = {k: rs.tpe_photon_statistics(c, "gauss", th, gt, n_steps=int(k * n0))["g2"]
         for k in (0.25, 0.5, 1.0, 2.0)}
    d1, d2 = abs(g[0.25] - g[0.5]), abs(g[0.5] - g[1.0])
    ratio = d1 / d2
    rich = g[2.0] + (g[2.0] - g[1.0]) / 3.0
    err = abs(g[1.0] - rich) / rich
    print(f"    n0 = {n0}: successive differences ratio = {ratio:.2f}; production-grid error vs Richardson = {err:.1e}")
    return 3.0 < ratio < 5.0 and err < 1e-4


# ------------------------------------------------------------------ F3

def _erfcx_closed(Gamma, s):
    """Closed form of the Gaussian-averaged Lorentzian overlap (gamma* = 0), s = std of the detuning difference."""
    return float(np.sqrt(np.pi / 2.0) * (Gamma / s) * erfcx(Gamma / (np.sqrt(2.0) * s)))


@check("(s1) gamma* = 0: Gaussian-averaged cross-HOM (hom_remote, hom_gauss_average, hom_vs_separation at dt >> tau_SD) == "
       "sqrt(pi/2)(Gamma/s) erfcx(Gamma/(sqrt2 s)) to rel 1e-8 for Gamma/s in {0.1, 1, 10}")
def _s1():
    ok = True
    G = 1.7
    for r in (0.1, 1.0, 10.0):
        s = G / r
        cf = _erfcx_closed(G, s)
        v1 = rs.hom_remote(G, 0.0, s / np.sqrt(2.0))
        v2 = rs.hom_gauss_average(G, 0.0, 0.0, s)
        v3 = rs.hom_vs_separation(1e6, G, 0.0, s / np.sqrt(2.0), 1.0)
        e = max(relerr(v1, cf), relerr(v2, cf), relerr(v3, cf))
        print(f"    Gamma/s = {r:g}: closed form {cf:.10f}, quadrature {v1:.10f} (worst rel err {e:.1e})")
        ok &= e < 1e-8
    # with fast dephasing: Gamma -> Gamma + gamma* inside, prefactor Gamma/(Gamma+gamma*)
    gs = 0.9
    Gt = G + gs
    for r in (0.1, 1.0, 10.0):
        s = Gt / r
        cf = G / Gt * _erfcx_closed(Gt, s)
        ok &= relerr(rs.hom_remote(G, gs, s / np.sqrt(2.0)), cf) < 1e-8
    return bool(ok)


@check("(s2) Monte Carlo over sampled detuning pairs (200k iid pairs; 200k OU-correlated pairs at dt = tau_SD) reproduces (s1) "
       "and hom_vs_separation within 3 standard errors")
def _s2():
    rng = np.random.default_rng(20261008)
    ok = True
    n = 200_000
    G, sigma = 1.0, 1.0 / np.sqrt(2.0)              # Gamma/s = 1 with s = sqrt2 sigma
    d1 = sigma * rng.standard_normal(n)
    d2 = sigma * rng.standard_normal(n)
    I = rs._hom_cross_batch(G, 0.0, d1, d2)         # per-pair NUMERICAL overlap (uses both detunings)
    m, se = float(I.mean()), float(I.std(ddof=1) / np.sqrt(n))
    cf = _erfcx_closed(G, np.sqrt(2.0) * sigma)
    z = (m - cf) / se
    print(f"    remote: MC {m:.5f} +- {se:.5f} vs closed form {cf:.5f} (z = {z:+.2f})")
    ok &= abs(z) < 3.0
    # OU-correlated pairs at dt = tau_SD, with fast dephasing
    tau, gs = 3.0, 0.4
    rho = np.exp(-1.0)
    x1 = sigma * rng.standard_normal(n)
    x2 = rho * x1 + np.sqrt(1.0 - rho ** 2) * sigma * rng.standard_normal(n)
    I2 = rs._hom_cross_batch(G, gs, x1, x2)
    m2, se2 = float(I2.mean()), float(I2.std(ddof=1) / np.sqrt(n))
    q = rs.hom_vs_separation(tau, G, gs, sigma, tau)
    z2 = (m2 - q) / se2
    print(f"    OU dt = tau_SD: MC {m2:.5f} +- {se2:.5f} vs hom_vs_separation {q:.5f} (z = {z2:+.2f})")
    ok &= abs(z2) < 3.0
    return bool(ok)


@check("(s3) limits: dt -> 0 recovers (d1) Gamma/(Gamma+gamma*) (rel 1e-12); dt >> tau_SD recovers hom_remote(offset = 0) (rel 1e-9)")
def _s3():
    ok = True
    for G, gs, sig, tau in [(1.0, 0.0, 0.3, 5.0), (1.0, 2.0, 0.7, 2.0), (3.0, 0.5, 1.5, 40.0)]:
        v0 = rs.hom_vs_separation(0.0, G, gs, sig, tau)
        vsmall = rs.hom_vs_separation(1e-12 * tau, G, gs, sig, tau)
        vlong = rs.hom_vs_separation(60.0 * tau, G, gs, sig, tau)
        rem = rs.hom_remote(G, gs, sig)
        ok &= relerr(v0, G / (G + gs)) < 1e-12 and relerr(vsmall, G / (G + gs)) < 1e-9
        ok &= relerr(vlong, rem) < 1e-9
        print(f"    G={G}, g*={gs}: dt=0 {v0:.8f} (d1 {G / (G + gs):.8f}); dt=60 tau {vlong:.8f} vs remote {rem:.8f}")
    # a static offset lowers the remote value
    ok &= rs.hom_remote(1.0, 0.0, 0.3, offset=0.5) < rs.hom_remote(1.0, 0.0, 0.3, offset=0.0)
    return bool(ok)


@check("(s4) hom_vs_separation is non-increasing in dt (60 points from 1e-3 tau to 30 tau, three parameter sets)")
def _s4():
    ok = True
    for G, gs, sig, tau in [(1.0, 0.0, 0.3, 5.0), (1.0, 2.0, 0.7, 2.0), (3.0, 0.5, 1.5, 40.0)]:
        dts = np.geomspace(1e-3 * tau, 30.0 * tau, 60)
        v = np.array([rs.hom_vs_separation(d, G, gs, sig, tau) for d in dts])
        ok &= bool(np.all(np.diff(v) <= 1e-13))
    return ok


@check("(s5) numerical cross-HOM == closed form Gamma(Gamma+gamma*)/((Gamma+gamma*)^2+d^2) (rel 1e-10 over a grid) and depends on d1 - d2 only (1e-12)")
def _s5():
    ok = True
    for G in (0.5, 2.0):
        for gs in (0.0, 0.7, 5.0):
            d = np.array([0.0, 0.1, 1.0, 3.0, 25.0])
            ok &= relerr(rs.hom_pair(G, gs, d), G * (G + gs) / ((G + gs) ** 2 + d ** 2)) < 1e-10
    a = rs._hom_cross_batch(1.0, 0.5, np.array([3.0, -2.0]), np.array([5.0, 0.3]))
    b = rs._hom_cross_batch(1.0, 0.5, np.array([0.0, -2.3]), np.array([2.0, 0.0]))
    ok &= relerr(a, b) < 1e-12
    return bool(ok)


@check("(s6) independent quadrature: plain Gauss-Hermite (200 nodes) of the closed-form overlap == the panelled rule for Gamma/s in {1, 10} (rel 1e-8)")
def _s7():
    x, w = np.polynomial.hermite.hermgauss(200)
    ok = True
    G = 1.0
    for r in (1.0, 10.0):
        s = G / r
        gh = float(np.sum(w / np.sqrt(np.pi) * G ** 2 / (G ** 2 + (np.sqrt(2.0) * s * x) ** 2)))
        ok &= relerr(rs.hom_gauss_average(G, 0.0, 0.0, s), gh) < 1e-8
    return bool(ok)


# ------------------------------------------------------------------ F4

CARD_PATH = Path(__file__).resolve().parents[1] / "cards" / "fock-source-4K-design.yaml"
SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_fock_source.py"


def _card():
    return rs.load_fock_card(CARD_PATH)


def _L_closed(src):
    """Self-consistent Purcell overlap at zero cavity detuning from the QUADRATIC
    L (w + kappa) = kappa with w = hbar G0 (m + b L) [ueV], b = m F_P (enhanced) or
    F_P (dot): hbar G0 b L^2 + (hbar G0 m + kappa) L - kappa = 0 (closed form; the
    module solves the same equation with Brent's method)."""
    h = lb.HBAR_MEV_NS * 1e3 * src.gamma0_ns
    m, F, k = src.hom_width_multiple, src.F_P, src.kappa_ueV
    b = m * F if src.dephasing_mode == "enhanced" else F
    if b == 0.0:
        return k / (h * m + k)
    B = h * m + k
    return (-B + np.sqrt(B * B + 4.0 * h * b * k)) / (2.0 * h * b)


@check("(m1) R_N = (f_rep/N) eta^N: closed form (exp/log evaluation, rel 1e-12) and bs_figure_of_merit's R_N_Hz for N = 1, 4, 8, 20")
def _m1():
    card = _card()
    src = rs.source_from_card(card)
    ok = True
    for N in (1, 4, 8, 20):
        out = rs.bs_figure_of_merit(src, N)
        eta = out["eta"]
        ref = src.f_rep_MHz * 1.0e6 / N * np.exp(N * np.log(eta))     # independent evaluation of eta^N
        ok &= relerr(out["R_N_Hz"], ref) < 1e-12 and relerr(rs.coincidence_rate_Hz(src.f_rep_MHz, N, eta), ref) < 1e-12
    # eta is the product of the documented factors (independent recomputation from the card numbers)
    c = out["eta_components"]
    ok &= relerr(out["eta"], c["n_ph"] * c["beta"] * c["eta_out"] * c["Z"] * c["demux"] * c["circuit"] * c["det"]) < 1e-14
    L = _L_closed(src)                      # self-consistent overlap, line = hbar (Gamma + gamma*) (finding 8)
    ok &= relerr(c["beta"], (src.beta0 + src.F_P * L) / (1.0 + src.F_P * L)) < 1e-12
    for mode in ("enhanced", "dot"):
        sm = rs.source_from_card(card, dephasing_mode=mode, F_P=7.0, kappa_ueV=50.0)
        ok &= relerr(rs.source_rates(sm)["L"], _L_closed(sm)) < 1e-12
        ok &= relerr(rs.source_rates(sm)["Gamma"], sm.gamma0_ns * (1.0 + sm.F_P * _L_closed(sm))) < 1e-12
    ok &= 0.0 < out["R_N_Hz"] < src.f_rep_MHz * 1e6
    return bool(ok)


@check("(m2) hom_matrix: symmetric, unit diagonal, off-diagonals == F3 hom_vs_separation (rel 1e-12) and == erfcx closed form at gamma* = 0 (rel 1e-8)")
def _m2():
    card = _card()
    src = rs.source_from_card(card, hom_width_multiple=1.0, sigma_sd_ueV=0.4, tau_sd_ns=60.0, pulse_fwhm_ps=2.0,
                              zpl_output=True)               # filtered branch: no Z^2 penalty, H == F3 exactly
    N = 6
    out = rs.bs_figure_of_merit(src, N)
    H = np.array(out["hom_matrix"])
    r = rs.source_rates(src)
    sig = rs._rate_from_ueV(src.sigma_sd_ueV)
    dt_unit = 1.0e3 / src.f_rep_MHz
    ok = np.array_equal(H, H.T) and np.all(np.diag(H) == 1.0)
    worst_f3, worst_cf = 0.0, 0.0
    for i in range(N):
        for j in range(i + 1, N):
            dt = (j - i) * dt_unit
            f3 = rs.hom_vs_separation(dt, r["Gamma"], r["gamma_star"], sig, src.tau_sd_ns)
            s = np.sqrt(2.0 * sig ** 2 * (1.0 - np.exp(-dt / src.tau_sd_ns)))
            cf = _erfcx_closed(r["Gamma"], s)
            worst_f3 = max(worst_f3, relerr(H[i, j], f3))
            worst_cf = max(worst_cf, relerr(H[i, j], cf))
    print(f"    worst rel diff vs F3 = {worst_f3:.1e}, vs erfcx closed form = {worst_cf:.1e}")
    # farther apart -> lower I (more spectral diffusion accumulated)
    ok &= bool(H[0, 1] >= H[0, 2] >= H[0, 5])
    ok &= worst_f3 < 1e-12 and worst_cf < 1e-8
    Hc = np.array(out["hom_matrix_pulse_corrected"])
    ok &= bool(np.allclose(Hc, Hc.T, rtol=0, atol=0) and np.all(np.diag(Hc) == 1.0))
    ok &= relerr(Hc[0, 1], H[0, 1] * out["pulse_factor"]) < 1e-14
    return bool(ok)


@check("(m2b) TPE pulse factor == cascade cap Gamma_XX/(Gamma_XX+Gamma_X) with the BARE gamma_XX = 2 Gamma_0: 2/(2+6) = 0.25 at F_P L = 5 (rel 1e-12), 2/(2+1+F_P L) from the closed-form L (rel 1e-9); resonant factor in (0.9, 1]")
def _m2b():
    G0 = 1.75
    # F_P = 5, L = 1 exactly at the _pulse_physics interface (Gamma = 6 Gamma_0, gamma_XX = 2 Gamma_0): cap = 2/(2+6)
    pf5 = rs._pulse_physics("tpe", "gauss", 6.0 * G0, 0.0, 0.02, 80.0, -2.0, 2.0, 1e-6, G0)[2]
    pf10 = rs._pulse_physics("tpe", "gauss", 11.0 * G0, 0.0, 0.02, 80.0, -2.0, 2.0, 1e-6, G0)[2]
    print(f"    cap(F_P L = 5) = {pf5:.15f} (0.25), cap(F_P L = 10) = {pf10:.12f} ({2.0 / 13.0:.12f})")
    ok = relerr(pf5, 2.0 / 8.0) < 1e-12 and relerr(pf10, 2.0 / 13.0) < 1e-12
    card = _card()
    tpe = rs.source_from_card(card, transition="X", excitation="tpe", pulse_fwhm_ps=20.0, F_P=2.0, hom_width_multiple=1.0)
    res = rs.source_from_card(card, pulse_fwhm_ps=5.0, F_P=2.0, hom_width_multiple=1.5)
    pf_t = rs.bs_figure_of_merit(tpe, 4)["pulse_factor"]
    pf_r = rs.bs_figure_of_merit(res, 4)["pulse_factor"]
    cap = 2.0 / (2.0 + 1.0 + tpe.F_P * _L_closed(tpe))
    print(f"    source-level TPE pulse factor = {pf_t:.12f} (cap {cap:.12f}); resonant pulse factor = {pf_r:.6f}")
    return bool(ok and relerr(pf_t, cap) < 1e-9 and 0.9 < pf_r <= 1.0)


@check("(m2c) ZPL-filtered output multiplies eta by Z(T) = exp(-S_total) of qd_gf (rel 1e-12); Z(4 K) in (0.5, 1)")
def _m2c():
    from fsim_core.qd_gf import PhononParams, huang_rhys
    card = _card()
    ok = True
    for alpha in (0.0181, 0.027):
        off = rs.bs_figure_of_merit(rs.source_from_card(card, zpl_output=False, phonon_alpha_ps2=alpha), 4)["eta"]
        on = rs.bs_figure_of_merit(rs.source_from_card(card, zpl_output=True, phonon_alpha_ps2=alpha), 4)["eta"]
        Z = float(np.exp(-huang_rhys(PhononParams(alpha_ps2=alpha), 4.0)))
        print(f"    alpha = {alpha}: eta_on/eta_off = {on / off:.9f}, exp(-S) = {Z:.9f}")
        ok &= relerr(on / off, Z) < 1e-12 and 0.5 < Z < 1.0
    return bool(ok)


@check("(m3) card: loads, every leaf carries unit/tag/source (tags A or E only), ranges lo < hi, round-trips through dump/load; malformed entries are rejected")
def _m3():
    import tempfile
    card = _card()
    ok = len(card.params) > 0
    for name, p in card.params.items():
        ok &= bool(p.unit) and bool(p.source) and p.tag in ("A", "E")
        if p.range is not None:
            ok &= p.range[0] < p.range[1]
    ok &= card["transition"].default == "X+" and card["excitation"].default == "resonant_pi"
    ok &= card["hom_width_multiple"].range == (1.0, 3.0)                  # work-order range
    with tempfile.TemporaryDirectory() as td:
        p2 = Path(td) / "roundtrip.yaml"
        rs.dump_fock_card(card, p2)
        back = rs.load_fock_card(p2)
    ok &= back.to_dict() == card.to_dict()
    bad = [{"unit": "x", "tag": "A", "value": 1.0},                           # no source
           {"unit": "x", "source": "s", "value": 1.0},                        # no tag
           {"unit": "x", "tag": "A", "source": "s", "value": 1.0, "range": [0, 1]},   # two kinds
           {"unit": "x", "tag": "A", "source": "s", "range": [2.0, 1.0]},     # lo >= hi
           {"unit": "x", "tag": "Q", "source": "s", "value": 1.0}]            # bad tag
    for raw in bad:
        try:
            rs.fock_card_from_dict({"meta": {"name": "t"}, "params": {"p": raw}})
            ok = False
        except ValueError:
            pass
    # TPE on a charged exciton is rejected (needs the neutral X-XX cascade)
    try:
        rs.source_from_card(card, excitation="tpe", transition="X+")
        ok = False
    except ValueError:
        pass
    return bool(ok)


@check("(m4) script: runs end to end (--quick, temp dir) twice and reproduces every output file byte for byte; JSON schema and CSV R_N closed form hold")
def _m4():
    import csv
    import hashlib
    import json
    import subprocess
    import tempfile
    hashes = []
    with tempfile.TemporaryDirectory() as td:
        for k in (1, 2):
            od = Path(td) / f"run{k}"
            r = subprocess.run([sys.executable, str(SCRIPT_PATH), "--quick", "--out-dir", str(od)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stderr[-800:])
                return False
            hashes.append({f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(od.iterdir())})
        ok = hashes[0] == hashes[1] and "sweep.csv" in hashes[0] and "results.md" in hashes[0]
        od = Path(td) / "run1"
        with open(od / "sweep.csv", newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        card = _card()
        f_rep = card["f_rep_MHz"].central()
        worst = 0.0
        for row in rows:
            for N in (4, 8, 20):
                eta = float(row["eta"])
                ref = f_rep * 1e6 / N * np.exp(N * np.log(eta))
                worst = max(worst, relerr(float(row[f"R_N{N}_Hz"]), ref))
        ok &= worst < 1e-12
        js = json.loads((od / "bs_best_N8.json").read_text(encoding="utf-8"))
        for key in ("N", "f_rep_MHz", "g2", "eta", "hom_matrix", "provenance"):
            ok &= key in js
        ok &= js["N"] == 8 and len(js["g2"]) == 8 and len(js["hom_matrix"]) == 8
        ok &= all(v in ("A", "E") for v in js["provenance"].values())
        text = (od / "results.md").read_text(encoding="utf-8")
        ok &= "VERDICT:" in text
        print(f"    {len(rows)} quick rows; worst CSV R_N rel err = {worst:.1e}; files: {sorted(hashes[0])}")
    return bool(ok)


@check("(f1) review 1: unfiltered off-diagonal HOM (hom_matrix and pulse-corrected) == Z^2 x the ZPL-filtered-branch HOM (rel 1e-12), Z = exp(-S_total) independent of the module; filtered branch has no penalty; eta unchanged by the HOM factor")
def _f1():
    from fsim_core.qd_gf import PhononParams, huang_rhys
    card = _card()
    ok = True
    for alpha in (0.0181, 0.027):
        Z = float(np.exp(-huang_rhys(PhononParams(alpha_ps2=alpha), 4.0)))
        un = rs.bs_figure_of_merit(rs.source_from_card(card, zpl_output=False, phonon_alpha_ps2=alpha), 5)
        fi = rs.bs_figure_of_merit(rs.source_from_card(card, zpl_output=True, phonon_alpha_ps2=alpha), 5)
        off = ~np.eye(5, dtype=bool)
        Hu, Hf = np.array(un["hom_matrix"]), np.array(fi["hom_matrix"])
        Cu, Cf = np.array(un["hom_matrix_pulse_corrected"]), np.array(fi["hom_matrix_pulse_corrected"])
        r1 = relerr(Hu[off], Z ** 2 * Hf[off])
        r2 = relerr(Cu[off], Z ** 2 * Cf[off])
        print(f"    alpha = {alpha}: Z^2 = {Z ** 2:.6f}; rel diff HOM {r1:.1e}, corrected {r2:.1e}; "
              f"mean I unfiltered {Hu[off].mean():.4f} vs filtered {Hf[off].mean():.4f}")
        ok &= r1 < 1e-12 and r2 < 1e-12 and np.all(np.diag(Hu) == 1.0)
        ok &= relerr(un["eta"] * Z, fi["eta"]) < 1e-12
        ok &= relerr(un["hom_phonon_factor"], Z ** 2) < 1e-12 and fi["hom_phonon_factor"] == 1.0
    return bool(ok)


@check("(f2) review 3: dephasing_mode 'enhanced' -> I = Gamma/(Gamma+gamma*) = 1/m independent of F_P (rel 1e-12; F3 numerics hom_pair rel 1e-9); 'dot' -> gamma* = (m-1) Gamma_0 fixed and I strictly increasing with F_P; modes agree at F_P = 0")
def _f2():
    card = _card()
    m = 1.7
    Fs = (0.0, 2.0, 5.0, 10.0)
    ok = True
    Id = []
    for F in Fs:
        se = rs.source_from_card(card, dephasing_mode="enhanced", F_P=F, hom_width_multiple=m)
        sd = rs.source_from_card(card, dephasing_mode="dot", F_P=F, hom_width_multiple=m)
        re_, rd = rs.source_rates(se), rs.source_rates(sd)
        Ie = rs.hom_pair(re_["Gamma"], re_["gamma_star"], 0.0)
        ok &= relerr(Ie, 1.0 / m) < 1e-9 and relerr(re_["Gamma"] / (re_["Gamma"] + re_["gamma_star"]), 1.0 / m) < 1e-12
        ok &= relerr(rd["gamma_star"], (m - 1.0) * sd.gamma0_ns) < 1e-12
        Id.append(rs.hom_pair(rd["Gamma"], rd["gamma_star"], 0.0))
        ok &= relerr(Id[-1], rd["Gamma"] / (rd["Gamma"] + rd["gamma_star"])) < 1e-9
        if F == 0.0:
            ok &= relerr(Id[-1], Ie) < 1e-12
    print(f"    I(F_P) enhanced = 1/m = {1.0 / m:.6f} for all F_P; dot = {[round(v, 5) for v in Id]}")
    ok &= all(b > a for a, b in zip(Id, Id[1:]))
    return bool(ok)


@check("(f3) review 2: TPE gamma_XX is the bare rate gamma_ratio x Gamma_0, so the cascade factor falls with F_P as 2/(3 + F_P L) from the closed-form L (rel 1e-9 at F_P = 0, 2, 5, 10)")
def _f3():
    card = _card()
    ok = True
    for F in (0.0, 2.0, 5.0, 10.0):
        t = rs.source_from_card(card, transition="X", excitation="tpe", pulse_fwhm_ps=20.0, F_P=F,
                                hom_width_multiple=1.0, kappa_ueV=300.0)
        pf = rs.bs_figure_of_merit(t, 4)["pulse_factor"]
        ref = 2.0 / (3.0 + F * _L_closed(t))
        print(f"    F_P = {F}: pulse factor {pf:.10f} vs 2/(3+F_P L) = {ref:.10f}")
        ok &= relerr(pf, ref) < 1e-9
    return bool(ok)


@check("(f4) review 4: the script verdict uses I - g2: summary/CSV I_cons == I_corr - g2 (rel 1e-12) and qualifying_rows recomputed from the CSV with I_corr - g2 >= 0.9 and g2 <= 0.05 matches the script; mean_pairwise_I_conservative == corrected - g2")
def _f4():
    import csv
    import json
    import subprocess
    import tempfile
    ok = True
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run([sys.executable, str(SCRIPT_PATH), "--quick", "--out-dir", td],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-800:])
            return False
        with open(Path(td) / "sweep.csv", newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        summ = json.loads((Path(td) / "summary.json").read_text(encoding="utf-8"))
    for N in (4, 8, 20):
        worst = max(relerr(float(x[f"I_cons_N{N}"]), float(x[f"I_mean_corr_N{N}"]) - float(x["g2"])) for x in rows)
        nq = sum(1 for x in rows if float(x["g2"]) <= 0.05 and float(x[f"I_mean_corr_N{N}"]) - float(x["g2"]) >= 0.9)
        nq_old = sum(1 for x in rows if float(x["g2"]) <= 0.05 and float(x[f"I_mean_corr_N{N}"]) >= 0.9)
        print(f"    N={N}: worst I_cons rel err {worst:.1e}; qualifying (I-g2) {nq} vs script {summ[f'N{N}']['qualifying_rows']} (uncorrected rule would give {nq_old})")
        ok &= worst < 1e-12 and nq == summ[f"N{N}"]["qualifying_rows"]
    out = rs.bs_figure_of_merit(rs.source_from_card(_card()), 4)
    ok &= abs(out["mean_pairwise_I_conservative"] - (out["mean_pairwise_I_corrected"] - out["g2"][0])) < 1e-15
    return bool(ok)


def main():
    t0 = time.time()
    n_pass = 0
    for name, fn in CHECKS:
        try:
            ok = bool(fn())
        except Exception as ex:          # a crash is a failure
            import traceback
            traceback.print_exc()
            print(f"    EXCEPTION {type(ex).__name__}: {ex}")
            ok = False
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
        n_pass += ok
    print(f"\n({time.time() - t0:.1f} s)")
    print(f"{n_pass}/{len(CHECKS)} verify_resonant_source checks passed")
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
