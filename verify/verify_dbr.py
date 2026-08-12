"""Phase-T3 regression suite: fsim_core.dbr (1-D waveguide/DBR tier).

Yuki standard: every closed form is checked against an independent second
method coded fresh in this file (Airy multiple-beam partial sum, direct
1-D Helmholtz RK4 integration, quarter-wave admittance closed forms,
numeric reflection-phase slopes) or against the quarantined S-1 EM oracle
(tests/oracles -- importable HERE because verify scripts sit outside the
oracle-only quarantine, which bans fsim_core/fsim_viz/fsim_gui only; see
the S-1 check in verify_fsim.py).

Run: python verify/verify_dbr.py   (exit code 0 iff all pass)
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.dbr import (
    HC_MEV_NM,
    Layer,
    cavity_mode,
    cavity_stack,
    dbr_reflectivity,
    dEdT_cav,
    energy_meV,
    invert_kappa,
    penetration_depth,
    planar_purcell,
    power_RT,
    quarter_wave_stack,
    qw_peak_reflectivity,
    rt_amplitudes,
    stopband_edges_analytic,
    transfer_matrix,
)
from tests.oracles.em_oracles import fresnel_slab  # verify-side use is allowed

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# Reference material system for the suite: GaAs/AlAs-class contrast at the
# 660 nm design point (values are generic H/L indices, not a material claim).
N_H, N_L, N_C, LAM0 = 3.5, 3.0, 3.5, 660.0
E0 = HC_MEV_NM / LAM0


# ------------------------------------------------- independent second methods

def airy_slab_sum(n1, n2, n3, lam, d, n_terms=400):
    """Fresh multiple-beam (Airy/Fabry-Perot) PARTIAL SUM for a single
    slab: the geometric series of multiple internal reflections summed
    term by term (no closed form, no transfer matrix) -- algebraically
    independent of both fsim_core.dbr and the S-1 oracle's closed-form
    geometric sum."""
    beta = 2.0 * math.pi * n2 * d / lam
    r12 = (n1 - n2) / (n1 + n2)
    r21 = -r12
    r23 = (n2 - n3) / (n2 + n3)
    t12, t21, t23 = 1.0 + r12, 1.0 + r21, 1.0 + r23
    ph = cmath.exp(2j * beta)
    q = r21 * r23 * ph
    r = complex(r12)
    t = 0j
    term_r = t12 * t21 * r23 * ph
    term_t = t12 * t23 * cmath.exp(1j * beta)
    qm = 1.0 + 0j
    for _ in range(n_terms):
        r += term_r * qm
        t += term_t * qm
        qm *= q
    return r, t, abs(r) ** 2, (n3 / n1) * abs(t) ** 2


def helmholtz_rt(layers, lam, n_in, n_out, steps_per_layer=300):
    """Direct numerical 1-D Helmholtz solve, E'' + k0^2 n(z)^2 E = 0:
    fixed-step RK4 integration of (E, E') backward from the exit medium
    (E = 1, E' = i k0 n_out at the back interface) to the front, then
    r = B/A, t = 1/A from the plane-wave decomposition at z = 0. No
    transfer matrices anywhere -- the T3 plan's 'direct numerical 1-D
    Helmholtz solve' second method."""
    k0 = 2.0 * math.pi / lam
    y = np.array([1.0 + 0j, 1j * k0 * n_out])
    for L in reversed(layers):
        ksq = (k0 * L.n) ** 2
        h = -L.d_nm / steps_per_layer

        def f(v):
            return np.array([v[1], -ksq * v[0]])

        for _ in range(steps_per_layer):
            k1 = f(y)
            k2 = f(y + 0.5 * h * k1)
            k3 = f(y + 0.5 * h * k2)
            k4 = f(y + h * k3)
            y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    kin = k0 * n_in
    E, dE = y
    A = 0.5 * (E + dE / (1j * kin))
    B = 0.5 * (E - dE / (1j * kin))
    return B / A, 1.0 / A


def mirror_phase_slope_Lpen(layers, n_c, n_sub, dlam=0.05):
    """Numeric phase-penetration depth from the reflection phase slope:
    L_pen = |dphi_r/dlambda| * lambda0^2 / (4 pi n_c)."""
    r_m, _ = rt_amplitudes(layers, LAM0 - dlam, n_in=n_c, n_out=n_sub)
    r_p, _ = rt_amplitudes(layers, LAM0 + dlam, n_in=n_c, n_out=n_sub)
    dphi = np.angle(r_p / r_m)  # phase difference, branch-safe for small dlam
    return abs(dphi) / (2.0 * dlam) * LAM0**2 / (4.0 * math.pi * n_c)


# ------------------------------------------------------------ 1: Fresnel

@check("single interface: r from transfer_matrix([]) equals Fresnel "
       "(n1-n2)/(n1+n2) exactly; R+T=1")
def _():
    for n1, n2 in ((1.0, 3.5), (3.5, 1.0), (3.0, 3.5), (2.0, 1.5), (1.0, 1.0)):
        M = transfer_matrix([], 660.0, n_in=n1, n_out=n2)
        r = M[1, 0] / M[0, 0]
        assert abs(r - (n1 - n2) / (n1 + n2)) < 1e-15, (n1, n2, r)
        R, T = power_RT([], 660.0, n1, n2)
        assert abs(R + T - 1.0) < 1e-15, (n1, n2, R + T)


# ------------------------------------------------------------ 2: Airy slab

@check("single slab: r,t,R,T match a fresh term-by-term Airy multiple-"
       "reflection sum to 1e-12; R+T=1 to 1e-12 (incl. a random 12-layer stack)")
def _():
    cases = [(1.0, 3.5, 3.0, 660.0, 80.0), (3.0, 3.5, 1.0, 700.0, 47.14),
             (1.0, 2.0, 1.5, 500.0, 300.0), (1.5, 1.0, 1.5, 800.0, 123.4)]
    for n1, n2, n3, lam, d in cases:
        r, t = rt_amplitudes([Layer(n2, d)], lam, n_in=n1, n_out=n3)
        R, T = power_RT([Layer(n2, d)], lam, n_in=n1, n_out=n3)
        r2, t2, R2, T2 = airy_slab_sum(n1, n2, n3, lam, d)
        assert abs(r - r2) < 1e-12 and abs(t - t2) < 1e-12, (n1, n2, n3)
        assert abs(R - R2) < 1e-12 and abs(T - T2) < 1e-12, (n1, n2, n3)
        assert abs(R + T - 1.0) < 1e-12, (R, T)
    rng = np.random.default_rng(0)
    stack = [Layer(float(n), float(d)) for n, d in
             zip(rng.uniform(1.4, 3.6, 12), rng.uniform(30.0, 300.0, 12))]
    R, T = power_RT(stack, 660.0, n_in=1.0, n_out=3.2)
    assert abs(R + T - 1.0) < 1e-12, R + T


# ------------------------------------------------------------ 3: S-1 oracle

@check("3-layer case: r,t,R,T agree with the quarantined S-1 oracle "
       "(tests.oracles.em_oracles.fresnel_slab) to 1e-12")
def _():
    for n1, n2, n3, lam, d in ((1.0, 3.5, 3.0, 660.0, 100.0),
                               (3.0, 3.5, 3.0, 700.0, 50.0),
                               (1.0, 1.5, 2.5, 532.0, 210.0)):
        o = fresnel_slab(n1, n2, n3, lam, d)
        r, t = rt_amplitudes([Layer(n2, d)], lam, n_in=n1, n_out=n3)
        R, T = power_RT([Layer(n2, d)], lam, n_in=n1, n_out=n3)
        assert abs(r - o["r"]) < 1e-12 and abs(t - o["t"]) < 1e-12, (n1, n2, n3)
        assert abs(R - o["R"]) < 1e-12 and abs(T - o["T"]) < 1e-12, (n1, n2, n3)


# ------------------------------------------------------- 4: Helmholtz ODE

@check("N-layer second method: 8-pair DBR r,t match a direct 1-D Helmholtz "
       "RK4 integration (no transfer matrices) to 1e-6 at three wavelengths")
def _():
    stack = quarter_wave_stack(N_H, N_L, 8, LAM0)
    for lam in (LAM0, 640.0, 703.0):
        r_tm, t_tm = rt_amplitudes(stack, lam, n_in=1.0, n_out=N_H)
        r_od, t_od = helmholtz_rt(stack, lam, 1.0, N_H)
        assert abs(r_tm - r_od) < 1e-6, (lam, r_tm, r_od)
        assert abs(t_tm - t_od) < 1e-6, (lam, t_tm, t_od)


# ----------------------------------------------------- 5: quarter-wave peak

@check("quarter-wave mirror: peak R at lambda0 matches the textbook closed "
       "form ((1-(n_out/n_in)(n_h/n_l)^2N)/(1+...))^2 for 5/10/20 pairs to 1e-10")
def _():
    for n_in, n_out in ((1.0, N_H), (N_H, 1.0), (1.0, 1.0)):
        for N in (5, 10, 20):
            stack = quarter_wave_stack(N_H, N_L, N, LAM0)
            R, _ = power_RT(stack, LAM0, n_in=n_in, n_out=n_out)
            R_ana = qw_peak_reflectivity(n_in, n_out, N_H, N_L, N)
            assert abs(R - R_ana) < 1e-10, (n_in, n_out, N, R, R_ana)


# --------------------------------------------------------------- 6: stopband

@check("stopband: R(lambda0) increases monotonically with pair count; "
       "numeric half-peak band is centered on lambda0 (frequency midpoint "
       "within 1%) and converges from above onto the analytic edges as N "
       "grows; |trace| of the period matrix = 2 exactly at the analytic edges")
def _():
    R_prev = -1.0
    for N in range(1, 13):
        stack = quarter_wave_stack(N_H, N_L, N, LAM0)
        R, _ = power_RT(stack, LAM0, n_in=1.0, n_out=N_H)
        assert R > R_prev, (N, R, R_prev)
        R_prev = R
    d = dbr_reflectivity(N_H, N_L, 15, LAM0, n_in=1.0, n_out=N_H)
    E_mid = 0.5 * (energy_meV(d["stopband_lo_nm"]) + energy_meV(d["stopband_hi_nm"]))
    assert abs(E_mid - E0) < 0.01 * E0, (E_mid, E0)
    assert abs(d["R_at_lambda0"] - qw_peak_reflectivity(1.0, N_H, N_H, N_L, 15)) < 1e-9
    assert d["R_peak"] >= d["R_at_lambda0"] > 0.98, (d["R_peak"], d["R_at_lambda0"])
    # second method for the analytic edges: at the infinite-crystal band edge
    # the Bloch condition |cos(K Lambda)| = |tr M_pair|/2 = 1 holds exactly
    lo_a, hi_a = stopband_edges_analytic(N_H, N_L, LAM0)
    pair = quarter_wave_stack(N_H, N_L, 1, LAM0)
    for lam_edge in (lo_a, hi_a):
        M = transfer_matrix(pair, lam_edge, n_in=N_H, n_out=N_H)
        assert abs(abs(np.trace(M)) - 2.0) < 1e-9, (lam_edge, np.trace(M))
    # the finite-N half-peak band is wider than the infinite-stack band and
    # converges onto it from above (fringe cliff sharpens ~1/N)
    w_ana = energy_meV(lo_a) - energy_meV(hi_a)
    ratios = []
    for N in (15, 30, 60):
        dN = dbr_reflectivity(N_H, N_L, N, LAM0, n_in=1.0, n_out=N_H,
                              n_points=8001)
        w_num = energy_meV(dN["stopband_lo_nm"]) - energy_meV(dN["stopband_hi_nm"])
        ratios.append(w_num / w_ana)
    assert ratios[0] > ratios[1] > ratios[2] > 1.0, ratios
    assert ratios[2] < 1.05, ratios


# ------------------------------------------------------ 7: penetration depth

@check("penetration depth: numeric phase slope dphi_r/dlambda of a 40-pair "
       "mirror reproduces the stated analytic convention within 2% for both "
       "the H-first (lambda0/4dn) and L-first (x n_h n_l/n_c^2) geometries")
def _():
    # H-first mirror seen from a low-index spacer n_c = n_l
    stack_h = quarter_wave_stack(N_H, N_L, 40, LAM0)
    L_num = mirror_phase_slope_Lpen(stack_h, n_c=N_L, n_sub=N_H)
    L_ana = penetration_depth(N_H, N_L, LAM0, first_layer="H")
    assert abs(L_num / L_ana - 1.0) < 0.02, (L_num, L_ana)
    # L-first mirror seen from a high-index spacer n_c = n_h
    stack_l = []
    for _ in range(40):
        stack_l.append(Layer(N_L, LAM0 / (4.0 * N_L)))
        stack_l.append(Layer(N_H, LAM0 / (4.0 * N_H)))
    L_num2 = mirror_phase_slope_Lpen(stack_l, n_c=N_H, n_sub=N_H)
    L_ana2 = penetration_depth(N_H, N_L, LAM0, n_c=N_C, first_layer="L")
    assert abs(L_num2 / L_ana2 - 1.0) < 0.02, (L_num2, L_ana2)
    # and the two conventions really differ (the factor is not cosmetic)
    assert abs(L_ana2 / L_ana - (N_H * N_L) / N_C**2) < 1e-12


# --------------------------------------------------------- 8: lambda cavity

@check("lambda cavity (symmetric 8/8, lossless): resonance within a linewidth "
       "of the design energy, on-resonance transmission -> 1, and Q = E/kappa "
       "self-consistent")
def _():
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, 8, 8, LAM0)
    m = cavity_mode(top + [spacer] + bottom, LAM0, n_in=1.0, n_out=1.0)
    assert abs(m["E_cav_meV"] - E0) < m["kappa_meV"], (m["E_cav_meV"], E0)
    assert m["T_peak"] > 0.999, m["T_peak"]  # symmetric + lossless -> unity
    assert abs(m["Q"] - m["E_cav_meV"] / m["kappa_meV"]) < 1e-9 * m["Q"], m
    assert m["kappa_meV"] > 0.0 and np.isfinite(m["kappa_meV"]), m


@check("lambda cavity: kappa decreases (Q increases) monotonically as top "
       "pairs increase at fixed bottom-mirror margin")
def _():
    prev_kappa, prev_Q = np.inf, 0.0
    for n_top in (4, 6, 8, 10):
        top, spacer, bottom = cavity_stack(N_H, N_L, N_C, n_top, n_top + 4, LAM0)
        m = cavity_mode(top + [spacer] + bottom, LAM0, n_in=1.0, n_out=N_H)
        assert m["kappa_meV"] < prev_kappa, (n_top, m["kappa_meV"], prev_kappa)
        assert m["Q"] > prev_Q, (n_top, m["Q"], prev_Q)
        prev_kappa, prev_Q = m["kappa_meV"], m["Q"]


# -------------------------------------------------------- 9: planar Purcell

@check("planar_purcell: both mirrors removed -> F = 1 exactly; dipole at the "
       "antinode of a symmetric high-R cavity -> F >> 1; at a node -> F < 1")
def _():
    d_c = LAM0 / N_C
    # no mirrors: half-'stacks' are bare n_c|n_c interfaces, r = 0 exactly
    F_free = planar_purcell([], [], N_C, d_c, 0.3 * d_c, LAM0,
                            n_in=N_C, n_out=N_C)
    assert abs(F_free - 1.0) < 1e-14, F_free
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, 10, 10, LAM0)
    F_anti = planar_purcell(top, bottom, N_C, spacer.d_nm, 0.5 * d_c, LAM0)
    F_node = planar_purcell(top, bottom, N_C, spacer.d_nm, 0.25 * d_c, LAM0)
    assert F_anti > 10.0, F_anti
    assert F_node < 1.0, F_node


@check("planar_purcell high-finesse limit: F(antinode) agrees with the "
       "standard symmetric-cavity estimate (1+|r|)/(1-|r|) ~ 4/(1-R) within 10%")
def _():
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, 10, 10, LAM0)
    F_anti = planar_purcell(top, bottom, N_C, spacer.d_nm,
                            0.5 * spacer.d_nm, LAM0)
    r_mag = abs(rt_amplitudes(list(reversed(top)), LAM0, n_in=N_C, n_out=1.0)[0])
    F_est = (1.0 + r_mag) / (1.0 - r_mag)
    assert abs(F_anti / F_est - 1.0) < 0.10, (F_anti, F_est)
    assert abs(F_est * (1.0 - r_mag**2) / (1.0 + r_mag) ** 2 - 1.0) < 1e-12


# ------------------------------------------------------------ 10: inversion

@check("invert_kappa: returned pair count meets the target under a direct "
       "cavity_mode recompute, and one fewer top pair violates it")
def _():
    target = 0.5  # meV
    sol = invert_kappa(target, N_H, N_L, N_C, LAM0, n_bottom_extra=6)
    assert sol["n_top"] >= 2, sol
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, sol["n_top"],
                                       sol["n_bottom"], LAM0)
    m = cavity_mode(top + [spacer] + bottom, LAM0, n_in=1.0, n_out=N_H)
    assert m["kappa_meV"] <= target, (m["kappa_meV"], target)
    assert abs(m["kappa_meV"] - sol["kappa_meV"]) < 1e-9 * target, (m, sol)
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, sol["n_top"] - 1,
                                       sol["n_top"] - 1 + 6, LAM0)
    m_less = cavity_mode(top + [spacer] + bottom, LAM0, n_in=1.0, n_out=N_H)
    assert m_less["kappa_meV"] > target, (m_less["kappa_meV"], target)
    assert sol["L_pen_nm"] == penetration_depth(N_H, N_L, LAM0, n_c=N_C,
                                                first_layer="L")


# ----------------------------------------------------------------- 11: dEdT

@check("dEdT_cav: uniform dn/dT < 0 gives dE_cav/dT > 0 with magnitude "
       "matching the analytic estimate dE/E = -dn/n within a factor of 2")
def _():
    dndt = -2.0e-4  # 1/K, uniform across all layers (sign per the T3 brief)
    top, spacer, bottom = cavity_stack(N_H, N_L, N_C, 8, 8, LAM0,
                                       dn_dT_h=dndt, dn_dT_l=dndt, dn_dT_c=dndt)
    slope = dEdT_cav(top + [spacer] + bottom, LAM0, n_in=1.0, n_out=1.0)
    assert slope > 0.0, slope
    est = -E0 * dndt / N_C  # dE/E = -dn/n with n ~ n_c (mode-weighted)
    assert 0.5 < slope / est < 2.0, (slope, est)


# -------------------------------------------------------------- 12: hygiene

@check("three-layer/oracle hygiene: fsim_core/dbr.py has no GUI imports and "
       "never references the quarantined oracle package")
def _():
    src = (Path(__file__).resolve().parents[1] / "fsim_core" / "dbr.py").read_text(
        encoding="utf-8")
    for banned in ("tkinter", "PyQt", "PySide", "wx", "matplotlib",
                   "streamlit", "plotly", "dearpygui"):
        assert banned not in src, banned
    for quarantined in ("tests" + ".oracles", "tests" + "/oracles"):
        assert quarantined not in src, quarantined


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
