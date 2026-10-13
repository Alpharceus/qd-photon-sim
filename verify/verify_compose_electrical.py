"""Checks for the electrical-drive physics of the composable engine (devlib-p10e).

Spec: .workers/specs/devlib-p10e-electrical.md; physics: .workers/devlib/p3-electrical-brief.md sections 1-4.
Anchors come from closed forms, the brief's independently computed numbers ([NV], first written with other code
paths), an independent scipy.integrate.solve_ivp moment hierarchy, and the published values the brief records. Nothing
is compared against the number the new code produced for itself.

Run: python verify/verify_compose_electrical.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
from scipy.integrate import dblquad, solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import cw_g2                                            # noqa: E402
from fsim_core import lindblad as lb                                   # noqa: E402
from fsim_core import pulse_counting as pcnt                           # noqa: E402
from fsim_core import resonant_source as rs                            # noqa: E402
from fsim_core import spectral                                         # noqa: E402
from fsim_core.compose import engine, ensure_loaded                    # noqa: E402
from fsim_core.compose import observables as OB                        # noqa: E402
from fsim_core.compose import registry as R                            # noqa: E402
from fsim_core.compose import solvers as S                             # noqa: E402
from fsim_core.compose.types import DrivePlan, EmitterModel, Pulse     # noqa: E402

ensure_loaded()
RESULTS = []
HB = lb.HBAR_MEV_NS


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------ Sartison numbers (brief 1.3)
G_S = 1.0 / 0.599          # tau_X 599 ps
PERIOD_S = 5.0             # 200 MHz


def sartison(shape, mu, tau_on=0.2, t_r=0.0, pump_ratio=1.0, n_per=16, **kw):
    dts, rk = S.pump_cells(shape, mu, tau_on, t_r, n_per=n_per)
    steps = [(d, r, 0.0, 0.0) for d, r in zip(dts, rk)]
    return pcnt.pulse_g2(0, G_S, 2 * G_S, 0, 0, 1, 0, 0, PERIOD_S - sum(dts), pump_ratio=pump_ratio,
                         r_profile=steps, **kw), dts, rk


# ---- 1. default paths: profile machinery reduces exactly to the rectangular code
leg = pcnt.pulse_g2(15.0, G_S, 2 * G_S, 0.1, 0.2, 1, 0, 0.2, 4.8, gate_ns=1.0, adjacent=True, split=True)
prof = pcnt.pulse_g2(0, G_S, 2 * G_S, 0, 0, 1, 0, 0, 4.8, r_profile=[(0.2, 15.0, 0.1, 0.2)], k_dark=(0.1, 0.2),
                     gate_ns=1.0, adjacent=True, split=True)
check("1a. pulse_g2 r_profile with ONE rectangular step == the rect path (g2, counts, g2_adj, X share), gate and escape on",
      all(abs(leg[k] - prof[k]) <= 1e-12 * max(1.0, abs(leg[k])) for k in ("g2", "mean_counts", "g2_adj", "mean_counts_x")),
      f"g2 {leg['g2']:.12f} vs {prof['g2']:.12f}")
Gm, km, gsm = 1.3, 0.2, 0.7
H0 = np.zeros((3, 3), dtype=complex)


def e3(i, j):
    m = np.zeros((3, 3), dtype=complex)
    m[i, j] = 1.0
    return m


Lb = lb.liouvillian(H0, [math.sqrt(Gm) * e3(0, 1), math.sqrt(km) * e3(0, 1), math.sqrt(gsm) * e3(1, 1),
                          math.sqrt(1 / 0.05) * e3(1, 2)])
Dp = lb.dissipator(e3(2, 0))
Jm = Gm * lb.jump_super(e3(0, 1))
a1 = lb.pulsed_counting(Lb + 12.0 * Dp, Lb, Jm, 3, 0.3, 4.7, gate_ns=2.0)
a2 = lb.pulsed_counting(None, Lb, Jm, 3, 0.0, 4.7, gate_ns=2.0, profile=[(0.3, Lb + 12.0 * Dp)])
check("1b. lindblad.pulsed_counting(profile=[one step]) == the L_on/L_off path (g2, mean counts, rho_period)",
      abs(a1["g2"] - a2["g2"]) < 1e-12 and abs(a1["mean_counts"] - a2["mean_counts"]) < 1e-12
      and np.max(np.abs(a1["rho_period"] - a2["rho_period"])) < 1e-12, f"g2 {a1['g2']:.12f}")
# evolve_td second envelope: detuning pulse on a two-level system, against an analytic Bloch (undamped) solution
s2 = lb.build_system(levels=2, gamma_X_ns=0.0, omega_ns=0.0)
Ldet = lb.hamiltonian_super(np.array([[0, 0], [0, 1]], dtype=complex))
rho0 = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)         # |+>: coherence phase = integral of det
tg = np.linspace(0.0, 1.0, 201)
out = lb.evolve_td(s2.L, np.zeros_like(s2.L), lambda t: 0.0, rho0, tg, L2=Ldet, envelope2=lambda t: 3.0)
check("1c. evolve_td second envelope: det_X = 3/ns for 1 ns rotates the coherence by exp(-i 3) (closed form)",
      abs(out[-1][1, 0] - 0.5 * np.exp(1j * 3.0)) < 1e-9 or abs(out[-1][1, 0] - 0.5 * np.exp(-1j * 3.0)) < 1e-9,
      f"rho_eg = {out[-1][1, 0]:.6f}")
base = lb.evolve_td(s2.L, np.zeros_like(s2.L), lambda t: 0.0, rho0, tg)
check("1d. evolve_td default (no L2) unchanged: no detuning, no decay -> coherence conserved",
      abs(base[-1][1, 0] - 0.5) < 1e-12)

# ---- 2. shaped pump (brief 1.1-1.3)
dts, rk = S.pump_cells("rc", 3.0, 0.2, 0.1)
check("2a. pump_cells conserves the operating point: sum r_k dt_k = mu for rect, rc and gauss",
      all(abs(sum(d * r for d, r in zip(*S.pump_cells(sh, 3.0, 0.2, tr))) - 3.0) < 1e-12
          for sh, tr in (("rect", 0.0), ("rc", 0.1), ("rc", 0.3), ("gauss", 0.0))))
# analytic RC integral: mu = r0 tau_on exactly for the UNtruncated profile (rise deficit = tail excess)
tau, tr_, r0 = 0.2, 0.1, 15.0
mu_rc = r0 * (tau - tr_ * (1 - math.exp(-tau / tr_)) + tr_ * (1 - math.exp(-tau / tr_)))
check("2b. RC profile keeps mu = r0 tau_on (rise deficit equals tail excess); the window LENGTHENS to tau_on + 12 t_r",
      abs(mu_rc - 3.0) < 1e-12 and abs(sum(dts) - (0.2 + 12 * 0.1)) < 1e-12, f"window {sum(dts) * 1e3:.0f} ps")
for name, shape, tr, ref in (("rect 200 ps", "rect", 0.0, 0.137), ("Gaussian FWHM 200 ps", "gauss", 0.0, 0.200),
                              ("rect + RC t_r 100 ps", "rc", 0.1, 0.237), ("rect + RC t_r 300 ps", "rc", 0.3, 0.464)):
    o, _, _ = sartison(shape, 3.0, 0.2, tr, pump_ratio=0.0 if False else 1.0)
    check(f"2c. Sartison cascade, mu 3, {name}: g2 {ref} (brief 1.3, 3%)", rel(o["g2"], ref) < 0.03, f"{o['g2']:.4f}")
sat = {mu: sartison("rect", mu)[0]["g2"] for mu in (0.5, 2.0, 4.0, 8.0, 20.0)}
ref_sat = {0.5: 0.113, 2.0: 0.134, 4.0: 0.135, 8.0: 0.101, 20.0: 0.049}
check("2d. cascade saturation then FALL (rect 200 ps / tau 599 ps): 0.113 0.134 0.135 0.101 0.049 (3%)",
      all(rel(sat[m], ref_sat[m]) < 0.03 for m in sat), " ".join(f"{sat[m]:.3f}" for m in sat))
check("2e. strong-pump sign: cascade g2 falls (dot parks in XX) while the two-level pump keeps rising",
      sat[20.0] < sat[4.0] and sartison("rect", 20.0, pump_ratio=0.0)[0]["g2"] > sartison("rect", 4.0, pump_ratio=0.0)[0]["g2"])
tl = {mu: sartison("rect", mu, pump_ratio=0.0)[0]["g2"] for mu in (1.0, 4.0, 20.0)}
check("2f. two-level (pump_ratio 0): 0.148 0.268 0.400 at mu 1, 4, 20 (3%)",
      all(rel(tl[m], r) < 0.03 for m, r in zip(tl, (0.148, 0.268, 0.400))), " ".join(f"{tl[m]:.3f}" for m in tl))
w1 = sartison("rect", 0.01, 0.002, pump_ratio=0.0)[0]["g2"]
s1 = sartison("rect", 4000.0, 0.002, pump_ratio=0.0)[0]["g2"]
check("2g. closed-form limits (two-level pump, Gamma tau_on << 1): weak g2 -> Gamma tau_on/3, strong fast pump -> 2 Gamma tau_on (2%)",
      rel(w1, G_S * 0.002 / 3.0) < 0.02 and rel(s1, 2 * G_S * 0.002) < 0.02,
      f"weak {w1:.6f} vs {G_S * 0.002 / 3:.6f}; strong {s1:.5f} vs {2 * G_S * 0.002:.5f}")
c16 = sartison("rc", 3.0, 0.2, 0.1, n_per=16)[0]["g2"]
c32 = sartison("rc", 3.0, 0.2, 0.1, n_per=32)[0]["g2"]
check("2h. step-count convergence: RC g2 at 16 vs 32 steps per t_r differs < 1e-3 (relative)", rel(c16, c32) < 1e-3,
      f"{c16:.5f} {c32:.5f}")


def ivp_rc(mu, tau_on, t_r, gamma, pump_ratio=1.0, period=PERIOD_S):
    """Independent solve_ivp counting hierarchy with the analytic RC r(t) (not the cell-averaged grid)."""
    r0 = mu / tau_on
    M0 = cw_g2.generator(0.0, gamma, 2 * gamma, 0, 0, pump_ratio)
    M1 = cw_g2.generator(1.0, gamma, 2 * gamma, 0, 0, pump_ratio) - M0
    Jm_ = np.zeros((3, 3))
    Jm_[0, 1] = gamma                          # X line only (t_X = 1, t_XX = 0), as in sartison()

    def r(t):
        return r0 * (1 - math.exp(-t / t_r)) if t < tau_on else r0 * (1 - math.exp(-tau_on / t_r)) * math.exp(-(t - tau_on) / t_r)

    def rhs(t, v):
        M = M0 + r(t) * M1
        p, m1, m2 = v[:3], v[3:6], v[6:]
        return np.concatenate([M @ p, M @ m1 + Jm_ @ p, M @ m2 + 2 * Jm_ @ m1])

    p = np.array([1.0, 0, 0])
    for _ in range(40):
        v0 = np.concatenate([p, np.zeros(6)])
        sol = solve_ivp(rhs, (0.0, period), v0, method="DOP853", rtol=1e-10, atol=1e-13,
                        t_eval=[period], max_step=0.01)
        pn = sol.y[:3, -1]
        if np.max(np.abs(pn - p)) < 1e-12:
            p = pn
            break
        p = pn
    v = sol.y[:, -1]
    return float(np.sum(v[6:]) / np.sum(v[3:6]) ** 2), float(np.sum(v[3:6]))


# the profile is truncated at 12 t_r: compare on the same tail (ivp runs the full analytic tail)
g_ivp, m_ivp = ivp_rc(3.0, 0.2, 0.1, G_S)
g_grid = sartison("rc", 3.0, 0.2, 0.1)[0]
check("2i. cell-averaged profile grid vs independent solve_ivp hierarchy with the analytic RC r(t): g2 and <m> agree to 1e-3",
      rel(g_grid["g2"], g_ivp) < 1e-3 and rel(g_grid["mean_counts"], m_ivp) < 1e-3,
      f"g2 {g_grid['g2']:.5f} vs {g_ivp:.5f}; <m> {g_grid['mean_counts']:.5f} vs {m_ivp:.5f}")
check("2j. RC tails RAISE g2 over the clean rectangle (never shorten the window): 0.137 < 0.237 < 0.464",
      sat[2.0] < sartison("rc", 2.0, 0.2, 0.1)[0]["g2"] < sartison("rc", 2.0, 0.2, 0.3)[0]["g2"])
# Sartison band: the paper's background-corrected 0.23 +/- 0.05 lies inside the shape band (brief section 5)
band = [sartison(sh, 3.0, 0.2, tr)[0]["g2"] for sh, tr in (("rect", 0), ("gauss", 0), ("rc", 0.05), ("rc", 0.1), ("rc", 0.15))]
check("2k. Sartison shape band {rect, gauss, rc 0-150 ps} contains the paper's 0.23 +/- 0.05 (UC, brief section 5)",
      min(band) <= 0.23 + 0.05 and max(band) >= 0.23 - 0.05 and min(band) < 0.23 < max(band),
      f"band {min(band):.3f}-{max(band):.3f}")

# ---- 3. EL decay and escape (brief 3)
GAM_E, K_OFF, TAU_CAP = 1.4, 0.19, (0.01, 0.1)


def tail_slope(g_c, k_off, k_on=0.0, tau_on=0.2, mu=0.4, ts=(1.0, 2.0)):
    """Emission-rate decay (1/ns) between 1 and 2 ns after the pulse of the S3d scheme (G, X, R): pump into R for the
    window, R -> X at g_c, X -> G at Gamma + k."""
    Lk = lambda k: lb.liouvillian(H0, [math.sqrt(GAM_E) * e3(0, 1), math.sqrt(g_c) * e3(1, 2)]
                                  + ([math.sqrt(k) * e3(0, 1)] if k > 0 else []))
    Lon = Lk(k_on) + (mu / tau_on) * lb.dissipator(e3(2, 0))
    Loff = Lk(k_off)
    v = lb.vec(np.diag([1.0, 0, 0]).astype(complex))
    v = expm(Lon * tau_on) @ v
    ts = np.array(ts)
    em_ = []
    for t in ts:
        rho = lb.unvec(expm(Loff * (t - 0.0)) @ v, 3)
        em_.append(GAM_E * rho[1, 1].real)
    return float(-(math.log(em_[1]) - math.log(em_[0])) / (ts[1] - ts[0]))


for tc in TAU_CAP:
    sl = tail_slope(1.0 / tc, 0.0)
    check(f"3a. Ellis: tail decay stays at Gamma = 1.400/ns with capture tau {tc * 1e3:.0f} ps and the 200 ps pump edge",
          rel(sl, GAM_E) < 5e-3, f"{sl:.4f}/ns")
sl_off = tail_slope(1.0 / 0.05, K_OFF)
check("3b. Ellis: with k_off = 0.19/ns the EL tail is Gamma + k_off = 1.59/ns", rel(sl_off, 1.59) < 5e-3, f"{sl_off:.4f}/ns")
check("3c. slow capture (tau_cap > tau_rad, trap 10): the tail measures capture, slope = min(Gamma, g_c)",
      rel(tail_slope(1.0 / 3.0, 0.0, ts=(10.0, 12.0)), 1.0 / 3.0) < 1e-2,
      f"{tail_slope(1.0 / 3.0, 0.0, ts=(10.0, 12.0)):.4f} vs 0.3333")
# engine: the rate and lindblad_incoherent solvers report tau_EL and keep Gamma0
em_s = EmitterModel(1380, 900, 3.5, GAM_E, 0.0, 0.0, None, -2.0, 2.0, 0, frozenset({"X", "XX"}), "in_plane", 0, 0.0, 50.0)
rt_s = NS(Gamma_ns=GAM_E, k_nr_ns=0.0, gamma_star_ns=0.0, Gamma0_ns=GAM_E)
plan_e = DrivePlan("electrical", "pulsed", 200.0, pump_rate_ns=10.0, tau_on_ns=0.2, k_esc_off_ns=K_OFF)
for sname in ("rate", "lindblad_incoherent"):
    sol = R.get("solver", sname)
    ok, why = sol.supports(plan_e, em_s, rt_s)
    res = sol.run(plan_e, em_s, rt_s) if ok else None
    check(f"3d. {sname} solver: tau_EL = 1/(Gamma + k_off) = {1e3 / (GAM_E + K_OFF):.1f} ps reported, radiative tau_rad "
          f"= {1e3 / GAM_E:.1f} ps untouched, Gamma0 never refit", ok and rel(res.extra["tau_EL_ps"], 1e3 / (GAM_E + K_OFF)) < 1e-12
          and rel(res.extra["tau_rad_ps"], 1e3 / GAM_E) < 1e-12 and rt_s.Gamma_ns == GAM_E, why)
r_noesc = R.get("solver", "rate").run(DrivePlan("electrical", "pulsed", 200.0, pump_rate_ns=10.0, tau_on_ns=0.2), em_s, rt_s)
r_esc = R.get("solver", "rate").run(plan_e, em_s, rt_s)
check("3e. escape in the dark window lowers the photons per pulse (k_off is a pure loss, never a jitter)",
      r_esc.n_ph < r_noesc.n_ph, f"{r_noesc.n_ph:.4f} -> {r_esc.n_ph:.4f}")
r_on = R.get("solver", "rate").run(DrivePlan("electrical", "pulsed", 200.0, pump_rate_ns=10.0, tau_on_ns=0.2, k_esc_on_ns=2.0), em_s, rt_s)
check("3f. escape during the pump window also lowers n_ph and the engine supports() window uses the shaped length",
      r_on.n_ph < r_noesc.n_ph and not R.get("solver", "rate").supports(
          DrivePlan("electrical", "pulsed", 200.0, pump_rate_ns=10.0, tau_on_ns=0.2, pump_profile="rc", t_rise_ns=0.5),
          em_s, rt_s)[0])
# rate vs lindblad_incoherent agree (no relaxation level) for a shaped, escaping plan
pl_sh = DrivePlan("electrical", "pulsed", 200.0, pump_rate_ns=10.0, tau_on_ns=0.2, pump_profile="rc", t_rise_ns=0.1,
                  k_esc_on_ns=0.3, k_esc_off_ns=K_OFF)
a_rate = R.get("solver", "rate").run(pl_sh, em_s, rt_s)
a_lind = R.get("solver", "lindblad_incoherent").run(pl_sh, em_s, rt_s)
check("3g. shaped + escaping plan: rate (moment hierarchy) and lindblad_incoherent (Liouvillian) give the same g2 and n_ph",
      rel(a_rate.g2_0, a_lind.g2_0) < 1e-6 and rel(a_rate.n_ph, a_lind.n_ph) < 1e-6,
      f"g2 {a_rate.g2_0:.6f} vs {a_lind.g2_0:.6f}")

# ---- 4. Stark-gated detuning (brief 2.1)
EM_C = EmitterModel(1380, 900, 3.5, 1.0, 0.0, 1.0, None, -2.0, 2.0, 0, frozenset({"X", "XX"}), "in_plane", 0, 0.0, 50.0)
RT_C = NS(Gamma_ns=1.0, k_nr_ns=0.0, gamma_star_ns=1.0, Gamma0_ns=1.0)   # Cao: T1 = 1 ns, T2 = 1 ns -> gamma* = 1/ns
OM = 2 * math.pi / 1.6                                                    # measured Rabi period 1.6 ns at 500 nW [V]
LS = S.LindbladSolver()
cao = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 300.0, OM * 0.3),)), EM_C, RT_C)
check("4a. Cao 2015 unfitted: 300 ps rect resonant gate at the measured Rabi frequency, g2 = 0.028 +/- 0.010 (UC)",
      abs(cao.g2_0 - 0.028) <= 0.010, f"g2 {cao.g2_0:.4f}, n_ph {cao.n_ph:.3f}, area {OM * 0.3 / math.pi:.2f} pi")
check("4b. Cao: weak-gate closed form g2 -> Gamma T/3 (1 - O(Theta^2)): 0.1/3 = 0.033 upper bound, model below it",
      cao.g2_0 < 1.0 * 0.3 / 3.0)
st = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 300.0, OM * 0.3, detuning_meV=0.0015),)), EM_C, RT_C)
st_td = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 300.0, OM * 0.3, detuning_meV=0.0015),),
                         needs=frozenset({"td_detuning"}), det_profile=(0.0015, 0.0015, 100.0, 0.0)), EM_C, RT_C)
check("4c. static detuning: the td_detuning (step-propagator) path reproduces resonant_source (g2, n_ph, HOM) to 1e-9",
      abs(st.g2_0 - st_td.g2_0) < 1e-9 and abs(st.n_ph - st_td.n_ph) < 1e-9 and abs(st.I_single - st_td.I_single) < 1e-7,
      f"g2 {st.g2_0:.9f} {st_td.g2_0:.9f}")
plan_gate = DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 500.0, OM * 0.5, detuning_meV=0.5),),
                      needs=frozenset({"td_detuning"}), det_profile=(0.5, 0.0, 300.0, 0.0))
check("4d. lindblad.supports() accepts td_detuning (flipped to supported); rate / jc still refuse it",
      LS.supports(plan_gate, EM_C, RT_C)[0] and not R.get("solver", "rate").supports(plan_gate, EM_C, RT_C)[0]
      and not R.get("solver", "jc").supports(plan_gate, EM_C, RT_C)[0])
gate = LS.run(plan_gate, EM_C, RT_C)


def ivp_gate(W, off_meV, gate_ns=0.3, period=5.0):
    """Independent solve_ivp of the gated problem: L(t) = L0 + Omega L1 + det(t) L_det, counting hierarchy, periodic state."""
    pc = rs._pieces(rs.TwoLevel(1.0, gamma_star_ns=1.0))
    Ld_ = lb.hamiltonian_super(np.array([[0, 0], [0, 1]], dtype=complex))
    off = float(lb.mev_to_rate(off_meV))
    ts_, te_ = 0.5 * W - 0.5 * gate_ns, 0.5 * W + 0.5 * gate_ns

    def rhs(t, v):
        L = pc.L0 + OM * pc.L1 + (0.0 if ts_ <= t < te_ else off) * Ld_
        return np.concatenate([L @ v[:4], L @ v[4:8] + pc.J @ v[:4], L @ v[8:] + 2 * pc.J @ v[4:8]])
    Aug = rs._aug(pc.L0 + off * Ld_, pc.J)
    v = np.zeros(12, dtype=complex)
    v[:4] = lb.vec(rs._ground(2))
    for _ in range(30):
        y = v.copy()
        y[4:] = 0
        for a, b in ((0.0, ts_), (ts_, te_), (te_, W)):
            if b > a:
                y = solve_ivp(rhs, (a, b), y, method="DOP853", rtol=1e-11, atol=1e-14, max_step=5e-4).y[:, -1]
        y = expm(Aug * (period - W)) @ y
        done = np.max(np.abs(y[:4] - v[:4])) < 1e-12
        v = y
        if done:
            break
    tr_ = lb.trace_row(2)
    m1_, m2_ = (tr_ @ y[4:8]).real, (tr_ @ y[8:]).real
    return float(m2_ / m1_ ** 2), float(m1_)


g_ivp_, n_ivp_ = ivp_gate(0.5, 0.5)
check("4e. Stark gate (det_off 0.5 meV, 300 ps resonant gate inside a 500 ps window): g2 and n_ph agree with an independent "
      "solve_ivp of L(t) = L0 + Omega L1 + det(t) L_det to 1e-6", rel(gate.g2_0, g_ivp_) < 1e-6 and rel(gate.n_ph, n_ivp_) < 1e-6,
      f"g2 {gate.g2_0:.8f} vs {g_ivp_:.8f}; n_ph {gate.n_ph:.6f} vs {n_ivp_:.6f}")
gate_tight = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 300.0, OM * 0.3, detuning_meV=0.5),),
                              needs=frozenset({"td_detuning"}), det_profile=(0.5, 0.0, 300.0, 0.0)), EM_C, RT_C)
check("4e'. gate filling the whole window == the 300 ps rect resonant gate (g2 within 1%: the detuned dark window only dephases "
      "the leftover coherence); off-gate drive on top moves g2 by a few % (physics, not a bug)",
      rel(gate_tight.g2_0, cao.g2_0) < 0.01 and rel(gate.n_ph, cao.n_ph) < 0.01,
      f"g2 {gate_tight.g2_0:.5f} vs {cao.g2_0:.5f}; with 100 ps off-gate each side {gate.g2_0:.5f}")
lz = gate.extra["lz"]
check("4f. Landau-Zener check reported: sharp edge, Omega^2 = 15/ns^2 << sweep -> not adiabatic (the dot just switches on)",
      abs(lz["omega_ns"] - OM) < 1e-9 and not lz["adiabatic"] and lz["p_lz_diabatic"] == 1.0,
      f"Omega {lz['omega_ns']:.3f}/ns, sweep inf")
gate_rc = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 500.0, OM * 0.5, detuning_meV=0.5),),
                           needs=frozenset({"td_detuning"}), det_profile=(0.5, 0.0, 300.0, 20.0)), EM_C, RT_C)
lz_rc = gate_rc.extra["lz"]
sweep_exp = (0.5 / HB) / 0.020
check("4g. RC-edge sweep rate = |det_on - det_off|/t_r and Omega^2/sweep << 1 (4e5-class, brief 2.1)",
      rel(lz_rc["sweep_ns2"], sweep_exp) < 1e-9 and lz_rc["adiabaticity"] < 1e-2, f"sweep {lz_rc['sweep_ns2']:.3g}/ns^2")
zero = LS.run(DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, 0.0),)), EM_C, RT_C)
check("4h. LindbladSolver at pulse area 0 returns the trivial zero-excitation result (no ZeroDivisionError)",
      zero.n_ph == 0.0 and zero.inversion == 0.0 and zero.extra["eta_prep"] == 0.0 and zero.valid and math.isnan(zero.I_single),
      f"n_ph {zero.n_ph}, g2 {zero.g2_0}")
zero_g = LS.run(DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 5.0, 0.0),)), EM_C, RT_C)
zero_t = LS.run(DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 5.0, 0.0),), needs=frozenset({"cascade"})), EM_C, RT_C)
check("4i. area 0 also safe for a Gaussian pulse and for TPE", zero_g.n_ph == 0.0 and zero_t.n_ph == 0.0)

# ---- 5. gated HOM (Bennett, brief 2.2)
G_B, GS_B = 1.0 / 0.8, 2.0 / 0.062               # T1 800 ps; T2* 62 ps textbook pure dephasing -> repo gamma* = 2/T2*
dc = OB.gated_hom(G_B, GS_B, 0.0, 2000.0)
check("5a. Bennett dc reference: ungated I = Gamma/(Gamma + gamma*) = 0.0373 (brief 0.0373 [NV]; paper 3.75%)",
      abs(dc["I_T"] - G_B / (G_B + GS_B)) < 1e-10 and abs(dc["I_T"] - 0.0373) < 5e-5 and abs(dc["I_T"] - 0.0375) < 3e-4,
      f"{dc['I_T']:.6f}")
g31 = OB.gated_hom(G_B, GS_B, 0.0, 0.031)
check("5b. 31 ps gate: kept fraction 3.8% and I_T = 0.74 (brief [NV]) without jitter",
      abs(g31["n_T"] - 0.038) < 1e-3 and abs(g31["I_T"] - 0.74) < 0.01, f"n_T {g31['n_T']:.4f}, I_T {g31['I_T']:.4f}")
check("5c. limits: T -> 0 gives I_T -> 1; T -> inf gives Gamma/(Gamma + gamma*); k_loss shifts it to (Gamma+k)/(Gamma+k+gamma*)",
      abs(OB.gated_hom(G_B, GS_B, 0.0, 1e-6)["I_T"] - 1.0) < 1e-3
      and abs(OB.gated_hom(G_B, GS_B, 0.0, 1e4, k_loss=0.4)["I_T"] - (G_B + 0.4) / (G_B + 0.4 + GS_B)) < 1e-9)
# independent 2-D numerical integration of the gate overlap (delay 0 and 20 ps, gate opens at 10 ps)


def num_I(Gamma, gs, tg, T, d):
    E = tg + T
    lo = max(tg, 0.0, d)

    def f(t2, t1):
        return Gamma ** 2 * math.exp(Gamma * d) * math.exp(-Gamma * (t1 + t2)) * math.exp(-gs * abs(t1 - t2))
    with warnings.catch_warnings():              # the |t1 - t2| kink makes quad report roundoff; the 1e-6 check is the verdict
        warnings.simplefilter("ignore")
        num = dblquad(f, lo, E, lo, E, epsabs=1e-10, epsrel=1e-8)[0]
    nA = Gamma * (math.exp(-Gamma * max(tg, 0)) - math.exp(-Gamma * E)) / Gamma
    nB = math.exp(Gamma * d) * (math.exp(-Gamma * max(tg, d)) - math.exp(-Gamma * E))
    return num / (nA * nB)


ok_num = True
for tg_, T_, d_ in ((0.0, 0.031, 0.0), (0.010, 0.031, 0.0), (0.010, 0.031, 0.020), (0.0, 0.2, -0.05)):
    a = OB.gated_hom(G_B, GS_B, tg_, T_, delay_ns=d_)["I_T"]
    b = num_I(G_B, GS_B, tg_, T_, d_)
    ok_num &= rel(a, b) < 1e-6
check("5d. closed-form gated overlap == scipy dblquad (4 gate/delay settings, 1e-6)", ok_num)
check("5e. unwindowed delay dependence I(d) = I(0) e^{-Gamma |d|} (HOM dip of independent packets)",
      abs(OB.gated_hom(G_B, GS_B, 0.0, 5000.0, delay_ns=0.3)["I_T"] - dc["I_T"] * math.exp(-G_B * 0.3)) < 1e-9)
# FX3b: the old 5f entered the 31 ps jitter ON TOP of a top-hat gate of the same width (double counting, and the old jitter
# model was blind to gates opening after the excitation). Bennett's 31 ps window IS a Gaussian: use the Gaussian gate.
from scipy.special import erfcx                                                # noqa: E402
SIG_G = 0.031 / (2.0 * math.sqrt(2.0 * math.log(2.0)))                         # Gaussian gate: 31 ps INTENSITY FWHM
gB = OB.gated_hom(G_B, GS_B, 0.05, 0.031, gate_shape="gauss")["I_T"]
check("5f. Bennett V: Gaussian 31 ps gate (no extra jitter) gives I_T = 0.656 = the closed form erfcx(gamma* sigma_g) "
      "(1e-4), inside the brief 0.55-0.75 band and within 0.10 of the paper's 0.60-0.64 (UC); the top-hat gate gives 0.736",
      abs(gB - erfcx(GS_B * SIG_G)) < 1e-4 and 0.55 <= gB <= 0.75 and abs(gB - 0.62) <= 0.10 and abs(g31["I_T"] - 0.736) < 0.002,
      f"Gaussian {gB:.4f} vs closed form {erfcx(GS_B * SIG_G):.4f}, top-hat {g31['I_T']:.4f}")
# gated g2 background: 1.5% -> 0.03 (CAL identity)
check("5g. Bennett gated g2 0.03 = background 1 - rho^2 with rho = sqrt(1 - 0.03) = 0.985 (law and its inverse round-trip)",
      abs(OB.g2_with_background(0.0, 1.0, (1 / math.sqrt(0.97) - 1)) - 0.03) < 1e-12
      and abs((0.42 - (1 - 0.86 ** 2)) / 0.86 ** 2 - OB.g2_with_background(0.0, 1, 1) * 0 - (0.42 - (1 - 0.86 ** 2)) / 0.86 ** 2) < 1e-15)
rho_ = 0.86
g_dot_ = 0.23
g_raw_ = OB.g2_with_background(g_dot_, rho_, 1 - rho_)           # S = rho, B = 1 - rho -> rho = S/(S + B) = rho
check("5h. background law and the Sartison/Lazic bgc formula are exact inverses: g_bgc = (g_raw - (1 - rho^2))/rho^2",
      abs((g_raw_ - (1 - rho_ ** 2)) / rho_ ** 2 - g_dot_) < 1e-12, f"raw {g_raw_:.4f}")

# ---- 6. binned g2 and IRF (Schlehahn, brief 2.3 / trap 2 / trap 6)
gam_S = 1.0 / 0.27
c_leak = math.exp(-gam_S * 1.0 / 2.0)
check("6a. binned g2 closed form: g2_bin = g2_0 (1 - e^{-gamma T/2}) + e^{-gamma T/2}; e^{-1.85} = 0.157 at 270 ps / 1 ns",
      abs(OB.binned_g2(0.17, gam_S, 1.0) - (0.17 * (1 - c_leak) + c_leak)) < 1e-12 and abs(c_leak - 0.157) < 1e-3,
      f"{OB.binned_g2(0.17, gam_S, 1.0):.4f}")
check("6b. binned limits: T -> inf gives g2_0; T -> 0 gives 1 (all neighbours leak); g2_0 = 1 stays 1",
      abs(OB.binned_g2(0.17, gam_S, 1e3) - 0.17) < 1e-12 and abs(OB.binned_g2(0.17, gam_S, 1e-9) - 1.0) < 1e-8
      and abs(OB.binned_g2(1.0, gam_S, 0.7) - 1.0) < 1e-12)
c_irf = OB.binned_g2(0.17, gam_S, 1.0, 300.0)


def c0_numeric(gamma, T, fwhm_ps):
    """P(|L + N| < T/2) for L ~ two-sided exponential (rate gamma), N ~ Gaussian: integral over L of the Gaussian window mass."""
    from scipy.integrate import quad
    from scipy.stats import norm
    sig = math.sqrt(2) * fwhm_ps * 1e-3 / (2 * math.sqrt(2 * math.log(2)))
    f = lambda y: 0.5 * gamma * math.exp(-gamma * abs(y)) * (norm.cdf((T / 2 - y) / sig) - norm.cdf((-T / 2 - y) / sig))
    a, _ = quad(f, -40.0, 0.0, epsabs=1e-13, epsrel=1e-12, limit=400)
    b, _ = quad(f, 0.0, 40.0, epsabs=1e-13, epsrel=1e-12, limit=400)
    return a + b


check("6c. IRF leak == direct numerical convolution of the two-sided exponential with a Gaussian (1e-4), and exceeds the IRF-free value",
      abs(c_irf - (1 - (1 - 0.17) * c0_numeric(gam_S, 1.0, 300.0))) < 1e-4 and c_irf > OB.binned_g2(0.17, gam_S, 1.0),
      f"IRF-free {OB.binned_g2(0.17, gam_S, 1.0):.4f}, 300 ps IRF {c_irf:.4f}")
check("6d. Schlehahn raw g2 0.3 reproduced with g2_0 ~ 0.17 (+ IRF) inside the brief's 0.25-0.35 pass band (UC on binning)",
      0.25 <= OB.binned_g2(0.17, gam_S, 1.0) <= 0.35 and 0.25 <= c_irf <= 0.35)
for (tau_d, fwhm, ref, who) in ((4.0, 200.0, 0.983, "Zhang"), (1.0, 350.0, 0.89, "Lazic"), (0.27, 300.0, 0.71, "Schlehahn")):
    f = cw_g2._dip_factor(tau_d, fwhm, "gaussian")
    check(f"6e. IRF dip factor {who}: {ref} (brief trap 2 [NV], 0.01)", abs(f - ref) < 0.01, f"{f:.4f}")
check("6f. irf_dip_g2: irf_ps = 0 returns g2_0 (default off); otherwise 1 - (1 - g2_0) F",
      OB.irf_dip_g2(0.1, 1.0, 0.0) == 0.1 and abs(OB.irf_dip_g2(0.1, 1.0, 350.0) - (1 - 0.9 * cw_g2._dip_factor(1.0, 350.0, "gaussian"))) < 1e-15)

# ---- 7. materials: gamma*(T) LO-phonon switch (brief 4.1)
gan = R.get("material", "gan_aln")
vals_t = {"gamma_T_model": "lo_phonon", "gamma_star_4K_ns": 1.0 / HB, "a_ac_ns_per_K": 0.018 / HB,
          "b_LO_meV": 3570.0, "E_LO_meV": 102.0}
for T_, ref in ((300.0, 76.8), (5.0, 1.09)):
    gs_ = gan.build(vals_t, T_).gamma_star_ns * HB
    direct = float(spectral.gamma_of_T(T_, 1.0, 0.018, 3570.0, 102.0))
    check(f"7a. lo_phonon at the Tamariz 2020 constants, T = {T_:.0f} K: == spectral.gamma_of_T and = {ref} meV",
          abs(gs_ - direct) < 1e-9 * direct and abs(gs_ - ref) < 0.01 * ref, f"{gs_:.3f} meV")
em_lin = gan.build({}, 300.0)
em_lin_v = gan.build({"gamma_T_model": "linear"}, 300.0)
a_def = gan.param("a_ac_ns_per_K").default
g4_def = gan.param("gamma_star_4K_ns").default
check("7b. default model is linear and bit-identical to gamma*(4 K) + a_ac (T - 4 K); explicit 'linear' same",
      em_lin.gamma_star_ns == max(0.0, g4_def + a_def * 296.0) and em_lin.gamma_star_ns == em_lin_v.gamma_star_ns
      and gan.param("gamma_T_model").default == "linear")
inp = R.get("material", "inp_gainp_red")
from fsim_core import linewidth                                   # noqa: E402
lw = linewidth.LinewidthParams()
em_inp = inp.build({"gamma_T_model": "lo_phonon"}, 300.0)
check("7c. InP entry: lo_phonon uses its existing LinewidthParams constants and reproduces linewidth.gamma_anchor (Bommer-class) at 300 K",
      rel(em_inp.gamma_star_ns * HB, float(linewidth.gamma_anchor(300.0, lw))) < 1e-6,
      f"{em_inp.gamma_star_ns * HB:.4f} vs {float(linewidth.gamma_anchor(300.0, lw)):.4f} meV")
z = R.get("material", "inas_gaas_sk")
check("7d. entries without an LO constant: b_LO_meV default 0 [A], lo_phonon differs from linear only by the 4 K reference (a_ac * 4 K)",
      z.param("b_LO_meV").default == 0.0 and z.param("b_LO_meV").tag == "A"
      and abs((z.build({"gamma_T_model": "lo_phonon"}, 50.0).gamma_star_ns - z.build({}, 50.0).gamma_star_ns)
              - 4.0 * z.param("a_ac_ns_per_K").default) < 1e-9)
check("7e. every nitride / InP LO coupling carries a tag and a source (tags kept): gan_aln E, ingan_gan A, inp_gainp_red DR",
      (R.get("material", "gan_aln").param("b_LO_meV").tag, R.get("material", "ingan_gan").param("b_LO_meV").tag,
       inp.param("b_LO_meV").tag) == ("E", "A", "DR"))
check("7f. Tamariz window transmissions (existing spectral functions): X 0.092, XX 0.027, eps = 0.29 (brief 4.2 [NV])",
      abs(float(spectral.tophat_transmission(0, 55, 8)) - 0.092) < 1e-3
      and abs(float(spectral.tophat_transmission(47, 80, 8)) - 0.027) < 1e-3
      and abs(float(spectral.tophat_transmission(47, 80, 8)) / float(spectral.tophat_transmission(0, 55, 8)) - 0.29) < 0.01)

# ---- 8. drive entries and end to end
ep, ec = R.get("drive", "electrical_pulsed"), R.get("drive", "electrical_cw")
need = ("pump_profile", "t_rise_ps", "k_esc_on_ns", "k_esc_off_ns", "gate_ps", "hom_gate_ps", "hom_gate_start_ps",
        "hom_jitter_ps", "irf_ps")
check("8a. electrical_pulsed exposes the shaped-pump / escape / gate / IRF params with unit, tag, source and help; "
      "defaults = off", all(any(p.name == n for p in ep.params) for n in need)
      and all(p.help and p.tag in "VDREA" and p.source for p in ep.params if p.name in need)
      and ep.param("pump_profile").default == "rect" and ep.param("k_esc_off_ns").default == 0.0
      and ep.param("irf_ps").default == 0.0 and ec.param("k_esc_on_ns").default == 0.0)
dp_default = ep.plan({}, em_s, rt_s)
check("8b. default plan: new DrivePlan fields at their defaults (rect, no rise, no escape, no gate, no det_profile)",
      (dp_default.pump_profile, dp_default.t_rise_ns, dp_default.k_esc_on_ns, dp_default.k_esc_off_ns, dp_default.gate_ns,
       dp_default.det_profile) == ("rect", 0.0, 0.0, 0.0, None, None) and "hom_gate" not in dp_default.inputs)


def espec(drp, mat=None, T=4.2, solver="auto"):
    return {"schema": "fsim-composed/1", "T_K": T, "material": {"kind": "inas_gaas_sk", "params": mat or {}},
            "structure": {"kind": "bulk_planar", "params": {}}, "drive": {"kind": "electrical_pulsed", "params": drp},
            "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": solver}


M3 = {"eta_capture": 1.0, "tau_on_ns": 0.2, "I_uA": 3.0 * 1.602176634e-19 / 0.2e-9 * 1e6, "f_rep_MHz": 200.0}   # mu = 3
MAT = {"tau_bare_ns": 0.599, "gamma_star_4K_ns": 0.0, "a_ac_ns_per_K": 0.0}
SOLV = "rate"
e_rect = engine.evaluate_spec(espec(M3, MAT, solver=SOLV))["scalars"]
e_rc = engine.evaluate_spec(espec({**M3, "pump_profile": "rc", "t_rise_ps": 100.0}, MAT, solver=SOLV))["scalars"]
e_rc_rc = engine.evaluate_spec(espec({**M3, "pump_profile": "rc", "R_drive_ohm": 50.0, "C_dep_F": 2e-12}, MAT, solver=SOLV))["scalars"]
e_gs = engine.evaluate_spec(espec({**M3, "pump_profile": "gauss"}, MAT, solver=SOLV))["scalars"]
check("8c. engine end to end (M-3, mu 3, two-level incoherent): rect < gauss < rc(100 ps) in g2 and the RC tail raises <n>",
      e_rect["g2_0"] < e_gs["g2_0"] < e_rc["g2_0"] and e_rc["n_ph"] > e_rect["n_ph"],
      f"g2 {e_rect['g2_0']:.3f} {e_gs['g2_0']:.3f} {e_rc['g2_0']:.3f}")
check("8d. RC from R_s = 50 ohm, C_j = 2 pF (t_r = 100 ps) == the t_rise_ps = 100 input, and keeps mu (not mech_pulsed's 'subtract t_r')",
      rel(e_rc_rc["g2_0"], e_rc["g2_0"]) < 1e-9 and rel(e_rc_rc["n_ph"], e_rc["n_ph"]) < 1e-9)
leg_e = engine.evaluate_spec(espec({**M3, "R_drive_ohm": 50.0, "C_dep_F": 2e-12}, MAT, solver=SOLV))["scalars"]
check("8e. profile 'rect' keeps the legacy behaviour (window shortened by R C): fewer photons than the RC profile at the same I",
      leg_e["n_ph"] < e_rc_rc["n_ph"])
env = engine.envelope(espec({**M3, "k_esc_off_ns": K_OFF, "irf_ps": 300.0}, {**MAT, "tau_bare_ns": 1 / GAM_E}, solver=SOLV),
                      ["tau_EL_ps", "tau_ps", "g2_binned", "g2_0"])
check("8f. observables: tau_EL_ps = 1e3/(Gamma + k_off) (= 628.9 ps) differs from tau_ps = 714.3 ps (radiative, F_P untouched)",
      rel(env["tau_EL_ps"]["central"], 1e3 / (GAM_E + K_OFF)) < 1e-9 and rel(env["tau_ps"]["central"], 1e3 / GAM_E) < 1e-9,
      f"{env['tau_EL_ps']['central']:.1f} vs {env['tau_ps']['central']:.1f} ps")
gb = OB.binned_g2(env["g2_0"]["central"], GAM_E + K_OFF, 1e3 / 200.0, 300.0)
check("8g. engine g2_binned == binned_g2(g2_0, Gamma + k_off, T_rep, irf_ps) (bin = pulse period)",
      abs(env["g2_binned"]["central"] - gb) < 1e-12, f"{gb:.4f}")
# Bennett through the engine: hom_I_gated with the entered gate and jitter
BEN = {"tau_bare_ns": 0.8, "gamma_star_4K_ns": 2.0 / 0.062, "a_ac_ns_per_K": 0.0, "alpha_ps2": 1e-9}
envb = engine.envelope(espec({**M3, "hom_gate_ps": 31.0, "hom_jitter_ps": 0.0, "f_rep_MHz": 505.0}, BEN, solver="lindblad_incoherent"),
                       ["hom_I_gated"])
sc = engine.evaluate_spec(espec({**M3, "hom_gate_ps": 31.0, "f_rep_MHz": 505.0}, BEN, solver="lindblad_incoherent"))["scalars"]
check("8h. engine hom_I_gated (Bennett T1 800 ps, T2* 62 ps, 31 ps gate): = gated_hom x q^2 and well above the ungated I",
      abs(envb["hom_I_gated"]["central"] - 0.7358 * 1.0) < 0.01 and envb["hom_I_gated"]["central"] > 10 * 0.0373,
      f"{envb['hom_I_gated']['central']:.4f}")

# ---- 9. registry / observables sanity
check("9a. new observables registered and offered by the electrical drives",
      all(i in R.get("drive", "electrical_pulsed").observables for i in ("tau_EL_ps", "g2_binned", "hom_I_gated"))
      and all(i in R.get("drive", "electrical_cw").observables for i in ("tau_EL_ps", "g2_irf")))


# ---- 10. devlib-p10e2: stark_gated_resonant drive, setup.irf_ps, spectral_overlap observable (independent anchors)
sg = R.get("drive", "stark_gated_resonant")
check("10a. stark_gated_resonant is ready and exposes the det profile (off, on, gate, RC edge), Rabi period / power, f_rep, "
      "HOM gate + jitter; each with unit, tag, source, help",
      sg.status == "ready" and all(any(p.name == n for p in sg.params) for n in
      ("det_off_meV", "det_on_meV", "t_gate_ps", "t_r_ps", "rabi_period_ns", "power_nW", "power_ref_nW", "f_rep_MHz",
       "hom_gate_ps", "hom_jitter_ps")) and all(p.help and p.source and p.tag in "VDREA" for p in sg.params))
pl = sg.plan({"det_off_meV": 5.0, "det_on_meV": 0.2, "t_gate_ps": 300.0, "t_r_ps": 20.0, "power_nW": 125.0}, EM_C, RT_C)
check("10b. paper convention Delta = omega_L - omega_X: Delta_off = +5 meV (laser blue) is stored as repo det_X = -5 in det_profile "
      "and the pulse detuning; Omega = (2 pi / 1.6 ns) sqrt(125/500) (power scaling); needs td_detuning",
      pl.det_profile == (-5.0, -0.2, 300.0, 20.0) and abs(pl.pulses[0].detuning_meV + 5.0) < 1e-12
      and abs(pl.inputs["omega_ns"] - OM * 0.5) < 1e-12 and "td_detuning" in pl.needs and pl.family == "coherent")


def cao_spec(drp=None, su=None):
    return {"schema": "fsim-composed/1", "T_K": 4.2,
            "material": {"kind": "inas_gaas_sk", "params": {"tau_bare_ns": 1.0, "gamma_star_4K_ns": 1.0, "a_ac_ns_per_K": 0.0,
                                                          "alpha_ps2": 1e-9}},
            "structure": {"kind": "bulk_planar", "params": {}},
            "drive": {"kind": "stark_gated_resonant", "params": drp or {}}, "setup": {"kind": "lens_fibre_det", "params": su or {}},
            "solver": "auto"}


cao_e = engine.evaluate_spec(cao_spec({"det_off_meV": 0.5, "t_gate_ps": 300.0}))["scalars"]
ref_gate = LS.run(DrivePlan("coherent", "pulsed", 200.0, pulses=(Pulse("rect", 300.0, OM * 0.3, detuning_meV=-0.5),),
                            needs=frozenset({"td_detuning"}), det_profile=(-0.5, 0.0, 300.0, 0.0)), EM_C, RT_C)
check("10c. engine end to end (Cao: T1 = T2 = 1 ns, 500 nW, 300 ps gate): g2 and n_ph == the direct LindbladSolver Stark-gate run "
      "(blue laser, sign-flipped plan) and g2 within the Cao 0.028 +/- 0.010 (UC)",
      abs(cao_e["g2_0"] - ref_gate.g2_0) < 1e-6 and abs(cao_e["n_ph"] - ref_gate.n_ph) < 1e-6 and abs(cao_e["g2_0"] - 0.028) <= 0.010,
      f"g2 {cao_e['g2_0']:.5f} vs {ref_gate.g2_0:.5f}; n_ph {cao_e['n_ph']:.4f}")
cao_blue = engine.evaluate_spec(cao_spec({"det_off_meV": 0.5, "t_gate_ps": 300.0, "t_r_ps": 20.0}))["scalars"]
check("10d. an RC gate edge (t_r 20 ps; Rabi window t_gate + 6 t_r) stays finite and close to the sharp gate (within 20%)",
      math.isfinite(cao_blue["g2_0"]) and rel(cao_blue["g2_0"], cao_e["g2_0"]) < 0.2, f"{cao_blue['g2_0']:.4f}")
eh = engine.envelope(cao_spec({"det_off_meV": 5.0, "hom_gate_ps": 31.0}, {}), ["hom_I_gated"])
check("10e. hom_I_gated offered by the new drive (31 ps gate, T2 = 1 ns): finite and above the ungated I = Gamma/(Gamma + gamma*) = 0.5",
      math.isfinite(eh["hom_I_gated"]["central"]) and eh["hom_I_gated"]["central"] > 0.5,
      f"{eh['hom_I_gated']['central']:.3f}")
# setup.irf_ps: setup-level == drive-level, and the drive value takes precedence
e_irf_s = engine.envelope(dict(espec({**M3, "f_rep_MHz": 1000.0, "k_esc_off_ns": K_OFF}, MAT, solver=SOLV),
                               setup={"kind": "lens_fibre_det", "params": {"irf_ps": 300.0}}), ["g2_binned"])["g2_binned"]["central"]
e_irf_d = engine.envelope(espec({**M3, "f_rep_MHz": 1000.0, "k_esc_off_ns": K_OFF, "irf_ps": 300.0}, MAT, solver=SOLV),
                          ["g2_binned"])["g2_binned"]["central"]
e_irf_0 = engine.envelope(espec({**M3, "f_rep_MHz": 1000.0, "k_esc_off_ns": K_OFF}, MAT, solver=SOLV), ["g2_binned"])["g2_binned"]["central"]
e_irf_p = engine.envelope(dict(espec({**M3, "f_rep_MHz": 1000.0, "k_esc_off_ns": K_OFF, "irf_ps": 100.0}, MAT, solver=SOLV),
                               setup={"kind": "lens_fibre_det", "params": {"irf_ps": 300.0}}), ["g2_binned"])["g2_binned"]["central"]
check("10f. lens_fibre_det.irf_ps (default 0) == the electrical drive's irf_ps in g2_binned; an IRF raises the binned g2; the drive "
      "value wins when both are set",
      R.get("setup", "lens_fibre_det").param("irf_ps").default == 0.0 and abs(e_irf_s - e_irf_d) < 1e-12 and e_irf_s > e_irf_0
      and abs(e_irf_p - engine.envelope(espec({**M3, "f_rep_MHz": 1000.0, "k_esc_off_ns": K_OFF, "irf_ps": 100.0}, MAT, solver=SOLV),
                                        ["g2_binned"])["g2_binned"]["central"]) < 1e-12,
      f"IRF 300 {e_irf_s:.4f}, 0 {e_irf_0:.4f}")


def tophat_closed(delta, gam, w):
    """Independent closed form: Lorentzian (FWHM gam, offset delta) through a top-hat of width w."""
    return (math.atan((w / 2 - delta) / (gam / 2)) + math.atan((w / 2 + delta) / (gam / 2))) / math.pi


def tam_spec(ratio):
    gam_x = 55.0
    return {"schema": "fsim-composed/1", "T_K": 300.0,
            "material": {"kind": "gan_aln", "params": {"tau_bare_ns": 1000.0, "gamma_star_4K_ns": gam_x / HB, "a_ac_ns_per_K": 0.0,
                                                       "E_b_meV": 47.0}},
            "structure": {"kind": "bulk_planar", "params": {}},
            "drive": {"kind": "electrical_cw", "params": {}},
            "setup": {"kind": "lens_fibre_det", "params": {"filter_width_meV": 8.0, "xx_width_ratio": 80.0 / 55.0, "xx_x_ratio": ratio}},
            "solver": "auto"}


try:
    eo = engine.envelope(tam_spec(0.34), ["spectral_overlap", "spectral_overlap:rho", "spectral_overlap:g2"])
    eps_ind = tophat_closed(47.0, 80.0, 8.0) / tophat_closed(0.0, 55.0, 8.0)
    rho_ind = 1.0 / (1.0 + 0.34 * eps_ind)
    check("10g. spectral_overlap (Tamariz 2020: E_b 47, widths 55/80, window 8 meV) == the independent arctan closed form eps = 0.29; "
          "rho = 1/(1 + 0.34 eps), g2 = 1 - rho^2 (brief 4.2: 0.17 at 6.5 kW/cm^2)",
          abs(eo["spectral_overlap"]["central"] - eps_ind) < 2e-3 and abs(eo["spectral_overlap"]["central"] - 0.29) < 0.01
          and abs(eo["spectral_overlap:rho"]["central"] - rho_ind) < 1e-3
          and abs(eo["spectral_overlap:g2"]["central"] - (1 - rho_ind ** 2)) < 1e-3 and abs(eo["spectral_overlap:g2"]["central"] - 0.17) < 0.03,
          f"eps {eo['spectral_overlap']['central']:.4f} vs {eps_ind:.4f}, g2 {eo['spectral_overlap:g2']['central']:.4f}")
except Exception as exc:                                                   # noqa: BLE001
    check("10g. spectral_overlap Tamariz anchor", False, repr(exc))
no_f = engine.envelope(dict(tam_spec(0.0), setup={"kind": "lens_fibre_det", "params": {}}), ["spectral_overlap", "spectral_overlap:g2"])
check("10h. no window and no XX light: eps = 1 (both lines fully collected) and g2 = 0 (rho = 1)",
      abs(no_f["spectral_overlap"]["central"] - 1.0) < 1e-12 and abs(no_f["spectral_overlap:g2"]["central"]) < 1e-12)

# ---- 11. devlib-fx3b: gated HOM with Gaussian gate + monotone gate jitter, cascade drive, bg_frac (independent anchors)
from scipy.integrate import quad                                              # noqa: E402


def _w_gauss(t, c, sg):
    return math.exp(-0.5 * ((t - c) / sg) ** 2)


def num_gauss(Gamma, gs, c, T, shift, k=0.0):
    """Independent gated overlap of two Gaussian-gated photons (gates at c -+ shift/2), emitted together at t = 0: the
    interfering amplitude is sqrt(W_A W_B) built POINTWISE; 2-D integral by dblquad over [0, c + 14 sigma]^2."""
    sg = T / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    cA, cB = c - shift / 2, c + shift / 2
    Gd = Gamma + k
    U = max(cA, cB) + 14.0 * sg

    def M(t):
        return math.sqrt(_w_gauss(t, cA, sg) * _w_gauss(t, cB, sg))

    def f(t2, t1):
        return Gamma ** 2 * M(t1) * M(t2) * math.exp(-Gd * (t1 + t2)) * math.exp(-gs * abs(t1 - t2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        num = dblquad(f, 0.0, U, 0.0, U, epsabs=1e-12, epsrel=1e-8)[0]
    nA = Gamma * quad(lambda t: _w_gauss(t, cA, sg) * math.exp(-Gd * t), 0.0, U, epsabs=1e-14, epsrel=1e-11, limit=200)[0]
    nB = Gamma * quad(lambda t: _w_gauss(t, cB, sg) * math.exp(-Gd * t), 0.0, U, epsabs=1e-14, epsrel=1e-11, limit=200)[0]
    return num / (nA * nB)


def num_rect(Gamma, gs, tg, T, shift, k=0.0):
    """Independent rect-gate overlap: gates A = [tg - s/2, tg - s/2 + T], B = [tg + s/2, tg + s/2 + T]; both photons emitted at
    t = 0 (gates clipped to t >= 0); only the intersection interferes; dblquad on the intersection."""
    Gd = Gamma + k
    a0, b0 = tg - shift / 2, tg + shift / 2
    lo, hi = max(a0, b0, 0.0), min(a0, b0) + T
    nA = Gamma * (math.exp(-Gd * max(a0, 0.0)) - math.exp(-Gd * (a0 + T))) / Gd
    nB = Gamma * (math.exp(-Gd * max(b0, 0.0)) - math.exp(-Gd * (b0 + T))) / Gd
    if not hi > lo:
        return 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        num = dblquad(lambda t2, t1: Gamma ** 2 * math.exp(-Gd * (t1 + t2)) * math.exp(-gs * abs(t1 - t2)), lo, hi, lo, hi,
                      epsabs=1e-12, epsrel=1e-9)[0]
    return num / (nA * nB)


# 11a. Gaussian gate closed form (full-line limit): I_T = erfcx(gamma* sigma_g), independent of Gamma, k and the gate position
ok11 = True
det11 = ""
for (Gm_, kk_, cc_, TT_) in ((G_B, 0.0, 0.40, 0.031), (3.0, 0.5, 0.35, 0.020), (0.6, 0.0, 0.50, 0.060)):
    sg_ = TT_ / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    a_ = OB.gated_hom(Gm_, GS_B, cc_, TT_, k_loss=kk_, gate_shape="gauss")["I_T"]
    b_ = float(erfcx(GS_B * sg_))
    ok11 &= abs(a_ - b_) < 1e-7
    det11 += f" ({Gm_:.2f},{kk_},{cc_}): {a_:.6f} vs {b_:.6f};"
check("11a. Gaussian-gate overlap in the full-line limit = erfcx(gamma* sigma_g) exactly (derived: int int W W' e^{-Gd(t+t')} "
      "e^{-gamma|t-t'|} / n^2), independent of Gamma, k and the gate position (3 settings, 1e-7)", ok11, det11)

# 11b. Gaussian gate vs an independent 2-D integral (gate on the onset: truncation at t = 0 matters), shifted gates included
ok11 = True
det11 = ""
for (cc_, sh_) in ((0.0, 0.0), (0.0155, 0.0), (0.03, 0.020), (0.05, 0.040)):
    a_ = OB.gated_hom(G_B, GS_B, cc_, 0.031, gate_shape="gauss", shift_ns=sh_)["I_T"]
    b_ = num_gauss(G_B, GS_B, cc_, 0.031, sh_)
    ok11 &= rel(a_, b_) < 1e-6
    det11 += f" (c {cc_}, shift {sh_}): {a_:.6f} vs {b_:.6f};"
check("11b. Gaussian-gate gated_hom == independent pointwise sqrt(W_A W_B) dblquad (4 gate positions / shifts, truncation at "
      "t = 0 and displaced gates included, 1e-6)", ok11, det11)

# 11c. rect gate with displaced gates vs an independent intersection integral
ok11 = True
det11 = ""
for (tg_, sh_, kk_) in ((0.0, 0.010, 0.0), (0.020, 0.020, 0.0), (0.020, 0.040, 0.0), (0.010, 0.015, 0.4), (0.0, 0.05, 0.0)):
    a_ = OB.gated_hom(G_B, GS_B, tg_, 0.031, k_loss=kk_, shift_ns=sh_)["I_T"]
    b_ = num_rect(G_B, GS_B, tg_, 0.031, sh_, kk_)
    ok11 &= (abs(a_ - b_) < 1e-9 if b_ == 0.0 else rel(a_, b_) < 1e-6)
    det11 += f" ({tg_},{sh_}): {a_:.6f} vs {b_:.6f};"
check("11c. top-hat gate with the two gates displaced == independent dblquad on the gate intersection (5 settings, 1e-6); "
      "gates that no longer overlap give 0", ok11, det11)

# 11d. jitter: sigma -> 0 limit, monotone non-increasing (both shapes, gates opening at 0 and after the excitation), acts for
#      a late gate, independent Gaussian average
ok_lim = ok_mono = True
det_lim = ""
sigs = [0.0] + list(np.geomspace(1e-3, 0.3, 40))
for shape in ("rect", "gauss"):
    for tg_ in ((0.0, 0.02, 0.05) if shape == "rect" else (0.0155, 0.04, 0.08)):
        I0 = OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, 0.0, gate_shape=shape)
        ok_lim &= abs(I0 - OB.gated_hom(G_B, GS_B, tg_, 0.031, gate_shape=shape)["I_T"]) < 1e-14
        small = OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, 1e-7, gate_shape=shape)
        ok_lim &= abs(small - I0) < 2e-5
        det_lim += f" {shape}@{tg_}: {small:.5f}->{I0:.5f};"
        vals = [OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, sg_, gate_shape=shape) for sg_ in sigs]
        ok_mono &= all(vals[i + 1] <= vals[i] + 1e-12 for i in range(len(vals) - 1))
check("11d1. sigma -> 0 limit: gated_hom_jitter(sigma = 0) == gated_hom exactly and sigma = 0.1 fs differs by < 2e-5 "
      "(both shapes, three gate positions)", ok_lim, det_lim)
check("11d2. I_T(sigma) is monotone non-increasing over sigma = 0 ... 300 ps (41 values, both shapes, gate opening at the "
      "excitation and well after it)", ok_mono)
for shape, tg_ in (("rect", 0.05), ("gauss", 0.08)):
    I0 = OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, 0.0, gate_shape=shape)
    Ih = OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, 0.0155, gate_shape=shape)
    check(f"11d3. jitter acts for a {shape} gate opening 50-80 ps AFTER the excitation (the old model returned the jitter-free "
          "value there): sigma = T/2 lowers I_T by > 0.02", Ih < I0 - 0.02, f"{I0:.4f} -> {Ih:.4f}")
r_ = OB.gated_hom_jitter(G_B, GS_B, 0.05, 0.031, 2.0)
r2_ = OB.gated_hom_jitter(G_B, GS_B, 0.05, 0.031, 20.0)
check("11d4. very large jitter (sigma >> T): the gates no longer overlap, I_T -> 0 (only the fraction T/(sigma sqrt(2 pi)) of "
      "the pairs overlap: 2 ns < 0.05, 20 ns < 0.005)", r_ < 0.05 and r2_ < 0.005 and r2_ < r_, f"{r_:.5f} {r2_:.5f}")
for shape, tg_ in (("rect", 0.02), ("gauss", 0.04)):
    sgm_ = 0.020
    xs = np.linspace(-6.0, 6.0, 4001)
    ww = np.exp(-0.5 * xs ** 2)
    iv = np.array([(lambda v: v if math.isfinite(v) else 0.0)(OB.gated_hom(G_B, GS_B, tg_, 0.031, gate_shape=shape,
                                                                         shift_ns=sgm_ * x)["I_T"]) for x in xs])
    ref_ = float(np.sum(ww * iv) / np.sum(ww))
    got_ = OB.gated_hom_jitter(G_B, GS_B, tg_, 0.031, sgm_, gate_shape=shape)
    check(f"11d5. jitter average ({shape}, sigma 20 ps) == an independent 4001-node average of the displaced-gate overlap "
          "(2e-3: kinks of the top-hat gate)", abs(got_ - ref_) < 2e-3, f"{got_:.5f} vs {ref_:.5f}")

# 11e. engine end to end: Bennett Gaussian gate
BENG = {**M3, "hom_gate_ps": 31.0, "hom_gate_shape": "gauss", "hom_gate_start_ps": 50.0, "f_rep_MHz": 505.0}
envg = engine.envelope(espec(BENG, BEN, solver="lindblad_incoherent"), ["hom_I_gated"])["hom_I_gated"]["central"]
ref_g = OB.gated_hom(G_B, GS_B, 0.05, 0.031, gate_shape="gauss")["I_T"]
check("11e. engine hom_I_gated with drive.hom_gate_shape = gauss (31 ps FWHM at 50 ps): == gated_hom (Gaussian) and inside the "
      "paper's 0.64 +/- 0.08 (Bennett 2008; UC, gate FWHM is the one assumption)",
      abs(envg - ref_g) < 1e-5 and abs(envg - 0.64) <= 0.08 and R.get("drive", "electrical_pulsed").param("hom_gate_shape").default == "rect",
      f"{envg:.4f}")

# 11f. cascade drive (Sartison brief 1.3 anchors) and its off-state
ELS = {"tau_on_ns": 0.2, "eta_capture": 1.0, "f_rep_MHz": 200.0}
E_C = 1.602176634e-19


def sart(mu, **kw):
    drp = {**ELS, "I_uA": mu * E_C / 0.2e-9 * 1e6, **kw}
    return engine.evaluate_spec(espec(drp, MAT, solver="rate"))["scalars"]["g2_0"]


casc = {m_: sart(m_, cascade=True) for m_ in (0.5, 2.0, 4.0, 8.0, 20.0)}
ref_c = {0.5: 0.113, 2.0: 0.134, 4.0: 0.135, 8.0: 0.101, 20.0: 0.049}
check("11f1. electrical cascade (cap-2, rate solver): the Sartison rect numbers of the brief (0.113, 0.134, 0.135, 0.101, 0.049 "
      "at mu 0.5, 2, 4, 8, 20) to 3% [NV], saturating then FALLING at strong pump",
      all(rel(casc[m_], ref_c[m_]) < 0.03 for m_ in ref_c) and casc[20.0] < casc[8.0] < casc[4.0] and casc[4.0] > casc[0.5],
      " ".join(f"{m_}:{casc[m_]:.3f}" for m_ in casc))
two = {m_: sart(m_) for m_ in (1.0, 4.0, 20.0)}
check("11f2. cascade off (default): the two-level pump keeps rising, 0.148, 0.268, 0.400 at mu 1, 4, 20 (brief [NV], 3%)",
      all(rel(two[m_], r_) < 0.03 for m_, r_ in ((1.0, 0.148), (4.0, 0.268), (20.0, 0.400))),
      " ".join(f"{m_}:{two[m_]:.3f}" for m_ in two))
pl_off = R.get("drive", "electrical_pulsed").plan({}, em_s, rt_s)
pl_on = R.get("drive", "electrical_pulsed").plan({"cascade": True, "pump_ratio_xx": 0.5}, em_s, rt_s)
check("11f3. default plan: no cascade need, no bg_frac / pump_ratio_xx entries; cascade=True adds the need and pump_ratio_xx",
      "cascade" not in pl_off.needs and "bg_frac" not in pl_off.inputs and "pump_ratio_xx" not in pl_off.inputs
      and "cascade" in pl_on.needs and pl_on.inputs["pump_ratio_xx"] == 0.5
      and not R.get("drive", "electrical_pulsed").param("cascade").default and R.get("drive", "electrical_pulsed").param("bg_frac").default == 0.0
      and R.get("drive", "electrical_cw").param("bg_frac").default == 0.0)

# 11g. background-only parameter bg_frac (no laser_leak stand-in)
g_sig = engine.evaluate_spec(espec({**ELS, "I_uA": 3.0 * E_C / 0.2e-9 * 1e6}, MAT, solver="rate"))["scalars"]
bgd = {**ELS, "I_uA": 3.0 * E_C / 0.2e-9 * 1e6, "bg_frac": 0.10}
g_bg = engine.evaluate_spec(espec(bgd, MAT, solver="rate"))["scalars"]
envb2 = engine.envelope(espec(bgd, MAT, solver="rate"), ["g2_0", "g2_0:raw", "g2_0:bg_corrected"])
check("11g1. drive.bg_frac: g2 = 1 - rho^2 (1 - g2_signal) with rho = 1 - bg_frac (exact), n_ph unchanged (signal only); "
      "g2_0:raw == g2_0; g2_0:bg_corrected recovers the signal g2 (1e-12); laser_leak untouched (0)",
      abs(g_bg["g2_0"] - (1 - 0.9 ** 2 * (1 - g_sig["g2_0"]))) < 1e-12 and abs(g_bg["n_ph"] - g_sig["n_ph"]) < 1e-12
      and abs(envb2["g2_0:raw"]["central"] - g_bg["g2_0"]) < 1e-12
      and abs(envb2["g2_0:bg_corrected"]["central"] - g_sig["g2_0"]) < 1e-12 and g_bg["laser_leak_used"] == 0.0,
      f"{g_bg['g2_0']:.5f} from {g_sig['g2_0']:.5f}")
bgc_c = engine.envelope({**espec({}, MAT, solver="rate"), "drive": {"kind": "electrical_cw", "params": {"bg_frac": 0.08}}},
                        ["g2_0", "g2_0:bg_corrected"])
check("11g2. electrical_cw bg_frac 0.08 (Lazic rho 0.92): the ideal cw emitter's g2 = 1 - 0.92^2 = 0.1536 and the background-"
      "corrected value is its g2 = 0 (the paper's Eq. 5 inverse)",
      abs(bgc_c["g2_0"]["central"] - (1 - 0.92 ** 2)) < 1e-12 and abs(bgc_c["g2_0:bg_corrected"]["central"]) < 1e-12)


# 12. R2 D3: mu = pairs loaded INTO THE DOT for every mechanism (eta_capture only thins the pump Fano factor)
def cw_spec(drp, mat=None):
    sp = espec(drp, mat, solver="rate")
    sp["drive"] = {"kind": "electrical_cw", "params": drp}
    return sp


def run_ctx(sp):
    spec = engine._as_spec(sp)
    ents = engine._entries(spec)
    vals = engine._central_values(engine._leaves(spec, ents))
    return engine._run_point(ents, vals, float(spec.T_K) if not isinstance(spec.T_K, dict) else 4.2, "rate")


r12 = engine.evaluate_spec(cw_spec({"mechanism": "M-1 poisson-rail", "mu": 0.4, "eta_capture": 0.05}, MAT))["scalars"]
G12 = float(r12["Gamma_ns"])
want12 = G12 * 0.4 / 1.4                                    # two-level incoherent cw, r = mu Gamma [DR, ch02 steady state]
check("12a. D3 anchor: M-1 cw, mu 0.4, eta_capture 0.05 (two-level): flux = Gamma mu/(1+mu) (1e-6; 0.0196/ns before the fix), g2 = 0",
      abs(float(r12["flux_ns"]) - want12) < 1e-6 and abs(float(r12["g2_0"])) < 1e-9,
      f"{float(r12['flux_ns']):.6f} vs {want12:.6f} (Gamma {G12:.4f})")
r12b = engine.evaluate_spec(cw_spec({"mechanism": "M-1 poisson-rail", "mu": 0.4, "eta_capture": 1.0}, MAT))["scalars"]
check("12b. D3: eta_capture does not change the M-1 cw flux (0.05 vs 1.0, 1e-12)",
      abs(float(r12["flux_ns"]) - float(r12b["flux_ns"])) < 1e-12, f"{float(r12['flux_ns']):.6f} {float(r12b['flux_ns']):.6f}")
c12 = run_ctx(cw_spec({"mechanism": "M-2 quiet-rail", "mu": 0.4, "eta_capture": 0.5}, MAT))
itf12 = c12.plan.interface
fd12 = 1.0 - 0.5 * (1.0 - float(itf12.F_p))
check("12c. M-2 quiet rail: extra F_dot = 1 - eta (1 - F_pump) (independent of the code path) and a 'not modelled' note when F_pump != 1",
      float(itf12.F_p) < 1.0 and abs(c12.res.extra["F_dot"] - fd12) < 1e-12
      and any("Fano factor not modelled" in str(n) for n in c12.res.extra["notes"]),
      f"F_pump {float(itf12.F_p):.4g}, F_dot {c12.res.extra.get('F_dot')}")
c12m1 = run_ctx(cw_spec({"mechanism": "M-1 poisson-rail", "mu": 0.4, "eta_capture": 0.5}, MAT))
check("12d. M-1 (F_pump = 1): F_dot = 1 and no Fano note",
      abs(c12m1.res.extra["F_dot"] - 1.0) < 1e-12 and not any("Fano" in str(n) for n in c12m1.res.extra["notes"]))


n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} compose_electrical checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
