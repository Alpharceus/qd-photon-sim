"""Checks for the polaron solver and the phonon-dependent drives (P9): fsim_core/compose/{polaron,drives_phonon}.py.

Anchors are numbers the code did not produce: Lindblad pulse (the existing solver the polaron solver reduces to
at alpha -> 0), closed forms (area theorem, detuned rect pulse, detailed balance, Rabi damping of the dressed-basis
Lindblad, spherical-dot spectral density), an independent scipy ODE of the same Hamiltonian (SUPER, chirp), and the
published Quilter et al. PRL 114, 137401 (2015) population. The Quilter comparison is PRE-REGISTERED: tolerance
+-0.12 (2 sigma of the measured 0.67 +- 0.06), phonon constants are the material-registry INPUTS (ingaas_gaas),
nothing is tuned.

Run: python verify/verify_compose_polaron.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import lindblad as lb                                      # noqa: E402
from fsim_core import qd_gf                                               # noqa: E402
from fsim_core import resonant_source as rs                               # noqa: E402
from fsim_core.compose import drives, engine, polaron as PL, solvers       # noqa: E402
from fsim_core.compose import papers as PP                                # noqa: E402
from fsim_core.compose import registry as R                               # noqa: E402
from fsim_core.compose.types import DrivePlan, Pulse, Rates               # noqa: E402
from fsim_core.spectral import KB                                         # noqa: E402

RESULTS = []
HB = lb.HBAR_MEV_NS
TO_RAD_PS = 1.0 / (HB * 1e3)          # meV -> rad/ps


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


# ------------------------------------------------------------------ fixtures
EM0 = R.get("material", "ingaas_gaas").build({}, 4.2)       # registry InGaAs class constants (INPUTS)
PH = EM0.phonon
RT = Rates(Gamma_ns=1.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
SOLVER = PL.PolaronSolver()


def em_with(alpha=None):
    return replace(EM0, phonon=PH if alpha is None else qd_gf.PhononParams(alpha_ps2=alpha))


def pplan(pulses, T=4.2, f=76.2, needs=("phonons",)):
    return DrivePlan("coherent", "pulsed", f, pulses=tuple(pulses), needs=frozenset(needs), inputs={"T_K": T})


def one(area_pi, fwhm, Delta_paper_meV=0.0, T=4.2, shape="gauss", f=76.2):
    return pplan([Pulse(shape, fwhm, area_pi * math.pi, detuning_meV=-Delta_paper_meV)], T, f)


# ------------------------------------------------------------------ 1. registration, planning, selection
check("polaron solver registered with handles_phonons", R.get("solver", "polaron").handles_phonons is True)
reg = R.REGISTRY
kinds = {k: reg.get("drive", k) for k in ("phonon_assisted", "super", "dichromatic", "arp")}
check("phonon_assisted / super / dichromatic / arp are all ready (FX2c flipped arp)",
      all(kinds[k].status == "ready" for k in ("phonon_assisted", "super", "dichromatic", "arp")))
pa = kinds["phonon_assisted"].plan({"detuning_meV": 0.83, "area_pi": 8.46, "fwhm_ps": 11.88}, EM0, RT)
check("phonon_assisted plan: paper Delta +0.83 -> repo det_X -0.83, needs {phonons}, one pulse",
      len(pa.pulses) == 1 and pa.pulses[0].detuning_meV == -0.83 and pa.needs == frozenset({"phonons"})
      and abs(pa.pulses[0].area_rad - 8.46 * math.pi) < 1e-12 and pa.inputs["Delta_paper_meV"] == 0.83)
sp = kinds["super"].plan({}, EM0, RT)
check("super plan: two pulses, repo det_X = +3.5 / +7.5 (paper -3.5 / -7.5), needs {phonons, two_colour}",
      len(sp.pulses) == 2 and (sp.pulses[0].detuning_meV, sp.pulses[1].detuning_meV) == (3.5, 7.5)
      and sp.needs == frozenset({"phonons", "two_colour"}))
dc = kinds["dichromatic"].plan({"asymmetry": 3.0, "delta_meV": 1.0, "area_pi": 4.0}, EM0, RT)
a_b, a_r = dc.pulses[0].area_rad, dc.pulses[1].area_rad
check("dichromatic plan: blue det_X -delta, red +delta, area split sqrt(a/(1+a)), sqrt(1/(1+a)); intensity kept",
      dc.pulses[0].detuning_meV == -1.0 and dc.pulses[1].detuning_meV == 1.0
      and abs(a_b ** 2 + a_r ** 2 - (4.0 * math.pi) ** 2) < 1e-9 and abs(a_b / a_r - math.sqrt(3.0)) < 1e-12)
p_ok = SOLVER.supports(pa, EM0, RT)[0]
sel, iss = solvers.select(replace(pa, inputs={**pa.inputs, "T_K": 4.2}), EM0, RT, "auto")
check("select auto: a phonon plan goes to the polaron solver", p_ok is False and sel is not None
      and sel.name == "polaron", f"supports-without-T_K={p_ok}, auto -> {None if sel is None else sel.name}")
sel2, iss2 = solvers.select(replace(pa, inputs={**pa.inputs, "T_K": 4.2}), EM0, RT, "lindblad")
check("select: lindblad refuses a phonon plan (explicit request)", sel2 is None and any("phonons" in i or "P9" in i for i in iss2),
      str(iss2))
sel3, iss3 = solvers.select(pa, EM0, RT, "auto")
check("select auto without T_K: no solver, issue lists 'needs P9 polaron solver' and the reason",
      sel3 is None and any("P9 polaron" in i for i in iss3) and any("T_K" in i for i in iss3), str(iss3))

# ------------------------------------------------------------------ 2. weak-phonon limit == lindblad
em_w = em_with(0.0)
pul = Pulse("gauss", 2.0, math.pi, detuning_meV=0.0)
rt_w = Rates(Gamma_ns=1.0, gamma_star_ns=0.05, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
rp = SOLVER.run(pplan([pul]), em_w, rt_w)
rl = solvers.LindbladSolver().run(DrivePlan("coherent", "pulsed", 76.2, pulses=(pul,)), em_w, rt_w)
d = (abs(rp.g2_0 - rl.g2_0), abs(rp.inversion - rl.inversion), abs(rp.I_single - rl.I_single),
     abs(rp.n_ph - rl.n_ph))
check("alpha -> 0: polaron == lindblad on a resonant pi pulse (g2, P_X, I, n_ph) to 1e-6",
      max(d) <= 1e-6, f"|dg2|={d[0]:.2e} |dP_X|={d[1]:.2e} |dI|={d[2]:.2e} |dn|={d[3]:.2e}; g2={rl.g2_0:.5f} I={rl.I_single:.5f}")
rp2 = SOLVER.run(pplan([pul]), em_with(1e-12), rt_w)
check("alpha = 1e-12 ps^2 (renormalisation and rates ~ 0) still equals lindblad to 1e-6",
      max(abs(rp2.g2_0 - rl.g2_0), abs(rp2.inversion - rl.inversion), abs(rp2.I_single - rl.I_single)) <= 1e-6)
tl_ref = rs.TwoLevel(1.0, gamma_star_ns=0.05)
ref = rs.pulsed_photon_statistics(tl_ref, "gauss", math.pi, 2e-3, 76.2, n_steps=4096)
check("polaron (alpha 0) vs a 4096-step resonant_source reference: P_X to 1e-6 (converged-grid cross-check)",
      abs(rp.inversion - ref["P_e_post"]) <= 1e-6, f"{rp.inversion:.8f} vs {ref['P_e_post']:.8f}")

# ------------------------------------------------------------------ 3. spectral density units, detailed balance
l_nm, alpha = 4.0, 0.03
sph = qd_gf.PhononParams(l_xy_nm=l_nm, l_z_nm=l_nm, alpha_ps2=alpha)
pm_s = PL.PhononModel(sph, 4.2)
w = 1.3                                                  # rad/ps
wc = math.sqrt(2.0) * sph.c_s_m_s * 1e-3 / l_nm          # rad/ps, spherical cutoff sqrt(2) c_s / l
J_exact = alpha * w ** 3 * math.exp(-(w / wc) ** 2)      # rad/ps
J_code = float(pm_s.J_rate_ns(w * 1e3)) * 1e-3
check("J(omega) rad/ps: solver units == alpha w^3 exp(-(w/w_c)^2) (spherical dot, independent form)",
      abs(J_code - J_exact) <= 1e-9 * J_exact, f"{J_code:.9e} vs {J_exact:.9e}")
worst = 0.0
for T in (0.5, 4.2, 10.0, 77.0):
    pmT = PL.PhononModel(sph, T)
    for lam in (0.2e3, 0.9e3, 1.8e3, 3.0e3):
        gd, gu = pmT.rates(lam, 0.3)
        db = math.exp(-HB * lam / (KB * T))
        worst = max(worst, abs(gu / gd - db) / db)
check("detailed balance Gamma_up/Gamma_down = exp(-hbar Lambda/kT) to 1e-12 (T 0.5-77 K, Lambda 0.2-3/ps)",
      worst <= 1e-12, f"worst rel dev {worst:.1e}")
# (pi/2) sin^2(2 theta) J (n+1) from the instantaneous eigenvectors of H
gen_t = PL.Generator(PL.PulseTrain([Pulse("rect", 10.0, 3.0 * math.pi, detuning_meV=-0.4)]),
                     rs.TwoLevel(1.0, delta_ns=float(lb.mev_to_rate(-0.4))), PL.PhononModel(sph, 4.2))
om_r, lam_r, mix_r, _ = gen_t.dressed(np.array([5e-3]))
Delta = float(gen_t.delta)
Om = abs(om_r[0])
sin2 = Om / math.sqrt(Delta ** 2 + Om ** 2)
gd_c, _ = gen_t.pm.rates(lam_r, mix_r)
n_b = 1.0 / math.expm1(HB * lam_r[0] / (KB * 4.2))
gd_ref = 0.5 * math.pi * sin2 ** 2 * float(gen_t.pm.J_rate_ns(lam_r[0])) * (n_b + 1.0)
check("Gamma_down = (pi/2) sin^2(2 theta) J(Lambda)(n+1) from the dressed eigenvectors (sin 2 theta = Omega_r/Lambda)",
      abs(gd_c[0] - gd_ref) <= 1e-12 * gd_ref, f"{gd_c[0]:.6e} vs {gd_ref:.6e}")

# ------------------------------------------------------------------ 4. Rabi damping (closed form) and detuned rect pulse
Om_r = 1.0 * 1e3                                              # rad/ns, resonant
dur = 0.040
B_s = PL.PhononModel(sph, 4.2).B
area = Om_r / B_s * dur                                       # Omega such that B Omega = Om_r
trn = PL.PulseTrain([Pulse("rect", dur * 1e3, area, detuning_meV=0.0)])
gen_r = PL.Generator(trn, rs.TwoLevel(1e-3), PL.PhononModel(sph, 4.2))
edges, Lr, Ur = gen_r.grid(40)
rho = np.array([[1, 0], [0, 0]], dtype=complex)
v = lb.vec(rho)
Jw = alpha * (Om_r * 1e-3) ** 3 * math.exp(-((Om_r * 1e-3) / wc) ** 2) * 1e3     # 1/ns
gam = 0.5 * 0.5 * math.pi * Jw / math.tanh(HB * Om_r / (2 * KB * 4.2))            # (Gd + Gu)/2 with Gd+Gu=(pi/2) J coth
dev = 0.0
for i in range(40):
    v = Ur[i] @ v
    t = edges[i + 1]
    px = 0.5 * (1.0 - math.exp(-gam * t) * math.cos(Om_r * t))
    dev = max(dev, abs(lb.unvec(v, 2)[1, 1].real - px))
check("resonant Rabi flopping: P_X(t) = (1 - e^{-g t} cos(Omega_r t))/2, g = (pi/4) J(Omega_r) coth(hbar Omega_r/2kT) "
      "(= (Gd+Gu)/2; (pi/2) J coth is the Gd+Gu of the brief/Ramsay)", dev <= 1e-4,
      f"max dev {dev:.1e}, g={gam:.3f}/ns, Omega_r={Om_r:.0f}/ns")
Om_d, Dd, tau = 0.8 * 1e3, 0.6 * TO_RAD_PS * 1e3, 0.012              # rad/ns, rad/ns, ns
trn_d = PL.PulseTrain([Pulse("rect", tau * 1e3, Om_d * tau, detuning_meV=0.6)])
gen_d = PL.Generator(trn_d, rs.TwoLevel(1e-3, delta_ns=trn_d.det0_ns), PL.PhononModel(None, 4.2))
_, Ld, Ud = gen_d.grid(24)
v = lb.vec(rho)
for u in Ud:
    v = u @ v
Lam = math.sqrt(Om_d ** 2 + Dd ** 2)
closed = (Om_d / Lam) ** 2 * math.sin(0.5 * Lam * tau) ** 2
check("detuned rect pulse (no phonons): P_X = (Omega/Lambda)^2 sin^2(Lambda tau/2) [brief 2.1]",
      abs(lb.unvec(v, 2)[1, 1].real - closed) <= 3e-5, f"{lb.unvec(v, 2)[1, 1].real:.6f} vs {closed:.6f}")

# ------------------------------------------------------------------ 5. renormalised Rabi frequency
Zq = qd_gf.zpl_weight(PH, 4.2)
Bq = math.sqrt(Zq)
areas = np.linspace(0.9, 1.25, 15)
pxs = np.array([SOLVER.run(one(a, 4.0, 0.0), em_with(), RT, with_hom=False).inversion for a in areas])
i0 = int(np.argmax(pxs))
c = np.polyfit(areas[max(i0 - 2, 0):i0 + 3], pxs[max(i0 - 2, 0):i0 + 3], 2)
a_max = -c[1] / (2 * c[0])
check("renormalised Rabi: first Rabi maximum (resonant, 4 ps) at area pi/B within 2%",
      abs(a_max * Bq - 1.0) <= 0.02, f"max at {a_max:.4f} pi, 1/B = {1 / Bq:.4f} (B = sqrt(Z) = {Bq:.4f}, Z = {Zq:.4f})")
check("the shift is real: the maximum sits > 3% above pi (so B is doing the work)", a_max > 1.03, f"{a_max:.4f}")

# ------------------------------------------------------------------ 6. sign of phonon-assisted excitation
blue = SOLVER.run(one(8.46, 11.88, +0.83), em_with(), RT, with_hom=False)
red = SOLVER.run(one(8.46, 11.88, -0.83), em_with(), RT, with_hom=False)
check("4.2 K: BLUE-detuned (paper +0.83 meV, 16.8 ps field FWHM) inversion >> RED-detuned (trap 7)",
      blue.inversion > 0.5 and red.inversion < 0.1 and blue.inversion > 10 * red.inversion,
      f"blue {blue.inversion:.3f} (P_X max {blue.extra['P_X_max']:.3f}), red {red.inversion:.3f} (max {red.extra['P_X_max']:.3f})")
red_cold = SOLVER.run(one(8.46, 11.88, -0.83, T=0.5), em_with(), RT, with_hom=False)
check("red detuning needs phonon ABSORPTION: at 0.5 K the red inversion stays small", red_cold.inversion < 0.1,
      f"{red_cold.inversion:.3f}")
lam_pk = math.sqrt(0.83 ** 2 + (8.46 * math.pi / (11.88 / (2 * math.sqrt(math.log(2))) * math.sqrt(2 * math.pi)) / TO_RAD_PS) ** 2)
bound = 1.0 / (1.0 + math.exp(-lam_pk / (KB * 4.2)))
check("thermalised bound P_X,max <= 1/(1+exp(-hbar Lambda/kT)) at the pulse peak",
      blue.extra["P_X_max"] <= bound, f"{blue.extra['P_X_max']:.3f} <= {bound:.3f} (Lambda_peak = {lam_pk:.2f} meV)")

# ------------------------------------------------------------------ 7. Quilter 2015 anchor (pre-registered)
quilter = PP.load_paper(ROOT / "cards" / "papers" / "quilter2015.yaml")
sc = engine.evaluate_spec(quilter.device)["scalars"]
Pq, Pmax = sc["inversion"], blue.extra["P_X_max"]
print(f"      Quilter 2015: model inversion at the window end {Pq:.3f}, P_X max over the pulse {Pmax:.3f} "
      f"(measured 0.67 +- 0.06, pre-registered tolerance +-0.12); solver {sc['solver']}")
check("Quilter 2015 (UC, InGaAs registry constants, nothing tuned): P_X within 0.67 +- 0.12 "
      "(window-end value and the maximum)", abs(Pq - 0.67) <= 0.12 and abs(Pmax - 0.67) <= 0.12)
check("engine routes the card through the polaron solver and flows T_K = 4.2 K",
      sc["solver"] == "polaron" and abs(sc["T_K"] - 4.2) < 1e-12)

# ------------------------------------------------------------------ 8. dichromatic closed check
dich = kinds["dichromatic"]
rt_id = replace(RT, gamma_star_ns=0.0)
pl = dich.plan({"asymmetry": 1.0, "delta_meV": 1.0, "area_pi": 5.0, "fwhm_ps": 10.0}, EM0, rt_id)
pl = replace(pl, inputs={**pl.inputs, "T_K": 4.2})
r_d0 = SOLVER.run(pl, em_with(0.0), rt_id, with_hom=False)
check("dichromatic, a = 1, no phonons, no dephasing: net area 0 -> P_X < 1e-3",
      r_d0.inversion < 1e-3, f"{r_d0.inversion:.2e} (P_X max during the pulse {r_d0.extra['P_X_max']:.2f})")
r_dp = SOLVER.run(pl, em_with(), rt_id, with_hom=False)
check("dichromatic with phonons (P1, approximate for a fast beat) creates population: P_X > 0.1",
      r_dp.inversion > 0.1, f"{r_dp.inversion:.3f} (brief: ~0.5 class value, not an anchor)")

# ------------------------------------------------------------------ 9. SUPER swing-up (phonons off), independent ODE
def ode_two_pulse(a1, a2, fw, d1_meV, d2_meV):
    D = -d1_meV * TO_RAD_PS
    dw = (d2_meV - d1_meV) * TO_RAD_PS
    s = fw / (2 * math.sqrt(math.log(2)))
    W = 16 * s
    tc = W / 2

    def Om(a, t):
        return a * math.pi / (s * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((t - tc) / s) ** 2)

    def f(t, y):
        g, x = y
        o = Om(a1, t) + Om(a2, t) * np.exp(-1j * dw * (t - tc))
        return [-0.5j * np.conj(o) * x, -1j * (D * x + 0.5 * o * g)]
    sol = solve_ivp(f, [0, W], [1 + 0j, 0j], rtol=1e-10, atol=1e-12, method="DOP853")
    return abs(sol.y[1, -1]) ** 2


rt_sup = Rates(Gamma_ns=5.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)   # Piccinini: (200 ps)^-1


def super_run(a1, a2, fw):
    pl_ = pplan([Pulse("gauss", fw, a1 * math.pi, detuning_meV=3.5), Pulse("gauss", fw, a2 * math.pi, detuning_meV=7.5)],
                T=4.0, f=80.0, needs=("phonons", "two_colour"))
    return SOLVER.run(pl_, em_with(0.0), rt_sup, with_hom=False)


rs_ode = []
for a1, a2, fw in ((5.4, 5.4, 3.8), (18.0, 20.0, 3.8)):
    rr = super_run(a1, a2, fw)
    o = ode_two_pulse(a1, a2, fw, -3.5, -7.5)
    rs_ode.append((rr.inversion, o, rr.extra["window_ps"] * 1e-3 * 5.0))
check("two-colour propagation vs independent scipy ODE of the same H (SUPER detunings): the Lindblad value is "
      "the unitary one reduced by at most the radiative decay over the window, ODE (1 - Gamma W) <= P_X <= ODE",
      all(o * (1.0 - gw) - 1e-3 <= x <= o + 1e-3 for x, o, gw in rs_ode),
      "; ".join(f"{x:.4f} vs ODE {o:.4f} (Gamma W = {gw:.2f})" for x, o, gw in rs_ode))
best = (0.0, None)
for a1 in (38.0, 41.0, 44.0):
    for a2 in (32.0, 35.0, 38.0):
        rr = super_run(a1, a2, 3.0)
        if rr.inversion > best[0]:
            best = (rr.inversion, (a1, a2, rr.extra["P_X_max"], rr.g2_0))
print(f"      SUPER (Piccinini detunings -3.5/-7.5 meV paper, 3 ps, gamma* 0, phonons off): best inversion {best[0]:.3f} "
      f"at areas {best[1][0]:.0f} pi / {best[1][1]:.0f} pi, P_X max {best[1][2]:.3f}, g2 {best[1][3]:.3f}")
check("SUPER swing-up with phonons off: inversion > 0.9 at some area pair of the scanned grid", best[0] > 0.9,
      f"best {best[0]:.3f} at {best[1][0]:.0f}/{best[1][1]:.0f} pi")
in_range = max(super_run(a1, a2, 3.8).inversion for a1, a2 in ((18.7, 18.6), (20.0, 20.0)))
print(f"      note: inside the card area ranges (<= 20 pi, 3.8 ps) the best inversion found is {in_range:.2f}; "
      "the swing-up needs ~2x larger pulse areas than the registry range")

# ------------------------------------------------------------------ 10. chirp support (the ARP drive is block 13)
def ode_chirp(a, fw, kappa_meV_ps, det_meV=0.0):
    s = fw / (2 * math.sqrt(math.log(2)))
    W = 16 * s
    tc = W / 2
    kap = kappa_meV_ps * TO_RAD_PS                       # rad/ps^2 (laser frequency rate)
    D = -det_meV * TO_RAD_PS

    def f(t, y):
        g, x = y
        o = a * math.pi / (s * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((t - tc) / s) ** 2) \
            * np.exp(-0.5j * kap * (t - tc) ** 2)
        return [-0.5j * np.conj(o) * x, -1j * (D * x + 0.5 * o * g)]
    sol = solve_ivp(f, [0, W], [1 + 0j, 0j], rtol=1e-10, atol=1e-12, method="DOP853")
    return abs(sol.y[1, -1]) ** 2


rt_c = Rates(Gamma_ns=1.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
cp = Pulse("gauss", 8.0, 12.0 * math.pi, detuning_meV=0.0, chirp_meV_ps=0.4)
rc = SOLVER.run(pplan([cp], f=80.0, needs=("phonons", "td_detuning")), em_with(0.0), rt_c, with_hom=False)
oc = ode_chirp(12.0, 8.0, 0.4)
gw_c = rc.extra["window_ps"] * 1e-3 * rt_c.Gamma_ns
check("chirped Gaussian (ARP-type, phonons off): ODE (1 - Gamma W) <= inversion <= ODE, and > 0.9",
      oc * (1.0 - gw_c) - 1e-3 <= rc.inversion <= oc + 1e-3 and rc.inversion > 0.9,
      f"{rc.inversion:.4f} vs ODE {oc:.4f} (Gamma W = {gw_c:.3f})")

# ------------------------------------------------------------------ 11. grid convergence, HOM, engine integration
base = SOLVER.run(one(8.46, 11.88, 0.83), em_with(), RT)
old = PL.DT_FACTOR
PL.DT_FACTOR = 0.04
fine = SOLVER.run(one(8.46, 11.88, 0.83), em_with(), RT)
PL.DT_FACTOR = old
check("fine-grid convergence (dt x 0.4): P_X, g2 and I change < 1e-3 relative",
      all(abs(a - b) <= 1e-3 * max(abs(b), 1e-3) for a, b in
          ((base.inversion, fine.inversion), (base.g2_0, fine.g2_0), (base.I_single, fine.I_single))),
      f"P_X {base.inversion:.5f}/{fine.inversion:.5f}, g2 {base.g2_0:.5f}/{fine.g2_0:.5f}, I {base.I_single:.5f}/{fine.I_single:.5f}")
check("phonon HOM overlap is physical (0 < I < I without phonons) and g2 in [0, 1)",
      0.0 < base.I_single < 1.0 and 0.0 <= base.g2_0 < 1.0, f"I {base.I_single:.4f}, g2 {base.g2_0:.4f}")
import copy
spec_sign = copy.deepcopy(quilter.device)
try:
    out_b = engine.evaluate_spec(spec_sign)["scalars"]
    spec_r = copy.deepcopy(spec_sign)
    spec_r["drive"]["params"]["detuning_meV"]["value"] = -0.83
    out_r = engine.evaluate_spec(spec_r)["scalars"]
    ok_eng = out_b["solver"] == "polaron" and out_b["inversion"] > 10 * out_r["inversion"]
    det_eng = f"engine inversion blue {out_b['inversion']:.3f} red {out_r['inversion']:.3f}"
except Exception as ex:                                                  # noqa: BLE001
    ok_eng, det_eng = False, f"{type(ex).__name__}: {ex}"
check("engine end to end (drive.plan -> select -> polaron -> observables): blue >> red", ok_eng, det_eng)
try:
    est = engine.estimate(spec_sign, 3)
    ok_est = math.isfinite(est) and est > 0
except Exception as ex:                                                  # noqa: BLE001
    ok_est, est = False, repr(ex)
check("engine.estimate works for a phonon spec", ok_est, f"{est}")

# ------------------------------------------------------------------ 12. FX2b: eta_prep, Thomas phonon-assisted HOM diagnosis
def bloch_prep(shape, fwhm_ps, area_pi, G, k, gs):
    """INDEPENDENT Bloch-equation integral of (P_X + emitted photons) over the resonant pulse window:
    dP = -(G+k) P + Om y, dy = -(G+k+gs)/2 y - (Om/2)(2P - 1), dE = G P (own code, no Lindblad module)."""
    dur = fwhm_ps * 1e-3
    if shape == "rect":
        win, Om = dur, (lambda t: area_pi * math.pi / dur)
    else:
        sg = dur / (2.0 * math.sqrt(math.log(2.0)))
        win = 2.0 * rs.GAUSS_WINDOW_SIGMA * sg
        Om = lambda t: area_pi * math.pi / (sg * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((t - win / 2) / sg) ** 2)
    gc = 0.5 * (G + k + gs)
    f = lambda t, y: [-(G + k) * y[0] + Om(t) * y[1], -gc * y[1] - 0.5 * Om(t) * (2 * y[0] - 1.0), G * y[0]]
    s = solve_ivp(f, [0.0, win], [0.0, 0.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, max_step=win / 200)
    return float(s.y[0, -1]), float(s.y[2, -1])


rt_e = Rates(Gamma_ns=6.0, gamma_star_ns=0.5, k_nr_ns=0.3, L=1.0, F_P_eff=0.0, beta=0.5)
ok, det = True, ""
for shape, fw in (("rect", 60.0), ("gauss", 20.0)):
    pul_e = Pulse(shape, fw, math.pi, detuning_meV=0.0)
    r_e = SOLVER.run(pplan([pul_e], f=76.2), em_with(1e-12), rt_e, with_hom=False)
    P_end, E = bloch_prep(shape, fw, 1.0, 6.0, 0.3, 0.5)
    tol = 2e-5 if shape == "rect" else 5e-5
    ok &= abs(r_e.extra["n_exc"] - (P_end + E)) <= tol and abs(r_e.extra["P_X_end"] - P_end) <= tol
    ok &= r_e.extra["n_exc"] > r_e.extra["P_X_end"] + 1e-3 and 0.0 <= r_e.extra["eta_prep"] <= 1.0
    det += f" {shape}{fw:g}ps: n_exc {r_e.extra['n_exc']:.6f}/{P_end + E:.6f}"
check("FX2b/FX2c polaron (phonons ~ 0) n_exc = P_X(end) + emitted (the expected count) equals an independent Bloch ODE "
      "(decay, dephasing, k_nr); eta_prep = P(n >= 1) lies in [0, 1]",
      ok, det)
rl_e = solvers.LindbladSolver().run(DrivePlan("coherent", "pulsed", 76.2, pulses=(pul,)), em_w, rt_w)
rp_e = SOLVER.run(pplan([pul]), em_w, rt_w, with_hom=False)
check("FX2b/FX2c alpha -> 0: polaron eta_prep and n_exc == lindblad (existing solver) to 1e-6",
      abs(rp_e.extra["eta_prep"] - rl_e.extra["eta_prep"]) <= 1e-6 and abs(rp_e.extra["n_exc"] - rl_e.extra["n_exc"]) <= 1e-6,
      f"eta_prep {rp_e.extra['eta_prep']:.8f} vs {rl_e.extra['eta_prep']:.8f}, n_exc {rp_e.extra['n_exc']:.8f} vs "
      f"{rl_e.extra['n_exc']:.8f}")
rq = SOLVER.run(one(8.46, 11.88, 0.83), em_with(), RT, with_hom=False)
check("FX2b/FX2c phonon-assisted: n_exc >= eta_prep >= P_X_end - 1e-9 and n_ph = n_exc when k_nr = 0 (every prepared "
      "excitation is emitted); eta_prep is the probability, n_exc the count",
      rq.extra["n_exc"] >= rq.extra["eta_prep"] - 1e-12 and rq.extra["eta_prep"] >= rq.extra["P_X_end"] - 1e-9
      and abs(rq.n_ph - rq.extra["n_exc"]) <= 1e-4,
      f"eta_prep {rq.extra['eta_prep']:.5f}, n_exc {rq.extra['n_exc']:.5f}, P_X_end {rq.extra['P_X_end']:.5f}, "
      f"n_ph {rq.n_ph:.5f}")

# Thomas 2021 diagnosis: Z is applied ONCE (q^2, q = Z/Z_eff) on top of the ZPL overlap of the polaron solver
th = PP.load_paper(ROOT / "cards" / "papers" / "thomas2021.yaml")
sp_th = engine._as_spec(th.device)
ent_th = engine._entries(sp_th)
c_th = engine._run_point(ent_th, engine._central_values(engine._leaves(sp_th, ent_th)), 8.0, "auto")
Zt, Zft = c_th.d["Z_phonon"], c_th.d["Z_filt"]
I1 = float(c_th.res.I_single)
check("Thomas 2021 diagnosis: engine hom_I = I_single(polaron ZPL overlap) x (Z/Z_eff)^2 exactly (Z once, no double "
      "counting) and exceeds I_single x Z^2",
      c_th.solver.name == "polaron" and abs(c_th.d["hom_I"] - I1 * (Zt / Zft) ** 2) <= 1e-12 and c_th.d["hom_I"] > I1 * Zt ** 2,
      f"hom_I {c_th.d['hom_I']:.4f}, I_single {I1:.4f}, Z {Zt:.4f}, Z_eff {Zft:.4f}, "
      f"Gamma/(Gamma+gamma*) = {c_th.d['Gamma_ns'] / (c_th.d['Gamma_ns'] + c_th.d['gamma_star_ns']):.3f} (gamma* {c_th.d['gamma_star_ns']:.2f}/ns "
      f"from the InAs/GaAs-class T slope at 8 K), g2 {c_th.d['g2_0']:.3f}")

def _try(fn):
    try:
        fn()
    except Exception as ex:                                                 # noqa: BLE001
        return ex
    return None


def drives_phonon_numbers(fw, phi2):
    from fsim_core.compose.drives_phonon import arp_pulse_numbers
    return arp_pulse_numbers(fw, phi2)


# ------------------------------------------------------------------ 13. FX2c: HOM double counting (solver side)
# Two instantaneous pi pulses a time T apart, Gamma only (no k_nr, no dephasing): n = 0 (no emission before the second
# pulse, which then returns the dot to |G>) or 2 (emission, re-excitation, second emission) with
# q = 1 - exp(-Gamma T), g2 = 1/(2q) and photons u (before the 2nd pulse) / v (after) in orthogonal modes, each a pure
# exponential wavepacket. G1 = |u><u| + q |v><v| => Int Int |G1|^2 / N^2 = 2 q^2/(2q)^2 = 1/2 although every photon is a
# perfect wavepacket (M_s = 1). The solver's I_single must be this whole-stream overlap (1/2), NOT M_s.
G_h, T_h, fw_h = 0.5, 3000.0, 6.0
q_h = 1.0 - math.exp(-G_h * T_h * 1e-3)
rt_h = Rates(Gamma_ns=G_h, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
pul_h = [Pulse("gauss", fw_h, math.pi, detuning_meV=0.0, t0_ps=-0.5 * T_h),
         Pulse("gauss", fw_h, math.pi, detuning_meV=0.0, t0_ps=+0.5 * T_h)]
rh = SOLVER.run(pplan(pul_h, f=40.0), em_with(0.0), rt_h)
g2_h, I_h = rh.g2_0, rh.I_single
check("FX2c double pi pulse (phonons off): g2 = 1/(2q) and n_ph = 2q (n = 0 or 2) within Gamma tau_p",
      abs(g2_h - 1.0 / (2.0 * q_h)) <= 0.01 * g2_h and abs(rh.n_ph - 2.0 * q_h) <= 0.01 * 2.0 * q_h,
      f"g2 {g2_h:.4f} vs {1.0 / (2.0 * q_h):.4f}, n_ph {rh.n_ph:.4f} vs {2.0 * q_h:.4f}")
check("FX2c solver I_single of that stream = 1/2 = Int Int |G1|^2/N^2 (whole-stream overlap: the multi-photon terms "
      "are already inside, M_s = 1) -> the old ((1 - g2) I - g2) double counts g2",
      abs(I_h - 0.5) <= 0.01, f"I_single {I_h:.4f}, conv {rh.extra['hom'].get('converged')}, "
      f"V = I - g2 = {I_h - g2_h:.4f} vs old ((1-g2) I - g2) = {(1 - g2_h) * I_h - g2_h:.4f}; exact = (q-1)/(2q) = "
      f"{(q_h - 1.0) / (2.0 * q_h):.4f}")

# ------------------------------------------------------------------ 14. FX2c: ARP drive
from scipy.integrate import quad                                            # noqa: E402

arp_d = reg.get("drive", "arp")
check("arp: ready, params chirp_ps2 / fwhm_ps / area_pi / detuning_meV / f_rep_MHz, observables incl. inversion",
      arp_d.status == "ready" and {"chirp_ps2", "fwhm_ps", "area_pi", "detuning_meV", "f_rep_MHz"} <= {x.name for x in arp_d.params}
      and "inversion" in arp_d.observables and arp_d.provides == frozenset({"coherent_pulse"}))


def arp_plan(**v):
    return arp_d.plan(v, EM0, RT)


# 14a. chirp mapping vs an independent FFT: spectral phase exp(+i phi2 w^2/2) on a transform-limited Gaussian
def fft_chirp(fwhm0, phi2):
    s0 = fwhm0 / (2.0 * math.sqrt(math.log(2.0)))
    n, dtp = 2 ** 16, 0.05
    w = 2.0 * math.pi * np.fft.fftfreq(n, dtp)
    spec = np.exp(-0.5 * (s0 * w) ** 2) * np.exp(0.5j * phi2 * w ** 2)
    E = np.fft.fft(spec)                          # sum_w spec exp(-i w t): the exp(-i w t) field convention
    E = np.fft.fftshift(E)
    t = (np.arange(n) - n // 2) * dtp
    I = np.abs(E) ** 2
    I /= I.max()
    tc = float(np.sum(t * I) / np.sum(I))
    var_I = float(np.sum((t - tc) ** 2 * I) / np.sum(I))            # intensity ~ exp(-t^2/sigma_c^2): var = sigma_c^2/2
    sig_c = math.sqrt(2.0 * var_I)
    sel = I > 0.05
    ph = np.unwrap(np.angle(E[sel]))
    cq = np.polyfit(t[sel] - tc, ph, 2)               # arg E = c2 t^2 + ... ; E ~ exp(-i phi(t)), phi'' = kappa
    kappa = -2.0 * cq[0]
    return sig_c, kappa, float(np.sum(I) * dtp), s0


ok, det = True, ""
for fw0, phi2 in ((3.0, 30.0), (3.0, -30.0), (5.0, 60.0), (2.0, 5.0)):
    sc_f, k_f, _, s0f = fft_chirp(fw0, phi2)
    sc_m, fw_c, k_m, s0m = drives_phonon_numbers(fw0, phi2)
    ok &= abs(sc_f - sc_m) <= 2e-3 * sc_m and abs(k_f - k_m) <= 2e-3 * abs(k_m) and k_f * phi2 > 0
    ok &= abs(fw_c - sc_m * 2.0 * math.sqrt(math.log(2.0))) < 1e-12
    det += f" ({fw0:g} ps, {phi2:+g} ps^2): sigma_c {sc_m:.3f}/{sc_f:.3f}, kappa {k_m:.5f}/{k_f:.5f};"
check("arp chirp mapping: sigma_c = sqrt(s0^4 + phi2^2)/s0 and kappa = phi2/(s0^4 + phi2^2) (rad/ps^2) match an "
      "independent FFT of the spectrally phased Gaussian, and chirp_ps2 > 0 gives a RISING laser frequency "
      "(Wei 2014 positive chirp)", ok, det)

pl_a = arp_plan(chirp_ps2=30.0, fwhm_ps=3.0, area_pi=2.0, detuning_meV=0.4)
pa_ = pl_a.pulses[0]
tr_a = PL.PulseTrain(pl_a.pulses)
E_chirped = quad(lambda t: abs(tr_a.omega(np.array([t]))[0]) ** 2, 0.0, tr_a.W, points=[0.5 * tr_a.W], limit=400)[0]
s0a = float(rs.gaussian_sigma(3.0)) * 1e-3
A0 = 2.0 * math.pi / (s0a * math.sqrt(2.0 * math.pi))
E_tl = A0 ** 2 * s0a * math.sqrt(math.pi)
check("arp plan: det_X = -Delta_paper, needs {phonons, td_detuning}, Int |Omega|^2 dt of the chirped pulse = that of the "
      "transform-limited pulse (energy-conserving stretch, area_pi = TL area) to 1e-4, chirp rate = hbar kappa",
      pa_.detuning_meV == -0.4 and pl_a.needs == frozenset({"phonons", "td_detuning"})
      and abs(E_chirped - E_tl) <= 1e-4 * E_tl and abs(pa_.chirp_meV_ps - qd_gf.HBAR * pl_a.inputs["chirp_rad_ps2"]) < 1e-15
      and pa_.chirp_meV_ps > 0,
      f"energy {E_chirped:.6f}/{E_tl:.6f} ns^-1.., fwhm_c {pa_.fwhm_ps:.2f} ps, area_c {pa_.area_rad / math.pi:.3f} pi, "
      f"chirp {pa_.chirp_meV_ps:.4f} meV/ps")
p0 = arp_plan(chirp_ps2=0.0, fwhm_ps=4.0, area_pi=1.3)
check("arp with zero chirp is the plain Gaussian pulse (fwhm, area, no chirp)",
      abs(p0.pulses[0].fwhm_ps - 4.0) < 1e-12 and abs(p0.pulses[0].area_rad - 1.3 * math.pi) < 1e-12
      and p0.pulses[0].chirp_meV_ps == 0.0)
rej = [_try(lambda: arp_plan(chirp_ps2=float("nan"))), _try(lambda: arp_plan(fwhm_ps=0.0))]
check("arp.plan rejects a non-finite chirp and a non-positive fwhm", all(isinstance(e, ValueError) for e in rej),
      repr(rej))

# 14b. phonons off: adiabatic plateau, against the independent ODE (ode_chirp) and the plateau criterion
RT_a = Rates(Gamma_ns=0.02, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)


def arp_run(area, chirp, T=4.2, alpha0=False, rt=RT_a, **kw):
    # f_rep 1 MHz: the period (1 us) is >> 1/Gamma, so every pulse starts from |G> (the 76 MHz default would carry the
    # excitation of the previous ARP pulse over at Gamma = 0.02/ns and invert it back)
    pl = arp_d.plan({"area_pi": area, "chirp_ps2": chirp, "f_rep_MHz": 1.0, **kw}, EM0, rt)
    pl.inputs["T_K"] = T
    em = em_with(0.0) if alpha0 else EM0
    return SOLVER.run(pl, em, rt, with_hom=False), pl


ok, det, plateau = True, "", {}
for chirp in (30.0, -30.0):
    for area in (0.5, 1.5, 2.0, 3.0):
        r_a, pl_ = arp_run(area, chirp, alpha0=True)
        pc_ = pl_.pulses[0]
        ref = ode_chirp(pc_.area_rad / math.pi, pc_.fwhm_ps, pc_.chirp_meV_ps)
        gw = r_a.extra["window_ps"] * 1e-3 * RT_a.Gamma_ns
        ok &= ref * (1.0 - gw) - 2e-3 <= r_a.inversion <= ref + 2e-3
        plateau[(chirp, area)] = r_a.inversion
    det += f" chirp {chirp:+g}: " + ", ".join(f"{a}pi {plateau[(chirp, a)]:.4f}" for a in (0.5, 1.5, 2.0, 3.0)) + ";"
check("arp, phonons off, |chirp| large: inversion matches the independent ODE of the chirped pulse (areas 0.5-3 pi, both "
      "signs; radiative decay during the window allowed)", ok, det)
plat = [plateau[(c, a)] for c in (30.0, -30.0) for a in (1.5, 2.0, 3.0)]
check("arp, phonons off: robust plateau, inversion -> 1 above the adiabatic threshold (> 0.97 and varies < 2% from 1.5 pi "
      "to 3 pi, both chirp signs; 0.5 pi is below threshold)",
      min(plat) > 0.97 and max(plat) - min(plat) < 0.02 and plateau[(30.0, 0.5)] < 0.8 * min(plat),
      f"plateau {min(plat):.4f}..{max(plat):.4f}, 0.5 pi {plateau[(30.0, 0.5)]:.3f}")

# 14b2. Wei-card-like case (Gamma 2.45/ns, chirp 32.3 ps^2, 1.6 pi, 82 MHz): the end-of-window inversion is LOW only because
# the excitation made in the second half of the 288 ps window decays (Gamma W/2 = 0.35); eta_prep = P(n >= 1) is ~1
RT_w = Rates(Gamma_ns=2.45, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
r_w, pl_w = arp_run(1.6, 32.3, alpha0=True, rt=RT_w, f_rep_MHz=82.0)
pc_w = pl_w.pulses[0]
ref_w = ode_chirp(pc_w.area_rad / math.pi, pc_w.fwhm_ps, pc_w.chirp_meV_ps)
check("arp card-like case (phonons off, Gamma 2.45/ns): eta_prep > 0.99 while the end-of-window inversion is < 0.8 "
      "(radiative decay inside the window, not an incomplete transfer); the lossless ODE of the same chirped pulse gives ~1",
      r_w.extra["eta_prep"] > 0.99 and r_w.inversion < 0.8 and ref_w > 0.99,
      f"eta_prep {r_w.extra['eta_prep']:.4f}, inversion {r_w.inversion:.3f}, ODE (no decay) {ref_w:.4f}, "
      f"window {r_w.extra['window_ps']:.0f} ps")

# 14c. phonons at 4.2 K: positive chirp robust, negative chirp degraded (Wei 2014: |+> -> |-> relaxation emits phonons)
RT_b = Rates(Gamma_ns=0.02, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
eta = {(c, a): arp_run(a, c, rt=RT_b)[0].extra["eta_prep"] for c in (30.0, -30.0) for a in (1.5, 2.0, 3.0)}
drop_pos, drop_neg = eta[(30.0, 1.5)] - eta[(30.0, 3.0)], eta[(-30.0, 1.5)] - eta[(-30.0, 3.0)]
check("arp with LA phonons at 4.2 K (Wei 2014 sign): the positive-chirp branch is robust, the negative one is degraded: "
      "eta_prep(+) - eta_prep(-) > 0.15 at 3 pi, and the area drop 1.5 -> 3 pi of the positive branch is < half of the "
      "negative one",
      eta[(30.0, 3.0)] - eta[(-30.0, 3.0)] > 0.15 and drop_pos < 0.5 * drop_neg and eta[(30.0, 3.0)] > 0.8,
      f"eta_prep(+): " + ", ".join(f"{a}pi {eta[(30.0, a)]:.3f}" for a in (1.5, 2.0, 3.0))
      + "; eta_prep(-): " + ", ".join(f"{a}pi {eta[(-30.0, a)]:.3f}" for a in (1.5, 2.0, 3.0)))
print("      note: the P1 (weak-coupling) model lets the positive branch lose "
      f"{drop_pos:.3f} between 1.5 and 3 pi at 4.2 K (thermal absorption |-> -> |+> and Z renormalisation); "
      "the paper's simulation (alpha 0.022 ps^2, w_c 2 ps^-1) shows a flat plateau")

# 14d. engine end to end
try:
    arp_spec = {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "ingaas_gaas", "params": {}},
                "structure": {"kind": "bulk_planar", "params": {}},
                "drive": {"kind": "arp", "params": {}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    out_arp = engine.evaluate_spec(arp_spec)["scalars"]
    ok_arp = (out_arp["solver"] == "polaron" and 0.0 < out_arp["eta_prep"] <= 1.0 and math.isfinite(out_arp["g2_0"])
              and out_arp["inversion"] > 0.5)
    det_arp = f"solver {out_arp['solver']}, eta_prep {out_arp['eta_prep']:.3f}, inversion {out_arp['inversion']:.3f}, g2 {out_arp['g2_0']:.4f}"
except Exception as ex:                                                  # noqa: BLE001
    ok_arp, det_arp = False, f"{type(ex).__name__}: {ex}"
check("engine end to end for arp defaults (drive.plan -> select -> polaron -> observables)", ok_arp, det_arp)


# ------------------------------------------------------------------ devlib-fx3b: excitation-induced dephasing and the dark-trion branch
from fsim_core.compose import observables as OBS                           # noqa: E402,F401
from fsim_core.compose.types import EmitterModel                           # noqa: E402,F401

RT_E = Rates(Gamma_ns=1e-6, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)    # Gamma T ~ 3e-7: negligible
EM_NP = em_with(0.0)                                                        # phonons off


def eid_plan(area_rad, dur_ps, c_ps, shape="rect"):
    return DrivePlan("coherent", "pulsed", 76.2, pulses=(Pulse(shape, dur_ps, area_rad),), needs=frozenset({"phonons"}),
                     inputs={"T_K": 0.0, **({"eid_ps": c_ps} if c_ps is not None else {})})


def damped_rabi(Omega, c_ns, T_ns):
    """Closed form for a resonant rect pulse with EID: z = rho_ee - rho_gg obeys z' = -Omega y, y' = Omega z - g y,
    g = coherence decay = gamma_EID/2 = c Omega^2/2 (Gamma = gamma* = 0): z(T) = -e^{-gT/2}[cos(Om' T) + (g/2Om') sin(Om' T)]."""
    g = 0.5 * c_ns * Omega ** 2
    Op = math.sqrt(Omega ** 2 - 0.25 * g * g)
    z = -math.exp(-0.5 * g * T_ns) * (math.cos(Op * T_ns) + g / (2.0 * Op) * math.sin(Op * T_ns))
    return 0.5 * (1.0 + z)


from scipy.linalg import expm as _expm                                      # noqa: E402

ok_e, det_e = True, ""
for (Om, dur, c_ps) in ((2 * math.pi / 0.2, 300.0, 15.0), (math.pi / 0.2, 300.0, 15.0), (2 * math.pi / 0.2, 300.0, 40.0)):
    # the Generator propagates |g><g| through ONE rect pulse (no periodic steady state, so Gamma = gamma* = 0 is fine)
    tl_e = rs.TwoLevel(0.0)
    train_e = PL.PulseTrain((Pulse("rect", dur, Om * dur * 1e-3),))
    gen_e = PL.Generator(train_e, tl_e, PL.PhononModel(None, 0.0), eid_ns=c_ps * 1e-3)
    L_e = gen_e.step_L(np.array([0.5 * train_e.W]))[0]
    v_e = _expm(L_e * dur * 1e-3) @ lb.vec(rs._ground(2))
    got = float(lb.unvec(v_e, 2)[1, 1].real)
    want = damped_rabi(Om, c_ps * 1e-3, dur * 1e-3)
    ok_e &= abs(got - want) < 1e-9
    det_e += f" Om {Om:.1f}, c {c_ps}: {got:.6f} vs {want:.6f};"
check("EID-a. excitation-induced dephasing gamma_EID = c |Omega|^2 on a resonant rect pulse (Generator propagator) == the "
      "closed-form damped Rabi oscillation z(T) = -e^{-gT/2}[cos Om'T + (g/2Om') sin Om'T], g = c Om^2/2 (3 settings incl. an "
      "Om^2 rescaling, 1e-9)", ok_e, det_e)
RT_B = Rates(Gamma_ns=1.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
base_e = SOLVER.run(eid_plan(3 * math.pi, 300.0, 15.0), EM_NP, RT_B, with_hom=False)
off_e = SOLVER.run(eid_plan(3 * math.pi, 300.0, None), EM_NP, RT_B, with_hom=False)
off0 = SOLVER.run(eid_plan(3 * math.pi, 300.0, 0.0), EM_NP, RT_B, with_hom=False)
check("EID-b. default OFF is bit-identical: eid_ps absent == eid_ps = 0 (g2, n_ph, inversion exactly equal), and the damping is "
      "a visible effect when on (inversion differs by > 0.05)",
      off_e.inversion == off0.inversion and off_e.n_ph == off0.n_ph and off_e.g2_0 == off0.g2_0
      and abs(base_e.inversion - off_e.inversion) > 0.05, f"{off_e.inversion:.5f} -> {base_e.inversion:.5f}")


def ivp_eid(area, fwhm_ps, c_ps, Gam, gs):
    """Independent 2x2 Lindblad solve_ivp for a Gaussian pulse with gamma_EID(t) = c Omega(t)^2 (no repo code)."""
    dur = fwhm_ps * 1e-3
    sig = rs.gaussian_sigma(dur)
    W = rs.pulse_window("gauss", dur)
    amp = area / (sig * math.sqrt(2 * math.pi))
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sm = np.array([[0, 1], [0, 0]], dtype=complex)          # |g><e| in the basis (g, e)
    P = np.array([[0, 0], [0, 1]], dtype=complex)

    def D(c_, r):
        return c_ @ r @ c_.conj().T - 0.5 * (c_.conj().T @ c_ @ r + r @ c_.conj().T @ c_)

    def rhs(t, y):
        r = y.reshape(2, 2)
        x = t - 0.5 * W
        Om = amp * math.exp(-0.5 * (x / sig) ** 2) if abs(x) <= rs.GAUSS_WINDOW_SIGMA * sig else 0.0
        H = 0.5 * Om * sx
        return (-1j * (H @ r - r @ H) + Gam * D(sm, r) + (gs + c_ps * 1e-3 * Om ** 2) * D(P, r)).reshape(-1)
    y0 = np.array([[1, 0], [0, 0]], dtype=complex).reshape(-1)
    sol = solve_ivp(rhs, (0.0, W), y0, method="DOP853", rtol=1e-11, atol=1e-13, max_step=W / 400)
    return float(sol.y[:, -1].reshape(2, 2)[1, 1].real)


RT_G = Rates(Gamma_ns=1.0, gamma_star_ns=0.5, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
ok_g, det_g = True, ""
for (ar, fw, c_ps) in ((3 * math.pi, 20.0, 5.0), (5 * math.pi, 10.0, 20.0)):
    r_g = SOLVER.run(eid_plan(ar, fw, c_ps, "gauss"), EM_NP, RT_G, with_hom=False)
    r_i = ivp_eid(ar, fw, c_ps, 1.0, 0.5)
    ok_g &= abs(r_g.inversion - r_i) < 2e-6
    det_g += f" {ar / math.pi:.0f}pi {fw} ps c {c_ps}: {r_g.inversion:.6f} vs {r_i:.6f};"
check("EID-c. Gaussian pulses with EID (Gamma 1/ns, gamma* 0.5/ns, area 3 pi / 5 pi): inversion == an independent 2x2 Lindblad "
      "solve_ivp with the time-dependent gamma_EID(t) (2e-6)", ok_g, det_g)

# plans, routing and refusals
d_res = R.get("drive", "resonant_pulse")
pl_e = d_res.plan({"eid_ps": 12.0, "p_trion": 0.1}, EM0, RT)
pl_0 = d_res.plan({}, EM0, RT)
check("EID-d. resonant_pulse / phonon drives expose eid_ps and p_trion (default 0 = absent from plan.inputs); the planner stores them",
      "eid_ps" not in pl_0.inputs and "p_trion" not in pl_0.inputs and pl_e.inputs["eid_ps"] == 12.0 and pl_e.inputs["p_trion"] == 0.1
      and all(R.get("drive", k).param("eid_ps").default == 0.0 and R.get("drive", k).param("p_trion").default == 0.0
              for k in ("resonant_pulse", "phonon_assisted", "super", "dichromatic", "arp"))
      and R.get("drive", "tpe").param("p_trion").default == 0.0 and not any(p.name == "eid_ps" for p in R.get("drive", "tpe").params))
LS_ = solvers.LindbladSolver()
em_l = replace(EM0, phonon=None)
RT_L = Rates(Gamma_ns=1.0, gamma_star_ns=0.5, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
pl_cw = DrivePlan("coherent", "cw", None, omega_cw_ns=1.0, inputs={"eid_ps": 10.0, "T_K": 4.2})
pl_cas = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 10.0, math.pi),), needs=frozenset({"cascade"}),
                   inputs={"eid_ps": 10.0, "T_K": 4.2})
pl_ok = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 10.0, math.pi),), inputs={"eid_ps": 10.0, "T_K": 4.2})
check("EID-e. lindblad accepts EID on pulsed two-level plans (runs on the polaron step propagators, phonons off) and refuses it for cw / TPE; "
      "jc refuses EID (p_trion is post-processed on the JC result since FX4; plain Rates without a cavity mode is refused)",
      LS_.supports(pl_ok, em_l, RT_L)[0] and not LS_.supports(pl_cw, em_l, RT_L)[0] and not LS_.supports(pl_cas, em_l, RT_L)[0]
      and not R.get("solver", "jc").supports(DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 10.0, math.pi),),
                                                       inputs={"eid_ps": 10.0}), em_l, RT_L)[0])
r_lin = LS_.run(pl_ok, em_l, RT_L)
r_pol = SOLVER.run(replace(eid_plan(math.pi, 10.0, 10.0, "gauss"), f_rep_MHz=80.0), EM_NP, RT_L)
check("EID-f. LindbladSolver with eid_ps == the polaron solver without phonons (inversion, g2, n_ph, HOM identical)",
      all(abs(a - b) < 1e-12 for a, b in ((r_lin.inversion, r_pol.inversion), (r_lin.g2_0, r_pol.g2_0), (r_lin.n_ph, r_pol.n_ph),
                                          (r_lin.I_single, r_pol.I_single))))

# dark-trion branch
pl_t0 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, math.pi),), inputs={"T_K": 4.2})
pl_t1 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, math.pi),), inputs={"T_K": 4.2, "p_trion": 0.2})
t0, t1 = LS_.run(pl_t0, em_l, RT_L), LS_.run(pl_t1, em_l, RT_L)
q_ = 0.8
check("TR-a. dark trion branch at a pi pulse (p_trion 0.2): n_ph, inversion, eta_prep x (1 - p), g2 / (1 - p), HOM unchanged "
      "(independent pulses: m1 -> (1-p) m1, m2 -> (1-p) m2, |G1|^2 and N^2 both ~ (1-p)^2)",
      abs(t1.n_ph - q_ * t0.n_ph) < 1e-12 and abs(t1.inversion - q_ * t0.inversion) < 1e-12
      and abs(t1.g2_0 - t0.g2_0 / q_) < 1e-12 and abs(t1.I_single - t0.I_single) < 1e-12
      and abs(t1.extra["eta_prep"] - q_ * t0.extra["eta_prep"]) < 1e-12,
      f"n_ph {t0.n_ph:.5f} -> {t1.n_ph:.5f}")
pl_t2 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, 2 * math.pi),), inputs={"T_K": 4.2, "p_trion": 0.2})
pl_t2b = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, 2 * math.pi),), inputs={"T_K": 4.2})
pl_t3 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, 1.5 * math.pi),), inputs={"T_K": 4.2, "p_trion": 0.2})
pl_t3b = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("rect", 10.0, 1.5 * math.pi),), inputs={"T_K": 4.2})
a2, a2b = LS_.run(pl_t2, em_l, RT_L), LS_.run(pl_t2b, em_l, RT_L)
a3, a3b = LS_.run(pl_t3, em_l, RT_L), LS_.run(pl_t3b, em_l, RT_L)
p3 = 0.2 * math.sin(0.75 * math.pi) ** 2
check("TR-b. p_trion(area) = p sin^2(area/2): a 2 pi pulse leaves the trion branch empty (unchanged to 1e-12), a 1.5 pi pulse has "
      "p = 0.2 sin^2(3 pi/4) = 0.1",
      abs(a2.n_ph - a2b.n_ph) < 1e-12 and abs(a3.n_ph - (1 - p3) * a3b.n_ph) < 1e-12 and abs(p3 - 0.1) < 1e-12)
# TPE (Vajner-like): the branch uses the nominal TPE area in units of the TPE pi pulse
em_t = replace(EM0, phonon=None, delta_xx_meV=-2.0, transitions=frozenset({"X", "XX"}))
pl_x0 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 10.0, math.pi),), needs=frozenset({"cascade"}), inputs={"T_K": 4.2})
pl_x1 = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 10.0, math.pi),), needs=frozenset({"cascade"}),
                  inputs={"T_K": 4.2, "p_trion": 0.16})
from types import SimpleNamespace as _NS                                    # noqa: E402
RT_X = _NS(Gamma_ns=1.0 / 1.21, Gamma0_ns=1.0 / 1.21, k_nr_ns=0.0, gamma_star_ns=0.0)
x0, x1 = LS_.run(pl_x0, em_t, RT_X), LS_.run(pl_x1, em_t, RT_X)
check("TR-c. TPE pi pulse with p_trion 0.16: XX inversion x 0.84 (the Vajner fidelity 0.82 vs 0.98 is a 16-17% loss), "
      "n_ph, n_XX x 0.84, g2 / 0.84, hom_XX unchanged",
      abs(x1.inversion - 0.84 * x0.inversion) < 1e-12 and abs(x1.n_ph - 0.84 * x0.n_ph) < 1e-12
      and abs(x1.extra["n_XX"] - 0.84 * x0.extra["n_XX"]) < 1e-12 and abs(x1.g2_0 - x0.g2_0 / 0.84) < 1e-12
      and abs(x1.extra["hom_XX"] - x0.extra["hom_XX"]) < 1e-12, f"inv {x0.inversion:.4f} -> {x1.inversion:.4f}")
# the trion branch through the polaron machinery (phonon drive) scales the same way
r_p0 = SOLVER.run(eid_plan(math.pi, 10.0, None, "gauss"), EM_NP, RT_L, with_hom=False)
pp1 = DrivePlan("coherent", "pulsed", 76.2, pulses=(Pulse("gauss", 10.0, math.pi),), needs=frozenset({"phonons"}),
                inputs={"T_K": 0.0, "p_trion": 0.3})
r_p1 = SOLVER.run(pp1, EM_NP, RT_L, with_hom=False)
check("TR-d. the polaron solver applies the same branch (p_trion 0.3 at a pi pulse): n_ph, inversion x 0.7, g2 / 0.7",
      abs(r_p1.n_ph - 0.7 * r_p0.n_ph) < 1e-12 and abs(r_p1.inversion - 0.7 * r_p0.inversion) < 1e-12
      and abs(r_p1.g2_0 - r_p0.g2_0 / 0.7) < 1e-12)
# engine end to end (real registry): EID through resonant_pulse
eng_sp = {"schema": "fsim-composed/1", "T_K": 4.2,
          "material": {"kind": "ingaas_gaas", "params": {"tau_bare_ns": 1.0, "gamma_star_4K_ns": 0.5, "a_ac_ns_per_K": 0.0}},
          "structure": {"kind": "bulk_planar", "params": {}},
          "drive": {"kind": "resonant_pulse", "params": {"shape": "gauss", "fwhm_ps": 20.0, "area_pi": 3.0, "eid_ps": 5.0}},
          "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
o_eid = engine.evaluate_spec(eng_sp)["scalars"]
want_eid = ivp_eid(3 * math.pi, 20.0, 5.0, float(o_eid["Gamma_ns"]), float(o_eid["gamma_star_ns"]))
check("EID-g. engine end to end (resonant_pulse, 3 pi Gaussian 20 ps, eid_ps 5, material gamma* 0.5/ns): the inversion equals the independent "
      "2x2 Lindblad solve_ivp with the engine's Gamma and gamma* (1e-5), and differs from the EID-off run by > 0.01",
      abs(o_eid["inversion"] - want_eid) < 1e-5
      and abs(engine.evaluate_spec({**eng_sp, "drive": {"kind": "resonant_pulse", "params": {**eng_sp["drive"]["params"], "eid_ps": 0.0}}})
              ["scalars"]["inversion"] - o_eid["inversion"]) > 0.01,
      f"{o_eid['inversion']:.6f} vs {want_eid:.6f}, solver {o_eid['solver']}")

n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} polaron checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
