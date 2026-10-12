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
check("phonon_assisted / super / dichromatic are ready; arp stays a placeholder",
      all(kinds[k].status == "ready" for k in ("phonon_assisted", "super", "dichromatic"))
      and kinds["arp"].status == "needs_physics")
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

# ------------------------------------------------------------------ 10. chirp support (ARP is not a drive yet)
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

n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} polaron checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
