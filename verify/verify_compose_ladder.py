"""Checks for the Piccinini Eq 4 pulse and the four-level ladder solver (devlib R3-LADDER):
fsim_core/compose/ladder.py, the `sinc_gauss` SUPER shape (drives_phonon.py, polaron.py) and the ladder optimum items of
cards/papers/piccinini2025.yaml (optimum.py).

Anchors are numbers the code did not produce: the spectral width the paper prints for its pulse (0.48 meV, arXiv:2510.21428v2),
Bracht et al., PRX Quantum 2, 040354 (2021) Table 1 row 4 (99.99 % without phonons, the same row verify_compose_polaron.py uses),
closed forms (erf normalisation of the sinc pulse) and the existing polaron solver (two-level limit). The scratch four-level solver
of r3-piccinini-answer.md is a cross-check by a second implementation, not a published number.

Run: OMP_NUM_THREADS=2 python verify/verify_compose_ladder.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import lindblad as lb                                      # noqa: E402
from fsim_core.compose import drives, engine as E, ladder as LAD         # noqa: E402
from fsim_core.compose import optimum as OPT, papers as PP, polaron as PL  # noqa: E402
from fsim_core.compose import registry as R                               # noqa: E402
from fsim_core.compose.types import DrivePlan, Pulse, Rates               # noqa: E402

RESULTS = []
HBAR_PS = lb.HBAR_MEV_NS * 1e3                  # meV ps
H_PLANCK = 2.0 * math.pi * HBAR_PS              # meV ps
CARD = ROOT / "cards" / "papers" / "piccinini2025.yaml"


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


# ------------------------------------------------------------------ (i) Eq 4 pulse: spectral width and normalisation
TP = 0.923 * H_PLANCK / 0.48                    # Methods: FWHM_lambda = 0.923 lambda0^2/(c t_p)  <=>  dE = 0.923 h / t_p [DR]
RR = 10.0                                       # Eq 4, r = 10 [V]


def spectral_fwhm_meV(tp, r, half=40.0, dt=0.01):
    """Independent FFT of s(t) = sinc(t/t_p) exp[-(t/(r t_p))^2] (unbounded in practice: +-40 t_p), FWHM of |FT|^2 in meV."""
    t = np.arange(-half * tp, half * tp, dt)
    s = np.sinc(t / tp) * np.exp(-(math.pi * t / (r * tp)) ** 2)
    F = np.fft.fftshift(np.fft.fft(s, n=1 << 20))
    E_ = np.fft.fftshift(np.fft.fftfreq(1 << 20, dt)) * H_PLANCK
    I = np.abs(F) ** 2
    above = E_[I >= 0.5 * I.max()]
    return float(above.max() - above.min())


fw = spectral_fwhm_meV(TP, RR)
check("(i) Eq 4 pulse at t_p = 0.923 h / 0.48 meV = 7.95 ps: numerical spectral FWHM within 10 % of the paper's 0.48 meV [V]",
      abs(TP - 7.95) < 0.01 and abs(fw - 0.48) <= 0.048, f"t_p {TP:.3f} ps, spectral FWHM {fw:.4f} meV")
tq = np.linspace(-60 * TP, 60 * TP, 600001)
num_area = float(np.trapezoid(np.sinc(tq / TP) * np.exp(-(math.pi * tq / (RR * TP)) ** 2), tq))
check("(i) Int s dt = t_p erf(r / 2) (normalised sinc x Gaussian; the paper's t_p' Erf(r/2) with t_p' = t_p/pi) against numerical quadrature to 1e-6",
      abs(num_area / PL.sinc_gauss_area(TP, RR) - 1.0) < 1e-6, f"{num_area:.6f} vs {PL.sinc_gauss_area(TP, RR):.6f} ps")

# the drive plans the same pulse (card leaves -> Pulse), records t_p and the numerical widths
spec_card = PP.load_paper(CARD)
sp = E._as_spec(spec_card.device)
entries = E._entries(sp)
vals = E._central_values(E._leaves(sp, entries))
em, setup, sr, rt, plan, _v, _r = E._prepare(entries, vals, float(vals["T_K"]))
inp = plan.inputs
check("(i) the SUPER sinc_gauss drive derives t_p = 0.923 h/dE = 7.95 ps and records the numerical spectral and main-lobe FWHM",
      inp["shape"] == "sinc_gauss" and abs(inp["t_p_ps"] - TP) < 1e-6 and abs(inp["spectral_fwhm_numeric_meV"] - fw) < 0.01
      and 6.5 < inp["main_lobe_fwhm_ps"] < 7.5 and plan.pulses[0].tp_ps == inp["t_p_ps"] and plan.pulses[0].r == RR,
      f"t_p {inp['t_p_ps']:.3f} ps, spectral {inp['spectral_fwhm_numeric_meV']:.3f} meV, main lobe {inp['main_lobe_fwhm_ps']:.2f} ps")
tr = PL.PulseTrain(plan.pulses)
_p0 = Pulse("sinc_gauss", inp["main_lobe_fwhm_ps"], math.pi, detuning_meV=0.0, tp_ps=TP, r=RR)
_tr0 = PL.PulseTrain((_p0,))
_tt = np.linspace(0.0, _tr0.W, 400001)
_area = float(np.trapezoid(_tr0.omega(_tt).real, _tt))
check("(i) the truncated sinc_gauss Omega(t) (window +-3 Gaussian widths = +-3 r t_p/pi ~ +-76 ps) integrates to the pulse area to 3e-3",
      abs(_area / math.pi - 1.0) < 3e-3 and abs(_tr0.W * 1e3 - 6 * RR * TP / math.pi) < 1e-6, f"area/pi {_area / math.pi:.5f}, window {_tr0.W * 1e3:.1f} ps")

# ------------------------------------------------------------------ (ii) XX coupling off -> Bracht Table 1 row 4
EM_B = R.get("material", "gaas_algaas_droplet").build({}, 4.0)
RT_B = Rates(Gamma_ns=1.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
d1, d2, sg1, sg2, a1_, a2_, tau_ = -8.0, -19.1630, 2.40, 3.04, 22.65, 19.29, -0.73     # Bracht 2021 Table 1 row 4 [V] (as in verify_compose_polaron.py)
tof = lambda sg: 2 * sg * math.sqrt(math.log(2))                                       # field sigma -> intensity FWHM
plan_b = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", tof(sg1), a1_ * math.pi, detuning_meV=-d1),
                                                        Pulse("gauss", tof(sg2), a2_ * math.pi, detuning_meV=-d2, t0_ps=tau_)),
                   needs=frozenset({"phonons", "two_colour"}), inputs={"T_K": 4.0})
S_b = LAD.prepare(plan_b, EM_B, RT_B, xx_coupling=False, gamma_ns=0.0, dt_ps=0.01)
r_b = LAD.scan(S_b, [np.array([a1_]), np.array([a2_])])
px_b = float(r_b["pX_end"][0])
check("(ii) ladder with the XX coupling off, no decay, Gaussian pulses of Bracht 2021 Table 1 row 4: P_X(end) reproduces the published 99.99 % to 5e-4 [V]",
      abs(px_b - 0.9999) < 5e-4 and float(r_b["pXX_end"][0]) < 1e-12, f"P_X {px_b:.5f}")
r_b2 = LAD.scan(LAD.prepare(plan_b, EM_B, RT_B, xx_coupling=True, gamma_ns=0.0, dt_ps=0.01), [np.array([a1_]), np.array([a2_])])
print(f"      info: the same Bracht row 4 pulses WITH the 4 meV biexciton (E_b of the GaAs droplet class) give P_X {abs(float(r_b2['pX_end'][0])):.4f}, "
      f"P_XX {abs(float(r_b2['pXX_end'][0])):.4f}: with the biexciton the swing-up of that row fails and the dot returns to |G>)")

# ------------------------------------------------------------------ (iv) two-level limit vs the polaron solver (rect 2.53 ps)
A1, A2 = 7.89195, 13.95303                       # the R2 SUPCMP argmax of the two-level inversion (card history)
plan_r = replace(plan, pulses=(Pulse("rect", 2.53, A1 * math.pi, detuning_meV=3.5), Pulse("rect", 2.53, A2 * math.pi, detuning_meV=7.5)))
res_pol = PL.PolaronSolver().run(plan_r, em, rt, with_hom=False)
g_two = LAD.g2_at(LAD.prepare(plan_r, em, rt, xx_coupling=False, dt_ps=0.005), [A1, A2])
check("(iv) two-level rect 2.53 ps limit (XX off) reproduces the polaron solver's g2 0.00425 [consistency] within 2 %",
      abs(g_two["g2_0"] / res_pol.g2_0 - 1.0) < 0.02 and abs(g_two["g2_0"] - 0.00425) < 1e-4 and abs(g_two["n_ph"] - res_pol.n_ph) < 0.01,
      f"ladder {g_two['g2_0']:.5f} (n_ph {g_two['n_ph']:.4f}) vs polaron {res_pol.g2_0:.5f} (n_ph {res_pol.n_ph:.4f})")
g_four_rect = LAD.g2_at(LAD.prepare(plan_r, em, rt, fss_ueV=12.0, dt_ps=0.005), [A1, A2])
print(f"      info: the same rect pulses with the four-level XX coupling: N_X {g_four_rect['N_X']:.4f}, g2 {g_four_rect['g2_0']:.5f} "
      f"(r3-piccinini-answer.md: 0.013)")

# ------------------------------------------------------------------ (v) independent dense 4 x 4 Lindblad (full G, XH, XV, XX; restart trajectories)
def dense_ladder(plan_, em_, gamma_ns, fss_ueV, a_pi, dt_ps=0.01, restart_ps=None):
    """Independent check: dense 4 x 4 density matrix with the explicit |XV> state and generic jump operators C rho C^+ - {C^+C, rho}/2,
    Omega(t) taken straight from PulseTrain.omega (not from the ladder's unit series), batch RK4; g2 by brute-force restart
    trajectories (quantum regression) instead of the ladder's backward counting operator."""
    train = PL.PulseTrain(tuple(Pulse(p.shape, p.fwhm_ps, a * math.pi, detuning_meV=p.detuning_meV, t0_ps=p.t0_ps) for p, a in zip(plan_.pulses, a_pi)))
    det = float(train.det0_ns)
    fs = float(lb.mev_to_rate(0.5 * fss_ueV * 1e-3))
    n = int(math.ceil(train.W / (dt_ps * 1e-3)))
    dt = train.W / n
    om = train.omega(np.arange(2 * n + 1) * dt / 2)
    H0 = np.diag([0.0, det + fs, det - fs, 2 * det + float(lb.mev_to_rate(em_.delta_xx_meV))]).astype(complex)

    def Hm(o):
        h = H0.copy()
        h[1, 0], h[3, 1] = 0.5 * o, 0.5 * o
        h[0, 1], h[1, 3] = np.conj(h[1, 0]), np.conj(h[3, 1])
        return h

    J = []
    for a_, b_, g_ in ((0, 1, gamma_ns), (0, 2, gamma_ns), (1, 3, gamma_ns * em_.gamma_xx_ratio), (2, 3, gamma_ns * em_.gamma_xx_ratio)):
        C = np.zeros((4, 4), complex)
        C[a_, b_] = math.sqrt(g_)
        J.append((C, C.conj().T @ C))

    def rhs(r, h):
        out = -1j * (h @ r - r @ h)
        for C, CdC in J:
            out = out + C @ r @ C.conj().T - 0.5 * (CdC @ r + r @ CdC)
        return out

    G = np.zeros((1, 4, 4), complex)
    G[0, 0, 0] = 1.0
    hs = [Hm(o) for o in om]

    def run(r, start=None):
        nem = np.zeros(r.shape[0])
        px = np.zeros((n + 1, r.shape[0]))
        px[0] = (r[:, 1, 1] + r[:, 2, 2]).real
        for i in range(n):
            k1 = rhs(r, hs[2 * i]); k2 = rhs(r + 0.5 * dt * k1, hs[2 * i + 1])
            k3 = rhs(r + 0.5 * dt * k2, hs[2 * i + 1]); k4 = rhs(r + dt * k3, hs[2 * i + 2])
            rn = r + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            nem += 0.5 * dt * gamma_ns * ((r[:, 1, 1] + r[:, 2, 2]).real + (rn[:, 1, 1] + rn[:, 2, 2]).real)
            r = rn
            if start is not None:
                held = start >= i + 1
                r[held] = G[0]
                nem[held] = 0.0
            px[i + 1] = (r[:, 1, 1] + r[:, 2, 2]).real
        return r, nem, px

    r, nem, px = run(G.copy())
    m1 = float(nem[0] + (r[0, 1, 1] + r[0, 2, 2] + r[0, 3, 3]).real)
    out = {"N_X": m1, "pXX_end": float(r[0, 3, 3].real)}
    if restart_ps:
        idx = np.unique(np.concatenate([np.arange(0, n + 1, max(1, int(round(restart_ps * 1e-3 / dt)))), [n]]))
        rB, nB, _ = run(np.repeat(G, idx.size, axis=0), start=idx)
        m2 = 2.0 * float(np.trapezoid(gamma_ns * px[idx, 0] * (nB + (rB[:, 1, 1] + rB[:, 2, 2] + rB[:, 3, 3]).real), idx * dt))
        out["g2_0"] = m2 / m1 ** 2
    return out


plan_d = DrivePlan("coherent", "pulsed", 80.0, pulses=(Pulse("gauss", 3.0, 12 * math.pi, detuning_meV=3.5), Pulse("gauss", 3.0, 14 * math.pi, detuning_meV=7.5)),
                   needs=frozenset(), inputs={"T_K": 4.0})
RT_D = Rates(Gamma_ns=5.0, gamma_star_ns=0.0, k_nr_ns=0.0, L=1.0, F_P_eff=0.0, beta=0.5)
d_ref = dense_ladder(plan_d, EM_B, 5.0, 12.0, [12.0, 14.0], restart_ps=0.1)
d_lad = LAD.g2_at(LAD.prepare(plan_d, EM_B, RT_D, fss_ueV=12.0, dt_ps=0.01), [12.0, 14.0])
check("(v) independent dense 4 x 4 Lindblad with explicit X_V and restart trajectories (Gaussian pulses 12 pi / 14 pi at -3.5 / -7.5 meV, FSS 12 ueV, E_b 4 meV) "
      "reproduces the ladder: N_X to 1e-4, P_XX to 1e-4, g2 to 1 %",
      abs(d_ref["N_X"] - d_lad["N_X"]) < 1e-4 and abs(d_ref["pXX_end"] - d_lad["pXX_end"]) < 1e-4 and abs(d_ref["g2_0"] / d_lad["g2_0"] - 1.0) < 0.01,
      f"N_X {d_ref['N_X']:.5f} / {d_lad['N_X']:.5f}, P_XX {d_ref['pXX_end']:.5f} / {d_lad['pXX_end']:.5f}, g2 {d_ref['g2_0']:.5f} / {d_lad['g2_0']:.5f}")

# ------------------------------------------------------------------ numerical convergence of the ladder g2 (step)
S_c = [LAD.prepare(plan, em, rt, fss_ueV=12.0, dt_ps=d) for d in (0.02, 0.01)]
g_c = [LAD.g2_at(S, [15.5, 16.5]) for S in S_c]
check("ladder g2 and N_X converged in the RK4 step: dt 0.02 ps vs 0.01 ps agree to 1 % (g2) and 1e-3 (N_X)",
      abs(g_c[0]["g2_0"] / g_c[1]["g2_0"] - 1.0) < 0.01 and abs(g_c[0]["N_X"] - g_c[1]["N_X"]) < 1e-3,
      f"g2 {g_c[0]['g2_0']:.5f} / {g_c[1]['g2_0']:.5f}, N_X {g_c[0]['N_X']:.5f} / {g_c[1]['N_X']:.5f}")

# ------------------------------------------------------------------ (iii) the card's ladder optimum and X-line g2
ent = {e["item"]: e for e in spec_card.reproduce}
t0 = time.time()
opt = OPT.optimise_ladder(spec_card.device, ent["NX_SUPER_ladder_opt"]["optimum"])
t_scan = time.time() - t0
a1o, a2o = opt["values"]["drive.area1_pi"], opt["values"]["drive.area2_pi"]
g_opt = LAD.g2_at(opt["setup"], [a1o, a2o])
print(f"      info: 81 x 81 ladder area scan (N_X) took {t_scan:.0f} s; argmax ({a1o:.2f}, {a2o:.2f}) pi, N_X {opt['value']:.4f}, "
      f"P_XX(end) {g_opt['pXX_end']:.4f}, X-line g2 {g_opt['g2_0']:.5f}")
check("(iii) ladder X-line g2 at the card's N_X optimum lies in 0.03 +/- 0.02 (2 sigma of the paper's 0.03 +/- 0.01 [V])",
      abs(g_opt["g2_0"] - 0.03) <= 0.02, f"g2 {g_opt['g2_0']:.5f}")
check("(iii) the N_X optimum is a point of the fixed 0.25 pi grid, interior to the 0-20 pi window, and N_X < 1 (XX takes part of the population)",
      not opt["info"]["on_boundary"] and abs(a1o / 0.25 - round(a1o / 0.25)) < 1e-9 and abs(a2o / 0.25 - round(a2o / 0.25)) < 1e-9
      and 0.4 < opt["value"] < 0.95, f"argmax ({a1o}, {a2o}) pi, N_X {opt['value']:.4f}")
check("cross-check (second implementation, NOT a published number): the scratch solver of the R3 physics answer re-run by the Opus reviewer with the "
      "corrected window (r_eff = 10/pi) gave N_X 0.652 and g2 0.0193; the ladder agrees within 0.005 (N_X) and 0.001 (g2) "
      "(the scratch g2 uses a coarse restart-trajectory quadrature, hence the looser g2 tolerance; the scratch's original 0.665 / 0.019 used the pi-too-wide window)",
      abs(opt["value"] - 0.652) < 0.005 and abs(g_opt["g2_0"] - 0.0193) < 0.001,
      f"N_X {opt['value']:.4f} vs 0.652, g2 {g_opt['g2_0']:.5f} vs 0.0193")

# ------------------------------------------------------------------ the card items through the paper engine
store: dict = {}
items = {}
for name in ("NX_SUPER_ladder_opt", "g2_SUPER_at_opt"):
    e = ent[name]
    ir = PP.ItemResult(id=name, cls=e["class"], kind="scalar" if e["class"] == "M" else "predict", observable="g2_0", convention="raw",
                       expect=e["expect"], tol=dict(e.get("tol") or {}))
    if e["class"] == "M":
        rep = next(x for x in spec_card.reported if x["id"] == e["compare_to"])
        ir.paper, ir.err = float(rep["value"]), float(rep["err"])
    OPT.run_item(spec_card, e, ir, store)
    items[name] = ir
check("card items: NX_SUPER_ladder_opt (P) and g2_SUPER_at_opt (M, ladder, tol 2 sigma | abs 0.02) run through the paper engine and the M item is judged 'pass' "
      "only because the computed value lands in tolerance",
      items["NX_SUPER_ladder_opt"].outcome == "pass" and items["g2_SUPER_at_opt"].outcome == "pass"
      and abs(items["g2_SUPER_at_opt"].model - g_opt["g2_0"]) < 1e-9 and abs(items["NX_SUPER_ladder_opt"].model - opt["value"]) < 1e-9,
      f"N_X {items['NX_SUPER_ladder_opt'].model:.4f}, g2 {items['g2_SUPER_at_opt'].model:.5f}")

# ------------------------------------------------------------------ result
n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} compose_ladder checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
