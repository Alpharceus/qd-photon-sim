"""Phase RT campaign: R1 requirement curve, R2 decision answers, R3 b_e
anchor, T-5 pulse-edge verdict, T-7 confirmation deltas.

All numbers inherit the V1-validated device tier (proposed x2 band, Raman to
confirm) and the published-dispersion inputs (data/dispersion_table_668nm.csv).
Quiet rail = mu 0.7, F_p 0.5 [A].

PLANAR CAVITY RATE (audit H5, fixed 2026-09-23): the planar lambda cavity
changes the dot's TOTAL emission rate by only ~1-2 % -- F_total_planar =
dbr.planar_total_rate at the centre antinode, ~1.013 at 695 nm / 10+16
pairs [E->analytic-1D-planar] -- so there is NO Purcell recovery from a
planar DBR cavity (Bjork et al., PRA 44, 669 (1991); Benisty et al., IEEE
JQE 34, 1612 (1998)). F_total_planar is a total-rate multiplier over ALL
channels (cavity cone, leaky, guided); it is NOT a cavity-mode F_P and is
therefore NOT passed into device.py as d.cavity.F_P (device.py reads F_P
as a single-Lorentzian mode-only factor). The purcell_wire rate path is
OFF for every planar cavity here; F_total_planar is carried as a reported
diagnostic only. Before H5 the 1-D on-axis LDOS planar_purcell (~1e2) was
wired in as F_P, which manufactured the earlier 'Purcell recovery'.
DEFERRED: applying the +1-2 % total-rate change inside the device model
needs a total-rate input in device.py (out of scope here); at this size
it cannot change any R1/R2 verdict.
CSVs -> out/rt_campaign/.
"""
import csv
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.cavity import emitter_energy
from fsim_core.dbr import cavity_mode, cavity_stack, planar_total_rate
from fsim_core.device import class_proxy_params, evaluate
from fsim_core.drive_mech import reexc_g2
from fsim_core.integrator import g2_from
from fsim_core.loading import f8_g2
from fsim_core.presets import preset_device
from fsim_core.spectral import gamma_of_T

OUT = Path(__file__).resolve().parents[1] / "out" / "rt_campaign"
OUT.mkdir(parents=True, exist_ok=True)
MU_Q, FP_Q = 0.7, 0.5
N_H, N_L, N_C = 3.495, 3.08, 3.4          # published x=0.45/AlAs design default
PROXY = class_proxy_params()
G300_PROXY = float(gamma_of_T(300.0, PROXY["gamma0"], PROXY["a_ac"],
                              PROXY["b_lo"], PROXY["E_lo"]))


def planar_total(top, sp, bot, lam0):
    """Diagnostic: TOTAL emission-rate multiplier (all channels) of an
    in-plane dipole at the centre antinode z = d_c/2 [A: standard QD
    placement] of a planar lambda cavity (dbr.planar_total_rate,
    spacer-propagating part incl. spacer-guided modes)
    [E->analytic-1D-planar, confirm: SIM-B]. NOT a cavity-mode F_P: never
    pass it through purcell_eff or into d.cavity.F_P."""
    return float(planar_total_rate(top, bot, N_C, sp.d_nm, 0.5 * sp.d_nm,
                                   lam0, n_out=N_H))


def rt_cavity(pairs=10):
    """300 K-tracked cavity at the published contrast (695 nm placement)."""
    lam0 = 1239.841984 / emitter_energy(300.0, 1.88)
    top, sp, bot = cavity_stack(N_H, N_L, N_C, pairs, pairs + 6, lam0)
    mode = cavity_mode(top + [sp] + bot, lam0, n_out=N_H)
    return mode["kappa_meV"], planar_total(top, sp, bot, lam0), lam0


def make(T, delta, gscale=1.0, cavity=None, track="hold"):
    d = preset_device("staged-inp-gaasp", "GaAs",
                      "planar-lambda" if cavity else "none", "cw-electrical")
    d.dot.lineshape = "ibm"
    d.dot.phonon = {"l_xy_nm": 7.0, "l_z_nm": 3.0}
    d.dot.delta_xx = delta
    d.dot.gamma_scale = gscale
    d.thermal.mesa_diameter_um = 1.0
    d.aperture.density_cm2 = 2.0e8
    d.thermal.T_hs = T
    d.filter.track = track
    if cavity:
        kap, _F_total_planar, _ = cavity   # diagnostic only, never F_P
        d.cavity.kappa = kap
        d.cavity.T_track = T
        d.cavity.dEdT_cav = -0.104        # V3.3 per-layer value [DR-based]
        # Planar DBR: rate wiring OFF (H5). d.cavity.F_P is left at the
        # preset value, which only feeds device.py's reported F_eff scalar
        # (not read here) while purcell_wire is False.
        d.cavity.purcell_wire = False
    s = evaluate(d, T_grid=[T])["scalars"]
    eps, rho = float(s["eps_op"]), float(s["rho_op"])
    g2q = float(g2_from(float(f8_g2(MU_Q, FP_Q, eps)), rho)) if eps < 1 else np.nan
    return {"g2_q": g2q, "g2_p": float(s["g2_op"]), "t_x": float(s["t_x_op"]),
            "eps": eps, "rho2": rho**2}


# ================================================================ R2 answers
print("=" * 74)
print("R2 -- planar-cavity Purcell answers [E->analytic-1D-planar, confirm: SIM-B]")
print("=" * 74)
cav_rt = rt_cavity(10)
print(f"RT cavity (10p, 695 nm): kappa={cav_rt[0]:.2f} meV, "
      f"F_total_planar={cav_rt[1]:.3f} (total-rate multiplier, all channels; "
      f"diagnostic, not wired)")
a_cav = make(300.0, 4.0, cavity=cav_rt)
print(f"300 K, Delta=4.0 (new prior central), planar cavity (spectral "
      f"filtering only): t_X {a_cav['t_x']:.3f}; g2_q {a_cav['g2_q']:.3f}")
print(f"  ANSWER (a): Purcell recovery of t_X >= 0.3 at 300 K: NO -- a planar "
      f"DBR cavity changes the total rate by {100 * (cav_rt[1] - 1):+.1f} % "
      f"(F_total_planar {cav_rt[1]:.3f}), about 1-2 %, not a Purcell "
      f"enhancement (t_X = {a_cav['t_x']:.3f} without rate wiring; applying "
      f"the ~1 % change is DEFERRED, see module docstring)")
lam120 = 668.0
top, sp, bot = cavity_stack(N_H, N_L, N_C, 10, 16, lam120)
kap120 = cavity_mode(top + [sp] + bot, lam120, n_out=N_H)["kappa_meV"]
F120 = planar_total(top, sp, bot, lam120)
b_slit = make(120.0, 4.0)
b_cav = make(120.0, 4.0, cavity=(kap120, F120, lam120))
print(f"120 K fork, Delta=4.0: slit-only g2_q={b_slit['g2_q']:.4f} "
      f"t_X={b_slit['t_x']:.3f} | planar cavity (no rate wiring, "
      f"F_total_planar {F120:.3f}) g2_q={b_cav['g2_q']:.4f} "
      f"t_X={b_cav['t_x']:.3f}")
print(f"  ANSWER (b): without Purcell (planar total-rate change ~1-2 %), the cavity route "
      f"{'BEATS' if (b_cav['g2_q'] < b_slit['g2_q'] and b_cav['t_x'] > b_slit['t_x']) else 'does NOT dominate'}"
      f" slit-only at the staged point (plus unmodeled G on top)")

# ================================================================ R1 curve
print()
print("=" * 74)
print("R1 -- Delta_XX*(Gamma_300): the RT closure threshold")
print("=" * 74)
rows = []
deltas = np.arange(1.0, 8.01, 0.5)
print(f"{'Gamma300':>8} | {'Dxx*(0.5) slit':>14} {'Dxx*(0.1) slit':>14} "
      f"{'Dxx*(0.5) cav':>15} {'Dxx*(0.1) cav':>15}   (quiet rail; planar cavity, no rate wiring)")
for g300 in (3.0, 5.0, 6.6, 8.0, 10.0):
    gs = g300 / G300_PROXY
    curves = {}
    for cfg, cav in (("slit", None), ("cav", cav_rt)):
        g2s = []
        for dd in deltas:
            c = make(300.0, float(dd), gscale=gs, cavity=cav)
            g2s.append(c["g2_q"])
            rows.append({"Gamma300_meV": g300, "delta_xx_meV": float(dd),
                         "config": cfg, **c})
        curves[cfg] = np.array(g2s)

    def cross(y, level):
        for i in range(len(y) - 1):
            if np.isfinite(y[i]) and np.isfinite(y[i + 1]) \
                    and y[i] > level >= y[i + 1]:
                return float(np.interp(level, [y[i + 1], y[i]],
                                       [deltas[i + 1], deltas[i]]))
        return np.nan
    xs = [cross(curves["slit"], 0.5), cross(curves["slit"], 0.1),
          cross(curves["cav"], 0.5), cross(curves["cav"], 0.1)]
    print(f"{g300:8.1f} | {xs[0]:14.2f} {xs[1]:14.2f} {xs[2]:15.2f} "
          f"{xs[3]:15.2f}")

with (OUT / "r1_delta_xx_requirement.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("anchors: family prior Dxx = +2..+6 meV (central ~4, Wimmer size-rising);"
      " piezo lever > 5 meV; Gamma300 UNMEASURED on platform (single-dot data"
      " stop at ~150 K)")

# ================================================================ R3 anchor
print()
print("=" * 74)
print("R3 -- b_e(300 K) provenance [DR-anchored range, replaces bare [A]]")
print("=" * 74)
KB = 0.08617333262
for E_act in (50.0, 100.0, 150.0):
    mult = np.exp(-E_act / (KB * 300.0)) / np.exp(-E_act / (KB * 80.0))
    print(f"  E_act={E_act:5.0f} meV: b_e(300)/b_e(80) = {mult:9.1f}x")
print("""  Anchor: family single-dot EL is background-limited at ~80 K
  (Reischle OE16: EL followed to 80 K, 'limit = background, not signal') ->
  b_e(80 K, plain family LED) = O(0.3-1) signal units [DR]. Extrapolation:
  b_e(300 K, unengineered) ~ 2e2..4e4 x that -- the RT b_e budget (<=0.01-0.1)
  therefore demands 4-6 ORDERS of background suppression vs the family
  baseline. b_e = 0.02 at RT is not an assumption but an ENGINEERING
  REQUIREMENT (aperture + spectral separation + injection selectivity);
  carried as such into R4.""")

# ================================================================ T-5 verdict
print("=" * 74)
print("T-5 -- pulse-edge RC and the M-3 verdict")
print("=" * 74)
eps0, epsr = 8.854e-12, 12.5
A = np.pi * (0.5e-6) ** 2
for d_i, pad_fF in ((0.2e-6, 100.0), (0.2e-6, 300.0), (0.2e-6, 1000.0)):
    Cj = eps0 * epsr * A / d_i
    Ctot = Cj + pad_fF * 1e-15
    R = 50.0 + 100.0                      # driver + series [E bracket]
    t_edge = 2.2 * R * Ctot
    print(f"  C_j={Cj*1e15:.2f} fF + pad {pad_fF:.0f} fF, R=150 Ohm: "
          f"10-90% edge = {t_edge*1e12:.0f} ps")
for w_ps in (100.0, 200.0, 500.0, 1000.0):
    r = reexc_g2(p_per_ps=0.94 / w_ps, tau_on_ps=w_ps, tau_x_ps=1000.0)
    print(f"  window {w_ps:6.0f} ps: g2_reexc = {r['g2']:.4f}"
          f"{'  <= 0.05 bar' if r['g2'] <= 0.05 else ''}")
print("  VERDICT: M-3 survives the staged 0.1 bar ONLY with ~100-200 ps"
      " windows, i.e. pad capacitance <= ~300 fF -- plausible but a real"
      " constraint; at 1 ns windows pulsed is excluded (g2_reexc ~ 0.3).")

# ================================================================ T-7 deltas
print()
print("=" * 74)
print("T-7 -- confirmation deltas vs tag tier-runs-2026-08-12 baselines")
print("=" * 74)
base_cells = {
    "120K corner cell (10p, quiet)": (0.0427, 0.407),
    "77K slit-only (Poisson)": (0.0486, 0.449),
    "77K-designed cavity 10p (Poisson)": (0.0242, 0.435),
}
new_cells = {}
lad10 = kap120
c = make(120.0, 4.0, cavity=(kap120, F120, lam120))
new_cells["120K corner cell (10p, quiet)"] = (c["g2_q"], c["t_x"])
c = make(77.0, 4.0)
new_cells["77K slit-only (Poisson)"] = (c["g2_p"], c["t_x"])
lam77 = 1239.841984 / emitter_energy(77.0, 1.88)
t7, s7, b7 = cavity_stack(N_H, N_L, N_C, 10, 16, lam77)
k77 = cavity_mode(t7 + [s7] + b7, lam77, n_out=N_H)["kappa_meV"]
F77 = planar_total(t7, s7, b7, lam77)
c = make(77.0, 4.0, cavity=(k77, F77, lam77))
new_cells["77K-designed cavity 10p (Poisson)"] = (c["g2_p"], c["t_x"])
print(f"{'cell':38s} {'g2 old':>8} {'g2 new':>8} {'t_X old':>8} {'t_X new':>8}")
for k in base_cells:
    print(f"{k:38s} {base_cells[k][0]:8.4f} {new_cells[k][0]:8.4f} "
          f"{base_cells[k][1]:8.3f} {new_cells[k][1]:8.3f}")
print("(new = published contrast + per-layer dEdT -0.104 + Delta_XX 4.0 "
      "central prior; old = tag baseline at assumed contrast, -0.1311, "
      "Delta 3.5)")
print(f"\nbundle -> {OUT}")
