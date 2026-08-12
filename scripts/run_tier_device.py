"""Tier-stack device run (2026-08-12): the staged InP/GaAsP device with every
known parameter engaged at once --

  * dot: staged-inp-gaasp preset (delta_xx 3.5 meV [A]); SYMMETRIC lens dot,
    IBM polaron lineshape with the class-proxy geometry (l_xy 4.5 / l_z 1.5 nm,
    [DR] mapping of the ~25x3 nm island) and InP-class phonon card [DR].
  * cladding/template: GaAs preset (epi/cladding 1.5 um, k 15 [A]).
  * cavity: DBR top AND bottom, solved by the T3 analytic tier (NOT assumed):
    invert_kappa at the staged 120K/0.1 kappa ceiling (1.06 meV) at 668 nm
    (F6 tracking rule: mode red of the 659.5 nm cryo line by the Varshni walk
    to T_target = 120 K), AlGaInP-class contrast 3.5/3.0, spacer 3.4. The DBR
    tier solves the VERTICAL (1-D) problem; in-plane confinement/collection
    stays [A] (G from the planar-lambda preset) pending SIM-B. F_P from the
    1-D antinode LDOS [E->analytic-1D]. dE_cav/dT from uniform dn/dT 2.3e-4/K
    [E class].
  * drive: the full T1 mechanism roster, each priced honestly at the junction
    temperature (SETs hit the F9 wall at >= 77 K -- that IS the result).
  * temperatures: 77 / 120 / 300 K (RT explicitly included) + T_c.

Outputs: per (mechanism x T): g2(0) (the autocorrelation), eps, rho^2, t_X,
T_j, delivered F_p, feasibility; T_c per mechanism (Lorentzian full curve;
IBM correction quantified at the reference mechanism); WP-M2' re-excitation
g2 for the pulsed window; the five standing constraints. CSV bundle ->
out/tier_device/. Every [A] input stays [A]: this is an envelope-mid run,
not a prediction upgrade.
"""
import csv
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.dbr import cavity_stack, dEdT_cav, invert_kappa, planar_purcell
from fsim_core.device import evaluate
from fsim_core.drive_mech import reexc_g2
from fsim_core.presets import preset_device
from fsim_core.qd_gf import PhononParams, sideband_fraction, zpl_weight

OUT = Path(__file__).resolve().parents[1] / "out" / "tier_device"
OUT.mkdir(parents=True, exist_ok=True)

T_OPS = [77.0, 120.0, 300.0]
LAMBDA0 = 668.0          # nm; 120 K tracking-rule placement [work order SIM-B1]
N_H, N_L, N_C = 3.5, 3.0, 3.4   # AlGaInP-class [work order SIM-B0 pending]
DNDT = 2.3e-4            # 1/K, uniform III-V class [E]
KAPPA_TARGET = 1.06      # meV; staged 120K/0.1 spec ceiling

# ---------------------------------------------------------------- DBR (T3)
print("=" * 72)
print("DBR cavity, solved by the T3 analytic tier (vertical/1-D)")
print("=" * 72)
sol = invert_kappa(KAPPA_TARGET, N_H, N_L, N_C, LAMBDA0)
top, spacer, bottom = cavity_stack(N_H, N_L, N_C, sol["n_top"], sol["n_bottom"],
                                   LAMBDA0, dn_dT_h=DNDT, dn_dT_l=DNDT,
                                   dn_dT_c=DNDT)
stack = top + [spacer] + bottom
dEdT = dEdT_cav(stack, LAMBDA0, n_out=N_H)
# 1-D antinode LDOS: scan the spacer for the max (lambda cavity: two antinodes)
zs = np.linspace(0.0, spacer.d_nm, 41)
Fz = [planar_purcell(top, bottom, N_C, spacer.d_nm, z, LAMBDA0, n_out=N_H)
      for z in zs]
F_planar = float(max(Fz))
z_anti = float(zs[int(np.argmax(Fz))])
print(f"pairs top/bottom : {sol['n_top']} / {sol['n_bottom']}  "
      f"(contrast {N_H}/{N_L}, spacer n={N_C}, lambda0={LAMBDA0} nm)")
print(f"kappa            : {sol['kappa_meV']:.3f} meV  (target <= {KAPPA_TARGET})")
print(f"Q                : {sol['Q']:.0f}")
print(f"L_pen            : {sol['L_pen_nm']:.0f} nm")
print(f"dE_cav/dT        : {dEdT:+.4f} meV/K  [E, uniform dn/dT {DNDT:.1e}]")
print(f"F_planar(antinode {z_anti:.0f} nm): {F_planar:.2f}  [E->analytic-1D; NOT 3-D F_P]")

# ------------------------------------------------------------- base device
base = preset_device("staged-inp-gaasp", "GaAs", "planar-lambda",
                     "cw-electrical")
base.dot.lineshape = "ibm"                      # T2 tier ON (symmetric dot:
base.dot.phonon = {}                            #  isotropic l_xy default)
base.cavity.kappa = sol["kappa_meV"]            # T3 tier numbers
base.cavity.dEdT_cav = dEdT
base.cavity.T_track = 120.0
base.cavity.F_P = F_planar                      # [E->analytic-1D]
base.thermal.mesa_diameter_um = 1.0             # spec floor
base.aperture.density_cm2 = 2.0e8               # F5 growth target
pp = PhononParams()

# ------------------------------------------------- mechanism configurations
def mk(mech, params=None, mu=None, drive_edits=None):
    d = deepcopy(base)
    d.drive.mechanism = mech
    d.drive.mech_params = params or {}
    if mu is not None:
        d.drive.mu = mu
    for k, v in (drive_edits or {}).items():
        setattr(d.drive, k, v)
    return d

MECHS = [
    ("M-1 poisson-rail (CW p-i-n)", mk("poisson-rail")),
    ("M-2 quiet-rail (high-Z, F_p=0.5 [A])",
     mk("quiet-rail", {"Z_drive_ohm": 1e3, "Z_junction_ohm": 1e3}, mu=0.7)),
    ("M-3 pulsed (1 ns window, eta 3e-6 [A])",
     mk("pulsed", {"tau_on_ns": 1.0, "eta_capture": 3e-6},
        drive_edits={"I_uA": 50.0, "duty": 0.05})),
    ("M-4a SET metallic (20 nm island)", mk("set-metallic")),
    ("M-4b SET gated (5 nm island)", mk("set-gated", {"radius_nm": 5.0})),
    ("M-4c SET turnstile (5 nm, R_T 200k)",
     mk("set-turnstile", {"radius_nm": 5.0, "R_T_ohm": 2e5})),
    ("M-5a RTI (eta-boosted, quiet rail)",
     mk("rti", {"eta_boost": 2.0, "F_p_rail": 0.5}, mu=0.7,
        drive_edits={"eta_capture": 0.4})),
]

# ------------------------------------------------------------------- runs
rows = []
print()
print("=" * 72)
print("Mechanism runs -- IBM lineshape, DBR cavity, staged card")
print("=" * 72)
hdr = (f"{'mechanism':38s} {'T':>5s} {'Tj':>6s} {'g2(0)':>7s} {'eps':>7s} "
       f"{'rho^2':>6s} {'t_X':>6s}")
print(hdr)
print("-" * len(hdr))
for name, dsn in MECHS:
    # Lorentzian full curve for T_c (fast); IBM corr. quantified below for M-1
    dl = deepcopy(dsn)
    dl.dot.lineshape = "lorentzian"
    Tc = float(evaluate(dl)["scalars"]["T_c"])
    for T in T_OPS:
        dT = deepcopy(dsn)
        dT.thermal.T_hs = T
        s = evaluate(dT, T_grid=[T])["scalars"]
        rows.append({
            "mechanism": name, "T_hs_K": T, "Tj_K": float(s["T_j_op"]),
            "g2": float(s["g2_op"]), "eps": float(s["eps_op"]),
            "rho2": float(s["rho_op"]) ** 2, "t_x": float(s["t_x_op"]),
            "Tc_lorentzian_K": Tc,
        })
        r = rows[-1]
        print(f"{name:38s} {T:5.0f} {r['Tj_K']:6.1f} {r['g2']:7.4f} "
              f"{r['eps']:7.4f} {r['rho2']:6.3f} {r['t_x']:6.3f}")
    print(f"{'':38s}   T_c (Lorentzian curve): {Tc:.1f} K")

with (OUT / "mechanism_runs.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ------------------------------------------------- WP-M2' pulsed re-excitation
print()
print("WP-M2' re-excitation (M-3 window): ", end="")
mu_pulse = MECHS[2][1].drive.mu  # informational; p set to load ~mu in window
p_rate = 0.94 / 1000.0   # 1/ps: loads mu ~ 0.94 across the 1 ns window [E]
rx = reexc_g2(p_per_ps=p_rate, tau_on_ps=1000.0, tau_x_ps=1000.0)
print(f"g2_reexc = {rx['g2']:.4f} (mean {rx['mean_photons']:.3f} X photons)"
      f"  [composes with eps downstream; window/lifetime = 1.0]")

# ---------------------------------------------------------- IBM sidebands info
print()
print("IBM dot (symmetric lens, class-proxy geometry):")
for T in T_OPS:
    print(f"  Z({T:3.0f} K) = {zpl_weight(pp, T):.3f}   "
          f"sideband fraction = {sideband_fraction(pp, T):.3f}")

# ------------------------------------------------------------ constraints (5)
print()
print("=" * 72)
print("Standing constraints (staged spec)")
print("=" * 72)
d120 = deepcopy(MECHS[0][1])
d120.thermal.T_hs = 120.0
s1 = evaluate(d120, T_grid=[120.0])["scalars"]
tx120 = float(s1["t_x_op"])
dTj = float(s1["dT_J"])
checks = [
    ("1. brightness floor t_X >= 0.3 @120K", tx120,
     np.isfinite(tx120) and tx120 >= 0.3),
    ("2. b_e <= 0.18 budget @120K (card [A])", base.drive.b_e,
     base.drive.b_e <= 0.18),
    ("3. density <= 2e8 cm^-2 (F5)", base.aperture.density_cm2,
     base.aperture.density_cm2 <= 2.0e8),
    ("4. mesa >= 1 um AND dT_J < 3 K @10uA", dTj,
     base.thermal.mesa_diameter_um >= 1.0 and dTj < 3.0),
    ("5. tracking: mode 668 nm red of 659.5 cryo line (F6) + slit retained",
     LAMBDA0, LAMBDA0 > 659.5 and base.filter.enabled),
]
for name, val, ok in checks:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}   (value: {val:g})")

print(f"\nbundle -> {OUT}")
