"""Dot-geometry x cavity-geometry sweep (2026-08-12): Purcell, purity, noise.

Everything not swept stays at the tier-device run's configuration (staged
card, GaAs template, 1 um mesa, 2e8 density, b_e 0.02, auto-w slit, five
standing constraints). Swept:

  * DOT GEOMETRY (IBM PhononParams confinement lengths, symmetric in-plane;
    [A] class-proxies until AFM/TEM #8): small / default lens / large lens /
    flat pancake / tall near-spherical.
  * CAVITY GEOMETRY (T3 DBR tier, contrast 3.5/3.0, spacer 3.4): top-pair
    ladder 10..18 (bottom = top+6) -> kappa ladder ~3 meV .. ~0.25 meV, plus
    the no-cavity slit-only reference. 120 K rows use the 120 K-tracked
    design (668 nm); 300 K rows RE-PLACE the mode per F6 at T_target=300
    (695 nm) -- a 120 K-tracked cavity at RT is 49 meV detuned and would
    measure the detuning, not the device.

Reported per cell: eps, g2(0) under Poisson (M-1) AND quiet rail (M-2,
F_p=0.5, mu=0.7 -- the noise axis, composed analytically from the same
eps/rho via F8), t_X, rho^2, ZPL weight Z(T), F_planar (1-D antinode LDOS
[E->analytic-1D]) and F_eff = F_P kappa/(kappa+Gamma) (F6 iii). NOTE: the
Purcell->retention coupling is deliberately NOT wired (roadmap item 13,
"free upside"), so F_eff is reported, not fed back into rho/T_c.

CSV bundle -> out/tier_geometry/.
"""
import csv
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.cavity import purcell_eff
from fsim_core.dbr import cavity_stack, cavity_mode, planar_purcell
from fsim_core.device import evaluate
from fsim_core.integrator import g2_from
from fsim_core.loading import f8_g2
from fsim_core.presets import preset_device
from fsim_core.qd_gf import PhononParams, zpl_weight

OUT = Path(__file__).resolve().parents[1] / "out" / "tier_geometry"
OUT.mkdir(parents=True, exist_ok=True)

N_H, N_L, N_C = 3.5, 3.0, 3.4
MU_Q, FP_Q = 0.7, 0.5          # quiet-rail noise point [A]

DOTS = {
    "small (l_xy 3.0, l_z 1.0)":        {"l_xy_nm": 3.0, "l_z_nm": 1.0},
    "default lens (4.5, 1.5)":          {},
    "large lens (7.0, 3.0)":            {"l_xy_nm": 7.0, "l_z_nm": 3.0},
    "flat pancake (6.0, 1.0)":          {"l_xy_nm": 6.0, "l_z_nm": 1.0},
    "tall near-sph (3.0, 3.0)":         {"l_xy_nm": 3.0, "l_z_nm": 3.0},
}

PAIR_LADDER = [10, 12, 14, 16, 18]


def dbr_solve(lambda0):
    """kappa + F_planar per top-pair count at this design wavelength."""
    out = {}
    for n_top in PAIR_LADDER:
        top, spacer, bottom = cavity_stack(N_H, N_L, N_C, n_top, n_top + 6,
                                           lambda0)
        mode = cavity_mode(top + [spacer] + bottom, lambda0, n_out=N_H)
        zs = np.linspace(0.0, spacer.d_nm, 21)
        F = max(planar_purcell(top, bottom, N_C, spacer.d_nm, z, lambda0,
                               n_out=N_H) for z in zs)
        out[n_top] = {"kappa": mode["kappa_meV"], "Q": mode["Q"],
                      "F_planar": float(F)}
    return out


def run_cells(T_hs, T_track, lambda0, dbr):
    rows = []
    for dot_name, geom in DOTS.items():
        for cav in ["none"] + PAIR_LADDER:
            d = preset_device("staged-inp-gaasp", "GaAs",
                              "planar-lambda" if cav != "none" else "none",
                              "cw-electrical")
            d.dot.lineshape = "ibm"
            d.dot.phonon = dict(geom)
            d.thermal.mesa_diameter_um = 1.0
            d.aperture.density_cm2 = 2.0e8
            d.thermal.T_hs = T_hs
            if cav != "none":
                d.cavity.kappa = dbr[cav]["kappa"]
                d.cavity.T_track = T_track
                d.cavity.F_P = dbr[cav]["F_planar"]
                d.cavity.dEdT_cav = -0.1311   # T3-computed [E]
            s = evaluate(d, T_grid=[T_hs])["scalars"]
            eps, rho = float(s["eps_op"]), float(s["rho_op"])
            gam = float(s["gamma_op"])
            g2_p = float(s["g2_op"])                       # Poisson M-1
            g2_q = float(g2_from(float(f8_g2(MU_Q, FP_Q, eps)), rho)) \
                if eps < 1.0 else np.nan                   # quiet M-2 (F8)
            kap = dbr[cav]["kappa"] if cav != "none" else np.nan
            Fpl = dbr[cav]["F_planar"] if cav != "none" else 1.0
            rows.append({
                "T_hs_K": T_hs, "dot": dot_name,
                "cavity_top_pairs": cav, "kappa_meV": kap,
                "eps": eps, "g2_poisson": g2_p, "g2_quiet": g2_q,
                "t_x": float(s["t_x_op"]), "rho2": rho**2,
                "Z_zpl": float(zpl_weight(PhononParams(**geom), float(s["T_j_op"]))),
                "F_planar": Fpl,
                "F_eff": float(purcell_eff(Fpl, kap, gam))
                if cav != "none" else 1.0,
                "brightness_ok": bool(s["t_x_op"] >= 0.3),
            })
    return rows


print("solving DBR ladders (T3 tier)...")
dbr120 = dbr_solve(668.0)
dbr300 = dbr_solve(695.0)
for n, v in dbr120.items():
    print(f"  668 nm  top={n:2d}: kappa={v['kappa']:.3f} meV  Q={v['Q']:.0f}  "
          f"F_planar={v['F_planar']:.0f}")

print("\nsweeping 120 K (tracked at 120, 668 nm)...")
rows = run_cells(120.0, 120.0, 668.0, dbr120)
print("sweeping 300 K (RE-TRACKED at 300, 695 nm)...")
rows += run_cells(300.0, 300.0, 695.0, dbr300)

with (OUT / "geometry_sweep.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ------------------------------------------------------------------- summary
for T in (120.0, 300.0):
    sub = [r for r in rows if r["T_hs_K"] == T]
    print(f"\n===== T_hs = {T:.0f} K "
          f"({'120-tracked/668' if T == 120 else '300-tracked/695'} nm) =====")
    hdr = (f"{'dot':28s} {'cav':>4s} {'kappa':>6s} {'eps':>7s} {'g2_P':>7s} "
           f"{'g2_Q':>7s} {'t_X':>6s} {'Z':>5s} {'F_eff':>6s} {'bright':>6s}")
    print(hdr)
    for r in sub:
        kap = f"{r['kappa_meV']:.2f}" if np.isfinite(r["kappa_meV"]) else "slit"
        print(f"{r['dot']:28s} {str(r['cavity_top_pairs']):>4s} {kap:>6s} "
              f"{r['eps']:7.4f} {r['g2_poisson']:7.4f} "
              f"{r['g2_quiet'] if np.isfinite(r['g2_quiet']) else float('nan'):7.4f} "
              f"{r['t_x']:6.3f} {r['Z_zpl']:5.3f} {r['F_eff']:6.1f} "
              f"{'ok' if r['brightness_ok'] else 'FAIL':>6s}")

# best purity-with-brightness cell per T
for T in (120.0, 300.0):
    ok = [r for r in rows if r["T_hs_K"] == T and r["brightness_ok"]
          and np.isfinite(r["g2_quiet"])]
    if ok:
        best = min(ok, key=lambda r: r["g2_quiet"])
        print(f"\nbest (t_X>=0.3) @ {T:.0f} K: {best['dot']} | "
              f"cav {best['cavity_top_pairs']} | g2_quiet={best['g2_quiet']:.4f} "
              f"g2_P={best['g2_poisson']:.4f} t_X={best['t_x']:.3f} "
              f"F_eff={best['F_eff']:.1f}")
    else:
        print(f"\nbest @ {T:.0f} K: NO cell meets the brightness floor")
print(f"\nbundle -> {OUT}")
