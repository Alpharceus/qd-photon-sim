"""Post-tier work order v2 gates: T-3b (transparency), T-R (corner
retarget), T-1 (slit-held staging comparison). 2026-08-12.

T-3b: Adachi-model transparency scan at 668 nm [E, verify vs Gehrsitz --
work-order pull #2]. Result: the transparent pair bracket is
Al0.40Ga0.60As/AlAs (3.62/3.10, aggressive) to Al0.45/AlAs (3.54/3.10,
conservative) vs the assumed 3.5/3.0 -- contrast bracket 0.52/0.44 vs 0.50.

T-R: the T3 objective retargeted from "saturate the kappa ceiling" to the
Pareto corner: MAXIMIZE t_X subject to g2 <= 0.1 (staged target, at the
design point, per rail) and dT_J <= 3 K, reporting the full frontier and
per-rung sensitivity rather than a single inversion. The old kappa <= 1.06
ceiling is dropped as a constraint: it was a Lorentzian/Poisson-era proxy
for the g2 target, and the sweep shows the direct g2 constraint admits
kappa ~ 3 meV under the quiet rail. Regression: on baseline inputs this
must reproduce the sweep's verdict (10-12 pairs = the admissible cavity
rungs; slit-only also passes) from the solver itself.

T-1: slit-held staging comparison at the corner-region cavity (12-pair per
the work order; 10-pair rows included). filter.track='hold': slit centered
ON X at each operating T (the lab monochromator convention -- the reading
of 'slit-fixed' under which the F6 inversion can actually be killed, which
is T-1's stated purpose; a slit frozen at the absolute design energy keeps
XX closer to center at 77 K and the inversion survives). Cavity walks per
dn/dT [T3-computed -0.1311 meV/K [E]].

CSVs -> out/tier_device/mechanism_runs_slit_held_12pair.csv and
out/tier_geometry/corner_frontier.csv.
"""
import csv
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.cavity import emitter_energy, tracking_detuning
from fsim_core.dbr import cavity_mode, cavity_stack
from fsim_core.device import evaluate
from fsim_core.integrator import g2_from
from fsim_core.loading import f8_g2
from fsim_core.presets import preset_device

OUT_G = Path(__file__).resolve().parents[1] / "out" / "tier_geometry"
OUT_D = Path(__file__).resolve().parents[1] / "out" / "tier_device"
MU_Q, FP_Q = 0.7, 0.5
G2_TARGET, TX_FLOOR = 0.1, 0.3
LARGE = {"l_xy_nm": 7.0, "l_z_nm": 3.0}   # the sweep's winning geometry

CONTRASTS = {
    "assumed 3.5/3.0 (baseline)":        (3.5, 3.0),
    "T-3b aggressive x=0.40 (3.62/3.10)": (3.62, 3.10),
    "T-3b conservative x=0.45 (3.54/3.10)": (3.54, 3.10),
}


def kappa_ladder(n_h, n_l, lambda0=668.0, n_c=3.4, pairs=(10, 12, 14, 16, 18)):
    out = {}
    for n_top in pairs:
        top, sp, bot = cavity_stack(n_h, n_l, n_c, n_top, n_top + 6, lambda0)
        out[n_top] = cavity_mode(top + [sp] + bot, lambda0, n_out=n_h)["kappa_meV"]
    return out


def cell(kappa, T_hs=120.0, track="mode", geom=LARGE):
    d = preset_device("staged-inp-gaasp", "GaAs",
                      "planar-lambda" if kappa is not None else "none",
                      "cw-electrical")
    d.dot.lineshape = "ibm"
    d.dot.phonon = dict(geom)
    d.thermal.mesa_diameter_um = 1.0
    d.aperture.density_cm2 = 2.0e8
    d.thermal.T_hs = T_hs
    d.filter.track = track
    if kappa is not None:
        d.cavity.kappa = kappa
        d.cavity.T_track = 120.0
        d.cavity.dEdT_cav = -0.1311
    s = evaluate(d, T_grid=[T_hs])["scalars"]
    eps, rho = float(s["eps_op"]), float(s["rho_op"])
    g2_q = float(g2_from(float(f8_g2(MU_Q, FP_Q, eps)), rho)) if eps < 1 else np.nan
    return {"eps": eps, "rho2": rho**2, "t_x": float(s["t_x_op"]),
            "g2_p": float(s["g2_op"]), "g2_q": g2_q,
            "dT_J": float(s["dT_J"]), "Tj": float(s["T_j_op"])}


# ------------------------------------------------------------ T-R frontier
print("=" * 74)
print("T-R -- Pareto frontier: max t_X s.t. g2 <= 0.1 (per rail), dT_J <= 3 K")
print("=" * 74)
rows = []
for cname, (n_h, n_l) in CONTRASTS.items():
    lad = kappa_ladder(n_h, n_l)
    print(f"\n[{cname}]  kappa ladder: " +
          ", ".join(f"{p}p={k:.2f}" for p, k in lad.items()))
    frontier = []
    for label, kap in [("slit", None)] + [(f"{p} pairs", lad[p]) for p in lad]:
        c = cell(kap)
        ok_q = np.isfinite(c["g2_q"]) and c["g2_q"] <= G2_TARGET and c["dT_J"] <= 3.0
        ok_p = c["g2_p"] <= G2_TARGET and c["dT_J"] <= 3.0
        frontier.append((label, kap, c, ok_q, ok_p))
        rows.append({"contrast": cname, "config": label,
                     "kappa_meV": kap if kap else np.nan, **c,
                     "passes_g2_quiet": ok_q, "passes_g2_poisson": ok_p,
                     "passes_brightness": c["t_x"] >= TX_FLOOR})
        print(f"  {label:9s} kappa={kap if kap else float('nan'):5.2f}  "
              f"g2_q={c['g2_q']:.4f}{'*' if ok_q else ' '} "
              f"g2_P={c['g2_p']:.4f}{'*' if ok_p else ' '} "
              f"t_X={c['t_x']:.3f}{'*' if c['t_x'] >= TX_FLOOR else ' '}")
    # corner: max t_X among quiet-rail-passing; cavity corner: same among cavity rows
    passing = [f for f in frontier if f[3]]
    cav_pass = [f for f in passing if f[1] is not None
                and f[2]["t_x"] >= TX_FLOOR]
    if passing:
        best = max(passing, key=lambda f: f[2]["t_x"])
        print(f"  CORNER (all configs): {best[0]}  t_X={best[2]['t_x']:.3f} "
              f"g2_q={best[2]['g2_q']:.4f}")
    if cav_pass:
        bc = max(cav_pass, key=lambda f: f[2]["t_x"])
        last = cav_pass[-1]
        print(f"  CORNER (cavity, both floors): {bc[0]}  t_X={bc[2]['t_x']:.3f}; "
              f"last admissible rung: {last[0]} (t_X={last[2]['t_x']:.3f})")
    n_pois = sum(1 for f in frontier if f[4])
    print(f"  Poisson rail passes g2<=0.1 in {n_pois}/{len(frontier)} configs")

with (OUT_G / "corner_frontier.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ------------------------------------------------------------- T-1 slit-held
print()
print("=" * 74)
print("T-1 -- slit-held staging comparison (large lens, 120 K-tracked design)")
print("=" * 74)
lad0 = kappa_ladder(3.5, 3.0)
t1_rows = []
for pairs in (12, 10):
    kap = lad0[pairs]
    for T in (77.0, 120.0):
        for track in ("hold", "mode"):
            c = cell(kap, T_hs=T, track=track)
            # detunings at the operating point
            dx_c = 1e3 * float(tracking_detuning(c["Tj"], 1.88, 120.0, "GaAs",
                                                 -0.1311e-3))
            t1_rows.append({"pairs": pairs, "kappa_meV": kap, "T_hs_K": T,
                            "track": track, **c,
                            "X_off_slit_meV": 0.0 if track == "hold" else dx_c,
                            "XX_off_slit_meV": (-3.5 if track == "hold"
                                                else dx_c - 3.5),
                            "X_off_mode_meV": dx_c})
            r = t1_rows[-1]
            print(f"  {pairs}p {T:5.0f} K {track:4s}: g2_q={r['g2_q']:.4f} "
                  f"g2_P={r['g2_p']:.4f} t_X={r['t_x']:.3f} eps={r['eps']:.4f} "
                  f"| X@slit {r['X_off_slit_meV']:+.2f}, XX@slit "
                  f"{r['XX_off_slit_meV']:+.2f} meV")

with (OUT_D / "mechanism_runs_slit_held_12pair.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(t1_rows[0].keys()))
    w.writeheader()
    w.writerows(t1_rows)

# staging verdict
h77 = next(r for r in t1_rows if r["pairs"] == 12 and r["T_hs_K"] == 77.0
           and r["track"] == "hold")
h120 = next(r for r in t1_rows if r["pairs"] == 12 and r["T_hs_K"] == 120.0
            and r["track"] == "hold")
print(f"\nVERDICT (12-pair, quiet rail): 77 K held (g2={h77['g2_q']:.4f}, "
      f"t_X={h77['t_x']:.3f}) vs 120 K (g2={h120['g2_q']:.4f}, "
      f"t_X={h120['t_x']:.3f})")
w77 = h77["g2_q"] < h120["g2_q"] and h77["t_x"] > h120["t_x"]
w120 = h120["g2_q"] < h77["g2_q"] and h120["t_x"] > h77["t_x"]
print("  -> " + ("77 K WINS BOTH" if w77 else
                 "120 K WINS BOTH" if w120 else "SPLIT -- margin check applies"))
print(f"\nbundles -> {OUT_G / 'corner_frontier.csv'}\n"
      f"           {OUT_D / 'mechanism_runs_slit_held_12pair.csv'}")
