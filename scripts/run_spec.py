"""Spec-mode driver (inverse design / requirements solver, planning-doc spec
mode): given an application target g2(0) <= t at an operating temperature
T_op, print and bundle the DESIGN-TARGET SPEC SHEETS -- what the cavity/
filter/background/mesa/growth MUST deliver, per the validated F-series chain
inverted in fsim_core/spec.py.

(i)  staged device: T_op in {77, 120} K x target in {0.1, 0.5} x Delta_XX in
     {2.5, 3.5, 5.0} meV -- one spec_sheet() per combination (class-proxy
     Gamma(T)/retention [A]).
(ii) 300 K route: Gamma(300 K) in {6.0, 7.0} meV [V anchor] x Delta_XX in
     {5.0, 5.4, 8.0} meV -- rho/G/kappa/V~ requirements against a fixed,
     measured-anchor linewidth instead of the T-swept class proxy (mirrors
     run_phase3.run_map's narrow-filter-bound convention).

HONESTY: every number in both tables is a NECESSARY CONDITION under the model
-- not proof a structure meeting it exists. See fsim_core/spec.py's module
docstring. Tag chain [A] throughout; INFEASIBLE cells are marked, not hidden.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import class_proxy_params
from fsim_core.integrator import retention
from fsim_core.spec import (
    G_required,
    eps_budget,
    kappa_ceiling,
    rho_required,
    spec_sheet,
    v_tilde_required,
)
from fsim_core.spectral import epsilon_narrow_filter
from fsim_viz.figures import spec_figure

FITP = class_proxy_params()

STAGED_T_OPS = (77.0, 120.0)
STAGED_TARGETS = (0.1, 0.5)
STAGED_DELTAS = (2.5, 3.5, 5.0)

GAMMA300 = (6.0, 7.0)
DELTAS_300 = (5.0, 5.4, 8.0)
TARGET_300 = 0.5
BE_300 = (0.01, 0.1)
MU_CW = 1e-6  # mu -> 0: recovers the CW/narrow-filter-bound F1 identity g2dot = eps


# ------------------------------------------------------------ (i) staged device

def run_staged():
    rows = []
    print("(i) staged-device design-target spec sheets (class-proxy Gamma(T)/retention "
          "[A]; T_op x target x Delta_XX; THREE independent routes per point -- baseline/"
          "slit/cavity -- verdict names exactly which route(s) close, per fsim_core.spec "
          "spec_sheet()'s per-route accounting):")
    for T_op in STAGED_T_OPS:
        for target in STAGED_TARGETS:
            for delta in STAGED_DELTAS:
                s = spec_sheet(T_op, target, delta)
                rows.append(s)
                print(f"    T_op={T_op:5.0f} K  target={target:4.2f}  Delta={delta:4.1f} meV: "
                      f"eps_bud={s['eps_budget']:.4f}  w_floor={s['w_floor']:.3f} meV  "
                      f"eps_slit={s['eps_slit']:.4f}  kappa_min={s['kappa_min']:.4f} meV  "
                      f"kappa_max={s['kappa_max']:.4f} meV  "
                      f"dens<={s['density_limit_cm2']:.2e} cm^-2  "
                      f"mesa_min={s['mesa_min_um']:.3f} um  -- {s['verdict']}")
    return rows


# --------------------------------------------------------------- (ii) 300 K route

def run_300k():
    T = 300.0
    S = float(retention(T, FITP["a_esc"], FITP["E_a"], FITP["b_p"], FITP["E_b"]))
    B0 = FITP["b0"] + FITP["beta"] * (1.0 - S)
    eb = eps_budget(MU_CW, TARGET_300)  # ~= target at mu -> 0

    rows = []
    print(f"\n(ii) 300 K route (S(300)={S:.4f}, base background B0(300)={B0:.4f} [A fit], "
          f"target g2(0) <= {TARGET_300}):")
    for gam300 in GAMMA300:
        for delta in DELTAS_300:
            eps_min = epsilon_narrow_filter(delta, gam300, gam300)
            rho_req = rho_required(eps_min, MU_CW, TARGET_300)
            G_lo = G_required(rho_req, S, B0 + BE_300[0])
            G_hi = G_required(rho_req, S, B0 + BE_300[1])
            kmax = kappa_ceiling(delta, gam300, gam300, eb)
            vtil = v_tilde_required(3.0, kmax / 2.0, gam300, 1.88) if np.isfinite(kmax) else float("nan")
            row = dict(gamma300=gam300, delta=delta, eps_min=eps_min, S=S, B0=B0,
                       rho_req=rho_req, G_req_be_lo=G_lo, G_req_be_hi=G_hi,
                       kappa_max=kmax, v_tilde_req=vtil)
            rows.append(row)
            flag = "  ** INFEASIBLE (rho alone) **" if not np.isfinite(rho_req) or rho_req >= 1.0 else ""
            print(f"    Gamma(300)={gam300:.1f} meV  Delta={delta:4.1f} meV: eps_min={eps_min:.4f}  "
                  f"rho_req={rho_req:.4f}  G_req(b_e=0.01/0.1)={G_lo:.3g}/{G_hi:.3g}  "
                  f"kappa_max={kmax:.4f} meV  V~_req={vtil:.3f}{flag}")
    return rows


# ------------------------------------------------------ G_req vs b_e curves (fig)

def build_be_curves(staged_rows, k300_rows):
    """G_required(b_e) sweep for a couple of representative design points -- the
    fsim_viz figure's panel B. Physics stays here (fsim_core-derived numbers
    only); the figure module just plots the arrays."""
    b_es = np.geomspace(1e-3, 0.3, 40)
    curves = {}
    for s in staged_rows:
        if s["T_op"] == 120.0 and s["target_g2"] == 0.1:
            label = f"staged 120K t=0.1 D={s['delta_xx']:.1f}"
            G = np.array([G_required(s["rho_required"], s["S_op"], s["B0_op"] + be)
                         for be in b_es])
            curves[label] = (b_es, G)
    for k in k300_rows:
        if k["delta"] == 5.4:
            label = f"300K G(300)={k['gamma300']:.1f} D=5.4"
            G = np.array([G_required(k["rho_req"], k["S"], k["B0"] + be) for be in b_es])
            curves[label] = (b_es, G)
    return curves


# -------------------------------------------------------------------------- I/O

def _write_csv(path: Path, rows: list) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    staged_rows = run_staged()
    k300_rows = run_300k()
    be_curves = build_be_curves(staged_rows, k300_rows)

    outdir = ROOT / "out" / "spec"
    outdir.mkdir(parents=True, exist_ok=True)
    _write_csv(outdir / "spec_staged.csv", staged_rows)
    _write_csv(outdir / "spec_300K.csv", k300_rows)
    spec_figure(staged_rows, k300_rows, be_curves, outdir)

    print(f"\nspec bundle -> {outdir}  (spec_staged.csv, spec_300K.csv, spec_figure.[pdf/svg/png])")
