"""Figure 03.11: sub-turn-on Boltzmann-tail loading suppression.

Physics review (2026-09-06, item 10): plot the dot's OWN suppression
f_qfl and the background's f_qfl_bg on the SAME dimensionless axis, not
the dot's full loading mu (which folds in eta_inj, f_QD and r_dot too, on
a current axis that hides which reservoir is actually being suppressed).
Both curves are evaluated via qfl_suppression(E, V_j, kT) --
fsim_core.transport.qfl_suppression -- at the SAME (V_j, T_j) sweep point:
f_qfl at the dot's own confinement-resolved E_X (evaluate() exports only
f_qfl_bg as a scalar, so f_qfl is recomputed here from the same resolved
E_X_eV and V_j via fsim_core.device._confinement_params, the identical
confinement solve device.py itself uses -- no independent physics). The
x-axis, (qV_j - E_X)/kT, is dimensionless and shared by both curves; the
background sits at E_X + dE_WL, so its own suppression always lags the
dot's on this same axis.

Sweeps drive.I_uA on the fallback 'gainp' design card through a fresh
fsim_core.device.evaluate() call at each point (300 K heat-sink).

Run: python presentation2/figures/fig_03_11_boltzmann_tail.py
Writes presentation2/figures/out/03_11_boltzmann_tail.png at 1600x900 px, dpi=150.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fsim_core import transport  # noqa: E402
from fsim_core.device import DeviceDesign, evaluate, _confinement_params  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_11_boltzmann_tail.png"
CARD = Path(__file__).resolve().parents[2] / "cards" / "edge-inp-gainp-design.yaml"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    base_design = DeviceDesign.load(CARD)
    I0_uA = base_design.drive.I_uA

    factors = np.logspace(-3, 1, 40)
    I_uA = I0_uA * factors
    x = np.full_like(I_uA, np.nan)          # (qV_j - E_X)/kT, dimensionless
    f_qfl = np.full_like(I_uA, np.nan)       # dot's own suppression at E_X
    f_qfl_bg = np.full_like(I_uA, np.nan)    # background's own suppression at E_X + dE_WL

    for i, I in enumerate(I_uA):
        design = DeviceDesign.load(CARD)
        design.drive.I_uA = float(I)
        sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
        Tj, Vj = sc["T_j_op"], sc["V_j_op"]
        E_X_eV = _confinement_params(design.ret, Tj)["E_X_eV"]
        kT_eV = transport.KB_EV * Tj
        x[i] = (Vj - E_X_eV) / kT_eV
        f_qfl[i] = transport.qfl_suppression(E_X_eV, Vj, kT_eV)
        f_qfl_bg[i] = sc["injection.f_qfl_bg"]

    # Same quantities at the card's own default operating current (I0_uA),
    # for the vertical marker and the "f_qfl = 1" annotation.
    sc0 = evaluate(base_design, T_grid=[base_design.thermal.T_hs])["scalars"]
    Tj0, Vj0 = sc0["T_j_op"], sc0["V_j_op"]
    E_X_eV0 = _confinement_params(base_design.ret, Tj0)["E_X_eV"]
    kT_eV0 = transport.KB_EV * Tj0
    x0 = (Vj0 - E_X_eV0) / kT_eV0
    f_qfl0 = transport.qfl_suppression(E_X_eV0, Vj0, kT_eV0)
    qVj_minus_EX_meV0 = (Vj0 - E_X_eV0) * 1e3

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.semilogy(x, f_qfl, color="#1769aa", lw=2.5, marker="o", markersize=4,
                label=r"dot suppression $f_{qfl}$ at its own $E_X$")
    ax.semilogy(x, f_qfl_bg, color="#c43d3d", lw=2.5, marker="s", markersize=4,
                label=r"background suppression $f_{qfl,bg}$ at $E_X + \Delta E_{WL}$")
    ax.axvline(x0, color="gray", ls=":", lw=1.5)
    ax.annotate(f"card's default operating point:\n"
                f"$f_{{qfl}}$ = {f_qfl0:.3g} (saturated -- "
                f"$qV_j$ exceeds $E_X$ by {qVj_minus_EX_meV0:.1f} meV)\n"
                f"$f_{{qfl,bg}}$ = {sc0['injection.f_qfl_bg']:.3g} "
                f"(only the background is suppressed here)",
                xy=(x0, f_qfl0), xytext=(x0 - 11.5, 3e-2), fontsize=11, color="#333333",
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.2))

    ax.set_xlabel(r"$(qV_j - E_X)\,/\,kT$  (dimensionless)", fontsize=16)
    ax.set_ylabel("suppressed fraction (dimensionless)", fontsize=16)
    ax.set_title("Sub-turn-on Boltzmann tail: the dot saturates ($f_{qfl}\\to 1$)\n"
                 "well before the background does (gainp card, 300 K)", fontsize=15)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="lower right", fontsize=12, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}  (f_qfl0={f_qfl0:.6g}, "
          f"qVj-EX={qVj_minus_EX_meV0:.6g} meV, f_qfl_bg0={sc0['injection.f_qfl_bg']:.6g})")


if __name__ == "__main__":
    main()
