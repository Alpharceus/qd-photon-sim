"""Figure 03.11: sub-turn-on Boltzmann-tail loading suppression.

Sweeps drive.I_uA on the fallback 'gainp' design card through a fresh
fsim_core.device.evaluate() call at each point (300 K heat-sink), reading
back scalars['mu_resolved'] (the dot's own per-pulse loading, which carries
its own f_qfl suppression at E_X) and scalars['injection.f_qfl_bg'] (the
background channel's OWN, independent suppression at E_X + dE_WL). Both are
genuine repository outputs, not a hand-rolled transport calculation -- see
the section-3 spec's requirement to always go through evaluate().

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

from fsim_core.device import DeviceDesign, evaluate  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_11_boltzmann_tail.png"
CARD = Path(__file__).resolve().parents[2] / "cards" / "edge-inp-gainp-design.yaml"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    base_design = DeviceDesign.load(CARD)
    I0_uA = base_design.drive.I_uA

    factors = np.logspace(-3, 1, 40)
    I_uA = I0_uA * factors
    mu = np.full_like(I_uA, np.nan)
    f_qfl_bg = np.full_like(I_uA, np.nan)

    for i, I in enumerate(I_uA):
        design = DeviceDesign.load(CARD)
        design.drive.I_uA = float(I)
        sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
        mu[i] = sc["mu_resolved"]
        f_qfl_bg[i] = sc["injection.f_qfl_bg"]

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.loglog(I_uA * 1e3, mu, color="#1769aa", lw=2.5, marker="o", markersize=4,
              label=r"dot loading $\mu$ = $r_{dot}\,\tau_{pulse}$ (own $f_{qfl}$ suppression)")
    ax.loglog(I_uA * 1e3, f_qfl_bg, color="#c43d3d", lw=2.5, marker="s", markersize=4,
              label=r"background suppression $f_{qfl,bg}$ at $E_X + \Delta E_{WL}$")
    ax.axvline(I0_uA * 1e3, color="gray", ls=":", lw=1.5)
    ax.text(I0_uA * 1e3 * 1.1, 1.5e-3, "card's default\noperating current",
            fontsize=12, color="gray")

    ax.set_xlabel("drive current (nA)", fontsize=16)
    ax.set_ylabel("suppressed fraction (dimensionless)", fontsize=16)
    ax.set_title("Sub-turn-on Boltzmann tail: dot loading rises faster\n"
                 "than background as current increases (gainp card, 300 K)", fontsize=15)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="upper left", fontsize=12, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
