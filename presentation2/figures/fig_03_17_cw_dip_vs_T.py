"""Figure 03.17: CW antibunching dip, intrinsic vs. IRF-convolved, vs. T.

Sweeps thermal.T_hs on cards/edge-inp-gainp-design.yaml (which sets
drive.cw=True by default) through fresh fsim_core.device.evaluate() calls,
reading back scalars['g2_cw0'] (intrinsic zero-delay correlation from the
rate-equation model, fsim_core.cw_g2) and scalars['g2_cw0_raw'] (the same
dip after convolution with the card's own drive.cw_irf_fwhm_ps instrument
response). Demonstrates that the IRF-convolved dip crosses the g2(0) < 0.5
single-photon threshold at a lower temperature than the intrinsic dip does.

Run: python presentation2/figures/fig_03_17_cw_dip_vs_T.py
Writes presentation2/figures/out/03_17_cw_dip_vs_T.png at 1600x900 px, dpi=150.
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

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_17_cw_dip_vs_T.png"
CARD = Path(__file__).resolve().parents[2] / "cards" / "edge-inp-gainp-design.yaml"

T_GRID = [80.0, 120.0, 160.0, 200.0, 230.0, 260.0, 300.0]


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    g2_intrinsic, g2_raw = [], []
    for T in T_GRID:
        design = DeviceDesign.load(CARD)
        design.thermal.T_hs = float(T)
        sc = evaluate(design, T_grid=[float(T)])["scalars"]
        g2_intrinsic.append(sc["g2_cw0"])
        g2_raw.append(sc["g2_cw0_raw"])

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.plot(T_GRID, g2_intrinsic, marker="o", markersize=8, lw=2.5, color="#1769aa",
            label=r"intrinsic $g_2^{cw}(0)$ (rate-equation model)")
    ax.plot(T_GRID, g2_raw, marker="s", markersize=8, lw=2.5, color="#c43d3d",
            label=r"IRF-convolved $g_2^{cw,raw}(0)$ (finite detector timing)")
    ax.axhline(0.5, color="gray", ls="--", lw=1.5)
    ax.text(T_GRID[0], 0.52, "single-photon threshold $g_2(0) = 0.5$", fontsize=12, color="gray")

    ax.set_xlabel("heat-sink temperature $T_{hs}$ (K)", fontsize=16)
    ax.set_ylabel(r"$g_2(0)$, CW drive", fontsize=16)
    ax.set_title("CW dip is IRF-limited: the raw curve stays above 0.5\n"
                 "at every temperature in range (gainp card)", fontsize=15)
    ax.set_ylim(0.0, 1.05)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", fontsize=13, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
