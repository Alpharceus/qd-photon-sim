"""Figure 03.05: homojunction built-in voltage V_bi vs. doping level.

Uses fsim_core.transport.homojunction_vbi (the textbook closed form
V_bi = (kT/q) ln(N_A N_D / n_i^2)) on a GaAs homojunction, sweeping
N_A = N_D from 1e15 to 1e19 cm^-3 at 300 K, and marks the 1e17/1e17 cm^-3
point cited on slide 03-05.

Run: python presentation2/figures/fig_03_05_vbi_doping.py
Writes presentation2/figures/out/03_05_vbi_doping.png at 1600x900 px, dpi=150.
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
from fsim_core.materials import binary  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_05_vbi_doping.png"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    gaas = binary("GaAs")
    N = np.logspace(15, 19, 200)
    vbi = np.array([transport.homojunction_vbi(gaas, n, n, 300.0) for n in N])

    N_mark = 1e17
    vbi_mark = transport.homojunction_vbi(gaas, N_mark, N_mark, 300.0)

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.plot(N, vbi, color="#1E2761", linewidth=2.5,
            label=r"$V_{bi} = (kT/q)\,\ln(N_A N_D / n_i^2)$, GaAs, 300 K")
    ax.plot(N_mark, vbi_mark, "o", color="#c43d3d", markersize=12, zorder=5)
    ax.annotate(f"$N_A = N_D = 10^{{17}}$ cm$^{{-3}}$\n$V_{{bi}}$ = {vbi_mark:.4f} V",
                xy=(N_mark, vbi_mark), xytext=(2e17, vbi_mark - 0.18),
                fontsize=14, fontweight="bold", color="#c43d3d",
                arrowprops=dict(arrowstyle="->", color="#c43d3d", lw=1.5))

    ax.set_xscale("log")
    ax.set_xlabel(r"doping $N_A = N_D$ (cm$^{-3}$)", fontsize=16)
    ax.set_ylabel(r"built-in voltage $V_{bi}$ (V)", fontsize=16)
    ax.set_title("GaAs homojunction built-in voltage vs. doping (300 K)", fontsize=16)
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
