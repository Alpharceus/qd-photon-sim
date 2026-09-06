"""Figure 03.09: diode I-V regimes (diffusion/SRH/high-injection + series
resistance rollover) and the mesa-vs-aperture current-crowding ratio.

Left panel: J(V_j) for fsim_core.transport.red_diode_preset() (a
representative repository diode preset, used here for its generic I-V shape
only -- see slide 03-09 notes), plotted on log-J vs V to show the
diffusion/SRH knee and the terminal I-V (including R_s) rolling over from
exponential toward linear at high current.

Right panel: the mesa (default Diode.area_um2) vs. current-aperture
(0.4 micron diameter, the edge-emitter cards' own aperture) area ratio that
was the source of the 2026-09-05 council-review current-density
inconsistency (module docstring, item 3 of transport.py / fig caption).

Run: python presentation2/figures/fig_03_09_iv_regimes.py
Writes presentation2/figures/out/03_09_iv_regimes.png at 1600x900 px, dpi=150.
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

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_09_iv_regimes.png"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    diode = transport.red_diode_preset()
    T = 300.0

    V_j = np.linspace(0.05, 2.2, 400)
    J = np.array([diode.j_of_vj(v, T) for v in V_j])
    V_applied = V_j + (J * diode.area_cm2) * diode.R_s_ohm

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")

    ax1.semilogy(V_j, J, color="#1769aa", lw=2.5, label=r"$J(V_j)$ -- junction voltage")
    ax1.semilogy(V_applied, J, color="#c43d3d", lw=2.5, ls="--",
                 label=r"$J(V_{applied})$ -- includes series resistance $R_s$")
    ax1.axvline(diode.vbi(T), color="gray", ls=":", lw=1.5)
    ax1.text(diode.vbi(T) + 0.02, J.min() * 3, r"$V_{bi}$", fontsize=13, color="gray")
    ax1.annotate("diffusion + SRH\n(exponential)", xy=(0.9, 1e0), fontsize=12,
                 color="#1769aa")
    ax1.annotate("$R_s$ rollover\n(linear-ish)", xy=(1.7, 3e2), fontsize=12,
                 color="#c43d3d")
    ax1.set_xlabel("voltage (V)", fontsize=15)
    ax1.set_ylabel(r"current density $J$ (A/cm$^2$)", fontsize=15)
    ax1.set_title("Two-diode I-V: from exponential to $R_s$-limited", fontsize=15)
    ax1.tick_params(labelsize=13)
    ax1.grid(True, alpha=0.25, which="both")
    ax1.legend(loc="lower right", fontsize=12, frameon=True)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    mesa_um2 = diode.area_um2
    aperture_um2 = float(np.pi * (0.4 / 2) ** 2)
    ratio = mesa_um2 / aperture_um2
    bars = ax2.bar(["mesa\n(1 um diam.)", "current aperture\n(0.4 um diam.)"],
                    [mesa_um2, aperture_um2], color=["#c43d3d", "#1769aa"], width=0.5)
    for b, v in zip(bars, [mesa_um2, aperture_um2]):
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f} um$^2$",
                  ha="center", fontsize=13, fontweight="bold")
    ax2.set_ylabel(r"area (um$^2$)", fontsize=15)
    ax2.set_title(f"Current-crowding fix: same area used everywhere\n"
                  f"(ratio = {ratio:.2f}x)", fontsize=15)
    ax2.tick_params(labelsize=13)
    ax2.grid(True, alpha=0.25, axis="y")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}  (mesa/aperture ratio = {ratio:.6f})")


if __name__ == "__main__":
    main()
