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

Physics review (2026-09-06, item 17): the left panel's R_s rollover is
real for transport.red_diode_preset() (used here only for its generic I-V
shape), but it appears at ~1e5-1e6 A/cm^2 -- five to six orders of
magnitude above the ~1.5 A/cm^2 current density either edge-inp design
card actually runs at (nanoamps over the 0.4 um aperture). The left panel
now marks that operating point honestly rather than implying the rollover
is reachable at this project's own currents.

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
from fsim_core.device import DeviceDesign  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_09_iv_regimes.png"
CARD = Path(__file__).resolve().parents[2] / "cards" / "edge-inp-gainp-design.yaml"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    diode = transport.red_diode_preset()
    T = 300.0

    # V_j is capped well below V_bi (2.10 V): beyond ~1.95 V, J grows so
    # steeply that I*R_s alone reaches kV-scale voltages no real mesa
    # survives (that failure regime is regimes IV, the previous slide, not
    # this one) -- this range covers J from ~1e-12 to ~4.5e6 A/cm^2, ample
    # to show both the exponential knee and the R_s rollover on one axis.
    V_j = np.linspace(0.05, 1.95, 400)
    J = np.array([diode.j_of_vj(v, T) for v in V_j])
    V_applied = V_j + (J * diode.area_cm2) * diode.R_s_ohm

    # Item 17: the card's own operating current density, for an honest
    # marker of how far this project's actual regime sits from the R_s
    # rollover shown above (fresh DeviceDesign.load, never a stale number).
    design = DeviceDesign.load(CARD)
    aperture_cm2 = np.pi * (design.aperture.diameter_um / 2.0) ** 2 * 1e-8
    J_op_Acm2 = design.drive.I_uA * 1e-6 / aperture_cm2

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")

    ax1.semilogy(V_j, J, color="#1769aa", lw=2.5, label=r"$J(V_j)$ -- junction voltage")
    ax1.semilogy(V_applied, J, color="#c43d3d", lw=2.5, ls="--",
                 label=r"$J(V_{applied})$ -- includes series resistance $R_s$")
    ax1.axhline(J_op_Acm2, color="#2e8b57", ls=":", lw=1.8)
    ax1.annotate(f"this project's own operating point:\n{J_op_Acm2:.2g} A/cm$^2$ "
                 f"(~{J[-1] / J_op_Acm2:.0e}x below the axis top)",
                 xy=(2.9, J_op_Acm2), xytext=(1.75, 3e-11), fontsize=10.5,
                 color="#2e8b57", arrowprops=dict(arrowstyle="->", color="#2e8b57", lw=1.2))
    ax1.annotate("diffusion + SRH\n(exponential)", xy=(0.62, 3e-7), xytext=(1.15, 1e-4),
                 fontsize=12, color="#1769aa",
                 arrowprops=dict(arrowstyle="->", color="#1769aa", lw=1.2))
    ax1.annotate("$R_s$ rollover (bends right):\nreal for this GENERIC preset,\n"
                 "but only at $10^5$-$10^6\\times$\nthis project's own current density",
                 xy=(2.05, 2e5), xytext=(2.15, 3e0), fontsize=10.5, color="#c43d3d",
                 arrowprops=dict(arrowstyle="->", color="#c43d3d", lw=1.2))
    ax1.set_xlim(0, 3.8)
    ax1.set_xlabel("voltage (V)", fontsize=15)
    ax1.set_ylabel(r"current density $J$ (A/cm$^2$)", fontsize=15)
    ax1.set_title("Two-diode I-V: from exponential to $R_s$-limited", fontsize=15)
    ax1.tick_params(labelsize=13)
    ax1.grid(True, alpha=0.25, which="both")
    ax1.legend(loc="upper left", fontsize=12, frameon=True)
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
    ax2.set_title(f"Current-crowding fix\n(same area everywhere; ratio = {ratio:.2f}x)",
                  fontsize=13)
    ax2.tick_params(labelsize=13)
    ax2.grid(True, alpha=0.25, axis="y")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}  (mesa/aperture ratio = {ratio:.6f}, "
          f"J_op = {J_op_Acm2:.6g} A/cm^2)")


if __name__ == "__main__":
    main()
