"""Figure 03.22: the electron-to-photon efficiency chain at the repository's
own favourable 300 K corner (out/rt_edge/verdict.md's best diagnostic-g2 row:
gainp card, dot.delta_xx=8 meV, dot.gamma300=6 meV, emission.NA=0.8,
emission.R_back=0.95, emission.L_um=250, pulsed drive at 80 MHz / 100 ps).

Obtains loading, t_X, S, eta_total and the reported collected flux from a
fresh fsim_core.device.evaluate() call (never from a possibly-stale sweep
file) and draws a waterfall of the multiplicative chain, then performs the
same self-check verdict.md reports: product * rep_rate reproduces
collected_flux_pulsed_s.

Run: python presentation2/figures/fig_03_22_efficiency_chain.py
Writes presentation2/figures/out/03_22_efficiency_chain.png at 1600x900 px, dpi=150.
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

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_22_efficiency_chain.png"
CARD = Path(__file__).resolve().parents[2] / "cards" / "edge-inp-gainp-design.yaml"

PULSE_WIDTH_NS = 0.1
REP_RATE_HZ = 80.0e6


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    design = DeviceDesign.load(CARD)
    design.dot.delta_xx = 8.0
    design.dot.gamma300 = 6.0
    design.emission.NA = 0.8
    design.emission.R_back = 0.95
    design.emission.L_um = 250.0
    design.drive.duty = PULSE_WIDTH_NS * 1e-9 * REP_RATE_HZ
    design.drive.cw = False
    design.drive.diode["tau_pulse_ns"] = PULSE_WIDTH_NS

    sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]

    mu = sc["mu_resolved"]
    loading = 1.0 - np.exp(-mu)
    t_x = sc["t_x_op"]
    S = sc["S_resolved"]
    eta_total = sc["edge_eta_total"]
    reported_flux = sc["collected_flux_pulsed_s"]

    factors = [("loading\n$1-e^{-\\mu}$", loading), ("$t_X$\n(spectral)", t_x),
               ("$S$\n(retention)", S), (r"$\eta_{total}$" + "\n(out-coupling)", eta_total)]

    cum = 1.0
    lefts, heights, labels = [], [], []
    for name, val in factors:
        lefts.append(cum)
        heights.append(cum * val - cum)
        labels.append(name)
        cum *= val
    product_flux = cum * REP_RATE_HZ
    rel_err = abs(product_flux / reported_flux - 1.0) if reported_flux else float("nan")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150,
                                   gridspec_kw={"width_ratios": [1.4, 1]})
    fig.patch.set_facecolor("white")

    colors = ["#1769aa", "#2e8b57", "#c43d3d", "#e08214"]
    running = 1.0
    xpos = np.arange(len(factors) + 1)
    values = [1.0]
    for (_, val) in factors:
        running *= val
        values.append(running)
    ax1.bar(xpos, values, color=["#888888"] + colors, width=0.55)
    for x, v in zip(xpos, values):
        ax1.text(x, v * 1.3 if v > 0 else 1e-6, f"{v:.3g}", ha="center", fontsize=12,
                 fontweight="bold")
    ax1.set_yscale("log")
    ax1.set_xticks(xpos)
    ax1.set_xticklabels(["start\n(per pulse)"] + [f[0] for f in factors], fontsize=11)
    ax1.tick_params(axis="x", pad=8)
    ax1.set_ylabel("running product (dimensionless)", fontsize=14)
    ax1.set_title("Electron-to-photon efficiency chain\n(favourable 300 K corner, gainp card)",
                  fontsize=14)
    ax1.grid(True, alpha=0.25, which="both", axis="y")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.axis("off")
    text = (
        f"rep rate = {REP_RATE_HZ:.3g} Hz\n\n"
        f"product x rep rate =\n{product_flux:,.1f} photons/s\n\n"
        f"reported collected_flux_pulsed_s =\n{reported_flux:,.1f} photons/s\n\n"
        f"self-check relative difference:\n{rel_err*100:.4f}%  "
        f"({'PASS' if rel_err < 0.01 else 'FAIL'}, within 1%)"
    )
    ax2.text(0.02, 0.5, text, fontsize=15, va="center", family="monospace",
             bbox=dict(boxstyle="round", facecolor="#f0f0f0", edgecolor="gray"))

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}  (self-check rel. diff. = {rel_err*100:.6f}%)")


if __name__ == "__main__":
    main()
