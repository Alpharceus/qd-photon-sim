"""Figure 03.12: injection efficiency eta_inj(T) for both RT edge-emitter
design cards.

Sweeps thermal.T_hs on each of cards/edge-inp-gaasp-design.yaml (primary,
HKUST-class transport diode) and cards/edge-inp-gainp-design.yaml (fallback,
Reischle-class transport diode) through fresh fsim_core.device.evaluate()
calls, reading back scalars['eta_inj'] at each temperature -- the actual
injection efficiency each card's own transport.Diode + confinement model
produces, not a hand-rolled preset call (which can silently diverge from the
card's plumbing; see the section-3 spec).

Run: python presentation2/figures/fig_03_12_eta_inj_vs_T.py
Writes presentation2/figures/out/03_12_eta_inj_vs_T.png at 1600x900 px, dpi=150.
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

OUT_PATH = Path(__file__).resolve().parent / "out" / "03_12_eta_inj_vs_T.png"
CARDS_DIR = Path(__file__).resolve().parents[2] / "cards"

CARDS = [
    ("edge-inp-gaasp-design.yaml", "primary (gaasp, HKUST-class)", "#c43d3d"),
    ("edge-inp-gainp-design.yaml", "fallback (gainp, Reischle-class)", "#1769aa"),
]

T_GRID = [80.0, 120.0, 160.0, 200.0, 230.0, 260.0, 300.0]


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for fname, label, color in CARDS:
        card_path = CARDS_DIR / fname
        eta = []
        for T in T_GRID:
            design = DeviceDesign.load(card_path)
            design.thermal.T_hs = float(T)
            sc = evaluate(design, T_grid=[float(T)])["scalars"]
            eta.append(sc["eta_inj"])
        ax.plot(T_GRID, eta, marker="o", markersize=7, lw=2.5, color=color, label=label)

    ax.axhline(1.0, color="gray", ls=":", lw=1.0)
    ax.set_xlabel("heat-sink temperature $T_{hs}$ (K)", fontsize=16)
    ax.set_ylabel(r"injection efficiency $\eta_{inj}$", fontsize=16)
    ax.set_title("Injection efficiency vs. temperature: barrier height decides\n"
                 "how fast thermionic leakage grows", fontsize=15)
    ax.set_ylim(0.0, 1.05)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower left", fontsize=13, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
