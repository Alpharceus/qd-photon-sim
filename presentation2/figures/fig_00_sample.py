"""presentation2/figures/fig_00_sample.py -- deterministic placeholder figure
for the schema-conformance sample section (00_sample.json). Not talk content:
exercises the "a figure script produces a real PNG the renderer embeds" path
that every real fig_0N_*.py script also follows.

Run: python presentation2/figures/fig_00_sample.py
Writes presentation2/figures/out/00_sample.png at 1600x900 px, dpi=150.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "out" / "00_sample.png"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    t = np.linspace(0.0, 4 * np.pi, 400)
    envelope = np.exp(-t / 8.0)
    signal = envelope * np.sin(t)

    fig, ax = plt.subplots(figsize=(1600 / 150, 900 / 150), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(t, signal, color="#1E2761", linewidth=2.5)
    ax.fill_between(t, signal, 0, color="#1E2761", alpha=0.12)
    ax.set_xlabel("time (arb. units)", fontsize=14)
    ax.set_ylabel("amplitude (arb. units)", fontsize=14)
    ax.set_title("Sample figure: schema-conformance fixture, not talk content", fontsize=14)
    ax.tick_params(labelsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
