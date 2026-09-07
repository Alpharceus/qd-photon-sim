"""Section 5 architecture map; deterministic schematic, not a data claim."""
from pathlib import Path
import matplotlib.pyplot as plt

out = Path(__file__).parent / "out" / "05_02_architecture.png"
out.parent.mkdir(exist_ok=True)
names = ["materials", "dot_levels", "linewidth", "transport", "waveguide", "cw_g2", "qd_gf", "spectral", "loading", "device", "cards", "run_rt_edge"]
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150, facecolor="white")
ax.set_axis_off()
centers = []
for i, name in enumerate(names):
    row = i // 4
    col = i % 4 if row % 2 == 0 else (3 - (i % 4))
    cx = 0.04 + col * 0.24 + 0.09
    cy = 0.76 - row * 0.29 + 0.09
    centers.append((cx, cy))
    ax.text(cx, cy, name, ha="center", va="center", fontsize=16, weight="bold",
            bbox=dict(boxstyle="round,pad=.45", fc="#e8f1fb", ec="#286090"), transform=ax.transAxes)

for i in range(len(names) - 1):
    cx1, cy1 = centers[i]
    cx2, cy2 = centers[i+1]
    if cy1 == cy2:
        if cx2 > cx1:
            ax.annotate("", (cx2 - 0.075, cy2), (cx1 + 0.075, cy1), xycoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color="#607080", lw=1.6))
        else:
            ax.annotate("", (cx2 + 0.075, cy2), (cx1 - 0.075, cy1), xycoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color="#607080", lw=1.6))
    else:
        ax.annotate("", (cx2, cy2 + 0.055), (cx1, cy1 - 0.055), xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="->", color="#607080", lw=1.6))
ax.text(.5, .04, "Inputs and models flow to card evaluation, sweep, evidence ledger, and verdict", ha="center", fontsize=16, transform=ax.transAxes)
fig.savefig(out, dpi=150, facecolor="white")
