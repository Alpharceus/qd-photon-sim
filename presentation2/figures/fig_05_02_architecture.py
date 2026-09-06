"""Section 5 architecture map; deterministic schematic, not a data claim."""
from pathlib import Path
import matplotlib.pyplot as plt

out = Path(__file__).parent / "out" / "05_02_architecture.png"
out.parent.mkdir(exist_ok=True)
names = ["materials", "dot_levels", "linewidth", "transport", "waveguide", "cw_g2", "qd_gf", "spectral", "loading", "device", "cards", "run_rt_edge"]
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150, facecolor="white")
ax.set_axis_off()
for i, name in enumerate(names):
    x, y = 0.04 + (i % 4)*.24, .76 - (i//4)*.29
    ax.text(x+.09, y+.09, name, ha="center", va="center", fontsize=16, weight="bold",
            bbox=dict(boxstyle="round,pad=.45", fc="#e8f1fb", ec="#286090"), transform=ax.transAxes)
    if i < len(names)-1:
        j=i+1; x2,y2=.04+(j%4)*.24+.09,.76-(j//4)*.29+.09
        ax.annotate("", (x2,y2), (x+.18,y+.09), xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="->", color="#607080", lw=1.4))
ax.text(.5,.04,"Inputs and models flow to card evaluation, sweep, evidence ledger, and verdict", ha="center", fontsize=16, transform=ax.transAxes)
fig.savefig(out, dpi=150, facecolor="white")
