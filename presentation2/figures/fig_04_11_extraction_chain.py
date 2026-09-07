"""Figure 04-11: Full Extraction Chain Decomposition and Acceptance Sweep Levers.

Data sourced directly from out/rt_edge/verdict.md and cards/edge-inp-gainp-design.yaml.
Contract reference: docs/rt_edge_contract.md.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_11_extraction_chain.png"

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "figure.titlesize": 18,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600/150, 900/150), dpi=150)
fig.patch.set_facecolor("white")

# Panel 1: Extraction Chain Sub-Factor Waterfall
# Sub-factors for the best diagnostic row from out/rt_edge/verdict.md:
# beta = 0.0292412 (2.924%)
# combined facet factor = 0.935012 (93.501%)
# propagation = 0.882497 (88.250%)
# NA collection = 0.857853 (85.785%)
# The cumulative product after all four factors IS eta_total (0.0206985) --
# a separate 5th "Total" bar at the same height as the 4th made their value
# labels sit on top of each other (round-2 review); dropped as redundant.
factors = [r"Waveguide $\beta$" + "\n(2.92%)",
           r"HR Mirror" + "\n" + r"$\eta_{facet}$" + "\n(93.5%)",
           r"Propagation" + "\n" + r"$\eta_{prop}$" + "\n(88.2%)",
           r"Lens Collect $\eta_{NA}$" + "\n" + r"$= \eta_{total}$" + "\n(85.8%)"]

# Values for cumulative product
vals = [0.0292412,
        0.0292412 * 0.935012,
        0.0292412 * 0.935012 * 0.882497,
        0.0292412 * 0.935012 * 0.882497 * 0.857853]

pcts = np.array(vals) * 100
x = np.arange(len(factors))

bars = ax1.bar(x, pcts, width=0.55, color=["#3182CE", "#4299E1", "#63B3ED", "#276749"], edgecolor="#2D3748", lw=1.5)

for i, (b, val) in enumerate(zip(bars, pcts)):
    ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.08,
             f"{val:.3f}%", ha="center", va="bottom", fontsize=12, fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels(factors, fontsize=12)
ax1.set_ylabel(r"Cumulative Out-Coupling Efficiency (%)")
ax1.set_title("Multiplicative Extraction Chain\n" + r"$\eta_{total}$ Decomposition", pad=10)
ax1.set_ylim(0, 3.5)
ax1.grid(True, axis="y", linestyle=":", alpha=0.5)

# Panel 2: Sweep Levers and Collected Photon Flux Progression (at 300 K)
# Base: NA=0.5, R_back=0, L=500 -> 552.6 photons/s (< 1000 floor)
# + HR Mirror (R_back=0.95): 552.6 * 2.600 = 1436.5 photons/s
# + Short Cavity (L=250 um): 1436.5 * 1.133 = 1627.8 photons/s
# + High NA (NA=0.80): 1627.8 * (0.857853 / 0.577843) = 2416.6 photons/s
# Short two-line category labels -- the full parameter values (R_back, L,
# NA) are already given in the bullets/notes and in the panel 1 breakdown;
# cramming them into the tick labels too made adjacent labels bleed into
# each other (round-2 review: "all nine x tick labels overlap").
steps = [
    "Baseline\n(uncoated)",
    "+ HR Mirror\n($2.60\\times$)",
    "+ Short Cavity\n($1.133\\times$)",
    "+ High-NA Lens\n($1.485\\times$)",
]
flux_vals = [552.6, 1436.5, 1627.8, 2416.6]

x2 = np.arange(len(steps))
colors = ["#E53E3E", "#DD6B20", "#D69E2E", "#38A169"]

bars2 = ax2.bar(x2, flux_vals, width=0.55, color=colors, edgecolor="#2D3748", lw=1.5)

# Add flux floor line at 1000 photons/s
ax2.axhline(1000.0, color="#C53030", ls="--", lw=2.5, label="Eligibility Floor [A] (1000 photons/s)")

for b, val in zip(bars2, flux_vals):
    lbl = f"{val:.0f} s$^{{-1}}$"
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 45,
             lbl, ha="center", va="bottom", fontsize=12, fontweight="bold")

# Annotate margin on the winning corner (kept fully inside the y-range and
# away from the legend -- the previous version ran past the bar it pointed
# to and past the axes into the title, and sat under the legend).
ax2.annotate("Passes flux floor:\nmargin 2.42x at 300 K",
             xy=(3, 2416.6), xytext=(0.85, 2750),
             arrowprops=dict(arrowstyle="->", color="#22543D", lw=2),
             fontsize=12, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

ax2.set_xticks(x2)
ax2.set_xticklabels(steps, fontsize=12)
ax2.set_ylabel("Collected Pulsed Flux (photons / s)")
ax2.set_title("Acceptance Sweep Levers at\n300 K vs Flux Floor", pad=10)
ax2.set_ylim(0, 3350)
ax2.grid(True, axis="y", linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9, fontsize=12)

fig.subplots_adjust(top=0.84, bottom=0.22, left=0.08, right=0.96, wspace=0.32)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
