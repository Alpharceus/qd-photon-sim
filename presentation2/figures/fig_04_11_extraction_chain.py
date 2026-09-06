"""Figure 04-11: Full Extraction Chain Decomposition and Acceptance Sweep Levers.

Data sourced directly from out/rt_edge/verdict.md and cards/edge-inp-gainp-design.yaml.
Contract reference: docs/rt_edge_contract.md.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

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
# Product eta_total = 0.0206985 (2.070%)
factors = [r"Waveguide $\beta$" + "\n(2.92%)",
           r"HR Mirror" + "\n" + r"$\eta_{facet}$" + "\n(93.5%)",
           r"Propagation" + "\n" + r"$\eta_{prop}$" + "\n(88.2%)",
           r"Lens Collect" + "\n" + r"$\eta_{NA}$" + "\n(85.8%)",
           r"Total Out-Coupling" + "\n" + r"$\eta_{total}$" + "\n(2.07%)"]

# Values for cumulative product
vals = [0.0292412,
        0.0292412 * 0.935012,
        0.0292412 * 0.935012 * 0.882497,
        0.0292412 * 0.935012 * 0.882497 * 0.857853,
        0.0206985]

pcts = np.array(vals) * 100
x = np.arange(len(factors))

bars = ax1.bar(x, pcts, width=0.55, color=["#3182CE", "#4299E1", "#63B3ED", "#90CDF4", "#276749"], edgecolor="#2D3748", lw=1.5)

for i, (b, val) in enumerate(zip(bars, pcts)):
    ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.08,
             f"{val:.3f}%", ha="center", va="bottom", fontsize=12, fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels(factors, fontsize=11)
ax1.set_ylabel(r"Cumulative Out-Coupling Efficiency (%)")
ax1.set_title(r"Multiplicative Extraction Chain $\eta_{total}$ Decomposition", pad=12)
ax1.set_ylim(0, 3.5)
ax1.grid(True, axis="y", linestyle=":", alpha=0.5)

# Panel 2: Sweep Levers and Collected Photon Flux Progression (at 300 K)
# Base: NA=0.5, R_back=0, L=500 -> 552.6 photons/s (< 1000 floor)
# + HR Mirror (R_back=0.95): 552.6 * 2.600 = 1436.5 photons/s
# + Short Cavity (L=250 um): 1436.5 * 1.133 = 1627.8 photons/s
# + High NA (NA=0.80): 1627.8 * (0.857853 / 0.577843) = 2416.6 photons/s
steps = [
    "Baseline Uncoated\n($R_b=0, L=500\\mu\\mathrm{m}, \\mathrm{NA}=0.5$)",
    "+ HR Back Mirror\n($R_{back} = 0.95$, $2.60\\times$)",
    "+ Shortened Cavity\n($L = 250\\mu\\mathrm{m}$, $1.133\\times$)",
    "+ High-NA Objective\n($\\mathrm{NA} = 0.80$, $1.485\\times$)",
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

# Annotate margin on the winning corner
ax2.annotate("Passes Flux Floor!\n(Margin = 2.42x at 300 K;\n230 K reaches 6890 s$^{-1}$)",
             xy=(3, 2416.6), xytext=(1.85, 2550),
             arrowprops=dict(arrowstyle="->", color="#22543D", lw=2),
             fontsize=11, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

ax2.set_xticks(x2)
ax2.set_xticklabels(steps, fontsize=10.5)
ax2.set_ylabel("Collected Pulsed Flux (photons / s)")
ax2.set_title("Acceptance Sweep Levers at 300 K vs Flux Floor", pad=12)
ax2.set_ylim(0, 3100)
ax2.grid(True, axis="y", linestyle=":", alpha=0.5)
ax2.legend(loc="upper left", framealpha=0.9, fontsize=11)

fig.subplots_adjust(top=0.91, bottom=0.14, left=0.08, right=0.95, wspace=0.28)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
