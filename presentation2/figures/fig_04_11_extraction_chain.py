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
# NA collection = 0.548244 (54.824%)
# Product eta_total = 0.0132281 (1.323%)
factors = [r"Waveguide $\beta$" + "\n(2.92%)",
           r"HR Mirror" + "\n" + r"$\eta_{facet}$" + "\n(93.5%)",
           r"Propagation" + "\n" + r"$\eta_{prop}$" + "\n(88.2%)",
           r"Lens Collect" + "\n" + r"$\eta_{NA}$" + "\n(54.8%)",
           r"Total Out-Coupling" + "\n" + r"$\eta_{total}$" + "\n(1.32%)"]

# Values for cumulative product
vals = [0.0292412,
        0.0292412 * 0.935012,
        0.0292412 * 0.935012 * 0.882497,
        0.0292412 * 0.935012 * 0.882497 * 0.548244,
        0.0132281]

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

# Panel 2: Sweep Levers and Collected Photon Flux Progression
# Base: NA=0.5, R_back=0, L=500 -> 286.3 photons/s
# + HR Mirror (R_back=0.95): 286.3 * (0.935012 / 0.359685) = 744.3 photons/s
# + Short Cavity (L=250 um): 744.3 * (0.882497 / 0.778801) = 843.4 photons/s
# + High NA (NA=0.80): 843.4 * (0.548244 / 0.299495) = 1544.4 photons/s
steps = [
    "Baseline Uncoated\n($R_b=0, L=500\\mu\\mathrm{m}, \\mathrm{NA}=0.5$)",
    "+ HR Back Mirror\n($R_{back} = 0.95$, $2.60\\times$)",
    "+ Shortened Cavity\n($L = 250\\mu\\mathrm{m}$, $1.13\\times$)",
    "+ High-NA Objective\n($\\mathrm{NA} = 0.80$, $1.83\\times$)",
]
flux_vals = [286.3, 744.3, 843.4, 1544.4]

x2 = np.arange(len(steps))
colors = ["#E53E3E", "#DD6B20", "#D69E2E", "#38A169"]

bars2 = ax2.bar(x2, flux_vals, width=0.55, color=colors, edgecolor="#2D3748", lw=1.5)

# Add flux floor line at 1000 photons/s
ax2.axhline(1000.0, color="#C53030", ls="--", lw=2.5, label="Eligibility Floor [A] (1000 photons/s)")

for b, val in zip(bars2, flux_vals):
    lbl = f"{val:.0f} s$^{{-1}}$"
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 35,
             lbl, ha="center", va="bottom", fontsize=12, fontweight="bold")

# Annotate margin on the winning corner
ax2.annotate("Passes Flux Floor!\n(Margin = 1.544x)",
             xy=(3, 1544.4), xytext=(2.1, 1650),
             arrowprops=dict(arrowstyle="->", color="#22543D", lw=2),
             fontsize=12, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

ax2.set_xticks(x2)
ax2.set_xticklabels(steps, fontsize=10.5)
ax2.set_ylabel("Collected Pulsed Flux (photons / s)")
ax2.set_title("Acceptance Sweep: Photonic Levers vs Flux Floor", pad=12)
ax2.set_ylim(0, 1950)
ax2.grid(True, axis="y", linestyle=":", alpha=0.5)
ax2.legend(loc="upper left", framealpha=0.9, fontsize=12)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
