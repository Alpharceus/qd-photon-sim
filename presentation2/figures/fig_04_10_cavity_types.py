"""Figure 04-10: Cavity Q vs V Parameter Space and Lemma 1 Collection Invariance.

Textbook reference: Fox, Quantum Optics (2006), ch. 8;
Michler (ed.), Quantum Dots for Quantum Information Technologies (2017);
Contract reference: docs/rt_edge_contract.md (Lemma 1).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_10_cavity_types.png"

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

# Panel 1: Cavity Architectures: Quality Factor Q vs Mode Volume V / (lambda/n)^3
# F_P = (3 / (4 * pi^2)) * (Q / V_norm) approx 0.076 * (Q / V_norm)
V_norm = np.logspace(-1, 2.5, 300)
Q_grid = np.logspace(1.5, 6.0, 300)
VV, QQ = np.meshgrid(V_norm, Q_grid)
FP = (3.0 / (4.0 * np.pi**2)) * (QQ / VV)

# Contour levels for Purcell Factor
levels = [1, 10, 100, 1000, 10000]
cs = ax1.contour(VV, QQ, FP, levels=levels, colors="#CBD5E0", linestyles="--", linewidths=1.5)
ax1.clabel(cs, inline=True, fmt=r"$F_P = %d$", fontsize=11)

# Architecture regions
# 1. Planar DBR Cavity
ax1.fill([10, 80, 80, 10], [200, 200, 1500, 1500], color="#BEE3F8", alpha=0.5)
ax1.text(30, 500, "Planar DBR\n" + r"($F_P \sim 1\mathrm{-}3$)", ha="center", va="center", fontsize=12, fontweight="bold", color="#2B6CB0")

# 2. Etched Micropillars
ax1.fill([2, 10, 10, 2], [1000, 1000, 15000, 15000], color="#C6F6D5", alpha=0.5)
ax1.text(4.5, 3500, "Micropillars\n" + r"($F_P \sim 5\mathrm{-}25$)", ha="center", va="center", fontsize=12, fontweight="bold", color="#22543D")

# 3. Photonic Crystal Defect Cavities (L3/H1)
ax1.fill([0.3, 1.5, 1.5, 0.3], [5000, 5000, 500000, 500000], color="#FED7D7", alpha=0.5)
ax1.text(0.65, 50000, "PhC Cavity\n" + r"($F_P > 100$)", ha="center", va="center", fontsize=12, fontweight="bold", color="#742A2A")

# 4. Ridge Waveguide (Continuum, unconfined along propagation)
ax1.scatter([30], [50], color="#D69E2E", s=180, zorder=6, marker="D")
ax1.annotate("Ridge Waveguide\n(1D Continuum,\n" + r"$\beta \approx 2.9\%$, Cavity OFF)",
             xy=(30, 50), xytext=(3, 40),
             arrowprops=dict(arrowstyle="->", color="#D69E2E", lw=2),
             fontsize=12, fontweight="bold", color="#744210",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FEFCBF", edgecolor="#D69E2E"))

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel(r"Normalized Mode Volume $V / (\lambda/n)^3$")
ax1.set_ylabel(r"Quality Factor $Q$")
ax1.set_title(r"Cavity Purcell Landscape: $F_P = \frac{3}{4\pi^2} \frac{Q}{V}$", pad=12)
ax1.set_xlim(0.15, 200)
ax1.set_ylim(30, 1e6)
ax1.grid(True, which="both", linestyle=":", alpha=0.4)

# Panel 2: Illustration of Lemma 1 (Collection Invariance of g2)
ax2.axis("off")

# Draw flowchart blocks
boxes = [
    (0.5, 0.82, "QD Emission Source\nPulsed / CW Electrical Drive\nExciton vs Biexciton Dynamics", "#EDF2F7", "#4A5568"),
    (0.5, 0.58, "Intrinsic Multi-Photon Purity\n" + r"$g^{(2)}_{\mathrm{intrinsic}}(0) = \frac{\langle n(n-1)\rangle}{\langle n \rangle^2} < 0.5$" + "\nDetermined SOLELY by Quantum Dot Loading", "#FEB2B2", "#9B2C2C"),
    (0.5, 0.32, "Linear Out-Coupling & Collection Chain\n" + r"$\eta_{\mathrm{total}} = \beta \times \eta_{\mathrm{facet}} \times \eta_{\mathrm{prop}} \times \eta_{\mathrm{NA}}$" + "\nRandom Bernoulli Partition (Lossy Beamsplitter)", "#EBF8FF", "#2B6CB0"),
    (0.5, 0.08, "Lemma 1 (docs/rt_edge_contract.md):\n" + r"$g^{(2)}_{\mathrm{detected}}(0) = \frac{\eta^2 \langle n(n-1)\rangle}{\eta^2 \langle n \rangle^2} \equiv g^{(2)}_{\mathrm{intrinsic}}(0)$" + "\nCollection Levers Scale FLUX ONLY, Never Purity!", "#C6F6D5", "#22543D"),
]

for bx, by, text, bg, border in boxes:
    ax2.text(bx, by, text, ha="center", va="center", fontsize=12,
             bbox=dict(boxstyle="round,pad=0.5", facecolor=bg, edgecolor=border, lw=2))

# Arrows
ax2.annotate("", xy=(0.5, 0.69), xytext=(0.5, 0.75), arrowprops=dict(arrowstyle="->", lw=2, color="#4A5568"))
ax2.annotate("", xy=(0.5, 0.44), xytext=(0.5, 0.50), arrowprops=dict(arrowstyle="->", lw=2, color="#9B2C2C"))
ax2.annotate("", xy=(0.5, 0.19), xytext=(0.5, 0.25), arrowprops=dict(arrowstyle="->", lw=2, color="#2B6CB0"))

ax2.set_title("Lemma 1: Out-Coupling vs Quantum Correlation", pad=12)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
