"""Figure 04-05: Spontaneous Emission Beta Factor vs Mode Area and Position.

Uses fsim_core.waveguide.beta_factor and edge_emission.
References:
- Lecamp, Lalanne, Hugonin, PRL 99, 023902 (2007)
- Arcari et al., PRL 113, 093603 (2014)
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.waveguide import beta_factor

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_05_beta_factor.png"

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

from fsim_core.materials import refractive_index

# Panel 1: Beta factor vs Mode Area for different group indices
# F_wg = (3 / (4 * pi)) * (lambda / n_dot)^2 * (n_g / n_dot) / A_mode
# beta = F_wg / (1 + F_wg)
lambda_nm = 770.0
n_dot = float(refractive_index("InP", lambda_nm))  # 3.418
A_grid = np.logspace(-2.5, 0.5, 300)  # um^2

ng_cases = [
    (4.355, r"Ridge Waveguide ($n_g \approx 4.36$)", "#2B6CB0", "-"),
    (15.0, r"Slow-Light Ridge ($n_g = 15$)", "#DD6B20", "--"),
    (40.0, r"PhC Waveguide ($n_g = 40$)", "#38A169", "-."),
]

for ng, label, col, ls in ng_cases:
    betas = [beta_factor(A, lambda_nm, n_dot, ng, 1.0)[1] for A in A_grid]
    ax1.plot(A_grid, np.array(betas) * 100, label=label, color=col, lw=2.5, linestyle=ls)

# Mark Ridge Waveguide nominal point
A_ridge = 0.5124  # um^2 (gainp fallback card)
F_ridge, b_ridge = beta_factor(A_ridge, lambda_nm, n_dot, 4.355, 1.0)
ax1.plot(A_ridge, b_ridge * 100, "o", color="#2B6CB0", markersize=10, zorder=5)
ax1.annotate(rf"Fallback $2\ \mu$m Ridge" + "\n" + rf"$\beta = {b_ridge*100:.2f}\%$" + "\n" + rf"($A = {A_ridge:.3f}\ \mu\mathrm{{m}}^2$)",
             xy=(A_ridge, b_ridge * 100), xytext=(0.08, 12),
             arrowprops=dict(arrowstyle="->", color="#2B6CB0", lw=2),
             fontsize=12, fontweight="bold", color="#1A365D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EBF8FF", edgecolor="#3182CE"))

# Mark Arcari et al. 2014 PhC Waveguide point (reported beta = 98.43%)
A_phc = 0.040  # um^2
beta_arcari = 98.43
ax1.plot(A_phc, beta_arcari, "s", color="#38A169", markersize=10, zorder=5)
ax1.annotate(r"Arcari et al. (PRL 2014)" + "\n" + rf"$\beta = {beta_arcari:.2f}\%$ (PhC W1)" + "\n" + r"($n_g \approx 30\mathrm{-}40,\ A \approx 0.04\ \mu\mathrm{m}^2$)",
             xy=(A_phc, beta_arcari), xytext=(0.005, 68),
             arrowprops=dict(arrowstyle="->", color="#38A169", lw=2),
             fontsize=12, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

ax1.set_xscale("log")
ax1.set_xlabel(r"Effective Mode Area $A_{mode}$ ($\mu\mathrm{m}^2$)")
ax1.set_ylabel(r"Guided Spontaneous Emission $\beta$ (%)")
ax1.set_title(r"$\beta$ Factor Scaling: Ridge vs Photonic Crystal (Lecamp 2007)", pad=12)
ax1.set_ylim(0, 105)
ax1.set_xlim(0.003, 3.0)
ax1.grid(True, which="both", linestyle=":", alpha=0.5)
ax1.legend(loc="center right", framealpha=0.9)

# Panel 2: Position Factor |E(z_dot)|^2 / max|E|^2 and beta vs QD Vertical Position
# Vertical profile approximated by cosine inside core (-148 to +148 nm) and exp decay in cladding
d_core_nm = 296.0
z_offsets = np.linspace(-300, 300, 200)
# Model normalized field
kz = np.pi / (d_core_nm * 1e-3)
gamma_clad = 8.0  # um^-1
field_z = []
for z in z_offsets:
    z_um = z * 1e-3
    if abs(z) <= d_core_nm / 2.0:
        field_z.append(np.cos(kz * z_um))
    else:
        dist = abs(z_um) - (d_core_nm / 2.0) * 1e-3
        field_z.append(np.cos(kz * (d_core_nm / 2.0) * 1e-3) * np.exp(-gamma_clad * dist))

field_z = np.array(field_z)
pos_factor = (field_z / np.max(field_z))**2
beta_vs_z = [beta_factor(A_ridge, lambda_nm, n_dot, 4.355, p)[1] * 100 for p in pos_factor]

ax2_twin = ax2.twinx()

p1, = ax2.plot(z_offsets, pos_factor, color="#805AD5", lw=2.5, label=r"Position Factor $\xi_{pos} = \frac{|E(z)|^2}{\max|E|^2}$")
p2, = ax2_twin.plot(z_offsets, beta_vs_z, color="#E53E3E", lw=2.5, linestyle="--", label=r"Coupling $\beta(z)$ (%)")

# Shade core region
ax2.axvspan(-d_core_nm / 2.0, d_core_nm / 2.0, color="#EDF2F7", alpha=0.6, label="Waveguide Core")
ax2.axvline(0, color="#D69E2E", ls=":", lw=2, label="Core Center (Antinode)")

ax2.set_xlabel(r"Dot Vertical Position Offset $\Delta z_{dot}$ (nm)")
ax2.set_ylabel(r"Position Overlap Factor $\xi_{pos}$", color="#805AD5")
ax2_twin.set_ylabel(r"Resulting $\beta$ Factor (%)", color="#E53E3E")
ax2.set_ylim(0, 1.05)
ax2_twin.set_ylim(0, 3.2)
ax2.set_xlim(-300, 300)
ax2.grid(True, linestyle=":", alpha=0.5)

ax2.set_title(r"Dipole Position Overlap $\xi_{pos}$ in Ridge Heterostructure", pad=12)

lines = [p1, p2]
labels = [l.get_label() for l in lines]
ax2.legend(lines, labels, loc="lower center", framealpha=0.9)

fig.subplots_adjust(top=0.91, bottom=0.13, left=0.08, right=0.92, wspace=0.32)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
