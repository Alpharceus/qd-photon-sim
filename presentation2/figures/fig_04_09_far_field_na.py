"""Figure 04-09: Gaussian Far-Field Emission and Objective NA Collection.

Uses fsim_core.waveguide.na_collection.
Textbook reference: Chuang, Physics of Photonic Devices, ch. 8;
Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated Circuits (2nd ed.), ch. 2.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.waveguide import na_collection

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_09_far_field_na.png"

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

# Parameters from fallback card (edge-inp-gainp-design)
wx_um = 0.82957155
wy_um = 0.19660935
lambda_nm = 770.0
lam_um = lambda_nm * 1e-3

# Divergence half-angles in radians and degrees (paraxial theta_div = lambda / (pi * w0))
tx_rad = lam_um / (np.pi * wx_um)
ty_rad = lam_um / (np.pi * wy_um)
tx_deg = np.degrees(tx_rad)
ty_deg = np.degrees(min(ty_rad, np.pi/2))

# Panel 1: Far-field angular cross-sections along x and y axes
theta_deg = np.linspace(-85, 85, 300)
theta_rad = np.radians(theta_deg)

# Gaussian far-field intensity: I(theta) ~ exp(-2 * (sin(theta)/sin(tx))^2)
Ix = np.exp(-2.0 * (np.sin(theta_rad) / np.sin(min(tx_rad, np.pi/2)))**2)
Iy = np.exp(-2.0 * (np.sin(theta_rad) / np.sin(min(ty_rad, np.pi/2)))**2)

ax1.plot(theta_deg, Ix, color="#3182CE", lw=2.8, label=rf"Lateral $\theta_x$ ($\theta_{{div}} \approx {tx_deg:.1f}^\circ$)")
ax1.plot(theta_deg, Iy, color="#E53E3E", lw=2.8, linestyle="--", label=rf"Vertical $\theta_y$ ($\theta_{{div}} \approx {np.degrees(ty_rad):.1f}^\circ$)")

# Mark acceptance angles for NA = 0.5, 0.75, 0.80
# theta_max = arcsin(NA)
na_cases = [(0.50, "#4A5568", ":"), (0.75, "#DD6B20", "-."), (0.80, "#38A169", "-")]

for na_val, col, ls in na_cases:
    th_acc = np.degrees(np.arcsin(na_val))
    ax1.axvline(-th_acc, color=col, ls=ls, lw=1.8)
    ax1.axvline(th_acc, color=col, ls=ls, lw=1.8, label=rf"$\mathrm{{NA}} = {na_val:.2f}$ Cone ($\pm {th_acc:.1f}^\circ$)")

ax1.set_xlabel(r"Emission Angle $\theta$ (degrees)")
ax1.set_ylabel("Normalized Far-Field Intensity")
ax1.set_title("Astigmatic Far-Field Angular Divergence", pad=12)
ax1.set_xlim(-85, 85)
ax1.set_ylim(0, 1.05)
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="upper right", framealpha=0.9, fontsize=11)

# Panel 2: NA Collection efficiency vs NA
na_grid = np.linspace(0.05, 0.95, 200)
eta_vals = np.array([na_collection(wx_um, wy_um, lambda_nm, na) for na in na_grid])

ax2.plot(na_grid, eta_vals * 100, color="#2B6CB0", lw=2.8, label=r"$\eta_{NA}(\mathrm{NA})$ (Circular Cone, Quad [DR])")

# Highlight NA = 0.5, 0.75, 0.8
pts = [
    (0.50, na_collection(wx_um, wy_um, lambda_nm, 0.50), "#4A5568", "Baseline Objective"),
    (0.75, na_collection(wx_um, wy_um, lambda_nm, 0.75), "#DD6B20", "Card Default NA"),
    (0.80, na_collection(wx_um, wy_um, lambda_nm, 0.80), "#38A169", "Optimal Sweep Corner"),
]

for na_pt, eta_pt, col, desc in pts:
    pct = eta_pt * 100
    ax2.plot(na_pt, pct, "o", color=col, markersize=10, zorder=5)

# Annotate NA 0.50 and 0.80
eta_05 = pts[0][1] * 100
eta_75 = pts[1][1] * 100
eta_80 = pts[2][1] * 100
gain_na = eta_80 / eta_05

ax2.annotate(rf"$\mathrm{{NA}} = 0.50$: $\eta = {eta_05:.1f}\%$",
             xy=(0.50, eta_05), xytext=(0.20, 20),
             arrowprops=dict(arrowstyle="->", color="#4A5568", lw=1.8),
             fontsize=11, fontweight="bold", color="#2D3748",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EDF2F7", edgecolor="#CBD5E0"))

msg_80 = f"$\\mathrm{{NA}} = 0.80$: $\\eta = {eta_80:.1f}\\%$\n({gain_na:.2f}$\\times$ Boost over NA 0.5)"
ax2.annotate(msg_80,
             xy=(0.80, eta_80), xytext=(0.38, 62),
             arrowprops=dict(arrowstyle="->", color="#38A169", lw=2),
             fontsize=11, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

ax2.set_xlabel(r"Objective Numerical Aperture ($\mathrm{NA}$)")
ax2.set_ylabel(r"Collected Fraction $\eta_{NA}$ (%)")
ax2.set_title(r"Collection Efficiency vs Objective $\mathrm{NA}$", pad=12)
ax2.set_xlim(0.1, 0.95)
ax2.set_ylim(0, 105)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper left", framealpha=0.9)

fig.subplots_adjust(top=0.91, bottom=0.13, left=0.08, right=0.95, wspace=0.28)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
