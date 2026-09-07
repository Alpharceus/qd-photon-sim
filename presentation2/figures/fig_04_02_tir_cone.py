"""Figure 04-02: Total Internal Reflection Cone and Surface vs Edge Emission.

Textbook reference: Chuang, Physics of Photonic Devices, ch. 7;
Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated Circuits, ch. 1.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_02_tir_cone.png"

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 18,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600/150, 900/150), dpi=150)
fig.patch.set_facecolor("white")

# Panel 1: Ray diagram of the escape cone vs trapped guided modes
n_semi = 3.5
theta_c_deg = np.degrees(np.arcsin(1.0 / n_semi))  # ~16.6 deg
theta_c_rad = np.radians(theta_c_deg)

# Draw semiconductor slab and air interface
ax1.fill_between([-4, 4], -3, 0, color="#E8EEF5", label="Semiconductor ($n=3.5$)")
ax1.fill_between([-4, 4], 0, 3, color="#FAFAFA", label="Air / Cladding ($n=1.0$)")
ax1.axhline(0, color="#1A365D", lw=2.5)

# QD dipole location at (0, -1.5)
qd_y = -1.5
ax1.scatter([0], [qd_y], color="#D9381E", s=180, zorder=5, label="Quantum Dot Dipole")

# Escape cone rays
angles_esc = np.linspace(-theta_c_deg, theta_c_deg, 5)
for ang in angles_esc:
    rad = np.radians(ang)
    x_surf = -qd_y * np.tan(rad)
    ax1.plot([0, x_surf], [qd_y, 0], color="#D9381E", lw=2, ls="-")
    sin_air = n_semi * np.sin(rad)
    sin_air = np.clip(sin_air, -1.0, 1.0)
    theta_air = np.arcsin(sin_air)
    dx_air = 2.0 * np.sin(theta_air)
    dy_air = 2.0 * np.cos(theta_air)
    ax1.plot([x_surf, x_surf + dx_air], [0, dy_air], color="#D9381E", lw=2, ls="-")

# Trapped / TIR rays
angles_tir = [-45, -30, 30, 45]
for ang in angles_tir:
    rad = np.radians(ang)
    x_surf = -qd_y * np.tan(rad)
    ax1.plot([0, x_surf], [qd_y, 0], color="#4A5568", lw=1.8, ls="--")
    dx_refl = 2.0 * np.sin(rad)
    dy_refl = -2.0 * np.cos(rad)
    ax1.plot([x_surf, x_surf + dx_refl], [0, dy_refl], color="#4A5568", lw=1.8, ls="--")

# Escape cone wedge highlight
cone_x = [-qd_y * np.tan(-theta_c_rad), 0, -qd_y * np.tan(theta_c_rad)]
cone_y = [0, qd_y, 0]
ax1.fill(cone_x, cone_y, color="#FBD38D", alpha=0.35, label=f"Escape Cone (±{theta_c_deg:.1f}°)")

ax1.text(0, -2.3, "Total Internal Reflection: ~98% trapped", ha="center", va="center",
         fontsize=13, fontweight="bold", color="#742A2A",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="#FED7D7", edgecolor="#FEB2B2"))
ax1.text(-2.55, 1.0, "Air Transmission\n$\\eta_{surf} = 2.08\\%$\n(paraxial $\\approx 2.04\\%$)", ha="center", va="center",
         fontsize=12, fontweight="bold", color="#276749",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="#C6F6D5", edgecolor="#9AE6B4"))

ax1.set_xlim(-3.5, 3.5)
ax1.set_ylim(-3, 3)
ax1.set_title("Surface Escape Cone at\nSemiconductor Interface", pad=10)
ax1.set_xlabel(r"Lateral Coordinate $x$ ($\mu$m)")
ax1.set_ylabel(r"Vertical Coordinate $z$ ($\mu$m)")
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="upper right", framealpha=0.9, fontsize=12)

# Panel 2: Extraction efficiency vs refractive index
n_vals = np.linspace(1.0, 4.0, 300)
theta_c_all = np.arcsin(1.0 / n_vals)
eta_exact = 0.5 * (1.0 - np.cos(theta_c_all))
eta_approx = 1.0 / (4.0 * n_vals**2)

ax2.plot(n_vals, eta_exact * 100, "b-", lw=2.5, label=r"Exact: $\eta = \frac{1 - \cos\theta_c}{2}$")
ax2.plot(n_vals, eta_approx * 100, "r--", lw=2, label=r"Paraxial: $\eta \approx \frac{1}{4 n^2}$")

materials = [
    ("Glass ($n=1.5$)", 1.5, "#3182CE"),
    ("GaN ($n=2.4$)", 2.4, "#805AD5"),
    ("GaAs / InP ($n=3.5$)", 3.5, "#D9381E")
]

for label, n_m, col in materials:
    tc = np.arcsin(1.0 / n_m)
    eta_m = 0.5 * (1.0 - np.cos(tc)) * 100
    ax2.plot(n_m, eta_m, "o", color=col, markersize=9)
    ax2.annotate(f"{label}\n$\\eta = {eta_m:.2f}\\%$",
                 xy=(n_m, eta_m), xytext=(n_m + 0.15, eta_m + 3),
                 arrowprops=dict(arrowstyle="->", color=col, lw=1.5),
                 fontsize=12, fontweight="bold", color=col)

ax2.set_title(r"Extraction Efficiency $\eta_{surf}$" + "\nvs Index $n$", pad=10)
ax2.set_xlabel("Substrate Refractive Index $n$")
ax2.set_ylabel("Extraction Fraction into Half-Space (%)")
ax2.set_ylim(0, 52)
ax2.set_xlim(1.0, 4.0)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9, fontsize=12)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
