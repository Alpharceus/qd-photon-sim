"""Figure 04-03: Dielectric Slab Waveguide Modes and Dispersion for HKUST Stack.

Uses fsim_core.waveguide.slab_modes and hkust_ridge_stack.
Textbook reference: Chuang, Physics of Photonic Devices, ch. 7;
Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated Circuits, ch. 7.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure root directory in sys.path
root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.waveguide import hkust_ridge_stack, slab_modes, Layer

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_03_slab_modes.png"

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

# 1. Evaluate HKUST stack at 668 nm
lambda_nm = 668.0
stack = hkust_ridge_stack(lambda_nm)
modes = slab_modes(stack, lambda_nm, pol="TE", n_modes=1)
m0 = modes[0]

# Extract z and index profile
z_nm = m0.z_nm
field = m0.field
intensity = field**2 / np.max(field**2)

# Reconstruct index step function along z
bounds = np.r_[0.0, np.cumsum([L.thickness_nm for L in stack])]
n_profile = np.zeros_like(z_nm)
for i, L in enumerate(stack):
    mask = (z_nm >= bounds[i]) & (z_nm <= bounds[i+1])
    n_profile[mask] = L.n

# Center z so core center is at z = 0
z_center = 0.5 * (bounds[0] + bounds[-1])
z_um = (z_nm - z_center) * 1e-3

# Panel 1: Vertical Index profile and Mode field / intensity
ax1_twin = ax1.twinx()

p1, = ax1.plot(z_um, n_profile, color="#2B6CB0", lw=2.5, label="Index $n(z)$")
ax1.axhline(m0.n_eff, color="#E53E3E", ls="--", lw=2, label=f"$n_{{eff}} = {m0.n_eff:.4f}$")

# Shade layers
core_half = 0.150  # ~150 nm
ax1.axvspan(-core_half, core_half, color="#EDF2F7", alpha=0.6, label="Core (AlGaInP)")
ax1.axvspan(-0.002, 0.002, color="#FEB2B2", alpha=0.9, label="QD Layer (InP, 2 nm)")

p2, = ax1_twin.plot(z_um, intensity, color="#DD6B20", lw=2.5, label=r"Intensity $|E(z)|^2$")
ax1_twin.fill_between(z_um, 0, intensity, color="#FBD38D", alpha=0.3)

ax1.set_xlabel(r"Vertical Coordinate $z$ ($\mu$m, relative to core center)")
ax1.set_ylabel("Refractive Index $n$", color="#2B6CB0")
ax1_twin.set_ylabel(r"Normalized Intensity $|E(z)|^2$", color="#DD6B20", fontsize=13, labelpad=10)
ax1.set_ylim(2.95, 3.65)
ax1_twin.set_ylim(0, 1.05)
ax1.set_xlim(-0.8, 0.8)
ax1.grid(True, linestyle=":", alpha=0.5)

dot_gamma = m0.gamma_layer("dot")
ax1.set_title(rf"HKUST Slab Mode: $n_{{eff}} = {m0.n_eff:.4f}$" + "\n" +
              rf"$\Gamma_{{dot}} = {dot_gamma*100:.2f}\%$", pad=10)

lines = [p1, ax1.get_lines()[1], p2]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc="upper right", framealpha=0.9)

# Panel 2: Dispersion vs Core Thickness d_core
# Core is split into two halves around the 4 nm well/dot stack
core_thicknesses_nm = np.linspace(60, 600, 45)
neff_vals = []
cutoff_index = 3.05  # cladding index

for d in core_thicknesses_nm:
    half_d = d / 2.0
    stk = [
        Layer("lower_cladding", 3.05, 1000),
        Layer("core_lower", 3.22, half_d),
        Layer("well_lower", 3.55, 1),
        Layer("dot", 3.50, 2, True),
        Layer("well_upper", 3.55, 1),
        Layer("core_upper", 3.22, half_d),
        Layer("upper_cladding", 3.05, 1000)
    ]
    try:
        md = slab_modes(stk, lambda_nm, "TE", 1)
        neff_vals.append(md[0].n_eff)
    except ValueError:
        neff_vals.append(np.nan)

ax2.plot(core_thicknesses_nm, neff_vals, color="#3182CE", lw=2.8, label=r"Fundamental $\mathrm{TE}_0$ Mode")
ax2.axhline(3.22, color="#4A5568", ls=":", lw=2, label=r"Core Index $n_{core} = 3.22$")
ax2.axhline(cutoff_index, color="#E53E3E", ls=":", lw=2, label=r"Cladding Index $n_{clad} = 3.05$ (Cutoff)")

# Mark HKUST nominal design (d_core = 296 nm)
ax2.plot(296, m0.n_eff, "o", color="#D69E2E", markersize=10, zorder=5)
ax2.annotate(rf"HKUST: $d=296$ nm, $n_{{eff}}={m0.n_eff:.4f}$" + "\n" +
             r"Symmetric 3.22/3.05: $n_{eff}=3.154888$" + "\n" +
             r"($k_z d - 2\arctan(\gamma/k_z) = 0$)",
             xy=(296, m0.n_eff), xytext=(200, 3.095),
             arrowprops=dict(arrowstyle="->", color="#D69E2E", lw=2),
             fontsize=10, fontweight="bold", color="#744210",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FEFCBF", edgecolor="#ECC94B"))

ax2.set_xlabel(r"Total Core Thickness $d_{core}$ (nm)")
ax2.set_ylabel(r"Effective Index $n_{eff}$", fontsize=13, labelpad=10)
ax2.set_title(r"Slab Dispersion Relation $n_{eff}(d_{core})$" + "\n" + r"at $\lambda = 668$ nm", pad=10)
ax2.set_xlim(50, 600)
ax2.set_ylim(3.03, 3.24)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="lower right", framealpha=0.9, fontsize=11)

fig.subplots_adjust(top=0.86, bottom=0.16, left=0.10, right=0.93, wspace=0.55)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
