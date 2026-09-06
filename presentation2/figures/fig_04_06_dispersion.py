"""Figure 04-06: Group Index and Waveguide Dispersion.

Uses fsim_core.waveguide and materials.refractive_index.
Textbook reference: Chuang, Physics of Photonic Devices, ch. 8;
Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated Circuits, ch. 2.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.materials import refractive_index
from fsim_core.waveguide import Layer, effective_index_ridge

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_06_dispersion.png"

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

# Wavelength range for GaInP/AlGaInP system (680 to 820 nm)
lambdas = np.linspace(680.0, 820.0, 35)

cl_mat = "(Al0.50Ga0.50)0.51In0.49P"
core_mat = "Ga0.51In0.49P"
dot_mat = "InP"

n_cl_list = []
n_core_list = []
n_dot_list = []
neff_list = []
ng_list = []

for lam in lambdas:
    n_cl = refractive_index(cl_mat, lam)
    n_core = refractive_index(core_mat, lam)
    n_dot = refractive_index(dot_mat, lam)
    n_cl_list.append(n_cl)
    n_core_list.append(n_core)
    n_dot_list.append(n_dot)
    
    stk = [
        Layer("barrier_lower", n_cl, 1000.0, material=cl_mat),
        Layer("matrix_lower", n_core, 148.0, material=core_mat),
        Layer("dot", n_dot, 3.0, True, material=dot_mat),
        Layer("matrix_upper", n_core, 148.0, material=core_mat),
        Layer("barrier_upper", n_cl, 1000.0, material=cl_mat),
    ]
    rm = effective_index_ridge(stk, lam, 2000.0, 1200.0)
    neff_list.append(rm.n_eff)

neff_arr = np.array(neff_list)
# Compute group index by numerical derivative ng = neff - lambda * dneff/dlambda
dneff_dlam = np.gradient(neff_arr, lambdas)
ng_arr = neff_arr - lambdas * dneff_dlam

# Panel 1: Material Index Dispersion
ax1.plot(lambdas, n_dot_list, color="#D69E2E", lw=2.5, label=f"Dot ({dot_mat})")
ax1.plot(lambdas, n_core_list, color="#E53E3E", lw=2.5, label=f"Core Matrix ({core_mat})")
ax1.plot(lambdas, n_cl_list, color="#3182CE", lw=2.5, label=f"Cladding Barrier (AlGaInP)")

ax1.set_xlabel(r"Wavelength $\lambda$ (nm)")
ax1.set_ylabel("Refractive Index $n$")
ax1.set_title("Constituent Material Dispersion", pad=12)
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="upper right", framealpha=0.9)
ax1.set_xlim(680, 820)

# Panel 2: Effective Index vs Group Index
ax2.plot(lambdas, neff_arr, color="#2B6CB0", lw=2.8, label=r"Phase Index $n_{eff}(\lambda)$")
ax2.plot(lambdas, ng_arr, color="#9B2C2C", lw=2.8, linestyle="--", label=r"Group Index $n_g(\lambda) = n_{eff} - \lambda \frac{dn_{eff}}{d\lambda}$")

# Highlight 770 nm design point (fallback card)
idx_770 = np.argmin(np.abs(lambdas - 770.0))
lam_pt = lambdas[idx_770]
neff_pt = neff_arr[idx_770]
ng_pt = ng_arr[idx_770]

ax2.plot(lam_pt, neff_pt, "o", color="#2B6CB0", markersize=9)
ax2.plot(lam_pt, ng_pt, "s", color="#9B2C2C", markersize=9)

enhancement = (ng_pt / neff_pt - 1.0) * 100
ax2.annotate(rf"$\lambda = 770$ nm:" + "\n" +
             rf"$n_{{eff}} = {neff_pt:.3f}$" + "\n" +
             rf"$n_g = {ng_pt:.3f}$" + "\n" +
             rf"(+ {enhancement:.1f}% LDOS boost)",
             xy=(lam_pt, ng_pt), xytext=(700, 4.35),
             arrowprops=dict(arrowstyle="->", color="#9B2C2C", lw=2),
             fontsize=12, fontweight="bold", color="#742A2A",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FED7D7", edgecolor="#E53E3E"))

# Fill the dispersion difference
ax2.fill_between(lambdas, neff_arr, ng_arr, color="#FEB2B2", alpha=0.3,
                 label=r"Dispersion Offset $-\lambda \frac{dn_{eff}}{d\lambda}$")

ax2.set_xlabel(r"Wavelength $\lambda$ (nm)")
ax2.set_ylabel("Index")
ax2.set_title(r"Ridge Waveguide Phase vs Group Index ($w = 2\ \mu\mathrm{m}$)", pad=12)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9)
ax2.set_xlim(680, 820)
ax2.set_ylim(3.15, 4.65)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
