"""Figure 04-04: Ridge Waveguide 2D Mode Profile from Effective-Index Method.

Uses fsim_core.waveguide.effective_index_ridge.
Textbook reference: Chuang, Physics of Photonic Devices, ch. 7;
Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated Circuits, ch. 7.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.waveguide import hkust_ridge_stack, effective_index_ridge

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_04_ridge_mode.png"

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

# Solve ridge mode for HKUST stack at 668 nm, width 2000 nm, etch depth 1200 nm
lambda_nm = 668.0
ridge_w_nm = 2000.0
etch_d_nm = 1200.0
stack = hkust_ridge_stack(lambda_nm)
rm = effective_index_ridge(stack, lambda_nm, ridge_w_nm, etch_d_nm)

# Lateral grid (x)
x_nm = rm.lateral.z_nm
x_center = 0.5 * (x_nm[0] + x_nm[-1])
x_um = (x_nm - x_center) * 1e-3
Ix = rm.lateral.field**2 / np.max(rm.lateral.field**2)

# Vertical grid (y)
y_nm = rm.vertical.z_nm
y_bounds = np.r_[0.0, np.cumsum([L.thickness_nm for L in stack])]
y_center = 0.5 * (y_bounds[0] + y_bounds[-1])
y_um = (y_nm - y_center) * 1e-3
Iy = rm.vertical.field**2 / np.max(rm.vertical.field**2)

# Create 2D meshgrid on cropped ROI (-2.5 to 2.5 um in x, -0.8 to 0.8 um in y)
roi_x = (x_um >= -2.5) & (x_um <= 2.5)
roi_y = (y_um >= -0.8) & (y_um <= 0.8)

x_sub = x_um[roi_x]
y_sub = y_um[roi_y]
Ix_sub = Ix[roi_x]
Iy_sub = Iy[roi_y]

I2D = np.outer(Iy_sub, Ix_sub)

# Panel 1: 2D Intensity Map with Ridge Outlines
im = ax1.imshow(I2D, extent=[x_sub[0], x_sub[-1], y_sub[0], y_sub[-1]],
                origin="lower", cmap="inferno", aspect="auto")

# Draw ridge boundary outlines
# Ridge is 2.0 um wide centered at x = 0
w_half = 0.5 * (ridge_w_nm * 1e-3)
# Upper cladding is at y > core_upper (+0.15 um) to +1.15 um
ax1.plot([-w_half, -w_half], [0.15, 0.8], color="cyan", lw=2, ls="--")
ax1.plot([w_half, w_half], [0.15, 0.8], color="cyan", lw=2, ls="--")
ax1.plot([-2.5, -w_half], [0.15, 0.15], color="cyan", lw=2, ls="--")
ax1.plot([w_half, 2.5], [0.15, 0.15], color="cyan", lw=2, ls="--")

# Colorbar
cbar = fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
cbar.set_label(r"Normalized Intensity $|E(x,y)|^2$", fontsize=12)
cbar.ax.tick_params(labelsize=11)

ax1.set_xlabel(r"Lateral Coordinate $x$ ($\mu$m)")
ax1.set_ylabel(r"Vertical Coordinate $y$ ($\mu$m)")
ax1.set_title(rf"2D Mode Intensity Profile" + "\n" + rf"($A_{{mode}} = {rm.A_mode_um2:.4f}\ \mu\mathrm{{m}}^2$)", pad=10)
ax1.set_xlim(-2.0, 2.0)
ax1.set_ylim(-0.6, 0.6)

# Annotate effective indices (two lines: the full sentence is wider than the
# panel and, since text is not clipped to its axes, one long line bled into
# the colorbar to the right)
ax1.text(0, -0.42, rf"Ridge: $n_{{ridge}} = {rm.n_ridge:.4f}$ | Etched: $n_{{outside}} = {rm.n_outside:.4f}$ [A]" + "\n" +
         rf"2D $n_{{eff}} = {rm.n_eff:.4f}$",
         color="white", ha="center", va="center", fontsize=10, fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.25", facecolor="black", alpha=0.65))

# Panel 2: 1D Cross Sections along x and y
ax2.plot(x_um, Ix, color="#3182CE", lw=2.5, label=rf"Lateral $I(x)$ ($w_x = {rm.wx_um:.3f}\ \mu\mathrm{{m}}$)")
ax2.plot(y_um, Iy, color="#E53E3E", lw=2.5, label=rf"Vertical $I(y)$ ($w_y = {rm.wy_um:.3f}\ \mu\mathrm{{m}}$)")

# 1/e^2 reference line. The plotted curves are the true (non-Gaussian) slab
# intensity profiles, not Gaussians, so they do NOT cross this level exactly
# at the Gaussian-equivalent waists w_x/w_y quoted in the legend above (those
# come from the separable effective-area definition A_mode = pi wx wy, not
# from where this raw profile crosses 1/e^2) -- label it as a reference only.
ax2.axhline(1.0 / np.e**2, color="#718096", ls=":", lw=1.8, label=r"$1/e^2 \approx 0.135$ reference level")

ax2.axvline(-w_half, color="#3182CE", ls="--", alpha=0.5)
ax2.axvline(w_half, color="#3182CE", ls="--", alpha=0.5)
ax2.text(0, 0.55, rf"Ridge Width" + "\n" + rf"$w = {ridge_w_nm*1e-3:.1f}\ \mu\mathrm{{m}}$", ha="center", va="center",
         fontsize=11, color="#2B6CB0",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="#EBF8FF", edgecolor="#BEE3F8"))
ax2.text(0, 0.28, rf"$A_{{mode}} = \pi w_x w_y = {rm.A_mode_um2:.4f}\ \mu\mathrm{{m}}^2$", ha="center", va="center",
         fontsize=11, fontweight="bold", color="#1A365D",
         bbox=dict(boxstyle="round,pad=0.2", facecolor="#EDF2F7", edgecolor="#CBD5E0"))

ax2.set_xlabel(r"Transverse Coordinate ($\mu$m)")
ax2.set_ylabel("Normalized Intensity")
ax2.set_title("Lateral vs Vertical Mode Profiles", pad=12)
ax2.set_xlim(-2.0, 2.0)
ax2.set_ylim(0, 1.05)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9, fontsize=10)

fig.subplots_adjust(top=0.86, bottom=0.13, left=0.095, right=0.95, wspace=0.42)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
