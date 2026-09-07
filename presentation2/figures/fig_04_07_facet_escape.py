"""Figure 04-07: Facet Fresnel Transmission and the Mid-Ridge Ray-Series Escape Fraction.

pr-pkg6-stale-text item B4: panel 2 used to plot the DELETED escape-rate switch
(eta_facet = T_f / (T_f + (1 - R_back)), R_back=None meaning an equal 50/50
split); it now plots the current model, facet_escape_fraction's mid-ridge
ray-probability series (fsim_core/waveguide.py), which folds single-pass
propagation in and resolves R_back=None to the bare cleaved facet
(R_back = R_front), not to an equal 50/50 split.

Textbook reference: Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated
Circuits (2nd ed.), ch. 2 (mirror loss and escape fraction).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from fsim_core.waveguide import facet_transmission, facet_escape_fraction
from fsim_core.device import DeviceDesign

# alpha_cm/L_um: this card's own ridge-loss/length levers, not hardcoded
# literals, so the figure tracks the card if either is ever revisited.
_gainp = DeviceDesign.load(root_dir / "cards" / "edge-inp-gainp-design.yaml")
ALPHA_CM = _gainp.emission.alpha_cm
L_UM = _gainp.emission.L_um

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_07_facet_escape.png"

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

# Panel 1: Fresnel Reflectivity and Transmission vs Index
n_grid = np.linspace(1.0, 4.2, 300)
T_grid = np.array([facet_transmission(n) for n in n_grid])
R_grid = 1.0 - T_grid

ax1.plot(n_grid, T_grid * 100, color="#2B6CB0", lw=2.8, label=r"Power Transmission $T_f = 1 - R_f$")
ax1.plot(n_grid, R_grid * 100, color="#E53E3E", lw=2.8, linestyle="--", label=r"Power Reflectivity $R_f = \left(\frac{n-1}{n+1}\right)^2$")

# Highlight semiconductor ridge point (neff = 3.253)
neff_nom = 3.2530
T_nom = facet_transmission(neff_nom)
R_nom = 1.0 - T_nom

ax1.plot(neff_nom, T_nom * 100, "o", color="#2B6CB0", markersize=9)
ax1.plot(neff_nom, R_nom * 100, "o", color="#E53E3E", markersize=9)

ax1.annotate(rf"Semiconductor Facet ($n_{{eff}} = {neff_nom:.3f}$):" + "\n" +
             rf"$T_f = {T_nom*100:.2f}\%$, $R_f = {R_nom*100:.2f}\%$",
             xy=(neff_nom, T_nom * 100), xycoords="data",
             xytext=(0.02, 0.98), textcoords="axes fraction",
             arrowprops=dict(arrowstyle="->", color="#2B6CB0", lw=2),
             fontsize=12, fontweight="bold", color="#1A365D", ha="left", va="top",
             wrap=True, annotation_clip=True,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EBF8FF", edgecolor="#3182CE"))

ax1.set_xlabel(r"Effective Index $n_{eff}$")
ax1.set_ylabel("Power Fraction (%)")
ax1.set_title("Cleaved Semiconductor\nFacet Transmission", pad=10)
ax1.set_xlim(1.0, 4.2)
ax1.set_ylim(0, 102)
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="lower left", framealpha=0.9, fontsize=12)

# Panel 2: mid-ridge ray-series escape fraction vs Back-Facet Reflectivity R_back
Rb_grid = np.linspace(0.001, 0.999, 300)
# Current model (fsim_core.waveguide.facet_escape_fraction): ray-probability
# series with the dot at mid-ridge, single-pass propagation folded in [DR].
eta_ray = np.array([facet_escape_fraction(T_nom, rb, ALPHA_CM, L_UM) for rb in Rb_grid])

# Coldren & Corzine front-facet fraction:
# F1 = (1-R1) sqrt(R2) / [(1-R1) sqrt(R2) + (1-R2) sqrt(R1)]
F1_grid = ((1.0 - R_nom) * np.sqrt(Rb_grid) /
           ((1.0 - R_nom) * np.sqrt(Rb_grid) + (1.0 - Rb_grid) * np.sqrt(R_nom)))

ax2.plot(Rb_grid, eta_ray * 100, color="#38A169", lw=2.8,
         label=r"Mid-Ridge Ray Series: $\eta_{facet}(R_{back})$ [DR]")
ax2.plot(Rb_grid, F1_grid * 100, color="#805AD5", lw=2.0, linestyle="--",
         label=r"Coldren & Corzine $F_1$ (cavity photon split)")

# Mark R_back = 0.95 nominal point (HR-coated back facet, this card's own value)
Rb_card = 0.95
eta_card = facet_escape_fraction(T_nom, Rb_card, ALPHA_CM, L_UM)
F1_card = ((1.0 - R_nom) * np.sqrt(Rb_card) /
           ((1.0 - R_nom) * np.sqrt(Rb_card) + (1.0 - Rb_card) * np.sqrt(R_nom)))

# Mark the bare, uncoated cleaved facet: R_back=None resolves to R_back=R_front
# (waveguide.edge_emission), NOT to a flat 0.5*T_f equal-split line -- it is
# one point on this SAME ray-series curve, at R_back = R_f.
eta_bare = facet_escape_fraction(T_nom, R_nom, ALPHA_CM, L_UM)
gain = eta_card / eta_bare

ax2.plot(Rb_card, eta_card * 100, "o", color="#38A169", markersize=10, zorder=5)
ax2.plot(Rb_card, F1_card * 100, "^", color="#805AD5", markersize=9, zorder=5)

ax2.annotate(rf"HR Mirror ($R_{{back}} = {Rb_card:.2f}$):" + "\n" +
             rf"Ray series: $\eta_{{facet}} = {eta_card*100:.2f}\%$" + "\n" +
             rf"Coldren $F_1 = {F1_card*100:.2f}\%$" + "\n" +
             rf"({gain:.2f}$\times$ Boost over Bare Facet)",
             xy=(Rb_card, eta_card * 100), xytext=(0.28, 62),
             arrowprops=dict(arrowstyle="->", color="#38A169", lw=2),
             fontsize=12, fontweight="bold", color="#22543D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FFF4", edgecolor="#38A169"))

# Mark R_back = R_front: the bare, uncoated cleaved back facet (R_back=None)
ax2.plot(R_nom, eta_bare * 100, "s", color="#E53E3E", markersize=9, zorder=5)
ax2.annotate(rf"Bare Facet ($R_{{back}}=R_f={R_nom:.3f}$)" + "\n" +
             r"$R_{back}=$None$\to R_{back}=R_{front}$" + "\n" +
             rf"$\eta_{{facet}} = {eta_bare*100:.2f}\%$",
             xy=(R_nom, eta_bare * 100), xytext=(0.36, 25),
             arrowprops=dict(arrowstyle="->", color="#E53E3E", lw=1.8),
             fontsize=12, fontweight="bold", color="#742A2A",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FED7D7", edgecolor="#FEB2B2"))

ax2.set_xlabel(r"Back Facet Reflectivity $R_{back}$")
ax2.set_ylabel(r"Front Escape Efficiency $\eta_{facet}$ (%)")
ax2.set_title("Front Facet Out-Coupling vs\nBack Mirror Reflectivity", pad=10)
ax2.set_xlim(-0.02, 1.02)
ax2.set_ylim(20, 102)
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper left", framealpha=0.9, fontsize=12)

fig.subplots_adjust(top=0.84, bottom=0.13, left=0.08, right=0.95, wspace=0.30)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
