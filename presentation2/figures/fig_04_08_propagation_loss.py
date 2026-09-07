"""Figure 04-08: Waveguide Propagation Loss and Cavity Length Optimization.

Textbook reference: Coldren, Corzine, Masanovic, Diode Lasers and Photonic Integrated
Circuits (2nd ed.), ch. 2; Chuang, Physics of Photonic Devices, ch. 8.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_08_propagation_loss.png"

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

# Panel 1: Transmission eta_prop = exp(-alpha * L) vs L
L_um = np.linspace(0, 1000, 300)
alphas = [
    (2.0, r"$\alpha = 2\ \mathrm{cm}^{-1}$ (Low-loss passive)", "#38A169", "-."),
    (5.0, r"$\alpha = 5\ \mathrm{cm}^{-1}$ (Class estimate [E])", "#2B6CB0", "-"),
    (10.0, r"$\alpha = 10\ \mathrm{cm}^{-1}$ (Doped cladding)", "#DD6B20", "--"),
    (20.0, r"$\alpha = 20\ \mathrm{cm}^{-1}$ (High scattering)", "#E53E3E", ":"),
]

for alpha, label, col, ls in alphas:
    # L in um -> L * 1e-4 in cm
    trans = np.exp(-alpha * (L_um * 1e-4))
    ax1.plot(L_um, trans * 100, color=col, lw=2.5, linestyle=ls, label=label)

# Highlight 250 um and 500 um points at alpha = 5 cm^-1
t_250 = np.exp(-5.0 * 250e-4) * 100
t_500 = np.exp(-5.0 * 500e-4) * 100

ax1.plot(250, t_250, "o", color="#2B6CB0", markersize=10, zorder=5)
ax1.plot(500, t_500, "s", color="#4A5568", markersize=10, zorder=5)

ax1.annotate(rf"$L = 250\ \mu\mathrm{{m}}$:" + "\n" + rf"$\eta_{{prop}} = {t_250:.1f}\%$",
             xy=(250, t_250), xytext=(120, 75),
             arrowprops=dict(arrowstyle="->", color="#2B6CB0", lw=2),
             fontsize=12, fontweight="bold", color="#1A365D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EBF8FF", edgecolor="#3182CE"))

ax1.annotate(rf"$L = 500\ \mu\mathrm{{m}}$:" + "\n" + rf"$\eta_{{prop}} = {t_500:.1f}\%$",
             xy=(500, t_500), xytext=(550, 70),
             arrowprops=dict(arrowstyle="->", color="#4A5568", lw=2),
             fontsize=12, fontweight="bold", color="#2D3748",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EDF2F7", edgecolor="#A0AEC0"))

ax1.set_xlabel(r"Cavity Length $L$ ($\mu$m)")
ax1.set_ylabel(r"Propagation Transmission $\eta_{prop} = e^{-\alpha L}$ (%)")
ax1.set_title("Internal Transmission vs\nWaveguide Length", pad=10)
ax1.set_xlim(0, 1000)
ax1.set_ylim(20, 102)
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="lower left", framealpha=0.9)

# Panel 2: Comparative Trade-Offs (250 um vs 500 um)
metrics = ["Optical Throughput\n(Transmission)", "Photon Escape\nSpeed (Round-trip)", "Bar Cleave &\nHandling Yield", "Contact Area &\nConductance"]
vals_250 = [88.2, 100.0, 70.0, 50.0]   # 250 um normalized
vals_500 = [77.9, 50.0, 95.0, 100.0]   # 500 um normalized

y = np.arange(len(metrics))
h = 0.35

ax2.barh(y + h/2, vals_250, height=h, color="#3182CE", label=r"$L = 250\ \mu\mathrm{m}$ (Fallback design)")
ax2.barh(y - h/2, vals_500, height=h, color="#A0AEC0", label=r"$L = 500\ \mu\mathrm{m}$ (Standard bar)")

# Annotate transmission gain
ax2.text(88.2 + 2, 0 + h/2, f"+13.3% Flux\n($\\eta={t_250:.1f}\\%$)", va="center", fontsize=11, fontweight="bold", color="#2B6CB0")
ax2.text(77.9 + 2, 0 - h/2, f"Baseline\n($\\eta={t_500:.1f}\\%$)", va="center", fontsize=11, color="#4A5568")

ax2.set_yticks(y)
ax2.set_yticklabels(metrics, fontsize=12)
ax2.set_xlabel("Illustrative / Qualitative Metric (%)")
ax2.set_title(r"[A] Illustrative Trade-Offs:" + "\n" + r"250 vs 500 $\mu$m Bars", pad=10)
ax2.set_xlim(0, 138)
ax2.grid(True, axis="x", linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9, fontsize=11)

fig.subplots_adjust(top=0.84, bottom=0.13, left=0.08, right=0.95, wspace=0.30)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
