"""Figure 01.20: Aperture spatial filtering and Poisson single-dot statistics.
Demonstrates why mean occupation N = n*A = 0.377 yields high single-dot fidelity (>82%).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_20_aperture_poisson.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Design card operating point: density = 3e8 cm^-2, diameter = 0.4 um
# Area = pi * (0.2 um)^2 = 0.12566 um^2 = 1.2566e-9 cm^2
# N = 3e8 * 1.2566e-9 = 0.377
N_design = 0.37699

# Left panel: Poisson distribution at N = 0.377
k_vals = np.arange(0, 5)
probs = [poisson.pmf(k, N_design) for k in k_vals]
colors = ['#757575', '#2e7d32', '#c43d3d', '#d32f2f', '#b71c1c']

bars = ax1.bar(k_vals, probs, color=colors, width=0.55, edgecolor='black', alpha=0.85)

for k, p in zip(k_vals, probs):
    ax1.text(k, p + 0.02, f'{p * 100:.1f}%', ha='center', fontsize=13, fontweight='bold')

ax1.set_xlabel('Number of QDs under Aperture ($k$)', fontsize=15)
ax1.set_ylabel('Probability $P(k)$', fontsize=15)
ax1.set_title(r'Poisson Distribution ($N = nA = 0.377$)', fontsize=16, fontweight='bold')
ax1.set_xticks(k_vals)
ax1.set_xticklabels(['0 (Dark)', '1 (Single QD)', '2 (Pair)', '3 (Cluster)', '4+'], fontsize=13)
ax1.set_ylim(0, 0.85)
ax1.grid(True, axis='y', alpha=0.25)

p_single_cond = probs[1] / (1.0 - probs[0]) * 100
ax1.annotate(rf'Conditional Single QD Purity:' + '\n' + rf'$P(k=1 \mid k \geq 1) = {p_single_cond:.1f}\%$',
             xy=(1, probs[1]), xytext=(1.8, 0.45),
             fontsize=13, fontweight='bold', color='#2e7d32',
             arrowprops=dict(arrowstyle='->', color='#2e7d32', lw=2))

# Right panel: Expected QD count vs aperture diameter
d_um = np.linspace(0.1, 1.0, 200)
A_cm2 = (np.pi * (d_um / 2)**2) * 1e-8

densities = [(1e8, '1e8 cm$^{-2}$', '#1769aa', ':'),
             (3e8, '3e8 cm$^{-2}$ (Design Card)', '#2e7d32', '-'),
             (1e9, '1e9 cm$^{-2}$', '#c43d3d', '--')]

for dens, label, col, ls in densities:
    N_curve = dens * A_cm2
    ax2.plot(d_um, N_curve, label=label, color=col, lw=2.5, linestyle=ls)

ax2.plot(0.4, N_design, 'o', color='#2e7d32', markersize=10, zorder=5)
ax2.annotate(r'$d = 0.4\ \mu\mathrm{m} \rightarrow N = 0.377$',
             xy=(0.4, N_design), xytext=(0.48, 0.45),
             fontsize=13, fontweight='bold', color='#2e7d32',
             arrowprops=dict(arrowstyle='->', color='#2e7d32', lw=2))

ax2.axhline(0.5, color='gray', linestyle='--', alpha=0.7, label=r'Single-dot threshold ($N < 0.5$)')
ax2.set_xlabel(r'Aperture Diameter $d$ ($\mu$m)', fontsize=15)
ax2.set_ylabel('Mean QD Count $N = n A$', fontsize=15)
ax2.set_title('Aperture Sizing vs Areal QD Density', fontsize=16, fontweight='bold')
ax2.set_xlim(0.1, 1.0)
ax2.set_ylim(0, 2.5)
ax2.grid(True, alpha=0.25)
ax2.legend(loc='upper left', fontsize=12, frameon=True)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
