"""Figure 01.05: Finite square well 1D solutions and BenDaniel-Duke boundary matching.
Uses fsim_core.dot_levels.finite_well_1d to compute eigenenergies and wavefunctions.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import dot_levels as D

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_05_finite_well.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# Parameters: V = 300 meV, width w = 5.0 nm, m_in = 0.08 m0, m_out = 0.12 m0
V_meV = 300.0
w_nm = 5.0
m_in = 0.08
m_out = 0.12
well = D.finite_well_1d(V_meV, w_nm, m_in, m_out)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Left panel: Transcendental root functions
a = w_nm / 2.0
r = m_in / m_out
Es = np.linspace(1.0, V_meV - 0.5, 1000)
k = np.sqrt(m_in * Es / D.HB2_2M0)
q = np.sqrt(m_out * (V_meV - Es) / D.HB2_2M0)

f_even = k * np.sin(k * a) - r * q * np.cos(k * a)
f_odd = k * np.cos(k * a) + r * q * np.sin(k * a)

ax1.plot(Es, f_even, color='#1769aa', lw=2.5, label=r'Even: $k\tan(ka) - \frac{m_{in}}{m_{out}}\kappa = 0$')
ax1.plot(Es, f_odd, color='#c43d3d', lw=2.5, label=r'Odd: $-k\cot(ka) - \frac{m_{in}}{m_{out}}\kappa = 0$')
ax1.axhline(0, color='gray', linestyle='--', lw=1.2)

for E_bound, p in zip(well.energies_meV, well.parity):
    color = '#1769aa' if p == 'even' else '#c43d3d'
    ax1.plot(E_bound, 0, 'o', color=color, markersize=10, zorder=5)
    ax1.annotate(f'{p.capitalize()}: {E_bound:.1f} meV',
                 xy=(E_bound, 0), xytext=(E_bound - 20, 0.4 if p == 'even' else -0.5),
                 fontsize=14, fontweight='bold', color=color,
                 arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

ax1.set_xlim(0, V_meV)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlabel('Energy E (meV)', fontsize=15)
ax1.set_ylabel('Transcendental Function Value', fontsize=15)
ax1.set_title('BenDaniel-Duke Secular Equations', fontsize=16, fontweight='bold')
ax1.grid(True, alpha=0.25)
ax1.legend(loc='upper right', fontsize=13, frameon=True)

# Right panel: Potential and wavefunctions
z = np.linspace(-6.0, 6.0, 800)
V_z = np.where(np.abs(z) <= a, 0.0, V_meV)

ax2.plot(z, V_z, color='black', lw=2.5, label='Potential V(z)')
colors = ['#1769aa', '#c43d3d']

for idx, (E_b, p) in enumerate(zip(well.energies_meV, well.parity)):
    k0 = np.sqrt(m_in * E_b / D.HB2_2M0)
    q0 = np.sqrt(m_out * (V_meV - E_b) / D.HB2_2M0)
    if p == 'even':
        psi = np.where(np.abs(z) <= a,
                       np.cos(k0 * z),
                       np.cos(k0 * a) * np.exp(-q0 * (np.abs(z) - a)))
    else:
        psi = np.where(np.abs(z) <= a,
                       np.sin(k0 * z),
                       np.sign(z) * np.sin(k0 * a) * np.exp(-q0 * (np.abs(z) - a)))
    scale = 40.0
    ax2.plot(z, E_b + scale * psi, color=colors[idx], lw=2.2,
             label=f'{p} state: E_{idx+1} = {E_b:.1f} meV')
    ax2.axhline(E_b, color=colors[idx], linestyle=':', lw=1.5)

ax2.axvline(-a, color='gray', linestyle='--', lw=1)
ax2.axvline(a, color='gray', linestyle='--', lw=1)
ax2.set_xlim(-6.0, 6.0)
ax2.set_ylim(-20, V_meV + 40)
ax2.set_xlabel('Position z (nm)', fontsize=15)
ax2.set_ylabel('Energy (meV)', fontsize=15)
ax2.set_title(f'Bound States in 1D Finite Well (w = {w_nm} nm)', fontsize=16, fontweight='bold')
ax2.grid(True, alpha=0.25)
ax2.legend(loc='upper right', fontsize=13, frameon=True)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
