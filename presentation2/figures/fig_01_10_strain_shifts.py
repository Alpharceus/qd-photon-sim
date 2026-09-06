"""Figure 01.10: Heterostructure strain shifts: hydrostatic shift and HH/LH splitting.
Uses fsim_core.materials.strain_shifts for InP on GaAs and InP on GaAs0.6P0.4.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import materials as M

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_10_strain_shifts.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

inp = M.binary('InP')
gaas = M.binary('GaAs')
gaasp = M.GaAsP(0.4)

# Mismatch values
f_gaas = M.mismatch(inp, gaas, 300.0)    # -0.03688 (-3.69%)
f_gaasp = M.mismatch(inp, gaasp, 300.0)  # -0.05069 (-5.07%)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Left panel: Band edge shifts vs strain
eps_range = np.linspace(-0.06, 0.02, 300)
dEc = []
dEhh = []
dElh = []

for ep in eps_range:
    s = M.strain_shifts(inp, gaas, eps_par=ep)
    dEc.append(s.dE_c * 1e3)    # meV
    dEhh.append(s.dE_hh * 1e3)  # meV
    dElh.append(s.dE_lh * 1e3)  # meV

ax1.plot(eps_range * 100, dEc, color='#1769aa', lw=3, label=r'CB shift $\Delta E_c = a_c \mathrm{tr}(\varepsilon)$')
ax1.plot(eps_range * 100, dEhh, color='#c43d3d', lw=3, label=r'HH shift $\Delta E_{hh} = a_v \mathrm{tr}(\varepsilon) - Q$')
ax1.plot(eps_range * 100, dElh, color='#e08e0b', lw=3, linestyle='--', label=r'LH shift $\Delta E_{lh}$')
ax1.axvline(0, color='gray', linestyle=':', lw=1.2)
ax1.axhline(0, color='gray', linestyle=':', lw=1.2)

# Mark InP on GaAs and InP on GaAs0.6P0.4
s_gaas = M.strain_shifts(inp, gaas)
s_gaasp = M.strain_shifts(inp, gaasp)

ax1.plot(f_gaas * 100, s_gaas.dE_c * 1e3, 'o', color='#1769aa', markersize=9)
ax1.plot(f_gaas * 100, s_gaas.dE_hh * 1e3, 's', color='#c43d3d', markersize=9)
ax1.axvline(f_gaas * 100, color='#2e7d32', linestyle='--', lw=1.5,
            label=f'InP on GaAs (f = {f_gaas*100:.2f}%)')

ax1.plot(f_gaasp * 100, s_gaasp.dE_c * 1e3, 'o', color='#1769aa', markersize=9)
ax1.plot(f_gaasp * 100, s_gaasp.dE_hh * 1e3, 's', color='#c43d3d', markersize=9)
ax1.axvline(f_gaasp * 100, color='#8e24aa', linestyle='--', lw=1.5,
            label=f'InP on GaAs0.6P0.4 (f = {f_gaasp*100:.2f}%)')

ax1.set_xlabel(r'In-plane misfit strain $\varepsilon_\parallel$ (%)', fontsize=15)
ax1.set_ylabel('Energy shift relative to unstrained (meV)', fontsize=15)
ax1.set_title('Strain Shifts from Pikus-Bir & Van de Walle', fontsize=16, fontweight='bold')
ax1.set_xlim(-6.0, 2.0)
ax1.grid(True, alpha=0.25)
ax1.legend(loc='lower right', fontsize=12, frameon=True)

# Right panel: Band structure schematic before and after compressive strain
x_unstrained = [0.5, 1.5]
x_compressive = [3.0, 4.0]

# Unstrained levels
E_c0 = 1.353  # eV
E_v0 = 0.000  # eV

# Compressive levels (InP on GaAs0.6P0.4)
E_c_str = E_c0 + s_gaasp.dE_c
E_hh_str = E_v0 + s_gaasp.dE_hh
E_lh_str = E_v0 + s_gaasp.dE_lh

# Plot Unstrained
ax2.plot(x_unstrained, [E_c0, E_c0], color='#1769aa', lw=4)
ax2.text(1.0, E_c0 + 0.03, r'CB ($E_{c0}$)', ha='center', fontsize=13, fontweight='bold', color='#1769aa')
ax2.plot(x_unstrained, [E_v0, E_v0], color='#c43d3d', lw=4)
ax2.text(1.0, E_v0 - 0.05, r'VB degenerate HH/LH', ha='center', fontsize=13, fontweight='bold', color='#c43d3d')

# Plot Compressive Strained
ax2.plot(x_compressive, [E_c_str, E_c_str], color='#1769aa', lw=4)
ax2.text(3.5, E_c_str + 0.03, rf'Strained CB (+{s_gaasp.dE_c*1e3:.0f} meV)', ha='center', fontsize=13, fontweight='bold', color='#1769aa')

ax2.plot(x_compressive, [E_hh_str, E_hh_str], color='#c43d3d', lw=4)
ax2.text(3.5, E_hh_str + 0.03, rf'HH (+{s_gaasp.dE_hh*1e3:.0f} meV)', ha='center', fontsize=13, fontweight='bold', color='#c43d3d')

ax2.plot(x_compressive, [E_lh_str, E_lh_str], color='#e08e0b', lw=4, linestyle='--')
ax2.text(3.5, E_lh_str - 0.06, rf'LH (+{s_gaasp.dE_lh*1e3:.0f} meV)', ha='center', fontsize=13, fontweight='bold', color='#e08e0b')

# Connect with dashed lines
ax2.plot([1.5, 3.0], [E_c0, E_c_str], 'k:', lw=1.5)
ax2.plot([1.5, 3.0], [E_v0, E_hh_str], 'k:', lw=1.5)
ax2.plot([1.5, 3.0], [E_v0, E_lh_str], 'k:', lw=1.5)

# Shear splitting bracket
ax2.annotate(rf'Shear splitting $\Delta E_{{HH-LH}} = {(s_gaasp.dE_hh - s_gaasp.dE_lh)*1e3:.0f}$ meV',
             xy=(4.1, (E_hh_str + E_lh_str) / 2),
             xytext=(4.3, (E_hh_str + E_lh_str) / 2),
             fontsize=12, fontweight='bold', color='#8e24aa',
             arrowprops=dict(arrowstyle='<->', color='#8e24aa', lw=1.8))

ax2.set_xlim(0, 5.8)
ax2.set_ylim(-0.2, 1.8)
ax2.set_xticks([1.0, 3.5])
ax2.set_xticklabels(['Unstrained InP', 'Compressive InP\n(on GaAs0.60P0.40)'], fontsize=14, fontweight='bold')
ax2.set_ylabel('Band Edge Energy (eV)', fontsize=15)
ax2.set_title('Compressive Strain: HH Tops Valence Band', fontsize=16, fontweight='bold')
ax2.grid(True, alpha=0.25)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
