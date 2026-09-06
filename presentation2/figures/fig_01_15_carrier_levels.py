"""Figure 01.15: Carrier states and thermal escape ceiling.
Uses fsim_core.dot_levels presets to compare confinement energies and retention curves S(T).
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import dot_levels as D

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_15_carrier_levels.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

presets = D.class_presets()
lv_gaasp = D.levels(presets['InP/GaAsP0.4/AlGaAs0.4 on GaAs']())
lv_gainp = D.levels(presets['InP/GaInP/AlGaInP0.55 on GaAs'](T=300.0))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Left panel: Potential well and bound state levels for InP/GaAsP0.4
# CB well
z_well = np.array([-10, -2, -2, 2, 2, 10])
V_e_well = np.array([lv_gaasp.V_e, lv_gaasp.V_e, 0, 0, lv_gaasp.V_e, lv_gaasp.V_e])

ax1.plot(z_well, V_e_well, color='#1769aa', lw=3, label=rf'CB Well ($V_e = {lv_gaasp.V_e:.0f}$ meV)')
ax1.axhline(lv_gaasp.E_e, color='#1769aa', linestyle='--', lw=2.5,
            label=rf'Electron $E_e = {lv_gaasp.E_e:.0f}$ meV ($\Delta E_{{e,\mathrm{{esc}}}} = {lv_gaasp.dE_e_matrix:.0f}$ meV)')

# VB well (inverted for hole energy view)
V_h_well = np.array([-lv_gaasp.V_h, -lv_gaasp.V_h, 0, 0, -lv_gaasp.V_h, -lv_gaasp.V_h])
ax1.plot(z_well, V_h_well, color='#c43d3d', lw=3, label=rf'VB Well ($V_h = {lv_gaasp.V_h:.0f}$ meV)')
ax1.axhline(-lv_gaasp.E_h, color='#c43d3d', linestyle='--', lw=2.5,
            label=rf'Hole $E_h = {lv_gaasp.E_h:.0f}$ meV ($\Delta E_{{h,\mathrm{{esc}}}} = {lv_gaasp.dE_h_matrix:.0f}$ meV)')

# Annotate escape arrows
ax1.annotate(rf'$\Delta E_{{e,\mathrm{{esc}}}} = {lv_gaasp.dE_e_matrix:.0f}$ meV',
             xy=(2.5, (lv_gaasp.E_e + lv_gaasp.V_e) / 2),
             xytext=(3.5, (lv_gaasp.E_e + lv_gaasp.V_e) / 2),
             fontsize=12, fontweight='bold', color='#1769aa',
             arrowprops=dict(arrowstyle='<->', color='#1769aa', lw=1.8))

ax1.annotate(rf'$\Delta E_{{h,\mathrm{{esc}}}} = {lv_gaasp.dE_h_matrix:.0f}$ meV',
             xy=(2.5, (-lv_gaasp.E_h - lv_gaasp.V_h) / 2),
             xytext=(3.5, (-lv_gaasp.E_h - lv_gaasp.V_h) / 2),
             fontsize=12, fontweight='bold', color='#c43d3d',
             arrowprops=dict(arrowstyle='<->', color='#c43d3d', lw=1.8))

ax1.set_xlim(-10, 10)
ax1.set_ylim(-130, 240)
ax1.set_xlabel('Position z (nm)', fontsize=14)
ax1.set_ylabel('Carrier Energy relative to Band Edges (meV)', fontsize=14)
ax1.set_title(r'InP/$\mathrm{GaAs}_{0.60}\mathrm{P}_{0.40}$ Confinement & Escape', fontsize=15, fontweight='bold')
ax1.grid(True, alpha=0.25)
ax1.legend(loc='upper right', fontsize=11, frameon=True)

# Right panel: Thermal retention S(T)
T_range = np.linspace(10, 320, 300)

p_gaasp = D.retention_params(lv_gaasp, 1.0, channel="pair_half", verbose=False)
p_gainp = D.retention_params(lv_gainp, 1.0, channel="pair_half", verbose=False)

def calc_S(T, p):
    return 1.0 / (1.0 + p['a_esc'] * np.exp(-p['E_a'] / (D.KB_MEV * T)) + p['b_p'] * np.exp(-p['E_b'] / (D.KB_MEV * T)))

S_gaasp = calc_S(T_range, p_gaasp)
S_gainp = calc_S(T_range, p_gainp)

ax2.semilogy(T_range, S_gainp, color='#2e7d32', lw=3,
             label=rf'InP/GaInP: $E_a = {p_gainp["E_a"]:.0f}$ meV, $S(300\mathrm{{K}}) \approx {calc_S(300, p_gainp):.3f}$')
ax2.semilogy(T_range, S_gaasp, color='#c43d3d', lw=3,
             label=rf'InP/GaAsP0.4: $E_a = {p_gaasp["E_a"]:.0f}$ meV, $S(300\mathrm{{K}}) \approx {calc_S(300, p_gaasp):.4f}$')

ax2.axvline(300, color='gray', linestyle=':', lw=1.5, label='Room Temp 300 K')
ax2.axhline(0.01, color='black', linestyle='--', alpha=0.5, label='1% Retention threshold')

ax2.set_xlim(0, 320)
ax2.set_ylim(1e-5, 1.5)
ax2.set_xlabel('Temperature T (K)', fontsize=14)
ax2.set_ylabel('Thermal Retention Factor S(T)', fontsize=14)
ax2.set_title('Arrhenius Retention S(T): The 300 K Ceiling', fontsize=15, fontweight='bold')
ax2.grid(True, which='both', alpha=0.25)
ax2.legend(loc='lower left', fontsize=11, frameon=True)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
