"""Figure 01.12: Band alignment and offsets: InP in GaAs0.6P0.4 vs InP in GaInP.
Uses fsim_core.materials.offsets to compute strained CB and VB alignments on GaAs substrate.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import materials as M

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_12_band_offsets.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

inp = M.binary('InP')
gaas = M.binary('GaAs')
gaasp = M.GaAsP(0.4)
gainp = M.GaInP()

off_gaasp = M.offsets(inp, gaasp, gaas, 300.0)
off_gainp = M.offsets(inp, gainp, gaas, 300.0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Spatial coordinates for band profile: Barrier (-15 to -5), Dot (-5 to 5), Barrier (5 to 15)
z = np.array([-15, -5, -5, 5, 5, 15])

# Panel 1: InP dot in GaAs0.60P0.40 matrix
w1 = off_gaasp.well
b1 = off_gaasp.barrier

Ec_1 = np.array([b1.E_c, b1.E_c, w1.E_c, w1.E_c, b1.E_c, b1.E_c])
Ev_hh_1 = np.array([b1.E_hh, b1.E_hh, w1.E_hh, w1.E_hh, b1.E_hh, b1.E_hh])
Ev_lh_1 = np.array([b1.E_lh, b1.E_lh, w1.E_lh, w1.E_lh, b1.E_lh, b1.E_lh])

ax1.plot(z, Ec_1, color='#1769aa', lw=3.5, label='Conduction Band $E_c$')
ax1.plot(z, Ev_hh_1, color='#c43d3d', lw=3.5, label='Heavy Hole $E_{hh}$')
ax1.plot(z, Ev_lh_1, color='#e08e0b', lw=2.5, linestyle='--', label='Light Hole $E_{lh}$')

# Annotate conduction offset
dEc_1_meV = off_gaasp.dE_c * 1e3
ax1.annotate(rf'$\Delta E_c = {dEc_1_meV:.0f}$ meV',
             xy=(5.0, (w1.E_c + b1.E_c) / 2),
             xytext=(6.5, (w1.E_c + b1.E_c) / 2),
             fontsize=13, fontweight='bold', color='#1769aa',
             arrowprops=dict(arrowstyle='<->', color='#1769aa', lw=1.8))

# Annotate valence offset (to top VB of barrier: barrier LH is top VB because barrier is tensile!)
V_h_meV = (w1.E_hh - b1.E_lh) * 1e3
ax1.annotate(rf'$\Delta E_v = {V_h_meV:.0f}$ meV' + '\n(to tensile LH)',
             xy=(-5.0, (w1.E_hh + b1.E_lh) / 2),
             xytext=(-14.0, (w1.E_hh + b1.E_lh) / 2 - 0.05),
             fontsize=12, fontweight='bold', color='#c43d3d',
             arrowprops=dict(arrowstyle='<->', color='#c43d3d', lw=1.8))

ax1.set_title(r'InP Dot in $\mathrm{GaAs}_{0.60}\mathrm{P}_{0.40}$ (Type I)', fontsize=15, fontweight='bold')
ax1.set_xlabel('Growth position z (nm)', fontsize=14)
ax1.set_ylabel('Absolute Energy (eV, V01 scale: InSb VB = 0)', fontsize=14)
ax1.set_xlim(-15, 15)
ax1.set_ylim(-1.3, 1.1)
ax1.grid(True, alpha=0.25)
ax1.legend(loc='upper right', fontsize=11, frameon=True)
ax1.text(0, w1.E_c - 0.15, 'InP Dot\n(compressive)', ha='center', fontsize=13, fontweight='bold', color='#2e7d32')
ax1.text(-10, b1.E_c - 0.15, 'GaAsP Matrix\n(tensile)', ha='center', fontsize=12, color='gray')
ax1.text(10, b1.E_c - 0.15, 'GaAsP Matrix\n(tensile)', ha='center', fontsize=12, color='gray')

# Panel 2: InP dot in Ga0.51In0.49P matrix
w2 = off_gainp.well
b2 = off_gainp.barrier

Ec_2 = np.array([b2.E_c, b2.E_c, w2.E_c, w2.E_c, b2.E_c, b2.E_c])
Ev_hh_2 = np.array([b2.E_hh, b2.E_hh, w2.E_hh, w2.E_hh, b2.E_hh, b2.E_hh])
Ev_lh_2 = np.array([b2.E_lh, b2.E_lh, w2.E_lh, w2.E_lh, b2.E_lh, b2.E_lh])

ax2.plot(z, Ec_2, color='#1769aa', lw=3.5, label='Conduction Band $E_c$')
ax2.plot(z, Ev_hh_2, color='#c43d3d', lw=3.5, label='Heavy Hole $E_{hh}$')
ax2.plot(z, Ev_lh_2, color='#e08e0b', lw=2.5, linestyle='--', label='Light Hole $E_{lh}$')

# Annotate conduction offset
dEc_2_meV = off_gainp.dE_c * 1e3
ax2.annotate(rf'$\Delta E_c = {dEc_2_meV:.0f}$ meV',
             xy=(5.0, (w2.E_c + b2.E_c) / 2),
             xytext=(6.5, (w2.E_c + b2.E_c) / 2),
             fontsize=13, fontweight='bold', color='#1769aa',
             arrowprops=dict(arrowstyle='<->', color='#1769aa', lw=1.8))

# Annotate valence offset
dEv_2_meV = off_gainp.dE_v_hh * 1e3
ax2.annotate(rf'$\Delta E_{{v,hh}} = {dEv_2_meV:.0f}$ meV',
             xy=(-5.0, (w2.E_hh + b2.E_hh) / 2),
             xytext=(-14.0, (w2.E_hh + b2.E_hh) / 2 - 0.05),
             fontsize=12, fontweight='bold', color='#c43d3d',
             arrowprops=dict(arrowstyle='<->', color='#c43d3d', lw=1.8))

ax2.set_title(r'InP Dot in $\mathrm{Ga}_{0.51}\mathrm{In}_{0.49}\mathrm{P}$ (Type I)', fontsize=15, fontweight='bold')
ax2.set_xlabel('Growth position z (nm)', fontsize=14)
ax2.set_ylabel('Absolute Energy (eV, V01 scale: InSb VB = 0)', fontsize=14)
ax2.set_xlim(-15, 15)
ax2.set_ylim(-1.3, 1.1)
ax2.grid(True, alpha=0.25)
ax2.legend(loc='upper right', fontsize=11, frameon=True)
ax2.text(0, w2.E_c - 0.15, 'InP Dot\n(compressive)', ha='center', fontsize=13, fontweight='bold', color='#2e7d32')
ax2.text(-10, b2.E_c - 0.15, 'GaInP Matrix\n(matched)', ha='center', fontsize=12, color='gray')
ax2.text(10, b2.E_c - 0.15, 'GaInP Matrix\n(matched)', ha='center', fontsize=12, color='gray')

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
