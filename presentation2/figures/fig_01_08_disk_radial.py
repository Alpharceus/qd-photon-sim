"""Figure 01.08: Separable 0D disk model and Bessel radial solutions.
Uses fsim_core.dot_levels.finite_disk_2d to compute radial wavefunctions and s/p shell levels.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, kve

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import dot_levels as D

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_08_disk_radial.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# Parameters for in-plane disk: V = 120 meV, R = 12.0 nm, m_in = 0.08, m_out = 0.08
V_meV = 120.0
R_nm = 12.0
m_in = 0.08
m_out = 0.08
disk = D.finite_disk_2d(V_meV, R_nm, m_in, m_out)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Radial wavefunctions
r_in = np.linspace(0.01, R_nm, 500)
r_out = np.linspace(R_nm, 25.0, 500)
r = np.concatenate([r_in, r_out])

# l = 0 (s shell)
k0 = np.sqrt(m_in * disk.E0_meV / D.HB2_2M0)
q0 = np.sqrt(m_out * (V_meV - disk.E0_meV) / D.HB2_2M0)
psi0_in = jv(0, k0 * r_in)
psi0_out = jv(0, k0 * R_nm) * kve(0, q0 * r_out) / kve(0, q0 * R_nm) * np.exp(-q0 * (r_out - R_nm))
psi0 = np.concatenate([psi0_in, psi0_out])
norm0 = np.sqrt(np.trapezoid(r * psi0**2, r))
psi0 /= norm0

# l = 1 (p shell)
k1 = np.sqrt(m_in * disk.E1_meV / D.HB2_2M0)
q1 = np.sqrt(m_out * (V_meV - disk.E1_meV) / D.HB2_2M0)
psi1_in = jv(1, k1 * r_in)
psi1_out = jv(1, k1 * R_nm) * kve(1, q1 * r_out) / kve(1, q1 * R_nm) * np.exp(-q1 * (r_out - R_nm))
psi1 = np.concatenate([psi1_in, psi1_out])
norm1 = np.sqrt(np.trapezoid(r * psi1**2, r))
psi1 /= norm1

# Left panel: Wavefunctions
ax1.plot(r, psi0, color='#1769aa', lw=3, label=rf's-shell ($l=0$): $J_0(kr) \leftrightarrow K_0(\kappa r)$')
ax1.plot(r, psi1, color='#c43d3d', lw=3, label=rf'p-shell ($l=1$): $J_1(kr) \leftrightarrow K_1(\kappa r)$')
ax1.axvline(R_nm, color='gray', linestyle='--', lw=1.5, label=f'Disk radius R = {R_nm} nm')
ax1.set_xlabel('Radial distance r (nm)', fontsize=15)
ax1.set_ylabel('Radial Wavefunction $R_{nl}(r)$ (a.u.)', fontsize=15)
ax1.set_title('In-Plane Bessel Eigenfunctions', fontsize=16, fontweight='bold')
ax1.set_xlim(0, 25)
ax1.grid(True, alpha=0.25)
ax1.legend(loc='upper right', fontsize=13, frameon=True)

# Right panel: Potential well and discrete levels
r_well = np.linspace(-25, 25, 1000)
V_r = np.where(np.abs(r_well) <= R_nm, 0.0, V_meV)
ax2.plot(r_well, V_r, color='black', lw=2.5, label='Radial well V(r)')

ax2.axhline(disk.E0_meV, color='#1769aa', lw=3, label=f's-shell: $E_0 = {disk.E0_meV:.1f}$ meV (degen=2)')
ax2.axhline(disk.E1_meV, color='#c43d3d', lw=3, label=f'p-shell: $E_1 = {disk.E1_meV:.1f}$ meV (degen=4)')

sp_split = disk.E1_meV - disk.E0_meV
ax2.annotate(rf'$\Delta E_{{sp}} = {sp_split:.1f}$ meV',
             xy=(0, (disk.E0_meV + disk.E1_meV) / 2),
             xytext=(3, (disk.E0_meV + disk.E1_meV) / 2),
             fontsize=14, fontweight='bold', color='#2e7d32',
             arrowprops=dict(arrowstyle='<->', color='#2e7d32', lw=2))

ax2.axvline(-R_nm, color='gray', linestyle='--', lw=1)
ax2.axvline(R_nm, color='gray', linestyle='--', lw=1)
ax2.set_xlim(-25, 25)
ax2.set_ylim(-10, V_meV + 25)
ax2.set_xlabel('In-plane coordinate r (nm)', fontsize=15)
ax2.set_ylabel('Energy (meV)', fontsize=15)
ax2.set_title(rf'Disk Shell Structure ($\Delta E_{{sp}} = {sp_split:.1f}$ meV)', fontsize=16, fontweight='bold')
ax2.grid(True, alpha=0.25)
ax2.legend(loc='upper right', fontsize=13, frameon=True)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
