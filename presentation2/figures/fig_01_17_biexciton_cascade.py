"""Figure 01.17: Biexciton binding energy, radiative cascade, and X-XX spectral splitting.
Illustrates the three-level ladder and the distinct emission lines used for single-photon filtering.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_17_biexciton_cascade.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

# Left panel: Three-level ladder
# Energies
E_0 = 0.0
E_X = 1.650      # eV
delta_xx = 0.005  # 5 meV
E_XX = 2 * E_X - delta_xx

# Plot levels
x_left = 1.0
x_right = 3.5

# |0> Ground State
ax1.plot([x_left, x_right], [E_0, E_0], color='black', lw=4)
ax1.text(x_right + 0.2, E_0, r'$|0\rangle$ (Empty Dot)', va='center', fontsize=14, fontweight='bold')

# |X> Exciton
ax1.plot([x_left, x_right], [E_X, E_X], color='#1769aa', lw=4)
ax1.text(x_right + 0.2, E_X, r'$|X\rangle$ (Exciton, $E_X$)', va='center', fontsize=14, fontweight='bold', color='#1769aa')

# |XX> Biexciton
ax1.plot([x_left, x_right], [E_XX, E_XX], color='#c43d3d', lw=4)
ax1.text(x_right + 0.2, E_XX, r'$|XX\rangle$ ($2E_X - \Delta_{XX}$)', va='center', fontsize=14, fontweight='bold', color='#c43d3d')

# Cascade arrows
# XX -> X
ax1.annotate('', xy=(2.0, E_X), xytext=(2.0, E_XX),
             arrowprops=dict(arrowstyle='->', color='#c43d3d', lw=3, mutation_scale=20))
ax1.text(1.3, (E_X + E_XX) / 2, r'$\gamma_{XX}$: $h\nu_{XX} = E_X - \Delta_{XX}$',
         fontsize=13, fontweight='bold', color='#c43d3d', va='center')

# X -> 0
ax1.annotate('', xy=(2.5, E_0), xytext=(2.5, E_X),
             arrowprops=dict(arrowstyle='->', color='#1769aa', lw=3, mutation_scale=20))
ax1.text(2.6, (E_0 + E_X) / 2, r'$\gamma_X$: $h\nu_X = E_X$',
         fontsize=13, fontweight='bold', color='#1769aa', va='center')

ax1.set_xlim(0.5, 5.5)
ax1.set_ylim(-0.3, 3.6)
ax1.set_ylabel('Total State Energy (eV)', fontsize=15)
ax1.set_title(r'Biexciton Cascade Ladder: $|XX\rangle \rightarrow |X\rangle \rightarrow |0\rangle$',
              fontsize=16, fontweight='bold')
ax1.set_xticks([])
ax1.grid(True, axis='y', alpha=0.25)

# Right panel: Resulting emission spectrum
detuning = np.linspace(-15, 10, 1000)  # meV relative to E_X
gamma_x = 1.0   # meV FWHM
gamma_xx = 1.2  # meV FWHM
delta_xx_meV = 5.0

# Lorentzian lineshapes
I_X = (gamma_x / 2)**2 / (detuning**2 + (gamma_x / 2)**2)
I_XX = 0.85 * (gamma_xx / 2)**2 / ((detuning + delta_xx_meV)**2 + (gamma_xx / 2)**2)

ax2.plot(detuning, I_X, color='#1769aa', lw=3, label=r'Exciton $X$ ($h\nu = E_X$)')
ax2.fill_between(detuning, I_X, color='#1769aa', alpha=0.15)

ax2.plot(detuning, I_XX, color='#c43d3d', lw=3, label=rf'Biexciton $XX$ ($h\nu = E_X - \Delta_{{XX}}$)')
ax2.fill_between(detuning, I_XX, color='#c43d3d', alpha=0.15)

# Bandpass filter transmission curve isolating X
w_filter = 3.0  # meV FWHM
T_filter = np.exp(-4 * np.log(2) * (detuning / w_filter)**2)
ax2.plot(detuning, T_filter, color='#2e7d32', lw=2.5, linestyle='--',
         label=r'Spectral Filter $T(\omega)$ (rejects $XX$)')

# Arrow for Delta_XX
ax2.annotate(rf'$\Delta_{{XX}} = {delta_xx_meV:.1f}$ meV',
             xy=(-delta_xx_meV, 0.85),
             xytext=(0.0, 0.85),
             fontsize=13, fontweight='bold', color='#8e24aa',
             arrowprops=dict(arrowstyle='<->', color='#8e24aa', lw=2))

ax2.set_xlim(-15, 10)
ax2.set_ylim(0, 1.15)
ax2.set_xlabel(r'Detuning from Exciton $E - E_X$ (meV)', fontsize=15)
ax2.set_ylabel('Emission Intensity (a.u.)', fontsize=15)
ax2.set_title(r'Spectral Splitting $\Delta_{XX}$ Enables Single-Photon Purity',
              fontsize=16, fontweight='bold')
ax2.grid(True, alpha=0.25)
ax2.legend(loc='upper left', fontsize=12, frameon=True)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
