"""Figure 01.06: Density of states across dimensionalities: 3D -> 2D -> 1D -> 0D.
Illustrates how spatial confinement transforms continuous spectra into discrete atomic-like levels.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

out_path = os.path.join(os.path.dirname(__file__), 'out', '01_06_dos.png')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(1600 / 150, 900 / 150), dpi=150)
plt.rcParams.update({'font.size': 14})

E = np.linspace(0.01, 100, 1000)

# (a) 3D Bulk
ax = axes[0, 0]
g_3d = np.sqrt(E)
ax.plot(E, g_3d, color='#1769aa', lw=3)
ax.fill_between(E, g_3d, color='#1769aa', alpha=0.15)
ax.set_title(r'3D Bulk: Continuum $g(E) \propto \sqrt{E}$', fontsize=15, fontweight='bold')
ax.set_xlabel('Energy E (meV)', fontsize=14)
ax.set_ylabel('DOS $g(E)$ (a.u.)', fontsize=14)
ax.set_xlim(0, 100)
ax.set_ylim(0, 12)
ax.grid(True, alpha=0.25)

# (b) 2D Quantum Well
ax = axes[0, 1]
subbands_2d = [15, 45, 85]
g_2d = np.zeros_like(E)
for En in subbands_2d:
    g_2d += np.where(E >= En, 3.2, 0.0)
ax.plot(E, g_2d, color='#2e7d32', lw=3)
ax.fill_between(E, g_2d, color='#2e7d32', alpha=0.15)
for En in subbands_2d:
    ax.axvline(En, color='gray', linestyle=':', alpha=0.7)
ax.set_title(r'2D Quantum Well: Staircase $g(E) \propto \sum \Theta(E-E_n)$', fontsize=15, fontweight='bold')
ax.set_xlabel('Energy E (meV)', fontsize=14)
ax.set_ylabel('DOS $g(E)$ (a.u.)', fontsize=14)
ax.set_xlim(0, 100)
ax.set_ylim(0, 12)
ax.grid(True, alpha=0.25)

# (c) 1D Quantum Wire
ax = axes[1, 0]
subbands_1d = [20, 55, 90]
g_1d = np.zeros_like(E)
for En in subbands_1d:
    dE = np.maximum(E - En, 1e-6)
    term = np.where(E >= En + 0.1, 4.0 / np.sqrt(dE), 0.0)
    term = np.clip(term, 0, 11)
    g_1d += term
ax.plot(E, g_1d, color='#e08e0b', lw=3)
ax.fill_between(E, g_1d, color='#e08e0b', alpha=0.15)
for En in subbands_1d:
    ax.axvline(En, color='gray', linestyle=':', alpha=0.7)
ax.set_title(r'1D Quantum Wire: Singularities $g(E) \propto \sum (E-E_n)^{-1/2}$', fontsize=15, fontweight='bold')
ax.set_xlabel('Energy E (meV)', fontsize=14)
ax.set_ylabel('DOS $g(E)$ (a.u.)', fontsize=14)
ax.set_xlim(0, 100)
ax.set_ylim(0, 12)
ax.grid(True, alpha=0.25)

# (d) 0D Quantum Dot ("Artificial Atom")
ax = axes[1, 1]
levels_0d = [(25, 's shell', '#c43d3d'), (50, 'p shell', '#8e24aa'), (80, 'd shell', '#00838f')]
for E0, name, col in levels_0d:
    ax.vlines(E0, 0, 9.5, colors=col, lw=4)
    ax.plot(E0, 9.5, 'o', color=col, markersize=8)
    ax.text(E0, 10.2, f'{name}\n({E0} meV)', ha='center', fontsize=12, fontweight='bold', color=col)
ax.set_title(r'0D Quantum Dot: Discrete $\delta$-Peaks ("Artificial Atom")', fontsize=15, fontweight='bold')
ax.set_xlabel('Energy E (meV)', fontsize=14)
ax.set_ylabel('DOS $g(E)$ (a.u.)', fontsize=14)
ax.set_xlim(0, 100)
ax.set_ylim(0, 12)
ax.grid(True, alpha=0.25)

plt.tight_layout()
fig.savefig(out_path, dpi=150, facecolor='white')
plt.close(fig)
print('Wrote', out_path)
