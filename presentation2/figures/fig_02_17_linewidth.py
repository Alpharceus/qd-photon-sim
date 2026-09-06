import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.linewidth import gamma_envelope

out = os.path.join(os.path.dirname(__file__), 'out', '02_17_linewidth.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
T = np.linspace(4, 300, 300); lo, hi = gamma_envelope(T); fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150)
ax.fill_between(T, lo, hi, color='#c43d3d', alpha=.25); ax.plot(T, lo, lw=3, label='lower class range'); ax.plot(T, hi, lw=3, label='upper class range')
ax.set(xlabel='temperature (K)', ylabel='FWHM $\\Gamma$ (meV)', title='Thermal broadening moves the linewidth into the meV range')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
