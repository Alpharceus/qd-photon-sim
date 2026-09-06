import os, csv
import numpy as np
import matplotlib.pyplot as plt

out = os.path.join(os.path.dirname(__file__), 'out', '02_22_landscape.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
src = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'out', 'rt_edge', 'sweep.csv'))
with open(src, encoding='utf-8') as h: rows = list(csv.DictReader(h))
g = sorted({float(x['gamma300_meV']) for x in rows}); d = sorted({float(x['delta_xx_meV']) for x in rows}); z = np.full((len(d), len(g)), np.nan)
for x in rows: z[d.index(float(x['delta_xx_meV'])), g.index(float(x['gamma300_meV']))] = float(x['g2_pulsed'])
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150); im = ax.imshow(z, origin='lower', aspect='auto', extent=[min(g), max(g), min(d), max(d)], vmin=0, vmax=1, cmap='magma_r')
fig.colorbar(im, ax=ax, label='pulsed $g^{(2)}(0)$'); ax.set(xlabel='$\\Gamma$ (meV)', ylabel='$\\Delta_{XX}$ (meV)', title='Repository sweep: linewidth versus biexciton separation')
fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
