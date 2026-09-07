import os
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14

out = os.path.join(os.path.dirname(__file__), 'out', '02_22_landscape.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
src = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'out', 'rt_edge', 'sweep.csv'))
df = pd.read_csv(src)
summary = df.groupby(['gamma300_meV', 'delta_xx_meV'])['g2_pulsed'].agg(['count', 'min', 'median', 'max']).reset_index()
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150)
labels = [f"$\\Gamma$={r.gamma300_meV:.0f}, $\\Delta_{{XX}}$={r.delta_xx_meV:.0f}" for r in summary.itertuples()]
ax.boxplot([df[(df.gamma300_meV == r.gamma300_meV) & (df.delta_xx_meV == r.delta_xx_meV)].g2_pulsed for r in summary.itertuples()], tick_labels=labels, patch_artist=True, boxprops={'facecolor': '#c43d3d', 'alpha': .55}, medianprops={'color': 'black', 'linewidth': 2})
for i, r in enumerate(summary.itertuples(), 1):
    ax.text(i, 0.97, f"n={r.count}\\nmin={r.min:.3f}", ha='center', va='top', fontsize=14)
ax.set(ylim=(0, 1), ylabel='pulsed $g^{(2)}(0)$', title='Repository sweep: distributions in its four sampled cells')
ax.grid(axis='y', alpha=.25)
fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
