"""Plot the favourable-corner temperature slices from the live acceptance verdict."""
from pathlib import Path
import matplotlib.pyplot as plt

out = Path(__file__).parent / "out" / "05_15_temperature.png"
out.parent.mkdir(exist_ok=True)
T=[230,250,273,300]; g=[.3214,.376,.3775,.3986]
fig, ax=plt.subplots(figsize=(1600/150,900/150),dpi=150,facecolor="white")
ax.plot(T,g,"o-",lw=3,color="#8b2c2c",label="minimum pulsed g²(0)")
ax.axhline(.5,color="#555",ls="--",label="headline gate")
ax.set(xlabel="heat-sink temperature (K)",ylabel="pulsed g²(0)",ylim=(0,.62))
ax.grid(alpha=.25); ax.legend(fontsize=14); ax.tick_params(labelsize=14)
fig.tight_layout(); fig.savefig(out,dpi=150,facecolor="white")
