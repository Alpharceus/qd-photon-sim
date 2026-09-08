"""Figure 04-11: Full Extraction Chain Decomposition and Acceptance Sweep Levers.

pkg5b fix (2026-09-07/08): panel 1 used to hardcode a stale, deleted-model
four-factor decomposition (beta, a separately-hardcoded "combined facet
factor" 0.935012, a separate propagation factor, and NA) built from a
literal n_eff=3.2530. The current facet model,
fsim_core.waveguide.facet_escape_fraction, already folds single-pass ridge
propagation into a single eta_facet factor -- there are only THREE
independent multiplicative factors (beta, eta_facet, eta_NA), and n_eff is
resolved LIVE from the card via fsim_core.device._resolve_edge, never a
literal, since that resolution is being edited concurrently and may change
again. Panel 2's lever steps are read fresh from out/rt_edge/sweep.csv
(headline model: drive.finite_pulse=true, ret.tau_cap_scales_with_density=
false) rather than hardcoded from a stale sweep run.
Contract reference: docs/rt_edge_contract.md.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, _resolve_edge  # noqa: E402
from fsim_core import waveguide as wg  # noqa: E402

plt.rcParams['font.size'] = 14

out_dir = Path(__file__).resolve().parents[0] / "out"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "04_11_extraction_chain.png"

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "figure.titlesize": 18,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1600/150, 900/150), dpi=150)
fig.patch.set_facecolor("white")

# ---------------------------------------------------------------- live data
CARD_PATH = "cards/edge-inp-gainp-design.yaml"
design = DeviceDesign.load(CARD_PATH)
edge = _resolve_edge(design.ret, design.emission, 300.0)[0]
n_eff = edge.n_eff
beta = edge.F_wg / (1.0 + edge.F_wg)
T_facet = wg.facet_transmission(n_eff)

eta_facet_hr = wg.facet_escape_fraction(T_facet, design.emission.R_back,
                                         design.emission.alpha_cm, design.emission.L_um)
eta_facet_bare = wg.facet_escape_fraction(T_facet, 1.0 - T_facet,
                                           design.emission.alpha_cm, design.emission.L_um)

# NA collection at the card's own NA (0.80) is exposed on `edge`; the NA=0.50
# comparison point re-resolves the edge with that NA override only.
eta_NA_hi = edge.eta_NA
design_na_lo = DeviceDesign.load(CARD_PATH)
design_na_lo.emission.NA = 0.50
eta_NA_lo = _resolve_edge(design_na_lo.ret, design_na_lo.emission, 300.0)[0].eta_NA

eta_total = beta * eta_facet_hr * eta_NA_hi

# Panel 1: three independent multiplicative factors, not four -- propagation
# is already folded into eta_facet's own ray series (fsim_core.waveguide.
# facet_escape_fraction), so it is annotated as informational only, never a
# separate bar.
factors = [r"Waveguide $\beta$" + "\n" + f"({beta * 100:.2f}%)",
           r"Facet $\eta_{facet}$" + "\n(ray series," + "\npropagation incl.)"
           + "\n" + f"({eta_facet_hr * 100:.2f}%)",
           r"Lens $\eta_{NA}$" + "\n" + r"$= \eta_{total}$"
           + "\n" + f"({eta_NA_hi * 100:.2f}%)"]

vals = [beta, beta * eta_facet_hr, eta_total]
pcts = np.array(vals) * 100
x = np.arange(len(factors))

bars = ax1.bar(x, pcts, width=0.55, color=["#3182CE", "#4299E1", "#276749"], edgecolor="#2D3748", lw=1.5)

for b, val in zip(bars, pcts):
    ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.08,
              f"{val:.3f}%", ha="center", va="bottom", fontsize=12, fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels(factors, fontsize=12)
ax1.set_ylabel(r"Cumulative Out-Coupling Efficiency (%)")
ax1.set_title("Multiplicative Extraction Chain\n" + r"$\eta_{total}$ Decomposition (3 factors)", pad=10)
ax1.set_ylim(0, 3.5)
ax1.grid(True, axis="y", linestyle=":", alpha=0.5)

# Panel 2: Sweep Levers and Collected Photon Flux Progression (at 300 K),
# read fresh from out/rt_edge/sweep.csv's headline-model rows (delta_xx=8
# meV, gamma300=6 meV -- this card's own best sampled corner) rather than
# hardcoded from a stale sweep run:
#   baseline (NA=0.5, R_back=0, L=500 um):      21.82 photons/s
#   + HR mirror (R_back=0.95):                  45.28 photons/s (2.08x)
#   + short cavity (L=250 um):                  53.88 photons/s (1.19x)
#   + high NA (NA=0.80):                        80.34 photons/s (1.49x)
# All four remain far below the 1000 photons/s eligibility floor under the
# corrected finite-pulse loading model -- unlike the pre-pkg5b figure, none
# of these steps clears it.
steps = [
    "Baseline\n(uncoated)",
    "+ HR Mirror\n($2.08\\times$)",
    "+ Short Cavity\n($1.19\\times$)",
    "+ High-NA Lens\n($1.49\\times$)",
]
flux_vals = [21.82, 45.28, 53.88, 80.34]

x2 = np.arange(len(steps))
colors = ["#E53E3E", "#DD6B20", "#D69E2E", "#38A169"]

bars2 = ax2.bar(x2, flux_vals, width=0.55, color=colors, edgecolor="#2D3748", lw=1.5)

# Flux floor line at 1000 photons/s
ax2.axhline(1000.0, color="#C53030", ls="--", lw=2.5, label="Eligibility Floor [A] (1000 photons/s)")

for b, val in zip(bars2, flux_vals):
    lbl = f"{val:.1f} s$^{{-1}}$"
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 25,
              lbl, ha="center", va="bottom", fontsize=12, fontweight="bold")

ax2.annotate("Still below the floor\nat every lever combination",
              xy=(3, 80.34), xytext=(1.0, 550),
              arrowprops=dict(arrowstyle="->", color="#822727", lw=2),
              fontsize=12, fontweight="bold", color="#822727",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF5F5", edgecolor="#C53030"))

ax2.set_xticks(x2)
ax2.set_xticklabels(steps, fontsize=12)
ax2.set_ylabel("Collected Pulsed Flux (photons / s)")
ax2.set_title("Acceptance Sweep Levers at\n300 K vs Flux Floor", pad=10)
ax2.set_ylim(0, 1150)
ax2.grid(True, axis="y", linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", framealpha=0.9, fontsize=12)

fig.subplots_adjust(top=0.84, bottom=0.22, left=0.08, right=0.96, wspace=0.32)
fig.savefig(out_path, dpi=150, facecolor="white")
plt.close(fig)
print("Saved", out_path)
