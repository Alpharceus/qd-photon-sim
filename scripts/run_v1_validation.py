"""V1 (work order v3): device-tier validation against the two measured
g2(T) series reaching toward/into room temperature.

Laferriere 2023 (InAsP/InP nanowire, O-band, optical): model ENVELOPE over
the 16 corners of the published-unknown [E] ranges (delta_xx, Gamma, w, mu)
per temperature, rho = 1 (authors background-correct), symmetric slit on X,
no cavity, IBM lineshape on. The 4 K point is excluded (re-excitation
mechanism, the authors' own attribution; WP-M2' tier, not spectral).

Chatzarakis 2023 ((211)B InAs/GaAs, PL): device-tier recomposition of the
V-a configuration at the fitted [DR] parameters with the published per-point
windows; Lorentzian pass = wiring check against the validated fit, IBM pass
= the new sideband content.

Pass criteria (v3): (a) ordering of g2 across T, both devices; (b) magnitudes
inside the PROPOSED x2 envelope band (band width is a Raman decision --
reported against the proposal, not stamped); (c) the 300 K Laferriere point
specifically. NO model iteration against misses (guardrail): outputs are
bounded and reported.

Artifacts: out/validation/gsq_T_validation.csv + gsq_T_validation.svg.
"""
import csv
import itertools
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.device import DeviceDesign, evaluate

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "validation"
OUT.mkdir(parents=True, exist_ok=True)


def base_design(lineshape, geom, alpha=None, rho_one=True):
    d = DeviceDesign()
    d.drive.V = 0.0
    d.drive.I_uA = 0.0
    d.drive.duty = 0.0
    d.drive.b_e = 0.0
    d.drive.mode = "PL"
    d.cavity.enabled = False
    d.filter.enabled = True
    d.filter.auto_w = False
    d.dot.lineshape = lineshape
    ph = dict(geom)
    if alpha is not None:
        ph["alpha_ps2"] = alpha
    d.dot.phonon = ph
    if rho_one:
        # rho -> 1 (Laferriere: authors background-correct their
        # correlations); tiny-but-truthy overrides dodge the 0->proxy
        # fallback in evaluate(). Chatzarakis keeps the FITTED retention/
        # background (rho(T) is part of the validated V-a configuration).
        d.ret.b0 = 1e-12
        d.ret.beta = 1e-12
    return d


def run_point(d, T, w, dx, mu, delta, gamma_target=None):
    d.thermal.T_hs = T
    d.filter.w = w
    d.filter.dx = dx
    d.drive.mu = mu
    d.dot.delta_xx = delta
    if gamma_target is not None:
        # per-point Gamma override via gamma_scale against the class proxy
        probe = DeviceDesign()
        probe.dot.gamma_scale = 1.0
        from fsim_core.device import class_proxy_params
        from fsim_core.spectral import gamma_of_T
        p = class_proxy_params()
        g_proxy = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        d.dot.gamma_scale = gamma_target / g_proxy
    s = evaluate(d, T_grid=[T])["scalars"]
    return float(s["g2_op"]), float(s["eps_op"]), float(s["t_x_op"])


rows = []

# ------------------------------------------------------------- Laferriere
card = yaml.safe_load((ROOT / "cards" / "laferriere2023.yaml").read_text())
dev = card["device"]
print("V1 -- Laferriere 2023 (envelope over 16 [E]-range corners per T)")
print(f"{'T':>5} {'meas':>6} {'model mid':>10} {'model band':>17} "
      f"{'x2 band':>13} {'in x2?':>6}")
for pt in card["data"]["g2_vs_T"]["points"]:
    T, g2m = float(pt["T"]), float(pt["g2"])
    key = str(int(T))
    gam = dev["gamma_meV"][key]
    win = dev["w_meV"][key]
    mu = dev["mu"]["mid_77"] if T < 100 else dev["mu"]["mid_high"]
    mu_rng = dev["mu"]["range_77"] if T < 100 else dev["mu"]["range_high"]
    d = base_design("ibm", dev["geometry"])
    g2_mid, eps_mid, _ = run_point(d, T, win["mid"], 0.0, mu,
                                   dev["delta_xx"]["mid"], gam["mid"])
    corners = []
    for dd, gg, ww, mm in itertools.product(dev["delta_xx"]["range"],
                                            gam["range"], win["range"],
                                            mu_rng):
        g2c, _, _ = run_point(d, T, ww, 0.0, mm, dd, gg)
        corners.append(g2c)
    lo, hi = float(np.min(corners)), float(np.max(corners))
    in_band = (g2m / 2.0) <= g2_mid <= (g2m * 2.0)
    overlap = not (hi < g2m / 2.0 or lo > g2m * 2.0)
    rows.append({"device": "laferriere", "T_K": T, "g2_measured": g2m,
                 "g2_model_mid": g2_mid, "g2_model_lo": lo, "g2_model_hi": hi,
                 "mid_in_x2_band": in_band, "envelope_overlaps_x2": overlap})
    print(f"{T:5.0f} {g2m:6.2f} {g2_mid:10.4f} [{lo:7.4f},{hi:7.4f}] "
          f"[{g2m/2:5.3f},{g2m*2:5.3f}] {'YES' if in_band else 'no':>6}")

laf = [r for r in rows if r["device"] == "laferriere"]
order_ok = all(laf[i]["g2_model_mid"] < laf[i + 1]["g2_model_mid"]
               for i in range(len(laf) - 1))
print(f"  ordering reproduced: {'YES' if order_ok else 'NO'}; "
      f"300 K point in x2 band: {'YES' if laf[-1]['mid_in_x2_band'] else 'NO'}")

# ------------------------------------------------------------ Chatzarakis
card_c = yaml.safe_load((ROOT / "cards" / "chatzarakis2023.yaml").read_text())
dev_c = card_c["device"]
print("\nV1 -- Chatzarakis 2023 (fitted [DR] params, published windows)")
print(f"{'T':>5} {'meas':>6} {'lorentzian':>10} {'ibm':>8} {'|d|_ibm':>8} "
      f"{'x2?':>4}")
for pt in card_c["data"]["g2_vs_T"]["points"]:
    T, g2m = float(pt["T"]), float(pt["g2"])
    res = {}
    for ls in ("lorentzian", "ibm"):
        d = base_design(ls, dev_c["geometry"], alpha=dev_c["phonon_alpha_ps2"],
                        rho_one=False)
        g2v, _, _ = run_point(d, T, float(pt["w"]), float(pt["dx"]),
                              dev_c["mu"], dev_c["delta_xx"])
        res[ls] = g2v
    in_band = (g2m / 2.0) <= res["ibm"] <= (g2m * 2.0) if g2m > 0 else True
    rows.append({"device": "chatzarakis", "T_K": T, "g2_measured": g2m,
                 "g2_model_mid": res["ibm"], "g2_model_lo": res["lorentzian"],
                 "g2_model_hi": res["ibm"], "mid_in_x2_band": in_band,
                 "envelope_overlaps_x2": in_band})
    print(f"{T:5.0f} {g2m:6.2f} {res['lorentzian']:10.4f} {res['ibm']:8.4f} "
          f"{abs(res['ibm'] - g2m):8.4f} {'YES' if in_band else 'no':>4}")

cha = [r for r in rows if r["device"] == "chatzarakis"]
order_c = all(cha[i]["g2_model_mid"] < cha[i + 1]["g2_model_mid"]
              for i in range(len(cha) - 1))
print(f"  ordering reproduced: {'YES' if order_c else 'NO'}")

with (OUT / "gsq_T_validation.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=False)
for ax, name, title in ((axes[0], "laferriere",
                         "Laferriere 2023 (nanowire, O-band)"),
                        (axes[1], "chatzarakis",
                         "Chatzarakis 2023 ((211)B InAs/GaAs)")):
    sub = [r for r in rows if r["device"] == name]
    Ts = [r["T_K"] for r in sub]
    ax.errorbar(Ts, [r["g2_measured"] for r in sub],
                yerr=[0.04] * len(sub), fmt="o", color="#B3403A",
                label="measured", zorder=3)
    ax.plot(Ts, [r["g2_model_mid"] for r in sub], "s-", color="#3E6DAE",
            label="model (IBM)" if name == "laferriere" else "model IBM")
    if name == "laferriere":
        ax.fill_between(Ts, [r["g2_model_lo"] for r in sub],
                        [r["g2_model_hi"] for r in sub], alpha=0.18,
                        color="#3E6DAE", label="[E]-range envelope")
    else:
        ax.plot(Ts, [r["g2_model_lo"] for r in sub], "^--", color="#557A2E",
                label="model Lorentzian")
    for r in sub:
        ax.fill_between([r["T_K"] - 3, r["T_K"] + 3],
                        [r["g2_measured"] / 2] * 2, [r["g2_measured"] * 2] * 2,
                        color="#8A7B7E", alpha=0.15)
    ax.set_xlabel("T (K)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
axes[0].set_ylabel("g$^{2}$(0)")
fig.suptitle("V1: device tier vs measured g$^2$(T)  --  grey = proposed x2 "
             "band (Raman to confirm)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "gsq_T_validation.svg")
fig.savefig(OUT / "gsq_T_validation.png", dpi=160)
print(f"\nbundle -> {OUT}")
