"""Phase 3 driver: prediction runs (planning doc section 4).

(i)   staged InP/GaAsP envelope at 77-120 K (the anchored deliverable)
(ii)  piezo-variant (Delta, rho) requirement map vs the 300 K target
(iii) aperture/density design rules from the F5 lemma
(iv)  measurement-priority ranking for the in-house queue

Every output is a requirement envelope over [A]-swept inputs (ranges swept,
never averaged); the tag chain of every panel is [A].
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card
from fsim_core.integrator import retention
from fsim_core.loading import aperture_g2, f1b_g2, n_window_competitors
from fsim_core.spectral import epsilon, gamma_of_T
from fsim_core.thermal import Layer, Stack, t_junction
from fsim_viz.figures import phase3_packet_figure

FITP = json.loads((ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]


def gamma_class(T, scale=1.0, anchor300=None):
    """Arsenide-class Gamma(T) from the V-a joint fit [A as proxy]; optionally
    rescaled so Gamma(300) hits a [V] anchor (removes the fit's known
    high-T overshoot -- the rescale itself is [A])."""
    g = gamma_of_T(T, FITP["gamma0"], FITP["a_ac"], FITP["b_lo"], FITP["E_lo"])
    if anchor300 is not None:
        g = g * anchor300 / gamma_of_T(300.0, FITP["gamma0"], FITP["a_ac"],
                                       FITP["b_lo"], FITP["E_lo"])
    return scale * g


def rho_class(T, b_extra=0.0):
    """Class retention/background rho from the V-a fit, plus extra (injection)
    background in signal units."""
    S = float(retention(T, FITP["a_esc"], FITP["E_a"], FITP["b_p"], FITP["E_b"]))
    B = FITP["b0"] + FITP["beta"] * (1.0 - S) + b_extra
    return S / (S + B)


def eps_min(T, delta, anchor300):
    g = gamma_class(T, anchor300=anchor300)
    return 1.0 / (1.0 + (2.0 * delta / g) ** 2)


# ---------------------------------------------------------- (ii) (Delta, rho) map

def tc_of_delta_rho(delta, rho, anchor300, T_hi=420.0):
    f = lambda T: rho**2 * (1.0 - eps_min(T, delta, anchor300)) - 0.5
    if f(4.0) < 0:
        return np.nan
    if f(T_hi) > 0:
        return T_hi
    lo, hi = 4.0, T_hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def run_map():
    card = load_card(ROOT / "cards" / "qcap-piezo-variant.yaml")
    g_lo, g_hi = card["gamma_300_lo"].fixed, card["gamma_300_hi"].fixed
    g_mid = 0.5 * (g_lo + g_hi)
    deltas = np.linspace(1.0, 13.0, 61)
    rhos = np.linspace(0.60, 1.00, 51)
    Tc = np.array([[tc_of_delta_rho(d, r, g_mid) for d in deltas] for r in rhos])

    # required rho for Tc >= 300 K, band over the [V] Gamma(300) interval
    rho_req = {}
    for g in (g_lo, g_hi):
        e300 = np.array([eps_min(300.0, d, g) for d in deltas])
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.sqrt(0.5 / (1.0 - e300))
        rho_req[g] = np.where(rr <= 1.0, rr, np.nan)

    pts = []
    for p in card.datasets["device_points"].rows:
        rho = p["rho"] if p["rho"] > 0 else rho_class(p["T"])
        pts.append({"name": p["name"], "delta": p["delta"], "rho": float(rho)})

    print("(ii) piezo-variant (Delta, rho) requirement map [A chain, Gamma(300) "
          f"anchored {g_lo}-{g_hi} meV [V]]:")
    for d in (4.0, 5.0, 5.4, 8.0):
        i = np.argmin(np.abs(deltas - d))
        rr_lo, rr_hi = rho_req[g_lo][i], rho_req[g_hi][i]
        print(f"    Delta = {d:4.1f} meV: rho >= {rr_hi:.3f}-{rr_lo:.3f} for 300 K "
              f"(narrow-filter bound; brightness cost extra)")
    for p in pts:
        print(f"    device point: {p['name']:20s} (Delta={p['delta']}, rho={p['rho']:.3f})")
    return {"deltas": deltas, "rhos": rhos, "Tc": Tc, "rho_req": rho_req,
            "points": pts, "g_band": (g_lo, g_hi)}


# ------------------------------------------------------- (i) staged 77-120 K envelope

def run_staged():
    card = load_card(ROOT / "cards" / "qcap-staged.yaml")
    d_lo, d_hi = card["delta_xx_staged"].bounds
    deltas = np.linspace(d_lo, d_hi, 15)
    scales = np.array(card["gamma_scale"].bounds + (1.0,))
    b_es = np.array([card["b_e_total"].bounds[0], 0.03, card["b_e_total"].bounds[1]])
    mus = np.array([card["mu_op_staged"].bounds[0], 0.5, card["mu_op_staged"].bounds[1]])

    # junction heating at the design point: 1 um mesa, 10 uA CW (worst epi-k)
    t_epi = card["epi_thickness"].fixed * 1e-6
    a_epi = card["epi_alpha"].fixed
    k_lo = card["epi_k"].bounds[0]
    st = Stack(layers=[Layer(t_epi, k_lo, a_epi)], k_sub300=55.0)
    P = 10e-6 * card["V_drive"].fixed

    out = {}
    for T_hs in (77.0, 120.0):
        Tj = t_junction(P, 0.5e-6, st, T_hs)
        band = []
        for d in deltas:
            g2s = []
            for sc in scales:
                gam = float(gamma_class(Tj, scale=sc))
                e = epsilon(d, gam, FITP.get("r_xx", 1.0) * gam, w=gam).eps  # w=Gamma convention [A]
                for b_e in b_es:
                    rho = rho_class(Tj, b_extra=b_e)
                    for mu in mus:
                        g2s.append(1.0 - rho**2 * (1.0 - float(f1b_g2(mu, e))))
            band.append((min(g2s), max(g2s)))
        out[T_hs] = {"Tj": Tj, "band": np.array(band)}

    print("\n(i) staged InP/GaAsP envelope (w = Gamma convention, class retention "
          "[A]; DeltaT_J at 1 um/10 uA worst epi-k):")
    for T_hs, r in out.items():
        print(f"    T_hs = {T_hs:.0f} K -> T_j = {r['Tj']:.1f} K")
        for target in (0.1, 0.5):
            ok = [d for d, (lo, hi) in zip(deltas, r["band"]) if hi <= target]
            req = f"Delta >= {min(ok):.1f} meV" if ok else "NOT REACHABLE in swept band"
            print(f"      g2 <= {target} guaranteed (worst-case over [A] sweep): {req}")
    return {"deltas": deltas, "envelopes": out}


# ------------------------------------------------------------ (iii) F5 aperture rules

def run_aperture():
    card = load_card(ROOT / "cards" / "qcap-staged.yaml")
    sigma = card["sigma_inh"].fixed
    r = card["comp_brightness"].fixed
    w = float(gamma_class(120.0))  # window ~ Gamma(T_target) [A]
    diams = np.linspace(0.2, 3.0, 40)
    dens = np.logspace(8, np.log10(2e10), 40)
    Nw = np.array([[n_window_competitors(n, np.pi * (d / 2) ** 2, w, sigma)
                    for d in diams] for n in dens])
    g2pen = 1.0 - (1.0 + Nw * r**2) / (1.0 + Nw * r) ** 2

    n0, A0 = 1e10, 0.1
    nw0 = n_window_competitors(n0, A0, 6.5, sigma)
    print(f"\n(iii) F5 aperture/density rules (w = Gamma(120 K) = {w:.1f} meV, "
          f"sigma_inh = {sigma:.0f} meV, competitor brightness r = {r}):")
    print(f"    foundation checkpoint: 1e10 cm^-2, 0.1 um^2, w=6.5 -> N_w = {nw0:.1f} "
          f"(geometric {n0 * A0 * 1e-8:.0f}) -- marginal at high T [reproduced]")
    for tgt in (0.01, 0.05):
        # largest density at 1 um aperture keeping worst-case penalty <= tgt
        i_d = np.argmin(np.abs(diams - 1.0))
        ok = [n for j, n in enumerate(dens) if g2pen[j, i_d] <= tgt]
        lim = f"density <= {max(ok):.1e} cm^-2" if ok else "no density in swept range"
        print(f"    1 um aperture, g2 penalty <= {tgt}: {lim}")
    print("    (note: one 30%-bright in-window competitor alone costs g2 = "
          f"{aperture_g2([1.0, 0.3]):.2f} -- the F5 lemma is unforgiving)")
    return {"diams": diams, "dens": dens, "Nw": Nw, "g2pen": g2pen, "w": w}


# ------------------------------------------- (iv) measurement-priority ranking

def run_ranking():
    card = load_card(ROOT / "cards" / "qcap-staged.yaml")
    grids = {
        "delta_xx": np.linspace(*card["delta_xx_staged"].bounds, 6),
        "gamma_scale": np.linspace(*card["gamma_scale"].bounds, 5),
        "b_e": np.geomspace(*card["b_e_total"].bounds, 5),
        "mu": np.linspace(*card["mu_op_staged"].bounds, 4),
    }

    def tc(delta, sc, b_e, mu, T_hi=400.0):
        def g2(T):
            gam = float(gamma_class(T, scale=sc))
            e = epsilon(delta, gam, FITP.get("r_xx", 1.0) * gam, w=gam).eps
            rho = rho_class(T, b_extra=b_e)
            return 1.0 - rho**2 * (1.0 - float(f1b_g2(mu, e)))
        if g2(4.0) >= 0.5:
            return 4.0
        if g2(T_hi) < 0.5:
            return T_hi
        lo, hi = 4.0, T_hi
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            lo, hi = (lo, mid) if g2(mid) >= 0.5 else (mid, hi)
        return 0.5 * (lo + hi)

    full = np.array([tc(d, s, b, m)
                     for d in grids["delta_xx"] for s in grids["gamma_scale"]
                     for b in grids["b_e"] for m in grids["mu"]])
    width_full = full.max() - full.min()

    mids = {k: g[len(g) // 2] for k, g in grids.items()}
    gains = {}
    for k in grids:
        gr = {kk: (np.array([mids[kk]]) if kk == k else g) for kk, g in grids.items()}
        sub = np.array([tc(d, s, b, m)
                        for d in gr["delta_xx"] for s in gr["gamma_scale"]
                        for b in gr["b_e"] for m in gr["mu"]])
        gains[k] = width_full - (sub.max() - sub.min())

    ranked = sorted(gains.items(), key=lambda kv: -kv[1])
    print(f"\n(iv) measurement-priority ranking (staged-device T_c envelope width "
          f"= {width_full:.0f} K over the [A] box):")
    labels = {"delta_xx": "Delta_XX distribution (in-house PL)",
              "gamma_scale": "Gamma(T) on InP/GaAsP (in-house PL)",
              "b_e": "injection background b_e(I,T) (EL spectra)",
              "mu": "operating loading mu (power series)"}
    for k, v in ranked:
        print(f"    measure {labels[k]:42s} -> envelope narrows by {v:5.1f} K")
    return {"width_full": width_full, "gains": gains, "labels": labels,
            "ranked": ranked}


if __name__ == "__main__":
    m = run_map()
    s = run_staged()
    a = run_aperture()
    r = run_ranking()
    phase3_packet_figure(m, s, a, r, ROOT / "out" / "phase3")
    print(f"\npacket bundle -> {ROOT / 'out' / 'phase3'}  (figure + CSV per panel)")
