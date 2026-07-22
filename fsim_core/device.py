"""Device assembly: one editable description of a complete device, evaluated
through the validated F-series chain (Modules A-E).

The device is NOT an arbitrary network -- the physics defines one signal chain
(dot -> barriers/retention -> cavity -> filter -> detector) with drive and
thermal attached. DeviceDesign holds the configurable blocks of that chain;
evaluate() runs the whole thing and returns curves + scalars.

Headless by construction (three-layer rule): fsim_gui/designer.py is a thin
client over this module. Designs round-trip to YAML in cards/ so every GUI
session stays reproducible from its file.

HONESTY: a design evaluated at point values yields conditional numbers, not
predictions -- unmeasured inputs keep their [A] tags and every result carries
the widest tag of its chain. The envelope treatment lives in run_phase3.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import yaml

from .cavity import purcell_eff, tracking_detuning
from .integrator import g2_from, retention
from .loading import (
    BackgroundChannel,
    aperture_g2,
    b_injection,
    f1b_g2,
    loading_probs,
    n_window_competitors,
)
from .spectral import epsilon, gamma_of_T
from .thermal import Layer, Stack, t_junction

_ROOT = Path(__file__).resolve().parents[1]


def class_proxy_params() -> dict:
    """Arsenide-class Gamma/retention parameters from the V-a joint fit [A as
    proxy for unmeasured platforms]."""
    p = json.loads((_ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]
    return {k: p[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo", "r_xx",
                              "a_esc", "E_a", "b_p", "E_b", "b0", "beta")}


@dataclass
class DotBlock:
    delta_xx: float = 3.5        # meV [A on InP/GaAsP]
    gamma_scale: float = 1.0     # multiplies the class Gamma(T) proxy [A]
    r_xx: float = 0.72           # Gamma_XX/Gamma_X [DR class]


@dataclass
class RetentionBlock:
    # class-proxy Arrhenius; override any field [A]
    a_esc: float = 0.0
    E_a: float = 0.0             # meV; 0 -> use class proxy
    b_p: float = 0.0
    E_b: float = 35.0
    b0: float = 0.0
    beta: float = 0.0


@dataclass
class DriveBlock:
    V: float = 1.9               # V
    I_uA: float = 10.0
    duty: float = 1.0
    mu: float = 0.5              # cap-2 loading per pulse [A]
    b_e: float = 0.02            # total injection background, signal units [A]
    b_e_m: float = 1.5           # current exponent of b_e [A]
    b_e_Eact: float = 100.0      # meV activation of the WL/barrier EL share [A]
    I_ref_uA: float = 10.0


@dataclass
class ThermalBlock:
    mesa_diameter_um: float = 1.0
    T_hs: float = 77.0
    layers: list = field(default_factory=lambda: [
        {"name": "epi/cladding", "t_um": 1.5, "k300": 15.0, "alpha": 0.5, "spread": False},
    ])
    substrate: dict = field(default_factory=lambda: {"name": "GaAs", "k300": 55.0, "alpha": 1.25})


@dataclass
class CavityBlock:
    enabled: bool = False
    type: str = "planar"         # "planar" (F6 resonant filter/gain) | "sin_waveguide"
    kappa: float = 1.0           # meV [A -> MEEP/COMSOL]; unused when type=sin_waveguide
    T_track: float = 120.0       # tracking-rule target (K)
    E_X0: float = 1.88           # eV cryogenic X
    dEdT_cav: float = -0.04      # meV/K
    F_P: float = 10.0            # [A]; unused when type=sin_waveguide
    G: float = 8.0               # collection gain on signal [A]; unused when type=sin_waveguide
    beta_sin: float = 1.0        # SiN evanescent coupling [A]; Lemma 1: brightness ONLY,
                                  # never eps/rho/g2/T_c (applied inline in evaluate())


@dataclass
class FilterBlock:
    enabled: bool = True
    auto_w: bool = True          # w = Gamma(T_j) operating convention [A]
    w: float = 2.0               # meV, used when auto_w = False
    dx: float = 0.0              # X offset from window center (meV)


@dataclass
class ApertureBlock:
    density_cm2: float = 7.0e8
    diameter_um: float = 1.0
    sigma_inh: float = 40.0
    comp_brightness: float = 0.3


@dataclass
class DeviceDesign:
    name: str = "my-device"
    dot: DotBlock = field(default_factory=DotBlock)
    ret: RetentionBlock = field(default_factory=RetentionBlock)
    drive: DriveBlock = field(default_factory=DriveBlock)
    thermal: ThermalBlock = field(default_factory=ThermalBlock)
    cavity: CavityBlock = field(default_factory=CavityBlock)
    filter: FilterBlock = field(default_factory=FilterBlock)
    aperture: ApertureBlock = field(default_factory=ApertureBlock)

    # ---- YAML round-trip (same reproducibility rule as the cards)
    def save(self, path):
        doc = {"meta": {"name": self.name, "role": "device-design",
                        "note": "editable device description; evaluate with "
                                "fsim_core.device.evaluate"},
               "design": asdict(self)}
        Path(path).write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    @staticmethod
    def load(path) -> "DeviceDesign":
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8"))["design"]
        return DeviceDesign(
            name=d.get("name", "my-device"),
            dot=DotBlock(**d["dot"]), ret=RetentionBlock(**d["ret"]),
            drive=DriveBlock(**d["drive"]), thermal=ThermalBlock(**d["thermal"]),
            cavity=CavityBlock(**d["cavity"]), filter=FilterBlock(**d["filter"]),
            aperture=ApertureBlock(**d["aperture"]),
        )


def _stack(th: ThermalBlock) -> Stack:
    return Stack(
        layers=[Layer(L["t_um"] * 1e-6, L["k300"], L.get("alpha", 1.25),
                      spread=bool(L.get("spread", False))) for L in th.layers],
        k_sub300=th.substrate["k300"], sub_alpha=th.substrate.get("alpha", 1.25),
    )


def evaluate(design: DeviceDesign, T_grid=None) -> dict:
    """Run the full chain. Returns {'curves': {...}, 'scalars': {...}}."""
    d = design
    proxy = class_proxy_params()
    gp = {k: proxy[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo")}
    rp = {k: (getattr(d.ret, k) or proxy[k]) for k in ("a_esc", "E_a", "b_p", "b0", "beta")}
    rp["E_b"] = d.ret.E_b or proxy["E_b"]

    a = 0.5 * d.thermal.mesa_diameter_um * 1e-6
    st = _stack(d.thermal)
    P = d.drive.duty * d.drive.I_uA * 1e-6 * d.drive.V
    chan = [BackgroundChannel("injection", A=d.drive.b_e, m=d.drive.b_e_m,
                              E_act=d.drive.b_e_Eact, I_ref=d.drive.I_ref_uA)]
    # SiN evanescent coupling: no resonant line (no kappa acceptance, no G boost);
    # beta_sin is applied to brightness only, below -- never to eps/rho/g2/T_c (Lemma 1).
    sin_mode = d.cavity.enabled and d.cavity.type == "sin_waveguide"

    def one(T_hs):
        Tj = t_junction(P, a, st, T_hs)
        if not np.isfinite(Tj):
            return dict(Tj=np.inf, gam=np.nan, eps=np.nan, rho=np.nan,
                        g2=np.nan, t_x=np.nan, runaway=True)
        gam = d.dot.gamma_scale * float(gamma_of_T(Tj, **gp))
        kappa = d.cavity.kappa if (d.cavity.enabled and not sin_mode) else None
        dx = d.filter.dx
        if d.cavity.enabled and not sin_mode:  # sin waveguide has no mode to track
            dx = 1e3 * float(tracking_detuning(Tj, d.cavity.E_X0, d.cavity.T_track,
                                               "GaAs", d.cavity.dEdT_cav * 1e-3))
        w = None
        if d.filter.enabled:
            w = gam if d.filter.auto_w else d.filter.w
        spec = epsilon(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=w, kappa=kappa, dx=dx)
        S = float(retention(Tj, rp["a_esc"], rp["E_a"], rp["b_p"], rp["E_b"]))
        B = rp["b0"] + rp["beta"] * (1.0 - S) + float(
            b_injection(chan, d.drive.I_uA, Tj))
        G = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
        rho = G * S / (G * S + B)
        g2_dot = float(f1b_g2(d.drive.mu, spec.eps)) if d.drive.mu > 0 else spec.eps
        return dict(Tj=Tj, gam=gam, eps=spec.eps, rho=rho,
                    g2=g2_from(g2_dot, rho), t_x=spec.t_x, runaway=False)

    Ts = np.asarray(T_grid if T_grid is not None else np.linspace(4.0, 350.0, 120))
    rows = [one(T) for T in Ts]
    curves = {
        "T_hs": Ts,
        "Tj": np.array([r["Tj"] for r in rows]),
        "g2": np.array([r["g2"] for r in rows]),
        "eps": np.array([r["eps"] for r in rows]),
        "rho2": np.array([r["rho"] for r in rows]) ** 2,
        "gamma": np.array([r["gam"] for r in rows]),
    }

    op = one(d.thermal.T_hs)
    # T_c: first heatsink temperature where g2 crosses 0.5 (above the g2 minimum)
    g2c = curves["g2"]
    Tc = np.nan
    if np.isfinite(g2c).any():
        imin = int(np.nanargmin(g2c))
        for i in range(imin, len(Ts) - 1):
            if np.isfinite(g2c[i]) and np.isfinite(g2c[i + 1]) and g2c[i] < 0.5 <= g2c[i + 1]:
                Tc = float(np.interp(0.5, [g2c[i], g2c[i + 1]], [Ts[i], Ts[i + 1]]))
                break

    _, P1, P2 = loading_probs(d.drive.mu)
    w_ap = op["gam"] if d.filter.auto_w or not d.filter.enabled else d.filter.w
    if np.isfinite(w_ap):
        Nw = float(n_window_competitors(
            d.aperture.density_cm2, np.pi * (d.aperture.diameter_um / 2) ** 2,
            w_ap, d.aperture.sigma_inh))
        n_comp = max(int(round(Nw)), 0)
        ap_pen = float(aperture_g2([1.0] + [d.aperture.comp_brightness] * n_comp)) \
            if n_comp else 0.0
    else:
        Nw, ap_pen = float("nan"), 0.0  # runaway: no meaningful operating window

    gain_factor = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
    beta_factor = d.cavity.beta_sin if sin_mode else 1.0
    scalars = {
        "T_j_op": op["Tj"], "dT_J": op["Tj"] - d.thermal.T_hs,
        "runaway": bool(op["runaway"]),
        "gamma_op": op["gam"], "eps_op": op["eps"], "rho_op": op["rho"],
        "g2_op": op["g2"], "t_x_op": op["t_x"],
        "brightness_per_pulse": float((P1 + P2)) * (op["t_x"] if np.isfinite(op["t_x"]) else 0.0)
        * gain_factor * beta_factor,
        "T_c": Tc,
        "F_eff": float(purcell_eff(d.cavity.F_P, d.cavity.kappa, op["gam"]))
        if (d.cavity.enabled and not sin_mode and np.isfinite(op["gam"])) else np.nan,
        "N_w": Nw, "aperture_g2_penalty": ap_pen,
        "tag_chain": "[A]",  # unmeasured inputs are always in the chain today
    }
    return {"curves": curves, "scalars": scalars}
