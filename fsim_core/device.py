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

import copy
import itertools
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
    f8_g2,
    f8b_thin_fano,
    gamma_eff,
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
    mode: str = "EL"             # "EL" (electrical: dg_inj + injection background
                                  # both active) | "PL" (optical excitation: both
                                  # disabled) [A]
    dg_inj: float = 0.0          # meV; F5' injection broadening amplitude [A]
    p_inj: float = 1.0           # F5' injection broadening current exponent [A]
    F_p: float = 1.0             # F8 pump Fano factor (1=Poisson) [A]
    eta_capture: float = 1.0     # F8b dot-capture fraction (thinning) [A]
    C_dep_pF: float = 10.0       # pF; F9 depletion capacitance, informational
                                  # (drives granularity_N/island_radius_nm only
                                  # -- never wired into g2) [E]


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
        # F5' injection broadening (v1.1): EL-mode-only, applied BEFORE the
        # auto-w filter width (below) and BEFORE epsilon.
        if d.drive.mode != "PL" and d.drive.dg_inj:
            gam = gamma_eff(gam, d.drive.dg_inj, d.drive.I_uA / d.drive.I_ref_uA, d.drive.p_inj)
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
        # PL mode: optical excitation, no injection-current background channel.
        inj_bg = float(b_injection(chan, d.drive.I_uA, Tj)) if d.drive.mode != "PL" else 0.0
        B = rp["b0"] + rp["beta"] * (1.0 - S) + inj_bg
        G = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
        rho = G * S / (G * S + B)
        if d.drive.mu > 0:
            if d.drive.F_p != 1.0:
                # F8b thinning: effective pump Fano factor AT the dot, after
                # dot-capture; F8: moment-matched (mu, Fano) cap-2 loading.
                F_eff = f8b_thin_fano(d.drive.eta_capture, d.drive.F_p)
                g2_dot = float(f8_g2(d.drive.mu, F_eff, spec.eps))
            else:
                # F_p == 1.0 (default): keep the original f1b_g2 path EXACTLY
                # -- f8_g2(mu, 1, eps) is not bit-identical to f1b_g2(mu, eps)
                # at finite mu (different, equally valid cap-2 conventions;
                # see loading.py module docstring), so every pre-v1.1 result
                # stays reproducible unless F_p is explicitly set != 1.0.
                g2_dot = float(f1b_g2(d.drive.mu, spec.eps))
        else:
            g2_dot = spec.eps
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


# ------------------------------------------------------------- envelope mode (D2)

_ENVELOPE_CAP = 4096
_ENV_CURVES = ("g2", "eps", "rho2", "Tj", "gamma")
_ENV_SCALARS = ("T_c", "g2_op", "eps_op", "rho_op", "dT_J")


def _set_path(design: DeviceDesign, path: str, value) -> None:
    block_name, field_name = path.split(".", 1)
    setattr(getattr(design, block_name), field_name, value)


def _value_set(spec, mode: str) -> tuple:
    """One parameter's sample set from a ranged spec: a 2-tuple (lo, hi) under
    mode="extremes" becomes the 2-level factorial [lo, hi]; anything else
    (len != 2, or an explicit 2-point grid the caller wants used verbatim) is
    taken as given."""
    seq = tuple(spec)
    if len(seq) == 2:
        if mode != "extremes":
            raise ValueError(f"evaluate_envelope: unsupported mode {mode!r} for a "
                             "(lo, hi) pair -- only mode='extremes' is implemented")
        return seq
    return seq


def _cartesian_eval(design: DeviceDesign, paths: list, value_sets: dict, T_grid,
                    raise_if_all_dropped=True):
    """Evaluate design at every point of the cartesian product of value_sets
    (one value set per path, in `paths` order); returns (curves_list,
    scalars_list), one entry per sample. raise_if_all_dropped=False lets a
    caller (the tornado's collapsed runs) receive an all-NaN result instead
    of an error when the whole collapsed box is outside the F8 domain."""
    combos = list(itertools.product(*(value_sets[p] for p in paths))) if paths else [()]
    curves_list, scalars_list = [], []
    n_dropped = 0
    for values in combos:
        d = copy.deepcopy(design)
        for path, v in zip(paths, values):
            _set_path(d, path, v)
        try:
            res = evaluate(d, T_grid=T_grid)
        except ValueError:
            # F8 domain violation (mu < 1 - F_eff): no cap-2 loading
            # distribution with these moments exists -- same category as
            # thermal runaway, so the sample contributes a gap, not a crash
            # (Opus v1.1 finding). The point-evaluation path keeps raising.
            n_T = len(T_grid) if T_grid is not None else 120
            nanarr = np.full(n_T, np.nan)
            res = {"curves": {k: nanarr.copy() for k in _ENV_CURVES},
                   "scalars": {k: float("nan") for k in _ENV_SCALARS}}
            res["scalars"]["runaway"] = False
            n_dropped += 1
        curves_list.append(res["curves"])
        scalars_list.append(res["scalars"])
    if n_dropped == len(combos) and raise_if_all_dropped:
        raise ValueError(
            "evaluate_envelope: every sample in the ranged box violates the F8 "
            "loading domain (mu >= 1 - F_eff). Narrow the mu range or raise "
            "F_p/eta_capture.")
    return curves_list, scalars_list


def _curve_bands(curves_list: list) -> dict:
    bands = {}
    for name in _ENV_CURVES:
        stacked = np.array([c[name] for c in curves_list], dtype=float)
        # nanmin/nanmax warn (via warnings.warn, not an FP flag) on all-NaN
        # slices -- e.g. a fully-runaway box; guard explicitly (Opus D2 finding)
        any_finite = np.isfinite(stacked).any(axis=0)
        lo = np.full(stacked.shape[1], np.nan)
        hi = np.full(stacked.shape[1], np.nan)
        if any_finite.any():
            lo[any_finite] = np.nanmin(stacked[:, any_finite], axis=0)
            hi[any_finite] = np.nanmax(stacked[:, any_finite], axis=0)
        bands[name] = (lo, hi)
    return bands


def _scalar_bands(scalars_list: list) -> dict:
    out = {}
    for name in _ENV_SCALARS:
        vals = np.array([s[name] for s in scalars_list], dtype=float)
        finite = vals[np.isfinite(vals)]
        if finite.size == 0:
            out[name] = (float("nan"), float("nan"))
        else:
            out[name] = (float(finite.min()), float(finite.max()))
    return out


def evaluate_envelope(design: DeviceDesign, ranged: dict, T_grid=None,
                      mode: str = "extremes") -> dict:
    """Headless envelope evaluation: sweep `ranged` inputs instead of pinning
    them to a point, per the honesty discipline (unmeasured [A] inputs are
    swept, never averaged into a single "prediction").

    ranged: {META path ("dot.delta_xx", ...): (lo, hi) | explicit value list}.
    A 2-tuple is expanded to a 2-level factorial [lo, hi] under the default
    mode="extremes" -- HONEST CAVEAT: this only bounds the true envelope when
    the response is monotone in that parameter across the box; a non-monotone
    response can have interior extrema the 2-level factorial misses. Anything
    else (an explicit list/tuple of length != 2) is used verbatim as that
    parameter's grid, e.g. to include known interior points.

    Samples the cartesian product over every ranged path's value set (capped
    at 4096 evaluations -- raises ValueError above that). Returns:
      bands:        {curve_name: (lo_array, hi_array)} pointwise NaN-safe
                    min/max over samples, for g2/eps/rho2/Tj/gamma.
      mid:          curves dict at the all-midpoint sample (every ranged path
                    pinned to the midpoint of its value set).
      scalar_bands: {scalar_name: (lo, hi)} NaN-safe (all-NaN -> (nan, nan))
                    for T_c/g2_op/eps_op/rho_op/dT_J.
      tornado:      {path: width_reduction} -- for each ranged path, the drop
                    in the baseline scalar's band width (T_c if any sample has
                    a finite T_c, else g2_op) when that one path alone is
                    collapsed to its midpoint (others stay ranged). Larger is
                    a higher measurement priority.
      n_samples:    number of cartesian-product evaluations.

    ranged={} degenerates exactly to evaluate(design, T_grid): one sample, so
    bands/mid/scalar_bands all collapse to that single evaluate() call.
    """
    Ts = np.asarray(T_grid if T_grid is not None else np.linspace(4.0, 350.0, 120))
    paths = list(ranged)
    value_sets = {p: _value_set(ranged[p], mode) for p in paths}

    n_samples = 1
    for p in paths:
        n_samples *= len(value_sets[p])
    if n_samples > _ENVELOPE_CAP:
        raise ValueError(
            f"evaluate_envelope: {len(paths)} ranged parameters -> {n_samples} "
            f"cartesian-product evaluations, over the cap of {_ENVELOPE_CAP}. "
            "Sweep fewer parameters at once, or use coarser/explicit grids."
        )

    curves_list, scalars_list = _cartesian_eval(design, paths, value_sets, Ts)
    bands = _curve_bands(curves_list)
    scalar_bands = _scalar_bands(scalars_list)

    mid_values = {p: 0.5 * (min(value_sets[p]) + max(value_sets[p])) for p in paths}
    mid_design = copy.deepcopy(design)
    for p in paths:
        _set_path(mid_design, p, mid_values[p])
    try:
        mid = evaluate(mid_design, T_grid=Ts)["curves"]
    except ValueError:
        # midpoint itself outside the F8 loading domain (possible for an
        # explicit user range straddling the floor): NaN reference curve,
        # bands remain valid (Opus v1.1 finding, mid-path leg)
        n_T = len(Ts) if Ts is not None else 120
        mid = {k: np.full(n_T, np.nan) for k in _ENV_CURVES}

    baseline = "T_c" if np.isfinite(scalar_bands["T_c"][0]) else "g2_op"
    full_lo, full_hi = scalar_bands[baseline]
    full_width = full_hi - full_lo
    tornado = {}
    for p in paths:
        collapsed = dict(value_sets)
        collapsed[p] = (mid_values[p],)
        _, sc_list = _cartesian_eval(design, paths, collapsed, Ts,
                                     raise_if_all_dropped=False)
        c_lo, c_hi = _scalar_bands(sc_list)[baseline]
        # collapsed midpoint may sit outside the F8 domain -> width undefined
        collapsed_width = c_hi - c_lo
        tornado[p] = full_width - collapsed_width

    return {"bands": bands, "mid": mid, "scalar_bands": scalar_bands,
            "tornado": tornado, "n_samples": n_samples}
