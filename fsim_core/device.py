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
from .drive_mech import mech_from_card
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
from .qd_gf import PhononParams, ibm_purcell_transmission, ibm_transmission
from .spectral import SpectralResult, epsilon, epsilon2, gamma_of_T
from .thermal import Layer, Stack, t_junction
from .linewidth import LinewidthParams, gamma_anchor
from . import cw_g2, dot_levels, materials, transport, waveguide

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
    lineshape: str = "lorentzian"  # T2 tier: "lorentzian" (legacy, bit-
                                   # identical) | "ibm" (polaron ZPL+sidebands;
                                   # the Lorentzian-optimism correction)
    phonon: dict = field(default_factory=dict)  # qd_gf.PhononParams overrides
                                   # (geometry l_xy_nm/l_z_nm, materials) [A/DR]
    linewidth: str = "class"    # "class" proxy | "anchored" [A; Matsuda, PRB 2001]
    gamma0: float = 0.25         # meV [V] Bommer et al., JAP 110, 063108 (2011)
    a_ac: float = 2.0e-3         # meV/K [V] Ortner et al., PRB 70, 201301 (2004)
    E_LO: float = 43.0           # meV [V] Ioffe NSM, InP LO phonon
    gamma300: float = 12.0       # meV [A] Matsuda et al., PRB 63, 121304 (2001) class anchor


@dataclass
class RetentionBlock:
    # class-proxy Arrhenius; override any field [A]
    a_esc: float = 0.0
    E_a: float = 0.0             # meV; 0 -> use class proxy
    b_p: float = 0.0
    E_b: float = 35.0
    b0: float = 0.0
    beta: float = 0.0
    mode: str = "proxy"         # "proxy" keeps class fit; "confinement" derives levels [E]
    preset: str = ""
    system: dict = field(default_factory=dict)  # inline DotSystem stack (mutually
                                  # exclusive with preset); never carries overrides
    tau_rad_ns: float = 1.0      # ns [E] dot_levels.retention_params convention
    channel: str = "pair_half"  # Gelinas et al., arXiv:0910.0480 [DR]
    overrides: dict = field(default_factory=dict)  # confinement mode only: explicit
                                  # a_esc/E_a/b_p/E_b/b0/beta replacements for the
                                  # retention_params()-derived values [A]; separate
                                  # from `system` so a preset stack can still be
                                  # combined with explicit overrides. An explicit
                                  # 0.0 here always wins (never a truthy-or fallback).


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
    mechanism: str = ""           # T1 drive-mechanism library: "" = legacy path
                                  # (bit-identical); else a drive_mech name
                                  # (poisson-rail/quiet-rail/pulsed/set-*/rti)
    mech_params: dict = field(default_factory=dict)  # per-mechanism params [A]
    diode: dict = field(default_factory=dict)  # transport preset plus Diode overrides [E/DR]
    n_dot_cm2: float = 0.0       # cm^-2; 0 uses aperture density [A]
    cw: bool = False              # opt-in CW g2(tau)/IRF diagnostics (docs/
                                  # rt_edge_contract.md); requires mode=
                                  # "EL-transport" (the CW pump rate is
                                  # transport.loading.r_dot, never inferred
                                  # from the dimensionless pulsed mu) [A]
    cw_irf_fwhm_ps: float = 500.0  # detector IRF FWHM, ps; Reischle 2008
                                  # class 0.5 ns [E]
    cw_irf_shape: str = "gaussian"  # cw_g2.IRF_SHAPES
    cw_pump_ratio: float = 1.0    # X->XX secondary CW pump ratio [A]
    cw_tau_max_ns: float = 10.0   # CW g2(tau) window half-width, ns [A]


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
    purcell_wire: bool = False   # R2 (roadmap 13), opt-in: wire F_eff into
                                  # (a) retention (radiative rate x rate_mult
                                  # -> escape ratios divided) and (b) the IBM
                                  # ZPL redistribution (Z_eff funneling).
                                  # [E->analytic-1D, confirm: SIM-B]. False =
                                  # legacy bit-identical (T_c conservative).


@dataclass
class FilterBlock:
    enabled: bool = True
    auto_w: bool = True          # w = Gamma(T_j) operating convention [A]
    w: float = 2.0               # meV, used when auto_w = False
    dx: float = 0.0              # X offset from window center (meV)
    track: str = "mode"          # "mode": slit follows the cavity mode (legacy,
                                 #   bit-identical) | "hold": F6 slit-held --
                                 #   slit centered ON X at each operating T
                                 #   (lab monochromator convention) while the
                                 #   cavity mode walks per dn/dT. The two
                                 #   coincide at T_track by construction.
    track_material: str = ""     # "" = legacy hard-coded GaAs Varshni tracking
                                 #   (cavity.tracking_detuning, bit-identical);
                                 #   "dot" | "matrix" opts into stack-based
                                 #   tracking via materials.bandgap() on the
                                 #   named layer of the shared inline stack
                                 #   (ret.preset / ret.system) -- never GaAs
                                 #   by default for a non-GaAs stack [DR].


@dataclass
class ApertureBlock:
    density_cm2: float = 7.0e8
    diameter_um: float = 1.0
    sigma_inh: float = 40.0
    comp_brightness: float = 0.3
    compose: bool = False        # False (legacy): aperture_g2_penalty stays an
                                 #   informational scalar (integer-rounded
                                 #   competitor count), never folded into g2.
                                 #   True: continuous [A] Poisson-competitor-
                                 #   bath composition into g2/g2_cw0 (docs/
                                 #   rt_edge_contract.md "Aperture assumptions").


@dataclass
class EmissionBlock:
    type: str = "none"           # "none" (legacy) | "edge" opts into
                                 #   waveguide.edge_emission out-coupling
                                 #   (Lemma 1: brightness only, never g2/eps/rho).
    ridge_width_nm: float = 2000.0
    etch_depth_nm: float = 1200.0
    L_um: float = 500.0
    NA: float = 0.5
    alpha_cm: float = 5.0        # ridge propagation loss, class estimate [E]
    R_back: float | None = None  # None: two facets share emission equally [A]
    coating: dict = field(default_factory=dict)  # facet_transmission() override [A]
    core_half_nm: float = 148.0  # matrix (well) layer thickness each side of
                                 #   the dot plane, hkust_ridge_stack class [E]
    cladding_nm: float = 1000.0  # barrier/cladding thickness each side [E]
    lambda_nm: float = 0.0       # 0 (default): the actual confinement-derived
                                 #   X transition, dot_levels.levels(system).
                                 #   lambda_nm (contract wording). MATERIAL_
                                 #   EXTRA's optical tables are sparse (a
                                 #   handful of tabulated points per alloy), so
                                 #   most physically-derived InP-dot
                                 #   wavelengths land outside any table and
                                 #   are correctly reported ineligible
                                 #   (missing optics evidence) rather than
                                 #   guessed; an explicit >0 value here is a
                                 #   stated [A] design wavelength (e.g. the
                                 #   668 nm point where the phosphide stack's
                                 #   indices are actually tabulated) used
                                 #   instead of the natural transition.


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
    emission: EmissionBlock = field(default_factory=EmissionBlock)
    provenance: dict = field(default_factory=dict)  # persisted card-side evidence notes

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
            emission=EmissionBlock(**d.get("emission", {})),
            provenance=dict(d.get("provenance", {})),
        )


def _stack(th: ThermalBlock) -> Stack:
    return Stack(
        layers=[Layer(L["t_um"] * 1e-6, L["k300"], L.get("alpha", 1.25),
                      spread=bool(L.get("spread", False))) for L in th.layers],
        k_sub300=th.substrate["k300"], sub_alpha=th.substrate.get("alpha", 1.25),
    )


def _material_from_mapping(value):
    """Decode a card material without changing the public material API [DR]."""
    if isinstance(value, str):
        return materials.binary(value)
    if not isinstance(value, dict):
        raise ValueError("ret.system material must be a name or mapping")
    kind = value.get("kind", value.get("alloy", "binary"))
    if kind == "binary":
        return materials.binary(value["name"])
    constructors = {"GaAsP": materials.GaAsP, "GaInP": materials.GaInP,
                    "AlGaAs": materials.AlGaAs, "AlGaInP": materials.AlGaInP,
                    "InGaAs": materials.InGaAs}
    if kind not in constructors:
        raise ValueError(f"ret.system unknown material kind {kind!r}")
    kw = {k: v for k, v in value.items() if k not in ("kind", "alloy", "name", "tag", "source")}
    return constructors[kind](**kw)


def _retention_system(ret: RetentionBlock):
    if ret.preset and ret.system:
        raise ValueError("ret.preset and ret.system cannot both be set")
    if ret.preset:
        presets = dot_levels.class_presets()
        if ret.preset not in presets:
            raise ValueError(f"unknown ret.preset {ret.preset!r}")
        # The actual junction temperature, rather than a preset's cryogenic
        # display temperature, is installed by _confinement_params below.
        return presets[ret.preset]()
    if not ret.system:
        raise ValueError("ret.mode='confinement' requires ret.preset or ret.system")
    raw = ret.system
    geometry = dot_levels.DotGeometry(**raw.get("geometry", {}))
    return dot_levels.DotSystem(
        _material_from_mapping(raw["dot"]), _material_from_mapping(raw["matrix"]),
        _material_from_mapping(raw["barrier"]), _material_from_mapping(raw["substrate"]),
        T=float(raw.get("T", 300.0)), geometry=geometry,
        vbo_override_eV=dict(raw.get("vbo_override_eV", {})),
        eps_r=raw.get("eps_r"), name=raw.get("name", "inline-system"),
    )


def _diode_from_drive(drive: DriveBlock):
    raw = dict(drive.diode)
    preset = raw.pop("preset", "")
    for key in ("tau_pulse_ns", "tau_rad_ns", "w_meV", "dE_WL_meV",
                "eta_rad_matrix", "E_urbach_meV", "eta_total"):
        raw.pop(key, None)
    # Inline material spellings are deliberately decoded here: no card is
    # allowed to silently replace its stated stack with a preset [DR].
    for key in ("p_cladding", "n_cladding", "active", "barrier", "substrate"):
        if key in raw:
            raw[key] = _material_from_mapping(raw[key])
    if preset == "hkust":
        return transport.hkust_preset(**raw)
    if preset in ("red", "red_diode"):
        return transport.red_diode_preset(**raw)
    raise ValueError("drive.mode='EL-transport' requires diode.preset 'hkust' or 'red'")


def _confinement_params(ret: RetentionBlock, Tj: float) -> dict:
    """Resolve confinement at the temperature consumed by escape [DR].  Also
    carries the level table's own E_X_eV (resolved transition energy) for
    the item-2 sub-turn-on loading suppression -- the SAME confinement
    solve, never a second one."""
    system = _retention_system(ret)
    # DotSystem is mutable, so copy rather than mutating the card-derived
    # object.  This is also why evaluate leaves its input design untouched.
    system = copy.copy(system)
    system.T = float(Tj)
    lv = dot_levels.levels(system)
    params = dot_levels.retention_params(lv, ret.tau_rad_ns, ret.channel, verbose=False)
    params["E_X_eV"] = lv.E_X_eV
    return params


def _tracked_material(ret: RetentionBlock, layer: str) -> materials.Material:
    """Resolve filter.track_material ("dot" | "matrix") against the SAME
    shared inline stack (ret.preset / ret.system) used by confinement
    retention -- never an independently selected material [DR]."""
    if layer not in ("dot", "matrix"):
        raise ValueError(f"unknown filter.track_material {layer!r} (use 'dot' or 'matrix')")
    system = _retention_system(ret)
    return system.dot if layer == "dot" else system.matrix


def _aperture_lambda(density_cm2, aperture_um2, w_full_meV, sigma_inh_meV,
                     comp_brightness) -> tuple:
    """Continuous (unrounded) competitor count N_w and composed bath mean
    lambda = N_w * comp_brightness (docs/rt_edge_contract.md "Aperture
    assumptions"). NaN in, NaN out (no meaningful operating window)."""
    if not np.isfinite(w_full_meV):
        return float("nan"), float("nan")
    Nw = float(n_window_competitors(density_cm2, aperture_um2, w_full_meV, sigma_inh_meV))
    return Nw, Nw * comp_brightness


def _compose_aperture_g2(g_target, lam):
    """Continuous aperture composition [A]: the competitor bath is an
    aggregate Poisson process of mean lam = N_w * comp_brightness (N_w
    continuous, never rounded), independent of the target and of any other
    loss already folded into g_target (filter/capture/loading are upstream of
    this call, composed exactly once). Factorial-moment derivation, target
    brightness normalized to 1:

        <n(n-1)> = g_target*1^2  +  lam^2  +  2*1*lam     (target + Poisson
                                                            bath + independent
                                                            cross term)
        <n>^2   = (1 + lam)^2

        g_mix = (g_target + 2 lam + lam^2) / (1 + lam)^2

    Continuous in lam (no integer rounding/clipping/additive penalty); lam=0
    (N_w=0 or comp_brightness=0) returns g_target exactly."""
    lam = np.asarray(lam, dtype=float)
    return (np.asarray(g_target, dtype=float) + 2.0 * lam + lam * lam) / (1.0 + lam) ** 2


def _resolve_edge(ret: RetentionBlock, emission: EmissionBlock, Tj: float):
    """Convert the shared inline stack (ret.preset / ret.system) into
    waveguide.Layer objects at the actual emission wavelength (the
    confinement-derived X transition, dot_levels.levels(system).lambda_nm)
    and run waveguide.edge_emission's documented public API. Never uses the
    hkust_ridge_stack() convenience builder (that stack's own composition is
    a separate literature fixture, not inferred from any card's stack) and
    never touches private waveguide internals."""
    if not (ret.preset or ret.system):
        raise ValueError("emission.type='edge' requires ret.preset or ret.system "
                          "(the shared inline stack; no independently selected preset)")
    system = copy.copy(_retention_system(ret))
    system.T = float(Tj)
    lv = dot_levels.levels(system)
    lambda_nm = emission.lambda_nm if emission.lambda_nm > 0 else lv.lambda_nm
    g = system.geometry

    def n_at(mat):
        return materials.refractive_index(mat.label, lambda_nm)

    layers = [
        waveguide.Layer("barrier_lower", n_at(system.barrier), emission.cladding_nm),
        waveguide.Layer("matrix_lower", n_at(system.matrix), emission.core_half_nm),
        waveguide.Layer("dot", n_at(system.dot), max(g.height_nm, 0.1), True),
        waveguide.Layer("matrix_upper", n_at(system.matrix), emission.core_half_nm),
        waveguide.Layer("barrier_upper", n_at(system.barrier), emission.cladding_nm),
    ]
    edge = waveguide.edge_emission(
        layers, emission.ridge_width_nm, emission.etch_depth_nm, lambda_nm,
        emission.L_um, emission.NA, alpha_cm=emission.alpha_cm,
        R_back=emission.R_back, coating=(emission.coating or None))
    return edge, lambda_nm


def _transport_options(drive: DriveBlock, w_meV: float) -> dict:
    """Card-only transport knobs; units are ns and meV [A/E]."""
    raw = drive.diode
    return dict(tau_pulse_ns=float(raw.get("tau_pulse_ns", 1.0)),
                tau_rad_ns=float(raw.get("tau_rad_ns", 1.0)),
                w_meV=float(raw.get("w_meV", w_meV)),
                dE_WL_meV=float(raw.get("dE_WL_meV", 100.0)),
                eta_rad_matrix=float(raw.get("eta_rad_matrix", 0.1)),
                E_urbach_meV=raw.get("E_urbach_meV"),
                eta_total=float(raw.get("eta_total", 0.01)))


def evaluate(design: DeviceDesign, T_grid=None) -> dict:
    """Run the full chain. Returns {'curves': {...}, 'scalars': {...}}."""
    d = design
    if d.drive.mode == "PL" and d.drive.diode:
        raise ValueError("drive.mode='PL' cannot be used with a non-empty diode")
    if d.drive.mode == "EL-transport" and d.drive.dg_inj:
        raise ValueError("drive.mode='EL-transport' cannot be used with dg_inj")
    if d.drive.mode == "EL-transport" and d.drive.mechanism:
        raise ValueError("drive.mode='EL-transport' cannot be combined with "
                         "drive.mechanism (conflicting loading resolutions)")
    if d.dot.linewidth not in ("class", "anchored"):
        raise ValueError(f"unknown dot.linewidth {d.dot.linewidth!r}")
    if d.ret.mode not in ("proxy", "confinement"):
        raise ValueError(f"unknown ret.mode {d.ret.mode!r}")
    if d.drive.mode not in ("EL", "PL", "EL-transport"):
        raise ValueError(f"unknown drive.mode {d.drive.mode!r}")
    if d.emission.type not in ("none", "edge"):
        raise ValueError(f"unknown emission.type {d.emission.type!r}")
    if d.filter.track_material not in ("", "dot", "matrix"):
        raise ValueError(f"unknown filter.track_material {d.filter.track_material!r}")
    if d.emission.type == "edge" and d.cavity.enabled:
        # Reject conflicting SiN/edge and resonant-cavity collection (docs/
        # rt_edge_contract.md): coexistence needs an explicit supported
        # leakage accounting this tier does not implement.
        raise ValueError("emission.type='edge' cannot be combined with "
                         "cavity.enabled (no supported cavity-leakage accounting)")
    if d.drive.cw and d.drive.mode != "EL-transport":
        # The CW pump rate is transport.loading.r_dot; never inferred from
        # the dimensionless pulsed mu (docs/rt_edge_contract.md).
        raise ValueError("drive.cw=True requires drive.mode='EL-transport'")
    proxy = class_proxy_params() if (d.dot.linewidth == "class" or d.ret.mode == "proxy") else {}
    gp = ({k: proxy[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo")}
          if d.dot.linewidth == "class" else None)
    rp = ({k: (getattr(d.ret, k) or proxy[k]) for k in ("a_esc", "E_a", "b_p", "b0", "beta")}
          if d.ret.mode == "proxy" else {})
    if d.ret.mode == "proxy":
        rp["E_b"] = d.ret.E_b or proxy["E_b"]
    retention_note = "class-proxy Arrhenius fit [A]"
    diode = _diode_from_drive(d.drive) if d.drive.mode == "EL-transport" else None

    a = 0.5 * d.thermal.mesa_diameter_um * 1e-6
    st = _stack(d.thermal)
    P = d.drive.duty * d.drive.I_uA * 1e-6 * d.drive.V
    chan = [BackgroundChannel("injection", A=d.drive.b_e, m=d.drive.b_e_m,
                              E_act=d.drive.b_e_Eact, I_ref=d.drive.I_ref_uA)]
    # SiN evanescent coupling: no resonant line (no kappa acceptance, no G boost);
    # beta_sin is applied to brightness only, below -- never to eps/rho/g2/T_c (Lemma 1).
    sin_mode = d.cavity.enabled and d.cavity.type == "sin_waveguide"

    def one(T_hs):
        # The legacy branch deliberately retains its one-shot P=IV calculation.
        # Transport uses the diode's junction power in a short fixed-point loop.
        Tj = t_junction(P, a, st, T_hs)
        if diode is not None:
            if d.drive.I_uA <= 0:
                return dict(Tj=np.nan, gam=np.nan, eps=np.nan, rho=np.nan,
                            g2=np.nan, t_x=np.nan, runaway=False, mu=np.nan,
                            eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                            g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                            cw_gamma_X_ns=np.nan, cw_rho=np.nan,
                            invalid_reason="EL-transport requires positive current")
            converged = False
            for _ in range(12):
                trial = transport.evaluate_injection(
                    diode=diode, I_uA=d.drive.I_uA, T=Tj,
                    n_dot_cm2=d.drive.n_dot_cm2 or d.aperture.density_cm2,
                    aperture_um2=np.pi * (d.aperture.diameter_um / 2) ** 2,
                    **_transport_options(d.drive, 1.0))
                next_Tj = t_junction(d.drive.duty * trial.P_junction_W, a, st, T_hs)
                if abs(next_Tj - Tj) < 1e-10:
                    converged = True
                    break
                Tj = next_Tj
            if not converged:
                return dict(Tj=np.nan, gam=np.nan, eps=np.nan, rho=np.nan,
                            g2=np.nan, t_x=np.nan, runaway=False, mu=np.nan,
                            eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                            g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                            cw_gamma_X_ns=np.nan, cw_rho=np.nan,
                            invalid_reason="transport self-heating did not converge")
        if not np.isfinite(Tj):
            return dict(Tj=np.inf, gam=np.nan, eps=np.nan, rho=np.nan,
                        g2=np.nan, t_x=np.nan, runaway=True, mu=np.nan,
                        eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                        g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                        cw_gamma_X_ns=np.nan, cw_rho=np.nan,
                        invalid_reason="thermal runaway")
        gam_base = (float(gamma_anchor(Tj, LinewidthParams(d.dot.gamma0, d.dot.a_ac,
                    d.dot.E_LO, d.dot.gamma300))) if d.dot.linewidth == "anchored"
                    else float(gamma_of_T(Tj, **gp)))
        gam = d.dot.gamma_scale * gam_base
        # F5' injection broadening (v1.1): EL-mode-only, applied BEFORE the
        # auto-w filter width (below) and BEFORE epsilon.
        if d.drive.mode != "PL" and d.drive.dg_inj:
            gam = gamma_eff(gam, d.drive.dg_inj, d.drive.I_uA / d.drive.I_ref_uA, d.drive.p_inj)
        kappa = d.cavity.kappa if (d.cavity.enabled and not sin_mode) else None
        dx = d.filter.dx
        if d.cavity.enabled and not sin_mode:  # sin waveguide has no mode to track
            if d.filter.track_material:
                # Stack-based tracking [DR]: the E_X0/T=0 anchor cancels
                # between the two bandgap differences below, so only the
                # named layer's OWN Varshni shape enters -- never GaAs unless
                # the stack's own dot/matrix material happens to be GaAs.
                mat = _tracked_material(d.ret, d.filter.track_material)
                dx = (1e3 * (materials.bandgap(mat, Tj) - materials.bandgap(mat, d.cavity.T_track))
                      - d.cavity.dEdT_cav * (Tj - d.cavity.T_track))
            else:
                dx = 1e3 * float(tracking_detuning(Tj, d.cavity.E_X0, d.cavity.T_track,
                                                   "GaAs", d.cavity.dEdT_cav * 1e-3))
        elif d.filter.track_material:
            # item 4 (council review 2026-09-05): track_material was gated
            # behind cavity.enabled, so it did nothing on any cavity-less
            # design (both edge-emitter cards) -- g2_op was bit-identical
            # for track_material in {dot, matrix, ""}.  With no cavity mode
            # to net against, the filter window still follows the dot's OWN
            # Varshni walk relative to cavity.T_track (the only reference
            # temperature this block carries); no dEdT_cav subtraction here
            # (there is no cavity mode drifting to subtract).  Legacy
            # default path (track_material == "") is untouched: dx stays
            # d.filter.dx exactly as before.
            mat = _tracked_material(d.ret, d.filter.track_material)
            dx = 1e3 * (materials.bandgap(mat, Tj) - materials.bandgap(mat, d.cavity.T_track))
        w = None
        if d.filter.enabled:
            w = gam if d.filter.auto_w else d.filter.w
        # Slit-held tracking (T-1): slit centered on X (dx_w = user offset
        # only), cavity at the physical mode walk (dx stays the tracking
        # detuning). Under "mode" (legacy) both centers share dx exactly.
        held = (d.filter.track == "hold" and d.filter.enabled
                and d.cavity.enabled and not sin_mode)
        dx_w = d.filter.dx if held else dx
        # R2 opt-in: overlap-penalized Purcell at this operating point,
        # scaled by the cavity Lorentzian at the X detuning (a mode far off
        # the line enhances nothing -- the factor that kills the 77 K +
        # 120 K-tracked wiring honestly).
        wire = (d.cavity.purcell_wire and d.cavity.enabled and not sin_mode)
        if wire:
            det = (kappa / 2.0) ** 2 / (dx**2 + (kappa / 2.0) ** 2)
            F_eff_op = 1.0 + (float(purcell_eff(d.cavity.F_P, kappa, gam))
                              - 1.0) * det
        else:
            F_eff_op = 1.0
        rate_mult = 1.0
        if d.dot.lineshape == "ibm":
            # T2 tier: sideband-aware transmissions. Same acceptance stack and
            # delta conventions as spectral.epsilon; the XX line reuses the X
            # exciton's PhononParams (single-form-factor approximation, [E] --
            # the XX phonon coupling differs at the O(1)-factor level, a
            # refinement that waits for data). ZPL widths carry the standing
            # Gamma(T) model unchanged (see qd_gf.effective_gamma_zpl note).
            pp = PhononParams(**d.dot.phonon)
            if wire:
                # R2: Purcell funnels emission into the ZPL (Z_eff) and
                # multiplies the radiative rate; X line only (the XX line's
                # own Purcell overlap differs -- kept unenhanced, [E]).
                t_x, _z_eff, rate_mult = ibm_purcell_transmission(
                    dx_w, pp, Tj, gam, F_eff_op, w_meV=w, kappa_meV=kappa,
                    delta_c_meV=dx if held else None)
            else:
                t_x = float(ibm_transmission(dx_w, pp, Tj, gam, w_meV=w,
                                             kappa_meV=kappa,
                                             delta_c_meV=dx if held else None))
            t_xx = float(ibm_transmission(dx_w - d.dot.delta_xx, pp, Tj,
                                          d.dot.r_xx * gam, w_meV=w,
                                          kappa_meV=kappa,
                                          delta_c_meV=(dx - d.dot.delta_xx)
                                          if held else None))
            spec = SpectralResult(eps=t_xx / t_x, t_x=t_x, t_xx=t_xx)
        elif held:
            spec = epsilon2(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=w,
                            kappa=kappa, dx_w=dx_w, dx_c=dx)
        else:
            spec = epsilon(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=w, kappa=kappa, dx=dx)
        if d.ret.mode == "confinement":
            derived = _confinement_params(d.ret, Tj)
            params = {k: derived[k] for k in ("a_esc", "E_a", "b_p", "E_b")}
            params["b0"], params["beta"] = d.ret.b0, d.ret.beta
            # ret.overrides is an explicit dict (distinct from ret.system, which
            # only ever carries the inline stack): its presence, not its
            # truthiness, decides whether a coefficient is replaced, so an
            # explicitly supplied zero always wins -- unlike legacy's
            # 0="use proxy" shorthand on the RetentionBlock fields themselves.
            for key, value in d.ret.overrides.items():
                if key not in ("a_esc", "E_a", "b_p", "E_b", "b0", "beta"):
                    raise ValueError(f"unsupported ret.overrides key {key!r}")
                params[key] = float(value)
            # item 2: the confinement-derived transition energy feeds
            # transport's sub-turn-on loading suppression below; proxy mode
            # has no confinement level structure, so E_X_eV stays None there
            # (transport.dot_loading's documented legacy-numerics default).
            E_X_eV = derived["E_X_eV"]
        else:
            params = rp
            E_X_eV = None
        # R2: Purcell speeds the radiative channel by rate_mult (photon-
        # weighted; Lorentzian path uses F_eff directly since there is no
        # sideband split), so the escape-to-radiative ratios divide by it.
        # rm == 1.0 whenever `wire` is False (rate_mult/F_eff_op are then
        # left at their 1.0 initializers) -- also reused below by the CW
        # gamma_X_ns "radiative enhancement" (docs/rt_edge_contract.md).
        rm = rate_mult if d.dot.lineshape == "ibm" else F_eff_op
        if wire:
            S = float(retention(Tj, params["a_esc"] / rm, params["E_a"],
                                params["b_p"] / rm, params["E_b"]))
        else:
            S = float(retention(Tj, params["a_esc"], params["E_a"], params["b_p"], params["E_b"]))
        # PL mode: optical excitation, no injection-current background channel.
        injection = None
        if diode is not None:
            opts = _transport_options(d.drive, w if w is not None else 1.0)
            injection = transport.evaluate_injection(
                diode=diode, I_uA=d.drive.I_uA, T=Tj,
                n_dot_cm2=d.drive.n_dot_cm2 or d.aperture.density_cm2,
                aperture_um2=np.pi * (d.aperture.diameter_um / 2) ** 2,
                S_dot=S, E_X_eV=E_X_eV, **opts)
            inj_bg = injection.b_e
        else:
            inj_bg = float(b_injection(chan, d.drive.I_uA, Tj)) if d.drive.mode != "PL" else 0.0
        G = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
        B = params["b0"] + params["beta"] * (1.0 - S) + (G * S * inj_bg if diode else inj_bg)
        rho = G * S / (G * S + B)
        # T1 mechanism library (opt-in): resolve the card's mechanism into the
        # DriveInterface at THIS junction temperature (SET pricing is T-honest)
        # and use its (mu, F_p, eta). mechanism == "" keeps the legacy fields.
        mu_use, Fp_use, eta_use = d.drive.mu, d.drive.F_p, d.drive.eta_capture
        if injection is not None:
            mu_use, eta_use = injection.mu or 0.0, injection.eta_capture
        if d.drive.mechanism:
            iface = mech_from_card(d.drive.mechanism, d.drive.mech_params,
                                   T_K=Tj, I_uA=d.drive.I_uA,
                                   mu=d.drive.mu, eta=d.drive.eta_capture)
            mu_use, Fp_use, eta_use = iface.mu, iface.F_p, iface.eta_capture
        if mu_use > 0:
            if Fp_use != 1.0:
                # F8b thinning: effective pump Fano factor AT the dot, after
                # dot-capture; F8: moment-matched (mu, Fano) cap-2 loading.
                F_eff = f8b_thin_fano(eta_use, Fp_use)
                if d.drive.mechanism:
                    # mechanism path: out-of-domain (mu, F) at this T is a
                    # gap, not a crash (a mechanism may drift out of the
                    # cap-2 domain as Tj moves along the curve)
                    try:
                        g2_dot = float(f8_g2(mu_use, F_eff, spec.eps))
                    except ValueError:
                        g2_dot = np.nan
                else:
                    # legacy F_p path: keep v1.1 semantics EXACTLY -- domain
                    # violations raise; the envelope legs (design_meta) do
                    # their own gap handling and verify4 asserts the raise
                    g2_dot = float(f8_g2(mu_use, F_eff, spec.eps))
            else:
                # F_p == 1.0 (default): keep the original f1b_g2 path EXACTLY
                # -- f8_g2(mu, 1, eps) is not bit-identical to f1b_g2(mu, eps)
                # at finite mu (different, equally valid cap-2 conventions;
                # see loading.py module docstring), so every pre-v1.1 result
                # stays reproducible unless F_p is explicitly set != 1.0.
                g2_dot = float(f1b_g2(mu_use, spec.eps))
        else:
            g2_dot = spec.eps
        # Continuous aperture composition [A] (docs/rt_edge_contract.md
        # "Aperture assumptions"): composed AFTER loading/capture/filter are
        # already folded into g2_dot (spec.eps upstream, mu/F_p above), and
        # BEFORE the background law below -- every loss is applied exactly
        # once. Legacy (compose=False) leaves g2_dot untouched; the rounded
        # aperture_g2_penalty scalar (below, at the operating point only)
        # stays informational-only in that case, unchanged from before.
        w_ap = gam if (d.filter.auto_w or not d.filter.enabled) else d.filter.w
        if d.aperture.compose:
            _, lam_row = _aperture_lambda(
                d.aperture.density_cm2, np.pi * (d.aperture.diameter_um / 2) ** 2,
                w_ap, d.aperture.sigma_inh, d.aperture.comp_brightness)
            g2_dot = float(_compose_aperture_g2(g2_dot, lam_row))
        else:
            lam_row = float("nan")

        # CW diagnostics (opt-in, drive.cw=True; requires EL-transport so the
        # pump rate comes from transport.loading.r_dot, never inferred from
        # the dimensionless pulsed mu). gamma_X_ns carries the same Purcell
        # "radiative enhancement" (rm) as the pulsed retention S above; k_X/
        # k_XX are the matching ABSOLUTE non-radiative rates for that same
        # enhanced gamma_X_ns (integrator.retention-consistent, see
        # cw_g2.escape_rates_from_retention). rho_cw is derived from the
        # TOTAL filtered X+XX signal, not 1/(1+b_e): b_e is a per-X ratio, so
        # naively applying it to the X+XX total would double-discount the XX
        # contribution.
        g2_cw0 = g2_cw0_raw = float("nan")
        cw_r_ns = cw_gamma_X_ns = cw_rho = float("nan")
        if d.drive.cw and injection is not None:
            gamma_X_ns = rm / d.ret.tau_rad_ns
            k_X, k_XX = cw_g2.escape_rates_from_retention(
                gamma_X_ns, params["a_esc"], params["E_a"], params["b_p"], params["E_b"], Tj)
            r_ns = injection.loading.r_dot / 1e9
            t_X, t_XX = spec.t_x, spec.eps * spec.t_x
            if r_ns > 0 and gamma_X_ns > 0 and t_X > 0:
                I_X, I_XX = cw_g2.photon_rates(r_ns, gamma_X_ns, 2.0 * gamma_X_ns,
                                               k_X, k_XX, d.drive.cw_pump_ratio)
                sig = t_X * I_X + t_XX * I_XX
                # item 3 fix (council review 2026-09-05): the old
                # `inj_bg * t_X * I_X` chained the PULSED ratio inj_bg
                # (denominator = transport's own min(r_dot,1/tau_rad)*S,
                # Background.rate_x) onto the CW rate equation's OWN X rate
                # I_X -- consistent only away from saturation, 13x too large
                # at the card operating point where the two X rates diverge.
                # Carry transport's absolute in-window background rate
                # (photons/s, the SAME physical operating point regardless
                # of pulsed/CW framing) into the CW calculation and divide
                # by the SAME I_X this rate-equation solution produces --
                # abs_bg_in_window / I_X * (t_X * I_X) -- so I_X cancels and
                # the ratio is formed once, consistently, instead of mixing
                # two different models' X rates across the multiplication.
                abs_bg_in_window_ns = injection.background.rate_bg_window / 1e9
                bg = abs_bg_in_window_ns * t_X
                rho_cw = sig / (sig + bg) if (sig + bg) > 0 else float("nan")
                report = cw_g2.cw_report(
                    r_ns, gamma_X_ns, 2.0 * gamma_X_ns, k_X, k_XX, t_X, t_XX, rho_cw,
                    d.drive.cw_irf_fwhm_ps, tau_max_ns=d.drive.cw_tau_max_ns,
                    pump_ratio=d.drive.cw_pump_ratio, irf_shape=d.drive.cw_irf_shape)
                if d.aperture.compose and np.isfinite(lam_row):
                    tau = report["curves"]["tau"]
                    g_mix = _compose_aperture_g2(report["curves"]["g2_dot"], lam_row)
                    g_meas = cw_g2.g2_with_background(g_mix, rho_cw)
                    g_raw = (cw_g2.convolve_irf(tau, g_meas, d.drive.cw_irf_fwhm_ps,
                                                d.drive.cw_irf_shape)
                             if d.drive.cw_irf_fwhm_ps > 0 else g_meas)
                    g2_cw0 = float(np.interp(0.0, tau, g_meas))
                    g2_cw0_raw = float(np.interp(0.0, tau, g_raw))
                else:
                    g2_cw0 = report["g2_meas0"]
                    g2_cw0_raw = report["g2_raw0"]
                cw_r_ns, cw_gamma_X_ns, cw_rho = r_ns, gamma_X_ns, rho_cw
        return dict(Tj=Tj, gam=gam, eps=spec.eps, rho=rho,
                    g2=g2_from(g2_dot, rho), t_x=spec.t_x, runaway=False,
                    mu=mu_use, eta_capture=eta_use, b_e=inj_bg,
                    injection=injection, S=S, g2_cw0=g2_cw0, g2_cw0_raw=g2_cw0_raw,
                    cw_r_ns=cw_r_ns, cw_gamma_X_ns=cw_gamma_X_ns, cw_rho=cw_rho)

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
    if d.drive.cw:
        # T-aligned CW diagnostics (docs/rt_edge_contract.md); g2_op/g2 above
        # remain the pulsed-intrinsic headline metric (Lemma 1 note).
        curves["g2_cw0"] = np.array([r["g2_cw0"] for r in rows])
        curves["g2_cw0_raw"] = np.array([r["g2_cw0_raw"] for r in rows])

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

    _, P1, P2 = loading_probs(op["mu"])
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
    lam_op = Nw * d.aperture.comp_brightness if d.aperture.compose and np.isfinite(Nw) else float("nan")

    gain_factor = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
    beta_factor = d.cavity.beta_sin if sin_mode else 1.0
    brightness = float((P1 + P2)) * (op["t_x"] if np.isfinite(op["t_x"]) else 0.0) \
        * gain_factor * beta_factor
    if diode is not None:
        # `mu` is captured-dot loading already: multiply retention and the
        # spectral window once, never eta_inj/eta_capture again.
        brightness *= op["S"]

    # emission.type="edge" (Lemma 1: brightness only -- eta_total already
    # contains beta and the chosen front-facet fraction exactly once; T_facet/
    # beta/cavity.beta_sin are never applied a second time here).
    edge, edge_lambda_nm, edge_err = None, float("nan"), None
    if d.emission.type == "edge" and np.isfinite(op["Tj"]):
        try:
            edge, edge_lambda_nm = _resolve_edge(d.ret, d.emission, op["Tj"])
            brightness *= edge.eta_total
        except ValueError as e:
            edge_err = str(e)
            brightness = float("nan")
    elif d.emission.type == "edge":
        edge_err = "no finite operating temperature for edge emission"
        brightness = float("nan")
    scalars = {
        "T_j_op": op["Tj"], "dT_J": op["Tj"] - d.thermal.T_hs,
        "runaway": bool(op["runaway"]),
        "gamma_op": op["gam"], "eps_op": op["eps"], "rho_op": op["rho"],
        "g2_op": op["g2"], "t_x_op": op["t_x"],
        "brightness_per_pulse": brightness,
        "T_c": Tc,
        "F_eff": float(purcell_eff(d.cavity.F_P, d.cavity.kappa, op["gam"]))
        if (d.cavity.enabled and not sin_mode and np.isfinite(op["gam"])) else np.nan,
        "N_w": Nw, "aperture_g2_penalty": ap_pen,
        "tag_chain": "[A]",  # unmeasured inputs are always in the chain today
        "linewidth_source": d.dot.linewidth,
        "retention_source": d.ret.mode,
        "drive_source": d.drive.mode,
        "mu_resolved": op["mu"],
        "eta_capture_resolved": op["eta_capture"],
        "b_e_resolved": op["b_e"],
        "V_j_op": op["injection"].V_j if op["injection"] is not None else d.drive.V,
        "V_j": op["injection"].V_j if op["injection"] is not None else d.drive.V,
        "V_applied": op["injection"].V_applied if op["injection"] is not None else d.drive.V,
        "V_bi": op["injection"].V_bi if op["injection"] is not None else np.nan,
        "P_junction_W": (op["injection"].P_junction_W
                         if op["injection"] is not None else d.drive.I_uA * 1e-6 * d.drive.V),
        "eta_inj": (op["injection"].leakage.eta_inj
                    if op["injection"] is not None else np.nan),
        "T_j_transport": op["injection"].T if op["injection"] is not None else np.nan,
        "S_resolved": op["S"],
        "loading.r_dot": (op["injection"].loading.r_dot
                          if op["injection"] is not None else np.nan),
        # Pre-collection photon rates staged for integration-b (transport.py
        # native units, photons/s -- not rescaled to ns^-1 here so that these
        # stay exact pass-throughs of the verified Background dataclass).
        "background.rate_bg_window": (op["injection"].background.rate_bg_window
                                      if op["injection"] is not None else np.nan),
        "background.rate_x": (op["injection"].background.rate_x
                              if op["injection"] is not None else np.nan),
        # emission.type="edge" (Lemma 1: reported, never re-multiplied into
        # g2/eps/rho -- see the brightness composition above).
        "emission_type": d.emission.type,
        "edge_lambda_nm": edge_lambda_nm,
        "edge_beta": edge.beta if edge is not None else np.nan,
        "edge_T_facet": edge.T_facet if edge is not None else np.nan,
        "edge_eta_prop": edge.eta_prop if edge is not None else np.nan,
        "edge_eta_NA": edge.eta_NA if edge is not None else np.nan,
        "edge_eta_total": edge.eta_total if edge is not None else np.nan,
        "edge_n_eff": edge.n_eff if edge is not None else np.nan,
        "edge_n_g": edge.n_g if edge is not None else np.nan,
        "edge_Gamma_dot": edge.Gamma_dot if edge is not None else np.nan,
        "edge_A_mode_um2": edge.A_mode_um2 if edge is not None else np.nan,
        # continuous aperture composition (docs/rt_edge_contract.md); N_w
        # above is already the unrounded competitor count.
        "aperture_compose": d.aperture.compose,
        "aperture_lambda_op": lam_op,
        # material tracking (filter.track_material)
        "track_material": d.filter.track_material,
        # CW diagnostics (drive.cw)
        "g2_cw0": op["g2_cw0"], "g2_cw0_raw": op["g2_cw0_raw"],
        "cw_r_ns": op["cw_r_ns"], "cw_gamma_X_ns": op["cw_gamma_X_ns"],
        "cw_rho_op": op["cw_rho"],
        "invalid_reasons": ([] if np.isfinite(op["g2"])
                            else [op.get("invalid_reason", "invalid operating point")])
        + ([edge_err] if edge_err else [])
        + (["CW diagnostics requested but not evaluable at this operating point"]
           if d.drive.cw and not np.isfinite(op["g2_cw0"]) else [])
        # item 2: sub-turn-on warning, informational -- the sweep can see the
        # dot is being loaded almost entirely off the Boltzmann tail of the
        # carrier reservoirs (qV_j far below E_X) without this making the
        # operating point itself invalid (g2_op stays finite).
        + (["sub_turn_on: f_qfl < 0.05 (qV_j well below the dot transition; "
            "loading exponentially suppressed, Sze & Ng ch. 12)"]
           if (op["injection"] is not None and op["injection"].f_qfl < 0.05) else []),
        "provenance": {
            "linewidth": {"tag": "A" if d.dot.linewidth == "anchored" else "A",
                          "note": "anchored linewidth model" if d.dot.linewidth == "anchored"
                                  else "arsenide class proxy"},
            "retention": {"tag": "E" if d.ret.mode == "confinement" else "A",
                          "note": retention_note},
            "drive": {"tag": "E" if diode is not None else "A",
                      "note": "transport Diode injection result" if diode is not None
                              else "legacy drive fields"},
            "mu_resolved": {"tag": "E" if diode is not None else "A",
                            "note": "transport dot loading" if diode is not None else "drive.mu"},
            "b_e_resolved": {"tag": "E" if diode is not None else "A",
                              "note": "in-window background per collected X photon"},
            "emission": {"tag": "E" if edge is not None else "A",
                        "note": ("waveguide.edge_emission on the shared inline stack"
                                 if edge is not None else "no edge collection (legacy)")},
            "aperture": {"tag": "A",
                        "note": ("continuous Poisson competitor-bath composition"
                                 if d.aperture.compose
                                 else "legacy informational aperture_g2_penalty only")},
            "tracking": {"tag": "DR" if d.filter.track_material else "A",
                        "note": (f"materials.bandgap on stack layer {d.filter.track_material!r}"
                                 if d.filter.track_material else "legacy hard-coded GaAs Varshni")},
            "cw": {"tag": "A", "note": ("cw_g2 rate-equation model, transport-derived pump rate"
                                        if d.drive.cw else "not requested")},
        },
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
