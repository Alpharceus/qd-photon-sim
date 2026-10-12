"""Material registry of the composable engine (design 2.2, P2).

Each material system is a `Material` entry registered on import. All numbers
live in `data/materials.yaml` as provenance-tagged leaves (unit, tag in
V/DR/E/A, source); this module turns them into ParamSpecs and a
`build(vals, T_K) -> EmitterModel`.

Content source: .workers/devlib/p1-physics-brief.md section 4. Repo numbers
(InP phonon constants, anchored linewidth model, DOT_PRESETS binding-energy
band, E_LO) are imported at load time, not copied.

Conventions: rates 1/ns, energies meV, hbar = 6.582e-4 meV ns; gamma* is the
repo convention (coherence decays at (Gamma + gamma*)/2, so a measured excess
Lorentzian FWHM of W gives gamma* = W / hbar). The gamma*(T) model is linear,
gamma*(4 K) + a_ac (T - 4 K), clipped at 0; the LO-phonon term is NOT included.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import ClassVar

import yaml

from .. import cavity as _cavity
from .. import linewidth as _linewidth
from .. import nitride_materials as _nitride
from .. import presets as _presets
from .. import qd_gf
from .params import Leaf, ParamSpec, central
from .registry import REGISTRY, Material, register
from .types import EmitterModel

HBAR_MEV_NS = 6.582119569e-4          # meV ns (repo value, lindblad.HBAR_MEV_NS)
HC_MEV_NM = 1239841.98                # meV nm
T_REF_K = 4.0
DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "materials.yaml")

_TAGS = ("V", "DR", "E", "A")
_CONVERT = {
    "ueV_to_rate": lambda x: x * 1.0e-3 / HBAR_MEV_NS,
    "meV_to_rate": lambda x: x / HBAR_MEV_NS,
}

_ADV = "advanced"
# name -> (ParamSpec unit, group, help). Order = display order.
_PARAMS = {
    "lambda_nm": ("nm", "main", "Zero-phonon exciton emission wavelength; E_X = 1239841.98 / lambda_nm meV."),
    "tau_bare_ns": ("ns", "main", "Bare (homogeneous-medium) exciton radiative lifetime; Gamma0 = 1 / tau_bare."),
    "fss_ueV": ("ueV", "main", "Exciton fine-structure splitting (a beat, not dephasing)."),
    "E_b_meV": ("meV", "main", "Biexciton binding energy, positive = binding (E_XX = 2 E_X - E_b)."),
    "gamma_xx_ratio": ("1", "main", "Radiative rate ratio Gamma_XX / Gamma_X of the biexciton to the exciton."),
    "gamma_star_4K_ns": ("1/ns", "main", "Pure dephasing rate gamma* at 4 K, repo convention (coherence decays at "
                                         "(Gamma + gamma*)/2; gamma* = excess FWHM / hbar). Linear T model; LO term not included."),
    "a_ac_ns_per_K": ("1/ns/K", "main", "Linear acoustic T-dephasing slope: gamma*(T) = gamma*(4 K) + a_ac (T - 4 K), clipped >= 0. "
                                        "LO term not included."),
    "alpha_ps2": ("ps^2", _ADV, "Direct phonon coupling alpha. 0 = compute from D_c, D_v, rho, c_LA, l_xy, l_z "
                                "(qd_gf.coupling_alpha_ps2); > 0 pins alpha."),
    "D_c_eV": ("eV", _ADV, "Conduction-band deformation potential (only D_c - D_v enters alpha)."),
    "D_v_eV": ("eV", _ADV, "Valence-band deformation potential."),
    "rho_kg_m3": ("kg/m^3", _ADV, "Matrix mass density (phonons are bulk LA of the matrix)."),
    "c_LA_m_s": ("m/s", _ADV, "Matrix LA sound velocity; cutoff E_c = hbar sqrt(2) c_s / l (qd_gf convention)."),
    "l_xy_nm": ("nm", _ADV, "In-plane Gaussian confinement length of the exciton wavefunction (not the island size)."),
    "l_z_nm": ("nm", _ADV, "Growth-axis Gaussian confinement length of the exciton wavefunction."),
    "n_host": ("1", "main", "Host refractive index at the emission wavelength."),
    "E_LO_meV": ("meV", "main", "LO phonon energy of the host (reference; the LO term is not used by build())."),
    "k_nr_ns": ("1/ns", "main", "Photonless escape rate (carrier loss without photon emission); 0 unless stated [A]."),
    "sigma_sd_ueV": ("ueV", _ADV, "Spectral-diffusion width (Gaussian sigma)."),
    "tau_sd_ns": ("ns", _ADV, "Spectral-diffusion correlation time."),
}
_ENUMS = (
    ParamSpec("charge_state", "", "X", choices=("X", "X+", "X-"), tag="A", sweepable=False,
              source="user choice of the emitting charge state [A]",
              help="Charge state of the emitter. X adds nothing; X+ / X- add the charged-exciton transition."),
    ParamSpec("dipole", "", "in_plane", choices=("in_plane", "axial"), tag="A", sweepable=False,
              source="user choice of the transition dipole orientation [A]",
              help="Transition dipole orientation: in_plane (planar dot, bright in-plane) or axial (along growth axis, e.g. nanowire)."),
)

# eta_blink is a uniform setup-like factor (not a per-material literature number): it lives here, not in the yaml.
_BLINK = ParamSpec(
    "eta_blink", "1", 1.0, lo=0.5, hi=1.0, tag="A", group="main",
    source="assumed: no blinking / full on-fraction [A]",
    help="Blinking duty factor (on-fraction of the time the dot is bright). Multiplies eta_rad and the collected "
         "brightness (reported as eta_rad_eff); it does NOT change the lifetime. Examples: Uppu et al., Sci. Adv. "
         "6, eabc8268 (2020) eta_rad 0.98 [V-from cards]; Wang et al., Phys. Rev. Lett. 122, 113602 (2019) 35 % "
         "blinking loss [V-from cards] (on-fraction 0.65).")

# EmitterModel field -> contributing parameter names (for EmitterModel.inputs)
_INPUT_FIELDS = {
    "E_X_meV": "lambda_nm", "lambda_nm": "lambda_nm", "n_host": "n_host", "gamma0_ns": "tau_bare_ns",
    "k_nr_ns": "k_nr_ns", "delta_xx_meV": "E_b_meV", "gamma_xx_ratio": "gamma_xx_ratio",
    "fss_ueV": "fss_ueV", "dipole": "dipole", "sigma_sd_ueV": "sigma_sd_ueV", "tau_sd_ns": "tau_sd_ns",
}


# --------------------------------------------------------------------------
# data loading
# --------------------------------------------------------------------------
def _repo_roots():
    return {"PhononParams": qd_gf.PhononParams(), "LinewidthParams": _linewidth.LinewidthParams(),
            "DOT_PRESETS": _presets.DOT_PRESETS, "E_LO_MEV": _linewidth.E_LO_MEV}


def _import_number(ref: str) -> float:
    root, *path = ref.split(".")
    cur = _repo_roots()[root]
    for key in path:
        cur = cur[key] if isinstance(cur, dict) else getattr(cur, key)
    return float(cur)


def load_data(path: str = DATA_PATH) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def leaf_numbers(leaf: dict) -> tuple:
    """-> (kind, value-or-range) of a yaml leaf with imports resolved and
    `convert` applied. kind is 'value' or 'range'."""
    conv = _CONVERT[leaf["convert"]] if "convert" in leaf else (lambda x: x)
    if "value" in leaf:
        return "value", conv(leaf["value"])
    if "value_from" in leaf:
        return "value", conv(_import_number(leaf["value_from"]))
    if "range" in leaf:
        lo, hi = leaf["range"]
        return "range", (conv(lo), conv(hi))
    if "range_from" in leaf:
        lo, hi = (_import_number(r) for r in leaf["range_from"])
        return "range", (conv(lo), conv(hi))
    raise ValueError(f"leaf has no value/range: {leaf}")


def _spec_from_leaf(name: str, leaf: dict) -> ParamSpec:
    unit, group, base_help = _PARAMS[name]
    kind, num = leaf_numbers(leaf)
    scale = leaf.get("scale", "lin")
    help_ = base_help
    if kind == "range":
        lo, hi = num
        default = central(Leaf(path=name, kind="range", unit=unit, tag=leaf["tag"], range=(lo, hi), scale=scale))
        help_ += " Default = midpoint of the class range, not a prediction."
    else:
        default = num
        lo, hi = leaf.get("band", (None, None)) if "band" in leaf else (None, None)
    help_ += f" Source: {leaf['source']}"
    return ParamSpec(name, unit, float(default), lo=None if lo is None else float(lo),
                     hi=None if hi is None else float(hi), scale=scale, tag=leaf["tag"],
                     source=leaf["source"], help=help_, group=group)


# --------------------------------------------------------------------------
# the Material entry
# --------------------------------------------------------------------------
def _dEdT_meV_K(varshni: str, T_K: float) -> float:
    """Slope of the Varshni gap E(T) = E0 - a T^2/(T + b) at T [meV/K]. GaAs and
    InP from fsim_core.cavity.VARSHNI; GaN from fsim_core.nitride_materials."""
    if varshni == "GaN":
        g = _nitride.binary("GaN")
        a, b = g.alpha_eVK, g.beta_K
    else:
        a, b = _cavity.VARSHNI[varshni]
    return -1.0e3 * a * T_K * (T_K + 2.0 * b) / (T_K + b) ** 2


@dataclass(frozen=True)
class DataMaterial(Material):
    axis: ClassVar[str] = "material"
    varshni: str = "GaAs"

    def build(self, vals: dict, T_K: float) -> EmitterModel:
        vals = vals or {}

        def g(n):
            return vals.get(n, self.param(n).default)

        T = float(T_K)
        lam = float(g("lambda_nm"))
        if not lam > 0:
            raise ValueError(f"{self.kind}: lambda_nm must be > 0, got {lam}")
        tau = float(g("tau_bare_ns"))
        if not tau > 0:
            raise ValueError(f"{self.kind}: tau_bare_ns must be > 0, got {tau}")
        alpha = float(g("alpha_ps2"))
        ph = qd_gf.PhononParams(D_e_eV=float(g("D_c_eV")), D_h_eV=float(g("D_v_eV")),
                                rho_m_kg_m3=float(g("rho_kg_m3")), c_s_m_s=float(g("c_LA_m_s")),
                                l_xy_nm=float(g("l_xy_nm")), l_z_nm=float(g("l_z_nm")),
                                alpha_ps2=alpha if alpha > 0 else None)
        gs = max(0.0, float(g("gamma_star_4K_ns")) + float(g("a_ac_ns_per_K")) * (T - T_REF_K))
        trans = {"X", "XX"}
        if g("charge_state") != "X":
            trans.add(str(g("charge_state")))
        inputs = {f: f"material.{p}" for f, p in _INPUT_FIELDS.items()}
        inputs["gamma_star_ns"] = "material.gamma_star_4K_ns+a_ac_ns_per_K"
        inputs["phonon"] = ("material.D_c_eV+D_v_eV+rho_kg_m3+c_LA_m_s+l_xy_nm+l_z_nm+alpha_ps2")
        # delta_xx_meV in the resonant_source.FockSource convention: NEGATIVE = bound (E_b > 0).
        return EmitterModel(
            E_X_meV=HC_MEV_NM / lam, lambda_nm=lam, n_host=float(g("n_host")),
            gamma0_ns=1.0 / tau, k_nr_ns=float(g("k_nr_ns")), gamma_star_ns=gs, phonon=ph,
            delta_xx_meV=-float(g("E_b_meV")), gamma_xx_ratio=float(g("gamma_xx_ratio")),
            fss_ueV=float(g("fss_ueV")), transitions=frozenset(trans), dipole=str(g("dipole")),
            dEdT_meV_K=_dEdT_meV_K(self.varshni, T), sigma_sd_ueV=float(g("sigma_sd_ueV")),
            tau_sd_ns=float(g("tau_sd_ns")), inputs=inputs)


def build_entries(path: str = DATA_PATH) -> list:
    """Fresh DataMaterial entries from the yaml (for tests that want a private Registry)."""
    data = load_data(path)
    out = []
    for kind, m in data["materials"].items():
        specs = [_spec_from_leaf(n, m["params"][n]) for n in _PARAMS]
        missing = set(m["params"]) - set(_PARAMS)
        if missing:
            raise ValueError(f"{kind}: unknown yaml params {sorted(missing)}")
        out.append(DataMaterial(kind=kind, label=m["label"], summary=m["summary"],
                                params=tuple(specs) + (_BLINK,) + _ENUMS, refs=tuple(m.get("refs", ())),
                                status="ready", varshni=m["varshni"],
                                provides=frozenset({"X", "XX", "X+", "X-"})))
    return out


MATERIALS = build_entries()
for _e in MATERIALS:
    if _e.kind in REGISTRY.kinds("material"):
        REGISTRY.unregister("material", _e.kind)   # re-import / reload friendly
    register(_e)
