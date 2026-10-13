"""Self-describing registries of the composable engine (design 2.1-2.6).

Axes (singular, matching the composed-spec keys): material, structure, drive,
setup, solver, observable. Each registered object carries ParamSpecs, so the
Studio renders forms generically from `registry_payload()`.

Registration is lazy-import friendly: later pieces add modules
`fsim_core.compose.{materials,structures,drives,setups,solvers,observables}`
that call `register(...)` at import; `ensure_loaded()` imports those that
exist and silently skips the ones that do not yet.
"""
from __future__ import annotations

import importlib
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar, Optional

from .params import ParamSpec
from .types import SetupChain

AXES = ("material", "structure", "drive", "setup", "solver", "observable")
_PLURAL = {"material": "materials", "structure": "structures", "drive": "drives",
           "setup": "setups", "solver": "solvers", "observable": "observables"}
_LAZY_MODULES = ("materials", "structures", "drives", "setups", "solvers", "observables")
PAYLOAD_SCHEMA = "fsim-compose-registry/1"


# --------------------------------------------------------------------------
# Entry base classes
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Entry(ABC):
    """Design 2.1. `axis` is a ClassVar set by the Material/Structure/Drive/
    Setup subclass (a dataclass field would break the default ordering)."""
    axis: ClassVar[str] = ""
    kind: str
    label: str
    summary: str
    params: tuple = ()               # tuple[ParamSpec, ...]
    refs: tuple = ()
    status: str = "ready"            # ready | needs_physics | tier
    provides: frozenset = frozenset()
    requires: frozenset = frozenset()
    viz: str = "fallback2d"
    eta_hint: tuple = (0.001, 0.01)

    def __post_init__(self):
        params = tuple(self.params)
        object.__setattr__(self, "params", params)
        for p in params:
            if not isinstance(p, ParamSpec):
                raise TypeError(f"{self.kind}: params must be ParamSpec, got {type(p).__name__}")
        names = [p.name for p in params]
        if len(set(names)) != len(names):
            raise ValueError(f"{self.kind}: duplicate param names in {names}")
        object.__setattr__(self, "provides", frozenset(self.provides))
        object.__setattr__(self, "requires", frozenset(self.requires))

    def param(self, name: str) -> ParamSpec:
        for p in self.params:
            if p.name == name:
                return p
        raise KeyError(f"{self.axis} {self.kind!r} has no param {name!r}; "
                       f"known: {[p.name for p in self.params]}")


@dataclass(frozen=True)
class Material(Entry):
    axis: ClassVar[str] = "material"

    @abstractmethod
    def build(self, vals: dict, T_K: float):
        """-> types.EmitterModel"""


@dataclass(frozen=True)
class Structure(Entry):
    axis: ClassVar[str] = "structure"

    @abstractmethod
    def response(self, vals: dict, em, T_K: float, NA: float):
        """-> types.StructureResponse"""


@dataclass(frozen=True)
class Drive(Entry):
    axis: ClassVar[str] = "drive"
    observables: frozenset = frozenset()

    @abstractmethod
    def plan(self, vals: dict, em, rt):
        """-> types.DrivePlan"""


@dataclass(frozen=True)
class Setup(Entry):
    axis: ClassVar[str] = "setup"

    @abstractmethod
    def chain(self, vals: dict) -> SetupChain:
        """Collection / detection chain from resolved parameter values."""


class Solver:
    """Design 2.5. Subclasses set `name` (the registry kind), `label`, `summary`."""
    axis: ClassVar[str] = "solver"
    name: str = ""
    label: str = ""
    summary: str = ""

    @property
    def kind(self) -> str:
        return self.name

    def supports(self, plan, em, rt) -> tuple:
        """-> (ok: bool, reason: str)"""
        raise NotImplementedError

    def run(self, plan, em, rt):
        """-> types.SolverResult"""
        raise NotImplementedError

    def eta_seconds(self, plan, n_points) -> float:
        raise NotImplementedError


@dataclass(frozen=True)
class Observable:
    """Design 2.6."""
    axis: ClassVar[str] = "observable"
    id: str
    label: str
    unit: str
    deps: tuple = ()
    conventions: tuple = ("intrinsic",)
    passthrough_of: Optional[str] = None
    tex: Optional[str] = None        # display symbol in KaTeX TeX (UI only)

    @property
    def kind(self) -> str:
        return self.id

    def compute(self, ctx, convention) -> float:
        raise NotImplementedError


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------
class Registry:
    def __init__(self):
        self._d = {a: {} for a in AXES}

    def register(self, entry):
        axis = getattr(entry, "axis", "")
        if axis not in self._d:
            raise ValueError(f"cannot register {entry!r}: unknown axis {axis!r}; axes: {list(AXES)}")
        kind = entry.kind
        if not kind or not isinstance(kind, str):
            raise ValueError(f"{axis} entry needs a non-empty string kind, got {kind!r}")
        if kind in self._d[axis]:
            raise ValueError(f"duplicate {axis} kind {kind!r}; registered: {sorted(self._d[axis])}")
        self._d[axis][kind] = entry
        return entry

    def unregister(self, axis: str, kind: str):
        self._d[axis].pop(kind, None)

    def get(self, axis: str, kind: str):
        if axis not in self._d:
            raise KeyError(f"unknown axis {axis!r}; axes: {list(AXES)}")
        try:
            return self._d[axis][kind]
        except KeyError:
            raise KeyError(f"unknown {axis} kind {kind!r}; registered: {sorted(self._d[axis])}") from None

    def kinds(self, axis: str) -> list:
        if axis not in self._d:
            raise KeyError(f"unknown axis {axis!r}; axes: {list(AXES)}")
        return sorted(self._d[axis])


REGISTRY = Registry()
_STATE = {"loaded": False, "loading": False}


def ensure_loaded(force: bool = False) -> None:
    """Import the optional registry modules that exist. A missing module is
    skipped; any other exception (including an ImportError raised INSIDE an
    existing module) propagates."""
    if (_STATE["loaded"] and not force) or _STATE["loading"]:
        return
    _STATE["loading"] = True
    try:
        for name in _LAZY_MODULES:
            full = f"{__package__}.{name}"
            try:
                importlib.import_module(full)
            except ModuleNotFoundError as e:
                if e.name != full:
                    raise
        _STATE["loaded"] = True
    finally:
        _STATE["loading"] = False


def _reg(registry):
    if registry is not None:
        return registry
    ensure_loaded()
    return REGISTRY


def register(entry, registry: Optional[Registry] = None):
    """Axis-keyed registration; a duplicate kind raises ValueError."""
    return (registry or REGISTRY).register(entry)


def unregister(axis: str, kind: str, registry: Optional[Registry] = None) -> None:
    (registry or REGISTRY).unregister(axis, kind)


def get(axis: str, kind: str, registry: Optional[Registry] = None):
    """Unknown kind -> KeyError whose message lists the registered kinds."""
    return _reg(registry).get(axis, kind)


def kinds(axis: str, registry: Optional[Registry] = None) -> list:
    return _reg(registry).kinds(axis)


def _param_dict(p: ParamSpec) -> dict:
    d = {"name": p.name, "default": p.default, "scale": p.scale, "help": p.help,
         "group": p.group, "visible_if": p.visible_if, "sweepable": p.sweepable,
         "lo": None, "hi": None, "choices": None}
    d.update(p.to_meta())            # unit, tag, source (+ lo, hi, choices when set)
    return d


def _entry_dict(e) -> dict:
    axis = e.axis
    if axis == "solver":
        return {"axis": axis, "kind": e.kind, "label": e.label, "summary": e.summary}
    if axis == "observable":
        d = {"axis": axis, "kind": e.kind, "id": e.id, "label": e.label, "unit": e.unit,
             "deps": list(e.deps), "conventions": list(e.conventions),
             "passthrough_of": e.passthrough_of}
        if e.tex is not None:
            d["tex"] = e.tex
        return d
    d = {"axis": axis, "kind": e.kind, "label": e.label, "summary": e.summary,
         "params": [_param_dict(p) for p in e.params], "refs": list(e.refs),
         "status": e.status, "provides": sorted(e.provides), "requires": sorted(e.requires),
         "viz": e.viz, "eta_hint": list(e.eta_hint)}
    if axis == "drive":
        d["observables"] = sorted(e.observables)
    return d


def registry_payload(registry: Optional[Registry] = None) -> dict:
    """JSON-serialisable description of everything registered (design 2.1).
    `tiers` and `compat` are filled by later pieces (tier passthroughs, explicit
    incompatibility reasons); empty lists here."""
    reg = _reg(registry)
    out = {"schema": PAYLOAD_SCHEMA}
    for axis in AXES:
        out[_PLURAL[axis]] = [_entry_dict(reg.get(axis, k)) for k in reg.kinds(axis)]
    out["tiers"] = []
    out["compat"] = []
    return out


# --------------------------------------------------------------------------
# Built-in setup entry
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class _LensFibreDet(Setup):
    def chain(self, vals: dict) -> SetupChain:
        g = lambda n: vals.get(n, self.param(n).default)  # noqa: E731
        return SetupChain(NA=float(g("NA")), eta_pol=float(g("eta_pol")),
                          eta_optics=float(g("eta_optics")), eta_fibre=float(g("eta_fibre")),
                          eta_det=float(g("eta_det")), laser_leak=float(g("laser_leak")),
                          dark_cps=float(g("dark_cps")), zpl_filter=bool(g("zpl_filter")))


_SRC = "assumed lab-setup placeholder; replace with the measured chain [A]"
LENS_FIBRE_DET = _LensFibreDet(
    kind="lens_fibre_det", label="Lens, fibre, detector",
    summary="First-lens collection into single-mode fibre and a single-photon detector; "
            "every efficiency is an entered number.",
    params=(
        ParamSpec("NA", "1", 0.65, tex=r"\mathrm{NA}", lo=0.05, hi=1.0, tag="A", source=_SRC,
                  help=r"Numerical aperture $\mathrm{NA}$ of the first collection lens."),
        ParamSpec("eta_pol", "1", 1.0, tex=r"\eta_{\mathrm{pol}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Polarisation-selection transmission $\eta_{\mathrm{pol}}$ (0.5 for a single linear analyser on unpolarised light)."),
        ParamSpec("eta_optics", "1", 0.5, tex=r"\eta_{\mathrm{optics}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Transmission $\eta_{\mathrm{optics}}$ of lenses, filters and mirrors between the first lens and the fibre."),
        ParamSpec("eta_fibre", "1", 1.0, tex=r"\eta_{\mathrm{fibre}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Fibre-coupling efficiency $\eta_{\mathrm{fibre}}$."),
        ParamSpec("eta_det", "1", 0.5, tex=r"\eta_{\mathrm{det}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Detector quantum efficiency $\eta_{\mathrm{det}}$."),
        ParamSpec("laser_leak", "1", 0.0, tex=r"f_{\mathrm{leak}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Residual laser light $f_{\mathrm{leak}}$ as a fraction of the signal rate; adds to the raw $g^{(2)}(0)$."),
        ParamSpec("dark_cps", "1/s", 0.0, tex=r"r_{\mathrm{dark}}", lo=0.0, hi=1e5, tag="A", source=_SRC,
                  help=r"Detector dark counts per second $r_{\mathrm{dark}}$; adds to the raw $g^{(2)}(0)$. Pulsed drives count only the darks inside "
                       r"$T_{\mathrm{coinc}}$ around each peak."),
        ParamSpec("coinc_window_ns", "ns", 2.0, tex=r"T_{\mathrm{coinc}}", lo=0.01, hi=1000.0, tag="A", source=_SRC,
                  help=r"Pulsed drives: histogram integration window $T_{\mathrm{coinc}}$ around each coincidence peak; the dark background per "
                       r"pulse and detector is $r_{\mathrm{dark}}\,T_{\mathrm{coinc}}$ (window capped at the pulse period). "
                       r"Ignored by continuous drives."),
        ParamSpec("eps_int", "1", 0.0, tex=r"\varepsilon_{\mathrm{int}}", lo=0.0, hi=1.0, tag="A", source=_SRC,
                  help=r"Interferometer imperfection $\varepsilon_{\mathrm{int}}$ (1 - fringe contrast of a perfect source) in the raw HOM visibility: "
                       r"$V_{\mathrm{raw}} = \left[(1 - g^{(2)})\,I - g^{(2)}\right](1 - \varepsilon_{\mathrm{int}})$."),
        ParamSpec("zpl_filter", "", False, tex=r"\mathrm{ZPL\ filter}", tag="A", source=_SRC,
                  help=r"Spectrally filter to the zero-phonon line (brightness $\times Z$, HOM unity)."),
        ParamSpec("irf_ps", "ps", 0.0, tex=r"\tau_{\mathrm{IRF}}", lo=0.0, hi=2000.0, tag="A", source=_SRC,
                  help=r"Detector timing response $\tau_{\mathrm{IRF}}$ (Gaussian FWHM); 0 = off. Feeds the IRF-convolved cw dip ($g^{(2)}_{\mathrm{IRF}}$) and "
                       r"the IRF leak of the binned GHz histogram ($g^{(2)}_{\mathrm{bin}}$). A drive-level $\tau_{\mathrm{IRF}}$ (electrical drives) "
                       r"takes precedence when non-zero."),
        ParamSpec("filter_width_meV", "meV", 0.0, tex=r"w_{\mathrm{filt}}", lo=0.0, hi=500.0, tag="A", source=_SRC, group="advanced",
                  help=r"Top-hat spectral window width $w_{\mathrm{filt}}$ for the X / XX leakage ($\varepsilon_{XX}$); 0 = no window "
                       r"(both lines fully collected)."),
        ParamSpec("filter_offset_meV", "meV", 0.0, tex=r"\delta_{\mathrm{filt}}", lo=-500.0, hi=500.0, tag="A", source=_SRC, group="advanced",
                  help=r"Offset $\delta_{\mathrm{filt}}$ of the X line from the window centre (negative = window blue of X)."),
        ParamSpec("xx_width_ratio", "1", 1.5, tex=r"w_{XX}/w_{X}", lo=0.1, hi=10.0, tag="E", group="advanced",
                  source="XX/X linewidth ratio ~1.5 (Tamariz et al. 2020: 80 vs 55 meV at 300 K) [E]",
                  help=r"Width of the XX line as a multiple of the X linewidth $\hbar(\Gamma + k + \gamma^{*})$."),
        ParamSpec("xx_x_ratio", "1", 0.0, tex=r"I_{XX}/I_{X}", lo=0.0, hi=100.0, tag="A", source=_SRC, group="advanced",
                  help=r"Unfiltered X-XX intensity ratio $I_{XX}/I_X$ at the operating point (from the power exponents); "
                       r"0 = no XX light. Sets $\rho$ and $g^{(2)} = 1 - \rho^2$ in the $\varepsilon_{XX}$ observable."),
    ),
)
REGISTRY.register(LENS_FIBRE_DET)
