"""Access-weighted sidewall nonradiative loss for the opt-in nitride
nanowire tier (piece 4; depends on piece 1, docs/nitride_nanowire_contract.md
and verify/data/nitride_nanowire_anchors.yaml).

This module prices ONE thing: a uniform cylindrical sidewall recombination
velocity turned into a per-nanosecond rate, split into two disjoint access
channels (a carrier sitting in the axial GaN reservoir, versus a carrier
already captured onto the occupied disc/dot).  It does not touch strain,
conducting radius, optical index, or barrier profile; those live in other
pieces.  It is not a fit: no parameter here is chosen by solving for a
measured ratio (see ``yield_ratio`` below and CLAUDE.md's provenance rule).

Provenance of every literal:
  - S_cm_s=1e3 [E]: Deshpande et al., Nat. Commun. 4, 1675 (2013), p.1,
    "the surface recombination velocity of GaN nanowires has been reported
    to be as low as 1e3 cm/s" -- a SECONDARY attribution to ref. 35 of that
    paper, not a value Deshpande et al. measured themselves.  See
    verify/data/nitride_nanowire_anchors.yaml:deshpande2013_thermal_and_pl,
    field surface_S_cm_s_secondary.
  - The uniform-cylinder rate k_side = 2S/R [DR]: standard surface/volume
    recombination arithmetic for a cylinder of radius R with a uniform bulk
    carrier density and a boundary condition set by S (Sze & Ng, Physics of
    Semiconductor Devices, 3rd ed., ch. 1).  Applying that BULK expression
    to a carrier sitting in a thin, LOCALIZED disc (occupied dot) or moving
    through the axial reservoir is an explicit MODEL TRANSFER [A], not
    itself verified by the Deshpande quotation -- the quotation only fixes
    S, not the geometry law.  A hard-wall envelope wavefunction that is
    exactly zero at the mathematical sidewall does not prove the true
    (diffusing, imperfectly confined) carrier never reaches a surface trap;
    reservoir_access and occupied_dot_access parametrize that missing
    diffusion/localization physics rather than assuming it away.
  - shell_multiplier=0.1 for AlGaN, 1.0 for none [A]: no passivation ratio
    is reported for this system; 0.1 is a round, explicitly assumed
    sensitivity factor, not a measured passivation improvement.
  - reservoir_access=1.0, occupied_dot_access=1.0 [A]: full, unweighted
    access is the conservative (maximal-loss) default in the absence of a
    measured localization length; both may be swept independently.
  - k_surface_XX = 2 * k_surface_X [A]: an independent-carrier cascade
    convention, matching the same doubling already used for other
    two-carrier (XX) rates elsewhere in the nitride tier (e.g.
    fsim_core/nitride_levels.py's k_XX_ns treatment) -- each of the two
    carriers in the biexciton is charged an independent copy of the
    single-carrier occupied-dot surface loss.

Reservoir vs occupied-dot channel semantics (contract-fixed, do not blur):
  - k_surface_reservoir_ns feeds piece 5's competition between reservoir
    surface loss, capture, and matrix recombination.  A reservoir surface
    death emits NO background photon.
  - k_surface_X_ns / k_surface_XX_ns feed piece 7's addition to the
    intrinsic k_nr and to escape AFTER capture.
  These act on disjoint carrier populations (uncaptured vs captured); this
  module never subtracts both from one supply, and never adjusts one
  channel to cancel the other.

Non-gating comparison, not a fit target: the measured ensemble PL ratio
I(300K)/I(10K) = 0.52 (Deshpande 2013, p.2; ledger key
deshpande2013_thermal_and_pl, field ensemble_PL_300K_over_10K) is an
ENSEMBLE, many-wire loss comparison.  It is not a single-dot internal
quantum efficiency of 0.52, and it is not proof that IQE at 10 K equals
unity.  Nothing in this module reads, imports, or solves for that number;
``yield_ratio`` below computes an independent MODEL proxy that a caller may
print next to 0.52 for context, never as a required match.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi


class NitrideNanowireSurfaceError(ValueError):
    """A surface-loss input is non-finite or outside this model's domain."""


_VALID_SHELLS = ("none", "AlGaN")


def _finite(name, value):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not isfinite(value):
        raise NitrideNanowireSurfaceError(f"{name} must be a finite number")
    return float(value)


def _unit_interval(name, value):
    v = _finite(name, value)
    if not 0.0 <= v <= 1.0:
        raise NitrideNanowireSurfaceError(f"{name} must be in [0, 1]")
    return v


@dataclass(frozen=True)
class NitrideNanowireSurfaceParams:
    """Explicit sidewall-loss design inputs; see module docstring for
    provenance of every default.

    ``shell='none'`` always uses a multiplier of 1.0 (unpassivated), even if
    ``shell_multiplier`` were changed from its AlGaN default; the AlGaN
    factor only applies when ``shell='AlGaN'``.
    """
    S_cm_s: float = 1e3
    shell: str = "none"
    shell_multiplier: float = 0.1
    reservoir_access: float = 1.0
    occupied_dot_access: float = 1.0

    def __post_init__(self):
        if self.shell not in _VALID_SHELLS:
            raise NitrideNanowireSurfaceError(f"shell must be one of {_VALID_SHELLS}")
        s = _finite("S_cm_s", self.S_cm_s)
        if s < 0.0:
            raise NitrideNanowireSurfaceError("S_cm_s must be non-negative")
        _unit_interval("shell_multiplier", self.shell_multiplier)
        _unit_interval("reservoir_access", self.reservoir_access)
        _unit_interval("occupied_dot_access", self.occupied_dot_access)

    def shell_factor(self) -> float:
        """The multiplier actually applied: 1.0 for shell='none' regardless
        of ``shell_multiplier``'s stored value; ``shell_multiplier`` only for
        shell='AlGaN'. [A]"""
        return 1.0 if self.shell == "none" else float(self.shell_multiplier)


def surface_rates(params: NitrideNanowireSurfaceParams, *, core_radius_nm: float, T_K: float) -> dict:
    """Uniform cylindrical sidewall loss, split into access-weighted channels.

    k_side_ns = 2 * S_cm_s / (core_radius_nm * 1e-7) * 1e-9  [DR unit
    conversion of the surface/volume expression 2S/R from cm/s and cm to a
    per-nanosecond rate: R_cm = core_radius_nm * 1e-7 (1 nm = 1e-7 cm); the
    trailing 1e-9 converts s^-1 to ns^-1 (1 ns = 1e-9 s)].

    T_K is accepted (and validated) for interface symmetry with the other
    nitride rate functions and for future callers that key provenance
    strings on temperature; this uniform-cylinder model has no explicit
    temperature dependence of its own -- no activation energy or
    temperature law is fit or implied here.
    """
    if not isinstance(params, NitrideNanowireSurfaceParams):
        raise NitrideNanowireSurfaceError("params must be a NitrideNanowireSurfaceParams")
    r_nm = _finite("core_radius_nm", core_radius_nm)
    if r_nm <= 0.0:
        raise NitrideNanowireSurfaceError("core_radius_nm must be positive")
    t_k = _finite("T_K", T_K)
    if t_k <= 0.0:
        raise NitrideNanowireSurfaceError("T_K must be positive")

    r_cm = r_nm * 1e-7
    k_side_ns = 2.0 * params.S_cm_s / r_cm * 1e-9

    shell_mult = params.shell_factor()
    k_surface_reservoir_ns = k_side_ns * shell_mult * params.reservoir_access
    k_surface_X_ns = k_side_ns * shell_mult * params.occupied_dot_access
    k_surface_XX_ns = 2.0 * k_surface_X_ns  # [A] independent-carrier cascade, see module docstring

    return {
        "k_side_ns": k_side_ns,
        "k_surface_reservoir_ns": k_surface_reservoir_ns,
        "k_surface_X_ns": k_surface_X_ns,
        "k_surface_XX_ns": k_surface_XX_ns,
        "shell_multiplier_used": shell_mult,
        "reservoir_access_used": params.reservoir_access,
        "occupied_dot_access_used": params.occupied_dot_access,
        "core_radius_nm": r_nm,
        "T_K": t_k,
        "provenance": {
            "S_cm_s": "[E] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.1, secondary attribution to ref. 35",
            "k_side_ns_formula": "[DR] 2S/R surface/volume unit conversion, Sze & Ng ch. 1",
            "geometry_transfer": "[A] uniform-cylinder bulk expression applied to reservoir/occupied-dot carriers; not itself verified by the S citation",
            "shell_multiplier": "[A] AlGaN passivation factor 0.1, or 1.0 for shell='none'",
            "reservoir_access / occupied_dot_access": "[A] conservative full access (1.0) absent a measured localization length",
            "k_surface_XX = 2 * k_surface_X": "[A] independent-carrier cascade convention",
        },
    }


def yield_ratio(*, gamma_300_ns: float, loss_300_ns: float, gamma_10_ns: float, loss_10_ns: float) -> dict:
    """Pure helper: the ratio of radiative yields (radiative rate / total
    rate) between two conditions, e.g. 300 K and 10 K.

    This is a MODEL PROXY.  Its baseline interpretation equals a measured PL
    intensity ratio ONLY under equal absorption, capture, and collection
    efficiency between the two conditions [A] -- an assumption this function
    neither checks nor needs to hold.  A caller comparing this proxy to the
    Deshpande 2013 ensemble ratio 0.52 (see module docstring) must record
    the model proxy and the measured ensemble ratio as two separate numbers;
    this function has no knowledge of 0.52 and never gates on it.

    Returns a dict with "ratio" set to a float when defined, or None with a
    "reason" string when the 10 K reference yield is zero (would divide by
    zero) -- e.g. no radiative channel at all at the reference condition.
    """
    g300 = _finite("gamma_300_ns", gamma_300_ns)
    l300 = _finite("loss_300_ns", loss_300_ns)
    g10 = _finite("gamma_10_ns", gamma_10_ns)
    l10 = _finite("loss_10_ns", loss_10_ns)
    for name, v in (("gamma_300_ns", g300), ("loss_300_ns", l300), ("gamma_10_ns", g10), ("loss_10_ns", l10)):
        if v < 0.0:
            raise NitrideNanowireSurfaceError(f"{name} must be non-negative")

    total_300 = g300 + l300
    total_10 = g10 + l10

    yield_300 = g300 / total_300 if total_300 > 0.0 else None
    yield_10 = g10 / total_10 if total_10 > 0.0 else None

    result = {
        "yield_300K": yield_300,
        "yield_10K": yield_10,
        "ratio": None,
        "defined": False,
        "reason": None,
        "provenance": "[A] baseline interpretation as a PL intensity ratio requires equal absorption/capture/collection between conditions; not itself verified here",
    }
    if total_300 <= 0.0:
        result["reason"] = "zero total rate at 300 K (undefined yield_300K)"
        return result
    if total_10 <= 0.0:
        result["reason"] = "zero total rate at 10 K (undefined yield_10K)"
        return result
    if yield_10 == 0.0:
        result["reason"] = "zero reference (10 K) radiative yield"
        return result

    result["ratio"] = yield_300 / yield_10
    result["defined"] = True
    return result
