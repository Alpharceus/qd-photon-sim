"""Shared frozen data types of the composable engine (design 2.2-2.5).

Units follow the repo: rates in 1/ns, energies in meV, times in ps unless a
name says otherwise. Fields may be added with defaults, never removed
(later pieces build against these names in parallel).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Optional


@dataclass(frozen=True)
class EmitterModel:
    """Material.build() output (design 2.2)."""
    E_X_meV: float
    lambda_nm: float
    n_host: float
    gamma0_ns: float                 # bare radiative rate Gamma0
    k_nr_ns: float                   # photonless loss rate
    gamma_star_ns: float             # pure dephasing, repo convention (coherence decays at (Gamma+gamma*)/2)
    phonon: Any                      # qd_gf.PhononParams
    delta_xx_meV: float              # E_XX - 2E_X in the resonant_source.FockSource convention: NEGATIVE = bound (binding energy E_b = -delta_xx_meV)
    gamma_xx_ratio: float
    fss_ueV: float
    transitions: frozenset
    dipole: str
    dEdT_meV_K: float
    sigma_sd_ueV: float              # spectral diffusion width
    tau_sd_ns: float
    inputs: dict = field(default_factory=dict)


@dataclass(frozen=True)
class StructureResponse:
    """Structure.response() output (design 2.3)."""
    F_P: float
    kappa_meV: Optional[float]
    E_mode_meV: Optional[float]
    g_meV: Optional[float]
    F_bg: float = 1.0
    beta0: float = 0.0
    eta_out: float = 1.0
    eta_out_fn: Optional[Callable[[float], float]] = None
    pol_factor: float = 1.0
    notes: tuple = ()
    inputs: dict = field(default_factory=dict)
    tag_floor: str = ""              # widest provenance tag of constants hidden inside response() ("" = none); tag_of applies it


@dataclass(frozen=True)
class Pulse:
    shape: str
    fwhm_ps: float                   # INTENSITY FWHM (physics brief pitfall 2)
    area_rad: float
    detuning_meV: float = 0.0        # laser minus dot is the paper convention; the solver converts (brief pitfall 1)
    chirp_meV_ps: float = 0.0
    t0_ps: float = 0.0
    phase_rad: float = 0.0


@dataclass(frozen=True)
class DrivePlan:
    """Drive.plan() output (design 2.4)."""
    family: str
    timing: str
    f_rep_MHz: Optional[float]
    target: str = "X"
    pulses: tuple = ()
    omega_cw_ns: float = 0.0
    pump_rate_ns: float = 0.0
    tau_on_ns: float = 0.0
    pump_ratio: float = 1.0
    relax_ns: float = 0.0
    interface: Any = None            # drive_mech.DriveInterface | None
    needs: frozenset = frozenset()   # {"phonons","cascade","td_detuning","two_colour"}
    inputs: dict = field(default_factory=dict)
    # devlib-p10e electrical fields (all default = the pre-existing behaviour)
    pump_profile: str = "rect"       # incoherent pump time profile: "rect" | "rc" (rect + RC rise/tail) | "gauss"
    t_rise_ns: float = 0.0           # RC time constant R_s C_j of the "rc" profile
    k_esc_on_ns: float = 0.0         # photonless escape X -> G during the pump window (bias on) [1/ns]
    k_esc_off_ns: float = 0.0        # ... and in the dark window (bias below turn-on; EL decay = Gamma + k_off) [1/ns]
    gate_ns: Optional[float] = None  # counting gate from the pump onset (None = whole period)
    det_profile: Optional[tuple] = None  # Stark gate (det_off_meV, det_on_meV, t_gate_ps, t_r_ps), REPO det_X = omega_X - omega_L


@dataclass(frozen=True)
class Rates:
    """rates.combine() output (design 1.4): Purcell-modified rates."""
    Gamma_ns: float
    gamma_star_ns: float
    k_nr_ns: float
    L: float                         # self-consistent cavity/linewidth overlap
    F_P_eff: float                   # F_P * L
    beta: float
    regime: str = "weak"             # weak | strong
    notes: tuple = ()


@dataclass(frozen=True)
class SolverResult:
    """Solver.run() output (design 2.5)."""
    g2_0: float
    n_ph: float
    I_single: float
    inversion: Optional[float] = None
    extra: dict = field(default_factory=dict)
    valid: bool = True
    reasons: tuple = ()


@dataclass(frozen=True)
class SetupChain:
    """Collection / detection chain of a Setup entry."""
    NA: float
    eta_pol: float
    eta_optics: float
    eta_fibre: float
    eta_det: float
    laser_leak: float
    dark_cps: float
    zpl_filter: bool
