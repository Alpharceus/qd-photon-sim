"""Drive registry entries of the composable engine (design 2.4, physics brief 2.3).

A Drive only PLANS: `plan(vals, em, rt) -> DrivePlan`; the solvers (P5) execute.
Importing this module registers every entry on the global REGISTRY.

Conventions (physics brief section 2, common rules):
  * User-facing `detuning_meV` is the PAPER convention Delta = omega_L - omega_X
    (positive = laser BLUE of the dot). The plan stores the REPO convention
    det_X = omega_X - omega_L = -Delta in `Pulse.detuning_meV` (and, for cw, in
    `plan.inputs["det_X_meV"]`). The paper value is kept in
    `plan.inputs["Delta_paper_meV"]`. NOTE: the comment in types.Pulse says the
    solver converts; here the conversion is already done, so P5 must NOT flip
    the sign again (see the report).
  * `fwhm_ps` is the INTENSITY FWHM of Omega^2; the amplitude sigma is
    resonant_source.gaussian_sigma (plan.inputs["sigma_ps"] for gauss).
  * Rates 1/ns, energies meV, times ps unless a name says otherwise.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

from .. import drive_mech
from ..lindblad import HBAR_MEV_NS, mev_to_rate
from ..resonant_source import gaussian_sigma
from .params import ParamSpec
from .registry import Drive, register
from .types import DrivePlan, Pulse

_P = ParamSpec

# observable id sets (design 2.6); cw drives exclude the per-pulse brightness (B_*) ids
_COMMON = frozenset({"g2_0", "rate_detected_Hz", "F_P_eff", "F_P_peak", "tau_ps", "beta", "eta_extract",
                     "gamma_star_ns", "Z_phonon", "linewidth_ueV", "lambda_nm", "detuning_ueV"})
_PULSED = frozenset({"hom_I", "hom_V", "B_first_lens", "B_fibre"})
OBS_COHERENT_PULSED = _COMMON | _PULSED | {"inversion", "hom_vs_dt"}
OBS_INCOHERENT_PULSED = _COMMON | _PULSED | {"hom_vs_dt"}
OBS_CW_RESONANT = _COMMON | {"g2_tau", "inversion"}
OBS_INCOHERENT_CW = _COMMON | {"hom_I", "hom_V"}
OBS_ELECTRICAL_PULSED = _COMMON | _PULSED
OBS_ELECTRICAL_CW = _COMMON | {"hom_I", "hom_V"}


@dataclass(frozen=True)
class PlannedDrive(Drive):
    """Drive whose plan() is a module-level function `planner(drive, vals, em, rt)`.
    `reason` is non-empty exactly for status == "needs_physics" entries."""
    planner: Optional[Callable] = None
    reason: str = ""

    def plan(self, vals: dict, em, rt):
        if self.status == "needs_physics":
            raise NotImplementedError(f"drive {self.kind!r}: awaiting physics: {self.reason}")
        return self.planner(self, vals, em, rt)


def _getter(drive, vals):
    def g(name):
        return vals[name] if name in vals else drive.param(name).default
    return g


def _positive(name, x):
    x = float(x)
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be a positive finite number, got {x!r}")
    return x


def _nonneg(name, x):
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError(f"{name} must be a non-negative finite number, got {x!r}")
    return x


def _pulse_inputs(shape, fwhm_ps):
    d = {"shape": shape, "fwhm_ps": fwhm_ps}
    if shape == "gauss":
        # amplitude sigma of Omega(t) ~ exp(-t^2/2 sigma^2); gaussian_sigma is unit-agnostic
        d["sigma_ps"] = float(gaussian_sigma(fwhm_ps))
    return d


_PAPER_CONV = ("Paper convention Delta = omega_L - omega_X, positive = laser BLUE of the dot; the plan "
               "stores the repo det_X = -Delta [DR].")
_SRC_A = "assumed default; set to the experiment [A]"
_REF_BRIEF = "physics brief 2.3 (.workers/devlib/p1-physics-brief.md)"


# ======================================================================
# resonant_pulse
# ======================================================================
def _plan_resonant_pulse(self, vals, em, rt):
    g = _getter(self, vals)
    shape = g("shape")
    fwhm = _positive("fwhm_ps", g("fwhm_ps"))
    area_pi = _nonneg("area_pi", g("area_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    delta = float(g("detuning_meV"))
    target = g("target")
    if target == "X+" and "X+" not in em.transitions:
        raise ValueError(f"resonant_pulse target X+ needs 'X+' in the material transitions, "
                         f"got {sorted(em.transitions)}")
    leak = _nonneg("laser_leak", g("laser_leak")) if g("cross_pol") else 0.0
    pulse = Pulse(shape=shape, fwhm_ps=fwhm, area_rad=area_pi * math.pi, detuning_meV=-delta)
    inputs = {**_pulse_inputs(shape, fwhm), "area_pi": area_pi, "Delta_paper_meV": delta,
              "det_X_meV": -delta, "detuning_convention": "repo det_X = omega_X - omega_L = -Delta_paper",
              "cross_pol": bool(g("cross_pol")), "laser_leak": leak}
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target=target,
                     pulses=(pulse,), needs=frozenset(), inputs=inputs)


RESONANT_PULSE = PlannedDrive(
    kind="resonant_pulse", label="Resonant pulse (rect / Gaussian)",
    summary="Coherent pulsed excitation of the exciton: area, detuning and intensity FWHM; "
            "optional cross-polarised collection with a laser-leak background.",
    params=(
        _P("shape", "", "gauss", choices=("gauss", "rect"), tag="A", source=_SRC_A, sweepable=False,
           help="Pulse envelope of the Rabi frequency."),
        _P("fwhm_ps", "ps", 3.0, lo=0.1, hi=200.0, scale="log", tag="A", source=_SRC_A,
           help="INTENSITY FWHM of Omega^2 (what papers quote as the laser autocorrelation width); "
                "the amplitude FWHM is sqrt(2) larger."),
        _P("area_pi", "pi", 1.0, lo=0.0, hi=20.0, tag="A", source=_SRC_A,
           help="Pulse area in units of pi (1 = pi pulse)."),
        _P("detuning_meV", "meV", 0.0, lo=-10.0, hi=10.0, tag="A", source=_SRC_A,
           help="Laser detuning, paper convention Delta = omega_L - omega_X (positive = laser BLUE "
                "of the dot). Stored in the plan as det_X = -Delta."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _P("target", "", "X", choices=("X", "X+"), tag="A", source=_SRC_A, sweepable=False,
           help="Transition driven; X+ needs a material that lists it."),
        _P("cross_pol", "", False, tag="A", source=_SRC_A, sweepable=False, group="advanced",
           help="Cross-polarised excitation/collection: laser light leaks in as a background."),
        _P("laser_leak", "1", 0.0, lo=0.0, hi=1.0, tag="A", source=_SRC_A, group="advanced",
           visible_if={"cross_pol": [True]},
           help="Residual laser as a fraction of the signal (used only when cross_pol is on; "
                "do not also enter it in the setup)."),
    ),
    refs=(_REF_BRIEF + ": strict resonant pi row; Rabi area theorem",),
    provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_resonant_pulse)


# ======================================================================
# tpe
# ======================================================================
def _plan_tpe(self, vals, em, rt):
    g = _getter(self, vals)
    if "XX" not in em.transitions:
        raise ValueError(f"tpe needs the biexciton: 'XX' not in the material transitions "
                         f"{sorted(em.transitions)}")
    shape = g("shape")
    fwhm = _positive("fwhm_ps", g("fwhm_ps"))
    area_pi = _nonneg("area_pi", g("area_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    delta = float(g("detuning_meV"))
    pulse = Pulse(shape=shape, fwhm_ps=fwhm, area_rad=area_pi * math.pi, detuning_meV=-delta)
    inputs = {**_pulse_inputs(shape, fwhm), "area_pi": area_pi,
              "area_relative_to_tpe_pi": True,
              "tpe_area_note": "area_rad = area_pi*pi where '1' is ONE TPE pi area; the solver multiplies "
                               "by resonant_source.tpe_pi_area(delta_ns, shape, duration)/pi",
              "Delta_paper_meV": delta, "det_X_meV": -delta,
              "detuning_convention": "offset from the two-photon resonance; repo sign det = -Delta_paper",
              "delta_xx_meV": float(em.delta_xx_meV)}
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="XX_tpe",
                     pulses=(pulse,), needs=frozenset({"cascade"}), inputs=inputs)


TPE = PlannedDrive(
    kind="tpe", label="Two-photon excitation of the biexciton",
    summary="Resonant two-photon excitation of XX via a virtual X level at half the binding energy; "
            "area is in units of the TPE pi area.",
    params=(
        _P("shape", "", "gauss", choices=("gauss", "rect"), tag="A", source=_SRC_A, sweepable=False,
           help="Pulse envelope of the Rabi frequency."),
        _P("fwhm_ps", "ps", 10.0, lo=1.0, hi=200.0, scale="log", tag="A", source=_SRC_A,
           help="INTENSITY FWHM of Omega^2."),
        _P("area_pi", "pi", 1.0, lo=0.0, hi=5.0, tag="A", source=_SRC_A,
           help="Area in units of the TPE pi area (computed by the solver via tpe_pi_area)."),
        _P("detuning_meV", "meV", 0.0, lo=-2.0, hi=2.0, tag="A", source=_SRC_A,
           help="Laser offset from the two-photon resonance, paper convention (positive = blue); "
                "0 sits at E_XX_total/2."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
    ),
    refs=(_REF_BRIEF + ": TPE row", "Mueller et al., Phys. Rev. Lett. 2014 (cited in the brief)"),
    requires=frozenset({"XX"}), provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_tpe)


# ======================================================================
# cw_resonant
# ======================================================================
def _plan_cw_resonant(self, vals, em, rt):
    g = _getter(self, vals)
    delta = float(g("detuning_meV"))
    s = _nonneg("saturation_s", g("saturation_s"))
    inputs = {"Delta_paper_meV": delta, "det_X_meV": -delta,
              "detuning_convention": "repo det_X = omega_X - omega_L = -Delta_paper"}
    if s > 0.0:
        # s = Omega^2 T1 T2 (brief 2.3 cw row), T1 = 1/(Gamma + k_nr), coherence decays at
        # (Gamma + k_nr + gamma*_repo)/2 -> T2 = 2/(Gamma + k_nr + gamma*) [DR]
        g1 = float(rt.Gamma_ns) + float(rt.k_nr_ns)
        g2c = 0.5 * (g1 + float(rt.gamma_star_ns))
        omega_ns = math.sqrt(s * g1 * g2c)
        inputs["omega_source"] = "saturation_s"
        inputs["saturation_s"] = s
    else:
        omega_ns = float(mev_to_rate(_nonneg("omega_ueV", g("omega_ueV")) * 1e-3))
        inputs["omega_source"] = "omega_ueV"
    inputs["omega_ueV"] = omega_ns * HBAR_MEV_NS * 1e3
    return DrivePlan(family="coherent", timing="cw", f_rep_MHz=None, target="X",
                     omega_cw_ns=float(omega_ns), needs=frozenset(), inputs=inputs)


def drive_mech_hbar_ueV() -> float:
    """hbar in ueV ns (6.582e-4 meV ns = 0.6582 ueV ns)."""
    return float(1e3 / mev_to_rate(1.0))


CW_RESONANT = PlannedDrive(
    kind="cw_resonant", label="CW resonant drive",
    summary="Continuous resonant laser (resonance fluorescence): set the Rabi energy or the "
            "saturation parameter s = Omega^2 T1 T2, plus a detuning.",
    params=(
        _P("omega_ueV", "ueV", 1.0, lo=0.0, hi=100.0, tag="A", source=_SRC_A,
           help="Rabi energy hbar*Omega; ignored when the saturation parameter is > 0."),
        _P("saturation_s", "1", 0.0, lo=0.0, hi=1000.0, tag="A", source=_SRC_A,
           help="Saturation parameter s = Omega^2 T1 T2; 0 means use omega_ueV instead."),
        _P("detuning_meV", "meV", 0.0, lo=-1.0, hi=1.0, tag="A", source=_SRC_A,
           help=_PAPER_CONV),
    ),
    refs=(_REF_BRIEF + ": cw resonance fluorescence row (Loudon / Mollow forms)",),
    provides=frozenset({"coherent_cw"}),
    observables=OBS_CW_RESONANT, planner=_plan_cw_resonant)


# ======================================================================
# above_band / p_shell (incoherent)
# ======================================================================
def _incoherent_plan(self, vals, em, rt, relax_name, shell):
    g = _getter(self, vals)
    pump = _nonneg("pump", g("pump"))
    cw = bool(g("cw"))
    relax_ps = _positive(relax_name, g(relax_name))
    cap2 = bool(g("cap2")) if not shell else False
    inputs = {"pump_P_over_Psat": pump, "cap2": cap2, "reexcitation": not shell,
              "relax_param": relax_name}
    if shell:
        inputs["b_dark"] = float(g("b_dark"))
        if not 0.0 <= inputs["b_dark"] <= 1.0:
            raise ValueError(f"b_dark must be in [0, 1], got {inputs['b_dark']!r}")
    # FX2b: trap refilling, uncorrelated background, thermal escape/recapture (all [A]/[E]; 0 = the original scheme)
    names = {p.name for p in self.params}
    for nm, lo, hi in (("p_refill", 0.0, 0.9), ("bg_frac", 0.0, 0.99), ("E_act_meV", 0.0, math.inf),
                       ("k0_ns", 0.0, math.inf)):
        if nm in names:
            v = float(g(nm))
            if not (math.isfinite(v) and lo <= v <= hi):
                raise ValueError(f"{nm} must be in [{lo}, {hi}], got {v!r}")
            inputs[nm] = v
    if "tau_refill_ps" in names:
        inputs["tau_refill_ps"] = _positive("tau_refill_ps", g("tau_refill_ps"))
    needs = frozenset({"cascade"}) if cap2 else frozenset()
    if cw:
        # P/P_sat = 1 means the pump rate equals the total X decay rate [A]
        rate = pump * (float(rt.Gamma_ns) + float(rt.k_nr_ns))
        return DrivePlan(family="incoherent", timing="cw", f_rep_MHz=None, target="X",
                         pump_rate_ns=rate, tau_on_ns=0.0, pump_ratio=pump,
                         relax_ns=relax_ps * 1e-3, needs=needs, inputs=inputs)
    window_ps = _positive("window_ps", g("window_ps"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    tau_on = window_ps * 1e-3
    # rect window: mean excitations per pulse p*tau_on = P/P_sat, so p_X = 1 - exp(-P/P_sat) [DR]
    inputs["window_ps"] = window_ps
    return DrivePlan(family="incoherent", timing="pulsed", f_rep_MHz=f_rep, target="X",
                     pump_rate_ns=pump / tau_on, tau_on_ns=tau_on, pump_ratio=pump,
                     relax_ns=relax_ps * 1e-3, needs=needs, inputs=inputs)


def _plan_above_band(self, vals, em, rt):
    return _incoherent_plan(self, vals, em, rt, "tau_cap_ps", shell=False)


def _plan_p_shell(self, vals, em, rt):
    return _incoherent_plan(self, vals, em, rt, "tau_relax_ps", shell=True)


_PUMP = _P("pump", "1", 1.0, lo=0.0, hi=100.0, tag="A", source=_SRC_A,
           help="Pump strength P/P_sat; brightness saturates as 1 - exp(-P/P_sat). In cw it is the pump "
                "rate in units of the total X decay rate [A].")
_CW = _P("cw", "", False, tag="A", source=_SRC_A, sweepable=False, help="Continuous pump instead of pulses.")
_WIN = _P("window_ps", "ps", 10.0, lo=0.1, hi=1000.0, scale="log", tag="A", source=_SRC_A,
          visible_if={"cw": [False]}, help="Length of the rectangular pump window.")
_REFILL = (
    _P("p_refill", "1", 0.0, lo=0.0, hi=0.9, tag="A", source=_SRC_A, group="advanced",
       help="Probability per emission of a second, uncorrelated capture after the first emission (trap "
            "refilling). Small-p limit: g2 = 2 p_refill / p_loaded (exactly, for a geometric refill chain) [DR]."),
    _P("tau_refill_ps", "ps", 100.0, lo=1.0, hi=10000.0, scale="log", tag="A", source=_SRC_A, group="advanced",
       visible_if={"cw": [False]}, help="Delay of the refilling capture after the emission."),
    _P("bg_frac", "1", 0.0, lo=0.0, hi=0.99, tag="A", source=_SRC_A, group="advanced",
       help="Uncorrelated background as a fraction of the detected light; g2 = 1 - rho^2 (1 - g2_signal) "
            "with rho = 1 - bg_frac. HOM and brightness are the signal-only values."),
)

_FREP = _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           visible_if={"cw": [False]}, help="Pulse repetition rate.")

ABOVE_BAND = PlannedDrive(
    kind="above_band", label="Above-band (incoherent) pump",
    summary="Non-resonant optical pumping into a carrier reservoir that relaxes into the dot with a "
            "capture time; the capture time is the timing jitter. Re-excitation while the pump is on.",
    params=(
        _PUMP, _CW, _WIN,
        _P("tau_cap_ps", "ps", 30.0, lo=1.0, hi=300.0, scale="log", tag="E",
           source="capture/relaxation jitter 10-100 ps class range [E], physics brief 2.3",
           help="Carrier capture / relaxation time (the HOM jitter)."),
        _FREP,
        _P("cap2", "", False, tag="A", source=_SRC_A, sweepable=False, group="advanced",
           help="Allow loading of the biexciton (capacity 2)."),
    ) + _REFILL,
    refs=(_REF_BRIEF + ": above-band incoherent row",),
    provides=frozenset({"incoherent_pump"}),
    observables=OBS_INCOHERENT_PULSED, planner=_plan_above_band)

P_SHELL = PlannedDrive(
    kind="p_shell", label="p-shell / quasi-resonant pump",
    summary="Quasi-resonant pump into the p shell; fast p->s relaxation sets the jitter and there is "
            "no re-excitation after relaxation. Optional branching into a dark state.",
    params=(
        _PUMP, _CW, _WIN,
        _P("tau_relax_ps", "ps", 10.0, lo=1.0, hi=100.0, scale="log", tag="E",
           source="p->s relaxation 5-30 ps class range [E], physics brief 2.3",
           help="p-shell to s-shell relaxation time (the HOM jitter)."),
        _FREP,
        _P("b_dark", "1", 0.0, lo=0.0, hi=1.0, tag="A", source=_SRC_A, group="advanced",
           help="Branching fraction of relaxation into a dark state."),
    ) + _REFILL + (
        _P("E_act_meV", "meV", 0.0, lo=0.0, hi=200.0, tag="E", source="Arrhenius activation of the thermal "
           "escape of s-shell carriers to the p shell / wetting layer; [E] class model, no paper activation", group="advanced",
           help="Activation energy of the thermal escape; k_esc(T) = k0 exp(-E_act / kT)."),
        _P("k0_ns", "1/ns", 0.0, lo=0.0, hi=1e6, tag="E", source="Arrhenius prefactor of the thermal escape; "
           "[E] class model; 0 = off (reproduces the T-independent result bit-for-bit)", group="advanced",
           help="Arrhenius prefactor. Escape X -> R (delays/jitters the emission) and recapture G -> R of "
                "thermally released carriers between pulses (uncorrelated extra photons: g2(T) grows)."),
    ),
    refs=(_REF_BRIEF + ": p-shell / quasi-resonant row",),
    provides=frozenset({"incoherent_pump"}),
    observables=OBS_INCOHERENT_PULSED, planner=_plan_p_shell)


# ======================================================================
# electrical_pulsed / electrical_cw: DriveInterface from drive_mech itself
# ======================================================================
# roster label -> mech_from_card name. The roster (drive_mech.roster()) is the source of the
# choice lists and of each mechanism's timing; this table only maps labels to the card factory.
_CARD_NAME = {
    "M-1 poisson-rail": "m1", "M-2 quiet-rail": "m2", "M-3 pulsed": "m3",
    "M-4 SET/metallic": "m4a", "M-4 SET/gated": "m4b", "M-4 SET/turnstile": "m4c",
    "M-5a RTI": "m5a",
}
_ROSTER = [i for i in drive_mech.roster() if i.mechanism in _CARD_NAME]
MECH_PULSED = tuple(i.mechanism for i in _ROSTER if i.timing in ("pulsed", "cycled"))
MECH_CW = tuple(i.mechanism for i in _ROSTER if i.timing == "cw")

# which parameters each mechanism (card name) consumes
_USES = {
    "m1": (),
    "m2": ("Z_drive_ohm", "Z_junction_ohm", "f_reg_Hz", "f_sample_Hz"),
    "m3": ("I_uA", "tau_on_ns", "F_p_rail", "R_drive_ohm", "C_dep_F"),
    "m4a": ("eps_cycle", "radius_nm", "R_T_ohm"),
    "m4b": ("eps_cycle", "radius_nm", "R_T_ohm"),
    "m4c": ("eps_cycle", "radius_nm", "R_T_ohm"),
    "m5a": ("eta_boost", "dg_inj_factor", "F_p_rail"),
}
_MU_IGNORED = ("M-3 pulsed", "M-4 SET/turnstile")      # mu is computed / fixed by these


def _labels_using(name):
    return [lab for lab, c in _CARD_NAME.items() if name in _USES[c]]


def _build_interface(self, vals, mech_choices):
    g = _getter(self, vals)
    mech = g("mechanism")
    if mech not in mech_choices:
        raise ValueError(f"{self.kind}: mechanism {mech!r} not in {list(mech_choices)}")
    card = _CARD_NAME[mech]
    p = {}
    for name in _USES[card]:
        v = g(name)
        # 0 stands for "not given" on the optional pairs (no unit-bearing None in a ParamSpec)
        if name in ("R_drive_ohm", "C_dep_F", "f_reg_Hz", "f_sample_Hz") and not v:
            continue
        p[name] = v
    for a, b in (("R_drive_ohm", "C_dep_F"), ("f_reg_Hz", "f_sample_Hz")):
        if (a in p) != (b in p):
            p.pop(a, None)
            p.pop(b, None)
    iface = drive_mech.mech_from_card(card, p, T_K=float(g("T_junction_K")),
                                      mu=float(g("mu")), eta=float(g("eta_capture")))
    return iface, g


def _plan_electrical(self, vals, em, rt, mech_choices, pulsed):
    iface, g = _build_interface(self, vals, mech_choices)
    inputs = {"mechanism": iface.mechanism, "feasible": bool(iface.feasible),
              "note": "generic electrical; use the RT-edge tier passthrough for the full transport chain"}
    if pulsed:
        f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
        tau_on = float(iface.notes["tau_eff_ns"]) if "tau_eff_ns" in iface.notes else 0.0
        timing = "pulsed"
    else:
        f_rep, tau_on, timing = None, 0.0, "cw"
    return DrivePlan(family="electrical", timing=timing, f_rep_MHz=f_rep, target="X",
                     tau_on_ns=tau_on, interface=iface, needs=frozenset(), inputs=inputs)


def _plan_electrical_pulsed(self, vals, em, rt):
    return _plan_electrical(self, vals, em, rt, MECH_PULSED, pulsed=True)


def _plan_electrical_cw(self, vals, em, rt):
    return _plan_electrical(self, vals, em, rt, MECH_CW, pulsed=False)


_ELEC_SUMMARY = ("generic electrical; use the RT-edge tier passthrough for the full transport chain. "
                 "Builds the drive_mech interface (pump Fano factor, loading, capture) for the chosen "
                 "mechanism.")
_SRC_DM = "drive_mech.mech_from_card default"


def _electrical_params(choices, default_mech, eta_default, pulsed):
    ps = [
        _P("mechanism", "", default_mech, choices=choices, tag="A", source="drive_mech.roster() names",
           sweepable=False, help="Injection mechanism (drive_mech roster)."),
        _P("mu", "1", 0.5, lo=0.0, hi=20.0, tag="A", source=_SRC_DM,
           visible_if={"mechanism": [m for m in choices if m not in _MU_IGNORED]},
           help="Operating point: mean carriers loaded per window."),
        _P("eta_capture", "1", eta_default, lo=0.0, hi=1.0, scale="lin", tag="A", source=_SRC_A,
           help="Delivery / capture efficiency into the dot."),
        _P("T_junction_K", "K", 300.0, lo=1.0, hi=500.0, tag="A", source=_SRC_DM,
           visible_if={"mechanism": [m for m in choices if _CARD_NAME[m] in ("m4a", "m4b", "m4c")]},
           help="Junction / island temperature used to price a Coulomb-blockade element."),
    ]
    if pulsed:
        ps.append(_P("f_rep_MHz", "MHz", 80.0, lo=0.01, hi=2000.0, scale="log", tag="A", source=_SRC_A,
                     help="Electrical pulse repetition rate."))
    spec = {
        "I_uA": _P("I_uA", "uA", 10.0, lo=0.0, hi=1e4, tag="A", source=_SRC_DM, help="Drive current."),
        "tau_on_ns": _P("tau_on_ns", "ns", 1.0, lo=0.0, hi=1e4, tag="A", source=_SRC_DM,
                        help="Current window length."),
        "F_p_rail": _P("F_p_rail", "1", 1.0, lo=0.0, hi=10.0, tag="A", source=_SRC_DM,
                       help="Pump Fano factor of the rail (1 = shot-noise limited)."),
        "R_drive_ohm": _P("R_drive_ohm", "ohm", 0.0, lo=0.0, hi=1e7, tag="A", source=_SRC_DM,
                          help="Rail resistance for the rise time; 0 = not given (needs C_dep_F too)."),
        "C_dep_F": _P("C_dep_F", "F", 0.0, lo=0.0, hi=1e-9, tag="A", source=_SRC_DM,
                      help="Depletion capacitance for the rise time; 0 = not given."),
        "Z_drive_ohm": _P("Z_drive_ohm", "ohm", 1e5, lo=1.0, hi=1e9, scale="log", tag="A", source=_SRC_DM,
                          help="Drive impedance of the quiet rail."),
        "Z_junction_ohm": _P("Z_junction_ohm", "ohm", 1e3, lo=1.0, hi=1e9, scale="log", tag="A",
                             source=_SRC_DM, help="Junction impedance (F9 impedance divider)."),
        "f_reg_Hz": _P("f_reg_Hz", "Hz", 0.0, lo=0.0, hi=1e12, tag="A", source=_SRC_DM,
                       help="Regulation bandwidth; 0 = ideal (needs f_sample_Hz too)."),
        "f_sample_Hz": _P("f_sample_Hz", "Hz", 0.0, lo=0.0, hi=1e12, tag="A", source=_SRC_DM,
                          help="Sampling rate of the regulated rail; 0 = not given."),
        "eps_cycle": _P("eps_cycle", "1", 0.01, lo=0.0, hi=1.0, tag="E", source=_SRC_DM + " [E]",
                        help="Per-cycle error of a single-carrier SET source (F_p ~ eps)."),
        "radius_nm": _P("radius_nm", "nm", 20.0, lo=1.0, hi=500.0, tag="A", source=_SRC_DM,
                        help="Island radius setting the charging energy."),
        "R_T_ohm": _P("R_T_ohm", "ohm", 1e6, lo=1e3, hi=1e9, scale="log", tag="A", source=_SRC_DM,
                      help="Tunnel resistance (needs R_T >> R_Q)."),
        "eta_boost": _P("eta_boost", "1", 2.0, lo=1.0, hi=20.0, tag="A", source=_SRC_DM,
                        help="Capture-efficiency boost of resonant-tunnelling injection."),
        "dg_inj_factor": _P("dg_inj_factor", "1", 0.5, lo=0.0, hi=1.0, tag="A", source=_SRC_DM,
                            help="Injection-broadening scale factor of resonant-tunnelling injection."),
    }
    for name, p in spec.items():
        labs = [m for m in _labels_using(name) if m in choices]
        if labs:
            ps.append(_P(p.name, p.unit, p.default, lo=p.lo, hi=p.hi, scale=p.scale, tag=p.tag,
                         source=p.source, help=p.help, group="advanced", visible_if={"mechanism": labs}))
    return tuple(ps)


ELECTRICAL_PULSED = PlannedDrive(
    kind="electrical_pulsed", label="Electrical, pulsed / cycled",
    summary=_ELEC_SUMMARY, status="ready",
    params=_electrical_params(MECH_PULSED, "M-3 pulsed", 1e-5, True),
    refs=("fsim_core/drive_mech.py (F8-F9 mechanism roster)",),
    provides=frozenset({"electrical_pulsed"}),
    observables=OBS_ELECTRICAL_PULSED, planner=_plan_electrical_pulsed)

ELECTRICAL_CW = PlannedDrive(
    kind="electrical_cw", label="Electrical, cw",
    summary=_ELEC_SUMMARY, status="ready",
    params=_electrical_params(MECH_CW, "M-1 poisson-rail", 0.05, False),
    refs=("fsim_core/drive_mech.py (F8-F9 mechanism roster)",),
    provides=frozenset({"electrical_cw"}),
    observables=OBS_ELECTRICAL_CW, planner=_plan_electrical_cw)


# ======================================================================
# needs_physics placeholders (full ParamSpecs so the UI can show them greyed)
# ======================================================================
_E_BRIEF = "physics brief 2.3 class range [E]"


def _placeholder(kind, label, summary, reason, params, refs, observables):
    return PlannedDrive(kind=kind, label=label,
                        summary=summary + " Awaiting physics: " + reason,
                        status="needs_physics", reason=reason, params=params, refs=refs,
                        observables=observables, planner=None)


_GAUSS = _P("shape", "", "gauss", choices=("gauss", "rect"), tag="A", source=_SRC_A, sweepable=False,
            help="Pulse envelope.")

# P9: the phonon-dependent coherent drives live in drives_phonon.py (real entries, status "ready"); they
# import the helpers above, so this import must stay below them. FX2c: ARP (chirped pulse; polaron solver with
# the PulseTrain chirp) is ready too; `_placeholder` stays for future needs_physics entries.
from . import drives_phonon as _dp  # noqa: E402

PHONON_ASSISTED, SUPER, DICHROMATIC, ARP = _dp.PHONON_ASSISTED, _dp.SUPER, _dp.DICHROMATIC, _dp.ARP

DRIVES = (RESONANT_PULSE, TPE, CW_RESONANT, ABOVE_BAND, P_SHELL, ELECTRICAL_PULSED, ELECTRICAL_CW,
          PHONON_ASSISTED, SUPER, DICHROMATIC, ARP)

for _d in DRIVES:
    register(_d)
