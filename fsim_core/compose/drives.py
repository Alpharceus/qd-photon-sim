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

# display symbols (KaTeX TeX; Scully-Zubairy / Metcalf-van der Straten notation), keyed by the registry id
_TEX = {
    "shape": r"\Omega(t)", "fwhm_ps": r"\tau_{\mathrm{FWHM}}", "area_pi": r"\Theta", "detuning_meV": r"\Delta", "f_rep_MHz": r"f_{\mathrm{rep}}",
    "target": r"\mathrm{target}", "cross_pol": r"\mathrm{cross\text{-}pol}", "laser_leak": r"f_{\mathrm{leak}}",
    "omega_ueV": r"\hbar\Omega", "saturation_s": r"s", "pump": r"P/P_{\mathrm{sat}}", "cw": r"\mathrm{cw}",
    "window_ps": r"T_{\mathrm{win}}", "p_refill": r"p_{\mathrm{refill}}", "tau_refill_ps": r"\tau_{\mathrm{refill}}",
    "bg_frac": r"f_{\mathrm{bg}}", "tau_cap_ps": r"\tau_{\mathrm{cap}}", "cap2": r"\mathrm{cap}_{2}",
    "tau_relax_ps": r"\tau_{\mathrm{relax}}", "b_dark": r"b_{\mathrm{dark}}", "E_act_meV": r"E_{\mathrm{act}}", "k0_ns": r"k_{0}",
    "pump_profile": r"r(t)", "t_rise_ps": r"t_{\mathrm r}", "k_esc_on_ns": r"k_{\mathrm{esc}}^{\mathrm{on}}",
    "k_esc_off_ns": r"k_{\mathrm{esc}}^{\mathrm{off}}", "gate_ps": r"T_{\mathrm{gate}}", "hom_gate_ps": r"T_{\mathrm{HOM}}",
    "hom_gate_shape": r"\mathrm{gate\ shape}", "hom_gate_start_ps": r"t_{\mathrm{HOM}}", "hom_jitter_ps": r"\sigma_{\mathrm{jit}}",
    "irf_ps": r"\tau_{\mathrm{IRF}}", "cascade": r"\mathrm{cascade}", "pump_ratio_xx": r"\rho_{XX}",
    "mechanism": r"\mathrm{mechanism}", "mu": r"\mu", "eta_capture": r"\eta_{\mathrm{cap}}", "T_junction_K": r"T_{\mathrm j}",
    "f_rep_MHz_el": r"f_{\mathrm{rep}}", "I_uA": r"I", "tau_on_ns": r"\tau_{\mathrm{on}}", "F_p_rail": r"F_{\mathrm{rail}}",
    "R_drive_ohm": r"R_{\mathrm{drive}}", "C_dep_F": r"C_{\mathrm{dep}}", "Z_drive_ohm": r"Z_{\mathrm{drive}}",
    "Z_junction_ohm": r"Z_{\mathrm j}", "f_reg_Hz": r"f_{\mathrm{reg}}", "f_sample_Hz": r"f_{\mathrm s}",
    "eps_cycle": r"\varepsilon_{\mathrm{cyc}}", "radius_nm": r"R_{\mathrm{island}}", "R_T_ohm": r"R_{\mathrm T}",
    "eta_boost": r"\eta_{\mathrm{boost}}", "dg_inj_factor": r"f_{\mathrm{inj}}", "det_off_meV": r"\Delta_{\mathrm{off}}",
    "det_on_meV": r"\Delta_{\mathrm{on}}", "t_gate_ps": r"t_{\mathrm{gate}}", "t_r_ps": r"t_{\mathrm r}",
    "rabi_period_ns": r"T_{\Omega}", "power_ref_nW": r"P_{\mathrm{ref}}", "power_nW": r"P",
    "eid_ps": r"\xi_{\mathrm{EID}}", "p_trion": r"p_{\mathrm{trion}}",
    "delta1_meV": r"\Delta_{1}", "delta2_meV": r"\Delta_{2}", "fwhm1_ps": r"\tau_{\mathrm{FWHM},1}", "fwhm2_ps": r"\tau_{\mathrm{FWHM},2}",
    "area1_pi": r"\Theta_{1}", "area2_pi": r"\Theta_{2}", "delta_meV": r"\delta", "asymmetry": r"a", "chirp_ps2": r"\varphi_{2}",
    "slm_T1": r"T_{1}", "slm_T2": r"T_{2}", "slm_alpha_pi": r"\tilde{\alpha}", "slm_bb_fwhm_meV": r"\hbar\Delta\omega_{\mathrm{bb}}",
    "slm_bb_centre_meV": r"\hbar(\omega_{\mathrm C}-\omega_{0})", "slm_mask_sigma_meV": r"\hbar\sigma", "slm_floor": r"t_{\mathrm{floor}}",
}


def _P(name, *args, **kw):
    """ParamSpec with the display symbol of this registry id (tex=None when the id has none)."""
    kw.setdefault("tex", _TEX.get(name))
    return ParamSpec(name, *args, **kw)

# observable id sets (design 2.6); cw drives exclude the per-pulse brightness (B_*) ids
_COMMON = frozenset({"g2_0", "rate_detected_Hz", "F_P_eff", "F_P_peak", "tau_ps", "tau_ratio", "beta", "eta_extract",
                     "gamma_star_ns", "Z_phonon", "linewidth_ueV", "lambda_nm", "detuning_ueV"})
_PULSED = frozenset({"hom_I", "hom_V", "B_first_lens", "B_fibre"})
OBS_COHERENT_PULSED = _COMMON | _PULSED | {"inversion", "hom_vs_dt"}
_CW_RATES = frozenset({"B_first_lens_Hz", "B_fibre_Hz"})     # devlib-fx3b: photon flux into the lens / fibre (cw)
OBS_INCOHERENT_PULSED = _COMMON | _PULSED | _CW_RATES | {"hom_vs_dt"}      # above_band / p_shell also run cw
# R2 F6: above_band / p_shell with cw = True are a cw stream: the per-pulse brightness B_*, the per-pulse HOM overlap
# (hom_I, hom_V) and its separation series hom_vs_dt are undefined there, so the cw variant does not advertise them
# (the photon fluxes B_*_Hz and rate_detected_Hz are the cw quantities). Decision: NOT ADVERTISED (r2-f56-decisions.md).
OBS_INCOHERENT_CW_VARIANT = OBS_INCOHERENT_PULSED - frozenset({"hom_I", "hom_V", "B_first_lens", "B_fibre", "hom_vs_dt"})
# R2 FIX-1a (audit 6): a cw stream has no per-pulse HOM overlap (it needs a timing-window definition), so the cw drives no
# longer advertise hom_I / hom_V; the engine returns NaN with a reason if one is requested anyway.
OBS_CW_RESONANT = _COMMON | _CW_RATES | {"g2_tau", "inversion", "s_0", "s_sat", "coherent_fraction"}
OBS_INCOHERENT_CW = _COMMON | _CW_RATES
# devlib-p10e: EL decay time, binned (GHz) histogram g2 and gated HOM for pulsed electrical drives; IRF-convolved g2 for cw
# R2 dead controls: F_dot = delivered pump Fano factor (mechanism, eta_cap, impedances, SET pricing), hom_gate_fraction = photons in
# the HOM gate (gate start / width)
OBS_ELECTRICAL_PULSED = _COMMON | _PULSED | {"tau_EL_ps", "g2_binned", "hom_I_gated", "hom_gate_fraction", "spectral_overlap", "F_dot"}
OBS_ELECTRICAL_CW = _COMMON | _CW_RATES | {"tau_EL_ps", "g2_irf", "spectral_overlap", "F_dot"}


@dataclass(frozen=True)
class PlannedDrive(Drive):
    """Drive whose plan() is a module-level function `planner(drive, vals, em, rt)`.
    `reason` is non-empty exactly for status == "needs_physics" entries."""
    planner: Optional[Callable] = None
    reason: str = ""
    # R2 F6: the advertised set when the drive runs with `cw` = True (None = the same as `observables`)
    observables_cw: Optional[frozenset] = None

    def observables_for(self, vals=None) -> frozenset:
        """Observables advertised for these parameter values (the timing variant decides: a cw stream drops the per-pulse ones)."""
        vals = vals or {}
        if self.observables_cw is not None and "cw" in {p.name for p in self.params}:
            if bool(vals["cw"] if "cw" in vals else self.param("cw").default):
                return self.observables_cw
        return self.observables

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


_PAPER_CONV = (r"Paper convention $\Delta = \omega_{\mathrm L} - \omega_X$, positive = laser BLUE of the dot; the plan "
               r"stores the repo $\mathrm{det}_X = -\Delta$ [DR].")
_SRC_A = "assumed default; set to the experiment [A]"

# devlib-fx3b, both default OFF (= bit-identical): excitation-induced dephasing and a dark-trion branch of coherent pulses
_SRC_FX3B = "Engine extension [E]: class model, the coefficient is NOT taken from a paper"
_EID = _P("eid_ps", "ps", 0.0, lo=0.0, hi=1e3, tag="E", source=_SRC_FX3B, group="advanced",
          help=r"Excitation-induced dephasing: pure dephasing rate $\gamma_{\mathrm{EID}} = \xi_{\mathrm{EID}}\,|\Omega(t)|^2$ ($\Omega$ in rad/ps) "
               r"added to the pure-dephasing linewidth $\Gamma_{\mathrm{ph}} = 2\gamma_{\mathrm{ph}}$ during the pulse (Koong-type HOM loss at large area). 0 = off. Runs on the polaron "
               r"step propagators (pulsed two-level drives, no TPE).")
_TRION = _P("p_trion", "1", 0.0, lo=0.0, hi=0.99, tag="E", source=_SRC_FX3B, group="advanced",
            help=r"Probability, at a $\pi$ pulse, that the excitation ends in a dark trion / excess-charge state "
                 r"(Vajner-type damping of the Rabi oscillation); the per-pulse probability is $p_{\mathrm{trion}}\sin^2(\Theta/2)$. "
                 r"Scales $n_{\mathrm{ph}}$ and the preparation by $(1 - p)$ and $g^{(2)}$ by $1/(1 - p)$. 0 = off.")


def _dark_inputs(g, inputs, eid=True):
    """Add the optional eid_ps / p_trion entries to a coherent plan's inputs (absent when 0)."""
    if eid:
        c = _nonneg("eid_ps", g("eid_ps"))
        if c > 0.0:
            inputs["eid_ps"] = c
    p = _nonneg("p_trion", g("p_trion"))
    if p >= 1.0:
        raise ValueError(f"p_trion must be < 1, got {p!r}")
    if p > 0.0:
        inputs["p_trion"] = p
    return inputs


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
    _dark_inputs(g, inputs)
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
           help=r"INTENSITY FWHM $\tau_{\mathrm{FWHM}}$ of $\Omega^2$ (what papers quote as the laser autocorrelation width); "
                r"the amplitude FWHM is $\sqrt{2}$ larger."),
        _P("area_pi", "pi", 1.0, lo=0.0, hi=20.0, tag="A", source=_SRC_A,
           help=r"Pulse area $\Theta$ in units of $\pi$ (1 = $\pi$ pulse)."),
        _P("detuning_meV", "meV", 0.0, lo=-10.0, hi=10.0, tag="A", source=_SRC_A,
           help=r"Laser detuning $\Delta$, paper convention $\Delta = \omega_{\mathrm L} - \omega_X$ (positive = laser BLUE "
                r"of the dot). Stored in the plan as $\mathrm{det}_X = -\Delta$."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _P("target", "", "X", choices=("X", "X+"), tag="A", source=_SRC_A, sweepable=False,
           help="Transition driven; X+ needs a material that lists it."),
        _P("cross_pol", "", False, tag="A", source=_SRC_A, sweepable=False, group="advanced",
           help="Cross-polarised excitation/collection: laser light leaks in as a background."),
        _P("laser_leak", "1", 0.0, lo=0.0, hi=1.0, tag="A", source=_SRC_A, group="advanced",
           visible_if={"cross_pol": [True]},
           help=r"Residual laser $f_{\mathrm{leak}}$ as a fraction of the signal (used only when cross-pol is on; "
                r"do not also enter it in the setup)."),
        _EID, _TRION,
    ),
    refs=("Rabi area theorem (Scully and Zubairy, Quantum Optics, 1997)",),
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
    _dark_inputs(g, inputs, eid=False)
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
           help=r"INTENSITY FWHM $\tau_{\mathrm{FWHM}}$ of $\Omega^2$."),
        _P("area_pi", "pi", 1.0, lo=0.0, hi=5.0, tag="A", source=_SRC_A,
           help=r"Area $\Theta$ in units of the TPE $\pi$ area (computed by the solver via tpe_pi_area)."),
        _P("detuning_meV", "meV", 0.0, lo=-2.0, hi=2.0, tag="A", source=_SRC_A,
           help=r"Laser offset $\Delta$ from the two-photon resonance, paper convention (positive = blue); "
                r"0 sits at $E_{XX,\mathrm{tot}}/2$."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _TRION,
    ),
    refs=("Two-photon biexciton excitation: Mueller et al., Phys. Rev. Lett. 2014",),
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
        # (Gamma + k_nr + Gamma_ph)/2, Gamma_ph = 2 gamma_ph = code gamma* (pure-dephasing linewidth) -> T2 = 2/(Gamma + k_nr + Gamma_ph) [DR]
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
           help=r"Rabi energy $\hbar\Omega$; ignored when the saturation parameter is $> 0$."),
        _P("saturation_s", "1", 0.0, lo=0.0, hi=1000.0, tag="A", source=_SRC_A,
           help=r"Saturation parameter $s = \Omega^2 T_1 T_2$; 0 means use $\hbar\Omega$ instead."),
        _P("detuning_meV", "meV", 0.0, lo=-1.0, hi=1.0, tag="A", source=_SRC_A,
           help=_PAPER_CONV),
    ),
    refs=("Resonance fluorescence: Loudon, The Quantum Theory of Light; Mollow, Phys. Rev. 188, 1969",),
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
    bright = 1.0
    if shell:
        inputs["b_dark"] = float(g("b_dark"))
        if not 0.0 <= inputs["b_dark"] <= 1.0:
            raise ValueError(f"b_dark must be in [0, 1], got {inputs['b_dark']!r}")
        # R2 F10c [A]: a fraction b_dark of the carriers relaxing p -> s ends in the optically dark exciton (spin flip in the
        # relaxation), so only (1 - b_dark) of the pump loads the bright state: rate_bright = (1 - b_dark) rate. The dark state
        # is not re-emitted (no dark -> bright flip in this model). b_dark = 0 multiplies by exactly 1 (bit-identical).
        bright = 1.0 - inputs["b_dark"]
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
        rate = pump * bright * (float(rt.Gamma_ns) + float(rt.k_nr_ns))
        return DrivePlan(family="incoherent", timing="cw", f_rep_MHz=None, target="X",
                         pump_rate_ns=rate, tau_on_ns=0.0, pump_ratio=pump,
                         relax_ns=relax_ps * 1e-3, needs=needs, inputs=inputs)
    window_ps = _positive("window_ps", g("window_ps"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    tau_on = window_ps * 1e-3
    # rect window: mean excitations per pulse p*tau_on = P/P_sat, so p_X = 1 - exp(-P/P_sat) [DR]
    inputs["window_ps"] = window_ps
    return DrivePlan(family="incoherent", timing="pulsed", f_rep_MHz=f_rep, target="X",
                     pump_rate_ns=pump * bright / tau_on, tau_on_ns=tau_on, pump_ratio=pump,
                     relax_ns=relax_ps * 1e-3, needs=needs, inputs=inputs)


def _plan_above_band(self, vals, em, rt):
    return _incoherent_plan(self, vals, em, rt, "tau_cap_ps", shell=False)


def _plan_p_shell(self, vals, em, rt):
    return _incoherent_plan(self, vals, em, rt, "tau_relax_ps", shell=True)


_PUMP = _P("pump", "1", 1.0, lo=0.0, hi=100.0, tag="A", source=_SRC_A,
           help=r"Pump strength $P/P_{\mathrm{sat}}$; brightness saturates as $1 - e^{-P/P_{\mathrm{sat}}}$. In cw it is the pump "
                r"rate in units of the total X decay rate [A].")
_CW = _P("cw", "", False, tag="A", source=_SRC_A, sweepable=False, help="Continuous pump instead of pulses.")
_WIN = _P("window_ps", "ps", 10.0, lo=0.1, hi=1000.0, scale="log", tag="A", source=_SRC_A,
          visible_if={"cw": [False]}, help="Length of the rectangular pump window.")
_REFILL = (
    _P("p_refill", "1", 0.0, lo=0.0, hi=0.9, tag="A", source=_SRC_A, group="advanced",
       help=r"Probability per emission of a second, uncorrelated capture after the first emission (trap "
            r"refilling). Small-$p$ limit: $g^{(2)} = 2\,p_{\mathrm{refill}}/p_{\mathrm{loaded}}$ (exactly, for a geometric refill chain) [DR]."),
    _P("tau_refill_ps", "ps", 100.0, lo=1.0, hi=10000.0, scale="log", tag="A", source=_SRC_A, group="advanced",
       visible_if={"cw": [False]}, help="Delay of the refilling capture after the emission."),
    _P("bg_frac", "1", 0.0, lo=0.0, hi=0.99, tag="A", source=_SRC_A, group="advanced",
       help=r"Uncorrelated background as a fraction $f_{\mathrm{bg}}$ of the detected light; $g^{(2)} = 1 - \rho^2\,(1 - g^{(2)}_{\mathrm{sig}})$ "
            r"with $\rho = 1 - f_{\mathrm{bg}}$. HOM and brightness are the signal-only values."),
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
           source="capture/relaxation jitter 10-100 ps class range [E]",
           help="Carrier capture / relaxation time (the HOM jitter)."),
        _FREP,
        _P("cap2", "", False, tag="A", source=_SRC_A, sweepable=False, group="advanced",
           help="Allow loading of the biexciton (capacity 2)."),
    ) + _REFILL,
    refs=("Above-band incoherent excitation: class model [E]",),
    provides=frozenset({"incoherent_pump"}),
    observables=OBS_INCOHERENT_PULSED, observables_cw=OBS_INCOHERENT_CW_VARIANT, planner=_plan_above_band)

P_SHELL = PlannedDrive(
    kind="p_shell", label="p-shell / quasi-resonant pump",
    summary="Quasi-resonant pump into the p shell; fast p->s relaxation sets the jitter and there is "
            "no re-excitation after relaxation. Optional branching into a dark state.",
    params=(
        _PUMP, _CW, _WIN,
        _P("tau_relax_ps", "ps", 10.0, lo=1.0, hi=100.0, scale="log", tag="E",
           source="p->s relaxation 5-30 ps class range [E]",
           help="p-shell to s-shell relaxation time (the HOM jitter)."),
        _FREP,
        _P("b_dark", "1", 0.0, lo=0.0, hi=1.0, tag="A", source=_SRC_A, group="advanced",
           help=r"Branching fraction $b_{\mathrm{dark}}$ of the p -> s relaxation into an optically dark state; only $(1 - b_{\mathrm{dark}})$ of the pump loads the bright exciton (the dark state is not re-emitted)."),
    ) + _REFILL + (
        _P("E_act_meV", "meV", 0.0, lo=0.0, hi=200.0, tag="E", source="Arrhenius activation of the thermal "
           "escape of s-shell carriers to the p shell / wetting layer; [E] class model, no paper activation", group="advanced",
           help=r"Activation energy of the thermal escape; $k_{\mathrm{esc}}(T) = k_0\,e^{-E_{\mathrm{act}}/k_{\mathrm B}T}$."),
        _P("k0_ns", "1/ns", 0.0, lo=0.0, hi=1e6, tag="E", source="Arrhenius prefactor of the thermal escape; "
           "[E] class model; 0 = off (reproduces the T-independent result bit-for-bit)", group="advanced",
           help=r"Arrhenius prefactor $k_0$. Escape X $\to$ R (delays/jitters the emission) and recapture G $\to$ R of "
                r"thermally released carriers between pulses (uncorrelated extra photons: $g^{(2)}(T)$ grows)."),
    ),
    refs=("p-shell quasi-resonant excitation: class model [E]",),
    provides=frozenset({"incoherent_pump"}),
    observables=OBS_INCOHERENT_PULSED, observables_cw=OBS_INCOHERENT_CW_VARIANT, planner=_plan_p_shell)


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
    "m5a": ("eta_boost", "F_p_rail"),
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
                                      mu=float(g("mu")) if _has(self, "mu") else 0.5, eta=float(g("eta_capture")))
    return iface, g


def _has(drive, name):
    return any(p.name == name for p in drive.params)


def _plan_electrical(self, vals, em, rt, mech_choices, pulsed):
    iface, g = _build_interface(self, vals, mech_choices)
    inputs = {"mechanism": iface.mechanism, "feasible": bool(iface.feasible),
              "note": "generic electrical drive; the room-temperature edge-emitter tier models the full transport chain"}
    extra = {}
    if pulsed:
        f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
        tau_on = float(iface.notes["tau_eff_ns"]) if "tau_eff_ns" in iface.notes else 0.0
        timing = "pulsed"
        profile = g("pump_profile")
        t_rise = _nonneg("t_rise_ps", g("t_rise_ps")) * 1e-3
        k_on = _nonneg("k_esc_on_ns", g("k_esc_on_ns"))
        k_off = _nonneg("k_esc_off_ns", g("k_esc_off_ns"))
        gate = _nonneg("gate_ps", g("gate_ps")) * 1e-3
        extra.update(pump_profile=profile, k_esc_on_ns=k_on, k_esc_off_ns=k_off, gate_ns=gate if gate > 0.0 else None)
        if profile != "rect":
            # devlib-p10e: an RC tail LENGTHENS the window (brief 1.1, trap 3): never subtract t_r. The nominal window
            # is the entered tau_on_ns (M-3) and mu = eta N_c is the UNSHORTENED operating point; t_r = R C when the
            # mechanism has it, else the t_rise_ps input.
            card = _CARD_NAME[iface.mechanism] if iface.mechanism in _CARD_NAME else ""
            tau_nom = float(g("tau_on_ns")) if (card == "m3" and _has(self, "tau_on_ns")) else tau_on
            if card == "m3" and tau_nom > 0.0:
                t_rise = float(iface.notes.get("t_rise_ns", 0.0)) or t_rise
                full = drive_mech.mech_pulsed(float(g("I_uA")), tau_nom, float(iface.eta_capture), float(iface.F_p))
                extra["pump_rate_ns"] = float(full.mu) / tau_nom           # mu = eta N_c already holds eta_capture (FX3b fix)
                tau_on = tau_nom
            extra["t_rise_ns"] = t_rise if profile == "rc" else 0.0
        hom_gate = _nonneg("hom_gate_ps", g("hom_gate_ps"))
        if hom_gate > 0.0:
            inputs["hom_gate"] = {"start_ns": _nonneg("hom_gate_start_ps", g("hom_gate_start_ps")) * 1e-3,
                                  "width_ns": hom_gate * 1e-3,
                                  "jitter_ns": _nonneg("hom_jitter_ps", g("hom_jitter_ps")) * 1e-3,
                                  "shape": g("hom_gate_shape")}
    else:
        f_rep, tau_on, timing = None, 0.0, "cw"
        extra["k_esc_on_ns"] = _nonneg("k_esc_on_ns", g("k_esc_on_ns"))
    irf = _nonneg("irf_ps", g("irf_ps"))
    if irf > 0.0:
        inputs["irf_ps"] = irf
    bg = _nonneg("bg_frac", g("bg_frac"))
    if bg >= 1.0:
        raise ValueError(f"bg_frac must be < 1, got {bg!r}")
    if bg > 0.0:
        inputs["bg_frac"] = bg
    needs = frozenset()
    if bool(g("cascade")):
        needs = frozenset({"cascade"})
        inputs["pump_ratio_xx"] = _nonneg("pump_ratio_xx", g("pump_ratio_xx"))
    return DrivePlan(family="electrical", timing=timing, f_rep_MHz=f_rep, target="X",
                     tau_on_ns=tau_on, interface=iface, needs=needs, inputs=inputs, **extra)


def _plan_electrical_pulsed(self, vals, em, rt):
    return _plan_electrical(self, vals, em, rt, MECH_PULSED, pulsed=True)


def _plan_electrical_cw(self, vals, em, rt):
    return _plan_electrical(self, vals, em, rt, MECH_CW, pulsed=False)


_ELEC_SUMMARY = ("Generic electrical injection. The chosen injection mechanism sets the pump noise (Fano factor), "
                 "how fully the dot is loaded and the carrier-capture efficiency; the full transport chain "
                 "(wire, contacts, heat) is modelled only in the room-temperature edge-emitter tier.")
_SRC_DM = "injection-mechanism default for the card [A]"


# Repo pointer: devlib-p10e electrical brief, .workers/devlib/p3-electrical-brief.md
_SRC_P10E = "Electrical-drive class model [A]"
# Repo pointer: p3-electrical-brief.md 1.3 (cascade saturation) / 1.5 (background law)
_SRC_FX3B_EL = "Cascade-saturation and background-law class model [A]"


def _electrical_shape_params():
    """Shaped pump, bias-dependent escape, counting gate, gated HOM and IRF inputs (all default = off)."""
    return [
        _P("pump_profile", "", "rect", choices=("rect", "rc", "gauss"), tag="A", source=_SRC_P10E, sweepable=False,
           group="advanced",
           help=r"Time profile of the incoherent electrical pump $r(t)$ at fixed operating point $\mu = \int r\,\mathrm dt$: rect, "
                r"rect + RC rise/tail (diode current through $R_{\mathrm s}C_{\mathrm j}$) or Gaussian (generator pulse). An RC tail "
                r"LENGTHENS the effective window ($\tau_{\mathrm{eff}} \approx \tau_{\mathrm{on}} + 1.7\,t_{\mathrm r}$) and raises the re-excitation $g^{(2)}$."),
        _P("t_rise_ps", "ps", 100.0, lo=0.0, hi=2000.0, tag="A", source="50 ohm x 2 pF class RC [E]; "
           "used only when the mechanism has no R and C of its own", group="advanced",
           visible_if={"pump_profile": ["rc"]},
           help=r"RC time constant $R_{\mathrm s}C_{\mathrm j}$ of the pump profile (the pulse tail $e^{-t/t_{\mathrm r}}$)."),
        _P("k_esc_on_ns", "1/ns", 0.0, lo=0.0, hi=50.0, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Photonless tunnelling escape rate $k_{\mathrm{esc}}^{\mathrm{on}}$, X $\to$ G, while the pump is on (flat band, small)."),
        _P("k_esc_off_ns", "1/ns", 0.0, lo=0.0, hi=50.0, tag="A", source="Ellis et al., Appl. Phys. Lett. 92, 2008: "
           "EL decay 1.59/ns vs TRPL 1.40/ns gives k_off = 0.19/ns [DR from V]; Schlehahn et al. 2015: 3.1/ns [DR from V]",
           group="advanced",
           help=r"Escape rate $k_{\mathrm{esc}}^{\mathrm{off}}$, X $\to$ G, in the dark window (bias below turn-on, higher field). The EL decay time is "
                r"$1/(\Gamma + k_{\mathrm{off}})$; it is reported ($\tau_{\mathrm{EL}}$) and never used as the radiative lifetime."),
        _P("gate_ps", "ps", 0.0, lo=0.0, hi=1e5, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Photon-counting gate $T_{\mathrm{gate}}$ from the pump onset ($g^{(2)}$ and brightness); 0 = count the whole period."),
        _P("hom_gate_ps", "ps", 0.0, lo=0.0, hi=1e4, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Detection-gate width $T_{\mathrm{HOM}}$ of a Stark-gated HOM (Bennett-like); 0 = off. Observable $I_{\mathrm{HOM}}^{\mathrm{gate}}$."),
        _P("hom_gate_shape", "", "rect", choices=("rect", "gauss"), tag="A", source=_SRC_P10E, sweepable=False,
           group="advanced",
           help=r"Shape of the HOM detection gate: rect (top-hat of width $T_{\mathrm{HOM}}$, opening at $t_{\mathrm{HOM}}$) or "
                r"gauss (Gaussian window of INTENSITY FWHM $T_{\mathrm{HOM}}$ CENTRED at $t_{\mathrm{HOM}}$; Bennett's 31 ps "
                r"window is a Gaussian)."),
        _P("hom_gate_start_ps", "ps", 0.0, lo=0.0, hi=1e4, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Time $t_{\mathrm{HOM}}$ after the excitation at which the HOM gate opens (rect) / its centre (gauss)."),
        _P("hom_jitter_ps", "ps", 0.0, lo=0.0, hi=1e3, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Standard deviation $\sigma_{\mathrm{jit}}$ of the relative timing jitter between the two photons' detection gates (Gaussian). The overlap "
                "falls monotonically with it; the displaced windows share only part of the photon."),
        _P("irf_ps", "ps", 0.0, lo=0.0, hi=2000.0, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Detector timing response $\tau_{\mathrm{IRF}}$ (Gaussian FWHM); 0 = off. Adds the IRF-convolved $g^{(2)}$ (cw dip) and the IRF "
                r"leak of the binned histogram ($g^{(2)}_{\mathrm{bin}}$)."),
    ]


def _electrical_cascade_params():
    """Cap-2 cascade (XX) loading and the background-only fraction (devlib-fx3b); default = off."""
    return [
        _P("cascade", "", False, tag="A", source=_SRC_FX3B_EL, sweepable=False, group="advanced",
           help=r"Allow loading of the biexciton (capacity 2, rate solver cascade G-X-XX). Strong pumping then PARKS the dot "
                r"in XX and the pulsed $g^{(2)}$ saturates and falls (Sartison anchor: rect cascade 0.135 at $\mu = 4$, 0.049 at "
                r"$\mu = 20$) instead of rising like the two-level pump."),
        _P("pump_ratio_xx", "1", 1.0, lo=0.0, hi=10.0, tag="A", source=_SRC_FX3B_EL, group="advanced",
           visible_if={"cascade": [True]},
           help=r"Rate $\rho_{XX}$ of loading X $\to$ XX relative to G $\to$ X (1 = equal; 0 = two-level pump)."),
        _P("bg_frac", "1", 0.0, lo=0.0, hi=0.99, tag="A", source=_SRC_FX3B_EL, group="advanced",
           help=r"Uncorrelated BACKGROUND as a fraction $f_{\mathrm{bg}}$ of the detected light (neighbouring dots, spectral leakage): "
                r"$g^{(2)} = 1 - \rho^2\,(1 - g^{(2)}_{\mathrm{sig}})$ with $\rho = 1 - f_{\mathrm{bg}}$, the same law as above_band. Use this instead of "
                r"the setup leak $f_{\mathrm{leak}}$ (which is a laser, not a background). HOM and brightness stay signal-only. "
                r"The background-corrected convention $g^{(2)}(0)$:bg_corrected recovers $g^{(2)}_{\mathrm{sig}}$."),
    ]


def _electrical_cw_params():
    return [
        _P("k_esc_on_ns", "1/ns", 0.0, lo=0.0, hi=50.0, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Photonless tunnelling escape rate $k_{\mathrm{esc}}^{\mathrm{on}}$, X $\to$ G, at the cw bias."),
        _P("irf_ps", "ps", 0.0, lo=0.0, hi=2000.0, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Detector timing response $\tau_{\mathrm{IRF}}$ (Gaussian FWHM); 0 = off. Adds the IRF-convolved dip $g^{(2)}_{\mathrm{IRF}}(0)$."),
    ]


def _electrical_params(choices, default_mech, eta_default, pulsed):
    ps = [
        _P("mechanism", "", default_mech, choices=choices, tag="A", source="injection-mechanism roster [A]",
           sweepable=False, help="Injection mechanism (sets pump noise, loading and capture)."),
    ]
    if not pulsed:
        # R2 F10b: mu is the cw operating point. Every pulsed mechanism ignores it (M-3 computes mu = eta I tau / e from the
        # current, the M-4 turnstile fixes mu = 1), so the pulsed kind does not show it (old saved specs still load: removed.py).
        ps.append(_P("mu", "1", 0.5, lo=0.0, hi=20.0, tag="A", source=_SRC_DM,
                     visible_if={"mechanism": [m for m in choices if m not in _MU_IGNORED]},
                     help=r"Operating point $\mu$: mean carriers loaded per window."))
    ps += [
        _P("eta_capture", "1", eta_default, lo=0.0, hi=1.0, scale="lin", tag="A", source=_SRC_A,
           help=r"Delivery / capture efficiency $\eta_{\mathrm{cap}}$ into the dot."),
        _P("T_junction_K", "K", 300.0, lo=1.0, hi=500.0, tag="A", source=_SRC_DM,
           visible_if={"mechanism": [m for m in choices if _CARD_NAME[m] in ("m4a", "m4b", "m4c")]},
           help=r"Junction / island temperature $T_{\mathrm j}$ used to price a Coulomb-blockade element."),
    ]
    if pulsed:
        ps.append(_P("f_rep_MHz", "MHz", 80.0, lo=0.01, hi=2000.0, scale="log", tag="A", source=_SRC_A,
                     help="Electrical pulse repetition rate."))
        ps += _electrical_shape_params()
    else:
        ps += _electrical_cw_params()
    ps += _electrical_cascade_params()
    spec = {
        "I_uA": _P("I_uA", "uA", 10.0, lo=0.0, hi=1e4, tag="A", source=_SRC_DM, help=r"Drive current $I$."),
        "tau_on_ns": _P("tau_on_ns", "ns", 1.0, lo=0.0, hi=1e4, tag="A", source=_SRC_DM,
                        help=r"Current window length $\tau_{\mathrm{on}}$."),
        "F_p_rail": _P("F_p_rail", "1", 1.0, lo=0.0, hi=10.0, tag="A", source=_SRC_DM,
                       help=r"Pump Fano factor $F_{\mathrm{rail}}$ of the rail (1 = shot-noise limited)."),
        "R_drive_ohm": _P("R_drive_ohm", "ohm", 0.0, lo=0.0, hi=1e7, tag="A", source=_SRC_DM,
                          help=r"Rail resistance $R_{\mathrm{drive}}$ for the rise time; 0 = not given (needs $C_{\mathrm{dep}}$ too)."),
        "C_dep_F": _P("C_dep_F", "F", 0.0, lo=0.0, hi=1e-9, tag="A", source=_SRC_DM,
                      help=r"Depletion capacitance $C_{\mathrm{dep}}$ for the rise time; 0 = not given."),
        "Z_drive_ohm": _P("Z_drive_ohm", "ohm", 1e5, lo=1.0, hi=1e9, scale="log", tag="A", source=_SRC_DM,
                          help=r"Drive impedance $Z_{\mathrm{drive}}$ of the quiet rail."),
        "Z_junction_ohm": _P("Z_junction_ohm", "ohm", 1e3, lo=1.0, hi=1e9, scale="log", tag="A",
                             source=_SRC_DM, help=r"Junction impedance $Z_{\mathrm j}$ (F9 impedance divider)."),
        "f_reg_Hz": _P("f_reg_Hz", "Hz", 0.0, lo=0.0, hi=1e12, tag="A", source=_SRC_DM,
                       help=r"Regulation bandwidth $f_{\mathrm{reg}}$; 0 = ideal (needs $f_{\mathrm s}$ too)."),
        "f_sample_Hz": _P("f_sample_Hz", "Hz", 0.0, lo=0.0, hi=1e12, tag="A", source=_SRC_DM,
                          help=r"Sampling rate $f_{\mathrm s}$ of the regulated rail; 0 = not given."),
        "eps_cycle": _P("eps_cycle", "1", 0.01, lo=0.0, hi=1.0, tag="E", source=_SRC_DM + " [E]",
                        help=r"Per-cycle error $\varepsilon_{\mathrm{cyc}}$ of a single-carrier SET source ($F_{\mathrm p} \approx \varepsilon_{\mathrm{cyc}}$)."),
        "radius_nm": _P("radius_nm", "nm", 20.0, lo=1.0, hi=500.0, tag="A", source=_SRC_DM,
                        help=r"Island radius $R$ setting the charging energy."),
        "R_T_ohm": _P("R_T_ohm", "ohm", 1e6, lo=1e3, hi=1e9, scale="log", tag="A", source=_SRC_DM,
                      help=r"Tunnel resistance $R_{\mathrm T}$ (needs $R_{\mathrm T} \gg R_{\mathrm Q}$)."),
        "eta_boost": _P("eta_boost", "1", 2.0, lo=1.0, hi=20.0, tag="A", source=_SRC_DM,
                        help=r"Capture-efficiency boost $\eta_{\mathrm{boost}}$ of resonant-tunnelling injection."),
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
    refs=("Electrical-injection mechanism class models [A]",),
    provides=frozenset({"electrical_pulsed"}),
    observables=OBS_ELECTRICAL_PULSED, planner=_plan_electrical_pulsed)

ELECTRICAL_CW = PlannedDrive(
    kind="electrical_cw", label="Electrical, cw",
    summary=_ELEC_SUMMARY, status="ready",
    params=_electrical_params(MECH_CW, "M-1 poisson-rail", 0.05, False),
    refs=("Electrical-injection mechanism class models [A]",),
    provides=frozenset({"electrical_cw"}),
    observables=OBS_ELECTRICAL_CW, planner=_plan_electrical_cw)


# ======================================================================
# stark_gated_resonant (devlib-p10e2): cw/periodic resonant laser, electrical Stark gate tunes the dot INTO resonance
# ======================================================================
_SRC_CAO = ("Cao et al., Nano Lett. 15, 2015 (electrically gated resonant emission): Rabi period 1.6 ns at 500 nW, "
            "300 ps gate at 200 MHz [V]; class defaults otherwise [A]")


def _plan_stark_gated(self, vals, em, rt):
    """Coherent rect Rabi window with a time-dependent detuning det_X(t) (plan.det_profile).

    Paper convention: Delta = omega_L - omega_X (positive = laser BLUE of the dot). det_off_meV is Delta OUTSIDE the
    gate (dot detuned by the bias), det_on_meV is Delta INSIDE the gate (0 = on resonance); the plan stores the repo
    det_X = omega_X - omega_L = -Delta (det_profile = (-det_off, -det_on, t_gate_ps, t_r_ps)). The laser is on
    throughout; the Rabi window of the solver is t_gate + 6 t_r (RC edges included), the rest of the period is dark
    (off-gate excitation ~ Omega^2 Gamma_2/(4 det_off^2) is negligible for |det_off| >> hbar Gamma_2, brief 2.1)."""
    g = _getter(self, vals)
    d_off, d_on = float(g("det_off_meV")), float(g("det_on_meV"))
    t_gate = _positive("t_gate_ps", g("t_gate_ps"))
    t_r = _nonneg("t_r_ps", g("t_r_ps"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    period = _positive("rabi_period_ns", g("rabi_period_ns"))
    p_ref = _positive("power_ref_nW", g("power_ref_nW"))
    p = _nonneg("power_nW", g("power_nW"))
    omega = 2.0 * math.pi / period * math.sqrt(p / p_ref)          # Omega ~ sqrt(P)
    width = t_gate + 6.0 * t_r
    pulse = Pulse(shape="rect", fwhm_ps=width, area_rad=omega * width * 1e-3, detuning_meV=-d_off)
    inputs = {**_pulse_inputs("rect", width), "omega_ns": omega, "rabi_period_ns": period, "power_nW": p,
              "Delta_off_paper_meV": d_off, "Delta_on_paper_meV": d_on, "t_gate_ps": t_gate, "t_r_ps": t_r,
              "detuning_convention": "repo det_X = omega_X - omega_L = -Delta_paper; det_profile in repo convention"}
    jit = _nonneg("hom_jitter_ps", g("hom_jitter_ps")) * 1e-3
    hg = _nonneg("hom_gate_ps", g("hom_gate_ps"))
    if hg > 0.0:
        inputs["hom_gate"] = {"start_ns": _nonneg("hom_gate_start_ps", g("hom_gate_start_ps")) * 1e-3,
                              "width_ns": hg * 1e-3, "jitter_ns": jit, "shape": g("hom_gate_shape")}
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(pulse,),
                     needs=frozenset({"td_detuning"}), inputs=inputs,
                     det_profile=(-d_off, -d_on, t_gate, t_r))


STARK_GATED_RESONANT = PlannedDrive(
    kind="stark_gated_resonant", label="Stark-gated resonant laser (electrical gate)",
    summary="Resonant laser held at a fixed frequency; an electrical gate (QCSE) tunes the dot into resonance for "
            "t_gate once per period (Cao-type). Time-dependent detuning (off, on, gate width, RC edge) and the Rabi "
            "frequency from the measured period and power; optional detection-gate HOM (Bennett-type).",
    params=(
        _P("det_off_meV", "meV", 5.0, lo=-30.0, hi=30.0, tag="A", source=_SRC_CAO,
           help=r"Laser detuning $\Delta_{\mathrm{off}}$ OUTSIDE the gate, paper convention $\Delta = \omega_{\mathrm L} - \omega_X$ (positive = laser BLUE of "
                r"the dot). Must be $\gg \hbar\Gamma_2$ so the off-gate excitation is negligible."),
        _P("det_on_meV", "meV", 0.0, lo=-30.0, hi=30.0, tag="A", source=_SRC_CAO,
           help=r"Laser detuning $\Delta_{\mathrm{on}}$ INSIDE the gate, same convention (0 = on resonance)."),
        _P("t_gate_ps", "ps", 300.0, lo=1.0, hi=1e4, scale="log", tag="A", source=_SRC_CAO,
           help=r"Width $t_{\mathrm{gate}}$ of the electrical gate (boxcar) during which the dot is tuned to $\Delta_{\mathrm{on}}$."),
        _P("t_r_ps", "ps", 0.0, lo=0.0, hi=2000.0, tag="A", source=_SRC_CAO,
           help=r"RC time constant $t_{\mathrm r}$ of the gate edges (0 = sharp). The Rabi window is $t_{\mathrm{gate}} + 6\,t_{\mathrm r}$."),
        _P("rabi_period_ns", "ns", 1.6, lo=0.05, hi=100.0, scale="log", tag="A", source=_SRC_CAO,
           help=r"Period $T_\Omega = 2\pi/\Omega$ of the measured Rabi oscillation at $P_{\mathrm{ref}}$."),
        _P("power_ref_nW", "nW", 500.0, lo=0.01, hi=1e6, scale="log", tag="A", source=_SRC_CAO,
           help=r"Laser power $P_{\mathrm{ref}}$ at which $T_\Omega$ was measured."),
        _P("power_nW", "nW", 500.0, lo=0.0, hi=1e6, tag="A", source=_SRC_CAO,
           help=r"Laser power $P$ used; the Rabi frequency scales as $\sqrt{P/P_{\mathrm{ref}}}$."),
        _P("f_rep_MHz", "MHz", 200.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_CAO,
           help="Gate repetition rate."),
        _P("hom_gate_ps", "ps", 0.0, lo=0.0, hi=1e4, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Detection-gate width $T_{\mathrm{HOM}}$ of a gated HOM; 0 = off. Observable $I_{\mathrm{HOM}}^{\mathrm{gate}}$."),
        _P("hom_gate_shape", "", "rect", choices=("rect", "gauss"), tag="A", source=_SRC_P10E, sweepable=False,
           group="advanced",
           help=r"Shape of the HOM detection gate: rect (top-hat, opens at $t_{\mathrm{HOM}}$) or gauss (INTENSITY FWHM "
                r"$T_{\mathrm{HOM}}$ centred at $t_{\mathrm{HOM}}$)."),
        _P("hom_gate_start_ps", "ps", 0.0, lo=0.0, hi=1e4, tag="A", source=_SRC_P10E, group="advanced",
           help="Time after the excitation at which the HOM gate opens (rect) / its centre (gauss)."),
        _P("hom_jitter_ps", "ps", 0.0, lo=0.0, hi=1e3, tag="A", source=_SRC_P10E, group="advanced",
           help=r"Standard deviation $\sigma_{\mathrm{jit}}$ of the relative timing jitter between the two photons' detection gates inside the HOM gate "
                r"(Gaussian); the overlap falls monotonically with it."),
    ),
    refs=("Cao et al., Nano Lett. 15, 2015 (gated resonance); Bennett et al. (gated emission)",),
    provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED | {"g2_binned", "hom_I_gated", "hom_gate_fraction", "spectral_overlap"}, planner=_plan_stark_gated)


# ======================================================================
# needs_physics placeholders (full ParamSpecs so the UI can show them greyed)
# ======================================================================
_E_BRIEF = "Class range [E]"


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

DRIVES = (RESONANT_PULSE, TPE, CW_RESONANT, ABOVE_BAND, P_SHELL, ELECTRICAL_PULSED, ELECTRICAL_CW, STARK_GATED_RESONANT,
          PHONON_ASSISTED, SUPER, DICHROMATIC, ARP)

for _d in DRIVES:
    register(_d)
