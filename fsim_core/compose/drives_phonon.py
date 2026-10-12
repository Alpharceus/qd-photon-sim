"""Phonon-dependent coherent drives of the composable engine (P9): phonon_assisted, super, dichromatic.

Imported at the END of drives.py (which supplies the helpers used here); import `drives`, never this module
first. Each drive only PLANS; the "polaron" solver (polaron.py) executes. Conventions as in drives.py:

  * user-facing detunings are the PAPER convention Delta = omega_L - omega_X (positive = laser BLUE of the
    dot); each Pulse stores the REPO det_X = omega_X - omega_L = -Delta (the solver does NOT flip again);
  * fwhm is the INTENSITY FWHM of Omega^2; areas are in units of pi;
  * plan.pulses[0] defines the rotating frame; the other pulses enter as Omega_k e^{-i(w_Lk - w_L1) t}
    (brief 2.3 SUPER row), so SUPER and dichromatic differ from phonon_assisted only in the pulse tuple.

Dichromatic [A: parametrisation]. Two spectral components at +-delta around the dot (paper Delta = +-delta),
Omega(t) = Omega_env(t) [sqrt(a/(1+a)) e^{-i delta t} + sqrt(1/(1+a)) e^{+i delta t}] with a = I_blue/I_red the
intensity asymmetry ratio; a = 1 is sqrt(2) Omega_env cos(delta t). The envelope area is split as
area_pi sqrt(a/(1+a)) and area_pi sqrt(1/(1+a)) so the time-averaged intensity of the pair is that of the
single envelope. Closed check (verify_compose_polaron): a = 1 with phonons off gives P_X ~ 0.
"""
from __future__ import annotations

import math

from .drives import (OBS_COHERENT_PULSED, _E_BRIEF, _GAUSS, _PAPER_CONV, _REF_BRIEF, _SRC_A, PlannedDrive, _P,
                     _getter, _nonneg, _positive, _pulse_inputs)
from .types import DrivePlan, Pulse

_PHONON_NOTE = ("Needs the polaron solver (tier P1 dressed-basis phonon rates, Omega_r = sqrt(Z) Omega); "
                "with a material whose alpha is 0 it reduces to the Lindblad pulse.")


# ======================================================================
# phonon_assisted
# ======================================================================
def _plan_phonon_assisted(self, vals, em, rt):
    g = _getter(self, vals)
    shape = g("shape")
    fwhm = _positive("fwhm_ps", g("fwhm_ps"))
    area_pi = _nonneg("area_pi", g("area_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    delta = float(g("detuning_meV"))
    pulse = Pulse(shape=shape, fwhm_ps=fwhm, area_rad=area_pi * math.pi, detuning_meV=-delta)
    inputs = {**_pulse_inputs(shape, fwhm), "area_pi": area_pi, "Delta_paper_meV": delta,
              "det_X_meV": -delta, "detuning_convention": "repo det_X = omega_X - omega_L = -Delta_paper",
              "mechanism": "LA-phonon-assisted (blue detuned: emission into the lower dressed state)"}
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(pulse,),
                     needs=frozenset({"phonons"}), inputs=inputs)


PHONON_ASSISTED = PlannedDrive(
    kind="phonon_assisted", label="Phonon-assisted excitation",
    summary="Detuned strong pulse driving the dot through LA-phonon relaxation between dressed states "
            "(blue detuning works at low T by phonon emission; red needs absorption). " + _PHONON_NOTE,
    params=(
        _GAUSS,
        _P("detuning_meV", "meV", 1.0, lo=-3.0, hi=3.0, tag="E", source=_E_BRIEF + " (blue, paper convention)",
           help=_PAPER_CONV + " Class values: +0.8 to +1.3 meV (blue); a red value is allowed so the "
                              "emission/absorption asymmetry can be explored."),
        _P("area_pi", "pi", 10.0, lo=5.0, hi=20.0, tag="E", source=_E_BRIEF, help="Pulse area in units of pi."),
        _P("fwhm_ps", "ps", 15.0, lo=10.0, hi=20.0, tag="E", source=_E_BRIEF, help="Intensity FWHM."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
    ),
    refs=(_REF_BRIEF + ": LA-phonon-assisted row (Quilter et al., PRL 114, 137401 (2015); Thomas et al., "
                       "PRL 126, 233601 (2021); Cosacchi et al. 2019 as cited there)",
          "Ramsay et al., PRL 104, 017402 (2010); Nazir & McCutcheon, J. Phys. Condens. Matter 28, 103002 (2016)"),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_phonon_assisted)


# ======================================================================
# super
# ======================================================================
def _plan_super(self, vals, em, rt):
    g = _getter(self, vals)
    shape = g("shape")
    d1, d2 = float(g("delta1_meV")), float(g("delta2_meV"))
    w1, w2 = _positive("fwhm1_ps", g("fwhm1_ps")), _positive("fwhm2_ps", g("fwhm2_ps"))
    a1, a2 = _nonneg("area1_pi", g("area1_pi")), _nonneg("area2_pi", g("area2_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    p1 = Pulse(shape=shape, fwhm_ps=w1, area_rad=a1 * math.pi, detuning_meV=-d1)
    p2 = Pulse(shape=shape, fwhm_ps=w2, area_rad=a2 * math.pi, detuning_meV=-d2)
    inputs = {"shape": shape, "fwhm1_ps": w1, "fwhm2_ps": w2, "area1_pi": a1, "area2_pi": a2,
              "Delta_paper_meV": d1, "Delta2_paper_meV": d2, "det_X_meV": -d1,
              "detuning_convention": "repo det_X = -Delta_paper; frame of pulse 1, pulse 2 beats at "
                                     "(Delta2 - Delta1)/hbar",
              "swing_up_splitting_meV": abs(d2 - d1)}
    if shape == "gauss":
        inputs["sigma1_ps"] = _pulse_inputs(shape, w1)["sigma_ps"]
        inputs["sigma2_ps"] = _pulse_inputs(shape, w2)["sigma_ps"]
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(p1, p2),
                     needs=frozenset({"phonons", "two_colour"}), inputs=inputs)


SUPER = PlannedDrive(
    kind="super", label="SUPER (two-colour red-detuned swing-up)",
    summary="Two red-detuned pulses whose frequency difference swings the population up "
            "(swing-up of the quantum emitter population). " + _PHONON_NOTE,
    params=(
        _GAUSS,
        _P("delta1_meV", "meV", -3.5, lo=-6.0, hi=-1.0, tag="E",
           source="Piccinini et al. value as quoted in physics brief 2.3 (not re-checked) [E]",
           help="Detuning of pulse 1, paper convention Delta = omega_L - omega_X (red = negative)."),
        _P("delta2_meV", "meV", -7.5, lo=-12.0, hi=-2.0, tag="E",
           source="Piccinini et al. value as quoted in physics brief 2.3 (not re-checked) [E]",
           help="Detuning of pulse 2, paper convention."),
        _P("fwhm1_ps", "ps", 5.0, lo=3.0, hi=10.0, tag="E", source=_E_BRIEF, help="Intensity FWHM of pulse 1."),
        _P("fwhm2_ps", "ps", 5.0, lo=3.0, hi=10.0, tag="E", source=_E_BRIEF, help="Intensity FWHM of pulse 2."),
        _P("area1_pi", "pi", 5.0, lo=0.0, hi=20.0, tag="A", source=_SRC_A, help="Area of pulse 1 in units of pi."),
        _P("area2_pi", "pi", 5.0, lo=0.0, hi=20.0, tag="A", source=_SRC_A, help="Area of pulse 2 in units of pi."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
    ),
    refs=(_REF_BRIEF + ": SUPER row (Bracht et al., PRX Quantum 2, 040354 (2021) as cited there)",
          "Piccinini et al., arXiv:2510.21428 (2025)"),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_super)


# ======================================================================
# dichromatic
# ======================================================================
def _plan_dichromatic(self, vals, em, rt):
    g = _getter(self, vals)
    shape = g("shape")
    delta = _nonneg("delta_meV", g("delta_meV"))
    a = _nonneg("asymmetry", g("asymmetry"))
    fwhm = _positive("fwhm_ps", g("fwhm_ps"))
    area_pi = _nonneg("area_pi", g("area_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    s_blue = math.sqrt(a / (1.0 + a))
    s_red = math.sqrt(1.0 / (1.0 + a))
    # blue component: paper Delta = +delta -> repo det_X = -delta (pulse 0 = the frame); red: det_X = +delta
    blue = Pulse(shape=shape, fwhm_ps=fwhm, area_rad=area_pi * s_blue * math.pi, detuning_meV=-delta)
    red = Pulse(shape=shape, fwhm_ps=fwhm, area_rad=area_pi * s_red * math.pi, detuning_meV=+delta)
    inputs = {**_pulse_inputs(shape, fwhm), "area_pi": area_pi, "asymmetry": a, "delta_meV": delta,
              "Delta_paper_meV": delta, "det_X_meV": -delta,
              "detuning_convention": "components at paper Delta = +-delta_meV around the dot; frame = blue "
                                     "component (repo det_X = -delta)",
              "component_area_pi": (area_pi * s_blue, area_pi * s_red)}
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(blue, red),
                     needs=frozenset({"phonons", "two_colour"}), inputs=inputs)


DICHROMATIC = PlannedDrive(
    kind="dichromatic", label="Dichromatic excitation",
    summary="Two spectral components at +-Delta around the line (Omega(t) ~ Omega_0(t) cos(Delta t) for equal "
            "intensities); an ideal two-level system returns to the ground state, phonons and an intensity "
            "asymmetry create the population. " + _PHONON_NOTE,
    params=(
        _GAUSS,
        _P("delta_meV", "meV", 1.0, lo=0.2, hi=5.0, tag="A", source=_SRC_A,
           help="Half spacing of the two spectral components."),
        _P("asymmetry", "1", 1.0, lo=0.0, hi=10.0, tag="A", source=_SRC_A,
           help="Intensity asymmetry ratio a = I_blue / I_red of the two components (1 = symmetric)."),
        _P("fwhm_ps", "ps", 10.0, lo=1.0, hi=100.0, scale="log", tag="A", source=_SRC_A, help="Intensity FWHM."),
        _P("area_pi", "pi", 5.0, lo=0.0, hi=20.0, tag="A", source=_SRC_A,
           help="Pulse area of the common envelope in units of pi (split between the components so the "
                "pair carries the intensity of one envelope)."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
    ),
    refs=(_REF_BRIEF + ": dichromatic row (Koong et al. 2021 as cited there)",),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_dichromatic)

DRIVES_PHONON = (PHONON_ASSISTED, SUPER, DICHROMATIC)
