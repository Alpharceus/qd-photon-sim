"""Phonon-dependent coherent drives of the composable engine (P9, FX2c): phonon_assisted, super, dichromatic, arp.

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

ARP (FX2c, adiabatic rapid passage). A transform-limited Gaussian (intensity FWHM fwhm_ps, amplitude sigma
s0 = fwhm/(2 sqrt ln 2), field exp(-t^2 / 2 s0^2)) is given the spectral phase exp(+i phi2 w^2 / 2) (group-delay
dispersion phi2 = `chirp_ps2`, ps^2, field convention exp(-i w t)). The Gaussian transform gives [DR]
    E(t) ~ exp(-t^2 / (2 (s0^2 - i phi2)))  ->  sigma_c^2 = (s0^4 + phi2^2) / s0^2,
    carrier phase  -phi2 t^2 / (2 (s0^4 + phi2^2)),   instantaneous frequency  w_L + kappa t,
    kappa = phi2 / (s0^4 + phi2^2)  (rad/ps^2; the Pulse chirp is hbar kappa in meV/ps).
So chirp_ps2 > 0 is the POSITIVE chirp of Wei et al., Phys. Rev. Lett. 113, 097401 (2014) ("the laser sweeps from
low to high frequency", d(omega_L)/dt > 0, the QD follows the lower dressed state |->); negative = the reverse.
Stretching conserves the pulse energy, so `area_pi` is the area of the TRANSFORM-LIMITED pulse and the chirped
envelope has area area_pi sqrt(sigma_c/s0) (peak amplitude A0 sqrt(s0/sigma_c), Int |Omega|^2 dt unchanged) [DR].
detuning_meV is the paper convention Delta = omega_L(t=0) - omega_X (centre of the sweep). Needs {phonons,
td_detuning}: the polaron solver runs it (alpha = 0 reduces to the chirped Lindblad pulse).
"""
from __future__ import annotations

import math

from .drives import (OBS_COHERENT_PULSED, _E_BRIEF, _EID, _GAUSS, _PAPER_CONV, _SRC_A, _TRION, PlannedDrive,
                     _P, _dark_inputs, _getter, _nonneg, _positive, _pulse_inputs)
from .. import qd_gf
from .slm_pulse import SlmField
from ..resonant_source import gaussian_sigma
from .types import DrivePlan, Pulse

_SRC_KOONG = ("Koong et al., Phys. Rev. Lett. 126, 047403 (2021) (arXiv:2009.02121): sidebands +-0.6 meV, 0.5 meV wide (3.65 ps [DR]), "
              "contrast C = -0.65, first maximum at 5 pi [V]; asymmetry a = (1 + C)/(1 - C) [DR]")
_PHONON_NOTE = (r"Needs the polaron solver (tier P1 dressed-basis phonon rates, $\Omega_{\mathrm r} = \sqrt{Z}\,\Omega$); "
                r"with a material whose $\alpha$ is 0 it reduces to the Lindblad pulse.")


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
    _dark_inputs(g, inputs)
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
        _P("area_pi", "pi", 10.0, lo=5.0, hi=20.0, tag="E", source=_E_BRIEF, help=r"Pulse area $\Theta$ in units of $\pi$."),
        _P("fwhm_ps", "ps", 15.0, lo=10.0, hi=20.0, tag="E", source=_E_BRIEF, help=r"Intensity FWHM $\tau_{\mathrm{FWHM}}$."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _EID, _TRION,
    ),
    refs=("LA-phonon-assisted excitation: Quilter et al., PRL 114, 137401 (2015); Thomas et al., "
                       "PRL 126, 233601 (2021); Cosacchi et al. 2019",
          "Ramsay et al., PRL 104, 017402 (2010); Nazir & McCutcheon, J. Phys. Condens. Matter 28, 103002 (2016)"),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_phonon_assisted)


# ======================================================================
# super
# ======================================================================
def _plan_super_slm(self, g, d1, d2, f_rep):
    """R2 SUPCMP: the pulse-shaper pair of Karli et al. 2022 Sec. V B (slm_pulse.py): a broadband Gaussian spectrum times a
    two-Gaussian field mask (with a floor between the peaks), inverse-transformed to Omega(t). The two Pulse entries only carry
    the frame (pulse 1), the window (the transform-limited width of the mask) and an area estimate; the field is in
    plan.inputs['slm_field'] (plain floats) and replaces the envelopes inside the polaron PulseTrain."""
    sig_m = _positive("slm_mask_sigma_meV", g("slm_mask_sigma_meV"))
    fld = {"d1_meV": d1, "d2_meV": d2, "T1": _nonneg("slm_T1", g("slm_T1")), "T2": _nonneg("slm_T2", g("slm_T2")),
           "alpha_rad": _positive("slm_alpha_pi", g("slm_alpha_pi")) * math.pi,
           "bb_fwhm_meV": _positive("slm_bb_fwhm_meV", g("slm_bb_fwhm_meV")), "bb_centre_meV": float(g("slm_bb_centre_meV")),
           "mask_sigma_meV": sig_m, "floor": _nonneg("slm_floor", g("slm_floor"))}
    if fld["T1"] > 1.0 or fld["T2"] > 1.0:
        raise ValueError("slm_T1 / slm_T2 are intensity transmissions in [0, 1]")
    sf = SlmField(**fld)
    fw = SlmField.mask_fwhm_ps(sig_m)
    p1 = Pulse(shape="gauss", fwhm_ps=fw, area_rad=math.sqrt(fld["T1"]) * sf.line_area_rad(1), detuning_meV=-d1)
    p2 = Pulse(shape="gauss", fwhm_ps=fw, area_rad=math.sqrt(fld["T2"]) * sf.line_area_rad(2), detuning_meV=-d2)
    inputs = {"shape": "slm_mask", "slm_field": fld, "Delta_paper_meV": d1, "Delta2_paper_meV": d2, "det_X_meV": -d1,
              "detuning_convention": "repo det_X = -Delta_paper; frame of pulse 1, mask line 2 beats at (Delta2 - Delta1)/hbar",
              "swing_up_splitting_meV": abs(d2 - d1), "mask_fwhm_ps": fw,
              "pulse_model": "Karli et al. 2022 Sec. V B: broadband Gaussian x two-Gaussian field mask (+ floor), inverse FT"}
    _dark_inputs(g, inputs)
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(p1, p2),
                     needs=frozenset({"phonons", "two_colour"}), inputs=inputs)


def _plan_super(self, vals, em, rt):
    g = _getter(self, vals)
    shape = g("shape")
    d1, d2 = float(g("delta1_meV")), float(g("delta2_meV"))
    if shape == "slm_mask":
        return _plan_super_slm(self, g, d1, d2, _positive("f_rep_MHz", g("f_rep_MHz")))
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
    _dark_inputs(g, inputs)
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(p1, p2),
                     needs=frozenset({"phonons", "two_colour"}), inputs=inputs)


_SLM_SRC = "Karli et al., Nano Lett. 22, 6567 (2022) (arXiv:2203.00712v1), Sec. V B (Eqs. 2-3)"
_SUPER_SHAPE = _P("shape", "", "gauss", choices=("gauss", "rect", "slm_mask"), tag="A", source=_SRC_A, sweepable=False,
                  help="Pulse envelope. slm_mask: the pulse-shaper pair of Karli et al. 2022 Sec. V B (broadband Gaussian spectrum "
                       "times a two-Gaussian field mask, inverse-transformed); then only delta1/delta2 and the slm_* "
                       "parameters act (areas and widths follow from the spectrum).")
_SLM_VIS = {"shape": ["slm_mask"]}

SUPER = PlannedDrive(
    kind="super", label="SUPER (two-colour red-detuned swing-up)",
    summary="Two red-detuned pulses whose frequency difference swings the population up "
            "(swing-up of the quantum emitter population). " + _PHONON_NOTE,
    params=(
        _SUPER_SHAPE,
        _P("delta1_meV", "meV", -5.0, lo=-12.0, hi=-1.0, tag="V",
           source="Bracht et al., PRX Quantum 2, 040354 (2021), Table 1 row 2: Delta1 = -5.0 meV [V]",
           help=r"Detuning $\Delta_1$ of pulse 1, paper convention $\Delta = \omega_{\mathrm L} - \omega_X$ (red = negative). "
                r"Default: Bracht et al. 2021 Table 1 row 2. Piccinini et al. 2025 use -3.5 / -7.5 meV, which is their XX-mediated "
                r"F2 point; their two-level X point is -3.7 / -6.9 meV (Fig. 3) [V]."),
        _P("delta2_meV", "meV", -13.56, lo=-30.0, hi=-2.0, tag="V",
           source="Bracht et al., PRX Quantum 2, 040354 (2021), Table 1 row 2: Delta2 = -13.56 meV [V] (Eq. 6 value -13.5624)",
           help=r"Detuning $\Delta_2$ of pulse 2, paper convention. Default: Bracht et al. 2021 Table 1 row 2. The swing-up "
                r"condition is Bracht Eq. 6, $\Delta_2 = \Delta_1 - \sqrt{\Omega_1(0)^2 + \Delta_1^2}$ with $\Omega_1(0)$ the "
                r"peak Rabi energy of pulse 1."),
        _P("fwhm1_ps", "ps", 2.58, lo=1.0, hi=15.0, tag="DR",
           source="Bracht et al. 2021 Table 1 row 2: field width sigma = 1.55 ps [V]; intensity FWHM = 2 sigma sqrt(ln 2) = 2.581 ps [DR]",
           help=r"Intensity FWHM $\tau_{\mathrm{FWHM},1}$ of pulse 1 ($=2\sigma\sqrt{\ln 2}$ for the field width $\sigma$ Bracht quotes). "
                r"At fixed peak Rabi energy the Gaussian area scales as $\tau_{\mathrm{FWHM}}$ ($\Theta = \Omega_{\mathrm{pk}}\sigma\sqrt{2\pi}$)."),
        _P("fwhm2_ps", "ps", 2.58, lo=1.0, hi=15.0, tag="DR",
           source="Bracht et al. 2021 Table 1 row 2: sigma = 1.55 ps [V]; intensity FWHM 2.581 ps [DR]",
           help=r"Intensity FWHM $\tau_{\mathrm{FWHM},2}$ of pulse 2."),
        _P("area1_pi", "pi", 13.06, lo=0.0, hi=60.0, tag="V",
           source="Bracht et al. 2021 Table 1 row 2: alpha1 = 13.06 pi (area = int Omega dt, H = Omega/2 sigma_x) [V]",
           help=r"Area $\Theta_1$ of pulse 1 in units of $\pi$. SUPER needs large areas (the peak Rabi energy must reach the "
                r"detuning scale) and the inversion is a narrow resonance in both areas, about $0.5\pi$ wide: use Bracht Eq. 6 "
                r"to set $\Delta_2$ and expect to tune the areas to that precision (Bracht et al., PRX Quantum 2, 040354 (2021))."),
        _P("area2_pi", "pi", 10.45, lo=0.0, hi=60.0, tag="V",
           source="Bracht et al. 2021 Table 1 row 2: alpha2 = 10.45 pi [V]",
           help=r"Area $\Theta_2$ of pulse 2 in units of $\pi$ (see pulse 1)."),
        _P("slm_T1", "1", 0.5, lo=0.0, hi=1.0, tag="V", source=_SLM_SRC + ": T1 = 0.5 (pulse 1 fixed, Table II / Fig. 2)", group="advanced",
           visible_if=_SLM_VIS, help=r"INTENSITY transmission $T_1$ of mask line 1 ($\tilde T_1 = \sqrt{T_1}$ on the field). slm_mask only."),
        _P("slm_T2", "1", 0.7, lo=0.0, hi=1.0, tag="V", source=_SLM_SRC + ": T2 is the y-axis of Fig. 2 (0 to 1)", group="advanced",
           visible_if=_SLM_VIS, help=r"INTENSITY transmission $T_2$ of mask line 2 ($\tilde T_2 = \sqrt{T_2}$ on the field). slm_mask only."),
        _P("slm_alpha_pi", "pi", 470.0, lo=1.0, hi=5000.0, scale="log", tag="V", source=_SLM_SRC + ": integrated pulse area alpha~ = 470 pi", group="advanced",
           visible_if=_SLM_VIS, help=r"Integrated area $\tilde{\alpha}$ of the broadband spectrum (Eq. 2; $\int \tilde{\Omega}\,d\omega = \tilde{\alpha}$)."),
        _P("slm_bb_fwhm_meV", "meV", 17.0, lo=1.0, hi=100.0, tag="V", source=_SLM_SRC + ": spectral FWHM 17 meV", group="advanced",
           visible_if=_SLM_VIS, help="Spectral FWHM of the broadband Gaussian (Eq. 2)."),
        _P("slm_bb_centre_meV", "meV", -8.4, lo=-30.0, hi=10.0, tag="V", source=_SLM_SRC + ": centred at hbar(w_C - w0) = -8.4 meV", group="advanced",
           visible_if=_SLM_VIS, help="Centre of the broadband spectrum relative to the dot (paper convention)."),
        _P("slm_mask_sigma_meV", "meV", 0.2, lo=0.02, hi=2.0, tag="V", source=_SLM_SRC + ": hbar sigma = 0.2 meV for all calculations", group="advanced",
           visible_if=_SLM_VIS, help=r"Width $\hbar\sigma$ of each Gaussian FIELD mask line (Eq. 3)."),
        _P("slm_floor", "1", 0.05, lo=0.0, hi=1.0, tag="DR", source=_SLM_SRC + ": 'a transmission of 5% between the two peaks'; read as a field weight [DR]", group="advanced",
           visible_if=_SLM_VIS, help="Mask floor between the two peak centres (non-zero SLM extinction), added to the field mask."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _EID, _TRION,
    ),
    refs=("SUPER excitation: Bracht et al., PRX Quantum 2, 040354 (2021)",
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
    _dark_inputs(g, inputs)
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(blue, red),
                     needs=frozenset({"phonons", "two_colour"}), inputs=inputs)


DICHROMATIC = PlannedDrive(
    kind="dichromatic", label="Dichromatic excitation",
    summary="Two spectral components at +-Delta around the line (Omega(t) ~ Omega_0(t) cos(Delta t) for equal "
            "intensities); an ideal two-level system returns to the ground state, phonons and an intensity "
            "asymmetry create the population. " + _PHONON_NOTE,
    params=(
        _GAUSS,
        _P("delta_meV", "meV", 0.6, lo=0.2, hi=5.0, tag="V", source=_SRC_KOONG,
           help=r"Half spacing $\delta$ of the two spectral components. Default: Koong et al. 2021 (sidebands $\pm 0.6$ meV)."),
        _P("asymmetry", "1", 0.212, lo=0.0, hi=10.0, tag="DR", source=_SRC_KOONG,
           help=r"Intensity asymmetry ratio $a = I_{\mathrm{blue}}/I_{\mathrm{red}}$ of the two components. Default: Koong's optimum "
                r"contrast $C = -0.65$, $a = (1 + C)/(1 - C) = 0.212$ (a symmetric pair, $a = 1$, returns an ideal two-level system to "
                r"the ground state; only phonons and the asymmetry create the population)."),
        _P("fwhm_ps", "ps", 3.65, lo=1.0, hi=100.0, scale="log", tag="DR", source=_SRC_KOONG,
           help=r"Intensity FWHM $\tau_{\mathrm{FWHM}}$. Default: transform limit of the 0.5 meV sidebands of Koong et al. 2021 (3.65 ps)."),
        _P("area_pi", "pi", 5.0, lo=0.0, hi=20.0, tag="V", source=_SRC_KOONG,
           help=r"Pulse area $\Theta$ of the common envelope in units of $\pi$ (split between the components so the "
                r"pair carries the intensity of one envelope). Default: first oscillation maximum, Koong et al. 2021 Fig 3d."),
        _P("f_rep_MHz", "MHz", 80.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate."),
        _EID, _TRION,
    ),
    refs=("Dichromatic excitation: Koong et al. 2021",),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_dichromatic)

# ======================================================================
# arp
# ======================================================================
def arp_pulse_numbers(fwhm_ps, chirp_ps2):
    """(sigma_c_ps, fwhm_c_ps, kappa_rad_ps2, sigma0_ps) of the chirped Gaussian (module doc) [DR]."""
    s0 = float(gaussian_sigma(fwhm_ps))
    phi2 = float(chirp_ps2)
    den = s0 ** 4 + phi2 ** 2
    sc = math.sqrt(den) / s0
    return sc, sc * 2.0 * math.sqrt(math.log(2.0)), phi2 / den, s0


def _plan_arp(self, vals, em, rt):
    g = _getter(self, vals)
    fwhm0 = _positive("fwhm_ps", g("fwhm_ps"))
    area_pi = _nonneg("area_pi", g("area_pi"))
    f_rep = _positive("f_rep_MHz", g("f_rep_MHz"))
    phi2 = float(g("chirp_ps2"))
    if not math.isfinite(phi2):
        raise ValueError(f"chirp_ps2 must be finite, got {phi2!r}")
    delta = float(g("detuning_meV"))
    sc, fwhm_c, kappa, s0 = arp_pulse_numbers(fwhm0, phi2)
    chirp_meV_ps = qd_gf.HBAR * kappa                   # meV ps x rad/ps^2 = meV/ps
    area_c = area_pi * math.pi * math.sqrt(sc / s0)     # energy-conserving stretch of the TL area
    pulse = Pulse(shape="gauss", fwhm_ps=fwhm_c, area_rad=area_c, detuning_meV=-delta, chirp_meV_ps=chirp_meV_ps)
    inputs = {**_pulse_inputs("gauss", fwhm_c), "area_pi": area_pi, "area_chirped_pi": area_c / math.pi,
              "fwhm_tl_ps": fwhm0, "fwhm_chirped_ps": fwhm_c, "chirp_ps2": phi2, "chirp_rad_ps2": kappa,
              "chirp_meV_ps": chirp_meV_ps, "Delta_paper_meV": delta, "det_X_meV": -delta,
              "detuning_convention": "repo det_X = omega_X - omega_L = -Delta_paper (centre of the sweep)",
              "chirp_sign": "chirp_ps2 > 0: laser frequency rises in time (Wei 2014 positive chirp)",
              "mechanism": "adiabatic rapid passage along the lower (positive chirp) / upper (negative chirp) "
                           "dressed state"}
    _dark_inputs(g, inputs)
    return DrivePlan(family="coherent", timing="pulsed", f_rep_MHz=f_rep, target="X", pulses=(pulse,),
                     needs=frozenset({"phonons", "td_detuning"}), inputs=inputs)


ARP = PlannedDrive(
    kind="arp", label="Adiabatic rapid passage (chirped pulse)",
    summary="Linearly chirped Gaussian sweeping through the resonance: robust inversion above the adiabatic "
            "threshold area. Positive chirp (frequency rising in time) follows the lower dressed state and stays "
            "robust at low T (the phonon relaxation |+> -> |-> emits); negative chirp is degraded by it "
            "(Wei et al. 2014). " + _PHONON_NOTE,
    params=(
        _P("shape", "", "gauss", choices=("gauss",), tag="A", source=_SRC_A, sweepable=False,
           help="Transform-limited envelope before the chirp (Gaussian only)."),
        _P("chirp_ps2", "ps^2", 30.0, lo=-200.0, hi=200.0, tag="DR",
           source="Wei et al., PRL 113, 097401 (2014): 3 ps transform-limited stretched to 30 ps gives "
                  "phi2 ~ s0 sigma_c ~ 30 ps^2 (derived from those two durations, not read from the paper) [DR]",
           help=r"Spectral chirp (group-delay dispersion $\varphi_2$) of the transform-limited pulse. Positive = laser "
                r"frequency rises in time (Wei positive chirp); the swept rate is $\varphi_2/(s_0^4 + \varphi_2^2)$."),
        _P("detuning_meV", "meV", 0.0, lo=-3.0, hi=3.0, tag="A", source=_SRC_A,
           help=_PAPER_CONV + " Detuning of the carrier at the pulse centre (centre of the sweep)."),
        _P("area_pi", "pi", 2.0, lo=0.0, hi=10.0, tag="A", source="Wei 2014 plateau 1.5 pi - 3 pi [A]",
           help=r"Pulse area $\Theta$ of the transform-limited pulse in units of $\pi$ (chirping conserves the pulse "
                r"energy, so the chirped envelope area is $\Theta\sqrt{\sigma_{\mathrm c}/\sigma_0}$)."),
        _P("fwhm_ps", "ps", 3.0, lo=0.5, hi=20.0, scale="log", tag="E",
           source="Wei et al., PRL 113, 097401 (2014): 3 ps transform-limited pulses (before stretching) [E]",
           help=r"Intensity FWHM $\tau_{\mathrm{FWHM}}$ of the transform-limited pulse before the chirp."),
        _P("f_rep_MHz", "MHz", 76.0, lo=1.0, hi=2000.0, scale="log", tag="A", source=_SRC_A,
           help="Pulse repetition rate (mode-locked Ti:sapphire class)."),
        _EID, _TRION,
    ),
    refs=("Wei et al., Phys. Rev. Lett. 113, 097401 (2014): positively chirped pulses robust, negatively "
          "chirped degraded by acoustic-phonon emission",),
    status="ready", provides=frozenset({"coherent_pulse"}),
    observables=OBS_COHERENT_PULSED, planner=_plan_arp)

DRIVES_PHONON = (PHONON_ASSISTED, SUPER, DICHROMATIC, ARP)
