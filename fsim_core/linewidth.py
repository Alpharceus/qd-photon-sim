"""Module C' -- literature-anchored linewidth model for the 300 K regime
(RT edge-emitter tier, 2026-09-02).

WHY. The class-proxy Gamma(T) that `device.evaluate` inherits from the V-a
joint fit (out/phase0/fit_params.json) was fitted to one (211)B InAs/GaAs
dot over 78-230 K. Run out to 300 K it gives 10.9 meV with an acoustic term
fitted to zero and an "LO" energy of 18 meV that matches no phonon of GaAs
(36 meV) or InP (43 meV) (code audit 2026-09-02, section 4.1). Nothing in
that fit knows about the InP-dot platform. This module provides a
transparent, material-aware alternative with explicit literature anchors and
an honest [E] range, so that every 300 K conclusion can be swept over the
range instead of resting on one extrapolated fit.

MODEL (standard exciton-phonon dephasing form, e.g. Bayer & Forchel, Phys.
Rev. B 65, 041308 (2002); Ortner et al., Phys. Rev. B 70, 201301 (2004);
Rudin, Reinecke, Segall, Phys. Rev. B 42, 11218 (1990)):

    Gamma(T) = Gamma_0 + a_ac * T + b_LO / (exp(E_LO / k_B T) - 1)

  Gamma_0  zero-temperature width (radiative + spectral diffusion) [meV]
  a_ac     acoustic-phonon (LA/TA) linear coefficient [meV/K]; single-dot
           literature 1-5 ueV/K for InGaAs-class dots (Ortner 2004, Bayer
           2002), larger for stronger confinement (GaN dots ~30 ueV/K,
           Holmes-class).
  b_LO     LO-phonon coupling [meV]; E_LO the LO phonon energy of the dot
           material: InP 43 meV, GaAs 36 meV, InAs 30 meV (Ioffe NSM tables).

ANCHORS for the 300 K width of III-V dots (FWHM, homogeneous or
single-dot):
  * In0.5Ga0.5As/GaAs single dot, NSOM: ~12 meV at 300 K -- Matsuda et al.,
    Phys. Rev. B 63, 121304(R) (2001) [V].
  * InAsP/InP nanowire dot: the 25 nm / ~16 meV figure at 1398 nm is the
    DETECTION BANDPASS FILTER width at 300 K, not a measured emission
    linewidth -- an instrument upper bound only [E] (Laferriere et al.,
    Nano Lett. 23, 962 (2023), p. 965). Not used as an anchor; see
    verify/data/rt_edge_anchors.yaml 'laferriere23-linewidth-class-proxy'.
  * InAs/InP single dot at 1.55 um: 17.5 meV homogeneous at 300 K (literature
    report cited in _goal/materials_research.md, section 6) [DR].
  * (211)B InAs/GaAs cavity dot: 6-7 meV at >= 250 K -- Chatzarakis et al.,
    PRApplied 20, 034011 (2023) [V] (lowest reported value; strong
    confinement).
  * InP/(Al,Ga)InP dots: 0.25 meV at low T with strong broadening above
    80 K; NO published 300 K single-dot width (Bommer et al., JAP 110,
    063108 (2011)) -- the InP-dot 300 K width is therefore UNMEASURED [A]
    and is carried as the class range below. The class range mixes
    homogeneous and spectrometer-FWHM values, and its 6 meV floor is based
    on one strongly confined (211)B dot.

Hence the class range used for envelopes [A]: Gamma(300 K) in [6, 20] meV,
default anchor 12 meV (Matsuda). The X-XX splittings of the same dot
families are 3-7 meV (InP/GaInP 4-7, InAs/InP 3.5, (211)B InAs 4-13), so at
300 K the lines overlap unless the dot sits at the low-width, high-splitting
corner -- this is the quantitative reason every 300 K single-photon result
of the literature saturates at g2 ~ 0.5 (materials_research.md, section 6).

The function `gamma_anchor` parameterizes the model by its 300 K value (the
measurable quantity) instead of b_LO: given Gamma_0, a_ac and E_LO, b_LO is
solved so that Gamma(300 K) hits the anchor. This keeps the extrapolation
between the two anchor points (low-T slope, 300 K width) rather than beyond
them.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .spectral import KB, gamma_of_T

E_LO_MEV = {"InP": 43.0, "GaAs": 36.0, "InAs": 30.0, "GaP": 51.0, "InGaAs": 33.0,
            "GaInP": 45.0}   # GaInP: InP-like mode dominant; [E] two-mode alloy

GAMMA300_CLASS_RANGE = (6.0, 20.0)   # meV, see anchors above
GAMMA300_DEFAULT = 12.0              # Matsuda 2001 In0.5Ga0.5As dot


@dataclass
class LinewidthParams:
    gamma0: float = 0.25        # meV; InP/(Al,Ga)InP low-T single-dot width (Bommer 2011)
    a_ac: float = 2.0e-3        # meV/K; InGaAs-class 1-5 ueV/K (Ortner 2004)
    E_LO: float = 43.0          # meV; InP
    gamma300: float = GAMMA300_DEFAULT  # meV anchor at 300 K
    # [A] InP/(Al,Ga)InP has no published 300 K single-dot width.
    tag: str = "A"

    @property
    def b_lo(self) -> float:
        return b_lo_from_anchor(self.gamma0, self.a_ac, self.E_LO, self.gamma300)


def b_lo_from_anchor(gamma0, a_ac, E_LO, gamma300, T_anchor=300.0):
    """LO coupling that makes Gamma(T_anchor) = gamma300 for the given
    Gamma_0, a_ac, E_LO. Raises if the anchor is below Gamma_0 + a_ac T."""
    resid = gamma300 - gamma0 - a_ac * T_anchor
    if resid < 0:
        raise ValueError("gamma300 anchor below Gamma_0 + a_ac*T: no positive b_LO")
    return resid * np.expm1(E_LO / (KB * T_anchor))


def gamma_anchor(T, p: LinewidthParams):
    """Gamma(T) [meV] of the anchored model (array-safe)."""
    return gamma_of_T(T, p.gamma0, p.a_ac, p.b_lo, p.E_LO)


def gamma_envelope(T, lo=GAMMA300_CLASS_RANGE[0], hi=GAMMA300_CLASS_RANGE[1],
                   gamma0=0.25, a_ac=2.0e-3, E_LO=43.0):
    """(Gamma_lo(T), Gamma_hi(T)) for the class range of 300 K widths."""
    pl = LinewidthParams(gamma0, a_ac, E_LO, lo)
    ph = LinewidthParams(gamma0, a_ac, E_LO, hi)
    return gamma_anchor(T, pl), gamma_anchor(T, ph)


def zpl_fraction_note():
    return ("The width above is the total exciton line (ZPL + acoustic sidebands) "
            "as measured by a spectrometer. When dot.lineshape='ibm' is used, the "
            "IBM sidebands are computed explicitly and the ZPL width should be the "
            "narrower dephasing-only width; using both at full value double-counts "
            "the acoustic contribution (qd_gf.effective_gamma_zpl documents the split).")
