"""Purcell-modified rates of one (emitter, structure) pair (design 1.4, 2.3).

    Gamma = Gamma0 (F_bg + F_P L)          radiative rate into the photon channel   [DR]
    beta  = (beta0 F_bg + F_P L)/(F_bg + F_P L)                                      [DR]
    k_nr                                   photonless loss, kept separate (not Purcell-enhanced)

L is the overlap of the emitter line with the cavity Lorentzian
(spectral.cavity_transmission). The emitter line FWHM depends on Gamma, which
depends on L, so L is solved self-consistently with the linewidth

    w(L) = hbar (Gamma(L) + gamma*)           [ueV, hbar in ueV ns]
    L    = cavity_transmission(detuning, w(L), kappa)

by Brent's method on (0, 1] (unique root: T decreases in w, w increases in L),
exactly as resonant_source.source_rates (checked against it in
verify/verify_compose_engine.py check 6 at F_bg = 1, k_nr = 0). Broadband
structures (kappa None) have L = 1.

Units and conventions (physics brief section 2): rates in 1/ns, energies in meV,
hbar = 6.582e-4 meV ns. gamma* is the REPO pure-dephasing rate = 2 x the
textbook 1/T2*: the e-g coherence decays at (Gamma + k + gamma*)/2, so the HOM
overlap is Gamma/(Gamma + gamma*_repo). F_P is the ADDITIONAL rate ratio
(total = F_bg + F_P L).

Strong-coupling flag [DR]. The weak-coupling (Purcell) model needs
g < (kappa + Gamma_bg + gamma*_repo)/4, i.e. the textbook condition
g > (kappa + gamma_tot)/4 with gamma_tot = Gamma_bg + 2 gamma*_textbook. Here
Gamma_bg = Gamma0 F_bg is the emitter decay into everything EXCEPT the cavity
mode (the cavity-mediated part is what the coupling competes with, so the
Purcell-enhanced Gamma must not enter). Regime "strong" when the inequality is
reversed. The JC flag of brief section 3 is jc_recommended = 4 g/kappa > 0.5
(adiabatic-elimination error ~ (2g/kappa)^2).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from scipy.optimize import brentq

from .. import lindblad as lb
from ..spectral import cavity_transmission
from .types import EmitterModel, Rates, StructureResponse

GAMMA_STAR_CONVENTION = ("repo: coherence decays at (Gamma + k + gamma*)/2, gamma*_repo = 2 x textbook "
                         "1/T2*; I = Gamma/(Gamma + gamma*_repo)")


@dataclass(frozen=True)
class RatesEx(Rates):
    """Rates plus the cavity quantities the solvers need (an extension of
    types.Rates, which stays untouched). All angular rates in 1/ns."""
    jc_recommended: bool = False
    g_ns: float = 0.0                       # JC coupling (0 = no cavity mode)
    kappa_ns: Optional[float] = None
    detuning_ueV: float = 0.0               # emitter minus cavity, ueV
    cav_minus_dot_ns: float = 0.0           # omega_cavity - omega_dot, 1/ns
    Gamma_bg_ns: float = 0.0                # Gamma0 F_bg: decay not into the cavity mode
    Gamma0_ns: float = 0.0                  # bare radiative rate
    linewidth_ueV: float = 0.0              # emitter FWHM used for L, hbar (Gamma + gamma*)


def combine(em: EmitterModel, sr: StructureResponse, T_K: float) -> RatesEx:
    """EmitterModel x StructureResponse -> Purcell-modified rates (module doc)."""
    G0 = float(em.gamma0_ns)
    F_bg = float(sr.F_bg)
    F_P = float(sr.F_P)
    gs = float(em.gamma_star_ns)
    k = float(em.k_nr_ns)
    h = lb.HBAR_MEV_NS * 1e3                      # ueV ns

    has_cav = sr.kappa_meV is not None and sr.kappa_meV > 0.0
    kappa_ueV = float(sr.kappa_meV) * 1e3 if has_cav else None
    if sr.E_mode_meV is not None:
        detuning_ueV = (float(em.E_X_meV) - float(sr.E_mode_meV)) * 1e3
    else:
        detuning_ueV = 0.0

    if has_cav:
        def resid(L):
            w = h * (G0 * (F_bg + F_P * L) + gs)
            return float(cavity_transmission(detuning_ueV, w, kappa_ueV)) - L
        L = float(brentq(resid, 0.0, 1.0, xtol=1e-15, rtol=8.9e-16))
    else:
        L = 1.0

    FPL = F_P * L
    Gamma = G0 * (F_bg + FPL)
    beta = (float(sr.beta0) * F_bg + FPL) / (F_bg + FPL)
    Gamma_bg = G0 * F_bg

    kappa_ns = float(lb.mev_to_rate(sr.kappa_meV)) if has_cav else None
    if sr.g_meV is not None:
        g_ns = float(lb.mev_to_rate(sr.g_meV))
    elif has_cav and F_P > 0.0 and G0 > 0.0:
        g_ns = math.sqrt(F_P * kappa_ns * G0 / 4.0)          # F_P = 4 g^2/(kappa Gamma0) [DR]
    else:
        g_ns = 0.0

    notes = [f"gamma* convention: {GAMMA_STAR_CONVENTION}"]
    regime = "weak"
    jc = False
    if has_cav and g_ns > 0.0:
        crit = (kappa_ns + Gamma_bg + gs) / 4.0
        if g_ns > crit:
            regime = "strong"
            notes.append(f"strong coupling: g = {g_ns:.4g}/ns > (kappa + Gamma_bg + gamma*)/4 = {crit:.4g}/ns")
        jc = 4.0 * g_ns / kappa_ns > 0.5
        if jc:
            notes.append(f"4g/kappa = {4.0 * g_ns / kappa_ns:.3g} > 0.5: adiabatic elimination (Purcell) "
                         "is unreliable, use the jc solver for lifetimes and coherent pulses")
    cav_minus_dot = (-detuning_ueV * 1e-3) / lb.HBAR_MEV_NS
    return RatesEx(Gamma_ns=Gamma, gamma_star_ns=gs, k_nr_ns=k, L=L, F_P_eff=FPL, beta=beta,
                   regime=regime, notes=tuple(notes), jc_recommended=bool(jc), g_ns=g_ns,
                   kappa_ns=kappa_ns, detuning_ueV=detuning_ueV, cav_minus_dot_ns=cav_minus_dot,
                   Gamma_bg_ns=Gamma_bg, Gamma0_ns=G0, linewidth_ueV=h * (Gamma + gs))
