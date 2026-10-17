"""Purcell-modified rates of one (emitter, structure) pair (design 1.4, 2.3).

    Gamma = Gamma0 (F_bg + F_P L)          radiative rate into the photon channel   [DR]
    beta  = (beta0 F_bg + F_P L)/(F_bg + F_P L)                                      [DR]
    k_nr                                   photonless loss, kept separate (not Purcell-enhanced)

L is the overlap of the DOT line with the cavity Lorentzian
(spectral.cavity_transmission):

    w_dot = hbar (Gamma_bg + k_nr + gamma*_repo)       [ueV, hbar in ueV ns]
    L     = cavity_transmission(detuning, w_dot, kappa)

with Gamma_bg = Gamma0 F_bg the decay into everything EXCEPT the cavity mode
[DR]. Adiabatic elimination of the cavity gives the cavity-channel rate
R = 4 g^2/(kappa + gamma_d) x Lorentzian(Delta; kappa + gamma_d) with
gamma_d = Gamma_bg + k_nr + gamma*_repo the dot's own dephasing/decay width;
the Purcell-enhanced rate is the OUTPUT of that channel, not part of the
width the cavity filter sees, so there is no self-consistent root solve
[DR; review finding 2]. (Tomm 2021: L 0.89 -> 0.987, tau 54.2 -> ~49 ps with
the old Purcell-broadened width, measured 47.5 ps.) The same expression is
used by resonant_source.source_rates. Broadband structures (kappa None)
have L = 1.

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
Purcell-enhanced Gamma must not enter) and the photonless loss k_nr is included
like everywhere else (g > (kappa + Gamma_bg + k_nr + gamma*_repo)/4). Regime
"strong" when the inequality is reversed. The JC flag of brief section 3 is jc_recommended = 4 g/kappa > 0.5
(adiabatic-elimination error ~ (2g/kappa)^2).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace as _dc_replace
from typing import Optional

from .. import lindblad as lb
from ..spectral import cavity_transmission
from .types import EmitterModel, Rates, StructureResponse

GAMMA_STAR_CONVENTION = (r"repo $\gamma^{*}$ is the pure-dephasing linewidth $\Gamma_{\mathrm{ph}} = 2\gamma_{\mathrm{ph}} = 2/T_2^{*}$ "
                         r"(Scully-Zubairy $\gamma_{\mathrm{ph}}$ = coherence decay rate): the $e$-$g$ coherence decays at "
                         r"$(\Gamma + k_{\mathrm{nr}} + \Gamma_{\mathrm{ph}})/2$ and $I = \Gamma/(\Gamma + \Gamma_{\mathrm{ph}})$")


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
    linewidth_ueV: float = 0.0              # emitter FWHM hbar (Gamma + k_nr + gamma*) (Purcell-broadened, informational)
    dot_linewidth_ueV: float = 0.0          # width used for L: hbar (Gamma_bg + k_nr + gamma*), no Purcell channel
    beta0: float = 0.0                      # background-channel collected fraction (the structure's beta0)


def apply_pump(rt: RatesEx, P_ns: float) -> RatesEx:
    """Rates of the same emitter/structure with the cw incoherent pump P = r (1 + p) added to the dot line width
    (R2, ch08 eq. dp-lor: FWHM = Gamma + k + P + Gamma_ph):

        w_dot = hbar (Gamma_bg + k_nr + gamma*_repo + P)     ->  L, Gamma, beta recomputed.

    Needs the cavity data kept in `rt` (kappa, detuning, F_P L); without a cavity (L = 1) nothing depends on w_dot and
    `rt` is returned unchanged. A pure function of `rt` so the solvers can call it once the pump rate is known
    (the planner needs Gamma to set r, so combine() itself cannot include P)."""
    P = float(P_ns)
    if not (P > 0.0) or rt.kappa_ns is None or rt.kappa_ns <= 0.0 or rt.Gamma0_ns <= 0.0 or rt.L <= 0.0:
        return rt
    h = lb.HBAR_MEV_NS * 1e3
    G0, k, gs = float(rt.Gamma0_ns), float(rt.k_nr_ns), float(rt.gamma_star_ns)
    F_bg = float(rt.Gamma_bg_ns) / G0
    F_P = float(rt.F_P_eff) / float(rt.L)
    w_dot = h * (float(rt.Gamma_bg_ns) + k + gs + P)
    L = float(cavity_transmission(float(rt.detuning_ueV), w_dot, float(rt.kappa_ns) * h))
    FPL = F_P * L
    Gamma = G0 * (F_bg + FPL)
    beta = (float(rt.beta0) * F_bg + FPL) / (F_bg + FPL)
    return _dc_replace(rt, Gamma_ns=Gamma, L=L, F_P_eff=FPL, beta=beta, dot_linewidth_ueV=w_dot,
                       linewidth_ueV=h * (Gamma + k + gs + P),
                       notes=tuple(rt.notes) + (f"cw pump P = {P:.4g}/ns added to the dot line width seen by the cavity (L = {L:.4g})",))


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

    w_dot = h * (G0 * F_bg + k + gs)              # dot line width without the Purcell channel [DR]
    if has_cav:
        L = float(cavity_transmission(detuning_ueV, w_dot, kappa_ueV))
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
        crit = (kappa_ns + Gamma_bg + k + gs) / 4.0
        if g_ns > crit:
            regime = "strong"
            notes.append(f"strong coupling: g = {g_ns:.4g}/ns > (kappa + Gamma_bg + k_nr + gamma*)/4 = {crit:.4g}/ns")
        jc = 4.0 * g_ns / kappa_ns > 0.5
        if jc:
            notes.append(f"4g/kappa = {4.0 * g_ns / kappa_ns:.3g} > 0.5: adiabatic elimination (Purcell) "
                         "is unreliable, use the jc solver for lifetimes and coherent pulses")
    cav_minus_dot = (-detuning_ueV * 1e-3) / lb.HBAR_MEV_NS
    return RatesEx(Gamma_ns=Gamma, gamma_star_ns=gs, k_nr_ns=k, L=L, F_P_eff=FPL, beta=beta,
                   regime=regime, notes=tuple(notes), jc_recommended=bool(jc), g_ns=g_ns,
                   kappa_ns=kappa_ns, detuning_ueV=detuning_ueV, cav_minus_dot_ns=cav_minus_dot,
                   Gamma_bg_ns=Gamma_bg, Gamma0_ns=G0, linewidth_ueV=h * (Gamma + k + gs),
                   dot_linewidth_ueV=w_dot, beta0=float(sr.beta0))
