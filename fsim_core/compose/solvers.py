"""Solvers of the composable engine (design 2.5): rate, lindblad,
lindblad_incoherent, jc, and select().

Each solver wraps an existing module and maps (DrivePlan, EmitterModel,
RatesEx) onto it:

  rate                 pulse_counting.pulse_g2 (pulsed) / cw_g2.g2_cw (cw): INCOHERENT pumping
                       only (pulse_counting validity guard D11). No coherence -> I_single NaN.
  lindblad             resonant_source: TwoLevel + pulsed_photon_statistics +
                       pulsed_indistinguishability (any area, detuning, rect/gauss); Cascade +
                       tpe_photon_statistics + tpe_cascade_hom for TPE; cw resonance
                       fluorescence from lindblad.build_system(levels=2, omega) + g2_tau.
  lindblad_incoherent  S3d scheme (G, X, R): pump G->R during the window, R->X at 1/tau_cap,
                       X->G at Gamma (+k), dephasing gamma*. g2 from lindblad.pulsed_counting;
                       HOM from lindblad.indistinguishability started in |R> (instantaneous
                       loading, window re-excitation ignored [A]). Anchor (brief 2.3):
                       I = [g_r/(g_r + Gamma_t)] Gamma_t/(Gamma_t + gamma*) with g_r = 1/tau_cap.
                       FX2b extensions (all off by default = the scheme above, bit-for-bit):
                       * p_refill, tau_refill (trap refilling [A]): each radiative emission is followed with
                         probability p by a dark state W that recaptures a second, uncorrelated pair at
                         1/tau_refill (X -> W branch, W -> X); photons count both branches. Geometric chain:
                         n_ph = p_loaded/(1 - p), g2 = 2 p / p_loaded (exact for a negligible pump re-excitation).
                       * bg_frac (uncorrelated background [A]): g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac
                         (Brouri form); HOM and n_ph stay signal-only.
                       * k0_ns, E_act_meV, T_K (thermal escape / recapture, Arrhenius [E] class model, p_shell):
                         k_esc = k0 exp(-E_act/kT); X -> R escape (delays the emission; in the HOM generator) and
                         G -> R recapture of thermally released carriers between pulses (uncorrelated extra photons,
                         NOT in the HOM generator: a steady background stream would make N diverge). Small
                         lam = k_esc T_rep: g2 = (2 p1 lam + lam^2)/(p1 + lam)^2. Note: a conserved single pair
                         escaping and returning cannot by itself raise g2; the growth comes from the recaptured
                         extra carriers, which is why the recapture channel is part of the model.
  jc                   dot + cavity mode (Fock cutoff 3) via lindblad.build_system(cavity=...).
                       Lifetime (slowest radiative-sector Liouvillian eigenvalue), photons per excitation (exact, projected
                       Laplace solve), pulsed g2 / mean photons through resonant_source._stats and
                       cavity-output HOM through pulsed_indistinguishability.

eta_prep (FX2c): the preparation PROBABILITY (papers define it so) = P(n >= 1) with n = the number of excitations
produced by one pulse from the ground state = photons emitted inside the window + the excitation (X; XX for TPE;
X or cavity photon for jc) still present at the window end. It is computed DIRECTLY as 1 - p0, p0 = the weight of the
no-jump (no counted photon) evolution that ends unexcited, so 0 <= eta_prep <= 1 exactly (_window_p0). The
expected excitation COUNT (FX2b's eta_prep = P_X(end) + emitted, which exceeds 1 under re-excitation) is kept as
extra["n_exc"]; for n <= 2 the counting hierarchy gives the same probability as m1 - m2/2 with the window moments
(truncation: that identity drops P(n >= 3), checked in verify_compose_engine). Nonradiative loss inside the window
is not added back. The old end-of-window population is kept as extra["P_X_end"] (informational; it includes the
decay during the window).

devlib-p10e (electrical drives; all off by default = bit-identical to the rectangular paths):
  * SHAPED incoherent pump r(t) (plan.pump_profile rect | rc | gauss, plan.t_rise_ns): piecewise-constant
    cell-averaged steps of `pump_cells` (mu = int r dt is the operating point, preserved exactly), handed to
    pulse_counting.pulse_g2(r_profile=...) (rate) and lindblad.pulsed_counting(profile=...) (lindblad_incoherent).
    An RC tail LENGTHENS the effective window (tau_eff ~ tau_on + 1.7 t_r); it is never subtracted.
  * bias-dependent photonless escape: plan.k_esc_on_ns (pump window) and plan.k_esc_off_ns (dark window) act as
    X -> G losses (XX scales with gamma_xx_ratio like k_nr); the EL decay time tau_EL = 1/(Gamma + k_nr + k_off) is
    REPORTED (extra["tau_EL_ps"]) and never fed back into Gamma0.
  * td_detuning is supported by the lindblad solver (and polaron): plans with a Stark gate (plan.det_profile) or a
    chirped pulse run on the polaron step-propagator machinery with the phonons switched off.

devlib-fx3b (all off by default = bit-identical):
  * M-3 pulsed pump: mu = eta_capture N_c already contains the capture efficiency, so _interp_pump does not multiply
    eta_capture a second time (R2 D3 extends this to every mechanism: mu = pairs loaded into the dot).
  * electrical cascade (plan.needs {"cascade"}, plan.inputs["pump_ratio_xx"]) runs on the rate solver (G-X-XX with
    pulse_counting.pulse_g2 and the shaped pump): strong pumping parks the dot in XX and g2 saturates then falls.
  * plan.inputs["bg_frac"] (electrical drives too): g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac (as above_band).
  * plan.inputs["eid_ps"] (excitation-induced dephasing gamma_EID = c |Omega|^2) routes a pulsed two-level plan through
    the polaron step propagators with the phonons off; plan.inputs["p_trion"] scales the result by the per-pulse dark-
    trion probability p sin^2(area/2) (polaron.apply_trion; HOM unchanged, g2 / (1 - p)).

Conventions (physics brief 2): rates 1/ns; Pulse.detuning_meV holds the REPO det_X =
omega_X - omega_L (P4 drives.py stores -Delta_paper); pulse fwhm is the INTENSITY FWHM;
gamma* is the repo value (2 x textbook). Phonons: only Z-factor bookkeeping in
observables.py (Z^2 on unfiltered HOM, Z on filtered brightness, never both). A plan that needs
"phonons" goes to the "polaron" solver (polaron.py, P9, tier P1 dressed-basis rates), registered at the
bottom of this module.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np
from dataclasses import replace as _dc_replace

from scipy.linalg import expm, solve

from .. import cw_g2
from .. import lindblad as lb
from .. import pulse_counting as pcnt
from .. import resonant_source as rs
from ..spectral import KB
from . import rates as _rates_mod
from .registry import Solver, get as _get, kinds as _kinds, register
from .types import SolverResult

NAN = float("nan")
_UNSUPPORTED_NEEDS = {
    "phonons": "needs P9 polaron solver",
    "td_detuning": "time-dependent detuning (Stark gate / chirp) is supported only by the lindblad and polaron solvers",
    "two_colour": "two-colour drives are not implemented",
}
_HOM_TOL = 1e-6                    # [A] numerical choice for the Richardson-converged HOM integral
_HOM_DOUBLINGS = 5


class SolverUnavailable(ValueError):
    """No registered solver can run the plan; the message lists the reasons."""


# ---------------------------------------------------------------- helpers
def _bad_needs(plan, allowed=()):
    out = []
    for n in sorted(plan.needs):
        if n in allowed:
            continue
        out.append(_UNSUPPORTED_NEEDS.get(n, f"need {n!r} is not supported by this solver"))
    return out


def _det_X_meV(plan) -> float:
    """Repo det_X (omega_X - omega_L) in meV: pulse value (already det_X, P4) or the cw input."""
    if plan.pulses:
        return float(plan.pulses[0].detuning_meV)
    return float(plan.inputs.get("det_X_meV", 0.0))


def _det_X_ns(plan) -> float:
    return float(lb.mev_to_rate(_det_X_meV(plan)))


def _f_rep(plan) -> Optional[float]:
    return None if plan.f_rep_MHz is None else float(plan.f_rep_MHz)


def _window_ns(pulse) -> float:
    return rs.pulse_window(pulse.shape, pulse.fwhm_ps * 1e-3)


def _pulsed_coherent_check(plan) -> Optional[str]:
    if plan.timing != "pulsed":
        return None
    if len(plan.pulses) != 1:
        return f"exactly one pulse supported, plan has {len(plan.pulses)}"
    p = plan.pulses[0]
    if p.shape not in ("rect", "gauss"):
        return f"pulse shape {p.shape!r} not supported (rect|gauss)"
    if abs(p.chirp_meV_ps) > 0.0 and "td_detuning" not in plan.needs:
        return "chirped pulses need time-dependent detuning (needs td_detuning)"
    f = _f_rep(plan)
    if f is None or f <= 0:
        return "pulsed plan needs f_rep_MHz"
    if 1e3 / f <= _window_ns(p):
        return "repetition period must exceed the pulse window"
    return None


def _interp_pump(plan, rt):
    """(pump rate 1/ns, window ns) of an incoherent / electrical plan.
    Electrical plans carry only a DriveInterface: pulsed r = mu / tau_on (mean loadings per window, tau_on =
    interface tau_eff_ns, default 1 ns); cw r = mu (Gamma + k_nr) (mean loadings per X lifetime, the same definition
    above_band cw uses: P/P_sat x (Gamma + k)) [A].
    R2 D3: mu is the number of pairs loaded INTO THE DOT per window for EVERY mechanism (ch05 eq. de-mu, ch10 eq.
    in-rdot: r_dot already contains eta_inj f_QD). eta_capture does NOT rescale it; it only thins the pump Fano factor
    (F_dot = 1 - eta (1 - F_pump), see _pump_fano) and, for M-3 only, converts the current to mu = eta I tau / e
    (drive_mech.mech_pulsed). Before R2 the M-1/M-2/M-4/M-5 paths multiplied mu by eta_capture a second time."""
    if plan.pump_rate_ns > 0.0:
        return float(plan.pump_rate_ns), float(plan.tau_on_ns)
    itf = plan.interface
    if itf is not None:
        eff = float(itf.mu)
        if plan.timing == "pulsed":
            tau = float(plan.tau_on_ns) if plan.tau_on_ns > 0 else 1.0
            return eff / tau, tau
        return eff * (float(rt.Gamma_ns) + float(rt.k_nr_ns)), 0.0
    return 0.0, float(plan.tau_on_ns)


def _pump_fano(plan):
    """(F_dot, notes) of an electrical plan: the DELIVERED pump Fano factor F_dot = 1 - eta (1 - F_pump)
    (loading.f8b_thin_fano). The compose solvers do not model sub-Poissonian loading in the counting hierarchy (separate
    spec), so when F_pump != 1 the value is only reported and a note says the effect is not in g2 (R2: M-2/M-4 were
    silent no-ops)."""
    itf = plan.interface
    if itf is None:
        return None, []
    F_p, eta = float(itf.F_p), float(itf.eta_capture)
    F_dot = 1.0 - eta * (1.0 - F_p)
    notes = []
    if abs(F_p - 1.0) > 1e-12:
        notes.append(f"pump Fano factor not modelled by this solver (F_pump = {F_p:.4g}, F_dot = {F_dot:.4g} reported only)")
    return F_dot, notes


def _eid(plan) -> float:
    """Excitation-induced dephasing coefficient c_EID in ps (gamma_EID = c |Omega|^2); 0 = off (devlib-fx3b)."""
    return float(plan.inputs.get("eid_ps", 0.0) or 0.0)


def _p_trion(plan) -> float:
    """Dark-trion probability at a pi pulse (p_trion(area) = p sin^2(area/2)); 0 = off (devlib-fx3b)."""
    return float(plan.inputs.get("p_trion", 0.0) or 0.0)


def _series_xy(x, y, unit_x="ns"):
    return {"x": [float(v) for v in x], "y": [float(v) for v in y], "x_unit": unit_x}


# ---------------------------------------------------------------- shaped pump profile (devlib-p10e)
_RC_TAIL_TR = 12.0                 # [A] RC tail truncated at tau_on + 12 t_r (the lost weight is < e^-12 of the tail)
_GAUSS_HALF_SIGMAS = 4.0           # [A] Gaussian window = +/- 4 sigma, renormalised to mu
_CELLS_PER_TR = 16                 # [A] steps per t_r (rc) / per sigma-pair (gauss); anchors converge < 1e-3


def pump_cells(shape, mu, tau_on_ns, t_rise_ns=0.0, n_per=_CELLS_PER_TR):
    """Piecewise-constant, CELL-AVERAGED incoherent pump r(t) [DR, brief 1.1]: -> (dts_ns, r_k_per_ns) with
    sum r_k dt_k = mu exactly (each cell carries the exact integral of r over it).
      rect : r0 = mu/tau_on on [0, tau_on)
      rc   : r0 (1 - e^{-t/t_r}) on [0, tau_on), r0 (1 - e^{-tau_on/t_r}) e^{-(t - tau_on)/t_r} after, r0 = mu/tau_on
             (the rise deficit equals the tail excess, so mu = r0 tau_on); window tau_on + 12 t_r
      gauss: r = mu N(t; t0, sigma), FWHM = tau_on (the pump pulse's own FWHM), window +/- 4 sigma, renormalised
    t_rise_ns <= 0 turns rc into rect. The step count grows with n_per; the shape anchors (Sartison, brief 1.3) are
    converged to < 1e-3 at the default."""
    mu, tau = float(mu), float(tau_on_ns)
    if shape == "rc" and not t_rise_ns > 0.0:
        shape = "rect"
    if shape == "rect":
        return [tau], [mu / tau]
    if shape == "rc":
        tr = float(t_rise_ns)
        r0 = mu / tau
        n1 = min(400, max(8, int(math.ceil(n_per * tau / tr))))
        W = tau + _RC_TAIL_TR * tr
        n2 = min(600, max(8, int(math.ceil(n_per * (W - tau) / tr))))
        nodes = np.concatenate([np.linspace(0.0, tau, n1 + 1), np.linspace(tau, W, n2 + 1)[1:]])
        Ft = r0 * (tau - tr * (1.0 - math.exp(-tau / tr)))
        amp = r0 * (1.0 - math.exp(-tau / tr)) * tr

        def F(t):
            return np.where(t <= tau, r0 * (t - tr * (1.0 - np.exp(-t / tr))),
                            Ft + amp * (1.0 - np.exp(-(t - tau) / tr)))
    elif shape == "gauss":
        sig = tau / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        W = 2.0 * _GAUSS_HALF_SIGMAS * sig
        nodes = np.linspace(0.0, W, int(math.ceil(W / (2.0 * sig / n_per))) + 1)
        from scipy.special import erf

        def F(t):
            return 0.5 * (1.0 + erf((t - 0.5 * W) / (sig * math.sqrt(2.0))))
    else:
        raise ValueError(f"pump_profile {shape!r} not supported (rect|rc|gauss)")
    Fn = F(nodes)
    dts = np.diff(nodes)
    dF = np.diff(Fn)
    if shape == "gauss":
        dF = dF / (Fn[-1] - Fn[0]) * mu
    else:
        dF = dF * (mu / (Fn[-1] - Fn[0]))        # renormalise the (tiny) truncated tail to mu
    return [float(d) for d in dts], [float(x) for x in dF / dts]


def _shaped(plan) -> bool:
    """True when the plan asks for anything beyond the rectangular pump with no bias-dependent escape."""
    shaped_pump = plan.pump_profile == "gauss" or (plan.pump_profile == "rc" and plan.t_rise_ns > 0.0)
    return bool(shaped_pump or plan.k_esc_on_ns > 0.0 or plan.k_esc_off_ns > 0.0)


def _plan_cells(plan, r, tau_on):
    """(dts, r_k) of a shaped plan: mu = r tau_on (the rectangular operating point) is preserved."""
    return pump_cells(plan.pump_profile, float(r) * float(tau_on), float(tau_on), float(plan.t_rise_ns))


def _pump_window_ns(plan, r, tau_on) -> float:
    """Length of the pump window including an RC tail (the rectangular tau_on when unshaped)."""
    if not _shaped(plan):
        return float(tau_on)
    return float(sum(_plan_cells(plan, r, tau_on)[0]))


def _el_extras(plan, rt) -> dict:
    """Escape / EL-decay report of a shaped plan (empty for the default rectangular, escape-free plan)."""
    if not _shaped(plan):
        return {}
    G, k, koff = float(rt.Gamma_ns), float(rt.k_nr_ns), float(plan.k_esc_off_ns)
    return {"pump_profile": plan.pump_profile, "t_rise_ns": float(plan.t_rise_ns),
            "k_esc_on_ns": float(plan.k_esc_on_ns), "k_esc_off_ns": koff,
            "tau_EL_ps": 1e3 / (G + k + koff), "tau_rad_ps": 1e3 / G,
            "eta_rad_esc": G / (G + k + koff)}


def _window_counts(pc, shape, area, duration, rho0=None, n_steps=None):
    """(rho at the end of the pulse window, photons counted inside the window) for the pieces `pc`,
    started in `rho0` (default |0><0|, the ground state). Same grid and midpoint propagation as
    resonant_source._stats; the photon count is the first factorial moment m1 of the counting hierarchy
    at the window end. Used for eta_prep = P_X(end) + emitted (FX2b)."""
    edges, om = rs._grid(shape, area, duration, n_steps, pc.max_dt)
    dt = np.diff(edges)
    N = pc.dim ** 2
    V = np.zeros(3 * N, dtype=complex)
    V[:N] = lb.vec(rs._ground(pc.dim) if rho0 is None else rho0)
    cache = {}
    for o, d in zip(om, dt):
        key = (round(float(o), 14), round(float(d), 14))
        if key not in cache:
            cache[key] = expm(rs._aug(pc.L0 + o * pc.L1, pc.J) * d)
        V = cache[key] @ V
    rho = lb.unvec(V[:N].copy(), pc.dim)
    m1 = float((lb.trace_row(pc.dim) @ V[N:2 * N]).real)
    return rho, m1


def _window_p0(pc, shape, area, duration, keep, n_steps=None):
    """FX2c: probability that NO counted photon is emitted inside the pulse window and the system is in the
    unexcited subspace `keep` (basis indices, populations only) at the window end, started from the ground
    state. The no-jump generator is L0 + Omega L1 - J (J = the jump superoperator the counting hierarchy
    uses), so the weight is exactly p0 of the photon-number distribution; eta_prep = 1 - p0. Same grid and
    midpoint propagation as _window_counts."""
    edges, om = rs._grid(shape, area, duration, n_steps, pc.max_dt)
    dt = np.diff(edges)
    v = lb.vec(rs._ground(pc.dim))
    cache = {}
    for o, d in zip(om, dt):
        key = (round(float(o), 14), round(float(d), 14))
        if key not in cache:
            cache[key] = expm((pc.L0 + o * pc.L1 - pc.J) * d)
        v = cache[key] @ v
    rho = lb.unvec(v, pc.dim)
    return float(sum(rho[i, i].real for i in keep))


# ---------------------------------------------------------------- rate
class RateSolver(Solver):
    name = "rate"
    label = "Rate equations"
    summary = ("Incoherent pump on the G-X-XX cascade: pulse_counting.pulse_g2 (pulsed) / cw_g2 (cw). "
               "No coherence: HOM and Rabi physics unavailable.")

    def supports(self, plan, em, rt):
        if plan.family == "coherent":
            return False, ("rate equations cannot describe Rabi oscillation / coherent drive "
                           "(pulse_counting validity guard); use lindblad")
        bad = _bad_needs(plan, allowed=("cascade",))
        if bad:
            return False, "; ".join(bad)
        r, tau = _interp_pump(plan, rt)
        if r <= 0.0:
            return False, "no incoherent pump rate in the plan"
        if plan.timing == "pulsed":
            f = _f_rep(plan)
            if f is None or tau <= 0.0:
                return False, "pulsed plan needs f_rep_MHz and a pump window"
            if 1e3 / f <= _pump_window_ns(plan, r, tau):
                return False, "repetition period must exceed the pump window"
        return True, ""

    def _rates(self, plan, em, rt):
        gX = float(rt.Gamma_ns)
        gXX = float(em.gamma_xx_ratio) * float(rt.Gamma0_ns)    # XX is not Purcell-enhanced [A]
        kX = float(rt.k_nr_ns)
        if plan.timing == "cw":
            kX = kX + float(plan.k_esc_on_ns)          # cw: the bias-on escape acts all the time (0 by default)
        kXX = kX * float(em.gamma_xx_ratio)
        p = float(plan.inputs.get("pump_ratio_xx", 1.0)) if "cascade" in plan.needs else 0.0
        return gX, gXX, kX, kXX, p

    def run(self, plan, em, rt):
        r, tau_on = _interp_pump(plan, rt)
        if plan.timing == "cw":
            # R2: the cw pump broadens the dot line the cavity filter sees, P = r (1 + p) (ch08 eq. dp-lor); r keeps
            # the planner's definition (unbroadened Gamma), the rates use the broadened L
            p_xx = float(plan.inputs.get("pump_ratio_xx", 1.0)) if "cascade" in plan.needs else 0.0
            if not any(str(n).startswith("cw pump P =") for n in rt.notes):      # the engine may have applied it (combine hook)
                rt = _rates_mod.apply_pump(rt, r * (1.0 + p_xx))
        gX, gXX, kX, kXX, pr = self._rates(plan, em, rt)
        reasons = ["rate solver: capture/relaxation jitter and coherence are not modelled"]
        F_dot, fano_notes = _pump_fano(plan)
        reasons += fano_notes
        fano_extra = {} if F_dot is None else {"F_dot": F_dot}
        bg = float(plan.inputs.get("bg_frac", 0.0) or 0.0)
        bgf = (lambda g: 1.0 - (1.0 - bg) ** 2 * (1.0 - g)) if bg > 0.0 else (lambda g: g)
        if float(plan.inputs.get("p_refill", 0.0) or 0.0) > 0.0 or float(plan.inputs.get("k0_ns", 0.0) or 0.0) > 0.0:
            reasons.append("p_refill / thermal recapture are not modelled by the rate solver (use lindblad_incoherent)")
        if bg > 0.0:
            reasons.append("bg_frac applied to g2 as 1 - rho^2 (1 - g2), rho = 1 - bg_frac [A]")
        if plan.timing == "pulsed":
            period = 1e3 / float(plan.f_rep_MHz)
            extra_el = _el_extras(plan, rt)
            if extra_el:
                ratio = float(em.gamma_xx_ratio)
                dts, rks = _plan_cells(plan, r, tau_on)
                kon, koff = float(plan.k_esc_on_ns), float(plan.k_esc_off_ns)
                steps = [(dt, rk, kX + kon, (kX + kon) * ratio) for dt, rk in zip(dts, rks)]
                W = float(sum(dts))
                out = pcnt.pulse_g2(0.0, gX, gXX, kX, kXX, 1.0, 0.0, 0.0, period - W, pump_ratio=pr,
                                    gate_ns=plan.gate_ns, r_profile=steps,
                                    k_dark=(kX + koff, (kX + koff) * ratio))
                extra_el.update(window_ns=W, mu=float(r * tau_on))
                reasons.append(f"shaped pump ({plan.pump_profile}, window {W * 1e3:.0f} ps) with bias-dependent escape "
                               "(k_on / k_off); EL decay is reported, never fed into Gamma0")
            else:
                out = pcnt.pulse_g2(r, gX, gXX, kX, kXX, 1.0, 0.0, tau_on, period - tau_on, pump_ratio=pr,
                                    gate_ns=plan.gate_ns)
            return SolverResult(g2_0=float(bgf(out["g2"])), n_ph=float(out["mean_counts"]), I_single=NAN,
                                inversion=None,
                                extra={"p_loaded": 1.0 - math.exp(-r * tau_on), "g2_signal": float(out["g2"]),
                                       "converged": bool(out["converged"]), "notes": reasons, **fano_extra, **extra_el},
                                valid=bool(out["converged"]),
                                reasons=() if out["converged"] else ("periodic steady state did not converge",))
        scale = max(gX + kX, 1e-12)
        taus = np.linspace(0.0, 12.0 / scale, 401)
        g2 = cw_g2.g2_cw(taus, r, gX, gXX, 1.0, 0.0, kX, kXX, pr)
        P = cw_g2.steady_state(r, gX, gXX, kX, kXX, pr)
        return SolverResult(g2_0=float(bgf(g2[0])), n_ph=NAN, I_single=NAN, inversion=None,
                            extra={"flux_ns": float(gX * P[1]), "tau_dip_ns": 1.0 / (gX + kX + r),
                                   "g2_signal": float(g2[0]),
                                   "series": {"g2_tau": _series_xy(taus, g2)}, "notes": reasons,
                                   **fano_extra, **_el_extras(plan, rt)})

    def eta_seconds(self, plan, n_points):
        return 0.005 * max(1, int(n_points))


# ---------------------------------------------------------------- lindblad (coherent)
class LindbladSolver(Solver):
    name = "lindblad"
    label = "Lindblad (coherent)"
    summary = ("Master equation of the driven two-level emitter (or X-XX cascade for TPE): pulsed g2, "
               "photon number, HOM overlap, inversion; cw resonance fluorescence g2(tau).")

    def supports(self, plan, em, rt):
        if plan.family != "coherent":
            return False, "coherent drives only; use lindblad_incoherent or rate"
        bad = _bad_needs(plan, allowed=("cascade", "td_detuning"))
        if bad:
            return False, "; ".join(bad)
        why = _pulsed_coherent_check(plan)
        if why:
            return False, why
        if "td_detuning" in plan.needs or _eid(plan) > 0.0:
            what = "td_detuning" if "td_detuning" in plan.needs else "excitation-induced dephasing (eid_ps)"
            if "cascade" in plan.needs or plan.timing != "pulsed":
                return False, f"{what}: pulsed two-level plans only (no TPE cascade, no cw)"
            ok, why = self._td_solver().supports(self._td_plan(plan), em, rt)
            if not ok:
                return False, f"{what} (polaron machinery, phonons off): " + why
        if _p_trion(plan) > 0.0 and plan.timing != "pulsed":
            return False, "p_trion: pulsed plans only (the dark-trion branch is a per-pulse probability)"
        if "cascade" in plan.needs:
            if plan.timing != "pulsed":
                return False, "TPE needs a pulsed plan"
            if "XX" not in em.transitions:
                return False, "TPE needs 'XX' in the material transitions"
            if abs(_det_X_meV(plan)) > 1e-12:
                return False, "TPE offset from the two-photon resonance is not supported"
        elif plan.timing == "cw" and plan.omega_cw_ns <= 0.0:
            return False, "cw drive with zero Rabi frequency has no emission"
        return True, ""

    @staticmethod
    def _td_solver():
        return _get("solver", "polaron")

    @staticmethod
    def _td_plan(plan):
        """The plan for the polaron step-propagator machinery (it needs plan.inputs["T_K"]; phonons are off)."""
        return plan if "T_K" in plan.inputs else _dc_replace(plan, inputs={**plan.inputs, "T_K": 0.0})

    def _td(self, plan, em, rt):
        """Time-dependent detuning (Stark gate det_profile / chirp): the polaron solver's piecewise-constant
        L(t) = L0 + Omega(t) L1 + (det(t) - det0) L_det propagators with the phonon coupling removed, so the same
        counting hierarchy, eta_prep and Grange HOM apply [DR]."""
        res = self._td_solver().run(self._td_plan(plan), _dc_replace(em, phonon=None), rt)
        extra = dict(res.extra)
        extra["notes"] = list(extra.get("notes", [])) + ["td_detuning / eid_ps: lindblad path = polaron step "
                                                         "propagators, phonons off"]
        return SolverResult(g2_0=res.g2_0, n_ph=res.n_ph, I_single=res.I_single, inversion=res.inversion,
                            extra=extra, valid=res.valid, reasons=res.reasons)

    def _zero_area(self, plan, kind="X"):
        """Trivial zero-excitation result (pulse area 0): no excitation, no photons; g2 and HOM are undefined."""
        extra = {"det_X_meV": _det_X_meV(plan), "eta_prep": 0.0, "n_exc": 0.0, "photons_in_window": 0.0,
                 "P_X_end": 0.0, "notes": ["zero pulse area: the emitter is never excited (n_ph = 0; g2 and HOM undefined)"]}
        if kind == "tpe":
            extra.update(n_XX=0.0)
        return SolverResult(g2_0=NAN, n_ph=0.0, I_single=NAN, inversion=0.0, extra=extra, valid=True,
                            reasons=("zero pulse area: no excitation, g2 and HOM undefined",))

    def run(self, plan, em, rt):
        res = self._run(plan, em, rt)
        pt = _p_trion(plan)
        if pt > 0.0 and not ("td_detuning" in plan.needs or _eid(plan) > 0.0) and plan.pulses:
            # dark-trion branch (devlib-fx3b); the polaron path applies it itself
            area = math.sqrt(sum(float(q.area_rad) ** 2 for q in plan.pulses))     # TPE: area_pi pi, 1 = the TPE pi pulse
            res = _polaron.apply_trion(res, _polaron.trion_probability(pt, area))
        return res

    def _run(self, plan, em, rt):
        if "cascade" in plan.needs:
            if plan.pulses and float(plan.pulses[0].area_rad) == 0.0:
                return self._zero_area(plan, "tpe")
            return self._tpe(plan, em, rt)
        if plan.timing == "cw":
            return self._cw(plan, em, rt)
        if plan.pulses and float(plan.pulses[0].area_rad) == 0.0:
            return self._zero_area(plan)
        if "td_detuning" in plan.needs or _eid(plan) > 0.0:
            return self._td(plan, em, rt)               # polaron step propagators (also apply the trion branch)
        p = plan.pulses[0]
        dur = p.fwhm_ps * 1e-3
        tl = rs.TwoLevel(float(rt.Gamma_ns), k_ns=float(rt.k_nr_ns), gamma_star_ns=float(rt.gamma_star_ns),
                         delta_ns=_det_X_ns(plan))
        st = rs.pulsed_photon_statistics(tl, p.shape, p.area_rad, dur, _f_rep(plan))
        I = rs.pulsed_indistinguishability(tl, p.shape, p.area_rad, dur, tol=_HOM_TOL,
                                           max_doublings=_HOM_DOUBLINGS)
        pcs = rs._pieces(tl)
        rho_w, emitted = _window_counts(pcs, p.shape, p.area_rad, dur)
        p0 = _window_p0(pcs, p.shape, p.area_rad, dur, (0,))
        return SolverResult(g2_0=float(st["g2"]), n_ph=float(st["mean_photons"]), I_single=float(I),
                            inversion=float(st["P_e_post"]),
                            extra={"det_X_meV": _det_X_meV(plan),
                                   "eta_prep": 1.0 - p0, "n_exc": float(rho_w[1, 1].real) + emitted,
                                   "photons_in_window": emitted, "P_X_end": float(st["P_e_post"])})

    def _tpe(self, plan, em, rt):
        p = plan.pulses[0]
        dur = p.fwhm_ps * 1e-3
        dxx = float(plan.inputs.get("delta_xx_meV", em.delta_xx_meV))   # card convention (negative = bound)
        delta = rs.tpe_detuning_ns(dxx)
        gs = float(rt.gamma_star_ns)
        kX = float(rt.k_nr_ns)
        casc = rs.Cascade(float(rt.Gamma_ns), float(em.gamma_xx_ratio) * float(rt.Gamma0_ns), k_X_ns=kX,
                          k_XX_ns=kX * float(em.gamma_xx_ratio),
                          deph=((gs, (0.0, 1.0, 0.0)),) if gs > 0 else (), delta_ns=delta)
        area_fn = getattr(rs, "_tpe_area_cached", None)
        pi_area = float(area_fn(delta, p.shape, dur)) if area_fn else float(rs.tpe_pi_area(delta, p.shape, dur))
        area = (p.area_rad / math.pi) * pi_area
        st = rs.tpe_photon_statistics(casc, p.shape, area, dur, _f_rep(plan))
        hom = rs.tpe_cascade_hom(casc, p.shape, area, dur)
        # FX2b: XX-line photon stream (the line a TPE experiment filters for the biexciton photon) and eta_prep
        st_xx = rs.tpe_photon_statistics(casc, p.shape, area, dur, _f_rep(plan), channel="XX")
        pcx = rs._pieces_cascade(casc, ("XX",))
        rho_w, emitted_xx = _window_counts(pcx, p.shape, area, dur)
        p0_xx = _window_p0(pcx, p.shape, area, dur, (0, 1))       # no XX photon and not in XX (G or X)
        s = casc.system()
        i_xx = 2 * s.dim + 2                                   # vec index of rho[XX, XX]
        tau_xx_ps = 1e3 / float(-s.L[i_xx, i_xx].real)         # XX population decay rate read off the Liouvillian
        prep = rs.tpe_prepare(casc, p.shape, area, dur)
        A_xx = math.sqrt(float(casc.gamma_XX_ns)) * s.ops["c_XX"]
        I_xx, N_xx = lb.indistinguishability(s.L, s.dim, prep["rho"], A_xx, return_counts=True)
        P_xx_end = float(st["P_XX_post"])
        return SolverResult(g2_0=float(st["g2"]), n_ph=float(st["mean_photons"]), I_single=float(hom["I"]),
                            inversion=P_xx_end,
                            extra={"tpe_pi_area_rad": pi_area, "tpe_area_rad": area,
                                   "P_prepared": list(hom["P"]),
                                   "eta_prep": 1.0 - p0_xx, "n_exc": float(rho_w[2, 2].real) + emitted_xx,
                                   "photons_in_window": emitted_xx, "P_X_end": P_xx_end,
                                   "tau_XX_ps": float(tau_xx_ps), "g2_XX": float(st_xx["g2"]),
                                   "n_XX": float(st_xx["mean_photons"]), "hom_XX": float(I_xx),
                                   "N_XX_prepared": float(N_xx)},
                            reasons=("TPE HOM prepares with dissipation off (optimistic by O(Gamma tau_p)); "
                                     "hom_XX likewise",))

    def _cw(self, plan, em, rt):
        gs = float(rt.gamma_star_ns)
        s = lb.build_system(levels=2, gamma_X_ns=float(rt.Gamma_ns), k_X=float(rt.k_nr_ns),
                            deph=[(gs, (0.0, 1.0))] if gs > 0 else [], det_X_ns=_det_X_ns(plan),
                            omega_ns=float(plan.omega_cw_ns))
        rate = float(rt.Gamma_ns) + float(rt.k_nr_ns) + gs
        taus = np.linspace(0.0, 24.0 / rate, 401)
        g2, rho, I = lb.g2_tau(s, taus)
        pxx = float(rho[1, 1].real)
        coh = float(abs(rho[0, 1]) ** 2 / pxx) if pxx > 0 else NAN
        return SolverResult(g2_0=float(g2[0]), n_ph=NAN, I_single=NAN, inversion=pxx,
                            extra={"flux_ns": float(I), "coherent_fraction": coh,
                                   "series": {"g2_tau": _series_xy(taus, g2)}},
                            reasons=("cw: no per-pulse photon number or HOM overlap",))

    def eta_seconds(self, plan, n_points):
        per = 0.01 if plan.timing == "cw" else (0.6 if plan.pulses and plan.pulses[0].shape == "gauss" else 0.25)
        if "cascade" in plan.needs:
            per = 2.0
        return per * max(1, int(n_points))


# ---------------------------------------------------------------- lindblad_incoherent
def _s3d_ops(with_res: bool, extra: int = 0):
    dim = (3 if with_res else 2) + int(extra)
    def e(i, j):
        m = np.zeros((dim, dim), dtype=complex)
        m[i, j] = 1.0
        return m
    return dim, e


class IncoherentLindbladSolver(Solver):
    name = "lindblad_incoherent"
    label = "Lindblad (incoherent pump, S3d)"
    summary = ("Three-level G-X-R scheme: incoherent pump into the reservoir R for the window, R->X at "
               "1/tau_cap, X->G radiative. Pulsed g2 by the counting hierarchy; HOM with capture jitter.")

    def supports(self, plan, em, rt):
        if plan.family not in ("incoherent", "electrical"):
            return False, "incoherent / electrical drives only; use lindblad for coherent drives"
        bad = _bad_needs(plan)
        if bad:
            return False, "; ".join(bad)
        if "cascade" in plan.needs:
            return False, "cap-2 (XX) loading needs the S4 cascade: use the rate solver"
        if plan.timing != "pulsed":
            return False, "pulsed plans only (cw: use rate)"
        r, tau = _interp_pump(plan, rt)
        f = _f_rep(plan)
        if r <= 0 or tau <= 0 or f is None:
            return False, "needs pump rate, pump window and f_rep_MHz"
        if 1e3 / f <= _pump_window_ns(plan, r, tau):
            return False, "repetition period must exceed the pump window"
        return True, ""

    @staticmethod
    def _extras(plan):
        """(p_refill, tau_refill_ns, bg_frac, k_esc_ns) of the plan (FX2b); all zero = the original S3d scheme.
        k_esc = k0 exp(-E_act/kT) at plan.inputs["T_K"] (Arrhenius, [E] class model, p_shell)."""
        inp = plan.inputs
        p_ref = float(inp.get("p_refill", 0.0) or 0.0)
        tau_ref = float(inp.get("tau_refill_ps", 100.0)) * 1e-3
        bg = float(inp.get("bg_frac", 0.0) or 0.0)
        k0 = float(inp.get("k0_ns", 0.0) or 0.0)
        k_esc = 0.0
        if k0 > 0.0:
            T = float(inp.get("T_K", 0.0) or 0.0)
            E = float(inp.get("E_act_meV", 0.0) or 0.0)
            k_esc = k0 * math.exp(-min(E / (KB * T), 700.0)) if T > 0.0 else (k0 if E <= 0.0 else 0.0)
        return p_ref, tau_ref, bg, k_esc

    def run(self, plan, em, rt):
        r, tau_on = _interp_pump(plan, rt)
        period = 1e3 / float(plan.f_rep_MHz)
        Gam, k, gs = float(rt.Gamma_ns), float(rt.k_nr_ns), float(rt.gamma_star_ns)
        relax = float(plan.relax_ns)
        with_res = relax > 0.0
        p_ref, tau_ref, bg, k_esc = self._extras(plan)
        n0 = 3 if with_res else 2
        iW = n0 if p_ref > 0.0 else None
        dim, e = _s3d_ops(with_res, 1 if iW is not None else 0)
        # basis: 0 = G, 1 = X, 2 = R (reservoir), iW = W (dot emptied after an emission, refill pending)
        iP = 2 if with_res else 1                                   # level the pump / thermal capture feeds
        c_pump = np.sqrt(r) * e(iP, 0)
        if iW is None:
            c_ops_off = [np.sqrt(Gam) * e(0, 1)]
            J = Gam * lb.jump_super(e(0, 1))
        else:
            # refilling: with probability p_refill an emission is followed by a second, uncorrelated capture
            # after ~tau_refill (X -> W branch, W -> X at 1/tau_refill); photons count either branch [A]
            c_ops_off = [np.sqrt(Gam * (1.0 - p_ref)) * e(0, 1), np.sqrt(Gam * p_ref) * e(iW, 1)]
            J = Gam * (1.0 - p_ref) * lb.jump_super(e(0, 1)) + Gam * p_ref * lb.jump_super(e(iW, 1))
            c_ops_off.append(np.sqrt(1.0 / tau_ref) * e(1, iW))
        if with_res:
            c_ops_off.append(np.sqrt(1.0 / relax) * e(1, 2))
        if k > 0:
            c_ops_off.append(np.sqrt(k) * e(0, 1))
        if gs > 0:
            c_ops_off.append(np.sqrt(gs) * e(1, 1))
        H0 = np.zeros((dim, dim), dtype=complex)
        # HOM generator: one prepared exciton. Thermal escape X -> R delays/jitters the emission (kept); the
        # thermal recapture G -> R below is an uncorrelated background stream and is NOT part of it.
        c_ops_hom = list(c_ops_off)
        if k_esc > 0.0 and with_res:
            c_ops_hom.append(np.sqrt(k_esc) * e(2, 1))
        extra_el = _el_extras(plan, rt)
        k_on, k_dark = float(plan.k_esc_on_ns), float(plan.k_esc_off_ns)
        c_ops_hom_base = list(c_ops_hom)
        if extra_el and k_dark > 0.0:
            c_ops_hom.append(np.sqrt(k_dark) * e(0, 1))            # dark-window bias-state escape (devlib-p10e)
        L_hom = lb.liouvillian(H0, c_ops_hom)
        c_ops_cnt = list(c_ops_hom)
        if k_esc > 0.0:
            c_ops_cnt.append(np.sqrt(k_esc) * e(iP, 0))            # recapture of thermally released carriers
        L_off = lb.liouvillian(H0, c_ops_cnt)
        if extra_el:
            # shaped pump window: L_k = L_on_base + r_k D[|P><G|], escape k_on on, k_off in the dark window
            c_on = list(c_ops_hom_base) + ([np.sqrt(k_esc) * e(iP, 0)] if k_esc > 0.0 else [])
            if k_on > 0.0:
                c_on.append(np.sqrt(k_on) * e(0, 1))
            L_base = lb.liouvillian(H0, c_on)
            D_pump = lb.dissipator(e(iP, 0))
            dts, rks = _plan_cells(plan, r, tau_on)
            W = float(sum(dts))
            profile = [(dt, L_base + rk * D_pump) for dt, rk in zip(dts, rks)]
            out = lb.pulsed_counting(None, L_off, J, dim, 0.0, period - W, gate_ns=plan.gate_ns, profile=profile)
            extra_el.update(window_ns=W, mu=float(r * tau_on))
        else:
            L_on = lb.liouvillian(H0, c_ops_cnt + [c_pump])
            out = lb.pulsed_counting(L_on, L_off, J, dim, tau_on, period - tau_on, gate_ns=plan.gate_ns)
        A = np.sqrt(Gam) * e(0, 1)
        rho0 = e(2, 2) if with_res else e(1, 1)
        I = lb.indistinguishability(L_hom, dim, rho0, A)
        g2_sig = float(out["g2"])
        rho_sig = 1.0 - bg
        g2 = 1.0 - rho_sig ** 2 * (1.0 - g2_sig) if bg > 0.0 else g2_sig   # uncorrelated background (Brouri form)
        reasons = ["HOM assumes instantaneous loading at t = 0 (window re-excitation ignored) [A]"]
        F_dot, fano_notes = _pump_fano(plan)
        reasons += fano_notes
        if iW is not None:
            reasons.append("p_refill: refill is a delayed second capture (trap refilling) [A]")
        if bg > 0.0:
            reasons.append("bg_frac: g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac; HOM and n_ph are the "
                           "signal-only values (background does not degrade the HOM here) [A]")
        if k_esc > 0.0:
            reasons.append(f"thermal escape/recapture (Arrhenius [E] class model): k_esc = {k_esc:.4g}/ns")
        if not out["converged"]:
            reasons.append("periodic steady state did not converge")
        extra = {"p_loaded": 1.0 - math.exp(-r * tau_on), "g_relax_ns": (1.0 / relax if with_res else None),
                 "g2_signal": g2_sig, "p_refill": p_ref, "bg_frac": bg, "k_esc_ns": k_esc, "notes": reasons,
                 **({} if F_dot is None else {"F_dot": F_dot}), **extra_el}
        return SolverResult(g2_0=float(g2), n_ph=float(out["mean_counts"]), I_single=float(I),
                            inversion=None, extra=extra,
                            valid=bool(out["converged"]),
                            reasons=(("periodic steady state did not converge",) if not out["converged"] else ()))

    def eta_seconds(self, plan, n_points):
        return 0.05 * max(1, int(n_points))


# ---------------------------------------------------------------- jc
class _CavOut:
    """TwoLevel-shaped adapter whose emitted field is the cavity output sqrt(kappa) a, so
    resonant_source.pulsed_indistinguishability yields the cavity-photon overlap."""

    def __init__(self, builder, kappa):
        self._b = builder
        self.gamma_ns = float(kappa)

    def system(self, omega_ns=0.0):
        s = self._b(omega_ns)
        a = s.ops["a"]
        s.ops = dict(s.ops, c_X=a)
        s.photon = {"X": (a, self.gamma_ns)}
        return s


JC_EID_REASON = "JC regime: EID not supported by the JC solver"


class _JCRegimeFallback:
    """auto-selection wrapper (FX4): the plan is in the JC regime (strong coupling or 4g/kappa > 0.5) but needs
    excitation-induced dephasing, which the JC solver does not model. The weak-coupling solver still runs, but its
    result is marked INVALID with JC_EID_REASON instead of being passed off as a valid number."""

    def __init__(self, inner):
        self._inner = inner

    def __getattr__(self, name):
        return getattr(self._inner, name)

    def run(self, plan, em, rt):
        res = self._inner.run(plan, em, rt)
        return SolverResult(g2_0=res.g2_0, n_ph=res.n_ph, I_single=res.I_single, inversion=res.inversion,
                            extra=dict(res.extra), valid=False, reasons=tuple(res.reasons) + (JC_EID_REASON,))


class JCSolver(Solver):
    name = "jc"
    label = "Jaynes-Cummings (dot + cavity)"
    summary = ("Dot coupled to one cavity mode (Fock cutoff 3), Lindblad with kappa, g. For 4g/kappa > 0.5 "
               "or strong coupling: lifetime, cavity-output photons per pulse, g2, HOM.")
    FOCK = 3

    _JC_NEEDS = {
        "cascade": "the biexciton ladder (two-photon excitation, cap-2 loading) is not included in the jc solver (two-level dot + cavity); "
                   "a TPE + cavity run needs a three-level ladder with the cavity, which is not implemented",
        "phonons": "phonon-dressed dynamics (polaron) with a cavity mode are not implemented in the jc solver",
        "two_colour": "two-colour drives with a cavity mode are not implemented in the jc solver",
        "td_detuning": "time-dependent detuning (chirp / Stark gate) with a cavity mode is not implemented in the jc solver",
    }

    def supports(self, plan, em, rt):
        if getattr(rt, "kappa_ns", None) is None or getattr(rt, "g_ns", 0.0) <= 0.0:
            return False, "jc needs a cavity mode (kappa and g)"
        for n in sorted(plan.needs):
            return False, self._JC_NEEDS.get(n, f"need {n!r} is not supported by the jc solver")
        if _eid(plan) > 0.0:
            return False, JC_EID_REASON
        if plan.family == "coherent":
            if plan.timing == "pulsed":
                why = _pulsed_coherent_check(plan)
                if why:
                    return False, why
            elif plan.timing == "cw":
                if plan.omega_cw_ns <= 0.0:
                    return False, "cw drive with zero Rabi frequency has no emission"
            else:
                return False, f"unknown timing {plan.timing!r}"
            return True, ""
        if plan.family in ("incoherent", "electrical"):
            if _shaped(plan):
                return False, ("shaped pump profile / bias-dependent escape (pump_profile, k_esc_on/off) is not implemented in the "
                               "jc solver; use the default rectangular pump with no escape")
            inp = plan.inputs
            if float(inp.get("p_refill", 0.0) or 0.0) > 0.0 or float(inp.get("k0_ns", 0.0) or 0.0) > 0.0:
                return False, "p_refill / thermal escape-recapture are not implemented in the jc solver"
            r, tau = _interp_pump(plan, rt)
            if r <= 0.0:
                return False, "no incoherent pump rate in the plan"
            if plan.timing == "pulsed":
                f = _f_rep(plan)
                if f is None or tau <= 0.0:
                    return False, "pulsed plan needs f_rep_MHz and a pump window"
                if 1e3 / f <= tau:
                    return False, "repetition period must exceed the pump window"
            return True, ""
        return False, f"drive family {plan.family!r} is not supported by the jc solver"

    def _builder(self, plan, em, rt):
        gs = float(rt.gamma_star_ns)
        det = _det_X_ns(plan)
        cav = dict(g_ns=float(rt.g_ns), kappa_ns=float(rt.kappa_ns),
                   det_ns=det + float(rt.cav_minus_dot_ns), n_max=self.FOCK)

        def build(omega_ns=0.0, r_ns=0.0, k_extra=0.0):
            return lb.build_system(levels=2, gamma_X_ns=float(rt.Gamma_bg_ns), k_X=float(rt.k_nr_ns) + float(k_extra),
                                   deph=[(gs, (0.0, 1.0))] if gs > 0 else [], det_X_ns=det,
                                   omega_ns=omega_ns, r_ns=float(r_ns), cavity=cav)
        return build

    def lifetime(self, plan, em, rt):
        """(tau_ps, n_cav, n_leak, tau_1e_ps): tau_ps = 1/|Re lambda| of the slowest decay eigenvalue of the
        Liouvillian restricted to the one-excitation population/coherence sector span{|X,0>, |G,1>} x its dual
        (the radiative sector reached by an instantaneous pi pulse); for a weakly coupled dot it is the Purcell
        lifetime, for strong coupling the slow polariton population decay (= 1/(2 |Re lambda_M|) of the
        non-Hermitian 2x2 amplitude matrix) [DR]. tau_1e_ps is the old 1/e time of the emitted flux after its
        peak (informational: it mixes both eigenmodes). Photons per excitation: exact projected Laplace solve."""
        build = self._builder(plan, em, rt)
        s = build(0.0)
        dim, nc = s.dim, self.FOCK
        rho0 = np.zeros((dim, dim), dtype=complex)
        ix = 1 * nc + 0                               # |X, n=0>
        rho0[ix, ix] = 1.0
        kap, gl = float(rt.kappa_ns), float(rt.Gamma_bg_ns)
        a = s.ops["a"]
        row = gl * lb.vec(s.ops["P_X"].T) + kap * lb.vec((a.conj().T @ a).T)
        lam = np.linalg.eigvals(s.L)
        nz = np.abs(lam.real)[np.abs(lam.real) > 1e-9 * np.abs(lam.real).max()]
        tmax = 10.0 / nz.min()
        times = np.linspace(0.0, tmax, 6001)
        V = lb.Propagator(s.L).apply(lb.vec(rho0), times)
        flux = (row @ V).real
        ip = int(np.argmax(flux))
        target = flux[ip] / math.e
        tau = NAN
        # radiative sector: vec(|a><b|) with a, b in {|X,0>, |G,1>} (the Lindbladian is closed on it)
        ia, ib = 1 * nc + 0, 0 * nc + 1

        def _idx(i, j):
            E = np.zeros((dim, dim), dtype=complex)
            E[i, j] = 1.0
            return int(np.nonzero(lb.vec(E))[0][0])
        sec = [_idx(i, j) for i in (ia, ib) for j in (ia, ib)]
        lam_sec = np.linalg.eigvals(s.L[np.ix_(sec, sec)])
        slow = float(np.min(np.abs(lam_sec.real)))
        tau_slow = 1.0 / slow if slow > 0 else NAN
        for i in range(ip, len(times) - 1):
            if flux[i + 1] <= target <= flux[i]:
                f = (flux[i] - target) / (flux[i] - flux[i + 1])
                tau = times[i] + f * (times[i + 1] - times[i]) - times[ip]
                break
        P = lb._stationary_projector(s.L)
        Q = np.eye(s.L.shape[0]) - P
        x = solve(-(s.L - P), Q @ lb.vec(rho0))
        n_cav = float(kap * (lb.vec((a.conj().T @ a).T) @ x).real)
        n_leak = float(gl * (lb.vec(s.ops["P_X"].T) @ x).real)
        return tau_slow * 1e3, n_cav, n_leak, tau * 1e3

    def run(self, plan, em, rt):
        res = self._run(plan, em, rt)
        pt = _p_trion(plan)
        if pt > 0.0 and plan.pulses:
            # FX4: dark-trion branch as per-pulse post-processing on the JC result (lifetime, HOM unchanged; photons,
            # inversion and preparation x (1 - p); g2 / (1 - p)); the cavity occupation at the window end scales too [E]
            area = math.sqrt(sum(float(q.area_rad) ** 2 for q in plan.pulses))
            p_t = _polaron.trion_probability(pt, area)
            res = _polaron.apply_trion(res, p_t)
            if p_t > 0.0 and "n_cav_end" in res.extra:
                res.extra["n_cav_end"] = float(res.extra["n_cav_end"]) * (1.0 - p_t)
        return res

    def _run(self, plan, em, rt):
        if plan.family == "coherent" and plan.timing == "pulsed":
            return self._run_coherent_pulsed(plan, em, rt)
        if plan.family == "coherent":
            return self._run_cw(plan, em, rt, coherent=True)
        if plan.timing == "cw":
            return self._run_cw(plan, em, rt, coherent=False)
        return self._run_incoherent_pulsed(plan, em, rt)

    # ---- R2 FIX-1a (audit 4): every supported drive has a strong-coupling path -----------------------------------------
    def _cw_grid(self, s, rt):
        """Delay grid of a cw g2(tau): out to ~12 slowest decay times, fine enough to resolve the vacuum-Rabi oscillation."""
        lam = np.linalg.eigvals(s.L)
        nz = np.abs(lam.real)[np.abs(lam.real) > 1e-9 * np.abs(lam.real).max()]
        tmax = 12.0 / max(float(nz.min()), 1e-9)
        fast = max(float(rt.g_ns), float(rt.kappa_ns) / 2.0, float(np.abs(lam.imag).max()))
        n = int(min(max(401, math.ceil(tmax * fast / 0.25)), 6001))
        return np.linspace(0.0, tmax, n)

    @staticmethod
    def spectrum(s, rho, op, omegas):
        """Incoherent emission spectrum of the field `op` of the stationary state `rho`:
        S(w) = 2 Re int_0^inf e^{-i w tau} <da^+(tau) da(0)> dtau, da = op - <op>, via the Liouvillian resolvent
        -(L - i w)^-1 (Lax / QRT) [DR]. The elastically scattered (coherent) part is excluded by the subtraction."""
        L, dim = s.L, s.dim
        mean = lb.expect(op, rho)
        da = op - mean * np.eye(dim)
        v = lb.vec(da @ rho)
        row = lb.vec(da.conj().T.T)                    # Tr[da^+ X] = vec((da^+)^T) . vec(X)
        eye = np.eye(L.shape[0])
        out = np.empty(len(omegas))
        for i, w in enumerate(omegas):
            x = -np.linalg.solve(L - 1j * float(w) * eye, v)
            out[i] = 2.0 * float((row @ x).real)
        return out

    def _spectrum_series(self, s, rho, rt):
        span = 2.5 * max(float(rt.g_ns), float(rt.kappa_ns) / 2.0, 1e-3)
        w = np.linspace(-span, span, 801) + 1e-7
        return {"x": [float(v) for v in w], "y": [float(v) for v in self.spectrum(s, rho, s.ops["a"], w)],
                "x_unit": "rad/ns (cavity-output spectrum relative to the dot frame)"}

    def _check_fock(self, s, rho, reasons):
        nc = self.FOCK
        n_op = s.ops["a"].conj().T @ s.ops["a"]
        top = np.real(np.diag(n_op)) > (nc - 1.5)
        pop_top = float(np.real(np.diag(rho))[top].sum())
        if pop_top > 1e-2:
            reasons.append(f"Fock cutoff {nc}: population of the top Fock level {pop_top:.3g} (> 1e-2): the cavity occupation "
                           "is truncated and the numbers are unreliable")
        return pop_top

    def _run_cw(self, plan, em, rt, coherent):
        """cw steady state of dot + cavity: coherent drive (Mollow triplet in a cavity) or an incoherent pump G -> X at
        r on the emitter (above_band, p_shell, electrical cw; S-Z form, pump term r D[sigma^+]). Both photon channels
        (leaky Gamma_bg and cavity output kappa) are counted; g2(tau) by the QRT; flux = total photon flux [DR]."""
        kap, Gbg = float(rt.kappa_ns), float(rt.Gamma_bg_ns)
        build = self._builder(plan, em, rt)
        tau_ps, _nc, _nl, tau_1e = self.lifetime(plan, em, rt)
        g_pop = 1e3 / float(tau_ps) if float(tau_ps) > 0 else float(rt.Gamma_ns) + float(rt.k_nr_ns)    # JC total decay rate 1/ns
        g_weak = float(rt.Gamma_ns) + float(rt.k_nr_ns)
        scale_notes = []
        if coherent:
            omega = float(plan.omega_cw_ns)
            if plan.inputs.get("omega_source") == "saturation_s" and g_weak > 0.0:
                # the planner set Omega = sqrt(s T1^-1 T2^-1) with the WEAK-coupling Gamma (Purcell); in the JC regime that rate is
                # meaningless, so the saturation parameter s = Omega^2 T1 T2 is re-applied with T1 = tau_jc and T2 = 2/(1/T1 + gamma*) [A]
                s_par = float(plan.inputs["saturation_s"])
                omega = math.sqrt(s_par * g_pop * 0.5 * (g_pop + float(rt.gamma_star_ns)))
                scale_notes.append(f"saturation parameter s = {s_par:.4g} applied with the JC decay rate 1/tau = {g_pop:.4g}/ns")
            s = build(omega)
            r = 0.0
        else:
            r, _ = _interp_pump(plan, rt)
            if g_weak > 0.0:
                # cw pump strength is "P/P_sat or mu in units of the total X decay rate" [A]: re-expressed with the JC decay rate
                # (the planner used the weak-coupling Gamma, which overstates the decay rate by F_P in the strong regime)
                r = r / g_weak * g_pop
                scale_notes.append(f"cw pump strength r = {r:.4g}/ns = (pump in units of the total decay rate) x the JC rate {g_pop:.4g}/ns")
            s = build(0.0, r_ns=r, k_extra=float(plan.k_esc_on_ns))
        rho = lb.steady_state(s.L, s.dim)
        a, cX = s.ops["a"], s.ops["c_X"]
        channels = [(cX, Gbg), (a, kap)]
        taus = self._cw_grid(s, rt)
        g2, _, flux = lb.g2_tau(s, taus, channels=channels, rho_ss=rho)
        n_c = float(lb.expect(a.conj().T @ a, rho).real)
        p_x = float(lb.expect(s.ops["P_X"], rho).real)
        n_cav_rate, n_leak_rate = kap * n_c, Gbg * p_x
        tot = n_cav_rate + n_leak_rate
        b0 = float(rt.beta0)
        reasons = ["cw: no per-pulse photon number or HOM overlap"] + scale_notes
        pop_top = self._check_fock(s, rho, reasons)
        extra = {"flux_ns": float(flux), "n_cav_ss": n_c,
                 "beta": ((n_cav_rate + b0 * n_leak_rate) / tot) if tot > 0 else NAN,
                 "tau_jc_ps": float(tau_ps), "tau_1e_ps": float(tau_1e), "fock_top_pop": pop_top,
                 "series": {"g2_tau": _series_xy(taus, g2), "spectrum_cav": self._spectrum_series(s, rho, rt)},
                 "notes": reasons}
        bg = float(plan.inputs.get("bg_frac", 0.0) or 0.0)
        g2_0 = float(g2[0])
        if coherent:
            cx = complex(lb.expect(cX, rho))
            extra["coherent_fraction"] = float(abs(cx) ** 2 / p_x) if p_x > 0 else NAN
            reasons.append("coherent cw drive with a cavity mode: Lindblad steady state (Fock cutoff 3), both photon channels counted")
        else:
            extra["tau_dip_ns"] = 1.0 / (1.0 / (float(tau_ps) * 1e-3) + r)       # dip recovery 1/(1/tau_jc + r) [A]
            extra["g2_signal"] = g2_0
            if bg > 0.0:
                g2_0 = 1.0 - (1.0 - bg) ** 2 * (1.0 - g2_0)
                reasons.append("bg_frac: g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac [A]")
            F_dot, fano_notes = _pump_fano(plan)
            reasons += fano_notes
            if F_dot is not None:
                extra["F_dot"] = F_dot
            reasons.append("incoherent cw pump on the emitter (S-Z pump term r D[sigma^+]); the pump broadening is in the Lindbladian")
        ok = pop_top < 0.05
        return SolverResult(g2_0=g2_0, n_ph=NAN, I_single=NAN, inversion=p_x, extra=extra, valid=ok,
                            reasons=() if ok else ("jc: Fock cutoff exceeded (top-level population > 5%)",))

    def _run_incoherent_pulsed(self, plan, em, rt):
        """Rectangular incoherent pump window on the emitter in the cavity: periodic steady state of the counting hierarchy
        (cavity-output and leaky photons both counted), HOM = cavity-output overlap of one exciton loaded at t = 0 [A]."""
        r, tau_on = _interp_pump(plan, rt)
        period = 1e3 / float(plan.f_rep_MHz)
        kap, Gbg = float(rt.kappa_ns), float(rt.Gamma_bg_ns)
        build = self._builder(plan, em, rt)
        s_off, s_on = build(0.0), build(0.0, r_ns=r)
        J = Gbg * lb.jump_super(s_off.ops["c_X"]) + kap * lb.jump_super(s_off.ops["a"])
        out = lb.pulsed_counting(s_on.L, s_off.L, J, s_off.dim, tau_on, period - tau_on, gate_ns=plan.gate_ns)
        tau_ps, n_cav, n_leak, tau_1e = self.lifetime(plan, em, rt)
        dim, nc = s_off.dim, self.FOCK
        rho0 = np.zeros((dim, dim), dtype=complex)
        rho0[1 * nc + 0, 1 * nc + 0] = 1.0                               # |X, 0 photons>
        reasons = ["HOM assumes instantaneous loading at t = 0 (window re-excitation ignored) [A]",
                   "HOM is the cavity-output photon overlap (leaky-mode photons excluded)"]
        if float(plan.relax_ns) > 0.0:
            reasons.append("capture / relaxation jitter (tau_cap, tau_relax) is not modelled by the jc solver")
        try:
            I = float(lb.indistinguishability(s_off.L, dim, rho0, np.sqrt(kap) * s_off.ops["a"]))
        except Exception as ex:                                          # noqa: BLE001
            I = NAN
            reasons.append(f"jc HOM failed: {type(ex).__name__}: {ex}")
        g2_sig = float(out["g2"])
        bg = float(plan.inputs.get("bg_frac", 0.0) or 0.0)
        g2 = 1.0 - (1.0 - bg) ** 2 * (1.0 - g2_sig) if bg > 0.0 else g2_sig
        if bg > 0.0:
            reasons.append("bg_frac: g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac; HOM and n_ph are the signal-only values [A]")
        F_dot, fano_notes = _pump_fano(plan)
        reasons += fano_notes
        if not out["converged"]:
            reasons.append("periodic steady state did not converge")
        tot = n_cav + n_leak
        b0 = float(rt.beta0)
        extra = {"p_loaded": 1.0 - math.exp(-r * tau_on), "g2_signal": g2_sig, "tau_jc_ps": float(tau_ps),
                 "tau_1e_ps": float(tau_1e), "n_cav_inst": n_cav, "n_leak_inst": n_leak,
                 "beta": ((n_cav + b0 * n_leak) / tot) if tot > 0 else NAN, "notes": reasons,
                 **({} if F_dot is None else {"F_dot": F_dot})}
        return SolverResult(g2_0=float(g2), n_ph=float(out["mean_counts"]), I_single=I, inversion=None, extra=extra,
                            valid=bool(out["converged"]),
                            reasons=() if out["converged"] else ("periodic steady state did not converge",))

    def _run_coherent_pulsed(self, plan, em, rt):
        p = plan.pulses[0]
        dur = p.fwhm_ps * 1e-3
        build = self._builder(plan, em, rt)
        tau_ps, n_cav, n_leak, tau_1e_ps = self.lifetime(plan, em, rt)
        s0, sd = build(0.0), build(1.0)
        kap = float(rt.kappa_ns)
        J = float(rt.Gamma_bg_ns) * lb.jump_super(s0.ops["c_X"]) + kap * lb.jump_super(s0.ops["a"])
        pc = rs._Pieces(L0=s0.L, L1=lb.hamiltonian_super(0.5 * sd.ops["sigma_x_drive"]), J=J,
                        A=np.sqrt(kap) * s0.ops["a"], ops=s0.ops, dim=s0.dim)
        st = rs._stats(pc, p.shape, p.area_rad, dur, _f_rep(plan), None, True)
        pop_X = float(lb.expect(s0.ops["P_X"], st["rho_post_pulse"]).real)
        # FX2b: excitation probability before decay = P_X + <a+a> at the window end + photons already emitted
        rho_w, emitted = _window_counts(pc, p.shape, p.area_rad, dur)
        n_cav_end = float(lb.expect(s0.ops["a"].conj().T @ s0.ops["a"], rho_w).real)
        n_exc = float(lb.expect(s0.ops["P_X"], rho_w).real) + n_cav_end + emitted
        dX = np.real(np.diag(s0.ops["P_X"]))
        dn = np.real(np.diag(s0.ops["a"].conj().T @ s0.ops["a"]))
        keep = tuple(int(i) for i in range(len(dX)) if dX[i] < 1e-12 and dn[i] < 1e-12)     # |G, 0 photons>
        eta_prep = 1.0 - _window_p0(pc, p.shape, p.area_rad, dur, keep)
        reasons = []
        try:
            I = float(rs.pulsed_indistinguishability(_CavOut(build, kap), p.shape, p.area_rad, dur,
                                                     tol=_HOM_TOL, max_doublings=_HOM_DOUBLINGS))
            reasons.append("HOM is the cavity-output photon overlap (leaky-mode photons excluded)")
        except Exception as ex:                                   # noqa: BLE001
            I = NAN
            reasons.append(f"jc HOM failed: {type(ex).__name__}: {ex}")
        tot = n_cav + n_leak
        b0 = float(rt.beta0)
        return SolverResult(g2_0=float(st["g2"]), n_ph=float(st["mean_photons"]), I_single=I,
                            inversion=pop_X,
                            extra={"tau_jc_ps": float(tau_ps), "tau_1e_ps": float(tau_1e_ps), "n_cav_inst": n_cav, "n_leak_inst": n_leak,
                                   "eta_prep": eta_prep, "n_exc": n_exc, "photons_in_window": emitted, "n_cav_end": n_cav_end,
                                   "P_X_end": pop_X,
                                   "beta": ((n_cav + b0 * n_leak) / tot) if tot > 0 else NAN, "notes": reasons},
                            reasons=())

    def eta_seconds(self, plan, n_points):
        return (6.0 if plan.family == "coherent" and plan.timing == "pulsed" else 1.0) * max(1, int(n_points))


# ---------------------------------------------------------------- registry + select
for _s in (RateSolver(), LindbladSolver(), IncoherentLindbladSolver(), JCSolver()):
    register(_s)
from . import polaron as _polaron  # noqa: E402,F401  (registers the "polaron" solver, P9)


def select(plan, em, rt, requested="auto", want_hom=True):
    """-> (solver | None, issues). Explicit solver wins if it supports the plan, else an issue (None).
    auto: phonon-needing plans go to a registered solver with handles_phonons (polaron; if it cannot run
    the plan the issues start with "needs P9 polaron solver" and list its reason); coherent pulsed resonant plans go to jc when rt.jc_recommended (strong
    coupling or 4g/kappa > 0.5), else lindblad; incoherent / electrical pulsed plans go to
    lindblad_incoherent (rate when want_hom is False or lindblad_incoherent cannot run, e.g. cap-2);
    cw incoherent to rate."""
    names = _kinds("solver")
    issues = []
    if requested not in (None, "auto"):
        if requested not in names:
            return None, [f"solver {requested!r} is not registered (registered: {names})"
                          + ("; needs P9 polaron solver" if requested == "polaron" else "")]
        s = _get("solver", requested)
        ok, why = s.supports(plan, em, rt)
        return (s, issues) if ok else (None, [f"solver {requested!r} cannot run this plan: {why}"])
    if "phonons" in plan.needs:
        for n in names:
            s = _get("solver", n)
            if getattr(s, "handles_phonons", False):
                ok, why = s.supports(plan, em, rt)
                if ok:
                    return s, issues
                issues.append(f"{n}: {why}")
        return None, ["needs P9 polaron solver (phonon-assisted drive)"] + issues
    strong = rt.regime == "strong"
    if plan.family == "coherent":
        if "cascade" in plan.needs:
            order = ["lindblad"]                 # TPE: the jc solver has no biexciton (the engine states that when strong)
        elif plan.timing == "pulsed":
            order = ["jc", "lindblad"] if (rt.jc_recommended or strong) else ["lindblad"]
        else:                                    # cw resonant: the jc path only where the weak-coupling model fails
            order = ["jc", "lindblad"] if strong else ["lindblad"]
    elif plan.timing == "pulsed":
        order = ["lindblad_incoherent", "rate"] if want_hom else ["rate", "lindblad_incoherent"]
        if strong:
            order = ["jc"] + order
    else:
        order = ["jc", "rate"] if strong else ["rate"]
    for n in order:
        if n not in names:
            continue
        s = _get("solver", n)
        ok, why = s.supports(plan, em, rt)
        if ok:
            if plan.family == "coherent" and "jc" in order and n != "jc" and any(
                    i == f"jc: {JC_EID_REASON}" for i in issues):
                return _JCRegimeFallback(s), issues
            return s, issues
        issues.append(f"{n}: {why}")
    return None, issues or ["no solver registered for this plan"]
