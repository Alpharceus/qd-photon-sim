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

eta_prep (FX2b): the excitation probability BEFORE decay = P_X(end of the pulse window) + photons emitted
inside the window (+ <a+a>(end) for jc; P_XX(end) + XX photons for TPE), started from the ground state
(single-pulse preparation; nonradiative loss inside the window is not added back). The old end-of-window
population is kept as extra["P_X_end"] (informational; it includes the decay during the window).

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
from scipy.linalg import expm, solve

from .. import cw_g2
from .. import lindblad as lb
from .. import pulse_counting as pcnt
from .. import resonant_source as rs
from ..spectral import KB
from .registry import Solver, get as _get, kinds as _kinds, register
from .types import SolverResult

NAN = float("nan")
_UNSUPPORTED_NEEDS = {
    "phonons": "needs P9 polaron solver",
    "td_detuning": "time-dependent detuning (chirp) is not implemented",
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
    if abs(p.chirp_meV_ps) > 0.0:
        return "chirped pulses need time-dependent detuning (not implemented)"
    f = _f_rep(plan)
    if f is None or f <= 0:
        return "pulsed plan needs f_rep_MHz"
    if 1e3 / f <= _window_ns(p):
        return "repetition period must exceed the pulse window"
    return None


def _interp_pump(plan, rt):
    """(pump rate 1/ns, window ns) of an incoherent / electrical plan.
    Electrical plans carry only a DriveInterface: pulsed r = mu eta_capture / tau_on (mean loadings per
    window, tau_on = interface tau_eff_ns, default 1 ns); cw r = mu eta_capture Gamma (mean loadings per
    radiative lifetime) [A]."""
    if plan.pump_rate_ns > 0.0:
        return float(plan.pump_rate_ns), float(plan.tau_on_ns)
    itf = plan.interface
    if itf is not None:
        eff = float(itf.mu) * float(itf.eta_capture)
        if plan.timing == "pulsed":
            tau = float(plan.tau_on_ns) if plan.tau_on_ns > 0 else 1.0
            return eff / tau, tau
        return eff * float(rt.Gamma_ns), 0.0
    return 0.0, float(plan.tau_on_ns)


def _series_xy(x, y, unit_x="ns"):
    return {"x": [float(v) for v in x], "y": [float(v) for v in y], "x_unit": unit_x}


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
            if 1e3 / f <= tau:
                return False, "repetition period must exceed the pump window"
        return True, ""

    def _rates(self, plan, em, rt):
        gX = float(rt.Gamma_ns)
        gXX = float(em.gamma_xx_ratio) * float(rt.Gamma0_ns)    # XX is not Purcell-enhanced [A]
        kX = float(rt.k_nr_ns)
        kXX = kX * float(em.gamma_xx_ratio)
        p = float(plan.inputs.get("pump_ratio_xx", 1.0)) if "cascade" in plan.needs else 0.0
        return gX, gXX, kX, kXX, p

    def run(self, plan, em, rt):
        gX, gXX, kX, kXX, pr = self._rates(plan, em, rt)
        r, tau_on = _interp_pump(plan, rt)
        reasons = ["rate solver: capture/relaxation jitter and coherence are not modelled"]
        bg = float(plan.inputs.get("bg_frac", 0.0) or 0.0)
        bgf = (lambda g: 1.0 - (1.0 - bg) ** 2 * (1.0 - g)) if bg > 0.0 else (lambda g: g)
        if float(plan.inputs.get("p_refill", 0.0) or 0.0) > 0.0 or float(plan.inputs.get("k0_ns", 0.0) or 0.0) > 0.0:
            reasons.append("p_refill / thermal recapture are not modelled by the rate solver (use lindblad_incoherent)")
        if bg > 0.0:
            reasons.append("bg_frac applied to g2 as 1 - rho^2 (1 - g2), rho = 1 - bg_frac [A]")
        if plan.timing == "pulsed":
            period = 1e3 / float(plan.f_rep_MHz)
            out = pcnt.pulse_g2(r, gX, gXX, kX, kXX, 1.0, 0.0, tau_on, period - tau_on, pump_ratio=pr)
            return SolverResult(g2_0=float(bgf(out["g2"])), n_ph=float(out["mean_counts"]), I_single=NAN,
                                inversion=None,
                                extra={"p_loaded": 1.0 - math.exp(-r * tau_on),
                                       "converged": bool(out["converged"]), "notes": reasons},
                                valid=bool(out["converged"]),
                                reasons=() if out["converged"] else ("periodic steady state did not converge",))
        scale = max(gX + kX, 1e-12)
        taus = np.linspace(0.0, 12.0 / scale, 401)
        g2 = cw_g2.g2_cw(taus, r, gX, gXX, 1.0, 0.0, kX, kXX, pr)
        P = cw_g2.steady_state(r, gX, gXX, kX, kXX, pr)
        return SolverResult(g2_0=float(bgf(g2[0])), n_ph=NAN, I_single=NAN, inversion=None,
                            extra={"flux_ns": float(gX * P[1]),
                                   "series": {"g2_tau": _series_xy(taus, g2)}, "notes": reasons})

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
        bad = _bad_needs(plan, allowed=("cascade",))
        if bad:
            return False, "; ".join(bad)
        why = _pulsed_coherent_check(plan)
        if why:
            return False, why
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

    def run(self, plan, em, rt):
        if "cascade" in plan.needs:
            return self._tpe(plan, em, rt)
        if plan.timing == "cw":
            return self._cw(plan, em, rt)
        p = plan.pulses[0]
        dur = p.fwhm_ps * 1e-3
        tl = rs.TwoLevel(float(rt.Gamma_ns), k_ns=float(rt.k_nr_ns), gamma_star_ns=float(rt.gamma_star_ns),
                         delta_ns=_det_X_ns(plan))
        st = rs.pulsed_photon_statistics(tl, p.shape, p.area_rad, dur, _f_rep(plan))
        I = rs.pulsed_indistinguishability(tl, p.shape, p.area_rad, dur, tol=_HOM_TOL,
                                           max_doublings=_HOM_DOUBLINGS)
        rho_w, emitted = _window_counts(rs._pieces(tl), p.shape, p.area_rad, dur)
        return SolverResult(g2_0=float(st["g2"]), n_ph=float(st["mean_photons"]), I_single=float(I),
                            inversion=float(st["P_e_post"]),
                            extra={"det_X_meV": _det_X_meV(plan),
                                   "eta_prep": float(rho_w[1, 1].real) + emitted,
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
        rho_w, emitted_xx = _window_counts(rs._pieces_cascade(casc, ("XX",)), p.shape, area, dur)
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
                                   "eta_prep": float(rho_w[2, 2].real) + emitted_xx,
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
        if 1e3 / f <= tau:
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
        L_hom = lb.liouvillian(H0, c_ops_hom)
        c_ops_cnt = list(c_ops_hom)
        if k_esc > 0.0:
            c_ops_cnt.append(np.sqrt(k_esc) * e(iP, 0))            # recapture of thermally released carriers
        L_off = lb.liouvillian(H0, c_ops_cnt)
        L_on = lb.liouvillian(H0, c_ops_cnt + [c_pump])
        out = lb.pulsed_counting(L_on, L_off, J, dim, tau_on, period - tau_on)
        A = np.sqrt(Gam) * e(0, 1)
        rho0 = e(2, 2) if with_res else e(1, 1)
        I = lb.indistinguishability(L_hom, dim, rho0, A)
        g2_sig = float(out["g2"])
        rho_sig = 1.0 - bg
        g2 = 1.0 - rho_sig ** 2 * (1.0 - g2_sig) if bg > 0.0 else g2_sig   # uncorrelated background (Brouri form)
        reasons = ["HOM assumes instantaneous loading at t = 0 (window re-excitation ignored) [A]"]
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
                 "g2_signal": g2_sig, "p_refill": p_ref, "bg_frac": bg, "k_esc_ns": k_esc, "notes": reasons}
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


class JCSolver(Solver):
    name = "jc"
    label = "Jaynes-Cummings (dot + cavity)"
    summary = ("Dot coupled to one cavity mode (Fock cutoff 3), Lindblad with kappa, g. For 4g/kappa > 0.5 "
               "or strong coupling: lifetime, cavity-output photons per pulse, g2, HOM.")
    FOCK = 3

    def supports(self, plan, em, rt):
        if plan.family != "coherent" or plan.timing != "pulsed":
            return False, "jc: coherent pulsed drives only"
        bad = _bad_needs(plan)
        if bad:
            return False, "; ".join(bad)
        if "cascade" in plan.needs:
            return False, "jc: no biexciton"
        why = _pulsed_coherent_check(plan)
        if why:
            return False, why
        if rt.kappa_ns is None or rt.g_ns <= 0.0:
            return False, "jc needs a cavity mode (kappa and g)"
        return True, ""

    def _builder(self, plan, em, rt):
        gs = float(rt.gamma_star_ns)
        det = _det_X_ns(plan)
        cav = dict(g_ns=float(rt.g_ns), kappa_ns=float(rt.kappa_ns),
                   det_ns=det + float(rt.cav_minus_dot_ns), n_max=self.FOCK)

        def build(omega_ns=0.0):
            return lb.build_system(levels=2, gamma_X_ns=float(rt.Gamma_bg_ns), k_X=float(rt.k_nr_ns),
                                   deph=[(gs, (0.0, 1.0))] if gs > 0 else [], det_X_ns=det,
                                   omega_ns=omega_ns, cavity=cav)
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
        eta_prep = float(lb.expect(s0.ops["P_X"], rho_w).real) + n_cav_end + emitted
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
                                   "eta_prep": eta_prep, "photons_in_window": emitted, "n_cav_end": n_cav_end,
                                   "P_X_end": pop_X,
                                   "beta": ((n_cav + b0 * n_leak) / tot) if tot > 0 else NAN, "notes": reasons},
                            reasons=())

    def eta_seconds(self, plan, n_points):
        return 6.0 * max(1, int(n_points))


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
    if plan.family == "coherent":
        order = ["jc", "lindblad"] if (rt.jc_recommended or rt.regime == "strong") else ["lindblad"]
        if "cascade" in plan.needs or plan.timing != "pulsed":
            order = ["lindblad"]
    elif plan.timing == "pulsed":
        order = ["lindblad_incoherent", "rate"] if want_hom else ["rate", "lindblad_incoherent"]
    else:
        order = ["rate"]
    for n in order:
        if n not in names:
            continue
        s = _get("solver", n)
        ok, why = s.supports(plan, em, rt)
        if ok:
            return s, issues
        issues.append(f"{n}: {why}")
    return None, issues or ["no solver registered for this plan"]
