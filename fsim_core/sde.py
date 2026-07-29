"""Module D -- stochastic drive engine (update-plan Task 4).

A Langevin (Ito-Euler) simulation of the three-level QD-laser rate equations
that Zhao et al. (arXiv:2309.09703 / PRResearch 6, L032021 (2024), Eqs. (1)-
(4)) [V, read -- see cards/zhao.yaml] use to connect pump-current statistics
(normal, shot-noise-limited vs "quiet", sub-Poissonian high-impedance-drive)
to the output photon-number g2(0). The paper's Eqs. (1)-(4) give only the
SKELETON (which reservoirs feed which, which rates the photon field):

    dN_RS/dt = R_pump + R_ES->RS - R_RS->ES - N_RS/tau_RS_spon
    dN_ES/dt = R_RS->ES + R_GS->ES - R_ES->RS - R_ES->GS - N_ES/tau_ES_spon
    dN_GS/dt = R_ES->GS - R_GS->ES - R_stim - N_GS/tau_GS_spon
    dS/dt    = R_stim - S/tau_p + beta_sp*N_GS/tau_GS_spon

The exact closed forms for R_RS<->ES, R_ES<->GS, R_stim and every numeric
rate constant are in the paper's Supplementary Material, which is NOT
available in this workspace (see cards/zhao.yaml meta.source). Every
parameter below is therefore drawn from the standard 1.3 um InAs/GaAs
QD-laser rate-equation literature class and TUNED within the stated [E]
bounds to reproduce cards/zhao.yaml's g2_vs_pump validation points -- this is
the plan's Tier-2 fit, not a first-principles prediction. See
QDLaserParams's field comments for the literature band behind each entry,
and scripts/run_zhao_fit.py's module docstring for the final tuned table.

Rate closed forms used here [E, standard QD-laser rate-equation form]:

    rho_ES = N_ES / (g_ES * N_dots),  rho_GS = N_GS / (g_GS * N_dots)   (occupations, clipped to [0,1])
    R_RS->ES = N_RS / tau_cap   * (1 - rho_ES)          (Pauli-blocked capture)
    R_ES->GS = N_ES / tau_relax * (1 - rho_GS)          (Pauli-blocked relaxation)
    R_ES->RS = N_ES / tau_esc_ES,      tau_esc_ES = tau_cap   * exp(+E_a_ES/kT)   (Arrhenius escape, unblocked -- RS is an unbound continuum)
    R_GS->ES = N_GS / tau_esc_GS * (1 - rho_ES),  tau_esc_GS = tau_relax * exp(+E_a_GS/kT)   (Arrhenius escape, blocked by the ES it lands in)
    R_stim   = g0 * (2*rho_GS - 1) * S                  (linear gain in GS inversion; forward gain g0*rho_GS*S minus absorption g0*(1-rho_GS)*S)
    R_loss   = S / tau_p
    R_spon   = beta_sp * N_GS / tau_GS_spon             (spontaneous coupling into the lasing mode)

Stochastic part (the point of this module): every elementary process above is
an independent (or, where it moves carriers between two reservoirs, a SHARED)
Poisson point process; each gets an Ito-Euler diffusion-approximation noise
term of variance = rate * dt (shot noise), added once per elementary process
and applied with the correct sign to every reservoir it touches (so e.g. a
single N_RS->N_ES capture "hop" carries ONE random draw, not two independent
ones -- the hop is a single event). R_stim is the one process expressed above
as a NET rate (forward gain minus absorption); its total shot-noise power is
the SUM of the two, g0*S*dt, not (2*rho_GS-1)*S*dt -- this matters most near
threshold (rho_GS ~ 0.5) where the net rate is small but both directions are
individually large (the well-known near-threshold noise peak). The pump term
alone carries the tunable Fano factor F_pump: its noise variance is
F_pump * R_pump * dt (F_pump=1 recovers Poisson/shot-noise "normal" pump;
F_pump ~ 0.05-0.1 approximates Zhao's high-impedance "quiet" current source,
per cards/zhao.yaml / fsim_core.loading's F9 fano_pump). State variables are
clipped to >= 0 after every step (reflecting boundary -- a standard, coarse
but adequate [E] numerical treatment for these SDEs at the operating points
used here, where the clip triggers rarely).

g2(0) extraction follows the paper's stated relation exactly: Var(S) =
(g2(0)-1)<S>^2 + <S>, i.e. g2(0) = 1 + (Var(S) - <S>)/<S>^2, applied to the
steady-state (post-burn-in) photon-number time series.

TUNING NEAR-MISS -- read before trusting the quiet-vs-normal separation:
the deterministic/threshold physics and the "normal" (F_pump=1) pump-noise
baseline tune cleanly to cards/zhao.yaml (g2 ~ 1.020 vs the published
1.0224 +- 0.003 at I=4 I_th -- see scripts/run_zhao_fit.py). Separating
"quiet" (F_pump ~ 0.05-0.1) from "normal" at the OUTPUT photon field,
however, was NOT achievable to the published magnitude (or sign -- quiet
does not cross below g2=1) with ANY combination of the [E]-bounded
parameters tried (12 tuning experiments spanning N_dots 500-30000, escape
activation energies up to 200/150 meV, gain-compression strengths 1e-3 to
3e-2, and F_pump values down to 0 -- i.e. even a PERFECTLY silent pump).
The mechanism, confirmed numerically (and independently re-verified in the
Opus v1.1 review via an A/B of the noise-coupling conventions): each hop
(capture, escape, relax, stim, loss) carries its OWN full-shot (F=1)
partition noise that is UNCORRELATED WITH THE PUMP -- note this is about
pump-correlation, NOT about the coupling convention: every hop's draw IS
correctly shared/anti-correlated between its source and sink reservoirs (one
dW per microscopic event, +/- into the coupled equations; see the noise
section below). Each hop dilutes whatever fraction of the pump's regulation survives
by roughly half (measured: N_RS's Fano factor moves from ~1.0 at F_pump=1 to
~0.5 at F_pump=0, matching a clean "half pump signal, half capture's own
noise" split), and by the third hop (GS -> S, via R_stim/R_loss, which are
THEMSELVES F=1 and of order R_pump in magnitude) the surviving pump-
attributable fraction of Var(S) is below 1% -- an order of magnitude too
small to move g2(0) by the published ~0.02-0.04. This is very likely why
the real device needs the paper's own tuned Supplement parameters (not
available here) and/or a frequency-resolved (not flat per-step) treatment
of the pump regulation bandwidth (Eq. 5's tau_te) to show the effect at the
size reported -- both out of reach of the literature-bounds-only Tier-2 fit
this module is scoped to. See scripts/run_zhao_fit.py's module docstring
and verify/verify_sde.py's ordering check for how this is honestly reported
downstream (the ordering check verifies the CORRECT, statistically
significant DIRECTION of the F_pump dependence rather than a g2=1 crossing
the model cannot reach).

eta_ext / thin_stats note (read before trusting "eta_ext moves g2 toward 1"):
binomial thinning of ANY photon-number random variable by an independent
efficiency eta -- m ~ Binomial(n, eta) -- gives Var(m) = eta^2 Var(n) +
eta(1-eta)<n> and <m> = eta<n>, and substituting into the g2 relation above
gives g2(m) = g2(n) EXACTLY (the well-known loss-invariance of g2(0); it is
why g2(0) is the metric of choice for a lossy single-photon-source
measurement in the first place). thin_stats() implements this exact,
textbook-correct transform. eta_ext therefore does NOT move the true g2(0);
what it DOES do, and what verify_sde.py's eta_ext check actually verifies, is
that the THINNED FANO FACTOR obeys the F8b lemma (fsim_core.loading.
f8b_thin_fano) exactly -- F_thin = 1 - eta*(1-F_source) -> 1 as eta -> 0 --
while g2, which further normalizes by the (also-shrinking) mean squared,
stays invariant. This is a deliberate, documented departure from a literal
reading of "thinning moves g2 toward 1": that phrasing is true of the raw
Fano factor, not of g2(0) itself, and this module reports the physically
correct result rather than fudging the estimator to produce a movement that
would not survive a textbook check.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from .spectral import KB  # meV/K

BURNIN_FRAC = 0.1  # fraction of simulate()/g2_vs_pump() steps discarded as SDE transient


# ============================================================================ params

@dataclass
class QDLaserParams:
    """Tuned literature-class three-level QD-laser parameter set. Every field
    is [E] (estimated from the standard 1.3 um InAs/GaAs QD-laser rate-
    equation literature -- the paper's own numbers are in its unavailable
    Supplement) unless the comment says otherwise. Defaults are the Tier-2
    fit to cards/zhao.yaml (see scripts/run_zhao_fit.py docstring for the
    tuned-parameter table)."""

    # ---- geometry / degeneracy
    N_dots: float = 2.0e4
    # [E] number of QDs coupled to the lasing mode. Not itself a published
    # number for this device; it sets the absolute carrier/photon-number
    # scale (hence the *magnitude* of the g2-1 excess, which scales roughly
    # as 1/<S>) and is the single most important TUNED knob for matching
    # cards/zhao.yaml's normal-pump +0.0224 offset (this value gives
    # g2 = 1.020 at F_pump=1, I=4 I_th -- see scripts/run_zhao_fit.py's
    # tuned-parameter table). Literature single-mode QD DFB active regions
    # hold 1e3-1e5 dots under the ridge.
    g_ES: float = 4.0     # [E] excited-state degeneracy (2 spin x 2 orbital) -- standard InAs/GaAs QD level count
    g_GS: float = 2.0     # [E] ground-state degeneracy (spin only)

    # ---- capture / relaxation times
    tau_cap: float = 3.0e-12     # [E] RS->ES capture time; literature band 1-10 ps
    tau_relax: float = 3.0e-12   # [E] ES->GS relaxation time; literature band 1-10 ps

    # ---- spontaneous lifetimes
    tau_RS_spon: float = 1.0e-9    # [E] literature band 0.5-2 ns
    tau_ES_spon: float = 1.0e-9    # [E] literature band 0.5-2 ns
    tau_GS_spon: float = 0.8e-9    # [E] literature band 0.5-2 ns

    # ---- thermal escape activation energies (Arrhenius suppression at T)
    E_a_ES_meV: float = 200.0
    # [E] RS-ES separation. Pushed to the deep end of the literature spread
    # for QD confinement (most-cited band is 40-100 meV for shallow/small
    # dots; deep, low-strain InAs/GaAs dots have been reported with
    # WL-to-ES separations of 150-250 meV) -- TUNED here, deliberately near
    # the top of that wider band, to suppress ES->RS thermal backfilling at
    # T=300K as much as the literature honestly allows. This was necessary
    # but NOT sufficient to give F_pump a measurable effect on the output
    # g2(0); see the "near-miss" note in this module's docstring and
    # scripts/run_zhao_fit.py for the full account.
    E_a_GS_meV: float = 150.0   # [E] ES-GS separation; same reasoning/band as E_a_ES_meV, kept below it (GS is the shallower-bound level)

    # ---- photon field / gain
    tau_p: float = 2.0e-12
    # [E] photon lifetime. Zhao et al. state the QD laser's photon lifetime
    # is "~1 ps" [V, read, qualitative]; literature band for these cavities
    # is 1-5 ps, tuned within that band.
    beta_sp: float = 4.0e-4    # [E] spontaneous coupling into the lasing mode; literature band 1e-4 - 1e-3
    g0: float = 1.4e12
    # [E] modal gain coefficient (1/s). Tuned so g0*tau_p > 1 (lasing is
    # reachable) and the deterministic threshold sits at a physically
    # reasonable GS inversion rho_GS_th ~ 0.6-0.8 (not squeezed against the
    # rho_GS -> 0.5 transparency point nor against saturation at 1).
    eps_gain: float = 2.5e-3
    # [E] nonlinear gain-compression coefficient: R_stim uses g0*(2rho_GS-1)*S
    # / (1 + eps_gain*S) rather than the bare linear form. Standard addition
    # in semiconductor-laser rate-equation models (spectral hole burning /
    # carrier heating); without it this multi-reservoir model's relaxation
    # oscillation is only very weakly damped (large-amplitude, near-undamped
    # ringing in S(t) even far above threshold -- verified numerically during
    # tuning) and Var(S) is dominated by that ringing rather than by the
    # pump-statistics-dependent noise floor the Zhao fit needs. Literature
    # band for QD lasers: 1e-3 - 1e-1 (units: 1/photon-number, so the natural
    # scale depends on the S normalization -- tuned here against this
    # model's own S scale, not taken directly off a cited number).

    # ---- pump axis (NOT SI amps -- see module docstring / find_threshold)
    I_ref: float = 3.2490e13
    # [E] reference pump rate (carriers/s) defining "1.0" on the I_over_Ith
    # axis used throughout this module: R_pump = I_over_Ith * I_ref.
    # Calibrated (one-off, via find_threshold) so that find_threshold(params)
    # returns ~1.0 for the DEFAULT parameter set -- i.e. "I_over_Ith" can be
    # read directly as I/I_th for these defaults. find_threshold() itself is
    # fully generic and makes no such assumption for other parameter sets.

    # ---- drive statistics
    F_pump: float = 1.0
    # Fano factor (Var/mean) of the pump point process. 1.0 = normal
    # (Poisson, shot-noise-limited) drive. Zhao's quiet high-impedance
    # current source [V regime] is approximated by F_pump ~ 0.05-0.1 [E
    # magnitude -- the paper reports the resulting *photon* squeezing
    # (3.1 dB, cards/zhao.yaml) and the qualitative high-impedance argument
    # (F9 in fsim_core.loading), not a measured pump Fano factor].
    eta_ext: float = 0.2
    # [E] external collection/outcoupling efficiency, used only by
    # thin_stats() for the F8b binomial-thinning consistency check (see
    # module docstring "eta_ext / thin_stats note" -- it does not move the
    # true g2(0)).

    T: float = 300.0   # K [V -- Zhao operate at 20 C +- 0.005 C]
    seed: int = 0

    def kT_meV(self) -> float:
        return KB * self.T


# ============================================================================ rates

def _rates(N_RS, N_ES, N_GS, S, R_pump, p: QDLaserParams):
    """Deterministic rate dict for state (N_RS, N_ES, N_GS, S) at pump R_pump.
    Accepts scalars or same-shape arrays (vectorized across parallel runs)."""
    kT = p.kT_meV()
    rho_ES = np.clip(N_ES / (p.g_ES * p.N_dots), 0.0, 1.0)
    rho_GS = np.clip(N_GS / (p.g_GS * p.N_dots), 0.0, 1.0)

    tau_esc_ES = p.tau_cap * np.exp(p.E_a_ES_meV / kT)
    tau_esc_GS = p.tau_relax * np.exp(p.E_a_GS_meV / kT)

    R_cap = np.maximum(N_RS, 0.0) / p.tau_cap * (1.0 - rho_ES)
    R_esc1 = np.maximum(N_ES, 0.0) / tau_esc_ES
    R_rel = np.maximum(N_ES, 0.0) / p.tau_relax * (1.0 - rho_GS)
    R_esc2 = np.maximum(N_GS, 0.0) / tau_esc_GS * (1.0 - rho_ES)
    R_RS_spon = np.maximum(N_RS, 0.0) / p.tau_RS_spon
    R_ES_spon = np.maximum(N_ES, 0.0) / p.tau_ES_spon
    R_GS_spon = np.maximum(N_GS, 0.0) / p.tau_GS_spon
    S_pos = np.maximum(S, 0.0)
    gain_sat = 1.0 / (1.0 + p.eps_gain * S_pos)  # nonlinear gain compression (RO damping)
    R_stim = p.g0 * (2.0 * rho_GS - 1.0) * S_pos * gain_sat
    R_stim_total = p.g0 * S_pos * gain_sat  # forward + absorption, for the noise power
    R_loss = np.maximum(S, 0.0) / p.tau_p
    R_spon_mode = p.beta_sp * R_GS_spon

    return dict(rho_ES=rho_ES, rho_GS=rho_GS, R_cap=R_cap, R_esc1=R_esc1, R_rel=R_rel,
                R_esc2=R_esc2, R_RS_spon=R_RS_spon, R_ES_spon=R_ES_spon, R_GS_spon=R_GS_spon,
                R_stim=R_stim, R_stim_total=R_stim_total, R_loss=R_loss, R_spon_mode=R_spon_mode)


def _derivs(N_RS, N_ES, N_GS, S, R_pump, p: QDLaserParams, disable_stim=False):
    r = _rates(N_RS, N_ES, N_GS, S, R_pump, p)
    R_stim = 0.0 if disable_stim else r["R_stim"]
    dN_RS = R_pump + r["R_esc1"] - r["R_cap"] - r["R_RS_spon"]
    dN_ES = r["R_cap"] + r["R_esc2"] - r["R_esc1"] - r["R_rel"] - r["R_ES_spon"]
    dN_GS = r["R_rel"] - r["R_esc2"] - R_stim - r["R_GS_spon"]
    dS = R_stim - r["R_loss"] + r["R_spon_mode"]
    return dN_RS, dN_ES, dN_GS, dS


# ======================================================================= steady state

def _resid_scale(N_RS, N_ES, N_GS, S, R_pump, params):
    d = np.array(_derivs(N_RS, N_ES, N_GS, S, R_pump, params))
    scale = np.array([max(N_RS, 1.0), max(N_ES, 1.0), max(N_GS, 1.0), max(S, 1e-6)])
    return float(np.max(np.abs(d) / scale * max(params.tau_p, params.tau_cap)))


def _steady_state_from_guess(params: QDLaserParams, R_pump: float, y0: np.ndarray,
                              disable_stim: bool = False) -> tuple:
    """Newton-polish only (no relax) -- fast path used for continuation along
    an I-grid, where the previous grid point's solution is an excellent
    initial guess. Returns (y_star, rel_resid)."""
    def resid(y):
        return np.array(_derivs(y[0], y[1], y[2], y[3], R_pump, params, disable_stim=disable_stim))
    y_star, info, ier, msg = fsolve(resid, y0, full_output=True, xtol=1e-13)
    y_star = np.maximum(y_star, 0.0)
    rel_resid = _resid_scale(*y_star, R_pump, params)
    return y_star, rel_resid


def steady_state(params: QDLaserParams, I_over_Ith: float, disable_stim: bool = False) -> dict:
    """Deterministic fixed point at pump R_pump = I_over_Ith * params.I_ref.
    Method: relax the noise-free ODE forward with a stiff solver to a good
    initial guess, then Newton-polish with fsolve to a tight root (residual
    typically << 1e-9 relative -- see verify_sde.py). Returns a dict with the
    fixed-point state, occupations, and the raw derivative residual.
    disable_stim=True gives the (different) fixed point of the no-lasing
    spontaneous-only system -- used as the initial condition when _integrate
    runs with disable_stim=True, so that run doesn't spend its whole
    trajectory relaxing away from the (very different) lasing fixed point."""
    R_pump = I_over_Ith * params.I_ref

    def rhs(t, y):
        return np.array(_derivs(y[0], y[1], y[2], y[3], R_pump, params, disable_stim=disable_stim))

    y0 = np.array([R_pump * params.tau_RS_spon * 0.1, 10.0, 10.0, 1.0])
    t_relax = 25.0 * params.tau_RS_spon
    sol = solve_ivp(rhs, (0.0, t_relax), y0, method="Radau", rtol=1e-9, atol=1e-6)
    y_relaxed = np.maximum(sol.y[:, -1], 0.0)

    y_star, rel_resid = _steady_state_from_guess(params, R_pump, y_relaxed, disable_stim=disable_stim)
    N_RS, N_ES, N_GS, S = y_star
    r = _rates(N_RS, N_ES, N_GS, S, R_pump, params)
    return dict(N_RS=float(N_RS), N_ES=float(N_ES), N_GS=float(N_GS), S=float(S),
                rho_ES=float(r["rho_ES"]), rho_GS=float(r["rho_GS"]),
                R_pump=float(R_pump), residual=rel_resid)


def _steady_state_curve(params: QDLaserParams, I_grid: np.ndarray) -> np.ndarray:
    """S(I) along I_grid (ascending) via continuation: each point's fsolve is
    warm-started from the previous point's converged state, with a single
    full relax+polish steady_state() only at the first point (or wherever a
    warm start fails to converge). Far cheaper than calling steady_state()
    independently at every grid point -- this is what keeps find_threshold
    and the tuning/scan workflow fast."""
    I_grid = np.sort(np.asarray(I_grid, dtype=float))
    S_out = np.empty_like(I_grid)
    ss0 = steady_state(params, I_grid[0])
    y_prev = np.array([ss0["N_RS"], ss0["N_ES"], ss0["N_GS"], ss0["S"]])
    S_out[0] = ss0["S"]
    for i in range(1, len(I_grid)):
        R_pump = I_grid[i] * params.I_ref
        y_star, rel_resid = _steady_state_from_guess(params, R_pump, y_prev)
        if not np.isfinite(rel_resid) or rel_resid > 1e-6:
            ss = steady_state(params, I_grid[i])
            y_star = np.array([ss["N_RS"], ss["N_ES"], ss["N_GS"], ss["S"]])
        y_prev = y_star
        S_out[i] = y_star[3]
    return S_out


# ============================================================================ threshold

def find_threshold(params: QDLaserParams, I_lo: float = 0.05, I_hi: float = 8.0,
                    n_grid: int = 80, s_lo_frac: float = 0.15, s_hi_frac: float = 0.75) -> float:
    """Operational threshold on the I_over_Ith axis, by the standard
    experimental L-I convention: linear-fit extrapolation. Above threshold,
    gain clamping pins rho_GS and S grows roughly linearly in I (S ~
    slope*(I-I_th)) until Pauli blocking of the ES reservoir saturates the
    whole cascade at high I; fit that line to the points whose S lies in
    [s_lo_frac, s_hi_frac] of S(I_hi) -- above the near-threshold curvature,
    below the high-I saturation shoulder -- and extrapolate to S=0. This is
    the same convention used to read I_th off a real L-I curve, and it is
    what "I = k * I_th" means when comparing to cards/zhao.yaml's I/I_th =
    4.0 operating point. (An alternative -- max curvature of log S vs log I
    -- was tried first and rejected: on this model's S(I) curve the two knee
    definitions disagree by a factor ~1.6, and the linear-extrapolation
    definition is the one that matches how I_th is actually read off a
    measured L-I curve, so it is the one used everywhere downstream. I_hi
    must be chosen inside the linear-clamped region and below the ES
    saturation shoulder; the default was picked for this module's tuned
    parameter set -- see scripts/run_zhao_fit.py.)"""
    I_grid = np.geomspace(I_lo, I_hi, n_grid)
    S_grid = _steady_state_curve(params, I_grid)
    S_ref = S_grid[-1]
    mask = (S_grid >= s_lo_frac * S_ref) & (S_grid <= s_hi_frac * S_ref)
    if mask.sum() < 4:
        mask = S_grid >= s_lo_frac * S_ref  # fall back: everything above the near-threshold knee
    I_fit, S_fit = I_grid[mask], S_grid[mask]
    A = np.vstack([I_fit, np.ones_like(I_fit)]).T
    slope, intercept = np.linalg.lstsq(A, S_fit, rcond=None)[0]
    if slope <= 0:
        raise RuntimeError("find_threshold: non-lasing slope in the fit region; widen I_hi")
    I_th = -intercept / slope
    return float(max(I_th, I_lo))


# ============================================================================ SDE core

def _integrate(params: QDLaserParams, I_over_Ith: float, t_end: float, dt: float,
                n_runs: int = 1, seed: Optional[int] = None, disable_stim: bool = False,
                F_pump: Optional[float] = None, keep_series: bool = True) -> dict:
    """Vectorized Ito-Euler integration of n_runs independent trajectories in
    lockstep (one Python loop over time steps, numpy arrays over runs -- this
    is what keeps verify_sde.py / run_zhao_fit.py inside their runtime
    budgets). Starts from the deterministic steady state (minimizes burn-in).
    Returns per-run summary statistics always; full (n_runs, n_steps) time
    series only if keep_series (set False for large n_runs to save memory)."""
    seed = params.seed if seed is None else seed
    F_pump = params.F_pump if F_pump is None else F_pump
    R_pump = I_over_Ith * params.I_ref
    n_steps = int(round(t_end / dt))
    rng = np.random.default_rng(seed)

    ss = steady_state(params, I_over_Ith, disable_stim=disable_stim)
    N_RS = np.full(n_runs, ss["N_RS"])
    N_ES = np.full(n_runs, ss["N_ES"])
    N_GS = np.full(n_runs, ss["N_GS"])
    S = np.full(n_runs, ss["S"])

    if keep_series:
        S_hist = np.empty((n_runs, n_steps))
        t_hist = np.arange(n_steps) * dt
        N_RS_hist = np.empty((n_runs, n_steps))
        N_ES_hist = np.empty((n_runs, n_steps))
        N_GS_hist = np.empty((n_runs, n_steps))
    sqdt = np.sqrt(dt)

    for i in range(n_steps):
        r = _rates(N_RS, N_ES, N_GS, S, R_pump, params)

        xi_pump, xi_cap, xi_esc1, xi_rel, xi_esc2 = rng.standard_normal((5, n_runs))
        xi_rsspon, xi_esspon, xi_gsspon, xi_stim, xi_loss, xi_spon = rng.standard_normal((6, n_runs))

        dW_pump = np.sqrt(np.maximum(F_pump, 0.0) * R_pump * dt) * xi_pump if F_pump > 0 else 0.0
        dW_cap = np.sqrt(r["R_cap"] * dt) * xi_cap
        dW_esc1 = np.sqrt(r["R_esc1"] * dt) * xi_esc1
        dW_rel = np.sqrt(r["R_rel"] * dt) * xi_rel
        dW_esc2 = np.sqrt(r["R_esc2"] * dt) * xi_esc2
        dW_rsspon = np.sqrt(r["R_RS_spon"] * dt) * xi_rsspon
        dW_esspon = np.sqrt(r["R_ES_spon"] * dt) * xi_esspon
        dW_gsspon = np.sqrt(r["R_GS_spon"] * dt) * xi_gsspon
        dW_loss = np.sqrt(r["R_loss"] * dt) * xi_loss
        dW_spon = np.sqrt(r["R_spon_mode"] * dt) * xi_spon

        if disable_stim:
            stim_drift = 0.0
            dW_stim = 0.0
        else:
            stim_drift = r["R_stim"] * dt
            dW_stim = np.sqrt(r["R_stim_total"] * dt) * xi_stim  # total (fwd+abs) shot power

        d_pump = R_pump * dt + dW_pump
        d_cap = r["R_cap"] * dt + dW_cap
        d_esc1 = r["R_esc1"] * dt + dW_esc1
        d_rel = r["R_rel"] * dt + dW_rel
        d_esc2 = r["R_esc2"] * dt + dW_esc2
        d_rsspon = r["R_RS_spon"] * dt + dW_rsspon
        d_esspon = r["R_ES_spon"] * dt + dW_esspon
        d_gsspon = r["R_GS_spon"] * dt + dW_gsspon
        d_loss = r["R_loss"] * dt + dW_loss
        d_spon = r["R_spon_mode"] * dt + dW_spon
        d_stim = stim_drift + dW_stim

        N_RS = N_RS + d_pump + d_esc1 - d_cap - d_rsspon
        N_ES = N_ES + d_cap + d_esc2 - d_esc1 - d_rel - d_esspon
        N_GS = N_GS + d_rel - d_esc2 - d_stim - d_gsspon
        S = S + d_stim - d_loss + d_spon

        N_RS = np.maximum(N_RS, 0.0)
        N_ES = np.maximum(N_ES, 0.0)
        N_GS = np.maximum(N_GS, 0.0)
        S = np.maximum(S, 0.0)

        if keep_series:
            N_RS_hist[:, i] = N_RS
            N_ES_hist[:, i] = N_ES
            N_GS_hist[:, i] = N_GS
            S_hist[:, i] = S

    n_burn = int(round(n_steps * BURNIN_FRAC))
    if keep_series:
        S_tail = S_hist[:, n_burn:]
    else:
        S_tail = None

    out = dict(n_steps=n_steps, n_burn=n_burn, dt=dt, R_pump=R_pump,
                steady_state=ss, F_pump=F_pump)
    if keep_series:
        out.update(t=t_hist, N_RS=N_RS_hist, N_ES=N_ES_hist, N_GS=N_GS_hist, S=S_hist)
        mean_S = S_tail.mean(axis=1)
        var_S = S_tail.var(axis=1, ddof=1)
        g2 = 1.0 + (var_S - mean_S) / np.maximum(mean_S, 1e-300) ** 2
        out.update(mean_S_per_run=mean_S, var_S_per_run=var_S, g2_per_run=g2)
    return out


def thin_stats(mean: float, var: float, eta: float) -> tuple:
    """Exact binomial thinning of a photon-number distribution's (mean, var)
    by independent detection efficiency eta: mean_thin = eta*mean,
    var_thin = eta^2*var + eta*(1-eta)*mean. g2 computed from (mean_thin,
    var_thin) is algebraically IDENTICAL to g2 from (mean, var) -- see module
    docstring. Returns (mean_thin, var_thin)."""
    mean_thin = eta * mean
    var_thin = eta**2 * var + eta * (1.0 - eta) * mean
    return mean_thin, var_thin


def g2_from_stats(mean: float, var: float) -> float:
    return 1.0 + (var - mean) / max(mean, 1e-300) ** 2


# ============================================================================ public API

def simulate(params: QDLaserParams, I_over_Ith: float, t_end: float, dt: float,
             seed: Optional[int] = None, disable_stim: bool = False,
             F_pump: Optional[float] = None) -> dict:
    """Single-trajectory Ito-Euler simulation. Returns a dict with the full
    time series (t, N_RS, N_ES, N_GS, S -- each length n_steps) plus derived
    scalars: mean_S, var_S, g2_0 (post-burn-in, BURNIN_FRAC discarded), and
    the eta_ext-thinned (mean_S_ext, var_S_ext, g2_0_ext) diagnostic (see
    module docstring -- g2_0_ext is expected to equal g2_0 within MC noise,
    NOT to move toward 1)."""
    r = _integrate(params, I_over_Ith, t_end, dt, n_runs=1, seed=seed,
                    disable_stim=disable_stim, F_pump=F_pump, keep_series=True)
    mean_S = float(r["mean_S_per_run"][0])
    var_S = float(r["var_S_per_run"][0])
    g2_0 = float(r["g2_per_run"][0])
    mean_ext, var_ext = thin_stats(mean_S, var_S, params.eta_ext)
    g2_0_ext = g2_from_stats(mean_ext, var_ext)
    return dict(t=r["t"], N_RS=r["N_RS"][0], N_ES=r["N_ES"][0], N_GS=r["N_GS"][0], S=r["S"][0],
                n_burn=r["n_burn"], R_pump=r["R_pump"], steady_state=r["steady_state"],
                mean_S=mean_S, var_S=var_S, g2_0=g2_0,
                mean_S_ext=mean_ext, var_S_ext=var_ext, g2_0_ext=g2_0_ext)


def g2_vs_pump(params: QDLaserParams, I_over_Ith: float, n_runs: int,
               t_end: float = 3.0e-9, dt: float = 1.0e-13,
               seed: Optional[int] = None, disable_stim: bool = False,
               F_pump: Optional[float] = None) -> dict:
    """n_runs independent trajectories (vectorized in one integration call)
    at fixed I_over_Ith; returns per-run g2 array plus the mean and standard
    error across runs -- the quantity compared against cards/zhao.yaml."""
    r = _integrate(params, I_over_Ith, t_end, dt, n_runs=n_runs, seed=seed,
                    disable_stim=disable_stim, F_pump=F_pump, keep_series=True)
    g2s = r["g2_per_run"]
    means = r["mean_S_per_run"]
    g2_mean = float(np.mean(g2s))
    g2_se = float(np.std(g2s, ddof=1) / np.sqrt(n_runs)) if n_runs > 1 else float("nan")
    return dict(g2_mean=g2_mean, g2_se=g2_se, g2_per_run=g2s, mean_S_per_run=means,
                mean_S=float(np.mean(means)), R_pump=r["R_pump"], n_runs=n_runs)
