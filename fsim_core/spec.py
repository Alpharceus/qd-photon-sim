"""Module F -- Spec mode (inverse design / requirements solver).

WHAT THIS MODULE IS NOT: fsim_core has no measured/simulated cavity numbers
(kappa, F_P, G) and no measured thermal R_th -- finding the REQUIRED values of
those quantities is exactly the point of this module. Every forward module
(spectral/loading/integrator/cavity/thermal) answers "given a design, what
g2(0) does it reach?"; spec.py answers the inverse question: "given a target
g2(0) <= t at an operating temperature T_op, what must the cavity/filter/
background/mesa/growth DELIVER?"

HONESTY (read this before trusting a spec_sheet() number): every function here
inverts the SAME validated F-series chain used everywhere else in fsim_core
(F1/F1b/F2/F5/F6, the master ceiling, Module A thermal) -- it is not new
physics. A requirement it returns (a kappa_max, a V_tilde_req, a mesa_min...)
is a NECESSARY CONDITION under the model: the design must clear that bar to
reach the target, given the class-proxy Gamma(T)/retention inputs FSIM is
carrying (all inherit [A] until measured). FSIM cannot say whether a physical
structure achieving a given spec EXISTS -- that is the later MEEP (mode
volume/kappa) or COMSOL (thermal) job. A finite kappa_max is a ceiling to
design under, not a promise that a cavity with that kappa can be built at the
required V_tilde. Every output of this module carries tag "[A]".

F-series math implemented (do not re-derive elsewhere -- reuse spectral.py /
loading.py / integrator.py / thermal.py closed forms):

  1. eps_budget(mu, t):        largest eps with F1b g2(mu, eps) <= t
  2. rho_required(...):        rho needed to reach t given an achieved eps_op
  3. G_required(...):          cavity collection gain implementing that rho
  4. b_e_budget(...):          background budget at G=1 (no cavity)
  5. kappa_ceiling(...):       largest tracked-cavity kappa keeping cavity-only
                                eps under the eps budget (F6 i, on resonance)
  6. v_tilde_required(...):    mode volume the cavity must deliver for a given
                                effective Purcell factor (F6 iii inverted)
  7. w_floor_of() / kappa_min_of(): closed-form brightness-floor window/kappa
     (t_X >= t_x_floor) the slit/cavity routes must clear
  8. density_limit(...):       F5 aperture/growth density ceiling (continuous)
  9. mesa_min(...):            smallest mesa clearing a junction-heating budget

spec_sheet() assembles (1)-(9) into THREE INDEPENDENT design routes (baseline
w=Gamma window with G=1; narrow slit with G=1; tracked cavity with gain G) at
one (T_op, target, delta_xx) design point, each with its own closing
condition, plus the thermal check. The overall verdict is closes_baseline OR
closes_slit OR closes_cavity -- NOT an AND over unrelated routes (that
conflation was a review-flagged honesty defect: a design can print
"INFEASIBLE" on the AND while one of its own routes already closes). See
spec_sheet()'s docstring for the per-route accounting.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .loading import f1b_g2
from .spectral import epsilon, gamma_of_T
from .thermal import Layer, Stack, t_junction

C_PUR = 3.0 / (4.0 * np.pi**2)  # F6 iii dipole-cavity coupling constant


# --------------------------------------------------------------- 1. eps budget

def eps_budget(mu, target):
    """Largest eps (spectral leakage ratio t_XX/t_X) with F1b g2(mu, eps) <=
    target, at rho = 1 (perfect background). f1b_g2 is monotone non-decreasing
    in eps at fixed mu (verify_fsim.py), so this is the tight ceiling on
    spectral leakage alone -- any rho < 1 tightens it further (see
    rho_required). Returns 1.0 if even eps=1 clears the target (unconstrained
    in [0, 1]); raises ValueError for target < 0 (not a valid g2 target).
    """
    if target < 0:
        raise ValueError(f"eps_budget: target g2 must be >= 0, got {target}")
    if target == 0.0:
        return 0.0  # f1b_g2(mu, eps) > 0 for any eps > 0 (mu > 0)
    if f1b_g2(mu, 1.0) <= target:
        return 1.0
    return float(brentq(lambda e: f1b_g2(mu, e) - target, 0.0, 1.0, xtol=1e-14, rtol=1e-14))


# ------------------------------------------------------------- 2. rho budget

def rho_required(eps_op, mu, target):
    """Signal purity rho required to reach `target` given an ALREADY-ACHIEVED
    operating eps_op (e.g. from a fixed filter/cavity window, not necessarily
    the tight eps_budget). g2dot = f1b_g2(mu, eps_op) is the cascade-only g2
    at that eps; the background law g2 = 1 - rho^2(1-g2dot) then requires
        rho_req = sqrt((1-target) / (1-g2dot)).
    If g2dot >= target, no rho <= 1 can reach the target (the optics alone
    already violate it) -- INFEASIBLE via rho alone; returns np.inf as the
    explicit infeasibility signal (eps_op must be reduced first, e.g. via a
    tighter kappa/window from kappa_ceiling/optimal_width).
    """
    g2dot = float(f1b_g2(mu, eps_op))
    if g2dot >= target:
        return float("inf")
    return float(np.sqrt((1.0 - target) / (1.0 - g2dot)))


# --------------------------------------------------------- 3. cavity gain req

def G_required(rho_req, S, B):
    """Collection gain G (F6 ii, acts on signal only) implementing rho_req
    given retention S and background B (both in the rho = GS/(GS+B) sense
    used throughout Module E/device.py). Inverting for G:
        G_req = rho_req * B / (S * (1 - rho_req)).
    Guard: rho_req >= 1 (including the rho_required() inf signal) has no
    finite solution -- returns np.inf.
    """
    if not np.isfinite(rho_req) or rho_req >= 1.0:
        return float("inf")
    if rho_req <= 0.0:
        return 0.0
    return float(rho_req * B / (S * (1.0 - rho_req)))


# ------------------------------------------------------------ 4. b_e budget

def b_e_budget(rho_req, S, B0):
    """Injection-background budget at G = 1 (no cavity): the largest extra
    background b_e the design can tolerate and still reach rho_req, given
    retention S and the T-dependent base background B0 (escaped-carrier
    channel, no injection). From rho_req = S/(S+B):
        B_max   = S (1 - rho_req) / rho_req
        b_e_max = B_max - B0.
    Returns -inf (hard infeasible, same convention as G_required) when
    rho_req >= 1. A finite but NEGATIVE b_e_max means the base background B0
    alone already exceeds the budget -- no injection background is left to
    spend; the caller must flag this INFEASIBLE (b_e >= 0 physically).
    """
    if not np.isfinite(rho_req) or rho_req >= 1.0:
        return float("-inf")
    B_max = S * (1.0 - rho_req) / rho_req
    return float(B_max - B0)


# ---------------------------------------------------------- 5. kappa ceiling

def _eps_cav(delta, gamma_x, gamma_xx, kappa):
    """Cavity-only eps on resonance (detuning 0 at T_op, the tracked-cavity
    convention): reuses spectral.epsilon's kappa acceptance directly -- no
    re-derivation of the Lorentzian-Lorentzian overlap here."""
    return epsilon(delta, gamma_x, gamma_xx, kappa=kappa, dx=0.0).eps


def kappa_ceiling(delta, gamma_x, gamma_xx, eps_budget):
    """Largest tracked cavity linewidth kappa (meV) keeping the cavity-only
    eps (F6 i, resonance at T_op) under eps_budget. eps_cav(kappa) is monotone
    increasing in kappa (narrow cavity ~ narrow filter, wide cavity -> no
    filtering, eps -> 1), so this is a ceiling: kappa <= kappa_max satisfies
    the budget, kappa > kappa_max does not.

    Bisects kappa in [1e-4, 50] meV. Returns 50.0 if even that wide a cavity
    clears the budget (unconstrained within the searched range -- NOT a claim
    that kappa can be arbitrarily large; only that this scan didn't bind).
    Returns np.nan if even kappa=1e-4 meV (near-zero linewidth) violates the
    budget -- the spectral leakage floor (X-XX splitting vs Gamma) is too
    tight for ANY cavity to fix; the eps budget itself must be relaxed
    (bigger target, smaller mu, or -- out of this module's scope -- a
    platform with a larger Delta_XX/Gamma ratio).
    """
    if _eps_cav(delta, gamma_x, gamma_xx, 50.0) <= eps_budget:
        return 50.0
    if _eps_cav(delta, gamma_x, gamma_xx, 1e-4) > eps_budget:
        return float("nan")
    return float(brentq(lambda k: _eps_cav(delta, gamma_x, gamma_xx, k) - eps_budget,
                        1e-4, 50.0, xtol=1e-10, rtol=1e-12))


# ------------------------------------------------------- 6. mode volume req

def v_tilde_required(F_eff_target, kappa, gamma, E_eV):
    """Mode volume V_tilde (units of (lambda/n)^3) the cavity must deliver for
    a target spectral-overlap-limited effective Purcell factor F_eff (F6 iii),
    given kappa and the emitter linewidth gamma (meV) it must beat, and the
    emitter energy E_eV.

    Forward relation (Purcell formula in these units, C_pur = 3/(4 pi^2)):
        F_P    = C_pur * (E_meV / kappa) / V_tilde
        F_eff  = F_P * kappa / (kappa + gamma)
               = C_pur * E_meV / ((kappa + gamma) * V_tilde)
    Inverted for V_tilde:
        V_tilde_req = C_pur * E_meV / ((kappa + gamma) * F_eff_target).

    NOTE: this depends on kappa and gamma only through their sum (kappa +
    gamma) -- the kappa dependence that matters here is the spectral-overlap
    penalty, not kappa on its own. Whether a real cavity can be built at the
    required V_tilde (diffraction limit, fabrication) is the EM solver's
    question (MEEP/COMSOL), not this module's.
    """
    E_meV = E_eV * 1000.0
    return float(C_PUR * E_meV / ((kappa + gamma) * F_eff_target))


# --------------------------------------------------- 7. brightness-floor routes

def w_floor_of(t_x_floor, gamma_x):
    """Closed-form top-hat width w with on-center X transmission t_X(w) =
    t_x_floor (the slit route's brightness floor -- [A] operating convention:
    the slit must transmit at least t_x_floor of the X line, else eps -> 0 at
    zero brightness always trivially "wins" and the route is ill-posed).

    Inverts spectral.tophat_transmission at dx=0: t_X = (2/pi) atan(w /
    (2*(gamma_x/2))) = (2/pi) atan(w/gamma_x), so
        w_floor = gamma_x * tan(pi * t_x_floor / 2).
    """
    return float(gamma_x * np.tan(np.pi * t_x_floor / 2.0))


def kappa_min_of(t_x_floor, gamma_x):
    """Closed-form tracked-cavity kappa with on-resonance X transmission
    t_X_cav(kappa) = t_x_floor (the cavity route's brightness floor -- same
    [A] convention as w_floor_of, on the cavity side).

    Inverts spectral.cavity_transmission's on-resonance value kappa/(gamma_x +
    kappa) = t_x_floor for kappa:
        kappa_min = t_x_floor * gamma_x / (1 - t_x_floor).
    """
    return float(t_x_floor * gamma_x / (1.0 - t_x_floor))


# --------------------------------------------------------- 8. density limit

def density_limit(penalty_budget, aperture_um, w, sigma, r):
    """Largest QD areal density (cm^-2) keeping the F5 aperture penalty (extra
    in-window competitor dots) under penalty_budget, for a given aperture
    diameter (um), spectral window w (meV, full width), inhomogeneous width
    sigma (meV), and relative competitor brightness r.

    CONTINUOUS inversion (matches scripts/run_phase3.run_map's forward formula
    EXACTLY -- this used to be a silent ~5x disagreement with an integer/
    rounding model here; fixed): the F5 penalty of one signal dot plus an
    EXPECTED (not rounded) number of competitors N_w is
        g2_pen(N_w) = 1 - (1 + N_w r^2) / (1 + N_w r)^2
    (the aperture_g2 lemma evaluated at the continuous expected-value limit,
    same identity run_phase3.py plots as g2pen from Nw). g2_pen is monotone
    increasing in N_w, so bisect N_w in [0, 1e6] for g2_pen(N_w) = penalty_
    budget (N_w = 1e6 if even that many competitors clear the budget --
    unconstrained within the searched range), then invert
    N_w = n_qd_cm2 * (pi (aperture_um/2)^2) * 1e-8 * w / sigma for the density:
        density = N_w * sigma / (area_um2 * 1e-8 * w).

    NOTE: this is the EXPECTED-VALUE model, not the previous rounded-integer-
    competitor-count model -- it is neither conservative nor permissive
    relative to a single realization, it is the correct continuous limit the
    rest of fsim_core (run_phase3.py) already uses. penalty_budget <= 0 returns
    density 0 (no competitors tolerated).
    """
    if penalty_budget <= 0.0:
        return 0.0

    def g2_pen(n_w):
        return 1.0 - (1.0 + n_w * r**2) / (1.0 + n_w * r) ** 2

    n_w_hi = 1e6
    if g2_pen(n_w_hi) <= penalty_budget:
        n_w = n_w_hi  # unconstrained within the searched range
    else:
        n_w = brentq(lambda n: g2_pen(n) - penalty_budget, 0.0, n_w_hi, xtol=1e-10, rtol=1e-12)

    area_um2 = np.pi * (aperture_um / 2.0) ** 2
    return float(n_w * sigma / (area_um2 * 1e-8 * w))


# ------------------------------------------------------------- 9. mesa_min

def mesa_min(P, stack_worst: Stack, T_hs, dT_max=3.0):
    """Smallest mesa diameter (um) keeping the junction-heating budget dT_J =
    T_j - T_hs <= dT_max at dissipated power P (W), evaluated at the WORST-
    case (lowest-k) epi stack -- the design must clear the budget under the
    unfavourable end of the [A] conductivity range, not just its midpoint.

    dT_J(d) is monotone decreasing in mesa diameter (more area, lower R_th):
    bisects d in [0.1, 5] um. Returns np.nan if even d=5 um violates dT_max
    (no mesa in the searched range clears the budget -- thermal runaway or a
    genuinely too-tight budget); returns 0.1 if even the smallest searched
    mesa already clears it (unconstrained within the searched range).
    """
    def dT(d_um):
        a = 0.5 * d_um * 1e-6
        Tj = t_junction(P, a, stack_worst, T_hs)
        return float("inf") if not np.isfinite(Tj) else Tj - T_hs

    if dT(5.0) > dT_max:
        return float("nan")
    if dT(0.1) <= dT_max:
        return 0.1
    return float(brentq(lambda d: dT(d) - dT_max, 0.1, 5.0, xtol=1e-6))


# ------------------------------------------------------------------ assembler

def spec_sheet(T_op, target, delta, mu=0.5, b_e_range=(1e-3, 0.3), gamma_scale=1.0,
               E_X0=1.88, F_eff_target=3.0, aperture_um=1.0, penalty_budget=0.05,
               sigma_inh=40.0, comp_brightness=0.3, V=1.9, I_uA=10.0, duty=1.0,
               epi_thickness_um=1.5, epi_k_worst=8.0, epi_alpha=0.5,
               sub_k300=55.0, sub_alpha=1.25, dT_max=3.0, t_x_floor=0.3) -> dict:
    """Assemble one staged-device requirement (design-target) sheet at a
    single (T_op, target g2(0), Delta_XX) design point, from the validated
    F-series chain. THIS IS SPEC MODE: every value returned is a NECESSARY
    CONDITION the physical device must meet to reach `target` at T_op under
    the class-proxy Gamma(T)/retention inputs -- not a claim such a device can
    be built (kappa/V_tilde achievability is MEEP's question; mesa/R_th
    achievability is COMSOL's). Tag "[A]" throughout.

    THREE INDEPENDENT ROUTES (fixes the route-conflation honesty defect: the
    old single AND-of-everything "feasible" flag mixed the un-engineered
    baseline, the no-cavity b_e budget, and the engineered kappa/slit routes,
    so a design could print INFEASIBLE while one of its own routes already
    closed). Each route is evaluated on its own terms; the overall verdict is
    the OR of the three:

      1. BASELINE (no filter engineering, G=1): eps_op at the codebase's
         standard w = Gamma(T) operating window (device.py FilterBlock.auto_w
         / run_phase3.run_staged convention). g2dot_base = f1b_g2(mu, eps_op).
         Closes iff g2dot_base < target AND b_e_budget_base >= 0 (rho_required
         + b_e_budget at that eps_op).
      2. SLIT (narrow top-hat window, no cavity, G=1): the window is narrowed
         to the brightness floor w_floor = w_floor_of(t_x_floor, Gamma) --
         the widest engineering freedom a bare slit has is trivially "make it
         infinitely narrow", which drives eps -> 0 but also brightness -> 0;
         t_x_floor pins the slit to transmit >= t_x_floor of the X line so the
         route is well-posed. eps_slit at w_floor; g2dot_slit = f1b_g2(mu,
         eps_slit). Closes iff g2dot_slit < target AND b_e_budget_slit >= 0.
      3. CAVITY (tracked cavity, gain G available): same brightness-floor
         logic on the cavity side, kappa_min = kappa_min_of(t_x_floor, Gamma).
         eps_cav_min at kappa_min; g2dot_cav = f1b_g2(mu, eps_cav_min). Closes
         iff g2dot_cav < target (no G cap here -- gain fights the background,
         not the spectral leakage; the REQUIRED G/kappa/V~ are reported, their
         achievability is the EM solver's question). kappa_max is the largest
         kappa keeping f1b_g2(mu, eps_cav(kappa)) < target (reuses
         kappa_ceiling(delta, Gamma, Gamma_XX, eps_budget(mu, target)) -- by
         construction that root is >= kappa_min whenever the route closes, so
         no new bisection is needed); G_required at worst/best b_e and
         V_tilde_required are evaluated at kappa_min (the route's design
         point: brightest cavity that still meets the floor).

    Pipeline (reuses class_proxy_params() -- the same arsenide-class Gamma(T)/
    retention proxy every other fsim_core module uses):
      * Gamma(T_op), Gamma_XX(T_op) from the class proxy, scaled by
        `gamma_scale`.
      * eps_budget(mu, target): the eps ceiling at rho=1, used only to locate
        kappa_max for the cavity route.
      * density_limit at (penalty_budget, aperture_um) -- independent of the
        three g2 routes (growth/aperture design rule).
      * mesa_min at duty*I_uA*V, worst-case epi-k -- independent thermal check.

    Returns a flat dict. "feasible" = closes_baseline OR closes_slit OR
    closes_cavity; "verdict" is a human-readable string naming exactly the
    routes that close (with their headline requirement) or, if none do,
    "INFEASIBLE: ...". "tag" is "[A]" throughout.
    """
    from .device import class_proxy_params  # local import: avoids a spec<->device cycle

    proxy = class_proxy_params()
    gam = gamma_scale * float(gamma_of_T(T_op, proxy["gamma0"], proxy["a_ac"],
                                          proxy["b_lo"], proxy["E_lo"]))
    gam_xx = proxy["r_xx"] * gam

    from .integrator import retention
    S = float(retention(T_op, proxy["a_esc"], proxy["E_a"], proxy["b_p"], proxy["E_b"]))
    B0 = proxy["b0"] + proxy["beta"] * (1.0 - S)

    b_e_lo, b_e_hi = min(b_e_range), max(b_e_range)
    eps_bud = eps_budget(mu, target)  # used only to locate the cavity route's kappa_max

    # ------------------------------------------------------------ route 1: baseline
    eps_op = epsilon(delta, gam, gam_xx, w=gam, dx=0.0).eps  # Gamma-convention window
    g2dot_base = float(f1b_g2(mu, eps_op))
    rho_req = rho_required(eps_op, mu, target)
    rho_feasible = np.isfinite(rho_req) and rho_req < 1.0
    b_e_bud_base = b_e_budget(rho_req, S, B0)
    b_e_feasible = np.isfinite(b_e_bud_base) and b_e_bud_base >= 0.0
    closes_baseline = bool(g2dot_base < target and b_e_feasible)

    # --------------------------------------------------------------- route 2: slit
    w_floor = w_floor_of(t_x_floor, gam)
    eps_slit = epsilon(delta, gam, gam_xx, w=w_floor, dx=0.0).eps
    g2dot_slit = float(f1b_g2(mu, eps_slit))
    rho_req_slit = rho_required(eps_slit, mu, target)
    b_e_bud_slit = b_e_budget(rho_req_slit, S, B0)
    slit_feasible = np.isfinite(b_e_bud_slit) and b_e_bud_slit >= 0.0
    closes_slit = bool(g2dot_slit < target and slit_feasible)

    # ------------------------------------------------------------- route 3: cavity
    kappa_min = kappa_min_of(t_x_floor, gam)
    eps_cav_min = _eps_cav(delta, gam, gam_xx, kappa_min)
    g2dot_cav = float(f1b_g2(mu, eps_cav_min))
    closes_cavity = bool(g2dot_cav < target)
    if closes_cavity:
        # eps_bud makes f1b_g2(mu, eps_cav(kappa))=target at the root by
        # construction (eps_budget's definition); since g2dot_cav already
        # clears target at kappa_min and eps_cav is monotone increasing in
        # kappa, that root sits at kappa >= kappa_min -- kappa_ceiling can be
        # reused directly, no separate bisection needed.
        kappa_max = kappa_ceiling(delta, gam, gam_xx, eps_bud)
    else:
        kappa_max = float("nan")
    kappa_feasible = np.isfinite(kappa_max)

    rho_req_cav = rho_required(eps_cav_min, mu, target)
    G_req_best = G_required(rho_req_cav, S, B0 + b_e_lo)   # least background -> easiest
    G_req_worst = G_required(rho_req_cav, S, B0 + b_e_hi)  # most background -> hardest
    G_feasible = np.isfinite(G_req_worst)
    v_tilde_req = v_tilde_required(F_eff_target, kappa_min, gam, E_X0)

    # -------------------------------------------------------- independent checks
    dens_lim = density_limit(penalty_budget, aperture_um, gam, sigma_inh, comp_brightness)

    P = duty * I_uA * 1e-6 * V
    stack_worst = Stack(layers=[Layer(epi_thickness_um * 1e-6, epi_k_worst, epi_alpha,
                                       spread=False)], k_sub300=sub_k300, sub_alpha=sub_alpha)
    Tj_op = t_junction(P, 0.5e-6, stack_worst, T_op)  # 1 um mesa reference point
    dT_J_op = Tj_op - T_op if np.isfinite(Tj_op) else float("inf")
    mesa_um = mesa_min(P, stack_worst, T_op, dT_max=dT_max)
    mesa_feasible = np.isfinite(mesa_um)

    # ------------------------------------------------------------------- verdict
    feasible = bool(closes_baseline or closes_slit or closes_cavity)
    parts = []
    if closes_baseline:
        parts.append("baseline")
    if closes_slit:
        parts.append(f"slit (b_e<={b_e_bud_slit:.3g})")
    if closes_cavity:
        g_lo, g_hi = sorted((G_req_best, G_req_worst))
        parts.append(f"cavity (G>={g_lo:.3g}..{g_hi:.3g}, "
                     f"kappa {kappa_min:.3g}..{kappa_max:.3g} meV, V~<={v_tilde_req:.3g})")
    verdict = ("closes via: " + "; ".join(parts) if parts else
              f"INFEASIBLE: eps floor exceeds budget at any brightness >= {t_x_floor:.2g}")

    return {
        "T_op": T_op, "target_g2": target, "delta_xx": delta, "mu": mu,
        "gamma_op": gam, "gamma_xx_op": gam_xx, "S_op": S, "B0_op": B0,
        "eps_budget": eps_bud, "t_x_floor": t_x_floor,
        # -- route 1: baseline (w = Gamma window, G = 1)
        "eps_op_gamma_window": eps_op, "g2dot_base": g2dot_base,
        "rho_required": rho_req, "rho_feasible": bool(rho_feasible),
        "b_e_budget_base": b_e_bud_base, "b_e_feasible": bool(b_e_feasible),
        "closes_baseline": closes_baseline,
        # -- route 2: slit (narrow top-hat at the brightness floor, G = 1)
        "w_floor": w_floor, "eps_slit": eps_slit, "g2dot_slit": g2dot_slit,
        "b_e_budget_slit": b_e_bud_slit, "slit_feasible": bool(slit_feasible),
        "closes_slit": closes_slit,
        # -- route 3: tracked cavity at the brightness floor, gain G available
        "kappa_min": kappa_min, "eps_cav_min": eps_cav_min, "g2dot_cav": g2dot_cav,
        "kappa_max": kappa_max, "kappa_feasible": bool(kappa_feasible),
        "G_required_best": G_req_best, "G_required_worst": G_req_worst,
        "G_feasible": bool(G_feasible),
        "v_tilde_required": v_tilde_req, "F_eff_target": F_eff_target,
        "closes_cavity": closes_cavity,
        # -- independent checks
        "density_limit_cm2": dens_lim, "penalty_budget": penalty_budget,
        "aperture_um": aperture_um,
        "T_j_op": Tj_op, "dT_J_op": dT_J_op,
        "mesa_min_um": mesa_um, "mesa_feasible": bool(mesa_feasible), "dT_max": dT_max,
        # -- overall verdict
        "feasible": feasible, "verdict": verdict,
        "tag": "[A]",
    }
