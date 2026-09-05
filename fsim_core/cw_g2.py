"""Module D-CW -- continuous-wave (DC-driven) second-order correlation of the
spectrally filtered exciton-biexciton cascade, Poissonian background, and the
detector-IRF convolution that turns the intrinsic g2(tau) into the RAW
measured g2(0) of a CW HBT experiment (Reischle 2008-style dc EL).

Why this module exists. The rest of the chain (loading.f1b_g2, f8_g2,
drive_mech.reexc_g2) computes the PULSED peak-area g2(0) = area(tau=0 peak) /
area(adjacent peak). Real electrically driven devices are frequently run dc,
and a dc HBT measurement reports the depth of a dip at tau = 0 of width
~1/(pump + decay rate), smeared by the detector timing resolution (IRF).
Those are different observables; this module provides the CW one.

Model (classical rate equations -- coherences neglected; adequate for
incoherent electrical/above-band pumping where dephasing >> Rabi coupling):

    states |0>, |X>, |XX>              (cap-2 ladder, consistent with loading.py)
    0  -> X   pump           r
    X  -> XX  pump           r_2 = pump_ratio * r           (default ratio 1)
    X  -> 0   radiative      gamma_X     (emits an X photon)
    X  -> 0   non-radiative  k_X         (no photon; thermal escape)
    XX -> X   radiative      gamma_XX    (emits an XX photon)
    XX -> X   non-radiative  k_XX        (no photon)

The escape rates come from the retention model of integrator.retention:
S = gamma/(gamma + k), k = gamma (a_esc e^{-E_a/kT} + b_p e^{-E_b/kT}); see
escape_rates_from_retention().

Two-time correlations use the quantum regression theorem, which for a
classical Markov chain is just "propagate the post-jump state with the
generator": after a photon of type i at t = 0 the dot is in the post-jump
state (X photon -> |0>, XX photon -> |X>), and

    G_ij(tau) = I_i * gamma_j * [exp(M tau) e_post(i)]_(state emitting j),

with I_X = gamma_X P_X, I_XX = gamma_XX P_XX the steady-state photon rates.
Behind an X-centered filter with transmissions t_X, t_XX (spectral.py;
eps = t_XX/t_X) the detected stream has

    G(tau)   = sum_ij t_i t_j G_ij(tau),   I_det = t_X I_X + t_XX I_XX,
    g2_dot   = G(|tau|) / I_det^2.

Background: g2_meas(tau) = 1 - rho^2 (1 - g2_dot(tau)), rho = S/(S+B) -- the
F-series background law applied pointwise (integrator.g2_from). It is exact
for a Poissonian background statistically independent of the dot, at every
tau, because the cross terms of independent stationary streams factorize.

IRF: raw g2(0) = (K * g2_meas)(0) for a normalized kernel K -- Gaussian
(parameterized by its FWHM) or two-sided exponential (Reischle 2008 fitted
"1 - A exp(-|tau|/tau_c) convolved with C exp(-|tau|/0.5 ns)", i.e. an
exponential IRF of time constant 0.5 ns). Closed forms for the single-
exponential dip are in dip_convolved / dip_convolved_exp, inverted by
deconvolve_dip (what Reischle did to get 0.15/0.25 from raw 0.41/0.43).

THE CW/PULSED DIFFERENCE (exact, derived in cw_vs_pulsed_note):

    g2_dot(0) = eps * p * (S_XX/S_X) * (1 + a + a b) / (1 + eps p (r/gamma_X) S_XX)^2

    a = r/(gamma_X + k_X),  b = p r/(gamma_XX + k_XX),  p = pump_ratio,
    S_X = gamma_X/(gamma_X + k_X),  S_XX = gamma_XX/(gamma_XX + k_XX).

    r -> 0 : g2_dot(0) -> eps p S_XX/S_X   (F1 identity recovered; = eps for
             p = 1 and no escape).  NOTE: the a-priori expectation that the
             CW leak penalty is amplified by gamma_XX/r >> 1 at low pump does
             NOT hold -- both G(0) and I_det^2 scale as r^2 (see the note).
    r >> gamma : g2_dot(0) -> gamma_X/(eps gamma_XX)  (unbounded cascade
             BUNCHING of the filtered stream; the pulsed drive factor is
             bounded by 2/(1+eps)^2).
    eps = 0 : g2_dot(0) = 0 exactly at every r (perfect filter or trion).

Units: rates in 1/ns, times in ns, IRF widths in ps (argument names say so),
energies meV, temperature K. Every closed form here is checked in
verify/verify_cw_g2.py.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
from scipy.special import erfc

from .integrator import g2_from
from .spectral import KB

_GAUSS_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # sigma = FWHM * this
_EXP_FWHM_TO_TAU = 1.0 / (2.0 * np.log(2.0))                       # tau_i = FWHM * this

IRF_SHAPES = ("gaussian", "exponential")


# ------------------------------------------------------------------- escape rates

def escape_rates_from_retention(gamma_X_ns, a_esc, E_a, b_p, E_b, T,
                                gamma_XX_ns=None):
    """Non-radiative escape rates (1/ns) consistent with integrator.retention:

        k = gamma * (a_esc e^{-E_a/kT} + b_p e^{-E_b/kT})   =>   S = gamma/(gamma + k)

    Returns (k_X, k_XX). The same Arrhenius factor is applied to the XX
    level; gamma_XX_ns defaults to 2 gamma_X_ns (the standard cascade ratio
    tau_XX = tau_X/2 used by drive_mech.reexc_g2). Energies in meV, T in K."""
    kT = KB * float(T)
    esc = a_esc * np.exp(-E_a / kT) + b_p * np.exp(-E_b / kT)
    if gamma_XX_ns is None:
        gamma_XX_ns = 2.0 * gamma_X_ns
    return float(gamma_X_ns * esc), float(gamma_XX_ns * esc)


# ------------------------------------------------------------------ rate generator

def generator(r_ns, gamma_X_ns, gamma_XX_ns, k_X=0.0, k_XX=0.0, pump_ratio=1.0):
    """Generator M (3x3, columns = from-state, rows = to-state) of dP/dt = M P
    for the ordered basis (|0>, |X>, |XX>)."""
    r2 = pump_ratio * r_ns
    GX = gamma_X_ns + k_X
    GXX = gamma_XX_ns + k_XX
    return np.array([
        [-r_ns,  GX,          0.0],
        [ r_ns, -(r2 + GX),   GXX],
        [ 0.0,   r2,         -GXX],
    ], dtype=float)


def steady_state(r_ns, gamma_X_ns, gamma_XX_ns, k_X=0.0, k_XX=0.0, pump_ratio=1.0):
    """Stationary occupations (P_0, P_X, P_XX) of the chain. The chain is a
    birth-death process, so detailed balance holds: P_X = a P_0, P_XX = b P_X
    with a = r/(gamma_X + k_X), b = r_2/(gamma_XX + k_XX). (verify checks
    M P_ss = 0 against the numeric null vector.)"""
    a = r_ns / (gamma_X_ns + k_X)
    b = pump_ratio * r_ns / (gamma_XX_ns + k_XX)
    P0 = 1.0 / (1.0 + a + a * b)
    return np.array([P0, a * P0, a * b * P0])


def photon_rates(r_ns, gamma_X_ns, gamma_XX_ns, k_X=0.0, k_XX=0.0, pump_ratio=1.0):
    """(I_X, I_XX) steady-state photon emission rates (1/ns) before the filter."""
    P = steady_state(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    return gamma_X_ns * P[1], gamma_XX_ns * P[2]


# ------------------------------------------------------------ propagation helpers

def propagate(M, v0, taus, method="eig"):
    """P(tau) = exp(M tau) v0 for every tau in `taus` (array, shape (3, N)).

    method="eig": one eigendecomposition M = V diag(lam) V^-1 and a vectorized
    exponential (the fast path). method="expm": scipy.linalg.expm per tau
    (the slow reference used by verify to check the eigen path)."""
    taus = np.atleast_1d(np.asarray(taus, dtype=float))
    v0 = np.asarray(v0, dtype=float)
    if method == "expm":
        return np.column_stack([expm(M * t) @ v0 for t in taus])
    lam, V = np.linalg.eig(M)
    c = np.linalg.solve(V, v0)
    # P(tau) = V @ (exp(lam tau) * c); lam may come out complex-typed; the
    # generator of a birth-death chain has real spectrum, so drop the 0j.
    E = np.exp(np.outer(lam, taus)) * c[:, None]
    out = np.real(V @ E)
    # exp(M 0) = 1 exactly: pin tau = 0 to v0 so that identities like
    # g2_dot(0) = 0 at eps = 0 hold bit-exactly instead of to roundoff.
    out[:, taus == 0.0] = v0[:, None]
    return out


# ---------------------------------------------------------------- correlations

def g2_cw(tau_ns, r_ns, gamma_X_ns, gamma_XX_ns, t_X=1.0, t_XX=0.0,
          k_X=0.0, k_XX=0.0, pump_ratio=1.0, method="eig", return_parts=False):
    """Intrinsic CW g2_dot(tau) of the filtered stream (background-free).

    tau_ns may be any array (negative values use |tau|: the total
    autocorrelation of a stationary stream is even). Returns g2_dot(tau); with
    return_parts=True returns (g2_dot, dict(G=..., I_det=..., G_ij=...)).
    """
    tau = np.abs(np.atleast_1d(np.asarray(tau_ns, dtype=float)))
    M = generator(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    P = steady_state(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    I_X, I_XX = gamma_X_ns * P[1], gamma_XX_ns * P[2]
    e0 = np.array([1.0, 0.0, 0.0])   # post-jump state after an X photon
    eX = np.array([0.0, 1.0, 0.0])   # post-jump state after an XX photon
    P_after_X = propagate(M, e0, tau, method)    # (3, N)
    P_after_XX = propagate(M, eX, tau, method)
    G_XX = I_X * gamma_X_ns * P_after_X[1]        # X then X
    G_XB = I_X * gamma_XX_ns * P_after_X[2]       # X then XX
    G_BX = I_XX * gamma_X_ns * P_after_XX[1]      # XX then X  (the cascade partner)
    G_BB = I_XX * gamma_XX_ns * P_after_XX[2]     # XX then XX
    G = (t_X * t_X * G_XX + t_X * t_XX * G_XB
         + t_XX * t_X * G_BX + t_XX * t_XX * G_BB)
    I_det = t_X * I_X + t_XX * I_XX
    g2 = G / I_det**2 if I_det > 0 else np.zeros_like(G)
    if return_parts:
        return g2, {"G": G, "I_det": I_det, "I_X": I_X, "I_XX": I_XX,
                    "G_ij": {"X,X": G_XX, "X,XX": G_XB, "XX,X": G_BX, "XX,XX": G_BB}}
    return g2


def g2_cw_zero(r_ns, gamma_X_ns, gamma_XX_ns, eps, k_X=0.0, k_XX=0.0,
               pump_ratio=1.0):
    """Closed form for g2_dot(0) (derived in cw_vs_pulsed_note; depends on the
    filter only through eps = t_XX/t_X):

        g2_dot(0) = eps p (S_XX/S_X) (1 + a + a b) / (1 + eps p (r/gamma_X) S_XX)^2
    """
    p = pump_ratio
    S_X = gamma_X_ns / (gamma_X_ns + k_X)
    S_XX = gamma_XX_ns / (gamma_XX_ns + k_XX)
    a = r_ns / (gamma_X_ns + k_X)
    b = p * r_ns / (gamma_XX_ns + k_XX)
    return (eps * p * (S_XX / S_X) * (1.0 + a + a * b)
            / (1.0 + eps * p * (r_ns / gamma_X_ns) * S_XX) ** 2)


def g2_cw_zero_low_pump(gamma_X_ns, gamma_XX_ns, eps, k_X=0.0, k_XX=0.0,
                        pump_ratio=1.0):
    """r -> 0 limit of g2_cw_zero: eps * pump_ratio * S_XX / S_X."""
    S_X = gamma_X_ns / (gamma_X_ns + k_X)
    S_XX = gamma_XX_ns / (gamma_XX_ns + k_XX)
    return eps * pump_ratio * S_XX / S_X


def g2_two_level(tau_ns, r_ns, gamma_X_ns, k_X=0.0):
    """Analytic incoherently pumped two-level result (the pump_ratio = 0 limit
    of g2_cw, filter-independent): g2(tau) = 1 - exp(-(r + gamma_X + k_X)|tau|)."""
    tau = np.abs(np.asarray(tau_ns, dtype=float))
    return 1.0 - np.exp(-(r_ns + gamma_X_ns + k_X) * tau)


# ------------------------------------------------------------------- background

def g2_with_background(g2_dot, rho):
    """Background law applied pointwise in tau: 1 - rho^2 (1 - g2_dot(tau)).
    Exact for a Poissonian background uncorrelated with the dot (the same law
    integrator.g2_from applies at tau = 0)."""
    return g2_from(np.asarray(g2_dot, dtype=float), rho)


# -------------------------------------------------------------------------- IRF

def irf_width(fwhm_ps, shape="gaussian"):
    """Shape parameter (ns) from the IRF FWHM (ps): sigma for a Gaussian,
    time constant tau_i for a two-sided exponential (FWHM = 2 tau_i ln 2)."""
    if shape not in IRF_SHAPES:
        raise ValueError(f"irf shape must be one of {IRF_SHAPES}, got {shape!r}")
    f = fwhm_ps * 1e-3
    return f * (_GAUSS_FWHM_TO_SIGMA if shape == "gaussian" else _EXP_FWHM_TO_TAU)


def irf_kernel(dt_ns, fwhm_ps, shape="gaussian"):
    """Discrete, unit-sum, symmetric IRF kernel sampled at step dt_ns."""
    wpar = irf_width(fwhm_ps, shape)
    if wpar <= 0:
        return np.array([1.0])
    reach = 8.0 * wpar if shape == "gaussian" else 25.0 * wpar
    n = int(np.ceil(reach / dt_ns))
    s = np.arange(-n, n + 1) * dt_ns
    k = np.exp(-0.5 * (s / wpar) ** 2) if shape == "gaussian" else np.exp(-np.abs(s) / wpar)
    return k / k.sum()


def convolve_irf(tau_ns, g2, fwhm_ps, shape="gaussian"):
    """Convolve g2(tau) sampled on a UNIFORM tau grid with the IRF. Beyond the
    grid g2 is continued with its edge values (the curve has already reached
    its asymptote there if tau_max is chosen sensibly), so no wrap-around."""
    tau = np.asarray(tau_ns, dtype=float)
    g = np.asarray(g2, dtype=float)
    dts = np.diff(tau)
    if tau.size < 2 or not np.allclose(dts, dts[0], rtol=1e-6, atol=0):
        raise ValueError("convolve_irf needs a uniform tau grid")
    k = irf_kernel(dts[0], fwhm_ps, shape)
    h = k.size // 2
    if h == 0:
        return g.copy()
    gp = np.concatenate([np.full(h, g[0]), g, np.full(h, g[-1])])
    return np.convolve(gp, k, mode="valid")


def raw_g2_0(tau_ns, g2, fwhm_ps, shape="gaussian"):
    """IRF-convolved value at tau = 0 (linear interpolation on the grid)."""
    conv = convolve_irf(tau_ns, g2, fwhm_ps, shape)
    return float(np.interp(0.0, np.asarray(tau_ns, dtype=float), conv))


def dip_convolved(A, tau_d, sigma):
    """Gaussian-IRF closed form for the single-exponential dip
    g2(tau) = 1 - A exp(-|tau|/tau_d):

        g2_raw(0) = 1 - A exp(sigma^2 / (2 tau_d^2)) erfc(sigma / (sqrt 2 tau_d))

    (sigma = IRF standard deviation, same time unit as tau_d)."""
    x = sigma / tau_d
    return 1.0 - A * np.exp(0.5 * x * x) * erfc(x / np.sqrt(2.0))


def dip_convolved_exp(A, tau_d, tau_irf):
    """Two-sided-exponential-IRF closed form (Reischle 2008's fit model,
    IRF C exp(-|tau|/tau_irf)):  g2_raw(0) = 1 - A tau_d / (tau_d + tau_irf)."""
    return 1.0 - A * tau_d / (tau_d + tau_irf)


def _dip_factor(tau_d_ns, fwhm_ps, shape):
    """F such that g2_raw(0) = 1 - A F for the single-exponential dip."""
    wpar = irf_width(fwhm_ps, shape)
    if shape == "gaussian":
        return 1.0 - dip_convolved(1.0, tau_d_ns, wpar)
    return tau_d_ns / (tau_d_ns + wpar)


def deconvolve_dip(g2_raw0, tau_d_ns, fwhm_ps, shape="gaussian"):
    """Invert dip_convolved for the intrinsic dip depth A given the raw
    (IRF-convolved) g2(0), the dip recovery time tau_d (ns) and the IRF FWHM
    (ps). The intrinsic g2(0) is 1 - A. This is the Reischle 2008 procedure."""
    return (1.0 - g2_raw0) / _dip_factor(tau_d_ns, fwhm_ps, shape)


def intrinsic_from_raw(g2_raw0, tau_dip_ns, irf_fwhm_ps, shape="gaussian"):
    """Intrinsic (IRF-deconvolved, background still included) g2(0) = 1 - A."""
    return 1.0 - deconvolve_dip(g2_raw0, tau_dip_ns, irf_fwhm_ps, shape)


# ------------------------------------------------------------------ dip fitting

def fit_dip(tau_ns, g2):
    """Least-squares fit of 1 - A exp(-|tau|/tau_d) to a g2(tau) curve
    (Reischle's fit model). Returns (A, tau_d). Falls back to the 1/e
    recovery point of the dip if the fit does not converge."""
    tau = np.asarray(tau_ns, dtype=float)
    g = np.asarray(g2, dtype=float)
    f = lambda t, A, td: 1.0 - A * np.exp(-np.abs(t) / td)
    A0 = float(np.clip(1.0 - g[np.argmin(np.abs(tau))], 1e-3, 2.0))
    td0 = max(0.1 * (tau.max() - tau.min()), 1e-6)
    try:
        popt, _ = curve_fit(f, tau, g, p0=[A0, td0], maxfev=20000)
        return float(popt[0]), float(abs(popt[1]))
    except (RuntimeError, ValueError):
        pos = tau >= 0
        t, d = tau[pos], 1.0 - g[pos]
        if d[0] <= 0:
            return float(d[0]), float("nan")
        below = np.where(d <= d[0] / np.e)[0]
        return float(d[0]), (float(t[below[0]]) if below.size else float("nan"))


# ----------------------------------------------------------------------- report

def cw_report(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, t_X, t_XX, rho,
              irf_fwhm_ps, tau_max_ns=10.0, pump_ratio=1.0,
              irf_shape="gaussian", n_tau=None):
    """CW g2 report for one operating point.

    Returns a dict with P_ss, I_X, I_XX, I_det, eps, g2_dot0 (exact closed
    form), g2_meas0 (background law), g2_raw0 (IRF-convolved g2_meas at tau
    = 0), tau_dip and A_dip (single-exponential dip fitted to g2_meas),
    g2_intrinsic_from_raw (what deconvolving the raw value with the fitted
    tau_dip gives back), curves (tau, g2_dot, g2_meas, g2_raw), and notes."""
    eps = t_XX / t_X if t_X > 0 else float("nan")
    wpar = irf_width(irf_fwhm_ps, irf_shape) if irf_fwhm_ps > 0 else np.inf
    rates = [x for x in (r_ns, pump_ratio * r_ns, gamma_X_ns + k_X, gamma_XX_ns + k_XX) if x > 0]
    dt = min(0.05 / max(rates), wpar / 20.0, tau_max_ns / 500.0)
    if n_tau is None:
        n_tau = 2 * int(np.ceil(tau_max_ns / dt)) + 1
    tau = np.linspace(-tau_max_ns, tau_max_ns, n_tau)
    g2_dot = g2_cw(tau, r_ns, gamma_X_ns, gamma_XX_ns, t_X, t_XX, k_X, k_XX, pump_ratio)
    g2_meas = g2_with_background(g2_dot, rho)
    g2_raw = (convolve_irf(tau, g2_meas, irf_fwhm_ps, irf_shape)
              if irf_fwhm_ps > 0 else g2_meas.copy())
    P = steady_state(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    I_X, I_XX = gamma_X_ns * P[1], gamma_XX_ns * P[2]
    g2_dot0 = g2_cw_zero(r_ns, gamma_X_ns, gamma_XX_ns, eps, k_X, k_XX, pump_ratio)
    g2_meas0 = float(g2_from(g2_dot0, rho))
    g2_raw0 = float(np.interp(0.0, tau, g2_raw))
    A_dip, tau_dip = fit_dip(tau, g2_meas)
    notes = [
        "rate-equation model (coherences neglected); cap-2 ladder |0>,|X>,|XX>",
        "g2_dot(0) closed form: eps p (S_XX/S_X)(1+a+ab)/(1+eps p (r/gamma_X) S_XX)^2",
        "background law applied pointwise (exact for uncorrelated Poissonian background)",
        f"IRF: {irf_shape}, FWHM {irf_fwhm_ps:g} ps; raw g2(0) is the IRF-convolved g2_meas at tau=0",
        "tau_dip/A_dip: least-squares single-exponential dip fitted to g2_meas(tau)",
    ]
    if g2_dot0 > 1.0:
        notes.append("g2_dot(0) > 1: cascade bunching of the leaked XX partner dominates at this pump")
    return {
        "P_ss": P, "I_X": float(I_X), "I_XX": float(I_XX),
        "I_det": float(t_X * I_X + t_XX * I_XX), "eps": float(eps),
        "g2_dot0": float(g2_dot0), "g2_meas0": g2_meas0, "g2_raw0": g2_raw0,
        "tau_dip": tau_dip, "A_dip": A_dip,
        "g2_intrinsic_from_raw": (float(intrinsic_from_raw(g2_raw0, tau_dip, irf_fwhm_ps, irf_shape))
                                  if irf_fwhm_ps > 0 and np.isfinite(tau_dip) else g2_raw0),
        "curves": {"tau": tau, "g2_dot": g2_dot, "g2_meas": g2_meas, "g2_raw": g2_raw},
        "notes": notes,
    }


# --------------------------------------------------------------------- the note

def cw_vs_pulsed_note():
    """Why the CW filtered g2(0) differs from the pulsed peak-area g2(0), with
    the exact CW closed form and its limits."""
    return """\
CW vs PULSED g2(0) of the filtered X line (cap-2 ladder, rate equations)

Pulsed (loading.f1b_g2): g2 = 2 P2 eps / [P1 + P2 (1+eps)]^2 -- the peak-AREA
ratio lumps every photon of a pulse period together, so the leak penalty is
bounded: drive factor g2/eps in [1, 2/(1+eps)^2).

CW: the only way to get two detected photons at the SAME instant is the
cascade partner: an XX photon leaking through the filter (t_XX) followed by
the X photon it triggers (t_X). At tau = 0 the post-jump state after an XX
photon is |X>, which emits X at rate gamma_X, so

    G(0)   = t_X t_XX gamma_X gamma_XX P_XX          (only the XX->X term survives:
                                                     [e^{M 0} e_0]_X = [e^{M 0} e_0]_XX = 0,
                                                     [e^{M 0} e_X]_XX = 0)
    I_det  = t_X gamma_X P_X + t_XX gamma_XX P_XX
    g2_dot(0) = G(0)/I_det^2
              = eps gamma_X gamma_XX P_XX / (gamma_X P_X + eps gamma_XX P_XX)^2 .

With detailed balance P_X = a P_0, P_XX = b P_X, a = r/(gamma_X+k_X),
b = p r/(gamma_XX+k_XX) (p = pump_ratio), S_X = gamma_X/(gamma_X+k_X),
S_XX = gamma_XX/(gamma_XX+k_XX):

    g2_dot(0) = eps p (S_XX/S_X) (1 + a + a b) / (1 + eps p (r/gamma_X) S_XX)^2   [exact]

Limits:
  r -> 0     : g2_dot(0) -> eps p S_XX/S_X  (= eps for p = 1, no escape). Both
               G(0) ~ P_XX ~ r^2 and I_det^2 ~ P_X^2 ~ r^2, so the CW leak
               penalty is NOT amplified by gamma_XX/r at low pump: the F1
               identity g2 = eps is recovered, exactly as for the pulsed
               mu -> 0 limit. (The a-priori guess of a gamma_XX/r enhancement
               was wrong and is superseded by this derivation.)
  eps = 0    : g2_dot(0) = 0 at every pump (perfect filter / trion).
  r >> gamma : g2_dot(0) -> gamma_X/(eps gamma_XX): unbounded BUNCHING of the
               filtered stream. Physically, the cascade pairs (XX leak, X)
               arrive within ~1/(gamma_X+k_X+p r) of each other and form a
               bunching peak at |tau| < tau_X on top of the antibunching dip;
               the pulsed area ratio cannot exceed 2 eps/(1+eps)^2.
  finite r   : the CW drive factor D_CW = (1+a+ab)/(1+eps p r S_XX/gamma_X)^2
               is ~ 1 + (r/gamma_X)(1 - 2 eps) at small r versus the pulsed
               ~ 1 + mu (1 - eps), and keeps growing ~ r^2 while the pulsed
               factor saturates; hence at r ~ gamma_X the CW raw g2(0) exceeds
               the pulsed peak-area value (verify_cw_g2 check e).

Measurement: the dc HBT histogram is the IRF convolution of
g2_meas(tau) = 1 - rho^2 (1 - g2_dot(tau)); a 0.5 ns IRF on a ~1 ns dip
raises 0.15 to ~0.4 (Reischle 2008). deconvolve_dip inverts that.
"""
