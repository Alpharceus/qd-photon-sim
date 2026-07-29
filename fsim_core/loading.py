"""Module D -- Loading & drive statistics (cap-2 Poisson baseline, F1b,
injection background channels, and the F-series v1.1 electrical-penalty
channels F8/F8b/F9).

F1b (finite-mu operating point). Pulse loads n ~ Poisson(mu) excitons, capped
at 2 (cap-2 loading assumption: n>=2 relaxes to the XX ground state). Emission
through the X-centered filter: n=1 -> X photon (transmission t_X); n=2 ->
cascade XX->X (t_XX then t_X, independent). With the experimental peak-area
convention g2(0) = area(tau=0)/area(adjacent),

    <m(m-1)> = 2 P2 t_X t_XX,   <m>^2 = t_X^2 [P1 + P2 (1+eps)]^2
    g2(mu)   = 2 P2 eps / [P1 + P2 (1+eps)]^2,   eps = t_XX/t_X

Limits: mu->0 recovers F1 (g2 = eps); mu->inf saturates at 2 eps/(1+eps)^2.
The drive factor f = g2(mu)/eps in [1, 2) is the entire finite-mu penalty.

Background law exactness: for ANY dot statistics plus an independent Poissonian
background of per-pulse mean b, with rho = s/(s+b) (s = mean signal photons),
    g2_total = 1 - rho^2 (1 - g2_dot)
holds exactly -- this is the F-series background law; verified by MC.

Injection background b_e(I,T) [A -> measured]: parameterized channel set
(WL EL, barrier EL, leakage light), each  A (I/I_ref)^m exp(-E_act/kT).

F-series v1.1 additions (electrical-drive penalty channels):
  * F8  -- loading-statistics factorization: cap-2 loading generalized from a
    pure Poisson pump (F_p=1) to an arbitrary pump Fano factor F_p, moment-
    matched to (mu, F_p) rather than derived from a specific counting
    distribution. See f8_g2 / f8_g2_load below.
  * F8b -- thinning lemma: binomial dot-capture thinning of a sub-Poissonian
    pump stream washes its statistics back toward Poisson as the capture
    efficiency eta -> 0. See f8b_thin_fano / mc_thin_fano.
  * F9  -- junction granularity: the single-electron/SET scale (Coulomb
    blockade energy vs kT) set by the depletion capacitance, and a qualitative
    impedance-divider parameterization of achievable pump Fano factor. See
    granularity_N / c_dep_for_N / island_radius_nm / fano_pump.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .spectral import KB


# --------------------------------------------------------------------- cap-2 loading

def loading_probs(mu):
    """(P0, P1, P2) for cap-2 Poisson loading."""
    mu = np.asarray(mu, dtype=float)
    P0 = np.exp(-mu)
    P1 = mu * np.exp(-mu)
    return P0, P1, 1.0 - P0 - P1


def f1b_g2(mu, eps):
    """F1b: pulsed g2(0) of the filtered cascade at mean loading mu."""
    _, P1, P2 = loading_probs(mu)
    denom = P1 + P2 * (1.0 + eps)
    return np.where(denom > 0, 2.0 * P2 * eps / denom**2, 0.0)


def drive_factor(mu, eps):
    """g2(mu)/eps in [~1, 2): the finite-mu penalty on the F1 identity."""
    return f1b_g2(mu, eps) / eps if eps > 0 else 1.0


def brightness_per_pulse(mu, t_x):
    """Mean X photons per pulse through the filter (leakage excluded)."""
    _, P1, P2 = loading_probs(mu)
    return (P1 + P2) * t_x


# ------------------------------------------------------------- injection background

@dataclass
class BackgroundChannel:
    """One injection background channel: b = A (I/I_ref)^m exp(-E_act/kT).

    E_act = 0 for T-independent channels (stray/leakage light); m > 1 for
    superlinear EL (barrier/WL population under strong injection). All [A]
    until measured on the device.
    """

    name: str
    A: float
    m: float = 1.0
    E_act: float = 0.0   # meV
    I_ref: float = 1.0   # same unit as I

    def rate(self, I, T):
        act = np.exp(-self.E_act / (KB * np.asarray(T, dtype=float))) if self.E_act else 1.0
        return self.A * (np.asarray(I, dtype=float) / self.I_ref) ** self.m * act


def b_injection(channels, I, T):
    """Total injection background (per-pulse mean, in units of the signal mean
    at full transmission -- the same normalization as b0/beta in Module E)."""
    return sum(c.rate(I, T) for c in channels) if channels else 0.0


# ------------------------------------------------------- F5 aperture lemma (unequal emitters)

def aperture_g2(p_list):
    """F5 aperture lemma: N independent single-photon emitters with brightnesses
    p_i give g2 = 1 - sum(p_i^2)/(sum p_i)^2  (= 1 - 1/N for equal emitters)."""
    p = np.asarray(p_list, dtype=float)
    return 1.0 - np.sum(p**2) / np.sum(p) ** 2


def n_window_competitors(n_qd_cm2, aperture_um2, w_full_meV, sigma_inh_meV):
    """Expected in-window competitor dots: N_w ~ n_QD * A_ap * w/sigma_inh with
    w the FULL window width (foundation section 5 checkpoint: 1e10 cm^-2,
    0.1 um^2, w = 6.5, sigma = 40 -> ~10 geometric, ~1.6 spectral). Drives the
    low-density / small-aperture requirement."""
    n_geom = n_qd_cm2 * aperture_um2 * 1e-8
    return n_geom * w_full_meV / sigma_inh_meV


def mc_aperture_g2(p_list, n_pulses=1_000_000, seed=0):
    """MC second method for the aperture lemma."""
    rng = np.random.default_rng(seed)
    m = np.zeros(n_pulses, dtype=np.int64)
    for p in p_list:
        m += rng.random(n_pulses) < p
    mean = m.mean()
    g2 = (m * (m - 1)).mean() / mean**2
    se = (m * (m - 1)).std() / np.sqrt(n_pulses) / mean**2
    return float(g2), float(se)


# ------------------------------------------------------------------ MC second method

def mc_pulsed_g2(mu, t_x, t_xx, b=0.0, n_pulses=2_000_000, seed=0):
    """Monte-Carlo second method for F1b (+ Poissonian background): simulate
    photon counts per pulse, return (g2, stderr) with the all-pairs peak-area
    normalization g2 = <m(m-1)>/<m>^2."""
    rng = np.random.default_rng(seed)
    n = np.minimum(rng.poisson(mu, n_pulses), 2)
    m = np.zeros(n_pulses, dtype=np.int64)
    m += (rng.random(n_pulses) < t_x) & (n >= 1)          # X photon
    m += (rng.random(n_pulses) < t_xx) & (n >= 2)         # leaked XX photon
    if b > 0:
        m += rng.poisson(b, n_pulses)
    pairs = m * (m - 1)
    mean = m.mean()
    g2 = pairs.mean() / mean**2
    # crude stderr via pair-count fluctuation (adequate for 4-sigma checks)
    se = pairs.std() / np.sqrt(n_pulses) / mean**2
    return float(g2), float(se)


# ============================================================ F8: loading factorization
#
# F-series v1.1 [A unless noted]. Generalized cap-2 loading parameterized by
# the pump mean mu and pump Fano factor F_p (Var/mean of the per-pulse pump
# count), instead of assuming a specific counting distribution. The second
# factorial moment is moment-matched to (mu, F_p):
#
#     <k(k-1)> = mu^2 + mu (F_p - 1)          (exact for ANY distribution with
#                                               that mean and Fano factor, to
#                                               second moment order)
#     P2 = <k(k-1)> / 2
#     P1 = mu - 2 P2
#     P0 = 1 - P1 - P2
#
# Domain: P2 >= 0 requires mu >= 1 - F_p; P1 >= 0 requires F_p <= 2 - mu.
# Both are guarded; violating either raises ValueError (outside the physical
# cap-2 two-level truncation, not a silent clamp).
#
# F_p = 1 (Poisson) is the F8 family's OWN cap-2 truncation, moment-matched to
# order mu^2 with the exact Poisson-cap-2 truncation of f1b_g2 -- the two
# conventions are NOT identical at finite mu (F8's P1 = mu - mu^2 vs Poisson's
# P1 = mu e^-mu = mu - mu^2 + mu^3/2 - ...; they agree to O(mu^2), differ at
# O(mu^3)). Both are valid conventions for "the" cap-2 loading law; F8's is
# the v1.1 canonical choice because it is the only one that generalizes to
# non-Poisson (F_p != 1) pump statistics without inventing a distribution.
# The chain (integrator.g2_of_T / device.evaluate) therefore only switches
# from f1b_g2 to f8_g2 when F_p != 1.0 exactly, so every existing validated
# F_p=1 result is preserved bit-for-bit.

def _f8_probs(mu, F_p):
    """(P0, P1, P2) for the F8 generalized cap-2 loading; raises ValueError if
    (mu, F_p) sits outside the domain where P1, P2 >= 0."""
    mu_a = np.asarray(mu, dtype=float)
    Fp_a = np.asarray(F_p, dtype=float)
    tol = 1e-9
    if np.any(mu_a < 1.0 - Fp_a - tol):
        raise ValueError(
            f"f8_g2: domain violated (need mu >= 1 - F_p); mu={mu}, F_p={F_p}"
        )
    if np.any(Fp_a > 2.0 - mu_a + tol):
        raise ValueError(
            f"f8_g2: domain violated (need F_p <= 2 - mu); mu={mu}, F_p={F_p}"
        )
    m2 = mu_a**2 + mu_a * (Fp_a - 1.0)
    P2 = np.clip(0.5 * m2, 0.0, None)
    P1 = np.clip(mu_a - 2.0 * P2, 0.0, None)
    P0 = 1.0 - P1 - P2
    return P0, P1, P2


def f8_g2(mu, F_p, eps):
    """F8: pulsed g2(0) of the filtered cascade for a pump of mean mu and Fano
    factor F_p (moment-matched cap-2 loading; see module docstring). Same
    area-ratio convention as f1b_g2:

        g2 = 2 P2 eps / [P1 + P2 (1+eps)]^2

    Limits (all encoded as verify/verify4.py checks):
      (i)   small-mu factorization: g2 -> f8_g2_load(mu, F_p) * eps
      (ii)  F_p=1 agrees with f1b_g2 to O(mu^2) at small mu (differs at O(mu^3)
            -- see module docstring; NOT identical at finite mu)
      (iii) turnstile/number-state limit: F_p -> 0 at mu = 1 gives g2 -> 0
      (iv)  domain edge mu = 1 - F_p exactly: g2 = 0 (P2 = 0 there)
    """
    _, P1, P2 = _f8_probs(mu, F_p)
    denom = P1 + P2 * (1.0 + eps)
    return np.where(denom > 0, 2.0 * P2 * eps / denom**2, 0.0)


def f8_g2_load(mu, F_p):
    """F8 small-mu factorization limit: g2(mu, F_p, eps) -> g2_load * eps as
    mu -> 0 (at fixed leakage eps), with

        g2_load = 1 + (F_p - 1) / mu

    Same domain guard as f8_g2 (mu >= 1 - F_p and F_p <= 2 - mu); note the
    domain forces F_p -> 1 as mu -> 0, so this "small-mu" limit is really a
    statement about the mu -> 0 edge of the allowed (mu, F_p) band, not an
    unconstrained mu -> 0 with F_p held fixed away from 1."""
    _f8_probs(mu, F_p)  # domain guard (raises on violation); probs unused here
    mu_a = np.asarray(mu, dtype=float)
    Fp_a = np.asarray(F_p, dtype=float)
    return 1.0 + (Fp_a - 1.0) / mu_a


def mc_f8_g2(mu, F_p, t_x, t_xx, n_pulses=2_000_000, seed=0):
    """MC second method for F8 (mirrors mc_pulsed_g2's construction): build the
    explicit {P0, P1, P2} loading distribution from (mu, F_p) via _f8_probs,
    draw the per-pulse load with np.random.Generator.choice, then simulate
    the same X/cascade-XX emission + all-pairs peak-area g2 as mc_pulsed_g2.
    Returns (g2, stderr)."""
    P0, P1, P2 = (float(x) for x in _f8_probs(mu, F_p))
    rng = np.random.default_rng(seed)
    n = rng.choice(np.array([0, 1, 2]), size=n_pulses, p=[P0, P1, P2])
    m = np.zeros(n_pulses, dtype=np.int64)
    m += (rng.random(n_pulses) < t_x) & (n >= 1)          # X photon
    m += (rng.random(n_pulses) < t_xx) & (n >= 2)         # leaked XX photon
    pairs = m * (m - 1)
    mean = m.mean()
    g2 = pairs.mean() / mean**2
    se = pairs.std() / np.sqrt(n_pulses) / mean**2
    return float(g2), float(se)


# ================================================================= F8b: thinning lemma

def f8b_thin_fano(eta, F_pump):
    """F8b thinning lemma: binomial thinning of a pump stream (Fano factor
    F_pump) by dot-capture probability eta gives a thinned-stream Fano factor

        F_thinned = 1 - eta (1 - F_pump)

    eta -> 0 washes any pump statistics back toward Poisson (F_thinned -> 1);
    eta -> 1 preserves the pump's own statistics exactly (F_thinned -> F_pump).
    A sub-Poissonian pump (F_pump < 1, e.g. Zhao's room-temperature quiet
    pump [V]) stays sub-Poissonian for any eta in (0, 1]."""
    return 1.0 - eta * (1.0 - F_pump)


def mc_thin_fano(n_source, eta, n_pulses=2_000_000, seed=0):
    """F8b MC second method: binomially thin an explicit integer per-pulse
    source stream by capture probability eta and measure both the source and
    thinned Fano factors (Var/mean) directly from samples -- independent of
    the f8b_thin_fano closed form.

    n_source: an array of per-pulse source counts, or a callable
    (rng, n_pulses) -> counts array (e.g. a fixed N=1-per-window stream giving
    F=0, or a two-point mixture giving an intermediate sub/super-Poissonian
    F). Returns (F_source, F_thin, se_thin); se_thin is a crude normal-approx
    stderr on F_thin, adequate for a coarse (few-sigma) comparison."""
    rng = np.random.default_rng(seed)
    n = n_source(rng, n_pulses) if callable(n_source) else np.asarray(n_source)
    thinned = rng.binomial(n, eta)
    F_source = float(n.var(ddof=1) / n.mean())
    m, v = thinned.mean(), thinned.var(ddof=1)
    F_thin = float(v / m)
    # crude delta-method stderr on Var/mean; generous factor for 4th-moment slack
    se_thin = float(3.0 * np.sqrt(v) / m / np.sqrt(n_pulses))
    return F_source, F_thin, se_thin


# ============================================================== F9: junction granularity
#
# SI constants (F9 is the one place in this module that works in SI, not
# meV/K -- the Coulomb-blockade energy scale e^2/2C needs e in coulombs).

KB_SI = 1.380649e-23        # J/K
E_SI = 1.602176634e-19      # C
EPS0_SI = 8.8541878128e-12  # F/m


def granularity_N(C_dep_F, T_K):
    """F9: junction granularity N = k_B T C_dep / e^2 [E] -- the number of
    thermal quanta k_B T sitting in the single-electron charging energy
    e^2/C_dep. N >> 1: the junction is thermally smooth, many electrons'
    worth of charging energy is washed out by kT (collective/statistical
    regulation regime, the regime Zhao's quiet-pump demonstration operates
    in). N ~ 1: single-electron/Coulomb-blockade (SET) regime, set by the
    device's own depletion capacitance -- this is a wall set by C_dep, not by
    engineering choice. The connection between N and the achievable pump
    Fano factor F_p is QUALITATIVE ONLY (see fano_pump below); this function
    computes N, it does not derive F_p from N."""
    return KB_SI * T_K * C_dep_F / E_SI**2


def c_dep_for_N(N, T_K):
    """F9 inverse of granularity_N: the depletion capacitance giving
    granularity N at temperature T_K."""
    return N * E_SI**2 / (KB_SI * T_K)


def island_radius_nm(C_F, eps_r=13.0):
    """F9: self-capacitance-sphere estimate of an island radius (nm) with
    total capacitance C_F, in a dielectric of relative permittivity eps_r --
    R = C / (4 pi eps0 eps_r). [E]: an idealized isolated-sphere model; a real
    gated island (image charges, substrate, leads) shifts the true C(R)
    relation by an O(1) factor, so treat this as an order-of-magnitude size
    estimate, not a device dimension."""
    return C_F / (4.0 * np.pi * EPS0_SI * eps_r) * 1e9


def fano_pump(Z_drive, Z_junction):
    """F9: impedance-divider parameterization [E] of the pump Fano factor
    delivered to the junction -- NOT a derived noise theory, just a
    monotone-limits parameterization consistent with two known regimes:

        Z_drive >> Z_junction -> F_p -> 0   (quiet/high-impedance current
                                              drive; Zhao's room-temperature
                                              demonstration [V] sits in this
                                              direction)
        Z_drive -> 0          -> F_p -> 1   (shot-noise-limited voltage
                                              drive, Poisson pump)

    F_p = Z_junction / (Z_junction + Z_drive)."""
    return Z_junction / (Z_junction + Z_drive)
