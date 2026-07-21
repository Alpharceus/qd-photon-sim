"""Module D -- Loading & drive statistics (cap-2 Poisson baseline, F1b,
injection background channels).

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
