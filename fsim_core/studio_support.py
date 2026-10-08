"""Small closed forms the FSIM Studio shows next to the literature cards and the CW explorer.

They live here, in fsim_core, so that the Flask layer (fsim_studio) carries no physics (CONTRACT
rule 1). Each function is a one-line closed form with its derivation and provenance tag in the
docstring, and verify/verify_studio_support.py checks it against an independent closed form or the
committed out/ CSVs. Nothing here is imported by the legacy paths; verify_fsim.py and
audit_physics.py are unaffected.

Provenance (CLAUDE.md): [V] verified against the cited paper, [DR] derived, [E] estimate,
[A] assumption.
"""
from __future__ import annotations

from .integrator import retention


def tau_anchored(T, T0, tau0_ns, a_esc, E_a, b_p, E_b):
    """X decay time anchored to one measured point: tau(T) = tau0 * S(T) / S(T0).  [DR]

    Derivation. The V-a retention model has a single lifetime-like factor, the retention
    S(T) = 1 / (1 + a_esc exp(-E_a/kT) + b_p exp(-E_b/kT)) (integrator.retention). The V-a fit
    constrains only how the decay time changes with T, i.e. the ratio tau(T)/tau(T0) = S(T)/S(T0);
    the absolute scale is not predicted. Anchoring the ratio to the measured tau0 at T0 gives the
    formula above. By construction tau(T0) = tau0.  Tag: [DR] (derived from the fitted [DR]
    retention parameters, out/phase0/fit_params.json; the anchor tau0 keeps its card tag).
    """
    return tau0_ns * (retention(T, a_esc, E_a, b_p, E_b) / retention(T0, a_esc, E_a, b_p, E_b))


def inverted_background_g2(g2_b, rho):
    """Intrinsic g2 recovered from a background-diluted one: g2_s = (g2_b - (1 - rho^2)) / rho^2.  [DR]

    Derivation. The F-series background law (integrator.g2_from) is
        g2_b = 1 - rho^2 (1 - g2_s),
    with rho = S/(S + B) the signal fraction. Solving for g2_s:
        1 - g2_b = rho^2 (1 - g2_s)  ->  g2_s = 1 - (1 - g2_b)/rho^2 = (g2_b - 1 + rho^2)/rho^2.
    This is the algebraic inverse Reischle et al., Opt. Express 16, 12771 (2008), Eq. (1) uses.
    Tag: [DR].
    """
    if not rho > 0.0:
        raise ValueError("rho must be positive")
    return (g2_b - (1.0 - rho * rho)) / (rho * rho)


def g2_cw_high_pump_limit(gamma_X, gamma_XX, eps, pump_ratio):
    """r -> infinity limit of the CW dot g2(0): gamma_X / (eps gamma_XX); None if undefined.  [DR]

    Derivation. cw_g2's exact result is
        g2_dot(0) = eps p (S_XX/S_X) (1 + a + a b) / (1 + eps p (r/gamma_X) S_XX)^2,
        a = r/(gamma_X + k_X), b = p r/(gamma_XX + k_XX), S_X = gamma_X/(gamma_X + k_X),
        S_XX = gamma_XX/(gamma_XX + k_XX).
    For r -> infinity the numerator grows as a b = p r^2/((gamma_X + k_X)(gamma_XX + k_XX)) and the
    denominator as (eps p r S_XX/gamma_X)^2, so
        g2_dot(0) -> gamma_X^2 / (eps S_X S_XX (gamma_X + k_X)(gamma_XX + k_XX))
                   = gamma_X^2 / (eps gamma_X gamma_XX) = gamma_X / (eps gamma_XX),
    because S_X (gamma_X + k_X) = gamma_X and S_XX (gamma_XX + k_XX) = gamma_XX: independent of
    the escape rates and of p. Undefined (None) for eps = 0 (g2_dot(0) = 0 at every r) or p = 0
    (no XX population). Matches cw_g2.py:65-77. Tag: [DR].
    """
    if eps > 0.0 and pump_ratio > 0.0:
        return gamma_X / (eps * gamma_XX)
    return None


def implied_dip_recovery_time(g2_raw0, g2_intrinsic, tau_irf_ns):
    """Dip recovery time implied by a published raw/intrinsic pair through an exponential IRF.  [DR]

    Derivation. cw_g2.dip_convolved_exp gives g2_raw(0) = 1 - A tau_d/(tau_d + tau_irf) for the
    dip 1 - A exp(-|tau|/tau_d) seen through a two-sided exponential IRF of time constant tau_irf.
    With the published intrinsic g2(0) = 1 - A (Reischle, Opt. Express 16, 12771 (2008): raw 0.41
    -> 0.15 with the 0.5 ns IRF), solving for tau_d:
        (1 - g2_raw)(tau_d + tau_irf) = A tau_d
        tau_d = tau_irf (1 - g2_raw) / (A - (1 - g2_raw)),   A = 1 - g2_intrinsic.
    The paper does not print tau_d; this is what the two published numbers imply. Raises ValueError
    when A <= 1 - g2_raw (no finite tau_d). Tag: [DR].
    """
    a = 1.0 - g2_intrinsic
    d = 1.0 - g2_raw0
    if not (a > d > 0.0):
        raise ValueError("no finite dip recovery time for this raw/intrinsic pair")
    return tau_irf_ns * d / (a - d)


__all__ = ["tau_anchored", "inverted_background_g2", "g2_cw_high_pump_limit", "implied_dip_recovery_time"]
