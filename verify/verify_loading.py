"""Verification for fsim_core.loading's cap-2 Poisson loading at small mu.

Regression pin: P2 = 1 - P0 - P1 cancelled catastrophically for mu < ~1e-3,
so f1b_g2 returned garbage (g2_op = 111 at mu = 5.4e-10) and 0 instead of its
F1 limit eps at mu = 0.

Expectations are NOT produced by the code under test:
  * P(n>=2) of a Poisson(mu) law, evaluated in 50-digit arithmetic (stdlib decimal,
    prec = 50) as the exact series e^-mu sum_{k>=2} mu^k/k!  [DR, Poisson identity].
  * lim_{mu->0} f1b_g2 = eps (F1 identity) and the small-mu expansion
    g2/eps = 1 + mu (1/3 - eps) + O(mu^2) (module docstring)  [DR].
  * Bit-identity with the legacy expression 1 - e^-mu - mu e^-mu (and the
    legacy f1b_g2) for mu >= 1e-3 (protects verify_fsim / audit_physics).
"""
import sys
import warnings
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.loading import f1b_g2, loading_probs  # noqa: E402

getcontext().prec = 50
checks = []


def ok(name, cond, detail=""):
    checks.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def p2_exact(mu):
    """Exact Poisson P(n>=2) in 50-digit arithmetic, summed as the series."""
    # Decimal(float) is the exact binary value of mu; the tail
    # sum_{k>=2} mu^k/k! is summed term by term until a term no longer
    # changes the 50-digit sum (mu <= 0.3 here, so convergence is fast).
    m = Decimal(mu)
    term = m * m / 2
    tail = Decimal(0)
    k = 2
    while True:
        new = tail + term
        if new == tail:
            break
        tail = new
        k += 1
        term = term * m / k
    return (-m).exp() * tail


def legacy_probs(mu):
    mu = np.asarray(mu, dtype=float)
    P0 = np.exp(-mu)
    P1 = mu * np.exp(-mu)
    return P0, P1, 1.0 - P0 - P1


def legacy_f1b(mu, eps):
    _, P1, P2 = legacy_probs(mu)
    denom = P1 + P2 * (1.0 + eps)
    return np.where(denom > 0, 2.0 * P2 * eps / denom**2, 0.0)


def main():
    # 1. P2 vs exact series at small mu, 1e-10 relative.
    for mu in (1e-12, 1e-9, 1e-7, 1e-5):
        p2 = float(loading_probs(mu)[2])
        ref = p2_exact(mu)
        rel = float(abs((Decimal(p2) - ref) / ref))
        ok(f"P2(mu={mu:g}) matches exact Poisson tail to 1e-10 rel", rel < 1e-10,
           f"P2={p2:.12e}, exact={float(ref):.12e}, rel={rel:.1e}")
    # also near and above the branch point (legacy side: loose, cancellation-limited)
    for mu in (9.99e-4, 1e-3, 0.3):
        p2 = float(loading_probs(mu)[2])
        rel = float(abs((Decimal(p2) - p2_exact(mu)) / p2_exact(mu)))
        ok(f"P2(mu={mu:g}) continuous across branch (1e-8 rel)", rel < 1e-8, f"rel={rel:.1e}")

    # 2. f1b_g2 mu -> 0 limit = eps.
    for eps in (0.01, 0.05, 0.3, 1.0):
        with warnings.catch_warnings():
            warnings.simplefilter("error")          # no 0/0 warnings at mu = 0
            g0 = float(f1b_g2(0.0, eps))
        g_tiny = float(f1b_g2(1e-15, eps))
        ok(f"f1b_g2(mu=0, eps={eps}) == eps to 1e-9", abs(g0 - eps) < 1e-9, f"{g0!r}")
        ok(f"f1b_g2(mu=1e-15, eps={eps}) == eps to 1e-9", abs(g_tiny - eps) < 1e-9, f"{g_tiny!r}")

    # 3. f1b_g2 in [0, 1+eps] for mu in [1e-12, 1e-2], and follows the
    #    analytic small-mu expansion 1 + mu(1/3 - eps) (error O(mu^2)).
    mus = np.logspace(-12, -2, 201)
    for eps in (0.01, 0.05, 0.3, 1.0):
        g = f1b_g2(mus, eps)
        ok(f"f1b_g2 in [0, 1+eps] for mu in [1e-12, 1e-2], eps={eps}",
           np.all(np.isfinite(g)) and np.all(g >= 0) and np.all(g <= 1 + eps),
           f"min={g.min():.6g}, max={g.max():.6g}")
        sel = mus <= 1e-4
        dev = g[sel] / eps - 1.0 - mus[sel] * (1.0 / 3.0 - eps)
        bound = 2.0 * mus[sel] ** 2 + 1e-14          # O(mu^2) remainder + rounding
        ok(f"f1b_g2/eps = 1 + mu(1/3-eps) + O(mu^2) for mu <= 1e-4, eps={eps}",
           np.all(np.abs(dev) <= bound), f"max|dev|={np.abs(dev).max():.1e}")

    # 4. Bit-identity with the legacy formula for mu >= 1e-3.
    rng = np.random.default_rng(12345)
    mu50 = np.sort(np.concatenate([[1e-3, 20.0], 10 ** rng.uniform(-3, np.log10(20), 48)]))
    same_scalar = all(
        all(float(a) == float(b) for a, b in zip(loading_probs(m), legacy_probs(m)))
        for m in mu50)
    ok("loading_probs == legacy exactly for 50 scalar mu in [1e-3, 20]", same_scalar)
    same_arr = all(np.array_equal(a, b) for a, b in zip(loading_probs(mu50), legacy_probs(mu50)))
    ok("loading_probs == legacy exactly for the 50-point array", same_arr)
    same_g2 = all(
        np.array_equal(f1b_g2(mu50, e), legacy_f1b(mu50, e))
        and all(float(f1b_g2(m, e)) == float(legacy_f1b(m, e)) for m in mu50)
        for e in (0.01, 0.05, 0.3, 1.0))
    ok("f1b_g2 == legacy exactly for mu in [1e-3, 20] (scalar and array)", same_g2)
    # mixed array: large-mu entries still legacy bit-for-bit
    mixed = np.concatenate([[0.0, 1e-9, 5e-4], mu50])
    ok("mixed small/large array keeps legacy entries bit-identical",
       np.array_equal(loading_probs(mixed)[2][3:], legacy_probs(mu50)[2])
       and np.array_equal(f1b_g2(mixed, 0.05)[3:], legacy_f1b(mu50, 0.05)))

    print(f"{sum(checks)}/{len(checks)} loading checks passed")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
