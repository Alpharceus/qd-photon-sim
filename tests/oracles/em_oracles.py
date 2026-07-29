"""Independent EM oracles for cross-checking FSIM's owned physics (S-1).

FSIM owns its filter/cavity physics in fsim_core/spectral.py (lorentzian,
tophat_transmission, cavity_transmission -- a peak-normalized Lorentzian
cavity acceptance). The formulas ported here duplicate that physics, so
under the program's audit adoption they are allowed to exist ONLY as
independent test oracles under tests/oracles/ -- they must NEVER be
imported by fsim_core, fsim_viz, or fsim_gui (the "oracle-only rule").
verify/verify_fsim.py enforces this quarantine with a source grep.

Ported from the teammate `cavity_pr_opt` package per AUDIT.md S-1:
  - benchmarks/plane_wave_barrier.py (solve_slab_coefficients): the exact
    normal-incidence 3-region dielectric-slab transfer-matrix solution
    (tangential E/H continuity at both interfaces, R + T = 1).
  - kernels/coupled_mode.py (coupled_mode_amplitude / coupled_mode_transmission):
    the standard coupled-mode-theory (CMT) single-resonance Lorentzian
    response.

Both are textbook results (transfer-matrix / Fresnel thin-film optics --
e.g. Hecht, *Optics*, or Born & Wolf, *Principles of Optics*; and CMT --
e.g. Haus, *Waves and Fields in Optoelectronics*). The source package
ships no license file and targets Python >=3.10 (PEP 604 `X | Y` unions),
so this is a clean-room reimplementation of the underlying math for plain
Python 3.9 -- written fresh, not copied comment-for-comment. Its purpose
is to be an INDEPENDENT second implementation: it should only ever agree
with fsim_core by virtue of both being correct, never because one was
copy-pasted from the other.
"""
import cmath
import math


# --------------------------------------------------------------------- slab

def fresnel_slab(n1, n2, n3, wavelength, thickness):
    """Exact normal-incidence transfer-matrix solution for a 3-region
    dielectric slab (isotropic, lossless, nonmagnetic media), obtained via
    the multiple-beam (Airy / Fabry-Perot) summation -- a textbook route
    that is algebraically independent of a direct 4x4 boundary-condition
    solve, useful as a cross-check of either method.

    Geometry: propagation along +x, slab occupies 0 <= x <= thickness.
    Fields are incident-normalized, phasor convention exp(-i*omega*t):

        region 1 (x<0):            E = exp(i k1 x) + r exp(-i k1 x)
        region 2 (0<=x<=thickness): E = a exp(i k2 x) + b exp(-i k2 x)
        region 3 (x>thickness):     E = t exp(i k3 (x-thickness))

    with k_j = 2*pi*n_j/wavelength. For nonmagnetic media at normal
    incidence the tangential H (normalized to the region-1 admittance) is
    H = (n_j/n1) * (forward - backward) in every region, so continuity of
    tangential E and H at each interface gives four scalar equations; see
    the residuals returned below.

    Uses the standard Fresnel/Airy construction: single-interface amplitude
    coefficients r_ij = (n_i-n_j)/(n_i+n_j), t_ij = 1+r_ij, and the
    multiple-reflection sum inside the slab with round-trip phase
    exp(2i*beta), beta = k2*thickness.

    Returns a dict with:
      r, a, b, t   -- incident-normalized complex amplitudes (reflection;
                       internal forward/backward; transmission)
      R, T         -- reflected/transmitted power fractions (R=|r|^2,
                       T=(n3/n1)|t|^2); R+T=1 for lossless media
      residual_E_left, residual_H_left, residual_E_right, residual_H_right
                   -- |LHS-RHS| of the four tangential-continuity equations,
                      evaluated independently of the derivation above; these
                      should vanish at machine precision for a correct
                      solve, and are used by oracle_checks.py both as a
                      self-consistency check and as the target of its
                      negative control.
      k1, k2, k3, beta -- wavenumbers and slab phase thickness.
    """
    if n1 <= 0 or n2 <= 0 or n3 <= 0:
        raise ValueError("refractive indices must be positive")
    if wavelength <= 0:
        raise ValueError("wavelength must be positive")
    if thickness < 0:
        raise ValueError("thickness must be non-negative")

    k1 = 2.0 * math.pi * n1 / wavelength
    k2 = 2.0 * math.pi * n2 / wavelength
    k3 = 2.0 * math.pi * n3 / wavelength
    beta = k2 * thickness  # one-way phase thickness of the slab

    r12 = (n1 - n2) / (n1 + n2)
    r23 = (n2 - n3) / (n2 + n3)
    t12 = 1.0 + r12
    t23 = 1.0 + r23

    phase_rt = cmath.exp(2j * beta)   # round-trip phase inside the slab
    denom = 1.0 + r12 * r23 * phase_rt

    r = (r12 + r23 * phase_rt) / denom
    a = t12 / denom
    b = a * r23 * phase_rt
    t = t12 * t23 * cmath.exp(1j * beta) / denom

    R = abs(r) ** 2
    T = (n3 / n1) * abs(t) ** 2

    phase_fwd = cmath.exp(1j * beta)
    phase_bwd = cmath.exp(-1j * beta)
    residual_E_left = abs((1.0 + r) - (a + b))
    residual_H_left = abs((1.0 - r) - (n2 / n1) * (a - b))
    residual_E_right = abs((a * phase_fwd + b * phase_bwd) - t)
    residual_H_right = abs((n2 / n1) * (a * phase_fwd - b * phase_bwd) - (n3 / n1) * t)

    return {
        "r": r, "a": a, "b": b, "t": t,
        "R": R, "T": T,
        "residual_E_left": residual_E_left,
        "residual_H_left": residual_H_left,
        "residual_E_right": residual_E_right,
        "residual_H_right": residual_H_right,
        "k1": k1, "k2": k2, "k3": k3, "beta": beta,
    }


# ---------------------------------------------------------------------- CMT

def cmt_intracavity(omega, omega0=1.0, kappa_i=0.08, kappa_e=0.16, s_in=1.0):
    """Standard coupled-mode-theory (CMT) intracavity field amplitude for a
    single resonant mode side-coupled to one input port:

        a(omega) = sqrt(kappa_e) * s_in / [i*(omega0 - omega) + (kappa_i+kappa_e)/2]

    kappa_i is the intrinsic (loss) decay rate, kappa_e the external
    coupling rate (both angular-frequency units matching omega and omega0).
    |a(omega)|^2 is a Lorentzian in omega centered at omega0 with
    FWHM = kappa_i + kappa_e.

    `omega` may be a scalar or a numpy array; the expression is fully
    vectorized via ordinary arithmetic.
    """
    return (kappa_e ** 0.5) * s_in / (1j * (omega0 - omega) + 0.5 * (kappa_i + kappa_e))


def cmt_transmission(omega, omega0=1.0, kappa_i=0.08, kappa_e=0.16):
    """Standard CMT through-port power transmission for the same
    single-resonance system (the familiar notch/all-pass filter response):

        T(omega) = |1 - kappa_e / [i*(omega - omega0) + (kappa_i+kappa_e)/2]|^2

    On resonance (omega=omega0): T = [(kappa_i-kappa_e)/(kappa_i+kappa_e)]^2
    (the under/over/critically-coupled notch depth).
    """
    denom = 1j * (omega - omega0) + 0.5 * (kappa_i + kappa_e)
    return abs(1.0 - kappa_e / denom) ** 2
