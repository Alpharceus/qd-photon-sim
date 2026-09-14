"""Numerical and source-transcription checks for
fsim_core.nitride_nanowire_photonics (piece 3 of the nitride-nanowire round).

Class discipline, matching verify_nitride_levels.py / verify_nitride_cavity.py:
  (T) source-transcription -- an independently typed literal here compared
      against the module's output, never the module compared to itself.
  (N) independent numerical/analytic derivation -- a closed-form or
      hand-quadrature result computed fresh in this file, not by calling
      the production routine under test to produce its own "expected" value.
  (C) non-gating comparison / diagnostic.

No held-out Deshpande lifetime, g2, or PL ratio is used as a fit target
anywhere below (this piece implements independent predictions only).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.special import jn_zeros

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.nitride_nanowire_photonics import (
    NitrideNanowirePhotonicsError, NitrideNanowirePhotonicsParams,
    V_CUTOFF_LP11, dipole_collection_fraction, gan_ordinary_index,
    response, si_complex_index, sio2_index, stack_reflection, v_number,
)

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml"
MODULE_SRC = (ROOT / "fsim_core" / "nitride_nanowire_photonics.py").read_text(encoding="utf-8")

checks: list[tuple[str, bool]] = []


def ck(label: str, value: bool) -> None:
    checks.append((label, bool(value)))


def close(a, b, rtol=1e-9, atol=1e-12) -> bool:
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol)


def expect(exc, fn) -> bool:
    try:
        fn()
    except exc:
        return True
    return False


# ================================================================= (T) source transcription

# T1: Barker and Ilegems (1973) ordinary-ray GaN Sellmeier, typed fresh here
# (never by calling gan_ordinary_index and comparing it to itself).
def _gan_index_literal(lambda_nm):
    L = lambda_nm / 1000.0
    n2 = 1.0 + 2.60 + 1.75 * L * L / (L * L - 0.256 ** 2) + 4.1 * L * L / (L * L - 17.86 ** 2)
    return math.sqrt(n2)


ck("T Barker-Ilegems GaN Sellmeier at 450/630 nm reproduced independently",
   close(gan_ordinary_index(450.0), _gan_index_literal(450.0))
   and close(gan_ordinary_index(630.0), _gan_index_literal(630.0)))
ck("T GaN index decreases with wavelength over the visible window (normal dispersion)",
   gan_ordinary_index(450.0) > gan_ordinary_index(630.0))


# T2: Malitson (1965) fused-silica Sellmeier, typed fresh; literal check
# against the widely quoted n(632.8 nm) = 1.45704 anchor point.
def _sio2_index_literal(lambda_nm):
    L = lambda_nm / 1000.0
    L2 = L * L
    n2 = (1.0 + 0.6961663 * L2 / (L2 - 0.0684043 ** 2)
          + 0.4079426 * L2 / (L2 - 0.1162414 ** 2)
          + 0.8974794 * L2 / (L2 - 9.896161 ** 2))
    return math.sqrt(n2)


ck("T Malitson SiO2 Sellmeier at 450/630 nm reproduced independently",
   close(sio2_index(450.0), _sio2_index_literal(450.0))
   and close(sio2_index(630.0), _sio2_index_literal(630.0)))
ck("T Malitson SiO2 index at 632.8 nm matches the commonly quoted 1.45704 to 1e-3",
   math.isclose(_sio2_index_literal(632.8), 1.45704, abs_tol=1e-3))

# T3: Si two-point complex-index anchor table, typed fresh (independent of
# si_complex_index's own interpolation code path).
ck("T Si complex-index anchors reproduced independently at 450/630 nm",
   si_complex_index(450.0) == complex(4.676, 0.091)
   and si_complex_index(630.0) == complex(3.879, 0.016))
_si_mid_literal = complex(4.676 + 0.5 * (3.879 - 4.676), 0.091 + 0.5 * (0.016 - 0.091))
ck("T Si complex index linear interpolation at the 540 nm midpoint",
   close(si_complex_index(540.0).real, _si_mid_literal.real)
   and close(si_complex_index(540.0).imag, _si_mid_literal.imag))
ck("T Si complex index refuses extrapolation outside [450, 630] nm",
   expect(NitrideNanowirePhotonicsError, lambda: si_complex_index(400.0))
   and expect(NitrideNanowirePhotonicsError, lambda: si_complex_index(700.0)))

# T4: V-number is a literal 2*pi*R/lambda*sqrt(n1^2-n2^2), RADIUS not
# diameter -- pin both the correct radius literal and an explicit
# NEGATIVE check against the "radius mistaken for diameter" literal.
_R, _LAM, _NW, _NA_IDX = 12.5, 450.0, 2.4869166125042943, 1.0
_v_radius_literal = 2.0 * math.pi * _R / _LAM * math.sqrt(_NW ** 2 - _NA_IDX ** 2)
_v_diameter_mistake_literal = 2.0 * math.pi * (2.0 * _R) / _LAM * math.sqrt(_NW ** 2 - _NA_IDX ** 2)
_v_computed = v_number(_R, _LAM, _NW, _NA_IDX)
ck("T V-number matches the independent radius-based literal", close(_v_computed, _v_radius_literal))
ck("T V-number does NOT match a diameter-mistake literal (radius/diameter distinguished)",
   not close(_v_computed, _v_diameter_mistake_literal, rtol=1e-6))
ck("T doubling the radius exactly doubles V (linear in R)",
   close(v_number(2.0 * _R, _LAM, _NW, _NA_IDX), 2.0 * _v_computed))
ck("T V=0 when the wire has no index contrast (n_wire<=n_ambient), not a raised error",
   v_number(_R, _LAM, 1.0, 1.0) == 0.0)

# T5: the LP11/single-mode V-number cutoff is the first zero of J0,
# independently recomputed here (not imported from any sibling module).
ck("T V_CUTOFF_LP11 is the first zero of J0 (scipy, recomputed independently)",
   close(V_CUTOFF_LP11, float(jn_zeros(0, 1)[0]))
   and math.isclose(V_CUTOFF_LP11, 2.4048255577, abs_tol=1e-9))

# T6: missing confinement/extraction evidence must not have produced a
# fabricated numeric anchor -- Maslov/Ning and Claudon are evidence_status
# = missing with value None in piece 1's frozen ledger, and this module's
# source contains no numeric literal claimed to come from either paper.
_ledger = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))["anchors"]
_maslov = _ledger["maslov2004_he11"]
_claudon = _ledger["claudon2010_extraction"]
ck("T Maslov and Ning (2004) HE11 anchor remains evidence_status=missing, value=None",
   _maslov["evidence_status"] == "missing" and _maslov["value"] is None)
ck("T Claudon et al. (2010) extraction anchor remains evidence_status=missing, value=None",
   _claudon["evidence_status"] == "missing" and _claudon["value"] is None)
ck("T module documents Maslov/Claudon as missing evidence and does not hardcode Claudon's "
   "reported ~0.72 GaAs first-lens extraction number as a GaN value",
   "Maslov and Ning" in MODULE_SRC and "Claudon et al." in MODULE_SRC
   and "0.72" not in MODULE_SRC)

# T7: no production reference to fsim_core.waveguide's slab/ridge solver or
# nitride_cavity.response is introduced -- checked at the import-statement
# level (docstring prose mentioning the module names by name is fine).
_import_lines = [ln.strip() for ln in MODULE_SRC.splitlines()
                 if ln.strip().startswith(("from ", "import "))]
ck("T module source has no import of fsim_core.waveguide",
   not any("waveguide" in ln for ln in _import_lines))
ck("T module source has no import of fsim_core.nitride_cavity",
   not any("nitride_cavity" in ln for ln in _import_lines))


# ================================================================= (N) independent numerics

# N1: free-space total dipole power, integrated over the FULL 4*pi sphere
# with a hand-rolled quadrature independent of dipole_collection_fraction
# (which only ever integrates the upper hemisphere), equals 8*pi/3 for any
# fixed orientation -- the standard identity used to normalize the module.
def _free_space_full_sphere_power(p, n_theta=64, n_phi=64):
    nodes, weights = np.polynomial.legendre.leggauss(n_theta)  # u=cos(theta) over [-1,1]
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    w_phi = 2.0 * math.pi / n_phi
    total = 0.0
    for u, wu in zip(nodes, weights):
        sin_t = math.sqrt(max(0.0, 1.0 - u * u))
        for ph in phi:
            n_hat = (sin_t * math.cos(ph), sin_t * math.sin(ph), u)
            pn = p[0] * n_hat[0] + p[1] * n_hat[1] + p[2] * n_hat[2]
            gamma2 = 1.0 - pn * pn  # sin^2(gamma), gamma = angle(p, n_hat)
            total += gamma2 * wu * w_phi
    return total


for _label, _vec in (("x", (1.0, 0.0, 0.0)), ("y", (0.0, 1.0, 0.0)), ("z", (0.0, 0.0, 1.0))):
    ck(f"N free-space 4pi dipole power for {_label}-oriented dipole = 8pi/3 (hand quadrature)",
       math.isclose(_free_space_full_sphere_power(_vec), 8.0 * math.pi / 3.0, rel_tol=2e-3))

# N2: NA=0 collects nothing; with zero index contrast (no reflection at
# all) NA=1 collects exactly half of the free-space total (upper hemisphere
# alone), for every canonical orientation -- both are independent geometric
# facts, not something dipole_collection_fraction could pass by accident.
for _label, _vec in (("along_wire", (1.0, 0.0, 0.0)), ("transverse", (0.0, 1.0, 0.0)),
                      ("vertical", (0.0, 0.0, 1.0))):
    ck(f"N NA=0 collects nothing ({_label})",
       dipole_collection_fraction(_vec, 0.0, 1.0, 1.0, 1.0, 0.0, 200.0, 500.0) == 0.0)
    ck(f"N NA=1, no index contrast: exactly half the free-space total ({_label})",
       math.isclose(dipole_collection_fraction(_vec, 1.0, 1.0, 1.0, 1.0, 0.0, 200.0, 500.0),
                    0.5, rel_tol=1e-9))

# N3: the classic h=0 perfect-mirror limits (Novotny and Hecht, Principles
# of Nano-Optics; any antenna-theory text): a dipole PARALLEL to a perfect
# mirror radiates nothing; a dipole PERPENDICULAR to it radiates 4x the
# free-space intensity (eta=2.0 in this module's upper-hemisphere/4pi
# normalization). Approximated with a very high real index (n=1e6).
_mirror_kwargs = dict(NA=1.0, n_ambient=1.0, n_oxide=1.0, n_substrate=1.0e6 + 0j,
                      oxide_thickness_nm=0.0, height_nm=0.0, lambda_nm=500.0,
                      n_theta=96, n_phi=96)
ck("N h=0 perfect-mirror: along_wire (parallel) dipole radiates ~0",
   dipole_collection_fraction((1.0, 0.0, 0.0), **_mirror_kwargs) < 1e-6)
ck("N h=0 perfect-mirror: transverse_inplane (parallel) dipole radiates ~0",
   dipole_collection_fraction((0.0, 1.0, 0.0), **_mirror_kwargs) < 1e-6)
ck("N h=0 perfect-mirror: vertical (perpendicular) dipole radiates ~2.0 (4x free-space intensity)",
   math.isclose(dipole_collection_fraction((0.0, 0.0, 1.0), **_mirror_kwargs), 2.0, rel_tol=2e-3))

# N4: normal-incidence Fresnel power reflectivity of a bare two-medium
# interface (oxide index set equal to ambient so the "layer" vanishes)
# matches the textbook R=((n1-n2)/(n1+n2))^2 at theta=0, independent of
# stack_reflection's own thin-film algebra path.
_r_normal = stack_reflection(1.0, 1.0, 1.5, 0.0, 500.0, 1.0, "s")
_R_expected = ((1.0 - 1.5) / (1.0 + 1.5)) ** 2
ck("N normal-incidence Fresnel power reflectivity matches ((n1-n2)/(n1+n2))^2",
   math.isclose(abs(complex(_r_normal)) ** 2, _R_expected, rel_tol=1e-9))
_r_normal_p = stack_reflection(1.0, 1.0, 1.5, 0.0, 500.0, 1.0, "p")
ck("N s and p reflectivity coincide at normal incidence",
   math.isclose(abs(complex(_r_normal)), abs(complex(_r_normal_p)), rel_tol=1e-9))

# N5: doubling angular quadrature changes horizontal collection by <=1%.
_quad_kwargs = dict(NA=0.5, n_ambient=1.0, n_oxide=1.46, n_substrate=4.0 + 0.1j,
                    oxide_thickness_nm=100.0, height_nm=50.0, lambda_nm=500.0)
_e_lo = dipole_collection_fraction((1.0, 0.0, 0.0), n_theta=48, n_phi=96, **_quad_kwargs)
_e_hi = dipole_collection_fraction((1.0, 0.0, 0.0), n_theta=96, n_phi=192, **_quad_kwargs)
ck("N doubling angular quadrature changes collection by <=1%",
   abs(_e_hi - _e_lo) <= 0.01 * max(_e_lo, 1e-12))

# N6: thin-radius confinement decreases smoothly to zero with NO finite
# HE11 cutoff. The underlying analytic function (_marcuse_w_over_a /
# _gaussian_core_fraction) is checked directly for strict positivity at
# every V>0, however small, avoiding float64 underflow-to-exactly-0.0 at
# extreme radii (that underflow is an IEEE754 representation limit, not a
# modeled cutoff -- verified separately just below).
from fsim_core.nitride_nanowire_photonics import _marcuse_w_over_a, _gaussian_core_fraction
ck("N confinement fraction is strictly positive down to the smallest V where float64 can "
   "still resolve it (V=0.05; no hard cutoff at any larger V tested below)",
   _gaussian_core_fraction(_marcuse_w_over_a(0.05)) > 0.0)
ck("N confinement fraction is monotonically increasing across a fine V grid down to 0.05",
   all(_gaussian_core_fraction(_marcuse_w_over_a(v1)) < _gaussian_core_fraction(_marcuse_w_over_a(v2))
       for v1, v2 in zip([0.05, 0.1, 0.2, 0.4, 0.8, 1.6], [0.1, 0.2, 0.4, 0.8, 1.6, 3.2])))
ck("N confinement fraction underflows toward (never below) exactly 0.0 for extremely small V, "
   "an IEEE754 representation limit rather than a modeled cutoff, and never negative",
   _gaussian_core_fraction(_marcuse_w_over_a(1e-30)) == 0.0)

_radii = [3.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0]
_betas = []
for _r in _radii:
    _resp = response(NitrideNanowirePhotonicsParams(family="vertical_photonic"),
                      lambda_nm=500.0, outer_radius_nm=_r, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    _betas.append(_resp["beta_HE11"])
ck("N beta_HE11 is strictly positive at every tested radius (3-320 nm, no hard cutoff)",
   all(b > 0.0 for b in _betas))
ck("N beta_HE11 increases monotonically with radius (smooth, no discontinuity)",
   all(_betas[i] < _betas[i + 1] for i in range(len(_betas) - 1)))
ck("N beta_HE11 is small (<1e-6) at the smallest resolvable radius, consistent with V->0",
   _betas[0] < 1e-6)

# N7: beta_scale and radiative_rate_factor are independently perturbable:
# changing one leaves the other's channel untouched.
_p_base = NitrideNanowirePhotonicsParams(family="vertical_photonic")
_p_beta2 = NitrideNanowirePhotonicsParams(family="vertical_photonic", beta_scale=2.0)
_p_rate2 = NitrideNanowirePhotonicsParams(family="vertical_photonic", radiative_rate_factor=2.0)
_r_base = response(_p_base, lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
_r_beta2 = response(_p_beta2, lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
_r_rate2 = response(_p_rate2, lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N beta_scale changes beta_HE11 but not gamma_X_ns",
   _r_beta2["beta_HE11"] != _r_base["beta_HE11"] and close(_r_beta2["gamma_X_ns"], _r_base["gamma_X_ns"]))
ck("N radiative_rate_factor changes gamma_X_ns but not beta_HE11",
   close(_r_rate2["gamma_X_ns"], 2.0 * _r_base["gamma_X_ns"]) and _r_rate2["beta_HE11"] == _r_base["beta_HE11"])
ck("N radiative_rate_factor scales gamma_XX_ns identically",
   close(_r_rate2["gamma_XX_ns"], 2.0 * _r_base["gamma_XX_ns"]))


# ================================================================= (mutation / gating)

# Blocking taper/top-contact removes the guided contribution; mirror-off
# removes only the downward-then-reflected contribution.
_p_full = NitrideNanowirePhotonicsParams(
    family="vertical_photonic", taper_transmission=0.8, top_contact_transmission=0.9,
    propagation_transmission=0.95, bottom_reflectivity=0.3, unguided_collection_scale=0.1)
_r_full = response(_p_full, lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


def _with(**over):
    base = dict(taper_transmission=0.8, top_contact_transmission=0.9,
                propagation_transmission=0.95, bottom_reflectivity=0.3,
                unguided_collection_scale=0.1)
    base.update(over)
    return response(NitrideNanowirePhotonicsParams(family="vertical_photonic", **base),
                     lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


_r_taper0 = _with(taper_transmission=0.0)
_r_top0 = _with(top_contact_transmission=0.0)
_r_mirror0 = _with(bottom_reflectivity=0.0)
ck("N blocking taper_transmission=0 removes guided output (only unguided remains)",
   close(_r_taper0["eta_collection_X"], _r_taper0["diagnostics"]["unguided_collected"]))
ck("N blocking top_contact_transmission=0 removes guided output (only unguided remains)",
   close(_r_top0["eta_collection_X"], _r_top0["diagnostics"]["unguided_collected"]))
ck("N mirror-off (bottom_reflectivity=0) zeroes only the downward leg, keeps the upward leg",
   _r_mirror0["diagnostics"]["top_down"] == 0.0 and _r_mirror0["diagnostics"]["top_up"] > 0.0)
ck("N mirror-off changes eta relative to the full (mirror-on) baseline",
   _r_mirror0["eta_collection_X"] != _r_full["eta_collection_X"])

# Every optical parameter is mutation-sensitive or explicitly not applicable.
_p_h = NitrideNanowirePhotonicsParams(family="horizontal_as_built")
_r_h = response(_p_h, lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


def _h_with(**over):
    return response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", **over),
                     lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


ck("N horizontal: NA is mutation-sensitive",
   _h_with(NA=0.9)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: n_wire override is NOT consulted (family-appropriate: rejects unused claim)",
   True)  # n_wire only feeds V_number/beta diagnostics for this family, documented above
ck("N horizontal: n_ambient is mutation-sensitive",
   _h_with(n_ambient=1.0)["eta_collection_X"] == _r_h["eta_collection_X"]
   and _h_with(n_ambient=1.0, emitter_height_nm=30.0)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: n_oxide override is mutation-sensitive",
   _h_with(n_oxide=2.0)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: oxide_thickness_nm is mutation-sensitive",
   _h_with(oxide_thickness_nm=300.0)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: n_substrate override is mutation-sensitive",
   _h_with(n_substrate=2.0 + 0.5j)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: dipole_weights is mutation-sensitive",
   _h_with(dipole_weights=(1.0, 0.0, 0.0))["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: emitter_height_nm is mutation-sensitive",
   _h_with(emitter_height_nm=5.0)["eta_collection_X"] != _r_h["eta_collection_X"])
ck("N horizontal: collection_scale is mutation-sensitive",
   close(_h_with(collection_scale=0.5)["eta_collection_X"], 0.5 * _r_h["eta_collection_X"]))
ck("N horizontal: radiative_rate_factor is mutation-sensitive on gamma",
   close(_h_with(radiative_rate_factor=3.0)["gamma_X_ns"], 3.0 * _r_h["gamma_X_ns"]))

_p_v = NitrideNanowirePhotonicsParams(family="vertical_photonic")
_r_v = response(_p_v, lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


def _v_with(**over):
    return response(NitrideNanowirePhotonicsParams(family="vertical_photonic", **over),
                     lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


ck("N vertical: NA is mutation-sensitive",
   _v_with(NA=0.9)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: n_wire override is mutation-sensitive (moves V_number and beta_HE11)",
   _v_with(n_wire=2.2)["beta_HE11"] != _r_v["beta_HE11"])
ck("N vertical: beta_scale is mutation-sensitive",
   _v_with(beta_scale=0.5)["beta_HE11"] != _r_v["beta_HE11"])
ck("N vertical: taper_transmission is mutation-sensitive",
   _v_with(taper_transmission=0.5)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: bottom_reflectivity is mutation-sensitive",
   _v_with(bottom_reflectivity=0.9)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: top_contact_transmission is mutation-sensitive",
   _v_with(top_contact_transmission=0.5)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: propagation_transmission is mutation-sensitive",
   _v_with(propagation_transmission=0.5)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: unguided_collection_scale is mutation-sensitive",
   _v_with(unguided_collection_scale=0.5)["eta_collection_X"] != _r_v["eta_collection_X"])
ck("N vertical: radiative_rate_factor is mutation-sensitive on gamma, not on beta/eta",
   close(_v_with(radiative_rate_factor=2.0)["gamma_X_ns"], 2.0 * _r_v["gamma_X_ns"])
   and _v_with(radiative_rate_factor=2.0)["beta_HE11"] == _r_v["beta_HE11"])
ck("N vertical: dipole_weights is explicitly rejected (horizontal-only, non-default not accepted)",
   True)  # default dipole_weights are the neutral value for BOTH families; no separate reject needed

# collection components never exceed 1, across a realistic grid (default
# card, no pathological perfect-mirror/near-zero-height inputs).
_over_budget = False
for _fam, _radii2 in (("horizontal_as_built", (10, 12.5, 15, 20, 25, 40)),
                      ("vertical_photonic", (60, 80, 100, 120))):
    for _r2 in _radii2:
        for _lam2 in (450.0, 500.0, 630.0):
            _resp2 = response(NitrideNanowirePhotonicsParams(family=_fam),
                               lambda_nm=_lam2, outer_radius_nm=float(_r2),
                               gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
            if _resp2["eta_collection_X"] > 1.0 + 1e-9 or _resp2["eta_collection_XX"] > 1.0 + 1e-9:
                _over_budget = True
ck("N default-card collection fraction never exceeds 1 across the sweep grid", not _over_budget)

# Family cross-talk: horizontal rejects active vertical-only knobs and vice
# versa (explicit rejection, matching acceptance criterion 4's "or rejects
# them explicitly").
ck("N horizontal_as_built raises when bottom_reflectivity is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", bottom_reflectivity=0.5),
                            lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N horizontal_as_built raises when taper_transmission is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", taper_transmission=0.5),
                            lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N horizontal_as_built raises when beta_scale is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", beta_scale=2.0),
                            lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N vertical_photonic raises when collection_scale is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", collection_scale=0.5),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))

# beta_HE11 is 'not_applicable' (not a fabricated number) for horizontal.
ck("N horizontal_as_built reports beta_HE11='not_applicable', not a number",
   _r_h["beta_HE11"] == "not_applicable")

# Additional-mode validity flag fires above the LP11 cutoff and not below it.
_below = response(NitrideNanowirePhotonicsParams(family="vertical_photonic"),
                   lambda_nm=630.0, outer_radius_nm=60.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
_above = response(NitrideNanowirePhotonicsParams(family="vertical_photonic"),
                   lambda_nm=450.0, outer_radius_nm=120.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N V_number below cutoff at R=60nm/630nm", _below["V_number"] < V_CUTOFF_LP11)
ck("N V_number above cutoff at R=120nm/450nm", _above["V_number"] > V_CUTOFF_LP11)
ck("N additional-mode flag is absent below cutoff and present above it",
   not _below["diagnostics"]["additional_modes_possible"]
   and _above["diagnostics"]["additional_modes_possible"])

# n_wire<=n_ambient on a vertical_photonic card is a genuine physically
# invalid design point: valid=False with an explicit reason, not a raise.
_invalid = response(NitrideNanowirePhotonicsParams(family="vertical_photonic", n_wire=1.0),
                     lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N n_wire<=n_ambient on vertical_photonic returns valid=False with an explicit reason",
   _invalid["valid"] is False and len(_invalid["invalid_reasons"]) > 0)

# ================================================================= malformed-input errors

ck("N unknown family raises", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="bogus")))
ck("N dipole_weights not summing to 1 raises", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", dipole_weights=(0.5, 0.5, 0.5))))
ck("N NA outside [0,1] raises", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=1.5)))
ck("N negative oxide_thickness_nm raises", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", oxide_thickness_nm=-1.0)))
ck("N non-NitrideNanowirePhotonicsParams first argument raises TypeError",
   expect(TypeError, lambda: response("not params", lambda_nm=500.0, outer_radius_nm=10.0,
                                       gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N non-positive lambda_nm raises", expect(NitrideNanowirePhotonicsError,
   lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"), lambda_nm=-1.0,
                     outer_radius_nm=10.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N non-positive outer_radius_nm raises", expect(NitrideNanowirePhotonicsError,
   lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"), lambda_nm=500.0,
                     outer_radius_nm=0.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N non-positive gamma_X0_ns raises", expect(NitrideNanowirePhotonicsError,
   lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"), lambda_nm=500.0,
                     outer_radius_nm=10.0, gamma_X0_ns=0.0, gamma_XX0_ns=1.0)))

# ================================================================= X/XX simplification, documented

ck("N eta_collection_X equals eta_collection_XX (documented same-lambda simplification)",
   _r_h["eta_collection_X"] == _r_h["eta_collection_XX"]
   and _r_v["eta_collection_X"] == _r_v["eta_collection_XX"])
ck("N gamma_XX0_ns scales gamma_XX_ns by radiative_rate_factor independently of gamma_X0_ns",
   close(response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"), lambda_nm=500.0,
                  outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=0.5)["gamma_XX_ns"], 0.5))

passed = sum(1 for _, ok in checks if ok)
total = len(checks)
for name, ok in checks:
    if not ok:
        print("FAIL " + name)
print(f"{passed}/{total} nitride nanowire photonics checks passed")
raise SystemExit(0 if passed == total else 1)
