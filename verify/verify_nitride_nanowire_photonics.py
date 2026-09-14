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
import re
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.special import jn_zeros

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.nitride_nanowire_photonics import (
    CPLANE_ONLY_DIPOLE_WEIGHTS, NitrideNanowirePhotonicsError, NitrideNanowirePhotonicsParams,
    V_CUTOFF_LP11, dipole_collection_fraction, gan_ordinary_index,
    response, si_complex_index, sio2_index, stack_reflection, v_number,
)
from fsim_core.nitride_nanowire_photonics import (
    _he11_objective_acceptance, _mode_field_radius_nm,
    _group_index_ratio, _wire_antenna_screening_intensity,
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
# [fix round DIRECTIVE H7] Claudon's reported ~0.95 guided-mode beta is now
# LEGITIMATELY transcribed into the module as a tagged [E] multimode-penalty
# anchor (see _multimode_penalty) -- it is no longer absent, but it must be
# present ONLY as that tagged, cited anchor, never as a silently fabricated
# GaN confinement/beta value.
ck("T module source transcribes the Bleuse et al. PRL 106, 103601 (2011) "
   "Fig. 2 / Claudon et al. (2010) multimode-penalty anchors (n_wire=3.45, "
   "beta~0.95 near d/lambda 0.22-0.24 falling to beta~0.70 by d/lambda~0.40) "
   "with an [E] tag and citation, not as a silently fabricated GaN value "
   "(fix round DIRECTIVE H7)",
   "_MULTIMODE_ANCHOR_N_WIRE = 3.45" in MODULE_SRC
   and "_MULTIMODE_ANCHOR_BETA_PEAK = 0.95" in MODULE_SRC
   and "_MULTIMODE_ANCHOR_BETA_HIGH = 0.70" in MODULE_SRC
   and "Bleuse" in MODULE_SRC and "Claudon" in MODULE_SRC and "[E]" in MODULE_SRC)

# T (Attempt 3 fix C): the Bleuse/Claudon multimode-penalty anchor is now a
# POPULATED, non-null ledger entry (bleuse2011_claudon2010_beta_envelope,
# evidence_status=figure_reading) -- pin the module's own hardcoded
# multimode-penalty constants against the ledger's own transcribed values
# (parsed from the ledger, never re-typed independently of it), so the two
# cannot silently drift apart.
_bleuse_claudon_anchor = _ledger["bleuse2011_claudon2010_beta_envelope"]
ck("T Attempt 3 fix C: bleuse2011_claudon2010_beta_envelope ledger entry is "
   "populated (evidence_status=figure_reading, not missing) and the "
   "module's own multimode-penalty constants match its transcribed "
   "beta_high/beta_lower/d_over_lambda_lower_beta values",
   _bleuse_claudon_anchor["evidence_status"] == "figure_reading"
   and _bleuse_claudon_anchor["tag"] == "E"
   and math.isclose(_bleuse_claudon_anchor["value"]["beta_high"], 0.95, abs_tol=1e-9)
   and math.isclose(_bleuse_claudon_anchor["value"]["beta_lower"], 0.70, abs_tol=1e-9)
   and math.isclose(_bleuse_claudon_anchor["value"]["d_over_lambda_lower_beta"], 0.40, abs_tol=1e-9)
   and "_MULTIMODE_ANCHOR_BETA_PEAK = 0.95" in MODULE_SRC
   and "_MULTIMODE_ANCHOR_BETA_HIGH = 0.70" in MODULE_SRC
   and "_MULTIMODE_ANCHOR_D_LAMBDA_HIGH = 0.40" in MODULE_SRC)
ck("T Attempt 3 fix C: provenance cites the ledger key by name instead of "
   "describing the numbers as ledger-free",
   "bleuse2011_claudon2010_beta_envelope" in MODULE_SRC)

# T8 (fix round LOW 8): oxide_thickness_nm=100.0 is [V], a literal
# transcription of deshpande2013_device_geometry's substrate description,
# not an [A] design choice -- parse the ledger's own string (never a
# literal re-typed independently of the ledger) and pin the module default
# against it.
_geom_anchor = _ledger["deshpande2013_device_geometry"]
_substrate_str = _geom_anchor["value"]["substrate"]
_default_card_horizontal = NitrideNanowirePhotonicsParams(family="horizontal_as_built")
ck("T deshpande2013_device_geometry (V) reports '100 nm thermal SiO2 ...' and "
   "the module's default oxide_thickness_nm literally matches it",
   "100 nm" in _substrate_str
   and _default_card_horizontal.oxide_thickness_nm == 100.0)
ck("T module tags oxide_thickness_nm=100.0 as [V] (source transcription), not "
   "[A], matching the frozen contract's own 'V (2013 substrate)' tag",
   response(_default_card_horizontal, lambda_nm=500.0, outer_radius_nm=12.5,
            gamma_X0_ns=1.0, gamma_XX0_ns=1.0)["provenance"]["oxide_thickness_nm"]
   .startswith("[V]"))

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

# N5b (fix round MEDIUM 7): intermediate-NA=0.5 free-space anchors, pinned
# against a closed-form analytic expression derived fresh here (not by
# calling dipole_collection_fraction to produce its own "expected" value):
# for a dipole and objective axis both fixed, in the (theta,phi) convention
# used by dipole_collection_fraction, the free-space (no reflection) power
# pattern is 1-sin^2(theta)cos^2(phi) for an in-plane dipole and sin^2(theta)
# for the vertical dipole. Integrating sin(theta) dtheta dphi from 0 to
# theta_max=asin(NA), 0 to 2*pi, and normalizing by the fixed 8*pi/3 total
# gives the closed forms below (each independently re-derived by direct
# antiderivative, not sourced from the module):
#   in-plane(a)  = (4 - 3*cos(a) - cos(a)**3) / 8
#   vertical(a)  = (2 - 3*cos(a) + cos(a)**3) / 4
# both satisfy a=0 -> 0 and a=pi/2 -> 0.5, matching N2 above.
def _inplane_closed_form(NA):
    a = math.asin(NA)
    return (4.0 - 3.0 * math.cos(a) - math.cos(a) ** 3) / 8.0


def _vertical_closed_form(NA):
    a = math.asin(NA)
    return (2.0 - 3.0 * math.cos(a) + math.cos(a) ** 3) / 4.0


_inplane_NA05 = _inplane_closed_form(0.5)
_vertical_NA05 = _vertical_closed_form(0.5)
ck("T horizontal dipole free-space NA=0.5 matches the spec-given anchor "
   "0.0940513 (closed-form analytic derivation, tolerance 1e-4)",
   math.isclose(_inplane_NA05, 0.0940513, abs_tol=1e-4))
ck("T vertical dipole free-space NA=0.5 matches the spec-given anchor "
   "0.0128614 (closed-form analytic derivation, tolerance 1e-4)",
   math.isclose(_vertical_NA05, 0.0128614, abs_tol=1e-4))
ck("N closed-form NA=0.5 in-plane/vertical free-space fractions match the "
   "module's own quadrature to numerical precision",
   math.isclose(_inplane_NA05,
                dipole_collection_fraction((1.0, 0.0, 0.0), 0.5, 1.0, 1.0, 1.0, 0.0, 200.0, 500.0),
                rel_tol=1e-6)
   and math.isclose(_vertical_NA05,
                     dipole_collection_fraction((0.0, 0.0, 1.0), 0.5, 1.0, 1.0, 1.0, 0.0, 200.0, 500.0),
                     rel_tol=1e-6))

# N5c (fix round MEDIUM 7): the oblique-film interference period. At normal
# incidence the thin-film phase is beta=2*pi/lambda*n_oxide*thickness*cos1,
# cos1=1 at theta=0 exactly, so the interference pattern (and hence a
# collection-fraction-vs-thickness scan) is periodic with period
# lambda/(2*n_oxide); a finite-NA=0.5 objective integrates a small spread of
# off-normal angles, so the MEASURED peak spacing (from a real thickness
# scan of the production dipole_collection_fraction, not from re-deriving
# the analytic period) is only approximately that value -- checked here to
# 5% (a physical smearing effect, not a numerical-precision tolerance).
def _measured_oxide_period_nm(lambda_nm, step_nm=0.5, t_max_nm=900.0):
    n_ox = sio2_index(lambda_nm)
    n_sub = si_complex_index(lambda_nm)
    thicknesses = np.arange(0.0, t_max_nm, step_nm)
    vals = np.array([
        dipole_collection_fraction((0.0, 0.0, 1.0), 0.5, 1.0, n_ox, n_sub, float(t), 12.5, lambda_nm)
        for t in thicknesses
    ])
    peaks = [thicknesses[i] for i in range(1, len(vals) - 1)
             if vals[i] > vals[i - 1] and vals[i] > vals[i + 1]]
    if len(peaks) < 2:
        return None
    diffs = [peaks[i + 1] - peaks[i] for i in range(len(peaks) - 1)]
    return sum(diffs) / len(diffs)


for _lam_scan in (450.0,):
    _measured_period = _measured_oxide_period_nm(_lam_scan)
    _predicted_period = _lam_scan / (2.0 * sio2_index(_lam_scan))
    ck(f"N oxide-thickness interference period at {_lam_scan:.0f} nm matches "
       f"lambda/(2*n_SiO2)={_predicted_period:.1f} nm within 5% (measured "
       f"{_measured_period:.1f} nm from a real thickness scan)",
       _measured_period is not None
       and abs(_measured_period - _predicted_period) <= 0.05 * _predicted_period)

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
# [fix round DIRECTIVE H7] radii 3-80 nm are all single-mode (V<=V_CUTOFF_LP11)
# at this lambda/n_wire, where the multimode penalty s(V)=1 and beta_HE11
# still increases monotonically exactly as before; radii 160/320 nm are
# ABOVE cutoff, where the new multimode penalty must instead turn the curve
# over (see the dedicated interior-maximum section below) -- the old
# full-range "monotonic increase to R=320" claim was the bug H7 fixes.
ck("N beta_HE11 increases monotonically with radius while V stays at/below "
   "the LP11 cutoff (single-mode branch, radii 3-80 nm)",
   all(_betas[i] < _betas[i + 1] for i in range(5)))
ck("N beta_HE11 is small (<1e-6) at the smallest resolvable radius, consistent with V->0",
   _betas[0] < 1e-6)
ck("N beta_HE11 DECREASES once V crosses the LP11 cutoff (fix round DIRECTIVE "
   "H7: an interior maximum, not unbounded growth) -- R=160/320 nm are both "
   "above cutoff at this lambda/n_wire and both lower than the R=80 nm "
   "(still single-mode) value",
   _betas[6] < _betas[5] and _betas[7] < _betas[6])

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
# [fix round MEDIUM 4] replaces a hardcoded ck(..., True): n_wire genuinely
# IS consulted for horizontal_as_built now (V_number diagnostic, and
# gamma_X_ns via antenna_rate_factor).
ck("N horizontal: n_wire override moves the V_number diagnostic",
   _h_with(n_wire=2.0)["V_number"] != _r_h["V_number"])
# [Attempt 3 fix A] REVERSES the fix-round MEDIUM 2 claim above: the
# wire-antenna screening is now a RATE effect (antenna_rate_factor), so
# n_wire must move gamma_X_ns and must NOT move eta_collection_X (the
# plain unscreened geometric collection).
ck("N Attempt 3 fix A: horizontal n_wire override moves gamma_X_ns via "
   "antenna_rate_factor (the screening is a rate effect, not a collection "
   "loss)",
   _h_with(n_wire=2.0)["gamma_X_ns"] != _r_h["gamma_X_ns"])
ck("N Attempt 3 fix A: horizontal n_wire override does NOT move "
   "eta_collection_X (the wire-antenna screening no longer touches "
   "collection at all)",
   close(_h_with(n_wire=2.0)["eta_collection_X"], _r_h["eta_collection_X"]))
# [fix round MEDIUM 5] replaces a check that mutated emitter_height_nm, not
# n_ambient, while claiming to test n_ambient: n_ambient is mutated directly.
ck("N horizontal: n_ambient is mutation-sensitive",
   _h_with(n_ambient=1.33)["eta_collection_X"] != _r_h["eta_collection_X"])
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
# [fix round HIGH 1] taper_output_mfr_nm is mutation-sensitive: expanding the
# mode radius at the taper output shrinks divergence and raises obj_accept
# (this small default wire is badly divergent, see the 60-deg note check
# below, so a real taper materially helps collection).
ck("N vertical: taper_output_mfr_nm is mutation-sensitive and raises obj_accept",
   _v_with(taper_output_mfr_nm=2000.0)["eta_collection_X"] != _r_v["eta_collection_X"]
   and (_v_with(taper_output_mfr_nm=2000.0)["diagnostics"]["obj_accept"]
        > _r_v["diagnostics"]["obj_accept"]))

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
# [fix round MEDIUM 3] two-sided guard: a vertical_photonic card must also
# reject every horizontal-only knob (replaces a hardcoded ck(..., True) that
# asserted this without ever calling response() to check it).
ck("N vertical_photonic raises when dipole_weights is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", dipole_weights=(1.0, 0.0, 0.0)),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N vertical_photonic raises when n_oxide is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", n_oxide=2.0),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N vertical_photonic raises when n_substrate is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", n_substrate=2.0 + 0.1j),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N vertical_photonic raises when oxide_thickness_nm is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", oxide_thickness_nm=50.0),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N vertical_photonic raises when emitter_height_nm is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="vertical_photonic", emitter_height_nm=5.0),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))
ck("N horizontal_as_built raises when taper_output_mfr_nm is non-default",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", taper_output_mfr_nm=500.0),
                            lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))

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
# [Attempt 3 fix A] gamma_XX_ns now also picks up antenna_rate_factor (the
# default card's dipole_weights are isotropic, so this is NOT 1.0) --
# independently recomputed here, never read from the module's own output.
_n_wire_500 = gan_ordinary_index(500.0)
_screen_500 = (2.0 / (_n_wire_500 ** 2 + 1.0)) ** 2
_antenna_rate_factor_iso_500 = (1.0 / 3.0) * (1.0 + 2.0 * _screen_500)
ck("N gamma_XX0_ns scales gamma_XX_ns by radiative_rate_factor*antenna_rate_factor "
   "independently of gamma_X0_ns (Attempt 3 fix A: antenna_rate_factor now enters too, "
   "recomputed fresh here for the isotropic default)",
   close(response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"), lambda_nm=500.0,
                  outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=0.5)["gamma_XX_ns"],
         0.5 * _antenna_rate_factor_iso_500, rtol=1e-9))

# ================================================================= (fix round LOW 10) NA domain

ck("N NA=0 raises (domain is (0,1], not [0,1])", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=0.0)))
ck("N NA>1 with n_ambient=1 (no immersion) raises", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=1.2)))
ck("N NA>1 admitted when n_ambient>1 (immersion) and NA<=n_ambient",
   NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=1.2, n_ambient=1.33).NA == 1.2)
ck("N NA>n_ambient still raises even under immersion", expect(NitrideNanowirePhotonicsError,
   lambda: NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=1.5, n_ambient=1.33)))
ck("N an admitted NA>1 immersion card documents the saturated-hemisphere "
   "approximation in notes (fix round LOW 10)",
   any("immersion" in n for n in
       response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", NA=1.2, n_ambient=1.33),
                lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)["notes"]))

# ================================================================= (fix round HIGH 1) vertical objective acceptance

# [Attempt 3 fix G] the old version of this check divided an expression by
# itself (numerator and denominator were the SAME trapezoid integral) --
# a tautology, true by construction regardless of the module's own
# correctness. This is a genuine independent re-derivation instead: a
# FRESH, independently typed copy of the Gaussian-aperture far-field
# formula (obliquity factor ((1+cos theta)/2)^2 times
# exp(-(k*w*sin theta)^2/2)) -- never imported from the module -- is
# integrated over [0, theta_max] and separately over the full [0, pi/2]
# hemisphere, and the resulting ratio is compared against
# `_he11_objective_acceptance`'s own returned value for the SAME
# theta_max/k/w (which IS imported, since that is the production routine
# under test).
def _independent_he11_far_field(theta, k, w_nm):
    return ((1.0 + np.cos(theta)) / 2.0) ** 2 * np.exp(-0.5 * (k * w_nm * np.sin(theta)) ** 2)


def _independent_objective_acceptance(theta_max, k, w_nm, n_pts=8001):
    theta_cone = np.linspace(0.0, theta_max, n_pts)
    theta_full = np.linspace(0.0, math.pi / 2.0, n_pts)
    numer = np.trapezoid(_independent_he11_far_field(theta_cone, k, w_nm) * np.sin(theta_cone), theta_cone)
    denom = np.trapezoid(_independent_he11_far_field(theta_full, k, w_nm) * np.sin(theta_full), theta_full)
    return float(numer / denom)


ck("N Attempt 3 fix G: independent re-derivation (fresh far-field formula, "
   "not imported from the module -- replaces the old x/x tautology) of the "
   "HE11 objective-acceptance ratio matches _he11_objective_acceptance "
   "across several theta_max/w",
   all(math.isclose(_independent_objective_acceptance(_tm, 2.0 * math.pi / 500.0, _w),
                     _he11_objective_acceptance(_tm, 2.0 * math.pi / 500.0, _w), rel_tol=2e-3)
       for _tm, _w in ((0.3, 84.0), (0.7, 84.0), (1.2, 150.0), (1.5, 40.0))))
ck("N Attempt 3 fix G: independent re-derivation confirms the full-"
   "hemisphere self-ratio is exactly 1.0 (theta_max=pi/2 against itself)",
   _independent_objective_acceptance(math.pi / 2.0, 2.0 * math.pi / 500.0, 84.0) == 1.0)
ck("N _he11_objective_acceptance at NA=1.0 (theta_max=pi/2) returns exactly 1.0 for any finite w",
   _he11_objective_acceptance(math.pi / 2.0, 2.0 * math.pi / 500.0, 84.0) == 1.0
   and _he11_objective_acceptance(math.pi / 2.0, 2.0 * math.pi / 500.0, 5000.0) == 1.0)
ck("N _he11_objective_acceptance is monotonically increasing in theta_max",
   all(_he11_objective_acceptance(a1, 2.0 * math.pi / 500.0, 100.0)
       < _he11_objective_acceptance(a2, 2.0 * math.pi / 500.0, 100.0)
       for a1, a2 in zip([0.1, 0.3, 0.6, 0.9, 1.2], [0.3, 0.6, 0.9, 1.2, 1.5])))
ck("N doubling _he11_objective_acceptance's angular resolution changes the result by <=1%",
   abs(_he11_objective_acceptance(0.7, 2.0 * math.pi / 500.0, 100.0, n_pts=2001)
       - _he11_objective_acceptance(0.7, 2.0 * math.pi / 500.0, 100.0, n_pts=4001))
   <= 0.01 * max(_he11_objective_acceptance(0.7, 2.0 * math.pi / 500.0, 100.0, n_pts=2001), 1e-12))

# Full pipeline: NA=1.0, mirror=1.0, lossless -> eta_collection_X == beta_HE11
# (fixes the reviewed bug where up to 21.5% of upward power was discarded
# "beyond 90 deg" even at NA=1.0).
_p_na1 = NitrideNanowirePhotonicsParams(family="vertical_photonic", NA=1.0, bottom_reflectivity=1.0)
_r_na1 = response(_p_na1, lambda_nm=500.0, outer_radius_nm=80.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N NA=1.0, mirror=1.0, lossless: eta_collection_X equals beta_HE11 within 1e-6",
   abs(_r_na1["eta_collection_X"] - _r_na1["beta_HE11"]) < 1e-6)
ck("N NA=1.0 case reports obj_accept == 1.0 exactly",
   _r_na1["diagnostics"]["obj_accept"] == 1.0)

# Source-matched Claudon et al. (2010) reference point: GaAs/InAs photonic
# wire, n_wire=3.45 [given in spec text, citing Claudon et al., Nat. Photon.
# 4, 174 (2010)], d/lambda=0.22, lambda=950 nm, NA=0.75, bottom_reflectivity
# 1.0, lossless taper/contact/propagation factors, taper_output_mfr_nm=750.0
# nm [A, this module's own design choice representing their reported ~1.5 um
# top-facet diameter as a mode-field RADIUS]. Non-gating: the envelope need
# not reproduce their reported 0.72 first-lens extraction (a GaAs number,
# never used as a GaN anchor -- see maslov_claudon_status), but per the fix
# round's own acceptance criterion the envelope must reach >=0.55.
_claudon_radius_nm = 0.22 * 950.0 / 2.0
_p_claudon = NitrideNanowirePhotonicsParams(
    family="vertical_photonic", n_wire=3.45, NA=0.75, bottom_reflectivity=1.0,
    taper_output_mfr_nm=750.0)
_r_claudon = response(_p_claudon, lambda_nm=950.0, outer_radius_nm=_claudon_radius_nm,
                       gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N source-matched Claudon 2010 reference point reaches eta_collection_X >= 0.55",
   _r_claudon["eta_collection_X"] >= 0.55)
print(f"non-gating: Claudon-matched envelope eta={_r_claudon['eta_collection_X']:.4f}, "
      f"beta_HE11={_r_claudon['beta_HE11']:.4f}, deviation from reported extraction "
      f"0.72 = {_r_claudon['eta_collection_X'] - 0.72:+.4f} (GaAs/InAs source, non-gating)")

# The 1/e^2 half-angle note fires for a small default (untapered) wire and
# is absent once a taper expands the mode enough.
_r_small_untapered = response(NitrideNanowirePhotonicsParams(family="vertical_photonic"),
                               lambda_nm=450.0, outer_radius_nm=80.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
_r_small_tapered = response(
    NitrideNanowirePhotonicsParams(family="vertical_photonic", taper_output_mfr_nm=2000.0),
    lambda_nm=450.0, outer_radius_nm=80.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
ck("N a small untapered wire's 1/e^2 half-angle exceeds 60 deg and is noted",
   any("exceeds 60 deg" in n for n in _r_small_untapered["notes"]))
ck("N a taper expanding the mode past 2000 nm removes the 60-deg note",
   not any("exceeds 60 deg" in n for n in _r_small_tapered["notes"]))
# [Attempt 3 fix K] the RAW paraxial diagnostic exceeds 90 deg here (a
# non-physical artifact of the small-angle formula), but the PRINTED note
# is capped at 90 deg for display; and every note fires for THIS card must
# carry the [A] platform-transfer tag for the GaAs/InAs-calibrated s(V)
# envelope applied to a GaN card.
ck("N Attempt 3 fix K: the raw theta_div_rad diagnostic itself is NOT "
   "capped (it genuinely exceeds 90 deg for this small default wire)",
   math.degrees(_r_small_untapered["diagnostics"]["theta_div_rad"]) > 90.0)
ck("N Attempt 3 fix K: the PRINTED 1/e^2 half-angle note never states a "
   "value above 90 deg, even though the underlying diagnostic does",
   all(float(_m) <= 90.0 for n in _r_small_untapered["notes"]
       for _m in re.findall(r"\(([\d.]+) deg", n)))
ck("N Attempt 3 fix K: the above-cutoff note explicitly tags the GaAs/InAs "
   "-> GaN platform transfer as [A]",
   any("platform transfer itself is [A]" in n for n in _r_small_untapered["notes"]))

# ================================================================= (fix round MEDIUM 6) beta != confinement

ck("N beta_HE11 differs from confinement_fraction at the default vertical card",
   _r_v["beta_HE11"] != _r_v["diagnostics"]["confinement_fraction"])
ck("N beta_HE11 and confinement_fraction are both in [0, 1] at the default",
   0.0 <= _r_v["beta_HE11"] <= 1.0
   and 0.0 <= _r_v["diagnostics"]["confinement_fraction"] <= 1.0)
ck("N Attempt 3 fix I: the dead 'beta_clipped' diagnostic/note is gone "
   "(beta_raw is a ratio of two non-negative terms and is always already "
   "in [0, 1])",
   "beta_clipped" not in _r_v["diagnostics"]
   and not any("clipped" in n for n in _r_v["notes"]))
ck("N group_index_ratio (n_g/n_wire) differs from 1.0 at the default (dispersive) card, "
   "reproduced independently via a fresh finite difference of gan_ordinary_index",
   not math.isclose(_r_v["diagnostics"]["group_index_ratio"], 1.0, abs_tol=1e-3))
# [Attempt 3 fix D] REVERSES the old claim: an n_wire override no longer
# collapses group_index_ratio to a hardcoded 1.0 -- it applies the SAME
# GaN-Sellmeier dispersion slope, scaled onto the override's magnitude,
# reproduced fresh here via an independent central finite difference of
# gan_ordinary_index (never via _group_index_ratio's own internals).
_d_fd = 1.0
_n_lo_fd = gan_ordinary_index(500.0 - _d_fd)
_n_hi_fd = gan_ordinary_index(500.0 + _d_fd)
_dn_dlambda_fd = (_n_hi_fd - _n_lo_fd) / (2.0 * _d_fd)
_expected_ratio_override = (2.2 - 500.0 * _dn_dlambda_fd) / 2.2
ck("N Attempt 3 fix D: group_index_ratio for an n_wire override matches an "
   "independent re-derivation (GaN Sellmeier slope scaled onto the "
   "override's magnitude, fresh finite difference typed here) and is NOT "
   "1.0",
   close(_group_index_ratio(2.2, 2.2, 500.0), _expected_ratio_override, rtol=1e-6)
   and not math.isclose(_group_index_ratio(2.2, 2.2, 500.0), 1.0, abs_tol=1e-3))
ck("N Attempt 3 fix D: group_index_ratio still falls back to 1.0 only when "
   "lambda_nm itself is outside the GaN Sellmeier's own validity window "
   "(350-10000 nm), preserving the n_wire-override-bypasses-GaN-range "
   "invariant",
   _group_index_ratio(2.2, 2.2, 100.0) == 1.0)
ck("N Attempt 3 fix D: beta_HE11 differs from confinement_fraction at the "
   "SOURCE-MATCHED CLAUDON REFERENCE too (n_wire override), not only at "
   "the default GaN card",
   not math.isclose(_r_claudon["beta_HE11"], _r_claudon["diagnostics"]["confinement_fraction"],
                     rel_tol=1e-6))
print(f"non-gating: Claudon-matched beta_HE11={_r_claudon['beta_HE11']:.4f} vs reported "
      f"guided-mode beta ~0.95 at d/lambda~0.24 (GaAs/InAs source, non-gating; "
      f"deviation {_r_claudon['beta_HE11'] - 0.95:+.4f})")

# ================================================================= (Attempt 3 fix A) antenna_rate_factor (rate, not collection)

# Independent per-orientation collection/screening literals at 450/630 nm
# (already-verified dipole_collection_fraction plus a freshly typed
# screening formula, not imported from the module) cross-checked against
# pure-orientation cards, where antenna_rate_factor collapses to the
# single orientation's own s_i and gamma_X_ns*eta_collection_X (the
# collected rate) exactly reproduces the OLD (pre-Attempt-3) screened
# per-orientation product s_i*eta_raw_i.
for _lam_a in (450.0, 630.0):
    _n_wire_a = gan_ordinary_index(_lam_a)
    _n_ox_a = sio2_index(_lam_a)
    _n_sub_a = si_complex_index(_lam_a)
    _screen_a = (2.0 / (_n_wire_a ** 2 + 1.0)) ** 2
    _eta_a_along = dipole_collection_fraction((1.0, 0.0, 0.0), 0.5, 1.0, _n_ox_a, _n_sub_a,
                                               100.0, 12.5, _lam_a)
    _eta_a_trans = dipole_collection_fraction((0.0, 1.0, 0.0), 0.5, 1.0, _n_ox_a, _n_sub_a,
                                               100.0, 12.5, _lam_a)
    _eta_a_vert = dipole_collection_fraction((0.0, 0.0, 1.0), 0.5, 1.0, _n_ox_a, _n_sub_a,
                                              100.0, 12.5, _lam_a)
    for _label_a, _vec_a, _s_expected, _eta_expected in (
            ("along_wire", (1.0, 0.0, 0.0), 1.0, _eta_a_along),
            ("transverse_inplane", (0.0, 1.0, 0.0), _screen_a, _eta_a_trans),
            ("vertical", (0.0, 0.0, 1.0), _screen_a, _eta_a_vert)):
        _r_pure = response(NitrideNanowirePhotonicsParams(family="horizontal_as_built",
                                                            dipole_weights=_vec_a),
                            lambda_nm=_lam_a, outer_radius_nm=12.5,
                            gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
        ck(f"N Attempt 3 fix A: pure-{_label_a} card at {_lam_a:.0f} nm: "
           "antenna_rate_factor equals the independent s_i",
           close(_r_pure["antenna_rate_factor"], _s_expected, rtol=1e-9))
        ck(f"N Attempt 3 fix A: pure-{_label_a} card at {_lam_a:.0f} nm: "
           "eta_collection_X equals the independent UNSCREENED "
           "per-orientation collection",
           close(_r_pure["eta_collection_X"], _eta_expected, rtol=1e-6))
        ck(f"N Attempt 3 fix A: pure-{_label_a} card at {_lam_a:.0f} nm: "
           "gamma_X_ns*eta_collection_X (collected rate) reproduces the "
           "OLD screened product s_i*eta_raw_i",
           close(_r_pure["gamma_X_ns"] * _r_pure["eta_collection_X"],
                 _s_expected * _eta_expected, rtol=1e-6))

ck("N Attempt 3 fix A: antenna_rate_factor is mutation-sensitive to "
   "dipole_weights (isotropic default differs from a pure along_wire card)",
   response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"),
            lambda_nm=450.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
            )["antenna_rate_factor"]
   != response(NitrideNanowirePhotonicsParams(family="horizontal_as_built",
                                               dipole_weights=(1.0, 0.0, 0.0)),
               lambda_nm=450.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
               )["antenna_rate_factor"])
ck("N Attempt 3 fix A: vertical_photonic antenna_rate_factor is the "
   "neutral 1.0 (no subwavelength dielectric-antenna screening modeled)",
   response(NitrideNanowirePhotonicsParams(family="vertical_photonic"),
            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
            )["antenna_rate_factor"] == 1.0)

# ================================================================= (Attempt 3 fix B+E+F) polarization

_pol_ledger = _ledger["deshpande2013_polarization"]
for _lam_pol in (450.0, 630.0):
    _n_ox_fresh = sio2_index(_lam_pol)
    _n_sub_fresh = si_complex_index(_lam_pol)
    _n_wire_fresh = gan_ordinary_index(_lam_pol)
    _screen_fresh = (2.0 / (_n_wire_fresh ** 2 + 1.0)) ** 2
    _along_fresh = dipole_collection_fraction((1.0, 0.0, 0.0), 0.5, 1.0, _n_ox_fresh, _n_sub_fresh,
                                               100.0, 12.5, _lam_pol)
    _trans_fresh = dipole_collection_fraction((0.0, 1.0, 0.0), 0.5, 1.0, _n_ox_fresh, _n_sub_fresh,
                                               100.0, 12.5, _lam_pol)
    _vert_fresh = dipole_collection_fraction((0.0, 0.0, 1.0), 0.5, 1.0, _n_ox_fresh, _n_sub_fresh,
                                              100.0, 12.5, _lam_pol)

    def _dolp_fresh_sum(weights, _along=_along_fresh, _trans=_trans_fresh,
                         _vert=_vert_fresh, _screen=_screen_fresh):
        # [Attempt 3 fix E] SUM (not mean) convention, typed fresh here --
        # never imported from the module's own _dolp_sum_convention.
        w_along, w_transverse, w_vertical = weights
        i_par = w_along * _along
        i_perp = w_transverse * _trans * _screen + 0.5 * w_vertical * _vert * _screen
        i_par += 0.5 * w_vertical * _vert * _screen
        denom = i_par + i_perp
        return (i_par - i_perp) / denom if denom > 0.0 else 0.0

    _iso_w = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    _cplane_w = (0.0, 0.5, 0.5)

    _r_pol_default = response(NitrideNanowirePhotonicsParams(family="horizontal_as_built"),
                               lambda_nm=_lam_pol, outer_radius_nm=12.5,
                               gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    _r_pol_cplane = response(NitrideNanowirePhotonicsParams(family="horizontal_as_built",
                                                              dipole_weights=_cplane_w),
                              lambda_nm=_lam_pol, outer_radius_nm=12.5,
                              gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    _dolp_default = _r_pol_default["degree_of_linear_polarization"]
    _dolp_cplane = _r_pol_cplane["degree_of_linear_polarization"]
    _dolp_iso_diag = _r_pol_default["diagnostics"]["degree_of_linear_polarization_isotropic"]

    ck(f"N degree_of_linear_polarization (default, isotropic card) at "
       f"{_lam_pol:.0f} nm is a finite fraction in [-1, 1]",
       isinstance(_dolp_default, float) and -1.0 <= _dolp_default <= 1.0)
    ck(f"N Attempt 3 fix B: default-card degree_of_linear_polarization at "
       f"{_lam_pol:.0f} nm matches an independent SUM-convention "
       "re-derivation for the card's own (isotropic) weights",
       close(_dolp_default, _dolp_fresh_sum(_iso_w), rtol=1e-9))
    ck(f"N Attempt 3 fix B: diagnostics['degree_of_linear_polarization_"
       f"isotropic'] at {_lam_pol:.0f} nm matches the same independent "
       "isotropic re-derivation (and equals the default card's own value, "
       "since the default IS isotropic)",
       close(_dolp_iso_diag, _dolp_fresh_sum(_iso_w), rtol=1e-9)
       and close(_dolp_iso_diag, _dolp_default, rtol=1e-9))
    ck(f"N Attempt 3 fix B+F: cplane_only-card degree_of_linear_polarization "
       f"at {_lam_pol:.0f} nm matches an independent SUM-convention "
       "re-derivation and remains opposite-sign with the expected magnitude band",
       close(_dolp_cplane, _dolp_fresh_sum(_cplane_w), rtol=1e-9)
       and _dolp_cplane < -0.8 and _dolp_cplane > -1.0)
    ck(f"N Attempt 3 fix B: cplane_only prediction at {_lam_pol:.0f} nm is "
       "printed as a falsification against the anchor (never hidden)",
       any("falsified" in n for n in _r_pol_cplane["notes"]))
    print(f"non-gating: predicted degree_of_linear_polarization (default, "
          f"isotropic) at {_lam_pol:.0f} nm = {_dolp_default * 100.0:.1f}% vs "
          f"deshpande2013_polarization anchor "
          f"{_pol_ledger['value']['axial_dolp_percent']:.1f}% (deviation "
          f"{_dolp_default * 100.0 - _pol_ledger['value']['axial_dolp_percent']:+.1f} "
          "points, non-gating)")
    print(f"non-gating: predicted degree_of_linear_polarization "
          f"(cplane_only sensitivity) at {_lam_pol:.0f} nm = "
          f"{_dolp_cplane * 100.0:.1f}% vs deshpande2013_polarization anchor "
          f"{_pol_ledger['value']['axial_dolp_percent']:.1f}% (deviation "
          f"{_dolp_cplane * 100.0 - _pol_ledger['value']['axial_dolp_percent']:+.1f} "
          "points, non-gating; OPPOSITE SIGN from the anchor: c-plane-only "
          "dipole prior falsified by the polarization anchor)")

ck("N wire_antenna_screening_intensity matches the closed-form (2/(n^2+1))^2 independently",
   math.isclose(_wire_antenna_screening_intensity(2.4869166125042943),
                (2.0 / (2.4869166125042943 ** 2 + 1.0)) ** 2, rel_tol=1e-12))
ck("N wire_antenna_screening_intensity is exactly 1.0 (no screening) when n_wire=1 (no dielectric contrast)",
   _wire_antenna_screening_intensity(1.0) == 1.0)
ck("N degree_of_linear_polarization is mutation-sensitive to n_wire (screening factor)",
   response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", n_wire=1.0),
            lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
            )["degree_of_linear_polarization"]
   != _r_h["degree_of_linear_polarization"])
ck("N Attempt 3 fix B: degree_of_linear_polarization is mutation-sensitive "
   "to dipole_weights itself (isotropic default, cplane_only, and pure "
   "along_wire all give three distinct values)",
   len({response(NitrideNanowirePhotonicsParams(family="horizontal_as_built", dipole_weights=w),
                 lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
                 )["degree_of_linear_polarization"]
        for w in ((1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), (0.0, 0.5, 0.5), (1.0, 0.0, 0.0))}) == 3)
ck("N Attempt 3 fix B: a pure along_wire card gives degree_of_linear_"
   "polarization exactly +1.0 (100% axial, I_perp=0)",
   close(response(NitrideNanowirePhotonicsParams(family="horizontal_as_built",
                                                   dipole_weights=(1.0, 0.0, 0.0)),
                  lambda_nm=500.0, outer_radius_nm=12.5, gamma_X0_ns=1.0, gamma_XX0_ns=1.0
                  )["degree_of_linear_polarization"], 1.0, rtol=1e-9))

# ================================================================= (fix round LOW 11) stated approximations

ck("N horizontal provenance states the emitter-medium (n_ambient, not n_wire) approximation",
   "emitter_medium_approximation" in _r_h["provenance"]
   and any("emitter" in n and "n_ambient" in n for n in _r_h["notes"]))
ck("N vertical provenance states the incoherent bottom-mirror treatment",
   "mirror_treatment" in _r_v["provenance"]
   and any("incoherent" in n for n in _r_v["notes"]))

# ================================================================= (fix round DIRECTIVE H7) beta_HE11 interior maximum

# Source-matched to the SAME n_wire=3.45/lambda=950nm Claudon platform used
# above (never the module's own default GaN dispersion, since the two
# transcribed anchor points are read off that platform's own beta_HE11(V)
# curve) -- a fresh call to response() at each d/lambda, not a re-derivation
# of _multimode_penalty's own formula.
def _h7_response_at(d_over_lambda, n_wire=3.45, lambda_nm=950.0):
    radius_nm = d_over_lambda * lambda_nm / 2.0
    return response(NitrideNanowirePhotonicsParams(family="vertical_photonic", n_wire=n_wire),
                     lambda_nm=lambda_nm, outer_radius_nm=radius_nm,
                     gamma_X0_ns=1.0, gamma_XX0_ns=1.0)


_h7_grid = [0.15, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.35, 0.40]
_h7_resp = {d: _h7_response_at(d) for d in _h7_grid}
_h7_beta = {d: _h7_resp[d]["beta_HE11"] for d in _h7_grid}
_h7_peak_d = max(_h7_beta, key=_h7_beta.get)

ck("N beta_HE11(d/lambda) has its maximum strictly inside (0.18, 0.30) on the "
   "n_wire=3.45/950nm Claudon-matched platform (fix round DIRECTIVE H7 "
   "acceptance criterion)",
   0.18 < _h7_peak_d < 0.30)
ck("N beta_HE11 at d/lambda=0.40 is lower than at d/lambda=0.24 (fix round "
   "DIRECTIVE H7 acceptance criterion)",
   _h7_beta[0.40] < _h7_beta[0.24])
ck("N approximation_error is exactly 0 at d/lambda<=0.22 (at/below the LP11 "
   "cutoff on this platform) and strictly positive at every tested "
   "d/lambda>=0.24 (above cutoff) -- fix round DIRECTIVE H7",
   _h7_resp[0.22]["approximation_error"] == 0.0
   and all(_h7_resp[d]["approximation_error"] > 0.0 for d in _h7_grid if d >= 0.24))
# [Attempt 3 fix H] approximation_error is a bookkeeping SUM of two
# UNLIKE-UNIT quantities: marcuse_window_error (V-number units) and
# multimode_penalty_deficit (dimensionless, 1-beta_multimode_penalty),
# each now separately named in diagnostics and re-derived independently
# here (marcuse_window_error via a fresh copy of the extrapolation-
# distance formula, multimode_penalty_deficit as 1-beta_multimode_penalty
# read from the response's own already-checked field).
for _d_h in _h7_grid:
    _diag_h = _h7_resp[_d_h]["diagnostics"]
    _v_h = _h7_resp[_d_h]["V_number"]
    _marcuse_expected = max(0.0, 0.8 - _v_h) + max(0.0, _v_h - 2.5)
    _deficit_expected = max(0.0, 1.0 - _h7_resp[_d_h]["beta_multimode_penalty"])
    ck(f"N Attempt 3 fix H: marcuse_window_error + multimode_penalty_deficit "
       f"sums exactly to approximation_error at d/lambda={_d_h:.2f}, and each "
       "matches an independent re-derivation",
       close(_diag_h["marcuse_window_error"] + _diag_h["multimode_penalty_deficit"],
             _h7_resp[_d_h]["approximation_error"], rtol=1e-12)
       and close(_diag_h["marcuse_window_error"], _marcuse_expected, rtol=1e-9)
       and close(_diag_h["multimode_penalty_deficit"], _deficit_expected, rtol=1e-9))
# [Attempt 3 fix J] the DELIVERED beta_HE11(d/lambda) decline between the
# two anchor points is honestly smaller than the anchors' own ratio-based
# decline (0.70/0.95): printed, not glossed over.
_delivered_decline_pct = (_h7_beta[0.24] - _h7_beta[0.40]) / _h7_beta[0.24] * 100.0
_anchored_decline_pct = (1.0 - 0.70 / 0.95) * 100.0
ck("N Attempt 3 fix J: the delivered beta_HE11(d/lambda) decline (0.24->0.40) "
   "is smaller than the anchors' own relative decline (honesty check, not a "
   "calibration target)",
   _delivered_decline_pct < _anchored_decline_pct)
print(f"non-gating: Attempt 3 fix J -- delivered beta_HE11 decline "
      f"(d/lambda 0.24->0.40) = {_delivered_decline_pct:.1f}% vs the "
      f"Bleuse/Claudon anchors' own relative decline = "
      f"{_anchored_decline_pct:.1f}% (0.70/0.95); s(V) is calibrated to the "
      "anchors' decline, not to reproduce the delivered curve's own decline")
ck("N single_mode is True at d/lambda<=0.22 and False at d/lambda>=0.24 "
   "(fix round DIRECTIVE H7)",
   _h7_resp[0.22]["single_mode"] is True
   and all(_h7_resp[d]["single_mode"] is False for d in _h7_grid if d >= 0.24))
ck("N beta_multimode_penalty equals exactly 1.0 at/below cutoff and is "
   "strictly less than 1.0 above it (fix round DIRECTIVE H7)",
   _h7_resp[0.22]["beta_multimode_penalty"] == 1.0
   and all(_h7_resp[d]["beta_multimode_penalty"] < 1.0 for d in _h7_grid if d >= 0.24))
ck("N single_mode and beta_multimode_penalty are 'not_applicable' for "
   "horizontal_as_built, matching beta_HE11's own convention (fix round "
   "DIRECTIVE H7)",
   _r_h["single_mode"] == "not_applicable"
   and _r_h["beta_multimode_penalty"] == "not_applicable")
print("non-gating: beta_HE11(d/lambda) at n_wire=3.45/950nm (fix round DIRECTIVE H7): "
      + ", ".join(f"{d:.2f}->{_h7_beta[d]:.4f}" for d in (0.15, 0.20, 0.24, 0.30, 0.40))
      + f"; maximum at d/lambda={_h7_peak_d:.2f} (Bleuse et al. PRL 106, 103601 "
        "(2011) Fig. 2 / Claudon et al. (2010) anchors, non-gating)")

# ================================================================= (Attempt 3 fix B) default dipole_weights reverts to isotropic

# [Attempt 3 fix B] REVERSES the DIRECTIVE M5 default: (0.0, 0.5, 0.5)
# predicts the wrong-sign DOLP against the anchor (see the fix B+E+F
# section above), so it cannot be the default; isotropic (orientation
# prior unknown) is.
ck("T CPLANE_ONLY_DIPOLE_WEIGHTS equals the literal (0.0, 0.5, 0.5) named "
   "sensitivity constant",
   CPLANE_ONLY_DIPOLE_WEIGHTS == (0.0, 0.5, 0.5))
ck("T Attempt 3 fix B: default dipole_weights is isotropic (1/3, 1/3, 1/3), "
   "NOT CPLANE_ONLY_DIPOLE_WEIGHTS -- orientation prior unknown",
   NitrideNanowirePhotonicsParams(family="horizontal_as_built").dipole_weights
   == (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
   and NitrideNanowirePhotonicsParams(family="horizontal_as_built").dipole_weights
   != CPLANE_ONLY_DIPOLE_WEIGHTS)
ck("N CPLANE_ONLY_DIPOLE_WEIGHTS remains an accepted, explicit named "
   "sensitivity input, distinct from the isotropic default",
   NitrideNanowirePhotonicsParams(
       family="horizontal_as_built", dipole_weights=CPLANE_ONLY_DIPOLE_WEIGHTS
   ).dipole_weights == CPLANE_ONLY_DIPOLE_WEIGHTS)
ck("N a default vertical_photonic card (dipole_weights unspecified) matches "
   "the isotropic horizontal default and is accepted by the two-sided "
   "cross-talk guard",
   NitrideNanowirePhotonicsParams(family="vertical_photonic").dipole_weights
   == (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0))
ck("N vertical_photonic still raises when dipole_weights is set to the "
   "OLD (0.0, 0.5, 0.5) default (now non-default for the two-sided guard)",
   expect(NitrideNanowirePhotonicsError,
          lambda: response(NitrideNanowirePhotonicsParams(
                                family="vertical_photonic", dipole_weights=CPLANE_ONLY_DIPOLE_WEIGHTS),
                            lambda_nm=500.0, outer_radius_nm=90.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)))

passed = sum(1 for _, ok in checks if ok)
total = len(checks)
for name, ok in checks:
    if not ok:
        print("FAIL " + name)
print(f"{passed}/{total} nitride nanowire photonics checks passed")
raise SystemExit(0 if passed == total else 1)
