"""Source transcription, numerical and non-gating comparison checks for the
GaN/AlGaN tunnel injector (piece 6, fsim_core/nitride_nanowire_injector.py).

Every expected number here is either transcribed from a cited published
source (never produced by the module under test) or independently derived
by a SEPARATE method written in this file (a fresh transcendental-equation
root-solve, a fresh analytic slab-transmission formula, or a fresh
brute-force quadrature) -- never by calling the production transmission()/
injector_feasibility() routines to manufacture their own "expected" value.
"""
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import simpson
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import nitride_materials as NM
from fsim_core.nitride_nanowire_injector import (
    NitrideNanowireInjectorParams, transmission, reflection, injector_feasibility,
    M0_KG, EV_J, HBAR_JS, KB_EV, _GROWTH_STEP_NM_DEFAULT,
)

c = []


def check(label, value):
    c.append(bool(value))
    if not value:
        print("FAIL", label)


# ---------------------------------------------------------------------
# Independent carrier-statistics helpers (audit M1, 2026-09-23), written
# fresh here and NOT imported from the module under test: the normalized
# complete Fermi-Dirac integral F_1/2(eta) = (2/sqrt(pi)) int_0^inf
# sqrt(x)/(1+exp(x-eta)) dx by scipy.integrate.quad [DR; Blakemore,
# Solid-State Electron. 25, 1067 (1982)], N_c = 2 (2 pi m k T / h^2)^(3/2)
# [DR; Sze & Ng, "Physics of Semiconductor Devices" 3rd ed. (2007) Ch. 1],
# the ellipsoidal DOS mass (me_xy^2 me_z)^(1/3) from nitride_materials'
# [V] Rinke PRB 77, 075202 (2008) GaN masses, and the T=0 degenerate
# electron-gas level (hbar^2/2m)(3 pi^2 n)^(2/3) [DR; Ashcroft & Mermin
# (1976) Ch. 2 Eq. 2.33].
from scipy.integrate import quad as _quad_indep  # noqa: E402

_H_JS_INDEP = 6.62607015e-34        # [V] CODATA 2018, exact
_KB_J_INDEP = 1.380649e-23          # [V] CODATA 2018, exact


def _f_half_indep(eta):
    val = _quad_indep(lambda x: math.sqrt(x) / (1.0 + math.exp(min(x - eta, 700.0))),
                      0.0, max(eta, 0.0) + 80.0, points=[eta] if eta > 0 else None,
                      limit=500, epsabs=0.0, epsrel=1e-12)[0]
    return 2.0 / math.sqrt(math.pi) * val


def _nc_indep_m3(m_ratio, T_K):
    return 2.0 * (2.0 * math.pi * m_ratio * M0_KG * _KB_J_INDEP * T_K / _H_JS_INDEP ** 2) ** 1.5


def _mdos_e_indep(me_z=None):
    """Bulk n-GaN emitter DOS mass from independently typed Rinke et al.,
    PRB 77, 075202 (2008) Table V literals (m_e-perp 0.209 in the c plane,
    m_e-par 0.186 along c): (0.209^2 * 0.186)^(1/3) = 0.2010.  The argument
    is accepted and IGNORED on purpose (strain-mass audit, 2026-09-23): the
    emitter mass must not follow a me_well (well tunnelling mass) override."""
    return (0.209 ** 2 * 0.186) ** (1.0 / 3.0)


def _mu_t0_indep_eV(n_m3, m_ratio):
    k_f = (3.0 * math.pi ** 2 * n_m3) ** (1.0 / 3.0)
    return HBAR_JS ** 2 * k_f ** 2 / (2.0 * m_ratio * M0_KG) / EV_J


def _n_cm3_for_mu_fd_indep(mu_eV, me_z, T_K):
    """Forward Fermi-Dirac density (no root solve): n = N_c F_1/2(mu/kT)."""
    kT_eV = _KB_J_INDEP * T_K / EV_J
    return _nc_indep_m3(_mdos_e_indep(me_z), T_K) * _f_half_indep(mu_eV / kT_eV) / 1e6


def _mu_fd_indep_eV(n_cm3, me_z, T_K):
    kT_eV = _KB_J_INDEP * T_K / EV_J
    ratio = n_cm3 * 1e6 / _nc_indep_m3(_mdos_e_indep(me_z), T_K)
    eta = brentq(lambda e: _f_half_indep(e) - ratio, -60.0, 400.0, xtol=1e-13, rtol=1e-14)
    return eta * kT_eV


# =====================================================================
# 1. Published structure/transmission benchmark (source-matched conditions)
# =====================================================================
# Encomendero et al., "Broken Symmetry Effects due to Polarization on
# Resonant Tunneling Transport in Double-Barrier Nitride Heterostructures",
# arXiv:2303.08352 (2023), Sec. II [V, full text fetched 2026-09-14]: the
# tb=1.5 nm device has AlN barriers 1.5 nm thick, a 3-nm-wide GaN quantum
# well, and degenerately doped (~2e19 cm-3 Si) n-GaN contacts; the paper
# reports peak current density 2.55e4 A/cm2 at 5.66 V for this device.  This
# module does NOT attempt to reproduce that onset voltage or current density
# (a Poisson-coupled, doping- and bias-dependent band-bending calculation is
# out of scope for a flat-barrier transmission engine) -- consistent with
# the spec's "if no source-matched measured rate exists, say so" escape
# hatch, those two numbers are recorded here for provenance only and are
# NOT gated on.  What IS independently checkable at these exact published
# materials/thicknesses is the flat-band RESONANCE ENERGY of the structure,
# against a symmetric mass-mismatched finite-square-well transcendental
# equation solved fresh below (Griffiths, "Introduction to Quantum
# Mechanics", symmetric well, BenDaniel-Duke matching condition) -- a
# textbook-standard, independently-coded root-solve, not a call into the
# module's own resonance finder.

benchmark = NitrideNanowireInjectorParams(
    al_fraction=1.0, electron_topology="double_barrier",
    electron_barrier_thickness_nm=1.5, electron_well_width_nm=3.0,
    n_cm3=2.0e19, p_cm3=2.0e19, include_polarization=False,
)
_ENCOMENDERO_PEAK_V = 5.66          # [V] non-gating comparison only
_ENCOMENDERO_PEAK_A_CM2 = 2.55e4    # [V] non-gating comparison only
_ENCOMENDERO_BARRIER_NM = 1.5
_ENCOMENDERO_WELL_NM = 3.0

# MEDIUM 8 fix: cross-check this inline transcription against the actual
# evidence ledger anchor (verify/data/nitride_nanowire_anchors.yaml), not
# just isfinite() on the two literals defined two lines above (a tautology
# that could never fail regardless of what the ledger says or even whether
# it exists).
_ANCHORS_PATH = Path(__file__).resolve().parent / "data" / "nitride_nanowire_anchors.yaml"
_anchors_doc = yaml.safe_load(_ANCHORS_PATH.read_text(encoding="utf-8")) or {}
_anchors = _anchors_doc.get("anchors") or {}
_encomendero_anchor = _anchors.get("encomendero2023_resonant_tunneling")
check("benchmark_numbers_match_the_evidence_ledger_anchor",
      _encomendero_anchor is not None
      and _encomendero_anchor.get("value", {}).get("peak_voltage_V") == _ENCOMENDERO_PEAK_V
      and _encomendero_anchor.get("value", {}).get("peak_current_density_A_cm2") == _ENCOMENDERO_PEAK_A_CM2
      and _encomendero_anchor.get("value", {}).get("barrier_nm") == _ENCOMENDERO_BARRIER_NM
      and _encomendero_anchor.get("value", {}).get("well_nm") == _ENCOMENDERO_WELL_NM)

# Along-c (tunnelling-axis) electron masses, Rinke PRB 2008 Table V row
# m_e-par [V]; re-pinned 2026-09-23 (strain-mass audit) from 0.209 / 0.329,
# which were the in-plane (m_e-perp) values under the old swapped axes.
m_w = 0.186   # [V] GaN conduction mass along c, Rinke PRB 2008 (via nitride_materials)
m_b = 0.322   # [V] AlN conduction mass along c, Rinke PRB 2008 (via nitride_materials)


def _independent_ground_state_eV(m_w, m_b, V0_eV, L_m):
    """Freshly-coded symmetric finite-well transcendental root-solve
    (even parity), independent of the production module's own resonance
    finder: k*tan(k L/2) = (m_w/m_b)*kappa."""
    def k_of(E):
        return math.sqrt(2.0 * m_w * M0_KG * E * EV_J) / HBAR_JS

    def kappa_of(E):
        return math.sqrt(2.0 * m_b * M0_KG * (V0_eV - E) * EV_J) / HBAR_JS

    def g(E):
        k = k_of(E)
        return k * math.tan(k * L_m / 2.0) - (m_w / m_b) * kappa_of(E)

    Es = np.linspace(V0_eV * 1e-6, V0_eV * (1 - 1e-6), 20000)
    prev = g(Es[0])
    for j in range(1, len(Es)):
        cur = g(Es[j])
        if math.isfinite(prev) and math.isfinite(cur) and prev < 0.0 <= cur:
            return brentq(g, Es[j - 1], Es[j])
        prev = cur
    raise RuntimeError("no root found")


# Barrier height at al_fraction=1 (pure AlN) is read back from the module's
# own resolved path only to fix V0 for the INDEPENDENT solver below (using
# the same [V] material inputs, not the module's transmission physics).
from fsim_core.nitride_nanowire_injector import _resolve_path, _find_resonance, _resonance_width_eV  # noqa: E402

_bench_path = _resolve_path(benchmark, "electron")
V0_bench = _bench_path.barrier_height_eV
L_bench = _bench_path.well_nm * 1e-9

E_ground_independent = _independent_ground_state_eV(m_w, m_b, V0_bench, L_bench)

# Locate the same resonance from the PRODUCTION transmission() API by a
# plain local scan (mimicking what any external consumer would do) --
# deliberately not importing the module's internal resonance finder.
window = np.linspace(E_ground_independent * 0.9, E_ground_independent * 1.1, 4000)
T_window = transmission(benchmark, window, carrier="electron")
i_peak = int(np.argmax(T_window))
for _ in range(6):
    lo = max(window[max(i_peak - 2, 0)], 1e-9)
    hi = window[min(i_peak + 2, len(window) - 1)]
    if hi <= lo:
        break
    window = np.linspace(lo, hi, 4000)
    T_window = transmission(benchmark, window, carrier="electron")
    i_peak = int(np.argmax(T_window))
E_ground_numeric = float(window[i_peak])
T_ground_numeric = float(T_window[i_peak])

check("benchmark_resonance_matches_independent_finite_well_eqn_2pct",
      abs(E_ground_numeric / E_ground_independent - 1.0) < 0.02)
check("benchmark_resonance_is_near_unity_transmission",
      T_ground_numeric > 0.99)
# Compare against nitride_materials directly (the actual cited [V] source
# of these numbers), not against this file's own re-declared literals --
# the previous version of this check only compared m_w/m_b to themselves
# and could never detect drift in, or a bug transferring, the underlying
# materials module.
check("benchmark_masses_match_nitride_materials_source",
      abs(m_w - NM.binary("GaN").me_z) < 1e-9 and abs(m_b - NM.binary("AlN").me_z) < 1e-9)


# =====================================================================
# 2. Transmission-engine unitarity, limits, resonance, thick-barrier WKB
# =====================================================================

rng = np.random.default_rng(20260914)

# T + R = 1 (unitarity/flux conservation) across random energies, both
# carriers, both topologies, with and without bias/field.
all_unitary = True
for carrier in ("electron", "hole"):
    for topo_kwargs in (dict(electron_topology="single_barrier", hole_topology="single_barrier"),
                        dict(electron_topology="double_barrier", hole_topology="double_barrier")):
        p = NitrideNanowireInjectorParams(**topo_kwargs)
        for E in rng.uniform(1e-3, 1.4, 12):
            for bias, field in ((0.0, 0.0), (0.3, 5.0), (-0.2, -3.0)):
                T = transmission(p, float(E), bias_V=bias, field_kVcm=field, carrier=carrier)
                R = reflection(p, float(E), bias_V=bias, field_kVcm=field, carrier=carrier)
                if abs((T + R) - 1.0) > 1e-6:
                    all_unitary = False
check("transmission_reflection_unitarity_T_plus_R_1", all_unitary)

# Zero-barrier-height, matched-mass limit: T = 1 identically (no scattering).
p_zero = NitrideNanowireInjectorParams(
    electron_topology="single_barrier", electron_barrier_thickness_nm=5.0,
    dEc_eV_override=0.0, me_barrier_override=NitrideNanowireInjectorParams().me_well,
)
T_zero = transmission(p_zero, np.linspace(0.01, 1.0, 20), carrier="electron")
check("zero_barrier_limit_T_equals_1", np.allclose(T_zero, 1.0, atol=1e-9))

# Mass-mismatch-only slab (V=0 both sides, m1 outside / m2 inside), checked
# against the standard step-index transfer-matrix slab formula, derived
# fresh here (impedance form Z=k/m): T = 1 / (1 + ((Z1^2-Z2^2)/(2 Z1 Z2))^2
# sin^2(k2 d)) -- an independent closed-form cross-check of the engine's
# BenDaniel-Duke boundary condition, not a call into the module's own code.
m1, m2, d_nm, E_mm = 0.209, 0.35, 3.0, 0.15
p_mm = NitrideNanowireInjectorParams(
    electron_topology="single_barrier", electron_barrier_thickness_nm=d_nm,
    dEc_eV_override=0.0, me_barrier_override=m2, me_well=m1,
)
T_num = transmission(p_mm, E_mm, carrier="electron")
k1 = math.sqrt(2 * m1 * M0_KG * E_mm * EV_J) / HBAR_JS
k2 = math.sqrt(2 * m2 * M0_KG * E_mm * EV_J) / HBAR_JS
Z1, Z2 = k1 / m1, k2 / m2
T_formula = 1.0 / (1.0 + ((Z1 ** 2 - Z2 ** 2) ** 2 / (4 * Z1 ** 2 * Z2 ** 2)) * math.sin(k2 * d_nm * 1e-9) ** 2)
check("mass_mismatch_matches_independent_slab_formula",
      abs(T_num / T_formula - 1.0) < 1e-6)

# Thick-barrier limit: the engine's own asymptotic decay exponent must
# converge to the standard WKB rate -2*kappa (independently computed here),
# not merely "some" exponential decay -- a bare exp(-2 kappa d) alone
# cannot represent the double-barrier resonance checked above, but it MUST
# be this engine's correct single-barrier deep-thickness asymptote.
# This analytic flat-barrier asymptote intentionally excludes the interface
# polarization staircase.  With polarization enabled the relevant exponent is
# the tilted-barrier WKB integral, rather than -2*kappa for one flat height.
p_thick = NitrideNanowireInjectorParams(electron_topology="single_barrier",
                                        include_polarization=False)
path_thick = _resolve_path(p_thick, "electron")
V0_t, m_t = path_thick.barrier_height_eV, path_thick.m_barrier
E_t = 0.1
kappa_t = math.sqrt(2 * m_t * M0_KG * (V0_t - E_t) * EV_J) / HBAR_JS
d1, d2 = 15.0, 16.0
p1 = NitrideNanowireInjectorParams(electron_topology="single_barrier",
                                   electron_barrier_thickness_nm=d1,
                                   include_polarization=False)
p2t = NitrideNanowireInjectorParams(electron_topology="single_barrier",
                                    electron_barrier_thickness_nm=d2,
                                    include_polarization=False)
T1 = transmission(p1, E_t, carrier="electron")
T2 = transmission(p2t, E_t, carrier="electron")
slope = (math.log(T2) - math.log(T1)) / ((d2 - d1) * 1e-9)
check("thick_barrier_WKB_asymptotic_slope_matches_minus_2_kappa",
      abs(slope / (-2.0 * kappa_t) - 1.0) < 1e-3)

# Double-barrier resonance: symmetric structure at resonance must reach
# T essentially 1 (textbook RTD result for identical barriers, no loss).
p_res = NitrideNanowireInjectorParams(include_polarization=False)
path_res = _resolve_path(p_res, "electron")
E_res_indep = _independent_ground_state_eV(path_res.m_well, path_res.m_barrier,
                                            path_res.barrier_height_eV, path_res.well_nm * 1e-9)
window2 = np.linspace(E_res_indep * 0.95, E_res_indep * 1.05, 4000)
T2w = transmission(p_res, window2, carrier="electron")
check("double_barrier_resonance_reaches_near_unity_transmission",
      float(np.max(T2w)) > 0.999)
check("double_barrier_resonance_location_matches_independent_finite_well_eqn",
      abs(window2[int(np.argmax(T2w))] / E_res_indep - 1.0) < 0.02)


# =====================================================================
# 2a. HIGH 4: explicit delta_Ev_GaN_AlN_eV partition knob
# =====================================================================
# Independent hand-computed barrier heights at al_fraction=0.30 for BOTH
# partitions (this module's default 0.70 eV [Martin, Yu, Waldrop, APL 68,
# 2541 (1996)] and the alternative Tsai & Bayram 0.30 eV pinned inside
# nitride_materials), from nitride_materials' own [V] Wu et al. bandgap()
# directly -- never by calling the production _default_barrier_material().
from fsim_core.nitride_nanowire_injector import _TSAI_BAYRAM_DELTA_EV_EV  # noqa: E402

_eg_gan_300 = NM.bandgap(NM.binary("GaN"), 300.0)
_eg_aln_300 = NM.bandgap(NM.binary("AlN"), 300.0)
_dEg_x030 = (_eg_aln_300 - _eg_gan_300) * 0.30

_partition_ok = True
for _delta_ev in (0.70, _TSAI_BAYRAM_DELTA_EV_EV):
    _hole_expected_eV = 0.30 * _delta_ev
    _electron_expected_eV = _dEg_x030 - _hole_expected_eV
    _p_part = NitrideNanowireInjectorParams(al_fraction=0.30, delta_Ev_GaN_AlN_eV=_delta_ev)
    _e_path_part = _resolve_path(_p_part, "electron")
    _h_path_part = _resolve_path(_p_part, "hole")
    if abs(_e_path_part.barrier_height_eV - _electron_expected_eV) > 1e-9:
        _partition_ok = False
    if abs(_h_path_part.barrier_height_eV - _hole_expected_eV) > 1e-9:
        _partition_ok = False
check("barrier_heights_match_hand_computed_partition_both_0.70_and_0.30", _partition_ok)
check("tsai_partition_constant_is_the_pinned_0.30_eV_value",
      abs(_TSAI_BAYRAM_DELTA_EV_EV - 0.30) < 1e-12)

_p_default_partition = NitrideNanowireInjectorParams()
check("default_partition_is_070_not_tsai_030",
      abs(_p_default_partition.delta_Ev_GaN_AlN_eV - 0.70) < 1e-12)

# Coherence H4 / Attempt-4 HIGH 3-4-5: independently reproduce the
# pseudomorphic Al0.30Ga0.70N polarization electrostatics from the
# published B97 binary constants -- the UNDOPED 2/4/2 nm stack sits
# between two reservoirs that pin the potential at BOTH ends, so D
# (displacement) is a single constant across every interface and the net
# potential drop across the WHOLE stack is zero: 2*d_b*E_b + d_w*E_w = 0,
# E_i = (D - P_i)/(eps0 eps_i).  Independently solved here (fresh
# arithmetic, not the production _stack_polarization_fields_eV_per_m) and
# compared against the LITERAL hand-computed numbers from the fix-round
# directive (never against the module's own algebra, LOW 17):
#   E_b = +1.4755 MV/cm, +0.2951 eV per 2 nm barrier
#   E_w = -1.4755 MV/cm, -0.5902 eV across the 4 nm well
g, a = NM.binary("GaN"), NM.binary("AlN")
x = 0.30
av = g.a_A + x * (a.a_A - g.a_A)
ep = (g.a_A - av) / av
c13 = g.C13_GPa + x * (a.C13_GPa - g.C13_GPa)
c33 = g.C33_GPa + x * (a.C33_GPa - g.C33_GPa)
e31 = g.e31_Cm2 + x * (a.e31_Cm2 - g.e31_Cm2)
e33 = g.e33_Cm2 + x * (a.e33_Cm2 - g.e33_Cm2)
p_al = (g.Psp_Cm2 + x * (a.Psp_Cm2 - g.Psp_Cm2)
        + 2 * e31 * ep + e33 * (-2 * c13 / c33 * ep))
eps_al = g.eps_r + x * (a.eps_r - g.eps_r)
p_gan, eps_gan = g.Psp_Cm2, g.eps_r
d_b, d_w = 2.0e-9, 4.0e-9
_D_closure = (2 * d_b * p_al / eps_al + d_w * p_gan / eps_gan) / (2 * d_b / eps_al + d_w / eps_gan)
E_b_hand = (_D_closure - p_al) / (NM.EPS0_SI * eps_al)      # V/m, signed
E_w_hand = (_D_closure - p_gan) / (NM.EPS0_SI * eps_gan)    # V/m, signed
_E_B_LITERAL_MV_CM = 1.4755      # [DR] hand-derived, see module docstring
_BARRIER_DROP_LITERAL_EV = 0.2951
_WELL_DROP_LITERAL_EV = -0.5902
check("zero_net_drop_closure_field_matches_hand_literal_1.4755_MVcm",
      abs(E_b_hand / 1.0e8 - _E_B_LITERAL_MV_CM) < 1e-3)
check("zero_net_drop_closure_barrier_drop_matches_hand_literal_0.2951_eV",
      abs(E_b_hand * d_b - _BARRIER_DROP_LITERAL_EV) < 1e-3)
check("zero_net_drop_closure_well_drop_matches_hand_literal_minus_0.5902_eV",
      abs(E_w_hand * d_w - _WELL_DROP_LITERAL_EV) < 1e-3)
check("zero_net_drop_closure_self_consistent", abs(2 * d_b * E_b_hand + d_w * E_w_hand) < 1e-3)

_p_pol = NitrideNanowireInjectorParams(polarity="Ga")
_p_pol_N = NitrideNanowireInjectorParams(polarity="N")
_p_flat = NitrideNanowireInjectorParams(include_polarization=False)

# The module's own reported per-barrier drop (rti_barrier_polarization_tilt_eV)
# must match the hand literal too -- via injector_feasibility's public output,
# not by importing the internal closure function.
_r_pol = injector_feasibility(_p_pol, T_K=300.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
check("module_barrier_polarization_tilt_matches_hand_literal_both_carriers",
      abs(_r_pol["rti_barrier_polarization_tilt_eV"]["electron"] - _BARRIER_DROP_LITERAL_EV) < 1e-3
      and abs(_r_pol["rti_barrier_polarization_tilt_eV"]["hole"] - _BARRIER_DROP_LITERAL_EV) < 1e-3)
check("polarization_changes_default_stack_transmission_over_10pct",
      abs(transmission(_p_pol, 0.10, carrier="electron")
          - transmission(_p_flat, 0.10, carrier="electron"))
      > 0.10 * max(transmission(_p_flat, 0.10, carrier="electron"), 1e-20))
check("mg_ionization_and_channel_count_are_reported",
      0.0 < _r_pol["rti_p_free_cm3"] < _p_pol.p_cm3
      and _r_pol["rti_reservoir_state_count_e"] == _p_pol.reservoir_state_count_e
      and _r_pol["rti_numerics_ok"] is True)

# MEDIUM 3 fix (Opus re-review of 12b39cd, 2026-09-14): the continuity
# check must probe the MODULE's OWN staircase, not a parallel hand-only
# loop that never calls production code at all -- the previous version of
# this check built _profile_points from E_b_hand/E_w_hand (correct
# NUMBERS, independently hand-derived above) but then only checked ITS OWN
# arithmetic for self-consistency, so it could never catch a break in the
# module's actual running polarization offset
# (_stack_polarization_fields_eV_per_m / _tilted_faces_eV).  Probe the
# module's own internal profile builder directly instead: _tilted_faces_eV
# returns the SAME per-segment (V_eV, entry_face_eV, exit_face_eV) tuples
# that _tilted_window_eV/_tilted_profile_max_eV/_find_resonance all consume,
# so this is exactly the profile the transmission engine itself sees.
from fsim_core.nitride_nanowire_injector import _tilted_faces_eV  # noqa: E402

_probe_path_e = _resolve_path(_p_pol, "electron")
_faces_probe, _tilt_probe = _tilted_faces_eV(_p_pol, _probe_path_e, 0.0, 0.0)
_max_jump_eV = 0.0
for _i in range(len(_faces_probe) - 1):
    _flat_step = _faces_probe[_i + 1][0] - _faces_probe[_i][0]
    _gap = (_faces_probe[_i + 1][1] - _faces_probe[_i][2]) - _flat_step
    _max_jump_eV = max(_max_jump_eV, abs(_gap))
check("polarization_potential_continuous_at_every_interface",
      _max_jump_eV < 1e-9)
# Ga-polar: exit face ABOVE entry face for the electron barrier (never
# abs()), read from the module's own faces, not a hand-rolled loop.
check("ga_polar_electron_barrier_exit_above_entry",
      _faces_probe[0][2] > _faces_probe[0][1])
# Same closure/sign applies identically to the hole path (same material
# electrostatics, see module docstring "Growth polarity and transport
# direction"): the module's own reported hole tilt must be POSITIVE too.
check("ga_polar_hole_barrier_tilt_reported_positive",
      _r_pol["rti_barrier_polarization_tilt_eV"]["hole"] > 0.0)
# N-polar reverses the sign for both carriers.
_r_pol_N = injector_feasibility(_p_pol_N, T_K=300.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
check("n_polar_reverses_barrier_polarization_tilt_sign_both_carriers",
      _r_pol_N["rti_barrier_polarization_tilt_eV"]["electron"] < 0.0
      and _r_pol_N["rti_barrier_polarization_tilt_eV"]["hole"] < 0.0)

# MEDIUM 6: independent (fresh) Mg-acceptor charge-neutrality mass-action
# solve, p^2/(Na-p) = (N_V/g) exp(-E_A/kT), NOT calling
# _hole_quasi_fermi_eV -- must match the module's rti_p_free_cm3 and the
# reviewer's independently-derived numbers (~3.10e16 / 9.55e16 cm-3 at
# 230/300 K).
def _n_v_hand(T_K, mh=1.88):
    m = mh * 9.1093837015e-31
    kT_J = KB_EV * T_K * 1.602176634e-19
    h_js = 4.135667696e-15 * 1.602176634e-19
    return 2.0 * (2.0 * math.pi * m * kT_J / (h_js ** 2)) ** 1.5


def _p_free_hand_cm3(T_K, Na_cm3=5.0e17, E_A_eV=0.170, g_deg=4.0):
    kT_eV = KB_EV * T_K
    NV = _n_v_hand(T_K)
    Na_m3 = Na_cm3 * 1e6
    K = (NV / g_deg) * math.exp(-E_A_eV / kT_eV)
    p_m3 = (-K + math.sqrt(K * K + 4.0 * K * Na_m3)) / 2.0
    return p_m3 / 1e6


_p_free_230_hand = _p_free_hand_cm3(230.0)
_p_free_300_hand = _p_free_hand_cm3(300.0)
check("mg_neutrality_p_free_230K_matches_independent_mass_action_solve",
      abs(_p_free_230_hand / 3.10e16 - 1.0) < 0.02)
check("mg_neutrality_p_free_300K_matches_independent_mass_action_solve",
      abs(_p_free_300_hand / 9.55e16 - 1.0) < 0.02)
_r_pol_230 = injector_feasibility(_p_pol, T_K=230.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
check("module_rti_p_free_cm3_matches_independent_mass_action_solve_230K",
      abs(_r_pol_230["rti_p_free_cm3"] / _p_free_230_hand - 1.0) < 1e-3)

# MEDIUM 7: rti_numerics_ok is a REAL computation, not a hardcoded True --
# a nominal (fine) stack is True; a DELIBERATELY coarse staircase (few,
# large slices) is False.
_p_numerics_fine = NitrideNanowireInjectorParams(occupancy_control_known=True, second_pair_control_known=True)
_r_numerics_fine = injector_feasibility(_p_numerics_fine, T_K=300.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
# A deliberately coarse staircase: a thick, deep barrier forced into just
# ONE flat slice (min_slices_per_segment=1, slice_length_nm larger than the
# barrier itself) drives the complex exponential in the backward recursion
# to overflow (deep evanescent decay collapsed into a single huge step),
# producing non-finite T/R that a hardcoded True would never catch.
_p_numerics_coarse = NitrideNanowireInjectorParams(
    electron_topology="single_barrier", hole_topology="single_barrier",
    electron_barrier_thickness_nm=500.0, hole_barrier_thickness_nm=500.0,
    dEc_eV_override=5.0, dEv_eV_override=-5.0,
    me_barrier_override=2.0, mh_barrier_override=2.0,
    slice_length_nm=1000.0, min_slices_per_segment=1, include_polarization=False,
    occupancy_control_known=True, second_pair_control_known=True)
import warnings as _warnings_numerics   # noqa: E402
with _warnings_numerics.catch_warnings():
    _warnings_numerics.simplefilter("ignore")
    _r_numerics_coarse = injector_feasibility(_p_numerics_coarse, T_K=300.0, rep_rate_hz=80e6,
        loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
        electron_spacing_meV=600.0, hole_spacing_meV=600.0,
        second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
check("rti_numerics_ok_true_for_fine_slicing", _r_numerics_fine["rti_numerics_ok"] is True)
check("rti_numerics_ok_false_for_deliberately_coarse_staircase",
      _r_numerics_coarse["rti_numerics_ok"] is False)

# MEDIUM 8: degeneracy is NOT dead -- doubling it (4 vs 2, channel count
# held fixed) must double the reported rate.
_p_deg2 = NitrideNanowireInjectorParams(degeneracy=2.0, reservoir_state_count_e=1.0, reservoir_state_count_h=1.0)
_p_deg4 = NitrideNanowireInjectorParams(degeneracy=4.0, reservoir_state_count_e=1.0, reservoir_state_count_h=1.0)
_r_deg2 = injector_feasibility(_p_deg2, T_K=300.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
_r_deg4 = injector_feasibility(_p_deg4, T_K=300.0, rep_rate_hz=80e6,
    loading_window_ns=0.1, electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=20.0, available_pair_rate_Hz=1e9)
check("degeneracy_4_doubles_rate_vs_degeneracy_2",
      abs(_r_deg4["rti_e_rate_Hz"] / _r_deg2["rti_e_rate_Hz"] - 2.0) < 1e-9
      and abs(_r_deg4["rti_h_rate_Hz"] / _r_deg2["rti_h_rate_Hz"] - 2.0) < 1e-9)


# =====================================================================
# 2b. HIGH 2: the resonance finder follows the field
# =====================================================================
# At each field, an independent dense (20000-point) full-window scan of
# the PRODUCTION transmission() finds the same lowest resonance
# (first local maximum scanning from the well floor) that _find_resonance
# reports -- not merely a value drawn from the same +/-10% neighborhood of
# the unbiased flat-band seed, which the review found the field can move
# the true peak well outside of.
from fsim_core.nitride_nanowire_injector import _tilted_window_eV  # noqa: E402

p_field = NitrideNanowireInjectorParams(include_polarization=False)
path_field = _resolve_path(p_field, "electron")
for field_kVcm in (0.0, 50.0, 200.0):
    well_floor_f, lower_top_f, _tilt_f = _tilted_window_eV(p_field, path_field, 0.0, field_kVcm)
    dense_grid = np.linspace(max(well_floor_f, 1e-9), lower_top_f * 0.999999, 20000)
    T_dense = transmission(p_field, dense_grid, field_kVcm=field_kVcm, carrier="electron")
    dense_peak_i = None
    # Strain-mass audit (2026-09-23): with the corrected along-c masses the
    # 200 kV/cm profile has a shallow non-resonant bump near 5 meV (T ~3e-6)
    # below the true ground resonance; a candidate counts only if it stands
    # at least 2x above BOTH flanking valleys (lowest sample between it and
    # the nearest higher sample, or the window edge) -- the same prominence
    # rule the production finder applies, re-coded here on the dense grid.
    for _j in range(1, len(T_dense) - 1):
        if T_dense[_j] > 1e-9 and T_dense[_j] >= T_dense[_j - 1] and T_dense[_j] >= T_dense[_j + 1]:
            _above_r = np.nonzero(T_dense[_j + 1:] > T_dense[_j])[0]
            _r_end = _j + 1 + (_above_r[0] if len(_above_r) else len(T_dense) - _j - 1)
            _above_l = np.nonzero(T_dense[:_j] > T_dense[_j])[0]
            _l_start = (_above_l[-1] + 1) if len(_above_l) else 0
            _valley_d = max(float(np.min(T_dense[_j:_r_end])), float(np.min(T_dense[_l_start:_j + 1])))
            if _valley_d <= 0.0 or T_dense[_j] >= 2.0 * _valley_d:
                dense_peak_i = _j
                break
    # The 20000-point coarse grid alone under-samples a resonance whose
    # FWHM (~0.02-0.04 meV here) is comparable to its own point spacing --
    # zoom in on the neighborhood it locates, using only the PUBLIC
    # transmission() API (same "mimic an external consumer" style as the
    # Encomendero benchmark scan above), so the dense-scan comparison value
    # is the actual peak, not an under-resolved sample of it.
    dense_E, dense_T = None, None
    if dense_peak_i is not None:
        lo_z = dense_grid[max(dense_peak_i - 1, 0)]
        hi_z = dense_grid[min(dense_peak_i + 1, len(dense_grid) - 1)]
        zoom = np.linspace(lo_z, hi_z, 4000)
        for _ in range(6):
            T_zoom = transmission(p_field, zoom, field_kVcm=field_kVcm, carrier="electron")
            k_peak = int(np.argmax(T_zoom))
            lo_z = zoom[max(k_peak - 2, 0)]
            hi_z = zoom[min(k_peak + 2, len(zoom) - 1)]
            if hi_z <= lo_z:
                break
            zoom = np.linspace(lo_z, hi_z, 4000)
        T_zoom = transmission(p_field, zoom, field_kVcm=field_kVcm, carrier="electron")
        k_peak = int(np.argmax(T_zoom))
        dense_E, dense_T = float(zoom[k_peak]), float(T_zoom[k_peak])
    found = _find_resonance(p_field, "electron", 0.0, field_kVcm)
    print(f"resonance vs field: {field_kVcm:.0f} kV/cm dense-scan E={dense_E*1000:.4f} "
          f"meV T={dense_T:.6f} | found E={found[0]*1000:.4f} meV T={found[1]:.6f}")
    check(f"resonance_follows_field_{int(field_kVcm)}kVcm_energy_within_1meV",
          dense_E is not None and found is not None
          and abs(found[0] - dense_E) * 1000.0 < 1.0)
    check(f"resonance_follows_field_{int(field_kVcm)}kVcm_transmission_within_1pct",
          dense_T is not None and found is not None
          and abs(found[1] - dense_T) < 0.01)


# =====================================================================
# 2c. HIGH 1 fix (Opus re-review of 12b39cd, 2026-09-14): the hole path
#     resolves real resonances once _tilted_window_eV uses the tilted
#     barrier's HIGHEST face (max(v0, v1)) as its top, not the lowest
#     (min(v0, v1)) -- the previous min() collapsed the hole window and
#     made _find_resonance bail entirely ("no hole resonance").
# =====================================================================
# Independent dense (50000-point) full-window scan of the PRODUCTION
# transmission() API for the DEFAULT designed stack's hole path (Ga
# polarity, 0.70 eV partition, al_fraction 0.30) -- same "mimic an
# external consumer" methodology as section 2b above, never the module's
# own resonance finder used to generate its own check.
p_hole_default = NitrideNanowireInjectorParams()
path_hole_default = _resolve_path(p_hole_default, "hole")
well_floor_h, lower_top_h, _tilt_h = _tilted_window_eV(p_hole_default, path_hole_default, 0.0, 0.0)
dense_grid_h = np.linspace(max(well_floor_h, 1e-9), lower_top_h * 0.999999, 50000)
T_dense_h = transmission(p_hole_default, dense_grid_h, carrier="hole")
_hole_peaks = []
for _j in range(1, len(T_dense_h) - 1):
    if (T_dense_h[_j] > 1e-12 and T_dense_h[_j] >= T_dense_h[_j - 1]
            and T_dense_h[_j] >= T_dense_h[_j + 1]):
        _hole_peaks.append((float(dense_grid_h[_j]), float(T_dense_h[_j])))
print("hole resonances at defaults (dense scan):",
      [(f"{e*1000:.3f}meV", f"T={t:.4e}") for e, t in _hole_peaks[:3]])
check("hole_path_resolves_resonances_once_window_uses_true_barrier_tops",
      len(_hole_peaks) >= 3)
_hole_e0_meV = _hole_peaks[0][0] * 1000.0 if _hole_peaks else float("nan")
_hole_T0 = _hole_peaks[0][1] if _hole_peaks else float("nan")
# Pinned literals from the reviewer's independent 200k-point dense scan
# (Opus re-review of 12b39cd, 2026-09-14, "fourth section"): 48.66 / 101.25
# / 166.05 meV.
check("hole_resonance_energies_match_reviewer_dense_scan_pinned_literals",
      len(_hole_peaks) >= 3
      and abs(_hole_peaks[0][0] * 1000.0 - 48.66) < 0.1
      and abs(_hole_peaks[1][0] * 1000.0 - 101.25) < 0.1
      and abs(_hole_peaks[2][0] * 1000.0 - 166.05) < 0.1)
check("hole_first_resonance_transmission_matches_reviewer_pinned_order_of_magnitude",
      math.isfinite(_hole_T0) and abs(_hole_T0 / 4.6553e-9 - 1.0) < 0.05)

found_hole_default = _find_resonance(p_hole_default, "hole", 0.0, 0.0)
check("production_find_resonance_matches_dense_scan_for_hole_at_defaults",
      found_hole_default is not None and math.isfinite(_hole_e0_meV)
      and abs(found_hole_default[0] * 1000.0 - _hole_e0_meV) < 0.1
      and abs(found_hole_default[1] / _hole_T0 - 1.0) < 0.05)

# The electron path is UNCHANGED by the HIGH 1 fix (its window's upper
# bound was already the correct barrier top): 268.920 meV with the old
# masses.  Re-pinned 2026-09-23 (strain-mass audit): with the along-c masses
# (well 0.186, barrier 0.2268) the ground resonance is at 291.104 meV
# (T 0.170).  The same transmission also has a shallow NON-resonant bump
# at ~7 meV (T 1.28e-7 over a 1.17e-7 valley); the production finder must
# not report it (resonance-contrast rule, _RESONANCE_MIN_CONTRAST).
found_electron_default = _find_resonance(p_hole_default, "electron", 0.0, 0.0)
check("electron_resonance_at_defaults_291.104meV_after_mass_axis_fix",
      found_electron_default is not None
      and abs(found_electron_default[0] * 1000.0 - 291.104) < 0.1)
_E_bump = np.linspace(0.001, 0.04, 400)
_T_bump = transmission(p_hole_default, _E_bump, carrier="electron")
_i_bump = int(np.argmax(_T_bump[:200]))
check("electron_finder_skips_the_7meV_background_bump_(dense_scan_shows_a_local_max_with_contrast_below_2)",
      0 < _i_bump < 199 and _T_bump[_i_bump] < 2.0 * float(np.min(_T_bump[_i_bump:]))
      and found_electron_default is not None and found_electron_default[0] > 0.1)

# MEDIUM 4 fix: a non-finite E_center must return a NaN rate, never a
# silent fallback to mu (direct unit-style probe of the private
# _forward_rate_hz, independent of any resonance being found at all).
from fsim_core.nitride_nanowire_injector import _forward_rate_hz  # noqa: E402

_rate_nan, _rate_above_nan, _warned_nan = _forward_rate_hz(
    p_hole_default, "electron", 0.0, 0.0, float("nan"), 0.01, 0.02, 0.05)
check("forward_rate_is_nan_for_unresolved_center_never_silently_mu",
      math.isnan(_rate_nan) and math.isnan(_rate_above_nan))

# LOW 8 fix: field -6469 kV/cm and bias 0.5 V verifier cases (module
# correct but previously unpinned) -- an extreme field pushes BOTH
# carriers' resonance out of the field-tilted window entirely (via
# injector_feasibility's field_kVcm), and a moderate bias alone (probed
# directly, since injector_feasibility itself always evaluates carriers at
# bias_V=0.0, see _carrier_screen) already leaves the hole path
# unresolved -- in both cases the module must report NaN, never a spurious
# zero alignment error.
p_extreme_field = NitrideNanowireInjectorParams(occupancy_control_known=True, second_pair_control_known=True)
r_extreme_field = injector_feasibility(
    p_extreme_field, T_K=230.0, rep_rate_hz=80e6, loading_window_ns=2.0,
    electron_level_eV=-0.6995, hole_level_eV=-0.4164,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=12.126, available_pair_rate_Hz=1e9,
    field_kVcm=-6469.0)
check("field_minus6469kVcm_gives_unresolved_reason_never_zero_alignment_error",
      r_extreme_field["valid"] is False
      and math.isnan(r_extreme_field["rti_alignment_error_e_meV"])
      and math.isnan(r_extreme_field["rti_alignment_error_h_meV"])
      and "electron_resonance_unresolved" in r_extreme_field["rti_failed_checks"]
      and "hole_resonance_unresolved" in r_extreme_field["rti_failed_checks"])

_found_bias_electron = _find_resonance(p_hole_default, "electron", 0.5, 0.0)
_found_bias_hole = _find_resonance(p_hole_default, "hole", 0.5, 0.0)
check("bias_0.5V_hole_resonance_unresolved_never_defaults_to_zero",
      _found_bias_hole is None and _found_bias_electron is not None)


# =====================================================================
# 3. Energy-grid refinement: rate/resonance integrals converge within 2%
# =====================================================================
# Independent brute-force Simpson quadrature of T(E) around the DESIGNED
# injector's own (narrow but not pathological) ground resonance, at
# increasing grid resolution, using only the PUBLIC transmission() API:
# a coarse grid must NOT be converged, while refining to a fine grid must
# converge to within 2% of a still-finer grid -- i.e. refinement is doing
# real work, not a no-op, and once done it is trustworthy.
E_design_res, T_design_res = _find_resonance(p_res, "electron", 0.0, 0.0)
fwhm_design = _resonance_width_eV(p_res, "electron", 0.0, 0.0, E_design_res,
                                   T_design_res, path_res.barrier_height_eV)
half_width_guess = fwhm_design if math.isfinite(fwhm_design) else max(E_design_res * 1e-4, 1e-9)
lo_i = E_design_res - 30 * half_width_guess
hi_i = E_design_res + 30 * half_width_guess
n_coarse, n_fine = 50, 4000

Ec = np.linspace(lo_i, hi_i, n_coarse)
Ef = np.linspace(lo_i, hi_i, n_fine)
Ic = simpson(transmission(p_res, Ec, carrier="electron"), x=Ec)
If = simpson(transmission(p_res, Ef, carrier="electron"), x=Ef)
Eff = np.linspace(lo_i, hi_i, 2 * n_fine)
Iff = simpson(transmission(p_res, Eff, carrier="electron"), x=Eff)

check("energy_grid_refinement_converged_within_2pct",
      If > 0 and Iff > 0 and abs(Iff / If - 1.0) < 0.02)
check("coarse_grid_is_not_converged_refinement_matters",
      abs(Ic / Iff - 1.0) > 0.05)


# =====================================================================
# 4. injector_feasibility: sensitivity, conjunction, falsification, safety
# =====================================================================

step = _GROWTH_STEP_NM_DEFAULT
thick_nom = 6 * step
hole_thick_nom = 2 * step
# MEDIUM 5/6 fix note: base_kwargs/call_kwargs are re-tuned for this round.
# With the hole reservoir's quasi-Fermi level now coming from the honest
# Mg-acceptor charge-neutrality solve (mu_h ~ -144 meV at 230 K, a dilute
# non-degenerate hole gas) instead of the old always-positive degenerate
# estimate, the hole path's capture rate is far smaller than the previous
# fixture assumed -- the hole barrier here is deliberately shallower/
# thinner than the electron barrier so the synthetic PASSING case still
# clears rate/bypass/margin under the corrected physics; single_barrier +
# include_polarization=False keeps this fixture's numbers independent of
# the polarization closure (exercised separately in section 2a).
#
# MEDIUM 5 fix (Opus re-review of 12b39cd, 2026-09-14): growth_tolerance
# is now an honest per-check REGRESSION diagnostic (a check that passes at
# nominal but fails under +/- growth_tolerance_steps perturbation), not a
# blanket AND of the perturbed screen -- the hole barrier here is only 2
# growth steps (~0.52 nm) thick, so a +/-1-step perturbation is a +/-50%
# thickness swing, and dEv_eV_override=-0.29 (the previous round's value)
# left NO headroom against that: the "+1 step" perturbed missed-load
# probability blew through the 1% ceiling even though the design's own
# nominal margin/bypass/missed-load all cleared comfortably.  A shallower
# hole barrier (-0.24 eV, still comfortably above the >=10 kT margin floor
# given hole_level_eV=0.005 eV and the 15 meV alignment uncertainty) and an
# explicit reservoir_state_count_h=2.0 (a caller-supplied transverse-
# channel count -- see the module's `reservoir_state_count_e/h` docstring;
# this synthetic fixture's own choice, not the disc-in-wire device's) give
# the hole rate enough headroom that BOTH the nominal design and its
# +/-1-growth-step perturbation clear every screen -- a genuinely robust,
# not merely round-number, synthetic passing case.
base_kwargs = dict(
    electron_topology="single_barrier", hole_topology="single_barrier",
    electron_barrier_thickness_nm=thick_nom, hole_barrier_thickness_nm=hole_thick_nom,
    me_barrier_override=0.25, mh_barrier_override=0.20,
    dEc_eV_override=0.30, dEv_eV_override=-0.24,
    occupancy_control_known=True, second_pair_control_known=True,
    growth_tolerance_steps=1.0, include_polarization=False,
    reservoir_state_count_h=2.0,
)
call_kwargs = dict(
    T_K=230.0, rep_rate_hz=80e6, loading_window_ns=12.4,
    electron_level_eV=0.02, hole_level_eV=0.005,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=600.0, available_pair_rate_Hz=1.0,
)


def run(pkw=None, ckw=None):
    p = NitrideNanowireInjectorParams(**{**base_kwargs, **(pkw or {})})
    kw = {**call_kwargs, **(ckw or {})}
    return injector_feasibility(p, **kw)


REQUIRED_KEYS = {
    "rti_feasible", "rti_status", "rti_transport_feasible", "rti_level_margin_kT",
    "rti_alignment_error_meV", "rti_linewidth_meV", "rti_rate_Hz", "rti_e_rate_Hz",
    "rti_h_rate_Hz", "rti_bypass_fraction", "rti_missed_load_probability",
    "rti_second_pair_probability", "rti_growth_feasible", "rti_failed_checks",
    "rti_evidence_status", "valid", "provenance",
}

baseline = run()
check("required_output_keys_present", REQUIRED_KEYS <= set(baseline))
check("orbital_margin_key_present_for_partial_selectivity", "rti_orbital_margin_kT" in baseline)

# --- 5. Synthetic fully-specified conditional passing case -----------------
check("synthetic_passing_case_is_feasible",
      baseline["rti_feasible"] is True and baseline["rti_status"] == "feasible")
check("synthetic_passing_case_no_failed_checks", baseline["rti_failed_checks"] == [])
check("synthetic_passing_case_all_thresholds_actually_cleared",
      baseline["rti_level_margin_kT"] >= 10.0
      and baseline["rti_bypass_fraction"] <= 0.01
      and baseline["rti_missed_load_probability"] <= 0.01
      and baseline["rti_second_pair_probability"] <= 0.01
      and baseline["rti_growth_feasible"] is True)

# --- Matched single-fault failing cases (criterion 4) -----------------------
r_margin = run(ckw=dict(electron_spacing_meV=5.0, hole_spacing_meV=5.0))
check("fault_thermal_margin_fails", r_margin["rti_feasible"] is False
      and r_margin["rti_level_margin_kT"] < 10.0)

r_window = run(ckw=dict(loading_window_ns=1e-5))
check("fault_rate_window_fails", r_window["rti_feasible"] is False
      and r_window["rti_missed_load_probability"] > 0.01
      and "pair_missed_load_window" in r_window["rti_failed_checks"])

r_bypass = run(pkw=dict(dEc_eV_override=0.02, dEv_eV_override=-0.02))
check("fault_thermionic_bypass_fails", r_bypass["rti_feasible"] is False
      and r_bypass["rti_bypass_fraction"] > 0.01)

# HIGH 1 fix: growth_tolerance_steps is a TOLERANCE (see module docstring),
# not a commensurability test -- an offset of a few growth steps from the
# baseline (itself an exact multiple of the growth step) still lands
# within tolerance of ITS OWN nearest commensurate thickness, so the
# fixture below is chosen where the NOMINAL transport screen still passes
# but perturbing by +growth_tolerance_steps growth steps (thicker barrier,
# lower rate) pushes the missed-load probability over threshold -- fragile
# under tolerance, not merely off a round number.
_growth_fault_offset = 1 * step
r_growth = run(pkw=dict(electron_barrier_thickness_nm=thick_nom + _growth_fault_offset,
                         hole_barrier_thickness_nm=hole_thick_nom + _growth_fault_offset))
check("fault_growth_tolerance_fails_but_nominal_transport_ok",
      r_growth["rti_feasible"] is False and r_growth["rti_growth_feasible"] is False
      and r_growth["rti_transport_feasible"] is True
      and "growth_tolerance" in r_growth["rti_failed_checks"])
check("growth_diagnostics_report_nearest_commensurate_and_perturbed_margins",
      "rti_growth_nearest_commensurate_nm" in r_growth
      and "electron" in r_growth["rti_growth_nearest_commensurate_nm"]
      and "rti_growth_perturbed_margins_kT" in r_growth
      and set(r_growth["rti_growth_perturbed_margins_kT"]) == {"minus", "plus"})

# The designed stack itself (barriers at exact multiples of the growth
# step) must be growth-feasible whenever its perturbed screens both pass --
# the bug this fix replaces flagged every round-number thickness as
# "submonolayer" regardless of physics. Same (post-MEDIUM-5/6) passing
# composition as base_kwargs above (see the comment there for why
# dEv_eV_override=-0.24 and reservoir_state_count_h=2.0 replace the
# previous round's -0.29/1.0, which left no perturbation headroom).
p_designed_growth = NitrideNanowireInjectorParams(**base_kwargs)
r_designed_growth = injector_feasibility(p_designed_growth, **call_kwargs)
check("designed_2nm_stack_growth_feasible_when_perturbed_screens_pass",
      r_designed_growth["rti_growth_perturbed_margins_kT"]["minus"]["ok"]
      and r_designed_growth["rti_growth_perturbed_margins_kT"]["plus"]["ok"]
      and r_designed_growth["rti_growth_feasible"] is True)

r_no_hole = run(ckw=dict(hole_level_eV=float("nan")))
check("fault_missing_hole_path_fails_safely_not_raises",
      r_no_hole["rti_feasible"] is False
      and r_no_hole["rti_status"] == "missing_or_invalid_input"
      and r_no_hole["valid"] is False)

r_no_second_pair_ctrl = run(pkw=dict(second_pair_control_known=False))
check("fault_missing_second_pair_exclusion_fails",
      r_no_second_pair_ctrl["rti_feasible"] is False
      and r_no_second_pair_ctrl["rti_status"] == "unknown_incomplete"
      and "second_pair_control_unspecified" in r_no_second_pair_ctrl["rti_failed_checks"])

# High-transmission passive injector with NO occupation control cannot pass,
# even though its transport necessary-conditions all clear.
r_no_occ = run(pkw=dict(occupancy_control_known=False, second_pair_control_known=False))
check("passive_high_transmission_injector_without_occupancy_control_never_feasible",
      r_no_occ["rti_transport_feasible"] is True
      and r_no_occ["rti_growth_feasible"] is True
      and r_no_occ["rti_feasible"] is False
      and r_no_occ["rti_status"] == "unknown_incomplete")

r_second_pair = run(ckw=dict(second_pair_addition_meV=1.0, available_pair_rate_Hz=1e12))
check("fault_second_pair_reload_fails",
      r_second_pair["rti_feasible"] is False
      and r_second_pair["rti_second_pair_probability"] > 0.01)

# --- Independence/sensitivity of diagnostics (criterion 3) ------------------
r_T_hot = run(ckw=dict(T_K=300.0))
check("temperature_changes_margin_and_bypass",
      r_T_hot["rti_level_margin_kT"] != baseline["rti_level_margin_kT"]
      and r_T_hot["rti_bypass_fraction"] >= baseline["rti_bypass_fraction"])

r_wide_e = run(pkw=dict(electron_barrier_thickness_nm=thick_nom * 3))
check("electron_barrier_thickness_changes_electron_rate_only",
      r_wide_e["rti_e_rate_Hz"] < baseline["rti_e_rate_Hz"]
      and r_wide_e["rti_h_rate_Hz"] == baseline["rti_h_rate_Hz"])

r_wide_h = run(pkw=dict(hole_barrier_thickness_nm=thick_nom * 3))
check("hole_barrier_thickness_changes_hole_rate_only",
      r_wide_h["rti_h_rate_Hz"] < baseline["rti_h_rate_Hz"]
      and r_wide_h["rti_e_rate_Hz"] == baseline["rti_e_rate_Hz"])

r_short_window = run(ckw=dict(loading_window_ns=0.5))
check("loading_window_changes_missed_load_probability",
      r_short_window["rti_missed_load_probability"] > baseline["rti_missed_load_probability"])

# LOW 14 fix: keep the rep-rate check HONEST about what it exercises.
# Use the disc charging scale [DR] 12.126 meV throughout.
#   (a) a gate SHORTER than BOTH the 80 MHz (12.5 ns) and 200 MHz (5 ns)
#       periods must give the IDENTICAL second-pair probability at both
#       rates -- reload is priced over the explicit counting gate, not the
#       electrical period, so with an unclipped gate the rep rate must be
#       irrelevant (this is correct physics, not a bug: a previous version
#       of this check used a 10 ns gate that happens to be clipped at 200
#       MHz only, and so "passed" for testing min() rather than reload
#       physics -- see the failed-attempt notes in this module's spec).
#   (b) a gate LONGER than the 200 MHz period (so it is clipped there but
#       not at 80 MHz) must then differ between the two rates.
r_rep_unclipped_a = run(ckw=dict(rep_rate_hz=80e6, gate_ns=2.0, loading_window_ns=2.0,
                                  second_pair_addition_meV=12.126, available_pair_rate_Hz=1e6))
r_rep_unclipped_b = run(ckw=dict(rep_rate_hz=200e6, gate_ns=2.0, loading_window_ns=2.0,
                                  second_pair_addition_meV=12.126, available_pair_rate_Hz=1e6))
check("rep_rate_independent_when_gate_shorter_than_both_periods",
      abs(r_rep_unclipped_a["rti_second_pair_probability"]
          - r_rep_unclipped_b["rti_second_pair_probability"]) < 1e-12)

r_rep_clipped_a = run(ckw=dict(rep_rate_hz=80e6, gate_ns=10.0, loading_window_ns=2.0,
                                second_pair_addition_meV=12.126, available_pair_rate_Hz=1e6))
r_rep_clipped_b = run(ckw=dict(rep_rate_hz=200e6, gate_ns=10.0, loading_window_ns=2.0,
                                second_pair_addition_meV=12.126, available_pair_rate_Hz=1e6))
check("rep_rate_changes_second_pair_probability_when_gate_clipped",
      r_rep_clipped_a["rti_second_pair_probability"] != r_rep_clipped_b["rti_second_pair_probability"])

# HIGH 2 fix: allowed-state ALIGNMENT error is versus the carrier's own
# EMITTER quasi-Fermi level (mu_e from n_cm3), NOT the dot level -- run on
# a genuine double-barrier resonance (base_kwargs' single-barrier fixture
# has alignment error identically 0 by construction and cannot exercise
# this gate; LOW 9 fix, replaces the previous vacuous
# "alignment_error_reduces_orbital_margin" check on that fixture).
# electron_level_eV no longer controls alignment at all (HIGH 2) -- it is
# used ONLY for the informational rti_well_to_dot_drop_meV.  Misalignment
# is instead constructed by choosing n_cm3 (electron reservoir doping) so
# that mu_e sits a KNOWN distance from the zero-field resonance E_res0,
# via the finite-T Fermi-Dirac forward formula n = N_c F_1/2(mu/kT) (DOS
# mass) evaluated fresh here at the call's own T_K (audit M1, 2026-09-23:
# the production emitter level is now that Fermi-Dirac inversion, no
# longer the T=0 formula this fixture used to invert) -- never by calling
# the production _electron_quasi_fermi_eV to manufacture the target.
p_align_probe = NitrideNanowireInjectorParams(
    hole_topology="single_barrier", hole_barrier_thickness_nm=thick_nom,
    mh_barrier_override=0.30, dEv_eV_override=-0.30,
    occupancy_control_known=True, second_pair_control_known=True,
    include_polarization=False,
)
E_res_probe, _T_res_probe = _find_resonance(p_align_probe, "electron", 0.0, 0.0)


def _n_cm3_for_mu_hand(mu_eV, m_ratio):
    k_f = math.sqrt(2.0 * m_ratio * M0_KG * mu_eV * EV_J) / HBAR_JS
    n_m3 = k_f ** 3 / (3.0 * math.pi ** 2)
    return n_m3 / 1e6


_me_well_probe = p_align_probe.me_well
_n_cm3_big = _n_cm3_for_mu_fd_indep(E_res_probe + 0.146, _me_well_probe, 230.0)
_n_cm3_small = _n_cm3_for_mu_fd_indep(E_res_probe + 0.002, _me_well_probe, 230.0)

align_call = dict(T_K=230.0, rep_rate_hz=80e6, loading_window_ns=12.4,
                   hole_level_eV=0.005, electron_level_eV=0.02,
                   electron_spacing_meV=600.0, hole_spacing_meV=600.0,
                   second_pair_addition_meV=600.0, available_pair_rate_Hz=1.0)

p_misaligned = replace(p_align_probe, n_cm3=_n_cm3_big)
r_misaligned = injector_feasibility(p_misaligned, **align_call)
check("146meV_misalignment_fails_feasibility",
      r_misaligned["rti_feasible"] is False
      and abs(r_misaligned["rti_alignment_error_e_meV"] - 146.0) < 1e-3
      and "electron_alignment" in r_misaligned["rti_failed_checks"])

p_small_misalign = replace(p_align_probe, n_cm3=_n_cm3_small)
r_small_misalign = injector_feasibility(p_small_misalign, **align_call)
check("2meV_misalignment_passes_alignment_gate",
      abs(r_small_misalign["rti_alignment_error_e_meV"] - 2.0) < 1e-3
      and "electron_alignment" not in r_small_misalign["rti_failed_checks"])

# rti_well_to_dot_drop_meV is purely informational and independent of the
# alignment gate: changing electron_level_eV must not change the reported
# alignment error or the alignment_ok gate above.
r_misaligned_other_dot = injector_feasibility(replace(p_align_probe, n_cm3=_n_cm3_big),
                                               **{**align_call, "electron_level_eV": 0.30})
check("well_to_dot_drop_is_informational_only",
      abs(r_misaligned_other_dot["rti_alignment_error_e_meV"]
          - r_misaligned["rti_alignment_error_e_meV"]) < 1e-6
      and r_misaligned_other_dot["rti_well_to_dot_drop_meV"]["electron"]
      != r_misaligned["rti_well_to_dot_drop_meV"]["electron"])

# MEDIUM 10 fix: alignment_tunable is gated on the ENGINE's own bias scan
# (a genuine sweep of the production resonance finder over bias, see
# _scan_required_bias_shift_meV), not an assumed 1:1 lever arm -- a
# generous range finds a bias that clears the gate; a narrow range does
# not.  The exact required shift is whatever the engine's sweep finds
# (not necessarily 146 meV 1:1), so this checks the GATING behavior, not a
# hardcoded lever arm.
p_align_tunable_ok = replace(p_align_probe, n_cm3=_n_cm3_big,
                              alignment_tunable=True, bias_tuning_range_meV=200.0)
r_tunable_ok = injector_feasibility(p_align_tunable_ok, **align_call)
check("alignment_tunable_within_bias_range_clears_alignment_gate",
      "electron_alignment" not in r_tunable_ok["rti_failed_checks"]
      and math.isfinite(r_tunable_ok["rti_required_bias_shift_meV"]))

p_align_tunable_bad = replace(p_align_probe, n_cm3=_n_cm3_big,
                               alignment_tunable=True, bias_tuning_range_meV=50.0)
r_tunable_bad = injector_feasibility(p_align_tunable_bad, **align_call)
check("alignment_tunable_beyond_bias_range_still_fails_alignment_gate",
      r_tunable_bad["rti_feasible"] is False
      and "electron_alignment" in r_tunable_bad["rti_failed_checks"])

# --- Rate normalization: independent analytic case (fully transparent,
#     zero-temperature ballistic supply saturates at g_s*mu/h) -------------
p_ballistic = NitrideNanowireInjectorParams(
    electron_topology="single_barrier", electron_barrier_thickness_nm=1e-6,
    dEc_eV_override=0.0, me_barrier_override=NitrideNanowireInjectorParams().me_well,
    n_cm3=3.0e18,
)
# audit M1 (2026-09-23): the production emitter level is the Fermi-Dirac
# inversion with the DOS mass; at 4 K it equals the T=0 degenerate level
# with that mass (checked separately below), computed fresh here.
mu_b = _mu_t0_indep_eV(p_ballistic.n_cm3 * 1e6, _mdos_e_indep(p_ballistic.me_well))
kT_cold = KB_EV * 4.0  # near-T=0 so Fermi-Dirac approximates a hard step
from fsim_core.nitride_nanowire_injector import H_EVS as _H_EVS  # noqa: E402


def _fd_fresh(E, mu, kT):
    x = (E - mu) / kT
    return 1.0 / (1.0 + math.exp(min(max(x, -60), 60)))


Es_b = np.linspace(1e-6, mu_b + 20 * kT_cold, 20000)
Ts_b = transmission(p_ballistic, Es_b, carrier="electron")
fs_b = np.array([_fd_fresh(E, mu_b, kT_cold) for E in Es_b])
rate_independent = p_ballistic.degeneracy * simpson(Ts_b * fs_b, x=Es_b) / _H_EVS
rate_expected_T0 = p_ballistic.degeneracy * mu_b / _H_EVS
check("rate_normalization_ballistic_case_matches_independent_quadrature",
      abs(rate_independent / rate_expected_T0 - 1.0) < 0.02)

r_ballistic = run(pkw=dict(electron_topology="single_barrier",
                            electron_barrier_thickness_nm=1e-6,
                            dEc_eV_override=0.0,
                            me_barrier_override=NitrideNanowireInjectorParams().me_well,
                            n_cm3=3.0e18),
                   ckw=dict(T_K=4.0, electron_level_eV=1e-6))
# MEDIUM 8 fix: the production rate is g_s (degeneracy) * channels
# (reservoir_state_count_e) * Landauer integral -- the channel count
# defaults to 2.0 (unset by base_kwargs/p_ballistic), so it must be
# included in the independent T=0 expectation too.
_p_ballistic_from_run = NitrideNanowireInjectorParams(**{**base_kwargs,
    "electron_topology": "single_barrier", "electron_barrier_thickness_nm": 1e-6,
    "dEc_eV_override": 0.0, "me_barrier_override": NitrideNanowireInjectorParams().me_well,
    "n_cm3": 3.0e18})
rate_expected_T0_production = _p_ballistic_from_run.reservoir_state_count_e * rate_expected_T0
check("rate_normalization_matches_production_within_5pct",
      abs(r_ballistic["rti_e_rate_Hz"] / rate_expected_T0_production - 1.0) < 0.05)

# --- Out-of-range / nonfinite inputs fail safely (never raise) -------------
try:
    r_bad1 = run(ckw=dict(T_K=float("nan")))
    r_bad2 = run(ckw=dict(T_K=-5.0))
    r_bad3 = run(ckw=dict(loading_window_ns=5.0, rep_rate_hz=1e9))  # window > period
    safe = (r_bad1["rti_feasible"] is False and r_bad2["rti_feasible"] is False
            and r_bad3["rti_feasible"] is False and all(
                not v["valid"] for v in (r_bad1, r_bad2, r_bad3)))
except Exception as exc:  # pragma: no cover - failure of the safety contract itself
    print("FAIL out_of_range_inputs_raised:", exc)
    safe = False
check("out_of_range_and_nonfinite_inputs_fail_safely_not_raise", safe)

# --- No hardcoded feasible=True switch: params defaults must not pass -----
r_defaults = injector_feasibility(
    NitrideNanowireInjectorParams(occupancy_control_known=True, second_pair_control_known=True),
    T_K=300.0, rep_rate_hz=200e6, loading_window_ns=0.1,
    electron_level_eV=0.05, hole_level_eV=0.01,
    electron_spacing_meV=float("nan"), hole_spacing_meV=float("nan"),
    second_pair_addition_meV=10.0, available_pair_rate_Hz=1e9,
)
check("default_designed_stack_is_not_hardcoded_feasible_true",
      r_defaults["rti_feasible"] is False)

# --- MEDIUM 6/7: NaN linewidth handling -------------------------------------
# Single-barrier path: no resonance exists, so its linewidth is NOT
# APPLICABLE (reported NaN) rather than a failure, and the combined width
# used for margin/rate broadening is the alignment uncertainty ALONE (no
# quadrature double-count against itself).
check("single_barrier_linewidth_not_applicable_but_screen_still_valid",
      math.isnan(baseline["rti_linewidth_meV"]) and baseline["valid"] is True)

# Double-barrier path: a resonance IS found, so its linewidth must be a
# finite number gating validity if it cannot be resolved.  HIGH 1 fix
# (Opus re-review of 12b39cd, 2026-09-14): the DEFAULT (Ga-polar, 0.70 eV
# partition) hole path DOES resolve real resonances once the window uses
# the tilted barrier's correct (highest-face) top -- three of them, at
# 48.66 / 101.25 / 166.05 meV (see section 2c above) -- so "no resolvable
# hole resonance at defaults" was an artefact of the previous (wrong)
# window, not a real device consequence.  This fixture still deliberately
# deepens the hole barrier further (dEv_eV_override=-0.5) for an unrelated
# reason: it gives a single, well-separated, generously-transmitting
# resonance for BOTH carriers so this specific check exercises the clean
# "found, finite linewidth" path rather than the default stack's much
# narrower (T ~4.7e-9) ground state.
p_db_linewidth = NitrideNanowireInjectorParams(
    dEv_eV_override=-0.5, mh_barrier_override=0.30,
    occupancy_control_known=True, second_pair_control_known=True)
r_db_linewidth = injector_feasibility(
    p_db_linewidth, T_K=230.0, rep_rate_hz=80e6, loading_window_ns=2.0,
    electron_level_eV=0.269, hole_level_eV=0.0004,
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=600.0, available_pair_rate_Hz=1.0,
)
check("double_barrier_resonance_reports_a_finite_linewidth",
      math.isfinite(r_db_linewidth["rti_linewidth_meV"]) and r_db_linewidth["rti_linewidth_meV"] > 0.0)

# --- LOW 10: honestly-named alias for the second-carrier upper bound -------
check("second_carrier_probability_alias_matches_documented_upper_bound",
      baseline["rti_second_carrier_probability"] == baseline["rti_second_pair_probability"])

# --- MEDIUM 8 / LOW 13: exercise injector_feasibility with the actual
#     DOUBLE-BARRIER designed stack (both carriers) at T 230/300 K and rep
#     80/200 MHz, E_C 12.126 meV, gate 2 ns, under both valence partitions
#     AND both polarities, printing every rti_* column -- the base_kwargs
#     fixture above is single_barrier for both carriers and never exercised
#     a real double-barrier feasibility call at all.  LOW 13 fix: uses the
#     REAL dot levels (from the levels module's dE_e/dE_h at the Deshpande
#     geometry, transcribed here as literals since this module does not
#     import the levels solver) instead of feeding the found resonance back
#     as the level -- the previous version made alignment trivially exact
#     by construction; the printed verdict is now the device's, not a
#     tautology.  Under HIGH 2 the electron_level_eV/hole_level_eV inputs
#     no longer gate alignment at all (that is now versus the emitter mu),
#     so this is safe regardless of whether a resonance is found.
_ELECTRON_DOT_LEVEL_EV = -0.6995   # [V] nitride_nanowire_levels dE_e at the Deshpande disc-in-wire geometry
_HOLE_DOT_LEVEL_EV = -0.4164       # [V] nitride_nanowire_levels dE_h at the Deshpande disc-in-wire geometry
# MEDIUM 7 fix (Opus re-review of 12b39cd, 2026-09-14): available_pair_rate_Hz
# must be a PHYSICALLY PLAUSIBLE supply rate, not 1.0 -- at 1.0 Hz the
# second-pair screen never fires (P ~3e-9) regardless of E_C, hiding the
# directive's honest failure line entirely.  1e9 Hz [A, the transport
# module's own r_captured order of magnitude at 1 nA bias current, see
# fsim_core.nitride_nanowire_transport] makes the screen fire at the
# disc's own E_C/kT (0.61 / 0.47 at 230 / 300 K), producing the
# "second_pair_reload (E_C below kT)" tag the fix round requires.
_designed_call_common = dict(
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=12.126, available_pair_rate_Hz=1e9,
    electron_level_eV=_ELECTRON_DOT_LEVEL_EV, hole_level_eV=_HOLE_DOT_LEVEL_EV,
    gate_ns=2.0,
)
_designed_results = {}
for _polarity in ("Ga", "N"):
    for _delta_ev in (0.70, _TSAI_BAYRAM_DELTA_EV_EV):
        _p_designed = NitrideNanowireInjectorParams(
            delta_Ev_GaN_AlN_eV=_delta_ev, polarity=_polarity,
            occupancy_control_known=True, second_pair_control_known=True,
        )
        for _T in (230.0, 300.0):
            for _rep in (80e6, 200e6):
                _r = injector_feasibility(
                    _p_designed, T_K=_T, rep_rate_hz=_rep, loading_window_ns=2.0,
                    **_designed_call_common,
                )
                _designed_results[(_polarity, _delta_ev, _T, _rep)] = _r
                _cols = " ".join(f"{k}={_r[k]}" for k in sorted(_r)
                                  if k.startswith("rti_") or k == "valid")
                print(f"designed stack polarity={_polarity} partition={_delta_ev:.2f}eV "
                      f"T={_T:.0f}K rep={_rep/1e6:.0f}MHz | {_cols}")
# Strain-mass audit (2026-09-23): under the Tsai 0.30 eV partition the
# hole window tops out at 90 meV and its transmission has no resolvable
# resonance (the 54.5 meV "peak" the old first-local-max rule accepted is a
# 1.06x ripple, T 7.7e-9 over a 7.1e-9 valley), so the Tsai-partition
# bypass is honestly NaN (unknown), never a fabricated number.  The check
# therefore requires finiteness OR a hole path the production finder
# reports as unresolved under that partition.
check("designed_stack_double_barrier_both_carriers_both_partitions_finite_diagnostics_or_explicitly_unresolved",
      all(math.isfinite(_r["rti_bypass_fraction_tsai_partition"])
          or _find_resonance(replace(NitrideNanowireInjectorParams(polarity=_k[0]),
                                     delta_Ev_GaN_AlN_eV=_TSAI_BAYRAM_DELTA_EV_EV), "hole", 0.0, 0.0) is None
          for _k, _r in _designed_results.items())
      and all(math.isfinite(_r["rti_bypass_fraction"]) for _k, _r in _designed_results.items() if _k[1] == 0.70))
check("designed_stack_rate_or_resonance_unresolved_is_explicit_never_silently_zero_error",
      all(math.isfinite(_r["rti_e_rate_Hz"]) or "electron_resonance_unresolved" in _r["rti_failed_checks"]
          for _r in _designed_results.values())
      and all(math.isfinite(_r["rti_h_rate_Hz"]) or "hole_resonance_unresolved" in _r["rti_failed_checks"]
              for _r in _designed_results.values()))
check("designed_stack_no_hardcoded_feasible_true_and_rti_feasible_false_is_diagnosed",
      all((not _r["rti_feasible"]) or (len(_r["rti_failed_checks"]) == 0)
          for _r in _designed_results.values())
      and any(_r["rti_failed_checks"] for _r in _designed_results.values()))
check("designed_stack_second_pair_reload_tag_fires_with_plausible_supply_rate",
      any("second_pair_reload (E_C below kT)" in _r["rti_failed_checks"]
          for _r in _designed_results.values()))

# LOW 9 fix (Opus re-review of 12b39cd, 2026-09-14): pin the designed
# stack's electron alignment error -- |E_res - mu_e| at the default 0.70 eV
# partition, Ga polarity, 230 K -- against the reviewer's independently
# derived literal (232.63 meV), not merely printed and never asserted.
# audit M1 re-pin (2026-09-23): 232.63 meV was |E_res - mu_e| with the old
# T=0 me_z emitter level (36.29 meV at every T), so the reviewer's literal
# fixes E_res = 232.63 meV + mu_T0(me_z).  With the Fermi-Dirac emitter
# level (DOS mass) the expected error is E_res - mu_FD(T), both terms
# computed fresh in this file: 240.00 meV at 230 K and 247.51 meV at 300 K.
# Strain-mass audit re-pin (2026-09-23): the along-c electron masses move
# the designed stack's ground resonance to E_res = 291.104 meV (pinned above,
# was 268.92) and the bulk-GaN DOS mass (0.2010, was 0.1934) moves mu_FD to
# 26.98 / 19.21 meV, so the expected error is 264.12 / 271.90 meV at 230 /
# 300 K (was 240.00 / 247.51).
_r_designed_default = _designed_results[("Ga", 0.70, 230.0, 80e6)]
_E_res_designed_meV = 291.104
for _T_al in (230.0, 300.0):
    _exp_al = _E_res_designed_meV - _mu_fd_indep_eV(3.0e18, NitrideNanowireInjectorParams().me_well, _T_al) * 1e3
    check(f"designed_stack_electron_alignment_error_matches_E_res_minus_FD_mu_{_T_al:.0f}K",
          abs(_designed_results[("Ga", 0.70, _T_al, 80e6)]["rti_alignment_error_e_meV"] - _exp_al) < 0.01)

# HIGH 2 fix (Opus re-review of 12b39cd, 2026-09-14): pin the designed
# stack's overall thermionic-bypass fraction (0.70 eV partition, Ga
# polarity -- the hole path dominates this max(e, h) at both temperatures)
# at 230 / 300 K against the reviewer's independently derived literals
# (0.0021 / 0.0151, computed from the ACTUAL tilted profile maximum, not
# the flat-band barrier height that previously flipped this screen from a
# ~40-95x fail to a spurious pass).
check("designed_stack_bypass_matches_pinned_literal_230K_0.0021",
      abs(_designed_results[("Ga", 0.70, 230.0, 80e6)]["rti_bypass_fraction"] / 0.0021 - 1.0) < 0.02)
check("designed_stack_bypass_matches_pinned_literal_300K_0.0151",
      abs(_designed_results[("Ga", 0.70, 300.0, 80e6)]["rti_bypass_fraction"] / 0.0151 - 1.0) < 0.02)

# --- AUDIT M1 (2026-09-23): electron emitter quasi-Fermi level ------------
# The production level must satisfy n = N_c F_1/2(mu/kT) with the DOS mass
# (F_1/2 evaluated independently by quad above) to 1e-4 relative at 230 and
# 300 K, match the Fermi-Dirac values 26.98 / 19.21 meV (strain-mass audit
# re-pin, 2026-09-23: bulk GaN DOS mass (0.209^2*0.186)^(1/3) = 0.2010 from
# the corrected Rinke axes, computed independently by _mu_fd_indep_eV; the
# transport audit M1 values were 28.9 / 21.4 meV with the old 0.1934), and
# reduce to the T=0 degenerate formula
# within 1 % at 4 K for strongly degenerate densities.
from fsim_core.nitride_nanowire_injector import _electron_quasi_fermi_eV  # noqa: E402
_p_fd = NitrideNanowireInjectorParams()
_mdos_fd = _mdos_e_indep(_p_fd.me_well)
for _T_fd, _audit_meV in ((230.0, 26.98), (300.0, 19.21)):
    _mu_prod = _electron_quasi_fermi_eV(_p_fd, _T_fd)
    _kT_fd = _KB_J_INDEP * _T_fd / EV_J
    _n_back = _nc_indep_m3(_mdos_fd, _T_fd) * _f_half_indep(_mu_prod / _kT_fd)
    check(f"M1_electron_mu_satisfies_n_eq_Nc_F12_{_T_fd:.0f}K_1e-4",
          abs(_n_back / (_p_fd.n_cm3 * 1e6) - 1.0) < 1e-4)
    check(f"M1_electron_mu_matches_FD_value_{_audit_meV}meV_{_T_fd:.0f}K",
          abs(_mu_prod * 1e3 - _audit_meV) < 0.01
          and abs(_mu_fd_indep_eV(_p_fd.n_cm3, None, _T_fd) * 1e3 - _audit_meV) < 0.01)
    # Checker condition (2026-09-23): the emitter is bulk n-GaN, so a
    # me_well (well tunnelling mass) override must not move its level.
    check(f"emitter_mu_independent_of_me_well_override_{_T_fd:.0f}K",
          _electron_quasi_fermi_eV(replace(_p_fd, me_well=0.35), _T_fd) == _mu_prod
          and _electron_quasi_fermi_eV(replace(_p_fd, me_well=0.10), _T_fd) == _mu_prod)
for _n_deg in (3.0e18, 2.0e19):
    _p_deg = replace(_p_fd, n_cm3=_n_deg)
    _mu4 = _electron_quasi_fermi_eV(_p_deg, 4.0)
    _mu0 = _mu_t0_indep_eV(_n_deg * 1e6, _mdos_fd)
    check(f"M1_electron_mu_reduces_to_T0_formula_within_1pct_at_4K_n{_n_deg:.0e}",
          abs(_mu4 / _mu0 - 1.0) < 0.01)

# --- AUDIT M2 (2026-09-23): unitarity scan covers the Landauer window -----
# Instrument the module's quad (records the rate integral's upper limit)
# and _transmission_scalar (records every energy the numerics scan visits)
# on the designed hole path at 230 K; the scan must reach exactly the
# integral's upper bound, which the audit computed independently as
# max(tilted top, mu_h) + 15 kT = 0.802 eV (flat-band bound was 0.507 eV).
import fsim_core.nitride_nanowire_injector as _mod_m2  # noqa: E402
_p_m2 = NitrideNanowireInjectorParams()
_kT_m2 = KB_EV * 230.0
_ph_m2 = _resolve_path(_p_m2, "hole")
_mu_h_m2 = _mod_m2._hole_quasi_fermi_eV(_p_m2, 230.0)[1]
_res_h_m2 = _find_resonance(_p_m2, "hole", 0.0, 0.0)
_Ec_h_m2 = _res_h_m2[0] if _res_h_m2 is not None else float("nan")
_w_m2 = _p_m2.alignment_uncertainty_meV / 1000.0
_quad_uppers = []
_orig_quad = _mod_m2.quad
_orig_ts = _mod_m2._transmission_scalar


def _rec_quad(f, a, b, *args, **kw):
    _quad_uppers.append(b)
    return _orig_quad(f, a, b, *args, **kw)


_scan_Es = []


def _rec_ts(params, E, *args, **kw):
    _scan_Es.append(E)
    return _orig_ts(params, E, *args, **kw)


try:
    _mod_m2.quad = _rec_quad
    _mod_m2._forward_rate_hz(_p_m2, "hole", 0.0, 0.0, _Ec_h_m2, _w_m2, _kT_m2, _mu_h_m2)
    _mod_m2.quad = _orig_quad
    _mod_m2._transmission_scalar = _rec_ts
    _scan_ok_m2 = _mod_m2._numerics_scan_ok(_p_m2, "hole", 0.0, _mu_h_m2, _w_m2, _kT_m2,
                                            E_center=_Ec_h_m2)
finally:
    _mod_m2.quad = _orig_quad
    _mod_m2._transmission_scalar = _orig_ts
_landauer_hi = _quad_uppers[0] if _quad_uppers else float("nan")
_scan_hi = max(_scan_Es) if _scan_Es else float("nan")
print(f"M2 hole 230 K: Landauer upper {_landauer_hi:.4f} eV, scan upper {_scan_hi:.4f} eV, "
      f"flat-band bound {max(_ph_m2.barrier_height_eV, _mu_h_m2) + 15 * _kT_m2:.4f} eV")
check("M2_numerics_scan_upper_bound_equals_landauer_upper_bound_hole_230K",
      math.isfinite(_landauer_hi) and abs(_scan_hi - _landauer_hi) < 1e-12 and _scan_ok_m2 is True)
check("M2_scan_upper_bound_matches_audit_0.802eV_not_flat_band_0.507eV",
      abs(_scan_hi - 0.802) < 0.001
      and _scan_hi > max(_ph_m2.barrier_height_eV, _mu_h_m2) + 15 * _kT_m2 + 0.2)

# --- AUDIT L6 (2026-09-23): tsai-partition max() is NaN-safe --------------
_src_l6 = Path(_mod_m2.__file__).read_text(encoding="utf-8")
check("L6_tsai_partition_bypass_uses_nan_safe_max",
      "float(max(e_tsai" not in _src_l6 and "_nan_safe_max(e_tsai[\"bypass_fraction\"]" in _src_l6)
# Behavioural L6 case: at the vertical-photonic depletion field (~79.8
# kV/cm, 300 K) the Tsai-partition HOLE resonance is unresolvable (no
# resonance, so its rate and bypass are NaN) while the electron path is
# finite; max(finite, nan) used to return the electron value and hide the
# missing hole number.  Expectation: NaN (a missing carrier is missing).
_p_l6 = NitrideNanowireInjectorParams(occupancy_control_known=True, second_pair_control_known=True)
_p_l6_tsai = replace(_p_l6, delta_Ev_GaN_AlN_eV=_TSAI_BAYRAM_DELTA_EV_EV)
_r_l6 = injector_feasibility(_p_l6, T_K=300.0, rep_rate_hz=80e6, loading_window_ns=0.1,
                             electron_level_eV=-0.7, hole_level_eV=-0.42,
                             electron_spacing_meV=600.0, hole_spacing_meV=600.0,
                             second_pair_addition_meV=12.126, available_pair_rate_Hz=1e9,
                             field_kVcm=79.8)
check("L6_tsai_partition_bypass_is_nan_when_tsai_hole_resonance_unresolved",
      _find_resonance(_p_l6_tsai, "hole", 0.0, 79.8) is None
      and _find_resonance(_p_l6_tsai, "electron", 0.0, 79.8) is not None
      and math.isnan(_r_l6["rti_bypass_fraction_tsai_partition"])
      and math.isfinite(_r_l6["rti_bypass_fraction"]))

# --- AUDIT C6 (2026-09-23, M3): loading window is an explicit parameter ----
# User decision Q6: the loading window is arbitrary (single photons in ANY
# window).  Expectation (internal audit): the second-carrier probability
# over a window t is 1 - exp(-r_2 t) with
# r_2 = r_e * exp(-E_C/kT), r_e = 4.85e7 Hz [E] (the sweep row's electron
# supply, audit), E_C = 12.126 meV [DR] (isolated-sphere e^2/C of the
# 12.5 nm GaN disc, drive_mech convention), T = 300 K -> 0.0030 / 0.0299 /
# 0.1408 at 0.1 / 1 / 5 ns (audit's printed values, pinned literally too).
# Fixture [A]: single-barrier paths with a negligible (1e-6 meV) alignment
# width so the charging energy is not width-subtracted (the audit's w = 0),
# available_pair_rate_Hz = r_e caps the (much faster) electron tunnel rate
# at the audit's r_e, and p_cm3 = 1e14 puts the hole supply at ~2.6e3 Hz
# (the audit's r_h 1.8e3 Hz class) so the second-HOLE term is < 1e-4.
_E_C_C6_MEV = 12.126017244990036
_R_E_C6_HZ = 4.85e7
_kT300_c6_eV = _KB_J_INDEP * 300.0 / 1.602176634e-19   # [V] CODATA 2018, exact
_r2_c6 = _R_E_C6_HZ * math.exp(-(_E_C_C6_MEV / 1000.0) / _kT300_c6_eV)
_AUDIT_C6 = {0.1: 0.0030, 1.0: 0.0299, 5.0: 0.1408}
_p_c6 = NitrideNanowireInjectorParams(electron_topology="single_barrier", hole_topology="single_barrier",
                                      alignment_uncertainty_meV=1e-6, growth_tolerance_steps=0.0,
                                      p_cm3=1e14)


def _c6(window_ns, gate_ns=None):
    kw = {} if gate_ns is None else {"gate_ns": gate_ns}
    return injector_feasibility(_p_c6, T_K=300.0, rep_rate_hz=80e6, loading_window_ns=window_ns,
                                electron_level_eV=-0.2, hole_level_eV=-0.2,
                                electron_spacing_meV=float("nan"), hole_spacing_meV=float("nan"),
                                second_pair_addition_meV=_E_C_C6_MEV, available_pair_rate_Hz=_R_E_C6_HZ,
                                **kw)


_c6_rows = {t: _c6(t) for t in (0.1, 1.0, 5.0)}
check("C6_audit_r2_closed_form_reproduces_audit_values_0.0030_0.0299_0.1408",
      all(abs((1.0 - math.exp(-_r2_c6 * t * 1e-9)) - v) < 1e-3 for t, v in _AUDIT_C6.items()))
check("C6_second_pair_probability_over_loading_window_equals_1_minus_exp_r2_t (0.1/1/5 ns, 1e-3 abs)",
      all(abs(_c6_rows[t]["rti_second_pair_probability"] - (1.0 - math.exp(-_r2_c6 * t * 1e-9))) < 1e-3
          and abs(_c6_rows[t]["rti_second_pair_probability"] - _AUDIT_C6[t]) < 1e-3
          and _c6_rows[t]["rti_gate_ns"] == t
          for t in _AUDIT_C6))
check("C6_second_pair_probability_strictly_monotone_in_window",
      _c6_rows[0.1]["rti_second_pair_probability"] < _c6_rows[1.0]["rti_second_pair_probability"]
      < _c6_rows[5.0]["rti_second_pair_probability"])
# Missed load is priced over the LOADING window, never the counting gate:
# at a fixed 0.1 ns window it is bit-identical whatever gate_ns is, and it
# equals the Poisson closed form 1-(1-e^{-r_e t})(1-e^{-r_h t}).  With the
# hole supply at ~kHz the pair-missed-load FAILURE is window-independent
# over 0.1..5 ns (a 99 % hole load needs -ln(0.01)/r_h ~ ms; audit M3).
_c6_gate = [_c6(0.1, gate_ns=g) for g in (0.1, 1.0, 5.0)]
check("C6_missed_load_independent_of_counting_gate_at_fixed_loading_window",
      len({r["rti_missed_load_probability"] for r in _c6_gate}) == 1
      and _c6_gate[0]["rti_second_pair_probability"] < _c6_gate[2]["rti_second_pair_probability"])
check("C6_missed_load_is_poisson_closed_form_and_fails_at_every_window (hole-limited, window-independent)",
      all(abs(r["rti_missed_load_probability"]
              - (1.0 - (1.0 - math.exp(-r["rti_e_rate_Hz"] * t * 1e-9))
                 * (1.0 - math.exp(-r["rti_h_rate_Hz"] * t * 1e-9)))) < 1e-12
          and r["rti_missed_load_probability"] > 0.99
          and "pair_missed_load_window" in r["rti_failed_checks"]
          and -math.log(0.01) / r["rti_h_rate_Hz"] > 1e-6
          for t, r in _c6_rows.items()))

# --- No production routine used to derive its own expected value ----------
import fsim_core.nitride_nanowire_injector as _mod  # noqa: E402
_src = Path(_mod.__file__).read_text(encoding="utf-8")
_import_lines = [ln for ln in _src.splitlines() if ln.strip().startswith(("import ", "from "))]
check("source_does_not_import_drive_mech_or_device_or_levels",
      not any(("drive_mech" in ln or "fsim_core.device" in ln or "nitride_nanowire_levels" in ln
               or "_photonics" in ln or "_surface" in ln)
              for ln in _import_lines))
check("source_ascii_only",
      all(ord(ch) < 128 for ch in _src))

print(f"{sum(c)}/{len(c)} nitride nanowire injector checks passed")
raise SystemExit(0 if all(c) else 1)
