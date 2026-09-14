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

m_w = 0.209   # [V] GaN conduction mass, Rinke PRB 2008 (via nitride_materials)
m_b = 0.329   # [V] AlN conduction mass, Rinke PRB 2008 (via nitride_materials)


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

# HIGH 3/4/5: the potential is CONTINUOUS at every interface (no reset) --
# probe the module's own transmission() profile indirectly is not possible
# (V(x) is internal), so reconstruct the SAME running-offset arithmetic
# fresh here (independent of _stack_polarization_fields_eV_per_m) and check
# it has no interface jump, using the SAME material segment list a
# double-barrier electron path has (barrier/well/barrier).
_pol_fields_hand = [E_b_hand, E_w_hand, E_b_hand]
_lengths_hand = [d_b, d_w, d_b]
_flat_heights_hand = [0.60615, 0.0, 0.60615]   # AlGaN barrier / GaN well / AlGaN barrier, x=0.30, 0.70 eV partition
_pos = 0.0
_pol_offset = 0.0
_profile_points = []   # (V_total_eV at entry, V_total_eV at exit) per segment
for _L, _Efield, _Vflat in zip(_lengths_hand, _pol_fields_hand, _flat_heights_hand):
    _entry = _Vflat + _pol_offset
    _pol_offset += _Efield * _L
    _exit = _Vflat + _pol_offset
    _profile_points.append((_entry, _exit))
_max_jump_eV = 0.0
for i in range(len(_profile_points) - 1):
    _flat_step = _flat_heights_hand[i + 1] - _flat_heights_hand[i]
    _gap = (_profile_points[i + 1][0] - _profile_points[i][1]) - _flat_step
    _max_jump_eV = max(_max_jump_eV, abs(_gap))
check("polarization_potential_continuous_at_every_interface",
      _max_jump_eV < 1e-9)
# Ga-polar: exit face ABOVE entry face for the electron barrier (never abs()).
check("ga_polar_electron_barrier_exit_above_entry",
      _profile_points[0][1] > _profile_points[0][0])
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
    for _j in range(1, len(T_dense) - 1):
        if T_dense[_j] > 1e-9 and T_dense[_j] >= T_dense[_j - 1] and T_dense[_j] >= T_dense[_j + 1]:
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
# MEDIUM 6 fix note: base_kwargs/call_kwargs are re-tuned for this round.
# With the hole reservoir's quasi-Fermi level now coming from the honest
# Mg-acceptor charge-neutrality solve (mu_h ~ -144 meV at 230 K, a dilute
# non-degenerate hole gas) instead of the old always-positive degenerate
# estimate, the hole path's capture rate is far smaller than the previous
# fixture assumed -- the hole barrier here is deliberately shallower/
# thinner than the electron barrier so the synthetic PASSING case still
# clears rate/bypass/margin under the corrected physics; single_barrier +
# include_polarization=False keeps this fixture's numbers independent of
# the polarization closure (exercised separately in section 2a).
base_kwargs = dict(
    electron_topology="single_barrier", hole_topology="single_barrier",
    electron_barrier_thickness_nm=thick_nom, hole_barrier_thickness_nm=hole_thick_nom,
    me_barrier_override=0.25, mh_barrier_override=0.20,
    dEc_eV_override=0.30, dEv_eV_override=-0.29,
    occupancy_control_known=True, second_pair_control_known=True,
    growth_tolerance_steps=1.0, include_polarization=False,
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
# "submonolayer" regardless of physics. Uses the same (post-MEDIUM-6)
# passing composition as base_kwargs above, since a symmetric literal-2.0nm
# hole barrier no longer clears rate/bypass under the corrected (much more
# dilute) hole quasi-Fermi level.
p_designed_growth = NitrideNanowireInjectorParams(
    electron_topology="single_barrier", hole_topology="single_barrier",
    electron_barrier_thickness_nm=thick_nom, hole_barrier_thickness_nm=hole_thick_nom,
    me_barrier_override=0.25, mh_barrier_override=0.20,
    dEc_eV_override=0.30, dEv_eV_override=-0.29,
    occupancy_control_known=True, second_pair_control_known=True,
    growth_tolerance_steps=1.0,
)
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
# inverting the SAME degenerate free-electron-gas formula fresh here
# (mu = hbar^2 k_f^2/2m, k_f=(3 pi^2 n)^(1/3)) -- never by calling the
# production _degenerate_mu_eV to manufacture the target.
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
_n_cm3_big = _n_cm3_for_mu_hand(E_res_probe + 0.146, _me_well_probe)
_n_cm3_small = _n_cm3_for_mu_hand(E_res_probe + 0.002, _me_well_probe)

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
path_b = _resolve_path(p_ballistic, "electron")
mu_b = path_b.mu_eV
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
# finite number gating validity if it cannot be resolved.  The DEFAULT
# (Ga-polar, 0.70 eV partition) hole barrier is only 0.21 eV tall -- less
# than the 0.2951 eV polarization swing across a single 2 nm barrier --
# so it has no resolvable resonance at all under Attempt-4 physics (a
# real, honest consequence of HIGH 3/4/5, not a bug: see
# hole_resonance_unresolved coverage below); deepen the hole barrier here
# (dEv_eV_override) so BOTH carriers resolve a resonance and this check
# exercises the genuine "found, finite linewidth" path for both.
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
_designed_call_common = dict(
    electron_spacing_meV=600.0, hole_spacing_meV=600.0,
    second_pair_addition_meV=12.126, available_pair_rate_Hz=1.0,
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
check("designed_stack_double_barrier_both_carriers_both_partitions_finite_diagnostics",
      all(math.isfinite(_r["rti_bypass_fraction_tsai_partition"])
          for _r in _designed_results.values()))
check("designed_stack_rate_or_resonance_unresolved_is_explicit_never_silently_zero_error",
      all(math.isfinite(_r["rti_e_rate_Hz"]) or "electron_resonance_unresolved" in _r["rti_failed_checks"]
          for _r in _designed_results.values())
      and all(math.isfinite(_r["rti_h_rate_Hz"]) or "hole_resonance_unresolved" in _r["rti_failed_checks"]
              for _r in _designed_results.values()))
check("designed_stack_no_hardcoded_feasible_true_and_rti_feasible_false_is_diagnosed",
      all((not _r["rti_feasible"]) or (len(_r["rti_failed_checks"]) == 0)
          for _r in _designed_results.values())
      and any(_r["rti_failed_checks"] for _r in _designed_results.values()))

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
