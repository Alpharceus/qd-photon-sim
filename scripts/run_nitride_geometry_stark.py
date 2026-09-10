"""Bounded geometry and Stark diagnostic sweep for the planar InGaN model.

Piece 6 (nitride round 2). Every predicted plotted point/cell traces back to
an evaluate() row through row_id; no numbers in results.md or the figures are
hand-entered. Geometry/orientation/shape axes and the illustrative
[-12,-8] meV/V compatibility window are [A], not measurements. Wang et al.,
Sci. Rep. 7, 12089 (2017) and Zhang et al., APL 108, 153102 (2016) are
retained only as non-gating literature guides (docs/nitride_geometry_stark_
contract.md); model predictions are never fitted to force agreement with
either.

Fix round 1 (2026-09-09, Opus review high finding): the opt-in nonpolar/QW
cards (cards/nitride-nonpolar-*.yaml, cards/nitride-qw-fluctuation-*.yaml)
were minimal provenance templates in round 2 piece 5's first pass, but piece
5's fix round (commits 2db3d60/999d969) rebuilt them as COMPLETE
one-variable copies of the c-plane parent headline cards -- every leaf the
parent carries (thermal stack, FilterBlock, aperture, drive.*, ret.*,
nitride.cavity/tau_rad0_ns/k_nr_ns/background_tau_ns) is now present
verbatim, with only the intended orientation/QW leaves changed. Every row in
this sweep therefore loads and evaluates the ACTUAL intended card file
directly (``_card_name`` below is the single source of truth for which file
that is); ``card_hash``/``card_file`` on every row are the sha256/name of
the file that was actually evaluated, never a stand-in.
"""
from __future__ import annotations
import argparse, copy, csv, hashlib, itertools, json, math, os, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

import fsim_core.device as device
from fsim_core.nitride_stark import screening_compatibility, stark_derivatives, ZHANG2016_SLOPE_MEV_PER_V, ZHANG2016_COMPARISON

REG = ("rectangular", "deterministic_pair")
ORI = ("c_plane", "semipolar_11_22", "m_plane", "a_plane")
NONPOLAR_ORI = ("m_plane", "a_plane")  # polarization_factor == 0 by construction [A]
CORE_H = (1., 2., 3., 4., 5., 7., 10.)
CORE_R = (5., 10., 15., 20., 30.)
TS = (230., 250., 273., 300.)
SCR = (0., .5, 1.)
G2, FLUX = .5, 1000.
STARK_H_FULL = (1., 2., 3., 5.)
STARK_V_FULL = tuple(round(.2 * i, 4) for i in range(17))          # 0..3.2 V, 17 pts
STARK_I_FULL = tuple(round(2.0 * 10 ** (-4 * (1 - i / 8.)), 8) for i in range(9))  # 2e-4..2 uA, 9 log pts
STARK_T_FULL = (230., 300.)
STARK_ORI = ("c_plane", "a_plane")
# Reference point shared by every one-at-a-time sensitivity axis [A, matches
# scripts/run_nitride_cavity.py's REF convention].
REF = {"height_nm": 3., "radius_nm": 10., "orientation": "c_plane", "x_in": .25,
       "screening_fraction": 0., "strain_fraction": 1., "Q": 2000., "current_uA": .02,
       "regime": "rectangular", "polarity": 1, "geometry_type": "isolated_dot",
       "shape": "disc", "top_radius_fraction": 1., "T_hs": 300.}
# Fixed reservation for the derivative-refinement convergence checks (halved
# voltage step around a representative point); a declared constant, not a
# post-hoc count, so the dry-run estimate matches the actual run exactly
# (acceptance criterion 1: no understated evaluate counts).
REFINEMENT_SPECS_QUICK = [
    dict(height_nm=3., orientation="c_plane", screening_fraction=0., regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="c_plane", screening_fraction=1., regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="a_plane", screening_fraction=0., regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="c_plane", screening_fraction=.5, regime="deterministic_pair", T_hs=300., polarity=1, target_V=1.0),
]
REFINEMENT_SPECS_FULL = REFINEMENT_SPECS_QUICK + [
    dict(height_nm=1., orientation="c_plane", screening_fraction=0., regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=5., orientation="c_plane", screening_fraction=1., regime="rectangular", T_hs=300., polarity=1, target_V=1.4),
    dict(height_nm=2., orientation="a_plane", screening_fraction=.5, regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="c_plane", screening_fraction=0., regime="rectangular", T_hs=230., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="a_plane", screening_fraction=1., regime="deterministic_pair", T_hs=300., polarity=1, target_V=1.0),
    dict(height_nm=3., orientation="c_plane", screening_fraction=1., regime="rectangular", T_hs=300., polarity=-1, target_V=1.0),
    dict(height_nm=5., orientation="c_plane", screening_fraction=0., regime="deterministic_pair", T_hs=300., polarity=1, target_V=1.6),
    dict(height_nm=1., orientation="a_plane", screening_fraction=1., regime="rectangular", T_hs=300., polarity=1, target_V=1.0),
]

# ---------------------------------------------------------------- utilities
def _safe(out):
    p = Path(out).resolve(); base = (ROOT / "out" / "nitride_geometry_stark").resolve()
    if p != base and base not in p.parents: raise ValueError("out-dir must stay under out/nitride_geometry_stark")
    return p

def _primitive(v):
    if isinstance(v, (str, int, float, bool)) or v is None: return v
    return json.dumps(v, sort_keys=True, default=str)

def _set(o, k, v):
    if isinstance(o, dict): o[k] = v
    else: setattr(o, k, v)

def _finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)

def _close(a, b, rtol=1e-8, atol=1e-10):
    """NaN/inf-safe closeness: both-NaN and both-(same-sign)-inf compare equal."""
    fa, fb = _finite(a), _finite(b)
    if not fa and not fb:
        if isinstance(a, float) and isinstance(b, float):
            if math.isnan(a) and math.isnan(b): return True
            if math.isinf(a) and math.isinf(b) and (a > 0) == (b > 0): return True
        return a == b
    if fa != fb: return False
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)

def eligible(valid, g2, flux):
    """Flux-floor eligibility gate, matching scripts/run_nitride_cavity.py's
    `eligible` EXACTLY (see run_nitride_cavity.py:61): valid + finite g2 +
    finite flux >= FLUX. g2<G2 is deliberately NOT folded in here (Opus
    fix-round finding: folding it in made this column incomparable with
    round 1) -- that gate lives in optical_pass below."""
    return bool(valid and _finite(g2) and _finite(flux) and flux >= FLUX)

def optical_pass(elig, g2, regime, one_pair_valid):
    """g2<G2 plus, for the SET regime, one_pair_valid -- matching
    scripts/run_nitride_cavity.py's optical_pass EXACTLY (run_nitride_cavity.
    py:69: ``r["eligible"] and r["g2"]<G2 and (reg!="deterministic_pair" or
    r["one_pair_valid"])``)."""
    return bool(elig and _finite(g2) and g2 < G2 and (regime != "deterministic_pair" or bool(one_pair_valid)))

def hardware_qualified(opt_pass, regime, hardware_feasible):
    """optical_pass plus, for the SET regime, hardware_feasible (=
    set_feasible AND pair_supply_possible) -- matching scripts/
    run_nitride_cavity.py's device_pass EXACTLY (run_nitride_cavity.py:70).
    Pure/importable so the verifier can unit-test the exact boundaries."""
    return bool(opt_pass and (regime != "deterministic_pair" or bool(hardware_feasible)))

# -------------------------------------------------------------- card + design
def _card_name(orientation, reg, geometry_type="isolated_dot"):
    """Pure filename resolution (no I/O) -- the single source of truth for
    which card a given (orientation, regime, geometry_type) evaluates, used
    both to load the card and to compute/report its hash."""
    if geometry_type == "qw_fluctuation":
        return "nitride-qw-fluctuation-set-design.yaml" if reg == "deterministic_pair" else "nitride-qw-fluctuation-pulse-design.yaml"
    if orientation == "a_plane":
        return "nitride-nonpolar-set-design.yaml" if reg == "deterministic_pair" else "nitride-nonpolar-pulse-design.yaml"
    return "nitride-cavity-set-design.yaml" if reg == "deterministic_pair" else "nitride-cavity-pulse-design.yaml"

_CARD_CACHE = {}
def _card(orientation, reg, geometry_type="isolated_dot"):
    """Load and evaluate the ACTUAL intended card file (Opus fix-round high
    finding). Parsed once per filename and deep-copied per row (spec: 'load
    each card once per (card, regime) and reuse it') -- a deepcopy of the
    already-parsed dataclass tree is far cheaper than re-reading and
    re-validating the YAML from disk on every row."""
    name = _card_name(orientation, reg, geometry_type)
    if name not in _CARD_CACHE:
        _CARD_CACHE[name] = device.DeviceDesign.load(ROOT / "cards" / name)
    return copy.deepcopy(_CARD_CACHE[name]), name

_DOT_KEYS = ("height_nm", "radius_nm", "x_in", "screening_fraction", "strain_fraction",
             "orientation", "polarization_factor", "shape", "top_radius_fraction",
             "geometry_type", "shape_height_fraction", "wl_thickness_nm")
_CAV_KEYS = ("Q", "detuning_offset_meV", "purcell_enabled")
_NIT_KEYS = ("tau_rad0_ns", "k_nr_ns", "background_tau_ns")

def _design(p):
    d, name = _card(p["orientation"], p["regime"], p.get("geometry_type", "isolated_dot"))
    dot, cav = d.nitride["dot"], d.nitride["cavity"]
    for k in _DOT_KEYS:
        if k in p: _set(dot, k, p[k])
    for k in _CAV_KEYS:
        if k in p: _set(cav, k, p[k])
    for k in _NIT_KEYS:
        if k in p: _set(d.nitride, k, p[k])
    _set(d.drive.diode, "x_in", p.get("x_in", .25))
    if p.get("geometry_type") == "qw_fluctuation":
        _set(d.drive.diode, "wl_thickness_nm", p["wl_thickness_nm"])
        _set(d.drive.diode, "d_i_nm", 24. + p["wl_thickness_nm"])
        _set(d.drive.diode, "d_active_nm", p["height_nm"])
    d.drive.cycle_loading = p["regime"]
    if "current_uA" in p: d.drive.I_uA = p["current_uA"]
    if "b_res" in p: d.drive.b_res = p["b_res"]
    if "n_dot_cm2" in p: d.drive.n_dot_cm2 = p["n_dot_cm2"]
    if "tau_cap_scales_with_density" in p: d.ret.tau_cap_scales_with_density = p["tau_cap_scales_with_density"]
    if "island_radius_nm" in p: _set(d.drive.set_params, "radius_nm", p["island_radius_nm"])
    # Fixed cavity reference (Opus fix-round high finding): every Stark row
    # (bias or current mode, either polarity) sets the SAME fixed
    # cavity_reference_V_j_V=0 so tau_rad_cavity/detuning/Fp_add are
    # comparable across the whole trace set -- not only when bias_mode was
    # junction_voltage or polarity was non-default, which left current-mode
    # polarity=+1 rows (the majority of the full run) on a moving,
    # non-comparable current-controlled reference.
    if "bias_mode" in p:
        if p["bias_mode"] == "junction_voltage":
            d.nitride["bias"] = {"mode": "junction_voltage", "field_polarity": p["polarity"],
                                  "V_j_V": p["V_j"], "T_j_K": p["T_j"], "cavity_reference_V_j_V": 0.}
        else:
            d.nitride["bias"] = {"mode": "current", "field_polarity": p["polarity"], "cavity_reference_V_j_V": 0.}
    return d, name

def _row(rid, kind, p, counter, cache):
    # Full input identity (incl. polarity/shape/temperature/sensitivity overrides)
    # so distinct settings never share a cache slot (spec: "never cache
    # different polarity/geometry/shape/temperature settings together").
    ident = json.dumps(p, sort_keys=True, separators=(",", ":"), default=str)
    hit = ident in cache
    if hit:
        s = cache[ident]
    else:
        d, card = _design(p)
        s = device.evaluate(d, T_grid=[p.get("T_j", p["T_hs"])])["scalars"]
        cache[ident] = s
        counter["evaluate_calls"] += 1
    # Card identity for hashing/provenance is a pure filename lookup (Opus
    # fix-round medium finding): no second _design()/YAML load needed just
    # to name the file that was (or, on a cache hit, previously was)
    # actually evaluated.
    card = _card_name(p["orientation"], p["regime"], p.get("geometry_type", "isolated_dot"))
    row = {"row_id": rid, "row_kind": kind, "card_file": card,
           "card_hash": hashlib.sha256((ROOT / "cards" / card).read_bytes()).hexdigest(),
           "cache_hit": hit, "cache_identity": ident, **p}
    for k, v in s.items(): row[k] = _primitive(v)
    row["E_X_eV"] = s.get("E_X_eV"); row["g2"] = s.get("g2_op"); row["signal_flux_s"] = s.get("collected_flux_pulsed_s")
    row["valid"] = bool(s.get("valid")); row["spectroscopy_valid"] = bool(s.get("spectroscopy_valid", s.get("valid")))
    row["invalid_reasons"] = json.dumps(s.get("invalid_reasons", []))
    one_pair_valid = bool(s.get("one_pair_valid", False))
    set_feasible = bool(s.get("set_feasible", False))
    pair_supply_possible = bool(s.get("pair_supply_possible", False))
    hardware_feasible = bool(set_feasible and pair_supply_possible)
    row["one_pair_valid"] = one_pair_valid
    row["set_feasible"] = set_feasible
    row["pair_supply_possible"] = pair_supply_possible
    # rectangular rows carry no SET hardware screen; leave blank (None), never
    # False, so they can never masquerade as hardware-screened SET data --
    # same convention as run_nitride_cavity.py.
    row["hardware_feasible"] = hardware_feasible if p["regime"] == "deterministic_pair" else None
    row["eligible"] = eligible(row["valid"], row["g2"], row["signal_flux_s"])
    row["optical_pass"] = optical_pass(row["eligible"], row["g2"], p["regime"], one_pair_valid)
    row["hardware_qualified"] = hardware_qualified(row["optical_pass"], p["regime"], hardware_feasible)
    return row

def _write(path, rows):
    keys = sorted({k for r in rows for k in r})
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

def _base(h, r, o, t, s, reg):
    return {"height_nm": h, "radius_nm": r, "orientation": o, "T_hs": t, "screening_fraction": s,
            "regime": reg, "x_in": .25, "Q": 2000., "current_uA": .02, "polarity": 1,
            "geometry_type": "isolated_dot", "shape": "disc", "top_radius_fraction": 1.}

# ------------------------------------------------------------- combo builders
# Pure (no evaluate() calls) so dry-run counts and the real run never drift
# apart (acceptance criterion 1).
def build_core(quick):
    hs = (1., 7.) if quick else CORE_H
    rs = (5., 30.) if quick else CORE_R
    temps = (230., 300.) if quick else TS
    return [_base(*z) for z in itertools.product(hs, rs, ORI, temps, SCR, REG)]

def build_shape(quick):
    if quick:
        sh_h, sh_r, sh_o, sh_t, sh_s = (3.,), (10.,), ("c_plane", "a_plane"), (300.,), (0.,)
    else:
        sh_h, sh_r, sh_o, sh_t, sh_s = (3., 7., 10.), (10., 20., 30.), ("c_plane", "a_plane"), (230., 300.), SCR
    out = []
    for h, r, o, t, s, reg, shape, frac in itertools.product(sh_h, sh_r, sh_o, sh_t, sh_s, REG, ("lens", "truncated_cone"), (None, 1.)):
        p = _base(h, r, o, t, s, reg)
        p.update(shape=shape, shape_height_fraction=frac, top_radius_fraction=.5 if shape == "truncated_cone" else 1.)
        out.append(p)
    return out

def build_qw(quick):
    if quick:
        ws, hs_off, rs, oo, tt, ss = (1.,), (.5,), (15.,), ("c_plane", "a_plane"), (300.,), (0.,)
    else:
        ws, hs_off, rs, oo, tt, ss = (1., 3.), (.5, 1.), (15., 20., 30.), ("c_plane", "a_plane"), (230., 300.), SCR
    out = []
    for w, dh, r, o, t, s, reg in itertools.product(ws, hs_off, rs, oo, tt, ss, REG):
        h = w + dh
        p = _base(h, r, o, t, s, reg); p.update(geometry_type="qw_fluctuation", wl_thickness_nm=w)
        out.append(p)
    return out

def build_stark(quick):
    """Returns (bias_params, current_params). Both lists carry bias_mode so
    _design routes them through nitride.bias correctly."""
    bias, current = [], []
    if quick:
        HH, tt, vv, ii = (1., 3.), (300.,), (0., .5, 1., 1.5, 2.), (.0002, .02, 2.)
    else:
        HH, tt, vv, ii = STARK_H_FULL, STARK_T_FULL, STARK_V_FULL, STARK_I_FULL
    for h, o, s, reg, t in itertools.product(HH, STARK_ORI, SCR, REG, tt):
        p0 = _base(h, 10., o, t, s, reg)
        for v in vv:
            p = dict(p0); p.update(bias_mode="junction_voltage", V_j=v, T_j=t, polarity=1); bias.append(p)
        for cur in ii:
            p = dict(p0); p.update(bias_mode="current", current_uA=cur, polarity=1); current.append(p)
    # opposite polarity at H=3 only, same remaining axes
    for o, s, reg, t in itertools.product(STARK_ORI, SCR, REG, tt):
        p0 = _base(3., 10., o, t, s, reg)
        for v in vv:
            p = dict(p0); p.update(bias_mode="junction_voltage", V_j=v, T_j=t, polarity=-1); bias.append(p)
        for cur in ii:
            p = dict(p0); p.update(bias_mode="current", current_uA=cur, polarity=-1); current.append(p)
    return bias, current

def build_sensitivities(quick):
    """One-at-a-time axes around REF; each axis holds every other input fixed
    (spec: 'Hold all other inputs fixed within each axis group'). Explicit
    null cases (nonpolar screening_fraction, m/a orientation identity) are
    NOT re-tested here -- they already fall out of the core grid, which
    the verifier checks directly against sweep.csv."""
    ts = (300.,) if quick else STARK_T_FULL
    out = []
    def add(axis, key, values, **extra):
        for v in values:
            for t in ts:
                p = dict(REF); p["T_hs"] = t; p[key] = v; p.update(extra)
                p["sensitivity_axis"] = axis; p["sensitivity_value"] = v
                out.append(p)
    add("semipolar_factor", "polarization_factor", (.1, .3), orientation="semipolar_11_22")
    add("x_in", "x_in", (.15, .4))
    add("strain_fraction", "strain_fraction", (0., .5))
    add("tau_rad0_ns", "tau_rad0_ns", (.5, 2.))
    add("k_nr_ns", "k_nr_ns", (0., 1.))
    add("background_tau_ns", "background_tau_ns", (0., 1.))
    add("b_res", "b_res", (0., .5))
    for q in (167., 500., 10000.):
        add("Q_purcell", "Q", (q,))
    add("Q_purcell", "Q", (500.,), purcell_enabled=False)
    add("detuning_offset_meV", "detuning_offset_meV", (0., 10.))
    for tcd in (False, True):
        add("tau_cap_density_convention", "tau_cap_scales_with_density", (tcd,), n_dot_cm2=1e9)
    for rad in ((.5, 5.) if quick else (.5, 1., 5.)):
        p = dict(REF); p["T_hs"] = ts[-1]; p["regime"] = "deterministic_pair"
        p["island_radius_nm"] = rad; p["sensitivity_axis"] = "island_radius_nm"; p["sensitivity_value"] = rad
        out.append(p)
    return out

def build_literature_aux():
    """Wang et al., Sci. Rep. 7, 12089 (2017): uncapped AFM ~7 nm height /
    ~35 nm diameter (radius ~17.5 nm) a-plane InGaN dots, 220 K [V]; a 300 K
    counterpart is added for comparison with this model's headline
    temperature. Composition/current/Q are this model's own assumptions
    [A] -- Wang's paper is optical excitation, not an electrical SPS, and
    is explicitly not a 300 K validation (docs/nitride_geometry_stark_
    contract.md)."""
    out = []
    for t, note in ((220., "comparison_only_reported_temperature"), (300., "comparison_only_model_headline_temperature")):
        p = _base(7., 17.5, "a_plane", t, 0., "rectangular")
        p["literature_source"] = "Wang et al., Sci. Rep. 7, 12089 (2017)"
        p["literature_note"] = note
        out.append(p)
    return out

def _deshpande_row(rid, counter, cache):
    """Replay the existing Deshpande comparison card exactly, as this piece's
    own independent evaluate() call (run_nitride_cavity.py's own comparison
    row is out of scope and untouched)."""
    path = ROOT / "cards" / "nitride-deshpande2014-comparison-design.yaml"
    ident = json.dumps({"card": "deshpande2014-replay", "T_hs": 300.}, sort_keys=True)
    hit = ident in cache
    if hit:
        s = cache[ident]
    else:
        d = device.DeviceDesign.load(path)
        s = device.evaluate(d, T_grid=[300.])["scalars"]
        cache[ident] = s; counter["evaluate_calls"] += 1
    row = {"row_id": rid, "row_kind": "literature_aux", "card_file": path.name,
           "card_hash": hashlib.sha256(path.read_bytes()).hexdigest(), "cache_hit": hit,
           "cache_identity": ident, "T_hs": 300., "regime": "rectangular", "orientation": "c_plane",
           "literature_source": "Deshpande et al., APL 105, 141109 (2014)",
           "literature_note": "abstract-only [V]; CONDITIONS INCOMPLETE", "measured_g2": .29}
    for k, v in s.items(): row[k] = _primitive(v)
    row["E_X_eV"] = s.get("E_X_eV"); row["g2"] = s.get("g2_op"); row["signal_flux_s"] = s.get("collected_flux_pulsed_s")
    row["valid"] = bool(s.get("valid")); row["spectroscopy_valid"] = bool(s.get("spectroscopy_valid", s.get("valid")))
    row["invalid_reasons"] = json.dumps(s.get("invalid_reasons", []))
    row["one_pair_valid"] = bool(s.get("one_pair_valid", False))
    row["set_feasible"] = bool(s.get("set_feasible", False))
    row["pair_supply_possible"] = bool(s.get("pair_supply_possible", False))
    row["hardware_feasible"] = None  # regime=="rectangular" here
    row["eligible"] = eligible(row["valid"], row["g2"], row["signal_flux_s"])
    row["optical_pass"] = optical_pass(row["eligible"], row["g2"], "rectangular", row["one_pair_valid"])
    row["hardware_qualified"] = False  # literature replay, never a headline hardware qualification
    return row

# ----------------------------------------------------------- Stark derivative wiring
def attach_derivatives(pool, key="V_j"):
    groups = {}
    for r in pool:
        groups.setdefault(tuple(r.get(k) for k in ("height_nm", "orientation", "screening_fraction", "regime", "T_hs", "polarity", "bias_mode")), []).append(r)
    for rr in groups.values():
        rr.sort(key=lambda x: float(x.get(key)) if isinstance(x.get(key), (int, float)) else -1e18)
        ds = stark_derivatives(rr, voltage_key=key)
        for r, dv in zip(rr, ds):
            r.update({k: dv[k] for k in ("dE_X_dV_meV_per_V", "derivative_valid", "derivative_row_ids", "derivative_policy")})

# ----------------------------------------------------------- convergence refinement
def run_refinement_checks(specs, bias_rows, counter, cache, rid_start, half_step):
    """half_step is HALF of the actual V_j trace step used to build
    bias_rows (0.1 V for the full 0.2 V grid, 0.25 V for the quick 0.5 V
    grid) -- a declared, mode-correct value, not a fixed constant that
    happens to coincide with the full grid only (Opus fix-round finding:
    'refinement probes at half the trace step'). Each probe pair must also
    share the SAME flat_band marker as the anchor point; a probe that
    crosses the flat-band kink invalidates the derivative instead of being
    silently averaged through it."""
    checks = []; rid = rid_start; probe_rows = []
    for spec in specs:
        target_v = spec["target_V"]
        group = [r for r in bias_rows if all(r.get(k) == spec[k] for k in
                 ("height_nm", "orientation", "screening_fraction", "regime", "T_hs", "polarity"))
                 and r.get("bias_mode") == "junction_voltage"]
        group = [r for r in group if r.get("derivative_valid") is True]
        if not group:
            checks.append({"spec": spec, "status": "no_valid_derivative_point", "converged": False}); continue
        anchor = min(group, key=lambda r: abs(float(r["V_j"]) - target_v))
        v0 = float(anchor["V_j"]); coarse = float(anchor["dE_X_dV_meV_per_V"])
        anchor_flat_band = anchor.get("flat_band")
        p0 = _base(spec["height_nm"], 10., spec["orientation"], spec["T_hs"], spec["screening_fraction"], spec["regime"])
        probes = {}
        for dv, tag in ((-half_step, "lo"), (half_step, "hi")):
            v = round(v0 + dv, 6)
            p = dict(p0); p.update(bias_mode="junction_voltage", V_j=v, T_j=spec["T_hs"], polarity=spec["polarity"])
            rid += 1
            row = _row(f"RF{rid:05d}", "refinement_probe", p, counter, cache)
            probes[tag] = row; probe_rows.append(row)
        e_lo, e_hi = probes["lo"]["E_X_eV"], probes["hi"]["E_X_eV"]
        fb_lo, fb_hi = probes["lo"].get("flat_band"), probes["hi"].get("flat_band")
        probe_ids = [probes["lo"]["row_id"], probes["hi"]["row_id"]]
        if not (_finite(e_lo) and _finite(e_hi)):
            checks.append({"spec": spec, "status": "probe_invalid", "converged": False,
                            "probe_row_ids": probe_ids}); continue
        if fb_lo != anchor_flat_band or fb_hi != anchor_flat_band:
            checks.append({"spec": spec, "status": "flat_band_marker_mismatch", "converged": False,
                            "probe_row_ids": probe_ids, "anchor_flat_band": anchor_flat_band,
                            "probe_flat_band": [fb_lo, fb_hi]}); continue
        fine = 1000. * (e_hi - e_lo) / (2. * half_step)
        abs_diff = abs(fine - coarse)
        rel_diff = abs_diff / abs(coarse) if coarse else float("inf")
        converged = abs_diff <= 0.5 or rel_diff <= 0.05
        checks.append({"spec": spec, "anchor_row_id": anchor["row_id"], "anchor_V_j": v0,
                        "coarse_slope_meV_per_V": coarse, "fine_slope_meV_per_V": fine,
                        "abs_diff_meV_per_V": abs_diff, "rel_diff": rel_diff,
                        "criterion": "abs<=0.5meV/V or rel<=5%", "converged": bool(converged),
                        "half_step_V": half_step, "flat_band_marker": anchor_flat_band,
                        "probe_row_ids": probe_ids})
    return checks, rid, probe_rows

# ----------------------------------------------------------- screening compatibility
def build_compatibility(bias_rows, slope_range, bias_window, slope_user_supplied, rid_start):
    compat = screening_compatibility(bias_rows, slope_range_meV_per_V=slope_range, voltage_window_V=tuple(bias_window))
    out = []; rid = rid_start
    for x in compat:
        rid += 1
        row = {"row_id": f"CT{rid:05d}",
               "slope_interval_kind": "user_supplied" if slope_user_supplied else "illustrative_not_measured",
               "zhang_guide_meV_per_V": ZHANG2016_SLOPE_MEV_PER_V,
               "source": ZHANG2016_COMPARISON["citation"] + " [V] non-gating"}
        gk = dict(x["group"])
        orientation = gk.get("orientation")
        if orientation in NONPOLAR_ORI and x["identification_status"] == "compatible":
            # Explicit expected null case: polarization_factor==0 makes
            # screening_fraction physically inert for m/a-plane -- the
            # fitted slope there tests only whether the model predicts a
            # near-flat nonpolar Stark response (falsifiable), never a
            # screening fraction. Override ONLY when the measured interval
            # actually ACCEPTS the shared curve (identification_status was
            # "compatible" before this override) -- an a-plane/m-plane slope
            # the window already REJECTS stays "incompatible" (Opus
            # fix-round finding: the previous unconditional override
            # relabelled genuinely out-of-window nonpolar slopes as
            # "unresolved" instead of excluding them).
            x = dict(x); x["compatible"] = None; x["identification_status"] = "screening_unidentifiable"
            x["nuisance_sensitivity"] = "nonpolar_orientation_zero_polarization_factor"
        for k, v in x.items(): row[k] = _primitive(v)
        out.append(row)
    return out, rid

# -------------------------------------------------------------------- plots
def _finite_num(v):
    try: return math.isfinite(float(v))
    except (TypeError, ValueError): return False

def _leg(ax, **kw): return ax.legend(fontsize=7, loc="best", **kw)

MAX_TRACES_PER_PANEL = 8  # spec: "split into panels if > 8 traces"

def _panel_plot(out, name, entries, xlabel, ylabel, suptitle, xlog=False, ylog=False):
    """entries: ordered [(label, xv, yv, row_ids, style), ...], ONE trace per
    fully-specified group -- never merged (Opus fix-round high finding).
    Splits into ceil(n/8) panels so no group is silently dropped (the
    previous [:14] truncation), each panel's legend placed OUTSIDE the axes.
    Returns the plot-construction contract {label: {row_ids, x, y}}."""
    import matplotlib.pyplot as plt
    contract = {}
    n = len(entries)
    n_panels = max(1, math.ceil(n / MAX_TRACES_PER_PANEL))
    ncols = min(3, n_panels); nrows = math.ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.6 * ncols, 4.6 * nrows), squeeze=False)
    for i in range(nrows * ncols):
        ax = axes[i // ncols][i % ncols]
        if i >= n_panels:
            ax.set_visible(False); continue
        chunk = entries[i * MAX_TRACES_PER_PANEL:(i + 1) * MAX_TRACES_PER_PANEL]
        for label, xv, yv, rids, style in chunk:
            ax.plot(xv, yv, style, ms=3, lw=1.2, label=label)
            contract[label] = {"row_ids": rids, "x": xv, "y": yv}
        if xlog: ax.set_xscale("log")
        if ylog: ax.set_yscale("log")
        ax.set(xlabel=xlabel, ylabel=ylabel)
        if n_panels > 1: ax.set_title(f"panel {i + 1}/{n_panels}", fontsize=8)
        if chunk: ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)
    fig.suptitle(suptitle, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / name, dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract

def plot_axis_response(out, name, rows, x, y, title, group_keys=("orientation", "screening_fraction", "regime")):
    groups = {}
    for r in rows:
        if _finite_num(r.get(x)) and _finite_num(r.get(y)):
            groups.setdefault(tuple(r.get(k) for k in group_keys), []).append(r)
    entries = []
    for key, rr in sorted(groups.items(), key=lambda kv: str(kv[0])):
        rr = sorted(rr, key=lambda z: float(z[x]))
        lab = "/".join(str(k) for k in key)
        xv = [float(z[x]) for z in rr]; yv = [float(z[y]) for z in rr]
        entries.append((lab, xv, yv, [z["row_id"] for z in rr], "o-"))
    return _panel_plot(out, name, entries, x, y, title)

def plot_stark(out, name, rows, x, y, title, ylog=False, zhang_anchor=False):
    """E_X/tau/overlap vs V_j or current. Grouped by the FULL coordinate set
    (orientation, height, regime, T_hs, polarity, screening) -- one trace
    per group, never merged (Opus fix-round high finding: the previous
    (orientation, screening, height, regime) grouping silently averaged
    opposite-polarity and multi-T_hs rows into one sawtooth trace). Screening
    0/0.5/1 use the same line-style/legend convention as the rest of the
    nitride pieces. Current traces (x=='current_uA') use a log x-axis."""
    label_for_s = {0.: "unscreened lower bound", .5: "midpoint 0.5", 1.: "screened upper bound"}
    style_for_s = {0.: "-o", .5: "--o", 1.: ":o"}
    groups = {}
    for r in rows:
        if _finite_num(r.get(x)) and _finite_num(r.get(y)):
            key = (r.get("orientation"), r.get("height_nm"), r.get("regime"), r.get("T_hs"),
                   r.get("polarity"), r.get("screening_fraction"))
            groups.setdefault(key, []).append(r)
    entries = []; zhang_entry = None
    for key, rr in sorted(groups.items(), key=lambda kv: str(kv[0])):
        o, h, reg, t, pol, s = key
        rr = sorted(rr, key=lambda z: float(z[x]))
        xv = [float(z[x]) for z in rr]; yv = [float(z[y]) for z in rr]
        lab = f"{o} h={h:g}nm {reg} T={t:g}K pol={pol:g} ({label_for_s.get(s, s)})"
        entries.append((lab, xv, yv, [z["row_id"] for z in rr], style_for_s.get(s, "-.o")))
        if zhang_anchor and zhang_entry is None and o == "c_plane" and s == 0. and reg == "rectangular" and pol == 1 and len(xv) >= 2:
            # Zhang guide anchored to this actual model row (not an invented
            # dataset); no row_ids of its own (it is a drawn annotation, not
            # a saved evaluate() trace).
            x0, y0 = xv[0], yv[0]
            zhang_entry = ("Zhang 2016 -10 meV/V guide (non-gating)",
                            [x0, x0 + 2.0], [y0, y0 + 2.0 * ZHANG2016_SLOPE_MEV_PER_V / 1000.], [], "k--")
    if zhang_entry is not None: entries.append(zhang_entry)
    return _panel_plot(out, name, entries, x, y, title, xlog=(x == "current_uA"), ylog=ylog)

def plot_envelope(out, name, core_rows, title):
    """Columns = orientation (c-plane AND a-plane -- Opus fix-round medium
    finding: both envelope figures previously filtered to c_plane only, no
    a-plane panel); rows = unscreened lower bound / screened upper bound.
    Cell value is the optical-pass fraction, matching run_nitride_cavity.
    py's envelope convention exactly (its `optical_pass`, not the flux-only
    `eligible` gate)."""
    import matplotlib.pyplot as plt
    orientations = [o for o in ("c_plane", "a_plane") if any(r["orientation"] == o for r in core_rows)] or ["c_plane"]
    fig, axes = plt.subplots(2, len(orientations), figsize=(6.2 * len(orientations), 8.6), squeeze=False)
    contract = {}
    hs = sorted({r["height_nm"] for r in core_rows}); ts = sorted({r["T_hs"] for r in core_rows})
    for col, o in enumerate(orientations):
        for row_i, (s, lab) in enumerate(((0., "unscreened lower bound"), (1., "screened upper bound"))):
            ax = axes[row_i][col]
            mat = []; cellmap = {}
            for t in ts:
                rowvals = []
                for h in hs:
                    q = [r for r in core_rows if r["orientation"] == o and r["height_nm"] == h and r["T_hs"] == t and r["screening_fraction"] == s]
                    rowvals.append(sum(r["optical_pass"] for r in q) / len(q) if q else float("nan"))
                    cellmap[f"h{h:g}_T{t:g}_s{s:g}"] = [r["row_id"] for r in q]
                mat.append(rowvals)
            im = ax.imshow(mat, vmin=0, vmax=1, aspect="auto", origin="lower", cmap="RdYlGn")
            fig.colorbar(im, ax=ax, label="optical-pass fraction")
            ax.set(xticks=range(len(hs)), xticklabels=[f"{v:g}" for v in hs],
                   yticks=range(len(ts)), yticklabels=[f"{v:g}" for v in ts],
                   xlabel="height nm", ylabel="T_hs K", title=f"{o} {lab}\n{title}")
            contract[f"{o} {lab}"] = cellmap
    fig.tight_layout(); fig.savefig(out / name, dpi=125); plt.close(fig)
    return contract

def plot_pulse_vs_set(out, core_rows):
    import matplotlib.pyplot as plt
    UNDERFLOW_FLOOR = 1e-9  # display floor for the LOG FLUX AXIS ONLY; a bar
    # drawn at this floor is annotated "underflow" (spec: "distinguish
    # numerical underflow from physical zero" -- Opus fix-round finding: the
    # previous max(fx,1e-12) silently drew a 3.6e-15/s bar as if it were a
    # real ~1e-12/s datapoint).
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.4)); contract = {}
    ts = sorted({r["T_hs"] for r in core_rows}); wdt = .18
    xs = range(len(ts))
    for i, (reg, s, lab, off) in enumerate([("rectangular", 0., "pulse unscreened", -1.5),
                                             ("rectangular", 1., "pulse screened", -.5),
                                             ("deterministic_pair", 0., "SET unscreened", .5),
                                             ("deterministic_pair", 1., "SET screened", 1.5)]):
        vals_g2, vals_fx, ids, underflow = [], [], [], []
        for t in ts:
            rr = [r for r in core_rows if r["regime"] == reg and r["screening_fraction"] == s and r["T_hs"] == t
                  and r["orientation"] == "c_plane" and r["height_nm"] == 1 and r["radius_nm"] == 5]
            g2v = rr[0]["g2"] if rr and _finite_num(rr[0]["g2"]) else float("nan")
            fx = rr[0]["signal_flux_s"] if rr and _finite_num(rr[0]["signal_flux_s"]) else float("nan")
            vals_g2.append(g2v)
            if _finite_num(fx) and fx < UNDERFLOW_FLOOR:
                vals_fx.append(UNDERFLOW_FLOOR); underflow.append(True)
            else:
                vals_fx.append(fx if _finite_num(fx) else float("nan")); underflow.append(False)
            ids += [r["row_id"] for r in rr]
        xpos = [xv + off * wdt for xv in xs]
        ax1.bar(xpos, vals_g2, width=wdt, label=lab)
        bars = ax2.bar(xpos, vals_fx, width=wdt, label=lab)
        for bar, flag in zip(bars, underflow):
            if flag:
                ax2.text(bar.get_x() + bar.get_width() / 2, UNDERFLOW_FLOOR, "underflow",
                          rotation=90, fontsize=5, va="bottom", ha="center")
        contract[lab] = ids
    ax1.axhline(G2, color="k", lw=.8); ax1.set(xticks=list(xs), xticklabels=[f"{t:g}" for t in ts], xlabel="T_hs K", ylabel="g2", title="pulse vs SET: g2")
    ax2.axhline(FLUX, color="gray", lw=.8, ls="--"); ax2.set_yscale("log"); ax2.set(xticks=list(xs), xticklabels=[f"{t:g}" for t in ts], xlabel="T_hs K", ylabel="flux /s (log)", title="pulse vs SET: flux\n('underflow' bars are display-floored, not a real value)")
    _leg(ax1); _leg(ax2); fig.tight_layout(); fig.savefig(out / "pulse_vs_set.png", dpi=125); plt.close(fig)
    return contract

def plot_set_island(out, sens_rows):
    import matplotlib.pyplot as plt
    z = [r for r in sens_rows if r.get("sensitivity_axis") == "island_radius_nm"]
    fig, ax = plt.subplots(figsize=(6.5, 4.2)); contract = {}
    rr = sorted(z, key=lambda r: float(r["sensitivity_value"]))
    xs = [float(r["sensitivity_value"]) for r in rr]
    ys = [r.get("set_EC_over_kT") for r in rr]
    ax.plot(xs, [float(y) if _finite_num(y) else float("nan") for y in ys], "o-", label="E_C / kT")
    margin = [10.] * len(xs)
    ax.plot(xs, margin, "k--", label="ec_margin threshold (10 kT)")
    ax.set(xlabel="assumed island radius nm", ylabel="E_C / kT", title="SET island charging-energy feasibility screen")
    _leg(ax); fig.tight_layout(); fig.savefig(out / "set_feasibility.png", dpi=125); plt.close(fig)
    contract["E_C_over_kT"] = {"row_ids": [r["row_id"] for r in rr], "x": xs, "y": [float(y) if _finite_num(y) else float("nan") for y in ys]}
    return contract

# ---------------------------------------------------------------- results.md
def _results_md(core, geo, bias, current, comp, refinement_checks, complete, runtime_s, evaluate_calls):
    lines = ["# Nitride geometry and Stark sweep results", "",
             "Model-only planar predictions; every number below is an evaluate() output "
             "or a traceable reduction of one (row_id in manifest.json's plot_row_mapping / "
             "the CSVs). Unscreened (screening_fraction=0) is labelled the conventional "
             "lower bound and screened (screening_fraction=1) the upper bound -- "
             "conventional polarization-screening SCENARIOS, not a rigorous ordered bound "
             "on every metric under every polarity; crossings are not suppressed. "
             "`eligible` is the flux-floor gate only (matches scripts/run_nitride_cavity.py's "
             "`eligible` exactly); `paired_optical_pass` additionally requires g2<0.5 and, "
             "for the SET regime, one_pair_valid; `hardware_qualified` further requires "
             "set_feasible AND pair_supply_possible for the SET regime.", ""]
    for family in ("c_plane", "a_plane"):
        for reg in REG:
            for s, label in ((0., "unscreened_lower"), (1., "screened_upper")):
                rs = [r for r in core if r["orientation"] == family and r["regime"] == reg and r["screening_fraction"] == s]
                elig = [r for r in rs if r["eligible"]]; opt = [r for r in rs if r["optical_pass"]]; hw = [r for r in rs if r["hardware_qualified"]]
                if reg == "deterministic_pair":
                    ideal = "pass_hardware_qualified" if hw else ("pass_hardware_infeasible" if opt else "no_idealized_pass")
                else:
                    ideal = "pass" if opt else "no_idealized_pass"
                lines.append(
                    "VERDICT: idealized_status=%s family=%s regime=%s screening=%s complete=%s "
                    "eligible=%d paired_optical_pass=%d hardware_qualified=%d coverage=%d/%d "
                    "invalid=%d flux_floor=1000/s" % (
                        ideal, family, reg, label, complete, len(elig), len(opt), len(hw),
                        len(rs), len(rs), sum(not r["valid"] for r in rs)))
    lines += ["",
        "## Supplementary orientation and geometry families (not headline coverage)", "",
        "| family | regime | screening | paired_optical_pass | hardware_qualified | total | invalid |",
        "|---|---|---|---|---|---|---|"]
    for family in ("semipolar_11_22", "m_plane"):
        for reg in REG:
            for s in SCR:
                rs = [r for r in core if r["orientation"] == family and r["regime"] == reg and r["screening_fraction"] == s]
                if not rs: continue
                opt = [r for r in rs if r["optical_pass"]]; hw = [r for r in rs if r["hardware_qualified"]]
                lines.append(f"| {family} | {reg} | {s:g} | {len(opt)} | {len(hw)} | {len(rs)} | {sum(not r['valid'] for r in rs)} |")
    # Maximum flux / eligible-g2 statistics (Opus fix-round finding: results.md
    # omitted both, along with reservoir choice and the round-1 coverage delta).
    all_rows = core + geo + bias + current
    flux_pairs = [(float(r["signal_flux_s"]), r["row_id"]) for r in all_rows
                  if r.get("valid") and _finite_num(r.get("signal_flux_s"))]
    opt_core = [r for r in core if r.get("optical_pass")]
    g2_opt = sorted(float(r["g2"]) for r in opt_core if _finite_num(r.get("g2")))
    def _median(xs):
        n = len(xs)
        if not n: return None
        return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])
    lines += ["", "## Maximum flux and eligible-g2 statistics", ""]
    if flux_pairs:
        mf, mf_rid = max(flux_pairs, key=lambda z: z[0])
        lines.append(f"Maximum valid signal_flux_s across all core/shape/QW/Stark rows: {mf:.6g}/s (row {mf_rid}).")
    else:
        lines.append("No valid rows with finite signal_flux_s.")
    if g2_opt:
        lines.append(f"Core-grid paired-optical-pass rows (g2<{G2:g}, flux>={FLUX:g}/s, plus "
                      f"one_pair_valid for SET): {len(opt_core)} of {len(core)}; "
                      f"g2_min={g2_opt[0]:.6g}, g2_median={_median(g2_opt):.6g}.")
    else:
        lines.append(f"No core-grid rows pass the paired-optical-pass gate (g2<{G2:g}, flux>={FLUX:g}/s).")
    lines.append("")
    # Per-card reservoir choice (device.py's reservoir_kind, one sample row
    # per distinct (card_file, orientation, geometry_type)).
    lines += ["## Per-card reservoir choice", "", "| card_file | orientation | geometry_type | reservoir_kind (sample) |", "|---|---|---|---|"]
    seen_cards = set()
    for r in core + geo:
        key = (r.get("card_file"), r.get("orientation"), r.get("geometry_type"))
        if key in seen_cards or "reservoir_kind" not in r: continue
        seen_cards.add(key)
        lines.append(f"| {r.get('card_file')} | {r.get('orientation')} | {r.get('geometry_type')} | {r.get('reservoir_kind')} |")
    lines.append("")
    lines += ["## Coverage and geometry change versus round 1 (scripts/run_nitride_cavity.py)", "",
        "Round 1's headline grid was c-plane only, radius fixed at 10 nm (height in {1,2,3,4,5} nm, "
        "x_in in {0.15,0.25,0.4}, Q in {500,2000,10000}, current in {0.002,0.02,0.2} uA), no "
        "orientation axis, and no Stark bias/current traces. This piece adds a "
        f"{len(ORI)}-orientation x {len(CORE_R)}-radius x {len(CORE_H)}-height x {len(TS)}-T_hs x "
        f"{len(SCR)}-screening x {len(REG)}-regime core grid ({len(core)} rows this run), "
        f"lens/truncated_cone shape ({len([r for r in geo if r['row_kind'] == 'shape'])} rows) and "
        f"QW-fluctuation ({len([r for r in geo if r['row_kind'] == 'qw'])} rows) geometry "
        f"supplements, and bias/current Stark diagnostic traces ({len(bias)}+{len(current)} rows) "
        "with derivative-refinement convergence checks and screening-compatibility fits. Round 1's "
        "own results.md and artifacts under out/nitride_cavity/ are unchanged and are not "
        "reinterpreted here (out of scope for this piece).", ""]
    shape_rows = [r for r in geo if r["row_kind"] == "shape"]; qw_rows = [r for r in geo if r["row_kind"] == "qw"]
    lines += ["", "## Shape supplement (mapping sensitivity, not a demonstrated shape accuracy)",
              f"{len(shape_rows)} rows (lens / truncated_cone, native and full-height variants). "
              "The numerical spread across shape variants at matched height/radius is a mapping-"
              "sensitivity diagnostic; systematic shape accuracy remains unquantified (no 3-D solver).",
              "", "## QW-fluctuation supplement",
              f"{len(qw_rows)} rows; diode.wl_thickness_nm/x_in and dot.wl_thickness_nm/x_in are matched "
              "per card contract; d_i_nm=24+w, d_active_nm=height [A]. No lens-in-QW cases were run "
              "(unsupported combination).", ""]
    wang = [r for r in geo if r["row_kind"] == "literature_aux" and "Wang" in str(r.get("literature_source", ""))]
    lines += ["## Literature comparisons (non-gating)", "",
              "Wang et al., Sci. Rep. 7, 12089 (2017): uncapped a-plane AFM ~7 nm/~35 nm dots, "
              "2.54 eV at 220 K, 19.0+/-0.4 meV linewidth, raw/corrected g2 0.47/0.21, optical "
              "excitation (76 MHz, 1 ps, 800 nm two-photon) -- NOT an electrical SPS or a 300 K "
              "validation. Model rows at the same nominal geometry (comparison-only, excluded "
              "from headline coverage):", ""]
    for r in wang:
        lines.append(f"- row {r['row_id']}: T_hs={r['T_hs']:g} K, E_X_eV={r.get('E_X_eV')}, g2={r.get('g2')}, valid={r['valid']} ({r.get('literature_note')})")
    desh = [r for r in geo if r["row_kind"] == "literature_aux" and "Deshpande" in str(r.get("literature_source", ""))]
    if desh:
        r = desh[0]
        lines.append(f"- Deshpande et al., APL 105, 141109 (2014) [V abstract-only; CONDITIONS INCOMPLETE], replay row {r['row_id']}: measured g2=0.29, model g2={r.get('g2')}.")
    lines += ["",
              f"Zhang et al., APL 108, 153102 (2016), Fig. 5: {ZHANG2016_SLOPE_MEV_PER_V:g} meV/V below 2 V, "
              "non-gating (different device, excitation and unspecified temperature); shown only as a "
              "labelled guide anchored to a model row on the Stark figures, never a fitted target.", ""]
    lines += ["## Screening-compatibility table", "",
              "The [-12,-8] meV/V window used when --slope-range is not supplied is EXPLICITLY "
              "illustrative around Zhang's approximate value, never a measured interval; the same "
              "--bias-window applies to any user-supplied slope. Nonpolar (m-plane/a-plane) rows "
              "are marked screening_unidentifiable: polarization_factor=0 makes screening_fraction "
              "physically inert for those orientations, so no compatibility test there can resolve "
              "a screening fraction -- this is a structural degeneracy, not a measurement result.", ""]
    n_unident = sum(1 for r in comp if r.get("identification_status") == "screening_unidentifiable")
    n_compat = sum(1 for r in comp if r.get("compatible") in (True, "True"))
    n_incompat = sum(1 for r in comp if r.get("identification_status") == "incompatible")
    lines.append(f"{len(comp)} compatibility rows: {n_compat} compatible, {n_incompat} incompatible, {n_unident} screening_unidentifiable (nonpolar or degenerate).")
    lines += ["", "## Derivative-refinement convergence checks", "",
              "| check | anchor V_j | coarse meV/V | fine (half-step) meV/V | abs diff | rel diff | converged |",
              "|---|---|---|---|---|---|---|"]
    for i, c in enumerate(refinement_checks):
        if "coarse_slope_meV_per_V" not in c:
            lines.append(f"| {i} | - | - | - | - | - | {c.get('status')} |"); continue
        lines.append(f"| {i} | {c['anchor_V_j']:.3g} | {c['coarse_slope_meV_per_V']:.4g} | {c['fine_slope_meV_per_V']:.4g} | {c['abs_diff_meV_per_V']:.4g} | {c['rel_diff']:.4g} | {c['converged']} |")
    lines += ["", "## Limitations",
              "Nonpolar strain, valence-band ordering and FSS are not modelled; shape mapping error "
              "is unquantified; finite-dot lateral fields, field-assisted escape and injection-"
              "dependent screening are not solved (screening_fraction is fixed along every current "
              "trace). m-plane and a-plane rows coincide under this scalar model (both "
              "polarization_factor=0) -- a model limitation, not independent evidence for either. "
              "SET hardware feasibility (island charging energy / RC bandwidth) is reported "
              "independently of the idealized optical pass and never substitutes for it. No "
              "field-free lifetime or screening fraction was ever fitted to force agreement with "
              "Wang or Zhang.",
              "", f"runtime_s={runtime_s:.1f} evaluate_calls={evaluate_calls} complete={complete}"]
    return "\n".join(lines) + "\n"

# -------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--stark", action="store_true")
    ap.add_argument("--out-dir", default=str(ROOT / "out" / "nitride_geometry_stark"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-evaluations", type=int, default=10000)
    ap.add_argument("--slope-range", nargs=2, type=float, default=None)
    ap.add_argument("--bias-window", nargs=2, type=float, default=(0., 2.))
    a = ap.parse_args(argv)

    core_p = build_core(a.quick); shape_p = build_shape(a.quick); qw_p = build_qw(a.quick)
    bias_p, current_p = build_stark(a.quick) if a.stark else ([], [])
    sens_p = build_sensitivities(a.quick); lit_p = build_literature_aux()
    refinement_specs = (REFINEMENT_SPECS_QUICK if a.quick else REFINEMENT_SPECS_FULL) if a.stark else []
    planned = (len(core_p) + len(shape_p) + len(qw_p) + len(bias_p) + len(current_p)
               + len(sens_p) + len(lit_p) + 1 + 2 * len(refinement_specs))
    if a.dry_run:
        print(json.dumps({"core_rows": len(core_p), "shape_rows": len(shape_p), "qw_rows": len(qw_p),
                           "stark_bias_rows": len(bias_p), "stark_current_rows": len(current_p),
                           "sensitivity_rows": len(sens_p), "literature_aux_rows": len(lit_p) + 1,
                           "refinement_evaluate_calls": 2 * len(refinement_specs),
                           "planned_evaluate_calls": planned, "max_evaluations": a.max_evaluations,
                           "within_cap": planned <= a.max_evaluations and planned <= 10000}, sort_keys=True))
        return 0 if planned <= a.max_evaluations and planned <= 10000 else 1
    if planned > a.max_evaluations or planned > 10000:
        raise SystemExit("planned evaluate calls (%d) exceed cap" % planned)

    out = _safe(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(out / "mplconfig")
    t0 = time.time(); cache = {}; counter = {"evaluate_calls": 0}; rid = 0

    def add_all(target, kind, plist):
        nonlocal rid
        for p in plist:
            rid += 1
            target.append(_row(f"{kind[:2].upper()}{rid:05d}", kind, p, counter, cache))

    core, geo, bias, current, sens = [], [], [], [], []
    add_all(core, "core", core_p)
    add_all(geo, "shape", shape_p)
    add_all(geo, "qw", qw_p)
    add_all(geo, "literature_aux", lit_p)
    rid += 1
    geo.append(_deshpande_row(f"LI{rid:05d}", counter, cache))
    add_all(bias, "stark_bias", bias_p)
    add_all(current, "stark_current", current_p)
    add_all(sens, "sensitivity", sens_p)

    if bias: attach_derivatives(bias, "V_j")
    if current: attach_derivatives(current, "V_j")

    # Half of the ACTUAL V_j trace step used to build `bias` in this mode
    # (0.1 V for the full 0.2 V grid, 0.25 V for the quick 0.5 V grid) --
    # spec: "refinement probes at half the trace step".
    stark_v_half_step = (0.5 if a.quick else 0.2) / 2.
    refinement_checks = []; probe_rows = []
    if refinement_specs and bias:
        refinement_checks, rid, probe_rows = run_refinement_checks(refinement_specs, bias, counter, cache, rid, stark_v_half_step)

    slope_range = tuple(a.slope_range) if a.slope_range else (-12., -8.)
    comp = []
    if bias:
        comp, rid = build_compatibility(bias, slope_range, a.bias_window, a.slope_range is not None, rid)

    _write(out / "sweep.csv", core)
    _write(out / "geometry_supplement.csv", geo)
    _write(out / "stark_bias.csv", bias)
    _write(out / "stark_current.csv", current)
    _write(out / "screening_compatibility.csv", comp)
    _write(out / "sensitivities.csv", sens)
    lit_csv = [
        {"row_id": "LC00001", "source": "Wang et al., Sci. Rep. 7, 12089 (2017)", "kind": "comparison_only",
         "temperature_K": 220, "height_nm": 7, "diameter_nm": 35, "note": "[V] not electrical SPS, not 300 K validation",
         "model_row_ids": json.dumps([r["row_id"] for r in geo if r["row_kind"] == "literature_aux" and "Wang" in str(r.get("literature_source", ""))])},
        {"row_id": "LC00002", "source": "Wang et al., APL 111, 053101 (2017)", "kind": "informational_only",
         "note": "[V] FSS 2-12 meV at 200 K (16 dots); omitted by this single-band model", "model_row_ids": "[]"},
        {"row_id": "LC00003", "source": ZHANG2016_COMPARISON["citation"], "kind": "non_gating_slope_guide",
         "slope_meV_per_V": ZHANG2016_SLOPE_MEV_PER_V, "note": "[V] Fig. 5, below 2 V; non-gating comparison only",
         "model_row_ids": json.dumps([r["row_id"] for r in comp])},
        {"row_id": "LC00004", "source": "Deshpande et al., APL 105, 141109 (2014)", "kind": "conditions_incomplete_replay",
         "note": "[V abstract-only] measured g2=0.29",
         "model_row_ids": json.dumps([r["row_id"] for r in geo if r["row_kind"] == "literature_aux" and "Deshpande" in str(r.get("literature_source", ""))])},
    ]
    _write(out / "literature_comparisons.csv", lit_csv)

    plots = {}
    # A fixed representative slice (radius=5nm, height=1nm) is used for the
    # "other axis" of these response plots -- both values are guaranteed
    # present in the reduced --quick core grid (H=[1,7], R=[5,30]) as well
    # as the full grid, so quick and full runs both produce populated
    # (non-empty) figures from the same slice logic.
    # height/radius response: group_keys include T_hs (Opus fix-round medium
    # finding: these previously selected T_hs in (230,300) but omitted it
    # from group_keys, mixing both temperatures into one curve at duplicated
    # x); orientation_flux keeps its default group_keys (T_hs is the x-axis
    # there, not a grouping key) -- panel-splitting (in plot_axis_response/
    # _panel_plot) now shows every orientation/screening/regime group
    # instead of truncating to the first 14 by sort key.
    plots["height_response.png"] = plot_axis_response(out, "height_response.png", [r for r in core if r["T_hs"] in (230., 300.) and r["radius_nm"] == 5.], "height_nm", "E_X_eV", "height response (r=5nm)", group_keys=("orientation", "screening_fraction", "regime", "T_hs"))
    plots["radius_response.png"] = plot_axis_response(out, "radius_response.png", [r for r in core if r["T_hs"] in (230., 300.) and r["height_nm"] == 1.], "radius_nm", "E_X_eV", "radius response (h=1nm)", group_keys=("orientation", "screening_fraction", "regime", "T_hs"))
    plots["temperature_response.png"] = plot_axis_response(out, "temperature_response.png", [r for r in core if r["radius_nm"] == 5. and r["height_nm"] == 1.], "T_hs", "g2", "temperature response (g2)")
    plots["orientation_flux.png"] = plot_axis_response(out, "orientation_flux.png", [r for r in core if r["radius_nm"] == 5. and r["height_nm"] == 1.], "T_hs", "signal_flux_s", "orientation flux comparison")
    plots["shape_comparison.png"] = plot_axis_response(out, "shape_comparison.png", [r for r in geo if r["row_kind"] == "shape"], "height_nm", "E_X_eV", "shape mapping sensitivity", group_keys=("shape", "orientation", "screening_fraction"))
    plots["qw_response.png"] = plot_axis_response(out, "qw_response.png", [r for r in geo if r["row_kind"] == "qw"], "wl_thickness_nm", "E_X_eV", "QW thickness response", group_keys=("orientation", "screening_fraction"))
    plots["pulse_vs_set.png"] = plot_pulse_vs_set(out, core)
    # Envelopes cover c-plane AND a-plane (Opus fix-round medium finding: both
    # envelope figures previously filtered to orientation=="c_plane" only).
    plots["envelope_pulse.png"] = plot_envelope(out, "envelope_pulse.png", [r for r in core if r["regime"] == "rectangular" and r["radius_nm"] == 5. and r["orientation"] in ("c_plane", "a_plane")], "pulse regime, r=5nm")
    plots["envelope_set.png"] = plot_envelope(out, "envelope_set.png", [r for r in core if r["regime"] == "deterministic_pair" and r["radius_nm"] == 5. and r["orientation"] in ("c_plane", "a_plane")], "SET regime, r=5nm")
    if sens:
        plots["set_feasibility.png"] = plot_set_island(out, sens)
    if bias:
        plots["stark_energy_bias.png"] = plot_stark(out, "stark_energy_bias.png", bias, "V_j", "E_X_eV", "E_X(V_j)", zhang_anchor=True)
        plots["stark_tau_bias.png"] = plot_stark(out, "stark_tau_bias.png", bias, "V_j", "tau_rad_bare_ns", "bare tau_rad(V_j)")
        plots["stark_tau_cavity_bias.png"] = plot_stark(out, "stark_tau_cavity_bias.png", bias, "V_j", "tau_rad_cavity_ns", "cavity tau_rad(V_j)")
        plots["stark_overlap_bias.png"] = plot_stark(out, "stark_overlap_bias.png", bias, "V_j", "overlap_sq", "overlap(V_j)")
    if current:
        plots["stark_energy_current.png"] = plot_stark(out, "stark_energy_current.png", current, "current_uA", "E_X_eV", "E_X(I), self-consistent T_j (heating)")
        plots["stark_tau_current.png"] = plot_stark(out, "stark_tau_current.png", current, "current_uA", "tau_rad_bare_ns", "tau_rad(I), self-consistent T_j (heating)", ylog=False)

    complete = (not a.quick) and counter["evaluate_calls"] <= a.max_evaluations and (time.time() - t0) < 1800
    runtime_s = time.time() - t0
    (out / "results.md").write_text(
        _results_md(core, geo, bias, current, comp, refinement_checks, complete, runtime_s, counter["evaluate_calls"]),
        encoding="utf-8")

    files = [p.name for p in out.iterdir() if p.is_file() and p.name != "manifest.json"]
    hashes = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}
    # Actual call counts by kind and invalid counts/reasons (Opus fix-round
    # medium finding: the manifest previously recorded no invalid/skipped
    # counts or reasons and no actual evaluate() call counts by kind, so the
    # 18 valid=False stark_bias rows in the fix-round-1 review appeared
    # nowhere). refinement_probe rows are included even though they are not
    # written to their own CSV (they are diagnostic-only, already referenced
    # by row_id inside convergence_checks).
    all_kind_rows = core + geo + bias + current + sens + probe_rows
    evaluate_calls_by_kind = {}; invalid_by_kind = {}; invalid_reasons_summary = {}
    for r in all_kind_rows:
        k = r.get("row_kind", "unknown")
        if not r.get("cache_hit"):
            evaluate_calls_by_kind[k] = evaluate_calls_by_kind.get(k, 0) + 1
        if not r.get("valid"):
            invalid_by_kind[k] = invalid_by_kind.get(k, 0) + 1
            try:
                row_reasons = json.loads(r.get("invalid_reasons") or "[]")
            except (TypeError, ValueError):
                row_reasons = []
            for reason in row_reasons:
                invalid_reasons_summary[reason] = invalid_reasons_summary.get(reason, 0) + 1
    manifest = {
        "quick": a.quick, "stark": a.stark, "complete": complete, "runtime_s": runtime_s,
        "evaluate_calls": counter["evaluate_calls"], "max_evaluations": a.max_evaluations,
        "actual_evaluate_calls_by_kind": evaluate_calls_by_kind,
        "invalid_counts_by_kind": invalid_by_kind,
        "invalid_reasons_summary": invalid_reasons_summary,
        "invalid_total": sum(invalid_by_kind.values()),
        "requested_rows": {"core": len(core), "shape": len([r for r in geo if r["row_kind"] == "shape"]),
                            "qw": len([r for r in geo if r["row_kind"] == "qw"]),
                            "literature_aux": len([r for r in geo if r["row_kind"] == "literature_aux"]),
                            "stark_bias": len(bias), "stark_current": len(current),
                            "sensitivities": len(sens), "compatibility": len(comp),
                            "refinement_probes": len(probe_rows)},
        "full_contract_rows": {"core": 3360, "shape": 864, "qw": 288, "stark": 3120, "total_before_caching": 7632},
        "axes": {"height_nm": list(CORE_H), "radius_nm": list(CORE_R), "orientation": list(ORI),
                 "T_hs": list(TS), "screening_fraction": list(SCR), "regime": list(REG)},
        "cache_hits": sum(1 for r in all_kind_rows if r.get("cache_hit")),
        "numeric_convergence_policy": "derivative refinement: halved voltage step must agree to <=0.5 meV/V absolute or <=5% relative [A convergence criteria]",
        "convergence_checks": refinement_checks,
        "slope_interval": {"range_meV_per_V": list(slope_range), "kind": "user_supplied" if a.slope_range else "illustrative_not_measured", "bias_window_V": list(a.bias_window)},
        "plot_row_mapping": plots,
        "output_hashes": hashes,
        "card_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT / "cards").glob("nitride-*.yaml")},
        "source_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT / "fsim_core").glob("*.py")},
        "literature_annotations": ["Wang2017SciRep", "Wang2017APL_FSS", "Zhang2016", "Deshpande2014"],
        "versions": {"python": sys.version.split()[0]},
        "resolved_out_dir": str(out),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print("evaluate_calls=%d runtime_s=%.2f complete=%s" % (counter["evaluate_calls"], runtime_s, complete))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
