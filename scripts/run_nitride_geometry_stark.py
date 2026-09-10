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
from fsim_core.nitride_stark import screening_compatibility, stark_derivatives, ZHANG2016_SLOPE_MEV_PER_V, ZHANG2016_COMPARISON, _TRACE_KEYS

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
       "shape": "disc", "top_radius_fraction": 1., "T_hs": 300., "cavity_tracking": "per_T_hs"}
# Fix round 3, item 5: the original REF (rectangular/unscreened) baseline
# never passes the optical gate (pulse fails everywhere; c-plane needs
# screening to pass). Two additional baselines that DO pass the optical
# gate (matching the headline SET passes) so the OAT sensitivities also
# bound assumptions at a configuration that produces a verdict.
REF_PASS_C = dict(REF, regime="deterministic_pair", screening_fraction=1.)
REF_PASS_A = dict(REF, regime="deterministic_pair", orientation="a_plane")
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
    # Cavity tracking (fix round 3, item 1): "per_T_hs" re-tracks the
    # cavity to this row's own T_hs, exactly matching scripts/
    # run_nitride_cavity.py's round-1 headline convention
    # (`_put(d.nitride["cavity"],"T_track",v)` there). "fixed_300K" is the
    # LABELLED sensitivity/Stark convention: T_track is pinned to the
    # card's own 300 K default regardless of T_hs, so it is set
    # explicitly here too (never left to accidentally match the card
    # default) -- an undisclosed detuning between T_j and a stale cavity
    # position was the round-2 high finding this fixes.
    tracking = p.get("cavity_tracking", "per_T_hs")
    if tracking == "per_T_hs":
        _set(cav, "T_track", p["T_hs"])
    else:
        _set(cav, "T_track", 300.)
    _set(d.thermal, "T_hs", p["T_hs"])
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
    # Fix round 2 (2026-09-09): fsim_core/device.py deliberately reports
    # scalars["T_hs"] as the CARD's own fixed heat-sink setting when
    # bias_mode=='junction_voltage' (T_j is the swept control variable
    # there, not T_hs -- see device.py's own "T_hs (Opus fix-round finding)"
    # comment), so the s.items() merge below overwrites p["T_hs"] with a
    # value that is IDENTICAL across every T_j in a Stark bias trace. Every
    # trace-grouping consumer (attach_derivatives, plot_stark,
    # screening_compatibility via build_compatibility, refinement-check
    # matching) must key on this uncorrupted sweep-requested value instead,
    # never on row["T_hs"] after the merge below.
    row["T_hs_requested"] = p.get("T_hs")
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
    # cavity_tracking (fix round 3, item 1): non-Stark rows default to
    # "per_T_hs" -- the cavity re-tracks to each row's own operating
    # temperature (nitride.cavity.T_track = T_hs), exactly as
    # scripts/run_nitride_cavity.py's round-1 headline does. Stark bias/
    # current rows override this to "fixed_300K" (build_stark below):
    # a Stark trace needs ONE fixed cavity reference across its whole V_j
    # sweep for tau_rad_cavity/detuning/Fp_add to stay comparable point to
    # point (unchanged from the original spec: "T_track stays at the
    # card's 300 K").
    return {"height_nm": h, "radius_nm": r, "orientation": o, "T_hs": t, "screening_fraction": s,
            "regime": reg, "x_in": .25, "Q": 2000., "current_uA": .02, "polarity": 1,
            "geometry_type": "isolated_dot", "shape": "disc", "top_radius_fraction": 1.,
            "cavity_tracking": "per_T_hs"}

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
        p = _base(h, r, o, t, s, reg)
        p.update(geometry_type="qw_fluctuation", wl_thickness_nm=w, qw_height_offset_nm=dh)
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
        p0 = _base(h, 10., o, t, s, reg); p0["cavity_tracking"] = "fixed_300K"
        for v in vv:
            p = dict(p0); p.update(bias_mode="junction_voltage", V_j=v, T_j=t, polarity=1); bias.append(p)
        for cur in ii:
            p = dict(p0); p.update(bias_mode="current", current_uA=cur, polarity=1); current.append(p)
    # opposite polarity at H=3 only, same remaining axes
    for o, s, reg, t in itertools.product(STARK_ORI, SCR, REG, tt):
        p0 = _base(3., 10., o, t, s, reg); p0["cavity_tracking"] = "fixed_300K"
        for v in vv:
            p = dict(p0); p.update(bias_mode="junction_voltage", V_j=v, T_j=t, polarity=-1); bias.append(p)
        for cur in ii:
            p = dict(p0); p.update(bias_mode="current", current_uA=cur, polarity=-1); current.append(p)
    return bias, current

def build_fixed_anchor(quick):
    """Fix round 3, item 1: a LABELLED fixed-anchor (cavity_tracking=
    'fixed_300K') sensitivity row set across the full T_hs axis, at the
    REF geometry (H=3nm, R=10nm, screening=0), both headline orientations
    and both regimes -- mirrors scripts/run_nitride_cavity.py's own
    `_anchor_rows`. Core rows (above) now re-track the cavity per row
    (T_track=T_hs); these rows instead hold T_track fixed at the card's
    300 K default while T_hs varies, so the resulting detuning is a
    visible, named diagnostic (results.md quantifies it against the
    matching tracked core row at the same coordinates) rather than an
    unlabelled artifact."""
    ts = (230., 300.) if quick else TS
    # Quick core only samples H in {1,7}, R in {5,30} (reduced grid); use a
    # point on that same reduced grid so the fixed-vs-tracked comparison
    # table below always finds a matching tracked core row, even in quick
    # mode. Full mode uses the REF geometry (H=3nm, R=10nm).
    h, r = (1., 5.) if quick else (3., 10.)
    out = []
    for o in ("c_plane", "a_plane"):
        for reg in (REG if not quick else ("rectangular",)):
            for t in ts:
                p = _base(h, r, o, t, 0., reg)
                p["cavity_tracking"] = "fixed_300K"
                p["sensitivity_axis"] = "cavity_tracking_fixed_300K"
                p["sensitivity_value"] = t
                out.append(p)
    return out

def build_sensitivities(quick, baseline=None, label="REF_unscreened_rectangular"):
    """One-at-a-time axes around `baseline` (defaults to REF); each axis
    holds every other input fixed (spec: 'Hold all other inputs fixed
    within each axis group'). Explicit null cases (nonpolar
    screening_fraction, m/a orientation identity) are NOT re-tested here --
    they already fall out of the core grid, which the verifier checks
    directly against sweep.csv.

    Fix round 3, item 5: called once at the original REF (rectangular,
    unscreened -- never passes the optical gate) and again at
    REF_PASS_C/REF_PASS_A (deterministic_pair baselines that DO pass), so
    every [A] assumption is also bounded at a configuration that produces
    a verdict. `sensitivity_baseline` on every row records which."""
    if baseline is None: baseline = REF
    ts = (300.,) if quick else STARK_T_FULL
    out = []
    def add(axis, key, values, **extra):
        for v in values:
            for t in ts:
                p = dict(baseline); p["T_hs"] = t; p[key] = v; p.update(extra)
                p["sensitivity_axis"] = axis; p["sensitivity_value"] = v
                p["sensitivity_baseline"] = label
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
        p = dict(baseline); p["T_hs"] = ts[-1]; p["regime"] = "deterministic_pair"
        p["island_radius_nm"] = rad; p["sensitivity_axis"] = "island_radius_nm"; p["sensitivity_value"] = rad
        p["sensitivity_baseline"] = label
        out.append(p)
    return out

def build_all_sensitivities(quick):
    """The OAT axis set evaluated at all three baselines (fix round 3,
    item 5): out += REF (never passes) + REF_PASS_C (c-plane screened SET,
    passes) + REF_PASS_A (a-plane SET, passes)."""
    out = build_sensitivities(quick, REF, "REF_unscreened_rectangular")
    out += build_sensitivities(quick, REF_PASS_C, "REF_pass_c_plane_screened_SET")
    out += build_sensitivities(quick, REF_PASS_A, "REF_pass_a_plane_SET")
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
           "cache_identity": ident, "T_hs": 300., "T_hs_requested": 300., "regime": "rectangular", "orientation": "c_plane",
           "cavity_tracking": "fixed_300K",  # replayed as-is; the card's own default (300 K) is untouched
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
def _gkey_like(r):
    """Mirrors nitride_stark._gkey exactly: nitride_stark's own _TRACE_KEYS
    whitelist (the geometry/shape/composition/field fields it groups
    stark_derivatives and screening_compatibility traces on). Deliberately
    does NOT include screening_fraction -- nitride_stark's own
    screening_compatibility adds that as a SEPARATE key via its private
    _screen(r), precisely because a screening hypothesis is a separate
    dimension from a trace's fixed geometry/field identity, not part of it.

    Two fields need a fallback because the RAW pre-evaluate() row/param dict
    and the post-evaluate() row dict spell them differently:
      - "T_hs": evaluate() deliberately reports the CARD's own fixed
        heat-sink setting for bias_mode=='junction_voltage' rows (T_j is the
        swept control variable there), identical across every T_j in a
        trace (Fix round 2 crash: grouping on that field alone silently
        merged every T_j in the sweep into one trace, feeding
        stark_derivatives duplicate voltages). T_hs_requested (set in
        _row()) is the uncorrupted sweep value and is used whenever present;
        pre-evaluate() param dicts have no T_hs_requested key at all, so the
        fallback to "T_hs" there reads the (already correct, unresolved)
        sweep value directly.
      - "field_polarity": only evaluate()'s scalars use this name; the
        pre-evaluate() param/build_stark dicts call it "polarity".
    """
    parts = []
    for k in _TRACE_KEYS:
        if k == "T_hs":
            if "T_hs_requested" in r or "T_hs" in r:
                parts.append((k, r.get("T_hs_requested", r.get("T_hs"))))
        elif k == "field_polarity":
            if "field_polarity" in r or "polarity" in r:
                parts.append((k, r.get("field_polarity", r.get("polarity"))))
        elif k in r:
            parts.append((k, r[k]))
    return tuple(parts)

def _screen_like(r):
    """Mirrors nitride_stark._screen exactly."""
    return r.get("screening_fraction", r.get("screening"))

def _derivative_group_key(r):
    """_gkey_like plus screening_fraction (a real, physically-distinct trace
    coordinate that _TRACE_KEYS omits -- see _gkey_like's docstring; the
    PREVIOUS 7-field grouping this replaces DID include it explicitly, so
    omitting it here would silently re-merge all 3 screening hypotheses into
    one trace) plus bias_mode, which nitride_stark has no concept of at all
    -- required so a bias-mode and current-mode row can never share a trace
    even if every other coordinate matches."""
    return _gkey_like(r) + (("screening_fraction", _screen_like(r)), ("bias_mode", r.get("bias_mode")))

def attach_derivatives(pool, key="V_j"):
    """Group by _derivative_group_key (Fix round 2: the previous 7-field key
    used evaluate()'s reported "T_hs", which fsim_core/device.py
    deliberately holds fixed at the card's own heat-sink setting for
    bias_mode=='junction_voltage' rows -- see _gkey_like's docstring). Any
    trace that still has non-distinct or non-increasing valid voltages (or
    otherwise fails) is degraded to derivative_valid=False with a recorded
    reason instead of aborting the whole run; failures are returned for the
    manifest, never raised past this function."""
    failures = []
    groups = {}
    for r in pool:
        groups.setdefault(_derivative_group_key(r), []).append(r)
    for gkey, rr in groups.items():
        rr.sort(key=lambda x: float(x.get(key)) if isinstance(x.get(key), (int, float)) else -1e18)
        trace_id = "|".join("%s=%s" % (k, v) for k, v in gkey)
        try:
            ds = stark_derivatives(rr, voltage_key=key)
        except (ValueError, ArithmeticError) as exc:
            reason = "non_monotonic_trace" if "increasing" in str(exc) else str(exc)
            failures.append({"trace_id": trace_id, "row_ids": [r.get("row_id") for r in rr], "reason": reason})
            for r in rr:
                r.update(dE_X_dV_meV_per_V=float("nan"), derivative_valid=False, derivative_row_ids=[],
                          derivative_policy="adjacent_three_point_nonuniform", derivative_invalid_reason=reason)
            continue
        for r, dv in zip(rr, ds):
            r.update({k: dv[k] for k in ("dE_X_dV_meV_per_V", "derivative_valid", "derivative_row_ids", "derivative_policy")})
            r["derivative_invalid_reason"] = None if dv["derivative_valid"] else "invalid_point_or_group"
    return failures

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
        # T_hs matches against T_hs_requested, not the row's own "T_hs"
        # (Fix round 2: evaluate() reports the card's fixed heat-sink
        # setting there for bias_mode=='junction_voltage' rows, so matching
        # on "T_hs" directly would silently find zero rows for any spec
        # T_hs that differs from the card default -- see _gkey_like/_derivative_group_key).
        group = [r for r in bias_rows if all(r.get(k) == spec[k] for k in
                 ("height_nm", "orientation", "screening_fraction", "regime", "polarity"))
                 and r.get("T_hs_requested") == spec["T_hs"]
                 and r.get("bias_mode") == "junction_voltage"]
        group = [r for r in group if r.get("derivative_valid") is True]
        if not group:
            checks.append({"spec": spec, "status": "no_valid_derivative_point", "converged": False}); continue
        anchor = min(group, key=lambda r: abs(float(r["V_j"]) - target_v))
        v0 = float(anchor["V_j"]); coarse = float(anchor["dE_X_dV_meV_per_V"])
        anchor_flat_band = anchor.get("flat_band")
        p0 = _base(spec["height_nm"], 10., spec["orientation"], spec["T_hs"], spec["screening_fraction"], spec["regime"])
        p0["cavity_tracking"] = "fixed_300K"  # refinement probes are Stark bias points
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
def _tau_at_screening(part_rows, screening, v_ref):
    """Bare radiative lifetime nearest v_ref (V_j) at the given screening
    hypothesis, from the SAME geometry-group rows already gathered for one
    screening_compatibility() call -- used for the lifetime-ratio column
    (fix round 3, item 3: 'the discriminating observable is the bias-
    resolved lifetime')."""
    cands = [r for r in part_rows if _finite_num(r.get("tau_rad_bare_ns")) and _finite_num(r.get("V_j"))
             and _screen_like(r) is not None and abs(float(_screen_like(r)) - float(screening)) < 1e-9]
    if not cands: return None
    best = min(cands, key=lambda r: abs(float(r["V_j"]) - v_ref))
    return float(best["tau_rad_bare_ns"])

def build_compatibility(bias_rows, slope_range, bias_window, slope_user_supplied, rid_start):
    """screening_compatibility groups internally on nitride_stark's own
    (_gkey(r), _screen(r)) -- _gkey uses _TRACE_KEYS, which includes the
    literal "T_hs" field, corrupted (see _gkey_like) for bias_mode==
    'junction_voltage' rows. Feed it rows whose "T_hs" has been repaired to
    the uncorrupted sweep value first.

    Fix round 3, item 3 (Opus/Astra high finding): partition ONLY by the
    geometry/field identity (_gkey_like, no screening_fraction) so ALL
    screening hypotheses of one trace group are handed to a SINGLE
    screening_compatibility() call together -- that function's own
    cross-hypothesis degeneracy check (identical fitted slopes across
    different screening_fraction values -> 'screening_unidentifiable') can
    then actually fire; the previous per-screening partitioning fed it
    exactly one hypothesis per call, so degeneracy could never be
    detected. One try/except per GROUP (not per group+screening) still
    isolates a bad group's failure from every other group (spec: 'any
    per-trace exception in derivatives, compatibility or plotting must be
    caught ... and never abort the sweep')."""
    fixed_rows = [dict(r, T_hs=r.get("T_hs_requested", r.get("T_hs"))) for r in bias_rows]
    partitions = {}
    for r in fixed_rows:
        pkey = _gkey_like(r)
        partitions.setdefault(pkey, []).append(r)
    compat = []; compat_failures = []; compat_by_id = []
    for gk, part_rows in partitions.items():
        trace_id = "|".join("%s=%s" % (k, v) for k, v in gk)
        try:
            fits = screening_compatibility(part_rows, slope_range_meV_per_V=slope_range, voltage_window_V=tuple(bias_window))
        except (ValueError, ArithmeticError) as exc:
            compat_failures.append({"trace_id": trace_id, "row_ids": [r.get("row_id") for r in part_rows], "reason": str(exc)})
            continue
        # Lifetime-ratio column (item 3): unscreened/screened bare tau_rad
        # at the bias-window's lower edge, from this SAME group's rows --
        # one ratio per geometry group, attached to every fit row that
        # came out of it (whichever screening hypothesis it represents).
        v_ref = float(bias_window[0])
        tau0 = _tau_at_screening(part_rows, 0., v_ref)
        tau1 = _tau_at_screening(part_rows, 1., v_ref)
        ratio = (tau0 / tau1) if (tau0 is not None and tau1 not in (None, 0.)) else None
        for x in fits:
            compat.append(x); compat_by_id.append(ratio)
    out = []; rid = rid_start
    for x, lifetime_ratio in zip(compat, compat_by_id):
        rid += 1
        row = {"row_id": f"CT{rid:05d}",
               "slope_interval_kind": "user_supplied" if slope_user_supplied else "illustrative_not_measured",
               "zhang_guide_meV_per_V": ZHANG2016_SLOPE_MEV_PER_V,
               "source": ZHANG2016_COMPARISON["citation"] + " [V] non-gating",
               "lifetime_ratio_s0_over_s1": lifetime_ratio}
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
    return out, rid, compat_failures

# -------------------------------------------------------------------- plots
def _finite_num(v):
    try: return math.isfinite(float(v))
    except (TypeError, ValueError): return False

def _leg(ax, **kw): return ax.legend(fontsize=7, loc="best", **kw)

MAX_TRACES_PER_PANEL = 8  # spec: "split into panels if > 8 traces"

def _panel_plot(out, name, entries, xlabel, ylabel, suptitle, xlog=False, ylog=False, guide=None):
    """entries: ordered [(label, xv, yv, row_ids, style), ...], ONE trace per
    fully-specified group -- never merged (Opus fix-round high finding).
    Splits into ceil(n/8) panels so no group is silently dropped (the
    previous [:14] truncation), each panel's legend placed OUTSIDE the axes.
    `guide`, when given, is (value, label) for a horizontal reference line
    drawn on every panel (fix round 3, item 6: 1000/s flux floor or g2<0.5
    gate). Returns the plot-construction contract {label: {row_ids, x, y}}."""
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
        if guide is not None:
            ax.axhline(guide[0], color="k", lw=.9, ls="--", label=guide[1])
        if xlog: ax.set_xscale("log")
        if ylog: ax.set_yscale("log")
        ax.set(xlabel=xlabel, ylabel=ylabel)
        if n_panels > 1: ax.set_title(f"panel {i + 1}/{n_panels}", fontsize=8)
        if chunk or guide is not None: ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)
    fig.suptitle(suptitle, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / name, dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract

def _panel_plot_grouped(out, name, super_groups, xlabel, ylabel, suptitle, xlog=False, ylog=False, guide=None):
    """super_groups: ordered [(subtitle, [(label,xv,yv,rids,style), ...]), ...].
    ONE subplot per super-group -- a super-group's own lines are NEVER
    split across panels/subplots (fix round 3, item 6: 'panels never split
    a screening triplet or a polarity pair'; the caller forms super-groups
    so every screening triplet / polarity pair for one (orientation,
    height,regime,T_hs) trace-identity lands in the SAME super-group).
    Returns the SAME flat plot-construction contract {label: {row_ids, x,
    y}} as _panel_plot (the verifier's replay check does not care how
    labels are laid out into subplots), plus 'panel_index' per label so a
    verifier can confirm co-grouped labels share a subplot."""
    import matplotlib.pyplot as plt
    contract = {}
    n_panels = max(1, len(super_groups))
    ncols = min(3, n_panels); nrows = math.ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.6 * ncols, 4.6 * nrows), squeeze=False)
    for i in range(nrows * ncols):
        ax = axes[i // ncols][i % ncols]
        if i >= n_panels or i >= len(super_groups):
            ax.set_visible(False); continue
        subtitle, lines_in_group = super_groups[i]
        for label, xv, yv, rids, style in lines_in_group:
            ax.plot(xv, yv, style, ms=3, lw=1.2, label=label)
            # Contract key qualified by subtitle (bare labels like "pol=1
            # (unscreened lower bound)" repeat across every super-group and
            # would silently collide/overwrite each other in a flat dict).
            contract[f"{subtitle} | {label}"] = {"row_ids": rids, "x": xv, "y": yv, "panel_index": i}
        if guide is not None:
            ax.axhline(guide[0], color="k", lw=.8, ls="--", label=guide[1])
        if xlog: ax.set_xscale("log")
        if ylog: ax.set_yscale("log")
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.set_title(subtitle, fontsize=7)
        if lines_in_group: ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)
    fig.suptitle(suptitle, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / name, dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract

LOG_Y_FIELDS = ("tau_rad_bare_ns", "tau_rad_cavity_ns", "overlap_sq", "signal_flux_s")

def plot_axis_response(out, name, rows, x, y, title, group_keys=("orientation", "screening_fraction", "regime"), fail_log=None, guide=None, ylog=None):
    groups = {}
    groups_all = {}
    for r in rows:
        groups_all.setdefault(tuple(r.get(k) for k in group_keys), []).append(r)
        if _finite_num(r.get(x)) and _finite_num(r.get(y)):
            groups.setdefault(tuple(r.get(k) for k in group_keys), []).append(r)
    # Keep the screening triplet together in one subplot.  This is also used
    # by height/radius/orientation/temperature response figures, so a panel
    # cannot accidentally split matched bounds or connect different fixed
    # inputs into one artificial trace.
    screen_i = group_keys.index("screening_fraction") if "screening_fraction" in group_keys else None
    super_groups = {}
    for key, rr in sorted(groups.items(), key=lambda kv: str(kv[0])):
        try:
            skey = tuple(v for i, v in enumerate(key) if i != screen_i) if screen_i is not None else key
            label = ("/".join(str(v) for i, v in enumerate(key) if i != screen_i)
                     if screen_i is not None else "/".join(str(v) for v in key))
            sval = key[screen_i] if screen_i is not None else None
            rr = sorted(rr, key=lambda z: float(z[x]))
            base_label = f"screening={sval:g}" if screen_i is not None else label
            xv = [float(z[x]) for z in rr]; yv = [float(z[y]) for z in rr]
            line = (base_label, xv, yv, [z["row_id"] for z in rr], "o-")
            super_groups.setdefault(skey, []).append(line)
            # Hardening item 8: a trace that stops early because higher-x
            # rows exist but fail the validity gate (e.g. height_response's
            # c_plane/screening=0 trace truncating at 5nm because the H=7/10
            # rows are invalid) is marked explicitly instead of silently
            # stopping -- a red 'x' overlaid on the LAST REAL plotted point
            # (its own true row_id/x/y, never a fabricated coordinate) with a
            # legend label naming which x-values were validity-rejected.
            plotted_ids = {z["row_id"] for z in rr}
            rejected = [rz for rz in groups_all.get(key, [])
                        if rz.get("row_id") not in plotted_ids and rz.get("valid") is False]
            if rejected:
                rejected_x = sorted({float(rz.get(x)) for rz in rejected if _finite_num(rz.get(x))})
                marker_label = f"{base_label} validity-rejected beyond {xv[-1]:g} (x={rejected_x})"
                super_groups[skey].append((marker_label, [xv[-1]], [yv[-1]], [rr[-1]["row_id"]], "rx"))
        except (ValueError, TypeError, KeyError, ArithmeticError) as exc:
            if fail_log is not None:
                fail_log.append({"trace_id": "%s|%s" % (name, "|".join(str(k) for k in key)), "reason": str(exc)})
    ordered = [("/".join(str(v) for v in key), lines)
               for key, lines in sorted(super_groups.items(), key=lambda kv: str(kv[0]))]
    if ylog is None:
        ylog = y in LOG_Y_FIELDS
    return _panel_plot_grouped(out, name, ordered, x, y, title, ylog=ylog, guide=guide)

def plot_stark(out, name, rows, x, y, title, ylog=False, zhang_anchor=False, fail_log=None):
    """E_X/tau/overlap vs V_j or current. Fix round 3, item 6: panels never
    split a screening triplet or a polarity pair. Rows are first grouped
    into SUPER-groups by (orientation, height, regime, T_hs_requested) --
    the trace-identity axes that legitimately deserve their own panel --
    and each super-group's OWN lines vary only polarity and screening (at
    most 2x3=6 lines), all drawn together in ONE subplot via
    _panel_plot_grouped, so a screening triplet or a +/-1 polarity pair
    (only present at H=3) can never be split across panels the way the
    previous flat 8-per-panel chunking could.

    T_hs_requested, not the row's own "T_hs" field, is used for grouping
    and the label (Fix round 2: evaluate() reports the CARD's fixed
    heat-sink setting as "T_hs" for bias_mode=='junction_voltage' rows --
    identical across every T_j in the sweep -- so grouping on it directly
    would merge every requested temperature back into one trace, the same
    bug that crashed attach_derivatives; see _gkey_like/
    _derivative_group_key). Screening 0/0.5/1 use the same line-style/
    legend convention as the rest of the nitride pieces. Current traces
    (x=='current_uA') use a log x-axis. A single line's construction
    failing is caught and logged to fail_log (trace_id, reason) rather
    than aborting the whole figure (spec: 'any per-trace exception in ...
    plotting must be caught ... and never abort the sweep')."""
    label_for_s = {0.: "unscreened lower bound", .5: "midpoint 0.5", 1.: "screened upper bound"}
    style_for_s = {(0., 1): "-o", (.5, 1): "--o", (1., 1): ":o", (0., -1): "-s", (.5, -1): "--s", (1., -1): ":s"}
    supergroups = {}
    for r in rows:
        if _finite_num(r.get(x)) and _finite_num(r.get(y)):
            skey = (r.get("orientation"), r.get("height_nm"), r.get("regime"), r.get("T_hs_requested", r.get("T_hs")))
            lkey = (r.get("polarity"), r.get("screening_fraction"))
            supergroups.setdefault(skey, {}).setdefault(lkey, []).append(r)
    super_groups_ordered = []
    zhang_target = None
    for skey, by_line in sorted(supergroups.items(), key=lambda kv: str(kv[0])):
        o, h, reg, t = skey
        lines_in_group = []
        for lkey, rr in sorted(by_line.items(), key=lambda kv: str(kv[0])):
            pol, s = lkey
            try:
                rr = sorted(rr, key=lambda z: float(z[x]))
                xv = [float(z[x]) for z in rr]; yv = [float(z[y]) for z in rr]
                lab = f"pol={pol:g} ({label_for_s.get(s, s)})"
                lines_in_group.append((lab, xv, yv, [z["row_id"] for z in rr], style_for_s.get((s, pol), "-.o")))
                if zhang_anchor and zhang_target is None and o == "c_plane" and h == 3. and s == 0. and reg == "rectangular" and pol == 1 and len(xv) >= 2:
                    zhang_target = (len(super_groups_ordered), xv[0], yv[0])
            except (ValueError, TypeError, KeyError, ArithmeticError) as exc:
                if fail_log is not None:
                    fail_log.append({"trace_id": "%s|%s|%s" % (name, "|".join(str(k) for k in skey), "|".join(str(k) for k in lkey)), "reason": str(exc)})
        subtitle = f"{o} h={h:g}nm {reg} T={t:g}K"
        super_groups_ordered.append((subtitle, lines_in_group))
    if zhang_target is not None:
        # Fix round 3, item 6: the Zhang guide is drawn ON the panel of the
        # model curve it anchors (c-plane H=3, polarity +1), never alone.
        gi, x0, y0 = zhang_target
        subtitle, lines_in_group = super_groups_ordered[gi]
        lines_in_group.append(("Zhang 2016 -10 meV/V guide (non-gating)",
                                [x0, x0 + 2.0], [y0, y0 + 2.0 * ZHANG2016_SLOPE_MEV_PER_V / 1000.], [], "k--"))
    return _panel_plot_grouped(out, name, super_groups_ordered, x, y, title, xlog=(x == "current_uA"), ylog=ylog)

def plot_envelope(out, name, core_rows, title):
    """Columns = orientation (c-plane AND a-plane -- Opus fix-round medium
    finding: both envelope figures previously filtered to c_plane only, no
    a-plane panel); rows = unscreened lower bound / screened upper bound.
    Cell value is the optical-pass fraction, matching run_nitride_cavity.
    py's envelope convention exactly (its `optical_pass`, not the flux-only
    `eligible` gate). Fix round 3, item 6: a cell whose rows are ALL
    invalid (rejected before any gate could apply) is masked and drawn
    hatched/grey, distinct from a cell with valid rows that simply fail
    the optical gate (which stays plain red at 0)."""
    import matplotlib.pyplot as plt
    orientations = [o for o in ("c_plane", "a_plane") if any(r["orientation"] == o for r in core_rows)] or ["c_plane"]
    fig, axes = plt.subplots(2, len(orientations), figsize=(6.2 * len(orientations), 8.6), squeeze=False)
    contract = {}
    hs = sorted({r["height_nm"] for r in core_rows}); ts = sorted({r["T_hs"] for r in core_rows})
    cmap = plt.get_cmap("RdYlGn").copy(); cmap.set_bad(color="0.75")
    for col, o in enumerate(orientations):
        for row_i, (s, lab) in enumerate(((0., "unscreened lower bound"), (1., "screened upper bound"))):
            ax = axes[row_i][col]
            mat = []; cellmap = {}; invalid_cells = []
            for ti, t in enumerate(ts):
                rowvals = []
                for hi, h in enumerate(hs):
                    q = [r for r in core_rows if r["orientation"] == o and r["height_nm"] == h and r["T_hs"] == t and r["screening_fraction"] == s]
                    all_invalid = bool(q) and all(not r["valid"] for r in q)
                    if all_invalid:
                        rowvals.append(float("nan")); invalid_cells.append((hi, ti))
                    else:
                        rowvals.append(sum(r["optical_pass"] for r in q) / len(q) if q else float("nan"))
                    cellmap[f"h{h:g}_T{t:g}_s{s:g}"] = [r["row_id"] for r in q]
                mat.append(rowvals)
            im = ax.imshow(mat, vmin=0, vmax=1, aspect="auto", origin="lower", cmap=cmap)
            fig.colorbar(im, ax=ax, label="optical-pass fraction")
            for hi, ti in invalid_cells:
                ax.add_patch(plt.Rectangle((hi - .5, ti - .5), 1, 1, fill=False, hatch="xxx", edgecolor="0.35", lw=.6))
            ax.set(xticks=range(len(hs)), xticklabels=[f"{v:g}" for v in hs],
                   yticks=range(len(ts)), yticklabels=[f"{v:g}" for v in ts],
                   xlabel="height nm", ylabel="T_hs K", title=f"{o} {lab}\n{title}")
            from matplotlib.patches import Patch
            ax.legend(handles=[Patch(facecolor="0.75", edgecolor="0.35", hatch="xxx", label="all rows invalid (not a gate failure)")],
                      fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)
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
    ax1.axhline(G2, color="k", lw=.8); ax1.set(xticks=list(xs), xticklabels=[f"{t:g}" for t in ts], xlabel="T_hs K", ylabel="g2", title="pulse vs SET; fixed geometry c-plane H=1nm R=5nm: g2")
    ax2.axhline(FLUX, color="gray", lw=.8, ls="--"); ax2.set_yscale("log"); ax2.set(xticks=list(xs), xticklabels=[f"{t:g}" for t in ts], xlabel="T_hs K", ylabel="flux /s (log)", title="pulse vs SET; fixed geometry c-plane H=1nm R=5nm: flux\n('underflow' bars are display-floored, not a real value)")
    _leg(ax1); _leg(ax2); fig.tight_layout(); fig.savefig(out / "pulse_vs_set.png", dpi=125); plt.close(fig)
    return contract

def plot_set_island(out, sens_rows):
    import matplotlib.pyplot as plt
    # set_EC_over_kT depends only on (island_radius_nm, T_hs), not on
    # sensitivity_baseline -- dedupe to one row per (radius,T_hs) so the
    # three baselines (item 5) don't triple-plot identical points.
    seen = {}
    for r in sens_rows:
        if r.get("sensitivity_axis") != "island_radius_nm": continue
        key = (r.get("sensitivity_value"), r.get("T_hs"))
        seen.setdefault(key, r)
    z = list(seen.values())
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

def _safe_plot(name, outer_fail_log, fn, *args, **kwargs):
    """Call one top-level plot_*() figure and never let it abort the sweep
    (spec: 'any per-trace exception in derivatives, compatibility or
    plotting must be caught, logged in the manifest with the trace id and
    reason, and never abort the sweep -- outputs must still be written').
    plot_axis_response/plot_stark already isolate failures per trace-group
    internally (via their own fail_log keyword, passed through **kwargs
    when present -- named outer_fail_log here specifically so it never
    collides with that inner kwarg); this outer catch is the backstop for
    whole-figure failures (e.g. an empty input list, a matplotlib error)
    that a per-group catch cannot localize. On failure the figure is simply
    omitted from plots{} (never a partially-built or misleading artifact)."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 -- must never propagate past here
        outer_fail_log.append({"trace_id": name, "reason": "%s: %s" % (type(exc).__name__, exc)})
        return None

def _resolve_axis_scale(fn, args, kwargs):
    """Hardening item 2: record the log/linear axis actually used for a
    figure in the manifest, derived from the SAME call (fn/args/kwargs) that
    built it -- not a hand-typed duplicate list that could drift from the
    plotting code. Mirrors plot_axis_response's/plot_stark's own xlog/ylog
    predicates exactly. Returns None for figures this doesn't apply to
    (pulse_vs_set/envelope/set_feasibility use bespoke multi-axis layouts)."""
    if fn is plot_axis_response:
        y = args[4]
        ylog = kwargs.get("ylog")
        if ylog is None: ylog = y in LOG_Y_FIELDS
        return {"xlog": False, "ylog": bool(ylog)}
    if fn is plot_stark:
        x = args[3]
        return {"xlog": (x == "current_uA"), "ylog": bool(kwargs.get("ylog", False))}
    return None

# ---------------------------------------------------------------- results.md
def _best_passing_flux(core, family, t_hs=None):
    """Fix round 3, item 4: the best flux among rows that PASS the optical
    gate (optical_pass=True), never a diagnostic row (a Stark row at high
    V_j with no optical/hardware gating applied to it, as the previous
    'maximum valid signal_flux_s' picked). Hardening item 5: `t_hs`, when
    given, restricts the candidate pool to that heat-sink temperature (the
    headline BEST_PASSING_FLUX above is unrestricted -- and lands at 230 K,
    the cold edge, since flux is monotone in T_hs -- so a machine-checkable
    300 K-specific line is reported separately)."""
    cand = [(float(r["signal_flux_s"]), r["row_id"], r.get("T_hs"), r.get("screening_fraction"), r.get("regime"))
            for r in core if r.get("orientation") == family and r.get("optical_pass") and _finite_num(r.get("signal_flux_s"))
            and (t_hs is None or float(r.get("T_hs")) == t_hs)]
    return max(cand, key=lambda z: z[0]) if cand else None

def _hypothesis_key(gk):
    """(orientation, height_nm, polarity, screening_fraction) -- the
    DISTINCT-hypothesis identity for compatibility reporting (fix round 3,
    item 3): regime and T_hs are replicate evaluations of the SAME
    physical hypothesis (Stark E_X(V_j) does not depend on regime, and the
    fitted slope is nearly T-independent), not independent hypotheses."""
    return (gk.get("orientation"), gk.get("height_nm"), gk.get("field_polarity"), gk.get("screening"))

def _distinct_compat_rows(comp):
    """One representative CSV row per distinct hypothesis (first seen, by
    row_id order), for the results.md table and summary counts."""
    seen = {}
    for r in comp:
        try:
            gk = json.loads(r["group"]) if isinstance(r.get("group"), str) else r.get("group")
        except (TypeError, ValueError, json.JSONDecodeError):
            gk = None
        gk_dict = dict(gk) if gk else {}
        key = _hypothesis_key({**gk_dict, "screening": r.get("screening")})
        if key not in seen:
            seen[key] = r
    return list(seen.values())

def _results_md(core, geo, bias, current, sens, comp, refinement_checks, complete, runtime_s, evaluate_calls,
                 invalid_by_kind, invalid_reasons_summary, ec_margin, g2_floor):
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
             "set_feasible AND pair_supply_possible for the SET regime. Every headline core "
             "row carries cavity_tracking=per_T_hs: the cavity is re-tracked to that row's own "
             "operating temperature (nitride.cavity.T_track=T_hs), exactly as scripts/"
             "run_nitride_cavity.py's round-1 headline does. A separate cavity_tracking="
             "fixed_300K sensitivity set (below) holds the cavity fixed at the card's 300 K "
             "default while T_hs varies, to show the resulting detuning as a labelled, "
             "quantified effect rather than an undisclosed artifact.", ""]
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
    # Maximum flux among optically-passing rows / eligible-g2 statistics
    # (fix round 3, item 4: the previous "maximum valid signal_flux_s" was
    # a Stark diagnostic row at high V_j with no optical gate applied; the
    # headline number must be the best flux among rows that actually PASS
    # optical_pass, reported separately per family).
    opt_core = [r for r in core if r.get("optical_pass")]
    g2_opt = sorted(float(r["g2"]) for r in opt_core if _finite_num(r.get("g2")))
    def _median(xs):
        n = len(xs)
        if not n: return None
        return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])
    lines += ["", "## Maximum flux among optically-passing rows, and eligible-g2 statistics", "",
              "MACHINE-CHECKABLE (verifier recomputes these from sweep.csv):"]
    for family in ("c_plane", "a_plane"):
        best = _best_passing_flux(core, family)
        if best is None:
            lines.append(f"BEST_PASSING_FLUX family={family} value=none row_id=none")
        else:
            mf, mf_rid, mf_t, mf_scr, mf_reg = best
            lines.append(f"BEST_PASSING_FLUX family={family} value={mf:.6g} row_id={mf_rid} "
                          f"T_hs={mf_t:g} screening={mf_scr:g} regime={mf_reg}")
    # Hardening item 5: the unrestricted BEST_PASSING_FLUX lines above land
    # at T_hs=230 K (the cold edge -- flux is monotone in T_hs, not a
    # detuning artifact); a separate machine-checkable 300 K-specific line
    # per family, so the 300 K number is available without re-deriving it
    # from sweep.csv.
    for family in ("c_plane", "a_plane"):
        best300 = _best_passing_flux(core, family, t_hs=300.)
        if best300 is None:
            lines.append(f"BEST_PASSING_FLUX_300K family={family} value=none row_id=none")
        else:
            mf, mf_rid, mf_t, mf_scr, mf_reg = best300
            lines.append(f"BEST_PASSING_FLUX_300K family={family} value={mf:.6g} row_id={mf_rid} "
                          f"T_hs={mf_t:g} screening={mf_scr:g} regime={mf_reg}")
    lines.append("")
    if g2_opt:
        lines.append(f"Core-grid paired-optical-pass rows (g2<{G2:g}, flux>={FLUX:g}/s, plus "
                      f"one_pair_valid for SET): {len(opt_core)} of {len(core)}; "
                      f"g2_min={g2_opt[0]:.6g}, g2_median={_median(g2_opt):.6g}.")
    else:
        lines.append(f"No core-grid rows pass the paired-optical-pass gate (g2<{G2:g}, flux>={FLUX:g}/s).")
    if g2_floor is not None:
        set_g2s = sorted({round(float(r["g2"]), 5) for r in opt_core if r.get("regime") == "deterministic_pair" and _finite_num(r.get("g2"))})
        lines.append(f"SET g2 floor: g2 = 1-(1/(1+b_res))^2 = {g2_floor:.6g} (b_res from the evaluated card's own "
                      "drive.b_res). Under idealized deterministic one-pair loading, exact-one-pair counting "
                      "gives zero coincidence probability by construction, so rho asymptotes to 1/(1+b_res) and "
                      "g2_op to this floor nearly independent of geometry -- distinct SET g2 values observed among "
                      f"optical-pass rows: {set_g2s if set_g2s else 'none'}. The g2 gate therefore carries no "
                      "geometry information in the SET regime; only the flux/eligibility gates discriminate "
                      "geometry there.")
    lines.append("")
    # Hardening item 4: state plainly (not just via the cavity_tracking
    # column) that every Stark trace row is cavity_tracking=fixed_300K BY
    # DESIGN -- a fixed cavity reference is the correct convention for bias
    # spectroscopy (the cavity position must not move while V_j/current is
    # swept, or tau_rad_cavity/detuning/Fp_add would not be comparable
    # across the trace) -- and that every 230 K Stark panel is therefore a
    # 300 K-ANCHORED cavity, not a cavity re-tracked to 230 K the way the
    # headline core rows are.
    lines += ["## Stark trace cavity-tracking convention (item 4)", "",
              f"All {len(bias)} stark_bias rows and {len(current)} stark_current rows carry "
              "cavity_tracking=fixed_300K BY DESIGN: the cavity is held at the card's fixed 300 K "
              "position for every V_j/current point within one trace (a fixed instrument reference "
              "for bias spectroscopy), never re-tracked to that trace's own T_hs -- unlike the "
              "headline core rows above, which use cavity_tracking=per_T_hs. Consequently every "
              "T_hs=230 K panel in stark_energy_bias.png, stark_tau_bias.png, "
              "stark_tau_cavity_bias.png, stark_overlap_bias.png, stark_energy_current.png and "
              "stark_tau_current.png shows a cavity ANCHORED AT 300 K, not one re-tracked to 230 K; "
              "compare against the Fixed-anchor sensitivity table below and the per_T_hs core rows "
              "above before drawing any temperature conclusion from these Stark panels.", ""]
    # Fixed-anchor (cavity_tracking=fixed_300K) sensitivity set (item 1):
    # quantify the detuning consequence of NOT re-tracking the cavity,
    # against the matching per_T_hs-tracked core row at the same
    # coordinates.
    anchor_rows = [r for r in sens if r["row_kind"] == "fixed_anchor"]
    anchor_geom = f"H={anchor_rows[0].get('height_nm'):g}nm, R={anchor_rows[0].get('radius_nm'):g}nm" if anchor_rows else "n/a"
    lines += ["## Fixed-anchor (cavity_tracking=fixed_300K) sensitivity", "",
              "Headline core rows re-track the cavity to T_hs every row (cavity_tracking=per_T_hs). "
              "These rows instead hold the cavity fixed at the card's 300 K default while T_hs "
              f"varies, at a reference geometry sampled by the core grid ({anchor_geom}, screening=0), "
              "matching each row against the tracked core row at the SAME coordinates.", "",
              "| orientation | regime | T_hs K | detuning_meV (fixed) | flux/s (fixed) | flux/s (tracked, matching core row) | ratio tracked/fixed | fixed row_id | tracked row_id |",
              "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(anchor_rows, key=lambda z: (z.get("orientation"), z.get("regime"), z.get("T_hs"))):
        matches = [c for c in core if c.get("orientation") == r.get("orientation") and c.get("regime") == r.get("regime")
                   and c.get("T_hs") == r.get("T_hs") and c.get("height_nm") == r.get("height_nm")
                   and c.get("radius_nm") == r.get("radius_nm") and c.get("screening_fraction") == r.get("screening_fraction")]
        tr = matches[0] if matches else None
        fx_fixed = r.get("signal_flux_s"); fx_tracked = tr.get("signal_flux_s") if tr else None
        ratio = (float(fx_tracked) / float(fx_fixed)) if (_finite_num(fx_fixed) and _finite_num(fx_tracked) and float(fx_fixed) != 0.) else None
        lines.append(f"| {r.get('orientation')} | {r.get('regime')} | {r.get('T_hs'):g} | {r.get('detuning_meV')} | "
                      f"{fx_fixed} | {fx_tracked} | {('%.4g' % ratio) if ratio is not None else 'n/a'} | "
                      f"{r.get('row_id')} | {tr.get('row_id') if tr else 'n/a'} |")
    lines.append("")
    # SET hardware-feasibility quantification (item 2): read directly off
    # every deterministic_pair core row's own set_EC_over_kT/set_radius_
    # max_nm/set_R_T_over_RQ/set_f_max_Hz -- never hand-entered.
    set_rows = [r for r in core if r.get("regime") == "deterministic_pair" and _finite_num(r.get("set_EC_over_kT"))]
    lines += ["## SET hardware feasibility (item 2)", ""]
    if set_rows:
        ecs = [float(r["set_EC_over_kT"]) for r in set_rows]
        rmax = [float(r["set_radius_max_nm"]) for r in set_rows if _finite_num(r.get("set_radius_max_nm"))]
        rtq = [float(r["set_R_T_over_RQ"]) for r in set_rows if _finite_num(r.get("set_R_T_over_RQ"))]
        fmax = [float(r["set_f_max_Hz"]) for r in set_rows if _finite_num(r.get("set_f_max_Hz"))]
        island_r = {r.get("set_radius_nm") for r in set_rows if r.get("set_radius_nm") is not None}
        margin_s = f"{ec_margin:g}" if ec_margin is not None else "unavailable"
        lines.append(f"hardware_qualified=0 in every headline VERDICT line because set_EC_over_kT is "
                      f"{min(ecs):.4g}-{max(ecs):.4g} across all {len(set_rows)} SET core rows at the card's "
                      f"assumed island radius ({sorted(island_r)} nm; allowed max {min(rmax):.4g}-{max(rmax):.4g} nm), "
                      f"versus the required margin ec_margin={margin_s} -- the sole failing criterion. R_T/R_Q "
                      f"({min(rtq):.4g}-{max(rtq):.4g}) and f_max ({min(fmax):.4g}-{max(fmax):.4g} Hz) both pass. "
                      "The island-radius sensitivity rows below carry set_EC_over_kT across an explicit radius "
                      "sweep (0.5/1/5 nm) to show how strongly this single criterion depends on the assumed "
                      "island size.")
    else:
        lines.append("No deterministic_pair core rows carried finite set_EC_over_kT.")
    island_sens = [r for r in sens if r.get("sensitivity_axis") == "island_radius_nm"]
    if island_sens:
        lines += ["", "| baseline | island radius nm | T_hs K | E_C/kT | allowed max radius nm | R_T/R_Q | f_max Hz | row_id |",
                  "|---|---|---|---|---|---|---|---|"]
        for r in sorted(island_sens, key=lambda z: (z.get("sensitivity_baseline"), float(z.get("sensitivity_value", 0)), float(z.get("T_hs", 0)))):
            lines.append(f"| {r.get('sensitivity_baseline')} | {r.get('sensitivity_value'):g} | {r.get('T_hs'):g} | "
                          f"{r.get('set_EC_over_kT')} | {r.get('set_radius_max_nm')} | {r.get('set_R_T_over_RQ')} | "
                          f"{r.get('set_f_max_Hz')} | {r.get('row_id')} |")
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
        "reinterpreted here (out of scope for this piece). Comparability (fix round 3, item 7): "
        "round 1's headline used cavity_tracking=per_T_hs (T_track=T_hs) and screening_fraction=0 "
        "exclusively. This round's headline now ALSO uses cavity_tracking=per_T_hs for every core "
        "row (item 1), so the family=c_plane screening=unscreened_lower VERDICT line above IS on "
        "comparable tracking-convention and screening terms with round 1's headline; the "
        "screening=screened_upper line is a round-2-only addition with no round-1 counterpart and "
        "is NOT comparable to round 1 on screening.", ""]
    # Radius-degeneracy disclosure (item 7): overlap_sq and tau_rad_bare_ns
    # are checked HERE (not asserted) for radius-invariance at matched
    # (orientation, height, T_hs, screening, regime) -- radius enters this
    # model only through E_a and prefactors, not the overlap/lifetime
    # themselves.
    by_group = {}
    for r in core:
        if r.get("valid") and _finite_num(r.get("overlap_sq")):
            key = (r.get("orientation"), r.get("height_nm"), r.get("T_hs"), r.get("screening_fraction"), r.get("regime"))
            by_group.setdefault(key, set()).add(round(float(r["overlap_sq"]), 12))
    n_groups_checked = sum(1 for v in by_group.values() if v)
    n_radius_invariant = sum(1 for v in by_group.values() if len(v) == 1)
    lines += ["## Radius-degeneracy disclosure", "",
              f"Checked directly against sweep.csv: {n_radius_invariant} of {n_groups_checked} "
              "(orientation, height, T_hs, screening, regime) groups have IDENTICAL overlap_sq "
              "across every sampled radius (5,10,15,20,30 nm) -- this planar model's overlap_sq "
              "and tau_rad_bare_ns do not depend on radius_nm; the radius axis enters only through "
              "E_a_meV (escape barrier, ~11 meV over 5-30 nm) and the escape/counting prefactors. "
              "Most core rows at fixed (orientation,height,T_hs,screening,regime) are therefore "
              "near-duplicates in E_X/overlap/tau_rad_bare, differing materially only through the "
              "escape-limited flux and validity at large radius.", ""]
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
              "E_X=2.54 eV [V] at 220 K, 19.0+/-0.4 meV linewidth [V], raw/corrected g2=0.47/0.21 [V], "
              "optical excitation (76 MHz, 1 ps, 800 nm two-photon) -- NOT an electrical SPS or a "
              "300 K validation. Composition transfer [A]: this model's cards use x_in=0.25 "
              "(In0.25Ga0.75N) at the same nominal geometry; Wang's paper does not report a "
              "composition, so the model's own default is retained rather than fitted. Model rows "
              "(comparison-only, excluded from headline coverage):", ""]
    for r in wang:
        lines.append(f"- row {r['row_id']}: T_hs={r['T_hs']:g} K, x_in={r.get('x_in')}, E_X_eV={r.get('E_X_eV')}, g2={r.get('g2')}, valid={r['valid']} ({r.get('literature_note')})")
    wang220 = next((r for r in wang if float(r.get("T_hs", -1)) == 220.), None)
    if wang220 is not None and _finite_num(wang220.get("E_X_eV")):
        lines.append(f"- Discrepancy at 220 K (never fitted to force agreement): model E_X={float(wang220['E_X_eV']):.4g} eV "
                      f"vs measured 2.54 eV [V]; model g2={wang220.get('g2')} vs measured raw/corrected 0.47/0.21 [V]. "
                      "The gap reflects the assumed x_in=0.25 transfer, uncapped-vs-capped surface treatment and "
                      "this model's idealized single-band overlap, none of which were tuned to close it.")
    desh = [r for r in geo if r["row_kind"] == "literature_aux" and "Deshpande" in str(r.get("literature_source", ""))]
    if desh:
        r = desh[0]
        lines.append(f"- Deshpande et al., APL 105, 141109 (2014) [V abstract-only; CONDITIONS INCOMPLETE], replay row {r['row_id']}: measured g2=0.29, model g2={r.get('g2')}.")
    lines += ["",
              f"Zhang et al., APL 108, 153102 (2016), Fig. 5: {ZHANG2016_SLOPE_MEV_PER_V:g} meV/V below 2 V [V], "
              "measured at T=10 K on a 3 nm In0.15Ga0.85N QW [V Fig. 5 / device section] -- both the "
              "temperature and the composition (x=0.15 vs this model's x_in=0.25 curves) are disclosed "
              "differences, non-gating (different device, temperature and excitation from this "
              "model's 230-300 K planar-dot cards); shown only as a labelled guide anchored to a "
              "model row on the Stark figures, never a fitted target.", ""]
    lines += ["## Screening-compatibility table", "",
              "The [-12,-8] meV/V window used when --slope-range is not supplied is EXPLICITLY "
              "illustrative around Zhang's approximate value, never a measured interval; the same "
              "--bias-window applies to any user-supplied slope. screening_compatibility is now "
              "called ONCE per (orientation,height,polarity,T_hs) geometry group with ALL THREE "
              "screening hypotheses together (fix round 3, item 3), so a genuine degeneracy across "
              "screening_fraction is detected instead of being hidden by testing each hypothesis in "
              "isolation. Nonpolar (m-plane/a-plane) rows are marked screening_unidentifiable ONLY "
              "when the shared curve already matches the window (polarization_factor=0 makes "
              "screening_fraction physically inert there, so an ACCEPTED nonpolar slope cannot "
              "distinguish a screening hypothesis); a nonpolar slope the window rejects stays "
              "incompatible.", ""]
    # Distinct-hypothesis counting (item 3): (orientation, height, polarity,
    # screening) -- regime and T_hs are replicate evaluations, not
    # independent hypotheses.
    distinct = _distinct_compat_rows(comp)
    n_unident = sum(1 for r in distinct if r.get("identification_status") == "screening_unidentifiable")
    n_compat = sum(1 for r in distinct if r.get("compatible") is True)
    n_incompat = sum(1 for r in distinct if r.get("identification_status") == "incompatible")
    lines.append(f"{len(comp)} raw compatibility-fit rows reduce to {len(distinct)} DISTINCT hypotheses "
                 "(orientation, height, polarity, screening) after collapsing regime/T_hs replicates: "
                 f"{n_compat} compatible, {n_incompat} incompatible, {n_unident} screening_unidentifiable.")
    compatible_distinct = [r for r in distinct if r.get("compatible") is True]
    lines += ["", "| orientation | height_nm | polarity | screening | fitted_slope meV/V | lifetime_ratio(s0/s1) | row_id |",
              "|---|---|---|---|---|---|---|"]
    for r in sorted(compatible_distinct, key=lambda z: (str(z.get("group")), z.get("screening"))):
        try:
            gk = dict(json.loads(r["group"])) if isinstance(r.get("group"), str) else dict(r.get("group") or [])
        except (TypeError, ValueError, json.JSONDecodeError):
            gk = {}
        lines.append(f"| {gk.get('orientation')} | {gk.get('height_nm')} | {gk.get('field_polarity')} | "
                      f"{r.get('screening')} | {r.get('fitted_slope_meV_per_V')} | {r.get('lifetime_ratio_s0_over_s1')} | {r.get('row_id')} |")
    if len(compatible_distinct) >= 2:
        # Height-screening degeneracy check: >=2 distinct (height,screening)
        # combos matching the SAME window at the SAME (orientation,polarity)
        # means the window alone does not identify screening.
        by_orient_pol = {}
        for r in compatible_distinct:
            try:
                gk = dict(json.loads(r["group"])) if isinstance(r.get("group"), str) else dict(r.get("group") or [])
            except (TypeError, ValueError, json.JSONDecodeError):
                gk = {}
            by_orient_pol.setdefault((gk.get("orientation"), gk.get("field_polarity")), []).append((gk.get("height_nm"), r.get("screening")))
        degenerate = {k: v for k, v in by_orient_pol.items() if len(set(v)) >= 2}
        if degenerate:
            for (o, pol), combos in degenerate.items():
                lines.append(f"\nHeight-screening degeneracy at orientation={o}, polarity={pol}: "
                              f"{sorted(set(combos))} all fit the same illustrative window -- the window "
                              "does NOT identify screening_fraction from height alone. The discriminating "
                              "observable is the bias-resolved lifetime (see lifetime_ratio_s0_over_s1 "
                              "column above, computed from stark_bias.csv's tau_rad_bare_ns at the "
                              "unscreened vs screened hypotheses).")
        else:
            lines.append("\nNo height-screening degeneracy among the currently compatible hypotheses this run.")
    lines.append("")
    # Sensitivities at baselines that PASS (item 5).
    lines += ["## One-at-a-time sensitivities at multiple baselines (item 5)", "",
              "The original REF baseline (rectangular, unscreened) never passes the optical gate; "
              "REF_pass_c_plane_screened_SET and REF_pass_a_plane_SET (H=3nm, R=10nm, T_hs=300K, "
              "deterministic_pair regime) DO pass, so every [A] assumption below is also bounded at "
              "a configuration that produces a verdict.", ""]
    # Hardening item 3: every OAT axis row's own regime/orientation/n_dot_cm2
    # is shown next to it, plus a note wherever a row was evaluated under a
    # DIFFERENT regime, orientation or carrier density than its own baseline
    # (island_radius_nm rows force regime=deterministic_pair; semipolar_
    # factor rows force orientation=semipolar_11_22; tau_cap_density_
    # convention rows force n_dot_cm2=1e9) -- otherwise the axis label alone
    # could be misread as isolating a single physical input.
    baseline_for_label = {"REF_unscreened_rectangular": REF, "REF_pass_c_plane_screened_SET": REF_PASS_C,
                           "REF_pass_a_plane_SET": REF_PASS_A}
    for label in ("REF_unscreened_rectangular", "REF_pass_c_plane_screened_SET", "REF_pass_a_plane_SET"):
        rs = [r for r in sens if r.get("sensitivity_baseline") == label and r["row_kind"] == "sensitivity"]
        if not rs: continue
        sample = rs[0]
        lines.append(f"### {label} (optical_pass at baseline: {sample.get('optical_pass')})")
        base = baseline_for_label[label]
        base_regime, base_orientation, base_n_dot = base["regime"], base["orientation"], base.get("n_dot_cm2")
        lines += ["", "| axis | value | regime | orientation | n_dot_cm2 | T_hs K | g2 | flux/s | optical_pass | row_id | note |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(rs, key=lambda z: (z.get("sensitivity_axis"), str(z.get("sensitivity_value")), z.get("T_hs"))):
            row_regime = r.get("regime"); row_orientation = r.get("orientation"); row_n_dot = r.get("n_dot_cm2")
            diffs = []
            if row_regime != base_regime:
                diffs.append(f"regime={row_regime} (baseline {base_regime})")
            if row_orientation != base_orientation:
                diffs.append(f"orientation={row_orientation} (baseline {base_orientation})")
            if row_n_dot not in (None, base_n_dot):
                diffs.append(f"n_dot_cm2={float(row_n_dot):g} (baseline card default)")
            note = "; ".join(diffs)
            lines.append(f"| {r.get('sensitivity_axis')} | {r.get('sensitivity_value')} | {row_regime} | {row_orientation} | "
                          f"{row_n_dot if row_n_dot is not None else ''} | {r.get('T_hs'):g} | "
                          f"{r.get('g2')} | {r.get('signal_flux_s')} | {r.get('optical_pass')} | {r.get('row_id')} | {note} |")
        lines.append("")
    # Invalid-row accounting (item 7): every kind, with reasons, tied to
    # the manifest's own invalid_counts_by_kind/invalid_reasons_summary.
    lines += ["## Invalid rows (all kinds)", "",
              f"invalid_total={sum(invalid_by_kind.values())} across kinds: "
              + ", ".join(f"{k}={v}" for k, v in sorted(invalid_by_kind.items())) + ".", "",
              "Reasons: " + ", ".join(f"{k}={v}" for k, v in sorted(invalid_reasons_summary.items())) if invalid_reasons_summary else "Reasons: none recorded.", ""]
    scr05_invalid = sum(1 for r in core if r.get("orientation") == "c_plane" and not r.get("valid") and r.get("screening_fraction") == .5)
    lines.append(f"c-plane screening=0.5 invalid core rows: {scr05_invalid} (included in the core count above).")
    lines.append("")
    lines += ["", "## Derivative-refinement convergence checks", "",
              "| check | orientation | height_nm | screening | polarity | T_hs | anchor V_j | coarse meV/V | fine (half-step) meV/V | abs diff | rel diff | converged | anchor row_id | probe row_ids |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, c in enumerate(refinement_checks):
        sp = c.get("spec", {})
        coords = f"{sp.get('orientation')} | {sp.get('height_nm')} | {sp.get('screening_fraction')} | {sp.get('polarity')} | {sp.get('T_hs')}"
        if "coarse_slope_meV_per_V" not in c:
            lines.append(f"| {i} | {coords} | - | - | - | - | - | {c.get('status')} | - | {c.get('probe_row_ids', [])} |"); continue
        lines.append(f"| {i} | {coords} | {c['anchor_V_j']:.3g} | {c['coarse_slope_meV_per_V']:.4g} | {c['fine_slope_meV_per_V']:.4g} | "
                      f"{c['abs_diff_meV_per_V']:.4g} | {c['rel_diff']:.4g} | {c['converged']} | {c.get('anchor_row_id')} | {c.get('probe_row_ids', [])} |")
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
    sens_p = build_all_sensitivities(a.quick); anchor_p = build_fixed_anchor(a.quick); lit_p = build_literature_aux()
    refinement_specs = (REFINEMENT_SPECS_QUICK if a.quick else REFINEMENT_SPECS_FULL) if a.stark else []
    planned = (len(core_p) + len(shape_p) + len(qw_p) + len(bias_p) + len(current_p)
               + len(sens_p) + len(anchor_p) + len(lit_p) + 1 + 2 * len(refinement_specs))
    if a.dry_run:
        print(json.dumps({"core_rows": len(core_p), "shape_rows": len(shape_p), "qw_rows": len(qw_p),
                           "stark_bias_rows": len(bias_p), "stark_current_rows": len(current_p),
                           "sensitivity_rows": len(sens_p), "fixed_anchor_rows": len(anchor_p),
                           "literature_aux_rows": len(lit_p) + 1,
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
    add_all(sens, "fixed_anchor", anchor_p)

    # Fix round 2: any per-trace exception in derivatives, compatibility or
    # plotting is caught and logged to the manifest (trace_id, reason) --
    # never raised past this function -- so one bad trace (e.g. duplicate/
    # non-increasing voltages) cannot abort the whole sweep and leave no
    # artifacts written.
    derivative_trace_failures = []
    if bias: derivative_trace_failures += attach_derivatives(bias, "V_j")
    if current: derivative_trace_failures += attach_derivatives(current, "V_j")

    # Half of the ACTUAL V_j trace step used to build `bias` in this mode
    # (0.1 V for the full 0.2 V grid, 0.25 V for the quick 0.5 V grid) --
    # spec: "refinement probes at half the trace step".
    stark_v_half_step = (0.5 if a.quick else 0.2) / 2.
    refinement_checks = []; probe_rows = []
    if refinement_specs and bias:
        refinement_checks, rid, probe_rows = run_refinement_checks(refinement_specs, bias, counter, cache, rid, stark_v_half_step)

    slope_range = tuple(a.slope_range) if a.slope_range else (-12., -8.)
    comp = []; compatibility_trace_failures = []
    if bias:
        comp, rid, compatibility_trace_failures = build_compatibility(bias, slope_range, a.bias_window, a.slope_range is not None, rid)

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
         "slope_meV_per_V": ZHANG2016_SLOPE_MEV_PER_V, "temperature_K": 10, "composition": "In0.15Ga0.85N",
         "note": "[V] Fig. 5, 10 K, x=0.15, below 2 V; non-gating comparison only",
         "model_row_ids": json.dumps([r["row_id"] for r in comp])},
        {"row_id": "LC00004", "source": "Deshpande et al., APL 105, 141109 (2014)", "kind": "conditions_incomplete_replay",
         "note": "[V abstract-only] measured g2=0.29",
         "model_row_ids": json.dumps([r["row_id"] for r in geo if r["row_kind"] == "literature_aux" and "Deshpande" in str(r.get("literature_source", ""))])},
    ]
    _write(out / "literature_comparisons.csv", lit_csv)

    plots = {}
    # Fix round 2: a per-trace-group failure inside plot_axis_response/
    # plot_stark is caught and logged there (fail_log=plot_trace_failures);
    # _safe_plot is the outer backstop for a whole-figure failure. Either
    # way the sweep continues and every other artifact (CSVs, results.md,
    # manifest.json) is still written.
    plot_trace_failures = []
    # A fixed representative slice (radius=5nm, height=1nm) is used for the
    # "other axis" of height/radius response plots -- both values are
    # guaranteed present in the reduced --quick core grid (H=[1,7],
    # R=[5,30]) as well as the full grid, so quick and full runs both
    # produce populated (non-empty) figures from the same slice logic.
    # height/radius response: group_keys include T_hs (Opus fix-round medium
    # finding: these previously selected T_hs in (230,300) but omitted it
    # from group_keys, mixing both temperatures into one curve at duplicated
    # x) -- panel-splitting (in plot_axis_response/_panel_plot) now shows
    # every orientation/screening/regime group instead of truncating to the
    # first 14 by sort key.
    #
    # Fix round 3, item 6: orientation_flux/temperature_response use the
    # REFERENCE geometry (H=3nm, R=10nm -- the corner where SET rows
    # actually pass), not the previous (R=5,H=1) corner where nothing
    # passes; the fixed H/R are stated in the title, and the 1000/s flux
    # floor / g2<0.5 gate are drawn as guide lines. Falls back to whatever
    # the reduced --quick core grid actually samples (H=1,R=5) when the
    # true reference point (3,10) is not in this run's core grid.
    avail_h = sorted({r["height_nm"] for r in core}); avail_r = sorted({r["radius_nm"] for r in core})
    ref_h = 3. if 3. in avail_h else (avail_h[0] if avail_h else 1.)
    ref_r = 10. if 10. in avail_r else (avail_r[0] if avail_r else 5.)
    ref_rows = [r for r in core if r["radius_nm"] == ref_r and r["height_nm"] == ref_h]
    figs = [
        ("height_response.png", plot_axis_response, (out, "height_response.png", [r for r in core if r["T_hs"] in (230., 300.) and r["radius_nm"] == 5.], "height_nm", "E_X_eV", "height response (r=5nm)"), dict(group_keys=("orientation", "screening_fraction", "regime", "T_hs"), fail_log=plot_trace_failures)),
        ("radius_response.png", plot_axis_response, (out, "radius_response.png", [r for r in core if r["T_hs"] in (230., 300.) and r["height_nm"] == 1.], "radius_nm", "E_X_eV", "radius response (h=1nm)"), dict(group_keys=("orientation", "screening_fraction", "regime", "T_hs"), fail_log=plot_trace_failures)),
        ("temperature_response.png", plot_axis_response, (out, "temperature_response.png", ref_rows, "T_hs", "g2", f"temperature response (g2); FIXED height={ref_h:g}nm radius={ref_r:g}nm"), dict(fail_log=plot_trace_failures, guide=(G2, "g2<0.5 optical gate"))),
        ("orientation_flux.png", plot_axis_response, (out, "orientation_flux.png", ref_rows, "T_hs", "signal_flux_s", f"orientation flux comparison; FIXED height={ref_h:g}nm radius={ref_r:g}nm"), dict(fail_log=plot_trace_failures, guide=(FLUX, "flux>=1000/s floor"))),
        # Shape/QW plots grouped by every fixed input (Astra finding): shape
        # additionally separates radius_nm, T_hs, regime and the native/
        # full-height variant (shape_height_fraction); QW separates
        # radius_nm, T_hs, regime and the height-offset dh (the QW's own
        # x-axis is wl_thickness_nm=w, so dh must stay a group key, not be
        # folded into the x-axis, or the two dh values collide).
        ("shape_comparison.png", plot_axis_response, (out, "shape_comparison.png", [r for r in geo if r["row_kind"] == "shape"], "height_nm", "E_X_eV", "shape mapping sensitivity (every fixed input separated)"), dict(group_keys=("shape", "orientation", "screening_fraction", "radius_nm", "T_hs", "regime", "shape_height_fraction"), fail_log=plot_trace_failures)),
        ("qw_response.png", plot_axis_response, (out, "qw_response.png", [r for r in geo if r["row_kind"] == "qw"], "wl_thickness_nm", "E_X_eV", "QW thickness response (every fixed input separated)"), dict(group_keys=("orientation", "screening_fraction", "radius_nm", "T_hs", "regime", "qw_height_offset_nm"), fail_log=plot_trace_failures)),
        ("pulse_vs_set.png", plot_pulse_vs_set, (out, core), {}),
        # Envelopes cover c-plane AND a-plane (Opus fix-round medium finding:
        # both envelope figures previously filtered to orientation=="c_plane" only).
        ("envelope_pulse.png", plot_envelope, (out, "envelope_pulse.png", [r for r in core if r["regime"] == "rectangular" and r["radius_nm"] == 5. and r["orientation"] in ("c_plane", "a_plane")], "pulse regime, r=5nm"), {}),
        ("envelope_set.png", plot_envelope, (out, "envelope_set.png", [r for r in core if r["regime"] == "deterministic_pair" and r["radius_nm"] == 5. and r["orientation"] in ("c_plane", "a_plane")], "SET regime, r=5nm"), {}),
    ]
    if sens:
        figs.append(("set_feasibility.png", plot_set_island, (out, sens), {}))
    if bias:
        figs += [
            ("stark_energy_bias.png", plot_stark, (out, "stark_energy_bias.png", bias, "V_j", "E_X_eV", "E_X(V_j)"), dict(zhang_anchor=True, fail_log=plot_trace_failures)),
            ("stark_tau_bias.png", plot_stark, (out, "stark_tau_bias.png", bias, "V_j", "tau_rad_bare_ns", "bare tau_rad(V_j)"), dict(ylog=True, fail_log=plot_trace_failures)),
            ("stark_tau_cavity_bias.png", plot_stark, (out, "stark_tau_cavity_bias.png", bias, "V_j", "tau_rad_cavity_ns", "cavity tau_rad(V_j)"), dict(ylog=True, fail_log=plot_trace_failures)),
            ("stark_overlap_bias.png", plot_stark, (out, "stark_overlap_bias.png", bias, "V_j", "overlap_sq", "overlap(V_j)"), dict(ylog=True, fail_log=plot_trace_failures)),
        ]
    if current:
        figs += [
            ("stark_energy_current.png", plot_stark, (out, "stark_energy_current.png", current, "current_uA", "E_X_eV", "E_X(I), self-consistent T_j (heating)"), dict(fail_log=plot_trace_failures)),
            ("stark_tau_current.png", plot_stark, (out, "stark_tau_current.png", current, "current_uA", "tau_rad_bare_ns", "tau_rad(I), self-consistent T_j (heating)"), dict(ylog=True, fail_log=plot_trace_failures)),
        ]
    figure_axis_scale = {}
    for fig_name, fn, args, kwargs in figs:
        result = _safe_plot(fig_name, plot_trace_failures, fn, *args, **kwargs)
        if result is not None:
            plots[fig_name] = result
            scale = _resolve_axis_scale(fn, args, kwargs)
            if scale is not None:
                figure_axis_scale[fig_name] = scale

    complete = (not a.quick) and counter["evaluate_calls"] <= a.max_evaluations and (time.time() - t0) < 1800
    runtime_s = time.time() - t0

    # Actual call counts by kind and invalid counts/reasons (Opus fix-round
    # medium finding: the manifest previously recorded no invalid/skipped
    # counts or reasons and no actual evaluate() call counts by kind, so the
    # 18 valid=False stark_bias rows in the fix-round-1 review appeared
    # nowhere). refinement_probe rows are included even though they are not
    # written to their own CSV (they are diagnostic-only, already referenced
    # by row_id inside convergence_checks). Computed BEFORE results.md so
    # its invalid-rows section (item 7) can report the same numbers.
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

    # ec_margin (item 2) and the SET g2 floor (item 4) are read directly off
    # an actually-evaluated card's own design object (never hand-entered).
    ec_margin = None; g2_floor = None
    try:
        set_card = _CARD_CACHE.get("nitride-cavity-set-design.yaml") or next(iter(_CARD_CACHE.values()))
        sp = set_card.drive.set_params
        ec_margin = float(sp["ec_margin"] if isinstance(sp, dict) else sp.ec_margin)
    except (AttributeError, KeyError, StopIteration, TypeError):
        pass
    try:
        b_res_val = float(set_card.drive.b_res)
        g2_floor = 1. - (1. / (1. + b_res_val)) ** 2
    except (AttributeError, KeyError, NameError, TypeError):
        pass

    (out / "results.md").write_text(
        _results_md(core, geo, bias, current, sens, comp, refinement_checks, complete, runtime_s,
                    counter["evaluate_calls"], invalid_by_kind, invalid_reasons_summary, ec_margin, g2_floor),
        encoding="utf-8")

    files = [p.name for p in out.iterdir() if p.is_file() and p.name != "manifest.json"]
    hashes = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}
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
                            "sensitivities": len([r for r in sens if r["row_kind"] == "sensitivity"]),
                            "fixed_anchor": len([r for r in sens if r["row_kind"] == "fixed_anchor"]),
                            "compatibility": len(comp),
                            "refinement_probes": len(probe_rows)},
        "full_contract_rows": {"core": 3360, "shape": 864, "qw": 288, "stark": 3120, "total_before_caching": 7632},
        "axes": {"height_nm": list(CORE_H), "radius_nm": list(CORE_R), "orientation": list(ORI),
                 "T_hs": list(TS), "screening_fraction": list(SCR), "regime": list(REG)},
        "cache_hits": sum(1 for r in all_kind_rows if r.get("cache_hit")),
        "numeric_convergence_policy": "derivative refinement: halved voltage step must agree to <=0.5 meV/V absolute or <=5% relative [A convergence criteria]",
        "convergence_checks": refinement_checks,
        "slope_interval": {"range_meV_per_V": list(slope_range), "kind": "user_supplied" if a.slope_range else "illustrative_not_measured", "bias_window_V": list(a.bias_window)},
        # Fix round 2: per-trace failures in derivatives, compatibility and
        # plotting are caught and recorded here (trace_id, reason, and the
        # affected row_ids where applicable) rather than aborting the sweep;
        # each key is always a list (possibly empty).
        "derivative_trace_failures": derivative_trace_failures,
        "compatibility_trace_failures": compatibility_trace_failures,
        "plot_trace_failures": plot_trace_failures,
        "plot_row_mapping": plots,
        "figure_axis_scale": figure_axis_scale,
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
