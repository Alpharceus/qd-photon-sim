"""Mechanical checker for the nonpolar a-plane and QW-thickness-fluctuation
geometry card pairs (piece 5, nitride-geometry-stark-contract):
cards/nitride-nonpolar-pulse-design.yaml, cards/nitride-nonpolar-set-design.
yaml, cards/nitride-qw-fluctuation-pulse-design.yaml, cards/
nitride-qw-fluctuation-set-design.yaml, docs/nitride_geometry_stark_contract
.md, and verify/data/nitride_geometry_stark_anchors.yaml.

Class discipline (README.md's five-way split): every check below is either
(N) numerical verification (DeviceDesign.load/evaluate wiring, the
one-variable structural diff against each card's c-plane parent, the
polarization_factor-contradiction scratch test, screening-degeneracy) or
(T) source transcription (the evidence ledger's own structure and literal
anchor values, checked against literals typed independently of the ledger
file itself). No check here is (C) parameter calibration or (P) held-out
prediction; evaluating a card through fsim_core.device.evaluate() is a (N)
wiring check on the shipped model (pieces 1-4, out of scope for this
piece), never a re-derivation of the underlying physics and never a claim
that a physically favorable result is itself evidence of a working device.

Import-side-effect-free: no module-level filesystem writes, no temp
directories. The round-trip check monkeypatches Path.read_text for the
duration of one DeviceDesign.load() call (mirroring verify/
verify_nitride_cards.py's own check_round_trip), never writing a scratch
file, so a sandboxed/read-only filesystem cannot make this check raise
PermissionError.
"""
from __future__ import annotations

import copy
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, evaluate
from fsim_core.nitride_levels import NitrideDotSystem, levels

CARDS_DIR = ROOT / "cards"
NONPOLAR_PULSE = CARDS_DIR / "nitride-nonpolar-pulse-design.yaml"
NONPOLAR_SET = CARDS_DIR / "nitride-nonpolar-set-design.yaml"
QW_PULSE = CARDS_DIR / "nitride-qw-fluctuation-pulse-design.yaml"
QW_SET = CARDS_DIR / "nitride-qw-fluctuation-set-design.yaml"
CPLANE_PULSE_PARENT = CARDS_DIR / "nitride-cavity-pulse-design.yaml"
CPLANE_SET_PARENT = CARDS_DIR / "nitride-cavity-set-design.yaml"

NEW_CARDS = (NONPOLAR_PULSE, NONPOLAR_SET, QW_PULSE, QW_SET)
# (new card, its c-plane parent)
PARENT_OF = {
    NONPOLAR_PULSE: CPLANE_PULSE_PARENT,
    NONPOLAR_SET: CPLANE_SET_PARENT,
    QW_PULSE: CPLANE_PULSE_PARENT,
    QW_SET: CPLANE_SET_PARENT,
}

CONTRACT_DOC = ROOT / "docs" / "nitride_geometry_stark_contract.md"
ANCHORS_PATH = ROOT / "verify" / "data" / "nitride_geometry_stark_anchors.yaml"
LEGACY_CARDS_VERIFIER = ROOT / "verify" / "verify_nitride_cards.py"

checks = []


def ok(name, value):
    checks.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


# ------------------------------------------------------------------ loading

def _load_yaml(path):
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def _raw_design(path):
    return _load_yaml(path)["design"]


def _encode(obj):
    """float.hex/array-exact encoding for bit-identical scalar comparison
    (mirrors verify/verify_nitride_cards.py's own _encode)."""
    if isinstance(obj, dict):
        return {k: _encode(v) for k, v in sorted(obj.items())}
    if isinstance(obj, np.ndarray):
        return [_encode(v) for v in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [_encode(v) for v in obj]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float.hex(float(obj))
    if obj is None or isinstance(obj, str):
        return obj
    raise TypeError(f"cannot encode {type(obj)!r}")


def _flatten(value, prefix=""):
    """Dotted-path flatten of a nested dict (mirrors _iter_leaves in
    verify/verify_nitride_cards.py but keyed as a dict of dotted-path ->
    leaf value; lists are treated as opaque leaves, not descended into, so
    a changed list element shows up as one differing key)."""
    out = {}
    if isinstance(value, dict):
        for k, v in value.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.update(_flatten(v, key))
    else:
        out[prefix] = value
    return out


# --------------------------------------------------------- structural diff

# The ONE-VARIABLE-DIFFERENCE allow-list per new card: every top-level
# design key except 'provenance' (documentation, not physics) is flattened
# and diffed against the c-plane parent; every dotted key that differs
# (changed value) or is added (present only on the new card) must appear
# here, or the check fails. No key may be REMOVED relative to the parent
# (a card that silently drops a parent leaf is never "verbatim copy plus
# intended fields").
ALLOWED_DIFF = {
    NONPOLAR_PULSE: {
        "name",
        "nitride.dot.orientation",
        "nitride.dot.polarization_factor",
    },
    NONPOLAR_SET: {
        "name",
        "nitride.dot.orientation",
        "nitride.dot.polarization_factor",
    },
    QW_PULSE: {
        "name",
        "nitride.dot.height_nm",
        "nitride.dot.radius_nm",
        "nitride.dot.wl_thickness_nm",
        "nitride.dot.geometry_type",
        "nitride.dot.shape",
        "drive.diode.d_i_nm",
        "drive.diode.d_active_nm",
        "drive.diode.wl_thickness_nm",
    },
    QW_SET: {
        "name",
        "nitride.dot.height_nm",
        "nitride.dot.radius_nm",
        "nitride.dot.wl_thickness_nm",
        "nitride.dot.geometry_type",
        "nitride.dot.shape",
        "drive.diode.d_i_nm",
        "drive.diode.d_active_nm",
        "drive.diode.wl_thickness_nm",
    },
}


def check_structural_diff():
    """AC2: one-variable-difference diff of the loaded design dicts,
    restricted to the intended keys. Proves the orientation/geometry
    change is isolated and that density, aperture, background, lifetime,
    and SET hardware parameters have not been opportunistically improved
    on top of the intended change."""
    for new_path, parent_path in PARENT_OF.items():
        label = new_path.stem
        d_new = asdict(DeviceDesign.load(str(new_path)))
        d_parent = asdict(DeviceDesign.load(str(parent_path)))
        d_new.pop("provenance", None)
        d_parent.pop("provenance", None)
        flat_new = _flatten(d_new)
        flat_parent = _flatten(d_parent)
        allowed = ALLOWED_DIFF[new_path]

        removed = set(flat_parent) - set(flat_new)
        ok(f"[{label}] no parent leaf is silently dropped", not removed)

        changed_or_added = set()
        for k, v in flat_new.items():
            if k not in flat_parent or flat_parent[k] != v:
                changed_or_added.add(k)
        unexpected = changed_or_added - allowed
        ok(f"[{label}] every differing/added leaf vs. {parent_path.stem} is on the intended-change allow-list",
           not unexpected)
        if unexpected:
            print("     unexpected diff keys: " + ", ".join(sorted(unexpected)))
        # And the allow-listed keys actually DO differ/exist (a stale
        # allow-list entry that matches nothing is a spec error, not a
        # passing check).
        ok(f"[{label}] every allow-listed key is actually present on the new card",
           allowed <= set(flat_new))


# ------------------------------------------------------------- card loading

def check_cards_load_and_evaluate():
    """AC1: every new card loads through fsim_core.device.DeviceDesign.load
    (never raw yaml.safe_load) and evaluates at 230 K and 300 K in both
    regimes, with explicit screening_fraction in {0, 0.5, 1}. Invalid rows
    are retained with non-empty invalid_reasons; a schema error (an
    exception from load()/evaluate()) is a hard failure, never relabelled
    as expected physics invalidity."""
    for path in NEW_CARDS:
        label = path.stem
        try:
            design = DeviceDesign.load(str(path))
            load_ok = True
        except Exception as exc:  # noqa: BLE001 - report, don't hide
            load_ok = False
            print(f"     [{label}] DeviceDesign.load raised: {exc!r}")
        ok(f"[{label}] loads through fsim_core.device.DeviceDesign.load", load_ok)
        if not load_ok:
            continue

        for screening in (0.0, 0.5, 1.0):
            d2 = copy.deepcopy(design)
            d2.nitride["dot"]["screening_fraction"] = screening
            try:
                result = evaluate(d2, T_grid=[230.0, 300.0])
                eval_ok = True
            except Exception as exc:  # noqa: BLE001
                eval_ok = False
                print(f"     [{label}] screening={screening} evaluate() raised: {exc!r}")
            ok(f"[{label}] screening={screening}: evaluate() at [230,300] K does not raise",
               eval_ok)
            if not eval_ok:
                continue
            s = result["scalars"]
            if s["valid"]:
                finite_keys = ("g2_op", "E_X_eV", "field_kVcm",
                               "collected_flux_pulsed_s", "S_X", "overlap_sq")
                all_finite = all(np.isfinite(s[k]) for k in finite_keys)
                ok(f"[{label}] screening={screening}: valid row has finite {finite_keys}",
                   all_finite)
            else:
                ok(f"[{label}] screening={screening}: invalid row retains non-empty invalid_reasons",
                   isinstance(s["invalid_reasons"], list) and len(s["invalid_reasons"]) > 0)


# --------------------------------------------------------- screening-degen

def check_nonpolar_screening_degenerate():
    """AC3: nonpolar (a_plane, polarization_factor=0) spectra are
    screening-degenerate at fixed bias -- the normal polarization field's
    intrinsic component is zero regardless of screening_fraction, so
    E_X_eV and field_kVcm must not move across screening in {0, 0.5, 1}."""
    for path in (NONPOLAR_PULSE, NONPOLAR_SET):
        label = path.stem
        design = DeviceDesign.load(str(path))
        e_x = []
        field = []
        for screening in (0.0, 0.5, 1.0):
            d2 = copy.deepcopy(design)
            d2.nitride["dot"]["screening_fraction"] = screening
            s = evaluate(d2, T_grid=[300.0])["scalars"]
            e_x.append(s["E_X_eV"])
            field.append(s["field_kVcm"])
        ok(f"[{label}] E_X_eV is screening-degenerate (0/0.5/1) at fixed bias",
           len(set(e_x)) == 1 and np.isfinite(e_x[0]))
        ok(f"[{label}] field_kVcm is screening-degenerate (0/0.5/1) at fixed bias",
           len(set(field)) == 1 and np.isfinite(field[0]))


# ------------------------------------------------- polarization-factor test

def check_polarization_factor_contradiction():
    """Required fix-round check: a card's polarization_factor must never
    silently contradict its orientation. Tried in a scratch copy (never
    the on-disk card): orientation=a_plane with an explicit
    polarization_factor=1.0 must make nitride_levels.levels() report
    valid=False with the module's own contradiction reason, not raise an
    uncaught exception and not silently succeed."""
    scratch = NitrideDotSystem(height_nm=3.0, radius_nm=10.0, x_in=0.25,
                                orientation="a_plane", polarization_factor=1.0)
    lv = levels(scratch, 300.0)
    ok("scratch orientation=a_plane, polarization_factor=1.0 is invalid",
       lv.valid is False)
    ok("scratch contradiction reports the module's own reason string",
       "polarization_factor contradicts the selected orientation" in lv.invalid_reasons)

    # The consistent case (explicit 0.0, matching a_plane's own factor)
    # must NOT raise the contradiction reason (it may still be valid=False
    # for unrelated physical-binding reasons; only the contradiction
    # reason is asserted absent).
    consistent = NitrideDotSystem(height_nm=3.0, radius_nm=10.0, x_in=0.25,
                                   orientation="a_plane", polarization_factor=0.0)
    lv_ok = levels(consistent, 300.0)
    ok("scratch orientation=a_plane, polarization_factor=0.0 raises no contradiction",
       "polarization_factor contradicts the selected orientation" not in lv_ok.invalid_reasons)


# --------------------------------------------------------------- QW wiring

def check_qw_reservoir_wiring():
    """AC3: QW cards expose the QW-fluctuation confinement model with a
    matched surrounding well (drive.diode.wl_thickness_nm ==
    nitride.dot.wl_thickness_nm, drive.diode.x_in == nitride.dot.x_in,
    drive.diode.d_i_nm >= nitride.dot.height_nm), no explicit
    nitride.background.reservoir_energy_eV override, and a reservoir_kind
    actually sourced from the qw_fluctuation confinement branch (not
    silently inherited from the isolated_dot default). nitride_levels'
    _qw_levels reports reservoir_kind per-carrier-channel-selection
    ('ingan_qw' when both carriers select the surrounding well,
    'gan_barrier' otherwise) -- either is a legitimate physically
    determined output for this card's own spec-mandated H=3.5/R=20/w=3 nm
    geometry; this check requires the QW branch was actually taken
    (geometry_type == 'qw_fluctuation') rather than requiring one specific
    channel-selection outcome."""
    for path in (QW_PULSE, QW_SET):
        label = path.stem
        design = DeviceDesign.load(str(path))
        dot = design.nitride["dot"]
        diode = design.drive.diode
        ok(f"[{label}] nitride.dot.geometry_type is qw_fluctuation",
           dot["geometry_type"] == "qw_fluctuation")
        ok(f"[{label}] nitride.dot.shape is disc (the only qw_fluctuation shape)",
           dot["shape"] == "disc")
        ok(f"[{label}] 0 < nitride.dot.wl_thickness_nm < nitride.dot.height_nm",
           0.0 < dot["wl_thickness_nm"] < dot["height_nm"])
        ok(f"[{label}] drive.diode.wl_thickness_nm == nitride.dot.wl_thickness_nm (matched well)",
           diode["wl_thickness_nm"] == dot["wl_thickness_nm"])
        ok(f"[{label}] drive.diode.x_in == nitride.dot.x_in (matched well)",
           diode["x_in"] == dot["x_in"])
        ok(f"[{label}] drive.diode.d_i_nm >= nitride.dot.height_nm",
           diode["d_i_nm"] >= dot["height_nm"])
        background = design.nitride.get("background", {})
        ok(f"[{label}] no nitride.background.reservoir_energy_eV override",
           "reservoir_energy_eV" not in background)

        s = evaluate(design, T_grid=[300.0])["scalars"]
        ok(f"[{label}] scalars['geometry_type'] is qw_fluctuation (QW branch actually ran)",
           s["geometry_type"] == "qw_fluctuation")
        ok(f"[{label}] scalars['reservoir_kind'] is a legitimate QW-branch value",
           s["reservoir_kind"] in ("ingan_qw", "gan_barrier"))
        print(f"     [{label}] reservoir_kind at 300 K, screening=0: {s['reservoir_kind']!r}")


def check_nonpolar_orientation_wiring():
    for path in (NONPOLAR_PULSE, NONPOLAR_SET):
        label = path.stem
        dot = DeviceDesign.load(str(path)).nitride["dot"]
        ok(f"[{label}] nitride.dot.orientation is a_plane", dot["orientation"] == "a_plane")
        ok(f"[{label}] nitride.dot.polarization_factor is explicit 0.0",
           dot["polarization_factor"] == 0.0)
        ok(f"[{label}] nitride.dot.geometry_type is isolated_dot (unchanged from parent)",
           dot.get("geometry_type", "isolated_dot") == "isolated_dot")


# ---------------------------------------------------- timing / T_track / SET

def check_timing_and_T_track():
    for path in NEW_CARDS:
        label = path.stem
        design = DeviceDesign.load(str(path))
        ok(f"[{label}] drive.diode.tau_pulse_ns is 100 ps",
           design.drive.diode["tau_pulse_ns"] == 0.1)
        ok(f"[{label}] drive.rep_rate_hz is 80 MHz",
           design.drive.rep_rate_hz == 8.0e7)
        ok(f"[{label}] nitride.cavity.T_track is fixed at 300 K",
           design.nitride["cavity"]["T_track"] == 300.0)

    ok("[nonpolar-pulse] cycle_loading is rectangular",
       DeviceDesign.load(str(NONPOLAR_PULSE)).drive.cycle_loading == "rectangular")
    ok("[nonpolar-set] cycle_loading is deterministic_pair",
       DeviceDesign.load(str(NONPOLAR_SET)).drive.cycle_loading == "deterministic_pair")
    ok("[qw-fluctuation-pulse] cycle_loading is rectangular",
       DeviceDesign.load(str(QW_PULSE)).drive.cycle_loading == "rectangular")
    ok("[qw-fluctuation-set] cycle_loading is deterministic_pair",
       DeviceDesign.load(str(QW_SET)).drive.cycle_loading == "deterministic_pair")


def check_independent_set_radius():
    """The SET island radius is independent of the optical dot radius and
    is never silently resized by a geometry sweep."""
    for path, dot_radius in ((NONPOLAR_SET, 10.0), (QW_SET, 20.0)):
        label = path.stem
        design = DeviceDesign.load(str(path))
        island_r = design.drive.set_params["radius_nm"]
        ok(f"[{label}] SET island radius (5.0) differs from the optical dot radius ({dot_radius})",
           island_r == 5.0 and dot_radius != island_r)
        ok(f"[{label}] SET island radius matches its c-plane parent's own 5.0 nm (unimproved)",
           island_r == DeviceDesign.load(str(PARENT_OF[path])).drive.set_params["radius_nm"] == 5.0)

        # A geometry sweep mutating nitride.dot.radius_nm must never
        # silently resize the SET island.
        d2 = copy.deepcopy(design)
        d2.nitride["dot"]["radius_nm"] = 30.0
        s = evaluate(d2, T_grid=[300.0])["scalars"]
        ok(f"[{label}] mutating nitride.dot.radius_nm does not change the reported set_radius_nm",
           s["set_radius_nm"] == 5.0)


# ----------------------------------------------------------- provenance

def check_provenance_completeness():
    for path in NEW_CARDS:
        label = path.stem
        raw = _raw_design(path)
        provenance = raw["provenance"]
        assumptions = provenance.get("assumptions")
        ok(f"[{label}] provenance.assumptions is a non-empty explicit list",
           isinstance(assumptions, list) and len(assumptions) > 0)
        sources = provenance.get("sources", {})
        missing_sources = [k for k in (assumptions or []) if k not in sources]
        ok(f"[{label}] every assumptions entry has a matching provenance.sources entry",
           not missing_sources)
        if missing_sources:
            print("     missing sources for: " + ", ".join(missing_sources))
        bad_tags = [k for k, v in sources.items() if v.get("tag") not in ("V", "DR", "E", "A")]
        ok(f"[{label}] every provenance.sources entry has tag in {{V,DR,E,A}}",
           not bad_tags)
        incomplete = [k for k, v in sources.items()
                      if not v.get("source") or "unit" not in v or "value" not in v]
        ok(f"[{label}] every provenance.sources entry has source/unit/value",
           not incomplete)
        ok(f"[{label}] provenance.card_id matches meta.name",
           provenance.get("card_id") == raw["name"] == _load_yaml(path)["meta"]["name"])

    # The two orientation-specific sources this fix round adds.
    for path in (NONPOLAR_PULSE, NONPOLAR_SET):
        sources = _raw_design(path)["provenance"]["sources"]
        ok(f"[{path.stem}] provenance.sources documents nitride.dot.orientation",
           "nitride.dot.orientation" in sources and sources["nitride.dot.orientation"]["tag"] == "A")
        ok(f"[{path.stem}] provenance.sources documents nitride.dot.polarization_factor",
           "nitride.dot.polarization_factor" in sources)
        ok(f"[{path.stem}] the nonpolar mass/strain/valence-ordering caveat is present",
           "Schade" in sources["nitride.dot.orientation"]["source"]
           and "valence" in sources["nitride.dot.orientation"]["source"].lower())

    for path in (QW_PULSE, QW_SET):
        sources = _raw_design(path)["provenance"]["sources"]
        for key in ("nitride.dot.geometry_type", "nitride.dot.shape",
                    "drive.diode.d_active_nm", "drive.diode.wl_thickness_nm"):
            ok(f"[{path.stem}] provenance.sources documents {key}", key in sources)


# ----------------------------------------------------------- round trip

def check_round_trip(path):
    """AC6: round-trip the new card through the actual schema without
    losing nitride geometry/bias/provenance blocks, using an in-memory
    YAML round trip (monkeypatch Path.read_text for one load() call, no
    temp files/dirs -- mirrors verify/verify_nitride_cards.py's own
    check_round_trip) -- never overwriting the on-disk card fixture."""
    label = path.stem
    d = DeviceDesign.load(str(path))
    doc_str = yaml.safe_dump({"meta": {"name": d.name}, "design": asdict(d)}, sort_keys=False)

    orig_read_text = Path.read_text
    Path.read_text = lambda self, *a, **kw: doc_str
    try:
        d2 = DeviceDesign.load("<in-memory round-trip of " + str(path) + ">")
    finally:
        Path.read_text = orig_read_text

    ok(f"[{label}] round-trip preserves platform", d2.platform == d.platform)
    ok(f"[{label}] round-trip preserves nitride block (geometry/orientation/bias)",
       d2.nitride == d.nitride)
    ok(f"[{label}] round-trip preserves provenance block", d2.provenance == d.provenance)
    ok(f"[{label}] round-trip preserves dot/ret/drive/thermal/cavity/emission/aperture blocks",
       d2.dot == d.dot and d2.ret == d.ret and d2.drive == d.drive and d2.thermal == d.thermal
       and d2.cavity == d.cavity and d2.emission == d.emission and d2.aperture == d.aperture)

    s1 = evaluate(d, [d.thermal.T_hs])["scalars"]
    s2 = evaluate(d2, [d.thermal.T_hs])["scalars"]
    ok(f"[{label}] round-trip evaluates bit-identically", _encode(s1) == _encode(s2))


# ----------------------------------------------------------------- ascii

def check_ascii():
    for path in NEW_CARDS + (CONTRACT_DOC, ANCHORS_PATH):
        text = path.read_text(encoding="utf-8")
        ok(f"{path.name} is pure ASCII", all(ord(c) < 128 for c in text))


# ------------------------------------------------------------- contract doc

REQUIRED_CONTRACT_LITERALS = (
    "3,360 rows",
    "<= 10,000",
    "semipolar_11_22",
    "m_plane",
    "a_plane",
    "H=3.5/R=20/w=3 nm",
    "height_nm=7.0", "radius_nm=17.5",
    "reservoir_kind",
    "non_gating_comparison",
    "source_transcription",
    "missing_evidence",
    "dE_X/dV",
    "screening_fraction",
    "verify/verify_nitride_geometry_cards.py",
    "T_track",
)


def check_contract_doc():
    text = CONTRACT_DOC.read_text(encoding="utf-8")
    for literal in REQUIRED_CONTRACT_LITERALS:
        ok(f"contract doc states '{literal}'", literal in text)
    ok("contract doc names its own verifier",
       "verify_nitride_geometry_cards.py" in text)
    ok("contract doc references the accepted field table (pieces 1-4)",
       "Accepted field table" in text)
    ok("contract doc has a per-card reservoir/SRH table",
       "Reservoir and SRH convention" in text
       and "nitride-qw-fluctuation-pulse-design.yaml" in text
       and "nitride-qw-fluctuation-set-design.yaml" in text
       and "nitride-nonpolar-pulse-design.yaml" in text
       and "nitride-nonpolar-set-design.yaml" in text)


# ----------------------------------------------------------------- ledger

REQUIRED_ANCHOR_FIELDS = (
    "citation", "doi_or_url", "location", "evidence_status", "evidence_kind",
    "tag", "platform", "excitation", "observable", "value", "unit",
    "tolerance", "transfer_notes",
)
ALLOWED_EVIDENCE_STATUS = {"full_text", "abstract_only", "figure_reading", "missing"}
ALLOWED_EVIDENCE_KIND = {"source_transcription", "non_gating_comparison", "missing_evidence"}

# Independently-typed literals (NOT read back from the ledger file) for the
# source-transcription cross-check -- ../_goal/nitride_digests.md Sec. 7d/
# 7e/8, transcribed by this verifier a second time, independently of
# whatever the ledger file itself says.
WANG_LINEWIDTH_220K_MEV = 19.0
WANG_LINEWIDTH_220K_TOL_MEV = 0.4
WANG_G2_RAW_220K = 0.47
WANG_G2_CORR_220K = 0.21
WANG_LIFETIME_FAST_220K_PS = 357.0
WANG_LIFETIME_FAST_TOL_PS = 20.0
WANG_LIFETIME_FAST_4P7K_PS = 480.0
WANG_LIFETIME_SLOW_NS_APPROX = 4.0
WANG_EMISSION_220K_EV = 2.54
WANG_GEOMETRY_HEIGHT_NM = 7.0
WANG_GEOMETRY_DIAMETER_NM = 35.0
WANG_FSS_200K_MIN_MEV = 2.0
WANG_FSS_200K_MAX_MEV = 12.0
WANG_FSS_N_DOTS = 16
ZHANG_SLOPE_MEV_PER_V = -10.0


def _load_anchors():
    return _load_yaml(ANCHORS_PATH)["anchors"]


def check_ledger_schema(anchors):
    ok("ledger has exactly 8 anchors", len(anchors) == 8)
    for anchor_id, entry in anchors.items():
        missing = [f for f in REQUIRED_ANCHOR_FIELDS if f not in entry]
        ok(f"[{anchor_id}] has every required ledger field", not missing)
        ok(f"[{anchor_id}] evidence_status is one of {sorted(ALLOWED_EVIDENCE_STATUS)}",
           entry.get("evidence_status") in ALLOWED_EVIDENCE_STATUS)
        ok(f"[{anchor_id}] evidence_kind is one of {sorted(ALLOWED_EVIDENCE_KIND)}",
           entry.get("evidence_kind") in ALLOWED_EVIDENCE_KIND)
        ok(f"[{anchor_id}] tag is one of V/DR/E/A", entry.get("tag") in ("V", "DR", "E", "A"))
        doi = entry.get("doi_or_url")
        ok(f"[{anchor_id}] doi_or_url is null or a non-empty string",
           doi is None or (isinstance(doi, str) and len(doi) > 0))

    ok("required entry wang2017_fss (APL 111, 053101) is present", "wang2017_fss" in anchors)
    ok("required entry zhang2016_stark_slope is present", "zhang2016_stark_slope" in anchors)
    ok("required entry schade2011_orientation is present", "schade2011_orientation" in anchors)

    # No fabricated DOI: only identifiers actually printed in the digest
    # (arXiv ids) for the Wang/Zhang entries; Schade's DOI is explicitly
    # sourced from the design brief, not the digest (recorded in its own
    # transfer_notes), and is the one allowed https://doi.org/ literal.
    for anchor_id, entry in anchors.items():
        if anchor_id == "schade2011_orientation":
            continue
        doi = entry.get("doi_or_url")
        ok(f"[{anchor_id}] doi_or_url is a digest-printed arXiv id or null, never a fabricated DOI",
           doi is None or doi.startswith("arXiv:"))
    ok("schade2011_orientation doi_or_url is the design-brief-given DOI",
       anchors["schade2011_orientation"]["doi_or_url"] == "https://doi.org/10.1002/pssb.201046350")


def check_ledger_literal_values(anchors):
    wang_lw = anchors["wang2017_linewidth"]["value"]
    ok("wang2017_linewidth 220K value matches an independently-typed literal (19.0 meV)",
       wang_lw["meV_220K"] == WANG_LINEWIDTH_220K_MEV)
    ok("wang2017_linewidth 220K tolerance matches an independently-typed literal (+/-0.4 meV)",
       anchors["wang2017_linewidth"]["tolerance"]["meV_220K"] == WANG_LINEWIDTH_220K_TOL_MEV)

    wang_g2 = anchors["wang2017_g2"]["value"]
    ok("wang2017_g2 raw (220K) matches an independently-typed literal (0.47)",
       wang_g2["g2_raw_220K"] == WANG_G2_RAW_220K)
    ok("wang2017_g2 corrected (220K) matches an independently-typed literal (0.21)",
       wang_g2["g2_corr_220K"] == WANG_G2_CORR_220K)
    ok("wang2017_g2 raw != corrected (never conflated)",
       wang_g2["g2_raw_220K"] != wang_g2["g2_corr_220K"])
    ok("wang2017_g2 excitation is optical, never electrical",
       "optical" in anchors["wang2017_g2"]["excitation"].lower()
       and "electr" not in anchors["wang2017_g2"]["excitation"].lower())

    wang_life = anchors["wang2017_lifetime"]["value"]
    ok("wang2017_lifetime fast component (220K) matches an independently-typed literal (357 ps)",
       wang_life["tau_fast_ps_220K"] == WANG_LIFETIME_FAST_220K_PS)
    ok("wang2017_lifetime fast component (4.7K) matches an independently-typed literal (480 ps)",
       wang_life["tau_fast_ps_4p7K"] == WANG_LIFETIME_FAST_4P7K_PS)
    ok("wang2017_lifetime slow component matches an independently-typed literal (~4 ns)",
       wang_life["tau_slow_ns_approx"] == WANG_LIFETIME_SLOW_NS_APPROX)
    ok("wang2017_lifetime transfer_notes explicitly disclaims fitting tau_rad0_ns to it",
       "never fit tau_rad0_ns" in str(anchors["wang2017_lifetime"]["transfer_notes"]).lower())

    wang_em = anchors["wang2017_emission"]["value"]
    ok("wang2017_emission 220K value matches an independently-typed literal (2.54 eV)",
       wang_em["E_eV_220K"] == WANG_EMISSION_220K_EV)
    ok("wang2017_emission never relabels 220 K as 300 K",
       "E_eV_300K" not in wang_em and "300" not in str(anchors["wang2017_emission"]["location"]))

    wang_geom = anchors["wang2017_geometry"]["value"]
    ok("wang2017_geometry height matches an independently-typed literal (7 nm)",
       wang_geom["height_nm"] == WANG_GEOMETRY_HEIGHT_NM)
    ok("wang2017_geometry diameter matches an independently-typed literal (35 nm)",
       wang_geom["diameter_nm"] == WANG_GEOMETRY_DIAMETER_NM)

    wang_fss = anchors["wang2017_fss"]["value"]
    ok("wang2017_fss range matches independently-typed literals (2-12 meV at 200 K)",
       wang_fss["FSS_meV_200K_min"] == WANG_FSS_200K_MIN_MEV
       and wang_fss["FSS_meV_200K_max"] == WANG_FSS_200K_MAX_MEV)
    ok("wang2017_fss dot count matches an independently-typed literal (16 dots)",
       wang_fss["n_dots"] == WANG_FSS_N_DOTS)
    ok("wang2017_fss evidence_kind is non_gating_comparison (omitted-physics caution)",
       anchors["wang2017_fss"]["evidence_kind"] == "non_gating_comparison")

    zhang = anchors["zhang2016_stark_slope"]
    ok("zhang2016_stark_slope matches an independently-typed literal (-10.0 meV/V)",
       zhang["value"]["slope_meV_per_V"] == ZHANG_SLOPE_MEV_PER_V)
    ok("zhang2016_stark_slope evidence_kind is non_gating_comparison",
       zhang["evidence_kind"] == "non_gating_comparison")
    ok("zhang2016_stark_slope excitation records optical-under-bias (not pure electrical EL)",
       "optical" in zhang["excitation"].lower())

    schade = anchors["schade2011_orientation"]
    ok("schade2011_orientation tag is A, never V (not a numeric measurement transcription)",
       schade["tag"] == "A")
    ok("schade2011_orientation value is null (no number is transcribed from it)",
       schade["value"] is None)
    ok("schade2011_orientation evidence_kind is missing_evidence",
       schade["evidence_kind"] == "missing_evidence")
    ok("schade2011_orientation transfer_notes disclaims a polarization-factor measurement",
       "not a numeric" in schade["transfer_notes"].lower()
       or "not" in schade["transfer_notes"].lower() and "polarization-factor measurement" in schade["transfer_notes"].lower())


# ------------------------------------------------------- existing suite

def check_existing_nitride_cards_suite_still_passes():
    """AC5: verify/verify_nitride_cards.py (>= 377 checks) keeps passing
    unchanged -- run as a fresh subprocess so this piece's own imports and
    monkeypatching cannot leak state into it."""
    result = subprocess.run(
        [sys.executable, str(LEGACY_CARDS_VERIFIER)],
        cwd=str(ROOT), capture_output=True, text=True,
    )
    ok("verify/verify_nitride_cards.py exits 0", result.returncode == 0)
    last_line = ""
    for line in result.stdout.splitlines():
        if "nitride card checks passed" in line:
            last_line = line
    passed_n = None
    if last_line:
        try:
            passed_n = int(last_line.split("/")[0].strip())
        except ValueError:
            passed_n = None
    ok("verify/verify_nitride_cards.py reports >= 377 passed checks",
       passed_n is not None and passed_n >= 377)
    print("     " + (last_line or "(no summary line found)"))
    if result.returncode != 0:
        print("     ---- verify_nitride_cards.py stdout (tail) ----")
        print("\n".join(result.stdout.splitlines()[-20:]))
        print("     ---- verify_nitride_cards.py stderr (tail) ----")
        print("\n".join(result.stderr.splitlines()[-20:]))


# --------------------------------------------------------------------- main

def main():
    check_cards_load_and_evaluate()
    check_structural_diff()
    check_nonpolar_screening_degenerate()
    check_polarization_factor_contradiction()
    check_nonpolar_orientation_wiring()
    check_qw_reservoir_wiring()
    check_timing_and_T_track()
    check_independent_set_radius()
    check_provenance_completeness()
    for path in NEW_CARDS:
        check_round_trip(path)
    check_ascii()
    check_contract_doc()

    anchors = _load_anchors()
    check_ledger_schema(anchors)
    check_ledger_literal_values(anchors)

    check_existing_nitride_cards_suite_still_passes()

    print(f"{sum(checks)}/{len(checks)} nitride geometry card checks passed")
    raise SystemExit(0 if all(checks) else 1)


if __name__ == "__main__":
    main()
