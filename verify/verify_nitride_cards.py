"""Mechanical checker for the planar InGaN/GaN cavity design cards (piece 6,
nitride-cards-contract): cards/nitride-cavity-pulse-design.yaml, cards/
nitride-cavity-set-design.yaml, cards/nitride-deshpande2014-comparison-
design.yaml, docs/nitride_cavity_contract.md, and verify/data/
nitride_cavity_anchors.yaml.

Class discipline (README.md's five-way split): every check below is either
(N) numerical verification (round-trip/evaluate wiring, internal
consistency of the frozen grid, the SET idealized-vs-hardware distinction)
or (T) source transcription (the evidence ledger's own structure and the
literal anchor values this piece cites). No check here is (C) parameter
calibration or (P) held-out prediction; running a card through
fsim_core.device.evaluate() is a (N) wiring check on the shipped model, not
a re-derivation or re-verification of the underlying physics modules
(pieces 1-5, out of scope for this piece) and never a claim that a
physically favorable result is itself evidence of a working device.

Import-side-effect-free: no module-level filesystem writes, no temp
directories. YAML round-trips are done in memory by monkeypatching
Path.read_text for the duration of one DeviceDesign.load() call (mirroring
verify/verify_nitride_device.py's own check_yaml_round_trip), never by
writing a scratch file.
"""
from __future__ import annotations

import math
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_core.device as devmod
from fsim_core.device import DeviceDesign, evaluate

CARDS_DIR = ROOT / "cards"
PULSE_CARD = CARDS_DIR / "nitride-cavity-pulse-design.yaml"
SET_CARD = CARDS_DIR / "nitride-cavity-set-design.yaml"
COMPARISON_CARD = CARDS_DIR / "nitride-deshpande2014-comparison-design.yaml"
CONTRACT_DOC = ROOT / "docs" / "nitride_cavity_contract.md"
ANCHORS_PATH = ROOT / "verify" / "data" / "nitride_cavity_anchors.yaml"

HEADLINE_CARDS = (PULSE_CARD, SET_CARD)
ALL_CARDS = (PULSE_CARD, SET_CARD, COMPARISON_CARD)

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
    (mirrors verify/verify_nitride_device.py's own _encode)."""
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


# ------------------------------------------------------------------- path walk

def _get_path(mapping, dotted):
    """Look up a dotted path (numeric segments index into lists) in a raw
    (yaml.safe_load'd) design mapping."""
    cur = mapping
    for seg in dotted.split("."):
        if isinstance(cur, list):
            cur = cur[int(seg)]
        else:
            cur = cur[seg]
    return cur


def _iter_leaves(mapping, prefix=""):
    """Yield (dotted_path, value) for every scalar/null leaf under `mapping`
    (dicts and lists are descended, scalars/null/bool/str are leaves)."""
    if isinstance(mapping, dict):
        for k, v in mapping.items():
            path = f"{prefix}.{k}" if prefix else str(k)
            yield from _iter_leaves(v, path)
    elif isinstance(mapping, list):
        for i, v in enumerate(mapping):
            path = f"{prefix}.{i}" if prefix else str(i)
            yield from _iter_leaves(v, path)
    else:
        yield prefix, mapping


# ---------------------------------------------------------------- anchors

def _load_anchors():
    doc = _load_yaml(ANCHORS_PATH)
    return doc["anchors"]


REQUIRED_ANCHOR_FIELDS = ("citation", "doi_or_url", "location", "evidence_status",
                          "tag", "platform", "excitation", "observable", "value",
                          "unit", "tolerance", "transfer_notes")
VALID_EVIDENCE_STATUS = {"full_text", "abstract_only", "figure_reading", "estimate", "missing"}
VALID_TAGS = {"V", "DR", "E", "A"}

REQUIRED_ANCHOR_IDS = (
    "bernardini1997-gan-polarization", "rinke2008-gan-electron-mass",
    "wu2003-ingan-bowing", "tsai2020-inn-gan-vbo",
    "deshpande2013-dot-thickness", "deshpande2013-dot-diameter",
    "deshpande2013-x-in", "deshpande2013-varshni-alpha",
    "deshpande2013-varshni-beta", "deshpande2013-xx-splitting",
    "deshpande2014-x-in", "deshpande2014-temperature", "deshpande2014-g2",
    "deshpande2014-lifetime", "deshpande2014-max-rep-rate",
    "deshpande2014-geometry-missing", "deshpande2014-waveform-missing",
    "deshpande2014-count-rate-missing", "deshpande2014-g2-measurement-rate-missing",
    "taylor2010-cavity-q", "taylor2010-no-purcell-observed",
    "zhang2016-qw-thickness", "zhang2016-x-in", "zhang2016-doping-na",
    "zhang2016-doping-nd", "zhang2016-g2-raw", "zhang2016-g2-corrected",
    "seguin2006-gan-lo-phonon", "kako2006-max-triggered-temperature",
    "holmes2014-g2-300k", "tamariz2020-g2-300k", "tamariz2020-linewidth-300k",
)


def check_anchor_ledger_structure():
    anchors = _load_anchors()
    ok("anchors ledger is a non-empty mapping", isinstance(anchors, dict) and len(anchors) > 0)
    for aid in REQUIRED_ANCHOR_IDS:
        ok(f"required anchor present: {aid}", aid in anchors)
    for aid, entry in anchors.items():
        missing_fields = [f for f in REQUIRED_ANCHOR_FIELDS if f not in entry]
        ok(f"anchor {aid} has all required fields", not missing_fields)
        if missing_fields:
            continue
        ok(f"anchor {aid} evidence_status is valid", entry["evidence_status"] in VALID_EVIDENCE_STATUS)
        ok(f"anchor {aid} tag is valid", entry["tag"] in VALID_TAGS)
        tol = entry["tolerance"]
        ok(f"anchor {aid} tolerance names its own meaning",
          isinstance(tol, dict) and "meaning" in tol and "value" in tol)
        if entry["evidence_status"] == "missing":
            ok(f"anchor {aid} missing status has null value", entry["value"] is None)
        else:
            ok(f"anchor {aid} non-missing status has a non-null value", entry["value"] is not None)
    return anchors


def check_anchor_literal_values(anchors):
    """AC3: independent literal anchor checks -- Deshpande2013 x0.25@10K vs
    Deshpande2014 x0.40@300K/g2=0.29/max200MHz, distinct raw/corrected
    Zhang values, optical-only GaN/AlN comparisons; missing 2014 fields
    stay null (never a fabricated number)."""
    d13 = anchors["deshpande2013-x-in"]
    ok("Deshpande2013 x_in literal value is 0.25", d13["value"] == 0.25)
    ok("Deshpande2013 x_in is full_text (10 K nanowire)", d13["evidence_status"] == "full_text")
    ok("Deshpande2013 x_in platform is the nanowire dot-in-wire device",
      "nanowire" in d13["platform"].lower())

    d14x = anchors["deshpande2014-x-in"]
    d14T = anchors["deshpande2014-temperature"]
    d14g2 = anchors["deshpande2014-g2"]
    d14rep = anchors["deshpande2014-max-rep-rate"]
    ok("Deshpande2014 x_in literal value is 0.40", d14x["value"] == 0.40)
    ok("Deshpande2014 temperature literal value is 300 K", d14T["value"] == 300.0)
    ok("Deshpande2014 g2 literal value is 0.29", d14g2["value"] == 0.29)
    ok("Deshpande2014 max repetition rate literal value is 200 MHz", d14rep["value"] == 2.0e8)
    for a in (d14x, d14T, d14g2, d14rep):
        ok(f"Deshpande2014 anchor {a['observable']} is abstract_only", a["evidence_status"] == "abstract_only")
    ok("Deshpande2013 x_in (0.25) != Deshpande2014 x_in (0.40): distinct compositions", d13["value"] != d14x["value"])

    # Missing 2014 fields stay null -- never a fabricated number.
    for aid in ("deshpande2014-geometry-missing", "deshpande2014-waveform-missing",
               "deshpande2014-count-rate-missing", "deshpande2014-g2-measurement-rate-missing"):
        a = anchors[aid]
        ok(f"{aid} is recorded missing with a null value (not fabricated)",
          a["evidence_status"] == "missing" and a["value"] is None)

    zraw = anchors["zhang2016-g2-raw"]; zcorr = anchors["zhang2016-g2-corrected"]
    ok("Zhang2016 raw g2 literal value is 0.42", zraw["value"] == 0.42)
    ok("Zhang2016 corrected g2 literal value is 0.38", zcorr["value"] == 0.38)
    ok("Zhang2016 raw and corrected g2 are distinct (never conflated)", zraw["value"] != zcorr["value"])

    holmes = anchors["holmes2014-g2-300k"]; tamariz = anchors["tamariz2020-g2-300k"]
    ok("Holmes2014 is an optical-pumping anchor (external context)", "optical" in holmes["excitation"].lower())
    ok("Holmes2014 platform is not the headline InGaN/GaN electrical platform",
      "InGaN" not in holmes["platform"] or "electrical" not in holmes["excitation"].lower())
    ok("Tamariz2020 is an optical-pumping anchor (external context)", "optical" in tamariz["excitation"].lower())
    ok("Tamariz2020 platform is GaN/AlN, not InGaN/GaN", "GaN" in tamariz["platform"] and "InGaN" not in tamariz["platform"])

    taylor_q = anchors["taylor2010-cavity-q"]; taylor_np = anchors["taylor2010-no-purcell-observed"]
    ok("Taylor2010 measured cavity Q literal value is 167", taylor_q["value"] == 167)
    ok("Taylor2010 reports no observed Purcell enhancement (value 0)", taylor_np["value"] == 0)


# ---------------------------------------------------------- provenance checks

def check_card_provenance(path, anchors):
    label = path.stem
    raw = _raw_design(path)
    provenance = raw.get("provenance", {})
    sources = provenance.get("sources", {})
    ok(f"[{label}] provenance.sources is present and non-empty", isinstance(sources, dict) and len(sources) > 0)

    # Every explicit leaf under `design` (excluding `provenance` and `name`)
    # must have a provenance.sources entry with the required keys.
    design_for_leaves = {k: v for k, v in raw.items() if k not in ("provenance", "name")}
    missing_sources = []
    value_mismatches = []
    for path_str, value in _iter_leaves(design_for_leaves):
        if path_str not in sources:
            missing_sources.append(path_str)
            continue
        entry = sources[path_str]
        if entry.get("value") != value:
            value_mismatches.append((path_str, entry.get("value"), value))
    ok(f"[{label}] every explicit design leaf has a provenance.sources entry", not missing_sources)
    if missing_sources:
        print("      missing:", missing_sources[:20])
    ok(f"[{label}] every provenance.sources value matches its card leaf", not value_mismatches)
    if value_mismatches:
        print("      mismatched:", value_mismatches[:10])

    required_keys = {"value", "unit", "tag", "source"}
    incomplete = [k for k, v in sources.items() if not required_keys <= set(v)]
    ok(f"[{label}] every provenance.sources entry has value/unit/tag/source", not incomplete)
    bad_tags = [k for k, v in sources.items() if v.get("tag") not in VALID_TAGS]
    ok(f"[{label}] every provenance.sources tag is one of V/DR/E/A", not bad_tags)

    # anchor_id cross-reference: every cited anchor must exist in the ledger.
    dangling = [k for k, v in sources.items() if v.get("anchor_id") and v["anchor_id"] not in anchors]
    ok(f"[{label}] every anchor_id resolves in the evidence ledger", not dangling)
    if dangling:
        print("      dangling anchor_ids:", dangling)

    # assumptions <-> A/E tag bidirectional consistency.
    assumptions = set(provenance.get("assumptions", []))
    ae_tagged = {k for k, v in sources.items() if v.get("tag") in ("A", "E")}
    ok(f"[{label}] design.provenance.assumptions lists exactly the A/E-tagged leaves",
      assumptions == ae_tagged)
    if assumptions != ae_tagged:
        print("      assumptions-only:", sorted(assumptions - ae_tagged)[:10])
        print("      A/E-only:", sorted(ae_tagged - assumptions)[:10])

    return raw, provenance


# ------------------------------------------------------------------ round trip

def check_round_trip(path):
    """AC1: load/save/reload without changing values, using an in-memory
    YAML round trip (monkeypatch Path.read_text for one load() call, no
    temp files/dirs) -- then evaluate both the original and the
    round-tripped design and require bit-identical scalars."""
    label = path.stem
    d = DeviceDesign.load(str(path))
    doc_str = yaml.safe_dump({"meta": {"name": d.name}, "design": asdict(d)}, sort_keys=False)
    reloaded_raw = yaml.safe_load(doc_str)["design"]

    orig_read_text = Path.read_text
    Path.read_text = lambda self, *a, **kw: doc_str
    try:
        d2 = DeviceDesign.load("<in-memory round-trip of " + str(path) + ">")
    finally:
        Path.read_text = orig_read_text

    ok(f"[{label}] round-trip preserves platform", d2.platform == d.platform)
    ok(f"[{label}] round-trip preserves nitride block", d2.nitride == d.nitride)
    ok(f"[{label}] round-trip preserves provenance block", d2.provenance == d.provenance)
    ok(f"[{label}] round-trip preserves dot/ret/drive/thermal/cavity/emission/aperture blocks",
      d2.dot == d.dot and d2.ret == d.ret and d2.drive == d.drive and d2.thermal == d.thermal
      and d2.cavity == d.cavity and d2.emission == d.emission and d2.aperture == d.aperture)

    s1 = evaluate(d, [d.thermal.T_hs])["scalars"]
    s2 = evaluate(d2, [d.thermal.T_hs])["scalars"]
    ok(f"[{label}] round-trip evaluates bit-identically", _encode(s1) == _encode(s2))

    # File-round-trip half of AC1 ("if a file round-trip is needed rewrite
    # only the same new card and restore identical bytes, never create a
    # temporary directory"): re-serializing this card's own already-loaded
    # design must reproduce the SAME field values the file on disk holds,
    # checked purely in memory against the raw YAML -- no file is written.
    raw_on_disk = _raw_design(path)
    ok(f"[{label}] in-memory re-serialization matches the on-disk nitride block",
      reloaded_raw.get("nitride") == raw_on_disk.get("nitride"))
    return d, s1


# ------------------------------------------------------------- common geometry

SET_ONLY_TOP_KEYS = ("name",)
SET_ONLY_DRIVE_KEYS = ("cycle_loading", "eta_load", "set_params")


def check_common_geometry():
    """AC1: exact common geometry/optics/material parameters between the two
    headline regimes -- only name, cycle_loading, and SET pricing metadata
    (drive.eta_load, drive.set_params) may differ."""
    pulse = _raw_design(PULSE_CARD)
    setc = _raw_design(SET_CARD)

    def strip(raw):
        raw = dict(raw)
        raw.pop("name", None)
        raw.pop("provenance", None)
        drive = dict(raw["drive"])
        for k in SET_ONLY_DRIVE_KEYS:
            drive.pop(k, None)
        raw["drive"] = drive
        return raw

    ok("headline cards share identical geometry/optics/material parameters "
      "(only name, cycle_loading, and SET pricing metadata differ)",
      strip(pulse) == strip(setc))
    ok("pulse card uses cycle_loading=rectangular", pulse["drive"]["cycle_loading"] == "rectangular")
    ok("SET card uses cycle_loading=deterministic_pair", setc["drive"]["cycle_loading"] == "deterministic_pair")
    ok("SET island radius is independent of the QD radius",
      setc["drive"]["set_params"]["radius_nm"] != setc["nitride"]["dot"]["radius_nm"])


# --------------------------------------------------------------------- timing

def check_timing_and_pumping():
    for path in ALL_CARDS:
        label = path.stem
        raw = _raw_design(path)
        drive = raw["drive"]
        ok(f"[{label}] finite_pulse is true", drive["finite_pulse"] is True)
        ok(f"[{label}] no optical pumping (drive.cw=false)", drive["cw"] is False)
        ok(f"[{label}] mechanism is empty (no drive_mech override)", not drive.get("mechanism", ""))

    for path in HEADLINE_CARDS:
        label = path.stem
        raw = _raw_design(path)
        drive = raw["drive"]
        ok(f"[{label}] 100 ps pump pulse", drive["diode"]["tau_pulse_ns"] == 0.1)
        ok(f"[{label}] 80 MHz repetition rate", drive["rep_rate_hz"] == 8.0e7)
        ok(f"[{label}] gate_ns not requested (resolves to the full period)", drive["gate_ns"] is None)
        period = 1e9 / drive["rep_rate_hz"]
        ok(f"[{label}] resolved period is 12.5 ns", abs(period - 12.5) < 1e-9)
        d = DeviceDesign.load(str(path))
        s = evaluate(d, [d.thermal.T_hs])["scalars"]
        ok(f"[{label}] counting gate resolves to the full 12.5 ns period, not the 100 ps pump pulse",
          abs(s["gate_ns_used"] - period) < 1e-9)

    comp_raw = _raw_design(COMPARISON_CARD)
    ok("[comparison] 200 MHz repetition rate (Deshpande 2014's reported maximum)",
      comp_raw["drive"]["rep_rate_hz"] == 2.0e8)
    ok("[comparison] rectangular (finite electrical) loading, distinct from the headline grid",
      comp_raw["drive"]["cycle_loading"] == "rectangular")


# ---------------------------------------------------------------- frozen grid

EXPECTED_GRID = {
    "nitride.dot.height_nm": [1.0, 2.0, 3.0, 4.0, 5.0],
    "nitride.dot.radius_nm": [5.0, 10.0, 15.0],
    "nitride.dot.x_in": [0.15, 0.25, 0.40],
    "thermal.T_hs": [230.0, 250.0, 273.0, 300.0],
    "nitride.cavity.Q": [500.0, 2000.0, 10000.0],
    "drive.I_uA": [0.002, 0.02, 0.2],
}
EXPECTED_HEADLINE_ROW_COUNT = 5 * 3 * 3 * 4 * 3 * 3 * 2  # 3240 (x2 loading regimes)

# path -> resolver against a raw design mapping ("cycle_loading" lives under drive)
_PATH_RESOLVERS = {
    "cycle_loading": lambda raw: raw["drive"]["cycle_loading"],
}


def _resolve(raw, path):
    if path in _PATH_RESOLVERS:
        return _PATH_RESOLVERS[path](raw)
    return _get_path(raw, path)


def check_frozen_grid():
    ok("frozen headline grid produces exactly 3240 rows (5x3x3x4x3x3x2)",
      EXPECTED_HEADLINE_ROW_COUNT == 3240)

    for path in HEADLINE_CARDS:
        label = path.stem
        raw = _raw_design(path)
        ranges = raw["provenance"]["ranges"]
        for axis, expected in EXPECTED_GRID.items():
            ok(f"[{label}] ranges['{axis}'] matches the frozen headline grid axis",
              ranges.get(axis, {}).get("values") == expected)
        ok(f"[{label}] provenance.headline_grid_row_count is 3240",
          raw["provenance"].get("headline_grid_row_count") == 3240)
        # Each card's own explicit value is a member of every headline axis
        # it participates in (the fixed-reference convention).
        for axis in EXPECTED_GRID:
            actual = _resolve(raw, axis)
            ok(f"[{label}] own value for '{axis}' is a member of the frozen grid axis",
              actual in EXPECTED_GRID[axis])
        ok(f"[{label}] own cycle_loading is a member of the frozen loading-regime axis",
          _resolve(raw, "cycle_loading") in ranges["cycle_loading"]["values"])

    comp_raw = _raw_design(COMPARISON_CARD)
    ok("[comparison] is explicitly NOT part of the 3240-row headline grid",
      comp_raw["provenance"].get("headline_grid_row_count") == 0)
    ok("[comparison] own radius_nm (12.5) is deliberately OUTSIDE the headline radius_nm axis",
      _resolve(comp_raw, "nitride.dot.radius_nm") not in EXPECTED_GRID["nitride.dot.radius_nm"])
    for axis in ("nitride.dot.height_nm", "nitride.dot.x_in", "thermal.T_hs",
                "nitride.cavity.Q", "drive.I_uA"):
        actual = _resolve(comp_raw, axis)
        ok(f"[comparison] own value for '{axis}' still coincides with a frozen grid point",
          actual in EXPECTED_GRID[axis])

    # SET island-radius sensitivity axis, independent of the headline grid.
    set_raw = _raw_design(SET_CARD)
    island_range = set_raw["provenance"]["ranges"]["drive.set_params.radius_nm"]["values"]
    ok("SET island-radius sensitivity axis is exactly [0.5, 1, 5] nm", island_range == [0.5, 1.0, 5.0])
    ok("SET island radius (5.0 nm) is a member of its own sensitivity axis",
      set_raw["drive"]["set_params"]["radius_nm"] in island_range)


REQUIRED_HEADINGS = (
    "## Platform and geometry", "## Electrical regimes", "## Cavity and tracking",
    "## Counting and background", "## Evidence status", "## Sweep grid and cost",
    "## Acceptance gates", "## Deshpande comparison", "## Artifacts and reproduction",
    "## Known limitations",
)


def check_contract_doc():
    text = CONTRACT_DOC.read_text(encoding="utf-8")
    for heading in REQUIRED_HEADINGS:
        ok(f"contract doc has heading '{heading}'", heading in text)
    for literal in ("height_nm=[1,2,3,4,5]", "radius_nm=[5,10,15]", "x_in=[0.15,0.25,0.40]",
                    "T_hs=[230,250,273,300]", "Q=[500,2000,10000]",
                    "current_uA=[0.002,0.02,0.2]", "3240"):
        ok(f"contract doc states grid literal '{literal}'", literal in text)
    ok("contract doc is ASCII", all(ord(c) < 128 for c in text))
    ok("contract doc names the exact test-command verifier",
      "verify_nitride_cards.py" in text)


# ------------------------------------------------------------------- evaluate

REQUIRED_SCALAR_KEYS = (
    "platform", "cycle_loading", "T_hs", "T_j", "g2_op", "collected_flux_pulsed_s",
    "rho_pulsed", "valid", "invalid_reasons", "gate_ns_used", "rep_rate_hz",
    "one_pair_valid", "set_feasible", "set_priced_F_p", "device_pass",
    "pair_supply_possible", "mean_counts",
)


def check_evaluate_each_card():
    """AC4: evaluate each card with T_grid=[thermal.T_hs]; scalar modes/
    units/assumptions are consistent. Bound/bright/pure results are not
    card-validity gates: a physical failure is valid output with reasons.
    The deterministic row exposes BOTH idealized output and hardware
    diagnostic, and cannot report a hardware pass when infeasible."""
    for path in ALL_CARDS:
        label = path.stem
        d = DeviceDesign.load(str(path))
        out = evaluate(d, [d.thermal.T_hs])
        s = out["scalars"]
        ok(f"[{label}] evaluate() returns the full required scalar contract",
          set(REQUIRED_SCALAR_KEYS) <= set(s))
        ok(f"[{label}] platform tag is ingan_gan_planar", s["platform"] == "ingan_gan_planar")
        ok(f"[{label}] T_hs threads through unchanged", s["T_hs"] == d.thermal.T_hs)
        ok(f"[{label}] valid is a bool", isinstance(s["valid"], (bool, np.bool_)))
        if s["valid"]:
            ok(f"[{label}] a valid row reports a finite g2_op", math.isfinite(s["g2_op"]))
            ok(f"[{label}] a valid row reports a finite collected flux",
              math.isfinite(s["collected_flux_pulsed_s"]))
        else:
            ok(f"[{label}] an invalid row still names its reasons (a physical failure is valid output)",
              len(s["invalid_reasons"]) > 0)
        # A physically unfavorable OR favorable result is never itself
        # treated as a pass/fail gate by this checker -- only the wiring
        # (finiteness/shape) is verified here; the acceptance sweep (piece
        # 7) owns the actual PASS/FAIL verdict.

    d_set = DeviceDesign.load(str(SET_CARD))
    s_set = evaluate(d_set, [d_set.thermal.T_hs])["scalars"]
    ok("[nitride-cavity-set-design] cycle_loading threads through as deterministic_pair",
      s_set["cycle_loading"] == "deterministic_pair")
    ok("[nitride-cavity-set-design] idealized counting (mean_counts) is exposed regardless of hardware feasibility",
      math.isfinite(s_set["mean_counts"]))
    ok("[nitride-cavity-set-design] a hardware-infeasible SET cannot report device_pass",
      not (s_set["device_pass"] and not s_set["set_feasible"]))
    ok("[nitride-cavity-set-design] device_pass requires one_pair_valid, set_feasible AND pair_supply_possible together",
      not s_set["device_pass"] or (s_set["one_pair_valid"] and s_set["set_feasible"] and s_set["pair_supply_possible"]))

    d_pulse = DeviceDesign.load(str(PULSE_CARD))
    s_pulse = evaluate(d_pulse, [d_pulse.thermal.T_hs])["scalars"]
    ok("[nitride-cavity-pulse-design] cycle_loading threads through as rectangular",
      s_pulse["cycle_loading"] == "rectangular")
    ok("[nitride-cavity-pulse-design] rectangular regime reports one_pair_valid unset (False, SET-only concept)",
      s_pulse["one_pair_valid"] is False)


def check_deshpande_comparison_card():
    label = COMPARISON_CARD.stem
    raw = _raw_design(COMPARISON_CARD)
    provenance = raw["provenance"]
    ok(f"[{label}] comparison_status is conditions_incomplete",
      provenance.get("comparison_status") == "conditions_incomplete")
    ok(f"[{label}] transfer_assumptions is a non-empty explicit list",
      isinstance(provenance.get("transfer_assumptions"), list) and len(provenance["transfer_assumptions"]) >= 5)

    anchors = _load_anchors()
    comp = provenance["deshpande_comparison"]
    ok(f"[{label}] deshpande_comparison.measured_g2 matches the ledger (0.29)",
      comp["measured_g2"]["value"] == anchors["deshpande2014-g2"]["value"] == 0.29)
    ok(f"[{label}] deshpande_comparison.measured_dot_size is null (unpublished, never fabricated)",
      comp["measured_dot_size"]["value"] is None)
    ok(f"[{label}] deshpande_comparison.measured_count_rate is null (unpublished, never fabricated as zero)",
      comp["measured_count_rate"]["value"] is None)
    ok(f"[{label}] deshpande_comparison.measured_g2_repetition_rate is null (max rate != g2 measurement rate)",
      comp["measured_g2_repetition_rate"]["value"] is None)

    # Never present 2013 (10 K) as the 2014 (300 K) device.
    ok(f"[{label}] nitride.dot.x_in (0.40) is the 2014 value, not 2013's 0.25",
      raw["nitride"]["dot"]["x_in"] == 0.40)
    ok(f"[{label}] thermal.T_hs (300 K) is the 2014 value, not 2013's 10 K",
      raw["thermal"]["T_hs"] == 300.0)
    ok(f"[{label}] height/radius placeholders are explicitly sourced from the 2013 geometry, tagged A",
      raw["provenance"]["sources"]["nitride.dot.height_nm"]["tag"] == "A"
      and raw["provenance"]["sources"]["nitride.dot.radius_nm"]["tag"] == "A")

    # This card must never gate the headline verdict: it carries no
    # headline grid membership and no acceptance-gate field of its own.
    ok(f"[{label}] never claims headline grid membership", provenance.get("headline_grid_row_count") == 0)


# ---------------------------------------------------------------------- ascii

def check_ascii_and_encoding():
    for path in ALL_CARDS + (ANCHORS_PATH,):
        text = path.read_text(encoding="utf-8")
        ok(f"{path.name} is pure ASCII", all(ord(c) < 128 for c in text))


def main():
    anchors = check_anchor_ledger_structure()
    check_anchor_literal_values(anchors)

    for path in ALL_CARDS:
        check_card_provenance(path, anchors)

    for path in ALL_CARDS:
        check_round_trip(path)

    check_common_geometry()
    check_timing_and_pumping()
    check_frozen_grid()
    check_contract_doc()
    check_evaluate_each_card()
    check_deshpande_comparison_card()
    check_ascii_and_encoding()

    print(f"{sum(checks)}/{len(checks)} nitride card checks passed")
    raise SystemExit(0 if all(checks) else 1)


if __name__ == "__main__":
    main()
