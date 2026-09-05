"""Mechanical checks for the room-temperature edge-emitter design cards
(cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml) against
docs/rt_edge_contract.md and verify/data/rt_edge_anchors.yaml.

This file does NOT judge whether either card's device is physically good
(a card whose corner is not eligible, or whose g2 is bad, is still a valid
card -- see docs/rt_edge_contract.md "Acceptance gates" vs. this file's own
job). It checks that the cards are self-describing and honest:

  1. both cards load/save/reload bit-for-bit through
     fsim_core.device.DeviceDesign (exact round trip);
  2. both explicitly opt into every new rt_edge_contract.md mode, and never
     hand-set a quantity the evaluator derives (drive.b_e, cavity.beta_sin);
  3. every scalar leaf the card actually sets carries a
     design.provenance.sources entry (value, tag, source) -- see
     "scalar leaf" below for exactly what counts;
  4. every anchor_id used cross-checks against
     verify/data/rt_edge_anchors.yaml (missing-status anchors force an
     E/A tag and an assumptions-list entry; verified numeric anchors force
     the card value inside the anchor's tolerance);
  5. the three required provenance.ranges entries exist with the exact
     endpoints the spec fixes, and each one keyed by an actual design path
     brackets the value the card sets;
  6. the single-dot-selection arithmetic (n_expected =
     density_cm2 * pi * (diameter_um/2)^2 * 1e-8) is reproduced
     independently and is < 0.5;
  7. fsim_core.device.evaluate() runs at the card's own thermal.T_hs
     without raising, returns curves/scalars, and the resolved scalars
     report every opted-in mode as active.

"Scalar leaf": every non-dict, non-list value found by recursing through
the card's `design` mapping (excluding `design.provenance` itself, which is
evidence bookkeeping, not a device input), EXCEPT values reached through a
key literally named "kind" (a materials-constructor routing tag, e.g.
{kind: GaAsP, x_p: 0.4} -- selects which constructor to call, not itself a
physical quantity) or "name" (a human-readable label, e.g. a thermal-layer
name or the top-level design name). Everything else -- every number, every
opt-in mode string, every boolean -- is a leaf and needs a sources entry.

Standalone, side-effect-free on import (all work happens under
`if __name__ == "__main__"`), importable without executing the checks.
Exits 0 iff every check passes; prints "N/N rt-edge card checks passed".
"""
from __future__ import annotations

import dataclasses
import math
import os
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, evaluate  # noqa: E402

CARDS = [
    ROOT / "cards" / "edge-inp-gaasp-design.yaml",
    ROOT / "cards" / "edge-inp-gainp-design.yaml",
]
ANCHORS_PATH = ROOT / "verify" / "data" / "rt_edge_anchors.yaml"

# Required opt-ins, common to every card (docs/rt_edge_contract.md).
REQUIRED_OPTINS = {
    "dot.linewidth": "anchored",
    "ret.mode": "confinement",
    "drive.mode": "EL-transport",
    "emission.type": "edge",
    "aperture.compose": True,
    "drive.cw": True,
    "cavity.enabled": False,
    "thermal.T_hs": 300.0,
}

# The three ranges the spec requires, with their exact endpoints.
REQUIRED_RANGES = {
    "dot.delta_xx": (4.0, 7.0, "meV"),
    "dot.gamma300": (6.0, 20.0, "meV"),
    "irf_ps": (50.0, 200.0, "ps"),
}

VALID_TAGS = {"V", "DR", "E", "A"}
STRUCTURAL_KEYS = {"kind", "name"}


# ------------------------------------------------------------------ plumbing

CHECKS = []


def ok(name: str, value: bool) -> None:
    CHECKS.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


def leaves(node, prefix: str = ""):
    """Yield (dotted_path, value) for every scalar leaf under `node` --
    see the module docstring's "Scalar leaf" definition."""
    if isinstance(node, dict):
        for k, v in node.items():
            if k in STRUCTURAL_KEYS:
                continue
            path = f"{prefix}.{k}" if prefix else str(k)
            yield from leaves(v, path)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from leaves(v, f"{prefix}.{i}")
    else:
        yield prefix, node


def values_equal(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)
    return a == b


def anchor_bound(anchor: dict) -> float:
    """Absolute tolerance bound: a plain number is absolute; a {rel: x} or
    {abs: x} mapping is relative-to-the-anchor-value or absolute resp."""
    tol = anchor["tolerance"]
    if isinstance(tol, dict):
        if "rel" in tol:
            return abs(float(anchor["value"])) * float(tol["rel"])
        if "abs" in tol:
            return float(tol["abs"])
        raise ValueError(f"unsupported tolerance mapping {tol!r}")
    return float(tol)


def resolve_design_path(design: DeviceDesign, path: str):
    """Best-effort getattr(getattr(design, block), field) for a two-segment
    dotted path; returns (True, value) if `path` names an actual
    DeviceDesign scalar field, else (False, None) -- used to decide whether
    a provenance.ranges key is a design path at all (e.g. "irf_ps" is not)."""
    if "." not in path:
        return False, None
    block_name, field_name = path.split(".", 1)
    try:
        block = getattr(design, block_name)
        value = getattr(block, field_name)
    except AttributeError:
        return False, None
    if isinstance(value, (dict, list)):
        return False, None
    return True, value


# --------------------------------------------------------------------- main

def check_card(path: Path, anchors: dict) -> None:
    tag = path.name
    raw_doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    raw_design = raw_doc["design"]
    provenance = raw_design.get("provenance", {})

    design = DeviceDesign.load(path)

    # ---- 1. exact round trip through save()/load()
    with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        design.save(tmp_path)
        reloaded = DeviceDesign.load(tmp_path)
    finally:
        os.unlink(tmp_path)
    ok(f"{tag}: exact DeviceDesign round trip", dataclasses.asdict(design) == dataclasses.asdict(reloaded))

    # ---- 2. required opt-ins, and the two forbidden derived-quantity keys
    for optin_path, expected in REQUIRED_OPTINS.items():
        resolvable, value = resolve_design_path(design, optin_path)
        ok(f"{tag}: opt-in {optin_path} == {expected!r}", resolvable and value == expected)
    diode_preset = raw_design.get("drive", {}).get("diode", {}).get("preset", "")
    ok(f"{tag}: opt-in drive.diode.preset is set", bool(diode_preset))
    track_material = raw_design.get("filter", {}).get("track_material", "")
    ok(f"{tag}: opt-in filter.track_material is non-empty", track_material in ("dot", "matrix"))
    ok(f"{tag}: raw YAML omits derived drive.b_e", "b_e" not in raw_design.get("drive", {}))
    ok(f"{tag}: raw YAML omits derived cavity.beta_sin", "beta_sin" not in raw_design.get("cavity", {}))

    # ---- 3. provenance coverage: every scalar leaf has a sources entry
    sources = provenance.get("sources", {})
    design_for_leaves = {k: v for k, v in raw_design.items() if k != "provenance"}
    leaf_paths = list(leaves(design_for_leaves))
    ok(f"{tag}: provenance.sources is non-empty", len(sources) > 0)
    for leaf_path, leaf_value in leaf_paths:
        entry = sources.get(leaf_path)
        if entry is None:
            ok(f"{tag}: provenance.sources[{leaf_path!r}] exists", False)
            continue
        value_ok = values_equal(entry.get("value"), leaf_value)
        tag_ok = entry.get("tag") in VALID_TAGS
        source_ok = isinstance(entry.get("source"), str) and entry.get("source").strip() != ""
        ok(f"{tag}: provenance.sources[{leaf_path!r}] value/tag/source",
           value_ok and tag_ok and source_ok)

    # ---- 4. anchor cross-checks
    assumptions = set(provenance.get("assumptions", []))
    for leaf_path, entry in sources.items():
        anchor_id = entry.get("anchor_id")
        if not anchor_id:
            continue
        anchor = anchors.get(anchor_id)
        if anchor is None:
            ok(f"{tag}: anchor_id {anchor_id!r} (from {leaf_path!r}) exists", False)
            continue
        if anchor.get("status") == "missing":
            ok(f"{tag}: {leaf_path!r} anchored to a missing-evidence anchor carries E/A + assumptions",
               entry.get("tag") in ("E", "A") and leaf_path in assumptions)
        elif anchor.get("status") == "verified" and anchor.get("value") is not None:
            bound = anchor_bound(anchor)
            ok(f"{tag}: {leaf_path!r} value lies within anchor {anchor_id!r} tolerance",
               abs(float(entry["value"]) - float(anchor["value"])) <= bound + 1e-12)

    # ---- 5. required ranges: exact endpoints, and bracketing where the key is a design path
    ranges = provenance.get("ranges", {})
    for range_path, (lo, hi, unit) in REQUIRED_RANGES.items():
        entry = ranges.get(range_path)
        if entry is None:
            ok(f"{tag}: provenance.ranges[{range_path!r}] exists", False)
            continue
        ok(f"{tag}: provenance.ranges[{range_path!r}] endpoints/unit",
           entry.get("lo") == lo and entry.get("hi") == hi and entry.get("unit") == unit)
        is_path, value = resolve_design_path(design, range_path)
        if is_path:
            ok(f"{tag}: range {range_path!r} brackets the card's own value",
               entry["lo"] <= value <= entry["hi"])

    # ---- 6. single-dot selection arithmetic
    sds = provenance.get("single_dot_selection", {})
    density_cm2 = design.aperture.density_cm2
    diameter_um = design.aperture.diameter_um
    n_expected_recomputed = density_cm2 * math.pi * (diameter_um / 2.0) ** 2 * 1e-8
    n_expected_card = sds.get("n_expected")
    ok(f"{tag}: single_dot_selection.n_expected matches the formula to 1e-6 relative",
       n_expected_card is not None
       and math.isclose(n_expected_recomputed, n_expected_card, rel_tol=1e-6))
    ok(f"{tag}: single_dot_selection.n_expected < 0.5", n_expected_recomputed < 0.5)
    ok(f"{tag}: single_dot_selection.density_cm2 matches design.aperture.density_cm2",
       values_equal(sds.get("density_cm2"), density_cm2))

    # ---- 7. evaluate() at the card's own thermal.T_hs: must not raise
    try:
        result = evaluate(design, T_grid=[design.thermal.T_hs])
        evaluated = True
    except Exception:
        result = None
        evaluated = False
    ok(f"{tag}: evaluate() at thermal.T_hs does not raise", evaluated)
    if evaluated:
        ok(f"{tag}: evaluate() returns curves and scalars",
           "curves" in result and "scalars" in result)
        sc = result["scalars"]
        ok(f"{tag}: resolved scalars report emission.type='edge' active",
           sc.get("emission_type") == "edge")
        ok(f"{tag}: resolved scalars report ret.mode='confinement' active",
           sc.get("retention_source") == "confinement")
        ok(f"{tag}: resolved scalars report drive.mode='EL-transport' active",
           sc.get("drive_source") == "EL-transport")
        ok(f"{tag}: resolved scalars report dot.linewidth='anchored' active",
           sc.get("linewidth_source") == "anchored")
        ok(f"{tag}: resolved scalars report aperture.compose active",
           sc.get("aperture_compose") is True)
        ok(f"{tag}: resolved scalars report the card's own filter.track_material",
           sc.get("track_material") == design.filter.track_material)


def main() -> int:
    anchors_doc = yaml.safe_load(ANCHORS_PATH.read_text(encoding="utf-8"))
    anchors = {a["id"]: a for a in anchors_doc["anchors"]}
    for path in CARDS:
        check_card(path, anchors)
    print(f"{sum(CHECKS)}/{len(CHECKS)} rt-edge card checks passed")
    return 0 if all(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
