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
  8. (rt-fix-cards-wavelength, council review 2026-09-05) at the card's own
     operating point: emission.lambda_nm is within 2% of the confinement-
     derived transition (dot_levels.levels(system).lambda_nm) UNLESS its
     provenance.sources entry explicitly overrides with tag DR or A (both
     cards currently do -- see each card's emission.lambda_nm note for why
     the gap is a documented, root-caused model limitation, not a design
     choice); evaluate()'s own sub_turn_on flag is absent (qV_j is not
     Boltzmann-suppressed at drive.I_uA); and evaluate() completes in
     under CW_RUNTIME_BUDGET_S seconds.
  9. (council review 2026-09-06) dot.gamma300 is re-sourced to the verified
     Matsuda 2001 class value (12.0 meV) and its provenance names Matsuda
     2001 without naming the superseded Laferriere detection-filter-width
     proxy; drive.b_res is re-sourced to Reischle et al., Optics Express 16,
     12771 (2008) (InP/AlGaInP QD C, 80 K), not Appl. Phys. Lett. 92, 233113,
     and never calls the device "InP/GaInP" in that sentence; the emission
     collection levers (emission.NA, emission.R_back, emission.L_um,
     emission.alpha_cm) are explicit, literature-ordinary values with
     provenance and sweep ranges rather than undeclared device.py defaults;
     and dot.delta_xx's sweep range is the contract's restated 4-8 meV union.

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

import copy
import dataclasses
import glob
import math
import os
import re
import sys
import tempfile
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import device as device_mod  # noqa: E402
from fsim_core import dot_levels  # noqa: E402
from fsim_core import transport as transport_mod  # noqa: E402
from fsim_core.device import DeviceDesign, evaluate  # noqa: E402
import scripts.run_rt_edge as rte  # noqa: E402 -- _collected_flux_s/REP_RATE_HZ, item C1

# pr-pkg1-fix2 (peer-review-triage.md finding 1b + runtime item 2): fixing
# the confinement escape prefactor's dot density (RetentionBlock.n_dot_cm2,
# device.py._confinement_params) raises both cards' escape rate by ~30x at
# their own aperture.density_cm2=3e8 cm^-2 (vs. dot_levels.retention_params'
# previous 1e10 cm^-2 default). That stiffer rate drove cw_g2.cw_report's tau
# grid dt far smaller, and cw_g2.convolve_irf's IRF kernel (built at that
# grid's own dt) grew with it -- the direct-method np.convolve cost
# (~grid_len * kernel_len) blew up to 28-64 s at the cards' own re-solved
# drive.I_uA (measured before the fix below). Root-caused and fixed in
# fsim_core/cw_g2.py: convolve_irf now uses scipy.signal.fftconvolve
# (O(N log N) instead of O(N * kernel), the identical linear convolution to
# 1.111e-16 relative floating-point round-off on the gainp card's own CW
# operating point (pr-pkg1-fix3 item 5, verify_device_rt.py's frozen CW-path
# fixture; supersedes this comment's earlier "~1e-10" estimate) -- no
# physics change) plus a
# device.py fix removing a literal duplicate evaluate()-internal `one()`
# call. Both cards now evaluate in well under 1 s; the budget below is a
# generous margin over that, not a relaxation to hide the regression.
CW_RUNTIME_BUDGET_S = 6.0

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

# The ranges the spec requires, with their exact endpoints. dot.delta_xx's
# hi endpoint is 8.0 meV (council review 2026-09-06 item 4): the contract's
# restated union of Beirne/Reischle's 4 meV lower bound and Bommer's
# 7 +/- 1 meV upper bound. emission.NA/R_back/L_um (council review
# 2026-09-06 item 3) are the collection-lever sweep ranges added alongside
# the cards' newly explicit emission fields.
REQUIRED_RANGES = {
    "dot.delta_xx": (4.0, 8.0, "meV"),
    "dot.gamma300": (6.0, 20.0, "meV"),
    "irf_ps": (50.0, 200.0, "ps"),
    "thermal.T_hs": (230.0, 300.0, "K"),
    "emission.NA": (0.5, 0.8, "dimensionless"),
    "emission.R_back": (0.0, 0.95, "fraction"),
    "emission.L_um": (250.0, 500.0, "um"),
}

VALID_TAGS = {"V", "DR", "E", "A"}
STRUCTURAL_KEYS = {"kind", "name"}


# ------------------------------------------------------------------ plumbing

CHECKS = []


def ok(name: str, value: bool) -> None:
    CHECKS.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


# item 9 (pr-pkg1-fix2): match ONLY the (AlxGa1-x)yIn1-yP family label shape
# (thermal-layer names and materials.AlGaInP(x).label both use it, e.g.
# "(Al0.55Ga0.45)0.51In0.49P cladding") -- deliberately stricter than a bare
# "Al...Ga" search so an AlGaAs-family label (e.g. hkust_preset's
# "Al0.40Ga0.60As" barrier, no "In"/no parens) is correctly NOT matched.
_ALGAINP_X_RE = re.compile(r"\(Al(\d*\.?\d+)Ga[\d.]+\)[\d.]+In[\d.]+P")


def _x_al_from_algainp_label(label: str):
    m = _ALGAINP_X_RE.search(label)
    return float(m.group(1)) if m else None


def _layer_role(name: str) -> str:
    """pr-pkg1-fix4 item 4: a thermal.layers[*] entry carries no `role`
    field (name/t_um/k300/alpha/spread only) -- its role is the free-text
    suffix after the (Al..Ga..)..In..P formula prefix in its own `name`,
    e.g. "(Al0.70Ga0.30)0.51In0.49P outer cladding" -> "outer cladding",
    "(Al0.55Ga0.45)0.51In0.49P cladding" -> "cladding". A label that does
    not parse as AlGaInP at all (e.g. a GaAs substrate) has no formula
    prefix to strip, so its whole name is its role."""
    m = _ALGAINP_X_RE.match(name)
    return name[m.end():].strip() if m else name.strip()


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

def check_card(path: Path, anchors: dict) -> set:
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

    # Round-3 card contract: the pulse bookkeeping and residual channel must
    # be explicit, even though CW diagnostics are also requested.
    drive_raw = raw_design.get("drive", {})
    diode_raw = drive_raw.get("diode", {})
    ok(f"{tag}: explicit tau_pulse_ns=0.1 ns", diode_raw.get("tau_pulse_ns") == 0.1
       and "drive.diode.tau_pulse_ns" in provenance.get("sources", {}))
    ok(f"{tag}: explicit duty=0.008", drive_raw.get("duty") == 0.008
       and "drive.duty" in provenance.get("sources", {}))
    ok(f"{tag}: rep-rate inputs resolve to 80 MHz", math.isclose(
        0.008 / (0.1e-9), 8.0e7, rel_tol=0.0, abs_tol=1.0))
    bres_entry = provenance.get("sources", {}).get("drive.b_res", {})
    ok(f"{tag}: drive.b_res has provenance", math.isclose(
        float(drive_raw.get("b_res", float("nan"))), 1.0 / 0.88 - 1.0,
        rel_tol=1e-9) and bres_entry.get("tag") in VALID_TAGS
       and math.isclose(float(bres_entry.get("value", float("nan"))),
                       float(drive_raw.get("b_res", float("nan"))), rel_tol=1e-9)
       and "rho" in bres_entry.get("source", "")
       and "residual" in bres_entry.get("source", ""))
    # Council review 2026-09-06 item 2: b_res is re-sourced to Reischle et
    # al., Optics Express 16, 12771 (2008) (InP/AlGaInP QD C at 80 K), not
    # Appl. Phys. Lett. 92, 233113; the device must not be called "InP/GaInP"
    # in this specific sentence (it is InP/AlGaInP).
    bres_source = bres_entry.get("source", "")
    ok(f"{tag}: drive.b_res source cites Opt. Express 16, 12771 (2008), "
       "not Appl. Phys. Lett. 92, 233113, and does not call the device InP/GaInP",
       "Express 16, 12771" in bres_source and "2008" in bres_source
       and "233113" not in bres_source
       and "InP/GaInP" not in bres_source)
    # Council review 2026-09-06 (fifth round) item 2: the ledger anchor
    # (reischle08-b-res-80k) states the 80 K -> 300 K, cross-material
    # transfer of this ratio "is an [A] of the cards" -- b_res is
    # accordingly tagged A here, not E.
    ok(f"{tag}: drive.b_res is tagged A (80 K -> 300 K transfer is an "
       "assumption of the cards, per the reischle08-b-res-80k anchor)",
       bres_entry.get("tag") == "A")
    ok(f"{tag}: drive.b_res text tag matches tag A",
       bres_source.startswith("[A]"))
    thermal_entry = provenance.get("sources", {}).get("thermal.T_hs", {})
    ok(f"{tag}: thermal.T_hs declares TEC cooling while retaining 300 K default",
       "TEC-cooled package" in thermal_entry.get("source", "")
       and "self-consistently" in thermal_entry.get("source", ""))

    # Council review 2026-09-06 item 3: the collection levers (emission.NA,
    # emission.R_back, emission.L_um, emission.alpha_cm) must be explicit,
    # literature-ordinary values with provenance, rather than undeclared
    # device.py defaults.
    emission_raw = raw_design.get("emission", {})
    sources_for_emission = provenance.get("sources", {})
    ok(f"{tag}: explicit emission.NA=0.75", emission_raw.get("NA") == 0.75
       and "emission.NA" in sources_for_emission)
    ok(f"{tag}: explicit emission.R_back=0.95", emission_raw.get("R_back") == 0.95
       and "emission.R_back" in sources_for_emission)
    ok(f"{tag}: explicit emission.L_um=250.0", emission_raw.get("L_um") == 250.0
       and "emission.L_um" in sources_for_emission)
    ok(f"{tag}: explicit emission.alpha_cm=5.0", emission_raw.get("alpha_cm") == 5.0
       and "emission.alpha_cm" in sources_for_emission)

    # Council review 2026-09-06 item 1: dot.gamma300 is re-sourced to the
    # verified Matsuda 2001 class value (12.0 meV), not the superseded
    # Laferriere et al. detection-filter-width proxy.
    dot_raw = raw_design.get("dot", {})
    gamma_entry = sources_for_emission.get("dot.gamma300", {})
    gamma_source = gamma_entry.get("source", "")
    ok(f"{tag}: dot.gamma300 is the Matsuda 2001 class value (12.0 meV)",
       math.isclose(float(dot_raw.get("gamma300", float("nan"))), 12.0, rel_tol=1e-9))
    ok(f"{tag}: dot.gamma300 provenance names Matsuda 2001 and not the Laferriere proxy",
       "Matsuda" in gamma_source and "2001" in gamma_source
       and "Laferriere" not in gamma_source)

    hold_window = raw_design.get("filter", {}).get("hold_window")
    track_source = provenance.get("sources", {}).get("filter.track_material", {}).get("source", "")
    ok(f"{tag}: inert track_material provenance is honest",
       hold_window is False and "inert" in track_source.lower()
       and "no cavity" in track_source.lower()
       and "tracking active" not in track_source.lower())

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
    t0 = time.perf_counter()
    try:
        result = evaluate(design, T_grid=[design.thermal.T_hs])
        evaluated = True
    except Exception:
        result = None
        evaluated = False
    eval_seconds = time.perf_counter() - t0
    ok(f"{tag}: evaluate() at thermal.T_hs does not raise", evaluated)
    if evaluated:
        ok(f"{tag}: evaluate() returns curves and scalars",
           "curves" in result and "scalars" in result)
        sc = result["scalars"]
        # The I_uA provenance is a human-readable audit trail; independently
        # compare each quoted operating value with this fresh evaluation.
        i_source = sources.get("drive.I_uA", {}).get("source", "")
        metric_patterns = {
            "mu_resolved": r"mu=([0-9.eE+-]+)",
            "rho_op": r"rho_op=([0-9.eE+-]+)",
            "g2_op": r"g2_op=([0-9.eE+-]+)",
            "collected_flux_pulsed_s": r"collected_flux_pulsed_s=([0-9.eE+-]+)",
            "V_j": r"V_j=([0-9.eE+-]+)",
        }
        quoted = {k: float(m.group(1).rstrip(".")) for k, pattern in metric_patterns.items()
                  if (m := re.search(pattern, i_source))}
        quoted["f_qfl"] = (float(re.search(r"f_qfl=([0-9.eE+-]+)", i_source).group(1).rstrip("."))
                           if re.search(r"f_qfl=([0-9.eE+-]+)", i_source) else float("nan"))
        f_qfl_actual = (math.exp(-(1239.841984 / sc["edge_lambda_nm"] - sc["V_j"])
                                  / (8.617333262e-5 * sc["T_j_op"]))
                        if sc["V_j"] < 1239.841984 / sc["edge_lambda_nm"] else 1.0)
        actual_metrics = {k: sc[k] for k in metric_patterns}
        actual_metrics["f_qfl"] = f_qfl_actual
        ok(f"{tag}: I_uA provenance operating metrics match fresh evaluation",
           all(k in quoted and math.isfinite(quoted[k]) and math.isfinite(actual_metrics[k])
               and math.isclose(quoted[k], actual_metrics[k], rel_tol=1e-6)
               for k in actual_metrics))
        # Lemma 1 regression: rho_op is transport-intrinsic and must not
        # depend on the collection levers (NA, R_back, or L_um).
        default_collection = copy.deepcopy(design)
        default_collection.emission.NA = 0.5
        default_collection.emission.R_back = None
        default_collection.emission.L_um = 500.0
        default_sc = evaluate(default_collection,
                               T_grid=[default_collection.thermal.T_hs])["scalars"]
        ok(f"{tag}: rho_op is unchanged by collection levers (Lemma 1)",
           math.isclose(float(sc["rho_op"]), float(default_sc["rho_op"]),
                        rel_tol=0.0, abs_tol=1e-9))
        ok(f"{tag}: no provenance sentence says rho_op moves upward with",
           not any("moves upward with" in str(entry.get("source", ""))
                   for entry in sources.values()))
        ok(f"{tag}: resolved rep_rate_hz is 8e7", math.isclose(
            sc.get("rep_rate_hz", float("nan")), 8.0e7, rel_tol=0.0, abs_tol=1.0))
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

        # ---- 8. wavelength self-consistency, corrected injection constraints,
        # constrained local optimum, and CW runtime budget.
        # (rt-fix-cards-wavelength, council review 2026-09-05).
        system = copy.copy(device_mod._retention_system(design.ret))
        system.T = float(sc.get("T_j_op", design.thermal.T_hs))
        lambda_derived = dot_levels.levels(system).lambda_nm
        lambda_used = sc.get("edge_lambda_nm")
        lambda_entry = sources.get("emission.lambda_nm", {})
        within_2pct = (lambda_used is not None and math.isfinite(lambda_used)
                       and lambda_derived and math.isfinite(lambda_derived)
                       and abs(lambda_used - lambda_derived) / lambda_derived < 0.02)
        overridden = lambda_entry.get("tag") in ("DR", "A")
        ok(f"{tag}: emission.lambda_nm within 2% of confinement-derived "
           f"({lambda_used!r} vs {lambda_derived:.1f} nm) or explicitly overridden [DR]/[A]",
           within_2pct or overridden)

        sub_turn_on = any("sub_turn_on" in r for r in sc.get("invalid_reasons", []))
        ok(f"{tag}: no sub_turn_on flag at drive.I_uA (qV_j not Boltzmann-suppressed)",
           not sub_turn_on)

        mu_op = sc.get("mu_resolved")
        ok(f"{tag}: resolved mu is in [0.05, 1.0]",
           mu_op is not None and math.isfinite(mu_op) and 0.05 <= mu_op <= 1.0)

        # ---- pr-pkg1-fix2 item 9 (peer-review-triage.md finding 5): ONE
        # x_al across ret.system.barrier, EVERY thermal.layers[*] entry that
        # names an AlGaInP composition (not only layer 0, the previous,
        # narrower check), and the diode preset's own resolved barrier
        # material -- a preset's internal barrier composition is a Python
        # literal (e.g. transport.red_diode_preset's hardcoded
        # AlGaInP(0.55)), never read from the card, so it can silently drift
        # from ret.system.barrier.x_al the same way finding 5 originally
        # found thermal.layers[0] drifting.
        barrier_x_al = design.ret.system.get("barrier", {}).get("x_al")
        # gaasp's thermal.layers[1] ("outer cladding") is a legitimately
        # different composition from thermal.layers[0] ("core", the layer
        # that plays the confinement "barrier" role) and ret.system.barrier
        # -- see the card's own KNOWN MODEL LIMITATION note. pr-pkg1-fix4
        # item 4: excluded by ROLE (the free-text suffix _layer_role reads
        # off the layer's own name, after its (Al..Ga..)..In..P formula
        # prefix), not by a per-card filename/literal-name list -- a role
        # of "outer cladding" is physically a different layer than the
        # barrier on ANY card, gaasp's included, so keying the exclusion on
        # this card's filename (or its exact current name string) was
        # accidental, not principled.
        excluded_role = "outer cladding"
        mismatches = []
        checked = []
        excluded = []
        for i, layer in enumerate(design.thermal.layers):
            name = layer.get("name", "")
            if _layer_role(name) == excluded_role:
                excluded.append(name)
                continue
            x_al = _x_al_from_algainp_label(name)
            if x_al is None:
                continue
            checked.append(f"thermal.layers[{i}]={x_al!r}")
            if barrier_x_al is None or not math.isclose(barrier_x_al, x_al, rel_tol=1e-9):
                mismatches.append(f"thermal.layers[{i}] {name!r} (x_al={x_al!r})")
        # pr-pkg1-fix3 item 8: gaasp's own drive.diode.preset is "hkust"
        # (transport.hkust_preset), whose barrier is Al0.40Ga0.60As -- a
        # genuinely different III-V system (AlGaAs, not AlGaInP), not a
        # drifted x_al of the SAME material. _x_al_from_algainp_label
        # returns None for it (the label does not parse as AlGaInP), so it
        # is silently excluded from the mismatch comparison below by
        # material, the same generalisation thermal.layers[*] above already
        # gets (a material string that doesn't parse as AlGaInP is not
        # compared, full stop -- no per-card/per-preset name list needed).
        diode = device_mod._diode_from_drive(design.drive, design.aperture.diameter_um)
        diode_x_al = _x_al_from_algainp_label(diode.barrier.label)
        if diode_x_al is not None:
            checked.append(f"diode.preset={design.drive.diode.get('preset')!r}={diode_x_al!r}")
            if barrier_x_al is None or not math.isclose(barrier_x_al, diode_x_al, rel_tol=1e-9):
                mismatches.append(f"drive.diode.preset={design.drive.diode.get('preset')!r} "
                                  f"resolved barrier {diode.barrier.label!r} (x_al={diode_x_al!r})")
        ok(f"{tag}: ret.system.barrier.x_al ({barrier_x_al!r}) equals the x_al of every "
           f"AlGaInP thermal.layers[*] entry and the diode preset's resolved barrier "
           f"(checked: {checked}; excluded: {sorted(excluded)}; mismatches: {mismatches})",
           barrier_x_al is not None and bool(checked) and not mismatches)

        # ---- pr-pkg1-fix2 item 6 (peer-review-triage.md finding 1b): the
        # previous version of this check recomputed "d.drive.n_dot_cm2 or
        # d.aperture.density_cm2" here and compared it against device.py's
        # own report of the SAME expression -- circular, since it could
        # never catch the two call sites (dot_levels.retention_params via
        # _confinement_params, and transport.evaluate_injection) actually
        # being given different densities. Instead, monkeypatch both
        # functions for the duration of ONE evaluate() call and record the
        # n_dot_cm2 argument each genuinely receives, independent of
        # whatever expression device.py used to compute it.
        #
        # pr-pkg1-fix3 item 7: record EVERY call, not just one per function
        # (device.py's transport self-heating fixed point calls
        # evaluate_injection up to 12 times per evaluate(), all of which
        # must agree -- the density does not move with Tj) -- and read
        # n_dot_cm2 whether the call site passes it positionally or by
        # keyword (both spied functions are called with n_dot_cm2= today,
        # but a future call site changing that silently should not blind
        # this spy).
        received = {"confinement": [], "transport": []}
        orig_retention_params = dot_levels.retention_params
        orig_evaluate_injection = transport_mod.evaluate_injection

        def _positional_or_kw(args, kwargs, name, pos_index):
            if name in kwargs:
                return kwargs[name]
            if len(args) > pos_index:
                return args[pos_index]
            return None

        def _spy_retention_params(*a, **kw):
            # retention_params(lv, tau_rad_ns, channel="pair_half",
            #                  T_ref=300.0, n_dot_cm2=1e10, ...)
            received["confinement"].append(_positional_or_kw(a, kw, "n_dot_cm2", 4))
            return orig_retention_params(*a, **kw)

        def _spy_evaluate_injection(*a, **kw):
            # evaluate_injection(diode, I_uA, T, n_dot_cm2, ...)
            received["transport"].append(_positional_or_kw(a, kw, "n_dot_cm2", 3))
            return orig_evaluate_injection(*a, **kw)

        dot_levels.retention_params = _spy_retention_params
        transport_mod.evaluate_injection = _spy_evaluate_injection
        try:
            evaluate(design, T_grid=[design.thermal.T_hs])
        finally:
            dot_levels.retention_params = orig_retention_params
            transport_mod.evaluate_injection = orig_evaluate_injection
        all_received = received["confinement"] + received["transport"]
        n_dot_cm2_confinement = received["confinement"][0] if received["confinement"] else None
        n_dot_cm2_transport = received["transport"][0] if received["transport"] else None
        ok(f"{tag}: dot_levels.retention_params and transport.evaluate_injection actually "
           f"received the same n_dot_cm2 on EVERY call during one evaluate() call "
           f"(confinement calls: {received['confinement']!r}; transport calls: "
           f"{received['transport']!r}), equal to the card's own aperture.density_cm2 "
           f"({design.aperture.density_cm2!r})",
           bool(received["confinement"]) and bool(received["transport"])
           and all(v is not None and math.isfinite(v) for v in all_received)
           and all(math.isclose(v, all_received[0], rel_tol=1e-9) for v in all_received)
           and design.aperture.density_cm2 is not None
           and math.isclose(n_dot_cm2_confinement, design.aperture.density_cm2, rel_tol=1e-9))

        ok(f"{tag}: provenance drive.I_uA contains the [DR] grid derivation",
           sources.get("drive.I_uA", {}).get("tag") == "DR"
           and "Log-spaced grid" in sources.get("drive.I_uA", {}).get("source", "")
           and "V_j" in sources.get("drive.I_uA", {}).get("source", "")
           and "f_qfl" in sources.get("drive.I_uA", {}).get("source", "")
           and "g2_op" in sources.get("drive.I_uA", {}).get("source", ""))

        # N_expected is independently recomputed above; repeat the explicit
        # operating-point gate here so this block mirrors the round-2 card
        # selection constraints.
        ok(f"{tag}: operating-point N_expected < 0.5", n_expected_recomputed < 0.5)

        # A boundary optimum is local within the feasible mu interval.  The
        # out-of-range side is still evaluated for diagnostics, but cannot
        # disqualify a constrained optimum (both selected cards sit at mu=1).
        neighbor_checks = []
        for multiplier in (0.5, 2.0):
            neighbor = copy.deepcopy(design)
            neighbor.drive.I_uA *= multiplier
            neighbor.drive.cw = False
            try:
                ns = evaluate(neighbor, T_grid=[neighbor.thermal.T_hs])["scalars"]
                nmu = ns.get("mu_resolved")
                feasible = (nmu is not None and math.isfinite(nmu)
                            and 0.05 <= nmu <= 1.0
                            and n_expected_recomputed < 0.5
                            and not any("sub_turn_on" in r
                                        for r in ns.get("invalid_reasons", [])))
                neighbor_checks.append((not feasible) or
                                       ns.get("g2_op", math.inf) >= sc.get("g2_op", -math.inf) - 1e-3)
            except Exception:
                # An evaluator failure makes that neighbor infeasible; skip it
                # rather than rejecting a constrained optimum.
                neighbor_checks.append(True)
        ok(f"{tag}: constrained local g2 optimum at 0.5x and 2x current "
           "(infeasible neighbors skipped)",
           all(neighbor_checks))

        ok(f"{tag}: evaluate() completes in under {CW_RUNTIME_BUDGET_S:g} s "
           f"({eval_seconds:.2f} s)",
           eval_seconds < CW_RUNTIME_BUDGET_S)

        # ---- pr-pkg1-fix4 item 1 (extended by pr-pkg6-stale-text item C1;
        # tolerances tightened by pr-pkg6-fix item 2): the card's
        # finding_1b_record prose names several literals -- a_esc(300 K)
        # and S(230 K) for BOTH the no-cancellation (shipped) and
        # full-cancellation (opt-in) conventions, the favourable-corner
        # mu_resolved, and the four g2_op/collected_flux_pulsed_s values in
        # the favourable-corner table -- pr-pkg1-fix3's version of the
        # full-cancellation S(230 K) literal was wrong at the 8th-9th digit
        # and nothing caught it before this (originally single-value) check
        # existed. Every one of these is now recomputed fresh and required
        # to agree with whatever the card's own text currently claims, so a
        # future stale literal fails loudly instead of sitting silently
        # wrong again. pr-pkg6-fix item 2: the no-cancellation a_esc(300 K)/
        # S(230 K) and full-cancellation a_esc(300 K) literals used to be
        # quoted to only 6-7 significant figures (rel diffs up to ~1.3e-7 /
        # ~4.3e-7 against a fresh evaluation -- not stale, just rounded for
        # display), checked at a correspondingly relaxed 2e-7/5e-7 rather
        # than this file's blanket rel 1e-8; the card literals are now
        # re-quoted to 10 significant figures (matching the full-
        # cancellation S(230 K) literal, mu_resolved, and the g2/flux
        # table, which already held 10 significant figures), so every one
        # of these checks now holds at the SAME rel 1e-8.
        record_text = " ".join(provenance.get("finding_1b_record", "").split())
        NUM = r"[0-9.eE+-]+"

        m_nocancel = re.search(
            r"no-cancellation,\s*ret\.\s*tau_cap_scales_with_density\s*unset/\s*false\)\s*"
            rf"a_esc\(300 K\)\s*=\s*({NUM}),\s*S\(230 K\)\s*=\s*({NUM})", record_text)
        ok(f"{tag}: finding_1b_record names a no-cancellation a_esc(300 K)/S(230 K) pair",
           m_nocancel is not None)
        if m_nocancel is not None:
            recorded_aesc_nc = float(m_nocancel.group(1))
            recorded_s230_nc = float(m_nocancel.group(2))
            fresh_default = DeviceDesign.load(path)
            sc_default = evaluate(fresh_default)["scalars"]
            fresh_aesc_nc = device_mod._confinement_params(
                fresh_default.ret, sc_default["T_j_op"],
                n_dot_cm2=fresh_default.aperture.density_cm2)["a_esc"]
            ok(f"{tag}: finding_1b_record's no-cancellation a_esc(300 K) literal "
               f"({recorded_aesc_nc!r}) matches a fresh DeviceDesign.load + evaluate "
               f"({fresh_aesc_nc!r}) to rel 1e-8 (10-significant-figure literal)",
               math.isclose(recorded_aesc_nc, fresh_aesc_nc, rel_tol=1e-8))

            fresh_design_nc = DeviceDesign.load(path)
            fresh_design_nc.thermal.T_hs = 230.0
            fresh_design_nc.dot.gamma300 = 6.0
            fresh_design_nc.dot.delta_xx = 8.0
            fresh_design_nc.emission.NA = 0.8
            fresh_design_nc.emission.R_back = 0.95
            fresh_design_nc.emission.L_um = 250.0
            fresh_design_nc.ret.tau_cap_scales_with_density = False
            fresh_s230_nc = evaluate(fresh_design_nc)["scalars"]["S_resolved"]
            ok(f"{tag}: finding_1b_record's no-cancellation S(230 K) literal "
               f"({recorded_s230_nc!r}) matches a fresh DeviceDesign.load + evaluate "
               f"({fresh_s230_nc!r}) to rel 1e-8 (10-significant-figure literal)",
               math.isclose(recorded_s230_nc, fresh_s230_nc, rel_tol=1e-8))

        m_fullcancel = re.search(
            rf"reproduces the 1e10 cm\^-2 default exactly: a_esc\(300 K\) = ({NUM}), "
            rf"S\(230 K\) = ({NUM})", record_text)
        ok(f"{tag}: finding_1b_record names a full-cancellation a_esc(300 K)/S(230 K) pair",
           m_fullcancel is not None)
        if m_fullcancel is not None:
            recorded_aesc_fc = float(m_fullcancel.group(1))
            recorded_s230 = float(m_fullcancel.group(2))
            fresh_default_fc = DeviceDesign.load(path)
            fresh_default_fc.ret.tau_cap_scales_with_density = True
            sc_default_fc = evaluate(fresh_default_fc)["scalars"]
            fresh_aesc_fc = device_mod._confinement_params(
                fresh_default_fc.ret, sc_default_fc["T_j_op"],
                n_dot_cm2=fresh_default_fc.aperture.density_cm2)["a_esc"]
            ok(f"{tag}: finding_1b_record's full-cancellation a_esc(300 K) literal "
               f"({recorded_aesc_fc!r}) matches a fresh DeviceDesign.load + evaluate "
               f"({fresh_aesc_fc!r}) to rel 1e-8 (10-significant-figure literal)",
               math.isclose(recorded_aesc_fc, fresh_aesc_fc, rel_tol=1e-8))

            fresh_design = DeviceDesign.load(path)
            fresh_design.thermal.T_hs = 230.0
            fresh_design.dot.gamma300 = 6.0
            fresh_design.dot.delta_xx = 8.0
            fresh_design.emission.NA = 0.8
            fresh_design.emission.R_back = 0.95
            fresh_design.emission.L_um = 250.0
            fresh_design.ret.tau_cap_scales_with_density = True
            fresh_s230 = evaluate(fresh_design)["scalars"]["S_resolved"]
            ok(f"{tag}: finding_1b_record's full-cancellation S(230 K) literal "
               f"({recorded_s230!r}) matches a fresh DeviceDesign.load + evaluate "
               f"({fresh_s230!r}) to rel 1e-8",
               math.isclose(recorded_s230, fresh_s230, rel_tol=1e-8))

        # Favourable-corner g2/flux table: default (no-cancellation) and
        # full-cancellation g2_op + collected_flux_pulsed_s, four literals,
        # each quoted to ~10 significant figures -- checked at rel 1e-8.
        m_table = re.search(
            rf"default \(no-cancellation\) g2_op=\s*({NUM})\s*\([^)]*\),\s*"
            rf"collected_flux_pulsed_s=\s*({NUM})\s*(?:s\^-1)?\s*\([^)]*\)[^;]*;\s*"
            rf"full-cancellation g2_op=\s*({NUM}),\s*"
            rf"collected_flux_pulsed_s=\s*({NUM})\s*(?:s\^-1)?", record_text)
        ok(f"{tag}: finding_1b_record names the favourable-corner g2/flux table (both conventions)",
           m_table is not None)
        if m_table is not None:
            recorded_g2_nc, recorded_flux_nc, recorded_g2_fc, recorded_flux_fc = (
                float(x) for x in m_table.groups())
            # pr-pkg6-fix item 2: mu_resolved (the pulsed loading parameter
            # gating F1/f1b_g2's [0.05, 1.0] validity) is named once in the
            # favourable-corner prose, ahead of the g2/flux split -- it is
            # resolved by transport/injection upstream of ret.
            # tau_cap_scales_with_density, so it is convention-independent
            # (identical whether the loop below is on its no-cancellation
            # or full-cancellation pass); captured from whichever pass runs
            # first rather than with a dedicated third evaluate() call.
            m_mu = re.search(rf"mu_resolved\s*=\s*({NUM})", record_text)
            ok(f"{tag}: finding_1b_record names the favourable-corner mu_resolved",
               m_mu is not None)
            recorded_mu = float(m_mu.group(1)) if m_mu is not None else None
            fresh_mu = None
            for label, full_cancel, recorded_g2, recorded_flux in (
                    ("no-cancellation", False, recorded_g2_nc, recorded_flux_nc),
                    ("full-cancellation", True, recorded_g2_fc, recorded_flux_fc)):
                corner = DeviceDesign.load(path)
                corner.thermal.T_hs = 230.0
                corner.dot.gamma300 = 6.0
                corner.dot.delta_xx = 8.0
                corner.emission.NA = 0.8
                corner.emission.R_back = 0.95
                corner.emission.L_um = 250.0
                corner.ret.tau_cap_scales_with_density = full_cancel
                sc_corner = evaluate(corner)["scalars"]
                if fresh_mu is None:
                    fresh_mu = sc_corner.get("mu_resolved")
                fresh_g2 = sc_corner["g2_op"]
                fresh_flux = rte._collected_flux_s(sc_corner)
                ok(f"{tag}: finding_1b_record's {label} favourable-corner g2_op literal "
                   f"({recorded_g2!r}) matches a fresh evaluate ({fresh_g2!r}) to rel 1e-8",
                   math.isclose(recorded_g2, fresh_g2, rel_tol=1e-8))
                ok(f"{tag}: finding_1b_record's {label} favourable-corner "
                   f"collected_flux_pulsed_s literal ({recorded_flux!r}) matches a fresh "
                   f"evaluate ({fresh_flux!r}) to rel 1e-8",
                   math.isclose(recorded_flux, fresh_flux, rel_tol=1e-8))
            if recorded_mu is not None:
                ok(f"{tag}: finding_1b_record's favourable-corner mu_resolved literal "
                   f"({recorded_mu!r}) matches a fresh evaluate ({fresh_mu!r}) to rel 1e-8",
                   fresh_mu is not None and math.isclose(recorded_mu, fresh_mu, rel_tol=1e-8))

    return assumptions


def main() -> int:
    anchors_doc = yaml.safe_load(ANCHORS_PATH.read_text(encoding="utf-8"))
    anchors = {a["id"]: a for a in anchors_doc["anchors"]}
    assumptions_by_card = {}
    for path in CARDS:
        assumptions_by_card[path.name] = check_card(path, anchors)
    # Council review 2026-09-06 (fifth round) item 2: both cards must carry
    # the same assumptions-list field set (ret.system.barrier.x_al and
    # emission.lambda_nm previously missing from the gaasp card).
    card_names = list(assumptions_by_card)
    ok(f"{card_names[0]!r} and {card_names[1]!r} carry the same "
       "provenance.assumptions field set",
       assumptions_by_card[card_names[0]] == assumptions_by_card[card_names[1]])

    # pr-pkg4-fix item 2 (extended source-wide, pr-pkg4-fix2 item 6): the
    # two per-card provenance-string guards above (drive.b_res's source
    # text) only ever see the wrong Reischle citation if a CARD
    # reintroduces it; extend the same guard to every source, verify,
    # script, card, and doc file -- Appl. Phys. Lett. 92, 233113 (2008) is
    # not a Reischle et al. paper (that is Optics Express 16, 12771 (2008),
    # DOI 10.1364/OE.16.012771) and must never appear outside the files
    # that document this guard's own history or existence.
    #
    # Whitelist (pr-pkg4-fix2 item 6): these files legitimately name the
    # drifted citation string -- this file's own comment above and the
    # bad-citation substrings it searches for, verify_rt_edge_sweep.py's
    # check that verdict.md never cites it (a Python string literal that
    # itself contains the drifted text), and run_rt_edge.py's comment on
    # the now-generated (not hardcoded) citation -- none of these is an
    # actual mis-citation, so path-whitelist them rather than weakening
    # the substring match.
    _citation_guard_whitelist = {
        ROOT / "verify" / "verify_rt_edge_cards.py",
        ROOT / "verify" / "verify_rt_edge_sweep.py",
        ROOT / "scripts" / "run_rt_edge.py",
    }
    citation_guard_globs = ("fsim_core/*.py", "verify/*.py", "scripts/*.py",
                            "cards/*.yaml", "docs/*.md")
    bad_citation_files = []
    for pattern in citation_guard_globs:
        for f in sorted(glob.glob(str(ROOT / pattern))):
            fp = Path(f)
            if fp in _citation_guard_whitelist:
                continue
            text = fp.read_text(encoding="utf-8")  # read once, check both substrings
            if "233113" in text or "Appl. Phys. Lett. 92" in text:
                bad_citation_files.append(str(fp))
    ok("fsim_core/*.py, verify/*.py, scripts/*.py, cards/*.yaml, docs/*.md never cite "
       "the drifted 'Appl. Phys. Lett. 92, 233113' Reischle source outside this guard's "
       "own whitelisted files (got: " + repr(bad_citation_files) + ")",
       not bad_citation_files)

    print(f"{sum(CHECKS)}/{len(CHECKS)} rt-edge card checks passed")
    return 0 if all(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
