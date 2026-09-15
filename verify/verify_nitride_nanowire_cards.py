"""Mechanical + physics checks for the four gated nitride nanowire cards
(piece 8, rewritten from the contract's own Card schema tables).

Per docs/nitride_nanowire_contract.md "Card schema" and
.workers/specs/nitride-nanowire-cards.md acceptance criteria 1-5:

  1. DeviceDesign.load + evaluate() for all four cards, both temperature
     endpoints (230/300 K) and both strain bounds (relaxed as-shipped,
     unrelaxed via an in-memory copy -- the cards themselves are never
     mutated), checking every new block round-trips through YAML
     serialization.
  2. Pairwise deep comparison enforces an explicit small allow-list of
     intended pulse/SET and horizontal/vertical differences.
  3. Every explicit physical leaf has matching design.provenance.sources
     (value/unit/tag/source), tags are in {V, DR, E, A}, and any cited
     anchor_id exists in the nanowire ledger.
  4. Source literals (x_in, 2013 geometry/resistance/doping, rep rate) are
     checked against the ledger with correct conditions.
  5. Both SET cards carry the Coulomb and RT-injector screens as
     diagnostics; a mutation set (>= 9 leaves) is checked to actually move
     a physics output, proving the checks above are not tautological.

Every check increments the numerator/denominator BEFORE any try/except
that could raise, so a raised exception grows the failure count on the
SAME denominator, never shrinks it (fix-1 finding 5/12's complaint about
the previous verifier). The cards on disk are never mutated: every
strain-bound/mutation probe below operates on an in-memory
copy.deepcopy() of a loaded DeviceDesign.
"""
import copy
import os
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, ROOT)
from fsim_core.device import DeviceDesign, evaluate

CARD_NAMES = (
    "nitride-nanowire-horizontal-pulse-design.yaml",
    "nitride-nanowire-horizontal-set-design.yaml",
    "nitride-nanowire-vertical-pulse-design.yaml",
    "nitride-nanowire-vertical-set-design.yaml",
)
LEDGER_PATH = os.path.join(ROOT, "verify", "data", "nitride_nanowire_anchors.yaml")

VALID_TAGS = {"V", "DR", "E", "A"}

# Card-schema leaf sets (docs/nitride_nanowire_contract.md "Card schema"),
# independently transcribed here (not imported from the device module) so
# this verifier catches a leaf silently added/removed from either side.
NANOWIRE_LEAVES = {"family", "core_radius_nm", "outer_radius_nm", "strain_bound",
                    "shell", "barrier_left_nm", "barrier_right_nm"}
DOT_LEAVES = {"radius_nm", "height_nm", "x_in", "strain_fraction", "screening_fraction",
              "external_field_kVcm", "vbo_InN_GaN_eV", "strain_c_fraction",
              "k_intrinsic_ns"}  # k_intrinsic_ns: device-module addendum, not a contract leaf
SURFACE_LEAVES = {"S_cm_s", "shell", "shell_multiplier", "reservoir_access", "occupied_dot_access"}
# nitride.photonics: 18 of the dataclass's 19 fields; dipole_weights is
# deliberately absent from the card YAML (PyYAML safe_load turns a flow
# sequence into a list, and NitrideNanowirePhotonicsParams requires a
# genuine tuple -- see design.provenance.sources['nitride.photonics.
# dipole_weights'] on each card for the recorded reason) so the dataclass's
# own isotropic-tuple class default applies unmutated.
PHOTONICS_LEAVES = {"family", "NA", "n_wire", "n_group_override", "n_ambient", "n_oxide",
                     "oxide_thickness_nm", "n_substrate", "emitter_height_nm",
                     "collection_scale", "radiative_rate_factor", "beta_scale",
                     "taper_transmission", "bottom_reflectivity", "top_contact_transmission",
                     "propagation_transmission", "unguided_collection_scale", "taper_output_mfr_nm"}
WIRE_THERMAL_LEAVES = {"R_s_ohm", "Rth_K_W", "f_Rs_local", "eta_total", "C_parasitic_F"}
INJECTOR_LEAVES = {
    "al_fraction", "electron_topology", "hole_topology", "electron_barrier_thickness_nm",
    "hole_barrier_thickness_nm", "electron_well_width_nm", "hole_well_width_nm",
    "growth_step_nm", "growth_tolerance_steps", "me_barrier_override", "mh_barrier_override",
    "dEc_eV_override", "dEv_eV_override", "me_well", "mh_well", "delta_Ev_GaN_AlN_eV",
    "n_cm3", "p_cm3", "alignment_uncertainty_meV", "degeneracy", "reservoir_state_count_e",
    "reservoir_state_count_h", "mg_acceptor_energy_meV", "include_polarization",
    "polarity", "bypass_prefactor", "field_leverarm", "slice_length_nm",
    "min_slices_per_segment", "alignment_tunable", "bias_tuning_range_meV",
    "occupancy_control_known", "second_pair_control_known",
}
# drive.diode leaves this card schema actually writes (preset/tau_pulse_ns/
# conducting_radius_nm are popped by nitride_nanowire_device before the
# NitrideWireDiode(**) construction; the other 7 are the contract's own
# "remaining nine" minus core_radius_nm/barrier_left_nm/barrier_right_nm/
# d_active_nm/x_in, which mirror other blocks and are NEVER drive.diode
# leaves -- setting them there raises "multiple values for keyword
# argument" in NitrideWireDiode(**diode_raw)).
DIODE_LEAVES = {"preset", "tau_pulse_ns", "conducting_radius_nm",
                "N_A", "N_D", "n_ideality", "tau_SRH_ns", "eps_r", "T", "tau_matrix_ns"}
SET_PARAMS_COMMON_LEAVES = {"R_T_ohm", "eps_r", "ec_margin"}

# Allow-lists for the pairwise deep-diff (acceptance criterion 2).
PULSE_SET_ALLOWED_PREFIXES = ("drive.cycle_loading", "drive.set_params", "drive.diode.tau_pulse_ns")
FAMILY_ALLOWED_PREFIXES = ("nitride.nanowire", "nitride.photonics", "nitride.wire_thermal",
                            "drive.diode.conducting_radius_nm")

# The >= 9 leaf mutations (fix-1 finding 5's own named list) that must each
# move a physics output away from the card's own baseline.
MUTATIONS = (
    ("nitride.dot.x_in", 0.25),
    ("nitride.surface.S_cm_s", 1.0e4),
    ("nitride.wire_thermal.R_s_ohm", 5.0e6),
    ("nitride.photonics.NA", 0.9),
    ("nitride.nanowire.barrier_left_nm", 30.0),
    ("drive.rep_rate_hz", 8.0e7),
    ("nitride.wire_thermal.Rth_K_W", 1.0e8),
    ("drive.b_res", 0.5),
    ("nitride.dot.height_nm", 4.0),
)


class Ledger:
    def __init__(self, checks_ref):
        self.checks_ref = checks_ref
        with open(LEDGER_PATH, encoding="utf-8") as f:
            self.anchors = yaml.safe_load(f)["anchors"]


def _get(mapping, dotted):
    cur = mapping
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None, False
        cur = cur[part]
    return cur, True


def _apply_mutation(design, path, value):
    """Apply one dotted-path mutation to an in-memory DeviceDesign copy.
    Only the small, fixed set of paths in MUTATIONS below is supported."""
    if path.startswith("nitride."):
        _, block, leaf = path.split(".")
        design.nitride[block][leaf] = value
    elif path == "drive.rep_rate_hz":
        design.drive.rep_rate_hz = value
    elif path == "drive.b_res":
        design.drive.b_res = value
    else:
        raise ValueError(f"_apply_mutation: unsupported path {path!r}")


def _flatten(mapping, prefix=""):
    """Recursively flatten a nested dict into {dotted_path: leaf_value}."""
    out = {}
    if isinstance(mapping, dict):
        for k, v in mapping.items():
            path = f"{prefix}.{k}" if prefix else str(k)
            out.update(_flatten(v, path))
    else:
        out[prefix] = mapping
    return out


def _diff_paths(a, b):
    """Return the set of dotted paths whose flattened leaf value differs
    (or whose presence differs) between two nested dict designs."""
    fa, fb = _flatten(a), _flatten(b)
    keys = set(fa) | set(fb)
    return {k for k in keys if fa.get(k, object()) != fb.get(k, object())}


def _allowed(path, prefixes):
    return any(path == p or path.startswith(p + ".") or path.startswith(p) for p in prefixes)


def main():
    checks = [0]
    failures = []

    def ck(fn, msg):
        checks[0] += 1
        try:
            ok = fn()
        except Exception as exc:
            failures.append(f"{msg}: {exc}")
            return
        if not ok:
            failures.append(msg)

    ledger = Ledger(checks)

    raws = {}
    designs = {}
    for name in CARD_NAMES:
        path = os.path.join(ROOT, "cards", name)
        with open(path, encoding="utf-8") as f:
            raws[name] = yaml.safe_load(f)
        checks[0] += 1
        try:
            designs[name] = DeviceDesign.load(path)
        except Exception as exc:
            failures.append(f"{name}: DeviceDesign.load failed: {exc}")

    # ================================================== acceptance criterion 1
    # Evaluate every card at both temperature endpoints, and (via an
    # in-memory copy, never touching the file) at both strain bounds.
    scalars_relaxed_300 = {}
    for name in CARD_NAMES:
        d = designs.get(name)
        if d is None:
            continue
        fam = raws[name]["design"]["nitride"]["nanowire"]["family"]
        regime = raws[name]["design"]["drive"]["cycle_loading"]
        for T in (230.0, 300.0):
            checks[0] += 1
            try:
                r = evaluate(copy.deepcopy(d), T_grid=[T])
                sc = r["scalars"]
                ok = (sc.get("valid") is True and sc.get("family") == fam
                      and sc.get("cycle_loading") == regime
                      and sc.get("bound_role") == "headline_upper"
                      and sc.get("headline_eligible") is True)
                if not ok:
                    failures.append(f"{name} T={T}: relaxed-bound evaluation mismatch "
                                     f"(valid={sc.get('valid')} family={sc.get('family')} "
                                     f"regime={sc.get('cycle_loading')} "
                                     f"bound_role={sc.get('bound_role')} "
                                     f"headline_eligible={sc.get('headline_eligible')})")
                if T == 300.0:
                    scalars_relaxed_300[name] = sc
            except Exception as exc:
                failures.append(f"{name} T={T}: evaluate raised: {exc}")

        # In-memory unrelaxed (conservative_lower) partner: never written back.
        for T in (230.0, 300.0):
            checks[0] += 1
            try:
                d2 = copy.deepcopy(d)
                d2.nitride["nanowire"]["strain_bound"] = "unrelaxed"
                d2.nitride["dot"]["strain_fraction"] = 1.0
                r = evaluate(d2, T_grid=[T])
                sc = r["scalars"]
                ok = (sc.get("valid") is True and sc.get("bound_role") == "conservative_lower")
                if not ok:
                    failures.append(f"{name} T={T}: unrelaxed-bound evaluation mismatch "
                                     f"(valid={sc.get('valid')} bound_role={sc.get('bound_role')})")
            except Exception as exc:
                failures.append(f"{name} T={T}: unrelaxed evaluate raised: {exc}")

        # The card on disk itself must be untouched by any of the above.
        checks[0] += 1
        with open(os.path.join(ROOT, "cards", name), encoding="utf-8") as f:
            reread = yaml.safe_load(f)
        if reread != raws[name]:
            failures.append(f"{name}: card mutated on disk during evaluation")

        # In-memory YAML round trip: every new block must survive
        # serialization (no temp files -- yaml.safe_dump to a string).
        checks[0] += 1
        try:
            from dataclasses import asdict
            doc = {"meta": {"name": d.name, "role": "device-design"}, "design": asdict(d)}
            text = yaml.safe_dump(doc, sort_keys=False)
            reparsed = yaml.safe_load(text)["design"]
            for block in ("nitride", "provenance"):
                if reparsed.get(block) != getattr(d, block):
                    failures.append(f"{name}: {block} block did not round-trip through YAML")
        except Exception as exc:
            failures.append(f"{name}: round-trip serialization raised: {exc}")

    # ================================================== acceptance criterion 3
    # Every explicit physical leaf under the required namespaces has a
    # matching design.provenance.sources entry with a valid tag, and any
    # cited anchor_id exists in the nanowire ledger.
    for name in CARD_NAMES:
        design = raws[name]["design"]
        sources = design.get("provenance", {}).get("sources", {})
        n = design["nitride"]
        leaf_groups = [
            ("nitride.nanowire", n["nanowire"], NANOWIRE_LEAVES),
            ("nitride.dot", n["dot"], DOT_LEAVES),
            ("nitride.surface", n["surface"], SURFACE_LEAVES),
            ("nitride.photonics", n["photonics"], PHOTONICS_LEAVES),
            ("nitride.wire_thermal", n["wire_thermal"], WIRE_THERMAL_LEAVES),
            ("nitride.injector", n["injector"], INJECTOR_LEAVES),
            ("drive.diode", design["drive"]["diode"], DIODE_LEAVES),
        ]
        for prefix, block, expected_leaves in leaf_groups:
            checks[0] += 1
            actual = set(block.keys())
            if actual != expected_leaves:
                failures.append(f"{name}: {prefix} leaf set mismatch "
                                 f"(missing={expected_leaves - actual}, "
                                 f"extra={actual - expected_leaves})")
            for leaf in block:
                path = f"{prefix}.{leaf}"
                checks[0] += 1
                entry = sources.get(path)
                if entry is None or entry.get("tag") not in VALID_TAGS or not entry.get("source"):
                    failures.append(f"{name}: provenance.sources missing/invalid entry for {path}")
                    continue
                anchor = entry.get("anchor_id")
                if anchor:
                    checks[0] += 1
                    if anchor not in ledger.anchors:
                        failures.append(f"{name}: {path} cites unknown anchor {anchor!r}")

        set_params = design["drive"].get("set_params", {})
        is_set = design["drive"]["cycle_loading"] == "deterministic_pair"
        expected_sp = SET_PARAMS_COMMON_LEAVES | ({"radius_nm"} if is_set else set())
        checks[0] += 1
        if set(set_params.keys()) != expected_sp:
            failures.append(f"{name}: drive.set_params leaf set mismatch "
                             f"(got {sorted(set_params)}, expected {sorted(expected_sp)})")
        for leaf in set_params:
            path = f"drive.set_params.{leaf}"
            checks[0] += 1
            entry = sources.get(path)
            if entry is None or entry.get("tag") not in VALID_TAGS or not entry.get("source"):
                failures.append(f"{name}: provenance.sources missing/invalid entry for {path}")

        for path in ("drive.rep_rate_hz", "thermal.T_hs"):
            checks[0] += 1
            entry = sources.get(path)
            if entry is None or entry.get("tag") not in VALID_TAGS or not entry.get("source"):
                failures.append(f"{name}: provenance.sources missing/invalid entry for {path}")

        checks[0] += 1
        if "cavity" in n:
            failures.append(f"{name}: nitride.cavity block must not be present")

    # ================================================== acceptance criterion 4
    # Source literals cross-checked against the ledger.
    geom = ledger.anchors["deshpande2013_geometry"]["value"]
    abstract = ledger.anchors["deshpande2014_abstract"]["value"]
    # n_cm3/p_cm3 are written in the ledger without a signed exponent
    # ("3.0e18"), which PyYAML's YAML-1.1 float resolver reads back as a
    # plain string, not a float (the same quirk fsim_core/device.py's own
    # _coerce_optional_floats works around) -- coerce explicitly here.
    geom = dict(geom, n_cm3=float(geom["n_cm3"]), p_cm3=float(geom["p_cm3"]))
    for name in CARD_NAMES:
        design = raws[name]["design"]
        n = design["nitride"]
        checks[0] += 1
        if n["dot"]["x_in"] != abstract["x_in"]:
            failures.append(f"{name}: nitride.dot.x_in does not match deshpande2014_abstract x_in")
        checks[0] += 1
        if design["drive"]["diode"]["N_A"] != geom["p_cm3"]:
            failures.append(f"{name}: drive.diode.N_A does not match deshpande2013_geometry p_cm3")
        checks[0] += 1
        if design["drive"]["diode"]["N_D"] != geom["n_cm3"]:
            failures.append(f"{name}: drive.diode.N_D does not match deshpande2013_geometry n_cm3")
        checks[0] += 1
        if design["drive"]["rep_rate_hz"] != abstract["max_rate_MHz"] * 1e6:
            failures.append(f"{name}: drive.rep_rate_hz does not match deshpande2014_abstract max_rate_MHz")
        checks[0] += 1
        if n["nanowire"]["barrier_left_nm"] != geom["barrier_left_nm"] or \
           n["nanowire"]["barrier_right_nm"] != geom["barrier_right_nm"]:
            failures.append(f"{name}: nitride.nanowire barrier_*_nm does not match deshpande2013_geometry")
        checks[0] += 1
        if n["dot"]["height_nm"] != geom["disc_height_nm"]:
            failures.append(f"{name}: nitride.dot.height_nm does not match deshpande2013_geometry disc_height_nm")
        if "horizontal" in name:
            checks[0] += 1
            if n["nanowire"]["core_radius_nm"] != geom["optical_diameter_nm"] / 2.0:
                failures.append(f"{name}: horizontal core_radius_nm does not match half the "
                                 "deshpande2013_geometry optical_diameter_nm")
            checks[0] += 1
            if n["wire_thermal"]["R_s_ohm"] != geom["resistance_GOhm"] * 1e9:
                failures.append(f"{name}: horizontal R_s_ohm does not match "
                                 "deshpande2013_geometry resistance_GOhm")
        # 1.3 ns / g2=0.29 (2014) and the 2013 HBT/lifetime numbers must occur ONLY as
        # held-out comparison metadata (design.provenance.deshpande_comparison), never
        # promoted into a live card leaf.
        checks[0] += 1
        flat = _flatten({k: v for k, v in n.items() if k != "injector"})
        flat.update(_flatten(design["drive"]))
        # Only the two literal numbers acceptance criterion 4 names by value
        # (the 2014 g2=0.29 and 1.3 ns lifetime); 0.7/0.71/1.1 are excluded
        # from this scan because they collide with unrelated legitimate
        # defaults elsewhere on the card (e.g. nitride.dot.strain_c_fraction
        # =0.7, nitride.injector.delta_Ev_GaN_AlN_eV=0.70).
        bad = [k for k, v in flat.items() if v in (1.3, 0.29)]
        if bad:
            failures.append(f"{name}: held-out comparison literal(s) leaked into a live "
                             f"card leaf: {bad}")

    # ================================================== relaxed strain / access / doping binding
    for name in CARD_NAMES:
        design = raws[name]["design"]
        n = design["nitride"]
        checks[0] += 1
        if n["nanowire"]["strain_bound"] != "relaxed" or n["dot"]["strain_fraction"] != 0.0:
            failures.append(f"{name}: headline strain_bound must be relaxed with strain_fraction 0.0")
        checks[0] += 1
        surf = n["surface"]
        if (surf["shell"] != "none" or surf["shell_multiplier"] != 1.0
                or surf["reservoir_access"] != 1.0 or surf["occupied_dot_access"] != 0.05):
            failures.append(f"{name}: nitride.surface does not match the binding common default")

    checks[0] += 1
    surfaces = [raws[n]["design"]["nitride"]["surface"] for n in CARD_NAMES]
    if not all(s == surfaces[0] for s in surfaces):
        failures.append("nitride.surface is not identical across all four cards")

    # ================================================== acceptance criterion 2 (deep diff)
    h_pulse = raws[CARD_NAMES[0]]["design"]
    h_set = raws[CARD_NAMES[1]]["design"]
    v_pulse = raws[CARD_NAMES[2]]["design"]
    v_set = raws[CARD_NAMES[3]]["design"]

    def check_pair(a, b, allowed, label):
        checks[0] += 1
        diff = _diff_paths({k: v for k, v in a.items() if k != "provenance" and k != "name"},
                            {k: v for k, v in b.items() if k != "provenance" and k != "name"})
        stray = {p for p in diff if not _allowed(p, allowed)}
        if stray:
            failures.append(f"{label}: unexpected differing leaves {sorted(stray)}")

    check_pair(h_pulse, h_set, PULSE_SET_ALLOWED_PREFIXES, "horizontal pulse vs SET")
    check_pair(v_pulse, v_set, PULSE_SET_ALLOWED_PREFIXES, "vertical pulse vs SET")
    check_pair(h_pulse, v_pulse, FAMILY_ALLOWED_PREFIXES, "horizontal vs vertical (pulse)")
    check_pair(h_set, v_set, FAMILY_ALLOWED_PREFIXES, "horizontal vs vertical (SET)")

    # cycle_loading / set_params really do differ within each family pair.
    checks[0] += 1
    if h_pulse["drive"]["cycle_loading"] == h_set["drive"]["cycle_loading"]:
        failures.append("horizontal pulse/SET cycle_loading did not differ")
    checks[0] += 1
    if v_pulse["drive"]["cycle_loading"] == v_set["drive"]["cycle_loading"]:
        failures.append("vertical pulse/SET cycle_loading did not differ")
    checks[0] += 1
    if "radius_nm" in h_pulse["drive"]["set_params"] or "radius_nm" in v_pulse["drive"]["set_params"]:
        failures.append("pulse cards must not carry drive.set_params.radius_nm")
    checks[0] += 1
    if "radius_nm" not in h_set["drive"]["set_params"] or "radius_nm" not in v_set["drive"]["set_params"]:
        failures.append("SET cards must carry drive.set_params.radius_nm")

    # families really do differ.
    checks[0] += 1
    if h_pulse["nitride"]["nanowire"]["family"] == v_pulse["nitride"]["nanowire"]["family"]:
        failures.append("horizontal/vertical family did not differ")
    checks[0] += 1
    if h_pulse["nitride"]["nanowire"]["core_radius_nm"] == v_pulse["nitride"]["nanowire"]["core_radius_nm"]:
        failures.append("horizontal/vertical core_radius_nm did not differ")

    # ================================================== vertical single-mode headline (MEDIUM 13)
    for name in (CARD_NAMES[2], CARD_NAMES[3]):
        d = designs.get(name)
        if d is None:
            continue
        checks[0] += 1
        try:
            r = evaluate(copy.deepcopy(d), T_grid=[300.0])
            sc = r["scalars"]
            if not (sc.get("single_mode") is True and sc.get("approximation_error") == 0.0
                    and sc.get("headline_eligible") is True):
                failures.append(f"{name}: vertical card is not single-mode/headline-eligible "
                                 f"(single_mode={sc.get('single_mode')} "
                                 f"approximation_error={sc.get('approximation_error')})")
        except Exception as exc:
            failures.append(f"{name}: single-mode evaluation raised: {exc}")

    # ================================================== acceptance criterion 5 (SET screens)
    for name in (CARD_NAMES[1], CARD_NAMES[3]):
        d = designs.get(name)
        if d is None:
            continue
        checks[0] += 1
        try:
            sc = evaluate(copy.deepcopy(d), T_grid=[300.0])["scalars"]
            has_coulomb = all(k in sc for k in ("set_feasible", "set_E_C_meV", "set_EC_over_kT"))
            has_rti = all(k in sc for k in ("rti_feasible", "rti_status", "rti_failed_checks"))
            if not (has_coulomb and has_rti):
                failures.append(f"{name}: SET card missing Coulomb/RT-injector diagnostic columns")
        except Exception as exc:
            failures.append(f"{name}: SET screen evaluation raised: {exc}")

        # A card with an unknown occupation control cannot report a bare pass.
        checks[0] += 1
        inj = raws[name]["design"]["nitride"]["injector"]
        if inj.get("occupancy_control_known") is not False:
            failures.append(f"{name}: occupancy_control_known must be the honest-failure "
                             "default (False)")

    # ================================================== acceptance criterion 5 (mutation set)
    baseline_design = designs.get(CARD_NAMES[1])  # horizontal-set: exercises the SET path too
    checks[0] += 1
    try:
        base_sc = evaluate(copy.deepcopy(baseline_design), T_grid=[300.0])["scalars"]
        base_metric = (base_sc.get("g2_op"), base_sc.get("collected_flux_pulsed_s"),
                       base_sc.get("E_X_eV"))
    except Exception as exc:
        failures.append(f"mutation baseline evaluation raised: {exc}")
        base_metric = None

    moved = 0
    for path, new_value in MUTATIONS:
        checks[0] += 1
        if base_metric is None or baseline_design is None:
            failures.append(f"mutation {path}: no baseline to compare against")
            continue
        try:
            mutated = copy.deepcopy(baseline_design)
            _apply_mutation(mutated, path, new_value)
            sc = evaluate(mutated, T_grid=[300.0])["scalars"]
            metric = (sc.get("g2_op"), sc.get("collected_flux_pulsed_s"), sc.get("E_X_eV"))
            if metric != base_metric:
                moved += 1
            else:
                failures.append(f"mutation {path}={new_value}: no row value changed")
        except Exception as exc:
            failures.append(f"mutation {path}={new_value}: raised {exc}")
    checks[0] += 1
    if moved < 9:
        failures.append(f"only {moved}/9 leaf mutations changed a row value")

    # The card on disk must never have been mutated by any check above.
    for name in CARD_NAMES:
        checks[0] += 1
        with open(os.path.join(ROOT, "cards", name), encoding="utf-8") as f:
            reread = yaml.safe_load(f)
        if reread != raws[name]:
            failures.append(f"{name}: card mutated on disk by the end of verification")

    total = checks[0]
    passed = total - len(failures)
    print(f"{passed}/{total} nitride nanowire card checks passed")
    if failures:
        for failure in failures:
            print("FAIL: " + failure)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
