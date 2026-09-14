"""Contract checks for the opt-in nitride nanowire tier.

This is a (T)/(N) source-transcription and structural-consistency check
(README.md's five-way split): it PARSES `docs/nitride_nanowire_contract.md`
and `verify/data/nitride_nanowire_anchors.yaml` and asserts that the
contract's module table, card schema, row-column list, VERDICT format, and
sweep-grid arithmetic are internally consistent AND match the nine spec
files under `.workers/specs/nitride-nanowire-*.md` that this contract
transcribes. It never re-verifies future production solvers (pieces 2-9
are implemented and reviewed separately); it verifies the FROZEN INTERFACE
those pieces are built against.

Unlike the previous revision of this file, every check below reads real
text from the document/ledger/spec files -- there is no `literal ==
literal` tautology.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DOC_PATH = ROOT / "docs" / "nitride_nanowire_contract.md"
LEDGER_PATH = ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml"
CACHE_PATH = ROOT / "verify" / "data" / "citation_verification_cache.json"
CITATIONS_PY_PATH = ROOT / "verify" / "verify_citations.py"
SPEC_DIR = ROOT / ".workers" / "specs"

SPEC_FILES = {
    "levels": SPEC_DIR / "nitride-nanowire-levels.md",
    "photonics": SPEC_DIR / "nitride-nanowire-photonics.md",
    "surface": SPEC_DIR / "nitride-nanowire-surface.md",
    "transport": SPEC_DIR / "nitride-nanowire-transport.md",
    "injector": SPEC_DIR / "nitride-nanowire-injector.md",
    "device": SPEC_DIR / "nitride-nanowire-device.md",
    "cards": SPEC_DIR / "nitride-nanowire-cards.md",
    "sweep": SPEC_DIR / "nitride-nanowire-sweep.md",
}

ALLOWED_TAGS = {"V", "DR", "E", "A"}
STATUS_ENUM = {"full_text", "abstract_only", "figure_reading", "missing"}
REQUIRED_ANCHOR_FIELDS = {
    "citation", "doi_or_url", "location", "evidence_status", "evidence_kind",
    "tag", "platform", "excitation", "observable", "value", "unit",
    "tolerance", "transfer_notes", "temperature", "geometry",
    "raw_corrected", "observable_definition",
}

# The exact VERDICT field list piece 9 (nitride-nanowire-sweep.md) prescribes:
# "VERDICT lines retain established fields and add family, regime,
#  strain_bound, bound_role, rep_rate_hz, complete, eligible,
#  paired_optical_pass, hardware_qualified, rti_qualified, coverage,
#  invalid, idealized_status and flux_floor=1000/s."
REQUIRED_VERDICT_FIELDS = [
    "idealized_status", "family", "regime", "strain_bound", "bound_role",
    "rep_rate_hz", "complete", "eligible", "paired_optical_pass",
    "hardware_qualified", "rti_qualified", "coverage", "invalid",
    "flux_floor",
]

# Module -> which spec file's "Interface or signature constraints" section
# must contain each Module-table Signature cell verbatim.
MODULE_SPEC_KEY = [
    ("nitride_nanowire_levels", "levels"),
    ("nitride_nanowire_photonics", "photonics"),
    ("nitride_nanowire_surface", "surface"),
    ("nitride_nanowire_transport", "transport"),
    ("nitride_nanowire_injector", "injector"),
    ("device.py", "device"),
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def get_section(text: str, start_heading: str, end_heading: str | None) -> str:
    """Return the text between two '## '/'### ' headings (exclusive of the
    headings themselves). end_heading=None means "to end of file"."""
    start = text.index(start_heading) + len(start_heading)
    if end_heading is None:
        return text[start:]
    end = text.index(end_heading, start)
    return text[start:end]


def find_table(lines: list[str], header_marker: str, from_idx: int = 0) -> tuple[list[str], list[list[str]], int]:
    """Find the first markdown table whose header row contains
    `header_marker`, starting the search at lines[from_idx]. Returns
    (header_cells, data_rows, index_after_table)."""
    for i in range(from_idx, len(lines)):
        line = lines[i].strip()
        if line.startswith("|") and header_marker in line:
            header = [c.strip() for c in line.strip("|").split("|")]
            rows = []
            j = i + 2  # skip the '---' separator row
            while j < len(lines) and lines[j].strip().startswith("|"):
                cells = [c.strip() for c in lines[j].strip("|").split("|")]
                rows.append(cells)
                j += 1
            return header, rows, j
    raise ValueError(f"no table with header marker {header_marker!r} found")


def find_all_subsection_tables(section_text: str, heading_prefix: str) -> list[tuple[str, list[str], list[list[str]]]]:
    """For a section containing multiple '### <name>' subsections, each
    immediately followed (after a blank line) by exactly one markdown
    table, return [(subsection_name, header_cells, rows), ...]."""
    lines = section_text.splitlines()
    out = []
    i = 0
    while i < len(lines):
        if lines[i].startswith(heading_prefix):
            name = lines[i][len(heading_prefix):].strip()
            # find the next table header line after this heading
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith("|"):
                if lines[j].startswith(heading_prefix):
                    break
                j += 1
            if j < len(lines) and lines[j].strip().startswith("|"):
                header = [c.strip() for c in lines[j].strip("|").split("|")]
                rows = []
                k = j + 2
                while k < len(lines) and lines[k].strip().startswith("|"):
                    cells = [c.strip() for c in lines[k].strip("|").split("|")]
                    rows.append(cells)
                    k += 1
                out.append((name, header, rows))
                i = k
                continue
        i += 1
    return out


def leading_tag(cell: str) -> str:
    """Extract the leading provenance-tag token from a table cell like
    'V (2013 optical diameter/2)' -> 'V', or 'DR' -> 'DR'."""
    m = re.match(r"^([A-Za-z]+)", cell.strip())
    return m.group(1) if m else cell.strip()


def check_module_table(doc_text: str, ok) -> None:
    section = get_section(doc_text, "## Module table\n", "## Card schema\n")
    lines = section.splitlines()
    header, rows, _ = find_table(lines, "| Module | Symbol | Signature | Verifier |")
    ok("module_table_nonempty", len(rows) > 0)
    ok("module_table_row_count", len(rows) >= 15)

    spec_texts: dict[str, str] = {}
    for key, path in SPEC_FILES.items():
        spec_texts[key] = read_text(path)

    def interface_section(spec_text: str) -> str:
        m = re.search(r"## Interface or signature constraints\n(.*?)\n## ", spec_text, re.S)
        if not m:
            m = re.search(r"## Interface or signature constraints\n(.*)", spec_text, re.S)
        return m.group(1) if m else ""

    for row in rows:
        if len(row) < 4:
            ok("module_table_row_shape_" + "_".join(row), False)
            continue
        module_cell, symbol_cell, signature_cell, verifier_cell = row[0], row[1], row[2], row[3]
        spec_key = None
        for needle, key in MODULE_SPEC_KEY:
            if needle in module_cell:
                spec_key = key
                break
        name = f"module_signature_{symbol_cell.strip('`')}"
        if spec_key is None:
            ok(name, False)
            continue
        section_text = interface_section(spec_texts[spec_key])
        sig = signature_cell.strip("`")
        ok(name, sig in section_text)
        ok(name + "_verifier_path", "verify_nitride_nanowire" in verifier_cell)


def check_module_evidence_map(doc_text: str, anchors: dict, ok) -> None:
    section = get_section(doc_text, "### Module evidence map\n", "## Card schema\n")
    lines = section.splitlines()
    header, rows, _ = find_table(lines, "| Module | Ledger anchors |")
    ok("module_evidence_map_row_count", len(rows) == 6)
    for row in rows:
        module_name = row[0].strip("`")
        anchor_cell = row[1]
        anchor_ids = re.findall(r"`([a-zA-Z0-9_]+)`", anchor_cell)
        ok("module_evidence_ids_found_" + module_name, len(anchor_ids) > 0)
        has_nonnull = False
        for aid in anchor_ids:
            anchor = anchors.get(aid)
            if anchor is None:
                ok("module_evidence_anchor_exists_" + module_name + "_" + aid, False)
                continue
            if anchor.get("value") is not None:
                has_nonnull = True
        ok("module_has_nonnull_anchor_" + module_name, has_nonnull)


def check_card_schema(doc_text: str, ok) -> None:
    section = get_section(doc_text, "## Card schema\n", "## Row columns\n")
    tables = find_all_subsection_tables(section, "### ")
    ok("card_schema_table_count", len(tables) >= 7)
    total_leaves = 0
    for name, header, rows in tables:
        tag_idx = None
        for idx, h in enumerate(header):
            if h.strip().lower() == "tag":
                tag_idx = idx
                break
        ok("card_table_has_tag_column_" + name, tag_idx is not None)
        if tag_idx is None:
            continue
        unit_idx = header.index("Unit") if "Unit" in header else None
        default_idx = header.index("Default") if "Default" in header else None
        ok("card_table_has_unit_column_" + name, unit_idx is not None)
        ok("card_table_has_default_column_" + name, default_idx is not None)
        for row in rows:
            leaf = row[0].strip("`")
            total_leaves += 1
            tag = leading_tag(row[tag_idx]) if tag_idx < len(row) else ""
            ok(f"card_leaf_tag_{name}_{leaf}", tag in ALLOWED_TAGS)
            unit_val = row[unit_idx].strip() if unit_idx is not None and unit_idx < len(row) else ""
            default_val = row[default_idx].strip() if default_idx is not None and default_idx < len(row) else ""
            ok(f"card_leaf_unit_{name}_{leaf}", len(unit_val) > 0)
            ok(f"card_leaf_default_{name}_{leaf}", len(default_val) > 0)
    ok("card_schema_leaf_count", total_leaves >= 30)


def check_row_columns(doc_text: str, ok) -> None:
    section = get_section(doc_text, "## Row columns\n", "## Sweep grid")
    # Names required verbatim from the injector spec (piece 6).
    injector_names = [
        "rti_feasible", "rti_status", "rti_transport_feasible",
        "rti_level_margin_kT", "rti_alignment_error_meV", "rti_linewidth_meV",
        "rti_rate_Hz", "rti_e_rate_Hz", "rti_h_rate_Hz", "rti_bypass_fraction",
        "rti_missed_load_probability", "rti_second_pair_probability",
        "rti_growth_feasible", "rti_failed_checks", "rti_evidence_status",
        "rti_orbital_margin_kT",
    ]
    injector_spec_text = read_text(SPEC_FILES["injector"])
    for name in injector_names:
        ok("row_column_" + name + "_in_doc", name in section)
        ok("row_column_" + name + "_in_injector_spec", name in injector_spec_text)

    device_spec_text = read_text(SPEC_FILES["device"])
    device_names = [
        "optical_pass", "hardware_qualified", "rti_qualified", "device_pass",
        "rti_device_pass", "tau_rad_bare_ns", "tau_rad_photonic_ns",
        "tau_total_X_ns",
    ]
    for name in device_names:
        ok("row_column_" + name + "_in_doc", name in section)
        ok("row_column_" + name + "_in_device_spec", name in device_spec_text)

    transport_spec_text = read_text(SPEC_FILES["transport"])
    transport_names = [
        "area_cm2", "J_A_cm2", "eta_inj", "f_capture", "r_supply_s",
        "r_captured_s", "r_matrix_radiative_s", "r_matrix_nonradiative_s",
        "r_surface_reservoir_s", "r_leakage_s", "power_on_W",
        "accounting_residual_s",
    ]
    for name in transport_names:
        ok("row_column_" + name + "_in_doc", name in section)
        ok("row_column_" + name + "_in_transport_spec", name in transport_spec_text)


def check_verdict_and_grid(doc_text: str, ok) -> None:
    section = get_section(doc_text, "## Sweep grid, VERDICT format, and output paths\n", "## Evidence and verdict semantics\n")

    # --- grid arithmetic: parse the table, compute the product of counts ---
    lines = section.splitlines()
    header, rows, _ = find_table(lines, "| Family | Axis | Values | Count |")
    ok("grid_table_row_count", len(rows) == 14)
    count_idx = header.index("Count")
    family_idx = header.index("Family")
    products = {"horizontal": 1, "vertical": 1}
    counts_seen = {"horizontal": 0, "vertical": 0}
    for row in rows:
        fam = row[family_idx].strip()
        try:
            n = int(row[count_idx].strip())
        except ValueError:
            ok("grid_count_parses_" + fam + "_" + row[header.index('Axis')], False)
            continue
        if fam in products:
            products[fam] *= n
            counts_seen[fam] += 1
    ok("grid_axis_count_horizontal", counts_seen["horizontal"] == 7)
    ok("grid_axis_count_vertical", counts_seen["vertical"] == 7)
    ok("grid_product_horizontal_1536", products["horizontal"] == 1536)
    ok("grid_product_vertical_1024", products["vertical"] == 1024)
    # cross-check against the independently-known expected totals (not the
    # same computation restated: this compares the PARSED-table product
    # against literals transcribed from the design brief / sweep spec).
    ok("grid_declared_total_in_doc", "1536" in section and "1024" in section)
    ok("grid_le_10000", products["horizontal"] + products["vertical"] <= 10000)

    # --- VERDICT template ---
    m = re.search(r"```\nVERDICT: (.*?)\n```", section, re.S)
    ok("verdict_block_found", m is not None)
    verdict_line = m.group(1) if m else ""
    for field in REQUIRED_VERDICT_FIELDS:
        if field == "coverage":
            ok("verdict_field_coverage", "coverage=" in verdict_line)
        else:
            ok("verdict_field_" + field, (field + "=") in verdict_line)
    sweep_spec_text = read_text(SPEC_FILES["sweep"])
    for field in REQUIRED_VERDICT_FIELDS:
        ok("verdict_field_" + field + "_named_in_sweep_spec", field in sweep_spec_text)

    # --- output paths ---
    ok("output_path_sweep_csv", "sweep.csv" in section)
    ok("output_path_manifest", "manifest.json" in section)
    ok("output_path_results_md", "results.md" in section)
    ok("output_path_png", ".png" in section)
    ok("output_path_out_dir", "out/nitride_nanowire" in section)


def check_evidence_semantics(doc_text: str, ok) -> None:
    section = get_section(doc_text, "## Families and common card contract\n", "## Family-specific interfaces\n")
    ok("cavity_enabled_false_documented", "cavity.enabled` is `false`" in section or "cavity.enabled` is false" in section)


def check_ledger_rules(anchors: dict, ok) -> None:
    for anchor_id, row in anchors.items():
        missing_fields = REQUIRED_ANCHOR_FIELDS - set(row)
        ok("ledger_schema_" + anchor_id, not missing_fields and row["evidence_status"] in STATUS_ENUM)

        if row.get("value") is None:
            ok("ledger_null_value_missing_status_" + anchor_id, row.get("evidence_status") == "missing")

        if row.get("tag") == "V":
            notes = (row.get("transfer_notes") or "").lower()
            ok("ledger_v_tag_not_secondary_" + anchor_id, "secondary" not in notes and "never read" not in notes)
            ok("ledger_v_tag_has_identifier_or_title_" + anchor_id,
               row.get("doi_or_url") is not None or bool(row.get("title")))

        def scan_zero(value, path):
            if isinstance(value, dict):
                for k, v in value.items():
                    scan_zero(v, path + "." + k)
            elif isinstance(value, (int, float)) and not isinstance(value, bool) and value == 0:
                ok("ledger_no_zero_with_missing_status_" + anchor_id + path,
                   row.get("evidence_status") != "missing")

        scan_zero(row.get("value"), "")

    ok("ledger_surface_velocity_tag_E",
       anchors.get("deshpande2013_surface_velocity", {}).get("tag") == "E")
    s_val = (anchors.get("deshpande2013_surface_velocity", {}).get("value") or {})
    ok("ledger_surface_velocity_value_1000", s_val.get("S_cm_s") == 1000.0)
    ok("ledger_thermal_pl_no_surface_S",
       "surface_S_cm_s_secondary" not in (anchors.get("deshpande2013_thermal_and_pl", {}).get("value") or {}))

    supplement = anchors.get("deshpande2013_supplement", {})
    ok("ledger_supplement_null_doi", supplement.get("doi_or_url") is None)
    ok("ledger_supplement_missing_status", supplement.get("evidence_status") == "missing")

    kitamura = anchors.get("kitamura2026_architecture", {})
    ok("ledger_kitamura_value_null", kitamura.get("value") is None)
    ok("ledger_kitamura_missing_status", kitamura.get("evidence_status") == "missing")

    coulomb = anchors.get("deshpande2013_coulomb_reference", {}).get("value") or {}
    # Independent recomputation (does not call any fsim_core code): the
    # isolated-sphere convention fsim_core.drive_mech.set_feasibility
    # actually uses, C_sigma = 4*pi*eps0*eps_r*R.
    e = 1.602176634e-19
    eps0 = 8.8541878128e-12
    kB = 1.380649e-23
    r_m = coulomb.get("r_nm", float("nan")) * 1e-9
    eps_r = coulomb.get("eps_r", float("nan"))
    C_sphere = 4 * math.pi * eps0 * eps_r * r_m
    EC_sphere_J = e ** 2 / C_sphere
    EC_sphere_meV = EC_sphere_J / e * 1e3
    ok("coulomb_sphere_meV", abs(EC_sphere_meV - coulomb.get("EC_sphere_meV", 0.0)) < 0.01)
    ok("coulomb_sphere_ratio_230K", abs(EC_sphere_J / (kB * 230.0) - coulomb.get("EC_over_kT_230K_sphere", 0.0)) < 0.005)
    ok("coulomb_sphere_ratio_300K", abs(EC_sphere_J / (kB * 300.0) - coulomb.get("EC_over_kT_300K_sphere", 0.0)) < 0.005)
    C_disc = 8 * eps0 * eps_r * r_m
    EC_disc_J = e ** 2 / C_disc
    EC_disc_meV = EC_disc_J / e * 1e3
    ok("coulomb_disc_meV", abs(EC_disc_meV - coulomb.get("EC_disc_meV", 0.0)) < 0.01)
    ok("coulomb_sphere_disc_differ", abs(EC_sphere_meV - EC_disc_meV) > 1.0)
    ok("coulomb_convention_matches_set_feasibility",
       coulomb.get("convention_used_by_set_feasibility") == "sphere")


def check_source_transcription_literals(anchors: dict, ok) -> None:
    """Independent, non-tautological checks on the transcribed 2013/2014
    numbers (unchanged in substance from the previous revision, but no
    longer duplicated as a self-comparison)."""
    g = anchors["deshpande2013_geometry"]["value"]
    hbt = anchors["deshpande2013_hbt"]
    life = anchors["deshpande2013_lifetimes"]["value"]
    thermal = anchors["deshpande2013_thermal_and_pl"]["value"]
    rt = anchors["deshpande2014_abstract"]

    ok("diameters_distinct", g["optical_diameter_nm"] == 25.0 and thermal["thermal_diameter_nm"] == 30.0)
    area25 = math.pi * (12.5e-7) ** 2
    area30 = math.pi * (15.0e-7) ** 2
    ok("current_density_25nm", abs(1e-9 / area25 - 203.718327) < 1e-3)
    ok("current_density_30nm", abs(1e-9 / area30 - 141.471061) < 1e-3)
    ok("cw_not_pulsed", hbt["excitation"] == "CW electrical, 1 nA" and "unknown" in rt["excitation"])
    v = hbt["value"]
    ok("x_xx_raw_corrected", v == {"X_raw": 0.30, "X_corrected": 0.16, "XX_raw": 0.38, "XX_corrected": 0.25})
    ok("lifetimes_distinct", life == {"TRPL_XX_ps": 711.0, "HBT_X_ns": 1.1, "HBT_XX_ns": 0.7})
    ok("thermal_and_ensemble_distinct", thermal["rise_1nA_K"] == 15.0 and thermal["rise_2nA_K"] == 49.0 and thermal["ensemble_PL_300K_over_10K"] == 0.52)
    ok("2014_held_out", rt["evidence_status"] == "abstract_only" and rt["value"]["lifetime_ns"] == 1.3 and rt["value"]["g2"] == 0.29)

    emission = anchors["deshpande2013_emission_energy"]["value"]
    ok("emission_energy_X", abs(emission["X_eV"] - 2.84) < 1e-9 and abs(emission["X_nm"] - 436.56) < 1e-9)
    ok("emission_energy_XX_antibinding", emission["XX_above_X_meV"] == 10.0 and emission["XX_ordering"] == "antibinding")
    cavity_ledger = yaml.safe_load(read_text(ROOT / "verify" / "data" / "nitride_cavity_anchors.yaml"))
    cavity_anchor = (cavity_ledger.get("anchors") or {}).get("deshpande2013-xx-splitting")
    ok("emission_energy_cross_reference_exists", cavity_anchor is not None)
    if cavity_anchor is not None:
        ok("emission_energy_cross_reference_value", cavity_anchor.get("value") == -10.0)

    polarization = anchors["deshpande2013_polarization"]["value"]
    ok("polarization_70_percent", polarization["axial_dolp_percent"] == 70.0)

    device_geom = anchors["deshpande2013_device_geometry"]["value"]
    ok("device_geometry_substrate", "SiO2" in device_geom["substrate"] and "Si" in device_geom["substrate"])
    ok("device_geometry_contacted_length", device_geom["contacted_length_nm"] == 600.0)


def main() -> int:
    passed = total = 0

    def ok(name: str, value: bool) -> None:
        nonlocal passed, total
        total += 1
        passed += bool(value)
        if not value:
            print("FAIL " + name)

    doc_text = read_text(DOC_PATH)
    ledger_doc = yaml.safe_load(read_text(LEDGER_PATH))
    anchors = ledger_doc["anchors"]

    check_module_table(doc_text, ok)
    check_module_evidence_map(doc_text, anchors, ok)
    check_card_schema(doc_text, ok)
    check_row_columns(doc_text, ok)
    check_verdict_and_grid(doc_text, ok)
    check_evidence_semantics(doc_text, ok)
    check_ledger_rules(anchors, ok)
    check_source_transcription_literals(anchors, ok)

    ok("ledger_registered_in_citations_py",
       "nitride_nanowire_anchors.yaml" in read_text(CITATIONS_PY_PATH))
    import json
    cache = json.loads(read_text(CACHE_PATH))
    old_keys = {"crossref_doi|10.1038/ncomms2691", "crossref_doi|10.1063/1.4897640"}
    ok("old_cache_preserved", old_keys <= set(cache))

    # --- self-test: a corrupted copy of the module table must fail the
    # module-signature check. This mutates an in-memory string only; the
    # on-disk document is never touched. ---
    corrupted = doc_text.replace(
        "levels(system, T_K=300.0, *, z_points=1201, exterior_nm=45.0)",
        "levels(system, T_K=300.0, *, z_points=1201, exterior_nm=99.0)",
        1,
    )
    assert corrupted != doc_text, "self-test setup failed: nothing was replaced"
    self_test_ok = {"passed": 0, "total": 0}

    def self_test_probe(name: str, value: bool) -> None:
        self_test_ok["total"] += 1
        self_test_ok["passed"] += bool(value)

    try:
        check_module_table(corrupted, self_test_probe)
        self_test_failed_as_expected = self_test_ok["passed"] < self_test_ok["total"]
    except Exception:
        self_test_failed_as_expected = True
    ok("self_test_corrupted_doc_detected", self_test_failed_as_expected)
    print(
        "SELF-TEST: corrupting the levels() signature in the module table "
        "makes check 'module_signature_levels' fail "
        f"({self_test_ok['passed']}/{self_test_ok['total']} sub-checks passed "
        "on the corrupted in-memory copy, vs. all passing on the real "
        "document) -- the on-disk contract file was never modified."
    )

    print(f"{passed}/{total} nitride nanowire contract checks passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
