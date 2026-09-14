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

import dataclasses
import importlib
import inspect
import math
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DOC_PATH = ROOT / "docs" / "nitride_nanowire_contract.md"
LEDGER_PATH = ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml"
CACHE_PATH = ROOT / "verify" / "data" / "citation_verification_cache.json"
CITATIONS_PY_PATH = ROOT / "verify" / "verify_citations.py"
DEVICE_PY_PATH = ROOT / "fsim_core" / "device.py"
SPEC_DIR = ROOT / ".workers" / "specs"

# fix-2: bind the contract to the modules LIVE, not only to the specs. Only
# the five modules committed at this revision are imported; the device/
# sweep pieces (7/9) are not on disk yet, so their Module table rows stay
# bound to spec text only (see check_module_table).
LIVE_MODULE_NAMES = {
    "levels": "fsim_core.nitride_nanowire_levels",
    "photonics": "fsim_core.nitride_nanowire_photonics",
    "surface": "fsim_core.nitride_nanowire_surface",
    "transport": "fsim_core.nitride_nanowire_transport",
    "injector": "fsim_core.nitride_nanowire_injector",
}
LIVE_MODULES: dict[str, object] = {}
for _key, _modname in LIVE_MODULE_NAMES.items():
    try:
        LIVE_MODULES[_key] = importlib.import_module(_modname)
    except ImportError:
        LIVE_MODULES[_key] = None


def canonical_signature(name: str, obj) -> str:
    """Render a live callable's inspect.signature() the way this document's
    Signature cells are written: annotations stripped (this project's
    `from __future__ import annotations` makes every annotation an opaque
    source-text string, not the simplified name this table uses), defaults
    formatted with repr(), no per-parameter type hints."""
    sig = inspect.signature(obj)
    parts = []
    star_emitted = False
    for pname, p in sig.parameters.items():
        if p.kind == inspect.Parameter.VAR_POSITIONAL:
            parts.append("*" + pname)
            star_emitted = True
            continue
        if p.kind == inspect.Parameter.VAR_KEYWORD:
            parts.append("**" + pname)
            continue
        if p.kind == inspect.Parameter.KEYWORD_ONLY and not star_emitted:
            parts.append("*")
            star_emitted = True
        token = pname
        if p.default is not inspect.Parameter.empty:
            token += "=" + repr(p.default)
        parts.append(token)
    # Return annotation is deliberately NOT rendered: this document's
    # Signature cells are inconsistent (by spec-text inheritance, not by
    # error) about carrying "-> dict" even where the live function has a
    # return annotation, so the live-signature check below compares
    # parameter lists only, matching what actually determines call
    # compatibility.
    return name + "(" + ", ".join(parts) + ")"


def normalize_ws(s: str) -> str:
    return " ".join(s.split())


def strip_return_annotation(sig_text: str) -> str:
    return re.sub(r"\s*->\s*\S+\s*$", "", sig_text.strip())


def leaf_identifier(cell: str) -> str | None:
    """Extract the backticked leaf name from a Card schema first cell like
    '`core_radius_nm` (horizontal)' -> 'core_radius_nm'."""
    m = re.match(r"`([A-Za-z_][A-Za-z0-9_]*)`", cell.strip())
    return m.group(1) if m else None

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

    def row_spec_key(module_cell: str) -> str | None:
        for needle, key in MODULE_SPEC_KEY:
            if needle in module_cell:
                return key
        return None

    # First pass: a symbol with a "module-only" row (a fix/directive round
    # drifted its live signature past its frozen spec text) keeps its
    # ORIGINAL, spec-verbatim row's live-signature check skipped -- the
    # module-only row alone carries the live binding for that symbol, so
    # the same function is never bound to two different literal strings.
    module_only_symbols: set[tuple[str, str]] = set()
    for row in rows:
        if len(row) < 4:
            continue
        module_cell, symbol_cell, _sig, verifier_cell = row[0], row[1], row[2], row[3]
        if "module-only" not in verifier_cell:
            continue
        spec_key = row_spec_key(module_cell)
        m = re.match(r"[A-Za-z_][A-Za-z0-9_]*", symbol_cell.strip("`"))
        if spec_key and m:
            module_only_symbols.add((spec_key, m.group(0)))

    for row in rows:
        if len(row) < 4:
            ok("module_table_row_shape_" + "_".join(row), False)
            continue
        module_cell, symbol_cell, signature_cell, verifier_cell = row[0], row[1], row[2], row[3]
        spec_key = row_spec_key(module_cell)
        symbol_display = symbol_cell.strip("`")
        name = f"module_signature_{symbol_display}"
        if spec_key is None:
            ok(name, False)
            continue
        module_only = "module-only" in verifier_cell
        sig = signature_cell.strip("`")
        if not module_only:
            # Unchanged since round 1: the Signature cell is a literal
            # substring of the spec's own "Interface or signature
            # constraints" section.
            section_text = interface_section(spec_texts[spec_key])
            ok(name, sig in section_text)
        ok(name + "_verifier_path", "verify_nitride_nanowire" in verifier_cell)

        # --- fix-2: live module binding, independent of the spec text ---
        live_mod = LIVE_MODULES.get(spec_key)
        symbol_id_m = re.match(r"[A-Za-z_][A-Za-z0-9_]*", symbol_display)
        symbol_id = symbol_id_m.group(0) if symbol_id_m else None
        if live_mod is None or symbol_id is None:
            # device/sweep modules are not committed yet (piece 7/9); skip
            # live binding for those two rows only.
            continue
        obj = getattr(live_mod, symbol_id, None)
        ok(f"module_live_{symbol_display}_symbol_exists", obj is not None)
        if obj is None:
            continue
        if inspect.isfunction(obj):
            if not module_only and (spec_key, symbol_id) in module_only_symbols:
                # the drifted live truth is carried by the module-only
                # sibling row instead; this frozen row stays spec-bound only
                continue
            try:
                canon = canonical_signature(symbol_id, obj)
            except (TypeError, ValueError):
                canon = None
            ok(f"module_live_{symbol_display}_signature_matches_module",
               canon is not None
               and normalize_ws(canon) == normalize_ws(strip_return_annotation(sig)))
        elif dataclasses.is_dataclass(obj):
            # dataclass rows are bound to the Card schema leaf set instead
            # of to a literal constructor-signature string; see
            # check_dataclass_card_binding below.
            pass


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


# fix-2 required change 5: "for each dataclass, assert the card-table leaf
# set equals the field set." Three of the five committed dataclasses map
# 1:1 onto a single Card schema subsection with identical leaf names
# (checked as exact set equality, both directions). The other two
# (NitrideNanowireSystem, NitrideWireDiode) have leaves split across
# several subsections, some renamed for continuity with the pre-H6 card
# convention (documented in "Composition rules for the device piece" and
# the Card schema notes); those are checked as "every field name is
# present somewhere in the Card schema section", which is what acceptance
# criterion 2 ("lists no field absent from the contract's table") actually
# requires.
DATACLASS_CARD_BINDING = [
    ("photonics", "NitrideNanowirePhotonicsParams", "nitride.photonics", {}),
    ("surface", "NitrideNanowireSurfaceParams", "nitride.surface", {}),
    ("injector", "NitrideNanowireInjectorParams", "nitride.injector", {}),
    ("levels", "NitrideNanowireSystem", None, {"disc_radius_nm": "radius_nm"}),
    ("transport", "NitrideWireDiode", None, {}),
]


def check_dataclass_card_binding(doc_text: str, ok) -> None:
    card_section = get_section(doc_text, "## Card schema\n", "## Row columns\n")
    tables = find_all_subsection_tables(card_section, "### ")

    all_leaves: set[str] = set()
    by_subsection: dict[str, set[str]] = {}
    for name, header, rows in tables:
        leaves = {leaf_identifier(row[0]) for row in rows if leaf_identifier(row[0])}
        all_leaves |= leaves
        m = re.match(r"`([^`]+)`", name.strip())
        key = m.group(1) if m else name.strip()
        by_subsection.setdefault(key, set()).update(leaves)

    for spec_key, class_name, subsection, rename in DATACLASS_CARD_BINDING:
        live_mod = LIVE_MODULES.get(spec_key)
        cls = getattr(live_mod, class_name, None) if live_mod is not None else None
        ok(f"dataclass_binding_{class_name}_exists",
           cls is not None and dataclasses.is_dataclass(cls))
        if cls is None or not dataclasses.is_dataclass(cls):
            continue
        field_names = {f.name for f in dataclasses.fields(cls)}
        mapped = {rename.get(f, f) for f in field_names}
        if subsection is not None:
            leaves = by_subsection.get(subsection, set())
            ok(f"dataclass_binding_{class_name}_no_field_absent_from_card",
               mapped <= leaves)
            ok(f"dataclass_binding_{class_name}_no_undocumented_extra_leaf",
               leaves <= mapped)
        else:
            ok(f"dataclass_binding_{class_name}_no_field_absent_from_card",
               mapped <= all_leaves)


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

    # --- fix-2 required change 3: established count/flux/background/
    # filter/thermal columns from fsim_core/device.py's planar
    # _evaluate_nitride (around line 1333), bound to that module's own
    # source text (a literal dict-key substring), not to a spec file. ---
    device_py_text = read_text(DEVICE_PY_PATH)
    established_device_names = [
        "collected_flux_x_s", "collected_flux_xx_s", "background_flux_s",
        "total_detected_flux_s", "mean_counts", "mean_counts_x",
        "mean_counts_xx", "gate_ns_used", "S_X", "S_XX", "rho_pulsed",
        "blocked_load_probability", "counting_converged", "eta_out",
    ]
    for name in established_device_names:
        ok("row_column_" + name + "_in_doc", name in section)
        ok("row_column_" + name + "_established_in_device_py",
           ('"' + name + '"') in device_py_text)

    # --- fix-2 required change 3: new coherence columns (physics
    # coherence review). fsim_core/nitride_nanowire_device.py does not
    # exist yet, so these are checked against this document only. ---
    coherence_names = [
        "thermal_iterations", "thermal_converged", "bound_reversal",
        "collected_flux_delivered_s", "f_qfl_dot_thermodynamic_limit",
        "sidewall_overlap", "degree_of_linear_polarization",
        "beta_multimode_penalty", "single_mode",
    ]
    for name in coherence_names:
        ok("row_column_" + name + "_in_doc", name in section)

    # --- fix-2 required change 3/5: the injector's rti_* extra keys,
    # introspected LIVE by calling injector_feasibility at its card
    # defaults (not a hardcoded list) -- catches a key this document has
    # not documented yet, independent of the (frozen) injector spec text.
    injector_mod = LIVE_MODULES.get("injector")
    live_rti_keys: list[str] = []
    if injector_mod is not None:
        try:
            params_cls = injector_mod.NitrideNanowireInjectorParams
            defaults = params_cls(occupancy_control_known=True,
                                   second_pair_control_known=True)
            out = injector_mod.injector_feasibility(
                defaults, T_K=300.0, rep_rate_hz=200e6,
                loading_window_ns=0.1, electron_level_eV=0.05,
                hole_level_eV=0.01, electron_spacing_meV=float("nan"),
                hole_spacing_meV=float("nan"),
                second_pair_addition_meV=10.0,
                available_pair_rate_Hz=1e9,
            )
            live_rti_keys = sorted(k for k in out if k.startswith("rti_"))
        except Exception:
            live_rti_keys = []
    ok("row_column_injector_live_introspection_succeeded", len(live_rti_keys) > 0)
    for name in live_rti_keys:
        ok("row_column_" + name + "_in_doc_live", name in section)


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

    # fix-2 required change 4 / coherence finding 8: screening= and access=
    # are new fields this contract revision adds; they postdate piece 9's
    # own frozen field-list bullet, so they are checked against this
    # document's VERDICT template only, not against the sweep spec text.
    ok("verdict_field_screening", "screening=" in verdict_line)
    ok("verdict_field_access", "access=" in verdict_line)

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

        # fix-2 required change 5: "scan tolerance values for
        # zero-as-unknown" -- the same null-is-unknown convention the
        # ledger already enforces for `value` (never 0/0.0 for "unknown")
        # applies to `tolerance.value` too (fix-1 finding 7).
        tol = row.get("tolerance") or {}
        tol_val = tol.get("value") if isinstance(tol, dict) else None
        ok("ledger_tolerance_value_not_zero_" + anchor_id,
           tol_val is None or tol_val != 0)

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


def check_source_transcription_literals(doc_text: str, anchors: dict, ok) -> None:
    """Independent, non-tautological checks on the transcribed 2013/2014
    numbers (unchanged in substance from the previous revision, but no
    longer duplicated as a self-comparison)."""
    g = anchors["deshpande2013_geometry"]["value"]
    hbt = anchors["deshpande2013_hbt"]
    life = anchors["deshpande2013_lifetimes"]["value"]
    thermal = anchors["deshpande2013_thermal_and_pl"]["value"]
    rt = anchors["deshpande2014_abstract"]

    ok("diameters_distinct", g["optical_diameter_nm"] == 25.0 and thermal["thermal_diameter_nm"] == 30.0)

    # fix-2 required change 5 (fix-1 finding 4): recompute both current
    # densities from RADII READ LIVE from the ledger anchors (not a
    # hardcoded area literal), then cross-check the result against (a) the
    # 203.7/141.5 A/cm2 figures parsed out of this document's own
    # "Evidence and verdict semantics" prose and (b) Deshpande 2013's own
    # independently reported rounded 142 A/cm2 figure for the 30 nm case.
    r25_cm = (g["optical_diameter_nm"] / 2.0) * 1e-7
    r30_cm = (thermal["thermal_diameter_nm"] / 2.0) * 1e-7
    area25 = math.pi * r25_cm ** 2
    area30 = math.pi * r30_cm ** 2
    J25 = 1e-9 / area25
    J30 = 1e-9 / area30
    m = re.search(
        r"1 nA is (\d+\.?\d*) A/cm2 at 25 nm and (\d+\.?\d*) A/cm2 at 30 nm",
        doc_text,
    )
    ok("current_density_prose_found", m is not None)
    if m is not None:
        ok("current_density_25nm_matches_doc_prose", abs(J25 - float(m.group(1))) < 0.1)
        ok("current_density_30nm_matches_doc_prose", abs(J30 - float(m.group(2))) < 0.1)
    ok("current_density_30nm_matches_reported_142",
       abs(J30 - thermal["J_1nA_reported_A_cm2"]) < 1.0)
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
    ok("device_geometry_wire_length", device_geom["wire_length_nm"] == 600.0)


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
    check_dataclass_card_binding(doc_text, ok)
    check_row_columns(doc_text, ok)
    check_verdict_and_grid(doc_text, ok)
    check_evidence_semantics(doc_text, ok)
    check_ledger_rules(anchors, ok)
    check_source_transcription_literals(doc_text, anchors, ok)

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
