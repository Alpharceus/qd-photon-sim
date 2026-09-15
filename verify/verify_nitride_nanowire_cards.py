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

Fix round (Opus review of 58eb5a2, HIGH 1 + MEDIUM 2/3/7/8): the previous
version of this file compared only leaf NAME sets and provenance
presence/anchor-existence against 8 hard-coded ledger literals, never a
leaf VALUE against the contract's own Card schema tables -- 8 of 12
Opus-planted single-leaf edits (b_res, I_uA, tau_rad0_ns, NA,
unguided_collection_scale, Rth_K_W, tau_cap_ps, set_params.eps_r) passed
556/556 unnoticed. This version parses
docs/nitride_nanowire_contract.md's "Card schema" markdown tables directly
(_parse_contract_card_schema below; see also DIRECT_DEFAULTS for the
handful of round-pinned leaves the contract states only in prose, never in
a table row: drive.b_res, drive.I_uA, nitride.tau_rad0_ns,
nitride.tau_cap_ps) and asserts every card leaf's value against the
contract's own default/range/tag for the card's own family, so a future
contract-table edit is picked up automatically instead of drifting out of
sync with a second, hand-maintained copy.

Fix round 3 (Opus re-review of 966d19a, HIGH verify:262 + 2 MEDIUM + 4 LOW):
(1) check_contract_leaf compared the contract's bare "horizontal"/
"vertical" leaf-cell suffix directly to the card's full family enum value
("horizontal_as_built"/"vertical_photonic"), so all six family-qualified
Card-schema rows always mismatched and silently skipped BOTH their value
and provenance-tag checks -- see _FAMILY_SUFFIX_TO_FAMILY. (2) drive.diode.
tau_pulse_ns and nitride.dot.k_intrinsic_ns had no value pin at all; both
are now in DIRECT_DEFAULTS. (3) the provenance-entry loop's anchor check
only fired when an anchor_id happened to be present (a deleted anchor_id
silently passed) and a missing/invalid entry `continue`d past it entirely
(shrinking the denominator); both are now unconditional, and a "V"-tagged
entry with no anchor_id is a failure unless the leaf cites a materials-
module constant, not a ledger anchor (LITERATURE_NO_ANCHOR_ALLOWLIST).
(4) nitride.photonics.dipole_weights is now written explicitly on all four
cards (isotropic, [A]; fsim_core/nitride_nanowire_device.py has coerced a
YAML list to the required tuple since 90ab4cf) and is in PHOTONICS_LEAVES;
check_contract_leaf special-cases it (a 3-tuple, not a scalar). (5)
dot.gamma300 and thermal.T_hs move no evaluate() scalar at this round's
T_grid, so DIRECT_DEFAULTS/DIRECT_DEFAULT_TAGS now pin their value and tag
explicitly rather than relying on the mutation-set check to catch a swap.
"""
import copy
import math
import os
import re
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

# Fix round 3 (Opus re-review of 966d19a, LOW verify:539-542): a "V"
# (verified/literature) provenance entry with no anchor_id is otherwise
# unaccountable -- the anchor check at verify:588-... only fires when an
# anchor_id string is PRESENT, so simply deleting one silently passed. The
# seven paths below are the only V-tagged leaves that cite a materials-
# module constant rather than a nanowire ledger anchor (Rinke 2008,
# Malitson 1965, Bernardini PRB 1997, Martin/Yu/Waldrop APL 1996, the two
# fsim_core.nitride_materials effective-mass constants, and the directive-
# round Mg acceptor ionization energy): a missing anchor_id on any OTHER
# V-tagged leaf is a failure, not a silent skip.
LITERATURE_NO_ANCHOR_ALLOWLIST = {
    "nitride.dot.vbo_InN_GaN_eV",             # Rinke et al., PRB 77, 075202 (2008)
    "nitride.photonics.n_oxide",              # Malitson (1965), sio2_index
    "nitride.injector.growth_step_nm",        # Bernardini, PRB 56, R10024 (1997)
    "nitride.injector.delta_Ev_GaN_AlN_eV",   # Martin, Yu, Waldrop, APL 68, 2541 (1996)
    "nitride.injector.me_well",               # fsim_core.nitride_materials
    "nitride.injector.mh_well",               # fsim_core.nitride_materials
    "nitride.injector.mg_acceptor_energy_meV",  # directive round e696bdd, Mg acceptor
}

# Card-schema leaf sets (docs/nitride_nanowire_contract.md "Card schema"),
# independently transcribed here (not imported from the device module) so
# this verifier catches a leaf silently added/removed from either side.
NANOWIRE_LEAVES = {"family", "core_radius_nm", "outer_radius_nm", "strain_bound",
                    "shell", "barrier_left_nm", "barrier_right_nm"}
DOT_LEAVES = {"radius_nm", "height_nm", "x_in", "strain_fraction", "screening_fraction",
              "external_field_kVcm", "vbo_InN_GaN_eV", "strain_c_fraction",
              "k_intrinsic_ns"}  # k_intrinsic_ns: device-module addendum, not a contract leaf
SURFACE_LEAVES = {"S_cm_s", "shell", "shell_multiplier", "reservoir_access", "occupied_dot_access"}
# nitride.photonics: all 19 of the dataclass's fields, including
# dipole_weights (fix round 3): fsim_core/nitride_nanowire_device.py has
# coerced a 3-element nitride.photonics.dipole_weights YAML list to a tuple
# of floats at the card boundary since 90ab4cf, so the leaf is reachable
# from a card and is now written explicitly (isotropic, [A]) instead of
# relying on the dataclass's own class default -- see
# design.provenance.sources['nitride.photonics.dipole_weights'] on each
# card.
PHOTONICS_LEAVES = {"family", "NA", "n_wire", "n_group_override", "n_ambient", "n_oxide",
                     "oxide_thickness_nm", "n_substrate", "dipole_weights", "emitter_height_nm",
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

# ============================================================ HIGH-1 fix
# Parse docs/nitride_nanowire_contract.md's own "Card schema" markdown
# tables (leaf, unit, default, range, tag[, notes]) instead of hand-copying
# a second schema that can silently drift from the doc. Only the eight
# named sub-tables below exist in that section; a leaf not covered by one
# of them (d_active_nm; drive.b_res/I_uA/nitride.tau_rad0_ns/tau_cap_ps,
# which the contract pins only in prose, never in a table row) is asserted
# separately -- see DIRECT_DEFAULTS and the b_res/I_uA/tau_rad0_ns/
# tau_cap_ps checks in main().
CONTRACT_PATH = os.path.join(ROOT, "docs", "nitride_nanowire_contract.md")
_SCHEMA_SECTION_HEADERS = {
    "`nitride.nanowire`": "nitride.nanowire",
    "`nitride.dot`": "nitride.dot",
    "`nitride.surface`": "nitride.surface",
    "`nitride.photonics`": "nitride.photonics",
    "`nitride.wire_thermal`": "nitride.wire_thermal",
    "`drive.diode`": "drive.diode",
    "`nitride.injector`": "nitride.injector",
    "`drive.set_params`": "drive.set_params",
}
_NUM_RE = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")


def _numbers_in(text):
    return [float(x) for x in _NUM_RE.findall(text)]


def _backtick_tokens(text):
    return re.findall(r"`([^`]+)`", text)


def _md_table_rows(block_text):
    lines = [l for l in block_text.splitlines() if l.strip().startswith("|")]
    if len(lines) < 2:
        return []
    return [[c.strip() for c in line.strip().strip("|").split("|")] for line in lines[2:]]


def parse_contract_card_schema(path=CONTRACT_PATH):
    """Return {section: [raw markdown table row cell-lists]} for every
    '### `section`' sub-table under docs/nitride_nanowire_contract.md's own
    '## Card schema' heading (independently re-parsed here every run, not
    imported from any cached copy)."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    section = text[text.index("## Card schema"):text.index("## Row columns")]
    out = {}
    for part in re.split(r"\n### ", section)[1:]:
        key = next((v for k, v in _SCHEMA_SECTION_HEADERS.items() if part.startswith(k)), None)
        if key:
            out[key] = _md_table_rows(part)
    return out


def _parse_leaf_cell(cell):
    m = re.match(r"`([^`]+)`(?:\s*\(([^)]+)\))?", cell)
    if not m:
        return cell.strip("` "), None
    return m.group(1), m.group(2)


def _parse_row_cells(section, cells):
    if section == "nitride.photonics":
        leaf_cell, _unit, default_cell, range_cell, tag_cell, _notes = cells
    else:
        leaf_cell, _unit, default_cell, range_cell, tag_cell = cells
    name, family_suffix = _parse_leaf_cell(leaf_cell)
    return name, family_suffix, default_cell, range_cell, tag_cell


def parse_contract_tags(tag_cell):
    """'V', 'V/E', or 'V (default)/A (override)' -> {'V'}, {'V','E'}, {'V','A'}."""
    toks = set()
    for part in tag_cell.split("/"):
        m = re.match(r"\s*([A-Z]{1,2})\b", part)
        if m:
            toks.add(m.group(1))
    return toks

# Leaf-name-keyed dynamic-default resolvers, transcribed one-for-one from
# the Default column's "equal to `X`" / "mirrors `X`" text (see the parsed
# cells printed by parse_contract_card_schema): the referenced leaf X is
# read from the SAME loaded card, never from a second hard-coded literal.
_EQUALITY_REFS = {
    ("nitride.dot", "radius_nm", "horizontal"): lambda nit: nit["nanowire"]["core_radius_nm"],
    ("drive.diode", "conducting_radius_nm", None): lambda nit: nit["nanowire"]["core_radius_nm"],
    ("drive.set_params", "radius_nm", None): lambda nit: nit["dot"]["radius_nm"],
    ("nitride.photonics", "family", None): lambda nit: nit["nanowire"]["family"],
}


def _outer_radius_expected(nit):
    core = nit["nanowire"]["core_radius_nm"]
    return core if nit["nanowire"]["shell"] == "none" else core + 3.0


def _strain_fraction_expected(nit, default_cell):
    m = re.search(r"([\d.]+)\s+unrelaxed.*?([\d.]+)\s+relaxed", default_cell)
    mapping = {"unrelaxed": float(m.group(1)), "relaxed": float(m.group(2))}
    return mapping[nit["nanowire"]["strain_bound"]]


def _shell_multiplier_expected(nit, default_cell):
    m_none = re.search(r"([\d.]+)\s*\(`none`\)", default_cell)
    m_alg = re.search(r"([\d.]+)\s*\(`AlGaN`\)", default_cell)
    return {"none": float(m_none.group(1)), "AlGaN": float(m_alg.group(1))}[nit["surface"]["shell"]]


def _family_split_numbers(cell):
    m_h = re.search(r"([-+0-9.eE]+)\s*\(horizontal\)", cell)
    m_v = re.search(r"([-+0-9.eE]+)\s*\(vertical\)", cell)
    return float(m_h.group(1)), float(m_v.group(1))


# The contract's own leaf-cell parenthetical spells the family qualifier
# "horizontal"/"vertical" (docs/nitride_nanowire_contract.md Card schema:
# "`core_radius_nm` (horizontal)"), while the card's own nitride.nanowire.
# family value is the full enum "horizontal_as_built"/"vertical_photonic".
# Fix round 3 (Opus re-review of 966d19a, HIGH verify:262): comparing the
# bare suffix to the full enum value directly always mismatched, so all six
# family-qualified Card-schema rows (core_radius_nm and R_s_ohm horizontal/
# vertical, dot.radius_nm horizontal/vertical) silently returned None and
# were skipped for BOTH the value check and the provenance-tag check below.
_FAMILY_SUFFIX_TO_FAMILY = {"horizontal": "horizontal_as_built", "vertical": "vertical_photonic"}


class SchemaParseError(Exception):
    """A Card-schema row this module has no resolver for (should never fire
    for the leaves the four cards actually write; a new leaf added to the
    contract without a matching resolver here raises loudly instead of
    silently skipping its value check)."""


def check_contract_leaf(section, name, family_suffix, default_cell, range_cell, *, nit, value, family):
    """Return None if this row's family_suffix does not apply to `family`,
    else (ok, note) asserting value equals the contract default (resolved
    for this card's own family/shell/strain_bound where the Default column
    is conditional) or lies inside an explicitly ENUMERATED contract range
    (a `{...}`/`sensitivity`/`sweep` cell). A continuous validity bound
    (`[0,1]`, `(0,1]`, `fixed`, ...) is not itself a licence to deviate from
    the default on a headline card -- only a discrete, explicitly named
    alternate value is."""
    dcell, rcell = default_cell.strip(), range_cell.strip()
    key = (section, name, family_suffix)
    if family_suffix is not None and _FAMILY_SUFFIX_TO_FAMILY.get(family_suffix, family_suffix) != family:
        return None

    # dipole_weights is a 3-tuple, not a scalar: none of the generic
    # numeric/enum parsing below applies to it (and the generic numeric
    # fallback would crash trying to hash an unhashable list into a set).
    # Both families share the same isotropic default and tag (contract Card
    # schema "both"), so this is a flat equality check, not a family split.
    if name == "dipole_weights":
        expected = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
        return tuple(value) == expected, "isotropic " + repr(expected)

    if dcell.startswith("equal to") or dcell.startswith("mirrors"):
        if name == "outer_radius_nm":
            expected = _outer_radius_expected(nit)
        elif key in _EQUALITY_REFS:
            expected = _EQUALITY_REFS[key](nit)
        elif (section, name, None) in _EQUALITY_REFS:
            expected = _EQUALITY_REFS[(section, name, None)](nit)
        else:
            raise SchemaParseError("no equality resolver for " + section + "." + name)
        return value == expected, "equal to " + repr(expected)
    if dcell.startswith("derived from"):
        if name == "strain_fraction":
            expected = _strain_fraction_expected(nit, dcell)
            return value == expected, "derived " + repr(expected)
        raise SchemaParseError("no derivation resolver for " + section + "." + name)
    if dcell == "required":
        allowed = set(_backtick_tokens(rcell))
        return value in allowed, "required, one of " + repr(allowed)
    if "(`none`)" in dcell and "(`AlGaN`)" in dcell:
        expected = _shell_multiplier_expected(nit, dcell)
        return value == expected, "shell-conditional default " + repr(expected)

    if "(horizontal)" in dcell and "(vertical)" in dcell:
        h, v = _family_split_numbers(dcell)
        default = h if family == "horizontal_as_built" else v
    elif dcell.startswith("None"):
        default = None
    elif re.match(r"^`(True|False)`", dcell):
        default = dcell.split("`")[1] == "True"
    elif re.match(r'^`"([^"]+)"`', dcell):
        default = re.match(r'^`"([^"]+)"`', dcell).group(1)
    elif re.match(r"^`([A-Za-z_][\w-]*)`", dcell):
        default = re.match(r"^`([A-Za-z_][\w-]*)`", dcell).group(1)
    else:
        nums = _numbers_in(dcell)
        default = nums[0] if nums else None

    m_fixed = re.search(r"fixed\s+(?:at\s+)?([-+0-9.eE]+)", rcell)
    if m_fixed:
        return value == float(m_fixed.group(1)), "fixed at " + m_fixed.group(1)
    if re.match(r"^fixed\b", rcell) or "fixed baseline" in rcell or "unchanged SET convention" in rcell:
        return value == default, "fixed, equal to default " + repr(default)
    if "{" in rcell:
        inner = re.search(r"\{([^}]*)\}", rcell).group(1)
        nums = _numbers_in(inner)
        if nums:
            allowed = set(nums) | ({default} if isinstance(default, (int, float)) else set())
            return value in allowed, "one of " + repr(allowed)
        toks = set(_backtick_tokens(inner)) | set(re.findall(r'"([^"]+)"', inner))
        toks = {True if t == "True" else False if t == "False" else t for t in toks}
        if default is not None:
            toks.add(default)
        return value in toks, "one of " + repr(toks)
    if ("sensitivity" in rcell or "sweep" in rcell) and not rcell.startswith("("):
        nums = _numbers_in(rcell)
        if nums:
            allowed = set(nums) | ({default} if isinstance(default, (int, float)) else set())
            return value in allowed, "one of " + repr(allowed)
        return value == default, "equal to default " + repr(default) + " (no enumerated sensitivity set)"
    if "one decade each way" in rcell:
        lo, hi = default / 10.0, default * 10.0
        return lo <= value <= hi, "within one decade of " + repr(default)
    if re.match(r"^,?\s*`", rcell) and "{" not in rcell:
        toks = {t for t in _backtick_tokens(rcell) if not re.search(r"[\[\](),]", t)}
        if default is not None:
            toks.add(default)
        return value in toks, "one of " + repr(toks) + " or equal to default " + repr(default)
    if rcell.startswith(">="):
        bound = _numbers_in(rcell)[0]
        return value >= bound, ">= " + repr(bound)
    if rcell.startswith(">"):
        head_nums = _numbers_in(rcell.split(",")[0])
        bound = head_nums[0] if head_nums else 0.0
        return value > bound, "> " + repr(bound)
    if rcell == "finite":
        return math.isfinite(value), "finite"
    if "positive if set" in rcell:
        return (value is None) or (value > 0), "None or > 0"
    return value == default, "equal to default " + repr(default)


SCHEMA_ROWS = parse_contract_card_schema()

# The handful of DriveBlock/nitride-top-level leaves this round pins to one
# exact value in prose (docs/nitride_nanowire_contract.md "Families and
# common card contract": "drive.b_res=0.1 [A]"; nitride.tau_rad0_ns=1.0;
# spec "Common explicit defaults": I_uA design choice = 2 nA = 2.0e-3 uA,
# tau_cap_ps=10) rather than a Card-schema table row, so parse_contract_
# card_schema above never sees them.
#
# Fix round 3 (Opus re-review of 966d19a, MEDIUM verify:107-108/73-75) adds
# two more: drive.diode.tau_pulse_ns (spec "Common explicit defaults":
# tau_pulse_ns=0.1) and nitride.dot.k_intrinsic_ns (a device-module
# addendum, not a Card-schema row -- see DOT_LEAVES above -- pinned to its
# neutral zero-extra-loss default); neither had a value pin before, so a
# mutated value applied to all four cards passed unnoticed.
DIRECT_DEFAULTS = {
    "drive.b_res": 0.1,
    "drive.I_uA": 2.0e-3,
    "nitride.tau_rad0_ns": 1.0,
    "nitride.tau_cap_ps": 10.0,
    "drive.diode.tau_pulse_ns": 0.1,
    "nitride.dot.k_intrinsic_ns": 0.0,
    "dot.gamma300": 35.0,
    "thermal.T_hs": 300.0,
}

# Fix round 3 (Opus re-review of 966d19a, LOW "gamma300 / thermal.T_hs
# unpinned (no scalar effect at T_grid)"): dot.gamma300 and thermal.T_hs
# move no evaluate() scalar the mutation-set check (acceptance criterion 5)
# samples, so DIRECT_DEFAULTS' value pin above is the only thing that can
# catch a mutated value; this dict additionally pins the provenance tag
# each of those two leaves carries on the card (dot.gamma300: A, design
# choice; thermal.T_hs: V, deshpande2014_abstract room-temperature
# condition), so a value-and-tag swap cannot pass either.
DIRECT_DEFAULT_TAGS = {
    "dot.gamma300": "A",
    "thermal.T_hs": "V",
}


class Ledger:
    def __init__(self):
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
    # LOW-11 fix: exact match or a real child path (dot-separated) only --
    # no bare startswith, which would also match e.g. "drive.set_paramsX"
    # against the allow-listed prefix "drive.set_params".
    return any(path == p or path.startswith(p + ".") for p in prefixes)


def main():
    checks = [0]
    failures = []

    ledger = Ledger()

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
                entry_ok = entry is not None and entry.get("tag") in VALID_TAGS and entry.get("source")
                if not entry_ok:
                    failures.append(f"{name}: provenance.sources missing/invalid entry for {path}")

                # MEDIUM-3 fix (kept): count this check unconditionally (not
                # only when anchor_id happens to be present) so deleting an
                # anchor_id line cannot shrink the failure count on the SAME
                # denominator (the previous verifier's "552/552" bug).
                # Fix round 3 (LOW verify:531-533): also run it
                # unconditionally when the entry itself is missing/invalid
                # (previously a bare `continue` here skipped this second
                # check entirely, shrinking the denominator by one relative
                # to a present-but-anchor-bad entry) and (LOW verify:539-
                # 542) require every "V"-tagged entry to cite a real
                # anchor_id unless the leaf is on the materials-constant
                # allow-list above.
                checks[0] += 1
                if not entry_ok:
                    failures.append(f"{name}: {path}: cannot check anchor rule, "
                                     "provenance entry missing/invalid")
                else:
                    anchor = entry.get("anchor_id")
                    if anchor:
                        if anchor not in ledger.anchors:
                            failures.append(f"{name}: {path} cites unknown anchor {anchor!r}")
                    elif entry.get("tag") == "V" and path not in LITERATURE_NO_ANCHOR_ALLOWLIST:
                        failures.append(f"{name}: {path} is V-tagged with no anchor_id and is "
                                         "not on the materials-constant allow-list")

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
    # MEDIUM-3 fix (the actual denominator-shrink case the finding names):
    # each of these ledger-governed leaves must ALSO cite the SAME anchor
    # its value is checked against above, so deleting the anchor_id line
    # (while leaving the value untouched) is its own separate failure, not
    # merely a value check that happens to still pass.
    def _anchor_ck(name, sources, path, expected_anchor):
        checks[0] += 1
        entry = sources.get(path)
        if entry is None or entry.get("anchor_id") != expected_anchor:
            failures.append(f"{name}: {path} does not cite anchor_id {expected_anchor!r}")

    for name in CARD_NAMES:
        design = raws[name]["design"]
        n = design["nitride"]
        sources = design.get("provenance", {}).get("sources", {})
        checks[0] += 1
        if n["dot"]["x_in"] != abstract["x_in"]:
            failures.append(f"{name}: nitride.dot.x_in does not match deshpande2014_abstract x_in")
        _anchor_ck(name, sources, "nitride.dot.x_in", "deshpande2014_abstract")
        checks[0] += 1
        if design["drive"]["diode"]["N_A"] != geom["p_cm3"]:
            failures.append(f"{name}: drive.diode.N_A does not match deshpande2013_geometry p_cm3")
        _anchor_ck(name, sources, "drive.diode.N_A", "deshpande2013_geometry")
        checks[0] += 1
        if design["drive"]["diode"]["N_D"] != geom["n_cm3"]:
            failures.append(f"{name}: drive.diode.N_D does not match deshpande2013_geometry n_cm3")
        _anchor_ck(name, sources, "drive.diode.N_D", "deshpande2013_geometry")
        checks[0] += 1
        if design["drive"]["rep_rate_hz"] != abstract["max_rate_MHz"] * 1e6:
            failures.append(f"{name}: drive.rep_rate_hz does not match deshpande2014_abstract max_rate_MHz")
        _anchor_ck(name, sources, "drive.rep_rate_hz", "deshpande2014_abstract")
        checks[0] += 1
        if n["nanowire"]["barrier_left_nm"] != geom["barrier_left_nm"] or \
           n["nanowire"]["barrier_right_nm"] != geom["barrier_right_nm"]:
            failures.append(f"{name}: nitride.nanowire barrier_*_nm does not match deshpande2013_geometry")
        _anchor_ck(name, sources, "nitride.nanowire.barrier_left_nm", "deshpande2013_geometry")
        _anchor_ck(name, sources, "nitride.nanowire.barrier_right_nm", "deshpande2013_geometry")
        checks[0] += 1
        if n["dot"]["height_nm"] != geom["disc_height_nm"]:
            failures.append(f"{name}: nitride.dot.height_nm does not match deshpande2013_geometry disc_height_nm")
        _anchor_ck(name, sources, "nitride.dot.height_nm", "deshpande2013_geometry")
        if "horizontal" in name:
            checks[0] += 1
            if n["nanowire"]["core_radius_nm"] != geom["optical_diameter_nm"] / 2.0:
                failures.append(f"{name}: horizontal core_radius_nm does not match half the "
                                 "deshpande2013_geometry optical_diameter_nm")
            _anchor_ck(name, sources, "nitride.nanowire.core_radius_nm", "deshpande2013_geometry")
            checks[0] += 1
            if n["wire_thermal"]["R_s_ohm"] != geom["resistance_GOhm"] * 1e9:
                failures.append(f"{name}: horizontal R_s_ohm does not match "
                                 "deshpande2013_geometry resistance_GOhm")
            _anchor_ck(name, sources, "nitride.wire_thermal.R_s_ohm", "deshpande2013_geometry")
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

    # ================================================== HIGH-1: contract schema value/range/tag
    # Parse docs/nitride_nanowire_contract.md's own Card schema tables and
    # assert every card leaf's value against the contract's own default (or
    # an explicitly enumerated sensitivity/sweep set), and that the card's
    # provenance tag is one the contract allows for that leaf. Opus review
    # of 58eb5a2 planted 8 single-leaf edits that the previous verifier
    # (leaf-name-set + provenance-presence only) never caught; this section
    # is the fix.
    _BLOCK_ACCESSOR = {
        "nitride.nanowire": lambda n, dr: n["nanowire"],
        "nitride.dot": lambda n, dr: n["dot"],
        "nitride.surface": lambda n, dr: n["surface"],
        "nitride.photonics": lambda n, dr: n["photonics"],
        "nitride.wire_thermal": lambda n, dr: n["wire_thermal"],
        "drive.diode": lambda n, dr: dr["diode"],
        "nitride.injector": lambda n, dr: n["injector"],
        "drive.set_params": lambda n, dr: dr.get("set_params", {}),
    }
    for name in CARD_NAMES:
        design = raws[name]["design"]
        n, dr = design["nitride"], design["drive"]
        fam = n["nanowire"]["family"]
        sources = design.get("provenance", {}).get("sources", {})
        for section, rows in SCHEMA_ROWS.items():
            block = _BLOCK_ACCESSOR[section](n, dr)
            for cells in rows:
                leaf, family_suffix, default_cell, range_cell, tag_cell = _parse_row_cells(section, cells)
                if leaf not in block:
                    continue
                value = block[leaf]
                checks[0] += 1
                try:
                    result = check_contract_leaf(section, leaf, family_suffix, default_cell, range_cell,
                                                  nit=n, value=value, family=fam)
                except SchemaParseError as exc:
                    failures.append(f"{name}: {section}.{leaf}: {exc}")
                    continue
                except Exception as exc:
                    failures.append(f"{name}: {section}.{leaf}: contract schema check raised: {exc}")
                    continue
                if result is None:
                    continue  # this table row is the other family's variant
                ok, note = result
                if not ok:
                    failures.append(f"{name}: {section}.{leaf}={value!r} violates the contract "
                                     f"Card schema (expected {note})")
                # provenance tag must be one the contract allows for this leaf.
                checks[0] += 1
                entry = sources.get(f"{section}.{leaf}")
                allowed_tags = parse_contract_tags(tag_cell)
                if entry is None or entry.get("tag") not in allowed_tags:
                    failures.append(f"{name}: {section}.{leaf} provenance tag "
                                     f"{entry.get('tag') if entry else None!r} not in contract tags {allowed_tags}")

    # The round's own explicit-default pins the contract states only in
    # prose, never in a Card-schema table row (see DIRECT_DEFAULTS above).
    for name in CARD_NAMES:
        design = raws[name]["design"]
        for path, expected in DIRECT_DEFAULTS.items():
            checks[0] += 1
            value, found = _get(design, path)
            if not found or value != expected:
                failures.append(f"{name}: {path}={value!r} must equal the round's pinned "
                                 f"default {expected!r}")

        sources = design.get("provenance", {}).get("sources", {})
        for path, expected_tag in DIRECT_DEFAULT_TAGS.items():
            checks[0] += 1
            entry = sources.get(path)
            tag = entry.get("tag") if entry else None
            if tag != expected_tag:
                failures.append(f"{name}: provenance.sources[{path!r}].tag={tag!r} must equal "
                                 f"the round's pinned tag {expected_tag!r}")

    # MEDIUM-2 fix: drive.duty is read directly by the nanowire device
    # (nitride_nanowire_device.py: duty = float(d.drive.duty)), not
    # recomputed from tau_pulse_ns/rep_rate_hz -- it must already equal
    # that product on every card, including a rep-rate sweep row.
    for name in CARD_NAMES:
        design = raws[name]["design"]
        checks[0] += 1
        expected_duty = design["drive"]["diode"]["tau_pulse_ns"] * 1e-9 * design["drive"]["rep_rate_hz"]
        if abs(design["drive"]["duty"] - expected_duty) > 1e-12:
            failures.append(f"{name}: drive.duty={design['drive']['duty']!r} does not equal "
                             f"tau_pulse_ns*1e-9*rep_rate_hz={expected_duty!r}")

    # acceptance-criterion-3 consistency checks (radii, d_i_nm, numeric types).
    for name in CARD_NAMES:
        design = raws[name]["design"]
        n = design["nitride"]
        nw, dot, diode = n["nanowire"], n["dot"], design["drive"]["diode"]
        fam = nw["family"]

        checks[0] += 1
        expected_outer = nw["core_radius_nm"] if nw["shell"] == "none" else nw["core_radius_nm"] + 3.0
        if nw["outer_radius_nm"] != expected_outer:
            failures.append(f"{name}: outer_radius_nm={nw['outer_radius_nm']!r} inconsistent with "
                             f"core_radius_nm/shell (expected {expected_outer!r})")

        checks[0] += 1
        if fam == "horizontal_as_built":
            radii_ok = dot["radius_nm"] == nw["core_radius_nm"]
        else:
            radii_ok = 0.0 < dot["radius_nm"] < nw["core_radius_nm"]
        if not radii_ok:
            failures.append(f"{name}: nitride.dot.radius_nm={dot['radius_nm']!r} inconsistent with "
                             f"family {fam!r}/core_radius_nm={nw['core_radius_nm']!r}")

        checks[0] += 1
        if not (0.0 < diode["conducting_radius_nm"] <= nw["core_radius_nm"]):
            failures.append(f"{name}: drive.diode.conducting_radius_nm={diode['conducting_radius_nm']!r} "
                             f"not in (0, core_radius_nm={nw['core_radius_nm']!r}]")

        set_params = design["drive"].get("set_params", {})
        if "radius_nm" in set_params:
            checks[0] += 1
            if set_params["radius_nm"] != dot["radius_nm"]:
                failures.append(f"{name}: drive.set_params.radius_nm={set_params['radius_nm']!r} "
                                 f"!= nitride.dot.radius_nm={dot['radius_nm']!r} (never core_radius_nm)")

        # d_i_nm = barrier_left_nm + height_nm + barrier_right_nm (contract
        # "Families and common card contract"): this round's cards all use
        # the V-tagged 2013-transferred geometry, so d_i_nm must reproduce
        # the ledger's own barrier/disc-height total.
        checks[0] += 1
        d_i_nm = nw["barrier_left_nm"] + dot["height_nm"] + nw["barrier_right_nm"]
        expected_d_i_nm = geom["barrier_left_nm"] + geom["disc_height_nm"] + geom["barrier_right_nm"]
        if d_i_nm != expected_d_i_nm:
            failures.append(f"{name}: d_i_nm (barrier_left+height+barrier_right)={d_i_nm!r} != "
                             f"the deshpande2013_geometry total {expected_d_i_nm!r}")

        # Numeric-typed leaves must not have drifted into a YAML string
        # (e.g. a quoted "12.5" instead of 12.5): every leaf this section's
        # Card-schema rows classify as numeric-default must still be a
        # Python int/float (and not bool, which is an int subclass) on disk.
        for section, rows in SCHEMA_ROWS.items():
            block = _BLOCK_ACCESSOR[section](n, design["drive"])
            for cells in rows:
                leaf, family_suffix, default_cell, range_cell, _tag = _parse_row_cells(section, cells)
                if leaf not in block or (family_suffix is not None and family_suffix != fam):
                    continue
                value = block[leaf]
                if value is None or isinstance(value, bool) or isinstance(value, str):
                    continue
                nums = _numbers_in(default_cell)
                looks_numeric = bool(nums) and not default_cell.strip().startswith(("`", "None", "required"))
                if looks_numeric:
                    checks[0] += 1
                    if not isinstance(value, (int, float)):
                        failures.append(f"{name}: {section}.{leaf}={value!r} is not numeric "
                                         f"(type {type(value).__name__})")

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
