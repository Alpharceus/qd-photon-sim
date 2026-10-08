"""fsim device designer (Dear PyGui): configure the device as a block diagram,
edit the fab stack, run, get graphs and numbers.

Multisim-style with one honest restriction: the physics defines ONE signal
chain (dot -> cavity -> filter -> detection) with drive and thermal attached,
so the canvas is that chain with configurable blocks -- not free-form wiring.

Three-layer rule: this file computes NO physics. It builds a DeviceDesign,
calls fsim_core.device.evaluate, and renders the results. Designs round-trip
to YAML in cards/ so every session is reproducible from its file. Phase D3
adds design comparison (Store A/B/C slots overlay the last RUN on the g2
plot + a delta table) and a "Report bundle" button that hands the current
run + stored slots to fsim_viz.report.designer_report -- imported lazily
inside that one callback so normal startup never pays for matplotlib.

Run:        python fsim_gui/designer.py
Smoke test: python fsim_gui/designer.py --frames 5
Selftest:   python fsim_gui/designer.py --selftest <outdir>
"""
import csv
import copy
import sys
from pathlib import Path

import dearpygui.dearpygui as dpg

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import design_meta, presets  # noqa: E402
from fsim_core.device import (DeviceDesign, evaluate, evaluate_envelope,  # noqa: E402
                              _legacy_density_cm2)
from fsim_core.loading import f8_g2_load, f8b_thin_fano, granularity_N  # noqa: E402
from fsim_theme import load_tokens  # noqa: E402
from fsim_theme.dpg import bind_theme  # noqa: E402

DESIGN_PATH = ROOT / "cards" / "staged-device-design.yaml"
LAYERS = []          # live fab-stack rows (list of dicts)
SUBSTRATE = {"name": "GaAs", "k300": 55.0, "alpha": 1.25}

# Preservation baseline (spec rt-edge-gui-preservation): the last design
# apply_design() was handed, kept verbatim. collect_design() deep-copies this
# and overwrites only the fields that have a live widget (WIDGET_TAG below),
# so any DeviceDesign field the GUI has no control for -- present today
# (ret.*, emission.*, drive.cw/diode/mechanism/..., aperture.compose,
# filter.track/track_material, ...) or added later -- round-trips untouched
# instead of silently reverting to its dataclass default.
_BASELINE_DESIGN = DeviceDesign()

# pr-pkg1-fix4 item 3: whether the aperture.density_cm2 widget has been
# touched by a real edit since it was last populated by apply_design().
# collect_design() used to decide "was this edited" by comparing the
# widget's CURRENT value against the constant unset_display = log10(7.0e8)
# -- but that display value is exactly what a None-density design's widget
# is ALSO populated with, so an edit that lands back on 7.0e8 (the same
# displayed default) was indistinguishable from no edit at all and was
# silently discarded (None written back instead of the user's 7.0e8). A
# value-only comparison can never fix this (the two cases produce the
# identical widget value); only tracking the edit event itself can, so
# this flag is set True by the widget's own callback (a real edit) and
# reset False by apply_design() (a re-population, not a user edit).
_DENSITY_EDITED = False

# Colours come from fsim_theme/tokens.json (dark), never hard-coded hexes.
# GREEN/AMBER/RED are the STATUS pair (verdicts / warnings only, always with a
# word); BLUE/PURPLE are categorical series slots. Provenance is NOT a hue:
# TAG_COLOR is one ink tone for every tag, the tag itself is carried by the
# chip text "[A]" (studio-02 / DIRECTION.md line-form grammar).
_TOK = load_tokens()
_C = _TOK["color"]["dark"]


def _rgb(hexv: str) -> tuple:
    h = hexv.lstrip("#")
    return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))


def _rgba_css(css: str) -> tuple:
    """'rgba(r,g,b,a)' -> DPG 0-255 RGBA."""
    r, g, b, a = (float(v) for v in css[css.index("(") + 1:css.index(")")].split(","))
    return (int(r), int(g), int(b), int(round(255 * a)))


SERIES = [_rgb(h) for h in _C["series"]]
INK1, INK2, MUTED, REF = _rgb(_C["ink-1"]), _rgb(_C["ink-2"]), _rgb(_C["muted"]), _rgb(_C["ref"])
REF_WASH = _rgba_css(_C["ref-wash"])
GREEN = _rgb(_C["status"]["pass"])
AMBER = _rgb(_C["status"]["warn"])
RED = _rgb(_C["status"]["fail"])
BLUE = SERIES[0]
PURPLE = SERIES[6]
TAG_COLOR = {"V": INK2, "DR": INK2, "E": INK2, "A": INK2}
FONTS = {"ui": None, "mono": None, "theme": None}   # filled by bind_theme() in main()

# ---- D3: comparison slots. LAST_RUN is the most recent successful RUN
# (point or envelope), same shape as a designer_report() entry minus
# "label": {design, name, curves, bands (or None), scalars, scalar_bands
# (or None), ranged (dict, possibly empty)}. SLOTS holds up to 3 stored
# copies of that shape, keyed "A"/"B"/"C".
LAST_RUN = None
SLOTS = {}
SLOT_LABELS = ("A", "B", "C")
SLOT_COLOR = {"A": SERIES[0], "B": SERIES[1], "C": SERIES[2]}   # comparison set, slots 1-3
DELTA_SCALARS = ["T_j_op", "dT_J", "eps_op", "rho_op", "g2_op", "brightness_per_pulse",
                "T_c", "F_eff", "N_w", "aperture_g2_penalty"]

# meta path -> live widget tag, for every parameter widget that carries a META
# entry. aperture.density_cm2 is edited as log10 in the widget (see transform
# in VIOLATIONS handling below).
WIDGET_TAG = {
    "dot.delta_xx": "dot.delta_xx", "dot.gamma_scale": "dot.gamma_scale",
    "dot.r_xx": "dot.r_xx",
    "drive.V": "drive.V", "drive.I_uA": "drive.I_uA", "drive.duty": "drive.duty",
    "drive.mu": "drive.mu", "drive.b_e": "drive.b_e", "drive.b_e_m": "drive.b_e_m",
    "drive.b_e_Eact": "drive.b_e_Eact",
    "drive.mode": "drive.mode", "drive.dg_inj": "drive.dg_inj",
    "drive.p_inj": "drive.p_inj", "drive.F_p": "drive.F_p",
    "drive.eta_capture": "drive.eta_capture", "drive.C_dep_pF": "drive.C_dep_pF",
    "thermal.mesa_diameter_um": "th.mesa", "thermal.T_hs": "th.T_hs",
    "cavity.enabled": "cav.enabled", "cavity.type": "cav.type",
    "cavity.kappa": "cav.kappa", "cavity.T_track": "cav.T_track",
    "cavity.E_X0": "cav.E_X0", "cavity.F_P": "cav.F_P", "cavity.G": "cav.G",
    "cavity.beta_sin": "cav.beta_sin",
    "filter.enabled": "fil.enabled", "filter.auto_w": "fil.auto_w",
    "filter.w": "fil.w", "filter.dx": "fil.dx",
    "aperture.density_cm2": "ap.log_density", "aperture.diameter_um": "ap.diam",
    "aperture.sigma_inh": "ap.sigma", "aperture.comp_brightness": "ap.r",
}

VIOLATIONS = set()   # meta paths currently outside their known-physical band


# ------------------------------------------------------------- tag bullets / validation

def _tag_bullet(path):
    """Colored provenance bullet + tooltip next to a parameter widget, keyed by
    design_meta.META. No-op for paths without a META entry (e.g. no widget)."""
    meta = design_meta.META.get(path)
    if meta is None:
        return
    color = TAG_COLOR.get(meta["tag"], INK2)
    b = dpg.add_text(f"[{meta['tag']}]", color=color)
    with dpg.tooltip(b):
        dpg.add_text(f"[{meta['tag']}] {meta['unit']}  |  band {meta['lo']:g}-{meta['hi']:g}"
                     f"  |  {meta['source']}")


def _validate(path, value):
    meta = design_meta.META.get(path)
    if meta is None or isinstance(value, str):
        return
    if meta["lo"] <= value <= meta["hi"]:
        VIOLATIONS.discard(path)
    else:
        VIOLATIONS.add(path)
    _refresh_warning()


def _refresh_warning():
    if VIOLATIONS and dpg.does_item_exist("param_warning"):
        names = ", ".join(sorted(VIOLATIONS))
        dpg.set_value("param_warning",
                      f"! {len(VIOLATIONS)} parameters outside known-physical bands: {names}")
    elif dpg.does_item_exist("param_warning"):
        dpg.set_value("param_warning", "")


# --------------------------------------------------------- envelope (D2) controls

def _toggle_range(path):
    on = dpg.get_value(f"rng.{path}")
    dpg.configure_item(f"rng.{path}.lo", show=on)
    dpg.configure_item(f"rng.{path}.hi", show=on)
    # the caption shortens to "range" beside the lo/hi fields (studio-08 R3:
    # the fields get the width; the checkbox tooltip keeps the full meaning)
    dpg.set_value(f"rng.{path}.cap", _range_caption(on))


def _range_caption(on):
    return "range" if on else "range (envelope)"


def _range_controls(path, extra=None):
    """Envelope toggle for a parameter with an ENV_DEFAULTS entry: checkbox
    '~' + lo/hi range inputs, shown only when checked. No-op for paths
    without an entry. Program decision: [A]-tagged params default RANGED on
    (see design_meta.default_ranged); apply_design() re-syncs the checked
    state and lo/hi values every time a design is loaded. `extra`, when
    given, also fires on the ~ toggle and on lo/hi edits (used by drive.mu
    to keep the F8 domain warning live -- see _update_drive_panel)."""
    bounds = design_meta.ENV_DEFAULTS.get(path)
    if bounds is None:
        return
    lo, hi = bounds
    on = path in design_meta.default_ranged(DeviceDesign())

    def _toggle_cb(s, v):
        _toggle_range(path)
        if extra:
            extra(s, v)

    with dpg.group(horizontal=True):
        dpg.add_checkbox(label="~", tag=f"rng.{path}", default_value=on,
                         callback=_toggle_cb)
        with dpg.tooltip(f"rng.{path}"):
            dpg.add_text("range (envelope): sweep this input between lo and hi")
        # 96 px + %.4g: lo/hi read in full next to the -/+ steppers
        # (1.5, 0.001, 2e+10); studio-08 R3
        dpg.add_input_float(tag=f"rng.{path}.lo", width=96, default_value=lo,
                            format="%.4g", show=on, callback=extra)
        dpg.add_input_float(tag=f"rng.{path}.hi", width=96, default_value=hi,
                            format="%.4g", show=on, callback=extra)
        dpg.add_text(_range_caption(on), tag=f"rng.{path}.cap",
                     color=(150, 150, 150))


def _reset_range_controls(d: DeviceDesign):
    """Re-derive the checked/unchecked state and lo/hi values for every
    envelope control from design_meta.default_ranged(d) -- called whenever a
    design is freshly loaded (apply_design), so a design's own [A] defaults
    always drive the initial envelope, not whatever the previous design left
    behind. Ranges are session state, not part of the saved DeviceDesign."""
    dr = design_meta.default_ranged(d)
    for path, bounds in design_meta.ENV_DEFAULTS.items():
        tag = f"rng.{path}"
        if not dpg.does_item_exist(tag):
            continue
        on = path in dr
        # use default_ranged's bounds where present (F8-narrowed mu floor for
        # quiet-pump designs -- Opus v1.1 finding), raw ENV_DEFAULTS otherwise
        lo, hi = dr.get(path, bounds)
        dpg.set_value(tag, on)
        dpg.set_value(f"{tag}.lo", lo)
        dpg.set_value(f"{tag}.hi", hi)
        dpg.configure_item(f"{tag}.lo", show=on)
        dpg.configure_item(f"{tag}.hi", show=on)
        dpg.set_value(f"{tag}.cap", _range_caption(on))


def _collect_ranged() -> dict:
    """Currently-checked envelope ranges, keyed by META path, as (lo, hi)."""
    out = {}
    for path in design_meta.ENV_DEFAULTS:
        tag = f"rng.{path}"
        if dpg.does_item_exist(tag) and dpg.get_value(tag):
            lo, hi = dpg.get_value(f"{tag}.lo"), dpg.get_value(f"{tag}.hi")
            out[path] = (min(lo, hi), max(lo, hi))
    return out


def _mk_cb(path, extra=None, transform=None):
    """Widget callback: validate against META, then run any pre-existing callback."""
    def cb(sender, value):
        _validate(path, value if transform is None else transform(value))
        if extra:
            extra(sender, value)
    return cb


def _mark_density_edited(sender=None, value=None):
    """pr-pkg1-fix4 item 3: the aperture.density_cm2 widget's own callback
    -- fires only on a real edit (never on apply_design()'s programmatic
    dpg.set_value(), which does not invoke callbacks), so collect_design()
    can tell "edited to the displayed default" apart from "never touched"
    without comparing widget values (see collect_design()'s docstring)."""
    global _DENSITY_EDITED
    _DENSITY_EDITED = True


def _update_drive_panel(sender=None, value=None):
    """Live refresh of the ELECTRICAL DRIVE node's v1.1 read-outs: the F9
    granularity info text, the F8 domain warning (amber, WARN-ONLY -- never
    blocks input, per the house rule), and the PL-mode note. Called from
    every v1.1 drive widget's callback and from apply_design(); a no-op
    before build_ui() has created the widgets. Every number here is computed
    by fsim_core.loading (granularity_N / f8b_thin_fano / f8_g2_load) --
    this function only reads widget values and formats text."""
    if not dpg.does_item_exist("drive.C_dep_pF"):
        return
    mode = dpg.get_value("drive.mode")
    F_p = dpg.get_value("drive.F_p")
    eta = dpg.get_value("drive.eta_capture")
    C_pF = dpg.get_value("drive.C_dep_pF")
    T_hs = dpg.get_value("th.T_hs") if dpg.does_item_exist("th.T_hs") else 77.0

    N = granularity_N(C_pF * 1e-12, T_hs)
    dpg.set_value("drive.N_info", f"N ~ {N:.1e} e-/quantum @ T_hs")

    dpg.set_value("drive.pl_note",
                  "PL: injection background + dg_inj disabled" if mode == "PL" else "")

    warn = ""
    if F_p != 1.0:
        floor = 1.0 - f8b_thin_fano(eta, F_p)
        if dpg.does_item_exist("rng.drive.mu") and dpg.get_value("rng.drive.mu"):
            mu_lo = dpg.get_value("rng.drive.mu.lo")
            mu_hi = dpg.get_value("rng.drive.mu.hi")
            mu_test = min(mu_lo, mu_hi)
        else:
            mu_test = dpg.get_value("drive.mu")
        if mu_test < floor:
            warn = (f"!  F8 domain: mu {mu_test:g} below floor {floor:.3f} "
                    f"(= 1 - F8b-thinned Fano) -- warn only, not blocked")
    dpg.set_value("drive.f8_warning", warn)


def revalidate_all():
    """Re-check every live widget against its META band (e.g. after a preset
    or a Load overwrites values without going through a widget callback)."""
    VIOLATIONS.clear()
    for path, tag in WIDGET_TAG.items():
        if not dpg.does_item_exist(tag):
            continue
        val = dpg.get_value(tag)
        if path == "aperture.density_cm2":
            val = 10.0 ** val
        _validate(path, val)
    _refresh_warning()


# ------------------------------------------------------------- design <-> widgets

def collect_design() -> DeviceDesign:
    """Deep-copy the preservation baseline (the last design apply_design()
    was handed -- see _BASELINE_DESIGN) and overwrite only the fields that
    have a live widget, via the same WIDGET_TAG map apply_design() writes
    from. Any field without a widget -- including nested dicts/lists on
    blocks the GUI has no editor for -- survives untouched; this is
    deliberately data-driven off WIDGET_TAG rather than a hand-written list
    of block/field names, so a field added later needs a WIDGET_TAG entry
    (or none, to keep it preserved) and nothing else here.

    aperture.density_cm2 (pr-pkg1-fix3 item 1, tracking fixed pr-pkg1-fix4
    item 3): the ONE WIDGET_TAG field that is float | None (device.py's
    ApertureBlock.density_cm2 -- None means "no explicit aperture density
    set", the legacy/unset state every consumer resolves its own way; see
    _legacy_density_cm2). The widget can only ever hold a float, so
    apply_design() displays the unset case as
    log10(_legacy_density_cm2(None)) (the same 7.0e8 legacy transport
    default other None-consumers fall back to). Whether to write None back
    is decided by the explicit _DENSITY_EDITED flag (set by the widget's
    own callback on a real edit, reset by apply_design() on repopulation)
    -- NOT by comparing the widget's current value against that displayed
    default, which a genuine edit landing back on 7.0e8 would be
    indistinguishable from (pr-pkg1-fix4 item 3: that value-only compare
    silently discarded such an edit and wrote None back).

    The widget itself is a Dear PyGui add_input_float, which stores its
    value as float32 (~7 significant figures), not the log10 value's own
    float64 precision -- so ``10.0 ** value`` on read-back does not
    reproduce whatever density was originally displayed bit-exactly, only
    to within float32 log-slider quantisation (rel ~1.05e-6 for the
    7.0e8 case verify/verify_designer_rt.py checks; see that file's
    WIDGET_FLOAT_PATHS/TOL for the general tolerance every widget-backed
    field is held to)."""
    d = copy.deepcopy(_BASELINE_DESIGN)
    d.name = dpg.get_value("design.name")
    for path, tag in WIDGET_TAG.items():
        if not dpg.does_item_exist(tag):
            continue
        value = dpg.get_value(tag)
        if path == "aperture.density_cm2":
            baseline_density = _BASELINE_DESIGN.aperture.density_cm2
            if baseline_density is None and not _DENSITY_EDITED:
                value = None
            else:
                value = 10.0 ** value
        block_name, field_name = path.split(".", 1)
        setattr(getattr(d, block_name), field_name, value)
    d.thermal.layers = [dict(L) for L in LAYERS]
    d.thermal.substrate = dict(SUBSTRATE)
    return d


def apply_design(d: DeviceDesign):
    """Push `d` onto every live widget via WIDGET_TAG (see collect_design)
    and record it verbatim as the new preservation baseline -- a deep copy,
    taken before anything below can mutate it, so a field with no widget
    (unknown to this GUI, today or in the future) still round-trips through
    the next collect_design().

    aperture.density_cm2 (pr-pkg1-fix3 item 1): the widget only holds a
    float (add_input_float), so an unset (None) density -- the legacy
    default every DeviceDesign() starts with -- is displayed at
    log10(_legacy_density_cm2(None)) = log10(7.0e8), the same transport
    legacy default. This used to be math.log10(value), which raised
    TypeError on None (verify_designer_rt.py dropped to 4/6, gate_v11_gui.py
    failed the --roundtrip-check subprocess -- a fresh DeviceDesign()'s
    density_cm2 is None). pr-pkg1-fix4 item 3: also clears _DENSITY_EDITED
    -- this is a repopulation from `d`, not a user edit, so collect_design()
    must not mistake the freshly-displayed default for one."""
    global LAYERS, SUBSTRATE, _BASELINE_DESIGN, _DENSITY_EDITED
    import math
    _BASELINE_DESIGN = copy.deepcopy(d)
    # pr-pkg6-stale-text item C4: reset unconditionally, BEFORE the
    # per-widget existence guard below -- this is a repopulation from `d`
    # regardless of whether the density widget happens to exist in the
    # current layout, so a stale True left over from a previous design
    # must not survive into this one just because the guard `continue`d
    # past the widget before the old reset (nested inside the loop body,
    # behind both the existence check and the path match) ever ran.
    _DENSITY_EDITED = False
    dpg.set_value("design.name", d.name)
    for path, tag in WIDGET_TAG.items():
        if not dpg.does_item_exist(tag):
            continue
        block_name, field_name = path.split(".", 1)
        value = getattr(getattr(d, block_name), field_name)
        if path == "aperture.density_cm2":
            value = math.log10(_legacy_density_cm2(value))
        dpg.set_value(tag, value)
    if dpg.does_item_exist("lemma1_note"):
        dpg.configure_item("lemma1_note", show=(d.cavity.type == "sin_waveguide"))
    LAYERS = [dict(L) for L in d.thermal.layers]
    SUBSTRATE = dict(d.thermal.substrate)
    rebuild_stack_table()
    draw_cross_section()
    revalidate_all()
    _reset_range_controls(d)
    _update_drive_panel()


def apply_presets():
    """Apply the four preset combos onto the live design, then push to widgets."""
    d = collect_design()
    presets.apply_dot_preset(d, dpg.get_value("preset.dot"))
    presets.apply_template_preset(d, dpg.get_value("preset.template"))
    presets.apply_cavity_preset(d, dpg.get_value("preset.cavity"))
    presets.apply_drive_preset(d, dpg.get_value("preset.drive"))
    presets.apply_injection_preset(d, dpg.get_value("preset.injection"))
    apply_design(d)
    dpg.set_value("status",
                  f"applied presets: dot={dpg.get_value('preset.dot')}  "
                  f"template={dpg.get_value('preset.template')}  "
                  f"cavity={dpg.get_value('preset.cavity')}  "
                  f"drive={dpg.get_value('preset.drive')}  "
                  f"injection={dpg.get_value('preset.injection')}")


# ------------------------------------------------------------------- fab stack UI

def rebuild_stack_table():
    dpg.delete_item("stack_table", children_only=True)
    dpg.add_table_column(label="layer", parent="stack_table")
    dpg.add_table_column(label="t (um)", parent="stack_table")
    dpg.add_table_column(label="k300", parent="stack_table")
    dpg.add_table_column(label="spread", parent="stack_table")
    dpg.add_table_column(label="", parent="stack_table")
    for i, L in enumerate(LAYERS):
        with dpg.table_row(parent="stack_table"):
            dpg.add_input_text(default_value=L["name"], width=90,
                               callback=_mk_edit(i, "name", str))
            dpg.add_input_float(default_value=L["t_um"], width=55, step=0,
                                callback=_mk_edit(i, "t_um", float))
            dpg.add_input_float(default_value=L["k300"], width=55, step=0,
                                callback=_mk_edit(i, "k300", float))
            dpg.add_checkbox(default_value=bool(L.get("spread", False)),
                             callback=_mk_edit(i, "spread", bool))
            dpg.add_button(label="x", callback=_mk_del(i), width=20)


def _mk_edit(i, key, cast):
    def cb(sender, val):
        LAYERS[i][key] = cast(val)
        draw_cross_section()
    return cb


def _mk_del(i):
    def cb():
        LAYERS.pop(i)
        rebuild_stack_table()
        draw_cross_section()
    return cb


def add_layer():
    LAYERS.append({"name": f"layer{len(LAYERS) + 1}", "t_um": 0.5, "k300": 30.0,
                   "alpha": 1.0, "spread": False})
    rebuild_stack_table()
    draw_cross_section()


def draw_cross_section():
    """Scale drawing of mesa + stack + substrate from the live values."""
    dpg.delete_item("xsec", children_only=True)
    W, H = 330, 250
    lateral_um = max(4.0, 2.0 * dpg.get_value("th.mesa"))
    px_per_um_x = W * 0.8 / lateral_um
    total_t = max(sum(L["t_um"] for L in LAYERS), 0.1)
    px_per_um_y = (H - 90) / total_t
    mesa_w = dpg.get_value("th.mesa") * px_per_um_x
    cx = W / 2

    y = 30
    dpg.draw_rectangle((cx - mesa_w / 2, y - 12), (cx + mesa_w / 2, y - 4),
                       fill=(200, 170, 60), parent="xsec")  # top contact
    dpg.draw_text((cx + mesa_w / 2 + 6, y - 14), "contact", size=12,
                  color=(180, 180, 180), parent="xsec")
    for i, L in enumerate(LAYERS):
        h = max(L["t_um"] * px_per_um_y, 8)
        wpx = mesa_w if not L.get("spread") else W * 0.8
        col = (70 + 35 * (i % 4), 95 + 20 * (i % 3), 160)
        dpg.draw_rectangle((cx - wpx / 2, y), (cx + wpx / 2, y + h), fill=col,
                           parent="xsec")
        dpg.draw_text((cx + wpx / 2 + 6, y + h / 2 - 7),
                      f"{L['name']} {L['t_um']:.2f} um", size=12,
                      color=(200, 200, 200), parent="xsec")
        if i == 0:  # QD plane sits under the first (top) layer
            dpg.draw_line((cx - wpx / 2, y + h - 2), (cx + wpx / 2, y + h - 2),
                          color=RED, thickness=2, parent="xsec")
            dpg.draw_text((cx - wpx / 2 - 62, y + h - 10), "QD layer", size=12,
                          color=RED, parent="xsec")
        y += h
    dpg.draw_rectangle((cx - W * 0.42, y), (cx + W * 0.42, y + 34),
                       fill=(60, 60, 68), parent="xsec")
    dpg.draw_text((cx - 40, y + 10), f"{SUBSTRATE['name']} substrate", size=12,
                  color=(200, 200, 200), parent="xsec")
    dpg.draw_text((10, y + 44),
                  f"mesa {dpg.get_value('th.mesa'):.2f} um  @  "
                  f"T_hs {dpg.get_value('th.T_hs'):.0f} K", size=12,
                  color=AMBER, parent="xsec")


def _f8_result_lines(d: DeviceDesign) -> list:
    """F8 results-panel lines for the given (point or mid-envelope) design:
    the g2_load factor (fsim_core.loading.f8_g2_load at the thinned F_eff --
    same F_eff device.evaluate/integrator.g2_of_T actually feed to f8_g2),
    "--" when mu/F_p sit outside the F8 domain, and a quiet-pump indicator
    when F_p < 1. Empty list when F_p == 1.0 (Poisson pump: F8 doesn't
    apply, nothing new to report)."""
    lines = []
    if d.drive.F_p != 1.0:
        F_eff = f8b_thin_fano(d.drive.eta_capture, d.drive.F_p)
        if d.drive.mu > 0:
            try:
                val = float(f8_g2_load(d.drive.mu, F_eff))
                lines.append(f"g2_load factor: {val:.3f} (F8)")
            except ValueError:
                lines.append("g2_load factor: -- (F8)")
        else:
            lines.append("g2_load factor: -- (F8)")
        if d.drive.F_p < 1.0:
            lines.append("[F8] quiet-pump")
    return lines


def _fmt_or_na(s: dict, key: str, fmt: str = "{:.4g}", unit: str = "") -> str:
    """One evaluate() scalar, formatted, or 'n/a' when the key is absent or
    NaN -- lets the results panel work both before and after a given
    evaluator key lands (spec rt-edge-gui-preservation). Every value here is
    a straight pass-through of an evaluate()/evaluate_envelope() scalar,
    never computed in this file (three-layer rule)."""
    v = s.get(key)
    if v is None or (isinstance(v, float) and v != v):
        return "n/a"
    suffix = f" {unit}" if unit else ""
    return fmt.format(v) + suffix


def _derived_result_lines(s: dict) -> list:
    """Read-only labelled lines for the transport/CW-derived scalars
    (integration-a/-b): junction voltage, applied voltage, built-in voltage,
    junction power, injection/capture efficiency, resolved loading mu, the
    background rate (per collected X photon), and the CW g2(0) pair --
    intrinsic g2_cw0 vs IRF-convolved g2_cw0_raw, kept distinct from the
    pulsed-intrinsic g2 headline number above. Labels only, never widgets:
    nothing here is ever written back into a design."""
    return [
        "",
        f"junction temperature (transport): {_fmt_or_na(s, 'T_j_transport', unit='K')}",
        f"junction voltage V_j: {_fmt_or_na(s, 'V_j_op', unit='V')}",
        f"applied voltage V_applied: {_fmt_or_na(s, 'V_applied', unit='V')}",
        f"built-in voltage V_bi: {_fmt_or_na(s, 'V_bi', unit='V')}",
        f"junction power P_junction: {_fmt_or_na(s, 'P_junction_W', unit='W')}",
        f"injection efficiency eta_inj: {_fmt_or_na(s, 'eta_inj')}",
        f"capture efficiency eta_capture: {_fmt_or_na(s, 'eta_capture_resolved')}",
        f"loading mu (resolved): {_fmt_or_na(s, 'mu_resolved')}",
        f"background per collected X photon: {_fmt_or_na(s, 'b_e_resolved')}",
        f"g2_cw0 (intrinsic CW): {_fmt_or_na(s, 'g2_cw0')}",
        f"g2_cw0_raw (IRF-convolved CW): {_fmt_or_na(s, 'g2_cw0_raw')}",
    ]


# ------------------------------------------- static (non-headline) labelling
#
# Spec audit-edge-cards-label (user decision Q2, 2026-09-23): the RT edge
# cards (cards/edge-inp-gainp-design.yaml, cards/edge-inp-gaasp-design.yaml)
# ship drive.finite_pulse=false, so a card-level evaluate() uses the static
# per-pulse loading that the post-peer-review verdict rejected as the
# headline model (headline = drive.finite_pulse=true with the tau_cap density
# law off, scripts/run_rt_edge.py HEADLINE_MODEL). The cards are NOT switched;
# every card-level g2/brightness this panel emits for such a design is labelled
# "static (non-headline)" and the headline finite-pulse g2 is named next to
# it. Labelling only: no emitted static number changes.
STATIC_LABEL = "static (non-headline)"


def _is_static_edge(d: DeviceDesign) -> bool:
    """True for an RT edge-emitter design evaluated with the static per-pulse
    loading (emission.type='edge', drive.finite_pulse false) -- the case
    whose g2/brightness is not the headline model's."""
    return (getattr(d.emission, "type", "none") == "edge"
            and not bool(getattr(d.drive, "finite_pulse", False)))


def _headline_g2(d: DeviceDesign) -> float:
    """Headline-model g2_op at this design's own thermal.T_hs: a straight
    pass-through of a fresh evaluate() on a copy with the headline switches
    (drive.finite_pulse=True, ret.tau_cap_scales_with_density=False; the CW
    diagnostic is skipped, it does not enter g2_op). Never changes `d` or
    the static numbers; NaN if the headline evaluation cannot run."""
    dh = copy.deepcopy(d)
    dh.drive.finite_pulse = True
    dh.ret.tau_cap_scales_with_density = False
    dh.drive.cw = False
    try:
        return float(evaluate(dh, T_grid=[dh.thermal.T_hs])["scalars"]["g2_op"])
    except Exception:  # noqa: BLE001 -- a label must never break the panel
        return float("nan")


def _static_label_lines(d: DeviceDesign, headline_g2: float) -> list:
    """The label block appended below the results for a static edge design:
    names the static loading, the headline model and its g2 at this point."""
    hv = "n/a" if headline_g2 != headline_g2 else f"{headline_g2:.4f}"
    return [
        "",
        f"NOTE: g2(0)/brightness above are {STATIC_LABEL}: drive.finite_pulse=false "
        "(static per-pulse loading).",
        "headline model = drive.finite_pulse=true (scripts/run_rt_edge.py HEADLINE_MODEL); "
        f"headline finite-pulse g2_op at T_hs {d.thermal.T_hs:.0f} K = {hv}",
    ]


def _f_eff_label(d: DeviceDesign) -> str:
    """Audit H5 (2026-09-23): the F_eff line is the SINGLE-MODE F_P input
    (device.py: F_P kappa/(kappa+Gamma)), never a planar-DBR total-rate
    enhancement (dbr.planar_total_rate ~1.0 for a planar lambda cavity;
    Bjork et al., PRA 44, 669 (1991)). Label only; no number changes."""
    pl = presets.CAVITY_PRESETS.get("planar-lambda", {})
    is_planar_preset = (bool(d.cavity.enabled) and d.cavity.type == pl.get("type")
                        and float(d.cavity.F_P) == float(pl.get("F_P", float("nan")))
                        and float(d.cavity.kappa) == float(pl.get("kappa", float("nan")))
                        and float(d.cavity.G) == float(pl.get("G", float("nan"))))
    if is_planar_preset:
        return "F_eff (mode-only F_P [A]; planar DBR total rate ~1.0)"
    return "F_eff (single-mode F_P input, not a planar total rate)"


# ------------------------------------------------------------------------- run

def run_device():
    """RUN dispatcher: any checked '~' range switches the whole run to
    envelope mode (bands + tornado); otherwise the original single-point
    path runs unchanged."""
    d = collect_design()
    ranged = _collect_ranged()
    if ranged:
        _run_envelope(d, ranged)
    else:
        _run_point(d)
    _refresh_open_expand_windows()


def _run_point(d: DeviceDesign):
    global LAST_RUN
    res = evaluate(d)
    c, s = res["curves"], res["scalars"]
    LAST_RUN = {"design": d, "name": d.name, "curves": c, "bands": None,
                "scalars": s, "scalar_bands": None, "ranged": {}}
    Ts = list(map(float, c["T_hs"]))
    dpg.set_value("s_g2", [Ts, list(map(float, c["g2"]))])
    dpg.set_value("s_eps", [Ts, list(map(float, c["eps"]))])
    dpg.set_value("s_rho2", [Ts, list(map(float, c["rho2"]))])
    dpg.set_value("s_half", [[Ts[0], Ts[-1]], [0.5, 0.5]])
    tc = s["T_c"]
    _set_tc_marks(tc, tc, Ts, point=True)
    dpg.set_value("s_tj", [Ts, [tj - t for tj, t in zip(map(float, c["Tj"]), Ts)]])
    dpg.set_value("s_gam", [Ts, list(map(float, c["gamma"]))])
    dpg.set_value("s_g2_band", [[], [], []])
    dpg.set_value("s_rho2_band", [[], [], []])
    dpg.set_value("s_tj_band", [[], [], []])
    dpg.set_value("envelope_header", "")
    _fit_result_axes()

    dpg.set_value("warn_runaway",
                  "!!!  THERMAL RUNAWAY -- no operating point  !!!" if s["runaway"] else "")
    dpg.set_value("warn_eps",
                  "!!!  CAVITY SELECTS XX (eps > 1) -- retune tracking/filter  !!!"
                  if s["eps_op"] > 1 else "")

    static_edge = _is_static_edge(d)
    static_tag = f"   [{STATIC_LABEL}]" if static_edge else ""
    lines = [
        f"tag chain {s['tag_chain']}  (unmeasured inputs -> conditional numbers;"
        f" envelopes: run_phase3)",
        "",
        f"T_j at operating point   {s['T_j_op']:.1f} K   (dT_J = {s['dT_J']:.2f} K)"
        + ("   ** THERMAL RUNAWAY **" if s["runaway"] else ""),
        f"Gamma(T_j)               {s['gamma_op']:.2f} meV",
        f"eps = t_XX/t_X           {s['eps_op']:.4f}",
        f"rho (signal purity)      {s['rho_op']:.3f}",
        f"g2(0) at operating point {s['g2_op']:.3f}" + static_tag,
        f"brightness/pulse (t_X)   {s['brightness_per_pulse']:.3f}" + static_tag,
        f"master ceiling T_c       "
        + (f"{s['T_c']:.0f} K" if s["T_c"] == s["T_c"] else "not crossed in range"),
        f"{_f_eff_label(d)}  "
        + (f"{s['F_eff']:.1f}" if s["F_eff"] == s["F_eff"] else "-- (cavity off)"),
        f"aperture: N_w = {s['N_w']:.2f}  ->  F5 g2 penalty {s['aperture_g2_penalty']:.3f}",
    ] + _derived_result_lines(s) + _f8_result_lines(d)
    if static_edge:
        lines += _static_label_lines(d, _headline_g2(d))
    dpg.set_value("results_text", "\n".join(lines))
    draw_cross_section()
    _refresh_delta_table()


def _interval(lo, hi, fmt="{:.3f}", unit=""):
    if lo != lo or hi != hi:  # NaN
        return "not reached in swept range"
    suffix = f" {unit}" if unit else ""
    return f"{fmt.format(lo)} - {fmt.format(hi)}{suffix}"


def _mid_design(d: DeviceDesign, ranged: dict) -> DeviceDesign:
    """Design with every ranged path pinned to its (lo+hi)/2 -- the same point
    evaluate_envelope's 'mid' curve already uses, so a point-scalars readout
    (T_j_op/brightness/F_eff/N_w/aperture penalty; not covered by
    scalar_bands) exists for envelope runs too, e.g. for slot storage and the
    report bundle. Bookkeeping only, not physics: evaluate() still does the
    actual chain."""
    dm = copy.deepcopy(d)
    for path, (lo, hi) in ranged.items():
        block_name, field_name = path.split(".", 1)
        setattr(getattr(dm, block_name), field_name, 0.5 * (lo + hi))
    return dm


def _run_envelope(d: DeviceDesign, ranged: dict):
    global LAST_RUN
    env = evaluate_envelope(d, ranged)
    mid, bands, sb = env["mid"], env["bands"], env["scalar_bands"]
    dm = _mid_design(d, ranged)
    mid_scalars = evaluate(dm, T_grid=[d.thermal.T_hs])["scalars"]
    LAST_RUN = {"design": d, "name": d.name, "curves": mid, "bands": bands,
                "scalars": mid_scalars, "scalar_bands": sb, "ranged": dict(ranged)}
    Ts = list(map(float, mid["T_hs"]))

    dpg.set_value("s_g2", [Ts, list(map(float, mid["g2"]))])
    dpg.set_value("s_eps", [Ts, list(map(float, mid["eps"]))])
    dpg.set_value("s_rho2", [Ts, list(map(float, mid["rho2"]))])
    dpg.set_value("s_half", [[Ts[0], Ts[-1]], [0.5, 0.5]])
    # a single marker can't honestly show a band: draw the T_c INTERVAL from
    # scalar_bands (or the "not crossed in range" note) instead
    _set_tc_marks(sb["T_c"][0], sb["T_c"][1], Ts, point=False)
    dTj_mid = [tj - t for tj, t in zip(map(float, mid["Tj"]), Ts)]
    dpg.set_value("s_tj", [Ts, dTj_mid])
    dpg.set_value("s_gam", [Ts, list(map(float, mid["gamma"]))])

    g2_lo, g2_hi = (list(map(float, a)) for a in bands["g2"])
    rho2_lo, rho2_hi = (list(map(float, a)) for a in bands["rho2"])
    tj_lo, tj_hi = bands["Tj"]
    dTj_lo = [tj - t for tj, t in zip(map(float, tj_lo), Ts)]
    dTj_hi = [tj - t for tj, t in zip(map(float, tj_hi), Ts)]
    dpg.set_value("s_g2_band", [Ts, g2_lo, g2_hi])
    dpg.set_value("s_rho2_band", [Ts, rho2_lo, rho2_hi])
    dpg.set_value("s_tj_band", [Ts, dTj_lo, dTj_hi])
    _fit_result_axes()

    dpg.set_value("warn_runaway", "")
    dpg.set_value("warn_eps", "")

    dpg.set_value("envelope_header",
                  f"ENVELOPE over {env['n_samples']} samples (honest mode):\n"
                  "[A] inputs swept; uncheck ~ for point designs")

    baseline = "T_c" if sb["T_c"][0] == sb["T_c"][0] else "g2_op"
    b_unit = "K" if baseline == "T_c" else ""
    ranked = sorted(env["tornado"].items(), key=lambda kv: kv[1], reverse=True)
    tornado_lines = []
    for path, red in ranked:
        if red != red:
            tornado_lines.append(f"    measure {path} first: ({baseline} band not finite)")
        else:
            tornado_lines.append(
                f"    measure {path} first: narrows {baseline} band by {red:.3g}"
                + (f" {b_unit}" if b_unit else ""))

    lines = [
        f"tag chain [A]  (envelope over {env['n_samples']} [A]-swept samples;"
        f" mid line = all ranges at their midpoint)",
        "",
        f"g2(0) at op:               {_interval(*sb['g2_op'])}"
        + (f"   [{STATIC_LABEL}]" if _is_static_edge(d) else ""),
        f"eps = t_XX/t_X at op:      {_interval(*sb['eps_op'])}",
        f"rho (signal purity) at op: {_interval(*sb['rho_op'])}",
        f"dT_J at op:                {_interval(*sb['dT_J'], fmt='{:.2f}', unit='K')}",
        f"master ceiling T_c:        {_interval(*sb['T_c'], fmt='{:.0f}', unit='K')}",
        "",
        "sensitivity (measurement priority, highest first):",
    ] + tornado_lines + _derived_result_lines(mid_scalars)
    f8_lines = _f8_result_lines(dm)
    if f8_lines:
        lines += [""] + f8_lines
    if _is_static_edge(d):
        # headline g2 at the mid design (all ranges at their midpoint), the
        # same point mid_scalars above reports.
        lines += _static_label_lines(dm, _headline_g2(dm))
    dpg.set_value("results_text", "\n".join(lines))
    draw_cross_section()
    _refresh_delta_table()


# --------------------------------------------------------- D3: comparison slots

def _fmt_delta(entry: dict, name: str) -> str:
    """Point value, or 'lo..hi' when `entry` carries a scalar_bands interval
    for this name -- same convention as fsim_viz.report._fmt_scalar, kept as
    a small local copy so the delta table (refreshed after every RUN) never
    has to import the viz package."""
    sb = entry.get("scalar_bands")
    if sb and name in sb:
        lo, hi = sb[name]
        return "nan" if (lo != lo or hi != hi) else f"{lo:.3g}..{hi:.3g}"
    v = entry["scalars"].get(name, float("nan"))
    return "nan" if (isinstance(v, float) and v != v) else f"{v:.3g}"


def _refresh_delta_table():
    if not dpg.does_item_exist("delta_text"):
        return
    if LAST_RUN is None or not SLOTS:
        dpg.set_value("delta_text", "")
        return
    cols = [("current", LAST_RUN)] + [(lbl, SLOTS[lbl]) for lbl in SLOT_LABELS if lbl in SLOTS]
    lines = ["", "comparison (current vs stored slots):",
             f"{'scalar':<22}" + "".join(f"{c:>17}" for c, _ in cols)]
    for name in DELTA_SCALARS:
        row = "".join(f"{_fmt_delta(entry, name):>17}" for _, entry in cols)
        lines.append(f"{name:<22}{row}")
    dpg.set_value("delta_text", "\n".join(lines))


def _slot_tags(label: str) -> tuple:
    return f"s_g2_{label}", f"s_g2_band_{label}"


def _build_slot_themes():
    """Bind each slot's line/band series to its fixed SLOT_COLOR (A amber, B
    blue, C purple) so overlays stay visually distinct from the auto-cycled
    current-run series. Called once from build_ui()."""
    for lbl in SLOT_LABELS:
        line_tag, band_tag = _slot_tags(lbl)
        r, g, b = SLOT_COLOR[lbl]
        with dpg.theme() as lt:
            with dpg.theme_component(dpg.mvLineSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Line, (r, g, b, 255),
                                    category=dpg.mvThemeCat_Plots)
        dpg.bind_item_theme(line_tag, lt)
        with dpg.theme() as bt:
            with dpg.theme_component(dpg.mvShadeSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Fill, (r, g, b, 70),
                                    category=dpg.mvThemeCat_Plots)
        dpg.bind_item_theme(band_tag, bt)


def _series_theme(item, rgb, *, alpha=1.0, weight=None, shade=False):
    """Bind one series' colour (line, or shade fill) from the token palette."""
    if not dpg.does_item_exist(item) or dpg.get_item_theme(item):
        return  # missing, or already themed (expand copies refresh every RUN)
    r, g, b = rgb[:3]
    a = int(round(255 * alpha)) if len(rgb) == 3 else rgb[3]
    with dpg.theme() as th:
        with dpg.theme_component(dpg.mvShadeSeries if shade else dpg.mvLineSeries):
            dpg.add_theme_color(dpg.mvPlotCol_Fill if shade else dpg.mvPlotCol_Line,
                                (r, g, b, a), category=dpg.mvThemeCat_Plots)
            if weight is not None and not shade:
                dpg.add_theme_style(dpg.mvPlotStyleVar_LineWeight, float(weight),
                                    category=dpg.mvThemeCat_Plots)
    dpg.bind_item_theme(item, th)


def _build_series_themes(prefix=""):
    """02-charts encodings for the live-run series: g2 slot 1, eps slot 2,
    rho^2 slot 3; the 0.5 ceiling and T_c are `ref` ink at 1 px (references,
    not series); bands at 12% of their series hue; the T_c interval is the
    neutral ref-wash. `prefix` themes the expand-window copies ("x_")."""
    band = _TOK["mark"]["band_opacity"]
    line_w = _TOK["mark"]["line"]
    _series_theme(prefix + "s_g2", SERIES[0], weight=line_w)
    _series_theme(prefix + "s_eps", SERIES[1], weight=line_w)
    _series_theme(prefix + "s_rho2", SERIES[2], weight=line_w)
    _series_theme(prefix + "s_half", REF, weight=_TOK["mark"]["ref"])
    _series_theme(prefix + "s_tc", REF, weight=_TOK["mark"]["ref"])
    _series_theme(prefix + "s_g2_band", SERIES[0], alpha=band, shade=True)
    _series_theme(prefix + "s_rho2_band", SERIES[2], alpha=band, shade=True)
    _series_theme(prefix + "s_tc_band", REF_WASH, shade=True)
    _series_theme(prefix + "s_tj", SERIES[0], weight=line_w)
    _series_theme(prefix + "s_tj_band", SERIES[0], alpha=band, shade=True)
    _series_theme(prefix + "s_gam", SERIES[0], weight=line_w)


def _g2_axis_limits(yaxis_tag="yax1"):
    """g2 axis fixed to 0..1 (the 0.5 ceiling always visible); only widened
    when a plotted value really exceeds 1 (e.g. eps > 1), never autoscaled
    down to a sliver."""
    top = 1.0
    for tag in ("s_g2", "s_eps", "s_rho2", "s_g2_band", "s_rho2_band"):
        if dpg.does_item_exist(tag):
            for arr in dpg.get_value(tag)[1:3]:
                vals = [v for v in arr if v == v]
                if vals:
                    top = max(top, max(vals))
    if dpg.does_item_exist(yaxis_tag):
        dpg.set_axis_limits(yaxis_tag, 0.0, top * 1.02 if top > 1.0 else 1.0)


def _set_tc_marks(tc_lo, tc_hi, Ts, point: bool):
    """T_c reference on the g2 plot. Point run: 1px ref rule labelled with
    the value. Envelope run: ref-wash span between scalar_bands['T_c'] lo/hi,
    or the 'not crossed in range' note (mirrors _interval())."""
    finite = tc_lo == tc_lo and tc_hi == tc_hi
    dpg.set_value("s_tc_band", [[], [], []])
    dpg.set_value("s_tc", [[], []])
    dpg.configure_item("tc_note", show=False)
    if point:
        if finite:
            dpg.set_value("s_tc", [[tc_lo, tc_lo], [0.0, 1.0]])
            dpg.configure_item("s_tc", label=f"T_c = {tc_lo:.0f} K")
        else:
            dpg.configure_item("s_tc", label="T_c")
    elif finite:
        dpg.set_value("s_tc_band", [[tc_lo, tc_hi], [0.0, 0.0], [1.0, 1.0]])
        dpg.configure_item("s_tc_band", label=f"T_c {tc_lo:.0f} - {tc_hi:.0f} K")
        dpg.configure_item("s_tc", label="T_c")
    if not finite and Ts:
        dpg.set_value("tc_note", (float(Ts[0]) + 0.02 * (float(Ts[-1]) - float(Ts[0])), 0.95))
        dpg.configure_item("tc_note", label="T_c: not crossed in range", show=True)


def _fit_result_axes():
    dpg.fit_axis_data("xax1")
    _g2_axis_limits("yax1")
    for ax in ("xax2", "yax2", "xax3", "yax3"):
        if dpg.does_item_exist(ax):
            dpg.fit_axis_data(ax)


def _refresh_slot_overlay():
    """Push every stored slot's g2 curve (+ band, if it has one) onto the
    pre-created hidden series on the main plot; hide a slot's series when it
    has no stored design."""
    for lbl in SLOT_LABELS:
        line_tag, band_tag = _slot_tags(lbl)
        if not dpg.does_item_exist(line_tag):
            continue
        entry = SLOTS.get(lbl)
        if entry is None:
            dpg.configure_item(line_tag, show=False)
            dpg.configure_item(band_tag, show=False)
            continue
        c = entry["curves"]
        Ts = list(map(float, c["T_hs"]))
        dpg.set_value(line_tag, [Ts, list(map(float, c["g2"]))])
        dpg.configure_item(line_tag, show=True, label=f"{lbl}: {entry['name']}")
        bands = entry.get("bands")
        if bands and "g2" in bands:
            lo, hi = (list(map(float, a)) for a in bands["g2"])
            dpg.set_value(band_tag, [Ts, lo, hi])
            dpg.configure_item(band_tag, show=True, label=f"{lbl} band")
        else:
            dpg.set_value(band_tag, [[], [], []])
            dpg.configure_item(band_tag, show=False)
    dpg.fit_axis_data("xax1")
    _g2_axis_limits("yax1")


def store_slot(label: str):
    """'Store <label>' button: capture LAST_RUN into a comparison slot. No-op
    (with a status warning) if nothing has been run yet."""
    if LAST_RUN is None:
        dpg.set_value("status", f"nothing run yet -- press RUN before storing slot {label}")
        return
    SLOTS[label] = LAST_RUN
    _refresh_slot_overlay()
    _refresh_delta_table()
    dpg.set_value("status", f"stored slot {label}: {LAST_RUN['name']}")


def clear_slots():
    SLOTS.clear()
    _refresh_slot_overlay()
    _refresh_delta_table()
    dpg.set_value("status", "cleared comparison slots")


# ------------------------------------------------------ expand-to-window plots

# Full-detail popouts for the two results plots. Each holds a lazily-created
# copy of the inline plot's series (own tags, prefixed "x_") that is
# refreshed from the live source series on every [ expand ] click and after
# every RUN (see _refresh_open_expand_windows, called from run_device()).
EXPAND_SPECS = {
    "main": {
        "window_tag": "expand_win_main",
        "title": "g2 / fractions vs heatsink T -- full detail",
        "xaxis_tag": "x_xax1", "yaxis_tag": "x_yax1",
        "xlabel": "heatsink T (K)", "ylabel": "g2 / fractions",
        "lines": ["s_g2", "s_eps", "s_rho2", "s_half", "s_tc"],
        "shades": ["s_tc_band", "s_g2_band", "s_rho2_band"],
        # D3 slot overlays: only copied once populated (see _copy_one_series)
        "slot_lines": [_slot_tags(lbl)[0] for lbl in SLOT_LABELS],
        "slot_shades": [_slot_tags(lbl)[1] for lbl in SLOT_LABELS],
    },
    "secondary": {
        "window_tag": "expand_win_secondary",
        "title": "dT_J / Gamma vs heatsink T -- full detail",
        "xaxis_tag": "x_xax2", "yaxis_tag": "x_yax2",
        "xlabel": "heatsink T (K)", "ylabel": "dT_J (K)",
        "lines": ["s_tj"],
        "shades": ["s_tj_band"],
        # second stacked panel (own units, linked x): Gamma in meV
        "panel2": {"xaxis_tag": "x_xax3", "yaxis_tag": "x_yax3", "ylabel": "Gamma (meV)",
                   "lines": ["s_gam"]},
        "slot_lines": [],
        "slot_shades": [],
    },
}


def _expand_tag(src_tag: str) -> str:
    return f"x_{src_tag}"


def _copy_one_series(src_tag: str, yaxis_tag: str, is_shade: bool, skip_if_empty: bool):
    """Push one source series' CURRENT data into its expanded-window twin,
    creating the twin lazily on first copy (distinct tag, see _expand_tag).
    `skip_if_empty` is used for the D3 slot overlays -- an unstored slot has
    no data and shouldn't clutter the expanded legend until it is populated
    (mirrors the inline plot's show=False convention, see
    _refresh_slot_overlay); once created, an emptied slot is hidden rather
    than deleted so it reappears cleanly if the slot is re-stored."""
    dst_tag = _expand_tag(src_tag)
    data = dpg.get_value(src_tag)
    has_data = len(data[0]) > 0
    if skip_if_empty and not has_data:
        if dpg.does_item_exist(dst_tag):
            dpg.configure_item(dst_tag, show=False)
        return
    label = dpg.get_item_label(src_tag)
    if not dpg.does_item_exist(dst_tag):
        if is_shade:
            dpg.add_shade_series(data[0], data[1], y2=data[2], label=label,
                                 tag=dst_tag, parent=yaxis_tag)
        else:
            dpg.add_line_series(data[0], data[1], label=label, tag=dst_tag,
                                parent=yaxis_tag)
    else:
        dpg.set_value(dst_tag, data)
        dpg.configure_item(dst_tag, label=label, show=True)


def _refresh_expand_window(kind: str):
    """Refresh every series in the expanded window for `kind` ("main" or
    "secondary") from its live source series. No-op if that window has never
    been opened."""
    spec = EXPAND_SPECS[kind]
    if not dpg.does_item_exist(spec["window_tag"]):
        return
    yax = spec["yaxis_tag"]
    for tag in spec["shades"]:
        _copy_one_series(tag, yax, is_shade=True, skip_if_empty=False)
    for tag in spec["lines"]:
        _copy_one_series(tag, yax, is_shade=False, skip_if_empty=False)
    for tag in spec["slot_shades"]:
        _copy_one_series(tag, yax, is_shade=True, skip_if_empty=True)
    for tag in spec["slot_lines"]:
        _copy_one_series(tag, yax, is_shade=False, skip_if_empty=True)
    p2 = spec.get("panel2")
    if p2:
        for tag in p2["lines"]:
            _copy_one_series(tag, p2["yaxis_tag"], is_shade=False, skip_if_empty=False)
        dpg.fit_axis_data(p2["xaxis_tag"])
        dpg.fit_axis_data(p2["yaxis_tag"])
    _build_series_themes(prefix="x_")
    dpg.fit_axis_data(spec["xaxis_tag"])
    if kind == "main":
        _g2_axis_limits(yax)
    else:
        dpg.fit_axis_data(yax)


def _open_expand_window(kind: str):
    """[ expand ] button callback: open a ~1200x700 resizable/closable window
    with a full-size copy of the given plot, or -- if it's already open --
    focus it and refresh its data. Native dpg zoom/pan (scroll, right-drag
    box-zoom, double-click autofit) is left at defaults on this plot too."""
    spec = EXPAND_SPECS[kind]
    win_tag = spec["window_tag"]
    if dpg.does_item_exist(win_tag):
        dpg.configure_item(win_tag, show=True)
        dpg.focus_item(win_tag)
    else:
        with dpg.window(tag=win_tag, label=spec["title"], width=1200, height=700,
                        pos=(120, 60)):
            p2 = spec.get("panel2")
            if p2 is None:
                with dpg.plot(height=-1, width=-1):
                    dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                    dpg.add_plot_axis(dpg.mvXAxis, label=spec["xlabel"],
                                      tag=spec["xaxis_tag"])
                    dpg.add_plot_axis(dpg.mvYAxis, label=spec["ylabel"],
                                      tag=spec["yaxis_tag"])
            else:  # two stacked panels, linked x, each with its own units
                with dpg.subplots(2, 1, height=-1, width=-1, link_all_x=True):
                    with dpg.plot():
                        dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                        dpg.add_plot_axis(dpg.mvXAxis, label="", tag=spec["xaxis_tag"],
                                          no_tick_labels=True)
                        dpg.add_plot_axis(dpg.mvYAxis, label=spec["ylabel"],
                                          tag=spec["yaxis_tag"])
                    with dpg.plot():
                        dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                        dpg.add_plot_axis(dpg.mvXAxis, label=spec["xlabel"],
                                          tag=p2["xaxis_tag"])
                        dpg.add_plot_axis(dpg.mvYAxis, label=p2["ylabel"],
                                          tag=p2["yaxis_tag"])
    _refresh_expand_window(kind)


def expand_main(sender=None, app_data=None):
    _open_expand_window("main")


def expand_secondary(sender=None, app_data=None):
    _open_expand_window("secondary")


def _refresh_open_expand_windows():
    """Called after every RUN (see run_device()): keep any already-open
    expanded window in sync with the new curves, without popping one open
    that the user hasn't asked for."""
    for kind, spec in EXPAND_SPECS.items():
        if dpg.does_item_exist(spec["window_tag"]) and dpg.is_item_shown(spec["window_tag"]):
            _refresh_expand_window(kind)


def report_bundle(outdir_override=None):
    """'Report bundle' button: hand the current run + every stored slot to
    fsim_viz.report.designer_report. Imported lazily here (not at module
    scope) so a normal designer launch never pays matplotlib's import cost."""
    if LAST_RUN is None:
        dpg.set_value("status", "nothing run yet -- press RUN before building a report")
        return
    from fsim_viz.report import designer_report

    entries = [dict(LAST_RUN, label="current")]
    entries += [dict(SLOTS[lbl], label=lbl) for lbl in SLOT_LABELS if lbl in SLOTS]
    outdir = Path(outdir_override) if outdir_override is not None else (
        ROOT / "out" / "designer" / f"{LAST_RUN['name']}-review")
    designer_report(entries, outdir, title=f"{LAST_RUN['name']} design review")
    try:
        shown = outdir.relative_to(ROOT)
    except ValueError:
        shown = outdir
    dpg.set_value("status", f"report bundle -> {shown}")
    return outdir


def save_design():
    d = collect_design()
    path = ROOT / "cards" / f"{d.name}-design.yaml"
    d.save(path)
    dpg.set_value("status", f"saved {path.name} (reproducible by anyone holding it)")


def load_design():
    name = dpg.get_value("design.name")
    path = ROOT / "cards" / f"{name}-design.yaml"
    if not path.exists():
        dpg.set_value("status", f"no such design: {path.name}")
        return
    apply_design(DeviceDesign.load(path))
    dpg.set_value("status", f"loaded {path.name}")


def _save_dialog_cb(sender, app_data):
    path = app_data["file_path_name"]
    if not path.lower().endswith(".yaml"):
        path += ".yaml"
    collect_design().save(path)
    dpg.set_value("status", f"saved {Path(path).name} (reproducible by anyone holding it)")


def _load_dialog_cb(sender, app_data):
    path = app_data["file_path_name"]
    apply_design(DeviceDesign.load(path))
    dpg.set_value("status", f"loaded {Path(path).name}")


def export_bundle():
    d = collect_design()
    res = evaluate(d)
    out = ROOT / "out" / "designer" / d.name
    out.mkdir(parents=True, exist_ok=True)
    d.save(out / "design.yaml")
    c = res["curves"]
    with open(out / "curves.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_hs_K", "T_j_K", "g2", "eps", "rho2", "gamma_meV"])
        w.writerows(zip(c["T_hs"], c["Tj"], c["g2"], c["eps"], c["rho2"], c["gamma"]))
    with open(out / "scalars.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for k, v in res["scalars"].items():
            w.writerow([k, v])
    dpg.set_value("status", f"bundle -> {out.relative_to(ROOT)} (design + CSVs)")


# --------------------------------------------------------------------------- UI

def build_ui():
    with dpg.window(tag="main"):
        with dpg.file_dialog(directory_selector=False, show=False, modal=True,
                             callback=_save_dialog_cb, tag="save_file_dialog",
                             default_path=str(ROOT / "cards"), width=700, height=400):
            dpg.add_file_extension(".yaml")
            dpg.add_file_extension(".*")
        with dpg.file_dialog(directory_selector=False, show=False, modal=True,
                             callback=_load_dialog_cb, tag="load_file_dialog",
                             default_path=str(ROOT / "cards"), width=700, height=400):
            dpg.add_file_extension(".yaml")
            dpg.add_file_extension(".*")

        with dpg.group(horizontal=True):
            dpg.add_input_text(tag="design.name", default_value="staged-device",
                               width=160)
            dpg.add_button(label="Load", callback=load_design)
            dpg.add_button(label="Save", callback=save_design)
            dpg.add_button(label="Save As...", callback=lambda: dpg.show_item("save_file_dialog"))
            dpg.add_button(label="Load...", callback=lambda: dpg.show_item("load_file_dialog"))
            dpg.add_button(label="Export bundle", callback=export_bundle)
            dpg.add_button(label="  RUN  ", callback=run_device, tag="run_button")
            dpg.add_text("tag chain [A] - every result inherits unmeasured inputs",
                         color=INK2)
        with dpg.group(horizontal=True):
            dpg.add_text("Compare:")
            dpg.add_button(label="Store A", callback=lambda: store_slot("A"))
            dpg.add_button(label="Store B", callback=lambda: store_slot("B"))
            dpg.add_button(label="Store C", callback=lambda: store_slot("C"))
            dpg.add_button(label="Clear slots", callback=clear_slots)
            dpg.add_button(label="Report bundle", callback=lambda: report_bundle())
        with dpg.group(horizontal=True):
            dpg.add_text("Presets:")
            dpg.add_combo(list(presets.DOT_PRESETS), default_value="chatzarakis-class",
                         tag="preset.dot", width=160)
            dpg.add_combo(list(presets.TEMPLATE_PRESETS), default_value="GaAs",
                         tag="preset.template", width=110)
            dpg.add_combo(list(presets.CAVITY_PRESETS), default_value="none",
                         tag="preset.cavity", width=140)
            dpg.add_combo(list(presets.DRIVE_PRESETS), default_value="cw-electrical",
                         tag="preset.drive", width=150)
            dpg.add_combo(list(presets.INJECTION_PRESETS), default_value="standard",
                         tag="preset.injection", width=150)
            dpg.add_button(label="Apply presets", callback=apply_presets)
        dpg.add_text("", tag="status", color=AMBER)
        dpg.add_text("", tag="param_warning", color=AMBER)

        with dpg.group(horizontal=True):
            # ---------------- left: block diagram (node columns 320 px apart so
            # the widest rows, "[E]" tags and "range (envelope)", are not
            # clipped by the next column; studio-07 F8)
            with dpg.child_window(width=900, height=640):
                with dpg.node_editor(width=-1, height=620, tag="editor",
                                     minimap=False):
                    with dpg.node(label="ELECTRICAL DRIVE", pos=(10, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="V", tag="drive.V", width=90,
                                                    default_value=1.9,
                                                    callback=_mk_cb("drive.V"))
                                _tag_bullet("drive.V")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="I (uA)", tag="drive.I_uA",
                                                    width=90, default_value=10.0,
                                                    callback=_mk_cb("drive.I_uA"))
                                _tag_bullet("drive.I_uA")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="duty", tag="drive.duty",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("drive.duty"))
                                _tag_bullet("drive.duty")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="mu/pulse", tag="drive.mu",
                                                    width=90, default_value=0.5,
                                                    callback=_mk_cb("drive.mu",
                                                                   extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.mu")
                            _range_controls("drive.mu", extra=lambda s, v: _update_drive_panel())
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="b_e", tag="drive.b_e",
                                                    width=90, default_value=0.02,
                                                    callback=_mk_cb("drive.b_e"))
                                _tag_bullet("drive.b_e")
                            _range_controls("drive.b_e")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="b_e exp m", tag="drive.b_e_m",
                                                    width=90, default_value=1.5,
                                                    callback=_mk_cb("drive.b_e_m"))
                                _tag_bullet("drive.b_e_m")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="b_e Eact", tag="drive.b_e_Eact",
                                                    width=90, default_value=100.0,
                                                    callback=_mk_cb("drive.b_e_Eact"))
                                _tag_bullet("drive.b_e_Eact")
                            dpg.add_separator()
                            with dpg.group(horizontal=True):
                                dpg.add_combo(("EL", "PL"), label="mode", tag="drive.mode",
                                             default_value="EL", width=60,
                                             callback=_mk_cb("drive.mode",
                                                            extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.mode")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="dGamma_inj (meV)", tag="drive.dg_inj",
                                                    width=90, default_value=0.0,
                                                    callback=_mk_cb("drive.dg_inj",
                                                                   extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.dg_inj")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="p_inj", tag="drive.p_inj",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("drive.p_inj"))
                                _tag_bullet("drive.p_inj")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="pump Fano F_p", tag="drive.F_p",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("drive.F_p",
                                                                   extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.F_p")
                            _range_controls("drive.F_p")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="eta_capture", tag="drive.eta_capture",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("drive.eta_capture",
                                                                   extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.eta_capture")
                            _range_controls("drive.eta_capture")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="C_dep (pF)", tag="drive.C_dep_pF",
                                                    width=90, default_value=10.0,
                                                    callback=_mk_cb("drive.C_dep_pF",
                                                                   extra=lambda s, v: _update_drive_panel()))
                                _tag_bullet("drive.C_dep_pF")
                            dpg.add_text("", tag="drive.N_info", color=(150, 150, 150))
                            dpg.add_text("", tag="drive.f8_warning", color=AMBER)
                            dpg.add_text("", tag="drive.pl_note", color=(150, 150, 150))
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="drive_out"):
                            dpg.add_text("I, V, heat")

                    with dpg.node(label="MESA / THERMAL", pos=(10, 700)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="th_in"):
                            dpg.add_text("P = duty * I * V")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="mesa d (um)", tag="th.mesa",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("thermal.mesa_diameter_um",
                                                                   extra=lambda s, v: draw_cross_section()))
                                _tag_bullet("thermal.mesa_diameter_um")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="T_hs (K)", tag="th.T_hs",
                                                    width=90, default_value=77.0,
                                                    callback=_mk_cb("thermal.T_hs",
                                                                   extra=lambda s, v: (
                                                                       draw_cross_section(),
                                                                       _update_drive_panel())))
                                _tag_bullet("thermal.T_hs")
                            dpg.add_text("stack: edit in Fab stack panel ->")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="th_out"):
                            dpg.add_text("T_j")

                    with dpg.node(label="QD EMITTER", pos=(330, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="dot_in"):
                            dpg.add_text("carriers @ T_j")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="Delta_XX (meV)",
                                                    tag="dot.delta_xx", width=90,
                                                    default_value=3.5,
                                                    callback=_mk_cb("dot.delta_xx"))
                                _tag_bullet("dot.delta_xx")
                            _range_controls("dot.delta_xx")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="Gamma scale",
                                                    tag="dot.gamma_scale", width=90,
                                                    default_value=1.0,
                                                    callback=_mk_cb("dot.gamma_scale"))
                                _tag_bullet("dot.gamma_scale")
                            _range_controls("dot.gamma_scale")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="r_XX", tag="dot.r_xx",
                                                    width=90, default_value=0.72,
                                                    callback=_mk_cb("dot.r_xx"))
                                _tag_bullet("dot.r_xx")
                            dpg.add_text("Gamma(T), retention: class proxy [A]",
                                         color=INK2)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="dot_out"):
                            dpg.add_text("X + XX photons")

                    with dpg.node(label="CAVITY (F6)", pos=(330, 330)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="cav_in"):
                            dpg.add_text("photons")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_checkbox(label="enabled", tag="cav.enabled",
                                                 default_value=False,
                                                 callback=_mk_cb("cavity.enabled"))
                                _tag_bullet("cavity.enabled")
                            with dpg.group(horizontal=True):
                                dpg.add_combo(("planar", "sin_waveguide"), label="type",
                                             tag="cav.type", default_value="planar",
                                             width=110,
                                             callback=_mk_cb("cavity.type",
                                                            extra=lambda s, v: dpg.configure_item(
                                                                "lemma1_note", show=(v == "sin_waveguide"))))
                                _tag_bullet("cavity.type")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="kappa (meV)", tag="cav.kappa",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("cavity.kappa"))
                                _tag_bullet("cavity.kappa")
                            _range_controls("cavity.kappa")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="track T (K)", tag="cav.T_track",
                                                    width=90, default_value=120.0,
                                                    callback=_mk_cb("cavity.T_track"))
                                _tag_bullet("cavity.T_track")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="E_X0 (eV)", tag="cav.E_X0",
                                                    width=90, default_value=1.88,
                                                    callback=_mk_cb("cavity.E_X0"))
                                _tag_bullet("cavity.E_X0")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="F_P", tag="cav.F_P",
                                                    width=90, default_value=10.0,
                                                    callback=_mk_cb("cavity.F_P"))
                                _tag_bullet("cavity.F_P")
                            _range_controls("cavity.F_P")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="gain G", tag="cav.G",
                                                    width=90, default_value=8.0,
                                                    callback=_mk_cb("cavity.G"))
                                _tag_bullet("cavity.G")
                            _range_controls("cavity.G")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="beta_sin", tag="cav.beta_sin",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("cavity.beta_sin"))
                                _tag_bullet("cavity.beta_sin")
                            dpg.add_text("Lemma 1: beta -> brightness only, never g2",
                                        tag="lemma1_note", color=RED, show=False)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="cav_out"):
                            dpg.add_text("filtered + boosted")

                    with dpg.node(label="SLIT FILTER", pos=(650, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="fil_in"):
                            dpg.add_text("spectrum")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_checkbox(label="enabled", tag="fil.enabled",
                                                 default_value=True,
                                                 callback=_mk_cb("filter.enabled"))
                                _tag_bullet("filter.enabled")
                            with dpg.group(horizontal=True):
                                dpg.add_checkbox(label="auto w = Gamma(T_j)",
                                                 tag="fil.auto_w", default_value=True,
                                                 callback=_mk_cb("filter.auto_w"))
                                _tag_bullet("filter.auto_w")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="w (meV)", tag="fil.w",
                                                    width=90, default_value=2.0,
                                                    callback=_mk_cb("filter.w"))
                                _tag_bullet("filter.w")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="dx (meV)", tag="fil.dx",
                                                    width=90, default_value=0.0,
                                                    callback=_mk_cb("filter.dx"))
                                _tag_bullet("filter.dx")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="fil_out"):
                            dpg.add_text("to detector")

                    with dpg.node(label="APERTURE / ENSEMBLE (F5)", pos=(650, 330)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="log10 density", tag="ap.log_density",
                                                    width=90, default_value=8.845,
                                                    callback=_mk_cb("aperture.density_cm2",
                                                                   extra=_mark_density_edited,
                                                                   transform=lambda v: 10.0 ** v))
                                _tag_bullet("aperture.density_cm2")
                            _range_controls("aperture.density_cm2")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="aperture (um)", tag="ap.diam",
                                                    width=90, default_value=1.0,
                                                    callback=_mk_cb("aperture.diameter_um"))
                                _tag_bullet("aperture.diameter_um")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="sigma_inh", tag="ap.sigma",
                                                    width=90, default_value=40.0,
                                                    callback=_mk_cb("aperture.sigma_inh"))
                                _tag_bullet("aperture.sigma_inh")
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="comp. bright r", tag="ap.r",
                                                    width=90, default_value=0.3,
                                                    callback=_mk_cb("aperture.comp_brightness"))
                                _tag_bullet("aperture.comp_brightness")

                dpg.add_node_link("drive_out", "th_in", parent="editor")
                dpg.add_node_link("th_out", "dot_in", parent="editor")
                dpg.add_node_link("dot_out", "cav_in", parent="editor")
                dpg.add_node_link("cav_out", "fil_in", parent="editor")

            # ---------------- middle: fab stack + cross-section
            with dpg.child_window(width=370, height=640):
                dpg.add_text("Fab stack (top -> substrate)")
                with dpg.table(tag="stack_table", header_row=True,
                               policy=dpg.mvTable_SizingFixedFit):
                    pass
                dpg.add_button(label="+ add layer", callback=lambda: add_layer())
                dpg.add_spacer(height=6)
                dpg.add_text("Cross-section (to scale in t)")
                dpg.add_drawlist(width=340, height=260, tag="xsec")

            # ---------------- right: results
            with dpg.child_window(width=-1, height=640):
                dpg.add_text("", tag="envelope_header", color=AMBER)
                dpg.add_button(label="[ expand ]", callback=expand_main)
                with dpg.plot(height=270, width=-1, tag="plot_main"):
                    dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                    dpg.add_plot_axis(dpg.mvXAxis, label="heatsink T (K)", tag="xax1")
                    with dpg.plot_axis(dpg.mvYAxis, label="g2 / fractions", tag="yax1"):
                        # T_c interval (envelope runs): ref-wash span lo..hi
                        dpg.add_shade_series([], [], y2=[], label="T_c interval",
                                             tag="s_tc_band")
                        dpg.add_shade_series([], [], y2=[], label="g2 band", tag="s_g2_band")
                        dpg.add_shade_series([], [], y2=[], label="rho^2 band", tag="s_rho2_band")
                        dpg.add_line_series([], [], label="g2(0)", tag="s_g2")
                        dpg.add_line_series([], [], label="eps", tag="s_eps")
                        dpg.add_line_series([], [], label="rho^2", tag="s_rho2")
                        dpg.add_line_series([], [], label="ceiling 0.5", tag="s_half")
                        dpg.add_line_series([], [], label="T_c", tag="s_tc")
                        # D3 comparison slots: pre-created hidden, shown/colored
                        # on store (see _refresh_slot_overlay / _build_slot_themes)
                        for _lbl in SLOT_LABELS:
                            _lt, _bt = _slot_tags(_lbl)
                            dpg.add_shade_series([], [], y2=[], label=f"{_lbl} band",
                                                 tag=_bt, show=False)
                            dpg.add_line_series([], [], label=_lbl, tag=_lt, show=False)
                    dpg.add_plot_annotation(label="", default_value=(0.0, 0.95), tag="tc_note",
                                            color=(0, 0, 0, 0), show=False)
                dpg.add_button(label="[ expand ]", callback=expand_secondary)
                # dT_J (K) and Gamma (meV) have different units: two stacked
                # plots on one linked x axis, never one shared y axis
                with dpg.subplots(2, 1, height=230, width=-1, link_all_x=True,
                                  row_ratios=[1.0, 1.0]):
                    with dpg.plot():
                        dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                        dpg.add_plot_axis(dpg.mvXAxis, label="", tag="xax2",
                                          no_tick_labels=True)
                        with dpg.plot_axis(dpg.mvYAxis, label="dT_J (K)", tag="yax2"):
                            dpg.add_shade_series([], [], y2=[], label="dT_J band",
                                                 tag="s_tj_band")
                            dpg.add_line_series([], [], label="dT_J", tag="s_tj")
                    with dpg.plot():
                        dpg.add_plot_legend(outside=True, location=dpg.mvPlot_Location_East)
                        dpg.add_plot_axis(dpg.mvXAxis, label="heatsink T (K)", tag="xax3")
                        with dpg.plot_axis(dpg.mvYAxis, label="Gamma (meV)", tag="yax3"):
                            dpg.add_line_series([], [], label="Gamma(T_j)", tag="s_gam")
                dpg.add_text("zoom: scroll | box: right-drag | reset: double-click",
                             color=(150, 150, 150))
                dpg.add_text("", tag="warn_runaway", color=RED)
                dpg.add_text("", tag="warn_eps", color=RED)
                dpg.add_text("press RUN", tag="results_text")
                dpg.add_text("", tag="delta_text", color=INK2)
    _build_slot_themes()
    _build_series_themes()
    _g2_axis_limits("yax1")
    # RUN is the one filled control: beam ink (DIRECTION.md, reserved colour)
    with dpg.theme() as run_theme:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button, _rgb(_C["beam"]))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, _rgb(_C["beam"]))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, _rgb(_C["beam"]))
            dpg.add_theme_color(dpg.mvThemeCol_Text, _rgb(_C["ground"]))
    dpg.bind_item_theme("run_button", run_theme)
    # numeric readouts in the mono face (columns line up, tabular digits)
    if FONTS.get("mono") is not None:
        for tag in ("results_text", "delta_text"):
            dpg.bind_item_font(tag, FONTS["mono"])


def _run_selftest(outdir: Path) -> bool:
    """D3 selftest: drive the real button callbacks (not a reimplementation)
    through one full compare-and-report cycle -- preset -> RUN (envelope
    defaults on) -> Store A -> tweak delta_xx -> RUN -> Store B -> Report
    bundle -- then check the bundle landed with the files a reviewer needs."""
    dpg.set_value("preset.dot", "staged-inp-gaasp")
    dpg.set_value("preset.template", "GaAs")
    dpg.set_value("preset.cavity", "none")
    dpg.set_value("preset.drive", "cw-electrical")
    apply_presets()
    run_device()
    store_slot("A")
    dpg.set_value("dot.delta_xx", 5.0)   # swept param: shows up in the YAML
    dpg.set_value("drive.I_uA", 50.0)    # NON-swept param: curves must diverge
    run_device()
    store_slot("B")
    report_bundle(outdir_override=outdir)
    for _ in range(5):
        dpg.render_dearpygui_frame()

    required = ["report.png", "report.pdf", "curves_A.csv", "curves_B.csv",
               "scalars_comparison.csv", "design_A.yaml", "design_B.yaml"]
    missing = [f for f in required if not (Path(outdir) / f).exists()]
    ok = not missing
    if ok:
        # the comparison must compare: A and B must not be byte-identical
        # (Opus D3 finding -- an all-swept edit yields identical curves)
        a = (Path(outdir) / "curves_A.csv").read_bytes()
        b = (Path(outdir) / "curves_B.csv").read_bytes()
        if a == b:
            ok = False
            print("selftest: FAIL  slots A and B are identical -- comparison "
                  "path not exercised")
    print(f"selftest: {'OK' if ok else 'FAIL'}  bundle -> {outdir}"
          + (f"  missing: {missing}" if missing else ""))
    return ok


def _run_roundtrip_check() -> bool:
    """v1.1 gate check (gate_v11_gui.py check c): apply_design() a design
    carrying the six new DriveBlock fields THROUGH THE REAL WIDGETS, then
    collect_design() and assert every field survives the round trip --
    exercises the widget wiring itself, not a reimplementation of it.

    pr-pkg1-fix3 item 1: `d` here is a bare DeviceDesign() -- its
    aperture.density_cm2 is None (the legacy/unset default), the exact
    case that used to crash apply_design()'s math.log10(None) before this
    fix (verify_designer_rt.py dropped to 4/6, gate_v11_gui.py's check (c)
    subprocess failed). Assert the untouched round trip preserves that
    None, not just that it no longer crashes.

    pr-pkg1-fix4 item 3: after the untouched check above passes, drive the
    SAME live aperture.density_cm2 widget through a real edit landing on
    7.0e8 (invoking its registered callback directly -- dpg.set_value()
    alone does not fire callbacks, so this is the one way to exercise the
    actual widget wiring, not a reimplementation of it, in a headless
    subprocess) and assert collect_design() now returns 7.0e8, not the
    None a value-only "still reads the display default" comparison used
    to silently fall back to."""
    import math
    d = DeviceDesign()
    d.drive.mode = "PL"
    d.drive.dg_inj = 2.0
    d.drive.F_p = 0.5
    d.drive.eta_capture = 0.7
    apply_design(d)
    for _ in range(2):
        dpg.render_dearpygui_frame()
    d2 = collect_design()
    tol = 1e-5  # add_input_float stores float32 internally (~7 sig figs)
    ok = (d2.drive.mode == d.drive.mode
          and abs(d2.drive.dg_inj - d.drive.dg_inj) < tol
          and abs(d2.drive.p_inj - d.drive.p_inj) < tol
          and abs(d2.drive.F_p - d.drive.F_p) < tol
          and abs(d2.drive.eta_capture - d.drive.eta_capture) < tol
          and abs(d2.drive.C_dep_pF - d.drive.C_dep_pF) < tol
          and d2.aperture.density_cm2 is None)

    edit_value = math.log10(7.0e8)
    dpg.set_value("ap.log_density", edit_value)
    dpg.get_item_configuration("ap.log_density")["callback"]("ap.log_density", edit_value)
    for _ in range(2):
        dpg.render_dearpygui_frame()
    d3 = collect_design()
    density_edit_ok = (d3.aperture.density_cm2 is not None
                       and abs(d3.aperture.density_cm2 - 7.0e8) / 7.0e8 < tol)
    ok = ok and density_edit_ok

    print(f"roundtrip-check: {'OK' if ok else 'FAIL'}  "
          f"mode={d2.drive.mode} dg_inj={d2.drive.dg_inj} p_inj={d2.drive.p_inj} "
          f"F_p={d2.drive.F_p} eta_capture={d2.drive.eta_capture} "
          f"C_dep_pF={d2.drive.C_dep_pF} aperture.density_cm2={d2.aperture.density_cm2!r}  "
          f"density_edit_to_7e8={d3.aperture.density_cm2!r}")
    return ok


def main(frames=None, selftest_outdir=None, roundtrip_check=False, *,
         design_path=None, screenshot=None, collect_dump=None):
    """`design_path`, `screenshot` and `collect_dump` are keyword-only (all
    None reproduces the pre-existing behavior exactly): `design_path`
    replaces the DESIGN_PATH default; `screenshot` writes a real PNG (plus
    the verbatim results-panel text at `screenshot + ".txt"`) once `frames`
    have rendered; `collect_dump` renders a couple of frames after
    apply_design() and saves collect_design() to that path, exiting before
    the event loop -- both without importing anything the GUI doesn't
    already use (three-layer rule)."""
    dpg.create_context()
    FONTS.update(bind_theme(dpg, "dark"))
    build_ui()
    dpath = Path(design_path) if design_path is not None else DESIGN_PATH
    default = DeviceDesign()
    if dpath.exists():
        default = DeviceDesign.load(dpath)
    apply_design(default)
    dpg.create_viewport(title="FSIM device designer", width=1720, height=760)
    dpg.setup_dearpygui()
    dpg.show_viewport()
    dpg.set_primary_window("main", True)
    if collect_dump is not None:
        # No run_device()/evaluate() here on purpose: a --collect-dump design
        # is only being checked for widget round-trip fidelity (spec
        # rt-edge-gui-preservation), and an arbitrary probe design may not be
        # a physically consistent one evaluate() would accept.
        for _ in range(2):
            dpg.render_dearpygui_frame()
        collect_design().save(collect_dump)
        dpg.destroy_context()
        sys.exit(0)
    if roundtrip_check:
        ok = _run_roundtrip_check()
        dpg.destroy_context()
        sys.exit(0 if ok else 1)
    if selftest_outdir is not None:
        ok = _run_selftest(Path(selftest_outdir))
        dpg.destroy_context()
        sys.exit(0 if ok else 1)
    if frames:
        run_device()  # exercise the full pipeline once
        for _ in range(frames):
            dpg.render_dearpygui_frame()
        if screenshot is not None:
            out_path = Path(screenshot)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            dpg.output_frame_buffer(str(out_path))  # async: needs further frames to flush
            for _ in range(10):
                dpg.render_dearpygui_frame()
            Path(str(out_path) + ".txt").write_text(
                dpg.get_value("results_text"), encoding="utf-8")
            print("screenshot ->", out_path)
        print("smoke: rendered", frames, "frames;",
              dpg.get_value("results_text").splitlines()[2].strip())
        # expand-to-window: drive both [ expand ] callbacks directly, then
        # confirm the popouts exist and their g2 series copied the full
        # source data (see EXPAND_SPECS / _refresh_expand_window).
        expand_main()
        expand_secondary()
        for _ in range(3):
            dpg.render_dearpygui_frame()
        src_len = len(dpg.get_value("s_g2")[0])
        exp_len = len(dpg.get_value("x_s_g2")[0])
        assert dpg.does_item_exist("expand_win_main"), "expand: main popout missing"
        assert dpg.does_item_exist("expand_win_secondary"), "expand: secondary popout missing"
        assert exp_len == src_len, f"expand: g2 length mismatch src={src_len} expanded={exp_len}"
        print(f"smoke: expand popouts OK; g2 series len {src_len} matches expanded copy")
    else:
        dpg.start_dearpygui()
    dpg.destroy_context()


if __name__ == "__main__":
    n = None
    selftest_outdir = None
    if "--frames" in sys.argv:
        n = int(sys.argv[sys.argv.index("--frames") + 1])
    if "--selftest" in sys.argv:
        selftest_outdir = sys.argv[sys.argv.index("--selftest") + 1]
    roundtrip_check = "--roundtrip-check" in sys.argv
    design_path = None
    if "--design" in sys.argv:
        design_path = sys.argv[sys.argv.index("--design") + 1]
    screenshot = None
    if "--screenshot" in sys.argv:
        screenshot = sys.argv[sys.argv.index("--screenshot") + 1]
    collect_dump = None
    if "--collect-dump" in sys.argv:
        collect_dump = sys.argv[sys.argv.index("--collect-dump") + 1]
    main(frames=n, selftest_outdir=selftest_outdir, roundtrip_check=roundtrip_check,
         design_path=design_path, screenshot=screenshot, collect_dump=collect_dump)
