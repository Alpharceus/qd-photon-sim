"""fsim device designer (Dear PyGui): configure the device as a block diagram,
edit the fab stack, run, get graphs and numbers.

Multisim-style with one honest restriction: the physics defines ONE signal
chain (dot -> cavity -> filter -> detection) with drive and thermal attached,
so the canvas is that chain with configurable blocks -- not free-form wiring.

Three-layer rule: this file computes NO physics. It builds a DeviceDesign,
calls fsim_core.device.evaluate, and renders the results. Designs round-trip
to YAML in cards/ so every session is reproducible from its file.

Run:        python fsim_gui/designer.py
Smoke test: python fsim_gui/designer.py --frames 5
"""
import csv
import sys
from pathlib import Path

import dearpygui.dearpygui as dpg

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import design_meta, presets  # noqa: E402
from fsim_core.device import DeviceDesign, evaluate, evaluate_envelope  # noqa: E402

DESIGN_PATH = ROOT / "cards" / "staged-device-design.yaml"
LAYERS = []          # live fab-stack rows (list of dicts)
SUBSTRATE = {"name": "GaAs", "k300": 55.0, "alpha": 1.25}

GREEN = (86, 166, 50)
AMBER = (227, 162, 26)
RED = (215, 25, 28)
BLUE = (60, 120, 216)
PURPLE = (120, 90, 200)
TAG_COLOR = {"V": GREEN, "DR": AMBER, "E": AMBER, "A": RED}

# meta path -> live widget tag, for every parameter widget that carries a META
# entry. aperture.density_cm2 is edited as log10 in the widget (see transform
# in VIOLATIONS handling below).
WIDGET_TAG = {
    "dot.delta_xx": "dot.delta_xx", "dot.gamma_scale": "dot.gamma_scale",
    "dot.r_xx": "dot.r_xx",
    "drive.V": "drive.V", "drive.I_uA": "drive.I_uA", "drive.duty": "drive.duty",
    "drive.mu": "drive.mu", "drive.b_e": "drive.b_e", "drive.b_e_m": "drive.b_e_m",
    "drive.b_e_Eact": "drive.b_e_Eact",
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
    color = TAG_COLOR.get(meta["tag"], AMBER)
    b = dpg.add_text("*", color=color)
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


def _range_controls(path):
    """Envelope toggle for a parameter with an ENV_DEFAULTS entry: checkbox
    '~' + lo/hi range inputs, shown only when checked. No-op for paths
    without an entry. Program decision: [A]-tagged params default RANGED on
    (see design_meta.default_ranged); apply_design() re-syncs the checked
    state and lo/hi values every time a design is loaded."""
    bounds = design_meta.ENV_DEFAULTS.get(path)
    if bounds is None:
        return
    lo, hi = bounds
    on = path in design_meta.default_ranged(DeviceDesign())
    with dpg.group(horizontal=True):
        dpg.add_checkbox(label="~", tag=f"rng.{path}", default_value=on,
                         callback=lambda s, v: _toggle_range(path))
        dpg.add_input_float(tag=f"rng.{path}.lo", width=68, default_value=lo,
                            show=on)
        dpg.add_input_float(tag=f"rng.{path}.hi", width=68, default_value=hi,
                            show=on)
        dpg.add_text("range (envelope)", color=(150, 150, 150))


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
        lo, hi = bounds
        dpg.set_value(tag, on)
        dpg.set_value(f"{tag}.lo", lo)
        dpg.set_value(f"{tag}.hi", hi)
        dpg.configure_item(f"{tag}.lo", show=on)
        dpg.configure_item(f"{tag}.hi", show=on)


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
    d = DeviceDesign()
    d.name = dpg.get_value("design.name")
    d.dot.delta_xx = dpg.get_value("dot.delta_xx")
    d.dot.gamma_scale = dpg.get_value("dot.gamma_scale")
    d.dot.r_xx = dpg.get_value("dot.r_xx")
    d.drive.V = dpg.get_value("drive.V")
    d.drive.I_uA = dpg.get_value("drive.I_uA")
    d.drive.duty = dpg.get_value("drive.duty")
    d.drive.mu = dpg.get_value("drive.mu")
    d.drive.b_e = dpg.get_value("drive.b_e")
    d.drive.b_e_m = dpg.get_value("drive.b_e_m")
    d.drive.b_e_Eact = dpg.get_value("drive.b_e_Eact")
    d.thermal.mesa_diameter_um = dpg.get_value("th.mesa")
    d.thermal.T_hs = dpg.get_value("th.T_hs")
    d.thermal.layers = [dict(L) for L in LAYERS]
    d.thermal.substrate = dict(SUBSTRATE)
    d.cavity.enabled = dpg.get_value("cav.enabled")
    d.cavity.type = dpg.get_value("cav.type")
    d.cavity.kappa = dpg.get_value("cav.kappa")
    d.cavity.T_track = dpg.get_value("cav.T_track")
    d.cavity.E_X0 = dpg.get_value("cav.E_X0")
    d.cavity.F_P = dpg.get_value("cav.F_P")
    d.cavity.G = dpg.get_value("cav.G")
    d.cavity.beta_sin = dpg.get_value("cav.beta_sin")
    d.filter.enabled = dpg.get_value("fil.enabled")
    d.filter.auto_w = dpg.get_value("fil.auto_w")
    d.filter.w = dpg.get_value("fil.w")
    d.filter.dx = dpg.get_value("fil.dx")
    d.aperture.density_cm2 = 10.0 ** dpg.get_value("ap.log_density")
    d.aperture.diameter_um = dpg.get_value("ap.diam")
    d.aperture.sigma_inh = dpg.get_value("ap.sigma")
    d.aperture.comp_brightness = dpg.get_value("ap.r")
    return d


def apply_design(d: DeviceDesign):
    global LAYERS, SUBSTRATE
    import math
    dpg.set_value("design.name", d.name)
    dpg.set_value("dot.delta_xx", d.dot.delta_xx)
    dpg.set_value("dot.gamma_scale", d.dot.gamma_scale)
    dpg.set_value("dot.r_xx", d.dot.r_xx)
    dpg.set_value("drive.V", d.drive.V)
    dpg.set_value("drive.I_uA", d.drive.I_uA)
    dpg.set_value("drive.duty", d.drive.duty)
    dpg.set_value("drive.mu", d.drive.mu)
    dpg.set_value("drive.b_e", d.drive.b_e)
    dpg.set_value("drive.b_e_m", d.drive.b_e_m)
    dpg.set_value("drive.b_e_Eact", d.drive.b_e_Eact)
    dpg.set_value("th.mesa", d.thermal.mesa_diameter_um)
    dpg.set_value("th.T_hs", d.thermal.T_hs)
    dpg.set_value("cav.enabled", d.cavity.enabled)
    dpg.set_value("cav.type", d.cavity.type)
    dpg.set_value("cav.kappa", d.cavity.kappa)
    dpg.set_value("cav.T_track", d.cavity.T_track)
    dpg.set_value("cav.E_X0", d.cavity.E_X0)
    dpg.set_value("cav.F_P", d.cavity.F_P)
    dpg.set_value("cav.G", d.cavity.G)
    dpg.set_value("cav.beta_sin", d.cavity.beta_sin)
    if dpg.does_item_exist("lemma1_note"):
        dpg.configure_item("lemma1_note", show=(d.cavity.type == "sin_waveguide"))
    dpg.set_value("fil.enabled", d.filter.enabled)
    dpg.set_value("fil.auto_w", d.filter.auto_w)
    dpg.set_value("fil.w", d.filter.w)
    dpg.set_value("fil.dx", d.filter.dx)
    dpg.set_value("ap.log_density", math.log10(d.aperture.density_cm2))
    dpg.set_value("ap.diam", d.aperture.diameter_um)
    dpg.set_value("ap.sigma", d.aperture.sigma_inh)
    dpg.set_value("ap.r", d.aperture.comp_brightness)
    LAYERS = [dict(L) for L in d.thermal.layers]
    SUBSTRATE = dict(d.thermal.substrate)
    rebuild_stack_table()
    draw_cross_section()
    revalidate_all()
    _reset_range_controls(d)


def apply_presets():
    """Apply the four preset combos onto the live design, then push to widgets."""
    d = collect_design()
    presets.apply_dot_preset(d, dpg.get_value("preset.dot"))
    presets.apply_template_preset(d, dpg.get_value("preset.template"))
    presets.apply_cavity_preset(d, dpg.get_value("preset.cavity"))
    presets.apply_drive_preset(d, dpg.get_value("preset.drive"))
    apply_design(d)
    dpg.set_value("status",
                  f"applied presets: dot={dpg.get_value('preset.dot')}  "
                  f"template={dpg.get_value('preset.template')}  "
                  f"cavity={dpg.get_value('preset.cavity')}  "
                  f"drive={dpg.get_value('preset.drive')}")


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


def _run_point(d: DeviceDesign):
    res = evaluate(d)
    c, s = res["curves"], res["scalars"]
    Ts = list(map(float, c["T_hs"]))
    dpg.set_value("s_g2", [Ts, list(map(float, c["g2"]))])
    dpg.set_value("s_eps", [Ts, list(map(float, c["eps"]))])
    dpg.set_value("s_rho2", [Ts, list(map(float, c["rho2"]))])
    dpg.set_value("s_half", [[Ts[0], Ts[-1]], [0.5, 0.5]])
    tc = s["T_c"]
    dpg.set_value("s_tc", [[tc, tc], [0.0, 1.0]] if tc == tc else [[], []])
    dpg.set_value("s_tj", [Ts, [tj - t for tj, t in zip(map(float, c["Tj"]), Ts)]])
    dpg.set_value("s_gam", [Ts, list(map(float, c["gamma"]))])
    dpg.set_value("s_g2_band", [[], [], []])
    dpg.set_value("s_rho2_band", [[], [], []])
    dpg.set_value("s_tj_band", [[], [], []])
    dpg.set_value("envelope_header", "")
    dpg.fit_axis_data("xax1"); dpg.fit_axis_data("yax1")
    dpg.fit_axis_data("xax2"); dpg.fit_axis_data("yax2")

    dpg.set_value("warn_runaway",
                  "!!!  THERMAL RUNAWAY -- no operating point  !!!" if s["runaway"] else "")
    dpg.set_value("warn_eps",
                  "!!!  CAVITY SELECTS XX (eps > 1) -- retune tracking/filter  !!!"
                  if s["eps_op"] > 1 else "")

    lines = [
        f"tag chain {s['tag_chain']}  (unmeasured inputs -> conditional numbers;"
        f" envelopes: run_phase3)",
        "",
        f"T_j at operating point   {s['T_j_op']:.1f} K   (dT_J = {s['dT_J']:.2f} K)"
        + ("   ** THERMAL RUNAWAY **" if s["runaway"] else ""),
        f"Gamma(T_j)               {s['gamma_op']:.2f} meV",
        f"eps = t_XX/t_X           {s['eps_op']:.4f}",
        f"rho (signal purity)      {s['rho_op']:.3f}",
        f"g2(0) at operating point {s['g2_op']:.3f}",
        f"brightness/pulse (t_X)   {s['brightness_per_pulse']:.3f}",
        f"master ceiling T_c       "
        + (f"{s['T_c']:.0f} K" if s["T_c"] == s["T_c"] else "not crossed in range"),
        f"F_eff (cavity)           "
        + (f"{s['F_eff']:.1f}" if s["F_eff"] == s["F_eff"] else "-- (cavity off)"),
        f"aperture: N_w = {s['N_w']:.2f}  ->  F5 g2 penalty {s['aperture_g2_penalty']:.3f}",
    ]
    dpg.set_value("results_text", "\n".join(lines))
    draw_cross_section()


def _interval(lo, hi, fmt="{:.3f}", unit=""):
    if lo != lo or hi != hi:  # NaN
        return "not reached in swept range"
    suffix = f" {unit}" if unit else ""
    return f"{fmt.format(lo)} - {fmt.format(hi)}{suffix}"


def _run_envelope(d: DeviceDesign, ranged: dict):
    env = evaluate_envelope(d, ranged)
    mid, bands, sb = env["mid"], env["bands"], env["scalar_bands"]
    Ts = list(map(float, mid["T_hs"]))

    dpg.set_value("s_g2", [Ts, list(map(float, mid["g2"]))])
    dpg.set_value("s_eps", [Ts, list(map(float, mid["eps"]))])
    dpg.set_value("s_rho2", [Ts, list(map(float, mid["rho2"]))])
    dpg.set_value("s_half", [[Ts[0], Ts[-1]], [0.5, 0.5]])
    dpg.set_value("s_tc", [[], []])  # a single marker can't honestly show a band -- see text
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
    dpg.fit_axis_data("xax1"); dpg.fit_axis_data("yax1")
    dpg.fit_axis_data("xax2"); dpg.fit_axis_data("yax2")

    dpg.set_value("warn_runaway", "")
    dpg.set_value("warn_eps", "")

    dpg.set_value("envelope_header",
                  f"ENVELOPE over {env['n_samples']} samples: [A] inputs swept "
                  "(honest mode); uncheck ~ to explore point designs")

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
        f"g2(0) at op:               {_interval(*sb['g2_op'])}",
        f"eps = t_XX/t_X at op:      {_interval(*sb['eps_op'])}",
        f"rho (signal purity) at op: {_interval(*sb['rho_op'])}",
        f"dT_J at op:                {_interval(*sb['dT_J'], fmt='{:.2f}', unit='K')}",
        f"master ceiling T_c:        {_interval(*sb['T_c'], fmt='{:.0f}', unit='K')}",
        "",
        "sensitivity (measurement priority, highest first):",
    ] + tornado_lines
    dpg.set_value("results_text", "\n".join(lines))
    draw_cross_section()


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
            dpg.add_button(label="  RUN  ", callback=run_device)
            dpg.add_text("tag chain [A] - every result inherits unmeasured inputs",
                         color=RED)
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
            dpg.add_button(label="Apply presets", callback=apply_presets)
        dpg.add_text("", tag="status", color=AMBER)
        dpg.add_text("", tag="param_warning", color=AMBER)

        with dpg.group(horizontal=True):
            # ---------------- left: block diagram
            with dpg.child_window(width=760, height=640):
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
                                                    callback=_mk_cb("drive.mu"))
                                _tag_bullet("drive.mu")
                            _range_controls("drive.mu")
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
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="drive_out"):
                            dpg.add_text("I, V, heat")

                    with dpg.node(label="MESA / THERMAL", pos=(10, 330)):
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
                                                                   extra=lambda s, v: draw_cross_section()))
                                _tag_bullet("thermal.T_hs")
                            dpg.add_text("stack: edit in Fab stack panel ->")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="th_out"):
                            dpg.add_text("T_j")

                    with dpg.node(label="QD EMITTER", pos=(270, 20)):
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
                                         color=RED)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="dot_out"):
                            dpg.add_text("X + XX photons")

                    with dpg.node(label="CAVITY (F6)", pos=(270, 330)):
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

                    with dpg.node(label="SLIT FILTER", pos=(530, 20)):
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

                    with dpg.node(label="APERTURE / ENSEMBLE (F5)", pos=(530, 330)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            with dpg.group(horizontal=True):
                                dpg.add_input_float(label="log10 density", tag="ap.log_density",
                                                    width=90, default_value=8.845,
                                                    callback=_mk_cb("aperture.density_cm2",
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
                with dpg.plot(height=270, width=-1):
                    dpg.add_plot_legend()
                    dpg.add_plot_axis(dpg.mvXAxis, label="heatsink T (K)", tag="xax1")
                    with dpg.plot_axis(dpg.mvYAxis, label="g2 / fractions", tag="yax1"):
                        dpg.add_shade_series([], [], y2=[], label="g2 band", tag="s_g2_band")
                        dpg.add_shade_series([], [], y2=[], label="rho^2 band", tag="s_rho2_band")
                        dpg.add_line_series([], [], label="g2(0)", tag="s_g2")
                        dpg.add_line_series([], [], label="eps", tag="s_eps")
                        dpg.add_line_series([], [], label="rho^2", tag="s_rho2")
                        dpg.add_line_series([], [], label="ceiling 0.5", tag="s_half")
                        dpg.add_line_series([], [], label="T_c", tag="s_tc")
                with dpg.plot(height=180, width=-1):
                    dpg.add_plot_legend()
                    dpg.add_plot_axis(dpg.mvXAxis, label="heatsink T (K)", tag="xax2")
                    with dpg.plot_axis(dpg.mvYAxis, label="dT_J (K) / Gamma (meV)",
                                       tag="yax2"):
                        dpg.add_shade_series([], [], y2=[], label="dT_J band", tag="s_tj_band")
                        dpg.add_line_series([], [], label="dT_J", tag="s_tj")
                        dpg.add_line_series([], [], label="Gamma(T_j)", tag="s_gam")
                dpg.add_text("", tag="warn_runaway", color=RED)
                dpg.add_text("", tag="warn_eps", color=RED)
                dpg.add_text("press RUN", tag="results_text")


def main(frames=None):
    dpg.create_context()
    build_ui()
    default = DeviceDesign()
    if DESIGN_PATH.exists():
        default = DeviceDesign.load(DESIGN_PATH)
    apply_design(default)
    dpg.create_viewport(title="FSIM device designer", width=1520, height=760)
    dpg.setup_dearpygui()
    dpg.show_viewport()
    dpg.set_primary_window("main", True)
    if frames:
        run_device()  # exercise the full pipeline once
        for _ in range(frames):
            dpg.render_dearpygui_frame()
        if "--screenshot" in sys.argv:
            out = sys.argv[sys.argv.index("--screenshot") + 1]
            dpg.output_frame_buffer(out)  # async: needs further frames to flush
            for _ in range(10):
                dpg.render_dearpygui_frame()
            print("screenshot ->", out)
        print("smoke: rendered", frames, "frames;",
              dpg.get_value("results_text").splitlines()[2].strip())
    else:
        dpg.start_dearpygui()
    dpg.destroy_context()


if __name__ == "__main__":
    n = None
    if "--frames" in sys.argv:
        n = int(sys.argv[sys.argv.index("--frames") + 1])
    main(frames=n)
