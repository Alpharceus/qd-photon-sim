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

from fsim_core.device import DeviceDesign, evaluate  # noqa: E402

DESIGN_PATH = ROOT / "cards" / "staged-device-design.yaml"
LAYERS = []          # live fab-stack rows (list of dicts)
SUBSTRATE = {"name": "GaAs", "k300": 55.0, "alpha": 1.25}

AMBER = (227, 162, 26)
RED = (215, 25, 28)
BLUE = (60, 120, 216)
PURPLE = (120, 90, 200)


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
    d.cavity.kappa = dpg.get_value("cav.kappa")
    d.cavity.T_track = dpg.get_value("cav.T_track")
    d.cavity.E_X0 = dpg.get_value("cav.E_X0")
    d.cavity.F_P = dpg.get_value("cav.F_P")
    d.cavity.G = dpg.get_value("cav.G")
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
    dpg.set_value("cav.kappa", d.cavity.kappa)
    dpg.set_value("cav.T_track", d.cavity.T_track)
    dpg.set_value("cav.E_X0", d.cavity.E_X0)
    dpg.set_value("cav.F_P", d.cavity.F_P)
    dpg.set_value("cav.G", d.cavity.G)
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
    d = collect_design()
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
    dpg.fit_axis_data("xax1"); dpg.fit_axis_data("yax1")
    dpg.fit_axis_data("xax2"); dpg.fit_axis_data("yax2")

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
        with dpg.group(horizontal=True):
            dpg.add_input_text(tag="design.name", default_value="staged-device",
                               width=160)
            dpg.add_button(label="Load", callback=load_design)
            dpg.add_button(label="Save", callback=save_design)
            dpg.add_button(label="Export bundle", callback=export_bundle)
            dpg.add_button(label="  RUN  ", callback=run_device)
            dpg.add_text("tag chain [A] - every result inherits unmeasured inputs",
                         color=RED)
        dpg.add_text("", tag="status", color=AMBER)

        with dpg.group(horizontal=True):
            # ---------------- left: block diagram
            with dpg.child_window(width=760, height=640):
                with dpg.node_editor(width=-1, height=620, tag="editor",
                                     minimap=False):
                    with dpg.node(label="ELECTRICAL DRIVE", pos=(10, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            dpg.add_input_float(label="V", tag="drive.V", width=90,
                                                default_value=1.9)
                            dpg.add_input_float(label="I (uA)", tag="drive.I_uA",
                                                width=90, default_value=10.0)
                            dpg.add_input_float(label="duty", tag="drive.duty",
                                                width=90, default_value=1.0)
                            dpg.add_input_float(label="mu/pulse", tag="drive.mu",
                                                width=90, default_value=0.5)
                            dpg.add_input_float(label="b_e", tag="drive.b_e",
                                                width=90, default_value=0.02)
                            dpg.add_input_float(label="b_e exp m", tag="drive.b_e_m",
                                                width=90, default_value=1.5)
                            dpg.add_input_float(label="b_e Eact", tag="drive.b_e_Eact",
                                                width=90, default_value=100.0)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="drive_out"):
                            dpg.add_text("I, V, heat")

                    with dpg.node(label="MESA / THERMAL", pos=(10, 330)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="th_in"):
                            dpg.add_text("P = duty * I * V")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            dpg.add_input_float(label="mesa d (um)", tag="th.mesa",
                                                width=90, default_value=1.0,
                                                callback=lambda s, v: draw_cross_section())
                            dpg.add_input_float(label="T_hs (K)", tag="th.T_hs",
                                                width=90, default_value=77.0,
                                                callback=lambda s, v: draw_cross_section())
                            dpg.add_text("stack: edit in Fab stack panel ->")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="th_out"):
                            dpg.add_text("T_j")

                    with dpg.node(label="QD EMITTER", pos=(270, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="dot_in"):
                            dpg.add_text("carriers @ T_j")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            dpg.add_input_float(label="Delta_XX (meV)",
                                                tag="dot.delta_xx", width=90,
                                                default_value=3.5)
                            dpg.add_input_float(label="Gamma scale",
                                                tag="dot.gamma_scale", width=90,
                                                default_value=1.0)
                            dpg.add_input_float(label="r_XX", tag="dot.r_xx",
                                                width=90, default_value=0.72)
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
                            dpg.add_checkbox(label="enabled", tag="cav.enabled",
                                             default_value=False)
                            dpg.add_input_float(label="kappa (meV)", tag="cav.kappa",
                                                width=90, default_value=1.0)
                            dpg.add_input_float(label="track T (K)", tag="cav.T_track",
                                                width=90, default_value=120.0)
                            dpg.add_input_float(label="E_X0 (eV)", tag="cav.E_X0",
                                                width=90, default_value=1.88)
                            dpg.add_input_float(label="F_P", tag="cav.F_P",
                                                width=90, default_value=10.0)
                            dpg.add_input_float(label="gain G", tag="cav.G",
                                                width=90, default_value=8.0)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="cav_out"):
                            dpg.add_text("filtered + boosted")

                    with dpg.node(label="SLIT FILTER", pos=(530, 20)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Input,
                                                tag="fil_in"):
                            dpg.add_text("spectrum")
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            dpg.add_checkbox(label="enabled", tag="fil.enabled",
                                             default_value=True)
                            dpg.add_checkbox(label="auto w = Gamma(T_j)",
                                             tag="fil.auto_w", default_value=True)
                            dpg.add_input_float(label="w (meV)", tag="fil.w",
                                                width=90, default_value=2.0)
                            dpg.add_input_float(label="dx (meV)", tag="fil.dx",
                                                width=90, default_value=0.0)
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Output,
                                                tag="fil_out"):
                            dpg.add_text("to detector")

                    with dpg.node(label="APERTURE / ENSEMBLE (F5)", pos=(530, 330)):
                        with dpg.node_attribute(attribute_type=dpg.mvNode_Attr_Static):
                            dpg.add_input_float(label="log10 density", tag="ap.log_density",
                                                width=90, default_value=8.845)
                            dpg.add_input_float(label="aperture (um)", tag="ap.diam",
                                                width=90, default_value=1.0)
                            dpg.add_input_float(label="sigma_inh", tag="ap.sigma",
                                                width=90, default_value=40.0)
                            dpg.add_input_float(label="comp. bright r", tag="ap.r",
                                                width=90, default_value=0.3)

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
                with dpg.plot(height=270, width=-1):
                    dpg.add_plot_legend()
                    dpg.add_plot_axis(dpg.mvXAxis, label="heatsink T (K)", tag="xax1")
                    with dpg.plot_axis(dpg.mvYAxis, label="g2 / fractions", tag="yax1"):
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
                        dpg.add_line_series([], [], label="dT_J", tag="s_tj")
                        dpg.add_line_series([], [], label="Gamma(T_j)", tag="s_gam")
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
