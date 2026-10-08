"""Dear PyGui theme binder driven by fsim_theme/tokens.json.

    import dearpygui.dearpygui as dpg
    from fsim_theme.dpg import bind_theme
    dpg.create_context()
    fonts = bind_theme(dpg, mode="dark")   # {"ui": id|None, "mono": id|None, "theme": id}

Safe before or after dpg.create_viewport (it only creates registry items and
binds them). A missing font file is skipped and its id comes back as None.
"""
from __future__ import annotations

from fsim_theme import FONT_DIR, load_tokens, _check_mode


def _rgba(hexv: str, alpha: float = 1.0) -> tuple[int, int, int, int]:
    h = hexv.lstrip("#")
    return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16), int(round(255 * alpha)))


def _font_file(tok: dict, family: str, weight: int) -> str | None:
    for face in tok["font"]["faces"]:
        if face["family"] == family and face["weight"] == weight and face.get("ttf"):
            p = FONT_DIR / face["ttf"]
            return str(p) if p.is_file() else None
    return None


def bind_theme(dpg, mode: str = "dark") -> dict:
    """Create + bind a global fsim theme and font registry. Returns item ids."""
    tok = load_tokens()
    c = tok["color"][_check_mode(mode)]
    rad = tok["radius"]
    sp = tok["space"]

    def col(name, value, cat=None):
        const = getattr(dpg, name, None)
        if const is None:
            return
        if cat is None:
            dpg.add_theme_color(const, value)
        else:
            dpg.add_theme_color(const, value, category=cat)

    def var(name, *vals):
        const = getattr(dpg, name, None)
        if const is not None:
            dpg.add_theme_style(const, *vals)

    def pvar(name, *vals):
        const = getattr(dpg, name, None)
        if const is not None:
            dpg.add_theme_style(const, *vals, category=dpg.mvThemeCat_Plots)

    ground, plate, rim = _rgba(c["ground"]), _rgba(c["plate"]), _rgba(c["rim"])
    surface, ink1, ink2 = _rgba(c["surface"]), _rgba(c["ink-1"]), _rgba(c["ink-2"])
    muted, grid, axis = _rgba(c["muted"]), _rgba(c["grid"]), _rgba(c["axis"])
    beam, focus = _rgba(c["beam"]), _rgba(c["focus"])
    clear = (0, 0, 0, 0)

    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvAll):
            col("mvThemeCol_WindowBg", ground)
            col("mvThemeCol_ChildBg", plate)
            col("mvThemeCol_PopupBg", plate)
            col("mvThemeCol_MenuBarBg", plate)
            col("mvThemeCol_Border", rim)
            col("mvThemeCol_BorderShadow", clear)
            col("mvThemeCol_Text", ink1)
            col("mvThemeCol_TextDisabled", muted)
            col("mvThemeCol_TitleBg", plate)
            col("mvThemeCol_TitleBgActive", plate)
            col("mvThemeCol_TitleBgCollapsed", plate)
            col("mvThemeCol_FrameBg", surface)
            col("mvThemeCol_FrameBgHovered", rim)
            col("mvThemeCol_FrameBgActive", rim)
            col("mvThemeCol_Button", plate)
            col("mvThemeCol_ButtonHovered", rim)
            col("mvThemeCol_ButtonActive", axis)
            col("mvThemeCol_Header", plate)
            col("mvThemeCol_HeaderHovered", rim)
            col("mvThemeCol_HeaderActive", axis)
            col("mvThemeCol_Tab", plate)
            col("mvThemeCol_TabHovered", rim)
            col("mvThemeCol_TabActive", rim)
            col("mvThemeCol_TabSelected", rim)
            col("mvThemeCol_SliderGrab", ink2)
            col("mvThemeCol_SliderGrabActive", beam)  # beam = the active control
            col("mvThemeCol_CheckMark", ink1)
            col("mvThemeCol_Separator", rim)
            col("mvThemeCol_ScrollbarBg", ground)
            col("mvThemeCol_ScrollbarGrab", rim)
            col("mvThemeCol_ScrollbarGrabHovered", axis)
            col("mvThemeCol_ScrollbarGrabActive", axis)
            col("mvThemeCol_TableHeaderBg", plate)
            col("mvThemeCol_TableBorderStrong", rim)
            col("mvThemeCol_TableBorderLight", grid)
            col("mvThemeCol_TableRowBg", clear)
            col("mvThemeCol_TableRowBgAlt", _rgba(c["grid"], 0.35))
            col("mvThemeCol_TextSelectedBg", _rgba(c["ref"], 0.25))
            col("mvThemeCol_NavHighlight", focus)
            col("mvThemeCol_NavCursor", focus)

            var("mvStyleVar_FrameRounding", rad)
            var("mvStyleVar_WindowRounding", rad)
            var("mvStyleVar_ChildRounding", rad)
            var("mvStyleVar_PopupRounding", rad)
            var("mvStyleVar_GrabRounding", rad)
            var("mvStyleVar_TabRounding", rad)
            var("mvStyleVar_ScrollbarRounding", rad)
            var("mvStyleVar_FrameBorderSize", 1)
            var("mvStyleVar_WindowPadding", sp * 1.5, sp * 1.5)
            var("mvStyleVar_FramePadding", sp, sp / 2)
            var("mvStyleVar_ItemSpacing", sp, sp * 0.75)
            var("mvStyleVar_ItemInnerSpacing", sp / 2, sp / 2)
            var("mvStyleVar_CellPadding", sp / 2, sp / 4)

            pc = dpg.mvThemeCat_Plots
            col("mvPlotCol_FrameBg", plate, pc)
            col("mvPlotCol_PlotBg", surface, pc)
            col("mvPlotCol_PlotBorder", rim, pc)
            col("mvPlotCol_LegendBg", _rgba(c["plate"], 0.85), pc)
            col("mvPlotCol_LegendBorder", clear, pc)
            col("mvPlotCol_LegendText", ink2, pc)
            col("mvPlotCol_TitleText", ink1, pc)
            col("mvPlotCol_InlayText", ink2, pc)
            col("mvPlotCol_AxisText", ink2, pc)
            col("mvPlotCol_AxisGrid", grid, pc)
            col("mvPlotCol_AxisTick", axis, pc)
            col("mvPlotCol_AxisBg", clear, pc)
            col("mvPlotCol_AxisBgHovered", _rgba(c["rim"], 0.5), pc)
            col("mvPlotCol_AxisBgActive", _rgba(c["rim"], 0.8), pc)
            col("mvPlotCol_Crosshairs", _rgba(c["ref"]), pc)
            col("mvPlotCol_Selection", _rgba(c["ref"], 0.6), pc)
            pvar("mvPlotStyleVar_LineWeight", float(tok["mark"]["line"]))
            pvar("mvPlotStyleVar_PlotPadding", sp, sp)

    fonts = {"ui": None, "mono": None}
    ui_file = _font_file(tok, tok["font"]["ui"], 400)
    mono_file = _font_file(tok, tok["font"]["num"], 500)  # readouts: Barlow Medium
    size = tok["size"]["dpg_font"]
    if ui_file or mono_file:
        try:
            with dpg.font_registry():
                if ui_file:
                    fonts["ui"] = dpg.add_font(ui_file, size)
                if mono_file:
                    fonts["mono"] = dpg.add_font(mono_file, size)
        except Exception:  # a bad font file must never kill the GUI
            fonts = {"ui": None, "mono": None}
    if fonts["ui"] is not None:
        dpg.bind_font(fonts["ui"])

    dpg.bind_theme(theme_id)
    fonts["theme"] = theme_id
    return fonts
