"""fsim-gui v1 (Phase V): Streamlit thin client over fsim-core.

Three-layer rule: this file computes NO physics -- it edits cards, calls
fsim_core, renders results. Every session reads/writes the same YAML cards the
CLI uses, so a live demo is reproducible by anyone holding the card file.

Run:  streamlit run fsim_gui/app.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import Tag, load_card, widest
from fsim_core.fitting import V_A_TOL, fit_phase0
from fsim_core.integrator import g2_from, g2_of_T, oat_sensitivity, solve_Tc
from fsim_core.loading import f1b_g2, f8_g2_load, f8b_thin_fano, gamma_eff, loading_probs
from fsim_core.spectral import epsilon, gamma_of_T, lorentzian
from fsim_theme import load_tokens, register_plotly
from fsim_viz.figures import phase0_bundle

# Theme (studio-02): one token source. The `fsim` Plotly template (dark) is the
# default for every figure here, and every st.plotly_chart passes theme=None so
# Streamlit does not restyle it on top.
register_plotly("dark", set_default=True)
TOK = load_tokens()
C = TOK["color"]["dark"]
S = C["series"]          # categorical slots 1-8 (S[0] = slot 1)
DIV = C["div"]           # diverging poles: blue (raises) / neutral / red (lowers)


def _rgba(hexv: str, alpha: float) -> str:
    h = hexv.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{alpha})"


# Provenance is graded by LINE FORM, not hue (DIRECTION.md): [V] solid chip,
# [DR] solid + inset rule, [E] outline, [A] dashed outline. The names stay
# defined for anything importing them, but both now carry one ink tone /
# form glyphs instead of the old green/amber/red traffic light.
TAG_HEX = {Tag.V: C["ink-1"], Tag.DR: C["ink-1"], Tag.E: C["ink-2"], Tag.A: C["ink-2"]}
TAG_DOT = {Tag.V: "■", Tag.DR: "▣", Tag.E: "□", Tag.A: "⬚"}
_CHIP_CLASS = {Tag.V: "v", Tag.DR: "dr", Tag.E: "e", Tag.A: "a"}

st.set_page_config(page_title="FSIM", layout="wide")

_MONO = "'B612 Mono', Consolas, ui-monospace, monospace"
_LABEL = "'Barlow Semi Condensed', Barlow, 'Segoe UI', sans-serif"
_NUM = "Barlow, 'Segoe UI', sans-serif"  # readouts: tabular Barlow (B612 Mono's decimal point is a full cell)
st.markdown(f"""<style>
@font-face {{font-family:'B612 Mono';font-weight:400;src:url('app/static/fonts/b612-mono-latin-400-normal.woff2') format('woff2');}}
@font-face {{font-family:'Barlow Semi Condensed';font-weight:600;src:url('app/static/fonts/barlow-semi-condensed-latin-600-normal.woff2') format('woff2');}}
.fsim-chip {{display:inline-block;font-family:{_MONO};font-size:0.74em;line-height:1.45;
  padding:1px 7px;margin:0 4px 2px 0;border-radius:2px;border:1px solid {C['ink-1']};
  color:{C['ink-1']};background:transparent;letter-spacing:0.02em;white-space:nowrap;}}
.fsim-chip-v {{background:{C['ink-1']};color:{C['ground']};}}
.fsim-chip-dr {{background:{C['ink-1']};color:{C['ground']};
  box-shadow:inset 0 0 0 2px {C['ink-1']}, inset 0 0 0 3px {C['ground']};}}
.fsim-chip-e {{border-color:{C['ink-2']};}}
.fsim-chip-a {{border:1px dashed {C['ink-2']};}}
.fsim-engraved {{font-family:{_LABEL};font-weight:600;text-transform:uppercase;
  letter-spacing:{TOK['font']['label_tracking']};font-size:0.78rem;color:{C['ink-2']};
  border-bottom:1px solid {C['rim']};padding-bottom:3px;margin:4px 0 8px 0;}}
h1, h2, h3 {{font-family:{_LABEL} !important;text-transform:uppercase;
  letter-spacing:{TOK['font']['label_tracking']};}}
[data-testid="stTab"] p {{font-family:{_LABEL};font-weight:600;text-transform:uppercase;
  letter-spacing:{TOK['font']['label_tracking']};}}
[data-testid="stMetricValue"], [data-testid="stMetricValue"] * {{font-family:{_NUM};font-variant-numeric:tabular-nums;font-weight:500;
  font-variant-numeric:tabular-nums;}}
[data-testid="stMetricLabel"] p {{font-family:{_LABEL};text-transform:uppercase;
  letter-spacing:{TOK['font']['label_tracking']};color:{C['ink-2']};}}
.fsim-readout {{font-family:{_NUM};font-variant-numeric:tabular-nums;}}
/* F8: theme primaryColor is neutral ink-2 (.streamlit/config.toml); the beam
   red marks only the primary action button. */
button[kind="primary"], [data-testid="stBaseButton-primary"] {{background-color:{C['beam']} !important;
  border-color:{C['beam']} !important;color:#ffffff !important;}}
button[kind="primary"]:hover, [data-testid="stBaseButton-primary"]:hover {{filter:brightness(1.08);}}
[data-testid="stSliderThumbValue"], [data-testid="stSliderTickBarMin"], [data-testid="stSliderTickBarMax"] {{
  font-family:{_NUM};font-variant-numeric:tabular-nums;color:{C['ink-1']} !important;}}
[data-testid="stTab"][aria-selected="true"] p {{color:{C['ink-1']};}}
[data-testid="stCheckbox"] label[data-selected="true"] svg polyline {{stroke:{C['ground']} !important;stroke-width:2px;}}
</style>""", unsafe_allow_html=True)


def engraved(text: str, where=st):
    """Laser-engraved section label (Barlow Semi Condensed caps, tracked)."""
    where.markdown(f"<div class='fsim-engraved'>{text}</div>", unsafe_allow_html=True)


# ------------------------------------------------------------------ card + params

def tag_badge(tag: Tag, label=""):
    return f"<span class='fsim-chip fsim-chip-{_CHIP_CLASS[tag]}'>{tag.label} {label}</span>"


@st.cache_data
def _card_raw(path_str: str, mtime: float):
    return yaml.safe_load(Path(path_str).read_text(encoding="utf-8"))


card_names = sorted(p.name for p in (ROOT / "cards").glob("*.yaml"))
with st.sidebar:
    st.title("FSIM")
    card_path = ROOT / "cards" / st.selectbox("parameter card", card_names)
    card = load_card(card_path)
    st.caption(card.meta.get("device", ""))
    ph = card.placeholders()
    if ph:
        st.error(f"PLACEHOLDER entries: {', '.join(ph)} — no gate can close on this card.")
    engraved("provenance")
    st.markdown(
        f"{tag_badge(Tag.V, 'measured')} {tag_badge(Tag.DR, 'derived')} "
        f"{tag_badge(Tag.E, 'estimated')} {tag_badge(Tag.A, 'assumed')}",
        unsafe_allow_html=True)
    st.caption("Every displayed number carries the widest tag in its input chain.")

has_va = card.name == "chatzarakis" and "g2_vs_T" in card.datasets

# fitted parameterization: session fit > saved bundle > range midpoints (exploration)
fit = st.session_state.get(f"fit_{card.name}")
fit_json = ROOT / "out" / "phase0" / "fit_params.json"
if fit is not None:
    p_model, p_src = dict(fit.params), "session fit"
elif has_va and fit_json.exists():
    p_model = json.loads(fit_json.read_text())["params"]
    p_src = f"saved bundle ({fit_json.relative_to(ROOT)})"
else:
    p_model = {}
    for n, prm in card.params.items():
        p_model[n] = prm.value if prm.value is not None else 0.5 * sum(prm.range)
    p_src = "range midpoints — EXPLORATION ONLY (ranges are swept, never averaged)"

if has_va:
    rows = sorted(card.datasets["g2_vs_T"].rows, key=lambda r: r["T"])
    Ts_d = np.array([r["T"] for r in rows])
    ws_d = np.array([r["w"] for r in rows])
    dxs_d = np.array([r["dx"] for r in rows])
    w_of_T = lambda T: float(np.interp(T, Ts_d, ws_d))
    dx_of_T = lambda T: float(np.interp(T, Ts_d, dxs_d))

spec_tab, casc_tab, dash_tab, card_tab = st.tabs(
    ["Spectral explainer", "Cascade", "Dashboard", "Card editor"])


# ------------------------------------------------------------- 1. spectral panel

with spec_tab:
    c1, c2 = st.columns([3, 1])
    with c2:
        T = st.slider("temperature T (K)", 4.0, 320.0, 78.0, 1.0)
        use_pub = has_va and st.checkbox("published window at this T", value=has_va)
        if use_pub:
            w, dx = w_of_T(T), dx_of_T(T)
            st.caption(f"w = {w:.2f} meV, dx = {dx:.2f} meV (interpolated from card)")
        else:
            w = st.slider("filter width w (meV)", 0.1, 15.0, 2.0, 0.1)
            dx = st.slider("X offset from window center dx (meV)", -6.0, 6.0, 0.0, 0.1)
        mu = st.slider("mean loading μ (F1b)", 0.0, 5.0,
                       float(p_model.get("mu", card.params.get("mu_op").value
                             if "mu_op" in card.params else 0.0) or 0.0), 0.05)
        drive_mode = st.radio("drive mode", ["PL", "EL"], horizontal=True,
                              key="spec_drive_mode")
        dg_inj, I_ratio = 0.0, 1.0
        p_inj = float(p_model.get("p_inj", 1.0) or 1.0)
        if drive_mode == "EL":
            dg_inj = st.number_input("dGamma_inj (meV)", 0.0, 20.0, 2.0, 0.5)
            I_ratio = st.number_input("I / I_ref", 0.0, 20.0, 1.0, 0.5)

    delta = p_model.get("delta_xx", 5.9)
    gam_pl = float(gamma_of_T(T, p_model.get("gamma0", 0.5), p_model.get("a_ac", 2e-3),
                              p_model.get("b_lo", 25.0), p_model.get("E_lo", 36.6)))
    # F5' injection broadening: THE ONE implementation lives in
    # fsim_core.loading.gamma_eff (also used by integrator.g2_of_T and
    # device.evaluate) -- the GUI never re-derives that broadening expression
    # itself (three-layer rule).
    gam_el = float(gamma_eff(gam_pl, dg_inj, I_ratio, p_inj))
    gam = gam_el if drive_mode == "EL" else gam_pl
    spec = epsilon(delta, gam, gam, w=w, dx=dx)
    g2_dot = float(f1b_g2(mu, spec.eps)) if mu > 0 else spec.eps

    x = np.linspace(-delta - 6 * gam - 4, 6 * gam + 4, 1200)
    LX = lorentzian(x, 0.0, gam)
    LXX = lorentzian(x, -delta, gam)
    peak = LX.max()
    cen = -dx  # window center relative to X
    fig = go.Figure()
    # filter window: neutral ref-wash with a direct label (not a series hue)
    fig.add_vrect(x0=cen - w / 2, x1=cen + w / 2, fillcolor=C["ref-wash"], opacity=1.0,
                  line_width=0, layer="below")
    fig.add_annotation(x=cen, y=1.0, yref="paper", yanchor="bottom", showarrow=False,
                       text=f"filter w = {w:.1f} meV", font=dict(color=C["ink-2"], size=12))
    fig.add_trace(go.Scatter(x=x, y=LX / peak, name="X", line=dict(color=S[0], width=2)))
    fig.add_trace(go.Scatter(x=x, y=LXX / peak, name="XX", line=dict(color=S[1], width=2)))
    inw = (x >= cen - w / 2) & (x <= cen + w / 2)
    fig.add_trace(go.Scatter(x=x[inw], y=(LXX / peak)[inw], fill="tozeroy", mode="none",
                             fillcolor=_rgba(S[1], 0.25),
                             name="t_XX leakage"))
    # direct labels at the peaks and at the leakage (the reported quantity)
    fig.add_annotation(x=0.0, y=1.0, text="X", showarrow=False, yshift=12,
                       font=dict(color=S[0], size=13))
    fig.add_annotation(x=-delta, y=float(LXX.max() / peak), text="XX", showarrow=False,
                       yshift=12, font=dict(color=S[1], size=13))
    if inw.any():
        k_leak = int(np.argmax(np.where(inw, LXX, -np.inf)))
        fig.add_annotation(x=float(x[k_leak]), y=float(LXX[k_leak] / peak),
                           text=f"ε = {spec.eps:.4f}", showarrow=True, arrowhead=0,
                           arrowcolor=C["ref"], ax=40, ay=-30,
                           font=dict(color=C["ink-1"], size=12))
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=40, b=10),
                      xaxis_title="energy − E_X (meV)", yaxis_title="peak-normalized intensity",
                      legend=dict(orientation="h", x=0, y=1.02, xanchor="left",
                                  yanchor="bottom"))
    with c1:
        engraved("spectral window: X / XX through the filter")
    c1.plotly_chart(fig, use_container_width=True, theme=None)

    ratio = (gam_el / gam_pl) if gam_pl > 0 else float("nan")
    ratio_style = (f"color:{C['status']['warn']};font-weight:bold" if ratio > 2
                   else "font-weight:bold")
    c2.markdown(
        f"PL Γ: {gam_pl:.2f} meV / EL Γ_eff: {gam_el:.2f} meV, ratio "
        f"<span style='{ratio_style}'>{ratio:.2f}x</span>"
        + ("  ← Kitamura ~4x" if ratio > 2 else ""),
        unsafe_allow_html=True)
    engraved("readouts", c2)

    def _ro(v):  # numeric readout in tabular Barlow
        return f"<span class='fsim-readout'>{v}</span>"

    c2.markdown(f"**Γ(T)** = {_ro(f'{gam:.2f}')} meV  (active: {drive_mode})",
                unsafe_allow_html=True)
    c2.markdown(f"**ε = t_XX/t_X** = {_ro(f'{spec.eps:.4f}')}", unsafe_allow_html=True)
    c2.markdown(f"**g²₀(μ)** = {_ro(f'{g2_dot:.4f}')}" + ("  (F1: ε)" if mu == 0 else "  (F1b)"),
                unsafe_allow_html=True)
    c2.markdown(f"**t_X (brightness)** = {_ro(f'{spec.t_x:.3f}')}", unsafe_allow_html=True)
    spec_tag = widest(card.params["delta_xx"].tag if "delta_xx" in card.params else Tag.A,
                      Tag.A)  # linewidth params are fitted [A]
    c2.markdown(tag_badge(spec_tag, "spectral chain"), unsafe_allow_html=True)
    c2.caption(f"params: {p_src}")


# ------------------------------------------------------------ 2. cascade diagram

with casc_tab:
    Tc2 = st.slider("temperature for rates (K)", 4.0, 320.0, 150.0, 1.0, key="Tcasc")
    pt = None
    if all(k in p_model for k in ("a_esc", "E_a", "b_p", "E_b", "b0", "beta")):
        p_full = {**p_model}
        if has_va:
            p_full["w"], p_full["dx"] = w_of_T(Tc2), dx_of_T(Tc2)
        p_full.setdefault("w", 2.0)
        pt = g2_of_T(Tc2, p_full)
    P0, P1, P2 = loading_probs(p_model.get("mu", 0.33) or 0.33)

    fig = go.Figure()
    for y, name in ((2.0, "|XX⟩"), (1.0, "|X⟩"), (0.0, "|0⟩")):
        fig.add_shape(type="line", x0=0.25, x1=0.75, y0=y, y1=y, line=dict(width=4, color=C["ink-1"]))
        fig.add_annotation(x=0.20, y=y, text=name, showarrow=False,
                           font=dict(size=18, color=C["ink-1"]))
    fig.add_annotation(x=0.5, y=1.5, ax=0.5, ay=2.0, axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowwidth=2, arrowcolor=S[1])
    fig.add_annotation(x=0.62, y=1.55, text=f"XX photon → filter: t_XX"
                       + (f" = {pt.eps * pt.t_x:.3f}" if pt else ""), showarrow=False)
    fig.add_annotation(x=0.5, y=0.5, ax=0.5, ay=1.0, axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowwidth=2, arrowcolor=S[0])
    fig.add_annotation(x=0.62, y=0.55, text=f"X photon → filter: t_X"
                       + (f" = {pt.t_x:.3f}" if pt else ""), showarrow=False)
    fig.add_annotation(x=0.32, y=2.35, text=f"cap-2 loading: P₁={P1:.2f}, P₂={P2:.2f}",
                       showarrow=False, font=dict(color=C["muted"]))
    if pt:
        fig.add_annotation(x=0.85, y=1.0,
                           text=f"background channels →<br>ρ({Tc2:.0f} K) = {pt.rho:.3f}",
                           showarrow=False, font=dict(color=C["ink-2"]))
        fig.add_annotation(x=0.5, y=-0.45,
                           text=(f"g²(0) = 1 − ρ²(1−g²₀) = {pt.g2:.3f}   "
                                 f"[ε = {pt.eps:.3f}, Γ = {pt.gamma:.2f} meV]"),
                           showarrow=False, font=dict(size=16, color=C["ink-1"]))
    fig.update_layout(height=480, xaxis=dict(visible=False, range=[0, 1.1]),
                      yaxis=dict(visible=False, range=[-0.7, 2.6]),
                      margin=dict(l=10, r=10, t=10, b=10), hovermode=False)
    engraved("cascade: XX → X → 0 through the filter")
    st.plotly_chart(fig, use_container_width=True, theme=None)
    st.caption("Redrawn from card values; rates and fractions are computed by fsim-core, "
               "never by the GUI (three-layer rule).")


# ---------------------------------------------------------------- 3. dashboard

with dash_tab:
    if not has_va:
        st.info("Dashboard requires a card with a g2_vs_T validation dataset "
                "(select chatzarakis.yaml).")
    else:
        cA, cB = st.columns([1, 3])
        with cA:
            n_starts = st.select_slider("fit multistarts", [8, 16, 24, 40], value=16)
            if st.button("Run V-a fit from this card", type="primary"):
                with st.spinner("fitting (ranges swept, anchors enforced)..."):
                    st.session_state[f"fit_{card.name}"] = fit_phase0(
                        card, n_starts=n_starts)
                st.rerun()
            if fit is not None and st.button("Write report bundle (out/phase0)"):
                phase0_bundle(fit, card_path, ROOT / "out" / "phase0")
                st.success("bundle written: figure + CSVs + params + card snapshot")

            st.markdown("---")
            F_p_dash = st.slider("F_p (pump Fano)", 0.0, 1.5, 1.0, 0.05, key="dash_Fp")
            mu_dash = st.slider("mu (loading, F8)", 0.0, 3.0,
                                float(p_model.get("mu", 0.5) or 0.5), 0.05, key="dash_mu")
            # F8 domain floor: mu >= 1 - F_eff, F_eff = f8b_thin_fano(eta, F_p) --
            # computed via the core (no eta_capture control here, so eta=1: no
            # dot-capture thinning, F_eff = F_p).
            F_eff_dash = f8b_thin_fano(1.0, F_p_dash)
            floor_dash = 1.0 - F_eff_dash
            if mu_dash < floor_dash:
                st.warning(f"F8 domain: mu {mu_dash:.2f} below floor "
                          f"{floor_dash:.2f} (= 1 - F_p) — g2_load undefined")
                st.caption("g2_load = —")
            elif mu_dash <= 0:
                st.caption("g2_load = —")
            else:
                try:
                    g2load_dash = float(f8_g2_load(mu_dash, F_eff_dash))
                    st.metric("g2_load", f"{g2load_dash:.3f}")
                except ValueError:
                    st.caption("g2_load = —")

        if not p_model or "a_esc" not in p_model:
            st.warning("No fitted parameterization available — run the fit.")
        else:
            p = dict(p_model)
            p["w"], p["dx"] = w_of_T, dx_of_T
            Ts = np.linspace(Ts_d.min() - 18, Ts_d.max() + 60, 250)
            pts = [g2_of_T(t, p) for t in Ts]
            Tc = solve_Tc(p)

            with cA:
                st.metric("master ceiling T_c", "—" if np.isnan(Tc) else f"{Tc:.1f} K")
                if fit is not None:
                    st.metric("max |residual|", f"{np.max(np.abs(fit.residuals)):.3f}",
                              delta=f"tol ±{V_A_TOL}", delta_color="off")
                    st.markdown(tag_badge(fit.tag, "fit chain"), unsafe_allow_html=True)
                st.caption(f"params source: {p_src}")

            # two stacked panels, shared x: g2 (model + data + refs) over the
            # eps / rho^2 decomposition on 0-1 (02-charts: no crowded single axis)
            fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.07,
                                row_heights=[0.62, 0.38])
            mark_ring = dict(width=TOK["mark"]["marker_ring"], color=C["surface"])
            is_ub = np.array([r.get("bound") == "upper" for r in rows])
            g2_d = np.array([r["g2"] for r in rows])
            err_d = np.array([r["err"] for r in rows])
            fig.add_trace(go.Scatter(
                x=Ts_d[~is_ub], y=g2_d[~is_ub],
                error_y=dict(type="data", array=err_d[~is_ub], color=C["ink-2"], thickness=1),
                mode="markers", name="data (Fig. 5) [V]",
                marker=dict(color=C["ink-1"], size=2 * TOK["mark"]["marker_radius"] + 1,
                            line=mark_ring)), row=1, col=1)
            if is_ub.any():
                fig.add_trace(go.Scatter(
                    x=Ts_d[is_ub], y=g2_d[is_ub], mode="markers+text", text=["≤"] * int(is_ub.sum()),
                    textposition="middle left", textfont=dict(color=C["ink-1"], size=13),
                    error_y=dict(type="data", array=err_d[is_ub], color=C["ink-2"], thickness=1),
                    name="data, upper bound (≤)",
                    marker=dict(color=C["ink-1"], symbol="triangle-down",
                                size=2 * TOK["mark"]["marker_radius"] + 3, line=mark_ring)),
                    row=1, col=1)
            fig.add_trace(go.Scatter(x=Ts, y=[q.g2 for q in pts], name="F-series model g²(0)",
                                     line=dict(color=S[0], width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(x=Ts, y=[q.eps for q in pts], name="ε(T)",
                                     line=dict(color=S[1], width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(x=Ts, y=[q.rho**2 for q in pts], name="ρ²(T)",
                                     line=dict(color=S[2], width=2)), row=2, col=1)
            # references (not verdicts): ref ink, 1px, direct labels
            fig.add_hline(y=0.5, line_width=TOK["mark"]["ref"], line_color=C["ref"],
                          row=1, col=1, annotation_text="g²(0) = 0.5 ceiling",
                          annotation_position="top right",
                          annotation_font=dict(color=C["ink-2"], size=12))
            if np.isfinite(Tc):
                fig.add_vline(x=Tc, line_width=TOK["mark"]["ref"], line_color=C["ref"],
                              row="all", col=1)
                fig.add_annotation(x=Tc, y=1.0, xref="x", yref="paper", yanchor="bottom",
                                   showarrow=False, text=f"T_c = {Tc:.0f} K",
                                   font=dict(color=C["ink-1"], size=12))
            fig.update_yaxes(range=[0, 1], title_text="g²(0)", row=1, col=1)
            fig.update_yaxes(range=[0, 1], title_text="ε, ρ²", row=2, col=1)
            fig.update_xaxes(title_text="T (K)", row=2, col=1)
            fig.update_xaxes(showspikes=True, spikemode="across", spikethickness=1,
                             spikecolor=C["ref"], spikedash="solid")
            fig.update_layout(height=560, margin=dict(l=10, r=10, t=40, b=10),
                              legend=dict(orientation="h", x=0, y=1.06, xanchor="left",
                                          yanchor="bottom"))
            with cB:
                engraved("g²(T): model vs Chatzarakis Fig. 5, and its decomposition")
            cB.plotly_chart(fig, use_container_width=True, theme=None)

            if fit is not None:
                st.dataframe(
                    [{"T (K)": r["T"],
                      "g2 data": ("≤" if r.get("bound") == "upper" else "") + f"{r['g2']:.3f}",
                      "g2 model": f"{m:.3f}", "residual": f"{res:+.3f}",
                      "within ±0.03": "ok" if abs(res) <= V_A_TOL else "FAIL"}
                     for r, m, res in zip(fit.data, fit.model_g2, fit.residuals)],
                    use_container_width=True, hide_index=True)

            sens = oat_sensitivity(p)
            # diverging bars sorted by |value| (largest on top), value at each tip;
            # blue pole raises T_c, red pole lowers it
            items = sorted(sens.items(), key=lambda kv: abs(kv[1]))
            figt = go.Figure(go.Bar(
                x=[v for _, v in items], y=[k for k, _ in items], orientation="h",
                marker=dict(color=[DIV[2] if v < 0 else DIV[0] for _, v in items],
                            cornerradius=4),
                text=[f"{v:+.1f} K" for _, v in items], textposition="outside",
                textfont=dict(family="Barlow, Segoe UI, sans-serif", color=C["ink-2"],
                              size=12),
                cliponaxis=False, hovertemplate="%{y}: ΔT_c %{x:+.2f} K<extra></extra>"))
            xmax = max([abs(v) for _, v in items] + [1e-9]) * 1.25
            figt.update_xaxes(range=[-xmax, xmax], zeroline=True, zerolinecolor=C["axis"],
                              zerolinewidth=1)
            figt.update_layout(height=max(260, 34 * len(items) + 90), bargap=0.45,
                               title="T_c sensitivity: ΔT_c for +5% of each parameter",
                               xaxis_title="ΔT_c (K)", margin=dict(l=10, r=10, t=40, b=10),
                               hovermode="closest")
            engraved("sensitivity")
            st.plotly_chart(figt, use_container_width=True, theme=None)


# --------------------------------------------------------------- 4. card editor

with card_tab:
    raw = _card_raw(str(card_path), card_path.stat().st_mtime)
    st.markdown(f"**{card.name}** — {card.meta.get('source', '')}")
    table = []
    for n, prm in card.params.items():
        table.append({
            "tag": f"{TAG_DOT[prm.tag]} {prm.tag.label}", "param": n,
            "value": prm.value, "lo": prm.range[0] if prm.range else None,
            "hi": prm.range[1] if prm.range else None,
            "unit": prm.unit, "source": prm.source,
        })
    edited = st.data_editor(
        table, use_container_width=True, hide_index=True,
        disabled=["tag", "param", "unit", "source"], key=f"editor_{card.name}")

    target = st.text_input("save card as (cards/…)", value=f"{card.name}-edited.yaml")
    if st.button("Save card"):
        out = dict(raw)
        for row in edited:
            n = row["param"]
            entry = dict(out["params"][n])
            if row["value"] is not None:
                entry["value"] = float(row["value"])
                entry.pop("range", None)
            elif row["lo"] is not None and row["hi"] is not None:
                entry["range"] = [float(row["lo"]), float(row["hi"])]
                entry.pop("value", None)
            out["params"][n] = entry
        dest = ROOT / "cards" / target
        dest.write_text(yaml.safe_dump(out, sort_keys=False, allow_unicode=True),
                        encoding="utf-8")
        st.success(f"saved {dest.relative_to(ROOT)} — same schema the CLI uses; "
                   "reproducible by anyone holding this file")
    for name, ds in card.datasets.items():
        st.markdown(f"**dataset `{name}`** {tag_badge(ds.tag)}", unsafe_allow_html=True)
        st.dataframe(ds.rows, use_container_width=True, hide_index=True)
