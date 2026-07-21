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
import streamlit as st
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import Tag, load_card, widest
from fsim_core.fitting import V_A_TOL, fit_phase0
from fsim_core.integrator import g2_from, g2_of_T, oat_sensitivity, solve_Tc
from fsim_core.loading import f1b_g2, loading_probs
from fsim_core.spectral import epsilon, gamma_of_T, lorentzian
from fsim_viz.figures import phase0_bundle

TAG_HEX = {Tag.V: "#1a9641", Tag.DR: "#e3a21a", Tag.E: "#e3a21a", Tag.A: "#d7191c"}
TAG_DOT = {Tag.V: "🟢", Tag.DR: "🟡", Tag.E: "🟡", Tag.A: "🔴"}

st.set_page_config(page_title="FSIM", layout="wide")


# ------------------------------------------------------------------ card + params

def tag_badge(tag: Tag, label=""):
    return (f"<span style='background:{TAG_HEX[tag]};color:white;padding:1px 8px;"
            f"border-radius:8px;font-size:0.8em'>{tag.label} {label}</span>")


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
    st.markdown(
        f"{tag_badge(Tag.V, 'measured')} {tag_badge(Tag.DR, 'derived/estimated')} "
        f"{tag_badge(Tag.A, 'assumed')}", unsafe_allow_html=True)
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

    delta = p_model.get("delta_xx", 5.9)
    gam = float(gamma_of_T(T, p_model.get("gamma0", 0.5), p_model.get("a_ac", 2e-3),
                           p_model.get("b_lo", 25.0), p_model.get("E_lo", 36.6)))
    spec = epsilon(delta, gam, gam, w=w, dx=dx)
    g2_dot = float(f1b_g2(mu, spec.eps)) if mu > 0 else spec.eps

    x = np.linspace(-delta - 6 * gam - 4, 6 * gam + 4, 1200)
    LX = lorentzian(x, 0.0, gam)
    LXX = lorentzian(x, -delta, gam)
    peak = LX.max()
    cen = -dx  # window center relative to X
    fig = go.Figure()
    fig.add_vrect(x0=cen - w / 2, x1=cen + w / 2, fillcolor="#bbbbbb", opacity=0.25,
                  line_width=0, annotation_text="filter w", annotation_position="top")
    fig.add_trace(go.Scatter(x=x, y=LX / peak, name="X", line=dict(color="#2166ac", width=3)))
    fig.add_trace(go.Scatter(x=x, y=LXX / peak, name="XX", line=dict(color="#5e3c99", width=3)))
    inw = (x >= cen - w / 2) & (x <= cen + w / 2)
    fig.add_trace(go.Scatter(x=x[inw], y=(LXX / peak)[inw], fill="tozeroy", mode="none",
                             fillcolor="rgba(215,25,28,0.55)",
                             name="t_XX leakage"))
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title="energy − E_X (meV)", yaxis_title="peak-normalized intensity",
                      legend=dict(orientation="h"))
    c1.plotly_chart(fig, use_container_width=True)

    c2.markdown(f"**Γ(T)** = {gam:.2f} meV")
    c2.markdown(f"**ε = t_XX/t_X** = {spec.eps:.4f}")
    c2.markdown(f"**g²₀(μ)** = {g2_dot:.4f}" + ("  (F1: ε)" if mu == 0 else "  (F1b)"))
    c2.markdown(f"**t_X (brightness)** = {spec.t_x:.3f}")
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
        fig.add_shape(type="line", x0=0.25, x1=0.75, y0=y, y1=y, line=dict(width=4))
        fig.add_annotation(x=0.20, y=y, text=name, showarrow=False, font=dict(size=18))
    fig.add_annotation(x=0.5, y=1.5, ax=0.5, ay=2.0, axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowwidth=2, arrowcolor="#5e3c99")
    fig.add_annotation(x=0.62, y=1.55, text=f"XX photon → filter: t_XX"
                       + (f" = {pt.eps * pt.t_x:.3f}" if pt else ""), showarrow=False)
    fig.add_annotation(x=0.5, y=0.5, ax=0.5, ay=1.0, axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowwidth=2, arrowcolor="#2166ac")
    fig.add_annotation(x=0.62, y=0.55, text=f"X photon → filter: t_X"
                       + (f" = {pt.t_x:.3f}" if pt else ""), showarrow=False)
    fig.add_annotation(x=0.32, y=2.35, text=f"cap-2 loading: P₁={P1:.2f}, P₂={P2:.2f}",
                       showarrow=False, font=dict(color="#666666"))
    if pt:
        fig.add_annotation(x=0.85, y=1.0,
                           text=f"background channels →<br>ρ({Tc2:.0f} K) = {pt.rho:.3f}",
                           showarrow=False, font=dict(color="#d7191c"))
        fig.add_annotation(x=0.5, y=-0.45,
                           text=(f"g²(0) = 1 − ρ²(1−g²₀) = {pt.g2:.3f}   "
                                 f"[ε = {pt.eps:.3f}, Γ = {pt.gamma:.2f} meV]"),
                           showarrow=False, font=dict(size=16))
    fig.update_layout(height=480, xaxis=dict(visible=False, range=[0, 1.1]),
                      yaxis=dict(visible=False, range=[-0.7, 2.6]),
                      margin=dict(l=10, r=10, t=10, b=10))
    st.plotly_chart(fig, use_container_width=True)
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

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=Ts_d, y=[r["g2"] for r in rows],
                error_y=dict(type="data", array=[r["err"] for r in rows]),
                mode="markers", name="data (Fig. 5)", marker=dict(color="#333333", size=9)))
            fig.add_trace(go.Scatter(x=Ts, y=[q.g2 for q in pts], name="F-series model",
                                     line=dict(color="#2166ac", width=3)))
            fig.add_trace(go.Scatter(x=Ts, y=[q.eps for q in pts], name="ε(T)",
                                     line=dict(color="#5e3c99", dash="dot")))
            fig.add_trace(go.Scatter(x=Ts, y=[q.rho**2 for q in pts], name="ρ²(T)",
                                     line=dict(color="#e66101", dash="dot")))
            fig.add_hline(y=0.5, line_dash="dot", line_color="#888888")
            if np.isfinite(Tc):
                fig.add_vline(x=Tc, line_dash="dash", line_color="#d7191c",
                              annotation_text=f"T_c = {Tc:.0f} K")
            fig.update_layout(height=460, xaxis_title="T (K)", yaxis_title="g²(0)",
                              margin=dict(l=10, r=10, t=30, b=10))
            cB.plotly_chart(fig, use_container_width=True)

            if fit is not None:
                st.dataframe(
                    [{"T (K)": r["T"],
                      "g2 data": ("≤" if r.get("bound") == "upper" else "") + f"{r['g2']:.3f}",
                      "g2 model": f"{m:.3f}", "residual": f"{res:+.3f}",
                      "within ±0.03": "ok" if abs(res) <= V_A_TOL else "FAIL"}
                     for r, m, res in zip(fit.data, fit.model_g2, fit.residuals)],
                    use_container_width=True, hide_index=True)

            sens = oat_sensitivity(p)
            items = sorted(sens.items(), key=lambda kv: abs(kv[1]))
            figt = go.Figure(go.Bar(
                x=[v for _, v in items], y=[k for k, _ in items], orientation="h",
                marker_color=["#d7191c" if v < 0 else "#2166ac" for _, v in items]))
            figt.update_layout(height=380, title="T_c sensitivity: ΔT_c for +5% of each parameter",
                               xaxis_title="ΔT_c (K)", margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(figt, use_container_width=True)


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
        st.markdown(f"**dataset `{name}`** {TAG_DOT[ds.tag]} {ds.tag.label}")
        st.dataframe(ds.rows, use_container_width=True, hide_index=True)
