"""SceneSpec builders for the FSIM Studio 3D views (fsim_studio/web/js/viz3d).

Three-layer rule: every number in a SceneSpec comes from an fsim_core call or
an out/ file.  This module only loads cards, calls fsim_core, reads committed
CSVs, resamples arrays for transport and attaches labels/provenance.  The
browser builds meshes from the spec and never computes physics.

All builders are pure functions (no globals mutated, no files written) that
return JSON-safe dicts (NaN/inf -> None, numpy -> python).

SceneSpec schema "fsim.scene/1"
-------------------------------
Common keys (every kind):
  schema      "fsim.scene/1"
  kind        "device" | "band" | "cascade" | "surface" | "lattice"
  variant     kind-specific layout selector (below)
  card        card name (or campaign id for lattice)
  title       short human title
  catalog_id  "FS-<card>-<hash6>" (sha256 of the scene request)
  units       {length: "um"|"nm", energy: "eV", time: "ns", ...}
  axes        {x|y|z: {label, unit, scale, badge}}; scale is the display
              factor relative to true scale.  scale != 1 REQUIRES a badge
              string (e.g. "z x20"); the front-end scale.js throws otherwise.
  layers      [{name, t_nm, n, material, tag, role, source}] bottom -> top
              (device views); tag is one of V | DR | E | A; role is one of
              substrate | cladding | core | active | mirror | contact |
              oxide | doped_p | doped_n | intrinsic | thermal
  overlays    kind-specific overlay data (mode texture, cones, curves ...)
  scalars     {key: {value, unit, tag, label, note?}}; value may be None
              with note giving the reason
  labels      mandatory on-screen labels (strings), e.g. "static (non-headline)",
              "NIR, false colour", "etch removes dot layer; sidewalls not modelled"
  provenance  {calls: [fsim_core calls / out files], card_file, ...}
  bookmarks   optional [{key, name}] camera bookmark names (1-4)

Variants
--------
device / edge_ridge (units um; origin: ridge centre x=0, stack bottom y=0,
  front facet z=0, propagation toward -z into the device):
  stacks.optical  {label "optical stack (waveguide.py)", layers, total_nm}
  stacks.thermal  {label "thermal stack (thermal.py)", layers, total_nm}
  geometry        ridge_width_nm, etch_depth_nm, etch_floor_nm (height of the
                  etched floor above the stack bottom), shoulder_nm (the mode
                  solver's lateral pad), L_um, front_len_um (true-scale drawn
                  length before the break), R_back, dot {height_nm, radius_nm,
                  y_nm (dot-plane centre height)}, dots {density_cm2,
                  per_um2, seed, region_um [w, l], xz_um [[x, z], ...]},
                  aperture {diameter_um, expected_dots, center_um [x, z],
                  placement}, substrate {material}
  overlays.mode   {nx, ny, x_um [x0, x1], y_um [y0, y1], encoding "u8-b64",
                  data (row-major, y from bottom), wx_um, wy_um, n_eff, tag,
                  label}: I(x, y) = lateral.field^2 * vertical.field^2 from
                  waveguide.effective_index_ridge, normalized to max 1
  overlays.na_cone {NA, half_angle_deg, n_ambient}
  overlays.emission {lambda_nm, false_colour}
device / nanowire_vertical | nanowire_horizontal (units um): geometry
  {core_radius_nm, outer_radius_nm, shell, dot {height_nm, radius_nm},
  barrier_left_nm, barrier_right_nm, oxide_thickness_nm, conducting_radius_nm,
  segments [{name, t_nm|null, role}], dipole_weights}, overlays {he11
  {beta_HE11, V_number, single_mode}, na_cone, emission}
device / planar_cavity (units um): stacks.thermal (thermal.layers),
  stacks.optical (DBR pairs from nitride_cavity.dbr_diagnostic), geometry
  {mesa_diameter_um, dot_depth_note}, overlays {mode_volume {V_norm,
  lambda_nm, n, radius_um}, standing_wave {lambda_nm, n, period_nm}, dbr}
device / thermal_mesa: stacks.thermal + geometry.mesa_diameter_um only.
band / nitride (units nm, eV): overlays {z_nm, cb_eV, vb_eV, psi2_e, psi2_h
  (normalized to peak 1), E_e_eV, E_h_eV (z-subband levels), dot_nm [lo, hi],
  field_kVcm, overlap_sq}
band / inp: overlays {regions [{name, z_nm [lo, hi], cb_eV, vb_eV, role}],
  levels {E_e_eV, E_h_eV}, arrows [{name, from_eV, to_eV, value_meV}]}
  (no wavefunction arrays: dot_levels does not return psi, none invented)
cascade (units ns): overlays {t_ns, P_G, P_X, P_XX, rate_X, rate_XX,
  rate_esc (1/ns instantaneous gamma*P and k*P), tau_on_ns, period_ns,
  playback [{t0_ns, t1_ns, s_per_ns}], sprite_boost {X, XX, esc},
  levels {E_X_eV, delta_xx_meV, gap_scale, gap_badge}, start_state}
surface (axes x = T_hs K, y = param, z = g2): overlays {x, y, lo [ny][nx],
  hi [ny][nx], iso 0.5, tc_points [{y, T_lo, T_hi}], sheet_note}; or a job
  descriptor {job: True, fn: "fsim_studio.scene:surface_compute", kwargs,
  eta_s, reason} when the computation is too slow for a request.
lattice: overlays {axes {x|y|z: {name, values}}, rows {xi, yi, zi, g2, flux,
  optical_pass, device_pass, eligible, reasons, model_finite_pulse,
  model_tau_cap_density (bool per row; null for campaigns without model axes),
  gates (bool: true only for rows of the gating model -- rt_edge: the
  headline model finite_pulse True + tau_cap_density False; other campaigns:
  every row), facet: {col: [idx]}}, facets {col: [values]}, facet_default
  {col: idx}, banner {text, style "pass"|"fail" (from GATING rows only), pass,
  gating_rows, rows, pass_column, non_gating_passes null | {count, rows,
  text ("... never gates")}}, counts}
"""
from __future__ import annotations

import base64
import copy
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from fsim_core import device as fdev
from fsim_core import (dot_levels, materials, nitride_cavity, nitride_levels, scene_support,
                       waveguide)
from fsim_core.design_meta import default_ranged
from fsim_core.nitride_nanowire_photonics import gan_ordinary_index

ROOT = Path(__file__).resolve().parents[1]
CARDS = ROOT / "cards"
OUT = ROOT / "out"
SCHEMA = "fsim.scene/1"

NANOWIRE_CSV = OUT / "nitride_nanowire" / "full" / "sweep.csv"
CAMPAIGNS = {
    "nitride_cavity": OUT / "nitride_cavity" / "sweep.csv",
    "nitride_nanowire": NANOWIRE_CSV,
    "rt_edge": OUT / "rt_edge" / "sweep.csv",
}

# Default tags for EmissionBlock schema defaults a card does not set, as
# documented in fsim_core/device.py EmissionBlock comments: core_half_nm and
# cladding_nm are hkust_ridge_stack class estimates [E]; etch_depth_nm is an
# unset schema default [A] (03-3d.md); alpha_cm class estimate [E].
_EMISSION_DEFAULT_TAG = {"core_half_nm": "E", "cladding_nm": "E", "etch_depth_nm": "A",
                         "ridge_width_nm": "A", "L_um": "A", "NA": "A", "R_back": "A",
                         "alpha_cm": "E", "lambda_nm": "A"}
_TAG_ORDER = {"V": 3, "DR": 2, "E": 1, "A": 0}
_HC_EV_NM = 1239.841984  # h c [V] CODATA 2018, eV nm (unit conversion only)


# ----------------------------------------------------------------- helpers

def _clean(obj):
    """JSON-safe copy: numpy -> python, NaN/inf -> None."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.ndarray):
        return [_clean(x) for x in obj.tolist()]
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(x) for x in obj]
    return str(obj)


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return f


def _load(card):
    """card: a card name (resolved only through fsim_studio.api_cards:
    _card_path's name check + load_design, so a path or anything outside
    cards/ is refused with CardError 400/404/422) or a DeviceDesign.
    Returns (design, name)."""
    if isinstance(card, fdev.DeviceDesign):
        return card, card.name
    from .api_cards import CardError, load_design
    if not isinstance(card, str) or not card:
        raise CardError(400, "card name required")
    return load_design(card), card


def _catalog(name, kind, extra=""):
    h = hashlib.sha256(f"{SCHEMA}|{kind}|{name}|{extra}".encode()).hexdigest()[:6]
    return f"FS-{name}-{h}"


def _sources(d):
    return (d.provenance or {}).get("sources", {}) or {}


def _tag_of(d, path, default="A"):
    src = _sources(d).get(path)
    if isinstance(src, dict) and src.get("tag"):
        t = str(src["tag"]).strip("[]").upper().split("/")[0]
        return t if t in _TAG_ORDER else default
    return default


def _source_text(d, path):
    src = _sources(d).get(path)
    if isinstance(src, dict):
        return " ".join(str(src.get("source", "")).split())[:240]
    return ""


def _chain_tag(d, scalars=None):
    """Weakest tag in the run's chain: evaluate's tag_chain when present,
    else the weakest tag among the card's provenance sources."""
    if scalars and isinstance(scalars.get("tag_chain"), str) and scalars["tag_chain"]:
        tags = [t.strip("[] ").upper() for t in scalars["tag_chain"].replace(",", " ").split()]
        tags = [t for t in tags if t in _TAG_ORDER]
        if tags:
            return min(tags, key=_TAG_ORDER.get)
    tags = []
    for v in _sources(d).values():
        if isinstance(v, dict) and v.get("tag"):
            t = str(v["tag"]).strip("[]").upper().split("/")[0]
            if t in _TAG_ORDER:
                tags.append(t)
    return min(tags, key=_TAG_ORDER.get) if tags else "A"


def _S(value, unit, tag, label, note=None):
    out = {"value": value, "unit": unit, "tag": tag, "label": label}
    if note:
        out["note"] = note
    return out


def _headline_copy(d):
    """RT-edge headline model switch (configuration, CONTRACT rule 5): finite
    pulse on, tau_cap density scaling off, CW off."""
    h = copy.deepcopy(d)
    h.drive.finite_pulse = True
    h.ret.tau_cap_scales_with_density = False
    h.drive.cw = False
    return h


def _is_edge(d):
    return d.platform == "legacy" and d.emission.type == "edge"


def _b64_u8(arr2d):
    a = np.clip(np.asarray(arr2d, dtype=float), 0.0, 1.0)
    return base64.b64encode((a * 255.0 + 0.5).astype(np.uint8).tobytes()).decode("ascii")


def _base(kind, variant, name, title, units, extra=""):
    return {"schema": SCHEMA, "kind": kind, "variant": variant, "card": name, "title": title,
            "catalog_id": _catalog(name, kind, extra), "units": units, "axes": {},
            "layers": [], "overlays": {}, "scalars": {}, "labels": [], "notes": [],
            "provenance": {"calls": []}}


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _flat_card_values(d):
    """Flattened scalar inputs used to locate a card's own row in a sweep."""
    nw = d.nitride.get("nanowire", {}) if isinstance(d.nitride, dict) else {}
    dot = d.nitride.get("dot", {}) if isinstance(d.nitride, dict) else {}
    return {"core_radius_nm": nw.get("core_radius_nm"), "height_nm": dot.get("height_nm"),
            "x_in": dot.get("x_in"), "T_hs": d.thermal.T_hs, "strain_bound": nw.get("strain_bound"),
            "rep_rate_hz": d.drive.rep_rate_hz, "family": nw.get("family")}


def _nanowire_row(d, name, T=None):
    """The card's own committed row in out/nitride_nanowire/full/sweep.csv
    (row_kind core, card_file == <name>.yaml, matching geometry/T/strain/rep)."""
    rows = [r for r in _read_csv(NANOWIRE_CSV)
            if r.get("card_file") == f"{name}.yaml" and r.get("row_kind") == "core"]
    want = _flat_card_values(d)
    if T is not None:
        want["T_hs"] = T
    def match(r):
        for k, v in want.items():
            if v is None:
                continue
            if isinstance(v, str):
                if r.get(k) != v:
                    return False
            elif not math.isclose(_num(r.get(k)), float(v), rel_tol=1e-9, abs_tol=1e-12):
                return False
        return True
    hits = [r for r in rows if match(r)]
    return (hits[0] if hits else None), len(hits), len(rows)


# ----------------------------------------------------------------- device

def device_scene(card, T=None):
    """3D device model SceneSpec for one card at heat-sink temperature T
    (default: the card's thermal.T_hs)."""
    d, name = _load(card)
    if _is_edge(d):
        return _clean(_edge_scene(d, name, T))
    if d.platform == "ingan_gan_nanowire":
        return _clean(_nanowire_scene(d, name, T))
    if d.platform == "ingan_gan_planar":
        return _clean(_planar_scene(d, name, T))
    return _clean(_mesa_scene(d, name, T))


def _thermal_layers(d):
    out = []
    for i, L in enumerate(d.thermal.layers or []):
        out.append({"name": L.get("name", f"layer {i}"), "t_nm": float(L.get("t_um", 0.0)) * 1000.0,
                    "n": None, "material": L.get("name", ""), "k300": L.get("k300"),
                    "tag": _tag_of(d, f"thermal.layers.{i}", _tag_of(d, "thermal.layers", "A")),
                    "role": "thermal", "source": _source_text(d, "thermal.layers")})
    return out


def _edge_scene(d, name, T):
    T = float(d.thermal.T_hs if T is None else T)
    em = d.emission
    system = copy.copy(fdev._retention_system(d.ret))
    system.T = T
    lv = dot_levels.levels(system)
    lam = em.lambda_nm if em.lambda_nm > 0 else lv.lambda_nm
    g = system.geometry

    def n_at(mat):
        return float(materials.refractive_index(mat.label, lam))

    def etag(field):
        return _tag_of(d, f"emission.{field}", _EMISSION_DEFAULT_TAG.get(field, "A"))

    # Same layer list as fsim_core.device._resolve_edge (bottom -> top).
    spec_layers = [
        ("barrier_lower", system.barrier, em.cladding_nm, False, "cladding", "cladding_nm"),
        ("matrix_lower", system.matrix, em.core_half_nm, False, "core", "core_half_nm"),
        ("dot", system.dot, max(g.height_nm, 0.1), True, "active", None),
        ("matrix_upper", system.matrix, em.core_half_nm, False, "core", "core_half_nm"),
        ("barrier_upper", system.barrier, em.cladding_nm, False, "cladding", "cladding_nm"),
    ]
    wl_layers = [waveguide.Layer(nm, n_at(mat), t, is_dot, material=mat.label)
                 for nm, mat, t, is_dot, _, _ in spec_layers]
    layers = []
    for (nm, mat, t, is_dot, role, field), wl in zip(spec_layers, wl_layers):
        tag = etag(field) if field else _tag_of(d, "ret.system.geometry", "A")
        layers.append({"name": nm, "t_nm": t, "n": wl.n, "material": mat.label, "tag": tag,
                       "role": role, "source": (f"emission.{field}" if field else "ret.system.geometry.height_nm")
                       + f"; n = materials.refractive_index('{mat.label}', {lam:.1f} nm)"})
    total = sum(L["t_nm"] for L in layers)

    rm = waveguide.effective_index_ridge(wl_layers, lam, em.ridge_width_nm, em.etch_depth_nm)
    pad = max(2000.0, em.ridge_width_nm)  # waveguide.effective_index_ridge lateral pad
    xv = np.asarray(rm.lateral.z_nm, float) - (pad + em.ridge_width_nm / 2.0)
    yv = np.asarray(rm.vertical.z_nm, float)
    nx = ny = 256
    xg = np.linspace(xv[0], xv[-1], nx)
    yg = np.linspace(yv[0], yv[-1], ny)
    fx = np.interp(xg, xv, np.asarray(rm.lateral.field, float) ** 2)
    fy = np.interp(yg, yv, np.asarray(rm.vertical.field, float) ** 2)
    I = np.outer(fy, fx)
    I = I / I.max()

    # Static evaluate (card as written) and the headline switch.
    stat = fdev.evaluate(d, T_grid=[T])["scalars"]
    head = fdev.evaluate(_headline_copy(d), T_grid=[T])["scalars"]
    chain = _chain_tag(d, stat)

    dens = float(d.aperture.density_cm2 or 0.0)
    per_um2 = dens * 1e-8
    front_len = 12.0
    rng = np.random.default_rng(20261007)
    n_dots = int(rng.poisson(per_um2 * (em.ridge_width_nm / 1000.0) * front_len))
    xs = rng.uniform(-em.ridge_width_nm / 2000.0, em.ridge_width_nm / 2000.0, n_dots)
    zs = rng.uniform(-front_len, 0.0, n_dots)
    ap_d = float(d.aperture.diameter_um or 0.0)
    expected = per_um2 * math.pi * (ap_d / 2.0) ** 2
    y_dot = em.cladding_nm + em.core_half_nm + max(g.height_nm, 0.1) / 2.0

    spec = _base("device", "edge_ridge", name, f"{name} - edge-emitting ridge", {"length": "um"}, f"T={T}")
    spec["axes"] = {"x": {"label": "lateral x", "unit": "um", "scale": 1, "badge": None},
                    "y": {"label": "growth y", "unit": "um", "scale": 1, "badge": None},
                    "z": {"label": "propagation z", "unit": "um", "scale": 1, "badge": None,
                          "break": {"after_um": front_len, "total_um": em.L_um}}}
    spec["layers"] = layers
    spec["stacks"] = {
        "optical": {"label": "optical stack (waveguide.py)", "layers": layers, "total_nm": total},
        "thermal": {"label": "thermal stack (thermal.py)", "layers": _thermal_layers(d),
                    "total_nm": sum(L["t_nm"] for L in _thermal_layers(d))},
    }
    etch_floor = total - em.etch_depth_nm
    spec["geometry"] = {
        "ridge_width_nm": em.ridge_width_nm, "etch_depth_nm": em.etch_depth_nm,
        "etch_floor_nm": etch_floor, "shoulder_nm": pad, "L_um": em.L_um,
        "front_len_um": front_len, "R_back": em.R_back,
        "etch_through_dot": bool(etch_floor < y_dot),
        "dot": {"height_nm": g.height_nm, "radius_nm": g.radius_nm, "y_nm": y_dot, "shape": g.shape,
                "wl_thickness_nm": g.wl_thickness_nm, "material": system.dot.label,
                "matrix": system.matrix.label},
        "dots": {"density_cm2": dens, "per_um2": per_um2, "seed": 20261007,
                 "region_um": [em.ridge_width_nm / 1000.0, front_len],
                 "xz_um": [[float(a), float(b)] for a, b in zip(xs, zs)],
                 "note": "positions are a seeded Poisson sample at the card density (illustrative placement)"},
        "aperture": {"diameter_um": ap_d, "expected_dots": expected, "center_um": [0.0, -2.0],
                     "placement": "aperture position is not a model input; drawn 2 um behind the facet"},
        "substrate": {"material": system.substrate.label if system.substrate else "substrate",
                      "note": "substrate thickness not modelled (break mark)"},
    }
    spec["overlays"] = {
        "mode": {"nx": nx, "ny": ny, "x_um": [xg[0] / 1000.0, xg[-1] / 1000.0],
                 "y_um": [yg[0] / 1000.0, yg[-1] / 1000.0], "encoding": "u8-b64", "data": _b64_u8(I),
                 "wx_um": rm.wx_um, "wy_um": rm.wy_um, "n_eff": rm.n_eff, "tag": "DR",
                 "label": "guided mode |E|^2 (waveguide.effective_index_ridge)"},
        "na_cone": {"NA": em.NA, "half_angle_deg": math.degrees(math.asin(min(em.NA, 1.0))),
                    "n_ambient": 1.0, "tag": etag("NA")},
        "emission": {"lambda_nm": lam, "false_colour": bool(lam > 700.0)},
    }
    st_lab = "static (non-headline)"
    spec["scalars"] = {
        "g2_headline": _S(head.get("g2_op"), "", chain, "g2(0) headline (finite-pulse)"),
        "flux_headline": _S(head.get("collected_flux_pulsed_s"), "1/s", chain, "collected flux, headline"),
        "g2_static": _S(stat.get("g2_op"), "", chain, f"g2(0) {st_lab}"),
        "edge_beta": _S(stat.get("edge_beta"), "", "DR", "guided beta"),
        "edge_T_facet": _S(stat.get("edge_T_facet"), "", "DR", "facet transmission"),
        "edge_eta_total": _S(stat.get("edge_eta_total"), "", "DR", "eta total (out-coupling)"),
        "edge_n_eff": _S(stat.get("edge_n_eff"), "", "DR", "n_eff"),
        "edge_lambda_nm": _S(lam, "nm", etag("lambda_nm"), "emission wavelength"),
        "mode_wx_um": _S(rm.wx_um, "um", "DR", "mode wx (1/e^2)"),
        "mode_wy_um": _S(rm.wy_um, "um", "DR", "mode wy (1/e^2)"),
        "expected_dots": _S(expected, "", _tag_of(d, "aperture.density_cm2", "A"), "expected dots in aperture"),
        "T_hs": _S(T, "K", _tag_of(d, "thermal.T_hs", "A"), "heat sink T"),
    }
    spec["labels"] = [st_lab + ": g2 static from the card as written; headline = finite-pulse switch",
                      "etch removes dot layer; sidewalls not modelled" if etch_floor < y_dot
                      else "etch stops above the dot layer",
                      "NIR, false colour" if lam > 700.0 else "visible, true colour",
                      f"optical stack {total/1000:.2f} um vs thermal stack "
                      f"{spec['stacks']['thermal']['total_nm']/1000:.2f} um: the two models disagree",
                      f"cavity length L = {em.L_um:g} um, drawn true-scale for {front_len:g} um then broken"]
    spec["bookmarks"] = [{"key": "1", "name": "Facet"}, {"key": "2", "name": "Side cut-away"},
                         {"key": "3", "name": "Top / aperture"}, {"key": "4", "name": "Dot magnifier"}]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": [
        "fsim_core.device._retention_system + dot_levels.levels (stack, as device._resolve_edge)",
        "fsim_core.materials.refractive_index", "fsim_core.waveguide.effective_index_ridge",
        "fsim_core.device.evaluate (card as written; static)",
        "fsim_core.device.evaluate (headline switch: finite_pulse, no tau_cap density, cw off)"]}
    return spec


def _nanowire_scene(d, name, T):
    nw = dict(d.nitride.get("nanowire", {}))
    dot = dict(d.nitride.get("dot", {}))
    ph = dict(d.nitride.get("photonics", {}))
    fam = nw.get("family", ph.get("family", "vertical_photonic"))
    row, n_hit, n_rows = _nanowire_row(d, name, T)
    chain = _chain_tag(d)
    Tq = float(d.thermal.T_hs if T is None else T)
    lam = _num(row.get("lambda_nm")) if row else float("nan")
    n_wire = gan_ordinary_index(lam) if math.isfinite(lam) else float("nan")
    variant = "nanowire_horizontal" if fam.startswith("horizontal") else "nanowire_vertical"
    spec = _base("device", variant, name, f"{name} - disc-in-nanowire ({fam})", {"length": "um"}, f"T={Tq}")
    spec["axes"] = {k: {"label": k, "unit": "um", "scale": 1, "badge": None} for k in "xyz"}
    bl, br = float(nw.get("barrier_left_nm", 0.0)), float(nw.get("barrier_right_nm", 0.0))
    h = float(dot.get("height_nm", 0.0))
    spec["layers"] = [
        {"name": "n-GaN", "t_nm": None, "n": n_wire, "material": "GaN", "tag": "A", "role": "doped_n",
         "source": "segment length not modelled"},
        {"name": "barrier (left)", "t_nm": bl, "n": n_wire, "material": "GaN", "role": "intrinsic",
         "tag": _tag_of(d, "nitride.nanowire.barrier_left_nm", "A"), "source": "nitride.nanowire.barrier_left_nm"},
        {"name": "InGaN disc", "t_nm": h, "n": None, "material": f"In{dot.get('x_in')}Ga N", "role": "active",
         "tag": _tag_of(d, "nitride.dot.height_nm", "A"), "source": "nitride.dot.height_nm"},
        {"name": "barrier (right)", "t_nm": br, "n": n_wire, "material": "GaN", "role": "intrinsic",
         "tag": _tag_of(d, "nitride.nanowire.barrier_right_nm", "A"), "source": "nitride.nanowire.barrier_right_nm"},
        {"name": "p-GaN", "t_nm": None, "n": n_wire, "material": "GaN", "tag": "A", "role": "doped_p",
         "source": "segment length not modelled"},
    ]
    spec["geometry"] = {
        "family": fam, "core_radius_nm": nw.get("core_radius_nm"), "outer_radius_nm": nw.get("outer_radius_nm"),
        "shell": nw.get("shell"), "dot": {"height_nm": h, "radius_nm": dot.get("radius_nm")},
        "barrier_left_nm": bl, "barrier_right_nm": br,
        "oxide_thickness_nm": ph.get("oxide_thickness_nm"),
        "conducting_radius_nm": (d.drive.diode or {}).get("conducting_radius_nm"),
        "dipole_weights": ph.get("dipole_weights"),
        "segments_note": "p/n segment lengths and wire length are not modelled (break marks); "
                         "p on top is a drawing convention",
    }
    num = (lambda k: _num(row.get(k)) if row else float("nan"))
    spec["overlays"] = {
        "he11": {"beta_HE11": num("beta_HE11"), "V_number": num("V_number"),
                 "single_mode": (row.get("single_mode") == "True") if row else None, "tag": "DR"},
        "na_cone": {"NA": ph.get("NA"), "half_angle_deg": math.degrees(math.asin(min(float(ph.get("NA", 0.5)), 1.0))),
                    "n_ambient": ph.get("n_ambient", 1.0), "tag": _tag_of(d, "nitride.photonics.NA", "A")},
        "emission": {"lambda_nm": lam, "false_colour": bool(math.isfinite(lam) and lam > 700.0)},
    }
    rid = row.get("row_id") if row else None
    spec["scalars"] = {
        "g2_op": _S(num("g2_op"), "", chain, "g2(0)"),
        "flux_commanded": _S(num("collected_flux_pulsed_s"), "1/s", chain, "collected flux, commanded"),
        "flux_delivered": _S(num("collected_flux_delivered_s"), "1/s", chain, "collected flux, delivered"),
        "beta_HE11": _S(num("beta_HE11"), "", "DR", "beta HE11"),
        "V_number": _S(num("V_number"), "", "DR", "V number"),
        "lambda_nm": _S(lam, "nm", chain, "emission wavelength"),
        "eta_collection_X": _S(num("eta_collection_X"), "", chain, "collection eta (X)"),
        "T_hs": _S(Tq, "K", _tag_of(d, "thermal.T_hs", "A"), "heat sink T"),
    }
    op = (row.get("optical_pass") == "True") if row else None
    dp = (row.get("device_pass") == "True") if row else None
    spec["verdict"] = {"optical_pass": op, "device_pass": dp,
                       "split": bool(op and not dp), "row_id": rid}
    spec["labels"] = [
        "wire length not modelled (break marks)",
        "nanowire numbers from the committed sweep row "
        + (f"{rid} (out/nitride_nanowire/full/sweep.csv)" if rid else "(no matching row)"),
        "flux shows commanded AND delivered",
    ]
    if not row:
        spec["notes"].append(f"no core row matched card {name} at T={Tq} ({n_rows} rows for this card)")
    elif n_hit > 1:
        spec["notes"].append(f"{n_hit} rows matched; first used")
    spec["bookmarks"] = [{"key": "1", "name": "Wire"}, {"key": "2", "name": "Side cut-away"},
                         {"key": "3", "name": "Top / cone"}, {"key": "4", "name": "Disc magnifier"}]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "row_id": rid,
                          "calls": ["out/nitride_nanowire/full/sweep.csv (core row for this card)",
                                    "fsim_core.nitride_nanowire_photonics.gan_ordinary_index"]}
    return spec


def _planar_scene(d, name, T):
    T = float(d.thermal.T_hs if T is None else T)
    sc = fdev.evaluate(d, T_grid=[T])["scalars"]
    chain = _chain_tag(d, sc)
    lam = float(sc.get("lambda_nm"))
    cav = dict(d.nitride.get("cavity", {}))
    vnorm = float(cav.get("mode_volume_norm", nitride_cavity.NitrideCavityParams.mode_volume_norm))
    opt = scene_support.planar_cavity_optics(lam, vnorm)
    n_gan, V_um3, dbr = opt["n_gan"], opt["V_um3"], opt["dbr"]
    pairs = dbr["pairs"]  # nitride_cavity.dbr_diagnostic's own default
    dbr_layers = []
    for i in range(pairs):
        dbr_layers.append({"name": f"H{i+1}", "t_nm": dbr["d_high_nm"], "n": dbr["n_high"],
                           "material": "high-n dielectric", "tag": "E", "role": "mirror",
                           "source": "nitride_cavity.dbr_diagnostic n_high [E/A]"})
        dbr_layers.append({"name": f"L{i+1}", "t_nm": dbr["d_low_nm"], "n": dbr["n_low"],
                           "material": "low-n dielectric", "tag": "E", "role": "mirror",
                           "source": "nitride_cavity.dbr_diagnostic n_low [E/A]"})
    th = _thermal_layers(d)
    spec = _base("device", "planar_cavity", name, f"{name} - planar InGaN/GaN cavity", {"length": "um"}, f"T={T}")
    spec["axes"] = {k: {"label": k, "unit": "um", "scale": 1, "badge": None} for k in "xyz"}
    spec["layers"] = th
    spec["stacks"] = {
        "thermal": {"label": "thermal stack (thermal.py)", "layers": th, "total_nm": sum(L["t_nm"] for L in th)},
        "optical": {"label": "DBR diagnostic (nitride_cavity.dbr_diagnostic)", "layers": dbr_layers,
                    "total_nm": dbr["total_nm"]},
    }
    spec["geometry"] = {"mesa_diameter_um": d.thermal.mesa_diameter_um,
                        "substrate": (d.thermal.substrate or {}).get("name", "substrate"),
                        "dot_depth_note": "[A] dot depth in the GaN film is not a model input; drawn mid-film",
                        "dot": {"height_nm": d.nitride.get("dot", {}).get("height_nm"),
                                "radius_nm": d.nitride.get("dot", {}).get("radius_nm")}}
    spec["overlays"] = {
        "mode_volume": {"V_norm": vnorm, "lambda_nm": lam, "n": n_gan, "V_um3": V_um3,
                        "radius_um": opt["radius_um"], "tag": "A",
                        "label": "[A] mode volume, no 3D geometry in model"},
        "standing_wave": {"lambda_nm": lam, "n": n_gan, "period_nm": opt["standing_wave_period_nm"],
                          "label": "standing-wave glyph: period lambda/2n (n_GaN); envelope illustrative"},
        "dbr": {"reflectance": dbr["reflectance"], "stopband_nm": list(dbr["stopband_nm"]),
                "d_high_nm": dbr["d_high_nm"], "d_low_nm": dbr["d_low_nm"], "pairs": pairs, "tag": "E"},
        "emission": {"lambda_nm": lam, "false_colour": bool(lam > 700.0)},
    }
    spec["scalars"] = {
        "g2_op": _S(sc.get("g2_op"), "", chain, "g2(0)"),
        "flux": _S(sc.get("collected_flux_pulsed_s"), "1/s", chain, "collected flux (pre-detector)"),
        "Q": _S(sc.get("Q"), "", _tag_of(d, "nitride.cavity.Q", "A"), "cavity Q"),
        "kappa_meV": _S(sc.get("kappa_meV"), "meV", "DR", "kappa"),
        "Fp_add": _S(sc.get("Fp_add"), "", "DR", "Purcell Fp (ideal)"),
        "lambda_nm": _S(lam, "nm", chain, "emission wavelength"),
        "dbr_R": _S(dbr["reflectance"], "", "E", "DBR R (10 pairs, free space)"),
        "T_hs": _S(T, "K", _tag_of(d, "thermal.T_hs", "A"), "heat sink T"),
    }
    spec["verdict"] = {"device_pass": sc.get("device_pass"), "optical_pass": None, "split": False}
    spec["labels"] = ["[A] mode volume, no 3D geometry in model",
                      f"thermal stack {spec['stacks']['thermal']['total_nm']/1000:.2f} um vs DBR diagnostic "
                      f"{spec['stacks']['optical']['total_nm']/1000:.2f} um: different models",
                      "dot depth not modelled"]
    spec["bookmarks"] = [{"key": "1", "name": "Mesa"}, {"key": "2", "name": "Side cut-away"},
                         {"key": "3", "name": "Top"}, {"key": "4", "name": "Dot magnifier"}]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": [
        "fsim_core.device.evaluate", "fsim_core.nitride_cavity.dbr_diagnostic(lambda_nm)",
        "fsim_core.nitride_nanowire_photonics.gan_ordinary_index"]}
    return spec


def _mesa_scene(d, name, T):
    T = float(d.thermal.T_hs if T is None else T)
    sc = fdev.evaluate(d, T_grid=[T])["scalars"]
    chain = _chain_tag(d, sc)
    th = _thermal_layers(d)
    spec = _base("device", "thermal_mesa", name, f"{name} - mesa (thermal stack)", {"length": "um"}, f"T={T}")
    spec["axes"] = {k: {"label": k, "unit": "um", "scale": 1, "badge": None} for k in "xyz"}
    spec["layers"] = th
    spec["stacks"] = {"thermal": {"label": "thermal stack (thermal.py)", "layers": th,
                                  "total_nm": sum(L["t_nm"] for L in th)}}
    spec["geometry"] = {"mesa_diameter_um": d.thermal.mesa_diameter_um,
                        "substrate": (d.thermal.substrate or {}).get("name", "substrate")}
    spec["overlays"] = {"emission": {"lambda_nm": None, "false_colour": False}}
    spec["scalars"] = {"g2_op": _S(sc.get("g2_op"), "", chain, "g2(0)"),
                       "T_j": _S(sc.get("T_j_op"), "K", chain, "junction T"),
                       "T_hs": _S(T, "K", _tag_of(d, "thermal.T_hs", "A"), "heat sink T")}
    spec["labels"] = ["only the thermal stack is modelled for this card (no optical geometry)"]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": ["fsim_core.device.evaluate"]}
    return spec


# ----------------------------------------------------------------- band

def band_scene(card, T=None):
    """Band-edge profile: nitride tilted CB/VB with psi_e/psi_h at the
    evaluated operating point, or InP/GaAsP piecewise-flat edges + levels."""
    d, name = _load(card)
    if d.platform in ("ingan_gan_planar", "ingan_gan_nanowire"):
        return _clean(_band_nitride(d, name, T))
    return _clean(_band_inp(d, name, T))


def _band_nitride(d, name, T):
    dot_kw = dict(d.nitride.get("dot", {}))
    dot_kw = {k: v for k, v in dot_kw.items()
              if k in {f.name for f in nitride_levels.NitrideDotSystem.__dataclass_fields__.values()}}
    calls = []
    if d.platform == "ingan_gan_planar":
        sc = fdev.evaluate(d, T_grid=[float(d.thermal.T_hs if T is None else T)])["scalars"]
        ext, Tj = float(sc["applied_field_kVcm"]), float(sc["T_j"])
        src_note = "operating point of fsim_core.device.evaluate (applied_field_kVcm, T_j)"
        calls.append("fsim_core.device.evaluate")
        chain = _chain_tag(d, sc)
    else:
        # The nanowire tier solves its levels with nitride_nanowire_levels
        # (barrier-clad wire model); the planar z_profile would not reproduce
        # its overlap_sq, so no tilted-band view is drawn for it.
        row, _, _ = _nanowire_row(d, name, T)
        spec = _base("band", "unavailable", name, f"{name} - band profile", {"length": "nm", "energy": "eV"})
        spec["labels"] = ["nanowire levels come from fsim_core.nitride_nanowire_levels; the planar "
                          "nitride_levels.z_profile is a different model, so no band view is drawn"]
        if row:
            spec["scalars"] = {"overlap_sq": _S(_num(row.get("overlap_sq")), "", _chain_tag(d), "overlap^2 (sweep row)"),
                               "field_kVcm": _S(_num(row.get("field_kVcm")), "kV/cm", _chain_tag(d), "field (sweep row)")}
        return spec
    system = nitride_levels.NitrideDotSystem(**{**dot_kw, "external_field_kVcm": ext})
    zp = nitride_levels.z_profile(system, Tj)
    lv = nitride_levels.levels(system, Tj)
    calls += ["fsim_core.nitride_levels.z_profile", "fsim_core.nitride_levels.levels"]
    z = np.asarray(zp["z_nm"]); h = zp["effective_height_nm"]
    win = 12.0 + h / 2.0
    zg = np.linspace(-win, win, 481)
    cb = np.interp(zg, z, zp["cb_eV"]); vb = np.interp(zg, z, zp["vb_eV"])
    pe = np.interp(zg, z, np.asarray(zp["psi_e"]) ** 2); ph_ = np.interp(zg, z, np.asarray(zp["psi_h"]) ** 2)
    spec = _base("band", "nitride", name, f"{name} - tilted band edges (QCSE)", {"length": "nm", "energy": "eV"},
                 f"Tj={Tj}")
    spec["axes"] = {"x": {"label": "growth z", "unit": "nm", "scale": 1, "badge": None},
                    "y": {"label": "energy", "unit": "eV", "scale": 1, "badge": "E axis: 1 eV drawn as 4 nm"},
                    "z": {"label": "depth (extrusion, no meaning)", "unit": "nm", "scale": 1, "badge": None}}
    spec["overlays"] = {
        "z_nm": zg, "cb_eV": cb, "vb_eV": vb,
        "psi2_e": pe / max(pe.max(), 1e-300), "psi2_h": ph_ / max(ph_.max(), 1e-300),
        "E_e_eV": zp["E_e_z_eV"], "E_h_eV": zp["E_h_z_eV"], "dot_nm": [-h / 2.0, h / 2.0],
        "field_kVcm": zp["field_kVcm"], "overlap_sq": zp["overlap_sq"],
        "energy_per_nm": 0.25, "level_note": "z-subband ground levels; radial confinement energy not drawn",
    }
    spec["scalars"] = {
        "field_kVcm": _S(zp["field_kVcm"], "kV/cm", chain, "field in dot"),
        "overlap_sq": _S(zp["overlap_sq"], "", chain, "e-h overlap^2"),
        "overlap_eval": _S(sc.get("overlap_sq"), "", chain, "overlap^2 (evaluate)"),
        "tilt_eV": _S(scene_support.band_tilt_eV(zp["field_kVcm"], h), "eV", chain, "tilt across dot"),
        "E_X_eV": _S(lv.E_X_eV, "eV", chain, "E_X"),
        "lambda_nm": _S(lv.lambda_nm, "nm", chain, "lambda"),
        "height_nm": _S(h, "nm", _tag_of(d, "nitride.dot.height_nm", "A"), "dot height (effective)"),
        "T_j": _S(Tj, "K", chain, "T_j"),
    }
    spec["labels"] = [src_note, "z at true scale; energy axis has its own aspect (badge)",
                      "|psi|^2 drawn at their z-subband energies, peak-normalized"]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": calls}
    return spec


def _band_inp(d, name, T):
    if not (d.ret.preset or d.ret.system):
        spec = _base("band", "unavailable", name, f"{name} - band profile", {"length": "nm", "energy": "eV"})
        spec["labels"] = ["no confinement stack on this card (ret.system/preset unset): no band view"]
        return spec
    system = copy.copy(fdev._retention_system(d.ret))
    T = float(d.thermal.T_hs if T is None else T)
    system.T = T
    lv = dot_levels.levels(system)
    g = system.geometry
    em = d.emission
    h = float(g.height_nm)
    core = float(em.core_half_nm)
    Eg = float(lv.Eg_dot_strained_eV)
    Ve, Vh, Vemb, Vhmb = lv.V_e / 1000, float(lv.V_h) / 1000, lv.V_e_mb / 1000, float(lv.V_h_mb) / 1000
    regions = [
        {"name": f"barrier {system.barrier.label}", "z_nm": [-(h / 2 + core + 40.0), -(h / 2 + core)],
         "cb_eV": Eg + Ve + Vemb, "vb_eV": -Vh - Vhmb, "role": "cladding", "open": "left"},
        {"name": f"matrix {system.matrix.label}", "z_nm": [-(h / 2 + core), -h / 2],
         "cb_eV": Eg + Ve, "vb_eV": -Vh, "role": "core"},
        {"name": f"dot {system.dot.label}", "z_nm": [-h / 2, h / 2], "cb_eV": Eg, "vb_eV": 0.0, "role": "active"},
        {"name": f"matrix {system.matrix.label}", "z_nm": [h / 2, h / 2 + core],
         "cb_eV": Eg + Ve, "vb_eV": -Vh, "role": "core"},
        {"name": f"barrier {system.barrier.label}", "z_nm": [h / 2 + core, h / 2 + core + 40.0],
         "cb_eV": Eg + Ve + Vemb, "vb_eV": -Vh - Vhmb, "role": "cladding", "open": "right"},
    ]
    E_e = Eg + lv.E_e / 1000.0
    E_h = -float(lv.E_h) / 1000.0
    spec = _base("band", "inp", name, f"{name} - band edges and confined levels", {"length": "nm", "energy": "eV"},
                 f"T={T}")
    spec["axes"] = {"x": {"label": "growth z", "unit": "nm", "scale": 1, "badge": None},
                    "y": {"label": "energy (dot VB top = 0)", "unit": "eV", "scale": 1,
                          "badge": "E axis: 1 eV drawn as 120 nm"},
                    "z": {"label": "depth (extrusion, no meaning)", "unit": "nm", "scale": 1, "badge": None}}
    spec["overlays"] = {
        "regions": regions, "levels": {"E_e_eV": E_e, "E_h_eV": E_h},
        "arrows": [
            {"name": "dE_e_matrix", "from_eV": E_e, "to_eV": Eg + Ve, "value_meV": lv.dE_e_matrix, "carrier": "e"},
            {"name": "dE_e_barrier", "from_eV": E_e, "to_eV": Eg + Ve + Vemb, "value_meV": lv.dE_e_barrier, "carrier": "e"},
            {"name": "dE_h_matrix", "from_eV": E_h, "to_eV": -Vh, "value_meV": float(lv.dE_h_matrix), "carrier": "h"},
            {"name": "dE_h_barrier", "from_eV": E_h, "to_eV": -Vh - Vhmb, "value_meV": float(lv.dE_h_barrier), "carrier": "h"},
        ],
        "energy_per_nm": 1.0 / 120.0, "type": lv.type,
        "barrier_note": f"barrier continues {em.cladding_nm:g} nm (drawn 40 nm, break mark)",
    }
    tg = _tag_of(d, "ret.system.geometry.height_nm", "A")
    spec["scalars"] = {
        "V_e": _S(lv.V_e, "meV", "DR", "V_e (dot->matrix CB)"),
        "V_h": _S(float(lv.V_h), "meV", "DR", "V_h (dot->matrix VB)"),
        "E_e": _S(lv.E_e, "meV", "DR", "E_e"), "E_h": _S(float(lv.E_h), "meV", "DR", "E_h"),
        "Eg_dot": _S(Eg, "eV", "DR", "Eg dot (strained)"),
        "Eg_matrix": _S(lv.Eg_matrix_eV, "eV", "DR", "Eg matrix"),
        "E_X_eV": _S(lv.E_X_eV, "eV", "DR", "E_X"), "lambda_nm": _S(lv.lambda_nm, "nm", "DR", "lambda"),
        "height_nm": _S(h, "nm", tg, "dot height"),
    }
    spec["labels"] = ["levels only: dot_levels returns no wavefunction arrays, none drawn",
                      f"type {lv.type}; energies relative to the dot VB top",
                      spec["overlays"]["barrier_note"]]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml",
                          "calls": ["fsim_core.dot_levels.levels(system)"]}
    return spec


# ----------------------------------------------------------------- cascade

def cascade_scene(card, T=None):
    """XX -> X -> 0 cascade over one drive period: populations from
    fsim_core.lindblad.evolve with the evaluated operating-point rates."""
    d, name = _load(card)
    T = float(d.thermal.T_hs if T is None else T)
    lab_head = None
    if _is_edge(d):
        sc = fdev.evaluate(_headline_copy(d), T_grid=[T])["scalars"]
        r_ns, gX, gXX = sc["finite_pulse_r_ns"], sc["finite_pulse_gamma_X_ns"], sc["finite_pulse_gamma_XX_ns"]
        kX, kXX = sc["finite_pulse_k_X"], sc["finite_pulse_k_XX"]
        tau_on, tau_dark = sc["finite_pulse_tau_on_ns"], sc["finite_pulse_tau_dark_ns"]
        g2, mean_x = sc["g2_op"], sc.get("finite_pulse_mean_counts_x")
        loading = "rectangular"
        system = copy.copy(fdev._retention_system(d.ret)); system.T = T
        E_X = dot_levels.levels(system).E_X_eV
        lab_head = "headline (finite-pulse) operating point"
        src = "fsim_core.device.evaluate (headline switch) finite_pulse_* rates"
        chain = _chain_tag(d, sc)
    elif d.platform == "ingan_gan_planar":
        sc = fdev.evaluate(d, T_grid=[T])["scalars"]
        r_ns = scene_support.per_s_to_per_ns(_num(sc.get("r_dot_s")))
        gX, gXX, kX, kXX = sc["gamma_X_ns"], sc["gamma_XX_ns"], sc["k_X_ns"], sc["k_XX_ns"]
        period, tau_on, tau_dark = scene_support.pulse_timing_ns(float(d.drive.rep_rate_hz),
                                                                 duty=float(d.drive.duty))
        g2, mean_x, E_X = sc["g2_op"], sc.get("mean_counts_x"), sc.get("E_X_eV")
        loading = d.drive.cycle_loading
        src = "fsim_core.device.evaluate scalars (r_dot_s, gamma_*, k_*)"
        chain = _chain_tag(d, sc)
    elif d.platform == "ingan_gan_nanowire":
        row, _, _ = _nanowire_row(d, name, T)
        if not row:
            spec = _base("cascade", "unavailable", name, f"{name} - cascade", {"time": "ns"})
            spec["labels"] = ["no committed sweep row for this card/T"]
            return _clean(spec)
        r_ns = scene_support.per_s_to_per_ns(_num(row["r_dot_s"]))
        gX, gXX, kX, kXX = (_num(row[k]) for k in ("gamma_X_ns", "gamma_XX_ns", "k_X_ns", "k_XX_ns"))
        period, tau_on, tau_dark = scene_support.pulse_timing_ns(_num(row["rep_rate_hz"]),
                                                                 tau_pulse_ns=_num(row["tau_pulse_ns"]))
        g2, mean_x, E_X = _num(row["g2_op"]), _num(row.get("mean_counts_x")), _num(row.get("E_X_eV"))
        loading = row.get("cycle_loading", "rectangular")
        src = f"out/nitride_nanowire/full/sweep.csv row {row.get('row_id')}"
        chain = _chain_tag(d)
    else:
        spec = _base("cascade", "unavailable", name, f"{name} - cascade", {"time": "ns"})
        spec["labels"] = ["no pulsed operating-point rates for this card (legacy CW path): cascade not drawn"]
        return _clean(spec)

    pump_ratio = float(d.drive.cw_pump_ratio)
    vals = [r_ns, gX, gXX, kX, kXX, tau_on, tau_dark]
    if not all(math.isfinite(_num(v)) for v in vals):
        spec = _base("cascade", "unavailable", name, f"{name} - cascade", {"time": "ns"})
        spec["labels"] = ["operating-point rates are not finite at this T (invalid row)"]
        return _clean(spec)
    tr = scene_support.cascade_trajectory(
        r_ns=r_ns, gamma_X_ns=gX, gamma_XX_ns=gXX, k_X=kX, k_XX=kXX, tau_on_ns=tau_on,
        tau_dark_ns=tau_dark, pump_ratio=pump_ratio, deterministic_pair=(loading == "deterministic_pair"))
    t = tr["t_ns"]
    P = {k: tr[k] for k in ("P_G", "P_X", "P_XX")}
    rate_X, rate_XX, rate_esc = tr["rate_X"], tr["rate_XX"], tr["rate_esc"]
    target = {"X": 4.0, "XX": 2.0, "esc": 14.0}  # sprites per loop (display choice)
    boost = {}
    for k in ("X", "XX", "esc"):
        n = tr["per_pulse"][k]
        boost[k] = (target[k] / n) if n > 0 else 0.0
    period = tr["period_ns"]
    playback = [{"t0_ns": 0.0, "t1_ns": float(tau_on), "s_per_ns": 2.0 / float(tau_on)},
                {"t0_ns": float(tau_on), "t1_ns": period, "s_per_ns": 6.0 / float(tau_dark)}]
    dxx = float(d.dot.delta_xx)
    E_X = _num(E_X)
    gap_scale = 20.0
    spec = _base("cascade", loading, name, f"{name} - XX -> X -> 0 cascade", {"time": "ns", "energy": "eV"},
                 f"T={T}")
    spec["axes"] = {"y": {"label": "energy (schematic)", "unit": "eV", "scale": gap_scale,
                          "badge": f"XX-X gap x{gap_scale:g}"},
                    "x": {"label": "time", "unit": "ns", "scale": 1, "badge": None}}
    spec["overlays"] = {
        "t_ns": t, "P_G": P["P_G"], "P_X": P["P_X"], "P_XX": P["P_XX"],
        "rate_X": rate_X, "rate_XX": rate_XX, "rate_esc": rate_esc,
        "tau_on_ns": float(tau_on), "period_ns": period, "playback": playback,
        "sprite_boost": boost, "start_state": tr["start_state"],
        "levels": {"E_X_eV": E_X, "delta_xx_meV": dxx, "gap_scale": gap_scale,
                   "gap_badge": f"XX-X gap x{gap_scale:g}"},
        "glow_scale": {"kind": "log10", "floor": 1e-6,
                       "label": "level glow = log10 population (floor 1e-6)"},
    }
    spec["scalars"] = {
        "g2_op": _S(g2, "", chain, "g2(0) evaluated" + (" (headline)" if lab_head else "")),
        "mean_counts_x": _S(mean_x, "", chain, "X counts per pulse (evaluated)"),
        "r_ns": _S(float(r_ns), "1/ns", chain, "pump r (on)"),
        "gamma_X_ns": _S(float(gX), "1/ns", chain, "gamma_X"),
        "gamma_XX_ns": _S(float(gXX), "1/ns", chain, "gamma_XX"),
        "k_X_ns": _S(float(kX), "1/ns", chain, "escape k_X"),
        "tau_on_ns": _S(float(tau_on), "ns", _tag_of(d, "drive.duty", "A"), "pulse on"),
        "period_ns": _S(period, "ns", _tag_of(d, "drive.rep_rate_hz", "A"), "period"),
        "delta_xx": _S(dxx, "meV", _tag_of(d, "dot.delta_xx", "A"), "XX binding delta_xx"),
    }
    spec["labels"] = ["illustrative sampling of computed rates",
                      "g2 shown is the evaluated g2_op, never counted from the animation",
                      "time is slowed: two playback speeds, labelled live",
                      f"sprite rates boosted: X x{boost['X']:.3g}, XX x{boost['XX']:.3g}, escape x{boost['esc']:.3g}"]
    if lab_head:
        spec["labels"].append(lab_head)
    spec["provenance"] = {"card_file": f"cards/{name}.yaml",
                          "calls": [src, "fsim_core.lindblad.build_system(levels=3)", "fsim_core.lindblad.evolve"]}
    return _clean(spec)


# ----------------------------------------------------------------- surface

_SURF_T = (4.0, 350.0, 30)
# Reduced grid for cards whose evaluate costs ~0.24 s per T point (RT edge):
# 5 param values x 5 T points, envelope over dot.delta_xx only, so the job
# finishes in well under a minute.
_EDGE_GRID = {"n_values": 5, "T": (230.0, 330.0, 5), "ranged_paths": ("dot.delta_xx",)}
_EDGE_EVAL_S = 0.4  # measured 2026-10-07: 144 evaluations in 56.6 s


def _param_grid(d, param, n=12):
    if param == "w":
        return "filter.w", list(np.geomspace(0.25, 8.0, n)), "filter window w", "meV"
    if param == "mu":
        return "drive.mu", list(np.linspace(0.1, 1.0, n)), "mean loading mu", ""
    raise ValueError("param must be 'w' or 'mu'")


def surface_compute(card, param="w", mode="envelope", progress=None, n_values=12, T=None,
                    ranged_paths=None):
    """The (possibly slow) g2(T, param) surface: for each param value, run
    evaluate_envelope over the card's default ranged [A] inputs (restricted
    to `ranged_paths` when given; mode 'point' ranges nothing) on a T grid
    (T = (lo, hi, n), default 4-350 K x 30).  Returns the surface SceneSpec.
    progress(k, n) is called after each param value."""
    d, name = _load(card)
    path, values, plabel, punit = _param_grid(d, param, int(n_values))
    Ts = np.linspace(*(T or _SURF_T))
    lo, hi, tcs = [], [], []
    base = copy.deepcopy(d)
    if path == "filter.w":
        base.filter.enabled = True
        base.filter.auto_w = False
    ranged = default_ranged(base) if mode == "envelope" else {}
    if ranged_paths is not None:
        ranged = {k: v for k, v in ranged.items() if k in set(ranged_paths)}
    ranged.pop(path, None)
    for k, v in enumerate(values):
        dd = copy.deepcopy(base)
        fdev._set_path(dd, path, float(v))
        env = fdev.evaluate_envelope(dd, ranged, T_grid=Ts)
        a, b = env["bands"]["g2"]
        lo.append(np.asarray(a, float)); hi.append(np.asarray(b, float))
        tc = env["scalar_bands"].get("T_c", (float("nan"), float("nan")))
        tcs.append({"y": float(v), "T_lo": tc[0], "T_hi": tc[1]})
        if progress:
            progress(k + 1, len(values))
    chain = _chain_tag(d)
    spec = _base("surface", "g2_T_" + param, name, f"{name} - g2(T, {param})",
                 {"x": "K", "y": punit, "z": ""}, f"{param}|{mode}|{len(values)}|{len(Ts)}")
    spec["axes"] = {"x": {"label": "heat-sink T", "unit": "K", "scale": 1, "badge": None},
                    "y": {"label": plabel, "unit": punit, "scale": 1, "badge": None,
                          "log": path == "filter.w"},
                    "z": {"label": "g2(0)", "unit": "", "scale": 1, "badge": None}}
    spec["overlays"] = {"x": Ts, "y": values, "lo": np.array(lo), "hi": np.array(hi), "iso": 0.5,
                        "tc_points": tcs,
                        "sheet_note": ("lo/hi sheets = pointwise min/max over ranged inputs "
                                       + ", ".join(sorted(ranged)) if ranged else
                                       "point evaluation: lo = hi (no ranged inputs)")}
    spec["scalars"] = {"n_values": _S(len(values), "", chain, f"{param} values"),
                       "n_T": _S(len(Ts), "", chain, "T points")}
    spec["labels"] = ["two sheets (lo, hi), never one averaged surface",
                      "g2 = 0.5 isoline is the T_c ridge; dots = evaluate_envelope T_c interval",
                      f"{param} grid is a display sweep choice ({len(values)} x {len(Ts)} T points)"]
    if _is_edge(d):
        spec["labels"].append("static (non-headline): card as written; reduced grid, envelope over "
                              + ", ".join(sorted(ranged)) + " only")
    spec["provenance"] = {"card_file": f"cards/{name}.yaml",
                          "calls": ["fsim_core.device.evaluate_envelope (default_ranged)"]}
    return _clean(spec)


def surface_scene(card, param="w", mode="envelope"):
    """g2(T, param) surface.  Fast cards (legacy, not edge) are computed
    directly.  The RT edge card (~0.24 s per T point) returns a job
    descriptor {job, fn, kwargs, eta_s} on a reduced grid (< 1 min).  The
    planar-nitride card reads out/nitride_cavity/sweep.csv (x = T_hs,
    y = Q) and nanowire cards out/nitride_nanowire/full/sweep.csv
    (x = T_hs, y = core_radius_nm); both are instant."""
    d, name = _load(card)
    if d.platform == "ingan_gan_nanowire":
        return _clean(_surface_nanowire(d, name))
    if d.platform == "ingan_gan_planar":
        return _clean(_surface_cavity_csv(d, name))
    if _is_edge(d):
        g = _EDGE_GRID
        ranged = {k: v for k, v in default_ranged(d).items() if k in g["ranged_paths"]}
        if mode == "envelope":
            per_env = 2 ** len(ranged) + 1 + len(ranged) * 2 ** max(len(ranged) - 1, 0)
        else:
            per_env = 1
        eta = _EDGE_EVAL_S * per_env * g["n_values"] * g["T"][2]
        return _clean({"job": True, "fn": "fsim_studio.scene:surface_compute",
                       "kwargs": {"card": name, "param": param, "mode": mode, "n_values": g["n_values"],
                                  "T": list(g["T"]), "ranged_paths": list(g["ranged_paths"])},
                       "eta_s": eta, "schema": SCHEMA, "kind": "surface", "card": name,
                       "reason": "evaluate costs ~0.24 s per T point on this card; "
                                 "reduced grid run as a background job"})
    return surface_compute(name, param, mode)


def _surface_cavity_csv(d, name):
    rows = [r for r in _read_csv(CAMPAIGNS["nitride_cavity"])
            if r.get("cycle_loading") == d.drive.cycle_loading
            and math.isclose(_num(r.get("current_uA")), float(d.drive.I_uA), rel_tol=1e-9)]
    Ts = sorted({_num(r["T_hs"]) for r in rows})
    Qs = sorted({_num(r["Q"]) for r in rows})
    lo = np.full((len(Qs), len(Ts)), np.nan)
    hi = np.full_like(lo, np.nan)
    for r in rows:
        g = _num(r.get("g2_op"))
        if not math.isfinite(g):
            continue
        i, j = Qs.index(_num(r["Q"])), Ts.index(_num(r["T_hs"]))
        lo[i, j] = g if not math.isfinite(lo[i, j]) else min(lo[i, j], g)
        hi[i, j] = g if not math.isfinite(hi[i, j]) else max(hi[i, j], g)
    spec = _base("surface", "g2_T_Q", name, f"{name} - g2(T, Q) from sweep",
                 {"x": "K", "y": "", "z": ""}, "csv")
    spec["axes"] = {"x": {"label": "heat-sink T", "unit": "K", "scale": 1, "badge": None},
                    "y": {"label": "cavity Q", "unit": "", "scale": 1, "badge": None, "log": True},
                    "z": {"label": "g2(0)", "unit": "", "scale": 1, "badge": None}}
    spec["overlays"] = {"x": Ts, "y": Qs, "lo": lo, "hi": hi, "iso": 0.5, "tc_points": [],
                        "sheet_note": "lo/hi = min/max over the sweep's height_nm x radius_nm x x_in "
                                      f"(cycle_loading {d.drive.cycle_loading}, I {d.drive.I_uA:g} uA)"}
    spec["scalars"] = {"rows": _S(len(rows), "", "A", "sweep rows used")}
    spec["labels"] = ["two sheets (lo, hi), never one averaged surface",
                      "from out/nitride_cavity/sweep.csv (committed), not recomputed; "
                      f"{len(Ts)} T x {len(Qs)} Q grid"]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": ["out/nitride_cavity/sweep.csv"]}
    return spec


def _surface_nanowire(d, name):
    rows = [r for r in _read_csv(NANOWIRE_CSV)
            if r.get("card_file") == f"{name}.yaml" and r.get("row_kind") == "core"]
    Ts = sorted({_num(r["T_hs"]) for r in rows})
    Rs = sorted({_num(r["core_radius_nm"]) for r in rows})
    lo = np.full((len(Rs), len(Ts)), np.nan); hi = np.full_like(lo, np.nan)
    for r in rows:
        g = _num(r.get("g2_op"))
        if not math.isfinite(g):
            continue
        i, j = Rs.index(_num(r["core_radius_nm"])), Ts.index(_num(r["T_hs"]))
        lo[i, j] = g if not math.isfinite(lo[i, j]) else min(lo[i, j], g)
        hi[i, j] = g if not math.isfinite(hi[i, j]) else max(hi[i, j], g)
    spec = _base("surface", "g2_T_core_radius", name, f"{name} - g2(T, core radius) from sweep",
                 {"x": "K", "y": "nm", "z": ""}, "csv")
    spec["axes"] = {"x": {"label": "heat-sink T", "unit": "K", "scale": 1, "badge": None},
                    "y": {"label": "core radius", "unit": "nm", "scale": 1, "badge": None, "log": True},
                    "z": {"label": "g2(0)", "unit": "", "scale": 1, "badge": None}}
    spec["overlays"] = {"x": Ts, "y": Rs, "lo": lo, "hi": hi, "iso": 0.5, "tc_points": [],
                        "sheet_note": "lo/hi = min/max over the sweep's other core factors "
                                      "(height_nm, x_in, strain_bound, rep_rate_hz)"}
    spec["scalars"] = {"rows": _S(len(rows), "", "A", "sweep rows used")}
    spec["labels"] = ["two sheets (lo, hi), never one averaged surface",
                      "from out/nitride_nanowire/full/sweep.csv (committed), not recomputed"]
    spec["provenance"] = {"card_file": f"cards/{name}.yaml", "calls": [str(NANOWIRE_CSV.relative_to(ROOT))]}
    return spec


# ----------------------------------------------------------------- lattice

_LATTICE = {
    "nitride_cavity": {"axes": ("height_nm", "radius_nm", "x_in"), "g2": "g2_op",
                       "flux": "collected_flux_pulsed_s", "pass": "device_pass",
                       "facets": ("Q", "T_hs", "current_uA", "cycle_loading"), "filter": None,
                       "prefer": {"Q": "2000.0", "T_hs": "300.0", "current_uA": "0.02",
                                  "cycle_loading": "rectangular"}},
    "nitride_nanowire": {"axes": ("T_hs", "core_radius_nm", "height_nm"), "g2": "g2_op",
                         "flux": "collected_flux_pulsed_s", "pass": "device_pass",
                         "facets": ("family", "x_in", "strain_bound", "rep_rate_hz", "card_file"),
                         "filter": ("row_kind", "core"),
                         "prefer": {"family": "vertical_photonic", "x_in": "0.4", "strain_bound": "relaxed",
                                    "rep_rate_hz": "200000000.0",
                                    "card_file": "nitride-nanowire-vertical-pulse-design.yaml"}},
    "rt_edge": {"axes": ("T_hs_K", "delta_xx_meV", "gamma300_meV"), "g2": "g2_pulsed",
                "flux": "collected_flux_pulsed_s", "pass": "headline_pass",
                "facets": ("card_id", "irf_ps", "I_uA", "emission_NA", "emission_R_back", "emission_L_um"),
                "filter": None, "prefer": {"card_id": "edge-inp-gaasp-design"}},
}


def lattice_scene(campaign, x=None, y=None, z=None, headline=True):
    """Factorial sweep lattice from a committed sweep CSV: one glyph per row
    (colour = g2, size = log flux), honest pass counts in the banner."""
    if campaign not in _LATTICE:
        raise ValueError(f"unknown campaign {campaign!r}; one of {sorted(_LATTICE)}")
    cfg = _LATTICE[campaign]
    path = CAMPAIGNS[campaign]
    rows_all = _read_csv(path)
    rows = rows_all
    labels = []
    if cfg["filter"]:
        k, v = cfg["filter"]
        rows = [r for r in rows if r.get(k) == v]
        labels.append(f"rows filtered to {k} = {v}")
    model_counts = None
    if campaign == "rt_edge":
        model_counts = []
        for fp in ("True", "False"):
            for tc in ("False", "True"):
                sub = [r for r in rows_all if r["model_finite_pulse"] == fp and r["model_tau_cap_density"] == tc]
                model_counts.append({"model_finite_pulse": fp == "True", "model_tau_cap_density": tc == "True",
                                     "headline_model": fp == "True" and tc == "False",
                                     "pass": sum(r["headline_pass"] == "True" for r in sub), "rows": len(sub)})
        if headline:
            rows = [r for r in rows if r["model_finite_pulse"] == "True" and r["model_tau_cap_density"] == "False"]
            labels.append("filtered to the headline model (model_finite_pulse True, model_tau_cap_density False)")
    ax = [x or cfg["axes"][0], y or cfg["axes"][1], z or cfg["axes"][2]]
    for a in ax:
        if a not in rows[0]:
            raise ValueError(f"column {a!r} not in {path.name}")

    def keyval(v):
        f = _num(v)
        return (0, f, "") if math.isfinite(f) else (1, 0.0, str(v))
    axvals = [sorted({r[a] for r in rows}, key=keyval) for a in ax]
    facets = [f for f in cfg["facets"] if f in rows[0] and f not in ax]
    fvals = {f: sorted({r[f] for r in rows}, key=keyval) for f in facets}
    fvals = {f: v for f, v in fvals.items() if 1 < len(v) <= 12}
    idx = {a: {v: i for i, v in enumerate(vs)} for a, vs in zip(ax, axvals)}
    out = {"xi": [], "yi": [], "zi": [], "g2": [], "flux": [], "optical_pass": [], "device_pass": [],
           "eligible": [], "reasons": [], "model_finite_pulse": [], "model_tau_cap_density": [],
           "gates": [], "facet": {f: [] for f in fvals}}
    has_model = "model_finite_pulse" in rows[0] and "model_tau_cap_density" in rows[0]
    for r in rows:
        out["xi"].append(idx[ax[0]][r[ax[0]]]); out["yi"].append(idx[ax[1]][r[ax[1]]])
        out["zi"].append(idx[ax[2]][r[ax[2]]])
        out["g2"].append(_num(r.get(cfg["g2"])))
        out["flux"].append(_num(r.get(cfg["flux"])))
        out["optical_pass"].append(r.get("optical_pass") == "True")
        out["device_pass"].append(r.get(cfg["pass"]) == "True")
        el = r.get("eligible", r.get("eligible_row", ""))
        out["eligible"].append(el == "True")
        rs = r.get("invalid_reasons", r.get("invalid_reasons_pulsed", "")) or ""
        out["reasons"].append(rs[:140] if rs not in ("[]", "()") else "")
        if has_model:
            fp, tc = r["model_finite_pulse"] == "True", r["model_tau_cap_density"] == "True"
            out["model_finite_pulse"].append(fp)
            out["model_tau_cap_density"].append(tc)
            out["gates"].append(fp and not tc)  # only the headline model gates PASS
        else:
            out["model_finite_pulse"].append(None)
            out["model_tau_cap_density"].append(None)
            out["gates"].append(True)
        for f in fvals:
            out["facet"][f].append(fvals[f].index(r[f]))
    n = len(rows)
    gate_idx = [i for i in range(n) if out["gates"][i]]
    n_gate = len(gate_idx)
    n_pass = sum(out["device_pass"][i] for i in gate_idx)
    n_opt = sum(out["optical_pass"][i] for i in gate_idx)
    ng_rows = n - n_gate
    ng_pass = sum(out["device_pass"][i] for i in range(n) if not out["gates"][i])
    pass_word = "headline_pass" if campaign == "rt_edge" else "the device gate (device_pass)"
    gate_rows_word = "headline-model rows" if has_model else "rows"
    text = f"{n_pass} / {n_gate} {gate_rows_word} pass {pass_word}"
    if campaign != "rt_edge":
        text += f"; optical_pass {n_opt} / {n_gate}"
    if cfg["filter"]:
        n_all = len(rows_all)
        p_all = sum(r.get(cfg["pass"]) == "True" for r in rows_all)
        text += f" (all {n_all} rows in the file: {p_all} pass)"
    non_gating = None
    if ng_rows:
        non_gating = {"count": ng_pass, "rows": ng_rows,
                      "text": f"{ng_pass} / {ng_rows} non-headline-model rows pass {cfg['pass']}: "
                              "never gates (shown greyed)"}
    banner = {"text": text, "style": "pass" if n_pass > 0 else "fail", "pass": n_pass,
              "gating_rows": n_gate, "rows": n, "pass_column": cfg["pass"],
              "non_gating_passes": non_gating}
    # Default facet selection: the preferred (card) value, else the first
    # value, so each lattice cell shows one row; None would mean "all".
    fdef = {}
    for f, vs in fvals.items():
        pv = cfg.get("prefer", {}).get(f)
        fdef[f] = vs.index(pv) if pv in vs else 0
    spec = _base("lattice", campaign, campaign, f"{campaign} sweep lattice", {}, "|".join(ax) + f"|{headline}")
    spec["axes"] = {k: {"label": a, "unit": "", "scale": 1, "badge": None, "values": vs}
                    for k, a, vs in zip("xyz", ax, axvals)}
    spec["overlays"] = {"axes": {k: {"name": a, "values": vs} for k, a, vs in zip("xyz", ax, axvals)},
                        "rows": out, "facets": fvals, "facet_default": fdef, "banner": banner,
                        "counts": {"rows": n, "gating_rows": n_gate, "pass": n_pass, "optical_pass": n_opt,
                                   "non_gating_rows": ng_rows, "non_gating_pass": ng_pass,
                                   "pass_column": cfg["pass"], "by_model": model_counts},
                        "encoding": {"colour": cfg["g2"] + " (perceptual)", "size": "log10 " + cfg["flux"],
                                     "outline": "optical_pass", "solid_shell": cfg["pass"]}}
    spec["scalars"] = {"rows": _S(n, "", "A", "rows"), "gating_rows": _S(n_gate, "", "A", "gating rows"),
                       "pass": _S(n_pass, "", "A", cfg["pass"] + " (gating rows)"),
                       "optical_pass": _S(n_opt, "", "A", "optical_pass (gating rows)")}
    spec["labels"] = [text] + labels + ["no pass region is drawn unless gating rows pass"]
    if non_gating:
        spec["labels"].append(non_gating["text"])
    if model_counts:
        spec["labels"].append("all models: " + ", ".join(
            f"{'headline' if m['headline_model'] else ('fp' if m['model_finite_pulse'] else 'static')}"
            f"{'+tc' if m['model_tau_cap_density'] else ''} {m['pass']}/{m['rows']}" for m in model_counts))
    spec["provenance"] = {"calls": [str(path.relative_to(ROOT))]}
    return _clean(spec)
