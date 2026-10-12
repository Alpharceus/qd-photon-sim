"""Geometry descriptors for composed devices (3D viewport of the Composer workspace).

`describe(spec) -> dict` turns a composed spec (schema "fsim-composed/1") into a SceneSpec-style
descriptor ("fsim.scene/1", kind "device") that fsim_studio/web/js/viz3d draws with one of the
builders micropillar | open_cavity | bullseye | phc_waveguide | phc_cavity | mesa.  An unknown or
unsupported structure kind, or a spec that does not resolve, gives the fallback descriptor
(variant/builder "fallback2d"): the front end then shows the existing 2D view.

Three-layer rule (as scene.py): this module only reads the registries and the structure response's `extra` geometry (FX4), resolves the spec's leaves
(`compose.params.resolve_leaf` / `central`) and does unit-only arithmetic; the browser never
computes physics.  Numbers in `geometry` are one of

  * a resolved registry parameter (central value of a range, first alternative of a choice),
    listed in `geometry_src[key]` with `display_default: false` and the leaf's own tag, or
  * a derived display number (quarter-wave thickness lambda/4n, Gaussian-beam waist, equal-volume
    sphere ...): tag DR, `display_default: false`, `from` lists its inputs, or
  * a DISPLAY DEFAULT (`display_default: true`, tag A) where the registry has no such parameter
    (ring count of a bullseye, lattice constant of a photonic crystal, mesa height ...).  These are
    illustrative; the builders print them as [A] and the label list says so.

Configuration, not physics: the kind -> builder map, the display defaults and the world layout
live here and in the builders.  The route GET/POST /api/compose/scene is wired by the orchestrator
(api_compose.py is not part of this piece): call `describe(spec)` and return the dict.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math

from fsim_core.compose import registry as _registry
from fsim_core.compose.params import central, resolve_leaf
from fsim_core.compose.spec import T_K_SPEC

SCHEMA = "fsim.scene/1"
FALLBACK = "fallback2d"
BUILDERS = ("micropillar", "open_cavity", "bullseye", "phc_waveguide", "phc_cavity", "mesa")

# structure kind -> builder id (used when the registry entry's own `viz` is not a known builder).
KIND_BUILDER = {
    "micropillar": "micropillar",
    "open_cavity": "open_cavity",
    "cbg": "bullseye",
    "phc_waveguide": "phc_waveguide",
    "phc_cavity": "phc_cavity",
    "sil_mesa": "mesa",
    "mesa_reflector": "mesa",
    "bulk_planar": "mesa",
}

# Display defaults: (key -> (value, why)).  Each is flagged display_default in `geometry_src`.
DISPLAY = {
    "micropillar": {
        "cavity_order": (1.0, "spacer = 1 lambda/n_high (the registry says the spacer is the high-index layer)"),
        "dot_frac": (0.5, "dot at the spacer centre (antinode); the registry has no dot depth"),
        "substrate_um": (1.0, "substrate thickness is not modelled; drawn as a broken plate"),
    },
    "open_cavity": {
        "fibre_d_um": (125.0, "standard single-mode fibre cladding diameter"),
        "fibre_core_um": (5.0, "single-mode core diameter, order of magnitude"),
        "fibre_len_um": (60.0, "drawn fibre length before the break"),
        "rim_clearance_frac": (0.5, "dimple radius chosen so the end-face rim sits L/2 above the planar mirror"),
        "dbr_pairs": (6, "planar DBR pairs drawn"),
        "n_high": (3.5, "GaAs-like high-index layer"),
        "n_low": (2.95, "AlAs-like low-index layer (same defaults as the micropillar entry)"),
        "membrane_n": (3.5, "membrane index for the lambda/n membrane that holds the dot"),
        "coating_um": (0.3, "fibre-side mirror coating thickness (illustrative)"),
        "chip_halfwidth_um": (80.0, "chip half-width drawn before the break"),
    },
    "bullseye": {
        "n_eff": (2.6, "effective index setting the Bragg period, period = lambda / n_eff"),
        "n_rings": (5, "number of grating rings"),
        "disk_over_period": (1.1, "centre-disk radius in units of the period"),
        "fill": (0.5, "ring (material) fraction of one period"),
        "h_over_lambda": (0.16, "membrane thickness in units of lambda"),
        "spacer_over_lambda": (0.3, "oxide spacer thickness in units of lambda"),
        "reflector": (True, "metal back reflector under the spacer"),
    },
    "phc": {
        "a_over_lambda": (0.26, "lattice constant in units of lambda (membrane PhC class value)"),
        "r_over_a": (0.3, "hole radius in units of the lattice constant"),
        "h_over_a": (0.6, "slab thickness in units of the lattice constant"),
        "dot_r_nm": (12.0, "dot drawn radius (the registry has no dot size)"),
        "n_slab": (3.4, "slab index used to turn V (lambda/n)^3 into a length"),
    },
    "phc_waveguide": {
        "n_cols_z": (17, "holes per row along the waveguide"),
        "n_rows": (4, "rows on each side of the missing row"),
    },
    "phc_cavity": {
        "cell": ("L3", "cavity type drawn: L3 (three missing holes) or H1 (one); override with describe(display=...)"),
        "n_cols": (15, "holes per row"),
        "n_rows": (5, "rows on each side of the cavity row"),
        "end_shift_over_a": (0.15, "outward shift of the two holes at the L3 ends (a common Q optimisation)"),
        "v_aspect": (2.0, "mode ellipsoid x : y aspect for the L3 mode (H1 uses 1.3)"),
    },
    "mesa": {
        "mesa_d_um": (5.0, "mesa diameter (the registry mesa entries hold only efficiency tables)"),
        "mesa_h_um": (2.0, "mesa height"),
        "slab_d_um": (16.0, "slab width drawn before the break (bulk_planar)"),
        "slab_h_um": (3.0, "slab thickness (bulk_planar)"),
        "sil_r_over_halfwidth": (1.2, "SIL radius in units of the mesa half-width (0.6 for the slab); a real SIL is ~ mm: not to scale"),
        "dot_depth_frac": (0.5, "dot depth, fraction of the mesa height from the top"),
        "reflector_um": (0.4, "metal reflector thickness (illustrative)"),
    },
}


class _Src:
    """Collects geometry numbers and their provenance."""

    def __init__(self):
        self.geometry, self.src, self.display = {}, {}, []

    def reg(self, key, value, leaf, path, ranged=None):
        self.geometry[key] = value
        self.src[key] = {"tag": leaf.tag, "display_default": False, "kind": "registry", "path": path,
                         "unit": leaf.unit, "source": leaf.source, "defaulted": bool(leaf.defaulted)}
        if leaf.kind == "range":
            self.src[key]["range"] = list(leaf.range)
        elif leaf.kind == "choice":
            self.src[key]["choice"] = list(leaf.choice)

    def derived(self, key, value, frm, unit, note=""):
        self.geometry[key] = value
        self.src[key] = {"tag": "DR", "display_default": False, "kind": "derived", "from": list(frm),
                         "unit": unit, "note": note}

    def disp(self, group, key, name=None, value=None):
        v, why = DISPLAY[group][key]
        v = v if value is None else value
        name = name or key
        self.geometry[name] = v
        self.src[name] = {"tag": "A", "display_default": True, "kind": "display", "why": why}
        self.display.append(name)
        return v


def _resolve_block(axis, blk, entry):
    """{name: (central value, Leaf)} for every ParamSpec of the entry (registry default when omitted)."""
    out = {}
    params = (blk or {}).get("params") or {}
    unknown = set(params) - {p.name for p in entry.params}
    if unknown:
        raise ValueError(f"{axis}.params: unknown {sorted(unknown)}")
    for p in entry.params:
        raw = params[p.name] if p.name in params else p.default
        leaf = resolve_leaf(raw, p, f"{axis}.{p.name}")
        out[p.name] = (central(leaf), leaf)
    return out


def _structure_extra(spec, sent, sblk, mblk):
    """The structure's own derived geometry numbers (StructureResponse.extra, FX4): built through the registry's
    material.build -> structure.response chain at the spec's central values, so the scene never re-derives them."""
    ment = _registry.get("material", mblk.get("kind"))
    T = spec.get("T_K")
    T_K = 5.0
    try:
        T_K = float(central(resolve_leaf(T, T_K_SPEC, "T_K")))
    except Exception:
        pass
    mv = {k: v[0] for k, v in _resolve_block("material", mblk, ment).items()}
    sv = {k: v[0] for k, v in _resolve_block("structure", sblk, sent).items()}
    em = ment.build(mv, T_K)
    return dict(getattr(sent.response(sv, em, T_K, 0.75), "extra", {}) or {})


def _false_colour(lam):
    return bool(lam is not None and lam > 700.0)


def _catalog(spec):
    blob = json.dumps({a: (spec.get(a) if isinstance(spec, dict) else None) for a in ("material", "structure")},
                      sort_keys=True, default=str)
    return "FS-composed-" + hashlib.sha256(blob.encode()).hexdigest()[:6]


def _base(variant, title, spec, lam):
    return {"schema": SCHEMA, "kind": "device", "variant": variant, "builder": variant, "card": "composed",
            "title": title, "catalog_id": _catalog(spec), "units": {"length": "um"},
            "axes": {k: {"label": k, "unit": "um", "scale": 1, "badge": None} for k in "xyz"},
            "layers": [], "overlays": {"emission": {"lambda_nm": lam, "false_colour": _false_colour(lam)}},
            "scalars": {}, "labels": [], "notes": [], "provenance": {"calls": []}, "fallback": False}


def _fallback(spec, reason, structure_kind=None, lam=None):
    d = _base(FALLBACK, f"composed device ({structure_kind or 'unknown structure'})", spec, lam)
    d["fallback"] = True
    d["reason"] = reason
    d["structure_kind"] = structure_kind
    d["labels"] = [f"no 3D view: {reason}"]
    return d


def _scalar(value, unit, leaf, label, note=None):
    s = {"value": value, "unit": unit, "tag": leaf.tag if leaf is not None else "A", "label": label}
    if note:
        s["note"] = note
    return s


def _layer(name, t_nm, role, tag, source, n=None, material=""):
    return {"name": name, "t_nm": t_nm, "role": role, "tag": tag, "source": source, "n": n, "material": material}


# ---------------------------------------------------------------- per-kind descriptors

def _micropillar(vals, lam, S, ex):
    d_um = vals["d_um"][0]
    n_top, n_bot = int(vals["N_top"][0]), int(vals["N_bot"][0])     # registry: fractional counts truncated
    nh, nl = vals["n_high"][0], vals["n_low"][0]
    S.reg("d_um", d_um, vals["d_um"][1], "structure.d_um")
    S.reg("N_top", n_top, vals["N_top"][1], "structure.N_top")
    S.reg("N_bot", n_bot, vals["N_bot"][1], "structure.N_bot")
    S.reg("n_high", nh, vals["n_high"][1], "structure.n_high")
    S.reg("n_low", nl, vals["n_low"][1], "structure.n_low")
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    t_h, t_l = float(ex["t_high_nm"]), float(ex["t_low_nm"])       # from the structure response (FX4)
    S.derived("t_high_nm", t_h, ["material.lambda_nm", "structure.n_high"], "nm", "quarter wave lambda/4n (structure response)")
    S.derived("t_low_nm", t_l, ["material.lambda_nm", "structure.n_low"], "nm", "quarter wave lambda/4n (structure response)")
    order = S.disp("micropillar", "cavity_order")
    spacer = order * float(ex["spacer_unit_nm"])
    S.derived("spacer_nm", spacer, ["material.lambda_nm", "structure.n_high"], "nm", "cavity order x lambda/n_high")
    S.disp("micropillar", "dot_frac")
    sub = S.disp("micropillar", "substrate_um")
    top_nm, bot_nm = n_top * (t_h + t_l), n_bot * (t_h + t_l)
    total = top_nm + bot_nm + spacer
    S.derived("height_um", total / 1000.0, ["structure.N_top", "structure.N_bot", "structure.n_high",
                                           "structure.n_low", "material.lambda_nm"], "um", "pillar height (no substrate)")
    layers = [_layer("substrate", sub * 1000.0, "substrate", "A", "thickness not modelled"),
              _layer(f"bottom DBR ({n_bot} pairs)", bot_nm, "mirror", "DR", "N_bot x (lambda/4n_h + lambda/4n_l)"),
              _layer("spacer", spacer, "core", "DR", "lambda/n_high"),
              _layer(f"top DBR ({n_top} pairs)", top_nm, "mirror", "DR", "N_top x (lambda/4n_h + lambda/4n_l)")]
    labels = [f"quarter-wave layers lambda/4n: {t_h:.1f} nm (high), {t_l:.1f} nm (low), true scale",
              "substrate thickness not modelled", "dot depth in the spacer is a display default [A]"]
    return layers, labels, ["d_um", "N_top", "N_bot", "n_high", "n_low", "lambda_nm"]


def _open_cavity(vals, lam, S, ex):
    roc, L = vals["ROC_um"][0], vals["L_um"][0]
    S.reg("ROC_um", roc, vals["ROC_um"][1], "structure.ROC_um")
    S.reg("L_um", L, vals["L_um"][1], "structure.L_um")
    S.reg("mode", vals["mode"][0], vals["mode"][1], "structure.mode")
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    fd = S.disp("open_cavity", "fibre_d_um")
    S.disp("open_cavity", "fibre_core_um")
    S.disp("open_cavity", "fibre_len_um")
    pairs = S.disp("open_cavity", "dbr_pairs")
    nh, nl = S.disp("open_cavity", "n_high"), S.disp("open_cavity", "n_low")
    mem_n = S.disp("open_cavity", "membrane_n")
    S.disp("open_cavity", "coating_um")
    S.disp("open_cavity", "chip_halfwidth_um")
    frac = S.disp("open_cavity", "rim_clearance_frac")
    lam_um = lam / 1000.0
    stable = roc > L > 0
    S.geometry["stable"] = bool(stable)
    S.src["stable"] = {"tag": "DR", "display_default": False, "kind": "derived", "from": ["structure.ROC_um", "structure.L_um"],
                       "unit": "", "note": "plano-concave stability 0 < L < ROC"}
    if stable and ex.get("waist_um") is not None:
        w0 = float(ex["waist_um"])                    # from the structure response (FX4)
        zr = math.pi * w0 * w0 / lam_um
        S.derived("waist_um", w0, ["material.lambda_nm", "structure.L_um", "structure.ROC_um"], "um",
                  "TEM00 waist on the planar mirror, w0^2 = (lambda/pi) sqrt(L (ROC - L)) (structure response)")
        S.derived("rayleigh_um", zr, ["waist_um", "material.lambda_nm"], "um", "z_R = pi w0^2 / lambda")
    else:
        S.geometry["waist_um"] = None
        S.geometry["rayleigh_um"] = None
        S.src["waist_um"] = {"tag": "DR", "display_default": False, "kind": "derived", "unit": "um",
                             "note": "no stable mode: L >= ROC"}
    # dimple: spherical cap of radius ROC whose rim is `frac` x L closer to the chip than the apex
    sag = min(frac * L, roc * 0.9)
    a = math.sqrt(max(roc * roc - (roc - sag) ** 2, 0.0))
    S.derived("sag_um", sag, ["structure.L_um", "structure.ROC_um"], "um", "cap depth = rim_clearance_frac x L (display rule)")
    S.derived("dimple_radius_um", a, ["structure.ROC_um", "sag_um"], "um", "a^2 = ROC^2 - (ROC - sag)^2")
    S.derived("t_high_nm", lam / (4 * nh), ["material.lambda_nm"], "nm", "quarter wave, display index")
    S.derived("t_low_nm", lam / (4 * nl), ["material.lambda_nm"], "nm", "quarter wave, display index")
    S.derived("membrane_nm", lam / mem_n, ["material.lambda_nm"], "nm", "lambda/n membrane holding the dot")
    layers = [_layer("planar DBR", pairs * lam * (1 / (4 * nh) + 1 / (4 * nl)), "mirror", "A", "display pairs"),
              _layer("membrane", lam / mem_n, "core", "A", "lambda/n display"),
              _layer("vacuum gap L", L * 1000.0, "cladding", vals["L_um"][1].tag, "structure.L_um")]
    labels = ["fibre diameter, core, drawn length, planar-DBR pairs and the dimple radius are display defaults [A]",
              f"ROC {roc:g} um, cavity length {L:g} um (to the dimple apex), true scale",
              "mode drawn: TEM00 Gaussian envelope w(z) in the gap" if stable else "unstable cavity (L >= ROC): no mode drawn"]
    return layers, labels, ["ROC_um", "L_um", "waist_um", "lambda_nm"]


def _bullseye(vals, lam, S):
    S.reg("F_P", vals["F_P"][0], vals["F_P"][1], "structure.F_P")
    S.reg("Q", vals["Q"][0], vals["Q"][1], "structure.Q")
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    n_eff = S.disp("bullseye", "n_eff")
    nr = S.disp("bullseye", "n_rings")
    dop = S.disp("bullseye", "disk_over_period")
    fill = S.disp("bullseye", "fill")
    h_l = S.disp("bullseye", "h_over_lambda")
    sp_l = S.disp("bullseye", "spacer_over_lambda")
    S.disp("bullseye", "reflector")
    period = lam / n_eff
    S.derived("period_nm", period, ["material.lambda_nm", "n_eff"], "nm", "Bragg period lambda / n_eff")
    S.derived("disk_radius_nm", dop * period, ["period_nm", "disk_over_period"], "nm")
    S.geometry["disk_over_period"] = dop
    S.derived("ring_width_nm", fill * period, ["period_nm", "fill"], "nm")
    S.derived("membrane_nm", h_l * lam, ["material.lambda_nm", "h_over_lambda"], "nm")
    S.derived("spacer_nm", sp_l * lam, ["material.lambda_nm", "spacer_over_lambda"], "nm")
    S.derived("trench_nm", (1.0 - fill) * period, ["period_nm", "fill"], "nm", "etched gap before each ring")
    S.derived("outer_radius_um", (dop * period + nr * period) / 1000.0,
              ["period_nm", "n_rings", "disk_over_period"], "um", "centre disk + N periods")
    layers = [_layer("back reflector", 100.0, "contact", "A", "display"),
              _layer("oxide spacer", sp_l * lam, "oxide", "A", "display"),
              _layer("membrane (rings)", h_l * lam, "core", "A", "display")]
    labels = ["ring count, period, centre disk, fill, membrane and spacer are display defaults [A]: the registry holds "
              "F_P, Q and an efficiency table only"]
    return layers, labels, ["F_P", "Q", "period_nm", "lambda_nm"]


def _phc_lattice(S, lam):
    a = S.disp("phc", "a_over_lambda") * lam
    S.geometry["a_nm"] = a
    S.src["a_nm"] = {"tag": "A", "display_default": True, "kind": "display", "why": "lattice constant = a/lambda x lambda"}
    S.display.append("a_nm")
    S.geometry["r_nm"] = S.disp("phc", "r_over_a", "r_over_a") * a
    S.src["r_nm"] = {"tag": "A", "display_default": True, "kind": "display", "why": "hole radius = r/a x a"}
    S.display.append("r_nm")
    S.geometry["h_nm"] = S.disp("phc", "h_over_a", "h_over_a") * a
    S.src["h_nm"] = {"tag": "A", "display_default": True, "kind": "display", "why": "slab thickness = h/a x a"}
    S.display.append("h_nm")
    S.disp("phc", "dot_r_nm")
    return a


def _phc_waveguide(vals, lam, S):
    ng, aeff = vals["n_g"][0], vals["A_eff_um2"][0]
    S.reg("n_g", ng, vals["n_g"][1], "structure.n_g")
    S.reg("A_eff_um2", aeff, vals["A_eff_um2"][1], "structure.A_eff_um2")
    S.reg("eta_coupler", vals["eta_coupler"][0], vals["eta_coupler"][1], "structure.eta_coupler")
    S.reg("eta_prop", vals["eta_prop"][0], vals["eta_prop"][1], "structure.eta_prop")
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    a = _phc_lattice(S, lam)
    S.disp("phc_waveguide", "n_cols_z")
    S.disp("phc_waveguide", "n_rows")
    S.derived("mode_radius_nm", 1000.0 * math.sqrt(aeff / math.pi), ["structure.A_eff_um2"], "nm",
              "radius of the circle with the effective mode area A_eff")
    S.geometry["missing_rows"] = 1
    S.src["missing_rows"] = {"tag": "A", "display_default": True, "kind": "display", "why": "W1: one missing row"}
    S.display.append("missing_rows")
    layers = [_layer("PhC slab", S.geometry["h_nm"], "core", "A", "display a, r/a, h/a")]
    labels = ["lattice constant a, hole radius, slab thickness, row/column counts and the dot size are display defaults [A]",
              "guided-mode tube radius from A_eff (energy-normalised area, equal-area circle)"]
    return layers, labels, ["n_g", "A_eff_um2", "eta_coupler", "eta_prop", "lambda_nm"]


def _phc_cavity(vals, lam, S, display):
    q, v = vals["Q"][0], vals["V_lambda3"][0]
    S.reg("Q", q, vals["Q"][1], "structure.Q")
    S.reg("V_lambda3", v, vals["V_lambda3"][1], "structure.V_lambda3")
    S.reg("g_ueV", vals["g_ueV"][0], vals["g_ueV"][1], "structure.g_ueV")
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    a = _phc_lattice(S, lam)
    n_slab = S.disp("phc", "n_slab")
    cell = S.disp("phc_cavity", "cell", value=display.get("cell"))
    if cell not in ("L3", "H1"):
        raise ValueError(f"display cell must be 'L3' or 'H1', got {cell!r}")
    S.disp("phc_cavity", "n_cols")
    S.disp("phc_cavity", "n_rows")
    S.disp("phc_cavity", "end_shift_over_a")
    # equal-volume ellipsoid of V (lambda/n)^3
    lam_n_um = lam / 1000.0 / n_slab
    vol = v * lam_n_um ** 3
    asp = S.disp("phc_cavity", "v_aspect", value=None if cell == "L3" else 1.3)
    # semi-axes (x, y, z) = (asp, 1, 0.6) k with 4/3 pi (asp*1*0.6) k^3 = V
    k = (vol * 3.0 / (4.0 * math.pi * asp * 0.6)) ** (1.0 / 3.0)
    S.derived("mode_semiaxes_um", [asp * k, k, 0.6 * k], ["structure.V_lambda3", "material.lambda_nm", "n_slab"], "um",
              "ellipsoid of volume V (lambda/n)^3 with display aspect (asp, 1, 0.6)")
    layers = [_layer("PhC slab", S.geometry["h_nm"], "core", "A", "display a, r/a, h/a")]
    labels = [f"{cell} cavity drawn (a display choice: the registry entry covers L3, H1 and nanobeam)",
              "lattice constant a, hole radius, slab thickness, counts and the dot size are display defaults [A]",
              "mode glow: ellipsoid of volume V (lambda/n)^3, n a display default"]
    return layers, labels, ["Q", "V_lambda3", "lambda_nm"]


def _mesa(kind, vals, lam, S):
    S.derived("lambda_nm", lam, ["material.lambda_nm"], "nm")
    if kind == "bulk_planar":
        n = vals["n"][0]
        S.reg("n", n, vals["n"][1], "structure.n")
        S.reg("sil", vals["sil"][0], vals["sil"][1], "structure.sil")
        S.reg("n_s", vals["n_s"][0], vals["n_s"][1], "structure.n_s")
        S.reg("mirror", bool(vals["mirror"][0]), vals["mirror"][1], "structure.mirror")
        mode, sil, mirror = "slab", vals["sil"][0], bool(vals["mirror"][0])
        d = S.disp("mesa", "slab_d_um")
        h = S.disp("mesa", "slab_h_um")
        S.geometry["d_um"], S.geometry["h_um"] = d, h
        S.src["d_um"] = {"tag": "A", "display_default": True, "kind": "display", "why": "slab width drawn before the break"}
        S.src["h_um"] = {"tag": "A", "display_default": True, "kind": "display", "why": "slab thickness"}
        S.display += ["d_um", "h_um"]
    else:
        mode = "mesa"
        sil = "hemi" if kind == "sil_mesa" else "none"
        mirror = kind == "mesa_reflector"
        d = S.disp("mesa", "mesa_d_um", "d_um")
        h = S.disp("mesa", "mesa_h_um", "h_um")
        S.geometry["n_s"] = 2.0
        S.src["n_s"] = {"tag": "A", "display_default": True, "kind": "display", "why": "SIL index used for the Weierstrass height"}
        S.display.append("n_s")
    S.geometry["mode"], S.geometry["sil"], S.geometry["mirror"] = mode, sil, mirror
    for key, why in (("mode", "mesa pillar or wide slab"), ("sil", "solid-immersion lens"), ("mirror", "metal back reflector")):
        if key not in S.src:
            S.src[key] = {"tag": "A", "display_default": True, "kind": "display", "why": f"{why} (set by the entry kind {kind})"}
            S.display.append(key)
    ratio = S.disp("mesa", "sil_r_over_halfwidth", value=0.6 if mode == "slab" else None)
    S.derived("sil_radius_um", ratio * d / 2.0 if sil != "none" else 0.0, ["d_um", "sil_r_over_halfwidth"], "um",
              "display rule: lens radius = ratio x half-width (a real SIL is ~mm)")
    ns = S.geometry["n_s"]
    S.derived("sil_height_um", (0.0 if sil == "none" else (S.geometry["sil_radius_um"] if sil == "hemi"
                                                          else S.geometry["sil_radius_um"] * (1.0 + 1.0 / ns))),
              ["sil_radius_um", "n_s"], "um", "hemisphere R; Weierstrass R (1 + 1/n_s)")
    S.disp("mesa", "dot_depth_frac")
    S.disp("mesa", "reflector_um")
    layers = []
    if mirror:
        layers.append(_layer("metal reflector", S.geometry["reflector_um"] * 1000.0, "contact", "A", "display"))
    layers.append(_layer("slab" if mode == "slab" else "mesa", h * 1000.0, "core", "A", "display"))
    labels = ["mesa/slab size, SIL radius, dot depth and reflector thickness are display defaults [A]",
              f"SIL not to scale: a real lens is ~mm across, drawn at {ratio:g} x the {mode} half-width"]
    return layers, labels, ["lambda_nm"]


# ---------------------------------------------------------------- entry point

def describe(spec, display=None) -> dict:
    """Geometry descriptor for a composed spec (dict, schema "fsim-composed/1").

    display: optional overrides of display defaults, {"cell": "H1"} (phc_cavity only).
    Never raises for an unknown structure kind or an unresolvable leaf: returns the fallback descriptor."""
    display = dict(display or {})
    try:
        _registry.ensure_loaded()
        if not isinstance(spec, dict):
            return _fallback({}, "spec is not a mapping")
        sblk = spec.get("structure") or {}
        skind = sblk.get("kind")
        try:
            sent = _registry.get("structure", skind)
        except Exception:
            return _fallback(spec, f"unknown structure kind {skind!r}", skind)
        builder = sent.viz if sent.viz in BUILDERS else KIND_BUILDER.get(skind, FALLBACK)
        # emission wavelength: the material's resolved lambda_nm
        lam = None
        mblk = spec.get("material") or {}
        try:
            ment = _registry.get("material", mblk.get("kind"))
            mv = _resolve_block("material", mblk, ment)
            lam = float(mv["lambda_nm"][0]) if "lambda_nm" in mv else None
        except Exception:
            lam = None
        if builder == FALLBACK:
            return _fallback(spec, f"no 3D builder for structure kind {skind!r}", skind, lam)
        if lam is None or not math.isfinite(lam) or lam <= 0:
            return _fallback(spec, "material has no resolvable emission wavelength", skind, lam)
        vals = _resolve_block("structure", sblk, sent)
        S = _Src()
        if builder == "micropillar":
            layers, labels, keys = _micropillar(vals, lam, S, _structure_extra(spec, sent, sblk, mblk))
        elif builder == "open_cavity":
            layers, labels, keys = _open_cavity(vals, lam, S, _structure_extra(spec, sent, sblk, mblk))
        elif builder == "bullseye":
            layers, labels, keys = _bullseye(vals, lam, S)
        elif builder == "phc_waveguide":
            layers, labels, keys = _phc_waveguide(vals, lam, S)
        elif builder == "phc_cavity":
            layers, labels, keys = _phc_cavity(vals, lam, S, display)
        else:
            layers, labels, keys = _mesa(skind, vals, lam, S)
    except Exception as e:          # an unresolved spec must not break the viewport
        return _fallback(spec if isinstance(spec, dict) else {}, f"{type(e).__name__}: {e}",
                         (spec.get("structure") or {}).get("kind") if isinstance(spec, dict) else None)
    d = _base(builder, f"{sent.label} - composed", spec, lam)
    d["structure_kind"] = skind
    d["material_kind"] = mblk.get("kind")
    d["geometry"] = copy.deepcopy(S.geometry)
    d["geometry_src"] = copy.deepcopy(S.src)
    d["display_defaults"] = sorted(set(S.display))
    d["layers"] = layers
    d["labels"] = labels + ["display defaults are illustrative [A]; registry and derived numbers are true scale"]
    scal = {}
    for k in keys:
        src = S.src.get(k, {})
        unit = src.get("unit", "")
        scal[k] = {"value": S.geometry.get(k), "unit": unit, "tag": src.get("tag", "A"), "label": k.replace("_", " "),
                   "display_default": bool(src.get("display_default"))}
    d["scalars"] = scal
    d["readout_keys"] = keys
    d["provenance"] = {"calls": ["fsim_core.compose.registry", "fsim_core.compose.params.resolve_leaf/central"],
                       "structure_kind": skind, "material_kind": mblk.get("kind")}
    return d
