"""Checks for fsim_studio/scene_compose.py (P13): descriptor numbers equal the composed spec's resolved values.

Reference values come from fsim_core.compose.spec.validate_spec (the validator's own `resolved` table) and
closed forms (quarter wave lambda/4n, Gaussian-beam waist, equal-area circle, escape cone, nanowire V number and
mirror relations, equal-volume ellipsoid) or from the physics modules the adapters wrap (dbr, waveguide,
nitride_cavity), not from scene_compose.  R2-3D adds the coverage gate: EVERY registered structure kind and EVERY
card in cards/papers (the _fixture included: it draws the lumped schematic) and cards/composed has a non-fallback
descriptor.  Run: python verify/verify_scene_compose.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import yaml                                                 # noqa: E402

from fsim_core import dbr, waveguide, nitride_cavity as _nc, nitride_nanowire_photonics as _nwp   # noqa: E402
from fsim_core.compose import registry as _reg              # noqa: E402
from fsim_core.compose.spec import validate_spec          # noqa: E402
from fsim_studio import api_compose as _api                # noqa: E402
from fsim_studio.scene_compose import describe, BUILDERS, KIND_BUILDER, FALLBACK   # noqa: E402

RES = []


def check(name, ok, detail=""):
    RES.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def close(a, b, rel=1e-12):
    return a is not None and b is not None and abs(a - b) <= rel * max(1.0, abs(b))


def spec(skind, params=None, mat="ingaas_gaas", mparams=None):
    return {"schema": "fsim-composed/1", "T_K": 4.2, "material": {"kind": mat, "params": mparams or {}},
            "structure": {"kind": skind, "params": params or {}}, "drive": {"kind": "x", "params": {}},
            "setup": {"kind": "y", "params": {}}}


def val(s, path):
    e = validate_spec(s)["resolved"][path]
    if "value" in e:
        return e["value"]
    if "range" in e:
        lo, hi = e["range"]
        return math.sqrt(lo * hi) if e.get("scale") == "log" else 0.5 * (lo + hi)
    return e["choice"][0]


# --- micropillar: resolved params, quarter-wave layers, height
s = spec("micropillar", {"d_um": 3.0, "N_top": {"value": 12, "unit": "1"}, "n_low": 3.0}, mparams={"lambda_nm": 900.0})
d = describe(s)
g = d["geometry"]
lam, nh, nl = val(s, "material.lambda_nm"), val(s, "structure.n_high"), val(s, "structure.n_low")
check("micropillar: builder id", d["variant"] == "micropillar" and not d["fallback"])
check("micropillar: d, N_top, N_bot, n equal resolved", g["d_um"] == val(s, "structure.d_um") == 3.0 and g["N_top"] == 12
      and g["N_bot"] == int(val(s, "structure.N_bot")) and g["n_high"] == nh and g["n_low"] == nl == 3.0)
check("micropillar: quarter-wave layers lambda/4n", close(g["t_high_nm"], lam / (4 * nh)) and close(g["t_low_nm"], lam / (4 * nl)),
      f"{g['t_high_nm']:.3f}/{g['t_low_nm']:.3f} nm")
tot = g["N_top"] * (lam / (4 * nh) + lam / (4 * nl)) + g["N_bot"] * (lam / (4 * nh) + lam / (4 * nl)) + lam / nh
check("micropillar: height = DBRs + lambda/n spacer", close(g["height_um"] * 1000, tot), f"{g['height_um']:.4f} um")
check("micropillar: emission wavelength is the material's", d["overlays"]["emission"]["lambda_nm"] == lam == 900.0)
check("micropillar: display defaults flagged, registry numbers not",
      g_src_ok := (d["geometry_src"]["dot_frac"]["display_default"] is True and d["geometry_src"]["d_um"]["display_default"] is False))
check("micropillar: JSON-safe", json.loads(json.dumps(d)) == d)

# --- range leaf -> central value (lin midpoint)
s = spec("micropillar", {"d_um": {"range": [2.0, 3.0], "unit": "um"}})
check("micropillar: range leaf drawn at its midpoint, range recorded", describe(s)["geometry"]["d_um"] == 2.5
      and describe(s)["geometry_src"]["d_um"]["range"] == [2.0, 3.0])

# --- open cavity: ROC, L, waist closed form, stability
s = spec("open_cavity", {"ROC_um": 25.0, "L_um": 3.0})
d = describe(s)
g = d["geometry"]
lam = val(s, "material.lambda_nm")
w0 = math.sqrt(lam / 1000 / math.pi * math.sqrt(3.0 * 22.0))
check("open_cavity: ROC, L equal resolved; waist closed form", g["ROC_um"] == 25.0 and g["L_um"] == 3.0 and close(g["waist_um"], w0),
      f"w0 {g['waist_um']:.4f} um")
check("open_cavity: Rayleigh range pi w0^2/lambda", close(g["rayleigh_um"], math.pi * w0 ** 2 / (lam / 1000)))
check("open_cavity: dimple fits the ROC (a^2 = ROC^2 - (ROC - sag)^2)",
      close(g["dimple_radius_um"] ** 2, 25.0 ** 2 - (25.0 - g["sag_um"]) ** 2, 1e-9) and g["sag_um"] == 0.5 * 3.0)
d2 = describe(spec("open_cavity", {"ROC_um": 2.0, "L_um": 3.0}))
check("open_cavity: unstable (L >= ROC) -> stable False, no waist", d2["geometry"]["stable"] is False and d2["geometry"]["waist_um"] is None)

# --- bullseye
s = spec("cbg", {"F_P": 12.0})
d = describe(s)
g = d["geometry"]
check("bullseye: F_P, Q equal resolved; builder", d["variant"] == "bullseye" and g["F_P"] == 12.0 and g["Q"] == val(s, "structure.Q"))
check("bullseye: period = lambda / n_eff, outer radius = disk + N periods",
      close(g["period_nm"], g["lambda_nm"] / g["n_eff"]) and
      close(g["outer_radius_um"] * 1000, g["disk_radius_nm"] + g["n_rings"] * g["period_nm"]))
check("bullseye: ring geometry all display defaults", all(d["geometry_src"][k]["display_default"] for k in ("n_eff", "n_rings", "fill")))

# --- PhC waveguide / cavity
s = spec("phc_waveguide", {"A_eff_um2": 0.05, "n_g": 20.0})
d = describe(s)
g = d["geometry"]
check("phc_waveguide: n_g, A_eff equal resolved; tube radius = equal-area circle",
      g["n_g"] == 20.0 and g["A_eff_um2"] == 0.05 and close(g["mode_radius_nm"], 1000 * math.sqrt(0.05 / math.pi)))
check("phc_waveguide: one missing row, r/a/h from display defaults",
      g["missing_rows"] == 1 and close(g["r_nm"], g["a_nm"] * g["r_over_a"]) and close(g["a_nm"], g["a_over_lambda"] * g["lambda_nm"]))
s = spec("phc_cavity", {"Q": 5e4, "V_lambda3": 1.5})
d = describe(s)
g = d["geometry"]
ax, ay, az = g["mode_semiaxes_um"]
vol = 4 / 3 * math.pi * ax * ay * az
check("phc_cavity: Q, V equal resolved; ellipsoid volume = V (lambda/n)^3",
      g["Q"] == 5e4 and g["V_lambda3"] == 1.5 and close(vol, 1.5 * (g["lambda_nm"] / 1000 / g["n_slab"]) ** 3, 1e-9))
check("phc_cavity: default L3, override H1", g["cell"] == "L3" and describe(s, display={"cell": "H1"})["geometry"]["cell"] == "H1")
check("phc_cavity: bad display cell -> fallback, no raise", describe(s, display={"cell": "X"})["fallback"] is True)

# --- mesa / SIL / slab
s = spec("bulk_planar", {"n": 3.2, "sil": "weierstrass", "n_s": 1.8, "mirror": True})
d = describe(s)
g = d["geometry"]
check("bulk_planar: slab with resolved n, SIL, n_s, mirror",
      d["variant"] == "mesa" and g["mode"] == "slab" and g["n"] == 3.2 and g["sil"] == "weierstrass" and g["n_s"] == 1.8 and g["mirror"] is True)
check("bulk_planar: Weierstrass height R (1 + 1/n_s)", close(g["sil_height_um"], g["sil_radius_um"] * (1 + 1 / 1.8)))
d = describe(spec("sil_mesa"))
check("sil_mesa: mesa + hemispherical SIL, height = radius", d["geometry"]["mode"] == "mesa" and d["geometry"]["sil"] == "hemi"
      and close(d["geometry"]["sil_height_um"], d["geometry"]["sil_radius_um"]))
d = describe(spec("mesa_reflector"))
check("mesa_reflector: reflector, no SIL", d["geometry"]["mirror"] is True and d["geometry"]["sil"] == "none" and d["geometry"]["sil_radius_um"] == 0.0)

# --- fallback descriptors
d = describe(spec("nope"))
check("nope: fallback descriptor (unknown structure)", d["fallback"] is True and d["variant"] == "fallback2d" and "unknown structure" in d["reason"], d.get("reason", ""))
d = describe(spec("edge_ridge"))
check("edge_ridge on a 950 nm dot: stack out of the index table -> fallback with the reason, no raise",
      d["fallback"] is True and d["variant"] == "fallback2d" and "tabulated" in d["reason"], d.get("reason", "")[:70])
check("non-mapping spec -> fallback, no raise", describe("x")["fallback"] is True)
check("bad leaf (unit mismatch) -> fallback, no raise",
      describe(spec("micropillar", {"d_um": {"value": 2.5, "unit": "nm"}}))["fallback"] is True)
check("descriptor layers present for the 2D fallback view",
      all(describe(spec(k))["layers"] for k in ("micropillar", "open_cavity", "cbg", "phc_waveguide", "phc_cavity", "sil_mesa")))

# ======================================================================= R2-3D
MAT = {"edge_ridge": "inp_gainp_red", "nitride_planar_cavity": "ingan_gan"}      # materials whose lambda the kind supports

# --- (a) every registered structure kind maps to a builder other than fallback2d
kinds = list(_reg.kinds("structure"))
check("coverage (a): the registry lists the 13 structure kinds", len(kinds) == 13, ",".join(kinds))
missing = [k for k in kinds if (_reg.get("structure", k).viz if _reg.get("structure", k).viz in BUILDERS else KIND_BUILDER.get(k, FALLBACK)) == FALLBACK]
check("coverage (a): every structure kind maps to a builder other than fallback2d", not missing, f"unmapped: {missing}")
bad = []
for k in kinds:
    dd = describe(spec(k, mat=MAT.get(k, "ingaas_gaas")))
    if dd["fallback"] or dd["variant"] == FALLBACK or dd["variant"] not in BUILDERS or not dd["layers"] or not dd["geometry"] or not dd["labels"]:
        bad.append((k, dd.get("reason", dd["variant"])))
check("coverage (a): describe() of every structure kind (default params) is a non-fallback descriptor with geometry and layers", not bad, str(bad))

# --- (b) every card produces a non-fallback descriptor (the _fixture is abstract: it draws the lumped schematic)
bad, n_cards = [], 0
for f in sorted((ROOT / "cards" / "papers").glob("*.yaml")):
    doc = yaml.safe_load(f.read_text(encoding="utf-8"))
    sp = dict(_api.central_design(doc["meta"]["id"], doc["device"])["design"]["composed"], schema="fsim-composed/1")
    dd = describe(sp, _api.display_hint(sp))
    n_cards += 1
    if dd["fallback"] or dd["variant"] not in BUILDERS:
        bad.append((f.stem, dd.get("reason", "")))
for f in sorted((ROOT / "cards" / "composed").glob("*.yaml")):
    doc = yaml.safe_load(f.read_text(encoding="utf-8"))
    sp = dict(doc["design"]["composed"], schema="fsim-composed/1")
    dd = describe(sp, _api.display_hint(sp))
    n_cards += 1
    if dd["fallback"] or dd["variant"] not in BUILDERS:
        bad.append((f.stem, dd.get("reason", "")))
check(f"coverage (b): all {n_cards} cards (papers incl. _fixture, composed) give a non-fallback descriptor", n_cards >= 46 and not bad, str(bad))

# --- (c) photonic_nanowire: numbers equal the resolved params and closed forms
s = spec("photonic_nanowire", {"d_nm": 250.0, "n_wire": 3.4, "mirror_R": 0.8, "mirror_phase_ok": True})
d = describe(s)
g = d["geometry"]
lam = val(s, "material.lambda_nm")
V = 2 * math.pi * 125.0 / lam * math.sqrt(3.4 ** 2 - 1.0)
check("nanowire (c): builder, d_nm/n_wire/mirror_R/mirror_phase_ok equal the spec", d["variant"] == "photonic_nanowire" and g["d_nm"] == 250.0
      and g["n_wire"] == 3.4 and g["mirror_R"] == 0.8 and g["mirror_phase_ok"] is True)
check("nanowire (c): V number = 2 pi a/lambda sqrt(n^2 - 1) (closed form)", close(g["V_number"], V, 1e-6), f"{g['V_number']:.5f} vs {V:.5f}")
check("nanowire (c): field |r| = sqrt(R)", close(g["mirror_field_r"], math.sqrt(0.8)))
bt, rb = g["beta_he11"], g["mirror_field_r"]
check("nanowire (c): upward fraction (beta/2)(1+|r|)^2/(1+beta|r|) and rate boost 1 + beta|r| (Munsch 2013 Eq. 1)",
      close(g["upward_fraction"], 0.5 * bt * (1 + rb) ** 2 / (1 + bt * rb)) and close(g["rate_boost"], 1 + bt * rb))
check("nanowire (c): needle: not a trumpet, top diameter = d, mode radius = the bare-wire Marcuse w", g["trumpet"] is False and g["top_d_nm"] == 250.0
      and g["flare_len_um"] == 0.0 and close(g["mode_radius_nm"], _nwp._marcuse_w_over_a(V) * 125.0, 1e-6), f"w {g['mode_radius_nm']:.2f} nm")
check("nanowire (c): antinode dot height lambda/(2 n_wire), flagged display default",
      close(g["dot_height_um"], 0.5 * lam / 3.4 / 1000.0) and d["geometry_src"]["dot_height_um"]["display_default"] is True)
s = spec("photonic_nanowire", {"taper_mfr_nm": 900.0})
d = describe(s)
g = d["geometry"]
want_flare = min(0.6 * 10.0, (900.0 - 110.0) / 1000.0 / math.tan(math.radians(5.0)))
r_top = 110.0 + want_flare * 1000.0 * math.tan(math.radians(5.0))
check("nanowire (c): trumpet: taper_mfr_nm = 900 -> trumpet, mode radius = the taper mfr, top radius / flare length from the display rule",
      g["trumpet"] is True and g["mode_radius_nm"] == 900.0 and close(g["flare_len_um"], want_flare) and close(g["top_d_nm"], 2 * r_top),
      f"flare {g['flare_len_um']:.3f} um, top Ø {g['top_d_nm']:.1f} nm")
check("nanowire (c): phase unknown -> dot height 0.3 of the wire, incoherent upward fraction (beta/2)(1+R)",
      close(g["dot_height_um"], 3.0) and close(g["upward_fraction"], 0.5 * g["beta_he11"] * (1 + 0.9)))

# --- (c) planar_dbr
s = spec("planar_dbr", {"N_top": 12, "N_bot": 18, "n_cavity": 3.3, "z_dot_frac": 0.25})
d = describe(s)
g = d["geometry"]
nh, nl = val(s, "structure.n_high"), val(s, "structure.n_low")
lam = val(s, "material.lambda_nm")
tot = 12 * (lam / (4 * nh) + lam / (4 * nl)) + 18 * (lam / (4 * nh) + lam / (4 * nl)) + lam / 3.3
check("planar_dbr (c): pair counts, n_cavity, z_dot_frac equal the spec", d["variant"] == "planar_dbr" and g["N_top"] == 12 and g["N_bot"] == 18
      and g["n_cavity"] == 3.3 and g["z_dot_frac"] == 0.25)
check("planar_dbr (c): lambda/4n layers, lambda/n_cavity spacer, dot depth, total height", close(g["t_high_nm"], lam / (4 * nh)) and close(g["t_low_nm"], lam / (4 * nl))
      and close(g["spacer_nm"], lam / 3.3) and close(g["z_dot_nm"], 0.25 * lam / 3.3) and close(g["height_um"] * 1000, tot), f"{g['height_um']:.4f} um")
Tt = 1.0 - dbr.qw_peak_reflectivity(1.0, 3.3, nh, nl, 12)
Tb = 1.0 - dbr.qw_peak_reflectivity(val(s, "structure.n_substrate"), 3.3, nh, nl, 18)
check("planar_dbr (c): exit asymmetry eta_out = T_top/(T_top + T_bottom) from dbr.qw_peak_reflectivity", close(g["eta_out"], Tt / (Tt + Tb), 1e-9), f"{g['eta_out']:.4f}")
check("planar_dbr (c): 1-D model flag, no mode volume", g["mode"] == "planar_1d")

# --- (c) nitride_planar_cavity
s = spec("nitride_planar_cavity", {"mode_volume_norm": 3.0, "Q": 4000.0}, mat="ingan_gan")
d = describe(s)
g = d["geometry"]
lam = val(s, "material.lambda_nm")
a, _, c = g["mode_semiaxes_um"]
check("nitride (c): builder, Q and V equal the spec", d["variant"] == "nitride_cavity" and g["Q"] == 4000.0 and g["mode_volume_norm"] == 3.0)
check("nitride (c): mode ellipsoid volume = V (lambda/n)^3 (closed form), pairs/indices flagged display defaults",
      close(4 / 3 * math.pi * a * a * c, 3.0 * (lam / 1000 / g["n_high"]) ** 3, 1e-9)
      and all(d["geometry_src"][k]["display_default"] for k in ("N_top", "N_bot", "n_high", "n_low")))
check("nitride (c): quarter-wave layers lambda/4n of the display indices; spacer lambda/n_high",
      close(g["t_high_nm"], lam / (4 * g["n_high"])) and close(g["t_low_nm"], lam / (4 * g["n_low"])) and close(g["spacer_nm"], lam / g["n_high"]))

# --- (c) edge_ridge (inp_gainp_red, 670 nm)
s = spec("edge_ridge", {"ridge_width_nm": 1800.0, "etch_depth_nm": 900.0, "R_back": -1.0}, mat="inp_gainp_red")
d = describe(s)
g = d["geometry"]
lam = val(s, "material.lambda_nm")
stack = waveguide.hkust_ridge_stack(lam)
mode = waveguide.effective_index_ridge(stack, lam, 1800.0, 900.0)
check("ridge (c): builder, width, etch depth, R_back equal the spec", d["variant"] == "ridge_waveguide" and g["ridge_width_nm"] == 1800.0
      and g["etch_depth_nm"] == 900.0 and g["R_back"] == -1.0)
check("ridge (c): stack = waveguide.hkust_ridge_stack(lambda): names, thicknesses, indices",
      [x["name"] for x in g["stack"]] == [L_.name for L_ in stack] and [x["t_nm"] for x in g["stack"]] == [L_.thickness_nm for L_ in stack]
      and all(close(x["n"], L_.n) for x, L_ in zip(g["stack"], stack)))
tot = sum(L_.thickness_nm for L_ in stack)
check("ridge (c): stack height, etch floor = height - etch, dot-layer centre", close(g["total_nm"], tot) and close(g["etch_floor_nm"], tot - 900.0)
      and close(g["dot_height_nm"], 1000 + 148 + 1 + 1.0))
check("ridge (c): mode radii wx, wy and n_eff = waveguide.effective_index_ridge", close(g["wx_um"], mode.wx_um, 1e-9) and close(g["wy_um"], mode.wy_um, 1e-9)
      and close(g["n_eff"], mode.n_eff, 1e-9), f"wx {g['wx_um']:.3f} wy {g['wy_um']:.3f} um")
e = waveguide.edge_emission(stack, 1800.0, 900.0, lam, g["L_um"], 0.75, alpha_cm=g["alpha_cm"], R_back=None, na_method="gaussian", dot_position=g["dot_position"])
check("ridge (c): F_P readout = waveguide.edge_emission F_wg; front-facet T = edge_emission T_facet", close(g["F_P"], e.F_wg, 1e-9) and close(g["T_facet"], e.T_facet, 1e-9))

# --- (c) generic_cavity: lumped schematic, and the illustrative micropillar of a paper card
s = spec("generic_cavity", {"F_P": 7.0, "kappa_ueV": 200.0})
d = describe(s)
g = d["geometry"]
check("lumped (c): builder lumped_cavity, geometry says lumped; F_P, kappa equal the spec", d["variant"] == "lumped_cavity" and g["lumped"] is True
      and g["F_P"] == 7.0 and close(g["kappa_ueV"], 200.0))
ment = _reg.get("material", "ingaas_gaas")
em = ment.build({p_.name: p_.default for p_ in ment.params}, 4.2)
check("lumped (c): Q_eff = E_mode/kappa", close(g["Q_eff"], em.E_X_meV / 0.2, 1e-6), f"{g['Q_eff']:.1f}")
check("lumped (c): the labels say there is no geometry", any("NO GEOMETRY" in x for x in d["labels"]))
d = describe(s, display={"display_as": "micropillar"})
g = d["geometry"]
check("lumped (c): display_as micropillar draws the micropillar builder with every geometry number flagged display_default + illustrative",
      d["variant"] == "micropillar" and g["illustrative"] is True and all(d["geometry_src"][k]["display_default"] for k in ("d_um", "N_top", "N_bot", "n_high", "spacer_nm"))
      and g["F_P"] == 7.0 and d["geometry_src"]["F_P"]["display_default"] is False)
check("lumped (c): the illustrative label is the first label", "illustrative geometry" in d["labels"][0])
check("lumped (c): display_as on a structure other than generic_cavity -> fallback, no raise", describe(spec("micropillar"), display={"display_as": "micropillar"})["fallback"] is True)
sp_s = dict(_api.central_design("somaschi2016", yaml.safe_load((ROOT / "cards/papers/somaschi2016.yaml").read_text(encoding="utf-8"))["device"])["design"]["composed"], schema="fsim-composed/1")
check("lumped (c): the somaschi/thomas cards carry meta.display_as and the API hint recovers it from the bare spec",
      _api.display_hint(sp_s) == {"display_as": "micropillar"}
      and all(yaml.safe_load((ROOT / f"cards/papers/{c}.yaml").read_text(encoding="utf-8"))["meta"]["display_as"] == "micropillar"
              for c in ("somaschi2016", "somaschi2016_qd1", "thomas2021", "thomas2021_res")))
sp_s["structure"]["params"]["F_P"] = {"value": 3.3, "unit": "1"}
check("lumped (c): editing a structure parameter drops the hint (lumped schematic again)", _api.display_hint(sp_s) is None
      and describe(sp_s, _api.display_hint(sp_s))["variant"] == "lumped_cavity")
check("new descriptors: JSON-safe and layers present for the 2D fallback",
      all(json.loads(json.dumps(describe(spec(k, mat=MAT.get(k, "ingaas_gaas"))))) == describe(spec(k, mat=MAT.get(k, "ingaas_gaas")))
          and describe(spec(k, mat=MAT.get(k, "ingaas_gaas")))["layers"] for k in ("photonic_nanowire", "planar_dbr", "nitride_planar_cavity", "edge_ridge", "generic_cavity")))

n = sum(RES)
print(f"{n}/{len(RES)} scene_compose checks passed")
sys.exit(0 if n == len(RES) else 1)
