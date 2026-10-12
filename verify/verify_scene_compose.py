"""Checks for fsim_studio/scene_compose.py (P13): descriptor numbers equal the composed spec's resolved values.

Reference values come from fsim_core.compose.spec.validate_spec (the validator's own `resolved` table) and
closed forms (quarter wave lambda/4n, Gaussian-beam waist, equal-area circle, escape cone), not from
scene_compose.  Run: python verify/verify_scene_compose.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.compose.spec import validate_spec          # noqa: E402
from fsim_studio.scene_compose import describe             # noqa: E402

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
for k, why in (("edge_ridge", "no 3D builder"), ("generic_cavity", "no 3D builder"), ("nope", "unknown structure")):
    d = describe(spec(k))
    check(f"{k}: fallback descriptor", d["fallback"] is True and d["variant"] == "fallback2d" and why in d["reason"], d.get("reason", ""))
check("non-mapping spec -> fallback, no raise", describe("x")["fallback"] is True)
check("bad leaf (unit mismatch) -> fallback, no raise",
      describe(spec("micropillar", {"d_um": {"value": 2.5, "unit": "nm"}}))["fallback"] is True)
check("descriptor layers present for the 2D fallback view",
      all(describe(spec(k))["layers"] for k in ("micropillar", "open_cavity", "cbg", "phc_waveguide", "phc_cavity", "sil_mesa")))

n = sum(RES)
print(f"{n}/{len(RES)} scene_compose checks passed")
sys.exit(0 if n == len(RES) else 1)
