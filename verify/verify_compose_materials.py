"""Checks for the material registry of the composable engine (fsim_core/compose/materials.py, P2).

Covers registration and the JSON payload, the provenance schema of data/materials.yaml, anchors the
code did not produce (qd_gf.PhononParams defaults, E_X = hc/lambda, Z(4 K) class range, Mueller E_b,
hbar conversions, Schweickert XX ratio) and finite EmitterModel output at 4, 40 and 300 K.

Run: python verify/verify_compose_materials.py   (exit code 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import cavity, linewidth, presets, qd_gf                  # noqa: E402
from fsim_core.compose import materials as M                             # noqa: E402
from fsim_core.compose import registry as R                              # noqa: E402
from fsim_core.compose import spec as S                                  # noqa: E402
from fsim_core.compose.types import EmitterModel                         # noqa: E402

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not ok else ""))


KINDS = ["inas_gaas_sk", "ingaas_gaas", "gaas_algaas_droplet", "inas_inp_cband",
         "inas_ingaas_metamorphic_cband", "inp_gainp_red", "ingan_gan", "gan_aln", "cdse_znse",
         "cspbbr3_colloidal"]
REQUIRED = ["lambda_nm", "tau_bare_ns", "fss_ueV", "E_b_meV", "gamma_xx_ratio", "gamma_star_4K_ns",
            "a_ac_ns_per_K", "alpha_ps2", "D_c_eV", "D_v_eV", "rho_kg_m3", "c_LA_m_s", "l_xy_nm", "l_z_nm",
            "n_host", "E_LO_meV", "k_nr_ns", "sigma_sd_ueV", "tau_sd_ns", "charge_state", "dipole"]
HBAR_MEV_NS = 6.582119569e-4          # independent CODATA value (meV ns)
HC_MEV_NM = 1239841.98                # hc in meV nm (spec value)

# 1. registration and payload ------------------------------------------------------------------
kinds = R.kinds("material")
check("all 10 material kinds registered on import", all(k in kinds for k in KINDS), str(kinds))
payload = R.registry_payload()
try:
    js = json.dumps(payload)
    ok_json = True
except (TypeError, ValueError) as e:
    ok_json, js = False, str(e)
check("registry_payload() is JSON-serialisable", ok_json, js[:80])
mats = {m["kind"]: m for m in payload["materials"]}
check("payload lists every material with all required params",
      all(k in mats and all(n in [p["name"] for p in mats[k]["params"]] for n in REQUIRED) for k in KINDS))
check("every param carries unit/tag/source/help and a default inside its band",
      all(p["tag"] in ("V", "DR", "E", "A") and p["source"] and p["help"] and p["unit"] is not None
          and (p["lo"] is None or p["hi"] is None or p["lo"] <= p["default"] <= p["hi"]
               or p["choices"] is not None)
          for m in mats.values() for p in m["params"]))
check("alpha_ps2, D_c, l_xy, l_z, sigma_sd, tau_sd are advanced; lambda_nm is main",
      all({p["name"]: p["group"] for p in m["params"]}[n] == "advanced"
          for m in mats.values() for n in ("alpha_ps2", "D_c_eV", "l_xy_nm", "l_z_nm", "sigma_sd_ueV", "tau_sd_ns"))
      and all({p["name"]: p["group"] for p in m["params"]}["lambda_nm"] == "main" for m in mats.values()))
check("class-only entries are labelled [E] class-only, not anchored in the summary",
      all("class-only, not anchored" in mats[k]["summary"] for k in
          ("cdse_znse", "cspbbr3_colloidal", "ingan_gan", "gan_aln", "inas_ingaas_metamorphic_cband")))

# 2. yaml provenance schema --------------------------------------------------------------------
data = yaml.safe_load(open(M.DATA_PATH, encoding="utf-8"))
leaf_kinds = ("value", "value_from", "range", "range_from")
bad = []
n_leaves = 0


def walk(prefix, d):
    global n_leaves
    n_leaves += 1
    present = [k for k in leaf_kinds if k in d]
    ok = (len(present) == 1 and isinstance(d.get("unit"), str) and d["unit"] != ""
          and d.get("tag") in ("V", "DR", "E", "A") and isinstance(d.get("source"), str) and d["source"].strip())
    if "value" in d:
        ok = ok and isinstance(d["value"], (int, float)) and math.isfinite(d["value"])
    if "range" in d:
        ok = ok and len(d["range"]) == 2 and d["range"][0] <= d["range"][1] and all(math.isfinite(x) for x in d["range"])
    if "range_from" in d:
        ok = ok and len(d["range_from"]) == 2
    if "convert" in d:
        ok = ok and d["convert"] in M._CONVERT
    if not ok:
        bad.append(prefix)


for name, leaf in data["_shared"].items():
    walk(f"_shared.{name}", leaf)
for kind, m in data["materials"].items():
    for name, leaf in m["params"].items():
        walk(f"{kind}.{name}", leaf)
check("every number in materials.yaml has unit, tag in V/DR/E/A and a non-empty source",
      not bad and n_leaves > 100, f"bad leaves: {bad[:6]} of {n_leaves}")
check("every material in the yaml defines exactly the known yaml params",
      all(set(m["params"]) == set(M._PARAMS) for m in data["materials"].values()))

# 3. anchors the code did not produce ------------------------------------------------------------
std = qd_gf.PhononParams()
for k in ("inp_gainp_red", "inas_inp_cband"):
    em = R.get("material", k).build({}, 4.0)
    ph = em.phonon
    same = (ph.D_e_eV == std.D_e_eV and ph.D_h_eV == std.D_h_eV and ph.rho_m_kg_m3 == std.rho_m_kg_m3
            and ph.c_s_m_s == std.c_s_m_s)
    check(f"{k}: InP phonon constants equal qd_gf.PhononParams defaults exactly", same)
em = R.get("material", "inp_gainp_red").build({}, 4.0)
check("inp_gainp_red: phonon card (alpha, cutoff geometry) == PhononParams() default, bit for bit",
      em.phonon == std and qd_gf.coupling_alpha_ps2(em.phonon) == qd_gf.coupling_alpha_ps2(std)
      and float(qd_gf.spectral_density(1.3, em.phonon)) == float(qd_gf.spectral_density(1.3, std)))
inp = R.get("material", "inp_gainp_red")
check("inp_gainp_red: gamma*(4 K), a_ac, E_LO, E_b band come from linewidth/presets (not copied)",
      abs(inp.param("gamma_star_4K_ns").default - 0.25 / HBAR_MEV_NS) < 1e-6
      and abs(inp.param("a_ac_ns_per_K").default - 2.0e-3 / HBAR_MEV_NS) < 1e-6
      and inp.param("E_LO_meV").default == 43.0
      and (inp.param("E_b_meV").lo, inp.param("E_b_meV").hi) == (3.5, 8.0)
      and abs(inp.param("gamma_star_4K_ns").default - linewidth.LinewidthParams().gamma0 / HBAR_MEV_NS) < 1e-9
      and presets.DOT_PRESETS["chatzarakis-class"]["delta_xx"] == 5.9
      and inp.param("E_b_meV").lo <= 5.9 <= inp.param("E_b_meV").hi)

for k in ("inas_gaas_sk", "ingaas_gaas", "gaas_algaas_droplet"):
    e = R.get("material", k)
    em = e.build({}, 4.0)
    gaas_consts = (em.phonon.D_e_eV, em.phonon.D_h_eV, em.phonon.rho_m_kg_m3, em.phonon.c_s_m_s) == (-14.6, -4.8, 5370.0, 5110.0)
    Z = qd_gf.zpl_weight(em.phonon, 4.0)
    check(f"{k}: GaAs constants (-14.6, -4.8, 5370, 5110) and Z(4 K) = {Z:.3f} in 0.85-0.97 "
          "(Thomas 2021 / Uppu ZPL 0.95 class [V])", gaas_consts and 0.85 <= Z <= 0.97, f"Z={Z}")
# explicit constants through the pure qd_gf API (independent of the registry build path)
Zc = qd_gf.zpl_weight(qd_gf.PhononParams(D_e_eV=-14.6, D_h_eV=-4.8, rho_m_kg_m3=5370.0, c_s_m_s=5110.0,
                                         l_xy_nm=4.0, l_z_nm=1.5), 4.0)
check("Z(4 K) from qd_gf.zpl_weight with explicit GaAs constants, l = (4, 1.5) nm, in 0.85-0.97",
      0.85 <= Zc <= 0.97, f"Z={Zc}")

worst = 0.0
for k in KINDS:
    e = R.get("material", k)
    for lam in (e.param("lambda_nm").default, 310.0, 1310.0, 1550.0):
        em = e.build({"lambda_nm": lam}, 4.0)
        worst = max(worst, abs(em.E_X_meV - HC_MEV_NM / lam))
check("E_X(lambda) = 1239841.98 / lambda_nm meV to 1e-9 for every kind", worst < 1e-9, f"{worst}")
check("inas_inp_cband O-band (1310 nm) E_X = 946.4 meV",
      abs(R.get("material", "inas_inp_cband").build({"lambda_nm": 1310.0}, 4.0).E_X_meV - 946.45) < 0.01)

mu = R.get("material", "inas_gaas_sk")
check("Mueller E_b = 2.3 meV [V] lies in the inas_gaas_sk range and is reproduced when overridden",
      mu.param("E_b_meV").lo <= 2.3 <= mu.param("E_b_meV").hi
      and mu.build({"E_b_meV": 2.3}, 4.0).delta_xx_meV == -2.3)  # EmitterModel: negative = bound
check("Somaschi tau_bulk = 1.3 ns [V] lies in the inas_gaas_sk tau range and gives Gamma0 = 1/1.3",
      mu.param("tau_bare_ns").lo <= 1.3 <= mu.param("tau_bare_ns").hi
      and abs(mu.build({"tau_bare_ns": 1.3}, 4.0).gamma0_ns - 1.0 / 1.3) < 1e-12)
dr = R.get("material", "gaas_algaas_droplet")
check("droplet: default X tau / Schweickert XX 0.125 ns [V] = gamma_xx_ratio; E_b range contains Piccinini 4.0",
      abs(dr.param("tau_bare_ns").default / 0.125 - dr.param("gamma_xx_ratio").default) < 1e-9
      and dr.param("E_b_meV").lo <= 4.0 <= dr.param("E_b_meV").hi)

for k in KINDS:
    e = R.get("material", k)
    g4 = e.param("gamma_star_4K_ns").default
    a = e.param("a_ac_ns_per_K").default
    ok4 = e.build({}, 4.0).gamma_star_ns == g4
    ok40 = abs(e.build({}, 40.0).gamma_star_ns - (g4 + a * 36.0)) < 1e-9
    ok300 = abs(e.build({}, 300.0).gamma_star_ns - (g4 + a * 296.0)) < 1e-9
    check(f"{k}: gamma*(4 K) equals the param; linear T model at 40 and 300 K", ok4 and ok40 and ok300)
check("gamma* is clipped at 0",
      mu.build({"gamma_star_4K_ns": 0.0, "a_ac_ns_per_K": 1.0}, 1.0).gamma_star_ns == 0.0)
check("1 ueV excess linewidth = 1.5193 / ns (gamma* = W / hbar)",
      abs(mu.param("gamma_star_4K_ns").hi - 1.0e-3 / HBAR_MEV_NS) < 1e-9
      and abs(mu.param("gamma_star_4K_ns").hi - 1.5193) < 1e-4)
check("alpha_ps2 override pins alpha (0.027 ps^2 Ramsay-class)",
      qd_gf.coupling_alpha_ps2(mu.build({"alpha_ps2": 0.027}, 4.0).phonon) == 0.027
      and qd_gf.coupling_alpha_ps2(mu.build({}, 4.0).phonon) > 0.0)
h = 1.0e-3
numeric = (cavity.varshni_shift(300.0 + h) - cavity.varshni_shift(300.0 - h)) / (2 * h) * 1e3
check("dE/dT of GaAs-class kinds at 300 K matches a numeric derivative of cavity.varshni_shift",
      abs(mu.build({}, 300.0).dEdT_meV_K - numeric) < 1e-6, f"{mu.build({}, 300.0).dEdT_meV_K} vs {numeric}")
check("charge_state X+ adds the charged transition; default has {X, XX}",
      mu.build({}, 4.0).transitions == frozenset({"X", "XX"})
      and mu.build({"charge_state": "X+"}, 4.0).transitions == frozenset({"X", "XX", "X+"}))

check("every material provides X, XX, X+ and X- (tpe drive requires XX)",
      all(R.get("material", k).provides == frozenset({"X", "XX", "X+", "X-"}) for k in KINDS))

# validate_spec accepts a material block with a tagged override leaf ------------------------------
vs = S.validate_spec({"schema": S.SCHEMA, "T_K": 4.2,
                      "material": {"kind": "inas_gaas_sk", "params": {
                          "E_b_meV": {"value": 2.3, "unit": "meV", "tag": "V", "source": "Mueller, Nat. Photon. 8, 224 (2014)"},
                          "tau_bare_ns": {"range": [0.8, 1.3], "unit": "ns", "tag": "E", "source": "class"}}},
                      "structure": {"kind": "x"}, "drive": {"kind": "x"}, "setup": {"kind": "lens_fibre_det"}})
mat_issues = [i for i in vs["issues"] if i["path"].startswith("material")]
check("validate_spec: a material block with tagged leaves has no material issues",
      not mat_issues and "material.E_b_meV" in vs["resolved"], str(mat_issues))

# 4. build() finite at 4, 40, 300 K ---------------------------------------------------------------
fin_ok = True
msg = ""
for k in KINDS:
    for T in (4.0, 40.0, 300.0):
        em = R.get("material", k).build({}, T)
        nums = [em.E_X_meV, em.lambda_nm, em.n_host, em.gamma0_ns, em.k_nr_ns, em.gamma_star_ns, em.delta_xx_meV,
                em.gamma_xx_ratio, em.fss_ueV, em.dEdT_meV_K, em.sigma_sd_ueV, em.tau_sd_ns]
        z = qd_gf.zpl_weight(em.phonon, T)
        if not (isinstance(em, EmitterModel) and all(math.isfinite(x) for x in nums) and em.gamma0_ns > 0
                and em.gamma_star_ns >= 0 and 0.0 < z <= 1.0 and isinstance(em.transitions, frozenset)):
            fin_ok, msg = False, f"{k} @ {T} K"
check("build() at 4, 40, 300 K returns a finite EmitterModel for every kind", fin_ok, msg)
check("EmitterModel.inputs maps fields to material.<param> paths",
      all(v.startswith("material.") for v in mu.build({}, 4.0).inputs.values()) and "E_X_meV" in mu.build({}, 4.0).inputs)

n, total = sum(RESULTS), len(RESULTS)
print(f"{n}/{total} compose materials checks passed")
sys.exit(0 if n == total else 1)
