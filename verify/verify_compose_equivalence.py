"""Checks for the tier passthroughs of the composable engine (fsim_core/compose/tiers.py, P6).

  1. every shipped design card through tiers.run == device.evaluate(DeviceDesign.load(card))
     exactly (scalars NaN-aware dict equality, curves np.array_equal); nanowire cards at T_grid=[300.0]
  2. overrides through tiers.run == the same META-path edit done by hand; the input design is not mutated
  3. tier table / classification / handoff descriptors / payload
  4. Fock card -> composed spec mapping round trip (against resonant_source's own FockSource)
  5. engine equivalence: engine.evaluate_spec(fock spec) vs resonant_source.bs_figure_of_merit,
     rtol 1e-9 (SKIPPED while the P5 engine is not merged)

Run: python verify/verify_compose_equivalence.py   (exit code 0 iff no check fails)
"""
from __future__ import annotations

import copy
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import resonant_source as rs                    # noqa: E402
from fsim_core.compose import tiers as T                       # noqa: E402
from fsim_core.compose.spec import ComposedSpec, validate_spec  # noqa: E402
from fsim_core.device import DeviceDesign, _set_path, evaluate  # noqa: E402

RESULTS = []
SKIPPED = []
T0 = time.time()


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not ok else ""))


def skip(name, reason):
    SKIPPED.append(name)
    print(f"[SKIPPED] {name}  ({reason})")


def same(a, b) -> bool:
    """Exact, NaN-aware structural equality (dicts, sequences, arrays, scalars)."""
    if isinstance(a, dict) or isinstance(b, dict):
        return (isinstance(a, dict) and isinstance(b, dict) and a.keys() == b.keys()
                and all(same(a[k], b[k]) for k in a))
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        return (isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)) and len(a) == len(b)
                and all(same(x, y) for x, y in zip(a, b)))
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    try:
        if isinstance(a, np.floating) or isinstance(b, np.floating):
            a, b = float(a), float(b)
            if math.isnan(a) and math.isnan(b):
                return True
    except Exception:
        pass
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    return bool(a == b)


def grid_for(tier_id):
    return [300.0] if tier_id == "nitride_nanowire" else None


# ---------------------------------------------------------------- 1. every design card, bit-identical
scan = T._scan_cards()
design_cards = sorted(c for c, (t, _v) in scan.items() if t != "fock")
all_design_files = sorted(f.stem for f in (ROOT / "cards").glob("*.yaml")
                          if isinstance(yaml.safe_load(f.read_text(encoding="utf-8")), dict)
                          and "design" in yaml.safe_load(f.read_text(encoding="utf-8")))
check("tier scan finds all 14 shipped design cards", design_cards == all_design_files and len(design_cards) == 14,
      f"{len(design_cards)} vs {len(all_design_files)}")

per_tier_time = {}
for card in design_cards:
    tier_id = scan[card][0]
    grid = grid_for(tier_id)
    t = time.time()
    got = T.TIERS[tier_id]["run"](card, None, grid)
    ref = evaluate(DeviceDesign.load(ROOT / "cards" / f"{card}.yaml"), grid)
    per_tier_time[tier_id] = per_tier_time.get(tier_id, 0.0) + time.time() - t
    check(f"passthrough {tier_id}: {card} == evaluate()",
          same(got["scalars"], ref["scalars"]) and same(got["curves"], ref["curves"])
          and set(got) == set(ref),
          f"scalars equal {same(got['scalars'], ref['scalars'])}, curves equal {same(got['curves'], ref['curves'])}")
print("   tier wall time (s):", {k: round(v, 1) for k, v in per_tier_time.items()})

# ---------------------------------------------------------------- 2. overrides == hand edit
STAGED = ROOT / "cards" / "staged-device-design.yaml"
base = DeviceDesign.load(STAGED)
ovr = {"dot.gamma_scale": 1.2, "thermal.T_hs": 280.0}
by_hand = DeviceDesign.load(STAGED)
for p, v in ovr.items():
    _set_path(by_hand, p, v)
ref_o = evaluate(by_hand)
got_o = T.TIERS["legacy"]["run"](STAGED, ovr)
check("legacy overrides (dot.gamma_scale, thermal.T_hs) == hand-edited DeviceDesign",
      same(got_o["scalars"], ref_o["scalars"]) and same(got_o["curves"], ref_o["curves"]))
ref_base = evaluate(DeviceDesign.load(STAGED))
check("the override actually changes the result (test is not vacuous)",
      not same(got_o["curves"], ref_base["curves"]) or not same(got_o["scalars"], ref_base["scalars"]))
snap = copy.deepcopy(base)
T.run_design(base, ovr)
check("run() does not mutate a DeviceDesign passed in", base == snap)
got_d = T.TIERS["legacy"]["run"](base, {"dot.gamma_scale": 1.2})
by_hand2 = DeviceDesign.load(STAGED)
_set_path(by_hand2, "dot.gamma_scale", 1.2)
check("a DeviceDesign object is accepted in place of a card path",
      same(got_d["scalars"], evaluate(by_hand2)["scalars"]))
# one override on a nitride tier (one T point for the nanowire)
PL = ROOT / "cards" / "nitride-cavity-pulse-design.yaml"
hand_pl = DeviceDesign.load(PL)
_set_path(hand_pl, "thermal.T_hs", 290.0)
got_pl = T.TIERS["nitride_planar"]["run"](PL, {"thermal.T_hs": 290.0})
ref_pl = evaluate(hand_pl)
check("nitride_planar override thermal.T_hs == hand edit",
      same(got_pl["scalars"], ref_pl["scalars"]) and same(got_pl["curves"], ref_pl["curves"]))
NW = ROOT / "cards" / "nitride-nanowire-vertical-pulse-design.yaml"
hand_nw = DeviceDesign.load(NW)
_set_path(hand_nw, "thermal.T_hs", 300.0)
got_nw = T.TIERS["nitride_nanowire"]["run"](NW, {"thermal.T_hs": 300.0}, [300.0])
ref_nw = evaluate(hand_nw, [300.0])
check("nitride_nanowire override thermal.T_hs == hand edit (one T point)",
      same(got_nw["scalars"], ref_nw["scalars"]) and same(got_nw["curves"], ref_nw["curves"]))


def raises(exc, fn, *a, contains=(), **kw):
    try:
        fn(*a, **kw)
    except exc as e:
        return all(c in str(e) for c in contains)
    except Exception:
        return False
    return False


check("unknown override path rejected with a clear error",
      raises(KeyError, T.run_design, STAGED, {"dot.no_such_field": 1.0}, contains=("dot.no_such_field",)))
check("malformed override path rejected", raises(KeyError, T.run_design, STAGED, {"gamma_scale": 1.0}))
check("a card of another tier is refused by a tier runner",
      raises(ValueError, T.TIERS["legacy"]["run"], PL, None, contains=("nitride_planar",)))
check("unknown card name raises FileNotFoundError", raises(FileNotFoundError, T.run_design, "no-such-card"))

# ---------------------------------------------------------------- 3. tier table, descriptors, payload
check("TIERS has the five tiers", list(T.TIERS) == ["legacy", "rt_edge", "nitride_planar", "nitride_nanowire", "fock"])
check("tier card counts legacy 1 / edge 2 / planar 7 / nanowire 4 / fock 1",
      [len(T.TIERS[k]["cards"]) for k in T.TIERS] == [1, 2, 7, 4, 1])
check("nanowire variants split 2 horizontal + 2 vertical",
      {k: len(v) for k, v in T.TIERS["nitride_nanowire"]["variants"].items()} == {"horizontal": 2, "vertical": 2})
check("every tier has a callable run", all(callable(t["run"]) for t in T.TIERS.values()))
check("spec_for_tier_card: edge, planar, nanowire H, legacy, fock",
      T.spec_for_tier_card("edge-inp-gainp-design") == {"tier": "rt_edge", "card": "edge-inp-gainp-design"}
      and T.spec_for_tier_card("nitride-nonpolar-set-design")["tier"] == "nitride_planar"
      and T.spec_for_tier_card("nitride-nanowire-horizontal-set-design")
      == {"tier": "nitride_nanowire", "card": "nitride-nanowire-horizontal-set-design", "variant": "horizontal"}
      and T.spec_for_tier_card(STAGED) == {"tier": "legacy", "card": "staged-device-design"}
      and T.spec_for_tier_card("fock-source-4K-design.yaml") == {"tier": "fock", "card": "fock-source-4K-design"})
check("spec_for_tier_card returns None for validation/prediction cards and unknown names",
      all(T.spec_for_tier_card(c) is None for c in ("zhao", "kitamura", "qcap-staged", "reischle", "nope")))
check("spec_for_tier_card accepts a loaded DeviceDesign",
      T.spec_for_tier_card(DeviceDesign.load(PL)) == {"tier": "nitride_planar", "card": "nitride-cavity-pulse-design"})
pay = T.tiers_payload()
check("tiers_payload is JSON-serialisable with label/summary/cards/meta_paths",
      len(json.dumps(pay)) > 100 and [p["id"] for p in pay] == list(T.TIERS)
      and all(p["label"] and p["summary"] and isinstance(p["cards"], list) and p["meta_paths"] > 0 for p in pay))
from fsim_core.design_meta import META                                   # noqa: E402
check("meta_paths == len(design_meta.META) for DeviceDesign tiers",
      all(p["meta_paths"] == len(META) for p in pay if p["id"] != "fock"))

# ---------------------------------------------------------------- 4. Fock card -> spec mapping round trip
FOCK = ROOT / "cards" / "fock-source-4K-design.yaml"
card = rs.load_fock_card(FOCK)
src = rs.source_from_card(card)
spec = T.fock_card_to_spec(FOCK)
res = validate_spec(spec)
check("fock spec validates against the registry with no issues", res["ok"] and not res["issues"], str(res["issues"]))
check("fock spec: ComposedSpec round trip (from_dict -> to_dict) is identical",
      ComposedSpec.from_dict(spec).to_dict() == spec)
check("fock spec: YAML dump/load round trip is identical", yaml.safe_load(yaml.safe_dump(spec)) == spec)
check("fock spec: kinds are inas_gaas_sk / generic_cavity / resonant_pulse / lens_fibre_det",
      [spec[a]["kind"] for a in ("material", "structure", "drive", "setup")]
      == ["inas_gaas_sk", "generic_cavity", "resonant_pulse", "lens_fibre_det"])
leaves = [("T_K", spec["T_K"])] + [(f"{a}.{k}", v) for a in ("material", "structure", "drive", "setup")
                                    for k, v in spec[a]["params"].items()]
check("every fock spec leaf carries tag in V/DR/E/A and a non-empty source",
      all(isinstance(v, dict) and v.get("tag") in ("V", "DR", "E", "A") and v.get("source") for _p, v in leaves))
back = T.fock_source_from_spec(spec)
direct = {"T_K": src.T_K, "transition": src.transition, "excitation": src.excitation, "F_P": src.F_P,
          "kappa_ueV": src.kappa_ueV, "cav_detune_ueV": src.cav_detune_ueV, "beta0": src.beta0,
          "eta_out": src.eta_out, "pulse_shape": src.pulse_shape, "pulse_fwhm_ps": src.pulse_fwhm_ps,
          "f_rep_MHz": src.f_rep_MHz, "gamma_xx_ratio": src.gamma_xx_ratio, "sigma_sd_ueV": src.sigma_sd_ueV,
          "tau_sd_ns": src.tau_sd_ns, "phonon_alpha_ps2": src.phonon_alpha_ps2, "eta_det": src.eta_det,
          "zpl_output": src.zpl_output}
bad = [k for k, v in direct.items() if back[k] != v]
check("fock round trip: directly mapped entries equal the FockSource at the card central point", not bad, str(bad))
check("fock round trip: gamma0 = 1/tau_bare and delta_xx = -E_b to rtol 1e-12",
      math.isclose(back["gamma0_ns"], src.gamma0_ns, rel_tol=1e-12)
      and math.isclose(back["delta_xx_meV"], src.delta_xx_meV, rel_tol=1e-12))
check("fock round trip: eta_optics == eta_demux x eta_circuit",
      math.isclose(back["eta_optics"], src.eta_demux * src.eta_circuit, rel_tol=1e-15))
rates = rs.source_rates(src)
check("fock round trip: gamma* leaf == source_rates gamma_star and I = Gamma/(Gamma+gamma*) == 1/m (enhanced mode)",
      back["gamma_star_ns"] == rates["gamma_star"]
      and math.isclose(rates["Gamma"] / (rates["Gamma"] + rates["gamma_star"]), 1.0 / src.hom_width_multiple,
                       rel_tol=1e-12))
src_dot = rs.source_from_card(card, dephasing_mode="dot")
spec_dot = T.fock_card_to_spec(FOCK, {"dephasing_mode": "dot"})
r_dot = rs.source_rates(src_dot)
check("fock mapping honours dephasing_mode=dot: gamma* = (m-1) Gamma0",
      T.fock_source_from_spec(spec_dot)["gamma_star_ns"] == r_dot["gamma_star"]
      and math.isclose(r_dot["gamma_star"], (src_dot.hom_width_multiple - 1.0) * src_dot.gamma0_ns, rel_tol=1e-15))
spec_tpe = T.fock_card_to_spec(FOCK, {"excitation": "tpe", "transition": "X"})
res_t = validate_spec(spec_tpe)
bt = T.fock_source_from_spec(spec_tpe)
check("fock tpe mapping: drive tpe, charge_state X, validates, E_b = -delta_xx",
      spec_tpe["drive"]["kind"] == "tpe" and spec_tpe["material"]["params"]["charge_state"]["value"] == "X"
      and res_t["ok"] and bt["excitation"] == "tpe"
      and math.isclose(bt["delta_xx_meV"], card["delta_xx_meV"].central(), rel_tol=1e-12), str(res_t["issues"]))
spec_ov = T.fock_card_to_spec(FOCK, {"F_P": 3.0, "pulse_fwhm_ps": 5.0})
check("fock mapping overrides are applied and recorded in the leaf source",
      spec_ov["structure"]["params"]["F_P"]["value"] == 3.0 and "overridden" in spec_ov["structure"]["params"]["F_P"]["source"]
      and spec_ov["drive"]["params"]["fwhm_ps"]["value"] == 5.0)
chk = T.TIERS["fock"]["run"](FOCK, None, None, 2)
ref_f = rs.bs_figure_of_merit(src, 2)
check("fock tier run == bs_figure_of_merit at the card central point",
      same(chk, ref_f))

# ---------------------------------------------------------------- 5. engine equivalence (needs P5)
try:
    from fsim_core.compose import engine as _engine
    out = _engine.evaluate_spec(spec)
except NotImplementedError as e:
    out = None
    skip("fock spec through engine.evaluate_spec == bs_figure_of_merit (g2, efficiency, I; rtol 1e-9)",
         f"P5 engine not merged yet: {e}")
except Exception as e:                                                   # noqa: BLE001
    out = None
    check("fock spec through engine.evaluate_spec runs", False, f"{type(e).__name__}: {e}")
if out is not None:
    sc = out["scalars"]
    fom = rs.bs_figure_of_merit(src, 2)

    def pick(*names):
        for n in names:
            if n in sc:
                return float(sc[n])
        raise KeyError(f"none of {names} in engine scalars {sorted(sc)}")

    rtol = 1e-9
    try:
        g2_e = pick("g2_0")
        eta_e = (pick("rate_detected_Hz") / (src.f_rep_MHz * 1e6)) if "rate_detected_Hz" in sc else pick("eta")
        I_e = pick("hom_I")
        # I reference: the intrinsic single-pulse overlap the card's own pulse physics gives,
        # I_p = pulse_factor x Gamma/(Gamma + gamma*) (resonant_source._pulse_physics), times the unfiltered
        # phonon factor Z^2 (hom_phonon_factor). No spectral diffusion (that is the pairwise hom_matrix).
        rr = rs.source_rates(src)
        I_ref = fom["pulse_factor"] * rr["Gamma"] / (rr["Gamma"] + rr["gamma_star"]) * fom["hom_phonon_factor"]
        check("engine g2 == bs_figure_of_merit g2 (rtol 1e-9)", math.isclose(g2_e, fom["g2"][0], rel_tol=rtol),
              f"{g2_e} vs {fom['g2'][0]}")
        check("engine efficiency == bs_figure_of_merit eta (rtol 1e-9)", math.isclose(eta_e, fom["eta"], rel_tol=rtol),
              f"{eta_e} vs {fom['eta']}")
        check("engine hom_I == pulse_factor x Gamma/(Gamma+gamma*) x Z^2 of bs_figure_of_merit (rtol 1e-9)",
              math.isclose(I_e, I_ref, rel_tol=rtol), f"{I_e} vs {I_ref}")
    except KeyError as e:
        check("engine scalars contain the Fock figures of merit", False, str(e))

if out is not None:
    # variants: ZPL-filtered output (Z enters eta) and two-photon excitation (cascade g2 and I)
    for label, ov in (("zpl_output", {"zpl_output": True}),
                      ("tpe", {"excitation": "tpe", "transition": "X"})):
        try:
            s_v = rs.source_from_card(card, **ov)
            fom_v = rs.bs_figure_of_merit(s_v, 2)
            sc_v = _engine.evaluate_spec(T.fock_card_to_spec(FOCK, ov))["scalars"]
            rr_v = rs.source_rates(s_v)
            I_v = (fom_v["pulse_factor"] * rr_v["Gamma"] / (rr_v["Gamma"] + rr_v["gamma_star"])
                   * fom_v["hom_phonon_factor"])
            eta_v = sc_v["rate_detected_Hz"] / (s_v.f_rep_MHz * 1e6)
            ok = (math.isclose(sc_v["g2_0"], fom_v["g2"][0], rel_tol=1e-9)
                  and math.isclose(eta_v, fom_v["eta"], rel_tol=1e-9)
                  and math.isclose(sc_v["hom_I"], I_v, rel_tol=1e-9))
            check(f"engine variant {label}: g2, eta, I == bs_figure_of_merit (rtol 1e-9)", ok,
                  f"g2 {sc_v['g2_0']} vs {fom_v['g2'][0]}; eta {eta_v} vs {fom_v['eta']}; I {sc_v['hom_I']} vs {I_v}")
        except Exception as e:                                           # noqa: BLE001
            check(f"engine variant {label} runs", False, f"{type(e).__name__}: {e}")

# ---------------------------------------------------------------- summary
dt = time.time() - T0
n_ok = sum(RESULTS)
print(f"runtime {dt:.1f} s" + (f"; SKIPPED: {len(SKIPPED)}" if SKIPPED else ""))
print(f"{n_ok}/{len(RESULTS)} compose_equivalence checks passed" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(0 if n_ok == len(RESULTS) else 1)
