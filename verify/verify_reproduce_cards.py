"""Gating reproduction of published measurements by the PRODUCTION model, against an
EXPECTATION LEDGER (verify/data/reproduce_ledger.yaml).  devlib-vf2.

Closes p0-verifiability defects H3/H4/H5 and M4/M5/M6: device-level model outputs for the
existing literature cards are compared with the published numbers under tolerances that
were written down in the ledger BEFORE the first run (class C calibration, U unfitted,
H held-out, KD known deviation).  Each item carries `expect: pass | fail`; a documented
model miss is a gated `expect: fail`, so an unexpected pass (the model changed) and an
unexpected fail (a regression) both turn the exit code to 1.

Evaluators call production code only (fsim_core.fitting.fit_phase0, fsim_core.device,
fsim_core.sde, fsim_core.dot_levels, fsim_core.nitride_levels, scripts.run_nitride_cavity,
fsim_core.nitride_nanowire_*); they do not re-implement the model.  scripts/run_v1_validation.py
runs its computation at import time (writes out/validation), so its CALL SEQUENCE is
replicated here instead of imported.

Usage:
    python verify/verify_reproduce_cards.py            # gating run (non-slow items)
    python verify/verify_reproduce_cards.py --all      # also the items marked `slow: true`
    python verify/verify_reproduce_cards.py --self-test

Exit 0 iff every item's outcome equals its ledger `expect`.
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

LEDGER_PATH = ROOT / "verify" / "data" / "reproduce_ledger.yaml"
HC_EV_NM = 1239.841984  # eV nm

# ----------------------------------------------------------------------------- evaluators
# Each returns {key: float}; a "<key>__se" entry carries a Monte-Carlo standard error.
# Groups listed in SLOW_GROUPS run only with --all.
SLOW_GROUPS: set[str] = set()


def _chatz_card():
    from fsim_core.card import load_card
    return load_card(ROOT / "cards" / "chatzarakis.yaml")


def ev_chatz_va():
    """V-a joint fit in memory (production fit_phase0), per-point model g2."""
    from fsim_core.fitting import fit_phase0
    fit = fit_phase0(_chatz_card())
    return {f"g2_T{int(r['T'])}": float(m) for r, m in zip(fit.data, fit.model_g2)}


def _ensure_class_proxy():
    """device.class_proxy_params() reads the gitignored out/phase0/fit_params.json.  On a
    clean clone it is absent: rebuild the same dict from the in-memory V-a fit and patch
    the module function (same values to ~1e-8, p0 audit), so the item still runs."""
    import fsim_core.device as dev
    if (ROOT / "out" / "phase0" / "fit_params.json").is_file():
        return "out/phase0/fit_params.json"
    from fsim_core.fitting import fit_phase0
    fit = fit_phase0(_chatz_card())
    keep = {k: float(fit.params[k]) for k in ("gamma0", "a_ac", "b_lo", "E_lo", "r_xx",
                                              "a_esc", "E_a", "b_p", "E_b", "b0", "beta")}
    dev.class_proxy_params = lambda: dict(keep)
    return "in-memory fit_phase0"


def _base_design(DeviceDesign, lineshape, geom, alpha=None, rho_one=True):
    # call sequence of scripts/run_v1_validation.py::base_design
    d = DeviceDesign()
    d.drive.V = 0.0
    d.drive.I_uA = 0.0
    d.drive.duty = 0.0
    d.drive.b_e = 0.0
    d.drive.mode = "PL"
    d.cavity.enabled = False
    d.filter.enabled = True
    d.filter.auto_w = False
    d.dot.lineshape = lineshape
    ph = dict(geom)
    if alpha is not None:
        ph["alpha_ps2"] = alpha
    d.dot.phonon = ph
    if rho_one:
        d.ret.b0 = 1e-12
        d.ret.beta = 1e-12
    return d


def _run_point(DeviceDesign, evaluate, d, T, w, dx, mu, delta, gamma_target=None):
    # call sequence of scripts/run_v1_validation.py::run_point
    d.thermal.T_hs = T
    d.filter.w = w
    d.filter.dx = dx
    d.drive.mu = mu
    d.dot.delta_xx = delta
    if gamma_target is not None:
        from fsim_core.device import class_proxy_params
        from fsim_core.spectral import gamma_of_T
        p = class_proxy_params()
        g_proxy = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        d.dot.gamma_scale = gamma_target / g_proxy
    return float(evaluate(d, T_grid=[T])["scalars"]["g2_op"])


def ev_chatz_dev_ibm():
    from fsim_core.device import DeviceDesign, evaluate
    _ensure_class_proxy()
    card = yaml.safe_load((ROOT / "cards" / "chatzarakis2023.yaml").read_text())
    dev = card["device"]
    out = {}
    for pt in card["data"]["g2_vs_T"]["points"]:
        d = _base_design(DeviceDesign, "ibm", dev["geometry"], alpha=dev["phonon_alpha_ps2"], rho_one=False)
        out[f"g2_T{int(pt['T'])}"] = _run_point(DeviceDesign, evaluate, d, float(pt["T"]), float(pt["w"]),
                                                float(pt["dx"]), dev["mu"], dev["delta_xx"])
    return out


def _centre(rng):
    return 0.5 * (float(rng[0]) + float(rng[1]))


def ev_laf_mid():
    """Pre-registered mid-point: centre of every [E] range (delta_xx, Gamma, w, mu)."""
    from fsim_core.device import DeviceDesign, evaluate
    _ensure_class_proxy()
    card = yaml.safe_load((ROOT / "cards" / "laferriere2023.yaml").read_text())
    dev = card["device"]
    out = {}
    for pt in card["data"]["g2_vs_T"]["points"]:
        T = float(pt["T"])
        key = str(int(T))
        mu_rng = dev["mu"]["range_77"] if T < 100 else dev["mu"]["range_high"]
        d = _base_design(DeviceDesign, "ibm", dev["geometry"])
        out[f"g2_T{int(T)}"] = _run_point(
            DeviceDesign, evaluate, d, T, _centre(dev["w_meV"][key]["range"]), 0.0, _centre(mu_rng),
            _centre(dev["delta_xx"]["range"]), _centre(dev["gamma_meV"][key]["range"]))
    return out


def ev_reischle():
    """Call-for-call copy of verify_rt_edge_papers._xcheck_reischle2008_g2_80K's fixture."""
    from fsim_core.device import DeviceDesign, evaluate
    d = DeviceDesign()
    d.ret.mode = "confinement"
    d.ret.preset = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
    d.drive.mode = "EL-transport"
    d.drive.diode = {"preset": "red", "tau_pulse_ns": 0.1}
    d.drive.duty = 0.008
    d.drive.n_dot_cm2 = 1e10
    d.thermal.T_hs = 80.0
    r = evaluate(d, [80.0])
    return {"g2_op": float(r["scalars"]["g2_op"])}


def ev_zhao():
    """scripts/run_zhao_fit.py main computation through production sde.py."""
    from fsim_core.sde import QDLaserParams, find_threshold, g2_vs_pump
    params = QDLaserParams()
    I_op = 4.0 * find_threshold(params)
    out = {}
    for label, F_pump in (("normal", 1.0), ("quiet", 0.08)):
        r = g2_vs_pump(params, I_op, n_runs=60, t_end=2.0e-9, dt=1.0e-13, seed=42, F_pump=F_pump)
        out[label] = float(r["g2_mean"])
        out[label + "__se"] = float(r["g2_se"])
    return out


def ev_dot_levels():
    from fsim_core import dot_levels as D
    P = D.class_presets()

    def lv(name, **kw):
        return D.levels(P[name](**kw))

    L = lv("InGaAs0.5/GaAs/GaAs on GaAs")
    Lc = lv("InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs")
    S = P["InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"]()
    S.geometry.matrix_thickness_nm = 2.0
    shift = D.levels(S).dE_pair_WL - L.dE_pair_WL
    Lr = lv("InP/GaInP/AlGaInP0.55 on GaAs")
    Lt = lv("InAs/InP/InP")
    Lg = lv("InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs", height_nm=5.5)
    return {"chatz_EX_eV": float(L.E_X_eV), "chatz_dEWL_meV": float(L.dE_pair_WL),
            "ssl0_dEWL_meV": float(Lc.dE_pair_WL), "ssl_shift_meV": float(shift),
            "reischle_EX_eV": float(Lr.E_X_eV), "pryor_Ve_meV": float(Lr.V_e),
            "laf_EX_eV": float(Lt.E_X_eV), "gu_EX_eV": float(Lg.E_X_eV)}


def ev_nitride():
    from fsim_core import nitride_levels as NL
    from fsim_core.nitride_levels import NitrideDotSystem, levels, rates
    desh = NitrideDotSystem(height_nm=NL.DESHPANDE2013_HEIGHT_NM,
                            radius_nm=NL.DESHPANDE2013_DIAMETER_NM / 2., x_in=NL.DESHPANDE2013_X_IN)
    zh = NitrideDotSystem(height_nm=NL.ZHANG2013_HEIGHT_NM,
                          radius_nm=NL.ZHANG2013_DIAMETER_NM / 2., x_in=NL.ZHANG2013_X_IN)
    lv_d, lv_z = levels(desh, 300.), levels(zh, 300.)
    if not (lv_d.valid and lv_z.valid):
        raise RuntimeError(f"planar nitride levels not bound: {lv_d.invalid_reasons} {lv_z.invalid_reasons}")
    r_z = rates(lv_z, 300., tau_rad0_ns=NL.TAU_RAD0_DEFAULT_NS)
    from scripts import run_nitride_cavity as sweep
    row = sweep._row("C000", "comparison", "rectangular", {"T_hs": 300.}, {}, True)
    return {"desh2013_lambda_nm": float(lv_d.lambda_nm), "zhang2013_tau_ns": float(1. / r_z["gamma_X0_ns"]),
            "desh2014_g2": float(row["g2"])}


def ev_nanowire():
    from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels
    from fsim_core.nitride_nanowire_photonics import NitrideNanowirePhotonicsParams, response
    desh_sys = NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.25, strain_bound="relaxed")
    lt = levels(desh_sys, 25., z_points=4807, exterior_nm=45.)  # as verify_nitride_nanowire_levels.py (M3)
    if not lt.valid:
        raise RuntimeError(f"nanowire levels not valid: {lt.invalid_reasons}")
    # Claudon-matched GaAs reference, parameters as verify_nitride_nanowire_photonics.py:821-826
    p = NitrideNanowirePhotonicsParams(family="vertical_photonic", n_wire=3.45, n_group_override=4.3,
                                       NA=0.75, bottom_reflectivity=1.0, taper_output_mfr_nm=750.0)
    r = response(p, lambda_nm=950.0, outer_radius_nm=0.22 * 950.0 / 2.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    return {"desh2013_lambda_nm": float(lt.lambda_nm), "claudon_eta": float(r["eta_collection_X"])}


EVALUATORS = {"chatz_va": ev_chatz_va, "chatz_dev_ibm": ev_chatz_dev_ibm, "laf_mid": ev_laf_mid,
              "reischle": ev_reischle, "zhao": ev_zhao, "dot_levels": ev_dot_levels,
              "nitride": ev_nitride, "nanowire": ev_nanowire}


# ----------------------------------------------------------------------------- judging
def judge(item, model, se=None):
    """Return (outcome 'pass'|'fail', paper-text, tol-text).  Tolerance logic only; no model-derived numbers."""
    paper, tol = item["paper"], item["tol"]
    kind = tol["kind"]
    if kind == "abs":
        ok = abs(model - paper["value"]) <= tol["value"] + 1e-12
        return ok, f"{paper['value']}", f"+-{tol['value']}"
    if kind == "upper_bound":
        ok = model <= paper["value"] + tol["value"] + 1e-12
        return ok, f"<{paper['value']}", f"<= {paper['value'] + tol['value']:.4g}"
    if kind == "range":
        ok = paper["lo"] - 1e-12 <= model <= paper["hi"] + 1e-12
        return ok, f"[{paper['lo']},{paper['hi']}]", "in range"
    if kind == "energy_meV":
        dE = abs(HC_EV_NM / model - HC_EV_NM / paper["value"]) * 1e3
        return dE <= tol["value"] + 1e-9, f"{paper['value']}", f"+-{tol['value']} meV (dE={dE:.0f})"
    if kind == "rel":
        ok = abs(model - paper["value"]) / abs(paper["value"]) <= tol["value"] + 1e-12
        return ok, f"{paper['value']}", f"+-{tol['value']:.0%}"
    if kind == "combined_sigma":
        if se is None:
            raise ValueError(f"{item['id']}: combined_sigma needs a model standard error")
        lim = tol.get("k", 1.0) * math.hypot(se, paper["err"])
        return abs(model - paper["value"]) <= lim, f"{paper['value']}+-{paper['err']}", f"+-{lim:.4f} ({tol.get('k', 1.0):g} sigma)"
    if kind == "regression":
        ok = abs(model - paper["pin"]) <= tol["value"] + 1e-12
        return ok, f"pin {paper['pin']}", f"+-{tol['value']}"
    raise ValueError(f"{item['id']}: unknown tol kind {kind!r}")


def load_ledger(path=LEDGER_PATH):
    led = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    ids = [it["id"] for it in led["items"]]
    assert len(ids) == len(set(ids)), "duplicate ledger ids"
    for it in led["items"]:
        for k in ("class", "tol", "tol_reason", "expect", "paper", "source"):
            assert k in it, f"{it['id']}: ledger item missing {k}"
        assert it["expect"] in ("pass", "fail"), it["id"]
        assert it["class"] in ("C", "U", "H", "KD"), it["id"]
        if it["expect"] == "fail" and not str(it.get("reason", "")).strip():
            raise AssertionError(f"{it['id']}: expect: fail needs a reason")
    return led


def evaluate_ledger(led, run_all=False, verbose=True):
    """Run evaluators lazily per group, judge every item.  Returns list of result dicts."""
    cache, errors = {}, {}
    results = []
    for it in led["items"]:
        group, key = it["source"].split(".", 1)
        res = {"id": it["id"], "class": it["class"], "expect": it["expect"], "item": it}
        if (it.get("slow") or group in SLOW_GROUPS) and not run_all:
            res.update(outcome="skipped", model=None, detail="slow (use --all)")
            results.append(res)
            continue
        if group not in cache and group not in errors:
            t0 = time.time()
            try:
                cache[group] = EVALUATORS[group]()
                if verbose:
                    print(f"  [group {group}: {time.time() - t0:.1f} s]", flush=True)
            except Exception as exc:  # an evaluator crash is a mismatch, never silently skipped
                errors[group] = f"{type(exc).__name__}: {exc}"
        if group in errors:
            res.update(outcome="error", model=None, detail=errors[group])
        else:
            try:
                model = cache[group][key]
                ok, ptxt, ttxt = judge(it, model, cache[group].get(key + "__se"))
                res.update(outcome="pass" if ok else "fail", model=model, paper_txt=ptxt, tol_txt=ttxt)
            except Exception as exc:
                res.update(outcome="error", model=None, detail=f"{type(exc).__name__}: {exc}")
        results.append(res)
    return results


def exit_code(results):
    return 0 if all(r["outcome"] == r["expect"] for r in results if r["outcome"] != "skipped") else 1


def report(results):
    print(f"{'id':26s} {'cl':>2s} {'paper':>16s} {'model':>10s} {'tol':>26s} {'outcome':>7s} {'expect':>6s}  match")
    for r in results:
        if r["outcome"] in ("skipped", "error"):
            m = "SKIP" if r["outcome"] == "skipped" else "MISMATCH"
            print(f"{r['id']:26s} {r['class']:>2s} {'-':>16s} {'-':>10s} {'-':>26s} {r['outcome']:>7s} {r['expect']:>6s}  {m}  ({r['detail']})")
            continue
        m = "MATCH" if r["outcome"] == r["expect"] else "MISMATCH"
        print(f"{r['id']:26s} {r['class']:>2s} {r['paper_txt']:>16s} {r['model']:10.4f} {r['tol_txt']:>26s} "
              f"{r['outcome']:>7s} {r['expect']:>6s}  {m}")
    ran = [r for r in results if r["outcome"] != "skipped"]
    skipped = len(results) - len(ran)
    matched = [r for r in ran if r["outcome"] == r["expect"]]
    P = sum(r["outcome"] == "pass" for r in matched)
    F = sum(r["outcome"] == "fail" for r in matched)
    print("\nper-class table (matched/total; expected-pass / expected-fail among matched):")
    for c, name in (("C", "calibration"), ("U", "unfitted"), ("H", "held-out"), ("KD", "known-deviation")):
        sub = [r for r in ran if r["class"] == c]
        if not sub:
            continue
        mt = [r for r in sub if r["outcome"] == r["expect"]]
        print(f"  {c:>2s} {name:16s} {len(mt)}/{len(sub)}   pass {sum(r['outcome'] == 'pass' for r in mt)}   "
              f"known-fail {sum(r['outcome'] == 'fail' for r in mt)}")
    if skipped:
        print(f"  ({skipped} slow items skipped, not counted as passed)")
    for r in ran:
        if r["outcome"] != r["expect"]:
            print(f"MISMATCH {r['id']}: outcome {r['outcome']} but ledger expects {r['expect']}")
    print(f"{len(matched)}/{len(ran)} reproduce checks matched ledger ({P} pass, {F} known-fail)")


def self_test(led):
    results = evaluate_ledger(led, run_all=False, verbose=False)
    base = exit_code(results)
    flipped = copy.deepcopy(results)
    victim = next(r for r in flipped if r["outcome"] != "skipped")
    victim["expect"] = "fail" if victim["expect"] == "pass" else "pass"
    code = exit_code(flipped)
    ok = base == 0 and code == 1
    print(f"self-test: baseline exit {base}; after flipping expect of {victim['id']} in memory exit {code}")
    print("self-test " + ("OK (flip detected)" if ok else "FAILED"))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    led = load_ledger()
    t0 = time.time()
    if "--self-test" in argv:
        rc = self_test(led)
        print(f"runtime {time.time() - t0:.1f} s")
        return rc
    results = evaluate_ledger(led, run_all="--all" in argv)
    report(results)
    print(f"runtime {time.time() - t0:.1f} s")
    return exit_code(results)


if __name__ == "__main__":
    sys.exit(main())
