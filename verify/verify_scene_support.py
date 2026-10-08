"""fsim_core.scene_support: the physics that used to live in the Flask layer
(fsim_studio/scene.py cascade / planar / band builders and the
fsim_studio/api_explain.py fallback constants) moved into fsim_core helpers.

The proof is output identity: verify/data/scene_support_baseline.json holds
the SceneSpecs / explain payloads the PRE-refactor code produced (captured
with --capture before the move, 2026-10-07). Every case must be identical
(canonical JSON, == on floats) to what the refactored code returns now. The
single documented exception is the mu=0 fix in api_explain.cascade (a fitted
mu of 0 was turned into 0.33 by `mu or 0.33`), checked separately on a
patched fit bundle.

Also checked: DBR indices / pairs in the planar SceneSpec come from
fsim_core.nitride_cavity.dbr_diagnostic's own defaults (no literals), and
api_explain answers 400 with the missing list when a fit parameter is absent.

Run: python verify/verify_scene_support.py            (exit 0 iff all pass)
     python verify/verify_scene_support.py --capture  (rewrite the baseline;
     only ever run on the pre-refactor code)
"""
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BASELINE = ROOT / "verify" / "data" / "scene_support_baseline.json"

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def canon(obj):
    return json.dumps(obj, sort_keys=True, allow_nan=True, separators=(",", ":"))


SCENE_CASES = [
    ("cascade", "edge-inp-gainp-design", None),
    ("cascade", "edge-inp-gainp-design", 230.0),
    ("cascade", "nitride-cavity-set-design", None),
    ("cascade", "nitride-cavity-pulse-design", None),
    ("cascade", "nitride-nonpolar-set-design", 230.0),
    ("cascade", "nitride-nanowire-vertical-set-design", None),
    ("cascade", "nitride-nanowire-horizontal-set-design", 273.0),
    ("device", "nitride-cavity-set-design", None),
    ("device", "nitride-cavity-pulse-design", 230.0),
    ("device", "nitride-nonpolar-set-design", None),
    ("band", "nitride-cavity-set-design", None),
    ("band", "nitride-nonpolar-set-design", 230.0),
    ("band", "nitride-cavity-pulse-design", None),
]
SPECTRAL_CASES = [
    {"T": 78.0, "pub": True},
    {"T": 230.0, "pub": True},
    {"T": 200.0, "w": 2.5, "dx": 0.4, "mu": 0.0, "mode": "EL", "dg_inj": 3.0, "I_ratio": 2.0, "pub": False},
    {"T": 150.0, "w": 1.0, "dx": -0.5, "mu": 0.8, "mode": "PL", "pub": False},
]
CASCADE_TS = [150.0, 230.0, 300.0]
ADDED_KEYS = {"param_sources"}  # documented additive explain key (not a changed value)


def collect():
    from fsim_studio import api_explain, scene
    out = {}
    for kind, card, T in SCENE_CASES:
        fn = getattr(scene, f"{kind}_scene")
        out[f"{kind}|{card}|{T}"] = fn(card, T)
    for i, kw in enumerate(SPECTRAL_CASES):
        out[f"spectral|{i}"] = api_explain.spectral(**kw)
    for T in CASCADE_TS:
        out[f"explain_cascade|{T}"] = api_explain.cascade(T=T)
    return json.loads(canon(out))


def main():
    os.environ.setdefault("FSIM_STUDIO_CACHE", tempfile.mkdtemp(prefix="fsim_scene_support_"))
    if "--capture" in sys.argv:
        data = collect()
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        BASELINE.write_text(canon(data), encoding="utf-8")
        print(f"captured {len(data)} cases -> {BASELINE.relative_to(ROOT)}")
        return 0

    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    now = collect()
    for key in base:
        a, b = base[key], now.get(key)
        if isinstance(b, dict) and key.startswith("spectral|"):
            # the one additive key: where each fit parameter came from (p_inj source stated)
            b = {k: v for k, v in b.items() if k not in ADDED_KEYS}
        same = canon(a) == canon(b)
        detail = ""
        if not same and isinstance(a, dict) and isinstance(b, dict):
            diff = sorted(k for k in set(a) | set(b) if canon(a.get(k)) != canon(b.get(k)))
            detail = "differs in " + ",".join(diff)
        check(f"identical to pre-refactor output: {key}", same,
              detail or hashlib.sha256(canon(b).encode()).hexdigest()[:10])

    # The moved code now lives in fsim_core.scene_support and the Flask layer calls it.
    import inspect

    from fsim_core import nitride_cavity, scene_support
    from fsim_studio import api_explain, scene
    src_scene = inspect.getsource(scene)
    src_explain = inspect.getsource(api_explain)
    check("scene.py builds no Lindblad drive itself (lindblad.* calls moved to fsim_core.scene_support)",
          "lindblad.build_system" not in src_scene.replace('"fsim_core.lindblad.build_system', "")
          and "lindblad.evolve(" not in src_scene and not re.search(r"\*\s*1e-9|1e9\s*/", src_scene)
          and "scene_support.cascade_trajectory" in src_scene)
    check("scene.py carries no DBR index / pair literals (2.1, 1.46, pairs = 10)",
          "2.1," not in src_scene and "1.46" not in src_scene and "pairs = 10" not in src_scene
          and "scene_support.planar_cavity_optics" in src_scene)
    check("scene.py band tilt via fsim_core.scene_support.band_tilt_eV (no 1e-4 field factor)",
          "1e-4" not in src_scene and "scene_support.band_tilt_eV" in src_scene)
    check("api_explain.py has no fallback physics constants (5.9, 0.5, 2e-3, 25.0, 36.6, 0.33)",
          not any(lit in src_explain for lit in ("5.9", "2e-3", "25.0", "36.6", "0.33", '"gamma0", 0.5'))
          and "scene_support.explain_fit_params" in src_explain)

    dflt = inspect.signature(nitride_cavity.dbr_diagnostic).parameters
    sp = now["device|nitride-cavity-set-design|None"]
    hl = [L for L in sp["stacks"]["optical"]["layers"] if L["name"].startswith("H")]
    ll = [L for L in sp["stacks"]["optical"]["layers"] if L["name"].startswith("L")]
    check("planar DBR n_high / n_low / pairs equal nitride_cavity.dbr_diagnostic defaults",
          all(L["n"] == dflt["n_high"].default for L in hl) and all(L["n"] == dflt["n_low"].default for L in ll)
          and len(hl) == dflt["pairs"].default == sp["overlays"]["dbr"]["pairs"],
          f"n_high {dflt['n_high'].default}, n_low {dflt['n_low'].default}, pairs {dflt['pairs'].default}")
    # single source: a different dbr default must flow into the SceneSpec (no literal in scene.py)
    orig_dd = scene_support.dbr_defaults
    scene_support.dbr_defaults = lambda: {"n_high": 2.3, "n_low": 1.5, "pairs": 7}
    try:
        sp2 = scene.device_scene("nitride-cavity-set-design")
    finally:
        scene_support.dbr_defaults = orig_dd
    L2 = sp2["stacks"]["optical"]["layers"]
    check("planar DBR stack follows fsim_core's DBR definition (patched 2.3 / 1.5 / 7 pairs flows through)",
          len(L2) == 14 and {L["n"] for L in L2 if L["name"].startswith("H")} == {2.3}
          and {L["n"] for L in L2 if L["name"].startswith("L")} == {1.5} and sp2["overlays"]["dbr"]["pairs"] == 7)

    # mu = 0 bug fix: a fitted mu of exactly 0 must stay 0 in the cascade explainer.
    orig_fit = api_explain._fit
    fit0 = json.loads(json.dumps(orig_fit()))
    fit0["params"]["mu"] = 0.0
    api_explain._fit = lambda: fit0
    try:
        d0 = api_explain.cascade(T=150.0)
    finally:
        api_explain._fit = orig_fit
    from fsim_core.loading import loading_probs
    check("cascade explainer keeps a fitted mu = 0 (was turned into 0.33)",
          d0["mu"] == 0.0 and (d0["P0"], d0["P1"], d0["P2"]) == tuple(float(x) for x in loading_probs(0.0)),
          f"mu {d0['mu']}, P0 {d0['P0']}")

    # Missing fit parameter -> 400 with the missing list (no silent physics fallback).
    from fsim_studio.server import create_app
    app = create_app()
    c = app.test_client()
    fitx = json.loads(json.dumps(orig_fit()))
    del fitx["params"]["gamma0"]
    del fitx["params"]["E_lo"]
    api_explain._fit = lambda: fitx
    try:
        r1 = c.get("/api/explain/spectral?T=100")
        r2 = c.get("/api/explain/cascade?T=100")
    finally:
        api_explain._fit = orig_fit
        app.config["JOB_MANAGER"].shutdown()
    j1 = r1.get_json()
    check("spectral explainer: missing fit params -> 400 with the missing list",
          r1.status_code == 400 and sorted(j1.get("missing") or []) == ["E_lo", "gamma0"], str(j1)[:120])
    check("cascade explainer: missing fit params -> 400",
          r2.status_code == 400 and "missing" in r2.get_json(), str(r2.get_json())[:120])
    sp0 = api_explain.spectral(T=200.0, mode="EL", pub=True)
    check("p_inj: absent from the fit bundle -> fsim_core.loading.gamma_eff default, source stated",
          sp0["p_inj"] == inspect.signature(__import__("fsim_core.loading", fromlist=["x"]).gamma_eff)
          .parameters["p_inj"].default and "gamma_eff" in sp0.get("param_sources", {}).get("p_inj", ""),
          str(sp0.get("param_sources", {}).get("p_inj")))

    # N4: the module docstring must not claim every helper "moved unchanged",
    # and must name the three helpers that are new in the move.
    import fsim_core.scene_support as ss
    doc = ss.__doc__ or ""
    flat = " ".join(doc.split())
    new_helpers = ("background_g2_floor", "is_background_floor", "explain_fit_params")
    check("scene_support docstring: no blanket 'moved ... unchanged' claim; names the new helpers as new",
          "unchanged" not in flat and "New in the move" in flat and all(n in flat for n in new_helpers),
          "missing: " + ",".join(n for n in new_helpers if n not in flat) if not all(n in flat for n in new_helpers) else "")
    # Independent values (hand arithmetic, not the helper): rho = 1/(1+b), g2 = 1 - rho^2.
    floor_ref = 0.17355371900826446  # 1 - (1/1.1)**2 for b_res = 0.1
    check("background_g2_floor(0.1) == 1-(1/1.1)^2 = 0.17355371900826446",
          abs(ss.background_g2_floor(0.1) - floor_ref) < 1e-12, repr(ss.background_g2_floor(0.1)))
    check("background_g2_floor(1.0) == 1-(1/2)^2 = 0.75", abs(ss.background_g2_floor(1.0) - 0.75) < 1e-12)
    check("is_background_floor: printed nitride floor 0.17355371900826433 at b_res 0.1 -> True",
          ss.is_background_floor(0.17355371900826433, 0.1) is True)
    check("is_background_floor: non-floor cases -> False (g2 0.5 / b_res 0.1; floor value at b_res 0.2; b_res 0; NaN)",
          ss.is_background_floor(0.5, 0.1) is False and ss.is_background_floor(floor_ref, 0.2) is False
          and ss.is_background_floor(0.0, 0.0) is False and ss.is_background_floor(float("nan"), 0.1) is False)

    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_pass}/{len(CHECKS)} scene_support checks passed")
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
