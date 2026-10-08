"""FSIM Studio Library (literature parameter cards) API checks (Flask test client, no network).

  1. all 10 shipped parameter cards load through GET /api/params/<name>;
  2. per-tag counts equal the counts of fsim_core.card.load_card (raw YAML for the two device-tier
     cards load_card cannot read; the API says so in `loader`);
  3. every param equals the card value / range exactly; a range is never collapsed to a value;
     dataset rows equal the card rows;
  4. every model overlay equals a direct fsim_core call made here, exactly (==, not approx), and
     overlays exist exactly where the rules file allows;
  5. PUT round-trip into a temporary copy of cards/ (the repository's cards/ is never written):
     409 on shipped names, value<->range swaps per the legacy editor schema, 403 without the header.

Run: python verify/verify_studio_library.py   (exit code 0 iff all pass)
"""
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CARDS = ["chatzarakis", "chatzarakis2023", "laferriere", "laferriere2023", "reischle", "kitamura", "zhao",
         "qcap-cavity", "qcap-staged", "qcap-piezo-variant"]
HDR = {"X-FSIM-Studio": "1"}
# The rules file (physics-brief.md section 8): class of each card and exactly where an overlay is allowed.
EXPECT_CLASS = {"chatzarakis": "C", "chatzarakis2023": "C", "laferriere": "M", "laferriere2023": "M",
                "reischle": "M", "kitamura": "T", "zhao": "C", "qcap-cavity": "P", "qcap-staged": "P",
                "qcap-piezo-variant": "P"}
EXPECT_OVERLAYS = {("chatzarakis", "g2_vs_T"), ("chatzarakis", "gamma_vs_T"), ("chatzarakis", "tau_vs_T"),
                   ("chatzarakis2023", "g2_vs_T"), ("laferriere", "g2_vs_T"), ("laferriere2023", "g2_vs_T"),
                   ("reischle", "rho_spectral"), ("reischle", "g2_oe2008"), ("zhao", "g2_vs_pump"),
                   ("qcap-piezo-variant", "device_points")}
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def counts_of(tags):
    out = {t: 0 for t in ("V", "DR", "E", "A")}
    for t in tags:
        out[t] += 1
    return out


def same_rows(a, b):
    if len(a) != len(b):
        return False
    for ra, rb in zip(a, b):
        if set(ra) != set(rb):
            return False
        for k in ra:
            va, vb = ra[k], rb[k]
            if isinstance(va, float) and isinstance(vb, float) and math.isnan(va) and math.isnan(vb):
                continue
            if va != vb:
                return False
    return True


def exact_equal(api_obj, direct_obj, path="overlay"):
    """Exact structural equality of the overlay JSON against the direct call (floats compared ==)."""
    if isinstance(direct_obj, dict):
        if not isinstance(api_obj, dict):
            return f"{path}: not an object"
        for k, v in direct_obj.items():
            if k not in api_obj:
                return f"{path}.{k}: missing"
            r = exact_equal(api_obj[k], v, f"{path}.{k}")
            if r:
                return r
        return None
    if hasattr(direct_obj, "tolist"):
        direct_obj = direct_obj.tolist()
    if isinstance(direct_obj, (list, tuple)):
        if not isinstance(api_obj, list) or len(api_obj) != len(direct_obj):
            return f"{path}: length {len(api_obj) if isinstance(api_obj, list) else 'n/a'} != {len(direct_obj)}"
        for i, (x, y) in enumerate(zip(api_obj, direct_obj)):
            r = exact_equal(x, y, f"{path}[{i}]")
            if r:
                return r
        return None
    if isinstance(direct_obj, float) and not math.isfinite(direct_obj):
        return None if api_obj is None else f"{path}: {api_obj!r} != null (non-finite)"
    if isinstance(direct_obj, float) or isinstance(api_obj, float):
        return None if float(api_obj) == float(direct_obj) else f"{path}: {api_obj!r} != {direct_obj!r}"
    return None if api_obj == direct_obj else f"{path}: {api_obj!r} != {direct_obj!r}"


def main():
    from fsim_core.card import load_card
    from fsim_studio import api_library
    from fsim_studio.server import create_app

    direct = direct_overlays()

    app = create_app()
    c = app.test_client()

    # ------------------------------------------------------------ list + load
    r = c.get("/api/params")
    lst = r.get_json() or []
    names = [x["name"] for x in lst]
    check("GET /api/params lists the 10 shipped parameter cards", r.status_code == 200 and all(n in names for n in CARDS),
          f"{len(lst)} listed")
    check("GET /api/params lists no design card", not any(n.endswith("-design") for n in names))
    for n in CARDS:
        r = c.get(f"/api/params/{n}")
        d = r.get_json() or {}
        check(f"{n}: loads (200)", r.status_code == 200, d.get("error", ""))
        if r.status_code != 200:
            continue
        path = ROOT / "cards" / f"{n}.yaml"
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        card = load_card(path)   # all 10 load now, the two device-tier cards included
        check(f"{n}: loaded by fsim_core.card.load_card (no raw fallback)", d["loader"] == "fsim_core.card.load_card",
              d["loader"])
        check(f"{n}: class {EXPECT_CLASS[n]} with a label and a reason (rules file)",
              d["cls"]["code"] == EXPECT_CLASS[n] and d["cls"]["label"] != "unclassified" and bool(d["cls"].get("reason")),
              str(d["cls"].get("label")))
        if d["has_device_block"]:
            dev_rng = [x for x in d["device_block"] if x.get("range")]
            check(f"{n}: device-tier ranges stay [lo, hi]; range centres are not shown",
                  all(x["value"] is None and len(x["range"]) == 2 for x in dev_rng)
                  and not any(x["path"].split(".")[-1].startswith("mid") for x in d["device_block"] if not x.get("hidden")))
            check(f"{n}: datasets without their own source say they use meta.source",
                  all(ds["source_note"] for ds in d["datasets"]))
        if card is not None:
            want = counts_of(p.tag.name for p in card.params.values())
            check(f"{n}: tag counts match fsim_core.card", d["tag_counts"] == want, f"{d['tag_counts']} vs {want}")
            ok = len(d["params"]) == len(card.params)
            for p in d["params"]:
                cp = card.params[p["name"]]
                if cp.range is not None:
                    ok &= p["range"] == [cp.range[0], cp.range[1]] and p["value"] is None
                else:
                    ok &= p["value"] == cp.value and p["range"] is None
                ok &= p["tag"] == cp.tag.name and p["unit"] == cp.unit and p["source"] == cp.source
            check(f"{n}: every param equals the card (ranges stay [lo, hi], never collapsed)", ok)
            ok = len(d["datasets"]) == len(card.datasets)
            for ds in d["datasets"]:
                cd = card.datasets[ds["name"]]
                ok &= ds["tag"] == cd.tag.name and same_rows(ds["rows"], cd.rows) and ds["n_rows"] == len(cd.rows)
            check(f"{n}: datasets, tags and row counts equal fsim_core.card", ok,
                  ", ".join(f"{x['name']}:{x['n_rows']}" for x in d["datasets"]))

    # ------------------------------------------------------------ overlays
    allowed = set(api_library.OVERLAYS)
    check("overlays exist exactly where the rules file allows (10 of them)", allowed == EXPECT_OVERLAYS,
          f"extra {sorted(allowed - EXPECT_OVERLAYS)} missing {sorted(EXPECT_OVERLAYS - allowed)}")
    for n in ("kitamura", "qcap-cavity", "qcap-staged"):
        st = (c.get(f"/api/params/{n}").get_json() or {}).get("overlays", {})
        check(f"{n}: no overlay and a card-level reason in the rules file",
              n in api_library.NO_OVERLAY and not any(k[0] == n for k in allowed) and not any(v["available"] for v in st.values()))
    check("reischle/g2_vs_ERR: refused with the rules-file reason (no g2_of_T overlay)",
          c.get("/api/params/reischle/overlay/g2_vs_ERR").status_code == 404
          and "temporal refilling" in (c.get("/api/params/reischle/overlay/g2_vs_ERR").get_json().get("reason") or ""))
    r = c.get("/api/params/chatzarakis-edited/overlay/g2_vs_T")
    check("an edited copy gets no overlay", r.status_code == 404)
    for n in CARDS:
        d = c.get(f"/api/params/{n}").get_json()
        for ds in d["datasets"]:
            key = (n, ds["name"])
            st = d["overlays"][ds["name"]]
            check(f"overlay status {n}/{ds['name']}: {'allowed' if key in allowed else 'none'}",
                  st["available"] == (key in allowed) and (st["available"] or bool(st.get("reason"))), st.get("reason", ""))
            r = c.get(f"/api/params/{n}/overlay/{ds['name']}")
            if key not in allowed:
                check(f"overlay {n}/{ds['name']}: refused 404 with a reason", r.status_code == 404 and r.get_json().get("reason"))
                continue
            ov = r.get_json()
            check(f"overlay {n}/{ds['name']}: 200", r.status_code == 200, (ov or {}).get("error", ""))
            if r.status_code != 200:
                continue
            fn = direct.get(key)
            if fn is None:
                check(f"overlay {n}/{ds['name']}: a direct fsim_core reference exists in this verify", False)
                continue
            diff = exact_equal(ov, fn())
            check(f"overlay {n}/{ds['name']}: equals the direct fsim_core call exactly", diff is None, diff or "")
            xs_data = [r[ds["plot"]["x"]] for r in ds["rows"] if isinstance(r.get(ds["plot"]["x"]), (int, float))] if ds.get("plot") else []
            if n in ("chatzarakis", "chatzarakis2023", "laferriere", "laferriere2023") and xs_data:
                ox = [x for s_ in ov.get("series", []) for x in s_["x"]] + [x for b in ov.get("bands", []) for x in b["x"]]
                check(f"overlay {n}/{ds['name']}: nothing extrapolates beyond the dataset's T range",
                      min(ox) >= min(xs_data) - 1e-9 and max(ox) <= max(xs_data) + 1e-9, f"{min(ox)}..{max(ox)} vs {min(xs_data)}..{max(xs_data)}")
            check(f"overlay {n}/{ds['name']}: bands are bands (lo <= hi) and series are finite",
                  all(all(l <= h for l, h in zip(b["lo"], b["hi"])) for b in ov.get("bands", []))
                  and all(all(v is not None for v in s_["y"]) for s_ in ov.get("series", [])))
            check(f"overlay {n}/{ds['name']}: labelled with its class and tag",
                  ov.get("cls", {}).get("code") == EXPECT_CLASS[n] and ov.get("tag") in ("V", "DR", "E", "A")
                  and ov.get("calls"), f"{ov.get('cls')} {ov.get('tag')}")
    z = c.get("/api/params/zhao/overlay/g2_vs_pump").get_json()
    check("zhao overlay: points with SE only (no curve through the quiet point), both rows verdict FAIL, caveat pinned",
          all(s_["mode"] == "markers" and len(s_["err"]) == len(s_["x"]) for s_ in z["series"])
          and z["readouts"]["verdict"] == {"normal": "FAIL", "quiet": "FAIL"}
          and "not a single-photon source" in z["label"])
    # review finding 3: the Zhao verdict word is the committed CSV's, never a re-threshold of the live sigma
    import csv as _csv
    with open(ROOT / "out/zhao/zhao_fit_comparison.csv", newline="", encoding="utf-8") as fh:
        zrows = {r["pump"]: r for r in _csv.DictReader(fh)}
    check("review 3: zhao overlay verdict equals the committed zhao_fit_comparison.csv verdict column, with its source named",
          z["readouts"]["verdict"] == {k: zrows[k]["verdict"] for k in ("normal", "quiet")}
          and z["readouts"]["verdict_source"] == "out/zhao/zhao_fit_comparison.csv"
          and all(abs(z["readouts"]["sigma_committed"][k] - float(zrows[k]["sigma_ratio"])) < 1e-12 for k in ("normal", "quiet")))
    zdir = Path(tempfile.mkdtemp(prefix="fsim_zhao_csv_"))
    try:
        fake = zdir / "zhao_fit_comparison.csv"
        fake.write_text((ROOT / "out/zhao/zhao_fit_comparison.csv").read_text(encoding="utf-8").replace("FAIL", "PASS"), encoding="utf-8")
        orig_csv = getattr(api_library, "ZHAO_COMPARISON_CSV", None)
        api_library.ZHAO_COMPARISON_CSV = fake
        try:
            zf = c.get("/api/params/zhao/overlay/g2_vs_pump").get_json()
        finally:
            api_library.ZHAO_COMPARISON_CSV = orig_csv
        check("review 3: the verdict follows the CSV (a CSV saying PASS gives PASS even though this live run is > 1 sigma off): not recomputed",
              zf["readouts"]["verdict"] == {"normal": "PASS", "quiet": "PASS"} and max(zf["readouts"]["sigma"].values()) > 1.0,
              str(zf["readouts"]["verdict"]))
    finally:
        shutil.rmtree(zdir, ignore_errors=True)
    # review finding 4: committed CSVs are the source of the rho overlays; tau_d is what the published pairs imply
    ro = c.get("/api/params/reischle/overlay/rho_spectral").get_json()
    po = c.get("/api/params/qcap-piezo-variant/overlay/device_points").get_json()
    oe = c.get("/api/params/reischle/overlay/g2_oe2008").get_json()
    imp = [0.5 * 0.59 / (0.85 - 0.59), 0.5 * 0.57 / (0.75 - 0.57)]
    check("review 4: reischle rho overlay is read from out/phase1/vc_rho_window.csv and the piezo band from out/phase3/map_rho_required_300K.csv",
          "out/phase1/vc_rho_window.csv" in ro["params_source"] and "out/phase3/map_rho_required_300K.csv" in po["params_source"]
          and "map_rho_required_300K.csv" in " ".join(po["calls"]), f"{ro['params_source'][:40]} | {po['params_source'][:50]}")
    check("review 4: reischle OE-2008 band uses the tau_d the published pairs imply (1.13 and 1.58 ns, [DR]), not an invented 1-3 ns [A] sweep",
          all(abs(a - b) < 1e-9 for a, b in zip(oe["readouts"]["tau_d_implied_ns"], imp)) and oe["tag"] == "DR"
          and "1.13" in oe["label"] and "swept [A]" not in oe["label"], oe["label"][-160:])
    src_lib = (ROOT / "fsim_studio/api_library.py").read_text(encoding="utf-8")
    check("review 4: no physics formula left in api_library.py (no sqrt, no retention ratio, no F2 inversion inline)",
          all(tok not in src_lib for tok in ("np.sqrt", "math.sqrt", "retention(", "** 2", "RHO_TARGET_G2", "epsilon_narrow_filter"))
          and "studio_support" in src_lib)
    # review finding 5: the card's `conflict:` note is passed through verbatim
    zc = c.get("/api/params/zhao").get_json()
    sq = next(p_ for p_ in zc["params"] if p_["name"] == "squeezing_dB")
    raw_conf = yaml.safe_load((ROOT / "cards/zhao.yaml").read_text(encoding="utf-8"))["params"]["squeezing_dB"]["conflict"]
    zl = next(x for x in c.get("/api/params").get_json() if x["name"] == "zhao")
    check("review 5: zhao squeezing_dB carries the card's conflict note verbatim (3.1 dB vs 0.9 dB); other params none; list counts it",
          sq["conflict"] == " ".join(raw_conf.split()) and "0.9" in sq["conflict"] and "3.1 dB" in sq["conflict"]
          and sum(1 for p_ in zc["params"] if p_.get("conflict")) == 1 and zl["n_conflicts"] == 1
          and all("conflict" in p_ for p_ in zc["params"]), sq["conflict"][:60])
    cz = c.get("/api/params/chatzarakis/overlay/g2_vs_T").get_json()
    check("chatzarakis overlay: labelled calibration, not held out; T_c 249.1 K (247.1 K fixed windows)",
          "not held out" in cz["label"] and abs(cz["readouts"]["T_c_fit_windows_K"] - 249.1) < 0.1
          and abs(cz["readouts"]["T_c_fixed_w3_dx-0.7_K"] - 247.1) < 0.1)
    lf = c.get("/api/params/laferriere/overlay/g2_vs_T").get_json()
    check("laferriere overlay: a band (no line) that excludes the 4 K point",
          not lf["series"] and lf["bands"] and 4.0 not in lf["bands"][0]["x"])
    check("overlay: unknown dataset is a 404", c.get("/api/params/zhao/overlay/nope").status_code == 404)
    check("params: unknown card 404, design card 422, bad name 400",
          c.get("/api/params/no-such-card").status_code == 404
          and c.get("/api/params/staged-device-design").status_code == 422
          and c.get("/api/params/..%5Cx").status_code in (400, 404))

    # ------------------------------------------------------------ PUT round-trip (temporary cards dir)
    before = sorted(p.name for p in (ROOT / "cards").glob("*.yaml"))
    tmp = Path(tempfile.mkdtemp(prefix="fsim_lib_"))
    try:
        for p in (ROOT / "cards").glob("*.yaml"):
            shutil.copy2(p, tmp / p.name)
        api_library.CARDS_DIR = tmp
        r = c.put("/api/params/chatzarakis", json={"params": {"delta_xx": {"value": 6.0}}}, headers=HDR)
        check("PUT shipped name -> 409", r.status_code == 409, (r.get_json() or {}).get("error", ""))
        r = c.put("/api/params/chatzarakis-edited", json={"params": {"delta_xx": {"value": 6.0}}})
        check("PUT without the Studio header -> 403", r.status_code == 403)
        # review finding 9: a Studio-saved card (cards/studio/<name>.yaml) with the same name would be hidden by this
        # copy (a cards/ file wins in _find_card), so the save is refused (409) and nothing is written
        from fsim_studio import api_cards
        uc = tmp / "user_cards"
        uc.mkdir()
        (uc / "chatzarakis-edited.yaml").write_text("design: {}\nmeta: {saved_by: fsim-studio}\n", encoding="utf-8")
        orig_uc = api_cards.USER_CARDS
        api_cards.USER_CARDS = uc
        try:
            r = c.put("/api/params/chatzarakis-edited", json={"params": {"delta_xx": {"value": 6.0}}}, headers=HDR)
        finally:
            api_cards.USER_CARDS = orig_uc
        check("review 9: PUT <name>-edited is refused (409) when a Studio-saved card already has that name; nothing is written",
              r.status_code == 409 and "Studio" in (r.get_json() or {}).get("error", "") and not (tmp / "chatzarakis-edited.yaml").exists(),
              str((r.get_json() or {}).get("error", ""))[:80])
        body = {"params": {
            "delta_xx": {"value": None, "lo": 5.5, "hi": 6.3},     # value -> range
            "E_lo": {"value": 30.0, "lo": 15.0, "hi": 37.0},       # range -> value (value wins, app.py rule)
            "gamma0": {"value": None, "lo": 0.1, "hi": 1.5},       # range -> new range
        }}
        r = c.put("/api/params/chatzarakis-edited", json=body, headers=HDR)
        check("PUT chatzarakis-edited -> 200", r.status_code == 200, str(r.get_json()))
        f = tmp / "chatzarakis-edited.yaml"
        ed = load_card(f) if f.is_file() else None
        orig = load_card(ROOT / "cards" / "chatzarakis.yaml")
        ok = ed is not None
        if ed is not None:
            ok &= ed["delta_xx"].range == (5.5, 6.3) and ed["delta_xx"].value is None
            ok &= ed["E_lo"].value == 30.0 and ed["E_lo"].range is None
            ok &= ed["gamma0"].range == (0.1, 1.5)
            ok &= all(ed[k].tag == orig[k].tag and ed[k].source == orig[k].source and ed[k].unit == orig[k].unit
                      for k in orig.params)
            ok &= all((ed[k].value, ed[k].range) == (orig[k].value, orig[k].range)
                      for k in orig.params if k not in body["params"])
            ok &= set(ed.datasets) == set(orig.datasets)
            rawe = yaml.safe_load(f.read_text(encoding="utf-8"))
            ok &= list(rawe) == list(yaml.safe_load((ROOT / "cards" / "chatzarakis.yaml").read_text(encoding="utf-8")))
        check("PUT round-trip: value replaces range, range replaces value; tags/units/sources/datasets kept; same schema", ok)
        g = c.get("/api/params/chatzarakis-edited").get_json()
        dx = next(p for p in g["params"] if p["name"] == "delta_xx")
        check("GET edited copy: the new range comes back as [lo, hi], not a value",
              dx["range"] == [5.5, 6.3] and dx["value"] is None and g["edited"] and g["base"] == "chatzarakis")
        r = c.put("/api/params/chatzarakis-edited", json={"params": {"delta_xx": {"value": None, "lo": 7.0, "hi": 6.0}}}, headers=HDR)
        check("PUT lo > hi -> 400", r.status_code == 400)
        r = c.put("/api/params/chatzarakis-edited", json={"params": {"nope": {"value": 1.0}}}, headers=HDR)
        check("PUT unknown parameter -> 400", r.status_code == 400)
        r = c.put("/api/params/chatzarakis2023-edited", json={"params": {"x": {"value": 1.0}}}, headers=HDR)
        check("PUT on a card without a params block -> 422", r.status_code == 422)
        r = c.put("/api/params/zhao-edited", json={"params": {"C_dep_pF": {"value": "abc"}}}, headers=HDR)
        check("PUT non-numeric value -> 400", r.status_code == 400)
    finally:
        api_library.CARDS_DIR = ROOT / "cards"
        shutil.rmtree(tmp, ignore_errors=True)
    after = sorted(p.name for p in (ROOT / "cards").glob("*.yaml"))
    check("the repository's cards/ was not written", before == after)

    n_pass = sum(ok for _, ok, _ in CHECKS)
    print(f"\n{n_pass}/{len(CHECKS)} studio library checks passed")
    return 0 if n_pass == len(CHECKS) else 1


def direct_overlays():
    """Direct fsim_core references for every overlay the rules file allows: {(card, dataset): fn}.
    Each fn re-derives the curve from the card file with public fsim_core calls only (the
    api_library builders are not imported), following the producing scripts named in the brief."""
    import itertools

    import numpy as np

    from fsim_core.card import load_card
    from fsim_core.cw_g2 import intrinsic_from_raw
    from fsim_core.device import DeviceDesign, class_proxy_params, evaluate
    from fsim_core.fitting import fit_phase0
    from fsim_core.integrator import g2_of_T, retention, rho_of_T, solve_Tc
    from fsim_core.loading import f1b_g2
    from fsim_core.sde import QDLaserParams, find_threshold, g2_vs_pump
    from fsim_core.spectral import epsilon, epsilon_narrow_filter, gamma_of_T

    cards = ROOT / "cards"
    memo = {}

    def fit():
        if "fit" not in memo:
            memo["fit"] = fit_phase0(load_card(cards / "chatzarakis.yaml"))
        return memo["fit"]

    def rows_of(name, ds):
        return yaml.safe_load((cards / f"{name}.yaml").read_text(encoding="utf-8"))["data"][ds]["points"]

    def chatz_g2():
        f = fit()
        Td = np.array([r["T"] for r in f.data])
        ws, dxs = np.array([r["w"] for r in f.data]), np.array([r["dx"] for r in f.data])
        p = dict(f.params)
        p["w"] = lambda T: np.interp(T, Td, ws)
        p["dx"] = lambda T: np.interp(T, Td, dxs)
        Ts = np.linspace(Td.min(), Td.max(), 121)
        return {"series": [{"x": Ts, "y": [g2_of_T(T, p).g2 for T in Ts]}], "tag": f.tag.name,
                "readouts": {"T_c_fit_windows_K": solve_Tc(p),
                             "T_c_fixed_w3_dx-0.7_K": solve_Tc(dict(f.params, w=3.0, dx=-0.7))}}

    def chatz_gamma():
        P = fit().params
        Td = [r["T"] for r in rows_of("chatzarakis", "gamma_vs_T")]
        Ts = np.linspace(min(Td), max(Td), 121)
        g = [float(gamma_of_T(T, P["gamma0"], P["a_ac"], P["b_lo"], P["E_lo"])) for T in Ts]
        return {"series": [{"x": Ts, "y": g}, {"x": Ts, "y": [P["r_xx"] * v for v in g]}]}

    def chatz_tau():
        P = fit().params
        rows = sorted(rows_of("chatzarakis", "tau_vs_T"), key=lambda r: r["T"])
        Ts = np.linspace(rows[0]["T"], rows[-1]["T"], 121)
        S = lambda T: retention(T, P["a_esc"], P["E_a"], P["b_p"], P["E_b"])
        return {"series": [{"x": Ts, "y": [rows[0]["tau_ns"] * float(S(T) / S(rows[0]["T"])) for T in Ts]}]}

    def design(ls, geom, alpha=None, rho_one=True):
        d = DeviceDesign()
        d.drive.V, d.drive.I_uA, d.drive.duty, d.drive.b_e, d.drive.mode = 0.0, 0.0, 0.0, 0.0, "PL"
        d.cavity.enabled = False
        d.filter.enabled, d.filter.auto_w = True, False
        d.dot.lineshape = ls
        d.dot.phonon = dict(geom, **({"alpha_ps2": alpha} if alpha is not None else {}))
        if rho_one:
            d.ret.b0 = d.ret.beta = 1e-12
        return d

    def point(d, T, w, dx, mu, delta, gt=None):
        d.thermal.T_hs, d.filter.w, d.filter.dx, d.drive.mu, d.dot.delta_xx = T, w, dx, mu, delta
        if gt is not None:
            p = class_proxy_params()
            d.dot.gamma_scale = gt / float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        return float(evaluate(d, T_grid=[T])["scalars"]["g2_op"])

    def chatz2023():
        raw = yaml.safe_load((cards / "chatzarakis2023.yaml").read_text(encoding="utf-8"))
        dev, pts = raw["device"], raw["data"]["g2_vs_T"]["points"]
        out = {"ibm": [], "lorentzian": []}
        for pt in pts:
            for ls in out:
                out[ls].append(point(design(ls, dev["geometry"], dev["phonon_alpha_ps2"], False), float(pt["T"]),
                                     float(pt["w"]), float(pt["dx"]), dev["mu"], dev["delta_xx"]))
        Ts = [float(p["T"]) for p in pts]
        return {"series": [{"x": Ts, "y": out["ibm"]}, {"x": Ts, "y": out["lorentzian"]}]}

    def laf2023():
        raw = yaml.safe_load((cards / "laferriere2023.yaml").read_text(encoding="utf-8"))
        dev, pts = raw["device"], raw["data"]["g2_vs_T"]["points"]
        mid, lo, hi = [], [], []
        for pt in pts:
            T = float(pt["T"])
            gam, win = dev["gamma_meV"][str(int(T))], dev["w_meV"][str(int(T))]
            mu, mr = ((dev["mu"]["mid_77"], dev["mu"]["range_77"]) if T < 100
                      else (dev["mu"]["mid_high"], dev["mu"]["range_high"]))
            d = design("ibm", dev["geometry"])
            mid.append(point(d, T, win["mid"], 0.0, mu, dev["delta_xx"]["mid"], gam["mid"]))
            c = [point(d, T, ww, 0.0, mm, dd, gg) for dd, gg, ww, mm in
                 itertools.product(dev["delta_xx"]["range"], gam["range"], win["range"], mr)]
            lo.append(float(np.min(c)))
            hi.append(float(np.max(c)))
        Ts = [float(p["T"]) for p in pts]
        return {"series": [{"x": Ts, "y": mid}], "bands": [{"x": Ts, "lo": lo, "hi": hi}]}

    def laf():
        card = load_card(cards / "laferriere.yaml")
        rows = sorted((r for r in card.datasets["g2_vs_T"].rows if "mechanism" not in r), key=lambda r: r["T"])
        G = lambda k, n: np.linspace(*card[k].bounds, n)
        lo, hi = [], []
        for r in rows:
            g = []
            for d in G("delta_xx", 9):
                if r["T"] < 100:
                    g += [float(f1b_g2(m, epsilon(d, ga, ga, w=w).eps)) for ga in G("gamma_77", 5)
                          for w in G("w_77", 5) for m in G("mu_77", 3)]
                elif r["T"] < 260:
                    g += [float(f1b_g2(m, epsilon(d, ga, ga, w=w).eps)) for ga in G("gamma_220", 6)
                          for w in G("w_220", 5) for m in G("mu_high", 6)]
                else:
                    g += [float(f1b_g2(m, epsilon(d, ga, ga, w=card["w_300"].fixed).eps)) for ga in G("gamma_300", 6)
                          for m in G("mu_high", 6)]
            lo.append(min(g))
            hi.append(max(g))
        return {"bands": [{"x": [r["T"] for r in rows], "lo": lo, "hi": hi}]}

    def reischle_rho():
        card = load_card(cards / "reischle.yaml")
        base = next(r for r in card.datasets["g2_vs_ERR"].rows
                    if r["device"] == 1 and r["position"] == 2 and r["ERR_MHz"] == 100.0)
        import csv
        with open(ROOT / "out/phase1/vc_rho_window.csv", newline="", encoding="utf-8") as fh:
            win = sorted(csv.DictReader(fh), key=lambda r: float(r["w_meV"]))
        req = [float(r["rho_required"]) for r in win]
        # independent closed form: rho_required = sqrt(1 - g2) of the 100 MHz point (eps = 0, g2 = 1 - rho^2)
        assert all(abs(v - float(np.sqrt(1.0 - base["g2"]))) < 1e-12 for v in req), "committed CSV disagrees with sqrt(1 - g2)"
        return {"series": [{"x": [float(r["w_meV"]) for r in win], "y": req}]}

    def reischle_oe():
        rows = load_card(cards / "reischle.yaml").datasets["g2_oe2008"].rows
        lo, hi = [], []
        # tau_d implied by each published (raw, intrinsic) pair, 0.5 ns exponential IRF: tau_irf (1 - raw) / ((1 - g2_b) - (1 - raw))
        implied = [0.5 * (1.0 - r["g2_raw"]) / ((1.0 - r["g2_b"]) - (1.0 - r["g2_raw"])) for r in rows]
        for r in rows:
            v = [float(intrinsic_from_raw(r["g2_raw"], td, 500.0, "exponential")) for td in np.linspace(min(implied), max(implied), 9)]
            lo.append(min(v))
            hi.append(max(v))
        T = [r["T"] for r in rows]
        return {"series": [{"x": T, "y": [(r["g2_b"] - (1 - r["rho"] ** 2)) / r["rho"] ** 2 for r in rows]}],
                "bands": [{"x": T, "lo": lo, "hi": hi}]}

    def zhao():
        params = QDLaserParams()
        I_th = find_threshold(params)
        res = [g2_vs_pump(params, 4.0 * I_th, n_runs=60, t_end=2.0e-9, dt=1.0e-13, seed=42, F_pump=F)
               for F in (1.0, 0.08)]
        return {"series": [{"x": ["normal", "quiet"], "y": [r["g2_mean"] for r in res],
                            "err": [r["g2_se"] for r in res]}]}

    def piezo():
        card = load_card(cards / "qcap-piezo-variant.yaml")
        import csv
        with open(ROOT / "out/phase3/map_rho_required_300K.csv", newline="", encoding="utf-8") as fh:
            rows_csv = list(csv.DictReader(fh))
        deltas = np.array([float(r["delta_meV"]) for r in rows_csv])
        req = [np.array([float(r[k]) for r in rows_csv]) for k in ("rho_req_gamma6.0", "rho_req_gamma7.0")]
        # independent recomputation of the committed map: rho >= sqrt(0.5/(1 - eps_min(300 K))), NaN where rho > 1
        glo, ghi = card["gamma_300_lo"].fixed, card["gamma_300_hi"].fixed
        for g, col in zip((glo, ghi), req):
            e = np.array([epsilon_narrow_filter(d, g, g) for d in deltas])
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = np.sqrt(0.5 / (1.0 - e))
            rr = np.where(rr <= 1.0, rr, np.nan)
            assert np.array_equal(np.isnan(rr), np.isnan(col)) and np.allclose(rr[~np.isnan(rr)], col[~np.isnan(col)], rtol=1e-9, atol=0), "CSV disagrees with the closed form"
        ok = np.isfinite(req[0]) & np.isfinite(req[1])
        P = fit().params
        pts = card.datasets["device_points"].rows
        ys = [p["rho"] if p["rho"] > 0 else
              float(rho_of_T(p["T"], P["a_esc"], P["E_a"], P["b_p"], P["E_b"], P["b0"], P["beta"])) for p in pts]
        return {"series": [{"x": [p["delta"] for p in pts], "y": ys}],
                "bands": [{"x": deltas[ok], "lo": np.minimum(req[0], req[1])[ok],
                           "hi": np.maximum(req[0], req[1])[ok]}]}

    return {
        ("chatzarakis", "g2_vs_T"): chatz_g2, ("chatzarakis", "gamma_vs_T"): chatz_gamma,
        ("chatzarakis", "tau_vs_T"): chatz_tau, ("chatzarakis2023", "g2_vs_T"): chatz2023,
        ("laferriere", "g2_vs_T"): laf, ("laferriere2023", "g2_vs_T"): laf2023,
        ("reischle", "rho_spectral"): reischle_rho, ("reischle", "g2_oe2008"): reischle_oe,
        ("zhao", "g2_vs_pump"): zhao, ("qcap-piezo-variant", "device_points"): piezo,
    }


if __name__ == "__main__":
    sys.exit(main())
