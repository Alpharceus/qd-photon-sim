"""FSIM Studio Composer + Papers API checks (devlib-p10; Flask test client, temp cache and card dirs).

  1. registry route == fsim_core.compose.registry_payload() (+ tiers_payload() when tiers.py exists);
  2. validate route == validate_spec(...) issues/resolved and eta_s == engine.estimate; bad specs -> issues;
  3. mutating routes without X-FSIM-Studio -> 403;
  4. run point / tsweep (g2, hom_I, B_lens, curves) == direct engine.evaluate_spec, bit for bit;
  5. series + reproduce are jobs: submit -> poll to done -> result == direct call; the second identical call is a
     cache hit with eta 0;
  6. /api/run with a composed design works (same numbers); envelope and headline -> 400;
  7. a composed design saves through PUT /api/cards/<slug> and round-trips;
  8. scene route == fsim_studio.scene_compose.describe (SKIP while that module is absent);
  9. papers routes (SKIP with a reason while fsim_core.compose.papers or a paper card is absent).

Run: python verify/verify_studio_compose.py   (exit code 0 iff every non-skipped check passes)
"""
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

HDR = {"X-FSIM-Studio": "1"}
CHECKS = []
SKIPS = []

SPEC = {"schema": "fsim-composed/1", "T_K": 4.2,
        "material": {"kind": "inas_gaas_sk", "params": {}},
        "structure": {"kind": "generic_cavity", "params": {}},
        "drive": {"kind": "resonant_pulse", "params": {}},
        "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
RATE_SPEC = dict(SPEC, solver="rate", drive={"kind": "above_band", "params": {}})


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def skip(name, reason):
    SKIPS.append((name, reason))
    print(f"[SKIP] {name}  (reason: {reason})", flush=True)


def have(modname):
    try:
        return importlib.util.find_spec(modname) is not None
    except ModuleNotFoundError:
        return False


def main():
    tmp = tempfile.TemporaryDirectory(prefix="fsim_studio_compose_")
    os.environ["FSIM_STUDIO_CACHE"] = str(Path(tmp.name) / "cache")
    os.environ["FSIM_STUDIO_USER_CARDS"] = str(Path(tmp.name) / "usercards")
    try:
        run_checks(tmp)
    finally:
        tmp.cleanup()
    n_ok = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_ok}/{len(CHECKS)} studio compose checks passed" + (f" ({len(SKIPS)} skipped)" if SKIPS else ""))
    return 0 if n_ok == len(CHECKS) else 1


def run_checks(tmp):
    from fsim_core.compose import engine, registry, validate_spec
    from fsim_studio import api_compose, cache, jobs
    from fsim_studio.serialize import to_jsonable
    from fsim_studio.server import create_app

    app = create_app()
    c = app.test_client()
    j = lambda r: json.loads(r.get_data(as_text=True))  # noqa: E731
    J = lambda o: json.loads(json.dumps(to_jsonable(o)))  # noqa: E731

    def post(path, body, hdr=HDR):
        return c.post(path, json=body, headers=hdr)

    def wait_job(job_id, timeout=240):
        t0 = time.time()
        while time.time() - t0 < timeout:
            snap = j(c.get(f"/api/jobs/{job_id}"))
            if snap["state"] in ("done", "error", "cancelled"):
                return snap
            time.sleep(0.2)
        return {"state": "timeout"}

    tiers_mod = importlib.util.find_spec("fsim_core.compose.tiers") if have("fsim_core.compose.tiers") else None
    papers_present = have("fsim_core.compose.papers")

    # ------------------------------------------------------------ 1. registry
    r = c.get("/api/compose/registry")
    registry.ensure_loaded()
    want = registry.registry_payload()
    if tiers_mod is not None:
        from fsim_core.compose import tiers
        want = dict(want, tiers=tiers.tiers_payload())
    check("1 GET /api/compose/registry == registry_payload()" + (" + tiers_payload()" if tiers_mod else ""),
          r.status_code == 200 and j(r) == J(want), f"{len(j(r).get('materials', []))} materials")
    if tiers_mod is None:
        skip("1b registry tiers list", "fsim_core.compose.tiers not present; registry route returns the empty tiers list")
    else:
        check("1b registry tiers list is non-empty and equals tiers_payload()", len(j(r)["tiers"]) > 0)
    check("1c second registry call is served from the per-hash cache (same object)",
          api_compose.registry_payload_cached() is api_compose.registry_payload_cached())

    # ------------------------------------------------------------ 2. validate
    r = post("/api/compose/validate", {"spec": SPEC})
    vd = validate_spec(SPEC)
    body = j(r)
    check("2 validate: ok, issues, resolved equal validate_spec", r.status_code == 200 and body["ok"] is True
          and body["issues"] == J(vd["issues"]) and body["resolved"] == J(vd["resolved"])
          and body["solver"] == vd["solver"])
    check("2b validate eta_s == engine.estimate(spec, 1)", body["eta_s"] == float(engine.estimate(SPEC, 1)),
          f"{body['eta_s']:.4g} s")
    bad = json.loads(json.dumps(SPEC))
    bad["structure"]["kind"] = "no_such_structure"
    r = post("/api/compose/validate", {"spec": bad})
    vb = validate_spec(bad)
    check("2c unknown kind: ok false, issues == validate_spec issues, eta null",
          r.status_code == 200 and j(r)["ok"] is False and j(r)["issues"] == J(vb["issues"]) and j(r)["eta_s"] is None)
    r = post("/api/compose/validate", {"spec": {"schema": "nope"}})
    check("2d structurally broken spec -> 200 ok false with the error issue",
          r.status_code == 200 and j(r)["ok"] is False and j(r)["issues"][0]["level"] == "error")
    check("2e validate without spec -> 400", post("/api/compose/validate", {}).status_code == 400)

    # ------------------------------------------------------------ 3. header rule
    pids = "x"
    for path in ("/api/compose/validate", "/api/compose/run", "/api/compose/series", "/api/compose/scene",
                 f"/api/papers/{pids}/reproduce"):
        check(f"3 POST {path} without X-FSIM-Studio -> 403", c.post(path, json={"spec": SPEC}).status_code == 403)

    # ------------------------------------------------------------ 4. run point / tsweep
    direct = engine.evaluate_spec(SPEC)
    r = post("/api/compose/run", {"spec": SPEC, "mode": "point"})
    resp = j(r)
    if "result" not in resp and "job_id" in resp:
        snap = wait_job(resp["job_id"])
        res = snap.get("result")
    else:
        res = resp.get("result")
    sc = res["scalars"] if res else {}
    check("4 run point: g2_0, hom_I, B_lens, rate_detected_Hz == direct evaluate_spec",
          r.status_code == 200 and res is not None
          and all(sc.get(k) == J(direct["scalars"][k]) for k in ("g2_0", "hom_I", "B_lens", "rate_detected_Hz")),
          f"g2={sc.get('g2_0')}")
    check("4b run point: curves equal direct (single T) and carry run_id/labels/series",
          res is not None and res["curves"]["g2"] == J(direct["curves"]["g2"]) and res["T_grid"] == [4.2]
          and res["run_id"].startswith("FS-composed-") and res["platform"] == "composed"
          and res.get("series") == J(direct["series"]) and isinstance(res["labels"], list))
    r = post("/api/compose/run", {"spec": RATE_SPEC, "mode": "point"})
    resp = j(r)
    dr = engine.evaluate_spec(RATE_SPEC)
    got = resp.get("result")
    if got is None and "job_id" in resp:
        got = wait_job(resp["job_id"]).get("result")
    check("4c rate-solver spec: g2_0 == direct; inline (result in the response) iff the engine estimate < 1 s",
          got is not None and got["scalars"]["g2_0"] == J(dr["scalars"]["g2_0"])
          and (("result" in resp) == (float(engine.estimate(RATE_SPEC, 1)) < jobs.COMPOSED_INLINE_MAX_ETA_S)),
          f"inline={'result' in resp}")
    grid = {"lo": 4.0, "hi": 20.0, "n": 3}
    r = post("/api/compose/run", {"spec": SPEC, "mode": "tsweep", "T_grid": grid})
    resp = j(r)
    got = resp.get("result") or wait_job(resp["job_id"]).get("result")
    dts = engine.evaluate_spec(SPEC, T_grid=jobs.resolve_T_grid(grid))
    check("4d tsweep n=3: g2 and hom_I curves == direct", got is not None and got["curves"]["g2"] == J(dts["curves"]["g2"])
          and got["curves"]["hom_I"] == J(dts["curves"]["hom_I"]) and got["curves"]["T_hs"] == [4.0, 12.0, 20.0])
    check("4e tsweep n=32 -> 400", post("/api/compose/run", {"spec": SPEC, "mode": "tsweep",
                                                              "T_grid": {"lo": 4, "hi": 40, "n": 32}}).status_code == 400)
    r = post("/api/compose/run", {"spec": bad})
    check("4f run with an invalid spec -> 400 carrying the validate issues",
          r.status_code == 400 and j(r).get("issues") == J(vb["issues"]))
    check("4g run with an unknown mode -> 400", post("/api/compose/run", {"spec": SPEC, "mode": "envelope"}).status_code == 400)

    # ------------------------------------------------------------ 5. series job + cache
    values = [1.0, 2.0, 5.0]
    obs = ["g2_0", "tau_ps"]
    body = {"spec": SPEC, "path": "structure.F_P", "values": values, "observables": obs}
    r = post("/api/compose/series", body)
    resp = j(r)
    check("5 series: 202 with a job id and a positive eta", r.status_code == 202 and "job_id" in resp
          and resp["cached"] is False and resp["eta_s"] > 0, f"eta {resp.get('eta_s')}")
    snap = wait_job(resp["job_id"])
    dser = engine.series(SPEC, "structure.F_P", values, obs)
    got = (snap.get("result") or {}).get("series") or {}
    check("5b series job done; x, y arrays, tags == direct engine.series", snap["state"] == "done"
          and got.get("x") == J(dser["x"]) and got.get("y") == J(dser["y"]) and got.get("tag") == dser["tag"]
          and got.get("failed") == [], snap["state"])
    # SSE: the finished job streams its done event
    sse = c.get(f"/api/jobs/{resp['job_id']}/events").get_data(as_text=True)
    check("5c SSE stream of the job ends with a done event", "event: done" in sse)
    r2 = post("/api/compose/series", body)
    r2j = j(r2)
    check("5d second identical series call is a cache hit with eta 0 and the same result",
          r2.status_code == 200 and r2j["cached"] is True and r2j["eta_s"] == 0.0 and r2j["result"] == snap["result"])
    check("5e series bad path / empty observables / unknown observable -> 400",
          post("/api/compose/series", {**body, "path": "structure.nope"}).status_code == 400
          and post("/api/compose/series", {**body, "observables": []}).status_code == 400
          and post("/api/compose/series", {**body, "observables": ["not_an_observable"]}).status_code == 400
          and post("/api/compose/series", {**body, "values": []}).status_code == 400)
    check("5f series_job is on the exact job allowlist and an arbitrary function is refused",
          "fsim_studio.api_compose:series_job" in jobs.ALLOWED_CALLS
          and "fsim_studio.api_compose:reproduce_job" in jobs.ALLOWED_CALLS
          and post("/api/compose/series", body).status_code == 200)

    # ------------------------------------------------------------ 6. /api/run composed
    design = {"name": "composed-probe", "platform": "composed", "composed": SPEC}
    r = post("/api/run", {"design": design, "mode": "point"})
    resp = j(r)
    res = resp.get("result") or wait_job(resp["job_id"]).get("result")
    check("6 /api/run with a composed design: g2_0 == direct, T_grid [T_K]",
          r.status_code == 200 and res is not None and res["scalars"]["g2_0"] == J(direct["scalars"]["g2_0"])
          and res["T_grid"] == [4.2] and res["platform"] == "composed")
    r = post("/api/run", {"design": design, "mode": "envelope"})
    check("6b /api/run composed + envelope -> 400 'envelope not supported for composed designs in v1'",
          r.status_code == 400 and j(r)["error"] == "envelope not supported for composed designs in v1"
          and j(r).get("kind") == "unsupported")
    check("6c /api/run composed + headline -> 400", post("/api/run", {"design": design, "mode": "headline"}).status_code == 400)
    badd = dict(design, composed=bad)
    check("6d /api/run composed with an unknown kind -> 400", post("/api/run", {"design": badd}).status_code == 400)
    check("6e platform_class(composed) == 'composed'; legacy classes unchanged",
          jobs.platform_class(__import__("fsim_studio.api_cards", fromlist=["x"]).design_from_dict(design)) == "composed")

    # ------------------------------------------------------------ 7. composed card round trip
    r = c.put("/api/cards/composed-probe-card", json={"design": design}, headers=HDR)
    saved = r.status_code == 200
    g = c.get("/api/cards/composed-probe-card") if saved else None
    check("7 PUT /api/cards/<slug> saves a composed design and GET returns it unchanged",
          saved and g.status_code == 200 and j(g)["design"]["composed"] == SPEC and j(g)["platform"] == "composed",
          f"PUT {r.status_code}: {r.get_data(as_text=True)[:200] if not saved else ''}")
    if saved:
        cards = j(c.get("/api/cards"))
        row = next((x for x in cards if x["name"] == "composed-probe-card"), None)
        check("7b /api/cards lists it with platform 'composed'", row is not None and row.get("platform") == "composed",
              str(row))
        # 7c. R2 finish (FIN1): the last run is stored beside the saved card, keyed by the design's spec hash
        fake_run = {"mode": "point", "scalars": {"g2_0": 0.123, "valid": True}, "run_id": "probe"}
        g0 = c.get("/api/compose/saved-run/composed-probe-card")
        check("7c saved-run GET before any run: run null, matches false",
              g0.status_code == 200 and j(g0)["run"] is None and j(g0)["matches"] is False, str(g0.status_code))
        rp = c.put("/api/compose/saved-run/composed-probe-card", json={"spec": SPEC, "run": fake_run}, headers=HDR)
        check("7d saved-run PUT stores the run and returns a 12-hex spec hash",
              rp.status_code == 200 and len(j(rp)["spec_hash"]) == 12, f"{rp.status_code} {rp.get_data(as_text=True)[:120]}")
        g1 = c.get("/api/compose/saved-run/composed-probe-card")
        check("7e saved-run GET returns the stored run, the same hash, and matches the card's design",
              g1.status_code == 200 and j(g1)["run"] == fake_run and j(g1)["spec_hash"] == j(rp)["spec_hash"] and j(g1)["matches"] is True,
              g1.get_data(as_text=True)[:200])
        c.put("/api/cards/composed-probe-card", json={"design": dict(design, composed=dict(SPEC, T_K=7.0))}, headers=HDR)
        g2 = c.get("/api/compose/saved-run/composed-probe-card")
        check("7f after the design changes, the stored run no longer matches its hash",
              g2.status_code == 200 and j(g2)["matches"] is False and j(g2)["spec_hash"] == j(rp)["spec_hash"] != j(g2)["current_hash"],
              g2.get_data(as_text=True)[:200])
        check("7g saved-run: GET on an unsaved card -> 200 run null; PUT -> 404; a bad name -> 400, a body without a run -> 400, no header -> 403",
              c.get("/api/compose/saved-run/no-such-card").status_code == 200 and j(c.get("/api/compose/saved-run/no-such-card"))["run"] is None
              and c.put("/api/compose/saved-run/no-such-card", json={"spec": SPEC, "run": fake_run}, headers=HDR).status_code == 404
              and c.put("/api/compose/saved-run/Bad_Name", json={"spec": SPEC, "run": fake_run}, headers=HDR).status_code == 400
              and c.put("/api/compose/saved-run/composed-probe-card", json={"spec": SPEC}, headers=HDR).status_code == 400
              and c.put("/api/compose/saved-run/composed-probe-card", json={"spec": SPEC, "run": fake_run}).status_code == 403)
        # 7h. an incompatible pair is not a device over HTTP: validate -> ok:false (an error), run -> 400
        inc_spec = dict(SPEC, material={"kind": "cdse_znse", "params": {}}, structure={"kind": "edge_ridge", "params": {}})
        vi = post("/api/compose/validate", {"spec": inc_spec})
        check("7h incompatible pair (edge_ridge + cdse_znse): /validate ok:false with an incompatible error, eta_s null",
              vi.status_code == 200 and j(vi)["ok"] is False and j(vi)["eta_s"] is None
              and any(i["level"] == "error" and "incompatible" in i["message"] for i in j(vi)["issues"]), vi.get_data(as_text=True)[:200])
        ri = post("/api/compose/run", {"spec": inc_spec, "mode": "point"})
        check("7i incompatible pair: /run is refused with a 4xx, not a crash", 400 <= ri.status_code < 500, str(ri.status_code))

    # ------------------------------------------------------------ 8. scene
    if have("fsim_studio.scene_compose"):
        from fsim_studio import scene_compose
        want = J(scene_compose.describe(SPEC))
        r = c.get("/api/compose/scene", query_string={"spec": json.dumps(SPEC)})
        r2 = post("/api/compose/scene", {"spec": SPEC})
        check("8 GET and POST /api/compose/scene == scene_compose.describe(spec)",
              r.status_code == 200 and j(r) == want and r2.status_code == 200 and j(r2) == want)
        check("8b scene without spec -> 400", c.get("/api/compose/scene").status_code == 400)
    else:
        skip("8 /api/compose/scene", "fsim_studio.scene_compose not present (P13)")

    # ------------------------------------------------------------ 9. papers
    if not papers_present:
        skip("9 /api/papers routes", "fsim_core.compose.papers not present (P7)")
        r = c.get("/api/papers")
        check("9a /api/papers while papers.py is absent -> 501 (not a crash)", r.status_code == 501)
        return
    from fsim_core.compose import papers
    pdir = Path(tmp.name) / "papers"
    pdir.mkdir()
    real = [p for p in sorted((ROOT / "cards" / "papers").glob("*.yaml")) if not p.name.startswith("_")]
    src = real[0] if real else None
    if src is None:
        fx = ROOT / "cards" / "papers" / "_fixture.yaml"
        if fx.is_file():
            import yaml
            doc = yaml.safe_load(fx.read_text(encoding="utf-8"))
            doc["meta"]["id"] = "fixturepaper"
            src = pdir / "fixturepaper.yaml"
            src.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    if src is None:
        skip("9 /api/papers routes", "no paper card (cards/papers/*.yaml) and no _fixture.yaml to copy")
        return
    if src.parent != pdir:
        shutil.copy(src, pdir / src.name)
    api_compose.PAPERS_DIR = pdir
    pid = src.stem
    try:
        paper = papers.load_paper(pdir / f"{pid}.yaml")
    except Exception as exc:  # noqa: BLE001
        check("9 paper card loads", False, f"{type(exc).__name__}: {exc}")
        return
    lst = j(c.get("/api/papers"))
    row = next((x for x in lst if x["id"] == pid), None)
    check("9 GET /api/papers lists the card with cite, read_level, kinds, class counts, n items",
          row is not None and row["cite"] == paper.meta["cite"] and row["read_level"] == paper.meta["read_level"]
          and row["kinds"]["material"] == paper.device["material"]["kind"]
          and row["n_items"] == len(paper.reproduce) and set(row["classes"]) == {"C", "M", "T", "P", "R"}
          and row["last_outcome"] is None)
    g = j(c.get(f"/api/papers/{pid}"))
    check("9b GET /api/papers/<id> == the loaded paper (device, reported, series, reproduce, inputs_from_paper)",
          g["device"] == J(paper.device) and g["reported"] == J(paper.reported) and g["series"] == J(paper.series)
          and g["reproduce"] == J(paper.reproduce) and g["inputs_from_paper"] == sorted(paper.inputs_from_paper)
          and g["resolved"] == J(validate_spec(paper.device)["resolved"]))
    check("9c unknown paper -> 404; traversal id -> 400/404",
          c.get("/api/papers/no_such_paper").status_code == 404 and c.get("/api/papers/_fixture").status_code in (400, 404)
          and c.get("/api/papers/..%2Fsecret").status_code in (400, 404))
    d = c.get(f"/api/papers/{pid}/design")
    dj = j(d)
    check("9d GET /api/papers/<id>/design is a valid composed design with no ranged leaf left",
          d.status_code == 200 and dj["design"]["platform"] == "composed"
          and validate_spec(dj["design"]["composed"])["ok"] is True
          and '"range":' not in json.dumps(dj["design"]["composed"]))
    # 9d2 (R2 finding 8): Open in Composer opens the SAME device Papers reproduced. For every card the collapsed design
    # evaluates to the central device (exact); a log-scale ParamSpec gives the geometric mean (params.central), not the mean.
    import yaml as _yaml
    from fsim_studio import api_compose as _ac
    bad_cards, n_ranged = [], 0
    for cf in real:
        dev = _yaml.safe_load(cf.read_text(encoding="utf-8")).get("device") or {}
        cd = _ac.central_design(cf.stem, dev)
        if not cd["collapsed"]:
            continue
        n_ranged += 1
        try:
            a = engine.evaluate_spec(dict(cd["design"]["composed"], schema="fsim-composed/1"))["scalars"]
            b = engine.evaluate_spec(dict(dev, schema="fsim-composed/1"))["scalars"]
        except Exception as exc:  # noqa: BLE001
            bad_cards.append(f"{cf.stem}: {type(exc).__name__}")
            continue
        for k in ("g2_0", "hom_I", "B_first_lens", "tau_ps"):
            if k in a or k in b:
                if a.get(k) != b.get(k) and not (a.get(k) != a.get(k) and b.get(k) != b.get(k)):
                    bad_cards.append(f"{cf.stem}:{k}")
    check(f"9d2 central_design == evaluate_spec(device) for every card with a ranged leaf ({n_ranged} cards)", not bad_cards and n_ranged > 0, str(bad_cards[:4]))
    if (ROOT / "cards" / "papers" / "somaschi2016_qd1.yaml").is_file():
        dev = _yaml.safe_load((ROOT / "cards" / "papers" / "somaschi2016_qd1.yaml").read_text(encoding="utf-8"))["device"]
        cd = _ac.central_design("somaschi2016_qd1", dev)["design"]["composed"]
        check("9d3 somaschi2016_qd1 hom_I is equal in Papers and Composer (log-scale tau_cap_ps centre)",
              engine.evaluate_spec(dict(cd, schema="fsim-composed/1"))["scalars"]["hom_I"]
              == engine.evaluate_spec(dict(dev, schema="fsim-composed/1"))["scalars"]["hom_I"])
    # 5e (R2 finding 5): observable requests are validated with their convention up front
    bad = [("hom_V", True), ("g2_0:nonsense", True), ("hom_V:raw", False), ("hom_V:corrected", False), ("g2_0", False)]
    got = []
    for req, must_fail in bad:
        try:
            _ac._check_observables([req])
            failed = False
        except _ac.ComposeError:
            failed = True
        got.append(failed == must_fail)
    check("5e _check_observables rejects a convention outside the observable's list ('hom_V', 'g2_0:nonsense') and accepts valid ones", all(got), str(got))
    rr = post("/api/compose/series", {"spec": SPEC, "path": "structure.F_P", "values": [1.0, 2.0], "observables": ["hom_V"]})
    check("5f series with a bare hom_V is a 400 at submit time with a clear message, not a failed job",
          rr.status_code == 400 and "convention" in json.dumps(j(rr)), str(rr.status_code))
    # R2 F56: the payload carries the cw-variant advertised set; an unavailable series request comes back with its reason
    pay56 = {d["kind"]: d for d in j(c.get("/api/compose/registry"))["drives"]}
    check("5g registry payload: above_band / p_shell carry observables_cw without the per-pulse ids, resonant_pulse does not carry it",
          all("observables_cw" in pay56[k] and not ({"hom_I", "hom_V", "B_first_lens", "B_fibre", "hom_vs_dt"} & set(pay56[k]["observables_cw"]))
              and "hom_I" in pay56[k]["observables"] for k in ("above_band", "p_shell"))
          and "observables_cw" not in pay56["resonant_pulse"])
    ep_spec = dict(SPEC, drive={"kind": "electrical_pulsed", "params": {}}, structure={"kind": "bulk_planar", "params": {}})
    rr = post("/api/compose/series", {"spec": ep_spec, "path": "T_K", "values": [4.0, 6.0], "observables": ["hom_I_gated", "g2_0"]})
    snap56 = wait_job(j(rr)["job_id"]) if rr.status_code == 202 else {"state": str(rr.status_code)}
    un56 = ((snap56.get("result") or {}).get("series") or {}).get("unavailable") or {}
    check("5h series of hom_I_gated with no HOM gate: done, 'unavailable' carries the reason, g2_0 is not listed",
          snap56["state"] == "done" and "no HOM gate set" in un56.get("hom_I_gated", "") and "g2_0" not in un56, json.dumps(un56)[:160])
    r = post(f"/api/papers/{pid}/reproduce", {})
    resp = j(r)
    check("9e reproduce: 202 job", r.status_code == 202 and "job_id" in resp and resp["cached"] is False)
    snap = wait_job(resp["job_id"], timeout=600)
    direct_rep = papers.reproduce(papers.load_paper(pdir / f"{pid}.yaml"))
    import dataclasses
    drep = dataclasses.asdict(direct_rep) if dataclasses.is_dataclass(direct_rep) else direct_rep
    check("9f reproduce job done; report items == direct papers.reproduce", snap["state"] == "done"
          and (snap.get("result") or {}).get("report") == J(drep), snap["state"])
    r2 = post(f"/api/papers/{pid}/reproduce", {})
    check("9g second reproduce is a cache hit with eta 0",
          r2.status_code == 200 and j(r2)["cached"] is True and j(r2)["eta_s"] == 0.0)
    row = next((x for x in j(c.get("/api/papers")) if x["id"] == pid), None)
    check("9h after reproduce, /api/papers carries last_outcome from the cache", row is not None and row["last_outcome"]
          is not None and row["last_outcome"]["cached"] is True, str(row and row["last_outcome"]))
    # the cache key follows the card file bytes
    k1 = api_compose._reproduce_key(pid, pdir / f"{pid}.yaml")
    with open(pdir / f"{pid}.yaml", "a", encoding="utf-8") as fh:
        fh.write("\n# edited\n")
    check("9i reproduce cache key changes when the card file changes",
          api_compose._reproduce_key(pid, pdir / f"{pid}.yaml") != k1)

    # ---------------------------------------------------------------- 10. examples + Fock tier (devlib-p11b)
    from fsim_core.device import DeviceDesign
    from fsim_core.compose import tiers as tiers_core
    edir = api_compose.EXAMPLES_DIR
    files = sorted(p.stem for p in edir.glob("*.yaml"))
    ex = j(c.get("/api/compose/examples"))
    check("10a /api/compose/examples == cards/composed/*.yaml (>= 6 files)",
          [e["name"] for e in ex] == files and len(files) >= 6, ",".join(files))
    check("10b every example has label, summary, four kinds and a composed design with schema",
          all(e["label"] and e["summary"] and set(e["kinds"]) == {"material", "structure", "drive", "setup"}
              and e["design"]["platform"] == "composed" and e["design"]["composed"]["schema"] == "fsim-composed/1"
              for e in ex))
    check("10c every example loads via DeviceDesign.load and validates",
          all(J(DeviceDesign.load(edir / f"{e['name']}.yaml").composed) == e["design"]["composed"]
              and validate_spec(e["design"]["composed"])["ok"] is True for e in ex))

    def leaves(node):
        if isinstance(node, dict):
            if "tag" in node and "value" in node:
                yield node
            else:
                for v in node.values():
                    yield from leaves(v)
    check("10d every example leaf carries a tag and a source; summaries say class example",
          all(lf.get("source") and lf["tag"] in ("V", "DR", "E", "A")
              for e in ex for lf in leaves(e["design"]["composed"])) and all("class example" in e["summary"].lower() for e in ex))
    for e in ex:
        r = post("/api/compose/run", {"spec": e["design"]["composed"], "mode": "point", "name": e["name"]})
        snap = wait_job(j(r)["job_id"], timeout=600) if r.status_code in (200, 202) else {"state": f"http {r.status_code}"}
        check(f"10e example {e['name']} runs (point) to done", snap["state"] == "done", snap["state"])
    fk = c.get("/api/compose/tier/fock")
    fj = j(fk)
    check("10f /api/compose/tier/fock returns a spec validate_spec accepts",
          fk.status_code == 200 and fj["design"]["platform"] == "composed"
          and validate_spec(fj["design"]["composed"])["ok"] is True, str(fk.status_code))
    check("10g that spec == tiers.fock_card_to_spec()", fj["design"]["composed"] == J(tiers_core.fock_card_to_spec()))


if __name__ == "__main__":
    sys.exit(main())
