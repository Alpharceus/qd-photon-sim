"""FSIM Studio backend API checks (Flask test client, no network).

Every number the API returns is compared against a direct fsim_core call or
the committed out/ file it claims to come from. The cache is redirected to a
temporary directory (FSIM_STUDIO_CACHE) so runs start cold and nothing is
written into the repository.

Run: python verify/verify_studio_api.py   (exit code 0 iff all pass)
"""
import concurrent.futures as cf
import copy
import json
import math
import os
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def wait_job(client, job_id, timeout=300.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        snap = client.get(f"/api/jobs/{job_id}").get_json()
        if snap["state"] in ("done", "error", "cancelled"):
            return snap
        time.sleep(0.05)
    raise TimeoutError(job_id)


def sse_events(client, job_id):
    resp = client.get(f"/api/jobs/{job_id}/events", buffered=False)
    text = "".join(c.decode() if isinstance(c, bytes) else c for c in resp.response)
    events = []
    for block in text.split("\n\n"):
        lines = [l for l in block.splitlines() if l and not l.startswith(":")]
        if not lines:
            continue
        etype = next((l[7:] for l in lines if l.startswith("event: ")), None)
        data = next((l[6:] for l in lines if l.startswith("data: ")), None)
        events.append((etype, json.loads(data) if data else None))
    return resp, events


def polish_checks(c, app, jm, tmp, direct_h):
    """studio-07-polish: H1 H3 H4 H5 H10-H14 H16 H17 H21 and interfaces I1-I6.
    Every check here fails on the pre-polish backend."""
    import re
    import shutil

    import fsim_studio
    from fsim_core.device import DeviceDesign
    from fsim_studio import cache, jobs, scene, stories
    from fsim_studio.api_cards import design_to_dict
    from fsim_studio.server import create_app
    blocks = []

    def _b0():  # I4 / H1 lattice: only headline-model rows gate; the rest is reported apart
        lat = c.get("/api/scene/lattice?campaign=rt_edge&headline=0").get_json()
        o = lat["overlays"]
        b = o["banner"]
        rows = o["rows"]
        gates = rows.get("gates") or []
        hl = [fp and not tc for fp, tc in zip(rows.get("model_finite_pulse") or [],
                                              rows.get("model_tau_cap_density") or [])]
        check("I4 lattice rt_edge headline=0: banner counts gating rows only (0/768), never pass-style",
              isinstance(b, dict) and b["pass"] == 0 and b["gating_rows"] == 768 and b["style"] == "fail"
              and b["rows"] == 3072 and o["counts"]["pass"] == 0, (b.get("text") if isinstance(b, dict) else b))
        ng = (b or {}).get("non_gating_passes") if isinstance(b, dict) else None
        check("I4 lattice: non-gating passes (232) reported separately with 'never gates'",
              ng is not None and ng["count"] == 232 and ng["rows"] == 2304 and "never gates" in ng["text"])
        check("I4 lattice rows carry model_finite_pulse / model_tau_cap_density / gates (true only headline)",
              len(gates) == 3072 and gates == hl and sum(gates) == 768
              and all(not rows["device_pass"][i] for i in range(3072) if gates[i]))
        lat_n = scene.lattice_scene("nitride_cavity")
        check("I4 lattice without model axes: every row gates, model_* null",
              all(lat_n["overlays"]["rows"].get("gates") or [False])
              and set(lat_n["overlays"]["rows"].get("model_finite_pulse") or [0]) == {None}
              and isinstance(lat_n["overlays"]["banner"], dict)
              and lat_n["overlays"]["banner"]["non_gating_passes"] is None)
    blocks.append(('I4 / H1 lattice: only headline-model rows gate; the rest is reported apart', _b0))

    def _b1():  # I1 / H3 tag_chain for nitride / nanowire results
        r = c.post("/api/run", json={"card": "nitride-cavity-set-design", "mode": "point",
                                     "T_grid": [300.0]}).get_json()
        snap = wait_job(c, r["job_id"]) if not r.get("result") else {"state": "done", "result": r["result"]}
        nres = snap.get("result") or {}
        check("I1 nitride run: tag_chain parsed from provenance (widest of [A/E/DR] -> A)",
              snap["state"] == "done" and nres.get("tag_chain") == "A"
              and nres.get("tag_chain_source") == "parsed from provenance"
              and nres["scalars"].get("tag_chain") is None, f"{nres.get('tag_chain')} / {nres.get('tag_chain_source')}")
        try:
            tags = __import__("importlib").import_module("fsim_studio.tags")
            helper_ok = (tags.chain_of({"provenance": {"a": "[V abstract-only] x", "b": {"c": "[DR, re-fit]"}}})
                         == ("DR", "parsed from provenance")
                         and tags.chain_of({}) == ("A", "no tag returned; assumed widest")
                         and tags.chain_of({"tag_chain": "[E]"}) == ("E", "fsim_core tag_chain"))
        except ImportError:
            helper_ok = False
        check("I1 tag helpers: [V abstract-only]/[DR, re-fit] parse; nothing -> A assumed", helper_ok)
        FLOOR = "structural floor set by b_res [A]"

        # ...... I2 labels: b_res floor on SET; non-headline model on edge re-runs (H17/H21)
        check("I2 SET g2 at the b_res floor carries 'structural floor set by b_res [A]'",
              getattr(jobs, "FLOOR_LABEL", None) == FLOOR and FLOOR in nres.get("labels", []))
        r = c.post("/api/run", json={"card": "nitride-cavity-pulse-design", "mode": "point",
                                     "T_grid": [300.0]}).get_json()
        snap = wait_job(c, r["job_id"]) if not r.get("result") else {"state": "done", "result": r["result"]}
        check("I2 pulse card (g2 above the floor): no floor label",
              snap["state"] == "done" and FLOOR not in snap["result"]["labels"]
              and snap["result"]["tag_chain"] == "A")
        d_edge = DeviceDesign.load(fsim_studio.CARDS / "edge-inp-gainp-design.yaml")
        d_ft = copy.deepcopy(d_edge)
        d_ft.drive.finite_pulse = True
        d_ft.ret.tau_cap_scales_with_density = True
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "design": design_to_dict(d_ft),
                                     "mode": "point", "T_grid": [230.0]}).get_json()
        snap = wait_job(c, r["job_id"])
        lab_ft = (snap.get("result") or {}).get("labels", [])
        check("I2/H17 re-run of a finite_pulse+tau_cap_density edge row: 'non-headline model', not static",
              snap["state"] == "done" and "non-headline model" in lab_ft and jobs.STATIC_LABEL not in lab_ft, str(lab_ft))
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "mode": "point", "T_grid": [300.0]}).get_json()
        snap = wait_job(c, r["job_id"])
        lab_st = (snap.get("result") or {}).get("labels", [])
        d_h = copy.deepcopy(d_edge)
        d_h.drive.finite_pulse = True
        d_h.ret.tau_cap_scales_with_density = False
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "design": design_to_dict(d_h),
                                     "mode": "point", "T_grid": [230.0]}).get_json()
        snap_h = wait_job(c, r["job_id"])
        check("I2/H21 static edge card point run: static AND non-headline model; headline switches: neither",
              jobs.STATIC_LABEL in lab_st and "non-headline model" in lab_st
              and "non-headline model" not in snap_h["result"]["labels"]
              and jobs.STATIC_LABEL not in snap_h["result"]["labels"], str(lab_st))
    blocks.append(('I1 / H3 tag_chain for nitride / nanowire results', _b1))

    def _b2():  # I3 gates
        g = c.get("/api/gates")
        gj = g.get_json() or {}
        check("I3 GET /api/gates: g2 ceiling 0.5 and flux floor 1000 photons/s [A], parsed with sources",
              g.status_code == 200 and gj["g2_ceiling"]["value"] == 0.5
              and gj["g2_ceiling"]["source"].startswith("out/rt_edge/verdict.md:")
              and gj["g2_ceiling"].get("contract_source", "").startswith("docs/rt_edge_contract.md:")
              and gj["flux_floor"]["value"] == 1000.0 and gj["flux_floor"]["tag"] == "A"
              and gj["flux_floor"]["unit"] == "photons/s", str(gj)[:160])
    blocks.append(('I3 gates', _b2))

    def _b3():  # I6 / H2 card -> verdict line(s)
        v1 = c.get("/api/cards/nitride-cavity-set-design/verdict").get_json() or {}
        check("I6 nitride-cavity-set -> nitride_cavity deterministic_pair line (exact; g2_min 0.17355...)",
              v1.get("match") == "exact" and v1["campaign"] == "nitride_cavity" and len(v1["lines"]) == 1
              and v1["lines"][0]["fields"]["regime"] == "deterministic_pair"
              and v1["lines"][0]["fields"]["g2_min"].startswith("0.17355"), str(v1.get("reason"))[:100])
        v2 = c.get("/api/cards/nitride-nonpolar-set-design/verdict").get_json() or {}
        check("I6 nonpolar SET -> a_plane deterministic_pair lines, ambiguous (2), not picked",
              v2.get("match") == "ambiguous" and v2["campaign"] == "nitride_geometry_stark"
              and len(v2["lines"]) == 2 and "ambiguous: 2 lines match" in v2["reason"]
              and all(l["fields"]["family"] == "a_plane" and l["fields"]["regime"] == "deterministic_pair"
                      and l["word"] == "pass_hardware_infeasible" for l in v2["lines"]))
        v3 = c.get("/api/cards/nitride-nanowire-vertical-set-design/verdict").get_json() or {}
        v4 = c.get("/api/cards/staged-device-design/verdict").get_json() or {}
        v5 = c.get("/api/cards/no-such-card/verdict")
        check("I6 nanowire exact on family/regime/strain/rep; legacy -> none; unknown card -> 404",
              v3.get("match") == "exact" and v3["lines"][0]["fields"]["rep_rate_hz"] == "2e+08"
              and v3["lines"][0]["fields"]["strain_bound"] == "relaxed" and v4.get("match") == "none"
              and v5.status_code == 404)
    blocks.append(('I6 / H2 card -> verdict line(s)', _b3))

    def _b4():  # H10 scene card resolution only through api_cards
        s1 = c.get("/api/scene/device?card=C:/Windows/win.ini")
        s2 = c.get("/api/scene/device?card=../cards/staged-device-design")
        s3 = c.get("/api/scene/device?card=no-such-card")
        s4 = c.get("/api/scene/band?card=chatzarakis")
        check("H10 scene card: paths refused 400 (invalid card name), unknown 404, param card 422",
              s1.status_code == 400 and "invalid card name" in s1.get_json()["error"]
              and s2.status_code == 400 and s3.status_code == 404 and s4.status_code == 422,
              f"{s1.status_code},{s2.status_code},{s3.status_code},{s4.status_code}")
    blocks.append(('H10 scene card resolution only through api_cards', _b4))

    def _b5():  # H11 set_state is a no-op after cancel / terminal
        j = jobs.Job("point", 1.0, "k" * 64, "FS-x-000000")
        jm._register(j)
        jm.cancel(j.id)
        revived = j.set_state("running")
        j2 = jobs.Job("point", 1.0, "k" * 64, "FS-x-000001")
        j2.finish({"ok": 1})
        j2.set_state("running")
        check("H11 set_state after cancel / done does not revive the job",
              j.state == "cancelled" and revived is False and j2.state == "done")
    blocks.append(('H11 set_state is a no-op after cancel / terminal', _b5))

    def _b6():  # H12 filter conflicting with the headline pin
        bad = c.get('/api/campaigns/rt_edge/sweep?headline=1&filter={"model_tau_cap_density":true}')
        ok = c.get('/api/campaigns/rt_edge/sweep?headline=1&filter={"model_finite_pulse":true}&limit=1')
        check("H12 headline=1 + conflicting filter key -> 400; agreeing key -> 200",
              bad.status_code == 400 and "headline" in bad.get_json()["error"] and ok.status_code == 200
              and ok.get_json()["filtered"] == 768)
    blocks.append(('H12 filter conflicting with the headline pin', _b6))

    def _b7():  # H13 Host allowlist, Studio header on mutating requests, bake + T_grid limits
        bare = app.test_client()
        h_evil = bare.get("/api/health", headers={"Host": "evil.example:8765"})
        p_nohdr = bare.post("/api/validate", json={"design": {}})
        check("H13 foreign Host -> 403; POST without X-FSIM-Studio -> 403",
              h_evil.status_code == 403 and p_nohdr.status_code == 403
              and "X-FSIM-Studio" in p_nohdr.get_json()["error"])
        tg = c.post("/api/run", json={"card": "staged-device-design", "mode": "point",
                                      "T_grid": {"lo": 4, "hi": 300, "n": 401}})
        check("H13 T_grid.n > 400 -> 400", tg.status_code == 400 and "400" in tg.get_json()["error"])
        errs = []
        for src in ({"kind": "md_regex", "file": ".git/HEAD", "pattern": "(ref)"},
                    {"kind": "md_regex", "file": "fsim_core/device.py", "pattern": "(def)"},
                    {"kind": "md_regex", "file": "out/rt_edge/verdict.md", "pattern": "((a+)+)$"},
                    {"kind": "md_regex", "file": "out/rt_edge/verdict.md", "pattern": "(" + "x" * 400 + ")"},
                    {"kind": "json", "file": "../outside.json", "key": "a"}):
            try:
                stories.RESOLVERS[src["kind"]](src)
                errs.append(None)
            except Exception as exc:  # noqa: BLE001 -- any refusal counts; a returned value does not
                errs.append(str(exc))
        check("H13 bake sources: only under out/; nested-quantifier / over-long patterns rejected",
              all(e is not None for e in errs) and "under out/" in errs[0] and "under out/" in errs[1]
              and "nested" in errs[2] and "longer" in errs[3], str(errs)[:160])
        app_p = create_app(job_manager=jm, port=8765)
        cp = app_p.test_client()
        check("H13 with a port: only 127.0.0.1:<port> / localhost:<port> accepted",
              cp.get("/api/health", headers={"Host": "127.0.0.1:8765"}).status_code == 200
              and cp.get("/api/health", headers={"Host": "localhost:8765"}).status_code == 200
              and cp.get("/api/health", headers={"Host": "127.0.0.1:9999"}).status_code == 403
              and cp.get("/api/health", headers={"Host": "localhost"}).status_code == 403)
    blocks.append(('H13 Host allowlist, Studio header on mutating requests, bake + T_grid limits', _b7))

    def _b8():  # H14 fsim_core hash follows edits; studio schema version in the cache key
        fake = Path(tmp.name) / "fakeroot"
        (fake / "fsim_core").mkdir(parents=True)
        f = fake / "fsim_core" / "m.py"
        f.write_text("x = 1\n", encoding="utf-8")
        orig_root = fsim_studio.ROOT
        fsim_studio.ROOT = fake
        try:
            h1 = fsim_studio.fsim_core_hash()
            f.write_text("x = 22\n", encoding="utf-8")
            st = f.stat()
            os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 5_000_000_000))
            h2 = fsim_studio.fsim_core_hash()
        finally:
            fsim_studio.ROOT = orig_root
        check("H14 fsim_core hash recomputed when a source changes (mtime/size signature)",
              h1 != h2 and fsim_studio.fsim_core_hash() == fsim_studio.fsim_core_hash())
        payload = {"x": 1}
        k1 = cache.make_key(payload)
        orig_schema = getattr(cache, "STUDIO_SCHEMA", None)
        cache.STUDIO_SCHEMA = "other-schema"
        try:
            k2 = cache.make_key(payload)
        finally:
            cache.STUDIO_SCHEMA = orig_schema
        check("H14 cache key includes the studio result-schema version", orig_schema is not None and k1 != k2)
    blocks.append(('H14 fsim_core hash follows edits; studio schema version in the cache key', _b8))

    def _b9():  # H16 exact job allowlist
        blocked = []
        for fn in ("fsim_studio.scene:lattice_scene", "fsim_studio.scene:_load", "os:system"):
            try:
                jm.submit_call(fn, {}, 0.0)
                blocked.append(False)
            except ValueError:
                blocked.append(True)
        try:
            jobs._w_call("fsim_studio.scene:device_scene", {"card": "staged-device-design"})
            w_blocked = False
        except ValueError:
            w_blocked = True
        # studio-p2d: the explorer modules add EXACT entries (no wildcard, no prefix match at run time);
        # the check pins the scene entry, the six api_explore_b entries, and that every other entry is a
        # named function of an api_explore_{a,b} module.
        allowed = getattr(jobs, "ALLOWED_CALLS", frozenset())
        b_entries = {f"fsim_studio.api_explore_b:{n}" for n in ("phonon_frames", "transport_sweep", "stark_trace", "sde_frame", "sde_li", "sde_trace")}
        check("H16 job allowlist is exact: the scene entry plus named api_explore_{a,b} functions",
              isinstance(allowed, frozenset) and "fsim_studio.scene:surface_compute" in allowed and b_entries <= allowed
              and all(fn == "fsim_studio.scene:surface_compute" or re.fullmatch(r"fsim_studio\.api_explore_[ab]:[a-z_0-9]+", fn) for fn in allowed)
              and all(blocked) and w_blocked)
    blocks.append(('H16 exact job allowlist', _b9))

    def _b10():  # H4 / H5 / I5 stories (on a temp copy; the shipped file is not rewritten)
        story_tmp = Path(tmp.name) / "stories_polish"
        story_tmp.mkdir()
        shutil.copy(stories.STORIES / "rt-single-photons.json", story_tmp)
        orig_dir = stories.STORIES
        stories.STORIES = story_tmp
        try:
            st = c.get("/api/stories/rt-single-photons").get_json()
            if st["stale"]:
                c.post("/api/stories/rt-single-photons/bake")
                st = c.get("/api/stories/rt-single-photons").get_json()
            scenes = {sc["id"]: sc for sc in st["scenes"]}
            prose_ok = all("{" not in (sc.get(f) or "") for sc in st["scenes"] for f in ("title", "claim", "caveat", "notes"))
            check("I5 story GET: prose interpolated (no {name} left), templates returned, no prose errors",
                  prose_ok and st.get("prose_errors") == []
                  and (scenes["s01"].get("claim_template") or "").count("{t_points[") == 2
                  and "from 78 to 230 K" in scenes["s01"]["claim"])
            check("I5 every scene has non-empty presenter notes (13/13)",
                  len(st["scenes"]) == 13 and all((sc.get("notes") or "").strip() for sc in st["scenes"]))
            hand = re.compile(r"\b(78|77|220|230|300)\b|100 ps|230-300")
            flagged = {sid: " ".join(PLACEHOLDER_STRIP.sub("", scenes[sid].get(f + "_template") or scenes[sid].get(f) or "")
                                     for f in ("title", "claim", "caveat"))
                       for sid in ("s01", "s02", "s04", "s06", "s11", "s12")}
            check("H5 s01/s02/s04/s06/s11/s12 title/claim/caveat carry no hand-typed numbers",
                  not any(hand.search(t) for t in flagged.values()),
                  ";".join(f"{k}:{hand.search(t).group(0)}" for k, t in flagged.items() if hand.search(t)))
            nums = {n["name"]: n for sc in st["scenes"] for n in sc["numbers"]}
            check("H5 tags resolved at bake from the source (bracket tag, live tag_chain, json tag_chain)",
                  nums["flux_floor"]["tag"] == "A" and "bracket tag" in nums["flux_floor"]["tag_source"]
                  and "live run tag_chain" in nums["live_g2"]["tag_source"]
                  and "tag_chain" in nums["tc_fit"]["tag_source"] and nums["live_g2"]["at_T"] == 230.0)
            # PUT cannot store typed values: resolved fields are stripped and restored
            body = json.loads(json.dumps(st))
            tgt = next(n for n in body["scenes"][4]["numbers"] if n["name"] == "g2_min")
            tgt["value"] = 0.1
            tgt["raw"] = "0.1"
            tgt["tag"] = "V"
            body["scenes"][0]["claim"] = "typed over the template"
            pr = c.put("/api/stories/rt-single-photons", json=body)
            stp = c.get("/api/stories/rt-single-photons").get_json()
            n2 = next(n for n in stp["scenes"][4]["numbers"] if n["name"] == "g2_min")
            check("H4 PUT strips typed values / tags / rendered prose; story stays fresh",
                  pr.status_code == 200 and n2["value"] == 0.9817 and n2["raw"] == "0.9817" and n2["tag"] != "V"
                  and stp["stale"] is False and "{t_points[0]}" in stp["scenes"][0]["claim_template"])
            # a value typed straight into the file marks the story stale
            raw_doc = json.loads((story_tmp / "rt-single-photons.json").read_text(encoding="utf-8"))
            next(n for n in raw_doc["scenes"][4]["numbers"] if n["name"] == "g2_min")["value"] = 0.1
            (story_tmp / "rt-single-photons.json").write_text(json.dumps(raw_doc), encoding="utf-8")
            check("H4 a value typed into the stored story -> stale",
                  c.get("/api/stories/rt-single-photons").get_json()["stale"] is True)
            raw_doc_ok = stories.load_story("rt-single-photons")
            next(n for n in raw_doc_ok["scenes"][4]["numbers"] if n["name"] == "g2_min")["value"] = 0.9817
            (story_tmp / "rt-single-photons.json").write_text(json.dumps(raw_doc_ok), encoding="utf-8")
            files = stories._referenced_files(raw_doc_ok)
            orig_fd = stories._file_digest
            stories._file_digest = lambda rel: (b"edited" if rel == "cards/edge-inp-gainp-design.yaml"
                                                else orig_fd(rel))
            try:
                stale_card = c.get("/api/stories/rt-single-photons").get_json()["stale"]
            finally:
                stories._file_digest = orig_fd
            check("H4 live-source card files are hashed: an edited card marks the story stale",
                  "cards/edge-inp-gainp-design.yaml" in files and stale_card is True
                  and c.get("/api/stories/rt-single-photons").get_json()["stale"] is False)
        finally:
            stories.STORIES = orig_dir
    blocks.append(('H4 / H5 / I5 stories (on a temp copy; the shipped file is not rewritten)', _b10))

    def _b11():  # studio-08-final: N1 / N2 / R4 (each fails on the pre-08 backend)
        import time as _time
        t0 = _time.perf_counter()
        r = c.post("/api/run", json={"card": "nitride-nanowire-vertical-pulse-design", "mode": "envelope"})
        dt = _time.perf_counter() - t0
        j = r.get_json() or {}
        check("N1 nanowire envelope run -> immediate 400 {kind: unsupported, message} (no job, no KeyError 'g2')",
              r.status_code == 400 and j.get("kind") == "unsupported" and "job_id" not in j
              and j.get("message") == jobs.NANOWIRE_ENVELOPE_UNSUPPORTED
              and "g2_op/T_j" in j["message"] and "Results" in j["message"] and dt < 5.0,
              f"{r.status_code} {str(j)[:120]} in {dt:.2f} s")
        d_pl = DeviceDesign.load(str(ROOT / "cards" / "nitride-cavity-set-design.yaml"))
        exc = ValueError("evaluate_envelope: every sample in the ranged box violates the F8 loading domain "
                         "(mu >= 1 - F_eff). Narrow the mu range or raise F_p/eta_capture.")
        e_pl = jobs.classify_error(d_pl, exc)
        check("N1 nitride (cycle_loading) all-dropped envelope: honest message, no drive.mu floor or remedy",
              e_pl.get("kind") == "envelope_all_dropped" and "f8_floor" not in e_pl and "suggestion" not in e_pl
              and "not a remedy" in e_pl.get("message", "") and "Narrow drive.mu" not in e_pl.get("message", ""),
              str(e_pl)[:160])
        e_lg = jobs.classify_error(DeviceDesign.load(str(ROOT / "cards" / "staged-device-design.yaml")), exc)
        check("N1 legacy F8 envelope error unchanged (f8_domain with floor + drive.mu suggestion)",
              e_lg.get("kind") == "f8_domain" and "f8_floor" in e_lg and "drive.mu" in e_lg.get("suggestion", ""))
        n_nw = jobs.NITRIDE_DISPLAY_GRID["ingan_gan_nanowire"][2]
        eta16 = jobs.eval_seconds("nanowire", n_nw) + jobs.eval_seconds("nanowire", 1)
        # re-pinned 2026-10-08: the vectorized injector (studio-p2c) cut the measured 16-point
        # nanowire grid from ~105 s to 33.9 s (vertical pulse, warm); estimate must sit near it
        check("N2 nanowire 16-point display grid estimate within the measured 25-60 s (33.9 s after p2c)",
              n_nw == 16 and 25.0 <= eta16 <= 60.0, f"{eta16:.1f} s")
        js = (ROOT / "fsim_studio" / "web" / "js" / "workspaces" / "designer.js").read_text(encoding="utf-8")
        m = re.search(r"nanowire:\s*\[([\d.]+),\s*([\d.]+)\]", js)
        check("N2 designer.js ETA mirrors the re-measured jobs.ETA_TABLE['nanowire'] (not the old 3.6/8.9)",
              m is not None and (float(m.group(1)), float(m.group(2))) == tuple(jobs.ETA_TABLE["nanowire"])
              and float(m.group(2)) > 100.0,
              m.group(0) if m else "no nanowire ETA in designer.js")
        st = stories.get_story("rt-single-photons")
        snake = re.compile(r"\b[a-z0-9]+(?:_[a-z0-9]+)+\b")
        strip = re.compile(r"\S*/\S*")
        bad = {sc["id"]: snake.findall(strip.sub(" ", (sc.get("claim") or "") + " " + (sc.get("notes") or "")))
               for sc in st["scenes"]}
        bad = {k: v for k, v in bad.items() if v}
        s13 = next(sc for sc in st["scenes"] if sc["id"] == "s13")
        check("R4 story claims + notes in room language: no raw snake_case tokens; s13 says "
              "'passes optically, fails in hardware'",
              not bad and "passes optically, fails in hardware" in s13["claim"]
              and "passes optically, fails in hardware" in s13["notes"] and not st.get("prose_errors"),
              str(bad)[:160])
        say = stories._render_number({"name": "x", "value": "no_idealized_pass", "raw": "no_idealized_pass",
                                      "source": {"kind": "verdict"}}, None, "say", None)
        check("R4 {name.say} maps printed status tokens to room language; the raw value is kept",
              say == "has no optical pass"
              and any(n.get("raw") == "pass_hardware_infeasible" for n in s13["numbers"]), say)
    blocks.append(("studio-08-final N1 / N2 / R4", _b11))

    for title, fn in blocks:
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 -- a crash in one block is a FAIL, not an abort
            check(f"{title}: block raised", False, f"{type(exc).__name__}: {exc}"[:200])


def p2a_checks(c, jm, tmp):
    """studio-p2a: preset builder (every axis value), Save as card (cards/studio/, 409 on shipped
    names, DeviceDesign.save round-trip), response animation (POST /api/animate: frames equal to
    direct evaluate(), per-frame cache shared with /api/run, frame SSE events). Every check here
    fails on the pre-p2a backend (no /api/animate, no meta on presets/apply, saves into cards/)."""
    import numpy as np

    import fsim_studio
    from fsim_core import presets
    from fsim_core.device import DeviceDesign, evaluate
    from fsim_studio import api_cards, jobs
    from fsim_studio.api_cards import design_from_dict, design_to_dict
    from fsim_studio.serialize import to_jsonable
    blocks = []
    SUB = ("g2_op", "collected_flux_pulsed_s", "collected_flux_delivered_s", "brightness_per_pulse",
           "T_j_op", "eps_op", "rho_op", "T_c")

    def _presets():
        base = {"dot": "chatzarakis-class", "template": "GaAs", "cavity": "none", "drive": "cw-electrical",
                "injection": "standard"}
        axes = {"dot": presets.DOT_PRESETS, "template": presets.TEMPLATE_PRESETS, "cavity": presets.CAVITY_PRESETS,
                "drive": presets.DRIVE_PRESETS, "injection": presets.INJECTION_PRESETS}
        bad, n = [], 0
        for ax, opts in axes.items():
            for key in opts:
                sel = dict(base, **{ax: key})
                r = c.post("/api/presets/apply", json=sel)
                pj = r.get_json() or {}
                d_ref = presets.preset_device(sel["dot"], sel["template"], sel["cavity"], sel["drive"])
                presets.apply_injection_preset(d_ref, sel["injection"])
                run = c.post("/api/run", json={"card": "preset-check", "design": pj.get("design"), "mode": "point"}).get_json() or {}
                s = (run.get("result") or {}).get("scalars") or {}
                direct = to_jsonable(evaluate(d_ref))["scalars"]
                meta = pj.get("meta") or {}
                ok = (r.status_code == 200 and pj["design"] == design_to_dict(d_ref)
                      and meta.get("dot.delta_xx", {}).get("tag") == "A" and "drive.mu" in meta
                      and isinstance(s.get("g2_op"), float) and s["g2_op"] == direct["g2_op"]
                      and (s.get("g2_op_valid") is not False or bool(s.get("g2_op_invalid_reason"))))
                n += 1
                if not ok:
                    bad.append(f"{ax}={key}")
        check("p2a presets: every axis value (15) composes like preset_device + injection, carries META "
              "tags, and evaluates clean (g2_op == direct evaluate; an invalid g2 carries its reason)",
              n == 15 and not bad, ",".join(bad) or f"{n} combinations")
    blocks.append(("p2a presets", _presets))

    def _save():
        user = Path(api_cards.USER_CARDS)
        pj = c.post("/api/presets/apply", json={"dot": "piezo-variant", "template": "GaAs-Si", "cavity": "micropillar",
                                                "drive": "pulsed-electrical", "injection": "dual-aperture-contact"}).get_json()
        r = c.put("/api/cards/p2a-verify-preset", json={"design": pj["design"], "presets": {"dot": "piezo-variant"}})
        rj = r.get_json() or {}
        f = user / "p2a-verify-preset.yaml"
        back = DeviceDesign.load(f) if f.is_file() else None
        want = design_from_dict(pj["design"])
        want.name = "p2a-verify-preset"
        d1 = DeviceDesign.load(f) if f.is_file() else None
        f2 = Path(tmp.name) / "roundtrip.yaml"
        if d1 is not None:
            d1.save(f2)
        check("p2a save as card: new slug saved to cards/studio (temp here) via DeviceDesign.save; YAML round-trips "
              "(load == sent design, save/load again identical), not into cards/",
              r.status_code == 200 and rj.get("saved") is True and back is not None
              and design_to_dict(back) == design_to_dict(want)
              and design_to_dict(DeviceDesign.load(f2)) == design_to_dict(back)
              and not (fsim_studio.CARDS / "p2a-verify-preset.yaml").exists(), str(rj)[:120])
        lst = c.get("/api/cards").get_json()
        got = c.get("/api/cards/p2a-verify-preset").get_json() or {}
        check("p2a saved card listed at once (saved: true) and served by GET /api/cards/<name>",
              any(x["name"] == "p2a-verify-preset" and x.get("saved") for x in lst)
              and got.get("design") == design_to_dict(want) and got.get("saved") is True)
        again = c.put("/api/cards/p2a-verify-preset", json={"design": pj["design"]})
        shipped = [c.put(f"/api/cards/{n}", json={"design": pj["design"]}).status_code
                   for n in ("staged-device-design", "edge-inp-gainp-design", "chatzarakis")]
        bad = [c.put(f"/api/cards/{n}", json={"design": pj["design"]}).status_code
               for n in ("Bad_Name", "unsaved", "x--y", "trailing-")]
        (user / "foreign-card.yaml").write_text((fsim_studio.CARDS / "staged-device-design.yaml").read_text(encoding="utf-8"),
                                               encoding="utf-8")
        foreign = c.put("/api/cards/foreign-card", json={"design": pj["design"]}).status_code
        bad_presets = [c.put("/api/cards/p2a-bad-presets", json={"design": pj["design"], "presets": x}).status_code
                       for x in ("piezo", ["dot", "x"], 5, True)]
        stray = sorted(q.name for q in user.glob("*.tmp"))
        check("review 7: a non-object presets body is a 400 (not a 500), saves nothing and leaves no <name>.<pid>.tmp behind",
              bad_presets == [400] * 4 and stray == [] and not (user / "p2a-bad-presets.yaml").exists(),
              f"{bad_presets} stray {stray}")
        # a failure after the temp file is written (round-trip / marker step) must also clean the temp file
        orig_read = api_cards._read_yaml
        api_cards._read_yaml = lambda pth: (_ for _ in ()).throw(RuntimeError("probe"))
        try:
            try:
                api_cards.save_card("p2a-tmp-probe", pj["design"], None)
            except RuntimeError:
                pass
        finally:
            api_cards._read_yaml = orig_read
        stray2 = sorted(q.name for q in user.glob("*.tmp"))
        check("review 7: an exception inside save_card unlinks the temp file (finally)", stray2 == [] and not (user / "p2a-tmp-probe.yaml").exists(), str(stray2))
        check("p2a save rules: own saved card may be re-saved; shipped names 409; non-slug / reserved 400; "
              "a non-Studio file in the user folder 409",
              again.status_code == 200 and shipped == [409, 409, 409] and bad == [400, 400, 400, 400] and foreign == 409,
              f"{again.status_code} {shipped} {bad} {foreign}")
    blocks.append(("p2a save as card", _save))

    def _animate():
        d_st = DeviceDesign.load(fsim_studio.CARDS / "staged-device-design.yaml")
        body = {"card": "staged-device-design", "param": "thermal.T_hs", "lo": 20.0, "hi": 300.0, "n": 8, "mode": "point"}
        t0 = time.time()
        r = c.post("/api/animate", json=body)
        rj = r.get_json() or {}
        snap = wait_job(c, rj["job_id"]) if r.status_code == 200 and not rj.get("result") else {"state": "done", "result": rj.get("result")}
        res = snap.get("result") or {}
        vals = [float(x) for x in np.linspace(20.0, 300.0, 8)]
        mism = []
        for i, v in enumerate(vals):
            d = copy.deepcopy(d_st)
            d.thermal.T_hs = v
            ds = to_jsonable(evaluate(d))["scalars"]
            fr = (res.get("frames") or [None] * 8)[i] or {}
            fs = fr.get("scalars") or {}
            if fr.get("value") != v or any(fs.get(k, "absent") != ds.get(k, "absent") for k in SUB) \
                    or fr.get("tag_chain") != "A" or "pre-retention" not in fr.get("labels", []):
                mism.append(i)
        check("p2a animate staged T_hs (n=8, full grid): 8 frames, values = linspace, frame scalars exactly equal "
              "direct evaluate() (g2_op ... T_c from the run), tag + labels per frame",
              r.status_code == 200 and snap["state"] == "done" and res.get("values") == vals
              and res.get("computed") == list(range(8)) and not mism and res.get("grid") == "full",
              f"mismatch {mism}; {time.time() - t0:.1f} s")
        t1 = time.time()
        r2 = c.post("/api/animate", json=body).get_json() or {}
        dt = time.time() - t1
        check("p2a animate: the repeated animation is served from the per-frame cache at once (cached, < 0.5 s, identical)",
              r2.get("cached") is True and dt < 0.5 and r2.get("result", {}).get("frames") == res.get("frames"), f"{dt * 1000:.0f} ms")
        d3 = copy.deepcopy(d_st)
        d3.thermal.T_hs = vals[3]
        rr = c.post("/api/run", json={"card": "staged-device-design", "design": design_to_dict(d3), "mode": "point"}).get_json() or {}
        check("p2a animate: a frame is the /api/run cache entry for the same design (shared key and run_id)",
              rr.get("cached") is True and rr.get("run_id") == res["frames"][3]["run_id"], f"{rr.get('cached')} {rr.get('run_id')}")

        # edge headline, operating point per frame, streamed frames
        d_e = DeviceDesign.load(fsim_studio.CARDS / "edge-inp-gainp-design.yaml")
        eb = {"card": "edge-inp-gainp-design", "param": "thermal.T_hs", "lo": 231.0, "hi": 297.0, "n": 3, "mode": "headline"}
        r = c.post("/api/animate", json=eb).get_json() or {}
        resp, events = sse_events(c, r["job_id"])
        kinds = [e[0] for e in events]
        fidx = sorted(e[1]["index"] for e in events if e[0] == "frame")
        done = next((e[1] for e in events if e[0] == "done"), None) or {}
        eres = done.get("result") or {}
        ok_e = []
        for i, v in enumerate([231.0, 264.0, 297.0]):
            dh = copy.deepcopy(d_e)
            dh.thermal.T_hs = v
            dh.drive.finite_pulse = True
            dh.ret.tau_cap_scales_with_density = False
            dh.drive.cw = False
            ds = to_jsonable(evaluate(dh, T_grid=[v]))["scalars"]
            fs = (eres.get("frames") or [{}] * 3)[i].get("scalars") or {}
            ok_e.append(all(fs.get(k, "absent") == ds.get(k, "absent") for k in SUB))
        check("p2a animate edge headline (n=3, op grid): frame SSE events stream (one per frame, before done); "
              "scalars equal direct evaluate of the headline switches at [T]; T_c null with its single-T reason path",
              r.get("cached") is False and resp.mimetype == "text/event-stream" and fidx == [0, 1, 2]
              and kinds[-1] == "done" and max(i for i, k in enumerate(kinds) if k == "frame") < kinds.index("done")
              and all(ok_e) and eres.get("grid") == "op" and jobs.HEADLINE_LABEL in eres["frames"][0]["labels"]
              and eres["frames"][0]["scalars"].get("T_c") is None, f"{kinds.count('frame')} frame events; {ok_e}")
        refused = [c.post("/api/animate", json=dict(body, **x)).status_code for x in
                   ({"param": "drive.nope"}, {"param": "drive.mode"}, {"n": 121}, {"n": 1}, {"lo": 5.0, "hi": 5.0},
                    {"mode": "envelope"})]
        check("p2a animate refuses non-META / non-numeric params, n outside 2..120, lo == hi, envelope mode (400)",
              refused == [400] * 6, str(refused))
        # review finding 1: the animated range must sit inside the META validity band; n must be an integer 2..120
        from fsim_core import design_meta
        band = design_meta.META["thermal.T_hs"]
        rb = [c.post("/api/animate", json=dict(body, **x)) for x in
              ({"lo": -50.0, "hi": 300.0}, {"lo": 20.0, "hi": 1000.0}, {"lo": band["lo"] - 1.0, "hi": band["hi"]},
               {"n": 0}, {"n": 2.5}, {"n": "7"}, {"n": True}, {"n": -3}, {"lo": "abc"})]
        edge = c.post("/api/animate", json=dict(body, lo=band["lo"], hi=band["hi"], n=2))
        ej = edge.get_json() or {}
        if ej.get("job_id") and not ej.get("cached"):
            c.delete(f"/api/jobs/{ej['job_id']}")
        check("review 1: animate refuses lo/hi outside the META band (T_hs -50..300, 20..1000, 3..350), n 0 / 2.5 / '7' / true / -3 "
              "and a non-number lo (400 with a reason); the band edges themselves (4..350, n=2) are accepted",
              [r.status_code for r in rb] == [400] * 9 and edge.status_code == 200 and ej.get("n") == 2
              and "validity band" in ((rb[0].get_json() or {}).get("error", "")) and "integer" in ((rb[4].get_json() or {}).get("error", "")) and "2..120" in ((rb[3].get_json() or {}).get("error", "")),
              f"{[r.status_code for r in rb]} {edge.status_code}")
        nw = c.post("/api/animate", json={"card": "nitride-nanowire-vertical-pulse-design", "param": "thermal.T_hs",
                                          "lo": 250.0, "hi": 300.0, "mode": "point"}).get_json() or {}
        if nw.get("job_id"):
            c.delete(f"/api/jobs/{nw['job_id']}")
        check("p2a animate default frame count: 16 for nanowire (48 otherwise), ETA from the measured table",
              nw.get("n") == 16 and nw.get("eta_s", 0) > 1.0, f"n={nw.get('n')} eta={nw.get('eta_s')}")
    blocks.append(("p2a animate", _animate))

    def _animate_honesty():
        """studio-p2d Coder A: physics brief section 7 (animation honesty rules 3-7) and the
        Studio-card verdict refusal. Every check fails on the pre-p2d backend."""
        from fsim_studio import animate as anim
        d_st = DeviceDesign.load(fsim_studio.CARDS / "staged-device-design.yaml")
        # rule 4: T_c from the frame's own full-grid evaluate(); never from an op frame
        body = {"card": "staged-device-design", "param": "thermal.T_hs", "lo": 20.0, "hi": 300.0, "n": 8, "mode": "point"}
        rj = c.post("/api/animate", json=body).get_json() or {}
        res = rj.get("result") or wait_job(c, rj["job_id"]).get("result") or {}
        fr = res.get("frames") or []
        tc_direct = to_jsonable(evaluate(d_st))["scalars"].get("T_c")
        check("p2d rule 4: full-grid frames carry grid 'full' and T_c equal to the direct evaluate() T_c",
              len(fr) == 8 and all(f.get("grid") == "full" for f in fr) and tc_direct is not None
              and all(f["scalars"].get("T_c") == tc_direct for f in fr), f"{[f.get('grid') for f in fr][:2]} {tc_direct}")
        eb = {"card": "edge-inp-gainp-design", "param": "thermal.T_hs", "lo": 231.0, "hi": 297.0, "n": 2, "mode": "headline"}
        er = c.post("/api/animate", json=eb).get_json() or {}
        eres = er.get("result") or wait_job(c, er["job_id"]).get("result") or {}
        ef = eres.get("frames") or []
        check("p2d rule 4: a grid 'op' frame is marked grid op and never carries a T_c (null + its reason, even if the "
              "result held one); a synthetic op frame built from a result WITH T_c still has T_c null",
              len(ef) == 2 and all(f.get("grid") == "op" and f["scalars"].get("T_c") is None
                                   and "never read from an op frame" in f["scalars"].get("T_c__nan_reason", "") for f in ef)
              and anim.frame_of(0, 1.0, {"scalars": {"T_c": 250.0, "g2_op": 0.1}}, "op")["scalars"]["T_c"] is None
              and anim.frame_of(0, 1.0, {"scalars": {"T_c": 250.0}}, "full")["scalars"]["T_c"] == 250.0)
        # rule 3: the flags the player must snap travel with the frame (copied, never derived here)
        ff = anim.frame_of(0, 1.0, {"scalars": {"g2_op": 0.1, "one_pair_valid": True, "set_feasible": False,
                                                 "flat_band": True, "depletion_regime": "flat_band"}}, "full")["scalars"]
        check("p2d rule 3: one_pair_valid / set_feasible / flat_band / depletion_regime ride along on a frame verbatim",
              ff.get("one_pair_valid") is True and ff.get("set_feasible") is False and ff.get("flat_band") is True
              and ff.get("depletion_regime") == "flat_band")
        # rule 5: the first runaway ends the displayable range
        fakes = [{"index": i, "value": float(i), "scalars": {"runaway": i >= 3}} for i in range(6)]
        an = anim.annotate_frames(fakes)
        check("p2d rule 5: runaway ends the displayable range: runaway_at is the first runaway frame and it and every "
              "later frame is flagged after_runaway (earlier frames and the cached originals untouched)",
              an["runaway_at"] == 3 and [bool(f.get("after_runaway")) for f in an["frames"]] == [False] * 3 + [True] * 3
              and "after_runaway" not in fakes[4] and anim.annotate_frames(fakes[:3])["runaway_at"] is None
              and res.get("runaway_at", "missing") is None)
        # rule 6: densities
        dd = res.get("density") or {}
        sug = dd.get("suggest") or {}
        ok5 = anim.density_report("thermal.T_hs", "legacy", [float(x) for x in range(150, 201, 5)],
                                  [{"grid": "full", "scalars": {"T_c": 175.5}}])
        check("p2d rule 6: T_hs frames 40 K apart around T_c (175.5 K) are reported too coarse (> 5 K within +-20 K) "
              "with a narrowed lo/hi/n suggestion of 9 frames at 5 K; a 5 K run passes",
              dd.get("rule") == "arrhenius" and dd.get("ok") is False and abs(dd.get("max_step", 0) - 40.0) < 1e-9
              and sug.get("n") == 9 and abs(sug["hi"] - sug["lo"] - 40.0) < 1e-9 and ok5["ok"] is True,
              str(dd)[:160])
        check("p2d rule 6: nitride Stark bias gets one frame per 0.2 V (0..5 V -> 26; 0..2 V -> 11) and a coarser run is reported",
              anim.default_frames("planar", "ingan_gan_planar", "drive.V", 0.0, 5.0) == 26
              and anim.default_frames("planar", "ingan_gan_planar", "drive.V", 0.0, 2.0) == 11
              and anim.default_frames("edge", "legacy", "drive.V", 0.0, 5.0) == 48
              and anim.default_frames("nanowire", "ingan_gan_nanowire", "thermal.T_hs", 250.0, 300.0) == 16
              and anim.density_report("drive.V", "ingan_gan_planar", [0.0, 0.5, 1.0], [])["ok"] is False
              and anim.density_report("drive.V", "ingan_gan_planar", [0.0, 0.2, 0.4], [])["ok"] is True)
        nz = c.post("/api/animate", json={"card": "nitride-cavity-pulse-design", "param": "drive.V", "lo": 0.0, "hi": 2.0,
                                          "mode": "point"}).get_json() or {}
        if nz.get("job_id"):
            c.delete(f"/api/jobs/{nz['job_id']}")
        check("p2d rule 6: POST /api/animate on a nitride card, drive.V 0..2 V, no n -> 11 frames", nz.get("n") == 11, str(nz.get("n")))
        # rule 7: an animated [A] input is a one-at-a-time slice; envelope runs are not animated
        pa = c.post("/api/animate", json={"card": "staged-device-design", "param": "dot.gamma_scale", "n": 2, "mode": "point"}).get_json() or {}
        pr = pa.get("result") or (wait_job(c, pa["job_id"]).get("result") if pa.get("job_id") else {}) or {}
        check("p2d rule 7: an animated [A] parameter's result carries the slice note 'one range at a time, not a prediction'; "
              "an [E] parameter's carries none; envelope mode is still refused (400)",
              pr.get("param_tag") == "A" and pr.get("slice_note") == "one range at a time, not a prediction"
              and res.get("slice_note") is None
              and c.post("/api/animate", json=dict(body, mode="envelope")).status_code == 400,
              f"{pr.get('param_tag')} {pr.get('slice_note')}")
        # saved-card verdicts
        vu = c.get("/api/cards/unsaved/verdict")
        vj = vu.get_json() or {}
        vs = c.get("/api/cards/p2a-verify-preset/verdict")
        vsj = vs.get_json() or {}
        vk = (c.get("/api/cards/staged-device-design/verdict").get_json() or {})
        check("p2d card_verdict refuses the unsaved draft and a Studio-saved card with a reason (200, match none, no lines, "
              "campaign null, refused: unsaved|saved); a shipped card is still matched",
              vu.status_code == 200 and vj.get("match") == "none" and vj.get("lines") == [] and vj.get("refused") == "unsaved"
              and "unsaved" in vj.get("reason", "").lower()
              and vs.status_code == 200 and vsj.get("match") == "none" and vsj.get("lines") == [] and vsj.get("refused") == "saved"
              and "saved" in vsj.get("reason", "").lower() and vsj.get("campaign") is None
              and vk.get("refused") is None and "match" in vk,
              f"{vu.status_code} {vj.get('refused')} / {vs.status_code} {vsj.get('refused')} / {vk.get('match')}")
    blocks.append(("p2d animate honesty", _animate_honesty))

    for title, fn in blocks:
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 -- a crash in one block is a FAIL, not an abort
            check(f"{title}: block raised", False, f"{type(exc).__name__}: {exc}"[:200])


PLACEHOLDER_STRIP = __import__("re").compile(r"\{[^{}]*\}")


def main():
    tmp = tempfile.TemporaryDirectory(prefix="fsim_studio_verify_")
    os.environ["FSIM_STUDIO_CACHE"] = tmp.name

    import numpy as np

    import fsim_studio
    from fsim_core import design_meta
    from fsim_core.device import DeviceDesign, evaluate, evaluate_envelope
    from fsim_studio import cache, jobs, stories
    from fsim_studio.api_cards import design_from_dict, design_to_dict
    from fsim_studio.serialize import dumps, to_jsonable
    from fsim_studio.server import create_app

    from fsim_studio import api_cards
    api_cards.USER_CARDS = Path(tmp.name) / "user_cards"  # Save as card writes here, never into cards/
    jm = jobs.JobManager(max_workers=4)
    app = create_app(job_manager=jm)
    c = app.test_client()
    c.environ_base["HTTP_X_FSIM_STUDIO"] = "1"  # mutating requests carry the Studio header (H13)

    try:
        # 1 health
        h = c.get("/api/health").get_json()
        check("health ok + fsim_core hash", h["ok"] is True and h["fsim_core_hash"] == fsim_studio.fsim_core_hash()
              and len(h["fsim_core_hash"]) == 64, h["fsim_core_hash"][:12])

        # 2 cards list
        cards = c.get("/api/cards").get_json()
        design_cards = [x for x in cards if x["kind"] == "design"]
        expected = sorted(p.stem for p in (ROOT / "cards").glob("*.yaml")
                          if "design:" in p.read_text(encoding="utf-8"))
        platforms = {x["name"]: x["platform"] for x in design_cards}
        check("cards: 14 design cards with platforms",
              len(design_cards) == 14 and sorted(platforms) == expected
              and all(v in ("legacy", "ingan_gan_planar", "ingan_gan_nanowire") for v in platforms.values()),
              f"{len(design_cards)} design, platforms {sorted(set(platforms.values()))}")
        check("cards: parameter cards listed as param",
              any(x["name"] == "chatzarakis" and x["kind"] == "param" for x in cards))

        # 3 GET card round-trip
        ok_rt, bad = True, []
        for name in expected:
            g = c.get(f"/api/cards/{name}").get_json()
            ref = DeviceDesign.load(ROOT / "cards" / f"{name}.yaml")
            back = design_from_dict(json.loads(json.dumps(g["design"])))
            if design_to_dict(back) != design_to_dict(ref) or g["platform"] != ref.platform:
                ok_rt = False
                bad.append(name)
        check("GET card round-trips through DeviceDesign.load (all 14)", ok_rt, ",".join(bad))
        g = c.get("/api/cards/staged-device-design").get_json()
        check("GET card meta carries tag/unit/band/source",
              g["meta"]["dot.delta_xx"]["tag"] == "A" and g["meta"]["dot.delta_xx"]["band"] == [2.0, 14.0])

        # 4 point run (legacy, inline)
        d_staged = DeviceDesign.load(ROOT / "cards" / "staged-device-design.yaml")
        r = c.post("/api/run", json={"card": "staged-device-design", "mode": "point"})
        raw_text = r.get_data(as_text=True)
        rj = r.get_json()
        res = rj["result"]
        direct = to_jsonable(evaluate(d_staged))
        g2 = res["curves"]["g2"]
        check("point run staged: inline result, g2 in [0,1], tag_chain",
              rj["cached"] is False and all(v is None or 0.0 <= v <= 1.0 for v in g2)
              and res["tag_chain"] == "A" and res["scalars"]["tag_chain"] == "[A]"
              and res["tag_chain_source"] == "fsim_core tag_chain",
              f"run_id {res['run_id']}")
        check("point run staged: g2 array identical to direct evaluate",
              g2 == direct["curves"]["g2"] and res["scalars"]["g2_op"] == direct["scalars"]["g2_op"])
        check("run_id format FS-<card>-<hex6>",
              res["run_id"] == "FS-staged-device-" + rj["run_id"].rsplit("-", 1)[-1]
              and len(res["run_id"].rsplit("-", 1)[-1]) == 6)
        check("NaN serialised as null (F_eff with cavity off)",
              math.isnan(evaluate(d_staged, T_grid=[77.0])["scalars"]["F_eff"])
              and res["scalars"]["F_eff"] is None and "NaN" not in raw_text and "Infinity" not in raw_text)
        ser = to_jsonable({"g2_op": float("nan"), "g2_op_invalid_reason": "runaway", "x": np.float64(np.inf)})
        check("NaN reason key <k>__nan_reason", ser == {"g2_op": None, "g2_op_invalid_reason": "runaway",
                                                         "x": None, "g2_op__nan_reason": "runaway"})
        check("pre-retention label from brightness_convention", "pre-retention" in res["labels"])

        # 5 cache hit + key change
        r2 = c.post("/api/run", json={"card": "staged-device-design", "mode": "point"}).get_json()
        check("cache hit on second identical run", r2["cached"] is True and r2["result"] == res
              and r2["run_id"] == rj["run_id"])
        payload = {"design": design_to_dict(d_staged), "mode": "point", "ranged": None, "T_grid": None}
        k1 = cache.make_key(payload)
        orig = cache.core_hash
        cache.core_hash = lambda: "0" * 64
        try:
            k2 = cache.make_key(payload)
        finally:
            cache.core_hash = orig
        check("cache key changes with the fsim_core hash", k1 != k2 and cache.make_key(payload) == k1)
        check("cache file written under FSIM_STUDIO_CACHE",
              (Path(tmp.name) / f"{k1}.json").is_file())

        # 6 headline run on edge-inp-gainp (science-agent point)
        d_edge = DeviceDesign.load(ROOT / "cards" / "edge-inp-gainp-design.yaml")
        d_pt = copy.deepcopy(d_edge)
        d_pt.thermal.T_hs = 230.0
        d_pt.emission.NA = 0.8
        d_pt.emission.R_back = 0.95
        d_pt.emission.L_um = 250.0
        d_pt.dot.delta_xx = 8.0
        d_pt.dot.gamma300 = 6.0
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "design": design_to_dict(d_pt),
                                     "mode": "headline", "T_grid": [230.0]}).get_json()
        snap = wait_job(c, r["job_id"])
        hres = snap.get("result") or {}
        dh = copy.deepcopy(d_pt)
        dh.drive.finite_pulse = True
        dh.ret.tau_cap_scales_with_density = False
        dh.drive.cw = False
        direct_h = to_jsonable(evaluate(dh, T_grid=[230.0]))["scalars"]
        check("headline run: g2_op equal to direct evaluate with switches",
              snap["state"] == "done" and hres["scalars"]["g2_op"] == direct_h["g2_op"]
              and abs(direct_h["g2_op"] - 0.98170) < 5e-6,
              f"g2_op {hres.get('scalars', {}).get('g2_op')}")
        check("headline run: no 'static (non-headline)' label",
              "static (non-headline)" not in hres.get("labels", []) and jobs.HEADLINE_LABEL in hres.get("labels", []))
        check("headline run: flux equal to direct",
              hres["scalars"]["collected_flux_pulsed_s"] == direct_h["collected_flux_pulsed_s"])

        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "mode": "point",
                                     "T_grid": [300.0]}).get_json()
        snap = wait_job(c, r["job_id"])
        sres = snap.get("result") or {}
        direct_s = to_jsonable(evaluate(d_edge, T_grid=[300.0]))["scalars"]
        check("plain point run of static edge card carries 'static (non-headline)'",
              snap["state"] == "done" and "static (non-headline)" in sres["labels"]
              and sres["scalars"]["g2_op"] == direct_s["g2_op"], f"g2_op {sres.get('scalars', {}).get('g2_op')}")

        # 7 progressive partial on a slow platform
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "mode": "headline",
                                     "T_grid": {"lo": 230, "hi": 300, "n": 3}}).get_json()
        resp, events = sse_events(c, r["job_id"])
        kinds = [e[0] for e in events]
        partial = next((e[1] for e in events if e[0] == "partial"), None)
        done = next((e[1] for e in events if e[0] == "done"), None)
        check("SSE stream: text/event-stream with state + done events",
              resp.mimetype == "text/event-stream" and "state" in kinds and kinds[-1] == "done", ",".join(kinds))
        check("progressive: partial (single T_hs) before full curve",
              partial is not None and len(partial["result"]["curves"]["g2"]) == 1
              and partial["result"]["partial"] is True and len(done["result"]["curves"]["g2"]) == 3
              and kinds.index("partial") < kinds.index("done"))

        # 8 envelope run vs direct evaluate_envelope (serial and process-pool map)
        r = c.post("/api/run", json={"card": "staged-device-design", "mode": "envelope"}).get_json()
        resp, events = sse_events(c, r["job_id"])
        snap = c.get(f"/api/jobs/{r['job_id']}").get_json()
        eres = snap["result"]
        ranged = design_meta.default_ranged(d_staged)
        env = evaluate_envelope(d_staged, ranged)
        env_j = to_jsonable(env)
        same = all(eres["bands"][k]["lo"] == env_j["bands"][k][0] and eres["bands"][k]["hi"] == env_j["bands"][k][1]
                   for k in env_j["bands"])
        check("envelope run equals direct evaluate_envelope bands exactly",
              snap["state"] == "done" and same
              and eres["scalar_bands"] == {k: list(v) for k, v in env_j["scalar_bands"].items()}
              and eres["tornado"] == env_j["tornado"] and eres["n_samples"] == env["n_samples"]
              and eres["curves"]["g2"] == env_j["mid"]["g2"], f"n_samples {eres['n_samples']}")
        prog = [e[1] for e in events if e[0] == "progress"]
        check("envelope emits progress k/n", len(prog) > 0 and prog[-1]["k"] == prog[-1]["n"],
              f"{len(prog)} progress events, n={prog[-1]['n'] if prog else None}")
        check("envelope labels: mid curve not a prediction",
              "all ranges at midpoint (not a prediction)" in eres["labels"])
        with cf.ProcessPoolExecutor(max_workers=4) as ex:
            ticks = []
            env_pool = evaluate_envelope(d_staged, ranged, map_fn=ex.map,
                                         progress=lambda k, n: ticks.append((k, n)))
        check("evaluate_envelope(map_fn=ProcessPool.map) == serial, exactly",
              dumps(env_pool) == dumps(env) and ticks[-1][0] == ticks[-1][1] == len(ticks))

        # 9 F8 domain error shape
        r = c.post("/api/run", json={"card": "nitride-cavity-pulse-design", "mode": "envelope",
                                     "T_grid": [300.0]}).get_json()
        snap = wait_job(c, r["job_id"]) if "job_id" in r else {}
        err = snap.get("error") or {}
        # studio-08 N1: this nitride (cycle_loading) card never reads drive.mu, so the all-dropped
        # envelope reports the honest message (raw F8 text kept as detail), not a drive.mu floor/remedy
        check("envelope all samples dropped on a nitride card -> {kind: envelope_all_dropped}, no drive.mu remedy",
              snap.get("state") == "error" and err.get("kind") == "envelope_all_dropped"
              and "f8_floor" not in err and "suggestion" not in err and "F8" in err.get("detail", ""),
              err.get("message", "")[:60])

        # 10 cancel
        r = c.post("/api/run", json={"card": "edge-inp-gainp-design", "mode": "headline",
                                     "T_grid": {"lo": 230, "hi": 300, "n": 5}}).get_json()
        cs = c.delete(f"/api/jobs/{r['job_id']}").get_json()
        time.sleep(0.2)
        snap = c.get(f"/api/jobs/{r['job_id']}").get_json()
        check("DELETE job cancels", cs["state"] == "cancelled" and snap["state"] == "cancelled"
              and "result" not in snap)

        # 11 campaigns
        camps = c.get("/api/campaigns").get_json()
        ids = {x["id"] for x in camps}
        check("campaigns list: rt_edge, nitride_cavity, nitride_geometry_stark, nitride_nanowire/full",
              {"rt_edge", "nitride_cavity", "nitride_geometry_stark", "nitride_nanowire/full"} <= ids,
              ",".join(sorted(ids)))
        check("campaigns list excludes c6_scratch_lw_quick",
              not any("c6_scratch_lw_quick" in i for i in ids))
        check("legacy F-series collections listed",
              all(x in ids for x in ("phase0", "phase3", "validation", "zhao"))
              and next(x for x in camps if x["id"] == "phase0")["kind"] == "legacy")
        check("pre-audit badge on tier_geometry / tier_device",
              next(x for x in camps if x["id"] == "tier_geometry")["stale"] is True
              and next(x for x in camps if x["id"] == "tier_device")["stale"] is True
              and next(x for x in camps if x["id"] == "rt_edge")["stale"] is False)
        rt = c.get("/api/campaigns/rt_edge").get_json()
        vline = next(l for l in (ROOT / "out/rt_edge/verdict.md").read_text(encoding="utf-8").splitlines()
                     if l.startswith("VERDICT:"))
        v0 = rt["verdicts"][0]
        check("rt_edge verdict word FAIL, eligible 16/768 (parsed from verdict.md)",
              v0["word"] == "FAIL" and v0["fields"]["eligible"] == "16/768" and v0["raw"] == vline
              and v0["fields"]["g2_min"] == "0.9817")
        cav = c.get("/api/campaigns/nitride_cavity").get_json()
        check("prose-tolerant parse (value with spaces)",
              cav["verdicts"][0]["fields"]["model"] == "finite electrical pulse"
              and cav["verdicts"][1]["fields"]["idealized_status"] == "pass_hardware_infeasible")
        nw = c.get("/api/campaigns/nitride_nanowire/full").get_json()
        check("nanowire: 16 VERDICT + 4 BEST lines",
              len(nw["verdicts"]) == 16 and len(nw["best"]) == 4
              and nw["best"][0]["fields"]["delivered_flux"] == "113102")
        check("campaign detail: columns with roles",
              {"name": "model_finite_pulse", "role": "input", "kind": "bool"} in rt["columns"]
              and any(col["role"] == "verdict" for col in rt["columns"]))

        sw = c.get("/api/campaigns/rt_edge/sweep?headline=1").get_json()
        check("sweep headline=1 on rt_edge returns 768 rows",
              sw["filtered"] == 768 and len(sw["rows"]) == 768 and sw["total"] == 3072)
        sw2 = c.get('/api/campaigns/rt_edge/sweep?headline=1&filter={"eligible_pulsed":true}'
                    "&cols=T_hs_K,g2_pulsed,collected_flux_pulsed_s").get_json()
        check("sweep filter + cols: 16 eligible headline rows, all at 230 K",
              sw2["filtered"] == 16 and all(row[0] == 230.0 for row in sw2["rows"])
              and min(row[1] for row in sw2["rows"]) == 0.9816951150623668)
        c.get("/api/campaigns/nitride_nanowire/full/sweep?limit=1")
        t0 = time.time()
        big = c.get("/api/campaigns/nitride_nanowire/full/sweep?limit=5000")
        dt = time.time() - t0
        bj = big.get_json()
        check("nanowire sweep (265 cols) answers < 2 s after first load",
              dt < 2.0 and len(bj["columns"]) == 265 and len(bj["rows"]) == 3216, f"{dt:.2f} s")
        rng = c.get('/api/campaigns/nitride_nanowire/full/sweep?filter={"T_hs":{"min":250,"max":300},'
                    '"row_kind":["core"]}&cols=T_hs').get_json()
        check("sweep range + list filters", rng["filtered"] > 0
              and all(250 <= row[0] <= 300 for row in rng["rows"]))

        # 12 cards misc
        put = c.put("/api/cards/staged-device-design", json={"design": g["design"]})
        check("PUT shipped card refused with 409", put.status_code == 409)
        val = c.post("/api/validate", json={"design": g["design"]}).get_json()
        check("validate: default_ranged + f8_floor",
              val["default_ranged"] == {k: list(v) for k, v in design_meta.default_ranged(d_staged).items()}
              and val["f8_floor"] == 0.0 and val["violations"] == [])
        pa = c.post("/api/presets/apply", json={"dot": "piezo-variant", "template": "GaAs",
                                                "cavity": "micropillar", "drive": "pulsed-electrical",
                                                "injection": "rti-quiet"}).get_json()
        check("presets/apply composes preset_device + injection",
              pa["design"]["dot"]["delta_xx"] == 8.0 and pa["design"]["drive"]["F_p"] == 0.5)
        meta = c.get("/api/meta").get_json()
        check("meta: META, ENV_DEFAULTS, presets", "dot.delta_xx" in meta["META"]
              and meta["ENV_DEFAULTS"]["dot.delta_xx"] == [1.5, 5.0] and "injection" in meta["presets"])
        th = c.get("/api/theme/plotly?mode=dark")
        check("theme endpoint returns a Plotly template", th.status_code == 200 and "layout" in th.get_json())

        # 13 compare report (bundle under the temp cache)
        rep = c.post("/api/compare/report", json={"entries": [
            {"label": "A", "run_id": rj["run_id"]}, {"label": "B", "run_id": eres["run_id"]}]})
        rp = rep.get_json()
        check("compare/report writes the designer_report bundle",
              rep.status_code == 200 and "report.png" in rp["files"] and "scalars_comparison.csv" in rp["files"]
              and Path(rp["path"]).resolve().is_relative_to(Path(tmp.name).resolve()))

        # 14 stories: bake a copy in a temp stories dir (the shipped file is never rewritten here)
        shipped = c.get("/api/stories/rt-single-photons").get_json()
        print(f"       (shipped story baked_at {shipped.get('baked_at')}, stale={shipped['stale']})")
        story_tmp = Path(tmp.name) / "stories"
        story_tmp.mkdir()
        (story_tmp / "rt-single-photons.json").write_text(
            (ROOT / "fsim_studio/stories/rt-single-photons.json").read_text(encoding="utf-8"), encoding="utf-8")
        orig_dir = stories.STORIES
        stories.STORIES = story_tmp
        try:
            bk = c.post("/api/stories/rt-single-photons/bake")
            st = c.get("/api/stories/rt-single-photons").get_json()
            nums = {(sc["id"], n["label"]): n for sc in st["scenes"] for n in sc["numbers"]}
            check("story bake: 13 scenes, every number resolved with a source",
                  bk.status_code == 200 and len(st["scenes"]) == 13 and st["bake_errors"] == []
                  and all(n.get("resolved_from") and "error" not in n for n in nums.values()),
                  f"{len(nums)} numbers, {bk.get_json().get('bake_seconds', 0):.1f} s")
            check("story numbers match out/ exactly (FAIL, 16/768, 0.9817)",
                  nums[("s05", "verdict")]["value"] == "FAIL"
                  and nums[("s04", "eligible (headline)")]["value"] == "16/768"
                  and nums[("s05", "pulsed g2 min (eligible)")]["raw"] == "0.9817"
                  and nums[("s05", "pulsed g2 min (eligible)")]["raw"] == v0["fields"]["g2_min"])
            check("story live number equals the headline API run",
                  nums[("s05", "live g2_op at 230 K (headline switches)")]["value"] == direct_h["g2_op"])
            cur = stories.current_hash(stories.load_story("rt-single-photons"))
            check("story baked_hash equals current hash (not stale)",
                  st["stale"] is False and st["baked_hash"] == cur)
            orig_h = stories.fsim_core_hash
            stories.fsim_core_hash = lambda: "f" * 64
            try:
                st2 = c.get("/api/stories/rt-single-photons").get_json()
            finally:
                stories.fsim_core_hash = orig_h
            check("story stale: true when the fsim_core hash moves", st2["stale"] is True)
        finally:
            stories.STORIES = orig_dir

        # 15 static + scene routing
        hz = c.get("/api/does-not-exist")
        check("unknown /api path -> JSON 404", hz.status_code == 404 and "error" in hz.get_json())
        sc = c.get("/api/scene/device?card=staged-device-design")
        check("scene route wired (200 SceneSpec, or 501 while scene.py is absent)",
              sc.status_code in (200, 501) and isinstance(sc.get_json(), dict), str(sc.status_code))
        v = c.get("/vendor/plotly/" + next(p.name for p in (ROOT / "fsim_studio/web/vendor/plotly").iterdir()))
        check("static /vendor served", v.status_code == 200)
        v.close()
        polish_checks(c, app, jm, tmp, direct_h)
        p2a_checks(c, jm, tmp)
    finally:
        jm.shutdown()

    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_pass}/{len(CHECKS)} studio api checks passed")
    try:
        tmp.cleanup()
    except OSError:
        pass
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
