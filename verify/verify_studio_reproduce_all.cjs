// FSIM Studio: every paper card reproduced through the API, compared with an independent engine run, then opened and run in
// the Composer; plus a full exercise of the Composer controls (R2 test audit 8, 9, 13, 5).
//
// Usage: node verify/verify_studio_reproduce_all.cjs [http://127.0.0.1:PORT]   (an already-running server; the server must
//        have been started with FSIM_STUDIO_USER_CARDS pointing at a temporary folder)
//        node verify/verify_studio_reproduce_all.cjs [--port 8840]            (starts its own private server, killed by PID)
// Prints "N/N reproduce-all checks passed" (KNOWN-OPEN checks count as passed and are listed); exits 0 iff all pass.
//
// 1. For EVERY card in cards/papers (not _fixture): POST /api/papers/<id>/reproduce (the call the Papers UI makes), poll the
//    job, and compare every item's `model` with a DIRECT engine run done by a separate python process that reads the card
//    YAML itself (engine.envelope central / engine.series rms) to 1e-9 relative. `model_series` must carry no error and
//    finite y where the card defines series.
// 2. For every card: Papers -> "Open in Composer" through the real UI, RUN there, and compare the Composer result with the
//    Papers result and with the direct engine run of the card device (finding 8: log-scale ranged leaves).
//    The readout text of every KPI must equal the page's own formatter applied to the API scalar.
// 3. Composer: swap every structure and every drive kind through the chip picker (a greyed entry must be greyed exactly when
//    validate says the swap adds an error), edit one enum, one bool and the tag selector (the posted spec carries them), RUN
//    each loaded tier, Series with HOM visibility (finding 5), KPI readouts equal the API values.
// KNOWN_OPEN lists audit findings being fixed elsewhere: their failures print KNOWN-OPEN and do not fail the exit code.

const path = require("path");
const fs = require("fs");
const os = require("os");
const { spawnSync } = require("child_process");
const { startStudioServer } = require("./studio_server.cjs");
const { assertPng } = require("./png_check.cjs");

const ROOT = path.resolve(__dirname, "..");
const argv = process.argv.slice(2);
let BASE = (argv.find((a) => /^https?:/.test(a)) || "").replace(/\/$/, "");
const ONLY_CARDS = argv.includes("--cards") ? argv[argv.indexOf("--cards") + 1].split(",") : null;   // debugging aid: restrict stages 1-2 to these cards
const SKIP_CARDS = argv.includes("--skip-cards");                                                      // debugging aid: stage 3 only
const PORT = Number(argv.includes("--port") ? argv[argv.indexOf("--port") + 1] : 8840);
const SHOTS = path.join(ROOT, ".workers", "devlib", "shots");
let OWN = null;

// ---- KNOWN_OPEN: {id, finding, desc, name: RegExp over the check name}
const KNOWN_OPEN = [
  { id: "F8", finding: 8, desc: "Open in Composer collapses a ranged log-scale leaf to the linear midpoint, the engine uses the geometric mean",
    name: /^composer==papers (ding2016|munsch2013_trumpet|somaschi2016_qd1|tomm2021)\b/ },
  { id: "F5", finding: 5, desc: "Composer Series with 'HOM visibility' sends the bare id hom_V, which has no 'intrinsic' convention",
    name: /Series with HOM visibility|API series accepts every observable the Composer can send/ },
  { id: "N2", finding: 0, desc: "NEW (found by this verifier): after a RUN that finished 'done', some KPI readouts of the shelf still show the empty dash although the API scalar is finite (repro: default example, swap the structure to planar_dbr, RUN, swap the drive to tpe, RUN: B, rate, F_P L and tau read '—')",
    name: /^swap (structure|drive) -> .*: KPI readouts equal the API scalar/, detail: /page '—'/ },
];

function loadPlaywright() {
  const tries = [path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright-core"),
    path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright"), "playwright-core", "playwright"];
  for (const t of tries) { try { return require(t); } catch { /* next */ } }
  console.error("playwright-core not found: npm i playwright-core under .workers/studio/node");
  process.exit(2);
}

const results = [];
const hits = {};
function check(name, ok, detail = "") {
  const kn = ok ? null : KNOWN_OPEN.find((k) => k.name.test(name) && (!k.detail || k.detail.test(String(detail))));
  results.push({ name, ok: !!ok || !!kn, known: kn ? kn.id : null });
  if (kn) hits[kn.id] = (hits[kn.id] || 0) + 1;
  console.log(`${ok ? "PASS" : kn ? `KNOWN-OPEN ${kn.id}` : "FAIL"}  ${name}${detail && !ok ? `  (${String(detail).slice(0, 500)})` : ""}`);
}

const HDR = { "X-FSIM-Studio": "1", "Content-Type": "application/json" };
async function getJSON(url) { const r = await fetch(`${BASE}${url}`); if (!r.ok) throw new Error(`${url} -> ${r.status}`); return r.json(); }
async function postJSON(url, body) { const r = await fetch(`${BASE}${url}`, { method: "POST", headers: HDR, body: JSON.stringify(body) }); return { status: r.status, body: await r.json().catch(() => ({})) }; }
async function awaitJob(resp, timeoutMs = 2400000) {
  if (resp.result) return { state: "done", result: resp.result };
  const t0 = Date.now();
  for (;;) {
    const snap = await getJSON(`/api/jobs/${resp.job_id}`);
    if (["done", "error", "cancelled"].includes(snap.state)) return snap;
    if (Date.now() - t0 > timeoutMs) return { state: "timeout" };
    await new Promise((r) => setTimeout(r, 400));
  }
}
const isNum = (v) => typeof v === "number" && Number.isFinite(v);
const relOK = (a, b, tol = 1e-9) => (!isNum(a) && !isNum(b)) || (isNum(a) && isNum(b) && Math.abs(a - b) <= tol * Math.max(Math.abs(a), Math.abs(b), 1e-300));

// ---- the page's own formatter (fsim_studio/web/js/ui/format.js `fmt`, 4 significant figures), copied verbatim
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const superscript = (n) => String(n).split("").map((c) => SUP[c] ?? c).join("");
function fmt(v, sig = 4) {
  if (v == null) return "n/a";
  if (typeof v !== "number" || !Number.isFinite(v)) return "n/a";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-3) { const [m, e] = v.toExponential(sig - 1).split("e"); return `${m}×10${superscript(Number(e))}`; }
  if (a >= 1000) return Math.round(v).toLocaleString("en-US");
  return Number(v.toPrecision(sig)).toString();
}

// ---- independent direct engine run (a separate python process; reads the YAML itself, never calls papers.reproduce)
const DIRECT = String.raw`
import sys, json, math, warnings
warnings.filterwarnings("ignore")
import yaml
from pathlib import Path
root, outp = sys.argv[1], sys.argv[2]
only = set(sys.argv[3].split(",")) if len(sys.argv) > 3 and sys.argv[3] else None
sys.path.insert(0, root)
from fsim_core.compose import engine
KEYS = ("g2_0", "hom_I", "B_lens", "tau_ps", "F_P_eff", "rate_detected_Hz")
def f(x):
    try:
        x = float(x)
    except Exception:
        return None
    return x if math.isfinite(x) else None
out = {}
for fp in sorted(Path(root, "cards", "papers").glob("*.yaml")):
    if fp.name.startswith("_") or (only and fp.stem not in only):
        continue
    raw = yaml.safe_load(fp.read_text(encoding="utf-8"))
    spec = {"schema": "fsim-composed/1", **raw["device"]}      # the card loader adds the schema key
    ids = {r["id"]: ("scalar", r) for r in raw.get("reported") or []}
    ids.update({s["id"]: ("series", s) for s in raw.get("series") or []})
    card = {"items": {}, "scalars": None, "error": None}
    try:
        sc = engine.evaluate_spec(spec)["scalars"]
        card["scalars"] = {k: f(sc.get(k)) for k in KEYS}
        card["valid"] = bool(sc["valid"])
    except Exception as e:
        card["error"] = type(e).__name__ + ": " + str(e)[:200]
    for e in raw.get("reproduce") or []:
        item = e["item"]
        try:
            if e["class"] == "P":
                kind, oid, conv = "predict", e["observable"], e.get("convention", "intrinsic")
            else:
                kind, ent = ids[item]
                if kind == "scalar":
                    oid, conv = ent["observable"], ent.get("convention", "intrinsic")
                else:
                    oid, conv = ent["y"]["observable"], ent["y"].get("convention", "intrinsic")
            req = oid if conv == "intrinsic" else oid + ":" + conv
            if kind in ("scalar", "predict"):
                env = engine.envelope(spec, [req], max_corners=0, n_lhs=0)
                card["items"][item] = {"kind": kind, "model": f(env[req]["central"])}
            else:
                pts = ent["points"]
                xs = [float(p[0]) for p in pts]
                py = [float(p[1]) for p in pts]
                sr = engine.series(spec, ent["x"]["path"], xs, [req])
                my = [float(v) for v in sr["y"][req]]
                r = [((m - p) / p) ** 2 for m, p in zip(my, py) if p != 0.0]
                card["items"][item] = {"kind": "series", "model": f(math.sqrt(sum(r) / len(r))) if r and all(math.isfinite(v) for v in my) else None}
        except Exception as ex:
            card["items"][item] = {"kind": "error", "model": None, "error": type(ex).__name__ + ": " + str(ex)[:160]}
    out[fp.stem] = card
Path(outp).write_text(json.dumps(out), encoding="utf-8")
`;

function runDirect() {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "fsim_direct_"));
  const script = path.join(tmp, "direct.py");
  const outp = path.join(tmp, "direct.json");
  fs.writeFileSync(script, DIRECT);
  const r = spawnSync(process.env.PYTHON || "python", ["-I", script, ROOT, outp, (ONLY_CARDS || []).join(",")], { cwd: tmp, env: { ...process.env, OMP_NUM_THREADS: "2", PYTHONPATH: "" }, encoding: "utf-8", maxBuffer: 1 << 26, timeout: 3600000 });
  if (r.status !== 0 || !fs.existsSync(outp)) throw new Error(`direct engine helper failed: ${(r.stderr || "").slice(-400)}`);
  const data = JSON.parse(fs.readFileSync(outp, "utf-8"));
  fs.rmSync(tmp, { recursive: true, force: true });
  return data;
}

async function settle(page, ms = 400) { await page.waitForTimeout(ms); }

(async () => {
  if (!BASE) { OWN = await startStudioServer(ROOT, PORT); BASE = OWN.base; }
  fs.mkdirSync(SHOTS, { recursive: true });
  const t0 = Date.now();
  let browser = null;
  try {
    // ================================================================ 1. API vs direct engine, every card
    let list = (await getJSON("/api/papers")).filter((p) => !String(p.id).startsWith("_"));
    if (ONLY_CARDS) list = list.filter((p) => ONLY_CARDS.includes(p.id));
    if (SKIP_CARDS) list = [];
    if (!ONLY_CARDS && !SKIP_CARDS) check("GET /api/papers lists the cards of cards/papers", list.length >= 30, `${list.length} cards`);
    const onDisk = fs.readdirSync(path.join(ROOT, "cards", "papers")).filter((f) => f.endsWith(".yaml") && !f.startsWith("_")).map((f) => f.slice(0, -5)).sort();
    if (!ONLY_CARDS && !SKIP_CARDS) check("the API lists exactly the cards on disk (none skipped)", JSON.stringify(list.map((p) => p.id).sort()) === JSON.stringify(onDisk), `${list.length} vs ${onDisk.length}`);
    console.log(`direct engine helper over ${onDisk.length} cards ...`);
    const direct = SKIP_CARDS ? {} : runDirect();
    const papersRes = {};
    for (const p of list) {
      const id = p.id;
      const post = await postJSON(`/api/papers/${id}/reproduce`, {});
      const snap = post.status < 300 ? await awaitJob(post.body) : { state: `http ${post.status}` };
      check(`API reproduce ${id}: job done`, snap.state === "done" && Array.isArray(snap.result?.items), snap.state + " " + JSON.stringify(snap.error || post.body.error || "").slice(0, 200));
      if (snap.state !== "done") continue;
      papersRes[id] = snap.result;
      const d = direct[id] || { items: {} };
      const bad = [];
      for (const it of snap.result.items) {
        const dd = d.items[it.id];
        if (!dd) { bad.push(`${it.id}: no direct value`); continue; }
        const apiModel = isNum(it.model) ? it.model : null;
        // an item the API marked unavailable/error carries no model; the direct run may still produce a number
        if (apiModel === null && ["unavailable", "error"].includes(it.outcome)) continue;
        if (!relOK(apiModel, dd.model)) bad.push(`${it.id}: API ${apiModel} vs direct ${dd.model}`);
      }
      check(`API model == direct engine run, every item of ${id} (${snap.result.items.length} items, 1e-9 rel)`, bad.length === 0, bad.join(" | "));
      const card = await getJSON(`/api/papers/${id}`);
      if ((card.series || []).length) {
        const ms = snap.result.model_series || [];
        const sbad = [];
        for (const s of card.series) {
          const m = ms.find((x) => x.id === s.id);
          if (!m) sbad.push(`${s.id}: missing`);
          else if (m.error) sbad.push(`${s.id}: ${m.error}`);
          else if (!(m.y || []).length || m.y.length !== s.points.length || !m.y.every(isNum)) sbad.push(`${s.id}: y not finite/complete`);
        }
        check(`model_series ${id}: every card series has a finite model curve, no error`, sbad.length === 0, sbad.join(" | "));
      }
    }

    // ================================================================ 2. Papers -> Open in Composer -> RUN, compare
    const { chromium } = loadPlaywright();
    browser = await chromium.launch({ args: ["--use-angle=d3d11", "--enable-gpu", "--ignore-gpu-blocklist"] });
    const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
    const page = await context.newPage();
    const errors = [];
    page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
    page.on("console", (m) => { if (m.type() === "error" && !/status of (4\d\d|5\d\d)/.test(m.text())) errors.push(m.text()); });
    const lastValidate = { body: null };
    page.on("request", (r) => { if (r.method() === "POST" && r.url().endsWith("/api/compose/validate")) { try { lastValidate.body = JSON.parse(r.postData()); } catch { /* ignore */ } } });

    const waitPlaque = () => page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status, null, { timeout: 90000 });
    const waitRunDone = (timeout = 240000) => page.waitForFunction(() => {
      const s = document.querySelector(".ws-compose .status-word");
      return s && ["done", "error", "cancelled"].includes(s.dataset.state);
    }, null, { timeout });
    /** click RUN, return {state, scalars?, resp} where scalars come from the job the page submitted. */
    async function runInComposer(urlPart = "/api/compose/run") {
      const respP = page.waitForResponse((r) => r.url().endsWith(urlPart) && r.request().method() === "POST", { timeout: 30000 });
      await page.click(".ws-compose .run-key");
      const resp = await respP;
      const body = await resp.json().catch(() => ({}));
      await waitRunDone();
      const state = await page.$eval(".ws-compose .status-word", (e) => e.dataset.state);
      let result = null;
      if (state === "done") {
        const snap = await awaitJob(body, 60000);
        result = snap.result || null;
      }
      return { state, result, status: resp.status(), body };
    }
    const readouts = () => page.$$eval(".ws-compose .readout", (els) => Object.fromEntries(els.map((e) => [e.dataset.readout, { text: (e.querySelector(".ro-value")?.textContent || "").replace(/\s+/g, ""), state: e.dataset.state }])));
    const KPI = ["g2_0", "hom_I", "B_first_lens", "rate_detected_Hz", "F_P_eff", "tau_ps"];
    const SC = { g2_0: "g2_0", hom_I: "hom_I", B_first_lens: "B_lens", rate_detected_Hz: "rate_detected_Hz", F_P_eff: "F_P_eff", tau_ps: "tau_ps" };
    // the page renders "10⁶" as 10<sup>6</sup>, so textContent holds plain digits: compare in that form
    const plainSup = (t) => t.replace(/⁻/g, "-").replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]/g, (c) => "⁰¹²³⁴⁵⁶⁷⁸⁹".indexOf(c)).replace(/\s+/g, "");
    async function readoutsMatch(scalars) {
      let bad = [];
      for (let attempt = 0; attempt < 8; attempt++) {          // RUN's "done" lands a tick before the shelf is painted
        const ro = await readouts();
        bad = [];
        for (const k of KPI) {
          const v = scalars[SC[k]];
          const want = plainSup(fmt(v));
          const got = ro[k]?.text ?? "(missing)";
          if (isNum(v) ? got !== want : !(ro[k]?.state === "nan" || ro[k]?.state === "empty")) bad.push(`${k}: page '${got}' vs fmt(API) '${want}'`);
        }
        if (!bad.length) break;
        await page.waitForTimeout(500);
      }
      return bad;
    }

    for (const p of list) {
      const id = p.id;
      const pr = papersRes[id];
      if (!pr) continue;
      try {
        await page.goto(`${BASE}/#/papers/${id}`, { waitUntil: "domcontentloaded" });
        await page.waitForSelector("#workspace .pp-composer", { timeout: 30000 });
        await page.waitForSelector('#workspace[data-ready="1"]', { timeout: 30000 });
        await page.locator(".pp-composer").click();
        await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
        await waitPlaque();
        await page.waitForFunction(() => !/Checking|VALIDATING/.test(document.querySelector(".ws-compose .plaque")?.textContent || ""), null, { timeout: 60000 });
        const plaque = await page.$eval(".ws-compose .plaque", (e) => ({ s: e.dataset.status, t: e.textContent.slice(0, 160) }));
        const seedNote = await page.$eval(".compose-seed", (e) => e.textContent).catch(() => "");
        const kinds = await page.$$eval(".ws-compose .kind-chip", (els) => Object.fromEntries(els.map((e) => [e.dataset.axis, e.dataset.kind])));
        const want = (await getJSON(`/api/papers/${id}`)).device;
        const kindsOK = ["material", "structure", "drive", "setup"].every((a) => kinds[a] === want[a].kind);
        check(`Open in Composer ${id}: the card's kinds load${seedNote ? "" : " (no seed note)"} and the plaque says ${plaque.s}`, kindsOK && !!seedNote && plaque.s === "pass", `${JSON.stringify(kinds)} ${plaque.t}`);
        if (plaque.s !== "pass") continue;
        const run = await runInComposer();
        if (run.state !== "done" || !run.result?.scalars) { check(`composer==papers ${id}: RUN completes`, false, `${run.state} ${JSON.stringify(run.body.error || "").slice(0, 200)}`); continue; }
        const sc = run.result.scalars;
        // (a) Composer scalars vs the Papers items on the same observable (intrinsic convention, scalar items only)
        const bad = [];
        for (const it of pr.items) {
          if (!["scalar", "predict"].includes(it.kind) || it.convention !== "intrinsic" || !isNum(it.model)) continue;
          const key = it.observable === "B_first_lens" ? "B_lens" : it.observable;
          if (!(key in sc) || !isNum(sc[key])) continue;
          if (!relOK(sc[key], it.model)) bad.push(`${it.id} (${it.observable}): Composer ${sc[key]} vs Papers ${it.model}`);
        }
        // (b) Composer scalars vs the direct engine run of the card device
        const d = direct[id]?.scalars;
        if (d) for (const k of Object.keys(d)) if (!relOK(sc[k], d[k])) bad.push(`${k}: Composer ${sc[k]} vs card device ${d[k]}`);
        check(`composer==papers ${id}: the Composer result equals the Papers result and the card device`, bad.length === 0, bad.join(" | "));
        const rb = await readoutsMatch(sc);
        check(`KPI readouts ${id} equal the API scalar through the page formatter`, rb.length === 0, rb.join(" | "));
      } catch (e) {
        check(`Open in Composer + RUN ${id}`, false, e.message.split("\n")[0]);
      }
    }

    // ================================================================ 3. Composer controls
    const reg = await getJSON("/api/compose/registry");
    const gotoComposer = async () => {
      await page.goto(`${BASE}/#/overview`, { waitUntil: "load" });
      await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
      await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
      await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "pass", null, { timeout: 90000 });
      await settle(page, 600);
    };
    await gotoComposer();
    const kindOf = (axis) => page.$eval(`.ws-compose .kind-chip[data-axis="${axis}"]`, (e) => e.dataset.kind);

    for (const [axis, regKey] of [["structure", "structures"], ["drive", "drives"]]) {
      for (const e of reg[regKey]) {
        const kind = e.kind;
        try {
          if (e.status === "needs_physics") continue;
          const cur = await kindOf(axis);
          const baseSpec = lastValidate.body?.spec;
          await page.click(`.ws-compose .kind-chip[data-axis="${axis}"]`);
          await page.waitForSelector(".cp-entry", { timeout: 15000 });
          await page.waitForFunction(() => !/Checking/.test(document.querySelector(".cp-note")?.textContent || ""), null, { timeout: 180000 });
          const info = await page.$eval(`.cp-entry[data-kind="${kind}"]`, (b) => ({ dis: b.getAttribute("aria-disabled") === "true" || b.dataset.incompatible === "true", compat: b.dataset.compat, why: (b.querySelector(".cp-why")?.textContent || "").trim() }));
          // independent expectation: validate the swap through the API; "fresh" errors = errors the current chain does not already have
          const baseV = await postJSON("/api/compose/validate", { spec: baseSpec });
          const baseErr = new Set((baseV.body.issues || []).filter((i) => i.level === "error").map((i) => `${i.path}|${i.message}`));
          const swapped = JSON.parse(JSON.stringify(baseSpec)); swapped[axis] = { kind, params: {} };
          const sv = await postJSON("/api/compose/validate", { spec: swapped });
          const fresh = (sv.body.issues || []).filter((i) => i.level === "error" && !baseErr.has(`${i.path}|${i.message}`));
          const expectGrey = fresh.length > 0 && kind !== cur;
          check(`picker ${axis} ${kind}: greyed exactly when validate adds an error${expectGrey ? " (reason shown)" : ""}`, info.dis === expectGrey && (!info.dis || info.why.length > 5), JSON.stringify({ ...info, fresh: fresh.length }));
          if (info.dis || kind === cur) { await page.keyboard.press("Escape"); await page.waitForSelector(".cp-entry", { state: "detached", timeout: 5000 }).catch(() => {}); continue; }
          const vP = page.waitForResponse((r) => r.url().endsWith("/api/compose/validate") && r.request().method() === "POST", { timeout: 30000 });
          await page.click(`.cp-entry[data-kind="${kind}"]`);
          await vP;
          await page.waitForFunction(([a, k]) => document.querySelector(`.ws-compose .kind-chip[data-axis="${a}"]`)?.dataset.kind === k, [axis, kind], { timeout: 10000 });
          await settle(page, 700);
          const spec = lastValidate.body.spec;
          const api = await postJSON("/api/compose/validate", { spec });
          const pl = await page.$eval(".ws-compose .plaque", (x) => x.dataset.status);
          check(`swap ${axis} -> ${kind}: the plaque (${pl}) equals /api/compose/validate (${api.body.ok ? "pass" : "fail"})`, spec[axis].kind === kind && (api.body.ok ? "pass" : "fail") === pl, JSON.stringify(api.body.issues || []).slice(0, 200));
          if (api.body.ok) {
            const run = await runInComposer();
            check(`swap ${axis} -> ${kind}: RUN completes for a chain validate accepted`, run.state === "done" && !!run.result?.scalars, `${run.state} ${JSON.stringify(run.body.error || "").slice(0, 200)}`);
            if (run.result?.scalars) {
              const rb = await readoutsMatch(run.result.scalars);
              check(`swap ${axis} -> ${kind}: KPI readouts equal the API scalar`, rb.length === 0, rb.join(" | "));
            }
          }
        } catch (err) {
          check(`swap ${axis} -> ${kind}`, false, err.message.split("\n")[0]);
          await page.keyboard.press("Escape").catch(() => {});
        }
      }
    }

    // ---- edit one enum, one bool, the tag selector: the posted spec carries them
    await gotoComposer();
    let enumDone = false, boolDone = false, tagDone = false;
    for (const axis of ["drive", "material", "structure", "setup"]) {
      if (enumDone && boolDone && tagDone) break;
      await page.click(`.ws-compose .mount-body[data-axis="${axis}"]`);
      await page.waitForSelector(".ws-compose .sheet .field-row", { timeout: 10000 }).catch(() => {});
      await settle(page, 400);
      if (!enumDone) {
        const sel = page.locator(".ws-compose .sheet .pf-item select.fr-select").first();
        if (await sel.count()) {
          const param = await sel.evaluate((e) => e.closest(".pf-item").dataset.param);
          const opts = await sel.evaluate((e) => [...e.options].map((o) => o.value));
          const cur = await sel.inputValue();
          const to = opts.find((o) => o !== cur);
          if (to) {
            const vP = page.waitForResponse((r) => r.url().endsWith("/api/compose/validate") && r.request().method() === "POST", { timeout: 20000 });
            await sel.selectOption(to);
            await vP;
            const leaf = lastValidate.body.spec[axis].params[param];
            check(`edit enum ${axis}.${param} -> ${to}: the posted spec carries it with the [A] tag`, leaf && leaf.value === to && leaf.tag === "A", JSON.stringify(leaf));
            enumDone = true;
          }
        }
      }
      if (!boolDone) {
        const sw = page.locator('.ws-compose .sheet .pf-item input[role="switch"]').first();
        if (await sw.count()) {
          const param = await sw.evaluate((e) => e.closest(".pf-item").dataset.param);
          const was = await sw.isChecked();
          const vP = page.waitForResponse((r) => r.url().endsWith("/api/compose/validate") && r.request().method() === "POST", { timeout: 20000 });
          await sw.evaluate((e) => e.click());
          await vP;
          const leaf = lastValidate.body.spec[axis].params[param];
          check(`edit bool ${axis}.${param}: the posted spec carries ${!was}`, leaf && leaf.value === !was, JSON.stringify(leaf));
          boolDone = true;
        }
      }
      if (!tagDone) {
        const row = page.locator(".ws-compose .sheet .pf-item").filter({ has: page.locator("input.fr-num") }).first();
        if (await row.count()) {
          const param = await row.evaluate((e) => e.dataset.param);
          const vP = page.waitForResponse((r) => r.url().endsWith("/api/compose/validate") && r.request().method() === "POST", { timeout: 20000 });
          await row.locator('.pf-tagsel .tag[data-tag="E"]').click();
          await vP;
          const leaf = lastValidate.body.spec[axis].params[param];
          check(`tag selector ${axis}.${param} -> E: the posted spec leaf carries tag E`, leaf && leaf.tag === "E", JSON.stringify(leaf));
          tagDone = true;
        }
      }
      await page.click(`.ws-compose .mount-body[data-axis="${axis}"]`).catch(() => {});
    }
    check("one enum, one bool and the tag selector were all exercised", enumDone && boolDone && tagDone, JSON.stringify({ enumDone, boolDone, tagDone }));

    // ---- each tier
    for (const t of reg.tiers || []) {
      try {
        await gotoComposer();
        await page.click('.ws-compose .start-btn[data-menu="tiers"]');
        await page.waitForSelector(".start-pop .start-item", { timeout: 5000 });
        await page.waitForTimeout(300);
        const chip = page.locator(`.start-pop .tier-chip[data-tier="${t.id}"]`);
        const off = await chip.isDisabled();
        if (off) { check(`tier ${t.id}: not mapped, the menu entry is disabled`, true); await page.keyboard.press("Escape"); continue; }
        await chip.click();
        if (t.handoff === "design") {
          await page.waitForFunction(() => location.hash.startsWith("#/design/"), null, { timeout: 10000 });
          check(`tier ${t.id}: hands off to the Designer (${await page.evaluate(() => location.hash)})`, true);
        } else {
          await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status, null, { timeout: 60000 });
          await settle(page, 900);
          const ps = await page.$eval(".ws-compose .plaque", (e) => e.dataset.status);
          check(`tier ${t.id}: opens here and validates`, ps === "pass", ps);
          if (ps === "pass") {
            const run = await runInComposer();
            check(`tier ${t.id}: RUN completes`, run.state === "done" && !!run.result?.scalars, `${run.state} ${JSON.stringify(run.body.error || "").slice(0, 200)}`);
            if (run.result?.scalars) check(`tier ${t.id}: KPI readouts equal the API scalar`, (await readoutsMatch(run.result.scalars)).length === 0);
          }
        }
      } catch (err) { check(`tier ${t.id}`, false, err.message.split("\n")[0]); }
    }

    // ---- Series with HOM visibility (finding 5)
    await gotoComposer();
    {
      // the UI sends, for the first convention of an observable, the bare id (composer.js: req = i === 0 ? o.id : `${o.id}:${c}`)
      const homObs = reg.observables.find((o) => o.id === "hom_V");
      const bare = homObs.id;
      await page.click('.ws-compose .seg[data-value="series"]');
      await page.waitForSelector(".ws-compose .rb-series .rb-path");
      await page.click(".ws-compose .rb-obs summary");
      for (const cb of await page.$$(".ws-compose .rb-obs input[type=checkbox]")) if (await cb.isChecked()) await cb.uncheck();
      // whatever request string the page now uses for HOM visibility (the first checkbox for hom_V)
      const request = await page.$eval('.ws-compose .rb-obs input[data-obs^="hom_V"]', (e) => e.dataset.obs);
      await page.check(`.ws-compose .rb-obs input[data-obs="${request}"]`);
      const run = await runInComposer("/api/compose/series");
      let ok = run.state === "done";
      let detail = `state ${run.state} ${JSON.stringify(run.body.error || "").slice(0, 160)}`;
      if (ok) {
        const y = run.result?.series?.y?.[request] || [];
        ok = y.length > 1 && y.every(isNum);
        detail = `y: ${JSON.stringify(y).slice(0, 120)}`;
        await page.waitForFunction(() => document.querySelector('.ws-compose [data-chart="series"] .chart-plot .main-svg'), null, { timeout: 20000 }).catch(() => {});
        const tr = await page.$eval('.ws-compose [data-chart="series"] .chart-plot', (gd) => (gd.data || []).filter((t) => (t.y || []).some((v) => v !== null)).length).catch(() => 0);
        ok = ok && tr >= 1;
      }
      check(`Series with HOM visibility (the request the Composer sends: ${request}) finishes with a finite curve drawn`, ok, detail);
      const api = await postJSON("/api/compose/series", { spec: lastValidate.body.spec, path: "T_K", values: [4, 10, 20], observables: [bare] });
      let apiOK = api.status === 400;
      let apiD = `status ${api.status}`;
      if (api.status < 300) {
        const snap = await awaitJob(api.body, 120000);
        apiOK = snap.state === "done" && (snap.result?.series?.y?.[bare] || []).every(isNum);
        apiD = `job ${snap.state} ${JSON.stringify(snap.error || "").slice(0, 160)}`;
      }
      check("API series accepts every observable the Composer can send, or rejects it at submit time (hom_V)", apiOK, apiD);
    }

    // ---- screenshots taken here are real images
    const shot = path.join(SHOTS, "reproduce-all-composer.png");
    await page.screenshot({ path: shot });
    const png = assertPng(shot);
    check("composer capture is a real image (pixel variance, colour share)", png.ok, png.reason);
    check("no console errors", errors.length === 0, errors.slice(0, 3).join(" | "));
  } catch (e) {
    check("verifier ran to the end", false, e.stack || e.message);
  } finally {
    if (browser) await browser.close();
    OWN?.stop();
  }
  const passed = results.filter((r) => r.ok).length;
  for (const k of KNOWN_OPEN) console.log(hits[k.id] ? `KNOWN-OPEN ${k.id} finding ${k.finding}: ${hits[k.id]} check(s) - ${k.desc}` : `[FIXED?] KNOWN_OPEN ${k.id} (finding ${k.finding}) matched nothing: delete the entry`);
  console.log(`${passed}/${results.length} reproduce-all checks passed (${results.filter((r) => r.known).length} KNOWN-OPEN, ${Math.round((Date.now() - t0) / 1000)} s)`);
  process.exit(passed === results.length ? 0 : 1);
})();
