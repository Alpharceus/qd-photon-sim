// FSIM Studio front-end smoke test + screenshots (Playwright, headless Chromium).
//
// Usage: node verify/verify_studio_web.cjs            (starts its own server on 8791 with a TEMP saved-card
//                                                       folder and cache, and stops it at the end)
//        node verify/verify_studio_web.cjs http://127.0.0.1:8791   (an already-running server: it must have been
//          started with FSIM_STUDIO_USER_CARDS=<temp dir>, and this process needs the same variable, otherwise
//          the Save-as-card checks would write into the real cards/studio/; the script refuses to run)
// Playwright is not a repo dependency: this resolves playwright-core from
// node_modules (temp install), then from NODE_PATH / global.
// Exits 0 iff every check passes; prints "N/N studio web checks passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const { startStudioServer, requireTempUserCards } = require("./studio_server.cjs");
let BASE = (process.argv[2] || "").replace(/\/$/, "");
let OWN = null;          // the private server this script started (none when a URL was given)
let USER_CARDS = "";     // the saved-card folder the server writes to (temp)
const SHOTS = path.join(ROOT, ".workers", "studio", "shots", "web");
const SHOTS_P2 = path.join(ROOT, ".workers", "studio", "shots", "p2");
const VIEW = { width: 1600, height: 1000 };

function loadPlaywright() {
  const tries = [
    path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright-core"),
    path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright"),
    "playwright-core",
    "playwright",
  ];
  for (const t of tries) {
    try { return require(t); } catch { /* next */ }
  }
  console.error("playwright-core not found: npm i playwright-core under .workers/studio/node");
  process.exit(2);
}

const results = [];
function check(name, ok, detail = "") {
  results.push({ name, ok: !!ok, detail });
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}${detail ? `  (${detail})` : ""}`);
}

async function waitRunDone(page, timeout = 120000) {
  await page.waitForFunction(() => {
    const s = document.querySelector(".status-word");
    return s && ["done", "error", "cancelled"].includes(s.dataset.state);
  }, null, { timeout });
  return page.$eval(".status-word", (el) => el.dataset.state);
}

async function settle(page, ms = 700) {
  await page.waitForTimeout(ms);
  await page.evaluate(() => document.fonts.ready);
}

async function kpiTagAudit(page) {
  return page.$$eval(".readouts .readout", (els) => els.filter((e) => !e.hidden).map((e) => ({
    key: e.dataset.readout,
    state: e.dataset.state,
    hasValue: e.dataset.state === "value",
    hasTag: !!e.querySelector(".ro-tag .tag"),
  })));
}

/** F3: no dotted Plotly line anywhere; every dashed ([A]/[E]) line >= 2.5px; no zero-width hover proxy (H20). */
async function traceAudit(page, scope = "body") {
  return page.evaluate((sel) => {
    const out = { traces: 0, dotted: [], thinDash: [], zeroWidth: [], widths: [] };
    for (const gd of document.querySelectorAll(`${sel} .chart-plot`)) {
      for (const t of gd.data || []) {
        out.traces++;
        const ln = t.line || {};
        if (ln.dash === "dot") out.dotted.push(t.name);
        if ((ln.dash === "dash" || ln.dash === "longdash") && !(ln.width >= 2.5)) out.thinDash.push(`${t.name}:${ln.width}`);
        if (t.mode && t.mode.includes("lines") && ln.width === 0) out.zeroWidth.push(t.name);
        if (t.mode === "lines" && ln.width) out.widths.push(ln.width);
      }
      for (const sh of gd.layout?.shapes || []) if (sh.line?.dash === "dot") out.dotted.push(`shape ${sh.type}`);
    }
    return out;
  }, scope);
}

/** F4: no numeral anywhere is set in B612 Mono (one numeral face: Barlow tabular). */
async function b612Audit(page) {
  return page.evaluate(() => {
    const bad = [];
    const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
    let n;
    while ((n = walker.nextNode())) {
      if (!/\d/.test(n.nodeValue)) continue;
      const el = n.parentElement;
      if (!el || !el.getClientRects().length) continue;
      const ff = getComputedStyle(el).fontFamily || "";
      if (/B612/i.test(ff.split(",")[0])) bad.push(`${el.className?.baseVal ?? el.className ?? el.tagName}:${n.nodeValue.trim().slice(0, 20)}`);
    }
    return bad;
  });
}

(async () => {
  if (BASE) USER_CARDS = requireTempUserCards(ROOT);
  else {
    OWN = await startStudioServer(ROOT, 8791);
    BASE = OWN.base;
    USER_CARDS = OWN.userCards;
  }
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch();
  const context = await browser.newContext({ viewport: VIEW, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const consoleErrors = [];
  page.on("console", (m) => { if (m.type() === "error") consoleErrors.push(m.text()); });
  page.on("pageerror", (e) => consoleErrors.push(`pageerror: ${e.message}`));
  page.on("requestfailed", (r) => {
    const u = r.url();
    // EventSource streams are closed by the client once a job finishes; not an error.
    if (!u.includes("/events")) consoleErrors.push(`requestfailed: ${u} ${r.failure()?.errorText}`);
  });

  const requests = [];
  page.on("request", (r) => requests.push(`${r.method()} ${r.url().replace(BASE, "")}`));
  let gates = null;

  try {
    // ---------------------------------------------------------------- landing route + gates
    await page.goto(`${BASE}/`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => location.hash === "#/overview", null, { timeout: 15000 }).catch(() => {});
    check("landing: an empty hash opens #/overview", (await page.evaluate(() => location.hash)) === "#/overview");
    gates = await page.evaluate(async () => (await fetch("/api/gates")).json()).catch(() => null);
    check("gates: GET /api/gates serves the g2 ceiling and the flux floor", gates && typeof gates.g2_ceiling?.value === "number" && typeof gates.flux_floor?.value === "number", JSON.stringify(gates)?.slice(0, 120));

    // ---------------------------------------------------------------- staged: point
    await page.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".designer:not(.is-loading)", { timeout: 30000 });
    let st = await waitRunDone(page);
    check("staged auto run finishes", st === "done", st);
    // explicit point run
    await page.click('.segmented [data-value="point"]');
    await page.click(".run-key");
    st = await waitRunDone(page);
    check("staged point run done", st === "done", st);
    await page.waitForSelector(".chart-hero .chart-plot .main-svg", { timeout: 20000 });
    let audit = await kpiTagAudit(page);
    check("staged point: every KPI with a value carries a tag chip", audit.length >= 5 && audit.every((a) => !a.hasValue || a.hasTag), JSON.stringify(audit.map((a) => `${a.key}:${a.state}:${a.hasTag}`)));
    check("staged point: g2_op readout has a value", audit.find((a) => a.key === "g2_op")?.hasValue);
    const runId = await page.$eval("#run-chip .id-v", (e) => e.textContent);
    check("run ID chip shows FS-<card>-<hash6>", /^FS-[\w-]+-[0-9a-f]{6}$/.test(runId), runId);
    const ceilTxt = await page.$$eval(".chart-hero .annotation-text", (e) => e.map((x) => x.textContent).join(" | "));
    check("H8: hero ceiling line is the /api/gates value, labelled with it", gates && ceilTxt.includes(`= ${gates.g2_ceiling.value} ceiling`), ceilTxt);
    const grat = await page.$eval(".chart-hero .chart-plot", (gd) => !!(gd.layout?.yaxis?.minor?.showgrid && gd.layout?.xaxis?.minor?.showgrid));
    check("ceiling: g2(T) screen carries a minor graticule", grat);
    // H19: a cancel issued before the POST returns still DELETEs the job once its id arrives
    const T = (120 + Math.random() * 100).toFixed(3);
    const cancelled = await page.evaluate(async (t) => {
      const api = await import("/js/api.js");
      const j = api.runJob({ card: "edge-inp-gainp-design", mode: "point", T_grid: [Number(t), Number(t) + 1, Number(t) + 2] });
      j.cancel();
      let kind = null;
      try { await j.promise; } catch (e) { kind = e.body?.error?.kind || e.message; }
      await new Promise((r) => setTimeout(r, 400));
      return { id: j.jobId, kind };
    }, T);
    check("H19: cancel before the job id -> DELETE sent when the id arrives", cancelled.id && requests.includes(`DELETE /api/jobs/${cancelled.id}`) && cancelled.kind === "cancelled", JSON.stringify(cancelled));
    // open a mount sheet to show the in-place editor in the hero shot
    await page.click('.mount[data-block="dot"]');
    await page.waitForSelector(".sheet .field-row");
    const rows = await page.$$eval(".sheet .field-row", (e) => e.length);
    check("Dot block sheet lists META fields", rows >= 5, `${rows} rows`);
    await page.click('.sheet .field-row[data-path="dot.delta_xx"] .fr-tools .tag');
    await page.waitForSelector("#inspector:not([hidden]) .insp-title");
    const inspTag = await page.$eval("#inspector .insp-tag strong", (e) => e.textContent);
    check("inspector shows the META tag of the selected field", inspTag.startsWith("[A]"), inspTag);
    await settle(page, 400);
    await page.screenshot({ path: path.join(SHOTS, "designer-sheet.png") });
    await page.keyboard.press("Escape");
    await page.click('.mount[data-block="dot"]'); // close again
    await settle(page);
    const bd = await b612Audit(page);
    check("F4: no B612 numerals in the Designer", bd.length === 0, bd.slice(0, 5).join(", "));
    await page.screenshot({ path: path.join(SHOTS, "designer-dark.png") });

    // pin A
    await page.click(".pin-main");
    await settle(page, 300);

    // ---------------------------------------------------------------- staged: envelope
    await page.click('.segmented [data-value="envelope"]');
    await page.click(".run-key");
    st = await waitRunDone(page, 180000);
    check("staged envelope run done", st === "done", st);
    const tcText = await page.$eval('[data-readout="T_c"] .ro-value', (e) => e.textContent);
    check("envelope: T_c rendered as an interval", tcText.includes("–"), tcText);
    const g2Text = await page.$eval('[data-readout="g2_op"] .ro-value', (e) => e.textContent);
    check("envelope: g2_op rendered as an interval", g2Text.includes("–"), g2Text);
    const legend = await page.$$eval(".chart-hero .legendtext", (e) => e.map((x) => x.textContent).join(" | "));
    check("envelope: mid curve labelled as not a prediction", legend.includes("all ranges at midpoint (not a prediction)"), legend);
    audit = await kpiTagAudit(page);
    check("envelope: every KPI with a value carries a tag chip", audit.every((a) => !a.hasValue || a.hasTag));
    const ta = await traceAudit(page, ".designer");
    check("F3: no dotted line; [A]/[E] dashed lines >= 2.5px (designer envelope)", ta.traces > 0 && !ta.dotted.length && !ta.thinDash.length, JSON.stringify({ d: ta.dotted, t: ta.thinDash }));
    check("H20: no zero-width (averaged) hover series in the envelope figure", ta.zeroWidth.length === 0, ta.zeroWidth.join(","));
    await settle(page);
    await page.screenshot({ path: path.join(SHOTS, "designer-envelope.png") });
    // pin B via the split menu
    await page.click(".pin-more");
    await page.click('.pin-opt[data-slot="B"]');
    await settle(page, 300);

    // ---------------------------------------------------------------- theme toggle
    await page.click("#theme-btn");
    await page.waitForFunction(() => document.documentElement.dataset.theme === "light");
    await settle(page, 900);
    check("theme toggles to light", (await page.evaluate(() => document.documentElement.dataset.theme)) === "light");
    const bg = await page.$eval(".chart-hero .main-svg", (e) => getComputedStyle(e).backgroundColor);
    check("charts re-templated on theme switch", bg && !bg.includes("21, 23, 28"), bg);
    await page.screenshot({ path: path.join(SHOTS, "designer-light.png") });
    await page.click("#theme-btn");
    await page.waitForFunction(() => document.documentElement.dataset.theme === "dark");

    // ---------------------------------------------------------------- edge: headline, 1-point grid
    await page.goto(`${BASE}/#/design/edge-inp-gainp-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("edge-inp-gainp"), null, { timeout: 30000 });
    st = await waitRunDone(page, 120000);
    check("edge auto run finishes", st === "done", st);
    await page.click('.segmented [data-value="headline"]');
    const opBtn = await page.$('.segmented [data-value="op"]');
    if (opBtn) await opBtn.click();
    await page.click(".run-key");
    st = await waitRunDone(page, 120000);
    check("edge headline run (1-point T grid) done", st === "done", st);
    const mode = await page.$eval(".status-mode", (e) => e.textContent).catch(() => "");
    check("edge result is the headline model", mode.includes("headline"), mode);
    audit = await kpiTagAudit(page);
    check("edge headline: every KPI with a value carries a tag chip", audit.every((a) => !a.hasValue || a.hasTag), JSON.stringify(audit.map((a) => `${a.key}:${a.state}:${a.hasTag}`)));
    // F2: one continuous beam behind the mounts, emission point -> Detection, >= 4px, visible at rest
    const beam = await page.evaluate(() => {
      const chain = document.querySelector(".chain").getBoundingClientRect();
      const b = document.querySelectorAll(".beam-svg .beam");
      const el = b[0];
      const m = /M\s*([\d.]+)\s+([\d.]+)\s+L\s*([\d.]+)\s+([\d.]+)/.exec(el?.getAttribute("d") || "");
      const c = (id) => { const r = document.querySelector(`.mount[data-block="${id}"]`).getBoundingClientRect(); return r.left - chain.left + r.width / 2; };
      const ap = document.querySelector('.mount[data-block="filter"] .mount-aperture').getBoundingClientRect();
      const mount = document.querySelector('.mount[data-block="filter"]');
      return {
        n: b.length, sw: el ? parseFloat(getComputedStyle(el).strokeWidth) : 0, op: el ? parseFloat(el.style.opacity) : 0,
        x0: m ? +m[1] : null, x1: m ? +m[3] : null, y: m ? +m[2] : null, dotX: c("dot"), detX: c("detection"),
        apY: ap.top - chain.top + ap.height / 2,
        behind: Number(getComputedStyle(document.querySelector(".beam-svg")).zIndex) < Number(getComputedStyle(document.querySelector(".mount-cell")).zIndex),
        windowed: (getComputedStyle(mount).backgroundImage.match(/gradient/g) || []).length >= 2,
        cap: document.querySelector(".beam-cap")?.textContent || "",
      };
    });
    check("F2: exactly one beam path, emission point (Dot) to Detection", beam.n === 1 && Math.abs(beam.x0 - beam.dotX) < 2 && Math.abs(beam.x1 - beam.detX) < 2, JSON.stringify(beam));
    check("F2: beam >= 4px, visible at rest (opacity >= 0.45), dims with log flux", beam.sw >= 4 && beam.op >= 0.45 && /log₁₀ collected flux/.test(beam.cap), `${beam.sw}px op ${beam.op}`);
    check("F2: beam runs behind the mounts through their windowed apertures", beam.behind && beam.windowed && Math.abs(beam.y - beam.apY) < 2);
    const posts = await page.$$eval(".mount-cell > .post", (ps) => ps.map((p) => { const r = p.getBoundingClientRect(); return ((r.left + r.width / 2) % 25 + 25) % 25; }));
    check("ceiling: mounts snapped to the 25px hole pitch (posts on hole centres)", posts.length === 6 && posts.every((x) => Math.abs(x - 12.5) < 0.6), posts.map((x) => x.toFixed(1)).join(","));
    const plq = await page.$eval(".plaque", (e) => ({ t: e.textContent, s: e.dataset.status }));
    check("I6: edge card plaque is the backend-matched rt_edge VERDICT line 6", /VERDICT line 6/.test(plq.t), plq.t.slice(0, 90));
    const engraved = await page.$eval(".mount .engrave", (e) => getComputedStyle(e).textShadow);
    check("ceiling: engraved labels use an inset (two-sided) text shadow", (engraved.match(/rgba?\(/g) || []).length >= 2, engraved);
    const bezel = await page.$eval(".readout", (e) => getComputedStyle(e).boxShadow);
    check("ceiling: shelf readouts sit in inset windows (bezel)", /inset/.test(bezel), bezel.slice(0, 60));
    await settle(page, 900);
    await page.screenshot({ path: path.join(SHOTS, "designer-edge-headline.png") });

    // signature interaction: heatsink slider re-runs a single-T point
    const before = await page.$eval('[data-readout="g2_op"] .ro-value', (e) => e.textContent);
    await page.$eval("#hs-slider", (el) => { el.value = "250"; el.dispatchEvent(new Event("input", { bubbles: true })); });
    await page.waitForFunction(() => document.querySelector(".status-op")?.textContent.includes("250"), null, { timeout: 60000 });
    const after = await page.$eval('[data-readout="g2_op"] .ro-value', (e) => e.textContent);
    const changed = await page.$eval('[data-readout="g2_op"]', (e) => e.classList.contains("is-changed"));
    check("heatsink drag re-runs the operating point", after !== before || changed, `${before} -> ${after}`);

    // static label on the plain point run of the edge card
    await page.click('.segmented [data-value="point"]');
    const op2 = await page.$('.segmented [data-value="op"]');
    if (op2) await op2.click();
    await page.click(".run-key");
    st = await waitRunDone(page, 120000);
    const qual = await page.$eval('[data-readout="g2_op"] .ro-qual', (e) => e.textContent);
    check("edge static point run carries 'static (non-headline)' in the KPI title", qual === "static (non-headline)", qual);
    // H21: the headline-model g2 at the same T_hs beside the static tiles
    await page.waitForSelector('[data-readout="g2_head"][data-state="value"]', { timeout: 60000 }).catch(() => {});
    const head = await page.$eval('[data-readout="g2_head"]', (e) => ({ hidden: e.hidden, state: e.dataset.state, tag: !!e.querySelector(".ro-tag .tag"), cap: e.querySelector(".ro-caption")?.textContent || "" }));
    check("H21: point run on a static edge card shows the headline g2 beside the static tiles", !head.hidden && head.state === "value" && head.tag && /headline model at T_hs/.test(head.cap), JSON.stringify(head));

    // ---------------------------------------------------------------- nitride: tags, b_res floor, verdict mapping
    await page.goto(`${BASE}/#/design/nitride-cavity-set-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("nitride-cavity-set"), null, { timeout: 30000 });
    st = await waitRunDone(page, 120000);
    check("nitride SET auto run finishes", st === "done", st);
    audit = await kpiTagAudit(page);
    const g2n = audit.find((a) => a.key === "g2_op");
    check("H3: nitride g2 readout renders a tagged value (backend tag_chain present)", g2n?.hasValue && g2n?.hasTag, JSON.stringify(g2n));
    await page.waitForSelector(".chart-hero .main-svg", { timeout: 20000 }).catch(() => {});
    const wh = await page.$eval(".chart-hero", (e) => ({ w: e.dataset.withheld || "0", svg: !!e.querySelector(".main-svg:not([hidden])") && !e.querySelector(".chart-plot").hidden, empty: e.querySelector(".chart-empty:not([hidden])")?.textContent || "" }));
    // Nitride results carry no curves.T_hs (backend gap, reported): the card must then say so, never "Not run yet".
    check("H3: nitride hero chart draws, or states that no T_hs axis was returned (never withheld)", wh.w === "0" && (wh.svg || /No temperature axis returned/.test(wh.empty)), JSON.stringify(wh));
    const refused = await page.evaluate(async () => {
      const { buildG2T } = await import("/js/charts/g2T.js");
      const r = buildG2T([{ name: "untagged", result: { tag_chain: null, curves: { T_hs: [200, 300], g2: [0.1, 0.2] }, scalars: {} } }]);
      return { n: r?.data?.length ?? -1, withheld: r?.withheld || [] };
    });
    check("H3: chart builder refuses a series with a null tag (withheld, reason listed)", refused.n === 0 && refused.withheld.includes("untagged"), JSON.stringify(refused));
    const mark = await page.$eval('[data-readout="g2_op"] .ro-marker', (e) => ({ hidden: e.hidden, t: e.textContent })).catch(() => ({ hidden: true, t: "" }));
    check("I2: SET g2 tile carries the b_res structural-floor marker", !mark.hidden && /structural floor set by b_res/.test(mark.t), mark.t);
    const plN = await page.$eval(".plaque", (e) => e.textContent);
    check("H2/I6: nitride-cavity-set plaque is VERDICT line 7 (deterministic_pair), not line 5", /VERDICT line 7/.test(plN) && !/line 5\b/.test(plN), plN.slice(0, 100));
    await page.screenshot({ path: path.join(SHOTS, "designer-nitride-set.png") });
    await page.goto(`${BASE}/#/design/nitride-nonpolar-set-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("nitride-nonpolar-set"), null, { timeout: 30000 });
    await page.waitForSelector(".plaque", { timeout: 30000 });
    const plP = await page.$eval(".plaque", (e) => ({ t: e.textContent, s: e.dataset.status }));
    check("H2/I6: nonpolar SET plaque maps to a_plane lines only (ambiguous shown, not picked)", plP.s === "ambiguous" && /line 11/.test(plP.t) && /line 12/.test(plP.t) && !/line 5\b/.test(plP.t), `${plP.s}: ${plP.t.slice(0, 120)}`);

    // ---------------------------------------------------------------- N1 / N2 nanowire designer
    await page.goto(`${BASE}/#/design/nitride-nanowire-vertical-pulse-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("nitride-nanowire-vertical-pulse"), null, { timeout: 30000 });
    await page.waitForSelector('.runbar button[data-value="envelope"]', { timeout: 30000 });
    const envSeg = await page.$eval('.runbar button[data-value="envelope"]', (b) => ({ dis: b.disabled, title: b.title, sub: b.querySelector(".seg-sub")?.textContent || "" }));
    check("N1: nanowire card disables the Envelope segment with the backend's reason",
      envSeg.dis && /not available for the nanowire tier/.test(envSeg.title) && /g2_op\/T_j/.test(envSeg.title) && /nanowire/.test(envSeg.sub), JSON.stringify(envSeg).slice(0, 160));
    await page.click('.runbar button[data-value="full"]');
    await page.waitForTimeout(200);
    const etaTxt = await page.$eval(".run-eta", (e) => e.textContent);
    const etaS = (() => { const m = /estimated (\d+(?:\.\d+)?) (min|s)\b/.exec(etaTxt); return m ? Number(m[1]) * (m[2] === "min" ? 60 : 1) : NaN; })();
    // jobs.ETA_TABLE nanowire (1.0 s single T, 295 s at 120 points, re-measured 2026-10-08 after the vectorized injector):
    // 16-point grid = 1.0 + 294 * 15/119 = 38.1 s, plus the 1.0 s single-T partial = 39 s
    check("N2: nanowire full-curve estimate is the 16-point display grid from the re-measured ETA table (about 39 s: 30-50 s)",
      /16-point T_hs grid/.test(etaTxt) && etaS >= 30 && etaS <= 50, `${etaTxt} -> ${etaS}`);
    await page.click('.runbar button[data-value="op"]');

    // ---------------------------------------------------------------- compare
    await page.goto(`${BASE}/#/compare`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".slot-plates");
    await page.waitForSelector(".chart-hero .main-svg", { timeout: 20000 });
    const nSlots = await page.$$eval(".slot-plate:not(.is-empty)", (e) => e.length);
    check("compare shows pinned slots A and B", nSlots >= 2, `${nSlots}`);
    const deltaRows = await page.$$eval(".delta-table tbody tr", (e) => e.length);
    check("compare delta table has rows", deltaRows >= 3, `${deltaRows}`);
    const cells = await page.$$eval(".delta-table td[data-metric]", (tds) => tds.map((t) => ({
      m: t.dataset.metric, slot: t.dataset.slot, v: !!t.querySelector(".cell-v"), tag: !!t.querySelector(".cell-line .tag"), mid: t.dataset.mid === "1",
      w: t.querySelector(".delta-w")?.textContent || "", conv: t.querySelector(".cell-conv")?.textContent || "", ref: !!t.querySelector(".delta.is-ref"),
    })));
    check("F5/H6: every delta-table value carries a tag chip", cells.filter((c) => c.v).length >= 3 && cells.filter((c) => c.v).every((c) => c.tag), JSON.stringify(cells.filter((c) => c.v && !c.tag)));
    const mids = cells.filter((c) => c.mid);
    check("F5/H6: envelope mid-design values are qualified and never compared better/worse", mids.length >= 1 && mids.every((c) => c.ref || c.w === "midpoint, not compared"), JSON.stringify(mids));
    const midMetrics = new Set(mids.map((c) => c.m));
    check("F5/H6: no better/worse word on a row where either side is mid-design", cells.filter((c) => midMetrics.has(c.m)).every((c) => !/better|worse/i.test(c.w)));
    check("F5/H6: brightness carries its convention caption", cells.some((c) => c.m === "brightness_per_pulse" && c.conv.length > 3), cells.filter((c) => c.m === "brightness_per_pulse").map((c) => c.conv).join(" | "));
    const tc = await traceAudit(page, ".compare");
    check("F3: compare overlay has no dotted line; dashed lines >= 2.5px", !tc.dotted.length && !tc.thinDash.length, JSON.stringify(tc.dotted.concat(tc.thinDash)));
    const bc = await b612Audit(page);
    check("F4: no B612 numerals on Compare", bc.length === 0, bc.slice(0, 5).join(", "));
    await settle(page, 900);
    await page.screenshot({ path: path.join(SHOTS, "compare.png") });

    // ---------------------------------------------------------------- present
    await page.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".designer:not(.is-loading)");
    await waitRunDone(page);
    await page.click("#present-btn");
    await page.waitForFunction(() => document.documentElement.hasAttribute("data-present"));
    await settle(page, 900);
    const railHidden = await page.$eval("#rail", (e) => getComputedStyle(e).display === "none");
    const pBeam = await page.$eval(".beam-svg .beam", (e) => parseFloat(getComputedStyle(e).strokeWidth));
    const pt = await traceAudit(page, ".designer");
    check("F3: Present widens the beam and the model curves", pBeam >= 6 && pt.widths.some((w) => w >= 2.8), `beam ${pBeam}px, widths ${pt.widths.join(",")}`);
    check("present: light theme, rail hidden", (await page.evaluate(() => document.documentElement.dataset.theme)) === "light" && railHidden);
    await page.screenshot({ path: path.join(SHOTS, "present.png") });
    await page.keyboard.press("Escape");
    await page.waitForFunction(() => !document.documentElement.hasAttribute("data-present"));
    check("present: Esc exits and restores the theme", (await page.evaluate(() => document.documentElement.dataset.theme)) === "dark");

    // ---------------------------------------------------------------- p2a: presets -> unsaved design -> save as card -> reopen
    fs.mkdirSync(SHOTS_P2, { recursive: true });
    await page.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("staged-device"), null, { timeout: 30000 });
    await waitRunDone(page);
    await page.click(".bh-presets");
    await page.waitForSelector(".preset-grid .sro", { timeout: 20000 });
    const PICK = [["dot", "piezo-variant"], ["template", "GaAs-Si"], ["cavity", "planar-lambda"], ["drive", "pulsed-electrical"], ["injection", "rti-quiet"]];
    for (const [ax, key] of PICK) await page.click(`.po[data-axis="${ax}"][data-key="${key}"]`);
    await page.waitForFunction(() => /piezo-variant \/ GaAs-Si \/ planar-lambda \/ pulsed-electrical \/ rti-quiet/.test(document.querySelector(".preset-combo")?.textContent || "")
      && !document.querySelector(".preset-summary.is-pending"), null, { timeout: 15000 });
    const pa = await page.evaluate(() => {
      const axes = [...document.querySelectorAll(".preset-axis")].map((f) => f.dataset.axis);
      const opts = [...document.querySelectorAll(".po")].map((o) => ({ k: `${o.dataset.axis}:${o.dataset.key}`, tag: !!o.querySelector(".tag"), desc: o.querySelector(".po-desc")?.textContent || "" }));
      const sum = [...document.querySelectorAll(".preset-grid .sro")].map((s) => ({ k: s.dataset.key, st: s.dataset.state, tag: !!s.querySelector(".tag") }));
      return { axes, opts, sum };
    });
    check("p2a presets sheet: five axes (dot, template, cavity, drive, injection), every option with a one-line description and its provenance chip (a no-change option has none)",
      pa.axes.join(",") === "dot,template,cavity,drive,injection" && pa.opts.length === 15 && pa.opts.every((o) => o.desc.length > 3 && (o.tag || o.k === "injection:standard")),
      JSON.stringify(pa.opts.filter((o) => !o.tag).map((o) => o.k)));
    check("p2a presets: live summary of the composed design, every key field a tagged value",
      pa.sum.length >= 8 && pa.sum.every((s) => s.st === "value" && s.tag) && pa.sum.some((s) => s.k === "drive.F_p"), `${pa.sum.length} fields`);
    await settle(page, 500);
    await page.screenshot({ path: path.join(SHOTS_P2, "designer-presets.png") });
    await page.click(".preset-open");
    await page.waitForFunction(() => location.hash === "#/design/unsaved", null, { timeout: 15000 });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("unsaved:"), null, { timeout: 15000 });
    st = await waitRunDone(page);
    const draft = await page.evaluate(() => ({
      title: document.querySelector(".board-title").textContent,
      run: document.querySelector("#run-chip .id-v").textContent,
    }));
    const chipTxt = await page.$eval("#unsaved-chip", (e) => (e.hidden ? "" : e.textContent)).catch(() => "");
    check("p2a Open in Designer: unsaved design titled 'unsaved: <combo>', runs (FS-unsaved-...), context bar says unsaved",
      st === "done" && draft.title === "unsaved: piezo-variant / GaAs-Si / planar-lambda / pulsed-electrical / rti-quiet"
      && /^FS-unsaved-[0-9a-f]{6}$/.test(draft.run) && /unsaved/.test(chipTxt), `${draft.title} | ${draft.run} | ${chipTxt}`);
    await page.evaluate(() => { location.hash = "#/design/staged-device-design"; });
    await page.waitForSelector(".leave-dialog", { timeout: 5000 });
    const nDlg = await page.$$eval(".leave-dialog", (e) => e.length);
    await page.click(".leave-stay");
    await page.waitForTimeout(200);
    const stayHash = await page.evaluate(() => location.hash);
    check("p2a leaving an unsaved design asks once (one dialog); Stay keeps it", nDlg === 1 && stayHash === "#/design/unsaved"
      && (await page.$$eval(".leave-dialog", (e) => e.length)) === 0, `${nDlg} dialog(s), ${stayHash}`);
    await page.click(".bh-save");
    await page.fill(".save-name", "staged-device-design");
    const shippedBlocked = await page.$eval(".save-go", (b) => b.disabled);
    await page.fill(".save-name", "Bad Name");
    const slugBlocked = await page.$eval(".save-go", (b) => b.disabled);
    // (from Node, not the page: a 4xx in the page would log a console error)
    const putShipped = (await page.request.put(`${BASE}/api/cards/staged-device-design`, { headers: { "X-FSIM-Studio": "1" }, data: { design: {} } })).status();
    check("p2a Save as card: shipped names and non-slugs blocked in the form; the server refuses a shipped name (409)",
      shippedBlocked && slugBlocked && putShipped === 409, `${shippedBlocked} ${slugBlocked} ${putShipped}`);
    await page.fill(".save-name", "p2a-verify-preset");
    await page.click(".save-go");
    await page.waitForFunction(() => location.hash === "#/design/p2a-verify-preset", null, { timeout: 15000 });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent === "p2a-verify-preset", null, { timeout: 15000 });
    st = await waitRunDone(page);
    const saved = await page.evaluate(() => ({
      picker: [...document.querySelectorAll("#card-picker option")].some((o) => o.value === "p2a-verify-preset" && o.selected),
      chip: document.querySelector("#unsaved-chip")?.hidden !== false,
    }));
    check("p2a saved card appears in the card picker at once (selected) and the unsaved indicator clears", saved.picker && saved.chip && st === "done", JSON.stringify(saved));
    check("p2a Save as card wrote into FSIM_STUDIO_USER_CARDS (a temp folder), not the user's cards/studio/",
      fs.existsSync(path.join(USER_CARDS, "p2a-verify-preset.yaml")) && path.resolve(USER_CARDS) !== path.resolve(ROOT, "cards", "studio"), USER_CARDS);
    await page.reload({ waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent === "p2a-verify-preset", null, { timeout: 30000 });
    st = await waitRunDone(page);
    const rt = await page.evaluate(async () => {
      const H = { "Content-Type": "application/json", "X-FSIM-Studio": "1" };
      const ap = await (await fetch("/api/presets/apply", { method: "POST", headers: H, body: JSON.stringify({ dot: "piezo-variant", template: "GaAs-Si", cavity: "planar-lambda", drive: "pulsed-electrical", injection: "rti-quiet" }) })).json();
      const card = await (await fetch("/api/cards/p2a-verify-preset")).json();
      const a = { ...ap.design, name: "p2a-verify-preset" };
      return { same: JSON.stringify(a) === JSON.stringify(card.design), saved: card.saved };
    });
    check("p2a reopen after reload: the saved card loads, runs, and equals the composed preset design", st === "done" && rt.same && rt.saved === true, JSON.stringify(rt));

    // ---------------------------------------------------------------- p2a: timeline honesty (deterministic, no server)
    const tlh = await page.evaluate(async () => {
      const { createTimeline } = await import("/js/ui/timeline.js");
      const seen = [];
      const tl = createTimeline({ n: 6, values: [0, 1, 2, 3, 4, 5], label: "x", onChange: (p, w) => seen.push({ p, k: w.kind, held: w.held, next: w.next }) });
      document.body.append(tl.el);
      for (const i of [0, 1, 4, 5]) tl.setComputed(i);
      tl.setJobOver(true);
      tl.seek(1.5);
      const w15 = tl.where().kind;
      tl.seek(0.5);
      const w05 = tl.where().kind;
      tl.seek(0);
      seen.length = 0;
      tl.play();
      await new Promise((r) => setTimeout(r, 1400));
      tl.pause();
      tl.el.remove();
      const bad = seen.filter((s) => s.k === "interp" && (s.held === 1 || s.held === 2 || s.held === 3));
      const sat = seen.filter((s) => s.p > 1 + 1e-9 && s.p < 4 - 1e-9);
      return { w15, w05, n: seen.length, bad: bad.length, sat: sat.length, reached4: seen.some((s) => s.p >= 4) };
    });
    check("p2a timeline honesty: interp only between computed neighbours; a gap is 'gap'; playback skips the uncomputed frames 2-3 (never sits between 1 and 4)",
      tlh.w15 === "gap" && tlh.w05 === "interp" && tlh.bad === 0 && tlh.sat === 0 && tlh.reached4, JSON.stringify(tlh));

    // ---------------------------------------------------------------- p2a: response animation, staged (n = 48)
    async function animAudit(pg, n) {
      return pg.evaluate(async (N) => {
        const A = document.querySelector(".anim-sheet").__anim;
        const val = () => document.querySelector('[data-readout="g2_op"] .ro-value').textContent;
        const qual = () => document.querySelector('[data-readout="g2_op"] .ro-qual').textContent;
        const cells = [...document.querySelectorAll(".tl-cell")].map((c) => c.dataset.state);
        A.seek(0); const v0 = val();
        A.seek(N - 1); const v1 = val();
        A.seek(3); const qFrame = qual(); const capFrame = document.querySelector('[data-readout="g2_op"] .ro-caption').textContent;
        A.seek(3.5); const qInterp = qual(); const vMid = val();
        const marker = document.querySelector(".anim-marker");
        const mk = { hidden: marker?.hidden, kind: marker?.dataset.kind };
        // verdict-like text (labels, flags) comes from the NEAREST computed frame (rule 3) and so
        // changes only where the nearest frame changes (at the midpoint between two frames)
        const comp = A.states();
        const samples = [];
        for (let p = 0; p <= N - 1 + 1e-9; p += 0.25) {
          A.seek(p);
          const lab = document.querySelector(".anim-labels");
          samples.push({ p, held: lab.dataset.frame, text: lab.textContent });
        }
        let badSwitch = 0, badHeld = 0;
        for (let i = 1; i < samples.length; i++) {
          const a = samples[i - 1], b = samples[i];
          if (b.text !== a.text && Math.round(b.p) === Math.round(a.p)) badSwitch++;
          const nr = Math.min(Math.round(b.p), N - 1);
          if (comp[nr] === "computed" && b.held !== String(nr)) badHeld++;
        }
        // Performance API: JS time per playback frame
        performance.clearMeasures("fsim-anim-frame");
        A.seek(0);
        A.play();
        await new Promise((r) => setTimeout(r, 2600));
        A.pause();
        const d = performance.getEntriesByName("fsim-anim-frame").map((e) => e.duration).sort((x, y) => x - y);
        const mean = d.reduce((x, y) => x + y, 0) / Math.max(1, d.length);
        return { cells, v0, v1, vMid, qFrame, qInterp, capFrame, mk, badSwitch, badHeld, nSamples: samples.length,
          perf: { n: d.length, mean, p95: d[Math.floor(d.length * 0.95)] ?? null, max: d[d.length - 1] ?? null } };
      }, n);
    }
    await page.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("staged-device"), null, { timeout: 30000 });
    await waitRunDone(page);
    await page.click(".hs-anim");
    await page.waitForSelector(".anim-sheet");
    const setup = await page.evaluate(() => ({ n: document.querySelector(".anim-n").value, eta: document.querySelector(".anim-eta").textContent }));
    check("p2a staged animate: defaults to 48 frames with the ETA shown before start", setup.n === "48" && /48 frames/.test(setup.eta) && /about/.test(setup.eta), setup.eta.slice(0, 90));
    await page.click(".anim-go");
    await page.waitForFunction(() => document.querySelector(".anim-sheet")?.dataset.phase === "ready", null, { timeout: 120000 });
    const sa = await animAudit(page, 48);
    check("p2a staged animate: the scrubber shows all 48 frames computed", sa.cells.length === 48 && sa.cells.every((c) => c === "computed"), sa.cells.filter((c) => c !== "computed").length + " not computed");
    check("p2a staged animate: readouts change with the playhead", sa.v0 !== sa.v1 && sa.vMid !== "", `${sa.v0} -> ${sa.v1}`);
    check("p2a staged animate: 'interp.' between computed frames, not on a computed frame; the chart marker says which",
      /interp\./.test(sa.qInterp) && !/interp\./.test(sa.qFrame) && /computed frame 4 of 48/.test(sa.capFrame) && sa.mk.hidden === false && sa.mk.kind === "interp", `${sa.qFrame} | ${sa.qInterp} | ${sa.capFrame}`);
    check("p2a staged animate: labels/verdict text switches only when the NEAREST computed frame changes and is held from it (rule 3)",
      sa.badSwitch === 0 && sa.badHeld === 0 && sa.nSamples > 180, `${sa.badSwitch} bad switches, ${sa.badHeld} bad holds, ${sa.nSamples} samples`);
    check("p2a staged animate: JS per playback frame < 2 ms (Performance API: mean and p95)", sa.perf.n > 30 && sa.perf.mean < 2 && sa.perf.p95 < 2, JSON.stringify(sa.perf));
    await page.evaluate(() => { document.querySelector(".anim-sheet").__anim.seek(20.5); document.querySelector(".anim-player").scrollIntoView({ block: "end" }); });
    await settle(page, 500);
    await page.screenshot({ path: path.join(SHOTS_P2, "designer-animate-staged.png") });

    // ---------------------------------------------------------------- p2d: animation honesty rules (staged, 48 frames)
    const hon = await page.evaluate(async () => {
      const A = document.querySelector(".anim-sheet").__anim;
      const ro = (k) => {
        const e = document.querySelector(`[data-readout="${k}"]`);
        return { v: e.querySelector(".ro-value").textContent, q: e.querySelector(".ro-qual").textContent, c: e.querySelector(".ro-caption").textContent, tag: e.querySelector(".ro-tag .tag")?.dataset.tag || "" };
      };
      const out = {};
      // rule 2: a computed frame wears the run's tag; a blend wears the dashed [A] form
      const tags0 = [A.frames[5].tag_chain, A.frames[6].tag_chain];
      A._patch(5, {}, null, "E"); A._patch(6, {}, null, "E");
      A.seek(5); out.tagFrame = ro("g2_op").tag;
      A.seek(5.5); out.tagInterp = ro("g2_op").tag; out.tcTagInterp = ro("T_c").tag;
      A.seek(6); out.tagBack = ro("g2_op").tag;
      A._patch(5, {}, null, tags0[0]); A._patch(6, {}, null, tags0[1]);
      out.tcValFrame = (A.seek(5), ro("T_c").v);
      // rule 1: invalid g2 at frame 10 -> no blend on either side of it; the frame labels differ at 20|21
      A._patch(10, { g2_op_valid: false });
      out.refuse = [A.refusal(9, 10), A.refusal(10, 11), A.refusal(11, 12)];
      A.seek(9.5); out.snapLo = { where: A.where().kind, ...ro("g2_op") };
      { // review finding 10: on a snap the note and the scrubber text show the value of the frame the readouts come from
        const w = A.where();
        out.snapV = { kind: w.kind, near: w.near, vn: A.values[w.near], va: (A.values[9] + A.values[10]) / 2, cap: ro("g2_op").c,
          pos: document.querySelector(".tl-pos").textContent };
      }
      A.seek(10.5); out.snapHi = { where: A.where().kind, ...ro("g2_op") };
      A.seek(11.5); out.okNear = { where: A.where().kind, ...ro("g2_op") };
      const lab20 = A.frames[20].labels.slice();
      A._patch(21, {}, lab20.concat(["synthetic label"]));
      out.labelRefusal = A.refusal(20, 21);
      A.seek(20.25); out.labSnap = { where: A.where().kind, txt: document.querySelector(".anim-labels").textContent };
      A.seek(20.75); out.labSnap2 = document.querySelector(".anim-labels").textContent;
      // rule 3/4: T_c is never blended: it is the nearest computed frame's own value; the change is "between frames"
      const tc0 = A.frames[14].scalars.T_c, tc1 = A.frames[15].scalars.T_c;
      A._patch(14, { T_c: 100 });
      const seen = new Set();
      let cross = "";
      for (let p = 14; p <= 15.0001; p += 0.125) { A.seek(p); seen.add(ro("T_c").v); if (p === 14.5) cross = document.querySelector(".anim-cross").textContent; }
      out.tcSeen = [...seen]; out.tcCross = cross; out.tcOrig = [tc0, tc1];
      A._patch(14, { T_c: tc0 });
      // rule 5: runaway at frame 30 ends the range
      A._patch(30, { runaway: true });
      out.runIdx = A.runawayIdx;
      A.seek(29.5); out.runBefore = { where: A.where().kind, q: ro("g2_op").q };
      A.seek(30); out.runAt = { q: ro("g2_op").q, c: ro("g2_op").c, qT: ro("T_j_op").q, v: ro("g2_op").v };
      A.seek(40.5); out.runAfter = { q: ro("g2_op").q, c: ro("g2_op").c };
      out.beyond = document.querySelectorAll('.tl-cell[data-beyond="1"]').length;
      A.seek(28); A.play();
      await new Promise((r) => setTimeout(r, 1800));
      out.runPlay = { pos: A.pos, playing: document.querySelector(".tl-play").getAttribute("aria-pressed") };
      A.pause();
      A._patch(30, { runaway: false }); A._patch(10, { g2_op_valid: true });
      out.restored = A.runawayIdx == null && A.refusal(9, 10) === null;
      A._patch(21, {}, lab20);
      return out;
    });
    check("p2d rule 2: on a computed frame the chip is the run's tag ([E] here), between two frames the dashed [A] form, back to the run's tag on the next frame; T_c (never blended) keeps the frame's tag",
      hon.tagFrame === "E" && hon.tagInterp === "A" && hon.tagBack === "E" && hon.tcTagInterp === "E", `${hon.tagFrame} -> ${hon.tagInterp} -> ${hon.tagBack}; T_c ${hon.tcTagInterp}`);
    {
      const num = (t, re) => { const m = re.exec(t); return m ? Number(m[1]) : NaN; };
      const sv = hon.snapV;
      const capN = num(sv.cap, /= (-?\d+(?:\.\d+)?)/), posN = num(sv.pos.replace(/^\S+\s/, ""), /(-?\d+(?:\.\d+)?)/);
      const near = (a, b) => Math.abs(a - b) <= 2e-3 * Math.abs(b) + 1e-9;
      check("review 10: on a snap the readout note and the scrubber text show the value of the frame snapped to, not the interpolated parameter value",
        sv.kind === "snap" && near(capN, sv.vn) && near(posN, sv.vn) && !near(capN, sv.va), `cap ${capN} pos ${posN} frame ${sv.vn} interp ${sv.va}`);
    }
    check("p2d rule 1: no blend next to a frame with g2_op_valid false; the playhead snaps ('snap'), says 'not interpolated' and loses the 'interp.' mark; a valid neighbour pair still blends",
      hon.refuse[0] && hon.refuse[1] && hon.refuse[2] === null && hon.snapLo.where === "snap" && hon.snapHi.where === "snap"
      && !/interp\./.test(hon.snapLo.q) && !/interp\./.test(hon.snapHi.q) && /not interpolated/.test(hon.snapLo.c + hon.snapHi.c)
      && hon.okNear.where === "interp" && /interp\./.test(hon.okNear.q), JSON.stringify([hon.refuse, hon.snapLo.where, hon.snapLo.q, hon.okNear.where]));
    check("p2d rule 1: frames whose label sets differ are not blended; the labels are the nearest computed frame's (switch at the midpoint)",
      /label/.test(hon.labelRefusal || "") && hon.labSnap.where === "snap" && !/synthetic label/.test(hon.labSnap.txt) && /synthetic label/.test(hon.labSnap2) && /nearest computed frame/.test(hon.labSnap2),
      JSON.stringify([hon.labelRefusal, hon.labSnap.where]));
    check("p2d rule 3/4: T_c is never blended: across a T_c change only the two computed values appear, and the crossing reads 'between frames 15 and 16'",
      hon.tcSeen.length === 2 && /between frames 15 and 16/.test(hon.tcCross), JSON.stringify([hon.tcSeen, hon.tcCross]));
    check("p2d rule 5: a runaway frame ends the displayable range: no blend into it, the readouts read 'runaway' at and after it, the cells past it are marked, and playback stops on it",
      hon.runIdx === 30 && hon.runBefore.where === "snap" && /^runaway$/.test(hon.runAt.q) && /displayable range ends/.test(hon.runAt.c) && hon.runAt.v === "n/a"
      && /^runaway$/.test(hon.runAfter.q) && hon.beyond === 17 && hon.runPlay.pos <= 30 + 1e-9 && hon.runPlay.playing === "false" && hon.restored,
      JSON.stringify([hon.runIdx, hon.runBefore, hon.runAt.q, hon.beyond, hon.runPlay]));
    const dens = await page.evaluate(() => ({ txt: document.querySelector(".anim-density")?.textContent || "", hidden: document.querySelector(".anim-density")?.hidden, btn: !!document.querySelector(".anim-refine"), leg: document.querySelector(".tl-legend").textContent }));
    check("p2d rule 6: the 48-frame T_hs run (7 K apart) near T_c is flagged too coarse with a one-click refine around T_c; the legend names the snap",
      dens.hidden === false && /within ±20 K of T_c/.test(dens.txt) && dens.btn && /snapped/.test(dens.leg), dens.txt.slice(0, 120));
    await page.screenshot({ path: path.join(SHOTS_P2, "designer-animate-staged.png") });
    const sliceT = await page.evaluate(() => !!document.querySelector(".anim-slice"));
    await page.selectOption(".anim-param", "dot.gamma_scale");
    await page.waitForFunction(() => document.querySelector(".anim-sheet")?.dataset.param === "dot.gamma_scale");
    const sliceA = await page.evaluate(() => ({ chip: document.querySelector(".anim-slice")?.textContent || "", tag: document.querySelector(".anim-field .tag")?.dataset.tag }));
    check("p2d rule 7: an animated [A] input is labelled 'one range at a time, not a prediction' on the sheet; the [E] T_hs sheet carries no such chip",
      !sliceT && sliceA.tag === "A" && sliceA.chip === "one range at a time, not a prediction", JSON.stringify([sliceT, sliceA]));
    await page.click(".anim-sheet .sheet-head .icon-btn");

    // ---------------------------------------------------------------- p2a: response animation, edge headline (n = 12)
    await page.goto(`${BASE}/#/design/edge-inp-gainp-design`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("edge-inp-gainp"), null, { timeout: 30000 });
    await waitRunDone(page, 120000);
    await page.click('.segmented [data-value="headline"]');
    await page.click(".hs-anim");
    await page.waitForSelector(".anim-sheet");
    {
      // review finding 1: From/To are constrained to the META validity band, and an out-of-band run is refused before any request
      const band = await page.evaluate(() => ({ min: document.querySelector(".anim-lo").min, max: document.querySelector(".anim-hi").max }));
      await page.fill(".anim-lo", "-50");
      await page.fill(".anim-hi", "1000");
      await page.click(".anim-go");
      const msg = await page.$eval(".anim-status", (e) => e.textContent);
      const phase = await page.$eval(".anim-sheet", (e) => e.dataset.phase);
      check("review 1: the Designer From/To inputs carry the band (min 4, max 350) and an out-of-band run is refused in place (validity band named, nothing computed)",
        band.min === "4" && band.max === "350" && /validity band/.test(msg) && phase === "setup", `${band.min}..${band.max} | ${msg.slice(0, 70)} | ${phase}`);
    }
    await page.fill(".anim-lo", "200");
    await page.fill(".anim-hi", "350");
    await page.fill(".anim-n", "12");
    const modeChip = await page.$eval(".anim-mode", (e) => e.textContent);
    await page.click(".anim-go");
    await page.waitForFunction(() => document.querySelector(".anim-sheet")?.dataset.phase === "ready", null, { timeout: 180000 });
    const ea = await animAudit(page, 12);
    const edgeLab = await page.evaluate(() => { const A = document.querySelector(".anim-sheet").__anim; A.seek(5); return { lab: document.querySelector(".anim-labels").textContent, tc: document.querySelector('[data-readout="T_c"] .ro-caption').textContent, q: document.querySelector('[data-readout="g2_op"] .ro-qual').textContent }; });
    check("p2a edge animate: headline model, 12 computed frames on the scrubber", modeChip === "headline model" && ea.cells.length === 12 && ea.cells.every((c) => c === "computed"), `${modeChip}; ${ea.cells.join(",")}`);
    check("p2a edge animate: readouts change with the playhead; 'headline · interp.' between frames, 'headline' on a frame",
      ea.v0 !== ea.v1 && /^headline · interp\.$/.test(ea.qInterp) && ea.qFrame === "headline", `${ea.v0} -> ${ea.v1}; ${ea.qFrame} | ${ea.qInterp}`);
    check("p2a edge animate: frame labels carry the headline-model label; T_c says it needs the full grid (never interpolated)",
      /headline model/.test(edgeLab.lab) && /full T_hs grid/.test(edgeLab.tc), `${edgeLab.lab.slice(0, 80)} | ${edgeLab.tc}`);
    check("p2a edge animate: labels switch only at computed frames", ea.badSwitch === 0 && ea.badHeld === 0, `${ea.badSwitch}/${ea.badHeld}`);
    check("p2a edge animate: JS per playback frame < 2 ms (Performance API: mean and p95)", ea.perf.n > 30 && ea.perf.mean < 2 && ea.perf.p95 < 2, JSON.stringify(ea.perf));
    await page.waitForFunction(() => document.querySelector(".v3-root")?.dataset.ready === "1", null, { timeout: 30000 }).catch(() => {});
    await settle(page, 2500); // let the rest-time scene rebuild from the audit's pause finish before counting contexts
    const gl = await page.evaluate(async () => {
      const A = document.querySelector(".anim-sheet").__anim;
      const root = document.querySelector(".v3-root");
      const api = root?.__viz3d;
      let ctxCalls = 0;
      const orig = HTMLCanvasElement.prototype.getContext;
      HTMLCanvasElement.prototype.getContext = function (t, ...r) { if (/webgl/.test(String(t))) ctxCalls++; return orig.call(this, t, ...r); };
      const nCanvas = document.querySelectorAll("canvas").length;
      const g = [];
      A.seek(0); g.push(api?._glow?.());
      A.seek(11); g.push(api?._glow?.());
      A.seek(1);
      A.play();
      await new Promise((r) => setTimeout(r, 1500));
      const during = { ctxCalls, nCanvas: document.querySelectorAll("canvas").length, playing: A.phase };
      A.pause();
      HTMLCanvasElement.prototype.getContext = orig;
      return { hasApi: !!api?.setPlayhead, mode: api?.mode, g, during };
    });
    check("p2d 3D hook: the scene follows the playhead through setPlayhead (a glow value is set, no rebuild) and no WebGL context or canvas is created during playback",
      gl.hasApi && (gl.mode !== "webgl" || gl.g.every((x) => typeof x === "number")) && gl.during.ctxCalls === 0, JSON.stringify(gl));
    await page.focus(".tl-track");
    await page.keyboard.press("ArrowRight");
    await page.waitForFunction(() => /3D at computed frame/.test(document.querySelector("#vm-note")?.textContent || ""), null, { timeout: 10000 }).catch(() => {});
    const vm = await page.$eval("#vm-note", (e) => e.textContent).catch(() => "");
    check("p2a edge animate: the 3D scene follows the playhead when it rests (render on demand, at a computed frame)", /3D at computed frame \d+: T_hs/.test(vm), vm);
    await page.evaluate(() => { document.querySelector(".anim-sheet").__anim.seek(4.5); document.querySelector(".anim-player").scrollIntoView({ block: "end" }); });
    await settle(page, 900);
    await page.screenshot({ path: path.join(SHOTS_P2, "designer-animate-edge.png") });
    await page.click(".anim-sheet .sheet-head .icon-btn");

    // reduced motion: playback off, scrubbing still works
    const rctx = await browser.newContext({ viewport: VIEW, reducedMotion: "reduce" });
    const rp = await rctx.newPage();
    rp.on("pageerror", (e) => consoleErrors.push(`reduced pageerror: ${e.message}`));
    await rp.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await rp.waitForFunction(() => document.querySelector(".board-title")?.textContent.startsWith("staged-device"), null, { timeout: 30000 });
    await waitRunDone(rp);
    await rp.click(".hs-anim");
    await rp.click(".anim-go");
    await rp.waitForFunction(() => document.querySelector(".anim-sheet")?.dataset.phase === "ready", null, { timeout: 60000 });
    const rm = await rp.evaluate(() => {
      const A = document.querySelector(".anim-sheet").__anim;
      const v = () => document.querySelector('[data-readout="g2_op"] .ro-value').textContent;
      A.seek(0); const a = v(); A.seek(40); const b = v();
      A.play();
      return { disabled: document.querySelector(".tl-play").disabled, playing: document.querySelector(".timeline").dataset.playing, a, b };
    });
    check("p2a prefers-reduced-motion: playback disabled, the scrubber still moves the readouts", rm.disabled && rm.playing === "0" && rm.a !== rm.b, JSON.stringify(rm));
    await rctx.close();

    // ---------------------------------------------------------------- mobile 390px
    const mctx = await browser.newContext({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 2 });
    const m = await mctx.newPage();
    m.on("console", (msg) => { if (msg.type() === "error") consoleErrors.push(`mobile: ${msg.text()}`); });
    m.on("pageerror", (e) => consoleErrors.push(`mobile pageerror: ${e.message}`));
    await m.goto(`${BASE}/#/design/staged-device-design`, { waitUntil: "domcontentloaded" });
    await m.waitForSelector(".designer:not(.is-loading)");
    await waitRunDone(m);
    await m.waitForTimeout(900);
    const hscroll = await m.evaluate(() => document.documentElement.scrollWidth > window.innerWidth + 1);
    check("mobile 390px: no horizontal page scroll", !hscroll);
    const chainA = await m.evaluate(() => {
      const w = document.querySelector(".chain-wrap"), sc = document.querySelector(".chain-scroll"), tr = document.querySelector(".chain-track");
      if (!w || !sc || !tr) return { missing: true };
      const th = tr.firstElementChild.getBoundingClientRect(), trr = tr.getBoundingClientRect();
      const mask = getComputedStyle(sc).maskImage || getComputedStyle(sc).webkitMaskImage || "none";
      return { overflow: w.dataset.overflow, edge: w.dataset.edge, mask, trackW: trr.width, trackH: trr.height, thumbW: th.width,
        sw: sc.scrollWidth, cw: sc.clientWidth, thin: getComputedStyle(sc).scrollbarWidth };
    });
    check("R1: mobile chain row overflows with a visible scroll affordance (edge fade mask + thin track, thumb < track)",
      !chainA.missing && chainA.overflow === "1" && chainA.edge === "start" && /gradient/.test(chainA.mask) && chainA.trackW > 0 && chainA.trackH >= 2
      && chainA.thumbW > 0 && chainA.thumbW < chainA.trackW - 1 && chainA.thin === "thin", JSON.stringify(chainA).slice(0, 220));
    await m.screenshot({ path: path.join(SHOTS, "designer-mobile.png"), fullPage: true });
    const chainEnd = await m.evaluate(async () => {
      const sc = document.querySelector(".chain-scroll");
      sc.scrollLeft = sc.scrollWidth;
      await new Promise((r) => setTimeout(r, 150));
      return document.querySelector(".chain-wrap").dataset.edge;
    });
    check("R1: scrolling the chain to its end flips the fade to the leading edge", chainEnd === "end", String(chainEnd));
    await mctx.close();
  } catch (e) {
    check("script ran to completion", false, e.message.split("\n")[0]);
    try { await page.screenshot({ path: path.join(SHOTS, "failure.png") }); } catch { /* ignore */ }
  }

  const relevant = consoleErrors.filter((t) => !/favicon/i.test(t));
  check("no console errors", relevant.length === 0, relevant.slice(0, 5).join(" || "));
  await browser.close();
  if (OWN) OWN.stop();

  const passed = results.filter((r) => r.ok).length;
  console.log(`${passed}/${results.length} studio web checks passed`);
  console.log(`screenshots: ${SHOTS}`);
  process.exit(passed === results.length ? 0 : 1);
})();
