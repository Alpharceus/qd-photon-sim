// FSIM Studio Composer workspace: Playwright smoke + screenshots (devlib-p11).
//
// Usage: node verify/verify_studio_composer.cjs http://127.0.0.1:8772     (an already-running server)
//        node verify/verify_studio_composer.cjs                           (starts its own on 8793)
// Saving a composed card is exercised against a MOCKED PUT /api/cards/<slug> (the request body is
// asserted), so a server without a temp card folder never gets written to. No physics number is
// pinned: the checks look at structure, provenance chips, network calls and chart traces.
// Playwright resolves from .workers/studio/node/node_modules (NODE_PATH also honoured).
// Prints "N/N composer checks passed"; exits 0 iff all pass.

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const { startStudioServer } = require("./studio_server.cjs");
let BASE = (process.argv[2] || "").replace(/\/$/, "");
let OWN = null;
const SHOTS = path.join(ROOT, ".workers", "devlib", "shots");

function loadPlaywright() {
  const tries = [
    path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright-core"),
    path.join(ROOT, ".workers", "studio", "node", "node_modules", "playwright"),
    "playwright-core", "playwright",
  ];
  for (const t of tries) { try { return require(t); } catch { /* next */ } }
  console.error("playwright-core not found: npm i playwright-core under .workers/studio/node");
  process.exit(2);
}

const results = [];
function check(name, ok, detail = "") {
  results.push({ name, ok: !!ok, detail });
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}${detail ? `  (${detail})` : ""}`);
}

async function settle(page, ms = 600) { await page.waitForTimeout(ms); await page.evaluate(() => document.fonts.ready); }

async function b612Audit(page) {
  return page.evaluate(() => {
    const bad = [];
    const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
    let n;
    while ((n = walker.nextNode())) {
      if (!/\d/.test(n.nodeValue)) continue;
      const el = n.parentElement;
      if (!el || !el.getClientRects().length) continue;
      if (/B612/i.test((getComputedStyle(el).fontFamily || "").split(",")[0])) bad.push(n.nodeValue.trim().slice(0, 20));
    }
    return bad;
  });
}

// Examples and Tiers are two menus on the "Start from:" row (popover lists), not chip rows.
async function openStart(page, kind) {
  await page.click(`.ws-compose .start-btn[data-menu="${kind}"]`);
  await page.waitForSelector(".start-pop .start-item", { timeout: 5000 });
  // popover() attaches its outside-click listener in a setTimeout(0): let it land before the menu is closed
  // (a menu closed within that tick leaves a stale listener that eats the next popover's first click).
  await page.waitForTimeout(300);
}
async function closeStart(page) {
  await page.keyboard.press("Escape");
  await page.waitForSelector(".start-pop", { state: "detached", timeout: 5000 }).catch(() => {});
}
async function pickExample(page, name) {
  await openStart(page, "examples");
  await page.click(`.start-pop .example-chip[data-example="${name}"]`);
  await page.waitForFunction((n) => document.querySelector(".ws-compose .start-examples")?.dataset.example === n, name, { timeout: 10000 });
}

async function waitRun(page, timeout = 120000) {
  await page.waitForFunction(() => {
    const s = document.querySelector(".ws-compose .status-word");
    return s && ["done", "error", "cancelled"].includes(s.dataset.state);
  }, null, { timeout });
  return page.$eval(".ws-compose .status-word", (el) => el.dataset.state);
}

(async () => {
  if (!BASE) { OWN = await startStudioServer(ROOT, 8793); BASE = OWN.base; }
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=d3d11", "--enable-gpu", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const errors = [];
  let mockedFailure = false;   // the failed-run capture answers 500 on purpose; its console line is not an app error
  page.on("console", (m) => { if (m.type() === "error" && !(mockedFailure && /status of 500/.test(m.text()))) errors.push(m.text()); });
  page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
  page.on("requestfailed", (r) => { if (!r.url().includes("/events")) errors.push(`requestfailed: ${r.url()} ${r.failure()?.errorText}`); });
  const validates = [];
  page.on("request", (r) => { if (r.url().endsWith("/api/compose/validate")) validates.push(r.postData()); });
  let putBody = null;
  await page.route("**/api/cards/*", async (route) => {
    const req = route.request();
    if (req.method() === "PUT") {
      putBody = JSON.parse(req.postData() || "{}");
      await route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ path: "cards/studio/mock.yaml", name: "mock" }) });
    } else if (req.method() === "GET" && putBody && req.url().endsWith("/composer-verify-card")) {
      await route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ name: "composer-verify-card", kind: "design", design: putBody.design, meta: {} }) });
    } else await route.continue();
  });

  try {
    const reg = await (await fetch(`${BASE}/api/compose/registry`)).json();
    const AX = [["drive", "drives"], ["material", "materials"], ["structure", "structures"], ["setup", "setups"]];

    // ---------------------------------------------------------------- load
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque"), null, { timeout: 60000 });
    await settle(page);
    const rail = await page.$$eval(".rail-item", (els) => els.map((e) => e.dataset.ws));
    check("rail: Composer follows Designer, Papers present", rail.indexOf("design") + 1 === rail.indexOf("compose") && rail.includes("papers"), rail.join(","));
    const chips = await page.$$eval(".ws-compose .kind-chip", (els) => els.map((e) => e.dataset.axis));
    check("a kind chip on every axis of the chain", ["drive", "material", "structure", "setup"].every((a) => chips.includes(a)), chips.join(","));
    check("registry has entries on every axis", AX.every(([, k]) => (reg[k] || []).length > 0), AX.map(([, k]) => `${k}:${(reg[k] || []).length}`).join(" "));
    const plaque0 = await page.$eval(".ws-compose .plaque", (e) => e.dataset.status);
    check("default spec validates (plaque VALID)", plaque0 === "pass", plaque0);
    check("one 'Start from:' row holds the Examples and Tiers menus (no chip rows)", (await page.$$(".ws-compose .start-row .start-btn")).length === 2 && (await page.$$(".ws-compose .tier-strip")).length === 0 && /^Start from:/.test(await page.$eval(".ws-compose .start-row", (e) => e.textContent.trim())));
    await openStart(page, "tiers");
    const tierItems = await page.$$eval(".start-pop .tier-chip", (els) => els.map((e) => e.textContent));
    check("the Tiers menu lists the registry tiers", tierItems.length === (reg.tiers || []).length && (reg.tiers || []).length > 0, `${tierItems.length}`);
    await closeStart(page);
    const startCaps = await page.$$eval(".ws-compose .start-btn", (els) => els.map((e) => getComputedStyle(e).textTransform));
    check("Start-from menu names are sentence case, not engraved caps", startCaps.every((t) => t === "none"), startCaps.join(","));
    const tp0 = await page.evaluate(() => { const t = document.querySelector(".ws-compose .temp-plate").getBoundingClientRect(); const rb = document.querySelector(".ws-compose .runbar").getBoundingClientRect(); return { top: t.top, bottom: t.bottom, vh: window.innerHeight, rb: rb.top }; });
    check("the temperature plate is fully visible in the 1440x900 first viewport", tp0.top >= 48 && tp0.bottom <= Math.min(tp0.vh, tp0.rb + 1), JSON.stringify(tp0));
    const pq0 = await page.$eval(".ws-compose .plaque", (e) => ({ line: e.classList.contains("plaque-line"), bg: getComputedStyle(e).backgroundColor, word: e.querySelector(".plaque-word").textContent }));
    check("before a run the validity is an outline status line, not a filled plaque", pq0.line && /VALID/.test(pq0.word) && (pq0.bg === "rgba(0, 0, 0, 0)" || pq0.bg === "transparent"), JSON.stringify(pq0));
    const hint0 = await page.$eval(".ws-compose .shelf-hint", (e) => ({ hidden: e.hidden, t: e.textContent }));
    check("the empty shelf says 'Press RUN to fill these (about ...)'", !hint0.hidden && /Press RUN to fill these/.test(hint0.t), hint0.t);
    const ids0 = await page.$$eval(".ws-compose .id-k", (els) => els.map((e) => e.textContent));
    check("the device ID is labelled DEVICE", ids0.includes("DEVICE"), ids0.join(","));
    check("the Composer subtitle is three steps", (await page.$$(".ws-compose .board-steps li")).length === 3);
    const ctxCard = await page.$eval(".ctx-card", (e) => ({ hidden: e.hidden, disp: getComputedStyle(e).display }));
    check("the Designer's card selector is hidden on #/compose", ctxCard.hidden && ctxCard.disp === "none", JSON.stringify(ctxCard));
    const kcH = await page.$$eval(".ws-compose .kind-chip", (els) => els.map((e) => Math.round(e.getBoundingClientRect().height)));
    check("kind chips are one line (short label, full label in the picker)", kcH.every((x) => x <= 34), kcH.join(","));
    const mvk = await page.$$eval(".ws-compose .mv-k", (els) => els.map((e) => e.textContent));
    check("mount summaries show no raw snake_case names", mvk.every((t) => !/_/.test(t)), mvk.join(","));

    // ---------------------------------------------------------------- pickers on every axis
    for (const [axis, regKey] of AX) {
      await page.click(`.ws-compose .kind-chip[data-axis="${axis}"]`);
      await page.waitForSelector(".cp-entry", { timeout: 10000 });
      await page.waitForFunction(() => !/Checking/.test(document.querySelector(".cp-note")?.textContent || ""), null, { timeout: 60000 });
      const info = await page.$$eval(".cp-entry", (els) => els.map((e) => ({ kind: e.dataset.kind, status: e.dataset.status, dis: e.getAttribute("aria-disabled") === "true", txt: e.textContent })));
      const ready = (reg[regKey] || []).filter((e) => e.status === "ready").length;
      check(`picker ${axis}: one card per registry entry, with summary text`, info.length === reg[regKey].length && info.every((i) => i.txt.length > 20), `${info.length}/${reg[regKey].length}`);
      const np = info.filter((i) => i.status === "needs_physics");
      check(`picker ${axis}: needs_physics entries are greyed, say "awaiting physics" and are not selectable`, np.every((i) => i.dis && /awaiting physics/.test(i.txt)), `${np.length} awaiting; ${ready} ready`);
      if (axis === "drive" && np.length) {
        const before = await page.$eval(`.ws-compose .kind-chip[data-axis="drive"]`, (e) => e.dataset.kind);
        await page.click(`.cp-entry[data-kind="${np[0].kind}"]`, { force: true });
        await page.waitForTimeout(200);
        const after = await page.$eval(`.ws-compose .kind-chip[data-axis="drive"]`, (e) => e.dataset.kind);
        check("clicking an awaiting-physics entry changes nothing", before === after && (await page.$(".cp-entry")) !== null);
      }
      const dashed = await page.$$eval(".cp-entry[aria-disabled='true']", (els) => els.every((e) => getComputedStyle(e).borderTopStyle === "dashed"));
      check(`picker ${axis}: greyed entries carry a dashed rim`, dashed);
      await page.keyboard.press("Escape");
      await page.waitForSelector(".cp-entry", { state: "detached", timeout: 5000 });
    }

    // ---------------------------------------------------------------- swap a structure
    await page.click(`.ws-compose .kind-chip[data-axis="structure"]`);
    await page.waitForSelector(".cp-entry");
    await page.waitForFunction(() => !/Checking/.test(document.querySelector(".cp-note")?.textContent || ""), null, { timeout: 60000 });
    const cur = await page.$eval(`.ws-compose .kind-chip[data-axis="structure"]`, (e) => e.dataset.kind);
    const alt = await page.$$eval(".cp-entry", (els, c) => els.filter((e) => e.dataset.kind !== c && e.getAttribute("aria-disabled") !== "true").map((e) => e.dataset.kind), cur);
    check("at least one other structure is selectable from the default chain", alt.length > 0, alt.slice(0, 4).join(","));
    const sceneBefore = await page.evaluate(() => document.querySelector(".ws-compose .vm-note")?.textContent || "");
    const nV0 = validates.length;
    let picked = null;
    for (const k of alt) {
      await page.click(`.cp-entry[data-kind="${k}"]`);
      await page.waitForFunction((kk) => document.querySelector(`.ws-compose .kind-chip[data-axis="structure"]`)?.dataset.kind === kk, k, { timeout: 5000 });
      await page.waitForTimeout(900);
      picked = k;
      break;
    }
    const swapped = await page.$eval(`.ws-compose .kind-chip[data-axis="structure"]`, (e) => e.dataset.kind);
    check("picking a structure swaps the chip and revalidates", swapped === picked && picked !== cur && validates.length > nV0, `${cur} -> ${swapped}`);
    await page.waitForFunction(() => { const v = document.querySelector(".ws-compose .viewport-slot"); return v && (v.dataset.ready === "1" || v.querySelector("[data-mode='fallback2d']") || v.querySelector(".v3-message")); }, null, { timeout: 30000 }).catch(() => {});
    const sceneAfter = await page.evaluate(() => document.querySelector(".ws-compose .vm-note")?.textContent || "");
    const vpState = await page.evaluate(() => { const v = document.querySelector(".ws-compose .viewport-slot"); return { mode: v?.dataset.mode || v?.querySelector("[data-mode]")?.dataset.mode || "", kids: v?.children.length || 0 }; });
    check("the viewport shows the chosen structure (3D builder or the flat fallback plate) with a catalog ID", sceneAfter.startsWith("FS-") && sceneAfter !== sceneBefore && vpState.kids > 0, `${sceneBefore} -> ${sceneAfter}; mode=${vpState.mode}`);
    // restore a runnable chain if the swap broke it
    let ok = await page.$eval(".ws-compose .plaque", (e) => e.dataset.status);
    if (ok !== "pass") {
      await page.click(`.ws-compose .kind-chip[data-axis="structure"]`);
      await page.waitForSelector(".cp-entry");
      await page.waitForFunction(() => !/Checking/.test(document.querySelector(".cp-note")?.textContent || ""), null, { timeout: 60000 });
      await page.click(`.cp-entry[data-kind="${cur}"]`);
      await page.waitForTimeout(900);
    }
    ok = await page.$eval(".ws-compose .plaque", (e) => e.dataset.status);
    check("chain valid again before the run", ok === "pass", ok);

    // ---------------------------------------------------------------- parameter sheet
    await page.click(`.ws-compose .mount-body[data-axis="structure"]`);
    await page.waitForSelector(".ws-compose .sheet .field-row");
    await page.waitForTimeout(500);
    const sheetBox = await page.evaluate(() => { const r = document.querySelector(".ws-compose .sheet").getBoundingClientRect(); return { top: r.top, bottom: r.bottom, vh: window.innerHeight }; });
    check("opening a mount's sheet scrolls it into view (its header is on screen)", sheetBox.top >= 40 && sheetBox.top < sheetBox.vh - 120, JSON.stringify(sheetBox));
    const labs = await page.$$eval(".ws-compose .sheet .fr-label", (els) => els.map((e) => e.textContent));
    check("sheet labels are symbols or plain words, not raw snake_case", labs.length > 0 && labs.every((t) => !/_/.test(t)), labs.slice(0, 6).join(" | "));
    const nField = await page.$$eval(".ws-compose .sheet .pf-item", (e) => e.length);
    if (nField > 12) {
      check("a sheet of more than 12 fields has a filter box", (await page.$(".ws-compose .sheet .pf-filter")) !== null, `${nField} fields`);
      await page.fill(".ws-compose .sheet .pf-filter", "zzzz-no-such-field");
      await page.waitForTimeout(150);
      const shown = await page.$$eval(".ws-compose .sheet .pf-item:not([hidden])", (e) => e.length);
      await page.fill(".ws-compose .sheet .pf-filter", "");
      await page.waitForTimeout(150);
      check("the filter hides non-matching fields and restores them", shown === 0 && (await page.$$eval(".ws-compose .sheet .pf-item:not([hidden])", (e) => e.length)) > 0);
    }
    await page.screenshot({ path: path.join(SHOTS, "finish-composer-sheet-open.png") });
    const nRows = await page.$$eval(".ws-compose .sheet .field-row", (e) => e.length);
    const nTagsel = await page.$$eval(".ws-compose .sheet .pf-tagsel", (e) => e.length);
    check("parameter sheet is generated from the registry (rows with tag selectors)", nRows > 0 && nTagsel === nRows, `${nRows} rows`);
    const numRow = await page.$(".ws-compose .sheet .pf-item .fr-point input.fr-num");
    const nV1 = validates.length;
    if (numRow) {
      const before = await numRow.inputValue();
      const bump = String(Number(before || 1) * 1.01 + 0.01);
      await numRow.fill(bump);
      await numRow.press("Tab");
      await page.waitForTimeout(700);
      check("editing a parameter triggers validate", validates.length > nV1);
      const edited = await page.$$eval(".ws-compose .sheet .pf-item.is-edited", (els) => els.map((e) => ({ style: getComputedStyle(e).borderTopStyle, tag: e.dataset.tag })));
      check("an edited leaf gets a dashed rim and the assumption tag [A]", edited.length >= 1 && edited[0].style === "dashed" && edited[0].tag === "A", JSON.stringify(edited[0] || {}));
      const mountVal = await page.$$eval(`.ws-compose .mount[data-axis="structure"] .mount-foot`, (e) => e[0]?.textContent || "");
      check("the mount foot counts the edit", /1 edited|edited/.test(mountVal), mountVal);
      await page.click(".ws-compose .sheet .pf-item.is-edited .pf-reset");
      await page.waitForTimeout(500);
    } else check("editing a parameter triggers validate", false, "no numeric row on this structure");
    // advanced group folded
    const adv = await page.$(".ws-compose .sheet .pf-adv");
    if (adv) check("advanced parameters are folded by default", !(await adv.evaluate((e) => e.open)));
    await page.click(`.ws-compose .mount-body[data-axis="structure"]`);

    // ---------------------------------------------------------------- run: point
    await page.waitForFunction(() => !document.querySelector(".ws-compose .run-key")?.disabled, null, { timeout: 30000 });
    await page.click(".ws-compose .run-key");
    const st = await waitRun(page);
    check("RUN (point) completes", st === "done", st);
    await settle(page, 500);
    const ro = await page.$$eval(".ws-compose .readouts .readout", (els) => els.map((e) => ({ k: e.dataset.readout, state: e.dataset.state, tag: !!e.querySelector(".tag"), cap: e.querySelector(".ro-caption")?.textContent || "" })));
    check("six readouts, each with a value", ro.length === 6 && ro.every((r) => r.state === "value"), JSON.stringify(ro.map((r) => `${r.k}:${r.state}`)));
    check("every readout number carries a provenance tag chip", ro.every((r) => r.tag));
    check("every readout names its convention", ro.every((r) => r.cap.length > 3));
    const roSym = await page.$$eval(".ws-compose .readouts .ro-label", (els) => els.map((e) => ({ t: e.textContent, k: [...e.querySelectorAll(".katex")].map((x) => getComputedStyle(x.closest(".math")).textTransform), sup: !!e.querySelector(".katex .msupsub"), n: e.querySelectorAll(".katex").length })));
    check("KPI symbols are KaTeX notation and not upper-cased (g(2)(0) with a superscript, F_P, I_HOM, B_lens, tau)", roSym.filter((r) => r.n).length >= 4 && roSym.every((r) => r.k.every((x) => x === "none")) && roSym.some((r) => r.sup) && !roSym.some((r) => /^T$/.test(r.t)), JSON.stringify(roSym.map((r) => r.t)));
    const caps = await page.$$eval(".ws-compose .readouts .ro-caption", (els) => els.map((e) => e.textContent.length));
    check("KPI captions are one short line; the formula sits behind 'why'", caps.every((n) => n <= 50) && (await page.$$(".ws-compose .readouts .ro-why:not([hidden])")).length >= 4, caps.join(","));
    await page.click('.ws-compose .readout[data-readout="g2_0"] .ro-why');
    await page.waitForSelector(".popover .why-lede", { timeout: 3000 });
    check("'why' opens the formula", (await page.$eval(".popover .why-lede", (e) => e.textContent.length)) > 20);
    await page.waitForTimeout(300);   // popover()'s outside-click listener lands in a setTimeout(0)
    await page.keyboard.press("Escape");
    const plaque1 = await page.$eval(".ws-compose .plaque", (e) => ({ s: e.dataset.status, t: e.textContent, line: e.classList.contains("plaque-line") }));
    check("after a run the plaque is the filled result plaque", !plaque1.line && /RESULT/.test(plaque1.t), plaque1.t.slice(0, 30));
    const hint1 = await page.$eval(".ws-compose .shelf-hint", (e) => e.hidden);
    check("the 'Press RUN' hint is gone once the shelf is filled", hint1 === true);
    const ctxBox = await page.evaluate(() => {
      const els = [...document.querySelectorAll(".ctxbar > *")].filter((e) => e.offsetParent !== null && e.getBoundingClientRect().width > 0).map((e) => ({ c: e.className || e.id, l: e.getBoundingClientRect().left, r: e.getBoundingClientRect().right })).sort((a, b) => a.l - b.l);
      const bad = []; for (let i = 1; i < els.length; i++) if (els[i].l < els[i - 1].r - 1) bad.push(`${els[i - 1].c}/${els[i].c}`);
      const chip = document.querySelector("#run-chip").getBoundingClientRect();
      return { bad, chipRight: chip.right, vw: window.innerWidth };
    });
    check("top bar: no two items overlap at 1440 after a run, the run chip stays inside", ctxBox.bad.length === 0 && ctxBox.chipRight <= ctxBox.vw, JSON.stringify(ctxBox));
    check("plaque shows solver, regime and the phonon note after a run", /solver/.test(plaque1.t) && /regime/.test(plaque1.t) && /phonon/i.test(plaque1.t), plaque1.t.slice(0, 120));
    const beamA = await page.$eval(".ws-compose .beam", (e) => Number(e.dataset.alpha));
    check("the beam follows the run (alpha above rest level)", beamA > 0.45, String(beamA));

    // ---------------------------------------------------------------- run: T sweep
    await page.click('.ws-compose .seg[data-value="tsweep"]');
    await page.waitForTimeout(200);
    await page.click(".ws-compose .run-key");
    const st2 = await waitRun(page);
    check("RUN (T sweep) completes", st2 === "done", st2);
    await page.waitForFunction(() => document.querySelector('.ws-compose [data-chart="tcurve"] .chart-plot .main-svg'), null, { timeout: 30000 }).catch(() => {});
    const tc = await page.$eval('.ws-compose [data-chart="tcurve"] .chart-plot', (gd) => ({ n: (gd.data || []).length, pts: (gd.data?.[0]?.x || []).length })).catch(() => ({ n: 0, pts: 0 }));
    check("T sweep draws curves against temperature", tc.n >= 1 && tc.pts > 1, JSON.stringify(tc));

    // ---------------------------------------------------------------- run: series
    await page.click('.ws-compose .seg[data-value="series"]');
    await page.waitForSelector(".ws-compose .rb-series .rb-path");
    const nPaths = await page.$$eval(".ws-compose .rb-path option", (e) => e.length);
    check("series x offers T_K and the sweepable parameters of the chain", nPaths >= 2, `${nPaths} paths`);
    await page.click(".ws-compose .run-key");
    const st3 = await waitRun(page);
    check("series run completes", st3 === "done", st3);
    await page.waitForFunction(() => document.querySelector('.ws-compose [data-chart="series"] .chart-plot .main-svg'), null, { timeout: 30000 }).catch(() => {});
    const sr = await page.$eval('.ws-compose [data-chart="series"] .chart-plot', (gd) => (gd.data || []).map((t) => ({ dash: t.line?.dash, w: t.line?.width, n: (t.x || []).length, name: t.name }))).catch(() => []);
    check("series chart is drawn with one trace per observable", sr.length >= 1 && sr.every((t) => t.n >= 2), JSON.stringify(sr));
    check("series traces follow the tag dash grammar (no dotted; dashed >= 2.5 px)", sr.every((t) => t.dash !== "dot" && (t.dash !== "dash" && t.dash !== "longdash" || t.w >= 2.5)), JSON.stringify(sr.map((t) => `${t.dash}:${t.w}`)));
    check("an [A] series trace is dashed at >= 2.5 px", sr.some((t) => t.dash === "dash" && t.w >= 2.5) || sr.every((t) => t.dash === "solid"), JSON.stringify(sr.map((t) => `${t.dash}:${t.w}`)));
    const aChip = await page.evaluate(() => { const t = document.querySelector(".ws-compose .tag-a"); if (!t) return null; const s = getComputedStyle(t); return { style: s.borderTopStyle, w: parseFloat(s.borderTopWidth) }; });
    check("an [A] tag chip is a dashed outline", aChip && aChip.style === "dashed" && aChip.w >= 1, JSON.stringify(aChip));

    // ---------------------------------------------------------------- visual contract
    check("no B612 numerals", (await b612Audit(page)).length === 0);
    const beamFill = await page.evaluate(() => {
      const probe = document.createElement("i"); probe.style.color = "var(--beam)"; document.body.append(probe);
      const beam = getComputedStyle(probe).color; probe.remove();
      return [...document.querySelectorAll(".ws-compose button, .ws-compose .btn, .ws-compose .seg, .ws-compose .chip, .ws-compose .kind-chip")]
        .filter((e) => e.offsetParent !== null && getComputedStyle(e).backgroundColor === beam).map((e) => e.className);
    });
    check("only RUN carries the beam fill", beamFill.length === 1 && /run-key/.test(beamFill[0]), JSON.stringify(beamFill));
    const radii = await page.evaluate(() => [...document.querySelectorAll(".ws-compose .mount, .ws-compose .kind-chip, .ws-compose .btn, .ws-compose .fr-num, .ws-compose .plaque, .ws-compose .chart-card")].filter((e) => e.offsetParent !== null).map((e) => parseFloat(getComputedStyle(e).borderTopLeftRadius)).filter((r) => r > 2));
    check("plates, chips, buttons and inputs keep the 2 px corner", radii.length === 0, JSON.stringify(radii.slice(0, 5)));
    const src = fs.readFileSync(path.join(ROOT, "fsim_studio", "web", "js", "workspaces", "composer.js"), "utf8");
    const kinds = AX.flatMap(([, k]) => reg[k].map((e) => e.kind));
    const hits = kinds.filter((k) => new RegExp(`(^|[^A-Za-z0-9_])${k}([^A-Za-z0-9_]|$)`).test(src));
    check("composer.js names no registry kind (no per-kind branching)", hits.length === 0, hits.join(","));
    const forms = ["paramForm.js", "chipPicker.js"].flatMap((f) => { const t = fs.readFileSync(path.join(ROOT, "fsim_studio", "web", "js", "ui", f), "utf8"); return kinds.filter((k) => new RegExp(`(^|[^A-Za-z0-9_])${k}([^A-Za-z0-9_]|$)`).test(t)); });
    check("paramForm.js and chipPicker.js name no registry kind", forms.length === 0, forms.join(","));

    // ---------------------------------------------------------------- save as composed card
    await page.click(".ws-compose .bh-save");
    await page.waitForSelector(".save-form");
    await page.fill(".save-name", "composer-verify-card");
    await page.waitForTimeout(150);
    await page.click(".save-go");
    await page.waitForFunction(() => location.hash.startsWith("#/compose/"), null, { timeout: 8000 }).catch(() => {});
    check("Save as composed card PUTs {name, platform: composed, composed}", putBody && putBody.design && putBody.design.platform === "composed" && putBody.design.name === "composer-verify-card" && putBody.design.composed?.schema === "fsim-composed/1" && putBody.design.composed?.structure?.kind, JSON.stringify(putBody || {}).slice(0, 160));
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 30000 });
    const reopened = await page.$eval('.ws-compose .kind-chip[data-axis="structure"]', (e) => e.dataset.kind);
    check("the saved composed card reopens at #/compose/<name> with the same chain", reopened === putBody.design.composed.structure.kind, reopened);
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip");

    // ---------------------------------------------------------------- papers seed + route
    const papers = await (await fetch(`${BASE}/api/papers`)).json().catch(() => []);
    if (Array.isArray(papers) && papers.length) {
      const d = await (await fetch(`${BASE}/api/papers/${papers[0].id}/design`)).json();
      await page.evaluate((seed) => { sessionStorage.setItem("fsim-studio:compose-seed", JSON.stringify(seed)); }, { design: d.design, paper: papers[0].id, collapsed: d.collapsed, note: d.note });
      await page.goto(`${BASE}/#/design/x`, { waitUntil: "load" });
      await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
      await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
      const note = await page.$eval(".compose-seed", (e) => e.textContent).catch(() => "");
      const left = await page.evaluate(() => sessionStorage.getItem("fsim-studio:compose-seed"));
      check("Open in Composer: the paper design loads, says where it came from, and the seed is consumed", new RegExp(`from paper ${papers[0].id}`).test(note) && left === null, note.slice(0, 80));
    } else check("Open in Composer: paper catalogue available", false, "no papers");
    await page.goto(`${BASE}/#/papers`, { waitUntil: "load" });
    await page.waitForSelector(".ws-papers .pp-plaque, .error-state", { timeout: 30000 }).catch(() => {});
    await page.waitForTimeout(500);
    check("#/papers routes to a workspace", !(await page.$(".error-state")) && (await page.$eval("#workspace", (e) => e.textContent.length)) > 50);

    // ---------------------------------------------------------------- examples row, default load, Fock tier (devlib-p11b)
    const exList = await (await fetch(`${BASE}/api/compose/examples`)).json();
    const viewportKind = () => page.evaluate(() => {
      const v = document.querySelector(".ws-compose .viewport-slot");
      return { canvas: !!v?.querySelector("canvas") && v?.dataset.mode === "webgl", fallback: !!v?.querySelector("[data-mode='fallback2d']") || v?.dataset.mode === "fallback", note: document.querySelector(".ws-compose .vm-note")?.textContent || "" };
    });
    const waitScene = () => page.waitForFunction(() => { const v = document.querySelector(".ws-compose .viewport-slot"); return v && (v.querySelector("canvas") || v.querySelector("[data-mode='fallback2d']")); }, null, { timeout: 30000 }).catch(() => {});
    await page.setViewportSize({ width: 1440, height: 900 });
    await page.goto(`${BASE}/#/design/x`, { waitUntil: "load" });
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "pass", null, { timeout: 60000 }).catch(() => {});
    await waitScene();
    await settle(page, 1200);
    const defKinds = await page.$$eval(".ws-compose .kind-chip", (els) => Object.fromEntries(els.map((e) => [e.dataset.axis, e.dataset.kind])));
    const defEx = exList.find((e) => e.name === "example-pillar-resonant");
    check("default load is example-pillar-resonant (kinds match the example)", !!defEx && Object.keys(defKinds).every((a) => defKinds[a] === defEx.kinds[a]), JSON.stringify(defKinds));
    check("default example is named on the Examples menu button", (await page.$eval(".ws-compose .start-examples", (e) => e.dataset.example)) === "example-pillar-resonant");
    await openStart(page, "examples");
    check("the default example is marked in the Examples menu", (await page.$eval('.start-pop .example-chip[data-example="example-pillar-resonant"]', (e) => e.getAttribute("aria-pressed"))) === "true");
    await closeStart(page);
    let vk = await viewportKind();
    check("default load shows a real 3D builder (canvas, not the fallback plate)", vk.canvas && !vk.fallback, JSON.stringify(vk));
    check("default example validates (plaque VALID)", (await page.$eval(".ws-compose .plaque", (e) => e.dataset.status)) === "pass");
    await page.screenshot({ path: path.join(SHOTS, "composer-micropillar.png") });
    const v3 = await page.$$eval(".ws-compose .v3-readout-label", (els) => els.map((e) => e.textContent));
    check("3D readout labels are notation, not code names ('d um', 'lambda nm')", v3.length > 0 && !v3.some((t) => /\b(um|nm|lambda)\b|_/.test(t)), v3.join(" | "));
    const v3c = await page.evaluate(() => { const i = document.querySelector(".ws-compose .v3-slider input"); return i ? getComputedStyle(i).accentColor : null; });
    const beamCol = await page.evaluate(() => { const probe = document.createElement("i"); probe.style.color = "var(--beam)"; document.body.append(probe); const c = getComputedStyle(probe).color; probe.remove(); return c; });
    check("the 3D cut-away slider is ink, not beam red", v3c && v3c !== beamCol && !/255, 59, 47|210, 40, 30/.test(v3c), String(v3c));
    const notesOpen = await page.evaluate(() => [...document.querySelectorAll(".ws-compose .v3-notes-box")].some((d) => d.open));
    check("the 3D notes panel is collapsed by default", !notesOpen);
    await openStart(page, "examples");
    check("the Examples menu has one entry per example", (await page.$$(".start-pop .example-chip")).length === exList.length && exList.length >= 6, `${exList.length}`);
    await closeStart(page);
    for (const ex of exList) {
      await pickExample(page, ex.name);
      await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "pass", null, { timeout: 60000 }).catch(() => {});
      await waitScene();
      const kinds = await page.$$eval(".ws-compose .kind-chip", (els) => Object.fromEntries(els.map((e) => [e.dataset.axis, e.dataset.kind])));
      const st = await page.$eval(".ws-compose .plaque", (e) => e.dataset.status);
      check(`example chip ${ex.name}: loads its kinds and validates`, Object.keys(ex.kinds).every((a) => kinds[a] === ex.kinds[a]) && st === "pass", `${st}`);
      if (ex.name === "example-open-cavity-jc") {
        await settle(page, 1200);
        vk = await viewportKind();
        check("open-cavity example shows a 3D builder", vk.canvas && !vk.fallback, JSON.stringify(vk));
        await page.screenshot({ path: path.join(SHOTS, "composer-open-cavity.png") });
      }
    }
    await openStart(page, "tiers");
    const fockChip = await page.$('.start-pop .tier-chip[data-tier="fock"]');
    check("Fock tier entry is enabled", !!fockChip && !(await fockChip.isDisabled()));
    if (fockChip) {
      await fockChip.click();
      await page.waitForFunction(() => document.querySelector(".ws-compose .start-examples")?.dataset.example === "", null, { timeout: 10000 }).catch(() => {});
      await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status, null, { timeout: 60000 });
      const fk = await (await fetch(`${BASE}/api/compose/tier/fock`)).json();
      const kinds = await page.$$eval(".ws-compose .kind-chip", (els) => Object.fromEntries(els.map((e) => [e.dataset.axis, e.dataset.kind])));
      check("Fock chip opens the Fock source spec in the Composer", ["drive", "material", "structure", "setup"].every((a) => kinds[a] === fk.design.composed[a].kind), JSON.stringify(kinds));
    }

    // ---------------------------------------------------------------- states for the finish review (captured, and the failed run asserted)
    await page.setViewportSize({ width: 1440, height: 900 });
    await page.goto(`${BASE}/#/overview`, { waitUntil: "load" });
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "pass", null, { timeout: 60000 }).catch(() => {});
    await settle(page, 800);
    await page.click(`.ws-compose .kind-chip[data-axis="structure"]`);
    await page.waitForSelector(".cp-entry");
    await page.waitForFunction(() => !/Checking/.test(document.querySelector(".cp-note")?.textContent || ""), null, { timeout: 60000 });
    const greyed = await page.$$eval(".cp-entry[data-compat='no']", (els) => els.map((e) => { const b = e.querySelector(".cp-badge"); const w = e.querySelector(".cp-why"); return { badge: b ? getComputedStyle(b).display !== "none" && !b.hidden : false, why: !!w && !w.hidden && w.textContent.trim().length > 0 }; }));
    check("chip picker: greyed 'doesn't fit' entries show only the reason (no READY badge)", greyed.length > 0 && greyed.every((g) => !g.badge && g.why), JSON.stringify(greyed));
    await page.screenshot({ path: path.join(SHOTS, "finish-composer-chip-picker.png") });
    await page.keyboard.press("Escape");
    await page.waitForSelector(".cp-entry", { state: "detached", timeout: 5000 });
    // a failed run: the server answers 500 to the run job (mocked), the shelf shows the error plate
    mockedFailure = true;
    await page.route("**/api/compose/run", async (route) => { await route.fulfill({ status: 500, contentType: "application/json", body: JSON.stringify({ error: { kind: "engine", message: "The solver stopped: the Liouvillian did not converge (mocked failure for the capture).", suggestion: "Try the other solver or a coarser pulse." } }) }); });
    await page.click(".ws-compose .run-key");
    await page.waitForSelector(".ws-compose .shelf-alert .error-state", { timeout: 15000 });
    check("a failed run shows the error plate on the shelf with a way to retry", /did not converge/.test(await page.$eval(".ws-compose .shelf-alert", (e) => e.textContent)) && (await page.$$(".ws-compose .shelf-alert .btn")).length >= 1);
    await settle(page, 400);
    await page.screenshot({ path: path.join(SHOTS, "finish-composer-failed-run.png") });
    await page.unroute("**/api/compose/run");
    mockedFailure = false;
    // an invalid chain: a mocked validate answer with an error issue (fresh page, so the failed run above is not on the shelf)
    await page.goto(`${BASE}/#/overview`, { waitUntil: "load" });
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "pass", null, { timeout: 60000 }).catch(() => {});
    await settle(page, 800);
    await page.route("**/api/compose/validate", async (route) => { await route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ ok: false, issues: [{ level: "error", path: "structure.d_um", message: "pillar diameter 0.2 um is below the cutoff of the fundamental mode (mocked invalid chain for the capture)" }, { level: "warn", path: "drive.fwhm_ps", message: "pulse is longer than the cavity lifetime" }], resolved: {}, solver: null, eta_s: null, regime: "unknown" }) }); });
    await page.evaluate(() => { const e = document.querySelector('.ws-compose .mount-body[data-axis="structure"]'); e?.click(); });
    await page.waitForSelector(".ws-compose .sheet .pf-item input.fr-num", { timeout: 5000 });
    const dd = await page.$(".ws-compose .sheet .pf-item .fr-point input.fr-num");
    const dv = await dd.inputValue();
    await dd.fill(String(Number(dv || 1) * 1.05));
    await dd.press("Tab");
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque")?.dataset.status === "fail", null, { timeout: 8000 });
    const badPl = await page.$eval(".ws-compose .plaque", (e) => ({ line: e.classList.contains("plaque-line"), t: e.textContent, run: document.querySelector(".ws-compose .run-key").disabled }));
    check("an invalid chain: INVALID status line lists the reasons and RUN is disabled", badPl.line && /INVALID/.test(badPl.t) && /cutoff/.test(badPl.t) && badPl.run, badPl.t.slice(0, 60));
    const plaqueGrid = () => page.$eval(".ws-compose .plaque", (pl) => {
      const r = (sel) => pl.querySelector(sel)?.getBoundingClientRect();
      const box = pl.getBoundingClientRect(), word = r(".plaque-word span"), den = r(".plaque-denom .v"), reg = r(".plaque-reasons .plaque-kv");
      const ta = (sel) => getComputedStyle(pl.querySelector(sel)).textAlign;
      return { wordL: Math.round(word.left - box.left), solverR: Math.round(box.right - den.right), regL: Math.round(reg.left - box.left), wordTA: ta(".plaque-word"), regTA: ta(".plaque-reasons"), denomTA: ta(".plaque-denom"), w: Math.round(box.width) };
    });
    const gridOk = (g) => g.wordL < 48 && g.solverR <= 16 && g.regL <= 16 && g.wordTA !== "center" && g.regTA !== "center";
    const gInv = await plaqueGrid();
    check("INVALID status line: word left, solver right on the plaque grid, regime left-aligned", gridOk(gInv), JSON.stringify(gInv));
    await settle(page, 400);
    await page.screenshot({ path: path.join(SHOTS, "finish-composer-invalid-chain.png") });
    await page.unroute("**/api/compose/validate");

    // ---------------------------------------------------------------- screenshots (dark + light, desktop + phone)
    for (const theme of ["dark", "light"]) {
      for (const [w, hgt] of [[1440, 900], [390, 844]]) {
        await page.setViewportSize({ width: w, height: hgt });
        await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
        await page.evaluate((t) => { localStorage.setItem("fsim-studio:prefs", JSON.stringify({ theme: t })); }, theme);
        await page.reload({ waitUntil: "load" });
        await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
        await page.waitForFunction(() => document.querySelector(".ws-compose .plaque"), null, { timeout: 60000 });
        await page.click(".ws-compose .run-key").catch(() => {});
        await waitRun(page).catch(() => {});
        await settle(page, 1200);
        if (theme === "dark" && w === 1440) {
          const gRes = await plaqueGrid();
          check("result plaque: word left, solver right, regime/estimate left-aligned (same grid as the status line)", gridOk(gRes) && Math.abs(gRes.solverR - gInv.solverR) <= 4 && Math.abs(gRes.regL - gInv.regL) <= 4, JSON.stringify(gRes));
          const un = await page.evaluate(() => ({
            sheet: [...document.querySelectorAll(".ws-compose .mv-u, .ws-compose .ro-unit")].map((e) => e.textContent),
            dot: !!document.querySelector('.ws-compose .ro-label .sym-dot, .ws-compose [data-readout="F_P_eff"] .sym-dot'),
            fp: document.querySelector('.ws-compose [data-readout="F_P_eff"] .ro-label')?.textContent || "",
            notes: [...document.querySelectorAll(".ws-compose .plaque .pi-t")].map((e) => e.textContent).join(" "),
          }));
          check("units: ueV -> µeV, no 'ueV'/'um' and no bare '1' unit in mounts or readouts", un.sheet.every((t) => !/^(ueV|um|um2|1)$/.test(t)) && un.sheet.some((t) => /µ/.test(t)), JSON.stringify(un.sheet));
          const v3u = await page.$$eval(".ws-compose .v3-readout-unit", (els) => els.map((e) => e.textContent));
          check("3D readouts: units are µm / µeV, no bare '1' unit", v3u.every((t) => !/^(ueV|um|um2|1)$/.test(t)), JSON.stringify(v3u));
          check("notation: F_P L is typeset (KaTeX, subscript P) in the tile label", await page.$eval('.ws-compose [data-readout="F_P_eff"] .ro-label', (e) => !!e.querySelector(".katex .msupsub")), un.fp);
          check("notation: phonon note reads Z² and ×Z (no 'Z^2', no ' x Z')", !/Z\^2| x Z/.test(un.notes) && (!un.notes || /Z²/.test(un.notes)), un.notes.slice(0, 120));
          await page.click(".ws-compose .mount-body[data-axis=\"material\"]", { force: true, timeout: 5000 }).catch((e) => console.log("click err", e.message.slice(0, 200)));
          await page.waitForSelector(".ws-compose .sheet .field-row", { timeout: 5000 }).catch(() => {});
          const su = await page.$$eval(".ws-compose .sheet .fr-unit, .ws-compose .sheet .pf-item .fr-body span:last-child", (els) => els.map((e) => e.textContent));
          check("units: the open sheet shows µeV / µm, never ueV or um", su.some((t) => t === "µeV") && su.every((t) => !/^(ueV|um|um2|1)$/.test(t)), JSON.stringify(su));
          await page.click(".ws-compose .mount-body[data-axis=\"material\"]", { force: true }).catch(() => {});
          await settle(page, 300);
        }
        const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
        if (w === 390) {
          check(`phone width ${theme}: no horizontal page overflow`, overflow <= 1, `${overflow}px`);
          check(`phone width ${theme}: document scrollWidth is exactly 390`, (await page.evaluate(() => document.scrollingElement.scrollWidth)) === 390);
          if (theme === "dark") {
            const grid = await page.$$eval(".ws-compose .mount-cell", (els) => els.map((e) => { const r = e.getBoundingClientRect(); return [Math.round(r.left), Math.round(r.top)]; }));
            check("phone: the chain is a 2 x 2 grid of mounts (no clipped scroller)", new Set(grid.map((g) => g[0])).size === 2 && new Set(grid.map((g) => g[1])).size === 2, JSON.stringify(grid));
            const beamD = await page.$eval(".ws-compose .beam", (e) => e.getAttribute("d"));
            check("phone: the beam is routed through both rows", (beamD.match(/L/g) || []).length >= 5, beamD);
            const hd = await page.evaluate(() => {
              const li = document.querySelectorAll(".board-steps li")[1];
              const bef = li ? getComputedStyle(li, "::before") : null;
              const a = document.querySelector(".start-row .start-btn:not([data-menu='tiers'])")?.getBoundingClientRect();
              const t = document.querySelector(".start-row #start-tiers, .start-row [data-menu='tiers'], .start-row .start-btn:last-child")?.getBoundingClientRect();
              return { dash: bef ? bef.display : "n/a", sameRow: !!(a && t) && Math.abs(a.top - t.top) < 6 };
            });
            check("phone header: step dashes hidden, Tiers menu on the 'Start from:' row", hd.dash === "none" && hd.sameRow, JSON.stringify(hd));
            const navLabels = await page.$$eval(".rail-label", (els) => els.filter((e) => e.offsetParent !== null && getComputedStyle(e).display !== "none").map((e) => e.textContent));
            check("phone: the bottom nav shows icons and only the active label", navLabels.length === 1 && /composer/i.test(navLabels[0]), navLabels.join(","));
            const navBoxes = await page.$$eval(".rail-item", (els) => els.map((e) => e.getBoundingClientRect()).map((r) => [r.left, r.right]));
            check("phone: bottom nav items do not overlap", navBoxes.every((b, i) => i === 0 || b[0] >= navBoxes[i - 1][1] - 1), "");
          }
        }
        await page.screenshot({ path: path.join(SHOTS, `composer-${theme}-${w}.png`), fullPage: w === 390 });
      }
    }
    check("screenshots written (dark + light, 1440x900 and 390x844)", ["dark", "light"].every((t) => [1440, 390].every((w) => fs.existsSync(path.join(SHOTS, `composer-${t}-${w}.png`)))));
    // R2 audit 14: every screenshot THIS run wrote must be a real image (pixel variance, colour share), not just exist
    {
      const { assertPng } = require("./png_check.cjs");
      const t0 = Date.now() - process.uptime() * 1000;
      const badPng = [];
      let nPng = 0;
      for (const dir of [SHOTS]) {
        if (!fs.existsSync(dir)) continue;
        for (const f of fs.readdirSync(dir)) {
          if (!f.endsWith(".png") || f === "failure.png" || !/^(composer|finish-composer)-/.test(f)) continue;
          const full = path.join(dir, f);
          if (fs.statSync(full).mtimeMs < t0) continue;
          nPng++;
          const r = assertPng(full);
          if (!r.ok) badPng.push(`${f}: ${r.reason}`);
        }
      }
      check(`screenshots of this run are non-blank (${nPng} PNGs)`, nPng > 0 && badPng.length === 0, badPng.join(" | "));
    }
    check("no console errors", errors.length === 0, errors.slice(0, 3).join(" | "));
  } catch (e) {
    check("verifier ran to the end", false, `${e.stack || e.message}${errors.length ? ` | page errors: ${errors.slice(0, 3).join(" | ")}` : ""}`);
  } finally {
    await browser.close();
    OWN?.stop();
  }
  const passed = results.filter((r) => r.ok).length;
  console.log(`${passed}/${results.length} composer checks passed`);
  process.exit(passed === results.length ? 0 : 1);
})();
