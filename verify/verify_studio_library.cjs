// FSIM Studio Library workspace checks + screenshots (Playwright, headless Chromium).
//
// Usage: node verify/verify_studio_library.cjs http://127.0.0.1:8772
// The server must already be running (python -m fsim_studio --no-window --port 8772).
// Playwright resolves from node_modules (temp install), then NODE_PATH / global.
// Screenshots saved at 1600x1000.
// Exits 0 iff every check passes; prints "N/N studio library checks passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BASE = (process.argv[2] || "http://127.0.0.1:8772").replace(/\/$/, "");
const SHOTS = path.join(ROOT, ".workers", "studio", "shots", "p2");
const VIEW = { width: 1600, height: 1000 };
const CARDS = ["chatzarakis", "chatzarakis2023", "laferriere", "laferriere2023", "reischle", "kitamura", "zhao",
  "qcap-cavity", "qcap-staged", "qcap-piezo-variant"];

const EXPECT_CLASS = { chatzarakis: "C", chatzarakis2023: "C", laferriere: "M", laferriere2023: "M", reischle: "M", kitamura: "T",
  zhao: "C", "qcap-cavity": "P", "qcap-staged": "P", "qcap-piezo-variant": "P" };

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

async function settle(page, ms = 700) {
  await page.waitForTimeout(ms);
  await page.evaluate(() => document.fonts.ready);
}

async function ready(page, timeout = 45000) {
  try { await page.waitForSelector('#workspace[data-ready="1"]', { timeout }); return true; } catch { return false; }
}

/** Every numeral the Library renders sits in a row / cell that carries a provenance chip. */
async function numberTagAudit(page) {
  return page.evaluate(() => {
    const bad = [];
    let n = 0;
    for (const el of document.querySelectorAll("#workspace .lib-num, #workspace .lib-count")) {
      if (!el.getClientRects().length) continue;
      n++;
      const scope = el.closest("tr, .lib-mix-cell, .lib-ds");
      if (!scope || !scope.querySelector(".tag")) bad.push(el.textContent.trim());
    }
    return { n, bad };
  });
}

(async () => {
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: VIEW, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const consoleErrors = [];
  page.on("console", (m) => { if (m.type() === "error") consoleErrors.push(`${page.url()} ${m.text()}`); });
  page.on("pageerror", (e) => consoleErrors.push(`pageerror: ${e.message}`));
  page.on("requestfailed", (r) => { const u = r.url(); if (!u.includes("/events")) consoleErrors.push(`requestfailed: ${u} ${r.failure()?.errorText}`); });
  const shot = (name) => page.screenshot({ path: path.join(SHOTS, `${name}.png`) });

  try {
    // ------------------------------------------------------------ index
    await page.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
    check("index: ready", await ready(page));
    const plaques = await page.$$eval("a.lib-plaque", (els) => els.map((e) => e.dataset.card));
    check("index: one plaque per parameter card (the 10 shipped cards present)", CARDS.every((c) => plaques.includes(c)), plaques.join(","));
    const groups = await page.$$eval(".lib-group", (els) => els.map((e) => e.getAttribute("aria-label")));
    check("index: cards grouped by paper", groups.length >= 5, groups.join(" | "));
    const mixOk = await page.$$eval("a.lib-plaque", (els) => els.every((e) => e.querySelectorAll(".lib-mix-cell .tag").length === 4));
    check("index: every plaque shows the [V]/[DR]/[E]/[A] mix", mixOk);
    let a = await numberTagAudit(page);
    check("index: every number has a tag chip", a.n > 20 && a.bad.length === 0, `${a.n} numbers; bad ${a.bad.slice(0, 5).join(",")}`);
    const rail = await page.$eval('.rail-item[data-ws="library"]', (e) => ({ cur: e.getAttribute("aria-current"), svg: !!e.querySelector("svg path") }));
    check("rail: Library item (7th) with its glyph, current on #/library", rail.cur === "page" && rail.svg);
    const railIdx = await page.$$eval(".rail-item", (els) => els.findIndex((e) => e.dataset.ws === "library"));
    check("rail: Library is item 7", railIdx === 6, `${railIdx + 1}`);
    const picker = await page.$$eval("#card-picker optgroup", (gs) => gs.map((g) => ({ label: g.label, n: g.querySelectorAll("option").length })));
    const lit = picker.find((g) => g.label === "Literature");
    check("picker: second group 'Literature' lists the parameter cards", picker[0]?.label === "Design" && lit && lit.n >= 10, JSON.stringify(picker));
    await settle(page, 500);
    await shot("library-index");

    // picker routes a literature card to the Library page (not the Designer)
    await page.selectOption("#card-picker", "lib:zhao");
    await page.waitForFunction(() => location.hash === "#/library/zhao", null, { timeout: 10000 }).catch(() => {});
    check("picker: a Literature card routes to #/library/<card>", (await page.evaluate(() => location.hash)) === "#/library/zhao");
    await ready(page);

    // ------------------------------------------------------------ card pages
    const table = [];
    for (const card of CARDS) {
      const api = await page.evaluate(async (n) => (await fetch(`/api/params/${n}`)).json(), card);
      await page.goto(`${BASE}/#/library/${card}`, { waitUntil: "domcontentloaded" });
      const ok = await ready(page);
      check(`${card}: page ready`, ok);
      const hdr = await page.$eval(".lib-head h1", (e) => e.textContent).catch(() => "");
      check(`${card}: header names the card`, hdr === card, hdr);
      a = await numberTagAudit(page);
      check(`${card}: every number has a tag chip`, a.n > 0 && a.bad.length === 0, `${a.n} numbers; bad ${a.bad.slice(0, 5).join(",")}`);
      // ranges never collapsed: each range param renders two numbers inside a bracket span
      const ranges = await page.$$eval(".lib-params tr[data-param]", (trs) => trs.map((tr) => ({ p: tr.dataset.param, nums: [...tr.querySelectorAll(".lib-val .lib-range .lib-num")].length, single: tr.querySelectorAll(".lib-val > .lib-num").length })));
      const want = api.params.filter((p) => Array.isArray(p.range)).map((p) => p.name);
      const rok = want.every((p) => { const r = ranges.find((x) => x.p === p); return r && r.nums === 2 && r.single === 0; });
      check(`${card}: ${want.length} ranges drawn as bracketed spans (two numbers, no midpoint)`, rok, want.join(","));
      // datasets plotted
      const plotted = await page.$$eval(".lib-ds-card", (els) => els.map((e) => ({ ds: e.dataset.dataset, svg: !!e.querySelector(".chart-plot .main-svg"), n: (e.querySelector(".chart-plot")?.data || []).length, ov: e.querySelector(".lib-ov")?.dataset.overlay, traces: (e.querySelector(".chart-plot")?.data || []).map((t) => t.name) })));
      check(`${card}: every dataset plotted (${api.datasets.length})`, plotted.length === api.datasets.length && plotted.every((p) => p.svg && p.n > 0), JSON.stringify(plotted.map((p) => `${p.ds}:${p.n}`)));
      const ovWanted = Object.entries(api.overlays || {}).filter(([, s]) => s.available).map(([d]) => d);
      const ovShown = plotted.filter((p) => p.ov === "shown").map((p) => p.ds);
      check(`${card}: overlays exactly where the rules allow`, ovWanted.length === ovShown.length && ovWanted.every((d) => ovShown.includes(d)), `wanted [${ovWanted}] shown [${ovShown}]`);
      if (ovShown.length) {
        const labelled = plotted.filter((p) => p.ov === "shown").every((p) => p.traces.some((t) => /\((C|M|T|P)\)/.test(t)));
        check(`${card}: overlay traces labelled with their class`, labelled);
      }
      const cls = await page.$eval(".lib-titleline .lib-cls", (e) => ({ code: e.dataset.cls, text: e.textContent })).catch(() => null);
      check(`${card}: class chip shown (${EXPECT_CLASS[card]})`, cls && cls.code === EXPECT_CLASS[card] && cls.text.length > 3, JSON.stringify(cls));
      const ovTexts = await page.$$eval(".lib-ov[data-overlay='shown']", (els) => els.map((e) => e.textContent));
      check(`${card}: each overlay line carries its class label and its calls`,
        ovTexts.length === ovShown.length && ovTexts.every((t) => /Overlay: .*\(([CMTP])\)\./.test(t) && /Calls: \S/.test(t)), ovTexts.map((t) => t.slice(0, 60)).join(" | "));
      if (card === "zhao") {
        const eb = await page.$eval(".lib-ds-card[data-dataset='g2_vs_pump'] .chart-plot", (e) => (e.data || []).some((t) => t.error_y && t.error_y.visible && t.error_y.array.length === 2)).catch(() => false);
        check("zhao: overlay points carry SE error bars", eb);
        // review finding 5: the card's `conflict:` note (3.1 dB [V] vs 0.9 dB in the digest) is rendered prominently
        const cf = await page.evaluate(() => ({
          banner: document.querySelector('[data-role="conflict-banner"]')?.textContent || "",
          row: document.querySelector('tr[data-param="squeezing_dB"] .lib-conflict')?.textContent || "",
          flagged: [...document.querySelectorAll("tr[data-param]")].filter((r) => r.dataset.conflict === "1").map((r) => r.dataset.param),
        }));
        check("zhao: the squeezing_dB source conflict (3.1 dB vs 0.9 dB) shows as a banner above the table and a note on its row; no other row",
          /Source conflict on 1 parameter/.test(cf.banner) && /squeezing_dB/.test(cf.banner) && /3\.1 dB/.test(cf.row) && /0\.9/.test(cf.row) && /source conflict/i.test(cf.row)
          && cf.flagged.length === 1 && cf.flagged[0] === "squeezing_dB", JSON.stringify(cf).slice(0, 160));
      }
      const design = await page.$('a.btn[href^="#/design/"]');
      check(`${card}: 'Open as design' offered only where a design card derives`, !!design === !!api.design, `${api.design || "none"}`);
      await settle(page, 900);
      await shot(`library-${card}`);
      for (const p of plotted) table.push(`${card}/${p.ds}: overlay ${p.ov}`);
    }
    console.log(table.join("\n"));

    // ------------------------------------------------------------ edit copy (UI only; no save from the browser)
    await page.goto(`${BASE}/#/library/chatzarakis`, { waitUntil: "domcontentloaded" });
    await ready(page);
    await page.click("#lib-edit");
    const nIn = await page.$$eval(".lib-params input.lib-in", (e) => e.length);
    const saveTxt = await page.$eval("#lib-save", (e) => e.textContent).catch(() => "");
    check("edit copy: read-only by default, inputs on 'Edit copy', saves to <name>-edited", nIn === 60 && /chatzarakis-edited/.test(saveTxt), `${nIn} inputs; ${saveTxt}`);
  } catch (e) {
    check("harness ran without throwing", false, e.stack || String(e));
  }

  const relevant = consoleErrors.filter((t) => !/favicon/i.test(t));
  check("no console errors", relevant.length === 0, relevant.slice(0, 6).join(" || "));
  await browser.close();
  const passed = results.filter((r) => r.ok).length;
  console.log(`\n${passed}/${results.length} studio library checks passed`);
  console.log(`screenshots: ${SHOTS}`);
  process.exit(passed === results.length ? 0 : 1);
})();
