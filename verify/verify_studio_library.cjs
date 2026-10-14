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

// LIBSCENE: count live WebGL contexts (a context whose canvas is attached and not lost); the page polls the maximum.
function glTrack() {
  window.__gl = { list: [], max: 0, created: 0 };
  const orig = HTMLCanvasElement.prototype.getContext;
  HTMLCanvasElement.prototype.getContext = function (type, ...a) {
    const c = orig.call(this, type, ...a);
    if (c && /webgl/.test(String(type)) && !window.__gl.list.includes(c)) { window.__gl.list.push(c); window.__gl.created++; }
    return c;
  };
  setInterval(() => {
    const n = window.__gl.list.filter((c) => c.canvas && c.canvas.isConnected && !c.isContextLost()).length;
    if (n > window.__gl.max) window.__gl.max = n;
  }, 30);
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
  const browser = await chromium.launch({ args: ["--use-angle=d3d11", "--enable-gpu", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: VIEW, deviceScaleFactor: 1 });
  await context.addInitScript(glTrack);
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
    check("rail: Library item (9th, after Compose + Papers) with its glyph, current on #/library", rail.cur === "page" && rail.svg);
    const railIdx = await page.$$eval(".rail-item", (els) => els.findIndex((e) => e.dataset.ws === "library"));
    check("rail: Library is item 9 (Compose + Papers inserted after Designer, 2026-10-12)", railIdx === 8, `${railIdx + 1}`);
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
      const hdrTip = await page.$eval(".lib-head h1", (e) => e.title).catch(() => "");
      check(`${card}: header names the card (a human title; the id is the tooltip)`, hdr.replace(/\s+/g, "").toLowerCase() === card.toLowerCase() && hdrTip === card, `${hdr} / ${hdrTip}`);
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
    // R2 FIN2: human names (no raw ids) for card titles and dataset names; the id stays in the tooltip
    await page.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
    await ready(page);
    const titles = await page.$$eval("a.lib-plaque", (as) => as.map((a) => ({ id: a.dataset.card, t: a.querySelector(".cp-title").textContent, tip: a.querySelector(".cp-title").title, ds: [...a.querySelectorAll(".lib-ds-name")].map((x) => x.textContent) })));
    check("library cards: titles are human names (no digit-glued ids), the id is the tooltip",
      titles.length >= 10 && titles.every((x) => x.t === x.t.replace(/^./, (c) => c.toUpperCase()) && !/[a-z]\d/.test(x.t) && x.tip === x.id), titles.slice(0, 3).map((x) => `${x.id} -> ${x.t}`).join("; "));
    check("library cards: dataset names are labels, not ids (no underscore, no g2_vs_T)", titles.every((x) => x.ds.every((n) => !/_/.test(n) && !/g2_vs/.test(n))), titles.flatMap((x) => x.ds).slice(0, 4).join(" | "));
    // FIN4: every Library card shows a device model or a stated reason (never an empty slot); at most one WebGL canvas
    await page.waitForFunction(() => [...document.querySelectorAll("a.lib-plaque .lib-thumb")].every((t) => t.dataset.state !== "pending"), null, { timeout: 40000 }).catch(() => {});
    const thumbs = await page.$$eval("a.lib-plaque", (as) => as.map((a) => {
      const t = a.querySelector(".lib-thumb");
      const img = t?.querySelector("img.lib-thumb-img"), why = t?.querySelector(".lib-thumb-why");
      return { id: a.dataset.card, state: t?.dataset.state, img: !!img && img.naturalWidth > 0, why: (why?.textContent || "").trim(), shown: !!t && t.getClientRects().length > 0 };
    }));
    check("FIN4: every library card shows a device model or a stated reason (no empty thumbnail slot)",
      thumbs.length >= 10 && thumbs.every((x) => x.shown && ((x.state === "ready" && x.img) || (x.state === "reason" && x.why.length > 20))),
      thumbs.map((x) => `${x.id}:${x.state}`).slice(0, 12).join(" "));
    check("FIN4: the library page holds at most one WebGL canvas", (await page.$$eval("canvas", (c) => c.length)) <= 1);
    await page.goto(`${BASE}/#/library/chatzarakis`, { waitUntil: "domcontentloaded" });
    await ready(page);
    const chartTitles = await page.$$eval(".lib-ds-card h3.chart-title", (e) => e.map((x) => x.textContent));
    check("library card page: dataset chart titles are typeset labels, not ids", chartTitles.length > 0 && chartTitles.every((t) => !/_/.test(t)), chartTitles.join(" | ").slice(0, 160));

    // ------------------------------------------------------------ LIBSCENE: every parameter card shows a 3D model
    {
      await page.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
      await ready(page);
      await page.evaluate(() => { window.__gl.max = 0; });      // count from the index load of this section
      await page.waitForFunction(() => [...document.querySelectorAll("a.lib-plaque .lib-thumb")].every((t) => t.dataset.state !== "pending"), null, { timeout: 60000 }).catch(() => {});
      const ts = await page.$$eval("a.lib-plaque", (as) => as.map((a) => {
        const t = a.querySelector(".lib-thumb"), img = t?.querySelector("img.lib-thumb-img");
        return { id: a.dataset.card, state: t?.dataset.state, ok: !!img && img.naturalWidth > 0, bytes: img ? img.src.length : 0 };
      }));
      const mine = ts.filter((x) => CARDS.includes(x.id));
      check("LIBSCENE: all 10 Library cards show a device model thumbnail (state ready, a decoded image, not blank), none a stated reason",
        mine.length === 10 && mine.every((x) => x.state === "ready" && x.ok && x.bytes > 8000), mine.map((x) => `${x.id}:${x.state}:${x.bytes}`).join(" "));
      const gl1 = await page.evaluate(() => ({ max: window.__gl.max, created: window.__gl.created, live: window.__gl.list.filter((c) => c.canvas && c.canvas.isConnected && !c.isContextLost()).length }));
      check("LIBSCENE: the index (10 thumbnails) held at most one live WebGL context at any time, and none remain",
        gl1.max <= 1 && gl1.live === 0 && gl1.created >= 10, JSON.stringify(gl1));

      // FIN5: scroll the whole list at 1440x900 (every card passes through the viewport), wait (bounded) for the thumbnails to settle:
      // afterwards EVERY card shows a rendered picture or a stated reason, and the page never held more than one live WebGL context.
      {
        await page.setViewportSize({ width: 1440, height: 900 });
        await page.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
        await ready(page);
        await page.evaluate(() => { window.__gl.max = 0; });
        const nCards = await page.$$eval("a.lib-plaque", (as) => as.length);
        for (let i = 0; i < nCards; i++) {
          await page.evaluate((k) => document.querySelectorAll("a.lib-plaque")[k]?.scrollIntoView({ block: "center" }), i);
          await page.waitForTimeout(120);
        }
        await page.waitForFunction(() => [...document.querySelectorAll("a.lib-plaque .lib-thumb")].every((t) => t.dataset.state !== "pending"), null, { timeout: 90000 }).catch(() => {});
        const seen = [];
        for (let i = 0; i < nCards; i++) {
          await page.evaluate((k) => document.querySelectorAll("a.lib-plaque")[k]?.scrollIntoView({ block: "center" }), i);
          await page.waitForTimeout(60);
          seen.push(await page.evaluate((k) => {
            const a = document.querySelectorAll("a.lib-plaque")[k], t = a.querySelector(".lib-thumb");
            const img = t?.querySelector("img.lib-thumb-img, canvas"), why = (t?.querySelector(".lib-thumb-why")?.textContent || "").trim();
            const r = (img || t?.querySelector(".lib-thumb-why"))?.getBoundingClientRect();
            const vh = window.innerHeight;
            return { id: a.dataset.card, state: t?.dataset.state, rendered: !!img && (img.tagName === "CANVAS" || (img.complete && img.naturalWidth > 0)),
              why: why.length, inView: !!r && r.width > 0 && r.height > 0 && r.bottom > 0 && r.top < vh };
          }, i));
        }
        const glF = await page.evaluate(() => ({ max: window.__gl.max, live: window.__gl.list.filter((c) => c.canvas && c.canvas.isConnected && !c.isContextLost()).length }));
        check("FIN5: after scrolling the whole Library list at 1440x900, every card shows a rendered thumbnail or a stated reason (none blank)",
          nCards >= 10 && seen.every((x) => x.inView && ((x.state === "ready" && x.rendered) || (x.state === "reason" && x.why > 20))), seen.map((x) => `${x.id}:${x.state}`).join(" "));
        check("FIN5: the Library list held at most one live WebGL context while its thumbnails rendered", glF.max <= 1, JSON.stringify(glF));
      }

      const SH = path.join(ROOT, ".workers", "devlib", "libscene");
      fs.mkdirSync(SH, { recursive: true });
      const VPS = [["dark", 1440, 900], ["light", 1440, 900], ["dark", 390, 844]];
      for (const [theme, w, hgt] of VPS) {
        const tag = `${theme}-${w}`;
        const ctx2 = await browser.newContext({ viewport: { width: w, height: hgt }, deviceScaleFactor: 1 });
        await ctx2.addInitScript(glTrack);
        await ctx2.addInitScript((th) => { try { localStorage.setItem("fsim-studio:prefs", JSON.stringify({ theme: th })); } catch { /* blocked */ } }, theme);
        const pg = await ctx2.newPage();
        const errs = [];
        pg.on("console", (m) => { if (m.type() === "error") errs.push(m.text()); });
        pg.on("pageerror", (e) => errs.push(`pageerror: ${e.message}`));
        try {
          await pg.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
          await ready(pg);
          await pg.waitForFunction(() => [...document.querySelectorAll("a.lib-plaque .lib-thumb")].every((t) => t.dataset.state !== "pending"), null, { timeout: 60000 }).catch(() => {});
          const th = await pg.$$eval("a.lib-plaque", (as) => as.filter((a) => a.querySelector(".lib-thumb")?.dataset.state === "ready").length);
          check(`LIBSCENE ${tag}: the index shows ${CARDS.length} model thumbnails`, th >= CARDS.length, String(th));
          await pg.screenshot({ path: path.join(SH, `index-${tag}.png`), fullPage: true });
          for (const card of CARDS) {
            const want = await pg.evaluate(async (n) => (await fetch(`/api/scene/device?card=${n}`)).json(), card);
            await pg.goto(`${BASE}/#/library/${card}`, { waitUntil: "domcontentloaded" });
            await ready(pg);
            const live = await pg.waitForFunction(() => document.querySelector('[data-role="device-model"]')?.dataset.state === "live", null, { timeout: 45000 }).then(() => true).catch(() => false);
            await pg.waitForFunction(() => document.querySelector('[data-role="device-model"] .dm-view')?.dataset.ready === "1", null, { timeout: 20000 }).catch(() => {});
            await settle(pg, 600);
            const st = await pg.evaluate(() => {
              const m = document.querySelector('[data-role="device-model"]'), v = m?.querySelector(".dm-view");
              const hudTxt = (m?.textContent || "");
              return { state: m?.dataset.state, variant: m?.dataset.variant, mode: v?.dataset.mode, ready: v?.dataset.ready,
                canvases: m ? m.querySelectorAll("canvas").length : 0, hud: hudTxt.includes("display-only scene"), tagline: hudTxt.includes("scene numbers taken from the card"),
                gl: window.__gl.list.filter((c) => c.canvas && c.canvas.isConnected && !c.isContextLost()).length,
                hscroll: document.documentElement.scrollWidth - document.documentElement.clientWidth, box: m ? m.getBoundingClientRect().height : 0 };
            });
            check(`LIBSCENE ${tag} ${card}: detail view mounts the live 3D model (WebGL, ready, builder ${want.builder}, one canvas, HUD names the card note and the tag summary)`,
              live && st.mode === "webgl" && st.ready === "1" && st.variant === want.variant && st.canvases === 1 && st.hud && st.tagline && st.box > 200,
              JSON.stringify(st));
            check(`LIBSCENE ${tag} ${card}: at most one live WebGL context on the card page, no horizontal scroll`, st.gl <= 1 && st.hscroll <= 1, `gl ${st.gl}, hscroll ${st.hscroll}`);
            const el = await pg.$('[data-role="device-model"]');
            const file = path.join(SH, `${card}-${tag}.png`);
            if (el) await el.screenshot({ path: file });
            const bytes = el ? fs.statSync(file).size : 0;
            check(`LIBSCENE ${tag} ${card}: the model screenshot is not blank (PNG > 15 KB)`, bytes > 15000, String(bytes));
            if (card === CARDS[0]) await pg.screenshot({ path: path.join(SH, `page-${card}-${tag}.png`) });
          }
          const mx = await pg.evaluate(() => window.__gl.max);
          check(`LIBSCENE ${tag}: across the index and all 10 card pages at most one WebGL context was live at any time`, mx <= 1, String(mx));
          check(`LIBSCENE ${tag}: no console errors`, errs.filter((t) => !/favicon/i.test(t)).length === 0, errs.slice(0, 3).join(" | "));
        } catch (e) {
          check(`LIBSCENE ${tag}: harness`, false, e.stack || String(e));
        } finally {
          await ctx2.close();
        }
      }
      await page.goto(`${BASE}/#/library`, { waitUntil: "domcontentloaded" });
    }
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
