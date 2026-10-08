// FSIM Studio workspaces (Overview, Results, Explain, Story): Playwright smoke + screenshots.
//
// Usage: node verify/verify_studio_workspaces.cjs http://127.0.0.1:8768
// The server must already be running (python -m fsim_studio --no-window --port 8768).
// Playwright is not a repo dependency: resolved from .workers/studio/node/node_modules first.
// Screenshots: .workers/studio/shots/ws/ at 1600x1000 (+ story scenes in the light theme,
// + one 390px mobile overview). Exits 0 iff every check passes; prints "N/N ... passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BASE = (process.argv[2] || "http://127.0.0.1:8768").replace(/\/$/, "");
const SHOTS = path.join(ROOT, ".workers", "studio", "shots", "ws");
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

async function settle(page, ms = 700) {
  await page.waitForTimeout(ms);
  await page.evaluate(() => document.fonts.ready);
}

/** Every number readout on screen (.sro / .readout with a value) must carry a tag chip. */
async function tagAudit(page, scope = "body") {
  return page.evaluate((sel) => {
    const root = document.querySelector(sel) || document.body;
    const ro = [...root.querySelectorAll(".sro, .readout")].filter((e) => e.offsetParent !== null || e.closest(".st-layer"));
    const withValue = ro.filter((e) => e.dataset.state === "value");
    const bad = withValue.filter((e) => !e.querySelector(".tag")).map((e) => e.dataset.key || e.dataset.readout || e.textContent.slice(0, 40));
    return { n: withValue.length, bad };
  }, scope);
}

/** F3: no dotted Plotly line; dashed ([A]/[E]) lines >= 2.5px. */
async function traceAudit(page, scope = "body") {
  return page.evaluate((sel) => {
    const out = { traces: 0, bad: [] };
    for (const gd of document.querySelectorAll(`${sel} .chart-plot`)) {
      for (const t of gd.data || []) {
        out.traces++;
        const ln = t.line || {};
        if (ln.dash === "dot") out.bad.push(`dot:${t.name}`);
        if ((ln.dash === "dash" || ln.dash === "longdash") && !(ln.width >= 2.5)) out.bad.push(`thin:${t.name}:${ln.width}`);
      }
      for (const sh of gd.layout?.shapes || []) {
        if (sh.line?.dash === "dot") out.bad.push(`dot shape ${sh.type}`);
        if ((sh.line?.dash === "dash" || sh.line?.dash === "longdash") && !(sh.line.width >= 2.5)) out.bad.push(`thin shape ${sh.type}`);
      }
    }
    return out;
  }, scope);
}

/** F4: no numeral anywhere is set in B612 Mono. */
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

async function waitReady(page, sel, timeout = 45000) {
  try {
    await page.waitForSelector(`${sel}[data-ready="1"]`, { timeout });
    return true;
  } catch { return false; }
}

(async () => {
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: VIEW, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const consoleErrors = [];
  const watch = (p, tag = "") => {
    p.on("console", (m) => { if (m.type() === "error") consoleErrors.push(`${tag}${m.text()}`); });
    p.on("pageerror", (e) => consoleErrors.push(`${tag}pageerror: ${e.message}`));
    p.on("requestfailed", (r) => {
      const u = r.url();
      if (!u.includes("/events")) consoleErrors.push(`${tag}requestfailed: ${u} ${r.failure()?.errorText}`);
    });
  };
  watch(page);
  const shot = (name) => page.screenshot({ path: path.join(SHOTS, `${name}.png`) });

  try {
    // ================================================================ OVERVIEW
    await page.goto(`${BASE}/#/overview`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector("a.cplaque", { timeout: 30000 });
    const nPlaques = await page.$$eval("a.cplaque", (e) => e.length);
    check("overview: one plaque per campaign (4)", nPlaques === 4, `${nPlaques}`);
    const words = await page.$$eval("a.cplaque", (els) => els.map((e) => `${e.dataset.campaign}:${e.querySelector(".cp-word .vb, .cp-word .vb-pair")?.textContent.trim()}`));
    check("overview: rt_edge plaque reads FAIL", words.some((w) => w.startsWith("rt_edge:") && /FAIL/.test(w)), words.join(" | "));
    const nw = await page.$eval('a.cplaque[data-campaign="nitride_nanowire/full"]', (e) => ({
      word: e.querySelector(".cp-word .vb")?.dataset.status, scope: e.querySelector(".cp-scope")?.textContent || "",
      ticks: [...e.querySelectorAll(".cp-tick")].map((t) => t.dataset.status), mix: e.querySelector(".cp-mix")?.textContent || "",
      counts: [...e.querySelectorAll(".cp-counts .sro")].map((x) => ({ k: x.dataset.key, v: x.querySelector(".sro-value")?.textContent })),
    }));
    check("H9: nanowire plaque leads with the conservative line (by its parsed word), not the most favourable", nw.word === "fail" && /least favourable/.test(nw.scope), JSON.stringify({ w: nw.word, s: nw.scope }));
    check("H9: nanowire line rail lists every line, split lines included", nw.ticks.length === 16 && nw.ticks.includes("split") && /8 × no idealized pass/.test(nw.mix), `${nw.ticks.length} ticks; ${nw.mix}`);
    check("H9: every nanowire count keeps its denominator and coverage is shown", nw.counts.length >= 3 && nw.counts.every((c) => /\//.test(c.v)) && nw.counts.some((c) => c.k === "count:coverage"), JSON.stringify(nw.counts));
    const statusFromWord = await page.evaluate(async () => {
      const { verdictPlaque } = await import("/js/ui/verdictPlaque.js");
      const el = verdictPlaque({ word: "no_idealized_pass", line: 1, fields: { paired_optical_pass: "5", hardware_qualified: "0" } });
      return el.dataset.status;
    });
    check("H9: plaque status comes from the parsed word only (no recomputed split from counts)", statusFromWord === "fail", statusFromWord);
    const plaqueMeta = await page.$$eval("a.cplaque", (els) => els.every((e) => /generated \d{4}-\d\d-\d\d/.test(e.textContent) && e.querySelector(".id-chip .id-v")?.textContent.length >= 7));
    check("overview: every plaque shows generated date + commit", plaqueMeta);
    const best = await page.$$eval(".best-row", (els) => els.map((e) => e.textContent));
    check("overview: BEST strip shows commanded and delivered side by side", best.length >= 2 && best.every((t) => t.includes("commanded") && t.includes("delivered")), `${best.length} rows`);
    check("overview: photons-per-cycle caveat inline", best.some((t) => t.includes("one delivered photon per 707.3 cycles")));
    const stale = await page.$$eval(".ov-coll .stale-badge", (e) => e.length);
    check("overview: pre-audit badges on stale collections", stale >= 2, `${stale}`);
    const heroReady = await waitReady(page, ".hero-viewport");
    check("overview: 3D device hero mounted", heroReady);
    let audit = await tagAudit(page, ".overview");
    check("overview: every number readout has a tag chip", audit.n > 10 && audit.bad.length === 0, `${audit.n} readouts; bad: ${audit.bad.slice(0, 4).join(", ")}`);
    const bo = await b612Audit(page);
    check("F4: no B612 numerals on the Overview", bo.length === 0, bo.slice(0, 5).join(", "));
    await settle(page, 1200);
    await shot("overview-dark");
    // plaque click -> Results
    await page.click('a.cplaque[data-campaign="rt_edge"]');
    await page.waitForSelector(".vt-panel", { timeout: 20000 });
    check("overview: plaque click opens Results", (await page.evaluate(() => location.hash)).startsWith("#/results/rt_edge"));

    // ================================================================ RESULTS
    for (const id of ["rt_edge", "nitride_cavity", "nitride_nanowire/full"]) {
      const slug = id.replace(/\//g, "_");
      // verdict
      await page.goto(`${BASE}/#/results/${id}/verdict`, { waitUntil: "domcontentloaded" });
      await page.waitForSelector(".vt-panel", { timeout: 30000 });
      const nBadges = await page.$$eval(".vt-panel .vb, .vt-panel .vb-pair", (e) => e.length);
      check(`results ${id}: verdict table with badges`, nBadges >= 1, `${nBadges}`);
      await settle(page, 400);
      await shot(`results-${slug}-verdict`);
      // explore
      await page.goto(`${BASE}/#/results/${id}/explore`, { waitUntil: "domcontentloaded" });
      await page.waitForSelector(".ex-chart .chart-plot .main-svg", { timeout: 45000 });
      await settle(page, 900);
      const cnt = await page.$eval(".ex-n", (e) => ({ rows: Number(e.dataset.rows), total: Number(e.dataset.total), text: e.parentElement.textContent }));
      check(`results ${id}: explore draws rows`, cnt.rows > 0, cnt.text);
      // N3: the flux floor is the campaign's own VERDICT flux_floor when it prints one, else /api/gates; source cited
      const fsrc = await page.$eval(".ex-floor-src", (e) => ({ t: e.textContent, o: e.dataset.origin })).catch(() => ({ t: "", o: "" }));
      if (id === "nitride_nanowire/full") {
        check("N3: nanowire explore flux floor comes from its own VERDICT flux_floor, cited file:line",
          fsrc.o === "campaign" && /out\/nitride_nanowire\/full\/results\.md:\d+ \(VERDICT flux_floor=1000\/s\)/.test(fsrc.t), JSON.stringify(fsrc));
      } else {
        check(`N3: ${id} explore flux floor falls back to /api/gates (no VERDICT flux_floor), cited`,
          fsrc.o === "gates" && /out\/rt_edge\/verdict\.md:\d+/.test(fsrc.t) && /this campaign prints no flux_floor/.test(fsrc.t), JSON.stringify(fsrc));
      }
      if (id === "rt_edge") {
        check("results rt_edge: explore shows 768 headline rows by default", cnt.rows === 768 && cnt.total === 768 && /768 headline rows/.test(cnt.text), cnt.text);
        check("results rt_edge: headline model chip visible", await page.$eval(".headline-chip", (e) => e.offsetParent !== null).catch(() => false));
        // H8: the gate box, the ceiling and the floor come from /api/gates
        const gates = await page.evaluate(async () => (await fetch("/api/gates")).json());
        await page.click('.ex-tools .seg[data-value="pareto"]');
        await page.waitForTimeout(1500);
        const refs = await page.$eval(".ex-chart .chart-plot", (gd) => ({
          shapes: (gd.layout.shapes || []).map((s) => ({ t: s.type, x0: s.x0, y0: s.y0, y1: s.y1 })),
          ann: (gd.layout.annotations || []).map((a) => a.text),
        }));
        const rect = refs.shapes.find((s) => s.t === "rect");
        check("H8: explore gate box and floor are the /api/gates values", rect && rect.x0 === gates.flux_floor.value && rect.y1 === gates.g2_ceiling.value && refs.ann.some((a) => a.includes(`flux floor`) && a.includes(`[${gates.flux_floor.tag}]`)), JSON.stringify({ rect, ann: refs.ann.slice(0, 3) }));
        // H18: g2 axes extend above 1 when the data does (CW g2 rows), nothing off-canvas
        await page.click('.ex-tools .seg[data-value="g2T"]');
        await page.waitForTimeout(1200);
        const ax = await page.$eval(".ex-chart .chart-plot", (gd) => ({ top: gd.layout.yaxis.range?.[1], max: Math.max(...gd.data.flatMap((t) => (t.y || []).filter((v) => typeof v === "number"))) }));
        const above = await page.$eval(".ex-count", (e) => e.textContent);
        check("H18: g2 axis extends to include data above 1 and reports it", ax.max > 1 && ax.top >= ax.max && /above g² = 1: axis extended/.test(above), `${JSON.stringify(ax)} · ${above.slice(0, 160)}`);
      }
      const tr = await traceAudit(page, ".res-main");
      check(`results ${id}: no dotted line; dashed lines >= 2.5px`, !tr.bad.length, tr.bad.slice(0, 4).join(","));
      if (id === "nitride_nanowire/full") {
        const leg = await page.$$eval(".ex-chart .legendtext", (e) => e.map((x) => x.textContent).join(" | "));
        // the geometry preset pairs commanded and delivered
        await page.click('.ex-tools .seg[data-value="geom"]');
        await page.waitForTimeout(1200);
        const leg2 = await page.$$eval(".ex-chart .legendtext", (e) => e.map((x) => x.textContent).join(" | "));
        check("results nanowire: commanded and delivered as paired series", /commanded/.test(leg2) && /delivered/.test(leg2), leg2);
        void leg;
      }
      audit = await tagAudit(page, ".res-main");
      check(`results ${id}: explore readouts carry tags`, audit.bad.length === 0, audit.bad.join(","));
      await shot(`results-${slug}-explore`);
      // 3D lattice
      await page.goto(`${BASE}/#/results/${id}/3d`, { waitUntil: "domcontentloaded" });
      const latReady = await waitReady(page, ".lat-viewport", 60000);
      check(`results ${id}: 3D lattice mounted`, latReady);
      if (id === "rt_edge") {
        // I4: the banner counts only gating (headline) rows; non-gating passes say "never gates"
        await page.waitForSelector(".lat-banner", { timeout: 20000 }).catch(() => {});
        const tog = !!(await page.$(".lat-bar label.switch"));
        if (tog) { await page.click(".lat-bar label.switch"); await page.waitForTimeout(2500); }
        const ban = await page.$eval(".lat-banner", (e) => ({ t: e.textContent, g: e.dataset.gating, s: e.dataset.style })).catch(() => ({ t: "", g: "" }));
        check("I4/H1: lattice banner counts gating rows only (0 headline passes), other models 'never gates'", ban.g === "0" && ban.s !== "pass" && (!/non-headline/.test(ban.t) || /never gates/.test(ban.t)), ban.t.slice(0, 160));
        if (tog) { await page.click(".lat-bar label.switch"); await page.waitForTimeout(1500); }
      }
      await settle(page, 1200);
      await shot(`results-${slug}-3d`);
      // figures
      await page.goto(`${BASE}/#/results/${id}/figures`, { waitUntil: "domcontentloaded" });
      await page.waitForSelector(".fig-tile img", { timeout: 20000 });
      await page.waitForFunction(() => [...document.querySelectorAll(".fig-tile img")].slice(0, 4).every((i) => i.complete), null, { timeout: 20000 }).catch(() => {});
      const nFig = await page.$$eval(".fig-tile", (e) => e.length);
      check(`results ${id}: figures gallery`, nFig >= 1, `${nFig}`);
      await settle(page, 500);
      await shot(`results-${slug}-figures`);
    }
    // lightbox
    await page.click(".fig-tile");
    await page.waitForSelector(".lightbox .lb-img", { timeout: 10000 });
    await settle(page, 600);
    await shot("results-lightbox");
    await page.keyboard.press("Escape");
    check("results: lightbox closes on Esc", !(await page.$(".lightbox")));
    // provenance + row drawer with live re-run (rt_edge)
    await page.goto(`${BASE}/#/results/rt_edge/provenance`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".prov-main", { timeout: 20000 });
    await settle(page, 400);
    await shot("results-rt_edge-provenance");
    await page.goto(`${BASE}/#/results/rt_edge/explore`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".ex-row", { timeout: 45000 });
    await page.click(".ex-row");
    await page.waitForSelector(".row-drawer .rd-dl", { timeout: 20000 });
    await page.click(".row-drawer .rd-actions .btn-outline");
    await page.waitForSelector(".row-drawer .cmp-table", { timeout: 60000 });
    const match = await page.$$eval(".row-drawer .cmp-table .match", (e) => e.map((x) => x.textContent));
    check("results rt_edge: live re-run reproduces the CSV row", match.length >= 2 && match.slice(0, 2).every((t) => t.includes("bit-identical")), match.join(" | "));
    const nanKey = await page.evaluate(async () => {
      const { liveNanReason, rtEdgePairs } = await import("/js/workspaces/rowSeed.js");
      const s = { g2_op: null, g2_op__nan_reason: "no crossing in range" };
      const key = rtEdgePairs({}, s)[0][4];
      return { key, reason: liveNanReason(s, key) };
    });
    check("H17: live n/a reason is looked up under the result key (g2_op), not the display label", nanKey.key === "g2_op" && nanKey.reason === "no crossing in range", JSON.stringify(nanKey));
    await settle(page, 400);
    await shot("results-rt_edge-drawer");
    // rerunSeed: "Open card in Designer" carries the row; the Designer shows CSV vs live side by side
    await page.click(".row-drawer .rd-actions .btn-quiet");
    await page.waitForFunction(() => location.hash.startsWith("#/design/edge-"), null, { timeout: 15000 });
    await page.waitForSelector('.seed[data-seeded="1"]', { timeout: 30000 }).catch(() => {});
    await page.waitForFunction(() => document.querySelectorAll(".seed .match .m-ok, .seed .match .m-diff").length >= 2, null, { timeout: 120000 }).catch(() => {});
    const seed = await page.$$eval(".seed .match", (e) => e.map((x) => x.textContent));
    check("rerunSeed: Designer applies the row and shows CSV vs live side by side (bit-identical)", seed.length >= 2 && seed.slice(0, 2).every((t) => t.includes("bit-identical")), seed.join(" | "));
    await settle(page, 600);
    await shot("designer-seeded");

    // ================================================================ EXPLAIN
    for (const t of ["spectral", "cascade", "band", "va-fit"]) {
      await page.goto(`${BASE}/#/explain/${t}`, { waitUntil: "domcontentloaded" });
      if (t === "spectral" || t === "va-fit") {
        await page.waitForSelector(".xp-main .chart-plot .main-svg", { timeout: 30000 });
        await page.waitForSelector('.xp-readouts .readout[data-state="value"]', { timeout: 30000 });
      } else {
        const ok = await waitReady(page, t === "band" ? ".xp-viewport.is-band" : ".xp-viewport", 60000);
        check(`explain ${t}: 3D scene mounted`, ok);
      }
      if (t === "cascade") await page.waitForSelector(".formula .f-res", { timeout: 30000 });
      if (t === "spectral") {
        // move T and confirm eps re-reads from the backend
        const before = await page.$eval('[data-readout="eps"] .ro-value', (e) => e.textContent);
        await page.$eval("#xs-T", (el) => { el.value = "250"; el.dispatchEvent(new Event("input", { bubbles: true })); });
        await page.waitForFunction((b) => document.querySelector('[data-readout="eps"] .ro-value')?.textContent !== b, before, { timeout: 15000 }).catch(() => {});
        const after = await page.$eval('[data-readout="eps"] .ro-value', (e) => e.textContent);
        check("explain spectral: eps follows the T slider (backend)", after !== before, `${before} -> ${after}`);
      }
      audit = await tagAudit(page, ".explain");
      check(`explain ${t}: every number readout has a tag chip`, audit.bad.length === 0, `${audit.n} readouts; bad ${audit.bad.join(",")}`);
      await settle(page, 1200);
      await shot(`explain-${t}`);
    }
    const calLabel = await page.$eval(".calib-label", (e) => e.textContent).catch(() => "");
    check("explain va-fit: labelled calibration (C), not held-out", calLabel.includes("calibration (C), not held-out"), calLabel.slice(0, 80));

    // ================================================================ STORY
    await page.goto(`${BASE}/#/overview`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector("a.cplaque", { timeout: 30000 });
    await page.goto(`${BASE}/#/story/rt-single-photons/1`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".st-title", { timeout: 30000 });
    check("story: Present default is the light theme", (await page.evaluate(() => document.documentElement.dataset.theme)) === "light");
    const nTicks = await page.$$eval(".st-tick", (e) => e.length);
    check("story: progress rail of 13 ticks", nTicks === 13, `${nTicks}`);
    for (let i = 1; i <= 13; i++) {
      if (i > 1) await page.keyboard.press(i % 2 ? "ArrowRight" : "PageDown");
      await page.waitForFunction((n) => location.hash.endsWith(`/${n}`) && document.querySelector(".st-counter")?.textContent.startsWith(String(n).padStart(2, "0")), i, { timeout: 15000 });
      // the visual: chart svg, a ready 3D view, or a board
      await page.waitForFunction(() => {
        const L = [...document.querySelectorAll(".st-layer")].pop();
        if (!L || L.classList.contains("is-in")) return false;
        if (L.querySelector(".board, .st-error")) return true;
        const v3 = [...L.querySelectorAll(".st-viewport")];
        const charts = [...L.querySelectorAll(".chart-plot")];
        return (v3.length ? v3.every((v) => v.dataset.ready === "1") : true) && (charts.length ? charts.every((c) => c.querySelector(".main-svg")) : true) && (v3.length + charts.length > 0);
      }, null, { timeout: 60000 }).catch(() => {});
      const s = await page.evaluate(() => ({
        title: document.querySelector(".st-title")?.textContent || "",
        caveat: document.querySelector(".st-caveat .cv-text")?.textContent || "",
        caveatVisible: (() => { const r = document.querySelector(".st-caveat")?.getBoundingClientRect(); return !!r && r.height > 20 && r.bottom <= window.innerHeight; })(),
        cat: document.querySelector(".st-cat .id-chip .id-v")?.textContent || "",
        commit: [...document.querySelectorAll(".st-cat .id-chip .id-k")].some((k) => k.textContent === "DATA" || k.textContent === "LIVE"),
        err: !!document.querySelector(".st-layer .st-error"),
      }));
      check(`story scene ${i}: claim, pinned caveat, catalog + data commit`, s.title && s.caveat && s.caveatVisible && s.cat && s.commit && !s.err, `${s.title.slice(0, 50)}${s.err ? " [visual error]" : ""}`);
      audit = await tagAudit(page, ".story");
      check(`story scene ${i}: every number readout has a tag chip`, audit.n > 0 && audit.bad.length === 0, `${audit.n} readouts; bad ${audit.bad.join(",")}`);
      const st3 = await traceAudit(page, ".story");
      if (st3.traces) check(`story scene ${i}: no dotted line; dashed lines >= 2.5px`, !st3.bad.length, st3.bad.slice(0, 3).join(","));
      if (i === 7 || i === 13) {
        const fill = await page.evaluate(() => {
          const L = [...document.querySelectorAll(".st-layer")].pop();
          const b = L?.querySelector(".board")?.getBoundingClientRect();
          const p = document.querySelector(".st-proof")?.getBoundingClientRect();
          return b && p ? b.height / p.height : 0;
        });
        check(`F7: scene ${i} grid fills the slide height`, fill >= 0.9, fill.toFixed(2));
        // R2: stretched readout cells carry numerals sized to the cell (value font >= 30% of the cell height)
        const kpi = await page.evaluate(() => {
          const L = [...document.querySelectorAll(".st-layer")].pop();
          return [...L.querySelectorAll(".m2-nums .sro, .sc-line-nums .sro")].map((e) => {
            const v = e.querySelector(".sro-value");
            return { k: e.dataset.key, h: Math.round(e.getBoundingClientRect().height), fs: parseFloat(getComputedStyle(v).fontSize), over: v.scrollWidth > v.clientWidth + 1 };
          });
        });
        const weak = kpi.filter((x) => x.fs / x.h < 0.3 || x.over);
        check(`R2: scene ${i} KPI numerals scale with their cells (font >= 0.3 x cell height, none clipped)`, kpi.length >= 3 && !weak.length,
          `${kpi.length} cells; weak ${weak.slice(0, 3).map((x) => `${x.k}:${x.fs}px/${x.h}px${x.over ? " clipped" : ""}`).join(", ")}`);
      }
      {
        // R4: no raw snake_case status token in the claim; on s13 none inside a readout value either (chips only)
        const tok = await page.evaluate(() => {
          const re = /\b[a-z0-9]+(?:_[a-z0-9]+)+\b/g;
          const strip = (t) => String(t || "").replace(/\S*\/\S*/g, " ");
          const lede = strip(document.querySelector(".st-lede")?.textContent).match(re) || [];
          const L = [...document.querySelectorAll(".st-layer")].pop();
          const vals = [...(L?.querySelectorAll(".sro-value") || []), ...document.querySelectorAll(".st-numbers .sro-value")].filter((v) => !v.querySelector(".vb")).flatMap((v) => v.textContent.match(re) || []);
          return { lede, vals };
        });
        check(`R4: scene ${i} claim and readout values carry no raw snake_case tokens`, !tok.lede.length && !tok.vals.length, JSON.stringify(tok));
      }
      if (i === 13) {
        const plates = await page.$$eval(".sc-plate", (ps) => ps.map((p) => ({ badges: p.querySelectorAll(".vb").length, pairs: p.querySelectorAll(".vb-pair").length })));
        check("F7: one status per platform card on the scorecard (no repeated or compound chip)", plates.length === 3 && plates.every((p) => p.badges === 1 && p.pairs === 0), JSON.stringify(plates));
      }
      if (i === 1) {
        const bs = await b612Audit(page);
        check("F4: no B612 numerals in Story", bs.length === 0, bs.slice(0, 5).join(", "));
      }
      await settle(page, 900);
      await shot(`story-${String(i).padStart(2, "0")}-light`);
    }
    // dark theme sample (T)
    await page.keyboard.press("t");
    await page.waitForFunction(() => document.documentElement.dataset.theme === "dark");
    await settle(page, 900);
    await shot("story-13-dark");
    await page.keyboard.press("t");
    // presenter view (B)
    const [pv] = await Promise.all([page.waitForEvent("popup", { timeout: 10000 }), page.keyboard.press("b")]);
    watch(pv, "presenter: ");
    await pv.waitForSelector("#pv-timer", { timeout: 10000 });
    const pvText = await pv.evaluate(() => document.body.textContent);
    check("story: presenter view opens with timer, notes and next-scene preview", pvText.includes("Notes") && pvText.includes("Next") && pvText.includes("Scene 13"), pvText.slice(0, 80));
    const pvNote = await pv.$eval(".pv-note", (e) => e.textContent).catch(() => "");
    check("F7: presenter view renders the scene's own notes (I5), with no unresolved {placeholder}", pvNote.length > 40 && !/^No presenter notes/.test(pvNote) && !/\{[a-z_0-9]+\}/i.test(pvNote), pvNote.slice(0, 100));
    await pv.setViewportSize({ width: 1180, height: 760 });
    await pv.screenshot({ path: path.join(SHOTS, "story-presenter.png") });
    // R4: every scene's rendered claim + presenter notes in room language (no raw snake_case token)
    const prose = await page.evaluate(async () => {
      const st = await (await fetch("/api/stories/rt-single-photons")).json();
      const re = /\b[a-z0-9]+(?:_[a-z0-9]+)+\b/g;
      const strip = (t) => String(t || "").replace(/\S*\/\S*/g, " ");
      return st.scenes.map((sc) => ({ id: sc.id, bad: [...(strip(sc.claim).match(re) || []), ...(strip(sc.notes).match(re) || [])] })).filter((x) => x.bad.length);
    });
    check("R4: all 13 scenes' claims + presenter notes carry no raw snake_case tokens", prose.length === 0, JSON.stringify(prose).slice(0, 200));
    const s13n = await pv.$eval(".pv-note", (e) => e.textContent).catch(() => "");
    check("R4: s13 presenter notes say 'passes optically, fails in hardware'", /passes optically, fails in hardware/.test(s13n) && !/pass_hardware_infeasible/.test(s13n), s13n.slice(0, 120));
    // presenter captures for every scene (notes visible), driven from the main window's hash
    let pvOk = 0;
    for (let i = 1; i <= 13; i++) {
      await page.evaluate((n) => { location.hash = `#/story/rt-single-photons/${n}`; }, i);
      const ok = await pv.waitForFunction((n) => document.querySelector(".pv-step")?.textContent === `Scene ${n} of 13`
        && (document.querySelector(".pv-note")?.textContent || "").length > 40, i, { timeout: 15000 }).then(() => true).catch(() => false);
      if (ok) pvOk++;
      await pv.waitForTimeout(150);
      await pv.screenshot({ path: path.join(SHOTS, `story-presenter-${String(i).padStart(2, "0")}.png`) });
    }
    check("presenter: notes rendered and captured for all 13 scenes", pvOk === 13, `${pvOk}/13`);
    const pvNums = await pv.evaluate(() => [...document.querySelectorAll(".pv-num-v")].map((e) => e.textContent));
    await page.evaluate(() => { location.hash = "#/story/rt-single-photons/9"; });
    await pv.waitForFunction(() => document.querySelector(".pv-step")?.textContent === "Scene 9 of 13", null, { timeout: 15000 }).catch(() => {});
    const pv9 = await pv.evaluate(() => [...document.querySelectorAll(".pv-num-v")].map((e) => e.textContent));
    check("R4: presenter 'numbers on screen' say status words in room language (s09, s13)",
      [...pvNums, ...pv9].every((t) => !/\b[a-z0-9]+(?:_[a-z0-9]+)+\b/.test(t)) && pv9.some((t) => /passes optically, fails in hardware/.test(t)),
      JSON.stringify(pv9).slice(0, 160));
    await pv.close();
    // E opens the live workspace; Esc returns
    await page.goto(`${BASE}/#/story/rt-single-photons/5`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => document.querySelector(".st-counter")?.textContent.startsWith("05"), null, { timeout: 15000 });
    await page.keyboard.press("e");
    await page.waitForFunction(() => location.hash.startsWith("#/results/rt_edge/explore"), null, { timeout: 15000 });
    await page.waitForSelector(".story-return", { timeout: 10000 });
    await page.waitForSelector(".ex-chart .main-svg", { timeout: 45000 });
    await page.keyboard.press("Escape");
    await page.waitForFunction(() => location.hash === "#/story/rt-single-photons/5", null, { timeout: 15000 }).catch(() => {});
    check("story: E opens the live workspace and Esc returns", (await page.evaluate(() => location.hash)) === "#/story/rt-single-photons/5");
    await page.waitForSelector(".st-title", { timeout: 15000 });
    await page.keyboard.press("Escape");
    await page.waitForFunction(() => !location.hash.startsWith("#/story"), null, { timeout: 10000 }).catch(() => {});
    check("story: Esc exits story mode", !(await page.evaluate(() => location.hash)).startsWith("#/story"));

    // ================================================================ mobile 390px overview
    const mctx = await browser.newContext({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 2 });
    const m = await mctx.newPage();
    watch(m, "mobile: ");
    await m.goto(`${BASE}/#/overview`, { waitUntil: "domcontentloaded" });
    await m.waitForSelector("a.cplaque", { timeout: 30000 });
    await m.waitForTimeout(1500);
    const hscroll = await m.evaluate(() => document.documentElement.scrollWidth > window.innerWidth + 1);
    check("mobile 390px overview: no horizontal page scroll", !hscroll);
    await m.screenshot({ path: path.join(SHOTS, "overview-mobile.png"), fullPage: true });
    await m.evaluate(() => window.scrollTo(0, document.documentElement.scrollHeight));
    await m.waitForTimeout(600);
    await m.screenshot({ path: path.join(SHOTS, "overview-mobile-end.png") });
    await mctx.close();
  } catch (e) {
    check("script ran to completion", false, e.message.split("\n")[0]);
    try { await page.screenshot({ path: path.join(SHOTS, "failure.png") }); } catch { /* ignore */ }
  }

  const relevant = consoleErrors.filter((t) => !/favicon/i.test(t));
  check("no console errors", relevant.length === 0, relevant.slice(0, 6).join(" || "));
  await browser.close();

  const passed = results.filter((r) => r.ok).length;
  console.log(`${passed}/${results.length} studio workspace checks passed`);
  console.log(`screenshots: ${SHOTS}`);
  process.exit(passed === results.length ? 0 : 1);
})();
