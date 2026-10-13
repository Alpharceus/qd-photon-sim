// FSIM Studio typeset mathematics: Playwright checks + screenshots (devlib-r2-math).
//
// Usage: node verify/verify_studio_math.cjs http://127.0.0.1:8822     (an already-running server)
//        node verify/verify_studio_math.cjs                           (starts its own on 8822, temp card/cache folders)
// Asserts: KaTeX is vendored and loads with no request leaving the Studio host (offline), `.katex` output on the Composer
// parameter sheet, the KPI tiles and the Papers model-vs-paper table, no raw `$` / backslash / snake_case parameter name in
// the visible text, richText() is XSS-safe (script-ish markup stays text, KaTeX trust is off), Plotly text gets sub/sup HTML.
// Screenshots go to .workers/devlib/shots/r2-tex-<name>-<desk|phone>-<dark|light>.png.
// Playwright resolves from .workers/studio/node/node_modules (NODE_PATH also honoured).
// Prints "N/N studio math checks passed"; exits 0 iff all pass.

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

const settle = async (page, ms = 700) => { await page.waitForTimeout(ms); await page.evaluate(() => document.fonts.ready); };
const setTheme = (page, t) => page.evaluate((x) => { document.documentElement.dataset.theme = x; }, t);

/** Ids of every registry parameter and observable (the snake_case names a user must never read as text). */
let REG_IDS = new Set();
async function loadIds() {
  const reg = await (await fetch(`${BASE}/api/compose/registry`)).json();
  const ids = new Set();
  for (const k of ["materials", "structures", "drives", "setups"]) for (const e of reg[k] || []) for (const p of e.params || []) ids.add(p.name);
  for (const o of reg.observables || []) ids.add(o.id);
  REG_IDS = new Set([...ids].filter((i) => i.includes("_")));
}

/** Visible text nodes under `sel` (not the hidden MathML copy, not inputs/code) that show a raw `$`, a backslash, or a registry id. */
async function rawText(page, sel) {
  return page.evaluate(([s, idList]) => {
    const ids = new Set(idList);
    const root = document.querySelector(s) || document.body;
    const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    const bad = [];
    let n;
    while ((n = w.nextNode())) {
      const v = n.nodeValue;
      if (!v || !v.trim()) continue;
      const el = n.parentElement;
      if (!el || el.closest(".katex-mathml, script, style, code, pre, textarea, option, .id-v, svg")) continue;
      if (!el.getClientRects().length) continue;
      // dotted code paths (fsim_core.compose.validate_spec, dbr.penetration_depth) are references, not parameter names
      const t = v.replace(/[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)+/g, " ");
      let hit = (t.match(/\$|\\/) || [])[0];
      if (!hit) for (const w2 of t.match(/\b[A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+\b/g) || []) if (ids.has(w2)) { hit = w2; break; }
      if (hit) bad.push(`${hit} in "${v.trim().slice(0, 50)}"`);
    }
    return bad;
  }, [sel, [...REG_IDS]]);
}

async function waitRun(page, timeout = 120000) {
  await page.waitForFunction(() => {
    const s = document.querySelector(".ws-compose .status-word");
    return s && ["done", "error", "cancelled"].includes(s.dataset.state);
  }, null, { timeout });
  return page.$eval(".ws-compose .status-word", (el) => el.dataset.state);
}

(async () => {
  if (!BASE) {
    // a Studio already listening on 8822 (the test command starts one for the URL-taking verifiers) is reused; this verifier writes nothing
    const up = await fetch("http://127.0.0.1:8822/api/health").then((r) => r.ok).catch(() => false);
    if (up) BASE = "http://127.0.0.1:8822";
    else { OWN = await startStudioServer(ROOT, 8822); BASE = OWN.base; }
  }
  fs.mkdirSync(SHOTS, { recursive: true });
  await loadIds();
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=d3d11", "--enable-gpu", "--ignore-gpu-blocklist"] });
  const host = new URL(BASE).host;
  const external = [];
  const errors = [];
  const newPage = async (vw = 1440, vh = 900) => {
    const ctx = await browser.newContext({ viewport: { width: vw, height: vh }, deviceScaleFactor: 1 });
    const page = await ctx.newPage();
    page.on("console", (m) => { if (m.type() === "error") errors.push(m.text()); });
    page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
    page.on("request", (r) => {
      const u = r.url();
      if (!u.startsWith("data:") && !u.startsWith("blob:") && new URL(u).host !== host) external.push(u);
    });
    return page;
  };
  const shot = (page, name, vp, theme, full = false) => page.screenshot({ path: path.join(SHOTS, `r2-tex-${name}-${vp}-${theme}.png`), fullPage: full });

  try {
    // ---------------------------------------------------------------- the vendored library
    const asset = async (p) => (await fetch(`${BASE}/${p}`)).status;
    check("vendor/katex: js, css, auto-render and the woff2 fonts are served by the Studio itself",
      (await asset("vendor/katex/katex.min.js")) === 200 && (await asset("vendor/katex/katex.min.css")) === 200
      && (await asset("vendor/katex/contrib/auto-render.min.js")) === 200 && (await asset("vendor/katex/fonts/KaTeX_Main-Regular.woff2")) === 200
      && (await asset("vendor/katex/fonts/KaTeX_Math-Italic.woff2")) === 200);
    const css = await (await fetch(`${BASE}/vendor/katex/katex.min.css`)).text();
    const fonts = [...new Set([...css.matchAll(/url\((fonts\/[^)]+?\.woff2)\)/g)].map((m) => m[1]))];
    const codes = await Promise.all(fonts.map((f) => asset(`vendor/katex/${f}`)));
    check(`every woff2 font the katex css names resolves (${fonts.length})`, fonts.length >= 19 && codes.every((c) => c === 200), `${codes.filter((c) => c !== 200).length} missing`);

    // ---------------------------------------------------------------- Composer
    let page = await newPage();
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.waitForFunction(() => document.querySelector(".ws-compose .plaque"), null, { timeout: 60000 });
    await settle(page);
    const has = await page.evaluate(() => ({ katex: typeof window.katex?.render === "function", plain: document.querySelectorAll(".math.math-plain").length }));
    check("window.katex is loaded and no formula fell back to plain text", has.katex && has.plain === 0, JSON.stringify(has));
    const nTile = await page.$$eval(".ws-compose .readout .ro-label .katex", (e) => e.length);
    check("KPI tiles: symbols are typeset (.katex in the tile labels)", nTile >= 4, `${nTile}`);
    const nMount = await page.$$eval(".ws-compose .mv-k .katex", (e) => e.length);
    check("mount summaries: parameter symbols are typeset", nMount >= 6, `${nMount}`);
    check("no raw $, backslash or registry id in the Composer first view", (await rawText(page, ".ws-compose")).length === 0, (await rawText(page, ".ws-compose")).slice(0, 4).join(" | "));
    await shot(page, "composer", "desk", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "composer", "desk", "light");
    await setTheme(page, "dark");

    await page.click('.ws-compose .mount-body[data-axis="drive"]');
    await page.waitForSelector(".ws-compose .sheet .pf-item", { timeout: 10000 });
    await settle(page, 400);
    const nSheet = await page.$$eval(".ws-compose .pf .fr-label .katex", (e) => e.length);
    const nHelp = await page.$$eval(".ws-compose .pf-help .katex", (e) => e.length);
    check("parameter sheet: every label is a typeset symbol and the help text carries formulas", nSheet >= 8 && nHelp >= 3, `${nSheet} labels, ${nHelp} help formulas`);
    const rs = await rawText(page, ".ws-compose .sheet");
    check("parameter sheet: no raw $, backslash or registry id is visible", rs.length === 0, rs.slice(0, 4).join(" | "));
    const pfItems = await page.$$eval(".ws-compose .pf-item", (els) => els.map((e) => e.querySelector(".fr-label")?.textContent || ""));
    check("parameter sheet: no label shows an underscore or a registry id", pfItems.every((t) => t && !/_/.test(t)), pfItems.slice(0, 3).join(","));
    const tt = await page.$eval(".ws-compose .pf-item .fr-label", (e) => e.title);
    check("a symbol's tooltip names the registry id (the id is findable, just not shown)", /registry id [a-z]/i.test(tt), tt);
    // the sheet search still finds a parameter by its id and by its symbol text
    await page.click('.ws-compose .mount-body[data-axis="material"]');
    await page.waitForSelector('.ws-compose .sheet[data-axis="material"] .pf-filter', { timeout: 10000 });
    await page.fill('.ws-compose .pf-filter', "tau_bare");
    await page.waitForTimeout(250);
    const found = await page.$$eval(".ws-compose .pf-item:not([hidden])", (e) => e.length);
    check("the field filter still matches a registry id", found >= 1, `${found}`);
    await page.fill('.ws-compose .pf-filter', "");
    await page.waitForTimeout(200);
    const sh = await page.$(".ws-compose .sheet");
    await sh.scrollIntoViewIfNeeded();
    await shot(page, "sheet", "desk", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "sheet", "desk", "light");
    await setTheme(page, "dark");

    // run a point so the KPI tiles hold numbers (and captions / why popovers hold formulas)
    await page.click(".ws-compose .run-key");
    const st = await waitRun(page);
    await settle(page, 400);
    check("a point run finishes (the tiles then hold numbers)", st === "done", st);
    const whyBtn = await page.$(".ws-compose .readout[data-readout='rate_detected_Hz'] .ro-why");
    if (whyBtn) {
      await whyBtn.click();
      await page.waitForSelector(".why .why-lede", { timeout: 5000 });
      const nWhy = await page.$$eval(".why .why-lede .katex", (e) => e.length);
      check("the 'why' popover of a KPI typesets its formula", nWhy >= 1, `${nWhy}`);
      await shot(page, "kpi-why", "desk", "dark");
      await page.keyboard.press("Escape");
    } else check("the 'why' popover of a KPI typesets its formula", false, "no why button on the detected-rate tile");
    const rr = await rawText(page, ".ws-compose .shelf");
    check("shelf after a run: no raw $, backslash or registry id", rr.length === 0, rr.slice(0, 4).join(" | "));
    const tiles = await page.$(".ws-compose .readouts");
    await tiles.scrollIntoViewIfNeeded();
    await shot(page, "kpi", "desk", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "kpi", "desk", "light");
    await setTheme(page, "dark");

    // ---------------------------------------------------------------- richText / XSS (the real module in the real page)
    const x = await page.evaluate(async () => {
      const m = await import("/js/ui/math.js");
      window.__xss = 0;
      const host = document.createElement("div");
      host.id = "xss-host";
      document.body.append(host);
      host.append(m.richText(String.raw`<img src=x onerror="window.__xss=1"> <script>window.__xss=2<\/script> $\href{javascript:window.__xss=3}{click}$ $\url{javascript:window.__xss=4}$ $\htmlClass{evil}{x}$ $\includegraphics{x.png}$`));
      const d = document.createElement("div");
      d.append(m.richText("cost \\$5 and $x^2$ and $unclosed"));
      const lit = d.textContent;
      const rich = document.createElement("div");
      rich.append(m.richText("a $\\eta_{\\mathrm{out}}$ b $$\\frac{a}{b}$$ c"));
      return {
        img: host.querySelectorAll("img").length, script: host.querySelectorAll("script").length,
        anchors: host.querySelectorAll("a[href]").length, evilClass: host.querySelectorAll(".evil").length, xss: window.__xss,
        text: host.textContent.includes("<img") && host.textContent.includes("<script>"),
        lit, katex: rich.querySelectorAll(".katex").length, display: rich.querySelectorAll(".katex-display, .math-display").length,
        symbols: [m.symbolFor({ name: "lambda_nm" }), m.symbolFor({ name: "fwhm_ps" }), m.symbolFor({ name: "n_high" }), m.symbolFor({ name: "mode_area" }), m.symbolFor({ name: "x", tex: "\\Omega" })],
        plain: m.plainSymbol("\\gamma^{*}(4\\,\\mathrm K) + \\Gamma_{XX}"),
        plotly: [m.plotlyText("Heatsink temperature T_hs (K)"), m.plotlyText("g²(0)"), m.plotlyText("$\\eta_{\\mathrm{out}}$ at $T_{\\mathrm c}$")],
      };
    });
    check("richText: markup in data stays text (no img, script, onerror, handlers)", x.img === 0 && x.script === 0 && x.xss === 0 && x.text, JSON.stringify({ img: x.img, script: x.script, xss: x.xss }));
    check("richText: KaTeX trust is off (\\href, \\url, \\htmlClass, \\includegraphics create nothing)", x.anchors === 0 && x.evilClass === 0, JSON.stringify({ a: x.anchors, c: x.evilClass }));
    check("richText: \\$ is a literal dollar and an unclosed $ stays text", x.lit === "cost $5 and x2 and $unclosed" || (x.lit.startsWith("cost $5 and ") && x.lit.endsWith("and $unclosed")), x.lit);
    check("richText: inline and display math both reach KaTeX", x.katex >= 2 && x.display >= 1, `${x.katex} katex, ${x.display} display`);
    check("symbolFor: tex wins; the derivation never returns raw snake_case", x.symbols[4] === "\\Omega" && x.symbols.every((s) => !/[A-Za-z]{2,}_[A-Za-z]{2,}/.test(s)), x.symbols.join(" ; "));
    check("plainSymbol gives a Unicode approximation", /γ\*\(4 K\)/.test(x.plain) && /Γ/.test(x.plain), x.plain);
    check("Plotly text: ASCII notation and TeX become <sub>/<sup>/<i> HTML", /T<sub>hs<\/sub>/.test(x.plotly[0]) && /<sup>\(2\)<\/sup>/.test(x.plotly[1]) && /<sub>out<\/sub>/.test(x.plotly[2]), x.plotly.join(" | "));
    await page.evaluate(() => document.getElementById("xss-host")?.remove());
    await page.context().close();

    // ---------------------------------------------------------------- Composer on a phone
    page = await newPage(390, 844);
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await settle(page);
    const hs = await page.evaluate(() => ({ sw: document.documentElement.scrollWidth, cw: document.documentElement.clientWidth }));
    check("phone 390: the Composer has no horizontal page scroll with typeset symbols", hs.sw <= hs.cw + 1, JSON.stringify(hs));
    await shot(page, "composer", "phone", "dark");
    await page.click('.ws-compose .mount-body[data-axis="drive"]');
    await page.waitForSelector(".ws-compose .sheet .pf-item", { timeout: 10000 });
    await settle(page, 400);
    await (await page.$(".ws-compose .sheet")).scrollIntoViewIfNeeded();
    await shot(page, "sheet", "phone", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "sheet", "phone", "light");
    const hs2 = await page.evaluate(() => ({ sw: document.documentElement.scrollWidth, cw: document.documentElement.clientWidth }));
    check("phone 390: the open parameter sheet has no horizontal page scroll", hs2.sw <= hs2.cw + 1, JSON.stringify(hs2));
    await page.context().close();

    // ---------------------------------------------------------------- Papers
    const papers = await (await fetch(`${BASE}/api/papers`)).json();
    check("the paper catalogue is not empty (the Papers checks need one)", papers.length > 0, `${papers.length}`);
    for (const [vw, vh, vp] of [[1440, 900, "desk"], [390, 844, "phone"]]) {
      page = await newPage(vw, vh);
      await page.goto(`${BASE}/#/papers/${encodeURIComponent(papers[0].id)}`, { waitUntil: "load" });
      await page.waitForSelector(".ws-papers .pp-table", { timeout: 60000 });
      await page.waitForFunction(() => document.querySelector(".ws-papers")?.dataset.ready === "1", null, { timeout: 60000 });
      await settle(page);
      const nObs = await page.$$eval(".ws-papers .pp-table .pp-item-obs .katex", (e) => e.length);
      const nRows = await page.$$eval(".ws-papers .pp-table tbody tr", (e) => e.length);
      check(`papers (${vp}): every ledger row names its observable in typeset notation`, nObs >= Math.min(nRows, 2) && nObs >= 1, `${nObs} formulas, ${nRows} rows`);
      const rp = await rawText(page, ".ws-papers .pp-table");
      check(`papers (${vp}): the model-vs-paper table shows no raw $, backslash or registry id`, rp.length === 0, rp.slice(0, 4).join(" | "));
      await (await page.$(".ws-papers .pp-table-panel")).scrollIntoViewIfNeeded();
      await shot(page, "papers", vp, "dark");
      if (vp === "desk") { await setTheme(page, "light"); await settle(page, 400); await shot(page, "papers", vp, "light"); }
      await page.context().close();
    }

    // ---------------------------------------------------------------- Designer + an explorer
    page = await newPage();
    await page.goto(`${BASE}/#/design/edge-inp-gainp-design`, { waitUntil: "load" });
    await page.waitForSelector(".designer .readout, .ws-design .readout, .readout", { timeout: 60000 });
    await settle(page, 1500);
    const nDes = await page.$$eval(".readout .ro-label .katex, .field-row .katex, .chart-card .katex", (e) => e.length);
    const rd = await rawText(page, "#workspace");
    check("designer: typeset symbols on the readouts / fields, no raw $ or backslash", nDes >= 3 && rd.every((t) => !/^[$\\]/.test(t)), `${nDes} formulas; ${rd.slice(0, 3).join(" | ")}`);
    await shot(page, "designer", "desk", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "designer", "desk", "light");
    await page.context().close();

    page = await newPage();
    await page.goto(`${BASE}/#/explain/phonon`, { waitUntil: "load" });
    await page.waitForSelector(".xs-row, .xp-readouts, .readout", { timeout: 60000 });
    await settle(page, 2500);
    const rx = await rawText(page, "#workspace");
    const nx = await page.$$eval("#workspace .katex", (e) => e.length);
    check("explorer (phonon): control labels and readouts are typeset, no raw $ or backslash", nx >= 3 && rx.every((t) => !/^[$\\]/.test(t)), `${nx} formulas; ${rx.slice(0, 3).join(" | ")}`);
    await shot(page, "explorer", "desk", "dark");
    await setTheme(page, "light"); await settle(page, 400);
    await shot(page, "explorer", "desk", "light");
    await page.context().close();

    // ---------------------------------------------------------------- fix2: 3D readouts, notes and numeric fields
    page = await newPage();
    await page.goto(`${BASE}/#/design`, { waitUntil: "load" });
    await page.waitForSelector(".viewport-slot[data-ready='1']", { timeout: 60000 });
    await page.waitForTimeout(1500);
    const v3 = await page.$$eval(".v3-readout-label", (els) => ({ n: els.length, typeset: els.filter((e) => e.querySelector(".katex")).length, raw: els.filter((e) => /[_\$]/.test(e.textContent)).map((e) => e.textContent) }));
    check("3D readout labels with a symbol are KaTeX, none shows an underscore or a backslash", v3.n > 0 && v3.typeset >= 1 && v3.raw.length === 0, JSON.stringify(v3));
    await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
    await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
    await page.click('.ws-compose .kind-chip[data-axis="structure"]');
    await page.waitForSelector(".cp-entry");
    const cp = await page.$$eval(".cp-sum", (els) => els.filter((e) => /(^|[\s(])[A-Za-z]_[A-Za-z0-9]{1,6}(?![\w.])/.test(e.textContent) && !e.querySelector(".katex")).map((e) => e.textContent.slice(0, 60)));
    check("picker summaries go through richText (no raw ASCII subscripts)", cp.length === 0, cp.slice(0, 3).join(" | "));
    await page.keyboard.press("Escape");
    await page.context().close();

    // ---------------------------------------------------------------- ui3: typeset identifiers in picker and sheet text (finding 17)
    {
      page = await newPage();
      await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
      await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
      await settle(page, 800);
      const RAW = String.raw`\b(E_b|tau_[a-z]+|lambda_nm|Omega_0|t_gate)\b`;
      const MORE = String.raw`\+-|\b(Omega|Delta|Lambda|kappa|theta)\b|\bOmega\^2`;     // ASCII stand-ins that must be symbols
      const CODE = String.raw`(qd_gf\.PhononParams|drive_mech|pulse_counting\.pulse_g2|nitride_cavity\.response|fsim_core\.[a-z_]+|DOT_PRESETS)`;
      // visible text nodes (not the MathML copy) of `sel` that show a raw identifier; code names found inside KaTeX output
      const scan = (sel) => page.evaluate(([s, raw, more, code]) => {
        const root = document.querySelector(s);
        const out = { raw: [], more: [], codeInKatex: [], codeSeen: 0, n: 0 };
        if (!root) return out;
        const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
        let n;
        while ((n = w.nextNode())) {
          const v = n.nodeValue, el = n.parentElement;
          if (!v.trim() || !el || el.closest(".katex-mathml") || !el.getClientRects().length) continue;
          out.n++;
          if (new RegExp(raw).test(v)) out.raw.push(v.trim().slice(0, 70));
          if (!el.closest("code") && new RegExp(more).test(v)) out.more.push(v.trim().match(new RegExp(more))[0] + " in " + v.trim().slice(0, 90));
          if (new RegExp(code).test(v)) { out.codeSeen++; if (el.closest(".katex")) out.codeInKatex.push(v.trim().slice(0, 50)); }
        }
        for (const m of root.querySelectorAll("[data-src]")) if (new RegExp(code).test(m.dataset.src)) out.codeInKatex.push(`data-src ${m.dataset.src}`);
        return out;
      }, [sel, RAW, MORE, CODE]);
      const tot = { raw: [], more: [], codeInKatex: [], codeSeen: 0, n: 0, kinds: 0, typeset: 0 };
      const add = (r) => { for (const k of ["raw", "more", "codeInKatex"]) tot[k].push(...r[k]); tot.codeSeen += r.codeSeen; tot.n += r.n; };
      for (const ax of ["drive", "material", "structure", "setup"]) {
        await page.click(`.ws-compose .kind-chip[data-axis="${ax}"]`);
        await page.waitForSelector(".cp-entry", { timeout: 10000 });
        add(await scan(".cp"));
        const kinds = await page.$$eval(".cp-entry", (els) => els.filter((e) => e.getAttribute("aria-disabled") !== "true").map((e) => e.dataset.kind));
        tot.typeset += await page.$$eval(".cp-sum .katex, .cp-label .katex", (e) => e.length);
        await page.keyboard.press("Escape");
        await page.waitForTimeout(250);
        // the open sheet, for every kind of the axis (a kind swap re-renders the sheet)
        for (const k of kinds) {
          await page.click(`.ws-compose .kind-chip[data-axis="${ax}"]`);
          await page.waitForSelector(".cp-entry", { timeout: 10000 });
          await page.click(`.cp-entry[data-kind="${k}"]`);
          await page.waitForFunction(([a, kk]) => document.querySelector(`.ws-compose .kind-chip[data-axis="${a}"]`)?.dataset.kind === kk, [ax, k], { timeout: 10000 });
          await page.waitForTimeout(500);
          // the mount carries .is-open while its sheet is the open one; click the body only when it is not
          if (!(await page.$eval(`.ws-compose .mount[data-axis="${ax}"]`, (m) => m.classList.contains("is-open")))) await page.click(`.ws-compose .mount-body[data-axis="${ax}"]`);
          await page.waitForSelector(`.ws-compose .sheet[data-axis="${ax}"]`, { timeout: 10000 });
          add(await scan(".ws-compose .sheet"));
          tot.kinds++;
        }
      }
      await page.goto(`${BASE}/#/compose`, { waitUntil: "load" });
      await page.waitForSelector(".ws-compose .mount .kind-chip", { timeout: 60000 });
      await page.click('.ws-compose .start-btn[data-menu="examples"]');
      await page.waitForSelector(".start-pop .start-item");
      add(await scan(".start-pop"));
      await page.keyboard.press("Escape");
      check("composer pickers and every kind's open sheet: no E_b / tau_x / lambda_nm / Omega_0 / t_gate as visible text", tot.raw.length === 0 && tot.kinds >= 20, `${tot.kinds} sheets, ${tot.n} text nodes; ${tot.raw.slice(0, 3).join(" | ")}`);
      check("picker and sheet text: no bare Omega/Delta/Lambda/kappa/theta, '+-' or 'Omega^2' left as ASCII", tot.more.length === 0, tot.more.slice(0, 3).join(" | "));
      check("code and module names stay code: none is inside .katex or a typeset data-src", tot.codeSeen > 0 && tot.codeInKatex.length === 0, `${tot.codeSeen} seen; ${tot.codeInKatex.slice(0, 3).join(" | ")}`);
      check("picker summaries are typeset (KaTeX present in the pickers)", tot.typeset >= 5, `${tot.typeset}`);
      await page.context().close();
    }

    // ---------------------------------------------------------------- offline + errors
    check("no request left the Studio host (KaTeX, fonts and Plotly are all local)", external.length === 0, [...new Set(external)].slice(0, 3).join(" "));
    const kerr = errors.filter((e) => !/favicon|ERR_ABORTED|net::ERR_FAILED/.test(e));
    check("no console error or page error during the run", kerr.length === 0, kerr.slice(0, 3).join(" | "));
  } catch (e) {
    check("verifier ran to the end without throwing", false, String(e && e.stack || e).slice(0, 400));
  } finally {
    await browser.close();
    if (OWN) OWN.stop();
  }

  const ok = results.filter((r) => r.ok).length;
  console.log(`${ok}/${results.length} studio math checks passed`);
  process.exit(ok === results.length ? 0 : 1);
})();
