// FSIM Studio Papers workspace checks + screenshots (Playwright, headless Chromium).
//
// Usage: node verify/verify_studio_papers.cjs http://127.0.0.1:8773
// The server must already be running (python -m fsim_studio --no-window --port 8773).
// Does not pin any physics number: it checks structure (counts, badges, fills, fonts, charts).
// Env PAPERS_INJECT=1: dev aid for before app.js routes #/papers; rewrites app.js in flight to add the workspace.
// Exits 0 iff every check passes; prints "N/N studio papers checks passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BASE = (process.argv[2] || "http://127.0.0.1:8773").replace(/\/$/, "");
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

async function getJSON(url) {
  const r = await fetch(`${BASE}${url}`);
  if (!r.ok) throw new Error(`${url} -> ${r.status}`);
  return r.json();
}

(async () => {
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const errors = [];
  let mockedFailure = false;   // the failed-reproduce capture answers 500 on purpose
  page.on("console", (m) => { if (m.type() === "error" && !(mockedFailure && /status of 500/.test(m.text()))) errors.push(m.text()); });
  page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
  page.on("requestfailed", (r) => { const u = r.url(); if (!u.includes("/events")) errors.push(`requestfailed: ${u}`); });
  if (process.env.PAPERS_INJECT === "1") {
    await page.route("**/js/app.js", async (route) => {
      const resp = await route.fetch();
      let body = await resp.text();
      if (!body.includes("workspaces/papers.js")) {
        body = body.replace('const WORKSPACES = [', 'const WORKSPACES = [\n  { id: "papers", label: "Papers", icon: "library", load: () => import("./workspaces/papers.js") },');
      }
      await route.fulfill({ response: resp, body });
    });
  }
  const shot = (name) => page.screenshot({ path: path.join(SHOTS, `papers-${name}.png`) });

  try {
    const list = await getJSON("/api/papers");
    check("GET /api/papers returns a list", Array.isArray(list), `${list.length} papers`);

    // pick a paper: prefer one with a series, else the first
    let target = list[0]?.id;
    let hasSeries = false;
    for (const p of list) {
      const d = await getJSON(`/api/papers/${p.id}`);
      if (d.series?.length) { target = p.id; hasSeries = true; break; }
    }
    const detail = target ? await getJSON(`/api/papers/${target}`) : null;

    await page.goto(`${BASE}/#/papers`, { waitUntil: "domcontentloaded" });
    let routed = true;
    try { await page.waitForSelector("#workspace.ws-papers", { timeout: 15000 }); } catch { routed = false; }
    check("route #/papers mounts the Papers workspace", routed, routed ? "" : "route #/papers is not wired in app.js yet");
    if (!routed) throw new Error("Papers workspace not mounted");
    await page.waitForSelector('#workspace[data-ready="1"]', { timeout: 30000 });

    const nPlaques = await page.locator("#workspace .pp-plaque").count();
    check("plaque count == GET /api/papers length", nPlaques === list.length, `${nPlaques} vs ${list.length}`);

    const readOk = await page.evaluate(() => [...document.querySelectorAll(".pp-plaque")].every((p) => p.querySelector(".pp-read") && p.querySelector(".pp-classes") && p.querySelector(".pp-last .ico") && /FS-paper-/.test(p.textContent)));
    check("each plaque has read-level badge, class counts, last-outcome badge with icon, catalog ID", readOk);
    check("the catalogue carries the M / C / T / P key", /unfitted/.test(await page.locator("#workspace .pp-key").innerText()) && /calibration/.test(await page.locator("#workspace .pp-key").innerText()));
    const listWords = await page.$$eval(".pp-plaque .pp-last", (els) => els.map((e) => e.textContent));
    check("no list badge says 'LEDGER MATCH' (ledger agreement is not model-equals-paper)", !listWords.some((t) => /LEDGER MATCH/i.test(t)), listWords.slice(0, 4).join(" | "));
    check("the Designer's card selector is hidden on #/papers", await page.$eval(".ctx-card", (e) => e.hidden && getComputedStyle(e).display === "none"));
    check("'Literature cards (legacy tiers)' is a quiet link in the Papers header", (await page.locator("#workspace .pp-head a.pp-legacy").count()) === 1 && (await page.locator("#workspace .pp-head a.btn").count()) === 0);

    // filter chips narrow the catalogue
    const chip = page.locator("#workspace .pp-filter").first();
    if (await chip.count()) {
      const axis = await chip.getAttribute("data-axis");
      const kind = await chip.getAttribute("data-kind");
      const expect = list.filter((p) => p.kinds[axis] === kind).length;
      await chip.click();
      const n = await page.locator("#workspace .pp-plaque").count();
      check("filter chip narrows the catalogue", n === expect, `${axis}=${kind}: ${n} vs ${expect}`);
      await page.locator("#workspace .pp-clear").click();
      check("clear filters restores the catalogue", (await page.locator("#workspace .pp-plaque").count()) === list.length);
    } else {
      check("filter chips (needs >= 2 kinds on an axis)", list.length < 2, "no chips rendered");
    }

    if (!target) {
      check("catalogue empty: how-to-add plaque shown", (await page.locator('#workspace [data-role="empty"]').count()) === 1);
    } else {
      await page.locator(`#workspace .pp-plaque[data-paper="${target}"]`).click();
      await page.waitForSelector("#workspace .pp-reproduce", { timeout: 20000 });
      await page.waitForSelector('#workspace[data-ready="1"]', { timeout: 30000 });
      check("opening a plaque shows the device header", (await page.locator("#workspace .pp-dev-title").count()) === 1, target);

      await page.locator("#pp-reproduce").click();
      let done = true;
      try { await page.waitForSelector('#workspace .pp-main[data-state="done"]', { timeout: 120000 }); } catch { done = false; }
      const errTxt = done ? "" : await page.locator("#workspace .pp-err").innerText().catch(() => "");
      check("Reproduce job completes", done, errTxt.slice(0, 120));
      if (done) {
        const rowN = await page.locator("#workspace .pp-table tbody tr").count();
        check("table rows == reproduce items", rowN === detail.reproduce.length, `${rowN} vs ${detail.reproduce.length}`);
        const badge = await page.evaluate(() => [...document.querySelectorAll(".pp-table tbody tr")].map((tr) => {
          const o = tr.querySelector(".pp-outcome");
          return !!(o && o.querySelector("svg") && /PASS|IN RANGE ONLY|FAIL|UNAVAILABLE|ERROR/.test(o.textContent));
        }));
        check("every outcome badge shows icon + word", badge.length > 0 && badge.every(Boolean));
        const cls = await page.evaluate(() => [...document.querySelectorAll(".pp-table tbody tr")].every((tr) => /^[CMTPR]/.test(tr.querySelector(".pp-cls")?.textContent || "") && tr.querySelector(".pp-cls")?.title));
        check("every row has a class badge with a tooltip", cls);
        const clsIn = await page.evaluate(() => [...document.querySelectorAll(".pp-table tbody tr")].every((tr) => tr.querySelector("th.pp-item .pp-cls") && !tr.querySelector("td .pp-cls")) && [...document.querySelectorAll(".pp-table tbody .pp-cls")].every((b) => { const bg = getComputedStyle(b).backgroundColor; return bg === "rgba(0, 0, 0, 0)" || bg === "transparent"; }));
        check("the class chip is an outline badge inside the Item cell", clsIn);
        const tagsPerModel = await page.$$eval(".pp-table td.pp-model", (tds) => tds.map((t) => t.querySelectorAll(".tag").length));
        check("one tag chip per model cell", tagsPerModel.length > 0 && tagsPerModel.every((n) => n <= 1), tagsPerModel.join(","));
        const tolTexts = await page.$$eval(".pp-table td.pp-tolc", (tds) => tds.map((t) => t.textContent));
        check("tolerances are plain words, not ledger codes", tolTexts.length > 0 && tolTexts.every((t) => !/\d(sigma)|rel\d|abs\d|\|/.test(t)), tolTexts.slice(0, 4).join(" | "));
        const deltaOk = await page.$$eval(".pp-table td.pp-model .pp-delta-v", (els) => els.map((e) => e.textContent));
        check("delta sits under the model value with a relative percentage", deltaOk.length > 0 && deltaOk.every((t) => /%\)/.test(t) || /^[+−]?[\d.]+$/.test(t)), deltaOk.slice(0, 3).join(" | "));
        const paperNotes = await page.$$eval(".pp-table .pp-where", (els) => els.map((e) => getComputedStyle(e).webkitLineClamp));
        check("the Paper source note wraps (no line clamp)", paperNotes.every((c) => c === "none"), paperNotes.slice(0, 3).join(","));
        const convs = await page.$$eval(".pp-table .pp-obsl", (els) => els.map((e) => e.textContent));
        check("conventions are in words ('raw g(2)(0)'), not ledger ids", convs.length > 0 && convs.every((t) => !/_/.test(t)), convs.slice(0, 4).join(" | "));
        const sumTxt = await page.locator("#workspace .pp-plaque-sum").innerText();
        const sumWord = await page.$eval(".pp-plaque-sum .plaque-word", (e) => e.textContent);
        check("the verdict plaque leads with what was tested: MODEL AGREES ON k OF N TESTED ITEMS", /^\s*MODEL AGREES ON \d+ OF \d+ TESTED ITEM/.test(sumWord) && !/LEDGER MATCHES/.test(sumWord), sumWord);
        const ledgerLine = await page.$eval(".pp-plaque-sum .pp-ledger-line", (e) => ({ t: e.textContent, bg: getComputedStyle(e).backgroundColor }));
        check("the ledger is a secondary outline line with its reason", /Ledger: \d+\/\d+ outcomes as the card expected/.test(ledgerLine.t) && /regression check|human review/.test(ledgerLine.t) && (ledgerLine.bg === "rgba(0, 0, 0, 0)" || ledgerLine.bg === "transparent"), ledgerLine.t.slice(0, 90));
        const plaqueFilled = await page.$eval(".pp-plaque-sum .plaque-word", (e) => ({ bg: getComputedStyle(e).backgroundColor }));
        const dSt = await page.$eval(".pp-plaque-sum", (e) => ({ st: e.dataset.status, ag: Number(e.dataset.agree), te: Number(e.dataset.tested) }));
        check("the plaque status follows the tested result (green only when every tested item agrees)", (dSt.te > 0 && dSt.ag === dSt.te) === (dSt.st === "pass") && (dSt.te === 0) === (dSt.st === "info"), JSON.stringify(dSt) + plaqueFilled.bg);
        await page.waitForSelector(`.pp-plaque[data-paper="${target}"] .pp-last-wrap`, { timeout: 8000 }).catch(() => {});
        const listAfter = await page.$eval(`.pp-plaque[data-paper="${target}"] .pp-last-wrap`, (e) => e.textContent).catch(() => "");
        check("the list badge follows: 'k/N AGREE' with the ledger as a quiet second part", /AGREE\s*\d+\/\d+/.test(listAfter) && /ledger \d+\/\d+/.test(listAfter), listAfter);
        const summaries = await page.$$eval("#workspace summary", (els) => els.filter((e) => e.offsetParent !== null).map((e) => getComputedStyle(e).listStyleType + "|" + getComputedStyle(e, "::after").content));
        check("disclosure markers are chevrons (no native triangles)", summaries.length > 0 && summaries.every((x) => x.startsWith("none|") && !/none$/.test(x)), summaries.slice(0, 4).join(" ; "));
        const mm = await page.evaluate(() => [...document.querySelectorAll(".pp-table tbody tr")].every((tr) => {
          const o = tr.querySelector(".pp-outcome");
          return !o || !o.classList.contains("pp-mismatch") || /LEDGER MISMATCH/.test(o.textContent);
        }));
        check("a ledger mismatch is a labelled split badge", mm);
        // FX4: pass-envelope is its own outline badge and is not counted in "MODEL AGREES ON k OF N"
        const envInfo = await page.evaluate(() => ({
          envRows: [...document.querySelectorAll(".pp-table tbody tr[data-outcome='pass-envelope']")].map((tr) => ({ cls: tr.dataset.cls, txt: tr.querySelector(".pp-outcome")?.textContent || "", filled: getComputedStyle(tr.querySelector(".pp-outcome")).backgroundColor })),
          mPass: [...document.querySelectorAll(".pp-table tbody tr[data-cls='M'][data-outcome='pass']")].length,
          mAll: [...document.querySelectorAll(".pp-table tbody tr[data-cls='M']")].length,
          mEnv: [...document.querySelectorAll(".pp-table tbody tr[data-cls='M'][data-outcome='pass-envelope']")].length,
          ag: Number(document.querySelector(".pp-plaque-sum")?.dataset.agree), te: Number(document.querySelector(".pp-plaque-sum")?.dataset.tested),
          inr: Number(document.querySelector(".pp-plaque-sum")?.dataset.inRange),
          rRows: [...document.querySelectorAll(".pp-table tbody tr[data-cls='R']")].map((tr) => tr.querySelector(".pp-item-kind")?.textContent),
        }));
        check("FX4: the plaque counts only central-within-tolerance M passes; in-range-only items are counted apart",
          envInfo.ag === envInfo.mPass && envInfo.te === envInfo.mAll && envInfo.inr === envInfo.mEnv, JSON.stringify({ ag: envInfo.ag, mPass: envInfo.mPass, te: envInfo.te, inr: envInfo.inr, mEnv: envInfo.mEnv }));
        check("FX4: a pass-envelope outcome shows the outline badge 'IN RANGE ONLY' (never a filled PASS)",
          envInfo.envRows.every((r) => /IN RANGE ONLY/.test(r.txt) && !/^PASS/.test(r.txt.trim()) && (r.filled === "rgba(0, 0, 0, 0)" || r.filled === "transparent")), JSON.stringify(envInfo.envRows));
        check("FX4: class R items are labelled 'paper reference'", envInfo.rRows.every((t) => /paper reference/.test(t || "")), JSON.stringify(envInfo.rRows));
        check("summary plaque shows n / N", /\d+ \/ \d+/.test(await page.locator("#workspace .pp-plaque-sum").innerText()));
        if (hasSeries) {
          await page.waitForSelector("#workspace .pp-ov .js-plotly-plot", { timeout: 30000 }).catch(() => {});
          check("at least one chart rendered (card has series)", (await page.locator("#workspace .pp-ov .js-plotly-plot").count()) >= 1);
          check("chart caption names the parameters source", (await page.locator('#workspace [data-role="params-source"]').first().innerText()).length > 20);
        } else {
          check("card without series shows the no-series note", (await page.locator('#workspace [data-role="no-series"]').count()) === 1);
        }
      }

      // beam fill: only Reproduce
      const beam = await page.evaluate(() => {
        const probe = document.createElement("i");
        probe.style.background = "var(--beam)";
        document.body.append(probe);
        const rgb = getComputedStyle(probe).backgroundColor;
        probe.remove();
        const hits = [];
        for (const el of document.querySelectorAll("#workspace *")) {
          if (!el.getClientRects().length) continue;
          const cs = getComputedStyle(el);
          if (cs.backgroundColor === rgb || /gradient/.test(cs.backgroundImage) && cs.backgroundImage.includes(rgb)) hits.push(el.id || el.className || el.tagName);
        }
        return hits;
      });
      check("only the Reproduce key uses beam fill", beam.length === 1 && String(beam[0]).includes("pp-reproduce"), beam.join(" | "));

      const b612 = await page.evaluate(() => {
        let bad = 0;
        for (const el of document.querySelectorAll("#workspace *")) {
          if (!el.getClientRects().length) continue;
          const own = [...el.childNodes].some((n) => n.nodeType === 3 && /\d/.test(n.textContent));
          if (own && /B612/i.test(getComputedStyle(el).fontFamily)) bad++;
        }
        return bad;
      });
      check("no B612 numerals", b612 === 0, `${b612} elements`);

      const corners = await page.evaluate(() => [...document.querySelectorAll("#workspace .pp-plaque, #workspace .panel, #workspace .run-key, #workspace .btn")].every((e) => parseFloat(getComputedStyle(e).borderTopLeftRadius) <= 2.01));
      check("2 px corners", corners);

      // finish round 2: item cell, header link, axis labels, units
      const cell = await page.$eval(".pp-table tbody tr th.pp-item", (th) => {
        const o = th.querySelector(".pp-item-obs"), i = th.querySelector(".pp-item-id");
        const co = getComputedStyle(o), ci = getComputedStyle(i);
        return { oFirst: th.firstElementChild === o, oW: parseInt(co.fontWeight, 10), iMono: /mono|consolas/i.test(ci.fontFamily), iSize: parseFloat(ci.fontSize), oSize: parseFloat(co.fontSize), id: i.textContent, label: o.textContent };
      });
      check("item cell: human label is the bold heading, the id is a small muted monospace line", cell.oFirst && cell.oW >= 600 && cell.iMono && cell.iSize < cell.oSize && cell.label.trim() && cell.label.trim() !== cell.id, JSON.stringify(cell));
      const leg = await page.$eval(".pp-head-links a", (a) => ({ cls: a.className, text: a.textContent }));
      check("header: 'Literature cards (legacy tiers)' is a quiet link, not a button", leg.text === "Literature cards (legacy tiers)" && !/btn/.test(leg.cls), JSON.stringify(leg));
      if (hasSeries) {
        const ax = await page.evaluate(() => [...document.querySelectorAll(".pp-ov .js-plotly-plot")].map((g) => ({ x: g._fullLayout?.xaxis?.title?.text || "", y: g._fullLayout?.yaxis?.title?.text || "" })));
        check("series charts: axis labels are human words with units, never raw ids", ax.length > 0 && ax.every((a) => a.x && a.y && !/[a-z0-9]_[A-Za-z0-9]/.test(a.x + a.y) && /\(.+\)/.test(a.x)), JSON.stringify(ax));
      }
      const txt = await page.evaluate(() => document.querySelector("#workspace").innerText);
      check("units: no 'ueV' or ' um' in the Papers text", !/ueV|\d\s?um/.test(txt));

      // screenshots: dark/light x desktop/phone
      const setTheme = async (mode) => {
        const cur = await page.evaluate(() => document.documentElement.dataset.theme);
        if (cur !== mode) { await page.click("#theme-btn"); await page.waitForTimeout(600); }
      };
      for (const [vp, size] of [["desktop", { width: 1440, height: 900 }], ["phone", { width: 390, height: 844 }]]) {
        await page.setViewportSize(size);
        for (const mode of ["dark", "light"]) {
          await setTheme(mode);
          await page.waitForTimeout(500);
          await shot(`${mode}-${vp}`);
          if (vp === "phone") {
            const sw = await page.evaluate(() => document.scrollingElement.scrollWidth);
            check(`phone ${mode}: Papers has no horizontal overflow (scrollWidth ${sw})`, sw === 390, String(sw));
          }
        }
        if (vp === "phone") {
          const ph = await page.evaluate(() => {
            const fold = document.querySelector(".pp-cat-fold");
            const main = document.querySelector(".pp-main").getBoundingClientRect();
            const cat = document.querySelector(".pp-cat").getBoundingClientRect();
            const sum = document.querySelector(".pp-cat-sum");
            const cs = getComputedStyle(document.querySelector(".pp-cat"));
            return { open: fold.open, mainTop: main.top, catTop: cat.top, sumShown: sum.offsetParent !== null, sumText: sum.textContent, maxH: cs.maxHeight, ov: cs.overflowY };
          });
          check("phone: 'Papers (n)' summary sits above the chosen paper, closed; the detail is first content below it", !ph.open && ph.catTop < ph.mainTop && ph.sumShown && /Papers \(\d+\)/.test(ph.sumText.replace(/\s+/g, " ")), JSON.stringify(ph));
          check("phone: the catalogue has no nested scroll area", ph.maxH === "none" && ph.ov !== "auto" && ph.ov !== "scroll", `${ph.maxH}/${ph.ov}`);
          const stacked = await page.$eval(".pp-table tbody tr", (tr) => getComputedStyle(tr).display);
          check("phone: the model-vs-paper table is one card per item", stacked === "block", stacked);
          const bottomNav = await page.$$eval(".rail-label", (els) => els.filter((e) => getComputedStyle(e).display !== "none").map((e) => e.textContent));
          check("phone: the bottom nav shows icons and the active label only", bottomNav.length === 1 && /papers/i.test(bottomNav[0]), bottomNav.join(","));
          await page.locator(".pp-table-panel").scrollIntoViewIfNeeded();
          await page.locator(".pp-table-panel").screenshot({ path: path.join(SHOTS, "finish-papers-table-390.png") });
          await page.click(".pp-cat-sum");
          await page.waitForTimeout(300);
          check("phone: tapping 'Papers (n)' opens the catalogue", await page.$eval(".pp-cat-fold", (e) => e.open));
          await page.evaluate(() => window.scrollTo(0, 0));
          await page.screenshot({ path: path.join(SHOTS, "finish-papers-catalogue-open-390.png") });
          await page.click(".pp-cat-sum");
        }
      }
      await setTheme("dark");
      await page.setViewportSize({ width: 1440, height: 900 });

      // a failed reproduce: the job endpoint answers 500 (mocked), the paper page shows the error plate
      mockedFailure = true;
      await page.route("**/api/papers/*/reproduce", async (route) => { await route.fulfill({ status: 500, contentType: "application/json", body: JSON.stringify({ error: { message: "The reproduce job failed: solver stopped (mocked failure for the capture)." } }) }); });
      await page.locator("#pp-reproduce").click();
      await page.waitForSelector("#workspace .pp-err .error-state", { timeout: 10000 });
      check("a failed reproduce shows the error plate with a retry", /mocked failure/.test(await page.locator("#workspace .pp-err").innerText()) && (await page.locator("#workspace .pp-err .btn").count()) >= 1);
      await page.waitForTimeout(300);
      await shot("finish-failed-reproduce");
      await page.unroute("**/api/papers/*/reproduce");
      mockedFailure = false;
      await page.locator("#pp-reproduce").click();
      await page.waitForSelector('#workspace .pp-main[data-state="done"]', { timeout: 120000 });

      // Open in Composer navigates to #/compose
      await page.locator(".pp-composer").click();
      await page.waitForTimeout(800);
      check("Open in Composer navigates to #/compose", page.url().endsWith("#/compose"), page.url());
    }

    // FX4: exercise the pass-envelope / weak / paper-reference rendering on cards that have them (ledger-driven, no physics pinned)
    const fx4 = async (id, want) => {
      const det = await getJSON(`/api/papers/${id}`);
      const led = det.reproduce || [];
      if (!led.some(want.ledger)) return check(`FX4 ${id}: card has the expected ledger entries`, false, "ledger lacks them");
      await page.goto(`${BASE}/#/papers/${id}`, { waitUntil: "domcontentloaded" });
      await page.waitForSelector("#workspace .pp-reproduce", { timeout: 20000 });
      await page.waitForSelector('#workspace[data-ready="1"]', { timeout: 30000 });
      await page.locator("#pp-reproduce").click();
      await page.waitForSelector('#workspace .pp-main[data-state="done"]', { timeout: 120000 });
      const got = await page.evaluate(() => ({
        env: [...document.querySelectorAll(".pp-table tbody tr[data-outcome='pass-envelope']")].map((tr) => ({ cls: tr.dataset.cls, txt: tr.querySelector(".pp-outcome")?.textContent || "", bg: getComputedStyle(tr.querySelector(".pp-outcome")).backgroundColor })),
        weak: [...document.querySelectorAll(".pp-table tbody tr .pp-weak")].map((e) => e.textContent.trim()),
        refs: [...document.querySelectorAll(".pp-table tbody tr[data-cls='R']")].map((tr) => tr.querySelector(".pp-item-kind")?.textContent),
        mPass: document.querySelectorAll(".pp-table tbody tr[data-cls='M'][data-outcome='pass']").length,
        ag: Number(document.querySelector(".pp-plaque-sum")?.dataset.agree), inr: Number(document.querySelector(".pp-plaque-sum")?.dataset.inRange),
        kv: [...document.querySelectorAll(".pp-plaque-sum .plaque-kv")].map((e) => e.textContent),
      }));
      want.assert(got, check, id);
      await shot(`fx4-${id}`);
    };
    await fx4("madsen2014", {
      ledger: (e) => e.expect === "pass-envelope" && e.weak === true || e.weak === true,
      assert: (g, ck, id) => {
        ck(`FX4 ${id}: pass-envelope row shows the outline badge IN RANGE ONLY, not counted in the plaque`, g.env.length >= 1 && g.env.every((r) => /IN RANGE ONLY/.test(r.txt) && (r.bg === "rgba(0, 0, 0, 0)" || r.bg === "transparent")) && g.ag === g.mPass && g.inr === g.env.filter((r) => r.cls === "M").length && g.kv.some((t) => /in range only/.test(t)), JSON.stringify({ env: g.env.length, ag: g.ag, mPass: g.mPass, inr: g.inr }));
        ck(`FX4 ${id}: a weak pass is flagged WEAK PASS`, g.weak.length >= 1 && g.weak.every((t) => /WEAK PASS/.test(t)), g.weak.join("|"));
      },
    });
    await fx4("munsch2013_trumpet", {
      ledger: (e) => e.class === "R",
      assert: (g, ck, id) => {
        ck(`FX4 ${id}: class R rows are labelled 'paper reference' and are not among the tested items`, g.refs.length >= 1 && g.refs.every((t) => /paper reference/.test(t || "")) && g.kv.some((t) => /paper reference/.test(t)), g.refs.join("|"));
      },
    });

    check("no console errors", errors.length === 0, errors.slice(0, 3).join(" | ").slice(0, 300));
  } catch (e) {
    if (!/not mounted/.test(e.message)) check("verifier ran to completion", false, e.message);
  } finally {
    await browser.close();
  }
  const pass = results.filter((r) => r.ok).length;
  console.log(`${pass}/${results.length} studio papers checks passed`);
  process.exit(pass === results.length ? 0 : 1);
})();
