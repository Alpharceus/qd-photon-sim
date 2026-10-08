// FSIM Studio model explorers A (cw-g2, pulse-counting, lindblad): Playwright checks + screenshots.
//
// Usage: node verify/verify_studio_explore_a.cjs http://127.0.0.1:8781
// The server must already be running (python -m fsim_studio --no-window --port 8781).
// Screenshots saved at 1600x1000.
// Exits 0 iff every check passes; prints "N/N studio explorer A web checks passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BASE = (process.argv[2] || "http://127.0.0.1:8781").replace(/\/$/, "");
const SHOTS = path.join(ROOT, ".workers", "studio", "shots", "p2");
const VIEW = { width: 1600, height: 1000 };

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

async function until(page, fn, arg, timeout = 30000) {
  try { await page.waitForFunction(fn, arg, { timeout }); return true; } catch { return false; }
}

/** Set a slider by data-ctl id; `value` is in real units (log sliders are mapped here, as the page does). */
async function setSlider(page, ctl, value, { log = null } = {}) {
  await page.evaluate(({ ctl, value, log }) => {
    const el = document.querySelector(`input.xs-range[data-ctl="${ctl}"]`);
    if (!el) throw new Error(`no slider ${ctl}`);
    el.value = log ? String(Math.log(value / log[0]) / Math.log(log[1] / log[0]) * 1000) : String(value);
    el.dispatchEvent(new Event("input", { bubbles: true }));
  }, { ctl, value, log });
}

const readoutText = (page, key) => page.evaluate((k) => {
  const el = document.querySelector(`.readout[data-readout="${k}"]`);
  return el ? { state: el.dataset.state, value: el.querySelector(".ro-value")?.textContent, qual: el.querySelector(".ro-qual")?.textContent, tag: el.querySelector(".tag")?.dataset.tag || null } : null;
}, key);

async function tagAudit(page) {
  return page.evaluate(() => {
    const bad = []; let n = 0;
    for (const el of document.querySelectorAll(".xp-body .readout")) {
      if (!el.getClientRects().length || el.dataset.state !== "value") continue;
      n++;
      if (!el.querySelector(".tag")) bad.push(el.dataset.readout);
    }
    let rows = 0;
    for (const el of document.querySelectorAll(".xp-body .xs-row[data-ctl]")) {
      if (!el.getClientRects().length) continue;
      rows++;
      if (!el.querySelector(".tag")) bad.push(`slider:${el.dataset.ctl}`);
    }
    return { n, rows, bad };
  });
}

async function chartReady(page, n = 1) {
  return until(page, (k) => document.querySelectorAll(".xp-body .chart-plot .main-svg").length >= k, n, 30000);
}

(async () => {
  fs.mkdirSync(SHOTS, { recursive: true });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const context = await browser.newContext({ viewport: VIEW, deviceScaleFactor: 1 });
  const page = await context.newPage();
  const errors = [];
  page.on("console", (m) => { if (m.type() === "error") errors.push(`${page.url()} ${m.text()}`); });
  page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
  page.on("requestfailed", (r) => { const u = r.url(); if (!u.includes("/events")) errors.push(`requestfailed: ${u} ${r.failure()?.errorText}`); });
  const shot = (name) => page.screenshot({ path: path.join(SHOTS, `explore-${name}.png`) });

  try {
    // ------------------------------------------------------------ the tabs list the new topics
    await page.goto(`${BASE}/#/explain/cw-g2`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".xp-tabs", { timeout: 30000 });
    const tabs = await page.$$eval(".xp-tabs a", (els) => els.map((e) => e.getAttribute("href")));
    check("explain tabs list the three explorers", ["#/explain/cw-g2", "#/explain/pulse-counting", "#/explain/lindblad"].every((h) => tabs.includes(h)), tabs.join(" "));

    // ------------------------------------------------------------ cw_g2
    check("cw-g2: both charts render", await chartReady(page, 2));
    check("cw-g2: readouts have values", await until(page, () => document.querySelectorAll('.xp-body .readout[data-state="value"]').length >= 6));
    const capt = await page.$eval('.xa-caveat[data-caveat="caveat"]', (e) => e.textContent).catch(() => "");
    check("cw-g2: pinned caveat visible (Rabi drive, 500 ps HBT)", capt.includes("not valid for resonant (Rabi) drive") && capt.includes("g2_raw is what a 500 ps HBT measures"), capt.slice(0, 60));
    const g0 = await readoutText(page, "g2_dot0");
    check("cw-g2: g2_dot(0) default reads 0.1474 with a tag chip", g0 && g0.value === "0.1474" && !!g0.tag, JSON.stringify(g0));
    const audit1 = await tagAudit(page);
    check("cw-g2: every readout and slider carries a tag chip", audit1.bad.length === 0 && audit1.n >= 6 && audit1.rows >= 7, `${audit1.n} readouts, ${audit1.rows} sliders, bad ${audit1.bad.join(",")}`);
    const rawBefore = (await readoutText(page, "g2_raw0")).value;
    await setSlider(page, "irf_fwhm_ps", 100);
    check("cw-g2: IRF slider updates g2_raw(0)", await until(page, (b) => document.querySelector('.readout[data-readout="g2_raw0"] .ro-value')?.textContent !== b, rawBefore), `${rawBefore} -> ${(await readoutText(page, "g2_raw0")).value}`);
    await setSlider(page, "r", 20, { log: [0.001, 100] });
    check("cw-g2: r = 20 flags cascade bunching and shows the note", await until(page, () => document.querySelector('.readout[data-readout="g2_dot0"] .ro-qual')?.textContent === "cascade bunching")
      && (await page.$eval('.xa-caveat[data-caveat="caveat"]', (e) => e.textContent)).includes("cascade bunching of the leaked XX partner"));
    await page.click('label.switch');
    check("cw-g2: thermal escape switch fills k_X", await until(page, () => document.querySelector('.readout[data-readout="k_X"]')?.dataset.state === "value"), JSON.stringify(await readoutText(page, "k_X")));
    await setSlider(page, "r", 0.5, { log: [0.001, 100] });
    await setSlider(page, "irf_fwhm_ps", 500);
    await page.click('label.switch');
    await until(page, () => document.querySelector('.readout[data-readout="g2_dot0"] .ro-value')?.textContent === "0.1474");
    await page.waitForTimeout(700);
    await shot("cw_g2");

    // ------------------------------------------------------------ pulse_counting
    await page.goto(`${BASE}/#/explain/pulse-counting`, { waitUntil: "domcontentloaded" });
    check("pulse-counting: three charts render", await chartReady(page, 3));
    check("pulse-counting: readouts have values", await until(page, () => document.querySelectorAll('.xp-body .readout[data-state="value"]').length >= 8));
    const capp = await page.$eval('.xa-caveat[data-caveat="caveat"]', (e) => e.textContent).catch(() => "");
    check("pulse-counting: pinned caveat visible (rectangular pump, long-delay peak, f1b gap, no background)",
      capp.includes("Rectangular pump [A]") && capp.includes("g2_adj is the adjacent-normalised") && capp.includes("re-excitation") && capp.includes("rho is a device.py addition"), capp.slice(0, 60));
    const gp = await readoutText(page, "g2");
    check("pulse-counting: g2 default reads a number with a tag chip", gp && gp.state === "value" && !!gp.tag, JSON.stringify(gp));
    const audit2 = await tagAudit(page);
    check("pulse-counting: every readout and slider carries a tag chip", audit2.bad.length === 0 && audit2.n >= 8 && audit2.rows >= 7, `${audit2.n} readouts, ${audit2.rows} sliders, bad ${audit2.bad.join(",")}`);
    const gBefore = gp.value;
    await setSlider(page, "tau_on", 1.0, { log: [0.01, 2] });
    check("pulse-counting: tau_on slider updates g2", await until(page, (b) => document.querySelector('.readout[data-readout="g2"] .ro-value')?.textContent !== b, gBefore), `${gBefore} -> ${(await readoutText(page, "g2")).value}`);
    await setSlider(page, "tau_on", 0.1, { log: [0.01, 2] });
    await until(page, (b) => document.querySelector('.readout[data-readout="g2"] .ro-value')?.textContent === b, gBefore);
    await page.click('.segmented[aria-label="Counting gate"] button[data-value="auto"]');
    check("pulse-counting: gate switch updates the counting window", await until(page, (b) => document.querySelector('.readout[data-readout="mean_counts"] .ro-value')?.textContent !== b, (await readoutText(page, "mean_counts")).value));
    await page.click('.segmented[aria-label="Counting gate"] button[data-value="none"]');
    const gate = await page.evaluate(() => (document.querySelector(".chart-plot")?.layout?.annotations || []).some((a) => /gate .*\/api\/gates/.test(a.text || "")));
    check("pulse-counting: the gate line is drawn with its /api/gates label", gate);
    await page.waitForTimeout(800);
    await shot("pulse_counting");

    // ------------------------------------------------------------ lindblad
    await page.goto(`${BASE}/#/explain/lindblad`, { waitUntil: "domcontentloaded" });
    check("lindblad: incoherent panel renders", await chartReady(page, 1));
    check("lindblad: readouts have values", await until(page, () => document.querySelectorAll('.xp-body .readout[data-state="value"]').length >= 3));
    const capl = await page.$eval('.xa-caveat[data-caveat="caveat"]', (e) => e.textContent).catch(() => "");
    check("lindblad: pinned caveat visible (Markovian baths, Lorentzian ZPL)", capl.includes("Markovian baths, Lorentzian ZPL"), capl.slice(0, 60));
    const md = await readoutText(page, "max_rel_diff");
    check("lindblad incoherent: the agreement readout is below 1e-10 and says agrees",
      md && md.state === "value" && /×10⁻[¹²³⁴⁵⁶⁷⁸⁹][⁰¹²³⁴⁵⁶⁷⁸⁹]/.test(md.value) && (await page.$eval('.readout[data-readout="max_rel_diff"] .ro-caption', (e) => e.textContent)).includes("agrees"), JSON.stringify(md));
    const st1 = await page.$eval(".xp-body .xa-check", (e) => e.textContent);
    check("lindblad incoherent: state check line is shown and physical", st1.includes("physical") && !st1.includes("NOT physical"), st1.slice(0, 80));
    const a3 = await tagAudit(page);
    check("lindblad incoherent: every readout and slider carries a tag chip", a3.bad.length === 0 && a3.n >= 4, `${a3.n} readouts, ${a3.rows} sliders, bad ${a3.bad.join(",")}`);
    await setSlider(page, "r", 20, { log: [0.001, 100] });
    check("lindblad incoherent: r slider updates the Lindblad g2(0)", await until(page, () => document.querySelector('.readout[data-readout="g2_0_lindblad"] .ro-value')?.textContent !== "0.1474"));
    await setSlider(page, "r", 0.5, { log: [0.001, 100] });
    await page.waitForTimeout(500);
    await shot("lindblad-incoherent");
    for (const [panel, ro1, extra] of [["rabi", "rho_ee", "rho_ee"], ["hom", "indistinguishability", null], ["filter", "g2_filtered_0", null]]) {
      await page.click(`.segmented[aria-label="Lindblad panel"] button[data-value="${panel}"]`);
      check(`lindblad ${panel}: panel renders with a tagged readout`, await until(page, (k) => document.querySelector(`.xp-split[data-panel="${k.p}"] .readout[data-readout="${k.r}"]`)?.dataset.state === "value", { p: panel, r: ro1 }, 30000)
        && !!(await readoutText(page, ro1)).tag, JSON.stringify(await readoutText(page, ro1)));
      const a = await tagAudit(page);
      check(`lindblad ${panel}: every readout and slider carries a tag chip`, a.bad.length === 0, `${a.n} readouts, ${a.rows} sliders, bad ${a.bad.join(",")}`);
      if (panel === "rabi") {
        const before = (await readoutText(page, "rho_ee")).value;
        await setSlider(page, "omega", 6);
        check("lindblad rabi: Omega slider updates rho_ee", await until(page, (b) => document.querySelector('.readout[data-readout="rho_ee"] .ro-value')?.textContent !== b, before));
        await setSlider(page, "omega", 3);
        await page.waitForTimeout(500);
      }
      if (panel === "hom") {
        const before = (await readoutText(page, "indistinguishability")).value;
        await setSlider(page, "gstar", 2);
        check("lindblad hom: gamma* slider updates I", await until(page, (b) => document.querySelector('.readout[data-readout="indistinguishability"] .ro-value')?.textContent !== b, before));
        await setSlider(page, "gstar", 10);
        await page.waitForTimeout(500);
      }
      await chartReady(page, 1);
      await page.waitForTimeout(400);
      await shot(`lindblad-${panel}`);
    }
    // mobile width: no horizontal page scroll
    await page.setViewportSize({ width: 390, height: 900 });
    await page.goto(`${BASE}/#/explain/cw-g2`, { waitUntil: "domcontentloaded" });
    await chartReady(page, 1);
    const sw = await page.evaluate(() => ({ sw: document.documentElement.scrollWidth, cw: document.documentElement.clientWidth }));
    check("cw-g2: no horizontal page scroll at 390 px", sw.sw <= sw.cw + 1, JSON.stringify(sw));
    await page.setViewportSize(VIEW);

    check("no console errors, page errors or failed requests", errors.length === 0, errors.slice(0, 3).join(" | "));
  } catch (e) {
    check("script ran to completion", false, e.stack || e.message);
  } finally {
    await browser.close();
  }
  const n = results.filter((r) => r.ok).length;
  console.log(`\n${n}/${results.length} studio explorer A web checks passed`);
  process.exit(n === results.length ? 0 : 1);
})();
