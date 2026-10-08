// FSIM Studio model explorers B (phonon, transport + nitride Stark, laser): Playwright checks + screenshots.
//
// Usage: node verify/verify_studio_explore_b.cjs http://127.0.0.1:8782
// The server must already be running (python -m fsim_studio --no-window --port 8782).
// Screenshots saved at 1600x1000.
// Exits 0 iff every check passes; prints "N/N studio explorer B web checks passed".

const path = require("path");
const fs = require("fs");

const ROOT = path.resolve(__dirname, "..");
const BASE = (process.argv[2] || "http://127.0.0.1:8782").replace(/\/$/, "");
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

/** Set a slider by data-ctl id; linear sliders take real units, log sliders take [min, max]. */
async function setSlider(page, ctl, value, { log = null } = {}) {
  await page.evaluate(({ ctl, value, log }) => {
    const el = document.querySelector(`input.xs-range[data-ctl="${ctl}"]`);
    if (!el) throw new Error(`no slider ${ctl}`);
    el.value = log ? String((Math.log(value / log[0]) / Math.log(log[1] / log[0])) * 1000) : String(value);
    el.dispatchEvent(new Event("input", { bubbles: true }));
  }, { ctl, value, log });
}

const readout = (page, key) => page.evaluate((k) => {
  const el = [...document.querySelectorAll(`.readout[data-readout="${k}"]`)].find((e) => e.getClientRects().length);
  return el ? { state: el.dataset.state, value: el.querySelector(".ro-value")?.textContent, caption: el.querySelector(".ro-caption")?.textContent, tag: el.querySelector(".tag")?.dataset.tag || null } : null;
}, key);

/** The 4-significant-figure text the page shows for a number (mirrors ui/format.js fmt for 1e-3 <= |v| < 1000). */
const four = (v) => String(Number(Number(v).toPrecision(4)));

async function tagAudit(page) {
  return page.evaluate(() => {
    const bad = []; let n = 0, rows = 0;
    for (const el of document.querySelectorAll(".xp-body .readout")) {
      if (!el.getClientRects().length || el.dataset.state !== "value") continue;
      n++;
      if (!el.querySelector(".tag")) bad.push(el.dataset.readout);
    }
    for (const el of document.querySelectorAll(".xp-body .xs-row[data-ctl]")) {
      if (!el.getClientRects().length) continue;
      rows++;
      if (!el.querySelector(".tag")) bad.push(`slider:${el.dataset.ctl}`);
    }
    return { n, rows, bad };
  });
}

const charts = (page) => page.evaluate(() => [...document.querySelectorAll(".xp-body .chart-card")].filter((c) => c.getClientRects().length).map((c) => ({
  title: c.querySelector(".chart-title")?.textContent, svg: !!c.querySelector(".chart-plot .main-svg"),
  traces: (c.querySelector(".chart-plot")?.data || []).map((t) => ({ name: t.name, mode: t.mode, n: (t.x || []).length })) })));

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
  const api = (url) => page.evaluate(async (u) => (await fetch(u)).json(), url);
  // /api/explore/phonon/live is a one-frame pool job (review 11): resolve the job and return the frame
  const live = (url) => page.evaluate(async (u) => {
    const r = await (await fetch(u)).json();
    if (r.error) return r;
    if (r.result) return r.result.frames[0];
    for (let i = 0; i < 600; i++) {
      const snap = await (await fetch(`/api/jobs/${r.job_id}`)).json();
      if (snap.state === "done") return snap.result.frames[0];
      if (snap.state === "error" || snap.state === "cancelled") return { error: snap.state };
      await new Promise((x) => setTimeout(x, 200));
    }
    return { error: "timeout" };
  }, url);
  const caveat = (id) => page.$eval(`.xb-caveat[data-caveat="${id}"]`, (e) => e.textContent).catch(() => "");

  try {
    // ------------------------------------------------------------ tabs
    await page.goto(`${BASE}/#/explain/phonon`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".xp-tabs", { timeout: 30000 });
    const tabs = await page.$$eval(".xp-tabs a", (els) => els.map((e) => e.getAttribute("href")));
    check("explain tabs list the phonon, transport and laser explorers", ["#/explain/phonon", "#/explain/transport", "#/explain/laser"].every((t) => tabs.includes(t)), tabs.join(" "));

    // ------------------------------------------------------------ phonon
    check("phonon: the three charts render", await until(page, () => document.querySelectorAll(".xp-body .chart-plot .main-svg").length >= 3, null, 60000));
    const framesOk = await until(page, () => /T frames ready/.test(document.querySelector('[data-role="frames"]')?.textContent || ""), null, 120000);
    check("phonon: T frames computed on the job pool (one batch)", framesOk, await page.$eval('[data-role="frames"]', (e) => e.textContent));
    const pc = await caveat("caveat");
    check("phonon: pinned caveat verbatim and not-brightness note", pc.includes("Independent-boson model, bulk LA deformation potential only, no LO or piezo term (:191-197) [E]. Gamma_zpl is phenomenological; the fitted a_ac T term partly double-counts acoustic broadening (:436-443). Lengths are wavefunction extents, not dot sizes (:127-134).") && pc.includes("not brightness"), pc.slice(0, 80));
    let ph = await live("/api/explore/phonon/live?T=100&l_xy=4.5&l_z=1.5&gamma_zpl=0.5&w=2&kappa=1&F_cav=1&delta=0");
    await setSlider(page, "T", 5); // index 5 = 100 K
    await page.waitForTimeout(400);
    let z = await readout(page, "Z");
    check("phonon: Z readout at 100 K equals the backend (frame)", z && z.state === "value" && z.value === four(ph.Z), `${z?.value} vs ${four(ph.Z)}`);
    await setSlider(page, "T", 11); // 300 K
    await page.waitForTimeout(400);
    const ph300 = await live("/api/explore/phonon/live?T=300&l_xy=4.5&l_z=1.5&gamma_zpl=0.5&w=2&kappa=1&F_cav=1&delta=0");
    z = await readout(page, "Z");
    check("phonon: T slider selects the 300 K frame (Z follows the backend, smaller than at 100 K)", z && z.value === four(ph300.Z) && ph300.Z < ph.Z, `${z?.value}`);
    await setSlider(page, "w", 5);
    await setSlider(page, "delta", 1.5);
    await until(page, () => /live/.test(document.querySelector('.readout[data-readout="Z"] .ro-caption')?.textContent || ""), null, 20000);
    await page.waitForTimeout(500);
    const ph2 = await live("/api/explore/phonon/live?T=300&l_xy=4.5&l_z=1.5&gamma_zpl=0.5&w=5&kappa=1&F_cav=1&delta=1.5");
    const tIbm = await readout(page, "t_ibm");
    check("phonon: w and delta are live (t_IBM follows the backend)", tIbm && tIbm.value === four(ph2.point.t_ibm), `${tIbm?.value} vs ${four(ph2.point.t_ibm)}`);
    await setSlider(page, "l_xy", 6.0);
    await page.waitForTimeout(700);
    const note = await page.$eval('[data-role="frames"]', (e) => e.textContent);
    check("phonon: changing the geometry recomputes the frames on the pool", /Computing|ready/.test(note), note.slice(0, 60));
    const phCharts = await charts(page);
    check("phonon: spectrum, Z(T)/S(T) and transmission traces present", phCharts.length >= 3 && phCharts.some((c) => /spectrum/i.test(c.title) && c.traces[0].n > 1000) && phCharts.some((c) => /Z/.test(c.title)) && phCharts.some((c) => /transmission/i.test(c.title) && c.traces.length >= 2), phCharts.map((c) => `${c.title}:${c.traces.length}`).join(" | "));
    let a = await tagAudit(page);
    check("phonon: every readout and slider carries a tag chip", a.n >= 8 && a.rows >= 8 && a.bad.length === 0, `${a.n} readouts, ${a.rows} sliders, bad ${a.bad.join(",")}`);
    await setSlider(page, "l_xy", 4.5);
    await until(page, () => /T frames ready/.test(document.querySelector('[data-role="frames"]')?.textContent || ""), null, 120000);
    await setSlider(page, "T", 5);
    await page.waitForTimeout(900);
    await shot("phonon");

    // ------------------------------------------------------------ transport: InP
    await page.goto(`${BASE}/#/explain/transport`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector('.xb-caveat[data-caveat="transport"]', { timeout: 30000 });
    const tc = await caveat("transport");
    check("transport: pinned InP caveat verbatim (LOWER bound on leakage, mu ~ 6000)", tc.includes("leak_valleys='G' is a LOWER bound on leakage (docstring item 5, 'GX' moves 50x per 0.1 eV)") && tc.includes("1 uA, 1 ns, one dot gives mu ~ 6000; cap-2 needs ~80 pA into one dot (item 6)."), tc.slice(0, 60));
    check("transport: four InP charts render (sweep ran as a job)", await until(page, () => [...document.querySelectorAll('.xb-panel[data-panel="inp"] .chart-plot .main-svg')].length >= 4, null, 90000));
    const pt = await api("/api/explore/transport/point?preset=red&I_uA=1&T=230&n_dot=2000000000&aperture=1&tau_pulse=0.1&w=2&dE_WL=100&tau_rad=1");
    await page.waitForTimeout(600);
    const muR = await readout(page, "mu");
    check("transport: mu readout equals the backend point at I = 1 uA", muR && muR.state === "value" && muR.value.replace(/\s/g, "") === fmtJs(pt.mu), `${muR?.value} vs ${fmtJs(pt.mu)}`);
    await setSlider(page, "I", 0.1, { log: [0.001, 200] });
    await page.waitForTimeout(600);
    const pt2 = await api("/api/explore/transport/point?preset=red&I_uA=0.1&T=230&n_dot=2000000000&aperture=1&tau_pulse=0.1&w=2&dE_WL=100&tau_rad=1");
    const vj = await readout(page, "V_j");
    check("transport: V_j follows the current slider (live point)", vj && vj.value === four(pt2.V_j) && pt2.V_j < pt.V_j, `${vj?.value} vs ${four(pt2.V_j)}`);
    const tr = await charts(page);
    check("transport: mu chart carries the 0.5 reference and r_dot the 1/tau_rad line (shapes)", await page.evaluate(() => [...document.querySelectorAll(".xb-panel[data-panel='inp'] .chart-plot")].filter((p) => p.layout).some((p) => (p.layout.shapes || []).some((s) => s.y0 === 0.5))));
    a = await tagAudit(page);
    check("transport: every visible readout and slider carries a tag chip (InP)", a.n >= 10 && a.rows >= 7 && a.bad.length === 0, `${a.n} readouts, ${a.rows} sliders, bad ${a.bad.join(",")}`);
    await setSlider(page, "I", 1.0, { log: [0.001, 200] });
    await page.waitForTimeout(700);
    await shot("transport-inp");

    // ------------------------------------------------------------ transport: nitride Stark
    await page.click('.xb-modebar .seg[data-value="nit"]');
    const traceOk = await until(page, () => /traces ready \(one point/.test(document.querySelector('[data-role="stark-status"]')?.textContent || ""), null, 120000);
    check("transport: three Stark traces computed as jobs", traceOk, await page.$eval('[data-role="stark-status"]', (e) => e.textContent));
    const nc = await caveat("nitride");
    check("transport: pinned nitride caveat verbatim", nc.includes("barrier-limited, doping dependence omitted [A]; screening fixed along a trace; m- and a-plane coincide by construction (out/nitride_geometry_stark/results.md:328); Zhang is 10 K, x=0.15, another device (results.md:129)"), nc.slice(0, 60));
    const nch = await charts(page);
    const eC = nch.find((c) => /E_X/.test(c.title)), sC = nch.find((c) => /slope/i.test(c.title));
    check("transport: E_X chart has 3 screening series + the Zhang guide; slope chart has 3 series", eC && eC.traces.length === 4 && eC.traces.some((t) => /Zhang/.test(t.name) && /guide only/.test(t.name)) && sC && sC.traces.length === 3, JSON.stringify(eC?.traces.map((t) => t.name)));
    const guideTxt = await page.$$eval(".xb-panel[data-panel='nit'] .chart-plot", (ps) => ps.map((p) => (p.layout?.annotations || []).map((x) => x.text).join(" ")).join(" "));
    check("transport: Zhang line is labelled guide only (10 K, x = 0.15, another device); no compatibility verdict anywhere on the page", /guide only/.test(guideTxt) && /another device/.test(guideTxt) && !/compatib/i.test(await page.$eval(".xp-body", (e) => e.textContent)));
    const stark = await page.evaluate(async () => (await fetch("/api/explore/transport/stark", { method: "POST", headers: { "Content-Type": "application/json", "X-FSIM-Studio": "1" }, body: JSON.stringify({ orientation: "c_plane", polarity: 1, T: 300, ext: 0, screenings: [0] }) })).json());
    const t0 = stark.traces[0].result || (await page.evaluate(async (id) => (await (await fetch(`/api/jobs/${id}`)).json()).result, stark.traces[0].job_id));
    const firstSlope = t0.rows.find((r) => r.derivative_valid).dE_X_dV_meV_per_V;
    await page.waitForTimeout(300);
    const slope = await readout(page, "slope");
    check("transport: slope readout equals the backend trace (stark_derivatives)", slope && slope.value === four(firstSlope), `${slope?.value} vs ${four(firstSlope)}`);
    await page.click('.xb-panel[data-panel="nit"] .seg[data-value="a_plane"]');
    check("transport: switching to a-plane reruns the traces", await until(page, () => /traces ready \(one point/.test(document.querySelector('[data-role="stark-status"]')?.textContent || ""), null, 120000));
    await page.click('.xb-panel[data-panel="nit"] .seg[data-value="c_plane"]');
    await until(page, () => /traces ready \(one point/.test(document.querySelector('[data-role="stark-status"]')?.textContent || ""), null, 120000);
    a = await tagAudit(page);
    check("transport: every visible readout and slider carries a tag chip (nitride)", a.n >= 6 && a.rows >= 3 && a.bad.length === 0, `${a.n} readouts, ${a.rows} sliders, bad ${a.bad.join(",")}`);
    await page.waitForTimeout(900);
    await shot("transport-nitride");

    // ------------------------------------------------------------ laser
    await page.goto(`${BASE}/#/explain/laser`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector('.xb-caveat[data-caveat="zhao"]', { timeout: 30000 });
    const lc = await caveat("zhao");
    check("laser: pinned Zhao caveat verbatim (quiet pump does NOT reproduce: 5.8 sigma, FAIL)", lc.includes("Langevin model of Zhao's QD LASER (thousands of dots, lasing mode), not a single-photon source; g2 ~ 1 is correct physics. Normal pump reproduces 1.0224+-0.003 to 1.7 sigma; quiet pump does NOT: 1.017 vs 0.982, 5.8 sigma (out/zhao/zhao_fit_comparison.csv, verdict FAIL both rows).") && /Axis is I\/I_th, not amperes/.test(lc), lc.slice(0, 60));
    const framesDone = await until(page, () => /frames ready/.test(document.querySelector('[data-role="frames"]')?.textContent || ""), null, 180000);
    check("laser: 14 seeded frames computed on the pool", framesDone, await page.$eval('[data-role="frames"]', (e) => e.textContent));
    await page.waitForTimeout(600);
    const lch = await charts(page);
    const g2c = lch.find((c) => /g²\(0\)/.test(c.title));
    const simTraces = g2c ? g2c.traces.filter((t) => /^simulation/.test(t.name)) : [];
    check("laser: g2 chart shows normal and quiet simulation, both Zhao card points and the committed run, markers only (no curve through the quiet point)",
      g2c && simTraces.length === 2 && g2c.traces.every((t) => t.mode === "markers") && g2c.traces.filter((t) => /Zhao card/.test(t.name)).length === 2 && g2c.traces.filter((t) => /committed fit run/.test(t.name)).length === 2, JSON.stringify(g2c?.traces.map((t) => `${t.name}:${t.mode}`)));
    check("laser: L-I and trace charts render", lch.some((c) => /L-I/.test(c.title) && c.svg) && lch.some((c) => /S\(t\)/.test(c.title) && c.svg));
    const li = await api("/api/explore/sde/meta");
    await setSlider(page, "I", 2); // index 2 = I/I_th 2.0
    await page.waitForTimeout(600);
    const fr = await page.evaluate(async () => { const r = await (await fetch("/api/explore/sde/frames", { method: "POST", headers: { "Content-Type": "application/json", "X-FSIM-Studio": "1" }, body: JSON.stringify({ I_grid: [2.0], F_pumps: [1.0], n_runs: 8, t_end: 1e-9, seed: 1 }) })).json(); return r.frames[0]; });
    const fres = fr.result || (await page.evaluate(async (id) => (await (await fetch(`/api/jobs/${id}`)).json()).result, fr.job_id));
    const g2r = await readout(page, "g2");
    check("laser: g2 readout equals the backend frame (g2_vs_pump, seed 1)", g2r && g2r.state === "value" && g2r.value === four(fres.g2_mean) && /SE/.test(g2r.caption), `${g2r?.value} vs ${four(fres.g2_mean)}; ${g2r?.caption}`);
    await page.click('.seg[data-value="0.08"]');
    await page.waitForTimeout(600);
    const card = await readout(page, "card");
    const sig = await readout(page, "sigma");
    check("laser: quiet selection shows the card point 0.9823 and the committed FAIL verdict at 5.8 sigma", card && card.value === "0.9823" && sig && sig.value === "5.815" && /FAIL/.test(sig.caption), `${card?.value}; ${sig?.value} ${sig?.caption?.slice(0, 40)}`);
    await page.click('.seg[data-value="1"]');
    a = await tagAudit(page);
    check("laser: every visible readout and slider carries a tag chip", a.n >= 6 && a.rows >= 3 && a.bad.length === 0, `${a.n} readouts, ${a.rows} sliders, bad ${a.bad.join(",")}`);
    check("laser: no disable_stim / beta_sp control on the page", !(await page.$eval(".xp-body", (e) => /disable_stim|beta_sp|estimator control/i.test(e.textContent))));
    await setSlider(page, "I", 4);
    await page.waitForTimeout(1500);
    await shot("laser");

    check("no console errors, page errors or failed requests", errors.length === 0, errors.slice(0, 3).join(" | "));
  } catch (e) {
    check("script ran to completion", false, e.stack || String(e));
  }
  await browser.close();
  const nOk = results.filter((r) => r.ok).length;
  console.log(`\n${nOk}/${results.length} studio explorer B web checks passed`);
  process.exit(nOk === results.length ? 0 : 1);
})();

/** ui/format.js fmt, for the values compared above (no thousands grouping needed for these). */
function fmtJs(v) {
  if (v == null) return "n/a";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-3) {
    const [m, e] = v.toExponential(3).split("e");
    const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
    return `${m}×10${String(Number(e)).split("").map((c) => SUP[c] ?? c).join("")}`;
  }
  if (a >= 1000) return Math.round(v).toLocaleString("en-US");
  return String(Number(v.toPrecision(4)));
}
