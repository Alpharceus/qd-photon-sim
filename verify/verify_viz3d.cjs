// Playwright check of the FSIM Studio 3D views (fsim_studio/web/viz3d-lab.html).
//
//   node verify/verify_viz3d.cjs http://127.0.0.1:8766
//
// Starts nothing: the server URL is argv[2].  For each of the 8 lab scenes:
// waits for data-ready, asserts no console/page errors, a non-blank WebGL
// canvas (sampled pixels), the draw-call / triangle budget, and writes a
// 1600x1000 screenshot.  Also checks the SVG fallback mounts from the same
// spec.  studio-07 adds: the rt_edge lattice (every row drawn, non-gating
// rows grey, banner from banner fields only), numbered pins + one control bar
// + zero overlaps at 870x420, 640x360 and 342x250, and the enclosure floor /
// beam.  studio-08 adds the tiny sizes (342x250, 390x280): plate <= 2 rows,
// the bar wraps with no control clipped, and the device box covers >= 35% of
// the canvas.  Playwright is resolved from node_modules (playwright@1.63.0,
// chromium-1243 in %LOCALAPPDATA%\ms-playwright); no package.json in the
// repo root.
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
let chromium;
try {
  ({ chromium } = require(path.join(ROOT, '.workers', 'studio', 'node', 'node_modules', 'playwright')));
} catch (e) {
  ({ chromium } = require('playwright'));
}

const BASE = (process.argv[2] || 'http://127.0.0.1:8766').replace(/\/$/, '');
const OUTDIR = path.join(ROOT, '.workers', 'studio', 'shots', '3d');
const SCENES = ['device-edge', 'device-nanowire', 'device-cavity', 'band-nitride', 'band-inp', 'cascade', 'surface', 'lattice', 'lattice-rt_edge'];
const MAX_CALLS = 120, MAX_TRIS = 300000;

const results = [];
function check(ok, name, detail = '') {
  results.push(!!ok);
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}${detail ? ' -- ' + detail : ''}`);
}

(async () => {
  fs.mkdirSync(OUTDIR, { recursive: true });
  const browser = await chromium.launch({ args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
  const ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 }, deviceScaleFactor: 1 });
  for (const id of SCENES) {
    const page = await ctx.newPage();
    const errors = [];
    page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('pageerror', (e) => errors.push(String(e)));
    try {
      await page.goto(`${BASE}/viz3d-lab.html?scene=${id}&only=1`, { waitUntil: 'load', timeout: 60000 });
      await page.waitForFunction(() => {
        const v = document.getElementById('view');
        return v && (v.dataset.ready === '1' || v.dataset.error);
      }, null, { timeout: 120000 });
      const st = await page.evaluate(() => {
        const v = document.getElementById('view');
        return { mode: v.dataset.mode, error: v.dataset.error || null };
      });
      check(!st.error, `${id}: scene loaded`, st.error || '');
      check(st.mode === 'webgl', `${id}: WebGL mode`, `mode=${st.mode}`);
      await page.waitForTimeout(id === 'cascade' ? 2600 : 900);
      const s = await page.evaluate(() => document.getElementById('view').__viz3d._sample());
      check(s && s.lit > s.samples * 0.04 && s.distinct >= 8, `${id}: canvas not blank`,
        s ? `lit ${s.lit}/${s.samples}, distinct colours ${s.distinct}` : 'no sample');
      const b = await page.evaluate(() => document.getElementById('view').__viz3d._stats());
      check(b && b.calls < MAX_CALLS && b.triangles < MAX_TRIS, `${id}: budget`, b ? `${b.calls} draw calls, ${b.triangles} triangles` : 'no stats');
      await page.screenshot({ path: path.join(OUTDIR, `${id}.png`) });
      // keyboard bookmarks tween (600 ms) without errors
      await page.focus('#view');
      for (const k of ['2', '1']) { await page.keyboard.press(k); await page.waitForTimeout(700); }
      check(fs.existsSync(path.join(OUTDIR, `${id}.png`)), `${id}: screenshot written`);
      check(errors.length === 0, `${id}: no console errors`, errors.slice(0, 3).join(' | '));
    } catch (err) {
      check(false, `${id}: run`, String(err).slice(0, 300));
    }
    await page.close();
  }
  // Ceiling: device scenes stand on the enclosure floor; the edge ridge
  // carries the breadboard beam into its facet with an emission glow.
  {
    const page = await ctx.newPage();
    for (const id of ['device-edge', 'device-staged', 'device-cavity']) {
      await page.goto(`${BASE}/viz3d-lab.html?scene=${id}&only=1`, { waitUntil: 'load' });
      await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 120000 });
      const f = await page.evaluate(() => { const a = document.getElementById('view').__viz3d; return { floor: a._has('enclosure-floor'), beam: a._has('breadboard-beam'), glow: a._has('emission-glow') }; });
      check(f.floor, `${id}: device posted on the enclosure floor (hole grid)`);
      if (id === 'device-edge') check(f.beam && f.glow, `${id}: beam enters the facet with an emission glow`);
    }
    await page.close();
  }

  // I4 / H1: rt_edge lattice, all models.  Every row is drawn; rows of the
  // non-headline models are grey ("never gates"); the banner text and its
  // style come only from banner fields (gating rows).
  {
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', (e) => errors.push(String(e)));
    for (const [id, rows, ng] of [['lattice-rt_edge', 3072, 2304], ['lattice-rt_edge-headline', 768, 0]]) {
      await page.goto(`${BASE}/viz3d-lab.html?scene=${id}&only=1`, { waitUntil: 'load' });
      await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 120000 });
      const r = await page.evaluate(() => {
        const v = document.getElementById('view');
        const b = v.querySelector('.v3-banner');
        return { info: v.__viz3d._info(), cls: b ? b.className : '', text: b ? b.textContent : '' };
      });
      const i = r.info || {};
      check(i.shown === rows && i.rows === rows, `${id}: all ${rows} rows drawn as glyphs`, `shown ${i.shown}/${i.rows}`);
      check(i.nonGating === ng && i.gating === rows - ng, `${id}: ${ng} non-gating rows greyed`, `gating ${i.gating}, non-gating ${i.nonGating}`);
      check(/v3-banner-fail/.test(r.cls) && !/v3-banner-pass/.test(r.cls) && /^0 \/ 768 headline-model rows/.test(i.banner || ''),
        `${id}: banner from banner fields, fail style (0 / 768 gating)`, `${r.cls} | ${i.banner}`);
      check(i.shells === 0, `${id}: no pass shell drawn for non-gating passes`, `shells ${i.shells}`);
      if (ng) check(/never gates/.test(r.text) && /never gates/.test(i.nonGatingText || ''), `${id}: non-gating passes reported separately ("never gates")`, i.nonGatingText || 'none');
    }
    check(errors.length === 0, 'lattice-rt_edge: no page errors', errors.slice(0, 2).join(' | '));
    await page.close();
  }

  // Extra views: the legacy mesa card full size, then the embedded sizes.
  // Compact (< 1000 x 500): callouts are numbered pins keyed in a list, the
  // readout plate shows <= 4 rows + "+N more", the control groups merge into
  // one bar under the canvas, and nothing overlaps (HUD panels and pins);
  // narrow (< 560 wide): the title never sits under the readout plate.
  const EMBED = [['device-staged', 1600, 1000]];
  for (const id of ['device-edge', 'device-staged']) for (const [w, h] of [[870, 420], [640, 360], [342, 250], [390, 280]]) EMBED.push([id, w, h]);
  for (const [id, w, h] of EMBED) {
    const name = w === 1600 ? id : `${id}-${w}x${h}`;
    const page = await ctx.newPage();
    const errors = [];
    page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('pageerror', (e) => errors.push(String(e)));
    try {
      const q = w === 1600 ? `scene=${id}` : `scene=${id}&w=${w}&h=${h}`;
      await page.goto(`${BASE}/viz3d-lab.html?${q}&only=1`, { waitUntil: 'load' });
      await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 120000 });
      await page.waitForTimeout(900);
      const st = await page.evaluate(() => {
        const v = document.getElementById('view');
        const R = v.getBoundingClientRect();
        const vis = (e) => { if (!e) return false; const r = e.getBoundingClientRect(); const cs = getComputedStyle(e); return r.width > 0 && r.height > 0 && cs.visibility !== 'hidden' && cs.display !== 'none'; };
        const named = [['title', v.querySelector('.v3-head')], ['readouts', v.querySelector('.v3-right')], ['notes', v.querySelector('.v3-notes')],
          ['controls', v.querySelector('.v3-controls')], ['banner', v.querySelector('.v3-banner')]].filter(([, e]) => vis(e));
        const chips = [...v.querySelectorAll('.v3-label-chip')].filter((e) => e.style.visibility !== 'hidden' && e.style.display !== 'none');
        const boxes = [...named.map(([n, e]) => [n, e.getBoundingClientRect()]), ...chips.map((e, i) => [`pin${i}`, e.firstChild.getBoundingClientRect()])];
        const over = [];
        for (let i = 0; i < boxes.length; i++) for (let j = i + 1; j < boxes.length; j++) {
          const [na, a] = boxes[i], [nb, b] = boxes[j];
          if (a.left < b.right - 0.5 && a.right > b.left + 0.5 && a.top < b.bottom - 0.5 && a.bottom > b.top + 0.5) over.push(`${na}/${nb}`);
        }
        const outside = boxes.filter(([, r]) => r.left < R.left - 0.5 || r.right > R.right + 0.5 || r.top < R.top - 0.5 || r.bottom > R.bottom + 0.5).map(([n]) => n);
        const compact = v.classList.contains('v3-compact');
        const fullText = chips.filter((e) => getComputedStyle(e.querySelector('.v3-label-text')).display !== 'none').length;
        const pinsShown = chips.filter((e) => getComputedStyle(e.querySelector('.v3-label-pin')).display !== 'none').length;
        const rows = [...v.querySelectorAll('.v3-readout')].filter((r) => getComputedStyle(r).display !== 'none').length;
        const more = v.querySelector('.v3-more-btn');
        const ctl = v.querySelector('.v3-controls').getBoundingClientRect(), cw = v.querySelector('.v3-canvas-wrap').getBoundingClientRect();
        const keyRows = v.querySelectorAll('.v3-key-row').length;
        const ctls = [v.querySelector('.v3-controls'), ...v.querySelectorAll('.v3-controls .v3-ctl-group, .v3-controls .v3-btn, .v3-controls .v3-seg-btn, .v3-controls .v3-slider')];
        const clipped = ctls.filter((e) => vis(e) && e.scrollWidth > e.clientWidth + 0.5)
          .map((e) => `${e.className}:${e.textContent.trim().slice(0, 24)} ${e.scrollWidth}>${e.clientWidth}`);
        const cvr = cw;
        const ctlOut = [...v.querySelectorAll('.v3-controls .v3-btn, .v3-controls .v3-seg-btn, .v3-controls .v3-slider')].filter(vis)
          .filter((e) => { const r = e.getBoundingClientRect(); return r.right > ctl.right + 0.5 || r.left < ctl.left - 0.5 || r.bottom > ctl.bottom + 0.5 || r.top < ctl.top - 0.5; }).length;
        const ctlRows = new Set([...v.querySelectorAll('.v3-controls .v3-ctl-group')].filter(vis).map((e) => Math.round(e.getBoundingClientRect().top))).size;
        return { tiny: v.classList.contains('v3-tiny'), clipped, ctlOut, ctlRows, cover: v.__viz3d._deviceCover(), canvasH: Math.round(cvr.height), compact, narrow: v.classList.contains('v3-narrow'), chips: chips.length, fullText, pinsShown, over, outside, rows,
          moreShown: !!more && getComputedStyle(more).display !== 'none', nReadouts: v.querySelectorAll('.v3-readout').length,
          barUnder: Math.abs(ctl.top - cw.bottom) < 1.5 && Math.abs(ctl.bottom - R.bottom) < 2 && ctl.width > R.width - 4, keyRows };
      });
      await page.screenshot({ path: path.join(OUTDIR, `${name}.png`), clip: w === 1600 ? undefined : { x: 0, y: 0, width: w, height: h } });
      check(st.over.length === 0 && st.outside.length === 0 && errors.length === 0, `${name}: renders, 0 overlaps (HUD panels + labels)`,
        `chips ${st.chips}, overlaps [${st.over.join(', ')}], outside [${st.outside.join(', ')}]${errors.length ? ', ' + errors[0] : ''}`);
      if (w < 1000) {
        check(st.compact && st.fullText === 0 && st.pinsShown === st.chips && st.keyRows > 0, `${name}: compact, callouts are numbered pins + key list`,
          `pins ${st.pinsShown}/${st.chips}, full callouts ${st.fullText}, key rows ${st.keyRows}`);
        check(st.rows <= 4 && (st.nReadouts <= 4 || st.moreShown), `${name}: readout plate <= 4 rows + "+N more"`, `${st.rows} of ${st.nReadouts} shown, more ${st.moreShown}`);
        check(st.barUnder, `${name}: VIEWS / STACK / CUT-AWAY in one bar under the canvas`);
        if (w < 560) check(st.narrow && !st.over.includes('title/readouts'), `${name}: narrow, title not under the readout plate`);
        if (w < 480 || h < 300) {
          check(st.tiny && st.rows <= 2 && (st.nReadouts <= 2 || st.moreShown), `${name}: tiny, readout plate <= 2 rows + "+N more"`, `${st.rows} of ${st.nReadouts} shown, more ${st.moreShown}`);
          check(st.clipped.length === 0 && st.ctlOut === 0 && st.ctlRows <= 2, `${name}: control bar wraps (<= 2 rows), no control clipped`,
            `rows ${st.ctlRows}, clipped [${st.clipped.join('; ')}], outside the bar ${st.ctlOut}`);
          check(st.cover !== null && st.cover >= 0.35, `${name}: device box covers >= 35% of the canvas`, `cover ${st.cover === null ? 'n/a' : (st.cover * 100).toFixed(1) + '%'}, canvas h ${st.canvasH}`);
        }
      } else {
        check(!st.compact && st.fullText === st.chips, `${name}: full-size callouts`, `chips ${st.chips}`);
      }
    } catch (err) {
      check(false, `${name}: run`, String(err).slice(0, 300));
    }
    await page.close();
  }

  // studio-p2d: the cheap playhead hook.  setGlow / setPlayhead scale the emission glow of a built device
  // scene without rebuilding it: no new WebGL context or canvas, one render, the pixels get dimmer at
  // glow 0 and return exactly at glow 1; kinds without a device scene report false.
  for (const id of ['device-edge', 'device-nanowire', 'device-cavity']) {
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', (e) => errors.push(String(e)));
    try {
      await page.goto(`${BASE}/viz3d-lab.html?scene=${id}&only=1`, { waitUntil: 'load' });
      await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 120000 });
      await page.waitForTimeout(700);
      const r = await page.evaluate(async () => {
        const v = document.getElementById('view');
        const api = v.__viz3d;
        const canvases = () => document.querySelectorAll('canvas').length;
        let ctxCalls = 0;
        const orig = HTMLCanvasElement.prototype.getContext;
        HTMLCanvasElement.prototype.getContext = function (t, ...a) { if (/webgl/.test(String(t))) ctxCalls++; return orig.call(this, t, ...a); };
        const nC = canvases();
        const frame = () => new Promise((res) => requestAnimationFrame(() => requestAnimationFrame(res)));
        const lit = () => { const s = api._sample(); return s ? s.lit : null; };
        await frame(); const base = lit();
        const okHi = api.setPlayhead({ T: 250, glow: 1 });
        await frame(); const hi = lit();
        const okLo = api.setGlow(0);
        await frame(); const lo = lit();
        const g0 = api._glow();
        const tLo = api._glowTargets();
        api.setPlayhead({ glow: 1 });
        await frame(); const back = lit();
        const tBack = api._glowTargets();
        const out = { okHi, okLo, g0, tLo, tBack, base, hi, lo, back, gBack: api._glow(), ctxCalls, nC, nC2: canvases(), T: v.dataset.playheadT, mode: api.mode };
        HTMLCanvasElement.prototype.getContext = orig;
        return out;
      });
      check(r.mode === 'webgl' && r.okHi === true && r.okLo === true && r.g0 === 0 && r.gBack === 1,
        `${id}: setGlow/setPlayhead accepted on a device scene`, JSON.stringify(r));
      check(r.ctxCalls === 0 && r.nC === r.nC2, `${id}: the playhead hook creates no WebGL context or canvas`, `contexts ${r.ctxCalls}, canvases ${r.nC}->${r.nC2}`);
      check(r.tLo.n > 0 && Math.abs(r.tLo.ratio - 0.2) < 1e-6 && Math.abs(r.tBack.ratio - 1) < 1e-9 && r.back === r.hi && r.lo <= r.hi && r.T === '250', `${id}: glow 0 scales ${r.tLo.n} emission materials to 20%, glow 1 restores them exactly; T recorded without a rebuild`, `ratio ${r.tLo.ratio} -> ${r.tBack.ratio}; lit ${r.hi} -> ${r.lo} -> ${r.back}`);
      check(errors.length === 0, `${id}: playhead hook, no page errors`, errors.slice(0, 2).join(' | '));
    } catch (err) {
      check(false, `${id}: playhead hook run`, String(err).slice(0, 300));
    }
    await page.close();
  }
  {
    const page = await ctx.newPage();
    await page.goto(`${BASE}/viz3d-lab.html?scene=cascade&only=1`, { waitUntil: 'load' });
    await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 120000 });
    const r = await page.evaluate(() => { const a = document.getElementById('view').__viz3d; return { g: a.setGlow(0.3), p: a.setPlayhead({ T: 1, glow: 0.3 }) }; });
    check(r.g === false && r.p === false, 'cascade: the playhead hook refuses a non-device scene (returns false, scene untouched)', JSON.stringify(r));
    await page.close();
  }

  // SVG fallback from the same SceneSpec.
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(String(e)));
  await page.goto(`${BASE}/viz3d-lab.html?scene=device-edge&only=1&fallback=1`, { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('view').dataset.ready === '1', null, { timeout: 60000 });
  const fb = await page.evaluate(() => ({ mode: document.getElementById('view').dataset.mode, rects: document.querySelectorAll('#view svg rect').length }));
  check(fb.mode === 'fallback' && fb.rects >= 5 && errors.length === 0, 'fallback: SVG cross-section from the same spec', `mode=${fb.mode}, rects=${fb.rects}`);
  await page.screenshot({ path: path.join(OUTDIR, 'fallback-device-edge.png') });
  await browser.close();
  const n = results.filter(Boolean).length;
  console.log(`${n}/${results.length} viz3d checks passed`);
  process.exit(n === results.length ? 0 : 1);
})().catch((e) => { console.error(e); process.exit(1); });
