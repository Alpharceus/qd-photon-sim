// Playwright check of the composed-device 3D builders (micropillar, openCavity, bullseye, phcWaveguide,
// phcCavity, mesa, photonicNanowire, planarDbr (+ nitride), ridgeWaveguide, lumpedCavity) in the FSIM Studio 3D lab page.
//
//   node verify/verify_viz3d_compose.cjs http://127.0.0.1:8774
//
// Starts nothing: the server URL is argv[2].  Descriptors come from fsim_studio.scene_compose.describe (run
// through python, no route needed); each is mounted into the lab page's #view with mountScene('device', ...).
// Per builder: data-ready, WebGL mode, no console / page errors, non-blank canvas (pixel variance), dark
// enclosure background (canvas corner pixels, in the light theme as well), the four camera bookmarks tween
// without errors, a 1600x1000 screenshot to .workers/devlib/shots/viz3d-<builder>.png.  Also: an unsupported
// structure kind mounts the message panel (no 3D view), not an error.
// R2-3D pass: every structure kind that has a builder, plus paper cards per kind (descriptor via api_compose
// central_design + display_hint), mounted in WebGL (dark and light) and in forced fallback: data-ready, data-mode, no
// console errors, non-blank canvas / SVG (pixel variance), screenshots .workers/devlib/shots/r2-3d-<name>[-light|-fallback].png.
const path = require('path');
const fs = require('fs');
const { execFileSync } = require('child_process');

const ROOT = path.resolve(__dirname, '..');
let chromium;
try {
  ({ chromium } = require(path.join(ROOT, '.workers', 'studio', 'node', 'node_modules', 'playwright')));
} catch (e) {
  ({ chromium } = require('playwright'));
}

const BASE = (process.argv[2] || 'http://127.0.0.1:8774').replace(/\/$/, '');
const OUTDIR = path.join(ROOT, '.workers', 'devlib', 'shots');

const results = [];
function check(ok, name, detail = '') {
  results.push(!!ok);
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}${detail ? ' -- ' + detail : ''}`);
}

// name -> [builder id expected, python expression for (structure kind, params, material, display)]
const CASES = [
  ['micropillar', 'micropillar', ['micropillar', {}, 'ingaas_gaas', {}]],
  ['openCavity', 'open_cavity', ['open_cavity', {}, 'ingaas_gaas', {}]],
  ['bullseye', 'bullseye', ['cbg', {}, 'ingaas_gaas', {}]],
  ['phcWaveguide', 'phc_waveguide', ['phc_waveguide', {}, 'ingaas_gaas', {}]],
  ['phcCavity', 'phc_cavity', ['phc_cavity', {}, 'ingaas_gaas', {}]],
  ['phcCavity-H1', 'phc_cavity', ['phc_cavity', {}, 'ingaas_gaas', { cell: 'H1' }]],
  ['mesa-sil', 'mesa', ['sil_mesa', {}, 'ingaas_gaas', {}]],
  ['mesa-reflector', 'mesa', ['mesa_reflector', {}, 'inas_inp_cband', {}]],
  ['mesa-slab', 'mesa', ['bulk_planar', { sil: 'weierstrass', mirror: true }, 'ingan_gan', {}]],
];

function describeAll() {
  const py = `
import json, sys
sys.path.insert(0, ${JSON.stringify(ROOT)})
from fsim_studio.scene_compose import describe
cases = json.loads(sys.stdin.read())
out = []
for kind, params, mat, disp in cases:
    spec = {"schema": "fsim-composed/1", "T_K": 4.2, "material": {"kind": mat, "params": {}},
            "structure": {"kind": kind, "params": params}, "drive": {"kind": "x", "params": {}}, "setup": {"kind": "y", "params": {}}}
    out.append(describe(spec, display=disp))
print(json.dumps(out))
`;
  const input = JSON.stringify([...CASES.map((c) => c[2]), ['edge_ridge', {}, 'cdse_znse', {}]]);
  const txt = execFileSync('python', ['-I', '-c', py], { input, cwd: ROOT, stdio: ['pipe', 'pipe', 'ignore'], maxBuffer: 1 << 26 }).toString();
  return JSON.parse(txt.trim().split('\n').pop());
}

// ---- R2-3D: new structure kinds and paper cards
const R2_KINDS = [
  ['photonic_nanowire', 'photonic_nanowire', ['photonic_nanowire', {}, 'ingaas_gaas', {}]],
  ['photonic_trumpet', 'photonic_nanowire', ['photonic_nanowire', { taper_mfr_nm: 900.0, mirror_phase_ok: true }, 'ingaas_gaas', {}]],
  ['planar_dbr', 'planar_dbr', ['planar_dbr', {}, 'ingaas_gaas', {}]],
  ['nitride_cavity', 'nitride_cavity', ['nitride_planar_cavity', {}, 'ingan_gan', {}]],
  ['edge_ridge', 'ridge_waveguide', ['edge_ridge', {}, 'inp_gainp_red', {}]],
  ['generic_cavity', 'lumped_cavity', ['generic_cavity', {}, 'ingaas_gaas', {}]],
  ['generic_cavity-as-micropillar', 'micropillar', ['generic_cavity', {}, 'ingaas_gaas', { display_as: 'micropillar' }]],
];
// [name, builder expected, source, id]: paper cards (cards/papers) and composed examples (cards/composed)
const R2_CARDS = [
  ['card-lazic2017', 'photonic_nanowire', 'paper', 'lazic2017'],
  ['card-munsch2013_trumpet', 'photonic_nanowire', 'paper', 'munsch2013_trumpet'],
  ['card-he2013', 'planar_dbr', 'paper', 'he2013'],
  ['card-somaschi2016', 'micropillar', 'paper', 'somaschi2016'],
  ['card-thomas2021_res', 'micropillar', 'paper', 'thomas2021_res'],
  ['card-example-nitride-rt', 'photonic_nanowire', 'composed', 'example-nitride-rt'],
];

function describeR2() {
  const py = `
import json, sys
sys.path.insert(0, ${JSON.stringify(ROOT)})
import yaml
from pathlib import Path
from fsim_studio.scene_compose import describe
from fsim_studio import api_compose as A
inp = json.loads(sys.stdin.read())
out = []
for kind, params, mat, disp in inp["kinds"]:
    spec = {"schema": "fsim-composed/1", "T_K": 4.2, "material": {"kind": mat, "params": {}},
            "structure": {"kind": kind, "params": params}, "drive": {"kind": "x", "params": {}}, "setup": {"kind": "y", "params": {}}}
    out.append(describe(spec, display=disp))
for src, cid in inp["cards"]:
    if src == "paper":
        doc = yaml.safe_load(Path("cards/papers/%s.yaml" % cid).read_text(encoding="utf-8"))
        spec = A.central_design(cid, doc["device"])["design"]["composed"]
    else:
        spec = yaml.safe_load(Path("cards/composed/%s.yaml" % cid).read_text(encoding="utf-8"))["design"]["composed"]
    spec = dict(spec, schema="fsim-composed/1")
    out.append(describe(spec, display=A.display_hint(spec)))
print(json.dumps(out))
`;
  const input = JSON.stringify({ kinds: R2_KINDS.map((c) => c[2]), cards: R2_CARDS.map((c) => [c[2], c[3]]) });
  const txt = execFileSync('python', ['-I', '-c', py], { input, cwd: ROOT, stdio: ['pipe', 'pipe', 'ignore'], maxBuffer: 1 << 26 }).toString();
  return JSON.parse(txt.trim().split('\n').pop());
}

(async () => {
  fs.mkdirSync(OUTDIR, { recursive: true });
  const descs = describeAll();
  const browser = await chromium.launch({ args: ['--use-angle=d3d11', '--enable-gpu', '--ignore-gpu-blocklist'] });
  const ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 }, deviceScaleFactor: 1 });

  // pixel luminance of the four canvas corners from a screenshot buffer, decoded in the page
  const cornerLum = async (page, buf) => page.evaluate(async (b64) => {
    const img = new Image();
    img.src = 'data:image/png;base64,' + b64;
    await img.decode();
    const c = document.createElement('canvas');
    c.width = img.width; c.height = img.height;
    const g = c.getContext('2d');
    g.drawImage(img, 0, 0);
    const pts = [[6, 60], [img.width - 6, img.height - 80], [6, img.height - 80], [img.width - 6, 60]];
    return pts.map(([x, y]) => { const d = g.getImageData(x, y, 1, 1).data; return (d[0] + d[1] + d[2]) / 3; });
  }, buf.toString('base64'));

  const mount = (page, desc, theme) => page.evaluate(async ({ desc, theme }) => {
    const { mountScene } = await import('/js/viz3d/index.js');
    // fresh container (disposing the lab's view in place would let its late context-loss event repaint the 2D fallback)
    const old = document.getElementById('view');
    const view = document.createElement('div');
    view.id = 'view'; view.dataset.ready = '0';
    view.style.cssText = old.style.cssText;
    old.replaceWith(view);
    try { window.__lab?.api?.dispose(); } catch (e) { /* old view is detached */ }
    try {
      window.__api = mountScene(view, 'device', desc, { theme });
    } catch (e) { view.dataset.error = String(e); }
  }, { desc, theme });

  const wait = (page) => page.waitForFunction(() => {
    const v = document.getElementById('view');
    return v && (v.dataset.ready === '1' || v.dataset.error);
  }, null, { timeout: 120000 });

  for (let i = 0; i < CASES.length; i++) {
    const [name, builder] = CASES[i];
    const desc = descs[i];
    const page = await ctx.newPage();
    const errors = [];
    page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('pageerror', (e) => errors.push(String(e)));
    try {
      check(desc.variant === builder && desc.fallback === false, `${name}: descriptor builder ${builder}`, `variant=${desc.variant}`);
      await page.goto(`${BASE}/viz3d-lab.html?scene=band-nitride&only=1`, { waitUntil: 'load', timeout: 60000 });
      await wait(page);
      errors.length = 0;
      await mount(page, desc, 'dark');
      await wait(page);
      const st = await page.evaluate(() => { const v = document.getElementById('view'); return { mode: v.dataset.mode, error: v.dataset.error || null }; });
      check(!st.error, `${name}: mounted`, st.error || '');
      check(st.mode === 'webgl', `${name}: WebGL mode`, `mode=${st.mode}`);
      await page.waitForTimeout(900);
      const s = await page.evaluate(() => document.getElementById('view').__viz3d._sample());
      check(s && s.lit > s.samples * 0.03 && s.distinct >= 8, `${name}: canvas not blank (pixel variance)`,
        s ? `lit ${s.lit}/${s.samples}, distinct colours ${s.distinct}` : 'no sample');
      const b = await page.evaluate(() => document.getElementById('view').__viz3d._stats());
      check(b && b.calls < 160 && b.triangles < 600000, `${name}: budget`, b ? `${b.calls} draw calls, ${b.triangles} triangles` : 'no stats');
      const shot = await page.screenshot({ path: path.join(OUTDIR, `viz3d-${name}.png`) });
      const lum = await cornerLum(page, shot);
      check(lum.every((v) => v < 45), `${name}: dark enclosure background (corners)`, lum.map((v) => v.toFixed(0)).join(','));
      await page.focus('#view');
      for (const k of ['2', '3', '4', '1']) { await page.keyboard.press(k); await page.waitForTimeout(700); }
      const labels = await page.evaluate(() => document.querySelectorAll('#view .v3-label').length);
      check(labels >= 4 && labels <= 30, `${name}: labelled (4..30 callouts)`, `${labels}`);
      check(errors.length === 0, `${name}: no console errors`, errors.slice(0, 3).join(' | '));
    } catch (err) {
      check(false, `${name}: run`, String(err).slice(0, 300));
    }
    await page.close();
  }

  // light theme keeps the dark enclosure (one micropillar + one open cavity)
  for (const idx of [0, 1]) {
    const [name] = CASES[idx];
    const page = await ctx.newPage();
    try {
      await page.goto(`${BASE}/viz3d-lab.html?scene=band-nitride&only=1`, { waitUntil: 'load', timeout: 60000 });
      await wait(page);
      await mount(page, descs[idx], 'light');
      await wait(page);
      await page.waitForTimeout(700);
      const shot = await page.screenshot();
      const lum = await cornerLum(page, shot);
      check(lum.every((v) => v < 45), `${name}: dark enclosure kept in the light theme`, lum.map((v) => v.toFixed(0)).join(','));
    } catch (err) {
      check(false, `${name}: light theme`, String(err).slice(0, 300));
    }
    await page.close();
  }

  // unsupported structure kind (edge_ridge on a 530 nm dot: no covered stack, a declared incompatibility): fallback descriptor -> message panel, no error
  {
    const desc = descs[CASES.length];
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', (e) => errors.push(String(e)));
    try {
      check(desc.fallback === true && desc.variant === 'fallback2d', 'edge_ridge: fallback descriptor', `variant=${desc.variant}`);
      await page.goto(`${BASE}/viz3d-lab.html?scene=band-nitride&only=1`, { waitUntil: 'load', timeout: 60000 });
      await wait(page);
      await mount(page, desc, 'dark');
      await wait(page);
      const st = await page.evaluate(() => { const v = document.getElementById('view'); return { mode: v.dataset.mode, msg: !!v.querySelector('.v3-message') }; });
      check(st.mode === 'unavailable' && st.msg && errors.length === 0, 'edge_ridge: message panel, no 3D, no error', `mode=${st.mode}`);
    } catch (err) {
      check(false, 'edge_ridge: run', String(err).slice(0, 300));
    }
    await page.close();
  }


  // ---- R2-3D: every structure kind + paper cards: WebGL dark, WebGL light, forced SVG fallback
  {
    const r2 = describeR2();
    const cases = [...R2_KINDS.map((c, i) => [c[0], c[1], r2[i]]), ...R2_CARDS.map((c, i) => [c[0], c[1], r2[R2_KINDS.length + i]])];
    const distinct = (page, buf) => page.evaluate(async (b64) => {
      const img = new Image();
      img.src = 'data:image/png;base64,' + b64;
      await img.decode();
      const c = document.createElement('canvas');
      c.width = img.width; c.height = img.height;
      const g2 = c.getContext('2d');
      g2.drawImage(img, 0, 0);
      const seen = new Set();
      const d = g2.getImageData(Math.floor(img.width * 0.25), Math.floor(img.height * 0.2), Math.floor(img.width * 0.5), Math.floor(img.height * 0.6)).data;
      for (let i = 0; i < d.length; i += 4 * 7) seen.add(((d[i] >> 4) << 8) | ((d[i + 1] >> 4) << 4) | (d[i + 2] >> 4));
      return seen.size;
    }, buf.toString('base64'));
    const fresh = async (page, desc, theme, opts) => {
      await page.evaluate(async ({ desc, theme, opts }) => {
        const { mountScene } = await import('/js/viz3d/index.js');
        const old = document.getElementById('view');
        const view = document.createElement('div');
        view.id = 'view'; view.dataset.ready = '0';
        view.style.cssText = old.style.cssText;
        old.replaceWith(view);
        try { window.__lab?.api?.dispose(); } catch (e) { /* detached */ }
        try { window.__api = mountScene(view, 'device', desc, { theme, ...opts }); } catch (e) { view.dataset.error = String(e); }
      }, { desc, theme, opts });
      await wait(page);
    };
    for (const [name, builder, desc] of cases) {
      const page = await ctx.newPage();
      const errors = [];
      page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
      page.on('pageerror', (e) => errors.push(String(e)));
      try {
        check(desc.variant === builder && desc.fallback === false, `r2 ${name}: descriptor builder ${builder}`, `variant=${desc.variant} ${desc.reason || ''}`);
        await page.goto(`${BASE}/viz3d-lab.html?scene=band-nitride&only=1`, { waitUntil: 'load', timeout: 60000 });
        await wait(page);
        errors.length = 0;
        await fresh(page, desc, 'dark', {});
        const st = await page.evaluate(() => { const v = document.getElementById('view'); return { mode: v.dataset.mode, ready: v.dataset.ready, error: v.dataset.error || null }; });
        check(!st.error && st.ready === '1' && st.mode === 'webgl', `r2 ${name}: WebGL data-ready, data-mode=webgl`, `mode=${st.mode} ready=${st.ready} ${st.error || ''}`);
        await page.waitForTimeout(900);
        const sm = await page.evaluate(() => document.getElementById('view').__viz3d._sample());
        check(sm && sm.lit >= 6 && sm.distinct >= 8, `r2 ${name}: WebGL canvas not blank (a 220 nm wire is a thin line at true scale)`, sm ? `lit ${sm.lit}/${sm.samples}, distinct ${sm.distinct}` : 'no sample');
        await page.screenshot({ path: path.join(OUTDIR, `r2-3d-${name}.png`) });
        const labs = await page.evaluate(() => document.querySelectorAll('#view .v3-label').length);
        check(labs >= 4 && labs <= 30, `r2 ${name}: labelled (4..30 callouts)`, `${labs}`);
        await page.focus('#view');
        for (const k of ['2', '3', '4', '1']) { await page.keyboard.press(k); await page.waitForTimeout(600); }
        check(errors.length === 0, `r2 ${name}: no console errors (WebGL dark)`, errors.slice(0, 3).join(' | '));
        await fresh(page, desc, 'light', {});
        await page.waitForTimeout(700);
        const shotL = await page.screenshot({ path: path.join(OUTDIR, `r2-3d-${name}-light.png`) });
        const lum = await cornerLum(page, shotL);
        check(lum.every((v) => v < 45), `r2 ${name}: dark enclosure kept in the light theme`, lum.map((v) => v.toFixed(0)).join(','));
        check(errors.length === 0, `r2 ${name}: no console errors (WebGL light)`, errors.slice(0, 3).join(' | '));
        await fresh(page, desc, 'dark', { forceFallback: true });
        const fb = await page.evaluate(() => {
          const v = document.getElementById('view');
          const svg = v.querySelector('.v3-fallback svg');
          return { mode: v.dataset.mode, ready: v.dataset.ready, shapes: svg ? svg.querySelectorAll('rect,polygon,ellipse,circle,line,path,polyline').length : 0, texts: svg ? svg.querySelectorAll('text').length : 0 };
        });
        check(fb.mode === 'fallback' && fb.ready === '1' && fb.shapes >= 4 && fb.texts >= 2, `r2 ${name}: forced fallback draws an SVG`, `mode=${fb.mode} shapes=${fb.shapes} texts=${fb.texts}`);
        const shotF = await page.screenshot({ path: path.join(OUTDIR, `r2-3d-${name}-fallback.png`) });
        const nd = await distinct(page, shotF);
        check(nd >= 8, `r2 ${name}: fallback not blank (pixel variance)`, `${nd} colour bins`);
        check(errors.length === 0, `r2 ${name}: no console errors (fallback)`, errors.slice(0, 3).join(' | '));
      } catch (err) {
        check(false, `r2 ${name}: run`, String(err).slice(0, 300));
      }
      await page.close();
    }
  }

  await browser.close();
  const n = results.filter(Boolean).length;
  console.log(`${n}/${results.length} viz3d_compose checks passed`);
  process.exit(n === results.length ? 0 : 1);
})().catch((e) => { console.error(e); process.exit(1); });
