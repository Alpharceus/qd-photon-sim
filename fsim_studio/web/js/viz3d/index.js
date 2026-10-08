// FSIM Studio 3D views: entry point.
//
//   mountScene(el, kind, spec, opts) -> {update(spec), setTheme(mode),
//                                        bookmark(n), exportPNG(scale), dispose(),
//                                        setGlow(0..1), setPlayhead({T, glow})}
//
// kind: 'device' | 'band' | 'cascade' | 'surface' | 'lattice'; spec is the
// SceneSpec JSON from GET /api/scene/<kind> (schema: fsim_studio/scene.py).
// opts: {theme: 'dark'|'light', reducedMotion?: bool (default: the
// prefers-reduced-motion media query), forceFallback?: bool, bloom?: bool,
// onReady?: fn(api)}.  The container gets data-ready="1" after the first
// frame and data-mode="webgl"|"fallback".  WebGL2 missing or a lost context
// switches to the SVG fallback drawn from the same spec.
import WebGL from 'three/addons/capabilities/WebGL.js';
import { Stage, BLOOM_LAYER } from './core/stage.js';
import { Hud } from './core/hud.js';
import { Labels } from './core/labels.js';
import { axesBadges } from './core/scale.js';
import { ScaleBar } from './core/scale.js';
import { updateLineResolution } from './builders/common.js';
import { buildRidge } from './builders/ridge.js';
import { buildNanowire } from './builders/nanowire.js';
import { buildPlanarCavity } from './builders/planarCavity.js';
import { buildBandProfile } from './builders/bandProfile.js';
import { buildCascade } from './builders/cascade.js';
import { buildSurface } from './builders/surface.js';
import { buildLattice } from './builders/lattice.js';
import { renderFallback } from './fallback2d.js';
import * as THREE from 'three';
import { addFloor, solidBox } from './builders/floor.js';

const KINDS = ['device', 'band', 'cascade', 'surface', 'lattice'];

let cssInjected = false;
function injectCss() {
  if (cssInjected) return;
  cssInjected = true;
  const link = document.createElement('link');
  link.rel = 'stylesheet';
  link.href = new URL('./viz3d.css', import.meta.url).href;
  document.head.appendChild(link);
}

function pickBuilder(kind, spec) {
  if (kind === 'device') {
    if (spec.variant === 'edge_ridge') return buildRidge;
    if (spec.variant?.startsWith('nanowire')) return buildNanowire;
    if (spec.variant === 'planar_cavity' || spec.variant === 'thermal_mesa') return buildPlanarCavity;
  }
  if (kind === 'band' && (spec.variant === 'nitride' || spec.variant === 'inp')) return buildBandProfile;
  if (kind === 'cascade' && spec.overlays?.t_ns) return buildCascade;
  if (kind === 'surface' && spec.overlays?.lo) return buildSurface;
  if (kind === 'lattice' && spec.overlays?.rows) return buildLattice;
  return null;
}

function reducedMotionDefault() {
  return !!window.matchMedia?.('(prefers-reduced-motion: reduce)').matches;
}

// Playhead glow (response animation): scales the emission glow of a built 'device' scene
// without rebuilding it. Only the opacity of additive glow materials and the emissive
// intensity of bloom-marked meshes are touched (base values are remembered once per
// material), plus the bloom strength; one requestRender, no new WebGL context. glow = 1
// is the scene exactly as built; 0 keeps GLOW_FLOOR of it so the device never goes dark.
const GLOW_FLOOR = 0.2;
function glowTargets(state) {
  if (state.glowMats) return state.glowMats;
  const seen = new Set();
  const out = [];
  state.built.root.traverse((o) => {
    if (!o.layers.isEnabled(BLOOM_LAYER)) return;
    for (const m of Array.isArray(o.material) ? o.material : [o.material]) {
      if (!m || seen.has(m)) continue;
      seen.add(m);
      if (m.isMeshBasicMaterial && m.blending === THREE.AdditiveBlending) out.push({ m, k: 'opacity', base: m.opacity });
      else if (m.emissive && m.emissiveIntensity !== undefined) out.push({ m, k: 'emissiveIntensity', base: m.emissiveIntensity });
    }
  });
  state.glowBloom = state.stage.bloomPass.strength;
  state.glowMats = out;
  return out;
}

function messagePanel(el, title, lines) {
  const box = document.createElement('div');
  box.className = 'v3-message';
  const h = document.createElement('div');
  h.className = 'v3-message-title';
  h.textContent = title;
  box.appendChild(h);
  for (const l of lines) {
    const p = document.createElement('div');
    p.className = 'v3-note';
    p.textContent = l;
    box.appendChild(p);
  }
  el.appendChild(box);
  return box;
}

export function mountScene(el, kind, spec, opts = {}) {
  if (!KINDS.includes(kind)) throw new Error(`mountScene: unknown kind ${kind}`);
  injectCss();
  const reducedMotion = opts.reducedMotion ?? reducedMotionDefault();
  el.classList.add('v3-root');
  el.dataset.theme = opts.theme || 'dark';
  if (!el.hasAttribute('tabindex')) el.tabIndex = 0;
  el.dataset.ready = '0';

  let state = null;   // {stage, hud, labels, built, scalebar} or {fallback}
  let current = spec;

  const teardown = () => {
    if (!state) return;
    el.removeEventListener('keydown', state.onKey);
    state.stage?.dispose();
    state.labels?.dispose();
    state.hud?.dispose();
    state.scalebar?.dispose();
    state.panel?.remove();
    state.fallback?.remove();
    state = null;
  };

  const useFallback = (s, why) => {
    teardown();
    const node = renderFallback(el, kind, s, { reason: why });
    state = { fallback: node, onKey: () => {} };
    el.dataset.mode = 'fallback';
    el.dataset.ready = '1';
  };

  const build = (s) => {
    teardown();
    if (s && s.job) {
      const panel = messagePanel(el, 'Computing in the background', [
        s.reason || 'this view needs a background job',
        `estimated ${Math.round(s.eta_s || 0)} s; the view fills when the job finishes`]);
      state = { panel, onKey: () => {} };
      el.dataset.mode = 'pending';
      el.dataset.ready = '1';
      return;
    }
    const builder = pickBuilder(kind, s || {});
    if (!builder) {
      const panel = messagePanel(el, s?.title || kind, (s?.labels && s.labels.length) ? s.labels : ['no 3D view for this card']);
      state = { panel, onKey: () => {} };
      el.dataset.mode = 'unavailable';
      el.dataset.ready = '1';
      return;
    }
    if (opts.forceFallback || !WebGL.isWebGL2Available()) { useFallback(s, opts.forceFallback ? 'forced' : 'WebGL2 unavailable'); return; }

    const stage = new Stage(el, { reducedMotion, bloom: opts.bloom ?? true, ortho: false });
    const hud = new Hud(el, s);
    const labels = new Labels();
    const badges = new Set(axesBadges(s));
    let built;
    try {
      built = builder(s, { stage, hud, labels, badges, reducedMotion, kind });
    } catch (err) {
      stage.dispose(); hud.dispose(); labels.dispose();
      throw err;
    }
    stage.scene.add(built.root);
    // The device's own extent (before the floor scenery joins the root):
    // tiny viewports fit the default view to it.
    const devBox = kind === 'device' ? solidBox(built.root) : null;
    if (kind === 'device' && !built.noFloor) {
      const fl = addFloor(stage, built.root);
      if (fl) built.extraNotes = [...(built.extraNotes || []), fl.note];
    }
    stage.contextLost = () => useFallback(current, 'WebGL context lost');
    hud.setReadouts(s.scalars, built.readoutKeys);
    // Compact (< ~1000 x 500): callouts become numbered pins keyed in a
    // list, the readout plate shows 4 rows + "+N more", and every control
    // group merges into one bar under the canvas.  Narrow (< 560 wide): the
    // readout plate stacks under the title so the two never overlap.  Tiny
    // (< 480 wide or < 300 high): the plate shows 2 rows, the bar wraps onto
    // two rows (no control clipped), and the default view is fitted to the
    // device so it gets most of the canvas.
    const rect0 = () => el.getBoundingClientRect();
    const compact = () => { const r = rect0(); return r.width < 1000 || r.height < 500; };
    const narrow = () => rect0().width < 560;
    const tiny = () => { const r = rect0(); return r.width < 480 || r.height < 300; };
    hud.setNotes([...(s.labels || []), ...(built.extraNotes || [])], [...badges], { open: !compact() });
    const pins = labels.pinned();
    hud.setKey(pins);
    stage.afterCss = () => {
      labels.layout(el, [hud.head, hud.right, hud.controls, hud.notes, hud.banner,
        el.querySelector('.v3-timeline'), el.querySelector('.v3-bigread'),
        el.querySelector('.v3-inset[style*="block"]')], 30);
      hud.syncKey();
    };
    if (s.verdict) hud.setVerdict(s.verdict);

    const scalebar = built.noScaleBar ? null : new ScaleBar(hud.notes, { unitToUm: built.unitToUm ?? 1 });
    if (scalebar) {
      scalebar.el.classList.add('v3-scalebar-main');
      hud.notes.prepend(scalebar.el);
      stage.onFrame.add(() => scalebar.update(stage.pxPerWorld(), built.worldPerUm ?? 1));
    }
    // Camera fit to the real canvas: the authored views are composed for a
    // 1600 x 1000 viewport; a smaller or narrower canvas (and the band the
    // stacked HUD covers when narrow) moves the camera back along the same
    // direction and shifts the projection centre into the free region.
    let distK = 1, lastView = null;
    // Tiny viewports: the default view aims at the device's centre and its
    // distance is scaled so the projected device box fills ~95% of the free
    // region (bookmarks keep their authored framing).
    const devFit = (v) => {
      if (!devBox || devBox.isEmpty() || v !== built.defaultView || !el.classList.contains('v3-tiny')) return null;
      const c = devBox.getCenter(new THREE.Vector3());
      const dir = new THREE.Vector3(...v.pos).sub(new THREE.Vector3(...v.target));
      const cam = stage.camera, corners = [];
      for (const x of [devBox.min.x, devBox.max.x]) for (const y of [devBox.min.y, devBox.max.y]) for (const z of [devBox.min.z, devBox.max.z]) corners.push(new THREE.Vector3(x, y, z));
      const span = (k) => {
        cam.position.copy(c).addScaledVector(dir, k); cam.lookAt(c); cam.updateMatrixWorld(true);
        let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
        for (const p of corners) { const q = p.clone().project(cam); x0 = Math.min(x0, q.x); x1 = Math.max(x1, q.x); y0 = Math.min(y0, q.y); y1 = Math.max(y1, q.y); }
        return Math.max((x1 - x0) / 2, (y1 - y0) / 2 / fitFree);
      };
      const p0 = cam.position.clone(), q0 = cam.quaternion.clone();
      let k = 1;
      for (let i = 0; i < 4; i++) k *= span(k) / 0.95;
      cam.position.copy(p0); cam.quaternion.copy(q0); cam.updateMatrixWorld(true);
      k = Math.min(Math.max(k, 0.05), 20);
      return { pos: [0, 1, 2].map((i) => c.getComponent(i) + dir.getComponent(i) * k), target: [c.x, c.y, c.z] };
    };
    let fitFree = 1, lastFree = 1;
    const applyView = (v, instant = false) => {
      lastView = v;
      const f = devFit(v);
      if (f) { stage.setView(f.pos, f.target, { instant }); return; }
      const t = v.target, p = v.pos;
      const q = [0, 1, 2].map((i) => t[i] + (p[i] - t[i]) * distK);
      stage.setView(q, t, { instant });
    };
    const fit = () => {
      const { w, h } = stage.size;
      let top = 0, bot = 0;
      if (el.classList.contains('v3-narrow')) {
        const wr = stage.wrap.getBoundingClientRect(), rr = hud.right.getBoundingClientRect(), hr = hud.head.getBoundingClientRect();
        const nr = hud.notes.getBoundingClientRect();
        top = Math.max(rr.height ? rr.bottom : 0, hr.bottom) - wr.top + 4;
        if (nr.height && !el.classList.contains('v3-tiny')) bot = wr.bottom - nr.top + 4;
        top = Math.min(Math.max(0, top), h * 0.55);
        bot = Math.min(Math.max(0, bot), h * 0.3);
      }
      const fh = Math.max(40, h - top - bot);
      fitFree = fh / h;
      const k = Math.max(1, h / Math.min(fh, w / 1.6));
      if (top || bot) stage.camera.setViewOffset(w, h, 0, -(top - bot) / 2, w, h);
      else stage.camera.clearViewOffset();
      stage.camera.updateProjectionMatrix();
      const changed = Math.abs(k - distK) > 1e-3 || Math.abs(fitFree - lastFree) > 1e-3;
      distK = k; lastFree = fitFree;
      if (changed && lastView) applyView(lastView, true);
      stage.requestRender();
    };
    const applyCompact = () => {
      const c = compact(), n = narrow(), t = c && tiny();
      el.classList.toggle('v3-compact', c);
      el.classList.toggle('v3-narrow', n);
      el.classList.toggle('v3-tiny', t);
      el.classList.toggle('v3-has-bar', c && hud.controls.children.length > 0);
      el.classList.toggle('v3-keep-legend', !!built.keepLegend);
      hud.setRowLimit(t ? 2 : 4);
      // the canvas ends where the (possibly two-row) bar starts
      const bh = c && hud.controls.children.length ? hud.controls.offsetHeight : 0;
      el.style.setProperty('--v3-bar-h', `${bh || 34}px`);
      // tiny: title and the 2-row plate sit side by side (shorter top band)
      hud.right.style.top = n && !t ? `${hud.head.offsetTop + hud.head.offsetHeight + 6}px` : '';
    };
    stage.onResize = [(w, h) => updateLineResolution(stage.scene, w, h), () => { applyCompact(); fit(); }];
    updateLineResolution(stage.scene, stage.size.w, stage.size.h);
    if (built.inset) updateLineResolution(built.inset, stage.size.w, stage.size.h);
    stage.controls.addEventListener('start', () => { lastView = null; });

    // camera
    if (built.controls) Object.assign(stage.controls, built.controls);

    // bookmarks 1-4 and the cut-away slider
    let cutSlider = null;
    const goBookmark = (n) => {
      const b = built.bookmarks?.[n - 1];
      if (!b) return;
      applyView(b);
      if (built.cut && b.cut !== undefined) {
        stage.clipPlane.constant = b.cut;
        if (cutSlider) cutSlider.value = b.cut;
      }
      if (built.magnifier && b.mag !== undefined) built.magnifier.set(b.mag);
      b.onEnter?.();
      stage.requestRender();
    };
    if (built.bookmarks?.length) {
      const g = hud.group('views');
      built.bookmarks.forEach((b, i) => hud.button(g, b.name, () => goBookmark(i + 1), { key: String(i + 1), title: `camera bookmark ${i + 1}` }));
      hud.controls.prepend(g);
    }
    if (built.cut) {
      stage.clipPlane.constant = built.cut.max + 1e-6;
      const g = hud.group('');
      cutSlider = hud.slider(g, {
        min: built.cut.min, max: built.cut.max + 1e-6, step: (built.cut.max - built.cut.min) / 400,
        value: built.cut.max + 1e-6, label: 'cut-away',
        onInput: (v) => { stage.clipPlane.constant = v; stage.requestRender(); },
      });
    } else {
      stage.clipPlane.constant = 1e9;
    }

    applyCompact();
    stage.resize();          // the bar may have shrunk the canvas
    fit();
    applyView(built.defaultView, true);

    // hover tooltips + click actions
    stage.onHover = (hit) => {
      if (!hit) { hud.showTip(null); return; }
      const info = typeof hit.info === 'function' ? hit.info(hit) : hit.info;
      hud.showTip(info, hit.x, hit.y);
    };
    let down = null;
    stage.css.domElement.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
    stage.css.domElement.addEventListener('pointerup', (e) => {
      if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) return;
      const rect = stage.css.domElement.getBoundingClientRect();
      stage.pointer.set(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
      stage.raycaster.setFromCamera(stage.pointer, stage.camera);
      const hit = stage.raycaster.intersectObjects(stage.pickables, false).find((h) => stage._clipOk(h.point));
      hit?.object.userData.onClick?.(hit);
    });
    stage.onBloomDrop = () => hud.addBadge('bloom off (frame budget)');

    const onKey = (e) => {
      if (e.target && /input|textarea|select/i.test(e.target.tagName)) return;
      const n = Number(e.key);
      if (n >= 1 && n <= 4) { goBookmark(n); e.preventDefault(); }
    };
    el.addEventListener('keydown', onKey);

    state = { stage, hud, labels, built, scalebar, onKey, goBookmark, devBox };
    el.dataset.mode = 'webgl';
    stage.render();
    requestAnimationFrame(() => {
      if (!state || state.stage !== stage) return;
      stage.render();
      el.dataset.ready = '1';
      opts.onReady?.(api);
    });
  };

  const api = {
    update(next) { current = next; build(next); },
    setTheme(mode) { el.dataset.theme = mode === 'light' ? 'light' : 'dark'; state?.stage?.requestRender(); },
    bookmark(n) { state?.goBookmark?.(n); },
    // Playhead hooks: cheap, non-rebuilding. Return false when there is nothing to scale
    // (fallback / message panels, non-device kinds).
    setGlow(g) {
      if (kind !== 'device' || !state?.stage || !state.built) return false;
      const x = Math.max(0, Math.min(1, Number.isFinite(g) ? g : 1));
      if (state.glow === x) return true;
      state.glow = x;
      const f = GLOW_FLOOR + (1 - GLOW_FLOOR) * x;
      for (const t of glowTargets(state)) t.m[t.k] = t.base * f;
      state.stage.bloomPass.strength = state.glowBloom * (0.4 + 0.6 * x);
      state.stage.requestRender();
      return true;
    },
    setPlayhead({ T, glow } = {}) {
      if (Number.isFinite(T)) el.dataset.playheadT = String(T);
      return glow == null ? false : api.setGlow(glow);
    },
    async exportPNG(scale = 2) {
      if (state?.stage) return state.stage.exportPNG(scale);
      return null;
    },
    dispose() { teardown(); el.classList.remove('v3-root'); delete el.dataset.ready; },
    // test hooks (not part of the shell contract)
    _sample() { return state?.stage ? state.stage.sample() : null; },
    _stats() { return state?.stage ? state.stage.stats() : null; },
    _has(name) { return !!state?.stage?.scene.getObjectByName(name); },
    _glow() { return state?.glow ?? null; },
    // number of materials the glow hook scales and their mean current/base ratio
    _glowTargets() {
      if (!state?.stage || !state.built) return null;
      const t = glowTargets(state);
      return { n: t.length, ratio: t.length ? t.reduce((x, y) => x + (y.base ? y.m[y.k] / y.base : 1), 0) / t.length : null };
    },
    _info() { return state?.built?.testInfo ? state.built.testInfo() : null; },
    // projected device box (solid meshes, floor excluded) as a fraction of
    // the canvas area, clipped to the canvas
    _deviceCover() {
      if (!state?.stage || !state.devBox || state.devBox.isEmpty()) return null;
      const cam = state.stage.camera, b = state.devBox;
      cam.updateMatrixWorld(true);
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (const x of [b.min.x, b.max.x]) for (const y of [b.min.y, b.max.y]) for (const z of [b.min.z, b.max.z]) {
        const q = new THREE.Vector3(x, y, z).project(cam);
        x0 = Math.min(x0, q.x); x1 = Math.max(x1, q.x); y0 = Math.min(y0, q.y); y1 = Math.max(y1, q.y);
      }
      const cl = (v) => Math.min(1, Math.max(-1, v));
      return ((cl(x1) - cl(x0)) / 2) * ((cl(y1) - cl(y0)) / 2);
    },
    get mode() { return el.dataset.mode; },
  };
  build(spec);
  el.__viz3d = api;
  return api;
}

export default mountScene;
