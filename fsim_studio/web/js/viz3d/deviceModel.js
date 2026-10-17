// Device model for the Papers detail panel and the Library cards (R2): the composed 3D scene of a stored design, built by
// the same server-side descriptor (POST /api/compose/scene) and the same viz3d builders as the Composer.
//
//   mountDeviceModel(host, source, { api, theme })  a live, orbitable view in `host`; loaded lazily (when the host
//                                                   scrolls into view); at most ONE live WebGL context at a time (a
//                                                   new mount disposes the previous one)
//   deviceThumbnail(source, { api, theme })         a static PNG data URL of the same scene, rendered one at a time
//                                                   on a short-lived offscreen context (the Library grid)
//
// source: {paper: "<id>"} (GET /api/papers/<id>/design, the central values of every ranged leaf), {spec: <composed spec>}
// or {card: "<parameter card>"} (GET /api/scene/device?card=: the display-only scene block of a Library parameter card;
// null when the card has no scene block).
// A design whose structure has no 3D builder gives the descriptor's own reason, shown as text: never an empty box.
import { h, clear } from '../ui/dom.js';

let live = null;            // {host, api}: the one live view
let queue = Promise.resolve();

let _papers = null;
// the ids of the paper cards (only those have a composed design); asked once
async function paperIds(api) {
  if (!_papers) _papers = api.getJSON('/api/papers').then((l) => new Set((Array.isArray(l) ? l : l?.papers || []).map((x) => x.id))).catch(() => new Set());
  return _papers;
}

async function sceneOf(api, source) {
  if (source.card) {
    try { return await api.getJSON(`/api/scene/device?card=${encodeURIComponent(source.card)}`); } catch (e) {
      if ([400, 404, 422].includes(e?.status)) return null;       // no scene block: the caller states the reason
      throw e;
    }
  }
  let spec = source.spec;
  if (!spec && source.paper) {
    const r = await api.getJSON(`/api/papers/${encodeURIComponent(source.paper)}/design`);
    spec = r?.design?.composed;
  }
  if (!spec) throw new Error('no composed design for this card');
  return api.composeScene({ schema: 'fsim-composed/1', ...spec });
}

const drawable = (s) => s && !s.fallback && s.builder && s.builder !== 'fallback2d' && s.variant !== 'fallback2d';

function note(host, text) {
  clear(host);
  host.append(h('p.dm-note', text));
  host.dataset.state = 'note';
}

export function disposeLive() {
  if (live) { try { live.api.dispose(); } catch { /* gone */ } live.host.dataset.state = 'idle'; live = null; }
}

export function mountDeviceModel(host, source, { api, theme = 'dark' } = {}) {
  host.dataset.state = 'pending';
  let dead = false;
  const start = async () => {
    try {
      const scene = await sceneOf(api, source);
      if (dead || !host.isConnected) return;
      if (!scene) { note(host, 'No 3D model: this card has no scene block.'); return; }
      if (!drawable(scene)) { note(host, scene.reason ? `No 3D view: ${scene.reason}` : 'No 3D view for this device'); return; }
      const mod = await import('./index.js');
      if (dead || !host.isConnected) return;
      disposeLive();
      clear(host);
      const view = h('div.dm-view');
      host.append(view);
      const a = mod.mountScene(view, 'device', scene, { theme });
      live = { host, api: a };
      host.dataset.state = 'live';
      host.dataset.variant = scene.variant || '';
      host.dataset.drive = scene.drive?.kind || '';
      if (view.dataset.mode === 'unavailable') { disposeLive(); note(host, 'No 3D view for this device'); }
    } catch (e) {
      note(host, `3D view unavailable: ${e.message || e}`);
    }
  };
  if ('IntersectionObserver' in window) {
    const io = new IntersectionObserver((es) => { if (es.some((x) => x.isIntersecting)) { io.disconnect(); start(); } });
    io.observe(host);
    return { dispose() { dead = true; io.disconnect(); if (live?.host === host) disposeLive(); } };
  }
  start();
  return { dispose() { dead = true; if (live?.host === host) disposeLive(); } };
}

// One thumbnail at a time, each on its own short-lived context (the grid never holds more than one).
export function deviceThumbnail(source, { api, theme = 'dark', width = 480, height = 300 } = {}) {
  const job = async () => {
    if (source.paper && !(await paperIds(api)).has(source.paper)) return null;     // a parameter card without a composed design: no picture
    const scene = await sceneOf(api, source);
    if (!scene || !drawable(scene)) return null;
    const mod = await import('./index.js');
    const box = h('div', { style: `position:fixed;left:-10000px;top:0;width:${width}px;height:${height}px;` });
    document.body.append(box);
    try {
      const a = mod.mountScene(box, 'device', scene, { theme, reducedMotion: true });
      await new Promise((res) => {
        const t0 = performance.now();
        const tick = () => (box.dataset.ready === '1' || performance.now() - t0 > 15000 ? res() : requestAnimationFrame(tick));
        tick();
      });
      await new Promise((r) => setTimeout(r, 250));
      const url = box.dataset.mode === 'webgl' ? await a.exportPNG(1) : null;
      a.dispose();
      return url;
    } finally {
      box.remove();
    }
  };
  const p = queue.then(job, job);
  queue = p.catch(() => null);
  return p;
}
