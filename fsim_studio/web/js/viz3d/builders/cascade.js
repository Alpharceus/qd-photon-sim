// XX -> X -> 0 cascade (cascade / rectangular | deterministic_pair).
// Level glow follows the populations P(t) from fsim_core.lindblad.evolve
// (log scale, stated).  Photon / escape sprites are an illustrative seeded
// sampling of the computed instantaneous rates gamma*P(t), k*P(t), boosted
// by stated factors.  The g2 shown is the evaluated g2_op.
import * as THREE from 'three';
import { glowMaterial, fatLine, nice } from './common.js';
import { lambdaToRGB, rgbCss } from '../core/colormap.js';
import { tagChip } from '../core/labels.js';
import { fmtNum } from '../core/hud.js';
import { exaggerate } from '../core/scale.js';

function mulberry32(a) {
  return function rnd() {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function interp(ts, ys, t) {
  let lo = 0, hi = ts.length - 1;
  if (t <= ts[0]) return ys[0];
  if (t >= ts[hi]) return ys[hi];
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ts[m] <= t) lo = m; else hi = m; }
  const f = (t - ts[lo]) / (ts[hi] - ts[lo] || 1);
  return ys[lo] + (ys[hi] - ys[lo]) * f;
}

export function buildCascade(spec, ctx) {
  const { stage, hud, labels, badges, reducedMotion } = ctx;
  const o = spec.overlays;
  const root = new THREE.Group();
  const lv = o.levels;
  const EX = Number.isFinite(lv.E_X_eV) ? lv.E_X_eV : 1.5;
  const u = 3.0;                                    // world units per E_X
  const gap = exaggerate(badges, 'XX-X gap', lv.gap_scale, lv.gap_badge);
  const yG = 0, yX = u, yXX = 2 * u - (lv.delta_xx_meV / 1000) * (u / EX) * gap;
  const lamX = 1239.841984 / EX;                    // unit conversion only (hc, CODATA)
  const cX = lambdaToRGB(lamX);
  const cXX = [0.62, 0.52, 0.98];                  // distinct hue for XX photons (not a lambda colour)
  const glowOf = (P) => Math.max(0, Math.min(1, (Math.log10(Math.max(P, 1e-30)) + 6) / 6));

  const levelDefs = [['G', '0', yG, [0.75, 0.78, 0.84], 'P_G'], ['X', 'X', yX, cX.rgb, 'P_X'], ['XX', 'XX', yXX, cXX, 'P_XX']];
  const levels = {};
  for (const [key, name, y, rgb, pk] of levelDefs) {
    const g = new THREE.Group();
    g.position.y = y;
    const plate = new THREE.Mesh(new THREE.CylinderGeometry(1.7, 1.7, 0.07, 96),
      new THREE.MeshStandardMaterial({ color: 0x1b1e24, emissive: new THREE.Color(...rgb), emissiveIntensity: 0, roughness: 0.35, metalness: 0.3 }));
    const halo = new THREE.Mesh(new THREE.CircleGeometry(2.1, 96), glowMaterial(rgb, { opacity: 0 }));
    halo.rotation.x = -Math.PI / 2; halo.position.y = 0.05;
    const rim = new THREE.Mesh(new THREE.TorusGeometry(1.72, 0.02, 8, 128), new THREE.MeshBasicMaterial({ color: 0x5a5e66 }));
    rim.rotation.x = Math.PI / 2;
    if (key !== 'G') { stage.markBloom(plate, 2.5); stage.markBloom(halo, 1.5); }
    g.add(plate, halo, rim);
    root.add(g);
    const lab = labels.add(root, [1.95, y, 0], `${name}`, { tag: spec.scalars.g2_op.tag, sub: 'P = ' });
    levels[key] = { plate, halo, lab, pk, rgb };
  }
  // Transition rails: XX -> X (XX photon), X -> G (X photon); pump from the left.
  root.add(fatLine([[0, yXX, 0], [0, yX, 0]], { color: new THREE.Color(...cXX).getHex(), width: 1.4, opacity: 0.6, dashed: true }));
  root.add(fatLine([[0, yX, 0], [0, yG, 0]], { color: new THREE.Color(...cX.rgb).getHex(), width: 1.4, opacity: 0.6, dashed: true }));
  const resX = -6.2;
  const res = new THREE.Mesh(new THREE.BoxGeometry(0.12, yXX + 1.2, 2.4), new THREE.MeshStandardMaterial({ color: 0x2a2e36, roughness: 0.6, transparent: true, opacity: 0.55 }));
  res.position.set(resX, yXX / 2, 0);
  root.add(res);
  labels.add(root, [resX, -0.9, 0], 'carrier reservoir', { kind: 'note', anchor: 'center', sub: 'pump r enters, escape k returns' });
  labels.add(root, [0, yXX + 1.2, 0], `XX-X gap drawn ×${gap}`, { kind: 'chip', tag: spec.scalars.delta_xx.tag, anchor: 'center', sub: `δxx ${nice(lv.delta_xx_meV)} meV; X at E_X ${nice(EX, 4)} eV` });
  const pump = fatLine([[resX + 0.3, yG + 0.4, 0], [-1.8, yX - 0.1, 0]], { color: 0xe9e8e3, width: 1.6, opacity: 0.0 });
  root.add(pump);

  // Sprites: one instanced pool, additive.
  const POOL = 160;
  const sprGeo = new THREE.SphereGeometry(0.075, 10, 8);
  const sprMat = new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, toneMapped: false });
  const spr = new THREE.InstancedMesh(sprGeo, sprMat, POOL);
  spr.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  for (let i = 0; i < POOL; i++) spr.setColorAt(i, new THREE.Color(1, 1, 1));
  spr.count = 0;
  stage.markBloom(spr);
  root.add(spr);
  const parts = [];
  const rnd = mulberry32(20261007);

  // HUD: big evaluated g2, live speed label, timeline strip.
  const big = document.createElement('div');
  big.className = 'v3-bigread';
  const g2 = spec.scalars.g2_op;
  big.innerHTML = '<div class="v3-bigread-label">g2(0), evaluated g2_op</div><div class="v3-bigread-value"></div><div class="v3-note">never counted from the animation</div><div class="v3-speed"></div>';
  big.querySelector('.v3-bigread-value').append(document.createTextNode(fmtNum(g2.value)), tagChip(g2.tag));
  hud.el.appendChild(big);
  const speedEl = big.querySelector('.v3-speed');
  const tl = document.createElement('div');
  tl.className = 'v3-timeline';
  const cv = document.createElement('canvas');
  tl.appendChild(cv);
  hud.el.appendChild(tl);

  const T = o.t_ns, period = o.period_ns, tauOn = o.tau_on_ns;
  const segFrac = 0.32;
  const xOfT = (t, W) => (t <= tauOn ? (t / tauOn) * segFrac * W : (segFrac + ((t - tauOn) / (period - tauOn)) * (1 - segFrac)) * W);
  const drawTimeline = (tNow) => {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const W = cv.clientWidth, H = cv.clientHeight;
    if (cv.width !== Math.round(W * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
    const c = cv.getContext('2d');
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, W, H);
    const pad = { l: 44, r: 8, t: 6, b: 18 };
    const w = W - pad.l - pad.r, h = H - pad.t - pad.b;
    const yOf = (P) => pad.t + h * (1 - (Math.log10(Math.max(P, 1e-12)) + 12) / 12);
    c.font = '500 11px Barlow, "Segoe UI", sans-serif';
    c.fillStyle = '#85837c';
    c.strokeStyle = '#2c2f35'; c.lineWidth = 1;
    for (const e of [0, -3, -6, -9, -12]) {
      const y = yOf(Math.pow(10, e));
      c.beginPath(); c.moveTo(pad.l, y); c.lineTo(pad.l + w, y); c.stroke();
      c.fillText(`1e${e}`, 4, y + 3);
    }
    const xs = pad.l + segFrac * w;
    c.strokeStyle = '#6b6e75'; c.setLineDash([3, 3]);
    c.beginPath(); c.moveTo(xs, pad.t); c.lineTo(xs, pad.t + h); c.stroke(); c.setLineDash([]);
    c.fillStyle = '#b9b8af';
    c.fillText(`pump on 0-${nice(tauOn)} ns`, pad.l + 4, H - 4);
    c.fillText(`dark ${nice(tauOn)}-${nice(period)} ns (own linear scale)`, xs + 6, H - 4);
    for (const [key, col] of [['P_G', [0.75, 0.78, 0.84]], ['P_X', cX.rgb], ['P_XX', cXX]]) {
      c.strokeStyle = rgbCss(col); c.lineWidth = 1.6;
      c.beginPath();
      for (let i = 0; i < T.length; i++) {
        const x = pad.l + xOfT(T[i], w), y = yOf(o[key][i]);
        if (i === 0) c.moveTo(x, y); else c.lineTo(x, y);
      }
      c.stroke();
    }
    const px = pad.l + xOfT(tNow, w);
    c.strokeStyle = '#ecebe6'; c.lineWidth = 1.2;
    c.beginPath(); c.moveTo(px, pad.t); c.lineTo(px, pad.t + h); c.stroke();
    c.fillStyle = '#ecebe6';
    c.fillText(`t = ${fmtNum(tNow)} ns`, Math.min(px + 4, pad.l + w - 90), pad.t + 10);
  };

  const setLevels = (t) => {
    for (const k of Object.keys(levels)) {
      const L = levels[k];
      const P = interp(T, o[L.pk], t);
      const gl = glowOf(P);
      // The ground state (empty dot) does not emit: it gets no glow.
      L.plate.material.emissiveIntensity = k === 'G' ? 0.02 : (0.05 + 1.6 * gl) * (L.plate.material.userData.hdr || 1);
      L.halo.material.opacity = k === 'G' ? 0 : 0.28 * gl * gl;
      const sub = L.lab?.element.querySelector('.v3-label-sub');
      if (sub) sub.textContent = `P = ${fmtNum(P)}  (glow log10, floor 1e-6)`;
    }
    pump.material.opacity = t <= tauOn ? 0.85 : 0.0;
  };

  let tSim = reducedMotion ? tauOn * 0.9 : 0;
  let playing = !reducedMotion;
  const speedAt = (t) => (o.playback.find((p) => t >= p.t0_ns && t < p.t1_ns) || o.playback[o.playback.length - 1]).s_per_ns;
  const spawn = (kind) => {
    if (parts.length >= POOL) return;
    if (kind === 'esc') {
      const fromXX = rnd() < 0.5 && interp(T, o.P_XX, tSim) * 2 > interp(T, o.P_X, tSim);
      const y0 = fromXX ? yXX : yX;
      parts.push({ x: -1.6, y: y0, z: (rnd() - 0.5) * 1.2, vx: -2.2 - rnd(), vy: -0.4 * rnd(), life: 1.6, age: 0, col: [0.55, 0.57, 0.6] });
    } else {
      const y0 = kind === 'XX' ? (yXX + yX) / 2 : (yX + yG) / 2;
      const ang = (rnd() - 0.5) * 0.5;
      parts.push({ x: 0.2, y: y0, z: (rnd() - 0.5) * 0.8, vx: 3.4 * Math.cos(ang), vy: 3.4 * Math.sin(ang), life: 2.4, age: 0, col: kind === 'XX' ? cXX : cX.rgb });
    }
  };
  const M = new THREE.Matrix4(), Q = new THREE.Quaternion(), V = new THREE.Vector3(), Sc = new THREE.Vector3(), C = new THREE.Color();
  const animate = (dt) => {
    if (!playing) return false;
    const sp = speedAt(tSim);                  // s of wall time per ns
    const dtNs = dt / sp;
    const t0 = tSim;
    tSim += dtNs;
    if (tSim >= period) { tSim -= period; }
    for (const [kind, key] of [['X', 'rate_X'], ['XX', 'rate_XX'], ['esc', 'rate_esc']]) {
      const lam = (o.sprite_boost[kind] || 0) * interp(T, o[key], t0) * dtNs;
      let n = Math.floor(lam);
      if (rnd() < lam - n) n += 1;
      for (let i = 0; i < Math.min(n, 12); i++) spawn(kind);
    }
    for (let i = parts.length - 1; i >= 0; i--) {
      const p = parts[i];
      p.age += dt; p.x += p.vx * dt; p.y += p.vy * dt;
      if (p.age > p.life) parts.splice(i, 1);
    }
    spr.count = parts.length;
    parts.forEach((p, i) => {
      const a = 1 - p.age / p.life;
      M.compose(V.set(p.x, p.y, p.z), Q, Sc.setScalar(0.6 + 0.8 * a));
      spr.setMatrixAt(i, M);
      spr.setColorAt(i, C.setRGB(p.col[0] * a * 1.6, p.col[1] * a * 1.6, p.col[2] * a * 1.6));
    });
    spr.instanceMatrix.needsUpdate = true;
    if (spr.instanceColor) spr.instanceColor.needsUpdate = true;
    setLevels(tSim);
    speedEl.textContent = `1 ns = ${fmtNum(sp)} s (${tSim <= tauOn ? 'pump on' : 'dark'})`;
    drawTimeline(tSim);
    return true;
  };
  setLevels(tSim);
  speedEl.textContent = reducedMotion ? 'animation off (prefers-reduced-motion); static frame inside the pump window'
    : `1 ns = ${fmtNum(speedAt(0))} s (pump on)`;
  requestAnimationFrame(() => drawTimeline(tSim));
  if (!reducedMotion) stage.addAnimator(animate);

  const gp = hud.group('playback');
  const btn = hud.button(gp, playing ? 'pause' : 'play', () => {
    playing = !playing;
    btn.textContent = playing ? 'pause' : 'play';
    if (playing && reducedMotion) stage.addAnimator(animate);
    stage.requestRender();
  });
  hud.button(gp, 'restart', () => { tSim = 0; parts.length = 0; spr.count = 0; setLevels(0); drawTimeline(0); stage.requestRender(); });
  stage.onResize = [...(stage.onResize || []), () => drawTimeline(tSim)];

  const cy = (yXX + yG) / 2;
  return {
    root,
    defaultView: { pos: [9.5, cy + 5.5, 19], target: [1.6, cy - 0.4, 0] },
    bookmarks: [
      { name: 'Cascade', pos: [9.5, cy + 5.5, 19], target: [1.6, cy - 0.4, 0] },
      { name: 'Side', pos: [15, cy + 0.5, 0.5], target: [0, cy, 0] },
      { name: 'Levels close-up', pos: [3.2, yXX + 1.2, 5.5], target: [0, (yXX + yX) / 2, 0] },
      { name: 'Top', pos: [0.5, cy + 14, 2], target: [0, cy, 0] },
    ],
    noScaleBar: true,
    readoutKeys: ['mean_counts_x', 'r_ns', 'gamma_X_ns', 'gamma_XX_ns', 'k_X_ns', 'tau_on_ns', 'period_ns', 'delta_xx'],
    extraNotes: [`glow: ${o.glow_scale.label}`, `start state: ${o.start_state}`],
  };
}
