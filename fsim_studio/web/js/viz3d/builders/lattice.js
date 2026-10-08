// Sweep lattice (lattice / <campaign>): one InstancedMesh glyph per sweep
// row on the factorial grid; every row is drawn (rows sharing a grid cell
// sit on a small sub-grid inside it).  Colour = g2 (viridis), size = log10
// flux, ring = optical_pass, shell = the device gate.  Rows of a model that
// never gates (rt_edge: anything but the headline model, I4 `gates` false)
// are drawn grey with no ring or shell, whatever their pass column says.
// The banner text and its pass/fail style come only from overlays.banner
// (computed by the backend from gating rows).
import * as THREE from 'three';
import { fatLine, nice } from './common.js';
import { viridis, viridisGradientCss } from '../core/colormap.js';
import { fmtNum } from '../core/hud.js';

const SP = 2.2;
const GREY = new THREE.Color(0x5d6067);

export function buildLattice(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const o = spec.overlays;
  const R = o.rows;
  const n = R.g2.length;
  const A = o.axes;
  const nx = A.x.values.length, ny = A.y.values.length, nz = A.z.values.length;
  const ox = ((nx - 1) * SP) / 2, oz = ((nz - 1) * SP) / 2;
  const pos = (i, j, k) => [i * SP - ox, j * SP, k * SP - oz];
  const root = new THREE.Group();
  const gates = (k) => (R.gates ? R.gates[k] !== false : true);
  const nNonGating = R.gates ? R.gates.filter((g) => g === false).length : 0;

  // Flux -> size on a fixed log scale across ALL rows.
  const lf = R.flux.map((f) => (f > 0 && Number.isFinite(f) ? Math.log10(f) : NaN));
  const fin = lf.filter(Number.isFinite);
  const fMin = fin.length ? Math.min(...fin) : 0, fMax = fin.length ? Math.max(...fin) : 1;
  const size = (k) => (Number.isFinite(lf[k]) ? 0.16 + 0.5 * (lf[k] - fMin) / ((fMax - fMin) || 1) : 0.12);

  const glyph = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(1, 0), new THREE.MeshStandardMaterial({ roughness: 0.35, metalness: 0.1, flatShading: true }), n);
  const bad = new THREE.InstancedMesh(new THREE.OctahedronGeometry(1, 0), new THREE.MeshStandardMaterial({ color: 0x4a4d55, roughness: 0.8 }), n);
  const ring = new THREE.InstancedMesh(new THREE.TorusGeometry(1, 0.07, 4, 20), new THREE.MeshBasicMaterial({ color: 0xf4f1e6, toneMapped: false }), n);
  const nShell = R.device_pass.filter((v, k) => v && gates(k)).length;
  const shell = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(1, 1), new THREE.MeshBasicMaterial({ color: 0xffffff, wireframe: true, transparent: true, opacity: 0.6 }), Math.max(1, nShell));
  for (let k = 0; k < n; k++) glyph.setColorAt(k, new THREE.Color(1, 1, 1));
  stage.markBloom(ring);
  root.add(glyph, bad, ring, shell);

  // Every row is drawn by default: facets start at "all" (null).
  const sel = Object.fromEntries(Object.keys(o.facets).map((f) => [f, null]));
  let gatingOnly = false;
  let shown = [], shownBad = [];
  const info = { rows: n, shown: 0, gating: 0, nonGating: 0, rings: 0, shells: 0, invalid: 0 };
  const M = new THREE.Matrix4(), Q = new THREE.Quaternion(), V = new THREE.Vector3(), S = new THREE.Vector3(), C = new THREE.Color();
  const qRing = new THREE.Quaternion().setFromEuler(new THREE.Euler(Math.PI / 2, 0, 0));
  const countEl = document.createElement('div');
  countEl.className = 'v3-key v3-lat-count';
  const refresh = () => {
    const cells = new Map();
    const pick = [];
    for (let k = 0; k < n; k++) {
      let ok = !(gatingOnly && !gates(k));
      if (ok) for (const [f, v] of Object.entries(sel)) if (v !== null && R.facet[f][k] !== v) { ok = false; break; }
      if (!ok) continue;
      const key = `${R.xi[k]},${R.yi[k]},${R.zi[k]}`;
      const slot = cells.get(key) || 0;
      cells.set(key, slot + 1);
      pick.push([k, slot]);
    }
    // Rows sharing a cell sit on an m x m x m sub-grid inside it.
    const multi = Math.max(...cells.values(), 1);
    const m = Math.ceil(Math.cbrt(multi) - 1e-9);
    const sub = m > 1 ? 1.4 / m : 0;
    const off = (slot) => {
      if (m <= 1) return [0, 0, 0];
      const a = slot % m, b = Math.floor(slot / m) % m, c = Math.floor(slot / (m * m));
      return [(a - (m - 1) / 2) * sub, (c - (m - 1) / 2) * sub, (b - (m - 1) / 2) * sub];
    };
    const shrink = m > 1 ? 0.85 / m : 1;
    shown = []; shownBad = [];
    let nr = 0, ns = 0, ng = 0, nn = 0;
    for (const [k, slot] of pick) {
      const [x, y, z] = pos(R.xi[k], R.yi[k], R.zi[k]);
      const [dx, dy, dz] = off(slot);
      V.set(x + dx, y + dy, z + dz);
      if (!Number.isFinite(R.g2[k])) {
        M.compose(V, Q, S.setScalar(0.1 * Math.max(shrink, 0.5)));
        bad.setMatrixAt(shownBad.length, M);
        shownBad.push(k);
        continue;
      }
      const s = size(k) * shrink;
      const i = shown.length;
      M.compose(V, Q, S.setScalar(s));
      glyph.setMatrixAt(i, M);
      if (gates(k)) {
        ng += 1;
        glyph.setColorAt(i, C.setRGB(...viridis(R.g2[k])));
        if (R.optical_pass[k]) { M.compose(V, qRing, S.setScalar(s * 1.45)); ring.setMatrixAt(nr++, M); }
        if (R.device_pass[k]) { M.compose(V, Q, S.setScalar(s * 1.7)); shell.setMatrixAt(ns++, M); }
      } else {
        nn += 1;
        glyph.setColorAt(i, GREY);
      }
      shown.push(k);
    }
    glyph.count = shown.length; bad.count = shownBad.length; ring.count = nr; shell.count = ns;
    for (const mm of [glyph, bad, ring, shell]) { mm.instanceMatrix.needsUpdate = true; mm.computeBoundingSphere(); }
    if (glyph.instanceColor) glyph.instanceColor.needsUpdate = true;
    Object.assign(info, { shown: shown.length + shownBad.length, gating: ng, nonGating: nn, rings: nr, shells: ns, invalid: shownBad.length });
    countEl.textContent = `${shown.length + shownBad.length} rows drawn · ${ng} gating · ${nn} never gate (grey) · ${nr} optical_pass (gating) · ${shownBad.length} invalid`;
    stage.requestRender();
  };
  const rowInfo = (k) => ({
    title: `row ${k}${gates(k) ? '' : ' · never gates'}`, tag: 'A',
    rows: [
      [A.x.name, A.x.values[R.xi[k]]], [A.y.name, A.y.values[R.yi[k]]], [A.z.name, A.z.values[R.zi[k]]],
      ...Object.keys(R.facet).map((f) => [f, o.facets[f][R.facet[f][k]]]),
      ...(R.model_finite_pulse && R.model_finite_pulse[k] !== null && R.model_finite_pulse[k] !== undefined
        ? [['model_finite_pulse', R.model_finite_pulse[k] ? 'yes' : 'no'], ['model_tau_cap_density', R.model_tau_cap_density[k] ? 'yes' : 'no']] : []),
      ['g2', fmtNum(R.g2[k])], ['flux', `${fmtNum(R.flux[k])} /s`],
      ['optical_pass', R.optical_pass[k] ? 'yes' : 'no'], [o.counts.pass_column, R.device_pass[k] ? 'yes' : 'no'],
      ['eligible', R.eligible[k] ? 'yes' : 'no'],
    ],
    note: R.reasons[k] ? `invalid: ${R.reasons[k]}`
      : gates(k) ? 'committed sweep row (gating model)'
        : 'non-headline model: this row never gates; its pass column is not a verdict',
  });
  stage.addPickable(glyph, (hit) => rowInfo(shown[hit.instanceId]));
  stage.addPickable(bad, (hit) => rowInfo(shownBad[hit.instanceId]));

  // Grid frame and axis ticks.
  const gcol = 0x2f333a;
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) root.add(fatLine([pos(i, j, 0), pos(i, j, nz - 1)], { color: gcol, width: 1 }));
    for (let k = 0; k < nz; k++) root.add(fatLine([pos(0, j, k), pos(nx - 1, j, k)], { color: gcol, width: 1 }));
  }
  const [x0, , z0] = pos(0, 0, 0);
  A.x.values.forEach((v, i) => labels.add(root, [pos(i, 0, 0)[0], -0.7, z0 - 0.9], `${nice(+v, 4)}`, { kind: 'tick', anchor: 'center' }));
  A.y.values.forEach((v, j) => labels.add(root, [x0 - 0.9, pos(0, j, 0)[1], z0 - 0.9], `${nice(+v, 4)}`, { kind: 'tick', anchor: 'right' }));
  A.z.values.forEach((v, k) => labels.add(root, [x0 - 0.9, -0.7, pos(0, 0, k)[2]], `${nice(+v, 4)}`, { kind: 'tick', anchor: 'right' }));
  labels.add(root, [0, -1.5, z0 - 0.9], A.x.name, { kind: 'axis', anchor: 'center' });
  labels.add(root, [x0 - 1.6, pos(0, ny - 1, 0)[1] + 0.9, z0 - 0.9], A.y.name, { kind: 'axis', anchor: 'right' });
  labels.add(root, [x0 - 1.6, -1.5, 0], A.z.name, { kind: 'axis', anchor: 'right' });

  const ban = typeof o.banner === 'string' ? { text: o.banner, style: 'fail' } : (o.banner || { text: '', style: 'fail' });
  hud.setBanner(ban.text, ban.style === 'pass' ? 'pass' : 'fail', ban.non_gating_passes ? ban.non_gating_passes.text : '');

  // Model filter (only when some rows never gate), then facet chips.
  if (nNonGating) {
    const g = hud.group('rows');
    hud.toggle(g, [[false, 'all models'], [true, 'gating model only']], false, (v) => { gatingOnly = v; refresh(); });
  }
  for (const [f, vals] of Object.entries(o.facets)) {
    const g = hud.group(f === 'Q' ? 'cavity Q' : f);
    const opts = vals.map((v, i) => [i, String(nice(+v, 4) === 'n/a' ? v : nice(+v, 4)).replace(/-design(\.yaml)?$/, '').slice(0, 22)]);
    opts.push([null, 'all']);
    hud.toggle(g, opts, sel[f], (v) => { sel[f] = v; refresh(); });
  }

  const bar = document.createElement('div');
  bar.className = 'v3-colorbar';
  bar.innerHTML = `<div class="v3-colorbar-title">colour: ${o.encoding.colour}</div><div class="v3-colorbar-bar" style="background:${viridisGradientCss()}"><span class="v3-colorbar-mark" style="left:50%"></span></div><div class="v3-colorbar-ticks"><span>0</span><span>0.5</span><span>1</span></div>`
    + `<div class="v3-key" style="margin-top:8px">size: ${o.encoding.size}, ${fmtNum(Math.pow(10, fMin))} to ${fmtNum(Math.pow(10, fMax))} /s</div>`
    + `<div class="v3-key">ring: optical_pass · shell: ${o.encoding.solid_shell} (gating rows only; none drawn if none pass)</div>`
    + (nNonGating ? '<div class="v3-key"><b class="v3-swatch-grey"></b>grey: non-headline model rows, never gate</div>' : '');
  bar.appendChild(countEl);
  hud.legend.appendChild(bar);

  refresh();
  const cy = ((ny - 1) * SP) / 2;
  const D = Math.max(nx, nz, ny) * SP;
  return {
    root,
    defaultView: { pos: [D * 1.1, cy + D * 0.75, D * 1.45], target: [0, cy, 0] },
    bookmarks: [
      { name: 'Overview', pos: [D * 1.1, cy + D * 0.75, D * 1.45], target: [0, cy, 0] },
      { name: `${A.x.name} x ${A.y.name}`, pos: [0, cy, D * 2.1], target: [0, cy, 0] },
      { name: `${A.x.name} x ${A.z.name}`, pos: [0.1, cy + D * 2.0, 0.4], target: [0, cy, 0] },
      { name: `${A.z.name} x ${A.y.name}`, pos: [D * 2.1, cy, 0], target: [0, cy, 0] },
    ],
    noScaleBar: true,
    keepLegend: true,
    testInfo: () => ({ ...info, banner: ban.text, bannerStyle: ban.style, nonGatingText: ban.non_gating_passes?.text || null }),
  };
}
