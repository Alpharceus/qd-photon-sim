// SVG fallback drawn from the SAME SceneSpec when WebGL2 is unavailable or
// the context is lost: true-scale cross-section (device), band diagram
// (band), population traces (cascade), lo/hi heat maps (surface) and
// small-multiple scatter (lattice).  Readouts keep their tag chips.
import { viridisCss, lambdaToRGB, ROLE_TINT } from './core/colormap.js';
import { tagChip, esc } from './core/labels.js';
import { fmtNum } from './core/hud.js';

const NS = 'http://www.w3.org/2000/svg';
const hex = (n) => `#${n.toString(16).padStart(6, '0')}`;

function svg(w, h) {
  const s = document.createElementNS(NS, 'svg');
  s.setAttribute('viewBox', `0 0 ${w} ${h}`);
  s.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  return s;
}
function add(parent, tag, attrs = {}, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  parent.appendChild(e);
  return e;
}
function polyline(parent, pts, attrs) {
  add(parent, 'polyline', { points: pts.map((p) => p.map((v) => v.toFixed(2)).join(',')).join(' '), fill: 'none', ...attrs });
}

function device(s, spec) {
  const W = 900, H = 520;
  const layers = spec.stacks?.optical?.layers || spec.layers || [];
  const total = layers.reduce((a, L) => a + (L.t_nm || 0), 0) || 1;
  const k = 380 / total;                      // px per nm, uniform (true scale)
  let y = 450;
  const x0 = 160, w = 420;
  for (const L of layers) {
    const h = (L.t_nm || 0) * k;
    y -= h;
    add(s, 'rect', { x: x0, y, width: w, height: Math.max(h, 1), fill: hex(ROLE_TINT[L.role] ?? 0x888888), 'fill-opacity': L.role === 'cladding' ? 0.5 : 0.9 });
    add(s, 'text', { x: x0 + w + 12, y: y + Math.max(h, 1) / 2 + 4, 'font-size': 12 }, `${L.name} · ${fmtNum(L.t_nm)} nm [${L.tag}]`);
  }
  const mode = spec.overlays?.mode;
  if (mode) {
    const raw = atob(mode.data);
    const c = document.createElement('canvas');
    c.width = mode.nx; c.height = mode.ny;
    const ctx = c.getContext('2d');
    const img = ctx.createImageData(mode.nx, mode.ny);
    const { rgb } = lambdaToRGB(spec.overlays.emission.lambda_nm);
    for (let i = 0; i < mode.nx * mode.ny; i++) {
      const v = raw.charCodeAt(i) / 255;
      const row = mode.ny - 1 - Math.floor(i / mode.nx), col = i % mode.nx;
      img.data.set([rgb[0] * 255, rgb[1] * 255, rgb[2] * 255, v * 255], (row * mode.nx + col) * 4);
    }
    ctx.putImageData(img, 0, 0);
    const mw = (mode.x_um[1] - mode.x_um[0]) * 1000 * k, mh = (mode.y_um[1] - mode.y_um[0]) * 1000 * k;
    add(s, 'image', { href: c.toDataURL(), x: x0 + w / 2 - mw / 2, y: 450 - mh, width: mw, height: mh, preserveAspectRatio: 'none' });
  }
  add(s, 'text', { x: x0, y: 480, 'font-size': 12 }, `cross-section at true scale (${fmtNum(total / 1000)} µm stack; lateral not to scale)`);
}

function band(s, spec) {
  const W = 900, H = 520, o = spec.overlays;
  if (spec.variant === 'nitride') {
    const zs = o.z_nm, all = [...o.cb_eV, ...o.vb_eV];
    const eMin = Math.min(...all) - 0.3, eMax = Math.max(...all) + 0.3;
    const X = (z) => 80 + ((z - zs[0]) / (zs[zs.length - 1] - zs[0])) * (W - 160);
    const Y = (e) => 40 + (1 - (e - eMin) / (eMax - eMin)) * (H - 100);
    polyline(s, zs.map((z, i) => [X(z), Y(o.cb_eV[i])]), { stroke: '#7fa6d6', 'stroke-width': 2 });
    polyline(s, zs.map((z, i) => [X(z), Y(o.vb_eV[i])]), { stroke: '#d99a6a', 'stroke-width': 2 });
    polyline(s, zs.map((z, i) => [X(z), Y(o.E_e_eV + 0.5 * o.psi2_e[i])]), { stroke: '#5c9eff', 'stroke-width': 1.5 });
    polyline(s, zs.map((z, i) => [X(z), Y(o.E_h_eV - 0.5 * o.psi2_h[i])]), { stroke: '#ff8c47', 'stroke-width': 1.5 });
    add(s, 'text', { x: 80, y: H - 20, 'font-size': 12 }, `z (nm), F = ${fmtNum(o.field_kVcm)} kV/cm, overlap² = ${fmtNum(o.overlap_sq)}`);
  } else {
    const regs = o.regions;
    const eMin = Math.min(...regs.map((r) => r.vb_eV)) - 0.2, eMax = Math.max(...regs.map((r) => r.cb_eV)) + 0.2;
    const z0 = regs[0].z_nm[0], z1 = regs[regs.length - 1].z_nm[1];
    const X = (z) => 80 + ((z - z0) / (z1 - z0)) * (W - 160);
    const Y = (e) => 40 + (1 - (e - eMin) / (eMax - eMin)) * (H - 100);
    for (const r of regs) {
      add(s, 'line', { x1: X(r.z_nm[0]), x2: X(r.z_nm[1]), y1: Y(r.cb_eV), y2: Y(r.cb_eV), stroke: '#7fa6d6', 'stroke-width': 2 });
      add(s, 'line', { x1: X(r.z_nm[0]), x2: X(r.z_nm[1]), y1: Y(r.vb_eV), y2: Y(r.vb_eV), stroke: '#d99a6a', 'stroke-width': 2 });
    }
    add(s, 'line', { x1: X(-10), x2: X(10), y1: Y(o.levels.E_e_eV), y2: Y(o.levels.E_e_eV), stroke: '#5c9eff', 'stroke-width': 2 });
    add(s, 'line', { x1: X(-10), x2: X(10), y1: Y(o.levels.E_h_eV), y2: Y(o.levels.E_h_eV), stroke: '#ff8c47', 'stroke-width': 2 });
    add(s, 'text', { x: 80, y: H - 20, 'font-size': 12 }, 'z (nm, true scale), energies relative to the dot VB top; levels only');
  }
}

function cascade(s, spec) {
  const W = 900, H = 520, o = spec.overlays;
  const X = (i) => 80 + (i / (o.t_ns.length - 1)) * (W - 160);
  const Y = (P) => 40 + (1 - (Math.log10(Math.max(P, 1e-12)) + 12) / 12) * (H - 100);
  for (const [k, c] of [['P_G', '#c0c4cc'], ['P_X', '#ff5a4a'], ['P_XX', '#9e85fa']]) polyline(s, o[k].map((p, i) => [X(i), Y(p)]), { stroke: c, 'stroke-width': 1.8 });
  add(s, 'text', { x: 80, y: H - 20, 'font-size': 12 }, `populations P(t), log scale, sample index axis (pump 0-${fmtNum(o.tau_on_ns)} ns, then dark to ${fmtNum(o.period_ns)} ns)`);
}

function heat(s, grid, x0, y0, w, h, title) {
  const ny = grid.length, nx = grid[0].length;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = grid[j][i];
    add(s, 'rect', { x: x0 + (i * w) / nx, y: y0 + h - ((j + 1) * h) / ny, width: w / nx + 0.5, height: h / ny + 0.5, fill: v === null ? '#2a2c32' : viridisCss(v) });
  }
  add(s, 'text', { x: x0, y: y0 - 8, 'font-size': 13 }, title);
}

function surface(s, spec) {
  const o = spec.overlays;
  heat(s, o.lo, 60, 60, 370, 380, 'lo sheet g2(T, param)');
  heat(s, o.hi, 470, 60, 370, 380, 'hi sheet g2(T, param)');
  add(s, 'text', { x: 60, y: 480, 'font-size': 12 }, `x: ${spec.axes.x.label} ${o.x[0]}-${o.x[o.x.length - 1]} K; y: ${spec.axes.y.label}; ${o.sheet_note}`);
}

function lattice(s, spec) {
  const o = spec.overlays, R = o.rows, A = o.axes;
  const nz = A.z.values.length, nx = A.x.values.length, ny = A.y.values.length;
  const pw = 820 / nz;
  for (let k = 0; k < nz; k++) {
    const x0 = 40 + k * pw;
    add(s, 'text', { x: x0, y: 40, 'font-size': 12 }, `${A.z.name} = ${A.z.values[k]}`);
    for (let r = 0; r < R.g2.length; r++) {
      if (R.zi[r] !== k) continue;
      let ok = true;
      for (const [f, v] of Object.entries(o.facet_default)) if (v !== null && R.facet[f][r] !== v) ok = false;
      if (!ok) continue;
      const cx = x0 + 20 + (R.xi[r] / Math.max(1, nx - 1)) * (pw - 60), cy = 420 - (R.yi[r] / Math.max(1, ny - 1)) * 340;
      const gating = !R.gates || R.gates[r] !== false;   // non-gating rows: grey, no pass ring
      add(s, 'circle', { cx, cy, r: 8, fill: !gating ? '#5d6067' : Number.isFinite(R.g2[r]) ? viridisCss(R.g2[r]) : '#4a4d55', stroke: gating && R.optical_pass[r] ? '#f4f1e6' : 'none', 'stroke-width': 2 });
    }
  }
  add(s, 'text', { x: 40, y: 480, 'font-size': 12 }, `x: ${A.x.name}, y: ${A.y.name}; ${typeof o.banner === 'string' ? o.banner : (o.banner?.text || '')}`);
}

export function renderFallback(el, kind, spec, { reason = '' } = {}) {
  const box = document.createElement('div');
  box.className = 'v3-fallback';
  const head = document.createElement('div');
  head.innerHTML = `<div class="v3-title">${esc(spec.title || kind)}</div><div class="v3-note">2D fallback (${esc(reason)}), drawn from the same SceneSpec</div>`;
  box.appendChild(head);
  const s = svg(900, 520);
  ({ device, band, cascade, surface, lattice })[kind]?.(s, spec);
  box.appendChild(s);
  const ro = document.createElement('div');
  ro.className = 'v3-note';
  for (const [k, v] of Object.entries(spec.scalars || {})) {
    if (!v?.tag) continue;
    const span = document.createElement('span');
    span.style.marginRight = '14px';
    span.textContent = `${v.label}: ${fmtNum(v.value)} ${v.unit || ''} `;
    span.appendChild(tagChip(v.tag));
    ro.appendChild(span);
    void k;
  }
  box.appendChild(ro);
  for (const l of spec.labels || []) {
    const n = document.createElement('div');
    n.className = 'v3-note';
    n.textContent = l;
    box.appendChild(n);
  }
  el.appendChild(box);
  return box;
}
