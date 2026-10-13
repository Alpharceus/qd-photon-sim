// SVG fallback drawn from the SAME SceneSpec when WebGL2 is unavailable or
// the context is lost: true-scale cross-section (device), band diagram
// (band), population traces (cascade), lo/hi heat maps (surface) and
// small-multiple scatter (lattice).  Readouts keep their tag chips.
import { viridisCss, lambdaToRGB, ROLE_TINT } from './core/colormap.js';
import { tagChip, esc } from './core/labels.js';
import { fmtNum } from './core/hud.js';
import { richText } from '../ui/math.js';
import { sceneTitle } from '../ui/labels.js';

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

// ---- R2-3D: the composed structure kinds with their own 2D cross-section (WebGL off / forced fallback)
const rgbCss = (lam) => { const { rgb } = lambdaToRGB(lam || 600); return `rgb(${rgb.map((v) => Math.round(v * 255)).join(',')})`; };

function nanowireSvg(s, spec) {
  const g = spec.geometry, H = g.length_um + g.mirror_um, k = 380 / H, cx = 360, y0 = 450;
  const r0 = (g.d_nm / 2000) * k, rT = (g.top_d_nm / 2000) * k, yM = y0 - g.mirror_um * k, yTop = y0 - H * k, yF = yTop + g.flare_len_um * k;
  const mw = g.mirror_halfwidth_um * k;
  add(s, 'rect', { x: cx - mw, y: yM, width: 2 * mw, height: g.mirror_um * k, fill: '#a8832f' });
  add(s, 'polygon', { points: `${cx - r0},${yM} ${cx - r0},${yF} ${cx - rT},${yTop} ${cx + rT},${yTop} ${cx + r0},${yF} ${cx + r0},${yM}`, fill: '#86a3c2', stroke: '#d8d6cf', 'stroke-width': 1 });
  const yd = yM - g.dot_height_um * k;
  add(s, 'circle', { cx, cy: yd, r: 5, fill: rgbCss(g.lambda_nm) });
  add(s, 'line', { x1: cx, x2: cx, y1: yTop - 4, y2: yTop - 40, stroke: '#ff3b2f', 'stroke-width': 3 });
  add(s, 'text', { x: cx + mw + 14, y: yM + 12, 'font-size': 12 }, `bottom mirror R ${fmtNum(g.mirror_R)} [${spec.geometry_src?.mirror_R?.tag || 'E'}]`);
  add(s, 'text', { x: cx + rT + 14, y: (yM + yF) / 2, 'font-size': 12 }, `wire Ø ${fmtNum(g.d_nm)} nm [${spec.geometry_src?.d_nm?.tag || 'A'}], ${fmtNum(g.length_um)} µm tall [A]`);
  add(s, 'text', { x: cx + rT + 14, y: yTop + 12, 'font-size': 12 }, g.trumpet ? `trumpet top Ø ${fmtNum(g.top_d_nm)} nm` : 'needle (no taper)');
  add(s, 'text', { x: cx - r0 - 14, y: yd + 4, 'font-size': 12, 'text-anchor': 'end' }, `dot ${fmtNum(g.lambda_nm)} nm [A]`);
  add(s, 'text', { x: 40, y: 500, 'font-size': 12 }, `cross-section at one true scale (${fmtNum(H)} µm; wire height and mirror size are display defaults [A])`);
}

function planarSvg(s, spec, k, x0, w, layers) {
  const g = spec.geometry;
  let off = 0;
  for (const L of layers) { if (L.name === 'spacer') break; off += L.t_nm || 0; }
  const yd = 450 - (off + g.z_dot_nm) * k;
  if (g.mode === 'ellipsoid') {
    const [a, , c] = g.mode_semiaxes_um;
    add(s, 'ellipse', { cx: x0 + w / 2, cy: yd, rx: Math.max(a * 1000 * k, 2), ry: Math.max(c * 1000 * k, 2), fill: rgbCss(g.lambda_nm), 'fill-opacity': 0.35 });
  } else {
    add(s, 'line', { x1: x0, x2: x0 + w, y1: yd, y2: yd, stroke: rgbCss(g.lambda_nm), 'stroke-dasharray': '6 4', 'stroke-width': 1.5 });
  }
  add(s, 'circle', { cx: x0 + w / 2, cy: yd, r: 5, fill: rgbCss(g.lambda_nm) });
  add(s, 'text', { x: x0 - 8, y: yd + 4, 'font-size': 12, 'text-anchor': 'end' }, 'dot');
}

function ridgeSvg(s, spec) {
  const g = spec.geometry, k = 110, cx = 330, y0 = 440;          // px per µm, cross-section at the front facet
  const nm = (v) => (v / 1000) * k, W = nm(g.ridge_width_nm), pad = nm(g.shoulder_nm), floorY = nm(g.etch_floor_nm);
  let y = 0;
  for (const L of g.stack) {
    const t = nm(L.t_nm), a = y, b = y + t;
    y = b;
    const fill = L.is_dot ? '#d9a845' : (L.name.includes('cladding') ? '#8794a6' : '#34587f');
    add(s, 'rect', { x: cx - W / 2, y: y0 - b, width: W, height: Math.max(t, 1), fill });
    if (a < floorY) {
      const top = Math.min(b, floorY);
      for (const sg of [-1, 1]) add(s, 'rect', { x: sg < 0 ? cx - W / 2 - pad : cx + W / 2, y: y0 - top, width: pad, height: Math.max(top - a, 1), fill });
    }
  }
  add(s, 'rect', { x: cx - W / 2 - pad, y: y0, width: W + 2 * pad, height: nm(g.substrate_um * 1000), fill: '#1c1f25' });
  const yDot = y0 - nm(g.dot_height_nm);
  add(s, 'ellipse', { cx, cy: yDot, rx: g.wx_um * k, ry: g.wy_um * k, fill: rgbCss(g.lambda_nm), 'fill-opacity': 0.35, stroke: rgbCss(g.lambda_nm) });
  add(s, 'line', { x1: cx - W / 2, x2: cx + W / 2, y1: yDot, y2: yDot, stroke: '#ffc457', 'stroke-width': 2 });
  add(s, 'text', { x: cx + W / 2 + pad + 14, y: yDot + 4, 'font-size': 12 }, `dot layer · mode ${fmtNum(2 * g.wx_um)} × ${fmtNum(2 * g.wy_um)} µm (1/e² intensity)`);
  add(s, 'text', { x: cx + W / 2 + pad + 14, y: y0 - nm(g.total_nm) / 2, 'font-size': 12 }, `ridge ${fmtNum(g.ridge_width_nm)} nm wide, etched ${fmtNum(g.etch_depth_nm)} nm`);
  add(s, 'text', { x: 40, y: 500, 'font-size': 12 }, `facet cross-section at true scale (${fmtNum(g.total_nm / 1000)} µm stack); ridge length ${fmtNum(g.L_um)} µm not drawn`);
}

function lumpedSvg(s, spec) {
  const g = spec.geometry, cx = 360;
  for (const yy of [150, 330]) add(s, 'rect', { x: cx - 140, y: yy, width: 280, height: 12, fill: '#a9b6c6', 'fill-opacity': 0.7 });
  add(s, 'rect', { x: cx - 110, y: 162, width: 220, height: 168, fill: 'none', stroke: '#d8d6cf', 'stroke-dasharray': '6 4', 'stroke-opacity': 0.5 });
  add(s, 'ellipse', { cx, cy: 246, rx: 70, ry: 60, fill: rgbCss(g.lambda_nm), 'fill-opacity': 0.3 });
  add(s, 'circle', { cx, cy: 246, r: 6, fill: rgbCss(g.lambda_nm) });
  add(s, 'line', { x1: cx, x2: cx, y1: 148, y2: 90, stroke: '#ff3b2f', 'stroke-width': 3 });
  add(s, 'text', { x: cx + 160, y: 156, 'font-size': 12 }, `output mirror · κ = ${g.kappa_ueV == null ? 'n/a' : fmtNum(g.kappa_ueV) + ' µeV'}`);
  add(s, 'text', { x: cx + 160, y: 246, 'font-size': 12 }, `F_P ${fmtNum(g.F_P)} · η_out ${fmtNum(g.eta_out)} · β0 ${fmtNum(g.beta0)}`);
  add(s, 'text', { x: 40, y: 440, 'font-size': 14, 'font-weight': 'bold' }, 'NO GEOMETRY IN THIS MODEL (lumped cavity)');
  add(s, 'text', { x: 40, y: 462, 'font-size': 12 }, 'schematic, not to scale: only F_P, κ, detuning, β and η_out are model quantities');
}

const OWN = { photonic_nanowire: nanowireSvg, ridge_waveguide: ridgeSvg, lumped_cavity: lumpedSvg };

function device(s, spec) {
  if (OWN[spec.variant] && spec.geometry) { OWN[spec.variant](s, spec); return; }
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
  if ((spec.variant === 'planar_dbr' || spec.variant === 'nitride_cavity') && spec.geometry) planarSvg(s, spec, k, x0, w, layers);
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
  head.innerHTML = `<div class="v3-title"></div><div class="v3-note">2D fallback (${esc(reason)}), drawn from the same SceneSpec</div>`;
  head.querySelector('.v3-title').replaceChildren(richText(sceneTitle({ ...spec, title: spec.title || kind })));
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
    n.replaceChildren(richText(l));
    box.appendChild(n);
  }
  el.appendChild(box);
  return box;
}
