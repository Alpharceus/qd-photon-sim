// Band-edge profile, 2.5D (band / nitride | inp).  x = growth z at true
// scale (nm), y = energy with its own aspect (badge from the spec), depth
// is an extrusion with no physical meaning.  Nitride: tilted CB/VB from
// nitride_levels.z_profile with |psi|^2 ribbons; InP: piecewise-flat edges
// and levels only (no wavefunctions exist in dot_levels, none invented).
import * as THREE from 'three';
import { fatLine, glowMaterial, breakMark, nice } from './common.js';

const CB_COL = 0x7fa6d6, VB_COL = 0xd99a6a, E_COL = [0.36, 0.62, 1.0], H_COL = [1.0, 0.55, 0.28];

function ribbon(xs, ys, depth, color, opacity = 0.9) {
  const n = xs.length;
  const pos = new Float32Array(n * 2 * 3);
  const idx = [];
  for (let i = 0; i < n; i++) {
    pos.set([xs[i], ys[i], -depth / 2], i * 6);
    pos.set([xs[i], ys[i], depth / 2], i * 6 + 3);
    if (i < n - 1) { const a = 2 * i; idx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2); }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return new THREE.Mesh(g, new THREE.MeshPhysicalMaterial({ color, roughness: 0.3, metalness: 0.2, clearcoat: 0.5, side: THREE.DoubleSide, transparent: opacity < 1, opacity }));
}

// Filled band region between a curve and a flat edge, on the back wall.
function wallFill(xs, ys, yEdge, z, color, opacity) {
  const shape = new THREE.Shape();
  shape.moveTo(xs[0], yEdge);
  for (let i = 0; i < xs.length; i++) shape.lineTo(xs[i], ys[i]);
  shape.lineTo(xs[xs.length - 1], yEdge);
  const m = new THREE.Mesh(new THREE.ShapeGeometry(shape), new THREE.MeshBasicMaterial({ color, transparent: true, opacity, depthWrite: false, side: THREE.DoubleSide }));
  m.position.z = z;
  return m;
}

function arrow(root, x, y0, y1, color, width = 1.6) {
  root.add(fatLine([[x, y0, 0], [x, y1, 0]], { color, width }));
  const s = Math.sign(y1 - y0) || 1;
  const head = new THREE.Mesh(new THREE.ConeGeometry(0.09, 0.26, 16), new THREE.MeshBasicMaterial({ color }));
  head.position.set(x, y1 - s * 0.13, 0);
  if (s < 0) head.rotation.z = Math.PI;
  root.add(head);
}

export function buildBandProfile(spec, ctx) {
  const { stage, labels } = ctx;
  const o = spec.overlays;
  const root = new THREE.Group();
  const ePerUnit = o.energy_per_nm;          // eV per nm of drawing (spec badge states it)
  const nitride = spec.variant === 'nitride';
  // Uniform display factor: the model spans ~18 world units in x.
  let zMin, zMax;
  if (nitride) { zMin = o.z_nm[0]; zMax = o.z_nm[o.z_nm.length - 1]; } else { zMin = o.regions[0].z_nm[0]; zMax = o.regions[o.regions.length - 1].z_nm[1]; }
  const f = 18 / (zMax - zMin);
  const X = (z) => z * f;
  const E0 = nitride ? (o.E_h_eV + o.E_e_eV) / 2 : 0.7;
  const Y = (E) => ((E - E0) / ePerUnit) * f;
  const depth = 5;
  let yTop, yBot;

  if (nitride) {
    const xs = o.z_nm.map(X), cb = o.cb_eV.map(Y), vb = o.vb_eV.map(Y);
    yTop = Math.max(...cb) + 1.6; yBot = Math.min(...vb) - 1.6;
    root.add(ribbon(xs, cb, depth, CB_COL), ribbon(xs, vb, depth, VB_COL));
    root.add(wallFill(xs, cb, yTop, -depth / 2 - 0.01, CB_COL, 0.16), wallFill(xs, vb, yBot, -depth / 2 - 0.01, VB_COL, 0.16));
    const [d0, d1] = o.dot_nm;
    for (const z of [d0, d1]) root.add(fatLine([[X(z), yBot, depth / 2], [X(z), yTop, depth / 2]], { color: 0x8a8f99, width: 1, dashed: true, dashSize: 0.25, gapSize: 0.18 }));
    const dotSlab = new THREE.Mesh(new THREE.PlaneGeometry(X(d1) - X(d0), yTop - yBot), new THREE.MeshBasicMaterial({ color: 0xd9b25c, transparent: true, opacity: 0.05, depthWrite: false }));
    dotSlab.position.set((X(d0) + X(d1)) / 2, (yTop + yBot) / 2, depth / 2 + 0.005);
    root.add(dotSlab);
    // |psi|^2 ribbons at their z-subband energies; amplitude is a drawing choice.
    const amp = 2.4;
    const eY = Y(o.E_e_eV), hY = Y(o.E_h_eV);
    const ePts = xs.map((x, i) => [x, eY + amp * o.psi2_e[i], depth / 2 + 0.02]);
    const hPts = xs.map((x, i) => [x, hY - amp * o.psi2_h[i], depth / 2 + 0.02]);
    const eCol = new THREE.Color(...E_COL).getHex(), hCol = new THREE.Color(...H_COL).getHex();
    const le = fatLine(ePts, { color: eCol, width: 2.4 }), lh = fatLine(hPts, { color: hCol, width: 2.4 });
    stage.markBloom(le, 1.6); stage.markBloom(lh, 1.6);
    root.add(le, lh);
    const fillE = wallFill(xs, ePts.map((p) => p[1]), eY, depth / 2 + 0.01, eCol, 0.22);
    const fillH = wallFill(xs, hPts.map((p) => p[1]), hY, depth / 2 + 0.01, hCol, 0.22);
    fillE.material.blending = fillH.material.blending = THREE.AdditiveBlending;
    root.add(fillE, fillH);
    root.add(fatLine([[xs[0], eY, depth / 2], [xs[xs.length - 1], eY, depth / 2]], { color: eCol, width: 1, opacity: 0.55, dashed: true }));
    root.add(fatLine([[xs[0], hY, depth / 2], [xs[xs.length - 1], hY, depth / 2]], { color: hCol, width: 1, opacity: 0.55, dashed: true }));
    const sc = spec.scalars;
    labels.add(root, [xs[xs.length - 1], cb[cb.length - 1] + 0.3, -depth / 2], 'conduction band edge', { tag: sc.field_kVcm.tag, anchor: 'left', sub: 'In(x)GaN dot / GaN, nitride_levels.z_profile' });
    labels.add(root, [xs[xs.length - 1], vb[vb.length - 1] - 0.3, -depth / 2], 'valence band edge', { tag: sc.field_kVcm.tag, anchor: 'left' });
    const iE = o.psi2_e.indexOf(Math.max(...o.psi2_e)), iH = o.psi2_h.indexOf(Math.max(...o.psi2_h));
    labels.add(root, [xs[iE], eY + amp + 0.2, depth / 2], `|ψe|² at E_e(z) ${nice(o.E_e_eV, 4)} eV`, { tag: sc.overlap_sq.tag, anchor: 'left' });
    labels.add(root, [xs[iH], hY - amp - 0.2, depth / 2], `|ψh|² at E_h(z) ${nice(o.E_h_eV, 4)} eV`, { tag: sc.overlap_sq.tag, anchor: 'right' });
    labels.add(root, [(X(d0) + X(d1)) / 2, yTop + 0.3, depth / 2], `dot ${nice(d1 - d0)} nm · F = ${nice(o.field_kVcm, 4)} kV/cm`, { tag: sc.field_kVcm.tag, anchor: 'center', sub: `tilt ${nice(sc.tilt_eV.value)} eV across the dot · overlap² ${nice(o.overlap_sq)}` });
  } else {
    const regs = o.regions;
    const cbPts = [], vbPts = [];
    for (const r of regs) {
      cbPts.push([X(r.z_nm[0]), Y(r.cb_eV)], [X(r.z_nm[1]), Y(r.cb_eV)]);
      vbPts.push([X(r.z_nm[0]), Y(r.vb_eV)], [X(r.z_nm[1]), Y(r.vb_eV)]);
    }
    const xs = cbPts.map((p) => p[0]);
    const cb = cbPts.map((p) => p[1]), vb = vbPts.map((p) => p[1]);
    yTop = Math.max(...cb) + 1.4; yBot = Math.min(...vb) - 1.4;
    root.add(ribbon(xs, cb, depth, CB_COL), ribbon(xs, vb, depth, VB_COL));
    root.add(wallFill(xs, cb, yTop, -depth / 2 - 0.01, CB_COL, 0.16), wallFill(xs, vb, yBot, -depth / 2 - 0.01, VB_COL, 0.16));
    for (const r of regs.filter((q) => q.open)) {
      const x = r.open === 'left' ? X(r.z_nm[0]) : X(r.z_nm[1]);
      root.add(breakMark([x, yBot, depth / 2], [x, yTop, depth / 2], 0.18, [1, 0, 0], 8));
    }
    const dot = regs.find((r) => r.role === 'active');
    const eY = Y(o.levels.E_e_eV), hY = Y(o.levels.E_h_eV);
    const lvE = fatLine([[X(dot.z_nm[0]), eY, depth / 2 + 0.02], [X(dot.z_nm[1]), eY, depth / 2 + 0.02]], { color: new THREE.Color(...E_COL).getHex(), width: 3 });
    const lvH = fatLine([[X(dot.z_nm[0]), hY, depth / 2 + 0.02], [X(dot.z_nm[1]), hY, depth / 2 + 0.02]], { color: new THREE.Color(...H_COL).getHex(), width: 3 });
    stage.markBloom(lvE, 1.6); stage.markBloom(lvH, 1.6);
    root.add(lvE, lvH);
    const dotSlab = new THREE.Mesh(new THREE.PlaneGeometry(X(dot.z_nm[1]) - X(dot.z_nm[0]), yTop - yBot), glowMaterial([0.85, 0.7, 0.36], { opacity: 0.12 }));
    dotSlab.position.set((X(dot.z_nm[0]) + X(dot.z_nm[1])) / 2, (yTop + yBot) / 2, depth / 2 + 0.005);
    root.add(dotSlab);
    // Escape arrows: matrix arrows at the dot edge, barrier arrows further out.
    const xA = { dE_e_matrix: X(dot.z_nm[1]) + 0.9, dE_e_barrier: X(dot.z_nm[1]) + 3.2, dE_h_matrix: X(dot.z_nm[0]) - 0.9, dE_h_barrier: X(dot.z_nm[0]) - 3.2 };
    for (const a of o.arrows) {
      const x = xA[a.name];
      const col = a.carrier === 'e' ? new THREE.Color(...E_COL).getHex() : new THREE.Color(...H_COL).getHex();
      arrow(root, x, Y(a.from_eV), Y(a.to_eV), col);
      labels.add(root, [x, (Y(a.from_eV) + Y(a.to_eV)) / 2, 0], `${a.name} ${nice(a.value_meV, 4)} meV`, { tag: 'DR', anchor: a.carrier === 'e' ? 'left' : 'right', kind: 'chip' });
    }
    const sc = spec.scalars;
    labels.add(root, [X(dot.z_nm[0]), eY + 0.3, depth / 2], `E_e ${nice(sc.E_e.value, 4)} meV`, { tag: sc.E_e.tag, anchor: 'right' });
    labels.add(root, [X(dot.z_nm[1]), hY - 0.3, depth / 2], `E_h ${nice(sc.E_h.value, 4)} meV`, { tag: sc.E_h.tag, anchor: 'left' });
    for (const r of [regs[1], regs[0]]) {
      labels.add(root, [(X(r.z_nm[0]) + X(r.z_nm[1])) / 2, Y(r.cb_eV) + 0.5, -depth / 2], r.name, { tag: 'DR', anchor: 'center', kind: 'note' });
    }
    labels.add(root, [0, yTop + 0.25, depth / 2], `${dot.name} · ${nice(dot.z_nm[1] - dot.z_nm[0])} nm`, { tag: sc.height_nm.tag, anchor: 'center', sub: 'true thickness on the z axis' });
  }

  // Axes: z ticks (nm) along the front bottom, energy ticks (eV) at left.
  const axisY = yBot - 0.4;
  root.add(fatLine([[X(zMin), axisY, depth / 2], [X(zMax), axisY, depth / 2]], { color: 0x6b6e75, width: 1 }));
  const zStep = nitride ? 5 : 50;
  for (let z = Math.ceil(zMin / zStep) * zStep; z <= zMax + 1e-9; z += zStep) {
    root.add(fatLine([[X(z), axisY, depth / 2], [X(z), axisY - 0.2, depth / 2]], { color: 0x6b6e75, width: 1 }));
    labels.add(root, [X(z), axisY - 0.5, depth / 2], `${z}`, { kind: 'tick', anchor: 'center' });
  }
  labels.add(root, [X(zMax), axisY - 1.05, depth / 2], 'growth z (nm, true scale)', { kind: 'axis', anchor: 'right' });
  const eMin = E0 + (yBot / f) * ePerUnit, eMax = E0 + (yTop / f) * ePerUnit;
  const eStep = (eMax - eMin) > 2.5 ? 1 : 0.5;
  const xAx = X(zMin) - 0.4;
  root.add(fatLine([[xAx, Y(Math.ceil(eMin / eStep) * eStep), -depth / 2], [xAx, Y(Math.floor(eMax / eStep) * eStep), -depth / 2]], { color: 0x6b6e75, width: 1 }));
  for (let e = Math.ceil(eMin / eStep) * eStep; e <= eMax + 1e-9; e += eStep) {
    labels.add(root, [xAx - 0.2, Y(e), -depth / 2], `${+e.toFixed(2)} eV`, { kind: 'tick', anchor: 'right' });
  }

  const cy = (yTop + yBot) / 2;
  const W = 18;
  return {
    root,
    defaultView: { pos: [W * 0.55, cy + W * 0.32, W * 2.3], target: [0.8, cy, 0] },
    bookmarks: [
      { name: 'Front', pos: [0.8, cy, W * 2.5], target: [0.8, cy, 0] },
      { name: 'Angled', pos: [W * 0.9, cy + W * 0.5, W * 1.9], target: [0.8, cy, 0] },
      { name: 'Dot close-up', pos: [W * 0.12, cy + 1, W * 1.0], target: [0, cy, 0] },
      { name: 'Top', pos: [0.2, cy + W * 2.0, W * 0.8], target: [0.8, cy, 0] },
    ],
    noScaleBar: true,
  };
}
