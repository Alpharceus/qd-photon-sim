// Planar DBR microcavity (composed device / structure planar_dbr) and the planar nitride cavity (nitride_planar_cavity).
// An unetched wafer piece drawn as a square coupon, broken laterally.  planar_dbr: pair counts, quarter-wave layer
// thicknesses lambda/4n, the lambda spacer and the dot depth in it all come from the structure response (true scale);
// the 1-D model has no lateral confinement, so the emission plane is a disc, not a mode volume.  nitride_cavity: the
// entry holds Q, V and an efficiency only, so the DBR pairs and indices are display defaults [A]; the mode is the
// equal-volume ellipsoid of V (lambda/n)^3 from the registry.
import * as THREE from 'three';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';
import { layerMaterial, glowMaterial, fatLine, breakMark, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeView } from './micropillar.js';

// Four side walls (open top and bottom) of a w x t x w block centred on x = z = 0, bottom at y0.
function walls(w, t, y0) {
  const h = w / 2, parts = [];
  for (const [rx, ry, px, pz] of [[0, 0, 0, h], [0, Math.PI, 0, -h], [0, Math.PI / 2, h, 0], [0, -Math.PI / 2, -h, 0]]) {
    const p = new THREE.PlaneGeometry(w, t);
    p.rotateX(rx); p.rotateY(ry); p.translate(px, y0 + t / 2, pz);
    parts.push(p);
  }
  return mergeGeometries(parts);
}

export function buildPlanarDbr(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const nitride = spec.variant === 'nitride_cavity';
  const H_um = g.height_um;
  const S = 11 / Math.max(H_um, g.lateral_um * 1.15);          // world units per um
  const nm = (v) => (v / 1000) * S;
  const W = g.lateral_um * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();

  const tH = nm(g.t_high_nm), tL = nm(g.t_low_nm), tS = nm(g.spacer_nm);
  const hiG = [], loG = [], edges = [];
  const edgeAt = (y) => {
    const h = W / 2 * 1.001;
    for (const [a, b] of [[[-h, -h], [h, -h]], [[h, -h], [h, h]], [[h, h], [-h, h]], [[-h, h], [-h, -h]]]) edges.push(a[0], y, a[1], b[0], y, b[1]);
  };
  const add = (y0, t, hi) => { const geo = walls(W, t, y0); (hi ? hiG : loG).push(geo); };
  const yBotTop = nm(g.N_bot * (g.t_high_nm + g.t_low_nm));
  let y = yBotTop;
  for (let i = 0; i < g.N_bot; i++) {                       // bottom DBR, built downward from the spacer: low layer touches it
    y -= tL; add(y, tL, false); edgeAt(y);
    y -= tH; add(y, tH, true); edgeAt(y);
  }
  const ySp1 = yBotTop + tS;
  y = ySp1;
  for (let i = 0; i < g.N_top; i++) {
    add(y, tL, false); y += tL; edgeAt(y);
    add(y, tH, true); y += tH; edgeAt(y);
  }
  const yTop = y;
  const hiMat = layerMaterial('mirror', { clip, opacity: 1 }); hiMat.color.setHex(0xb9c4d2);
  const loMat = layerMaterial('mirror', { clip, opacity: 1 }); loMat.color.setHex(0x566273);
  const hiMesh = new THREE.Mesh(mergeGeometries(hiG), hiMat);
  const loMesh = new THREE.Mesh(mergeGeometries(loG), loMat);
  const dispNote = nitride ? 'pairs and indices are display defaults [A]' : 'quarter wave at the emission wavelength; thickness derived [DR]';
  stage.addPickable(hiMesh, { title: 'high-index DBR layers', tag: nitride ? 'A' : tg('n_high'), rows: [['n', nice(nitride ? 2.45 : g.n_high, 4)], ['thickness λ/4n', `${nice(g.t_high_nm, 4)} nm`], ['pairs', `${g.N_bot} bottom + ${g.N_top} top`]], note: dispNote });
  stage.addPickable(loMesh, { title: 'low-index DBR layers', tag: nitride ? 'A' : tg('n_low'), rows: [['n', nice(nitride ? 2.15 : g.n_low, 4)], ['thickness λ/4n', `${nice(g.t_low_nm, 4)} nm`], ['pairs', `${g.N_bot} bottom + ${g.N_top} top`]], note: dispNote });
  root.add(hiMesh, loMesh);
  const topCap = new THREE.Mesh(new THREE.PlaneGeometry(W, W), hiMat);
  topCap.rotation.x = -Math.PI / 2; topCap.position.y = yTop;
  stage.addPickable(topCap, { title: 'top surface', tag: 'A', rows: [['lateral width drawn', `${nice(g.lateral_um)} µm`]], note: 'the wafer piece is broken laterally; the width is a display default [A]' });
  root.add(topCap);
  const rl = new THREE.LineSegments(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(edges, 3)),
    new THREE.LineBasicMaterial({ color: 0x0b0d10, transparent: true, opacity: 0.35, depthWrite: false }));
  rl.material.clippingPlanes = [clip];
  root.add(rl);

  // spacer (transparent, holds the dot)
  const spM = layerMaterial('core', { clip, opacity: 0.5 });
  spM.color.setHex(0x6f8fae);
  const spacer = new THREE.Mesh(new THREE.BoxGeometry(W, tS, W), spM);
  spacer.position.y = yBotTop + tS / 2;
  stage.addPickable(spacer, { title: 'cavity spacer', tag: 'DR', rows: [['thickness', `${nice(g.spacer_nm, 4)} nm`], ...(nitride ? [] : [['n', nice(g.n_cavity, 4)]])], note: nitride ? 'cavity order is a display default [A]' : 'λ / n_cavity (a one-wavelength cavity)' });
  root.add(spacer);

  // substrate: a broken box under the bottom DBR
  const subT = nm(g.substrate_um * 1000);
  const subMat = new THREE.MeshStandardMaterial({ color: 0x2a2f37, roughness: 0.95, side: THREE.DoubleSide, clippingPlanes: [clip] });
  const sub = new THREE.Mesh(new THREE.BoxGeometry(W, subT, W), subMat);
  sub.position.y = -subT / 2 - 0.004;
  stage.addPickable(sub, { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']], note: 'drawn as a broken plate under the bottom DBR' });
  root.add(sub);
  const h = W / 2;
  root.add(breakMark([-h, -subT / 2, h + 0.003], [h, -subT / 2, h + 0.003], 0.05, [0, 1, 0], 9));
  root.add(breakMark([h + 0.003, -subT / 2, -h], [h + 0.003, -subT / 2, h], 0.05, [0, 1, 0], 9));

  // dot, emission plane / mode
  const yd = yBotTop + nm(g.z_dot_nm);
  if (nitride) {
    const [a, , c] = g.mode_semiaxes_um;
    const m = new THREE.Mesh(new THREE.SphereGeometry(1, 40, 24), glowMaterial(rgb, { opacity: 0.2 }));
    m.scale.set(a * S, c * S, a * S);
    m.position.set(0, yd, 0);
    stage.markBloom(m);
    root.add(m);
  } else {
    const m = new THREE.Mesh(new THREE.CircleGeometry(W * 0.36, 64), glowMaterial(rgb, { opacity: 0.16 }));
    m.rotation.x = -Math.PI / 2; m.position.y = yd;
    stage.markBloom(m);
    root.add(m);
    root.add(fatLine(Array.from({ length: 65 }, (_, i) => [Math.cos((i / 64) * Math.PI * 2) * W * 0.36, yd, Math.sin((i / 64) * Math.PI * 2) * W * 0.36]), { color: 0xd8d6cf, width: 1.1, opacity: 0.55 }));
  }
  composeDot(stage, root, rgb, [0, yd, 0], Math.max(W * 0.012, 0.05),
    { title: 'quantum dot', tag: nitride ? 'A' : tg('z_dot_frac'), rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['height in the spacer', `${nice(g.z_dot_nm, 4)} nm`]], note: nitride ? 'dot depth is a display default [A]' : 'dot depth z_dot_frac x spacer from the spec' });

  // emission: up through the top DBR, and (planar_dbr) the leak down through the substrate, lengths ∝ the split
  const eta = g.eta_out;
  composeBeam(stage, root, [0, yTop + 0.02, 0], [0, 1, 0], 5, Math.max(W * 0.1, 0.2));
  if (!nitride && Number.isFinite(eta) && eta < 0.97) {
    const dn = fatLine([[0, -subT, 0], [0, -subT - Math.max(5 * (1 - eta), 0.5), 0]], { color: 0xff3b2f, width: 2.5, opacity: 0.55 });
    root.add(dn);
  }

  // annotations
  const lx = W / 2 + 0.3;
  labels.add(root, [lx, yTop - nm(g.N_top * (g.t_high_nm + g.t_low_nm)) / 2, W / 2], `top DBR · ${g.N_top} pairs · ${nice(g.N_top * (g.t_high_nm + g.t_low_nm) / 1000, 3)} µm`, { tag: nitride ? 'A' : tg('N_top'), priority: 70 });
  labels.add(root, [lx, yBotTop / 2, W / 2], `bottom DBR · ${g.N_bot} pairs · ${nice(g.N_bot * (g.t_high_nm + g.t_low_nm) / 1000, 3)} µm`, { tag: nitride ? 'A' : tg('N_bot'), priority: 70 });
  labels.add(root, [lx, yBotTop + tS / 2, W / 2], `spacer · ${nice(g.spacer_nm, 4)} nm`, { tag: 'DR', priority: 85, sub: nitride ? 'λ/n_high; cavity order is a display default' : 'λ/n_cavity' });
  labels.add(root, [-lx, yd, W * 0.3], `dot · ${nice(lam, 4)} nm`, { tag: nitride ? 'A' : tg('z_dot_frac'), anchor: 'right', priority: 90, sub: falseColour ? colourNote : `${nice(g.z_dot_nm, 3)} nm above the spacer bottom` });
  labels.add(root, [0, yTop + 0.5, W / 2], `${nitride ? 'nitride cavity' : 'planar cavity'} · ${nice(g.height_um, 3)} µm tall`, { tag: 'DR', anchor: 'center', priority: 80, sub: `${nice(g.lateral_um)} µm width is a display default [A]` });
  labels.add(root, [-lx, yTop - nm(g.t_high_nm) * 3, W / 2], `λ/4n · ${nice(g.t_high_nm, 3)} / ${nice(g.t_low_nm, 3)} nm`, { tag: nitride ? 'A' : 'DR', anchor: 'right', priority: 60, sub: nitride ? 'display indices' : `n ${nice(g.n_high, 3)} / ${nice(g.n_low, 3)}` });
  if (nitride) {
    const [a, , c] = g.mode_semiaxes_um;
    labels.add(root, [W * 0.18, yd + nm(g.spacer_nm) * 0.4, W * 0.18], `mode V = ${nice(g.mode_volume_norm, 3)} (λ/n)³`, { tag: tg('mode_volume_norm'), priority: 75, sub: `ellipsoid ${nice(2000 * a, 3)} × ${nice(2000 * c, 3)} nm; Q ${nice(g.Q, 3)}` });
  } else {
    labels.add(root, [W * 0.3, yd, W * 0.3], 'emission plane', { tag: 'DR', priority: 65, sub: '1-D model: no lateral confinement, no mode volume' });
    if (Number.isFinite(eta)) labels.add(root, [0, -subT - 1.3, 0], `top exit ${nice(eta, 3)} · substrate ${nice(1 - eta, 3)}`, { tag: 'DR', anchor: 'center', priority: 50, sub: 'T_top / (T_top + T_bottom)' });
  }
  labels.add(root, [-lx, -subT / 2, W * 0.4], 'substrate', { tag: 'A', anchor: 'right', priority: 25, sub: 'thickness not modelled' });

  const tc = [0, yTop * 0.5, 0];
  const main = composeView(tc, [0.62, 0.22, 0.86], yTop * 2.4 + W);
  const bm = (name, vw, cut) => ({ name, pos: vw.pos, target: vw.target, cut, mag: false });
  return {
    root,
    defaultView: main,
    bookmarks: [
      bm('Wafer piece', main, W * 0.6),
      bm('Side cut-away', composeView([0, yd, 0], [1, 0.04, 0.05], Math.max(yTop, W) * 1.9), 0),
      bm('Spacer close-up', composeView([0, yd, 0], [0.5, 0.12, 0.85], Math.max(tS * 9, 2.6)), 0),
      bm('Top / beam', composeView([0, yTop, 0], [0.45, 0.55, 0.7], W * 2.2), W * 0.6),
    ],
    cut: { min: -W / 2, max: W * 0.6 },
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
