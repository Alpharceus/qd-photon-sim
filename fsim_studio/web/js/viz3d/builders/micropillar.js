// DBR micropillar (composed device / structure micropillar).  True scale in
// every axis: the pillar diameter, the pair counts and the quarter-wave layer
// thicknesses lambda/4n all come from the descriptor (scene_compose.py), which
// takes them from the resolved spec.  Display defaults (dot depth, substrate
// thickness) carry an [A] tag.  Two merged meshes hold the high- and low-index
// layers (2 draw calls for ~130 layers).
//
// This file also exports the small helpers the other composed-device builders
// share (compose* below): tag lookup, emission dot, breadboard beam, frame.
import * as THREE from 'three';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';
import { layerMaterial, glowMaterial, fatLine, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';

export const BEAM_RGB = [1.0, 0.231, 0.184];       // the one reserved ink

// Provenance tag of a descriptor geometry key.
export function composeTag(spec, key, dflt = 'A') {
  return spec.geometry_src?.[key]?.tag || dflt;
}

// Emission dot: a hot core with a soft halo.  Not clipped, so it stays
// visible in the cut-away.
export function composeDot(stage, parent, rgb, pos, r, info) {
  const core = new THREE.Mesh(new THREE.SphereGeometry(r, 24, 16),
    new THREE.MeshStandardMaterial({ color: new THREE.Color(...rgb), emissive: new THREE.Color(...rgb), emissiveIntensity: 2.4, roughness: 0.4 }));
  core.position.set(...pos);
  stage.markBloom(core);
  if (info) stage.addPickable(core, info);
  const halo = new THREE.Mesh(new THREE.SphereGeometry(r * 3.2, 24, 16), glowMaterial(rgb, { opacity: 0.22 }));
  halo.position.set(...pos);
  stage.markBloom(halo, 2.0);
  parent.add(core, halo);
  stage.dotAnchor = { pos: new THREE.Vector3(...pos), r };      // the drive layer (drive.js) aims at the dot
  return core;
}

// The breadboard beam leaving the device along `dir` from `from`, with the
// emission glow at the exit point.  Beam ink, not the emission wavelength.
export function composeBeam(stage, parent, from, dir, len, glowR) {
  const d = new THREE.Vector3(...dir).normalize();
  const a = new THREE.Vector3(...from), b = a.clone().addScaledVector(d, len);
  const line = fatLine([[a.x, a.y, a.z], [b.x, b.y, b.z]], { color: 0xff3b2f, width: 4 });
  line.name = 'breadboard-beam';
  const glow = new THREE.Mesh(new THREE.SphereGeometry(glowR, 24, 16), glowMaterial(BEAM_RGB, { opacity: 0.6 }));
  glow.scale.set(1, 0.5, 1);
  glow.position.copy(a).addScaledVector(d, glowR * 0.25);
  glow.name = 'emission-glow';
  stage.markBloom(glow, 2.5);
  parent.add(line, glow);
}

// Zig-zag break mark running round a cylinder of radius r at height y (hugs the
// surface, so it reads from every angle).
export function composeRingBreak(r, y, amp, teeth = 30) {
  const pts = [];
  const n = teeth * 2;
  for (let i = 0; i <= n; i++) {
    const a = (i / n) * Math.PI * 2;
    pts.push([Math.cos(a) * r * 1.002, y + (i % 2 ? amp : -amp), Math.sin(a) * r * 1.002]);
  }
  return fatLine(pts, { color: 0xd8d6cf, width: 1.4 });
}

// Substrate plate: wall + top + a bottom cap lifted 4 mm-units off the plate's
// lowest plane (the enclosure floor's post top is coplanar with it and would
// z-fight in the cut-away).  top = y of the upper face.
export function composeSubstrate(stage, radius, thick, top, { color = 0x2a2f37, info = null } = {}) {
  const mat = new THREE.MeshStandardMaterial({ color, roughness: 0.95, metalness: 0, side: THREE.DoubleSide, clippingPlanes: [stage.clipPlane] });
  const grp = new THREE.Group();
  const wall = new THREE.Mesh(new THREE.CylinderGeometry(radius, radius, thick, 96, 1, true), mat);
  wall.position.y = top - thick / 2;
  const up = new THREE.Mesh(new THREE.CircleGeometry(radius, 96), mat);
  up.rotation.x = -Math.PI / 2; up.position.y = top - 0.004;
  const dn = new THREE.Mesh(new THREE.CircleGeometry(radius, 96), mat);
  dn.rotation.x = Math.PI / 2; dn.position.y = top - thick + 0.004;
  grp.add(wall, up, dn);
  if (info) stage.addPickable(wall, info);
  return grp;
}

// Camera position at distance d from target t along the direction dir.
export function composeView(t, dir, d) {
  const n = Math.hypot(...dir);
  return { pos: [0, 1, 2].map((i) => t[i] + (dir[i] / n) * d), target: t };
}

// Plain solid-colour lathe of a (r, y) profile.
export function composeLathe(profile, material, segs = 64) {
  return new THREE.Mesh(new THREE.LatheGeometry(profile.map(([r, y]) => new THREE.Vector2(r, y)), segs), material);
}

export function buildMicropillar(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const H_um = g.height_um;
  const S = 12 / Math.max(H_um, 2.5 * g.d_um);          // world units per um
  const nm = (v) => (v / 1000) * S;
  const R = (g.d_um / 2) * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();

  // --- layer stack, bottom -> top: bottom DBR, spacer, top DBR.  The layer next
  // to the spacer is the low-index one on both sides (alternating outward).
  const tH = nm(g.t_high_nm), tL = nm(g.t_low_nm), tS = nm(g.spacer_nm);
  const hiGeo = [], loGeo = [];
  const rings = [];
  const seg = 72;
  const addLayer = (y0, t, hi) => {
    const geo = new THREE.CylinderGeometry(R, R, t, seg, 1, true);     // walls only: coplanar caps would z-fight in the cut-away
    geo.translate(0, y0 + t / 2, 0);
    (hi ? hiGeo : loGeo).push(geo);
  };
  const ringAt = (y) => {
    for (let i = 0; i < seg; i++) {
      const a0 = (i / seg) * Math.PI * 2, a1 = ((i + 1) / seg) * Math.PI * 2;
      rings.push(Math.cos(a0) * R * 1.001, y, Math.sin(a0) * R * 1.001, Math.cos(a1) * R * 1.001, y, Math.sin(a1) * R * 1.001);
    }
  };
  let y = 0;
  // bottom DBR: build from the spacer downward so the low layer touches the spacer
  const yBotTop = nm(g.N_bot * (g.t_high_nm + g.t_low_nm));
  y = yBotTop;
  for (let i = 0; i < g.N_bot; i++) {           // low then high going down
    y -= tL; addLayer(y, tL, false); ringAt(y);
    y -= tH; addLayer(y, tH, true); ringAt(y);
  }
  const ySpacer0 = yBotTop, ySpacer1 = yBotTop + tS;
  y = ySpacer1;
  for (let i = 0; i < g.N_top; i++) {           // low then high going up
    addLayer(y, tL, false); y += tL; ringAt(y);
    addLayer(y, tH, true); y += tH; ringAt(y);
  }
  const yTop = y;
  const hiMat = layerMaterial('mirror', { clip, opacity: 1 }); hiMat.color.setHex(0xb9c4d2);
  const loMat = layerMaterial('mirror', { clip, opacity: 1 }); loMat.color.setHex(0x566273);
  const hiMesh = new THREE.Mesh(mergeGeometries(hiGeo), hiMat);
  const loMesh = new THREE.Mesh(mergeGeometries(loGeo), loMat);
  stage.addPickable(hiMesh, { title: 'high-index DBR layers', tag: tg('n_high'), rows: [['n', nice(g.n_high, 4)], ['thickness λ/4n', `${nice(g.t_high_nm, 4)} nm`], ['pairs', `${g.N_bot} bottom + ${g.N_top} top`]], note: 'quarter wave at the emission wavelength; thickness derived [DR]' });
  stage.addPickable(loMesh, { title: 'low-index DBR layers', tag: tg('n_low'), rows: [['n', nice(g.n_low, 4)], ['thickness λ/4n', `${nice(g.t_low_nm, 4)} nm`], ['pairs', `${g.N_bot} bottom + ${g.N_top} top`]], note: 'quarter wave at the emission wavelength; thickness derived [DR]' });
  root.add(hiMesh, loMesh);
  const topCap = new THREE.Mesh(new THREE.CircleGeometry(R, seg), hiMat);
  topCap.rotation.x = -Math.PI / 2; topCap.position.y = yTop;
  stage.addPickable(topCap, { title: 'top surface', tag: tg('d_um'), rows: [['pillar diameter', `${nice(g.d_um)} µm`]] });
  root.add(topCap);
  const rl = new THREE.LineSegments(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(rings, 3)),
    new THREE.LineBasicMaterial({ color: 0x0b0d10, transparent: true, opacity: 0.35, depthWrite: false }));
  rl.material.clippingPlanes = [clip];
  root.add(rl);

  // --- cavity spacer (transparent, holds the dot)
  const spM = layerMaterial('core', { clip, opacity: 0.5 });
  spM.color.setHex(0x6f8fae);
  const spacer = new THREE.Mesh(new THREE.CylinderGeometry(R, R, tS, seg, 1, false), spM);
  spacer.position.y = (ySpacer0 + ySpacer1) / 2;
  const spEdge = new THREE.LineSegments(new THREE.EdgesGeometry(spacer.geometry, 20), new THREE.LineBasicMaterial({ color: 0xd8d6cf, transparent: true, opacity: 0.55, depthWrite: false }));
  spEdge.material.clippingPlanes = [clip];
  spacer.add(spEdge);
  stage.addPickable(spacer, { title: 'cavity spacer', tag: 'DR', rows: [['thickness', `${nice(g.spacer_nm, 4)} nm`], ['order', `${nice(g.cavity_order)} λ/n_high`], ['n', nice(g.n_high, 4)]], note: 'lambda/n of the high-index layer; order is a display default [A]' });
  root.add(spacer);

  // --- substrate (thickness not modelled; broken)
  const subT = g.substrate_um * S, subR = R * 2.6;
  root.add(composeSubstrate(stage, subR, subT, 0, { info: { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']], note: 'drawn as a plate under the bottom DBR' } }));
  root.add(composeRingBreak(subR, -subT * 0.5, 0.06));

  // --- mode envelope (schematic) and dot
  const yd = ySpacer0 + tS * g.dot_frac;
  const w = R * 0.72, ym = (yTop - 0) * 0.26;
  const prof = [];
  for (let i = 0; i <= 40; i++) {
    const t = i / 40, yy = yd + (t - 0.5) * (yTop - 0) * 0.75;
    prof.push([Math.max(0.003, w * Math.exp(-Math.pow((yy - yd) / ym, 2))), yy]);
  }
  const mode = composeLathe(prof, glowMaterial(rgb, { opacity: 0.1 }), 48);
  stage.markBloom(mode);
  root.add(mode);
  composeDot(stage, root, rgb, [0, yd, 0], Math.max(R * 0.07, 0.03),
    { title: 'quantum dot', tag: tg('dot_frac'), rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['position', 'spacer centre (display default)']], note: 'dot depth is not a model input' });
  composeBeam(stage, root, [0, yTop + 0.02, 0], [0, 1, 0], 5, Math.max(R * 0.35, 0.12));

  // --- annotations
  const lx = R + 0.25;
  labels.add(root, [lx, yTop - nm(g.N_top * (g.t_high_nm + g.t_low_nm)) / 2, 0], `top DBR · ${g.N_top} pairs · ${nice(g.N_top * (g.t_high_nm + g.t_low_nm) / 1000, 3)} µm`, { tag: tg('N_top'), priority: 70 });
  labels.add(root, [lx, yBotTop / 2, 0], `bottom DBR · ${g.N_bot} pairs · ${nice(g.N_bot * (g.t_high_nm + g.t_low_nm) / 1000, 3)} µm`, { tag: tg('N_bot'), priority: 70 });
  labels.add(root, [lx, (ySpacer0 + ySpacer1) / 2, 0], `spacer · ${nice(g.spacer_nm, 4)} nm`, { tag: 'DR', priority: 85, sub: 'λ/n_high; cavity order is a display default' });
  labels.add(root, [-R - 0.2, yd, 0], `dot · ${nice(lam, 4)} nm`, { tag: tg('dot_frac'), anchor: 'right', priority: 90, sub: falseColour ? colourNote : 'depth not a model input' });
  labels.add(root, [0, yTop + 0.35, R], `pillar Ø ${nice(g.d_um)} µm · H ${nice(g.height_um, 3)} µm`, { tag: tg('d_um'), anchor: 'center', priority: 80 });
  labels.add(root, [-R - 0.2, yTop - tH * 3, 0], `λ/4n · ${nice(g.t_high_nm, 3)} / ${nice(g.t_low_nm, 3)} nm`, { tag: 'DR', anchor: 'right', priority: 60, sub: `n ${nice(g.n_high, 3)} / ${nice(g.n_low, 3)}` });
  labels.add(root, [subR * 0.75, -subT * 0.5, subR * 0.6], 'substrate', { tag: 'A', priority: 30, sub: 'thickness not modelled' });
  if (g.illustrative) labels.add(root, [0, yTop + 1.05, 0], 'illustrative geometry: the model is lumped', { tag: 'A', anchor: 'center', priority: 99, sub: 'only F_P, κ, η_out are model inputs; the pillar is a display default' });

  // --- magnifier: true-scale layers around the spacer
  const mag = makeMagnifier(stage, hud, {
    mag: 6, title: 'DBR layers at the spacer (true scale)',
    buildInset(scene) {
      const q = (wd, ht, x, yy, zz, c) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(wd, ht), new THREE.MeshBasicMaterial({ color: c, toneMapped: false })); o.position.set(x, yy, zz); scene.add(o); };
      const Wn = g.d_um * 1000;
      q(Wn, g.spacer_nm, 0, 0, 0, 0x6f8fae);
      const rgbC = new THREE.Color(...rgb).getHex();
      const dotC = new THREE.Mesh(new THREE.CircleGeometry(14, 24), new THREE.MeshBasicMaterial({ color: rgbC, toneMapped: false }));
      dotC.position.set(0, (g.dot_frac - 0.5) * g.spacer_nm, 0.2);
      scene.add(dotC);
      for (const s of [1, -1]) {
        let yy = s * g.spacer_nm / 2;
        for (let i = 0; i < 6; i++) {
          const hi = i % 2 === 1, t = hi ? g.t_high_nm : g.t_low_nm;
          q(Wn, t, 0, yy + s * t / 2, 0, hi ? 0xb9c4d2 : 0x566273);
          yy += s * t;
        }
      }
      return { center: [0, 0], insetUnitsPerWorld: 1000 / S, note: `spacer ${nice(g.spacer_nm, 4)} nm; λ/4n layers ${nice(g.t_high_nm, 3)} / ${nice(g.t_low_nm, 3)} nm` };
    },
  });
  const gv = hud.group('view');
  hud.button(gv, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the layer magnifier' });

  const tc = [0, yTop * 0.5, 0];
  const dMain = yTop * 2.15;
  const v = (t, dir, d) => composeView(t, dir, d);
  const main = v(tc, [0.62, 0.22, 0.86], dMain);
  const bm = (name, vw, cut, magOn = false) => ({ name, pos: vw.pos, target: vw.target, cut, mag: magOn });
  return {
    root,
    defaultView: main,
    bookmarks: [
      bm('Pillar', main, subR * 1.05),
      bm('Side cut-away', v([0, yd, 0], [1, 0.04, 0.05], R * 13), 0),
      bm('Top / beam', v([0, yTop, 0], [0.45, 0.55, 0.7], R * 8), subR * 1.05),
      bm('Spacer magnifier', v([0, yd, 0], [0.5, 0.12, 0.85], R * 12), subR * 1.05, true),
    ],
    cut: { min: -R, max: subR * 1.05 },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
