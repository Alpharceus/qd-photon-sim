// Mesa / slab with optional SIL and metal reflector (composed device /
// structures sil_mesa, mesa_reflector, bulk_planar).  The registry entries
// carry efficiency tables and (bulk_planar) n, mirror, sil, n_s only, so the
// mesa and slab sizes, SIL radius, dot depth and reflector thickness are
// display defaults [A] and labelled so; the SIL is not to scale against a real
// ~mm lens.  bulk_planar draws the escape cone asin(1/n) from its n (spec).
// Dot at the stated depth, beam up through the lens.
import * as THREE from 'three';
import { layerMaterial, coneMaterial, outline, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeRingBreak, composeSubstrate } from './micropillar.js';

export function buildMesa(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const slab = g.mode === 'slab';
  const sil = g.sil;
  const silR = g.sil_radius_um;
  const S = 9 / Math.max(g.d_um, (sil !== 'none' ? 2 * silR : 0), 2);          // world per um
  const um = (v) => v * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();
  const R = um(g.d_um / 2), h = um(g.h_um);
  const refT = g.mirror ? um(g.reflector_um) : 0;
  const plateR = R * (slab ? 1.0 : 2.4);
  const subT = um(1.2);

  root.add(composeSubstrate(stage, plateR, subT, 0, { info: { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']] } }));
  root.add(composeRingBreak(plateR, -subT / 2, 0.07));

  if (g.mirror) {
    const rp = slab ? plateR : R * 1.5;
    const rf = new THREE.Mesh(new THREE.CylinderGeometry(rp, rp, refT, 96), new THREE.MeshPhysicalMaterial({ color: 0xa8832f, metalness: 0.7, roughness: 0.55 }));
    rf.material.clippingPlanes = [clip];
    rf.position.y = refT / 2;
    stage.addPickable(rf, { title: 'metal reflector', tag: tg('mirror'), rows: [['thickness', `${nice(g.reflector_um * 1000, 3)} nm (illustrative)`]], note: spec.structure_kind === 'bulk_planar' ? 'back mirror flag from the spec' : 'the entry kind implies a reflector' });
    root.add(rf);
  }

  const body = new THREE.Mesh(new THREE.CylinderGeometry(R, R, h, 96), layerMaterial('core', { clip, opacity: 0.55 }));
  body.material.color.setHex(0x86a3c2);
  body.position.y = refT + h / 2;
  outline(body, { clip, opacity: 0.45 });
  stage.addPickable(body, { title: slab ? 'planar slab' : 'mesa', tag: 'A', rows: [['diameter', `${nice(g.d_um)} µm`], ['height', `${nice(g.h_um)} µm`], ...(slab ? [['n', nice(g.n, 4)]] : [])], note: 'size is a display default [A]' + (slab ? '; the slab is laterally extended (break mark)' : '') });
  root.add(body);
    const yTop = refT + h, yd = yTop - g.dot_depth_frac * h;

  // escape cone (bulk_planar: n from the spec)
  let escapeDeg = null;
  if (slab && g.n > 1) {
    const th = Math.asin(1 / g.n), hh = yTop - yd;
    const cone = new THREE.Mesh(new THREE.ConeGeometry(hh * Math.tan(th), hh, 64, 1, true), coneMaterial(rgb, 0.16));
    cone.rotation.x = Math.PI;
    cone.position.set(0, yd + hh / 2, 0);
    root.add(cone);
    escapeDeg = (th * 180) / Math.PI;
  }

  // SIL
  let yExit = yTop;
  if (sil !== 'none') {
    const Rs = um(silR);
    const ns = g.n_s;
    const thMax = sil === 'hemi' ? Math.PI / 2 : Math.acos(-1 / ns);
    const yc = sil === 'hemi' ? yTop : yTop + Rs / ns;
    const lens = new THREE.Mesh(new THREE.SphereGeometry(Rs, 96, 48, 0, Math.PI * 2, 0, thMax),
      layerMaterial('cladding', { clip, opacity: 0.3 }));
    lens.material.side = THREE.DoubleSide;
    lens.material.clearcoat = 0.6;
    lens.position.y = yc;
    stage.addPickable(lens, { title: `${sil === 'hemi' ? 'hemispherical' : 'Weierstrass'} SIL`, tag: tg('sil'), rows: [['type', sil], ['n_s', nice(ns, 3)], ['radius', `${nice(silR)} µm (display)`], ['height', `${nice(g.sil_height_um, 3)} µm`]], note: 'radius is a display rule; a real lens is ~mm across' });
    root.add(lens);
    const rimPts = [];
    const baseR = Rs * Math.sin(thMax);
    for (let i = 0; i <= 96; i++) { const t = (i / 96) * Math.PI * 2; rimPts.push(new THREE.Vector3(Math.cos(t) * baseR, yc + Rs * Math.cos(thMax), Math.sin(t) * baseR)); }
    const rim = new THREE.Line(new THREE.BufferGeometry().setFromPoints(rimPts), new THREE.LineBasicMaterial({ color: 0xaeb4bd, transparent: true, opacity: 0.6 }));
    rim.material.clippingPlanes = [clip];
    root.add(rim);
    yExit = yc + Rs;
  }

  composeDot(stage, root, rgb, [0, yd, 0], Math.max(um(0.05), 0.05),
    { title: 'quantum dot', tag: 'A', rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['depth', `${nice(g.dot_depth_frac * g.h_um)} µm below the top (display)`]], note: 'dot depth is not a model input' });
  composeBeam(stage, root, [0, yExit + 0.02, 0], [0, 1, 0], 5, Math.max(R * 0.07, 0.15));

  // annotations
  const kind = spec.structure_kind;
  labels.add(root, [R + 0.2, refT + h * 0.5, 0], `${slab ? 'slab' : 'mesa'} Ø ${nice(g.d_um)} µm · ${nice(g.h_um)} µm`, { tag: 'A', priority: 85, sub: 'display default size' });
  labels.add(root, [-R - 0.2, yd, 0], `dot · ${nice(lam, 4)} nm`, { tag: 'A', anchor: 'right', priority: 90, sub: falseColour ? colourNote : `${nice(g.dot_depth_frac * g.h_um)} µm below the top` });
  if (sil !== 'none') labels.add(root, [um(silR) * 0.7, yExit - um(silR) * 0.25, um(silR) * 0.5], `${sil === 'hemi' ? 'hemispherical' : 'Weierstrass'} SIL · n ${nice(g.n_s, 3)}`, { tag: kind === 'bulk_planar' ? tg('sil') : 'A', priority: 80, sub: 'radius not to scale (real lens ~mm)' });
  if (g.mirror) labels.add(root, [plateR * 0.8, refT / 2, plateR * 0.4], 'metal reflector', { tag: kind === 'bulk_planar' ? tg('mirror') : 'A', priority: 55 });
  if (escapeDeg) labels.add(root, [-R * 0.5, (yd + yTop) / 2, R * 0.5], `escape cone ${nice(escapeDeg, 3)}° (n ${nice(g.n, 3)})`, { tag: tg('n'), anchor: 'right', priority: 70, sub: 'asin(1/n); total internal reflection outside' });
  labels.add(root, [plateR * 0.75, -subT / 2, plateR * 0.5], 'substrate', { tag: 'A', priority: 25, sub: 'thickness not modelled' });

  const D = Math.max(R * 2, sil !== 'none' ? um(silR) * 2 : 0) * 2.2;
  const yc = yExit * 0.45;
  const cutMax = Math.max(plateR, sil !== 'none' ? um(silR) : 0) * 1.1;
  return {
    root,
    defaultView: { pos: [D * 0.78, yc + D * 0.36, D * 1.0], target: [0, yc * 0.9, 0] },
    bookmarks: [
      { name: slab ? 'Slab' : 'Mesa', pos: [D * 0.78, yc + D * 0.36, D * 1.0], target: [0, yc * 0.9, 0], cut: cutMax, mag: false },
      { name: 'Side cut-away', pos: [D * 1.1, yc, 0.2], target: [0, yc, 0], cut: 0, mag: false },
      { name: 'Top / beam', pos: [0.5, yExit + D * 0.9, D * 0.35], target: [0, yTop, 0], cut: cutMax, mag: false },
      { name: 'Dot close-up', pos: [R * 1.6, yd + R * 0.7, R * 2.0], target: [0, yd, 0], cut: 0, mag: false },
    ],
    cut: { min: -Math.max(R, sil !== 'none' ? um(silR) : 0), max: cutMax },
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
