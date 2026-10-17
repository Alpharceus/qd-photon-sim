// Circular Bragg grating / bullseye (composed device / structure cbg): a
// centre disk and concentric rings in a thin membrane over an oxide spacer
// and a metal back reflector.  The registry entry holds F_P, Q and an
// efficiency table only, so the ring period (lambda / n_eff), ring count,
// fill, membrane and spacer are display defaults [A]; the wavelength is the
// material's.  Everything is drawn at one uniform scale (true to those
// numbers).  The dot sits at the centre of the disk, mid-membrane.
import * as THREE from 'three';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';
import { layerMaterial, glowMaterial, outline, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeRingBreak, composeSubstrate } from './micropillar.js';

// Flat ring (inner radius ri, outer ro) standing on y0 with thickness t.
function ringGeo(ri, ro, y0, t) {
  const sh = new THREE.Shape();
  sh.absarc(0, 0, ro, 0, Math.PI * 2, false);
  if (ri > 0) { const h = new THREE.Path(); h.absarc(0, 0, ri, 0, Math.PI * 2, true); sh.holes.push(h); }
  const geo = new THREE.ExtrudeGeometry(sh, { depth: t, bevelEnabled: false, curveSegments: 96 });
  geo.rotateX(-Math.PI / 2);
  geo.translate(0, y0, 0);
  return geo;
}

export function buildBullseye(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const outer = g.outer_radius_um;
  const S = 6 / outer;                           // world units per um
  const nm = (v) => (v / 1000) * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();
  const tM = nm(g.membrane_nm), tSp = nm(g.spacer_nm), tRef = g.reflector ? nm(120) : 0;
  const R0 = nm(g.disk_radius_nm), P = nm(g.period_nm), W = nm(g.ring_width_nm), trench = nm(g.trench_nm);

  // membrane: centre disk + rings (one merged mesh)
  const parts = [ringGeo(0, R0, 0, tM)];
  const ringRadii = [];
  for (let k = 1; k <= g.n_rings; k++) {
    const ri = R0 + (k - 1) * P + trench, ro = ri + W;
    parts.push(ringGeo(ri, ro, 0, tM));
    ringRadii.push([ri, ro]);
  }
  const memMat = layerMaterial('core', { clip, opacity: 1 });
  memMat.color.setHex(0x8aa6c4);
  const mem = new THREE.Mesh(mergeGeometries(parts), memMat);
  mem.position.y = tSp + tRef;
  outline(mem, { clip, opacity: 0.5 });
  stage.addPickable(mem, { title: 'membrane: disk and rings', tag: 'A', rows: [['period', `${nice(g.period_nm, 4)} nm`], ['ring width', `${nice(g.ring_width_nm, 4)} nm`], ['disk radius', `${nice(g.disk_radius_nm, 4)} nm`], ['rings', `${g.n_rings}`], ['thickness', `${nice(g.membrane_nm, 4)} nm`]], note: 'ring geometry is a display default [A]; period = λ / n_eff' });
  root.add(mem);

  // oxide spacer + reflector + wafer, all as wide discs (laterally broken)
  const plateR = (outer * 1.35) * S;
  const sp = new THREE.Mesh(new THREE.CylinderGeometry(plateR, plateR, tSp, 96), layerMaterial('oxide', { clip, opacity: 0.55 }));
  sp.position.y = tRef + tSp / 2;
  stage.addPickable(sp, { title: 'oxide spacer', tag: 'A', rows: [['thickness', `${nice(g.spacer_nm, 4)} nm`]], note: 'display default' });
  root.add(sp);
  if (g.reflector) {
    const rf = new THREE.Mesh(new THREE.CylinderGeometry(plateR, plateR, tRef, 96), new THREE.MeshPhysicalMaterial({ color: 0xa8832f, metalness: 0.7, roughness: 0.55 }));
    rf.material.clippingPlanes = [clip];
    rf.position.y = tRef / 2;
    stage.addPickable(rf, { title: 'metal back reflector', tag: 'A', rows: [['thickness', 'illustrative']], note: 'display default' });
    root.add(rf);
  }
  const subT = 0.55;
  root.add(composeSubstrate(stage, plateR, subT, 0, { info: { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']] } }));
  root.add(composeRingBreak(plateR, -subT / 2, 0.05, 36));

  // dot and mode glow at the disk centre
  const yd = tSp + tRef + tM / 2;
  const modeG = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 16), glowMaterial(rgb, { opacity: 0.16 }));
  modeG.scale.set(R0 * 0.85, tM * 1.7, R0 * 0.85);
  modeG.position.set(0, yd, 0);
  stage.markBloom(modeG);
  root.add(modeG);
  composeDot(stage, root, rgb, [0, yd, 0], Math.max(R0 * 0.05, 0.03),
    { title: 'quantum dot', tag: 'A', rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['position', 'disk centre (display)']], note: 'dot position is not a model input' });
  composeBeam(stage, root, [0, tSp + tRef + tM + 0.02, 0], [0, 1, 0], 5, Math.max(R0 * 0.55, 0.2));

  // annotations
  const yTop = tSp + tRef + tM;
  labels.add(root, [R0 * 0.4, yTop + 0.18, -R0 * 0.7], `disk r ${nice(g.disk_radius_nm, 3)} nm`, { tag: 'A', priority: 85, anchor: 'left' });
  const [r1i, r1o] = ringRadii[0];
  labels.add(root, [(r1i + r1o) / 2, yTop + 0.2, 0], `period ${nice(g.period_nm, 4)} nm · ${g.n_rings} rings`, { tag: 'A', priority: 90, sub: 'Λ = λ / n_eff, n_eff a display default' });
  labels.add(root, [ringRadii[g.n_rings - 1][1] + 0.1, yTop, 0], `membrane ${nice(g.membrane_nm, 3)} nm`, { tag: 'A', priority: 60, sub: 'etched through' });
  labels.add(root, [plateR * 0.95, tRef + tSp / 2, plateR * 0.2], `oxide ${nice(g.spacer_nm, 3)} nm`, { tag: 'A', priority: 45 });
  if (g.reflector) labels.add(root, [plateR * 0.95, tRef / 2, -plateR * 0.2], 'metal reflector', { tag: 'A', priority: 40 });
  labels.add(root, [-plateR * 0.7, 0, plateR * 0.6], 'substrate', { tag: 'A', anchor: 'right', priority: 25, sub: 'thickness not modelled' });
  labels.add(root, [0, yTop + 1.5, 0], `λ ${nice(lam, 4)} nm`, { tag: 'DR', anchor: 'center', priority: 70, sub: falseColour ? colourNote : 'material emission wavelength' });

  // magnifier: radial section through disk and first rings (nm); the magnifier
  // camera sits at the origin, so the content is shifted to centre it
  const mag = makeMagnifier(stage, hud, {
    mag: 6, title: 'radial section (true scale)',
    buildInset(scene) {
      const grp = new THREE.Group();
      grp.position.y = -(g.membrane_nm / 2 - g.spacer_nm * 0.2);
      scene.add(grp);
      const q = (wd, ht, x, yy, zz, c) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(wd, ht), new THREE.MeshBasicMaterial({ color: c, toneMapped: false })); o.position.set(x, yy, zz); grp.add(o); };
      const rN = g.disk_radius_nm, total = g.disk_radius_nm + g.period_nm * Math.min(g.n_rings, 4) + 100;
      q(total * 2, g.spacer_nm, 0, -g.spacer_nm / 2, 0, 0x56626e);
      q(2 * rN, g.membrane_nm, 0, g.membrane_nm / 2, 0.1, 0x8aa6c4);
      for (let k = 1; k <= Math.min(g.n_rings, 4); k++) {
        const ri = rN + (k - 1) * g.period_nm + g.trench_nm, c = ri + g.ring_width_nm / 2;
        for (const s of [1, -1]) q(g.ring_width_nm, g.membrane_nm, s * c, g.membrane_nm / 2, 0.1, 0x8aa6c4);
      }
      q(40, 40, 0, g.membrane_nm / 2, 0.3, new THREE.Color(...rgb).getHex());
      return { center: [0, 0], insetUnitsPerWorld: 1000 / S, note: `Λ ${nice(g.period_nm, 4)} nm, ring ${nice(g.ring_width_nm, 3)} nm, trench ${nice(g.trench_nm, 3)} nm, membrane ${nice(g.membrane_nm, 3)} nm` };
    },
  });
  const gv = hud.group('view');
  hud.button(gv, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the radial-section magnifier' });

  const D = outer * S * 3.3;
  const cutMax = plateR * 1.1;
  return {
    root,
    defaultView: { pos: [D * 0.62, D * 0.55, D * 0.78], target: [0, 0.2, 0] },
    bookmarks: [
      { name: 'Grating', pos: [D * 0.62, D * 0.55, D * 0.78], target: [0, 0.2, 0], cut: cutMax, mag: false },
      { name: 'Side cut-away', pos: [D * 0.9, 0.5, 0.3], target: [0, 0.2, 0], cut: 0, mag: false },
      { name: 'Top', pos: [0.2, D * 1.1, 0.4], target: [0, 0, 0], cut: cutMax, mag: false },
      { name: 'Section magnifier', pos: [D * 0.3, D * 0.2, D * 0.5], target: [0, 0.2, 0], cut: 0, mag: true },
    ],
    cut: { min: -plateR, max: cutMax },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
