// Open / fibre micro-cavity (composed device / structure open_cavity): a
// fibre tip with a concave dimple (radius of curvature ROC) over a planar
// DBR chip, cavity length L measured from the membrane surface to the dimple
// apex.  True scale in every axis; ROC, L, the TEM00 waist (w0^2 =
// lambda/pi sqrt(L (ROC - L))) and the quarter-wave DBR come from the
// descriptor.  The fibre diameter / core, the dimple radius (rim L/2 below
// the apex), the DBR pair count and the chip width are display defaults [A].
// The mode is drawn only in the gap and the fibre core: it is the active path.
import * as THREE from 'three';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';
import { layerMaterial, glowMaterial, fatLine, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeLathe, composeView, composeRingBreak, composeSubstrate } from './micropillar.js';

export function buildOpenCavity(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const S = 24 / g.fibre_d_um;                     // world units per um
  const um = (v) => v * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();
  const Rf = um(g.fibre_d_um / 2), Rc = um(g.fibre_core_um / 2), chipR = um(g.chip_halfwidth_um);
  const L = um(g.L_um), roc = g.ROC_um, a = g.dimple_radius_um;
  const dimpleY = (r) => g.L_um - (roc - Math.sqrt(Math.max(roc * roc - r * r, 0)));     // um
  const yFace = um(g.L_um - g.sag_um);
  const yFibTop = yFace + um(g.fibre_len_um);

  // --- fibre: glass body (open cylinder + top cap), the flat end-face annulus and the dimple bowl
  const glass = layerMaterial('cladding', { clip, opacity: 0.22 });
  glass.side = THREE.DoubleSide;
  const fibre = new THREE.Group();
  const wall = new THREE.Mesh(new THREE.CylinderGeometry(Rf, Rf, yFibTop - yFace, 96, 1, true), glass);
  wall.position.y = (yFibTop + yFace) / 2;
  const topCap = new THREE.Mesh(new THREE.CircleGeometry(Rf, 96), glass);
  topCap.rotation.x = -Math.PI / 2; topCap.position.y = yFibTop;
  const face = new THREE.Mesh(new THREE.RingGeometry(um(a), Rf, 96, 1), glass);
  face.rotation.x = Math.PI / 2; face.position.y = yFace;
  const bowlP = [];
  for (let i = 0; i <= 24; i++) { const r = (a * i) / 24; bowlP.push([um(r), um(dimpleY(r))]); }
  const bowl = composeLathe(bowlP.reverse(), glass, 96);
  fibre.add(wall, topCap, face, bowl);
  const finfo = { title: 'fibre', tag: tg('fibre_d_um'), rows: [['diameter', `${nice(g.fibre_d_um)} µm`], ['dimple ROC', `${nice(roc)} µm`], ['dimple Ø', `${nice(2 * a)} µm`], ['sag', `${nice(g.sag_um)} µm`]], note: 'fibre diameter and dimple radius are display defaults [A]; ROC is the spec value' };
  stage.addPickable(wall, finfo);
  root.add(fibre);
  const rim = [];
  for (let i = 0; i <= 96; i++) { const t = (i / 96) * Math.PI * 2; rim.push([Math.cos(t) * Rf, yFace, Math.sin(t) * Rf]); }
  root.add(fatLine(rim, { color: 0xaeb4bd, width: 1, opacity: 0.55 }));
  const rimD = [];
  for (let i = 0; i <= 64; i++) { const t = (i / 64) * Math.PI * 2; rimD.push([Math.cos(t) * um(a), yFace, Math.sin(t) * um(a)]); }
  root.add(fatLine(rimD, { color: 0xe9e8e3, width: 1.4, opacity: 0.9 }));
  root.add(composeRingBreak(Rf, yFibTop - 0.6, 0.12, 36));

  // --- fibre-side mirror coating on the dimple (thin shell, thickness illustrative)
  const ct = um(g.coating_um);
  const cprof = [];
  for (let i = 0; i <= 24; i++) { const r = a * (i / 24); cprof.push([um(r), um(dimpleY(r))]); }
  for (let i = 24; i >= 0; i--) { const r = a * (i / 24); cprof.push([um(r), um(dimpleY(r)) + ct]); }
  const coatMat = layerMaterial('mirror', { clip, opacity: 1 }); coatMat.color.setHex(0xd7dde6); coatMat.side = THREE.DoubleSide;
  const coat = composeLathe(cprof, coatMat, 96);
  stage.addPickable(coat, { title: 'fibre-side mirror coating', tag: 'A', rows: [['thickness', `${nice(g.coating_um * 1000)} nm (illustrative)`]], note: 'coating thickness and stack are not model inputs' });
  root.add(coat);

  // --- single-mode core above the dimple (active path only)
  const core = new THREE.Mesh(new THREE.CylinderGeometry(Rc, Rc, yFibTop - L, 32, 1, true), glowMaterial(rgb, { opacity: 0.009 }));
  core.position.y = (yFibTop + L) / 2;
  stage.markBloom(core);
  root.add(core);

  // --- planar chip: membrane (holds the dot) over a DBR, then the substrate
  const tMem = um(g.membrane_nm / 1000), tHn = um(g.t_high_nm / 1000), tLn = um(g.t_low_nm / 1000);
  const hiG = [], loG = [];
  let y = -tMem;
  for (let i = 0; i < g.dbr_pairs; i++) {
    y -= tLn; const lo = new THREE.CylinderGeometry(chipR, chipR, tLn, 96, 1, true); lo.translate(0, y + tLn / 2, 0); loG.push(lo);
    y -= tHn; const hi = new THREE.CylinderGeometry(chipR, chipR, tHn, 96, 1, true); hi.translate(0, y + tHn / 2, 0); hiG.push(hi);
  }
  const yDbr = y;
  const hiM = layerMaterial('mirror', { clip, opacity: 1 }); hiM.color.setHex(0xb9c4d2);
  const loM = layerMaterial('mirror', { clip, opacity: 1 }); loM.color.setHex(0x566273);
  const hiMesh = new THREE.Mesh(mergeGeometries(hiG), hiM), loMesh = new THREE.Mesh(mergeGeometries(loG), loM);
  const dinfo = (t, n) => ({ title: 'planar DBR', tag: 'A', rows: [['pairs drawn', `${g.dbr_pairs}`], ['layer', `λ/4n = ${nice(t, 4)} nm (n ${n} display)`]], note: 'pair count and indices are display defaults [A]' });
  stage.addPickable(hiMesh, dinfo(g.t_high_nm, g.n_high)); stage.addPickable(loMesh, dinfo(g.t_low_nm, g.n_low));
  const mem = new THREE.Mesh(new THREE.CylinderGeometry(chipR, chipR, tMem, 96), layerMaterial('core', { clip, opacity: 0.85 }));
  mem.material.color.setHex(0x4f6c8c);
  mem.position.y = -tMem / 2;
  stage.addPickable(mem, { title: 'membrane', tag: 'A', rows: [['thickness λ/n', `${nice(g.membrane_nm, 4)} nm`]], note: 'holds the dot at an antinode; thickness derived from the display index' });
  const subT = um(6);
  const sub = composeSubstrate(stage, chipR, subT, yDbr, { info: { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']] } });
  root.add(hiMesh, loMesh, mem, sub);
  root.add(composeRingBreak(chipR, yDbr - subT * 0.5, 0.2, 60));

  // --- TEM00 mode in the gap (waist on the planar mirror) and the dot
  let modeW = Rc * 0.6;
  if (g.stable && g.waist_um) {
    const w0 = g.waist_um, zr = g.rayleigh_um;
    const mp = [];
    for (let i = 0; i <= 28; i++) { const z = (i / 28) * g.L_um; mp.push([um(w0 * Math.sqrt(1 + (z / zr) ** 2)), um(z)]); }
    modeW = mp[0][0];
    const m = composeLathe(mp, glowMaterial(rgb, { opacity: 0.2 }), 48);
    stage.markBloom(m);
    stage.addPickable(m, { title: 'TEM00 mode', tag: 'DR', rows: [['waist w0', `${nice(w0, 3)} µm`], ['Rayleigh z_R', `${nice(zr, 3)} µm`], ['ROC / L', `${nice(roc)} / ${nice(g.L_um)} µm`]], note: 'w0² = (λ/π) √(L (ROC − L)); envelope w(z) = w0 √(1 + z²/z_R²)' });
    root.add(m);
  }
  const yDot = -tMem / 2;
  composeDot(stage, root, rgb, [0, yDot, 0], Math.max(um(0.02), 0.04),
    { title: 'quantum dot', tag: 'A', rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['position', 'membrane centre (display)']], note: 'dot position is not a model input' });
  composeBeam(stage, root, [0, yFibTop, 0], [0, 1, 0], 5, Math.max(Rc * 1.2, 0.2));

  // --- section lines (shown when the cut passes through the axis): the dimple meridian, the coating, the
  // end-face, the fibre walls and core, and every chip interface, all in the plane x = 0
  const section = new THREE.Group();
  const SL = (pts, o = {}) => section.add(fatLine(pts, { color: 0xe9e8e3, width: 1.3, opacity: 0.85, ...o }));
  const mer = [], merC = [];
  for (let i = -24; i <= 24; i++) { const r = (a * i) / 24; mer.push([0, um(dimpleY(Math.abs(r))), um(r)]); merC.push([0, um(dimpleY(Math.abs(r))) + ct, um(r)]); }
  SL(mer, { width: 2, opacity: 1 }); SL(merC, { width: 1, opacity: 0.6 });
  for (const s of [-1, 1]) {
    SL([[0, yFace, s * um(a)], [0, yFace, s * Rf]]);
    SL([[0, yFace, s * Rf], [0, yFibTop, s * Rf]], { opacity: 0.5 });
    SL([[0, L, s * Rc], [0, yFibTop, s * Rc]], { opacity: 0.45, dashed: true, dashSize: 0.3, gapSize: 0.2 });
  }
  const ifc = [];
  { let yy = 0; ifc.push(0, yy, -chipR, 0, yy, chipR);
    yy = -tMem; ifc.push(0, yy, -chipR, 0, yy, chipR);
    for (let i = 0; i < g.dbr_pairs; i++) { yy -= tLn; ifc.push(0, yy, -chipR, 0, yy, chipR); yy -= tHn; ifc.push(0, yy, -chipR, 0, yy, chipR); } }
  const ifl = new THREE.LineSegments(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(ifc, 3)),
    new THREE.LineBasicMaterial({ color: 0xe9e8e3, transparent: true, opacity: 0.55 }));
  section.add(ifl);
  section.visible = false;
  root.add(section);
  stage.onFrame.add(() => { section.visible = Math.abs(stage.clipPlane.constant) < 0.05; });

  // --- annotations
  labels.add(root, [Rf + 0.3, (yFace + yFibTop) / 2, 0], `fibre Ø ${nice(g.fibre_d_um)} µm · core Ø ${nice(g.fibre_core_um)} µm`, { tag: 'A', priority: 60, sub: 'display defaults' });
  labels.add(root, [um(a) + 0.2, yFace + 0.05, um(a) * 0.4], `dimple Ø ${nice(2 * a, 3)} µm · ROC ${nice(roc)} µm`, { tag: tg('ROC_um'), priority: 85, sub: 'ROC is the spec value; the dimple radius is a display rule' });
  labels.add(root, [Rc * 3, L * 0.55, -um(a) * 0.3], `L = ${nice(g.L_um)} µm`, { tag: tg('L_um'), priority: 90, sub: 'membrane surface to dimple apex' });
  if (g.stable && g.waist_um) labels.add(root, [-modeW * 2.2, L * 0.15, 0], `waist w0 ${nice(g.waist_um, 3)} µm`, { tag: 'DR', anchor: 'right', priority: 80, sub: 'TEM00, on the planar mirror' });
  else labels.add(root, [-Rc * 3, L * 0.5, 0], 'unstable (L ≥ ROC): no mode drawn', { tag: 'DR', anchor: 'right', priority: 80 });
  labels.add(root, [chipR * 0.55, yDbr / 2, chipR * 0.55], `planar DBR · ${g.dbr_pairs} pairs`, { tag: 'A', priority: 50, sub: 'λ/4n layers, display pairs' });
  labels.add(root, [-chipR * 0.5, -tMem / 2, chipR * 0.6], `membrane λ/n ${nice(g.membrane_nm, 3)} nm · dot`, { tag: 'A', anchor: 'right', priority: 55 });
  labels.add(root, [chipR * 0.55, yDbr - subT * 0.5, chipR * 0.55], 'substrate', { tag: 'A', priority: 25, sub: 'thickness not modelled' });

  // --- magnifier: section through the gap in µm
  const mag = makeMagnifier(stage, hud, {
    mag: 2.5, title: 'section through the gap (true scale)', unitLabel: 'µm',
    buildInset(scene) {
      // the magnifier camera sits at the origin; content is shifted so the gap is centred
      const grp = new THREE.Group();
      grp.position.y = -(g.L_um / 2 - 0.2);
      scene.add(grp);
      const q = (wd, ht, x, yy, zz, c, op = 1) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(wd, ht), new THREE.MeshBasicMaterial({ color: c, toneMapped: false, transparent: op < 1, opacity: op })); o.position.set(x, yy, zz); grp.add(o); return o; };
      const Wd = 2 * a * 1.6;
      q(Wd, g.membrane_nm / 1000, 0, -g.membrane_nm / 2000, 0, 0x6f8fae);
      let yy = -g.membrane_nm / 1000;
      for (let i = 0; i < g.dbr_pairs; i++) {
        q(Wd, g.t_low_nm / 1000, 0, yy - g.t_low_nm / 2000, 0, 0x566273); yy -= g.t_low_nm / 1000;
        q(Wd, g.t_high_nm / 1000, 0, yy - g.t_high_nm / 2000, 0, 0xb9c4d2); yy -= g.t_high_nm / 1000;
      }
      const sh = new THREE.Shape();
      sh.moveTo(-Wd / 2, g.L_um - g.sag_um + 6);
      sh.lineTo(-Wd / 2, g.L_um - g.sag_um); sh.lineTo(-a, g.L_um - g.sag_um);
      for (let i = 1; i <= 48; i++) { const r = -a + (2 * a * i) / 48; sh.lineTo(r, dimpleY(Math.abs(r))); }
      sh.lineTo(Wd / 2, g.L_um - g.sag_um); sh.lineTo(Wd / 2, g.L_um - g.sag_um + 6);
      const gm = new THREE.Mesh(new THREE.ShapeGeometry(sh), new THREE.MeshBasicMaterial({ color: 0x8fa3b8, toneMapped: false }));
      gm.position.z = 0.05; grp.add(gm);
      if (g.stable && g.waist_um) {
        const s2 = new THREE.Shape(), n = 40;
        for (let i = 0; i <= n; i++) { const z = (i / n) * g.L_um; const w = g.waist_um * Math.sqrt(1 + (z / g.rayleigh_um) ** 2); i ? s2.lineTo(w, z) : s2.moveTo(w, z); }
        for (let i = n; i >= 0; i--) { const z = (i / n) * g.L_um; s2.lineTo(-g.waist_um * Math.sqrt(1 + (z / g.rayleigh_um) ** 2), z); }
        const mm = new THREE.Mesh(new THREE.ShapeGeometry(s2), new THREE.MeshBasicMaterial({ color: new THREE.Color(...rgb), toneMapped: false, transparent: true, opacity: 0.5 }));
        mm.position.z = 0.1; grp.add(mm);
      }
      q(0.12, 0.05, 0, -g.membrane_nm / 2000, 0.2, new THREE.Color(...rgb).getHex());
      return { center: [0, 0], insetUnitsPerWorld: 1 / S, note: `L ${nice(g.L_um)} µm, ROC ${nice(roc)} µm, ${g.stable ? 'w0 ' + nice(g.waist_um, 3) + ' µm' : 'unstable'}` };
    },
  });
  const gv = hud.group('view');
  hud.button(gv, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the gap magnifier' });

  const cutMax = chipR * 1.1;
  const v = composeView;
  const gapT = [0, L * 0.5, 0];
  const close = v(gapT, [1, 0.16, 0.38], 10);
  const over = v([0, 3.5, 0], [0.62, 0.36, 0.9], 52);
  const bm = (name, vw, cut, magOn = false) => ({ name, pos: vw.pos, target: vw.target, cut, mag: magOn });
  return {
    root,
    defaultView: over,
    bookmarks: [
      bm('Fibre over chip', over, cutMax),
      bm('Cavity cut-away', close, 0),
      bm('Gap, side on', v(gapT, [1, 0.1, 0.0001], 7), 0),
      bm('Gap magnifier', v(gapT, [1, 0.28, 0.5], 14), 0, true),
    ],
    cut: { min: -Rf, max: cutMax },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
