// Lumped cavity (composed device / structure generic_cavity).  This model has NO geometry: it takes the Purcell
// factor F_P, the linewidth kappa (or Q), the detuning, beta and the outcoupling eta_out.  The scene is therefore an
// honest abstract schematic, labelled as such: two mirror plates, a single mode glow and a dot, all NOT to scale
// (no scale bar, no floor).  The only numbers drawn from the model are printed on the callouts.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeView } from './micropillar.js';

export function buildLumpedCavity(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();
  const R = g.plate_um, gap = g.gap_um * 2.2;           // schematic units, not a length scale
  const tp = 0.18;

  // two mirror plates (glass-grey discs) and a dashed frame marking the single lumped mode
  const plateMat = () => layerMaterial('mirror', { clip, opacity: 0.55 });
  for (const s of [-1, 1]) {
    const p = new THREE.Mesh(new THREE.CylinderGeometry(R, R, tp, 72), plateMat());
    p.position.y = s * gap / 2;
    stage.addPickable(p, { title: s > 0 ? 'output mirror (schematic)' : 'back mirror (schematic)', tag: 'DR', rows: [['eta_out', nice(g.eta_out, 3)], ['κ', g.kappa_ueV == null ? 'n/a' : `${nice(g.kappa_ueV, 4)} µeV`]], note: 'not a geometry: the model is lumped' });
    root.add(p);
  }
  const frame = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(R * 1.5, gap, R * 1.5)),
    new THREE.LineDashedMaterial({ color: 0xd8d6cf, dashSize: 0.18, gapSize: 0.14, transparent: true, opacity: 0.4 }));
  frame.computeLineDistances();
  frame.material.clippingPlanes = [clip];
  root.add(frame);

  // the single mode: one glow ellipsoid, plus the dot at its centre
  const m = new THREE.Mesh(new THREE.SphereGeometry(1, 40, 24), glowMaterial(rgb, { opacity: 0.13 }));
  m.scale.set(R * 0.42, gap * 0.36, R * 0.42);
  stage.markBloom(m);
  root.add(m);
  composeDot(stage, root, rgb, [0, 0, 0], 0.09,
    { title: 'quantum dot', tag: 'DR', rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['F_P', nice(g.F_P, 4)], ['detuning', `${nice(g.detuning_ueV, 3)} µeV`]], note: 'schematic position: the model has no geometry' });
  const topY = gap / 2 + tp / 2;
  composeBeam(stage, root, [0, topY + 0.02, 0], [0, 1, 0], 3.4, 0.32);

  // annotations
  labels.add(root, [R + 0.35, gap / 2, 0], 'output mirror · κ', { tag: 'DR', priority: 75, sub: g.kappa_ueV == null ? 'κ n/a' : `κ = ${nice(g.kappa_ueV, 4)} µeV (schematic plate)` });
  labels.add(root, [R + 0.35, -gap / 2, 0], 'back mirror', { tag: 'DR', priority: 40, sub: 'schematic plate, not a layer stack' });
  labels.add(root, [-R * 0.55, 0.0, R * 0.2], `dot · ${nice(lam, 4)} nm`, { tag: 'DR', anchor: 'right', priority: 90, sub: `detuning ${nice(g.detuning_ueV, 3)} µeV` });
  labels.add(root, [0, topY + 1.7, 0], `F_P ${nice(g.F_P, 4)} · η_out ${nice(g.eta_out, 3)}`, { tag: tg('F_P'), anchor: 'center', priority: 95, sub: `Q_eff ${g.Q_eff == null ? 'n/a' : nice(g.Q_eff, 4)} · β0 ${nice(g.beta0, 3)}` });
  labels.add(root, [-R * 0.9, -gap / 2 - 0.6, R * 0.8], 'NO GEOMETRY IN THIS MODEL', { tag: 'DR', anchor: 'right', priority: 99, sub: 'a lumped cavity: only F_P, κ, detuning, β, η_out; no mode volume' });
  labels.add(root, [R * 0.6, gap * 0.62 - 0.1, R * 0.5], 'mode (schematic)', { tag: 'A', priority: 45, sub: 'gap and plate size are display values, not a scale' });

  const tc = [0, 0.4, 0];
  const main = composeView(tc, [0.62, 0.22, 0.86], 11);
  const bm = (name, vw, cut) => ({ name, pos: vw.pos, target: vw.target, cut, mag: false });
  return {
    root,
    defaultView: main,
    bookmarks: [
      bm('Schematic', main, R * 1.1),
      bm('Side cut-away', composeView([0, 0, 0], [1, 0.04, 0.05], 10), 0),
      bm('Dot and mode', composeView([0, 0, 0], [0.5, 0.15, 0.85], 4.8), R * 1.1),
      bm('Output', composeView([0, topY, 0], [0.45, 0.5, 0.75], 7), R * 1.1),
    ],
    cut: { min: -R, max: R * 1.1 },
    noFloor: true, noScaleBar: true,
    unitToUm: 1, worldPerUm: 1,
    readoutKeys: spec.readout_keys,
  };
}
