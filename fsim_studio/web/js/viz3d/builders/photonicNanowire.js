// Photonic nanowire / trumpet (composed device / structure photonic_nanowire).  True scale in every axis for what
// the registry holds: the wire diameter d_nm, the taper's top mode-field radius (trumpet: the top flares until its
// radius reaches it; needle: no flare) and the bottom mirror (its power reflectivity R and, through
// mirror_phase_ok, whether the dot sits on a standing-wave antinode).  The wire HEIGHT, the flare half-angle, the
// mirror size and the dot height are display defaults [A] (descriptor geometry_src says so).
import * as THREE from 'three';
import { layerMaterial, glowMaterial, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeRingBreak, composeSubstrate, composeLathe, composeView } from './micropillar.js';

export function buildPhotonicNanowire(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const tm = g.mirror_um;
  const H_um = g.length_um + tm;
  const S = 12 / H_um;                                  // world units per um
  const um = (v) => v * S;
  const nm = (v) => (v / 1000) * S;
  const lam = g.lambda_nm;
  const { rgb, falseColour } = lambdaToRGB(lam);
  const root = new THREE.Group();

  const r0 = nm(g.d_nm / 2), rTop = nm(g.top_d_nm / 2);
  const yM = um(tm), yTop = um(H_um), yFl0 = yTop - um(g.flare_len_um);
  const mirR = um(g.mirror_halfwidth_um);

  // --- substrate + bottom mirror (a coated plate under the wire)
  const subT = um(0.7);
  root.add(composeSubstrate(stage, mirR, subT, 0, { info: { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']], note: 'drawn as a broken plate' } }));
  root.add(composeRingBreak(mirR, -subT / 2, 0.05, 28));
  const mirMat = new THREE.MeshPhysicalMaterial({ color: 0xa8832f, metalness: 0.7, roughness: 0.5, clippingPlanes: [clip] });
  const mir = new THREE.Mesh(new THREE.CylinderGeometry(mirR, mirR, yM, 96), mirMat);
  mir.position.y = yM / 2;
  stage.addPickable(mir, {
    title: 'bottom mirror', tag: tg('mirror_R'),
    rows: [['power reflectivity R', nice(g.mirror_R, 4)], ['field |r| = √R', nice(g.mirror_field_r, 4)], ['dot on antinode', g.mirror_phase_ok ? 'yes (coherent boost)' : 'no / unknown (incoherent return)'], ['rate boost 1 + β|r|', nice(g.rate_boost, 4)]],
    note: 'thickness and plate size are display defaults [A]',
  });
  root.add(mir);

  // --- the wire: needle (constant radius) or trumpet (flare over the top flare_len_um)
  const prof = [[r0, yM], [r0, yFl0]];
  if (g.trumpet && yTop > yFl0 + 1e-6) {
    for (let i = 1; i <= 24; i++) {                       // a smooth flare: radius grows as t^2 (adiabatic-looking), ends at rTop
      const t = i / 24;
      prof.push([r0 + (rTop - r0) * t * t, yFl0 + (yTop - yFl0) * t]);
    }
  } else {
    prof.push([r0, yTop]);
  }
  const wireMat = layerMaterial('core', { clip, opacity: 1 });
  wireMat.color.setHex(0x86a3c2);
  const wire = composeLathe(prof, wireMat, 72);
  const cap = new THREE.Mesh(new THREE.CircleGeometry(prof[prof.length - 1][0], 72), wireMat);
  cap.rotation.x = -Math.PI / 2; cap.position.y = yTop;
  wire.add(cap);
  stage.addPickable(wire, {
    title: g.trumpet ? 'photonic trumpet' : 'photonic nanowire', tag: tg('d_nm'),
    rows: [['diameter', `${nice(g.d_nm, 4)} nm`], ['n', nice(g.n_wire, 4)], ['V number', nice(g.V_number, 3)], ['β (HE11)', nice(g.beta_he11, 3)], ['top diameter', `${nice(g.top_d_nm, 4)} nm`], ['height', `${nice(g.length_um)} µm (display)`]],
    note: 'diameter is a registry parameter; the height and the flare half-angle are display defaults [A]',
  });
  root.add(wire);

  // --- guided mode: a slim glow inside the wire, and the top-facet mode spot (mode-field radius)
  const modeProf = prof.map(([r, y]) => [Math.max(r * 0.8, 0.004), y]);
  const mode = composeLathe(modeProf, glowMaterial(rgb, { opacity: 0.2 }), 36);
  stage.markBloom(mode);
  root.add(mode);
  const spotR = Math.max(nm(g.mode_radius_nm), rTop * 0.5);
  const spot = new THREE.Mesh(new THREE.CircleGeometry(spotR, 48), glowMaterial(rgb, { opacity: 0.55 }));
  spot.rotation.x = -Math.PI / 2; spot.position.y = yTop + 0.012;
  stage.markBloom(spot, 2.0);
  root.add(spot);

  // --- dot
  const yd = yM + um(g.dot_height_um);
  composeDot(stage, root, rgb, [0, yd, 0], Math.max(r0 * 0.5, 0.05),
    { title: 'quantum dot', tag: tg('dot_height_um'), rows: [['emission', `${nice(lam, 4)} nm${falseColour ? ' (NIR, false colour)' : ''}`], ['height above the mirror', `${nice(g.dot_height_um, 3)} µm (display)`]], note: 'dot height is not a model input' });
  composeBeam(stage, root, [0, yTop + 0.03, 0], [0, 1, 0], 5, Math.max(rTop * 0.7, 0.16));

  // --- annotations
  const lx = Math.max(rTop, r0) + 0.3;
  labels.add(root, [lx, (yM + yFl0) / 2, 0], `wire Ø ${nice(g.d_nm, 4)} nm · height ${nice(g.length_um)} µm`, { tag: tg('d_nm'), priority: 85, sub: 'diameter from the spec; height is a display default' });
  labels.add(root, [lx, yTop - 0.15, 0], g.trumpet ? `trumpet top Ø ${nice(g.top_d_nm, 4)} nm` : 'needle top (no taper)', { tag: g.trumpet ? tg('taper_mfr_nm') : 'A', priority: 70, sub: g.trumpet ? 'radius = taper mode-field radius (flare angle is a display default)' : 'taper_mfr_nm = 0' });
  labels.add(root, [-lx, yd, 0], `dot · ${nice(lam, 4)} nm`, { tag: tg('dot_height_um'), anchor: 'right', priority: 90, sub: g.mirror_phase_ok ? 'antinode: λ/(2 n) above the mirror' : 'height not a model input' });
  labels.add(root, [mirR * 0.8, yM / 2, mirR * 0.5], `mirror R ${nice(g.mirror_R, 3)}`, { tag: tg('mirror_R'), priority: 80, sub: g.mirror_phase_ok ? `|r| ${nice(g.mirror_field_r, 3)}, coherent: rate × ${nice(g.rate_boost, 3)}` : 'phase unknown: incoherent return' });
  labels.add(root, [0, yTop + 0.55, spotR * 0.8], `β ${nice(g.beta_he11, 3)} · up ${nice(g.upward_fraction, 3)}`, { tag: 'DR', anchor: 'center', priority: 75, sub: `top mode-field radius ${nice(g.mode_radius_nm, 3)} nm` });
  labels.add(root, [-mirR * 0.7, -subT / 2, mirR * 0.5], 'substrate', { tag: 'A', anchor: 'right', priority: 25, sub: 'thickness not modelled' });
  labels.add(root, [0, yTop + 1.9, 0], `λ ${nice(lam, 4)} nm`, { tag: 'DR', anchor: 'center', priority: 50, sub: falseColour ? 'NIR, false colour' : 'material emission wavelength' });

  const tc = [0, yTop * 0.5, 0];
  const main = composeView(tc, [0.62, 0.16, 0.86], yTop * 2.1);
  const bm = (name, vw, cut) => ({ name, pos: vw.pos, target: vw.target, cut, mag: false });
  return {
    root,
    defaultView: main,
    bookmarks: [
      bm('Wire', main, mirR * 1.1),
      bm('Side cut-away', composeView([0, yTop * 0.5, 0], [1, 0.03, 0.05], yTop * 1.45), 0),
      bm('Dot and mirror', composeView([0, (yM + yd) / 2, 0], [0.55, 0.2, 0.85], Math.max(um(g.dot_height_um) * 3.2, 2.2)), 0),
      bm('Top facet / beam', composeView([0, yTop, 0], [0.45, 0.5, 0.75], Math.max(rTop * 14, 3)), mirR * 1.1),
    ],
    cut: { min: -mirR, max: mirR * 1.1 },
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
