// Photonic-crystal waveguide (composed device / structure phc_waveguide): a
// suspended membrane with a triangular hole lattice and one missing row (W1),
// the dot in the row, and the guided mode as a glow tube whose radius is the
// equal-area circle of A_eff (spec value).  n_g, A_eff and the coupler /
// propagation efficiencies are the registry numbers; the lattice constant, hole
// radius, slab thickness and row / hole counts are display defaults [A] (the
// registry has no lattice).  The beam leaves along the row, the active path.
//
// Also exports phcLattice() / phcSlab(), the hole-lattice slab shared with phcCavity.js.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, outline, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeView } from './micropillar.js';

const SQ3 = Math.sqrt(3);

// Triangular-lattice hole centres (nm): rows along x, row spacing a sqrt3 / 2.
export function phcLattice(a, nx, nrEach, skip) {
  const out = [];
  for (let j = -nrEach; j <= nrEach; j++) {
    for (let i = 0; i < nx; i++) {
      const x = (i - (nx - 1) / 2) * a + (Math.abs(j) % 2 ? a / 2 : 0);
      const z = j * a * SQ3 / 2;
      if (skip(i, j, x, z)) continue;
      out.push([x, z, i, j]);
    }
  }
  return out;
}

// Slab with real through-holes, centred on y = 0, thickness h; coordinates in
// nm scaled by k (world per nm).  halfX / halfZ are the slab half extents (nm).
export function phcSlab(stage, parent, { holes, rNm, hNm, halfX, halfZ, k, info }) {
  const sh = new THREE.Shape();
  sh.moveTo(-halfX * k, -halfZ * k); sh.lineTo(halfX * k, -halfZ * k); sh.lineTo(halfX * k, halfZ * k); sh.lineTo(-halfX * k, halfZ * k); sh.closePath();
  for (const [x, z] of holes) {
    const p = new THREE.Path();
    p.absarc(x * k, -z * k, rNm * k, 0, Math.PI * 2, true);
    sh.holes.push(p);
  }
  const geo = new THREE.ExtrudeGeometry(sh, { depth: hNm * k, bevelEnabled: false, curveSegments: 24 });
  geo.rotateX(-Math.PI / 2);
  geo.translate(0, -hNm * k / 2, 0);
  const mat = layerMaterial('core', { clip: stage.clipPlane, opacity: 1 });
  mat.color.setHex(0x86a3c2);
  const slab = new THREE.Mesh(geo, mat);
  outline(slab, { clip: stage.clipPlane, opacity: 0.4 });
  stage.addPickable(slab, info);
  parent.add(slab);
  return slab;
}

export function buildPhcWaveguide(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const tg = (k) => composeTag(spec, k);
  const a = g.a_nm, nx = g.n_cols_z, nr = g.n_rows;
  const lam = g.lambda_nm;
  const { rgb, falseColour } = lambdaToRGB(lam);
  const k = 12 / (nx * a);                       // world per nm (uniform)
  const S = k * 1000;                              // world per um
  const root = new THREE.Group();
  // Rows run along x (the guided mode leaves toward -x), so the cut-away plane x = c
  // is a section perpendicular to the waveguide.
  const holes = phcLattice(a, nx, nr, (i, j) => j === 0);
  const halfX = nx * a / 2 + a * 0.6, halfZ = (nr + 0.5) * a * SQ3 / 2 + a * 0.35;
  phcSlab(stage, root, {
    holes, rNm: g.r_nm, hNm: g.h_nm, halfX, halfZ, k,
    info: { title: 'PhC membrane', tag: 'A', rows: [['lattice a', `${nice(a, 4)} nm`], ['hole r', `${nice(g.r_nm, 3)} nm`], ['thickness', `${nice(g.h_nm, 3)} nm`], ['missing rows', `${g.missing_rows} (W1)`]], note: 'lattice constant, hole radius and thickness are display defaults [A]' },
  });

  // guided mode: equal-area circle of A_eff, along the row (x axis)
  const rm = g.mode_radius_nm * k;
  const Lw = 2 * halfX * k, hz = halfZ * k, hh = g.h_nm * k;
  const mk = (len, x0, op) => {
    const m = new THREE.Mesh(new THREE.CylinderGeometry(rm, rm, len, 32, 1, true), glowMaterial(rgb, { opacity: op }));
    m.rotation.z = Math.PI / 2;
    m.position.set(x0, 0, 0);
    m.material.clippingPlanes = [stage.clipPlane];
    stage.markBloom(m);
    return m;
  };
  const out = mk(Lw / 2, -Lw / 4, 0.2);            // output side (-x), toward the coupler
  const back = mk(Lw / 2, Lw / 4, 0.06);
  stage.addPickable(out, { title: 'guided mode', tag: 'E', rows: [['n_g', nice(g.n_g, 3)], ['A_eff', `${nice(g.A_eff_um2, 3)} µm²`], ['tube radius', `${nice(g.mode_radius_nm, 3)} nm`], ['η_coupler · η_prop', `${nice(g.eta_coupler, 3)} · ${nice(g.eta_prop, 3)}`]], note: 'tube = equal-area circle of A_eff; not a field profile' });
  root.add(out, back);
  composeDot(stage, root, rgb, [0, 0, 0], Math.max(g.dot_r_nm * k, 0.04),
    { title: 'quantum dot', tag: 'A', rows: [['emission', `${nice(lam, 4)} nm${falseColour ? ' (NIR, false colour)' : ''}`], ['position', 'in the row, mid-membrane (display)']], note: 'dot size and position are not model inputs' });
  composeBeam(stage, root, [-Lw / 2, 0, 0], [-1, 0, 0], 5, Math.max(rm * 1.6, 0.2));

  // annotations
  labels.add(root, [Lw * 0.12, hh / 2 + 0.05, hz + 0.1], 'W1 · one missing row', { tag: 'A', priority: 70, sub: 'row along the waveguide axis' });
  labels.add(root, [Lw * 0.3, hh / 2, hz + 0.25], `lattice a ${nice(a, 4)} nm · hole r ${nice(g.r_nm, 3)} nm`, { tag: 'A', priority: 85, sub: 'display defaults' });
  labels.add(root, [Lw * 0.3, 0, -hz - 0.2], `membrane ${nice(g.h_nm, 3)} nm`, { tag: 'A', anchor: 'right', priority: 55 });
  labels.add(root, [-Lw * 0.28, 0.7, 0.4], `guided mode · n_g ${nice(g.n_g, 3)}`, { tag: tg('n_g'), priority: 90, sub: 'slow light: group index from the spec' });
  labels.add(root, [-Lw * 0.5, 0.45, -0.4], `A_eff ${nice(g.A_eff_um2, 3)} µm² · tube r ${nice(g.mode_radius_nm, 3)} nm`, { tag: tg('A_eff_um2'), anchor: 'right', priority: 80 });
  labels.add(root, [0.1, -hh * 0.5 - 0.15, -0.2], `dot · ${nice(lam, 4)} nm`, { tag: 'A', anchor: 'right', priority: 88, sub: 'in the missing row, mid-membrane (display)' });
  labels.add(root, [Lw / 2, 0, 0.3], `λ ${nice(lam, 4)} nm`, { tag: 'DR', priority: 50, sub: falseColour ? 'NIR, false colour' : 'material emission wavelength' });
  labels.add(root, [-Lw / 2 - 0.2, 0.2, 0.35], `coupler η ${nice(g.eta_coupler, 3)} · propagation η ${nice(g.eta_prop, 3)}`, { tag: tg('eta_coupler'), anchor: 'right', priority: 60 });

  const D = Lw * 1.9;
  const v = composeView;
  const bm = (name, vw, cut) => ({ name, pos: vw.pos, target: vw.target, cut, mag: false });
  const lat = v([-Lw * 0.05, 0, 0], [-0.35, 0.42, 0.6], D);
  const cutMax = halfX * k * 1.1;
  return {
    root,
    defaultView: lat,
    bookmarks: [
      bm('Lattice', lat, cutMax),
      bm('Cross-section', v([0, 0, 0], [1, 0.12, 0.3], 6.5), 0),
      bm('Top', v([-Lw * 0.05, 0, 0], [-0.02, 1, 0.03], D * 0.8), cutMax),
      bm('Along the row', v([-Lw * 0.3, 0, 0], [1, 0.16, 0.12], Lw * 0.9), cutMax),
    ],
    cut: { min: -halfX * k, max: cutMax },
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
