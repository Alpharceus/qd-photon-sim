// Photonic-crystal cavity (composed device / structure phc_cavity): a
// suspended membrane with a triangular hole lattice and an L3 (three missing
// holes, end holes shifted outward) or H1 (one missing hole) defect.  Q,
// V (lambda/n)^3 and g are the registry numbers; the glow ellipsoid has the
// volume V (lambda/n)^3 with n and the aspect a display default.  The lattice
// constant, hole radius, thickness, counts and the end shift are display
// defaults [A].  The dot is at the defect centre; the beam leaves vertically.
import * as THREE from 'three';
import { glowMaterial, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam } from './micropillar.js';
import { phcLattice, phcSlab } from './phcWaveguide.js';

const SQ3 = Math.sqrt(3);

export function buildPhcCavity(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const tg = (k) => composeTag(spec, k);
  const a = g.a_nm, nx = g.n_cols, nr = g.n_rows;
  const lam = g.lambda_nm;
  const { rgb, falseColour } = lambdaToRGB(lam);
  const k = 12 / (nx * a);                          // world per nm
  const S = k * 1000;
  const root = new THREE.Group();
  const L3 = g.cell === 'L3';
  const mid = (nx - 1) / 2;
  const skip = (i, j) => j === 0 && (L3 ? Math.abs(i - mid) <= 1 : i === mid);
  const holes = phcLattice(a, nx, nr, skip).map(([x, z, i, j]) => {
    // L3: the two holes next to the defect move outward by 0.15 a
    const shifted = L3 && j === 0 && Math.abs(i - mid) === 2;
    return [shifted ? x + Math.sign(x) * g.end_shift_over_a * a : x, z, i, j];
  });
  const halfX = nx * a / 2 + a * 0.6, halfZ = (nr + 0.5) * a * SQ3 / 2 + a * 0.35;
  phcSlab(stage, root, {
    holes, rNm: g.r_nm, hNm: g.h_nm, halfX, halfZ, k,
    info: { title: `PhC membrane (${g.cell})`, tag: 'A', rows: [['lattice a', `${nice(a, 4)} nm`], ['hole r', `${nice(g.r_nm, 3)} nm`], ['thickness', `${nice(g.h_nm, 3)} nm`], ['defect', g.cell === 'L3' ? 'L3, end holes shifted 0.15 a' : 'H1']], note: 'lattice, thickness and defect type are display choices [A]' },
  });

  // mode glow: ellipsoid of volume V (lambda/n)^3 (semi-axes in um from the descriptor)
  const [ax, ay, az] = g.mode_semiaxes_um;
  const mode = new THREE.Mesh(new THREE.SphereGeometry(1, 40, 20), glowMaterial(rgb, { opacity: 0.2 }));
  mode.scale.set(ax * S, az * S, ay * S);
  stage.markBloom(mode);
  stage.addPickable(mode, { title: 'cavity mode (volume)', tag: tg('V_lambda3'), rows: [['V', `${nice(g.V_lambda3, 3)} (λ/n)³`], ['Q', nice(g.Q, 4)], ['g', g.g_ueV > 0 ? `${nice(g.g_ueV, 3)} µeV` : 'unused'], ['semi-axes', `${[ax, ay, az].map((v) => nice(v * 1000, 3)).join(' / ')} nm`]], note: 'ellipsoid with volume V (λ/n)³; n and the aspect are display defaults [A]; not a field profile' });
  root.add(mode);
  composeDot(stage, root, rgb, [0, 0, 0], Math.max(g.dot_r_nm * k, 0.04),
    { title: 'quantum dot', tag: 'A', rows: [['emission', `${nice(lam, 4)} nm${falseColour ? ' (NIR, false colour)' : ''}`], ['position', 'defect centre (display)']], note: 'dot size and position are not model inputs' });
  composeBeam(stage, root, [0, g.h_nm * k / 2 + 0.02, 0], [0, 1, 0], 5, Math.max(az * S * 0.9, 0.2));

  // annotations
  const hw = halfX * k, hz = halfZ * k;
  labels.add(root, [ax * S + 0.1, 0.5, 0], `${g.cell} cavity · Q ${nice(g.Q, 3)}`, { tag: tg('Q'), priority: 90, sub: g.cell === 'L3' ? 'three missing holes, ends shifted' : 'one missing hole' });
  labels.add(root, [-ax * S - 0.1, 0.5, -hz * 0.15], `V ${nice(g.V_lambda3, 3)} (λ/n)³`, { tag: tg('V_lambda3'), anchor: 'right', priority: 85, sub: 'glow ellipsoid has this volume' });
  labels.add(root, [hw * 0.55, 0, hz + 0.15], `lattice a ${nice(a, 4)} nm · r ${nice(g.r_nm, 3)} nm`, { tag: 'A', priority: 70, sub: 'display defaults' });
  labels.add(root, [-hw * 0.7, g.h_nm * k / 2, hz * 0.9], `membrane ${nice(g.h_nm, 3)} nm`, { tag: 'A', anchor: 'right', priority: 50 });
  labels.add(root, [0, 1.9, 0], `λ ${nice(lam, 4)} nm`, { tag: 'DR', anchor: 'center', priority: 60, sub: falseColour ? 'NIR, false colour' : 'material emission wavelength' });

  const mag = makeMagnifier(stage, hud, {
    mag: 1.2, title: 'cavity region (true scale)', unitLabel: 'nm',
    buildInset(scene) {
      const mkc = (x, z, r, c) => { const o = new THREE.Mesh(new THREE.CircleGeometry(r, 32), new THREE.MeshBasicMaterial({ color: c, toneMapped: false })); o.position.set(x, -z, 0.1); scene.add(o); };
      const bg = new THREE.Mesh(new THREE.PlaneGeometry(a * 9, a * 6), new THREE.MeshBasicMaterial({ color: 0x86a3c2, toneMapped: false }));
      scene.add(bg);
      for (const [x, z] of holes) if (Math.abs(x) < a * 4.6 && Math.abs(z) < a * 3) mkc(x, z, g.r_nm, 0x0b0d10);
      mkc(0, 0, g.dot_r_nm * 2, new THREE.Color(...rgb).getHex());
      return { center: [0, 0], insetUnitsPerWorld: 1 / k, note: `top view, a ${nice(a, 4)} nm, r ${nice(g.r_nm, 3)} nm` };
    },
  });
  const gv = hud.group('view');
  hud.button(gv, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the top-view magnifier' });

  const D = hw * 3.4;
  return {
    root,
    defaultView: { pos: [D * 0.45, D * 0.5, D * 0.62], target: [0, 0, 0] },
    bookmarks: [
      { name: 'Membrane', pos: [D * 0.45, D * 0.5, D * 0.62], target: [0, 0, 0], cut: hw * 1.3, mag: false },
      { name: 'Side cut-away', pos: [0.2, 0.7, D * 0.75], target: [0, 0, 0], cut: 0, mag: false },
      { name: 'Top', pos: [0.2, D * 0.95, 0.35], target: [0, 0, 0], cut: hw * 1.3, mag: false },
      { name: 'Defect close-up', pos: [hw * 0.55, hw * 0.55, hw * 0.75], target: [0, 0, 0], cut: hw * 1.3, mag: true },
    ],
    cut: { min: -hw, max: hw * 1.3 },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
