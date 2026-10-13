// Edge-emitting ridge waveguide (composed device / structure edge_ridge).  World axes: x across the ridge, y up,
// the ridge runs along -z and the front (cleaved) facet is at z = 0, facing the camera.  The layer stack is
// waveguide.hkust_ridge_stack(lambda) (the stack the mode solve uses [E]), the ridge width and etch depth come from
// the spec, and the guided mode is drawn with its Gaussian-equivalent 1/e^2 intensity radii wx, wy from the
// effective-index solve, all at one true scale.  The front window length, the back stub, the etched shoulder width and
// the substrate thickness are display defaults [A]; the rest of the cavity length L is hidden behind a break.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, outline, box, fatLine, breakMark, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';
import { composeTag, composeDot, composeBeam, composeView } from './micropillar.js';

export function buildRidgeWaveguide(spec, ctx) {
  const { stage, labels } = ctx;
  const g = spec.geometry;
  const clip = stage.clipPlane;
  const tg = (k) => composeTag(spec, k);
  const fl = g.front_len_um, bl = g.back_len_um, gapUm = 1.6;
  const S = 12 / (fl + gapUm + bl);                       // world units per um
  const nm = (v) => (v / 1000) * S;
  const um = (v) => v * S;
  const W = nm(g.ridge_width_nm), pad = nm(g.shoulder_nm);
  const Lf = um(fl), Lb = um(bl), Lg = um(gapUm);
  const floorY = nm(g.etch_floor_nm);
  const lam = g.lambda_nm;
  const { rgb, falseColour, colourNote, colourSuffix } = lambdaToRGB(lam);
  const root = new THREE.Group();

  // --- the stack, bottom -> top, over the front window and the back stub
  let y = 0;
  const yDot = nm(g.dot_height_nm);
  for (const Ly of g.stack) {
    const t = nm(Ly.t_nm), y0 = y, y1 = y + t;
    y = y1;
    const active = !!Ly.is_dot;
    const role = active ? 'active' : (Ly.name.includes('cladding') ? 'cladding' : 'core');
    const mat = layerMaterial(role, { clip, emissive: active ? 0x8a5a12 : 0x000000, emissiveIntensity: active ? 0.9 : 0 });
    const info = { title: Ly.name.replace(/_/g, ' '), tag: 'E', rows: [['material', Ly.material || 'n/a'], ['thickness', `${nice(Ly.t_nm, 4)} nm`], ['n at λ', nice(Ly.n, 4)]], note: 'waveguide.hkust_ridge_stack(λ), the stack of the mode solve [E]' };
    const segs = [[-Lf / 2, Lf], [-(Lf + Lg) - Lb / 2, Lb]];
    for (const [zc, len] of segs) {
      const r = box(W, t, len, mat, [0, (y0 + y1) / 2, zc]);
      if (!active) outline(r, { clip, opacity: role === 'cladding' ? 0.28 : 0.18 });
      stage.addPickable(r, info);
      root.add(r);
      if (y0 < floorY - 1e-9) {                                    // etched shoulders: the layer below the etch floor continues sideways
        const top = Math.min(y1, floorY), th = top - y0;
        for (const s of [-1, 1]) {
          const sh = box(pad, th, len, mat, [s * (W / 2 + pad / 2), (y0 + top) / 2, zc]);
          outline(sh, { clip, opacity: 0.22 });
          stage.addPickable(sh, { ...info, rows: [...info.rows, ['etched shoulder', `${nice(th * 1000 / S, 4)} nm remain`]] });
          root.add(sh);
        }
      }
    }
  }
  const H = y;

  // --- substrate (thickness not modelled) and the break between front window and back stub
  const X = W / 2 + pad;
  const subT = um(g.substrate_um);
  const zMid = -(Lf + Lg + Lb) / 2;
  const sub = box(2 * X, subT, Lf + Lg + Lb, layerMaterial('substrate', { clip }), [0, -subT / 2, zMid]);
  stage.addPickable(sub, { title: 'substrate', tag: 'A', rows: [['thickness', 'not modelled']], note: 'drawn as a broken plate' });
  root.add(sub);
  root.add(breakMark([-X, -subT, 0.01], [X, -subT, 0.01], 0.06, [0, 1, 0], 9));
  for (const zb of [-Lf, -(Lf + Lg)]) {
    root.add(breakMark([-X, -subT - 0.005, zb], [X, -subT - 0.005, zb], 0.06, [0, 1, 0], 9));
    root.add(breakMark([-W / 2, H + 0.005, zb], [W / 2, H + 0.005, zb], 0.05, [0, 1, 0], 6));
  }

  // --- dot layer line across the front facet, the guided mode along the ridge and at the facet
  const wx = um(g.wx_um), wy = um(g.wy_um);
  root.add(fatLine([[-W / 2, yDot, 0.012], [W / 2, yDot, 0.012]], { color: 0xffc457, width: 2.2 }));
  const spot = new THREE.Mesh(new THREE.CircleGeometry(1, 48), glowMaterial(rgb, { opacity: 0.32 }));
  spot.scale.set(wx, wy, 1);
  spot.position.set(0, yDot, 0.014);
  stage.markBloom(spot, 1.0);
  root.add(spot);
  const tube = new THREE.Mesh(new THREE.CylinderGeometry(1, 1, Lf, 48, 1, true), glowMaterial(rgb, { opacity: 0.09 }));
  tube.rotation.x = Math.PI / 2;
  tube.scale.set(wx, 1, wy);
  tube.position.set(0, yDot, -Lf / 2);
  stage.markBloom(tube);
  root.add(tube);
  const dotZ = -Math.min(Lf * 0.55, Lf - 0.4);
  composeDot(stage, root, rgb, [0, yDot, dotZ], Math.max(W * 0.05, 0.035),
    { title: 'quantum dot', tag: tg('dot_position'), rows: [['emission', `${nice(lam, 4)} nm${colourSuffix}`], ['position along the ridge', `${nice(g.dot_position, 3)} of L`]], note: `dot_position is a fraction of the ${nice(g.L_um)} µm cavity; drawn inside the front window` });
  composeBeam(stage, root, [0, yDot, 0.03], [0, 0, 1], 5, Math.max(Math.max(wx, wy) * 1.1, 0.2));

  // --- back facet: coated (R_back >= 0) or bare cleaved
  const yb = -(Lf + Lg + Lb);
  if (g.R_back >= 0) {
    const coat = box(W, H, 0.03, new THREE.MeshPhysicalMaterial({ color: 0xa8832f, metalness: 0.7, roughness: 0.5, clippingPlanes: [clip] }), [0, H / 2, yb - 0.015]);
    stage.addPickable(coat, { title: 'back-facet coating', tag: tg('R_back'), rows: [['reflectivity R_back', nice(g.R_back, 3)]] });
    root.add(coat);
  }

  // --- annotations
  const xr = W / 2 + pad + 0.25;
  labels.add(root, [W / 2 + 0.1, H + 0.15, -Lf * 0.3], `ridge ${nice(g.ridge_width_nm, 4)} nm wide`, { tag: tg('ridge_width_nm'), priority: 85, sub: `etched ${nice(g.etch_depth_nm, 4)} nm from the top` });
  labels.add(root, [-xr, yDot, -Lf * 0.1], `dot layer · ${nice(lam, 4)} nm`, { tag: 'E', anchor: 'right', priority: 90, sub: falseColour ? `${colourNote}; drawn as a line` : 'a few nm thick: drawn as a line' });
  labels.add(root, [xr, H * 0.5, -Lf * 0.65], `stack ${nice(g.total_nm / 1000, 3)} µm`, { tag: 'E', priority: 60, sub: 'AlInP / AlGaInP core, dot in a well' });
  labels.add(root, [wx * 0.5, yDot + wy + 0.12, 0.05], `mode ${nice(2 * g.wx_um, 3)} × ${nice(2 * g.wy_um, 3)} µm`, { tag: 'DR', priority: 80, sub: `1/e² intensity diameters; n_eff ${nice(g.n_eff, 4)}` });
  labels.add(root, [0, H + 0.45, -(Lf + Lg / 2)], `break: cavity L = ${nice(g.L_um)} µm`, { tag: tg('L_um'), anchor: 'center', priority: 70, sub: `front ${nice(g.front_len_um)} µm and back stub ${nice(g.back_len_um)} µm drawn` });
  labels.add(root, [xr, 0, yb * 0.5 - 0.0], g.R_back >= 0 ? `back facet R ${nice(g.R_back, 3)}` : 'back facet: bare cleave', { tag: tg('R_back'), priority: 50, sub: g.R_back >= 0 ? 'coated back facet' : 'R_back < 0 in the spec' });
  labels.add(root, [-xr, -subT / 2, -Lf * 0.5], 'substrate', { tag: 'A', anchor: 'right', priority: 25, sub: 'thickness not modelled' });
  labels.add(root, [0, yDot - wy - 0.25, 0.9], `front facet T ${nice(g.T_facet, 3)}`, { tag: 'DR', anchor: 'center', priority: 55, sub: 'power transmission of the cleaved facet' });

  const tc = [0, H * 0.5, zMid * 0.55];
  const main = composeView(tc, [0.7, 0.42, 0.9], 15);
  const bm = (name, vw, cut) => ({ name, pos: vw.pos, target: vw.target, cut, mag: false });
  return {
    root,
    defaultView: main,
    bookmarks: [
      bm('Ridge', main, X * 1.1),
      bm('Side cut-away', composeView([0, H * 0.5, -Lf * 0.5], [1, 0.05, 0.12], 11), 0),
      bm('Facet', composeView([0, yDot, 0], [0.32, 0.18, 1], Math.max(H * 3.6, 3.2)), X * 1.1),
      bm('Top', composeView([0, H, zMid * 0.5], [0.12, 1, 0.3], 14), X * 1.1),
    ],
    cut: { min: -X, max: X * 1.1 },
    unitToUm: 1, worldPerUm: S,
    readoutKeys: spec.readout_keys,
  };
}
