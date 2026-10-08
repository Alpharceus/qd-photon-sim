// Edge-emitting ridge (device / edge_ridge).  World unit = 1 um, true scale
// in every axis.  The front 12 um of the 250 um cavity is drawn, then a
// break, then the back facet.  The facet carries the real guided-mode
// intensity from waveguide.effective_index_ridge.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, coneMaterial, outline, box, fatLine, breakMark, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';

function decodeMode(mode, rgb) {
  const raw = atob(mode.data);
  const n = mode.nx * mode.ny;
  const px = new Uint8Array(n * 4);
  let peak = 0, peakRow = 0;
  for (let i = 0; i < n; i++) {
    const v = raw.charCodeAt(i);
    if (v > peak) { peak = v; peakRow = Math.floor(i / mode.nx); }
    const a = v / 255;
    px[i * 4] = Math.round(255 * Math.min(1, rgb[0] * (0.35 + 0.9 * a) + 0.25 * a * a));
    px[i * 4 + 1] = Math.round(255 * Math.min(1, rgb[1] * (0.35 + 0.9 * a) + 0.18 * a * a));
    px[i * 4 + 2] = Math.round(255 * Math.min(1, rgb[2] * (0.35 + 0.9 * a) + 0.18 * a * a));
    px[i * 4 + 3] = Math.round(255 * Math.pow(a, 0.75));
  }
  const tex = new THREE.DataTexture(px, mode.nx, mode.ny, THREE.RGBAFormat);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.magFilter = THREE.LinearFilter; tex.minFilter = THREE.LinearFilter;
  tex.needsUpdate = true;
  const yPeak = mode.y_um[0] + (peakRow / (mode.ny - 1)) * (mode.y_um[1] - mode.y_um[0]);
  return { tex, yPeak };
}

export function buildRidge(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const um = (nm) => nm / 1000;
  const W = um(g.ridge_width_nm), pad = um(g.shoulder_nm), L = g.front_len_um;
  const X = W / 2 + pad;
  const floorY = um(g.etch_floor_nm);
  const clip = stage.clipPlane;
  const root = new THREE.Group();
  const optical = new THREE.Group();
  const thermal = new THREE.Group();
  root.add(optical, thermal);

  // --- optical stack (waveguide.py), bottom -> top
  let y = 0;
  let dotMesh = null;
  for (const Ly of spec.stacks.optical.layers) {
    const t = um(Ly.t_nm), y0 = y, y1 = y + t;
    y = y1;
    const active = Ly.role === 'active';
    const mat = layerMaterial(Ly.role, { clip, emissive: active ? 0x8a5a12 : 0x000000, emissiveIntensity: active ? 0.9 : 0 });
    const info = {
      title: Ly.name, tag: Ly.tag,
      rows: [['material', Ly.material], ['thickness', `${nice(Ly.t_nm, 4)} nm`], ['n at λ', Ly.n?.toFixed(4)], ['stack', 'optical (waveguide.py)']],
      note: Ly.source,
    };
    const r = box(W, t, L, mat, [0, (y0 + y1) / 2, -L / 2]);
    if (!active) outline(r, { clip, opacity: Ly.role === 'cladding' ? 0.28 : 0.18 });
    stage.addPickable(r, info);
    optical.add(r);
    if (active) dotMesh = r;
    if (y0 < floorY - 1e-9) {
      const top = Math.min(y1, floorY), th = top - y0;
      for (const s of [-1, 1]) {
        const sh = box(pad, th, L, mat, [s * (W / 2 + pad / 2), (y0 + top) / 2, -L / 2]);
        outline(sh, { clip, opacity: 0.22 });
        stage.addPickable(sh, { ...info, rows: [...info.rows, ['etched shoulder', `${nice(th * 1000, 4)} nm remain`]] });
        optical.add(sh);
      }
    }
  }
  const H = y;

  // --- substrate (thickness not modelled)
  const subT = 0.9;
  const sub = box(2 * X, subT, L, layerMaterial('substrate', { clip }), [0, -subT / 2, -L / 2]);
  stage.addPickable(sub, { title: `${g.substrate.material} substrate`, tag: 'A', rows: [['thickness', 'not modelled']], note: g.substrate.note });
  root.add(sub);
  root.add(breakMark([-X, -subT, 0.01], [X, -subT, 0.01], 0.07, [0, 1, 0], 9));

  // --- thermal stack (thermal.py): drawn in place of the optical stack on toggle
  let ty = 0;
  for (const Ly of spec.stacks.thermal.layers) {
    const t = um(Ly.t_nm);
    const m = box(W, t, L, layerMaterial('thermal', { clip, opacity: 0.55 }), [0, ty + t / 2, -L / 2]);
    m.material.color.setHex(ty === 0 ? 0x7a8792 : 0x98a1aa);
    stage.addPickable(m, { title: Ly.name, tag: Ly.tag, rows: [['thickness', `${nice(Ly.t_nm, 4)} nm`], ['k300', `${Ly.k300} W/m/K`], ['stack', 'thermal (thermal.py)']], note: 'thermal.layers entry; listed order' });
    thermal.add(m);
    ty += t;
  }
  const ghost = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(W, H, L)),
    new THREE.LineDashedMaterial({ color: 0xd8d6cf, dashSize: 0.12, gapSize: 0.08, transparent: true, opacity: 0.7 }));
  ghost.position.set(0, H / 2, -L / 2);
  ghost.computeLineDistances();
  thermal.add(ghost);
  thermal.visible = false;

  // --- dots (seeded Poisson sample at the card density), true size
  const dot = g.dot;
  const dr = um(dot.radius_nm), dh = um(dot.height_nm), yd = um(dot.y_nm);
  const pts = g.dots.xz_um;
  const dotGeo = new THREE.SphereGeometry(1, 14, 8, 0, Math.PI * 2, 0, Math.PI / 2);
  const dotMat = new THREE.MeshStandardMaterial({ color: 0xf0c060, emissive: 0xffb440, emissiveIntensity: 1.4, roughness: 0.4 });
  dotMat.clippingPlanes = [clip];
  const dots = new THREE.InstancedMesh(dotGeo, dotMat, Math.max(1, pts.length));
  const M = new THREE.Matrix4();
  pts.forEach(([x, z], i) => { M.compose(new THREE.Vector3(x, yd - dh / 2, z), new THREE.Quaternion(), new THREE.Vector3(dr, dh, dr)); dots.setMatrixAt(i, M); });
  dots.count = pts.length;
  stage.markBloom(dots);
  optical.add(dots);

  // --- aperture ring in the dot plane
  const ap = g.aperture;
  const ring = [];
  for (let i = 0; i <= 64; i++) {
    const a = (i / 64) * Math.PI * 2;
    ring.push([ap.center_um[0] + Math.cos(a) * ap.diameter_um / 2, yd + 0.002, ap.center_um[1] + Math.sin(a) * ap.diameter_um / 2]);
  }
  optical.add(fatLine(ring, { color: 0xe9e8e3, width: 1.5 }));

  // --- facet mode texture (real |E|^2) + guided-mode envelope inside the ridge
  const mode = spec.overlays.mode;
  const lam = spec.overlays.emission.lambda_nm;
  const { rgb } = lambdaToRGB(lam);
  const { tex, yPeak } = decodeMode(mode, rgb);
  const mw = mode.x_um[1] - mode.x_um[0], mh = mode.y_um[1] - mode.y_um[0];
  const facet = new THREE.Mesh(new THREE.PlaneGeometry(mw, mh),
    new THREE.MeshBasicMaterial({ map: tex, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, toneMapped: false }));
  facet.position.set((mode.x_um[0] + mode.x_um[1]) / 2, (mode.y_um[0] + mode.y_um[1]) / 2, 0.004);
  stage.markBloom(facet);
  stage.addPickable(facet, { title: 'guided mode |E|² on the facet', tag: mode.tag, rows: [['wx (1/e²)', `${nice(mode.wx_um)} µm`], ['wy (1/e²)', `${nice(mode.wy_um)} µm`], ['n_eff', nice(mode.n_eff, 5)], ['λ', `${nice(lam, 4)} nm`]], note: mode.label });
  root.add(facet);
  const env = new THREE.Mesh(new THREE.CylinderGeometry(1, 1, L, 48, 1, true), glowMaterial(rgb, { opacity: 0.035 }));
  env.rotation.x = Math.PI / 2;
  env.scale.set(mode.wx_um, 1, mode.wy_um);
  env.position.set(0, yPeak, -L / 2);
  root.add(env);

  // --- NA cone
  const na = spec.overlays.na_cone;
  const coneLen = 1.6;
  const th = (na.half_angle_deg * Math.PI) / 180;
  const cone = new THREE.Mesh(new THREE.ConeGeometry(coneLen * Math.tan(th), coneLen, 64, 1, true), coneMaterial(rgb, 0.07));
  cone.rotation.x = -Math.PI / 2;
  cone.position.set(0, yPeak, coneLen / 2 + 0.01);
  root.add(cone);

  // --- the breadboard beam: the 2D table's one live beam enters the 3D
  // enclosure here, leaving the facet along the guided mode's axis; a glow
  // marks the emission point.  Beam ink, not the emission wavelength colour.
  const BEAM = [1.0, 0.231, 0.184];
  const beamLine = fatLine([[0, yPeak, 0.02], [0, yPeak, 9.0]], { color: 0xff3b2f, width: 4 });
  beamLine.name = 'breadboard-beam';
  root.add(beamLine);
  const glow = new THREE.Mesh(new THREE.SphereGeometry(Math.max(mode.wy_um, 0.12), 24, 16), glowMaterial(BEAM, { opacity: 0.85 }));
  glow.scale.set(mode.wx_um / Math.max(mode.wy_um, 0.12), 1, 0.5);
  glow.position.set(0, yPeak, 0.03);
  glow.name = 'emission-glow';
  stage.markBloom(glow, 2.5);
  root.add(glow);

  // --- break + back facet (L = 250 um is NOT compressed: it is broken)
  const gap = 1.6, seg = 2.2;
  const zb0 = -L, zb1 = -L - gap;
  for (const z of [zb0 - 0.05, zb1 + 0.05]) {
    root.add(breakMark([-W / 2, H + 0.02, z], [W / 2, H + 0.02, z], 0.12, [0, 0, 1], 4));
    for (const s of [-1, 1]) root.add(breakMark([s * W / 2, floorY + 0.02, z], [s * X, floorY + 0.02, z], 0.12, [0, 0, 1], 4));
  }
  const back = new THREE.Group();
  let by = 0;
  for (const Ly of spec.stacks.optical.layers) {
    const t = um(Ly.t_nm);
    const m = layerMaterial(Ly.role, { opacity: Ly.role === 'cladding' ? 0.3 : 0.9 });
    back.add(box(W, t, seg, m, [0, by + t / 2, zb1 - seg / 2]));
    if (by < floorY - 1e-9) {
      const top = Math.min(by + t, floorY);
      for (const s of [-1, 1]) back.add(box(pad, top - by, seg, m, [s * (W / 2 + pad / 2), (by + top) / 2, zb1 - seg / 2]));
    }
    by += t;
  }
  back.add(box(2 * X, subT, seg, layerMaterial('substrate'), [0, -subT / 2, zb1 - seg / 2]));
  const mirror = new THREE.Mesh(new THREE.PlaneGeometry(W * 1.02, H * 1.02),
    new THREE.MeshPhysicalMaterial({ color: 0xcfd6e0, metalness: 1, roughness: 0.12, side: THREE.DoubleSide }));
  mirror.position.set(0, H / 2, zb1 - seg - 0.01);
  back.add(mirror);
  root.add(back);

  // --- annotations (CSS2D, <= 30)
  const lx = X + 0.15;
  let cy = 0;
  for (const Ly of spec.stacks.optical.layers) {
    const t = um(Ly.t_nm);
    const mid = cy + t / 2;
    cy += t;
    if (Ly.role === 'active') {
      labels.add(root, [-W / 2, mid, -10.5], `dot plane · ${Ly.material} · ${nice(Ly.t_nm)} nm`, { priority: 85, tag: Ly.tag, anchor: 'left', sub: `${nice(100 * Ly.t_nm / spec.stacks.optical.total_nm, 2)} % of the stack, true thickness; click for the x500 magnifier`, anchor: 'right' });
      continue;
    }
    const lyY = Ly.name === 'matrix_lower' ? mid - 0.06 : Ly.name === 'matrix_upper' ? mid + 0.06 : mid;
    labels.add(optical, [lx, lyY, 0], `${Ly.name} · ${nice(Ly.t_nm, 4)} nm · n ${Ly.n.toFixed(3)}`, { tag: Ly.tag, sub: Ly.material });
  }
  labels.add(root, [W / 2 + pad * 0.55, floorY + 0.02, -L * 0.42], `etch ${nice(g.etch_depth_nm, 4)} nm · floor ${nice(g.etch_floor_nm, 4)} nm < dot plane ${nice(dot.y_nm, 4)} nm`,
    { priority: 80, tag: 'A', sub: g.etch_through_dot ? 'etch removes dot layer; sidewalls not modelled' : 'etch stops above the dot layer' });
  labels.add(root, [-X, H + 0.32, 0], `guided mode |E|² · wx ${nice(mode.wx_um)} µm · wy ${nice(mode.wy_um)} µm · n_eff ${nice(mode.n_eff, 4)}`, { priority: 90, tag: mode.tag, sub: 'waveguide.effective_index_ridge', anchor: 'left' });
  labels.add(root, [coneLen * Math.tan(th) * 0.72, yPeak + coneLen * Math.tan(th) * 0.72, coneLen], `collection NA ${nice(na.NA)} · ${nice(na.half_angle_deg, 3)}°`, { priority: 60, tag: na.tag });
  labels.add(root, [X, -subT * 0.6, -1.0], `${g.substrate.material} substrate`, { priority: 30, tag: 'A', sub: 'thickness not modelled' });
  labels.add(root, [-W / 2 - 0.1, H, zb1 - seg / 2], `back facet R ${nice(g.R_back)} · L = ${nice(g.L_um)} µm`, { priority: 55, tag: 'A', sub: `drawn true-scale for ${L} µm, then broken; length never compressed`, anchor: 'right' });
  labels.add(optical, [ap.center_um[0] - ap.diameter_um / 2 - 0.05, yd + 0.05, ap.center_um[1]], `aperture Ø ${nice(ap.diameter_um)} µm · expected dots ${nice(ap.expected_dots, 2)}`,
    { priority: 45, tag: spec.scalars?.expected_dots?.tag || 'A', anchor: 'right', sub: ap.placement });
  labels.add(optical, [-W / 2 + 0.05, yd + 0.04, -L + 0.6], `dots ${g.dots.xz_um.length} · ${nice(g.dots.per_um2)} µm⁻²`, { priority: 35, tag: 'A', anchor: 'left', sub: 'seeded sample at the card density; positions illustrative' });
  const thermalNote = labels.add(thermal, [lx, ty / 2, 0], `thermal stack ${nice(spec.stacks.thermal.total_nm / 1000)} µm`, { priority: 70, tag: 'A', sub: `optical stack ${nice(H)} µm (dashed): the two models disagree` });

  // --- controls
  const gStack = hud.group('stack');
  hud.toggle(gStack, [['optical', 'optical (waveguide.py)'], ['thermal', 'thermal (thermal.py)']], 'optical', (v) => {
    optical.visible = v === 'optical'; thermal.visible = v === 'thermal';
    facet.visible = env.visible = v === 'optical';
    stage.requestRender();
  });
  void thermalNote;

  // --- magnifier inset: true-scale section through one dot
  const mag = makeMagnifier(stage, hud, {
    mag: 500, title: 'section through one dot (true scale)',
    buildInset(scene) {
      const h = dot.height_nm, r = dot.radius_nm, wl = dot.wl_thickness_nm || 0, core = 148;
      const mat = new THREE.MeshStandardMaterial({ color: 0x5d7fa3, roughness: 0.6 });
      const q = (w, hh, x, yy, z, m) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(w, hh), m); o.position.set(x, yy, z); scene.add(o); return o; };
      q(2000, core, 0, -h / 2 - core / 2, 0, mat);
      q(2000, core, 0, h / 2 + core / 2, 0, mat);
      q(2000, h, 0, 0, 0, mat);
      const wlm = new THREE.MeshBasicMaterial({ color: 0xd9a845, toneMapped: false });
      if (wl > 0) q(2000, wl, 0, -h / 2 + wl / 2, 0.1, wlm);
      // spherical-cap profile, base radius r, height h (drawn section)
      const shape = new THREE.Shape();
      const R = (r * r + h * h) / (2 * h);
      shape.moveTo(-r, -h / 2);
      for (let i = 0; i <= 48; i++) {
        const xx = -r + (2 * r * i) / 48;
        const yy = Math.sqrt(Math.max(0, R * R - xx * xx)) - (R - h);
        shape.lineTo(xx, -h / 2 + yy);
      }
      shape.lineTo(r, -h / 2);
      const lens = new THREE.Mesh(new THREE.ShapeGeometry(shape), new THREE.MeshBasicMaterial({ color: 0xffc35c, toneMapped: false }));
      lens.position.z = 0.2;
      scene.add(lens);
      const box2 = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.PlaneGeometry(2000, h)), new THREE.LineDashedMaterial({ color: 0xe9e8e3, dashSize: 1.2, gapSize: 0.8 }));
      box2.position.z = 0.3; box2.computeLineDistances();
      scene.add(box2);
      return { center: [0, 0], insetUnitsPerWorld: 1000, note: `${dot.material} dot r ${r} nm h ${h} nm in ${dot.matrix}; dashed = the optical model's uniform ${h} nm layer` };
    },
  });
  const gView = hud.group('view');
  hud.button(gView, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the x500 dot magnifier (or click the dot plane)' });
  if (dotMesh) dotMesh.userData.onClick = () => mag.set(true);

  const bookmarks = [
    { name: 'Facet', pos: [1.6, 1.9, 7.2], target: [0, 1.05, 0], cut: X, mag: false },
    { name: 'Side cut-away', pos: [9.5, 2.6, -1.2], target: [0, 1.0, -5.5], cut: 0, mag: false },
    { name: 'Top / aperture', pos: [0.8, 7.6, 0.2], target: [0, 1.15, -2.4], cut: X, mag: false },
    { name: 'Dot magnifier', pos: [1.1, 1.5, 3.6], target: [0, 1.12, 0], cut: X, mag: true },
  ];
  return {
    root,
    defaultView: { pos: [7.6, 5.6, 9.8], target: [-0.6, 0.2, -5.6] },
    bookmarks,
    cut: { min: -X, max: X, axis: 'x' },
    magnifier: mag,
    unitToUm: 1, worldPerUm: 1,
  };
}
