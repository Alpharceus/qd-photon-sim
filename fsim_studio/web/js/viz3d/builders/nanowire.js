// Disc-in-nanowire (device / nanowire_vertical | nanowire_horizontal).
// True scale; wire and doped-segment lengths are not modelled, so they end
// in break marks.  Numbers come from the card's committed sweep row.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, coneMaterial, outline, fatLine, breakMark, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';

export function buildNanowire(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const horizontal = spec.variant === 'nanowire_horizontal';
  const S = horizontal ? 100 : 40;          // world units per um (uniform)
  const nm = (v) => (v / 1000) * S;
  const clip = stage.clipPlane;
  const root = new THREE.Group();
  const wire = new THREE.Group();
  root.add(wire);
  // Wire axis: vertical -> +y; horizontal -> +x (lying on the oxide).
  if (horizontal) wire.rotation.z = -Math.PI / 2;

  const rc = nm(g.core_radius_nm);
  const segDisp = horizontal ? 140 : 380;   // drawn length of the unmodelled doped segments (nm)
  const lam = spec.overlays.emission.lambda_nm;
  const { rgb, falseColour } = lambdaToRGB(lam);
  const L = spec.layers;
  const tOf = (Ly) => (Ly.t_nm == null ? segDisp : Ly.t_nm);
  const total = L.reduce((a, Ly) => a + tOf(Ly), 0);
  let y = -total / 2;
  let discY = 0;
  for (const Ly of L) {
    const t = tOf(Ly);
    const y0 = y, y1 = y + t;
    y = y1;
    const isDisc = Ly.role === 'active';
    if (isDisc) {
      discY = nm((y0 + y1) / 2);
      // GaN around the disc, then the InGaN dot disc itself (radius_nm).
      const ring = new THREE.Mesh(new THREE.CylinderGeometry(rc, rc, nm(t), 64, 1, true), layerMaterial('intrinsic', { clip, opacity: 0.4 }));
      ring.position.y = discY;
      wire.add(ring);
      const rd = nm(g.dot.radius_nm);
      const disc = new THREE.Mesh(new THREE.CylinderGeometry(rd, rd, nm(t), 48),
        new THREE.MeshStandardMaterial({ color: new THREE.Color(...rgb), emissive: new THREE.Color(...rgb), emissiveIntensity: 2.2, roughness: 0.4 }));
      disc.position.y = discY;
      stage.markBloom(disc);
      stage.addPickable(disc, { title: 'InGaN dot disc', tag: Ly.tag, rows: [['height', `${nice(t)} nm`], ['radius', `${nice(g.dot.radius_nm)} nm`], ['emission', `${nice(lam, 4)} nm`]], note: Ly.source });
      wire.add(disc);
      continue;
    }
    const mat = layerMaterial(Ly.role, { clip, opacity: Ly.t_nm == null ? 0.5 : 0.42 });
    const m = new THREE.Mesh(new THREE.CylinderGeometry(rc, rc, nm(t), 64, 1, Ly.t_nm != null), mat);
    m.position.y = nm((y0 + y1) / 2);
    outline(m, { clip, opacity: 0.25 });
    stage.addPickable(m, { title: Ly.name, tag: Ly.tag, rows: [['thickness', Ly.t_nm == null ? 'not modelled' : `${nice(Ly.t_nm)} nm`], ['material', Ly.material], ['n', Ly.n ? Ly.n.toFixed(3) : 'n/a']], note: Ly.source });
    wire.add(m);
  }
  const yTop = nm(total / 2), yBot = -yTop;
  for (const yy of [yTop, yBot]) wire.add(breakMark([-rc * 1.25, yy, 0], [rc * 1.25, yy, 0], rc * 0.12, [0, 1, 0], 6));

  // HE11 guided-mode glow: opacity follows beta_HE11 (no field profile drawn).
  const he = spec.overlays.he11;
  const beta = Number.isFinite(he.beta_HE11) ? he.beta_HE11 : 0;
  if (beta > 0) {
    for (const [rf, op] of [[0.92, 0.012], [0.5, 0.02]]) {
      const c = new THREE.Mesh(new THREE.CylinderGeometry(rc * rf, rc * rf, yTop - discY, 48, 1, true), glowMaterial(rgb, { opacity: op * (0.4 + beta) }));
      c.position.y = (yTop + discY) / 2;
      stage.markBloom(c);
      wire.add(c);
    }
  }

  // Collection cone above the top (vertical) / above the wire (horizontal).
  const na = spec.overlays.na_cone;
  const th = (na.half_angle_deg * Math.PI) / 180;
  const coneLen = horizontal ? 12 : 10;
  const cone = new THREE.Mesh(new THREE.ConeGeometry(coneLen * Math.tan(th), coneLen, 64, 1, true), coneMaterial(rgb, 0.05));
  cone.rotation.x = Math.PI;
  const apexY = horizontal ? rc : yTop + 0.4;
  cone.position.set(0, apexY + coneLen / 2, 0);
  root.add(cone);
  const rim = [];
  for (let i = 0; i <= 96; i++) { const a = (i / 96) * Math.PI * 2; rim.push([Math.cos(a) * coneLen * Math.tan(th), apexY + coneLen, Math.sin(a) * coneLen * Math.tan(th)]); }
  root.add(fatLine(rim, { color: new THREE.Color(...rgb).getHex(), width: 1, opacity: 0.35 }));

  // Horizontal family: oxide + substrate under the wire, dipole axes.
  if (horizontal) {
    const ox = nm(g.oxide_thickness_nm || 100);
    const w = 2 * nm(segDisp + 40), d = rc * 10;
    const oxide = new THREE.Mesh(new THREE.BoxGeometry(w, ox, d), layerMaterial('oxide', { clip }));
    oxide.position.y = -rc - ox / 2;
    stage.addPickable(oxide, { title: 'oxide', tag: 'A', rows: [['thickness', `${g.oxide_thickness_nm} nm`]], note: 'nitride.photonics.oxide_thickness_nm' });
    const sub = new THREE.Mesh(new THREE.BoxGeometry(w, ox * 0.8, d), layerMaterial('substrate', { clip }));
    sub.position.y = -rc - ox - ox * 0.4;
    root.add(oxide, sub);
    root.add(breakMark([-w / 2, -rc - ox * 1.8, d / 2], [w / 2, -rc - ox * 1.8, d / 2], 0.25, [0, 1, 0], 10));
    const wts = g.dipole_weights || [];
    const axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    wts.forEach((wt, i) => {
      const len = rc * 5 * wt;
      const a = axes[i];
      root.add(fatLine([[-a[0] * len, -a[1] * len, -a[2] * len], [a[0] * len, a[1] * len, a[2] * len]], { color: 0xe9e8e3, width: 2, opacity: 0.9 }));
    });
    labels.add(root, [rc * 2.2, rc * 1.4, 0], `dipole weights ${wts.map((v) => nice(v, 2)).join(' / ')} (x, y, z)`, { tag: 'A', sub: 'arrow length ∝ weight' });
    labels.add(root, [w / 2 - 2, -rc - ox / 2, d / 2], `oxide ${nice(g.oxide_thickness_nm)} nm`, { tag: 'A', anchor: 'left' });
  }

  // Annotations.
  const side = rc + (horizontal ? 1.2 : 0.6);
  const segLabelPos = (frac) => (horizontal ? [nm(frac), rc + 1.6, 0] : [side, nm(frac), 0]);
  let acc = -total / 2;
  for (const Ly of L) {
    const t = tOf(Ly);
    const mid = acc + t / 2;
    acc += t;
    if (Ly.role === 'active') {
      labels.add(root, horizontal ? [0, -rc - 0.6, rc] : [-side, discY, 0], `InGaN disc · h ${nice(Ly.t_nm)} nm · r ${nice(g.dot.radius_nm)} nm`, { tag: Ly.tag, anchor: horizontal ? 'left' : 'right', sub: `λ ${nice(lam, 4)} nm${falseColour ? ' (NIR, false colour)' : ''}` });
      continue;
    }
    if (Ly.t_nm == null) {
      labels.add(root, segLabelPos(mid), `${Ly.name}`, { tag: Ly.tag, sub: 'length not modelled (break mark)' });
    } else {
      labels.add(root, horizontal ? [nm(mid), -rc - 0.4, rc] : [side, nm(mid) + (Ly.name.includes('left') ? -0.25 : 0.25), 0], `${Ly.name} · ${nice(Ly.t_nm)} nm`, { tag: Ly.tag });
    }
  }
  labels.add(root, horizontal ? [-nm(segDisp), rc * 1.2, -rc] : [-rc * 0.9, (yTop + discY) / 2, rc], `HE11 β ${nice(he.beta_HE11)} · V ${nice(he.V_number)} · ${he.single_mode ? 'single mode' : 'multimode'}`,
    { tag: he.tag, anchor: 'right', sub: 'glow opacity follows β (not a field profile)' });
  labels.add(root, [coneLen * Math.tan(th) * 0.7, apexY + coneLen * 0.85, 0], `collection NA ${nice(na.NA)} · ${nice(na.half_angle_deg, 3)}°`, { tag: na.tag });
  labels.add(root, horizontal ? [0, -rc - 1, -rc * 3] : [side, yBot - 0.2, 0], `core Ø ${nice(2 * g.core_radius_nm)} nm · shell ${g.shell}`, { tag: 'A', sub: g.segments_note });

  const mag = makeMagnifier(stage, hud, {
    mag: 20, title: 'section through the disc (true scale)',
    buildInset(scene) {
      const r = g.core_radius_nm, h = g.dot.height_nm, rd = g.dot.radius_nm;
      const gan = new THREE.MeshBasicMaterial({ color: 0x5d8b84 });
      const bar = new THREE.MeshBasicMaterial({ color: 0x8c9199 });
      const q = (w, hh, x, yy, z, m) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(w, hh), m); o.position.set(x, yy, z); scene.add(o); };
      q(2 * r, g.barrier_left_nm, 0, -h / 2 - g.barrier_left_nm / 2, 0, bar);
      q(2 * r, g.barrier_right_nm, 0, h / 2 + g.barrier_right_nm / 2, 0, bar);
      q(2 * r, h, 0, 0, 0, bar);
      q(2 * r, 200, 0, -h / 2 - g.barrier_left_nm - 100, 0, gan);
      q(2 * r, 200, 0, h / 2 + g.barrier_right_nm + 100, 0, gan);
      q(2 * rd, h, 0, 0, 0.1, new THREE.MeshBasicMaterial({ color: new THREE.Color(...rgb), toneMapped: false }));
      return { center: [0, 0], insetUnitsPerWorld: 1000 / S, note: `disc r ${rd} nm h ${h} nm, barriers ${g.barrier_left_nm}/${g.barrier_right_nm} nm, core r ${r} nm` };
    },
  });
  const gView = hud.group('view');
  hud.button(gView, 'magnifier', () => mag.set(!mag.visible), { title: 'toggle the disc magnifier' });

  const R = rc * 7;
  const bookmarks = horizontal
    ? [
      { name: 'Wire', pos: [R * 0.9, R * 0.7, R * 1.6], target: [0, 0, 0], cut: rc * 6, mag: false },
      { name: 'Side cut-away', pos: [0.1, 0.3, R * 1.8], target: [0, 0, 0], cut: 0, mag: false },
      { name: 'Top / cone', pos: [R * 0.3, R * 2.2, R * 0.6], target: [0, rc, 0], cut: rc * 6, mag: false },
      { name: 'Disc magnifier', pos: [R * 0.4, R * 0.3, R * 0.9], target: [0, 0, 0], cut: rc * 6, mag: true },
    ]
    : [
      { name: 'Wire', pos: [R * 1.6, R * 0.6, R * 2.2], target: [0, discY + 2, 0], cut: rc * 1.3, mag: false },
      { name: 'Side cut-away', pos: [R * 2.4, discY + 1.5, 0.2], target: [0, discY, 0], cut: 0, mag: false },
      { name: 'Top / cone', pos: [R * 0.8, yTop + R * 2.2, R * 0.9], target: [0, yTop, 0], cut: rc * 1.3, mag: false },
      { name: 'Disc magnifier', pos: [R * 0.7, discY + 1.2, R * 1.1], target: [0, discY, 0], cut: rc * 1.3, mag: true },
    ];
  return {
    root,
    defaultView: horizontal ? { pos: [R * 1.3, R * 1.0, R * 2.0], target: [0, 0.6, 0] }
      : { pos: [R * 2.2, discY + R * 0.9, R * 3.2], target: [0, discY + 2.5, 0] },
    bookmarks,
    cut: horizontal ? { min: -rc * 6, max: rc * 6 } : { min: -rc * 1.3, max: rc * 1.3 },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: ['g2_op', 'flux_commanded', 'flux_delivered', 'beta_HE11', 'V_number', 'lambda_nm', 'eta_collection_X', 'T_hs'],
  };
}
