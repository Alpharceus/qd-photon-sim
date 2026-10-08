// Planar / micropillar nitride cavity (device / planar_cavity) and the
// thermal-only mesa (device / thermal_mesa).  True scale, uniform display
// factor S world units per um.  The card has no pillar geometry: only the
// thermal stack, the DBR diagnostic and a mode-volume PARAMETER exist, and
// each is labelled as such.
import * as THREE from 'three';
import { layerMaterial, glowMaterial, outline, fatLine, breakMark, makeMagnifier, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';

export function buildPlanarCavity(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const g = spec.geometry;
  const mesaD = g.mesa_diameter_um || 1;
  const S = 10 / Math.max(mesaD, 2);         // world units per um
  const R = (mesaD / 2) * S;
  const clip = stage.clipPlane;
  const root = new THREE.Group();
  const thermalG = new THREE.Group(), opticalG = new THREE.Group();
  root.add(thermalG, opticalG);
  const isCavity = spec.variant === 'planar_cavity';
  const lam = spec.overlays.emission?.lambda_nm;
  const { rgb } = lambdaToRGB(lam);

  const disc = (y0, t, role, opts = {}, rad = R) => {
    const m = new THREE.Mesh(new THREE.CylinderGeometry(rad, rad, t, 96, 1, false), layerMaterial(role, { clip, ...opts }));
    m.position.y = y0 + t / 2;
    outline(m, { clip, opacity: 0.3 });
    return m;
  };
  // Thermal stack, drawn top -> down in the listed order (top at y = 0).
  const th = spec.stacks.thermal.layers;
  let y = 0;
  const yOf = {};
  th.forEach((Ly, i) => {
    const t = (Ly.t_nm / 1000) * S;
    y -= t;
    const role = i === 0 ? 'core' : 'thermal';
    const m = disc(y, t, role, { opacity: isCavity ? (i === 0 ? 0.55 : 0.6) : (i === 0 ? 0.82 : 0.7) });
    if (i === 0 && isCavity) m.material.color.setHex(0x6f8fae);
    stage.addPickable(m, { title: Ly.name, tag: Ly.tag, rows: [['thickness', `${nice(Ly.t_nm / 1000)} µm`], ['k300', `${Ly.k300} W/m/K`], ['stack', 'thermal (thermal.py)']], note: Ly.source || 'thermal.layers entry' });
    (i === 0 ? root : thermalG).add(m);
    yOf[i] = [y, y + t];
  });
  const gaNTop = 0;
  const gaNBot = yOf[0] ? yOf[0][0] : 0;
  const stackBot = y;
  // Ground plate: the substrate wafer extends beyond the mesa (thickness
  // and width not modelled, break marked); a heat-sink plate sits under it.
  const subT = 0.9 * S;
  const subR = R * (isCavity ? 1.35 : 2.6);
  const sub = disc(stackBot - subT, subT, 'substrate', { opacity: 1 }, subR);
  sub.material.color.setHex(0x2a2f37);
  const hsT = 0.35 * S;
  const hs = new THREE.Mesh(new THREE.BoxGeometry(subR * 2.6, hsT, subR * 2.6),
    new THREE.MeshPhysicalMaterial({ color: 0x3a3e45, metalness: 0.85, roughness: 0.38 }));
  hs.position.y = stackBot - subT - hsT / 2 - 0.02;
  root.add(hs);
  labels.add(root, [subR * 1.3, stackBot - subT - hsT / 2, subR * 1.3], `heat sink \u00b7 T_hs ${spec.scalars?.T_hs?.value ?? ''} K`, { tag: spec.scalars?.T_hs?.tag || 'A', priority: 25, sub: 'thermal boundary of thermal.py' });
  stage.addPickable(sub, { title: `${g.substrate} substrate`, tag: 'A', rows: [['thickness', 'not modelled']] });
  root.add(sub);
  root.add(breakMark([-subR, stackBot - subT * 0.5, subR + 0.01], [subR, stackBot - subT * 0.5, subR + 0.01], 0.12, [0, 1, 0], 12));
  if (!isCavity) {
    // Dot plane: the card gives no depth for it; drawn mid-film and labelled.
    const yd = (yOf[0][0] + yOf[0][1]) / 2;
    const dp = new THREE.Mesh(new THREE.CylinderGeometry(R * 1.002, R * 1.002, Math.max(0.004 * S, 0.012), 96),
      new THREE.MeshStandardMaterial({ color: 0xd9a845, emissive: 0xffb440, emissiveIntensity: 0.6, roughness: 0.4, transparent: true, opacity: 0.4, depthWrite: false }));
    dp.material.clippingPlanes = [clip];
    dp.position.y = yd;
    stage.markBloom(dp, 2.0);
    stage.addPickable(dp, { title: 'dot plane', tag: 'A', rows: [['depth', 'not a model input']], note: 'drawn mid-film; this card has no optical geometry' });
    root.add(dp);
    labels.add(root, [-R, yd, 0], 'dot plane (depth not modelled)', { tag: 'A', anchor: 'right', priority: 80, sub: 'drawn mid-film' });
  }

  // Optical toggle: the DBR diagnostic replaces the effective DBR+contact layer.
  if (isCavity && spec.stacks.optical) {
    const L = spec.stacks.optical.layers;
    const geo = new THREE.CylinderGeometry(R, R, 1, 96);
    const hiM = layerMaterial('mirror', { clip, opacity: 0.9 }); hiM.color.setHex(0xc9d2de);
    const loM = layerMaterial('mirror', { clip, opacity: 0.55 }); loM.color.setHex(0x5b6573);
    const nH = L.filter((x) => x.name.startsWith('H')).length;
    const hi = new THREE.InstancedMesh(geo, hiM, nH), lo = new THREE.InstancedMesh(geo, loM, L.length - nH);
    let yy = gaNBot, ih = 0, il = 0;
    const M = new THREE.Matrix4();
    for (const Ly of L) {
      const t = (Ly.t_nm / 1000) * S;
      yy -= t;
      M.compose(new THREE.Vector3(0, yy + t / 2, 0), new THREE.Quaternion(), new THREE.Vector3(1, t, 1));
      if (Ly.name.startsWith('H')) hi.setMatrixAt(ih++, M); else lo.setMatrixAt(il++, M);
    }
    const dbr = spec.overlays.dbr;
    const info = { title: 'DBR diagnostic', tag: 'E', rows: [['pairs', dbr.pairs], ['d_H / d_L', `${nice(dbr.d_high_nm)} / ${nice(dbr.d_low_nm)} nm`], ['R at λ', nice(dbr.reflectance, 5)], ['stop band', `${dbr.stopband_nm.map((v) => nice(v, 4)).join(' - ')} nm`]], note: 'nitride_cavity.dbr_diagnostic: ideal lossless quarter-wave stack, free-space exit [E/A]' };
    stage.addPickable(hi, info); stage.addPickable(lo, info);
    opticalG.add(hi, lo);
    const bottom = yy;
    opticalG.add(fatLine([[R + 0.05, gaNBot, 0], [R + 0.05, bottom, 0]], { color: 0xe9e8e3, width: 1.4 }));
    labels.add(opticalG, [R + 0.1, (gaNBot + bottom) / 2, 0], `DBR diagnostic · ${dbr.pairs} pairs · ${nice(spec.stacks.optical.total_nm / 1000)} µm`, { tag: 'E', sub: `R ${nice(dbr.reflectance, 4)} at ${nice(lam, 4)} nm; replaces the 1 µm effective layer of the thermal model` });
    opticalG.visible = false;
    const gS = hud.group('stack');
    hud.toggle(gS, [['thermal', 'thermal (thermal.py)'], ['optical', 'DBR (nitride_cavity.py)']], 'thermal', (v) => {
      thermalG.visible = v === 'thermal'; opticalG.visible = v === 'optical';
      // the substrate moves below whichever stack is shown
      const bot = v === 'optical' ? bottom : stackBot;
      sub.position.y = bot - subT / 2;
      stage.requestRender();
    });
  }

  // Dot, mode-volume sphere and standing-wave glyph (cavity only).
  let dotY = (gaNTop + gaNBot) / 2;
  if (isCavity) {
    const mv = spec.overlays.mode_volume;
    const rS = mv.radius_um * S;
    const sphere = new THREE.Mesh(new THREE.SphereGeometry(rS, 48, 24), glowMaterial(rgb, { opacity: 0.35 }));
    sphere.position.y = dotY;
    stage.markBloom(sphere);
    stage.addPickable(sphere, { title: 'mode volume (parameter)', tag: 'A', rows: [['V', `${nice(mv.V_norm)} (λ/n)³ = ${nice(mv.V_um3)} µm³`], ['n (GaN)', nice(mv.n, 4)], ['equal-volume sphere r', `${nice(mv.radius_um * 1000)} nm`]], note: mv.label });
    root.add(sphere);
    const sw = spec.overlays.standing_wave;
    const per = (sw.period_nm / 1000) * S;
    const n = Math.min(40, Math.floor((gaNTop - gaNBot) / per));
    const geo = new THREE.CylinderGeometry(rS * 2.2, rS * 2.2, per * 0.18, 48);
    const fr = new THREE.InstancedMesh(geo, glowMaterial(rgb, { opacity: 0.22 }), n);
    const M = new THREE.Matrix4();
    for (let i = 0; i < n; i++) {
      const k = i - Math.floor(n / 2);
      const yy = dotY + k * per;
      const fade = Math.exp(-Math.pow(k / (n / 4), 2));
      M.compose(new THREE.Vector3(0, yy, 0), new THREE.Quaternion(), new THREE.Vector3(fade, 1, fade));
      fr.setMatrixAt(i, M);
    }
    stage.markBloom(fr);
    root.add(fr);
    labels.add(root, [rS * 2.4, dotY + per * 3, 0], `mode volume ${nice(mv.V_norm)} (λ/n)³`, { tag: 'A', sub: '[A] mode volume, no 3D geometry in model' });
    labels.add(root, [-rS * 2.4, dotY + per * 6, 0], `standing wave λ/2n = ${nice(sw.period_nm)} nm`, { tag: 'DR', anchor: 'right', sub: sw.label });
    labels.add(root, [R * 0.15, dotY - 0.9, R * 0.5], `InGaN dot · h ${nice(g.dot.height_nm)} nm · r ${nice(g.dot.radius_nm)} nm`, { tag: 'A', anchor: 'left', sub: g.dot_depth_note });
  }

  // Annotations: thermal layers.
  th.forEach((Ly, i) => {
    const [y0, y1] = yOf[i];
    labels.add(i === 0 ? root : thermalG, [R + 0.1, (y0 + y1) / 2 + (i === 0 ? 0 : 0), 0], `${Ly.name} · ${nice(Ly.t_nm / 1000)} µm`, { tag: Ly.tag, sub: 'thermal stack (thermal.py)' });
  });
  labels.add(root, [R, stackBot - subT * 0.5, 0.2], `${g.substrate} substrate`, { tag: 'A', sub: 'thickness not modelled' });
  labels.add(root, [0, gaNTop + 0.15, -R], `mesa Ø ${nice(mesaD)} µm`, { tag: 'A', anchor: 'center' });

  let mag = null;
  if (isCavity) {
    mag = makeMagnifier(stage, hud, {
      mag: 500, title: 'section through the dot (true scale)',
      buildInset(scene) {
        const h = g.dot.height_nm, r = g.dot.radius_nm;
        const q = (w, hh, x, yy, z, c) => { const o = new THREE.Mesh(new THREE.PlaneGeometry(w, hh), new THREE.MeshBasicMaterial({ color: c, toneMapped: false })); o.position.set(x, yy, z); scene.add(o); };
        q(4000, 4000, 0, 0, 0, 0x6f8fae);
        q(2 * r, h, 0, 0, 0.1, new THREE.Color(...rgb).getHex());
        return { center: [0, 0], insetUnitsPerWorld: 1000 / S, note: `In0.25GaN dot r ${r} nm h ${h} nm in GaN; depth not modelled` };
      },
    });
    const gv = hud.group('view');
    hud.button(gv, 'magnifier', () => mag.set(!mag.visible));
  }

  const D = R * 2.6;
  // The cut-away range spans the whole wafer (the substrate disc is wider
  // than the mesa); at its maximum nothing is cut.
  const cutMax = subR * 1.3;
  return {
    root,
    defaultView: isCavity ? { pos: [D * 0.8, D * 0.36, D * 1.05], target: [2.2, -0.9, 0] }
      : { pos: [R * 9.5, R * 3.2, R * 12.5], target: [R * 0.9, stackBot * 0.6, 0] },
    bookmarks: [
      { name: 'Mesa', pos: [D * 0.95, D * 0.55, D * 1.15], target: [0, -1.0, 0], cut: cutMax, mag: false },
      { name: 'Side cut-away', pos: [D * 1.3, -0.4, 0.2], target: [0, -1.2, 0], cut: 0, mag: false },
      { name: 'Top', pos: [0.4, D * 1.4, 0.6], target: [0, 0, 0], cut: cutMax, mag: false },
      { name: 'Dot magnifier', pos: [R * 0.25, dotY + R * 0.2, R * 0.9], target: [0, dotY, 0], cut: 0, mag: !!mag },
    ],
    cut: { min: -R, max: cutMax },
    magnifier: mag,
    unitToUm: 1, worldPerUm: S,
    readoutKeys: isCavity ? ['g2_op', 'flux', 'Q', 'kappa_meV', 'Fp_add', 'lambda_nm', 'dbr_R', 'T_hs'] : undefined,
  };
}
