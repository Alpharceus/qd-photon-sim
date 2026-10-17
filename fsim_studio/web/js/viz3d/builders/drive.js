// Drive layer of a composed device (R2): what excites the dot, drawn from the descriptor's `drive` block
// (fsim_studio/scene_compose.py _drive_block).  The structure builders know nothing about it: index.js adds the rig
// after the builder has run, anchored on the dot that composeDot() recorded (stage.dotAnchor).
//
//   optical drives   a laser head and a pump beam in its wavelength colour (a labelled false colour outside the
//                    visible: dashed core), a pulse train for pulsed drives and a continuous beam for cw, a spectral
//                    ruler with the dot line and the true detuning of every component (ARP: the swept band)
//   two-colour       two beams (SUPER, dichromatic), blob size = intensity weight where the block gives one
//   TPE              each pulse is drawn as two photons (a double core)
//   above-band /     a translucent absorption disc (barrier) or a ring (p shell) round the dot with relaxation arrows
//   p-shell
//   electrical       p and n contacts, a current path through the dot with arrows (dashed + pips when pulsed)
//   stark-gated      a lateral gate electrode beside the dot plus the cw laser
//
// Pulse trains are schematic in time (the pulse spacing is not f_rep to scale); the labels say so.  Detunings on the
// ruler are to scale.  Nothing here computes physics: every number comes from the descriptor.
import * as THREE from 'three';
import { fatLine, glowMaterial, nice } from './common.js';
import { lambdaToRGB } from '../core/colormap.js';

const UP = new THREE.Vector3(0, 1, 0);
const Z = new THREE.Vector3(0, 0, 1);
const SRC_DIR = new THREE.Vector3(-0.84, 0.3, 0.46).normalize();     // from the dot towards the laser head
const hex = (rgb) => new THREE.Color(rgb[0], rgb[1], rgb[2]).getHex();
const fmtSigned = (v) => (v === 0 ? '0' : `${v > 0 ? '+' : '-'}${nice(Math.abs(v), 3)}`);

function orientY(obj, dir) { obj.quaternion.setFromUnitVectors(UP, dir.clone().normalize()); return obj; }

function beamRgb(b) {
  const { rgb, falseColour } = lambdaToRGB(b.lambda_nm);
  return { rgb, falseColour };
}

function tempoText(d) {
  const p = d.pulse || {};
  if (d.temporal === 'cw') return 'cw';
  if (d.temporal === 'dc') return 'dc';
  return `pulsed${p.f_rep_MHz ? `, ${nice(p.f_rep_MHz, 3)} MHz` : ''}`;
}

export function addDrive(spec, { stage, labels, devBox, anchor }) {
  const d = spec.drive;
  if (!d || !d.family || d.family === 'none') return null;
  const box = devBox && !devBox.isEmpty() ? devBox : new THREE.Box3(new THREE.Vector3(-2, -2, -2), new THREE.Vector3(2, 2, 2));
  const size = box.getSize(new THREE.Vector3());
  // scale of the rig: the device's height (or 0.6 of its widest side for flat devices), so a wide substrate plate does not push the head out of the default view
  const U = Math.max(size.y, 0.6 * Math.max(size.x, size.z), 1);
  const dot = anchor?.pos ? anchor.pos.clone() : box.getCenter(new THREE.Vector3());
  const dotR = anchor?.r ?? 0.05;
  const rig = new THREE.Group();
  rig.name = 'drive-rig';
  stage.scene.add(rig);
  const probes = [];
  const info = { kind: d.kind, family: d.family, temporal: d.temporal, beams: 0, pulses: 0, hasGate: false, hasContacts: false };

  const rb = Math.min(Math.max(U * 0.012, 0.035), 0.14);               // beam radius (display)
  const L = Math.min(Math.max(U * 0.55, 1.6), 8);                       // beam length (display)
  const clear = Math.max(dotR * 3.6, 0.14);                             // the beam stops short of the dot halo
  const hw = Math.min(Math.max(U * 0.025, 0.1), 0.34);                  // laser head half-width
  const beams = d.beams || [];

  // ------------------------------------------------------------ optical beams
  const heads = [];
  beams.forEach((b, i) => {
    const n = beams.length;
    const s = SRC_DIR.clone().applyAxisAngle(UP, (i - (n - 1) / 2) * 0.18).normalize();   // dot -> head
    const toDot = s.clone().negate();
    const src = dot.clone().addScaledVector(s, L);
    const end = dot.clone().addScaledVector(s, clear);
    const len = src.distanceTo(end);
    const { rgb, falseColour } = beamRgb(b);
    const w = Math.max(0.55, Math.sqrt(b.weight ?? 1));
    // pulse size ~ sqrt(area) (schematic: the area is the pulse's integrated Rabi angle in units of pi)
    const areaK = d.pulse?.area_pi != null ? Math.min(Math.max(0.7 + 0.25 * Math.sqrt(d.pulse.area_pi), 0.7), 1.7) : 1;
    const cw = d.temporal === 'cw';
    const col = hex(rgb);
    // core line
    const core = fatLine([src.toArray(), end.toArray()], { color: col, width: cw ? 4 : 2, opacity: cw ? 1 : 0.6, dashed: falseColour, dashSize: 0.18, gapSize: 0.1 });
    core.name = 'drive-beam';
    rig.add(core);
    if (d.kind === 'tpe') {
      // two photons per pulse: a second, parallel core
      const off = new THREE.Vector3().crossVectors(s, UP).normalize().multiplyScalar(rb * 2.2);
      const c2 = fatLine([src.clone().add(off).toArray(), end.clone().add(off).toArray()], { color: col, width: 2, opacity: 0.6, dashed: falseColour, dashSize: 0.18, gapSize: 0.1 });
      c2.name = 'drive-beam-2';
      rig.add(c2);
    }
    // glow
    if (cw) {
      const cyl = new THREE.Mesh(new THREE.CylinderGeometry(rb * 1.5 * w, rb * 1.5 * w, len, 20, 1, true), glowMaterial(rgb, { opacity: falseColour ? 0.22 : 0.3 }));
      cyl.position.copy(src).add(end).multiplyScalar(0.5);
      orientY(cyl, toDot);
      cyl.name = 'drive-beam-glow';
      stage.markBloom(cyl, 2.0);
      rig.add(cyl);
      probes.push(src.clone().lerp(end, 0.75));
    } else {
      const N = 6;
      const chirp = d.kind === 'arp' && d.sweep;
      for (let k = 0; k < N; k++) {
        const t = (k + 0.5) / N;
        const c = src.clone().lerp(end, t * 0.94);
        const plen = rb * (chirp ? 14 : 5);
        if (chirp) {
          // a stretched pulse: the colour ramps along the pulse in the direction of the chirp
          const a = beamRgb({ lambda_nm: d.sweep.from_nm }), z = beamRgb({ lambda_nm: d.sweep.to_nm });
          const dark = (c0) => c0.falseColour ? c0.rgb.map((v) => v * 0.35) : c0.rgb;
          const c0 = dark(a), c1 = z.falseColour ? z.rgb : z.rgb;
          const p0 = c.clone().addScaledVector(toDot, -plen / 2), p1 = c.clone().addScaledVector(toDot, plen / 2);
          const g = fatLine([p0.toArray(), p1.toArray()], { color: 0xffffff, width: 7, opacity: 0.95 });
          g.geometry.setColors([...c0, ...c1]);
          g.material.vertexColors = true;
          g.material.needsUpdate = true;
          g.name = 'drive-pulse';
          rig.add(g);
          if (k === N - 1) probes.push(c.clone());
        } else {
          const blob = new THREE.Mesh(new THREE.SphereGeometry(1, 16, 10), glowMaterial(rgb, { opacity: falseColour ? 0.7 : 0.9 }));
          blob.scale.set(rb * 2.2 * w * areaK, plen / 2 * w * areaK, rb * 2.2 * w * areaK);
          blob.position.copy(c);
          orientY(blob, toDot);
          blob.name = 'drive-pulse';
          stage.markBloom(blob, 2.4);
          rig.add(blob);
          if (k === N - 1) probes.push(c.clone());
        }
        info.pulses += 1;
      }
    }
    // laser head: a dark body with an aperture
    const body = new THREE.Mesh(new THREE.BoxGeometry(hw * 2, hw * 2, hw * 3.2),
      new THREE.MeshStandardMaterial({ color: 0x2b3038, metalness: 0.55, roughness: 0.45 }));
    body.position.copy(src).addScaledVector(s, hw * 1.6);
    body.quaternion.setFromUnitVectors(Z, s);
    body.name = 'drive-head';
    const ap = new THREE.Mesh(new THREE.SphereGeometry(hw * 0.55, 16, 10), glowMaterial(rgb, { opacity: 0.95 }));
    ap.position.copy(src);
    stage.markBloom(ap, 2.2);
    rig.add(body, ap);
    stage.addPickable(body, {
      title: b.label || b.role, tag: b.tag || 'A',
      rows: [['wavelength', `${nice(b.lambda_nm, 5)} nm`], ['colour', b.colour_label || (falseColour ? 'false colour' : 'true colour')],
        ['detuning from the line', `${fmtSigned(b.detuning_meV)} meV`], ['timing', tempoText(d)]],
      note: (falseColour ? 'outside the visible: drawn in a stand-in colour, dashed' : 'drawn in the wavelength colour') + '; pulse size ~ sqrt(area), spacing not to time scale',
    });
    heads.push({ src, s, b, rgb, falseColour });
    info.beams += 1;
  });

  // ------------------------------------------------------------ labels for the beams (one per beam)
  heads.forEach((h, i) => {
    const pos = h.src.clone().addScaledVector(h.s, hw * 3.4).add(new THREE.Vector3(0, hw * 1.6 + i * hw * 2.6, 0));
    const nmTxt = `${nice(h.b.lambda_nm, 5)} nm`;
    const note = [h.b.colour_label, tempoText(d), h.b.display_default ? 'wavelength is a display default [A]' : null].filter(Boolean).join(' · ');
    labels.add(rig, pos.toArray(), `${h.b.label || d.label} · ${nmTxt}`,
      { tag: h.b.tag || 'A', anchor: 'right', priority: 80 - i * 2, sub: note });
  });

  // ------------------------------------------------------------ spectral ruler (to scale)
  const showRuler = heads.length && d.kind !== 'above_band' && d.kind !== 'p_shell';
  if (showRuler) {
    const deltas = heads.map((h) => h.b.delta_meV ?? 0);
    const sweep = d.sweep?.span_meV || 0;
    const hr = Math.max(1.25 * Math.max(...deltas.map(Math.abs), sweep / 2), 2);    // meV at half the ruler
    const W = Math.min(Math.max(U * 0.22, 1.2), 2.4);
    const mid = heads.reduce((a, h) => a.add(h.src), new THREE.Vector3()).multiplyScalar(1 / heads.length);
    const rc = mid.clone().add(new THREE.Vector3(0, hw * 5.2, 0));
    const seg = (x0, y0, x1, y1, o) => {
      const l = fatLine([[rc.x + x0, rc.y + y0, rc.z], [rc.x + x1, rc.y + y1, rc.z]], o);
      l.name = 'drive-ruler';
      rig.add(l);
      return l;
    };
    seg(-W, 0, W, 0, { color: 0x8a9099, width: 1.4, opacity: 0.9 });
    seg(0, -0.1, 0, 0.1, { color: 0xe9e8e3, width: 2.4 });                                   // the dot line
    if (sweep) seg(-sweep / 2 / hr * W + deltas[0] / hr * W, 0.02, sweep / 2 / hr * W + deltas[0] / hr * W, 0.02, { color: hex(heads[0].rgb), width: 6, opacity: 0.8 });
    heads.forEach((h, i) => {
      const x = (deltas[i] / hr) * W;
      if (Math.abs(x) > 1e-6 && !sweep) seg(0, 0.02, x, 0.02, { color: hex(h.rgb), width: 9, opacity: 0.55 });      // the gap the laser sits off the line
      seg(x, 0, x, 0.2 + 0.04 * i, { color: hex(h.rgb), width: 5 });
    });
    const txt = d.kind === 'super' ? `Δ1 = ${fmtSigned(deltas[0])} meV, Δ2 = ${fmtSigned(deltas[1])} meV`
      : d.kind === 'dichromatic' ? `Δ = ±${nice(Math.abs(deltas[0]), 3)} meV`
        : `Δ = ${fmtSigned(deltas[0])} meV`;
    const subTxt = d.kind === 'tpe' ? 'from the two-photon resonance; 2 photons per pulse'
      : d.kind === 'arp' ? `sweep ${nice(d.sweep.span_meV, 3)} meV, chirp ${d.sweep.direction}`
        : 'positive = laser blue of the line; ruler to scale';
    labels.add(rig, [rc.x + W, rc.y + 0.16, rc.z], txt, { tag: 'A', anchor: 'left', priority: 60, sub: subTxt });
  }

  // ------------------------------------------------------------ above-band / p-shell: absorption region
  if (d.kind === 'above_band' || d.kind === 'p_shell') {
    const rgb = heads[0]?.rgb || [1, 1, 1];
    if (d.kind === 'above_band') {
      const R = Math.max(U * 0.1, 0.5);
      const disc = new THREE.Mesh(new THREE.CylinderGeometry(R, R, Math.max(U * 0.004, 0.015), 48), glowMaterial(rgb, { opacity: 0.18 }));
      disc.position.copy(dot);
      disc.name = 'drive-absorber';
      rig.add(disc);
      labels.add(rig, [dot.x + R, dot.y - 0.1, dot.z + R * 0.2], 'barrier absorption', { tag: 'A', priority: 55, sub: 'carriers captured by the dot, then relax' });
    } else {
      const R = Math.max(dotR * 6, 0.24);
      const ring = new THREE.Mesh(new THREE.TorusGeometry(R, Math.max(R * 0.04, 0.012), 10, 48), glowMaterial(rgb, { opacity: 0.6 }));
      ring.position.copy(dot);
      ring.rotation.x = Math.PI / 2;
      ring.name = 'drive-absorber';
      stage.markBloom(ring, 1.8);
      rig.add(ring);
      labels.add(rig, [dot.x + R, dot.y - 0.05, dot.z + R * 0.2], 'p-shell', { tag: 'A', priority: 55, sub: 'relaxes to the s shell' });
    }
    for (let k = 0; k < 3; k++) {
      const a = (k / 3) * Math.PI * 2;
      const r = Math.max(dotR * 5, 0.22) * (d.kind === 'above_band' ? 3 : 1);
      const p = dot.clone().add(new THREE.Vector3(Math.cos(a) * r, 0, Math.sin(a) * r));
      const dir = dot.clone().sub(p).normalize();
      const cone = new THREE.Mesh(new THREE.ConeGeometry(Math.max(dotR * 0.8, 0.03), Math.max(dotR * 2.2, 0.09), 12), glowMaterial(rgb, { opacity: 0.85 }));
      cone.position.copy(p);
      orientY(cone, dir);
      cone.name = 'drive-relax';
      rig.add(cone);
    }
  }

  // ------------------------------------------------------------ electrical: contacts and the current path
  if (d.family === 'electrical' && d.electrical) {
    const e = d.electrical;
    const px = Math.max(Math.min(size.x / 2, U * 0.6) + U * 0.1, 0.9), pw = Math.max(U * 0.05, 0.16), ph = Math.max(U * 0.16, 0.7);
    const pm = new THREE.MeshStandardMaterial({ color: 0xb08a3c, metalness: 0.75, roughness: 0.4 });
    for (const [x, name] of [[-px, 'p contact'], [px, 'n contact']]) {
      const m = new THREE.Mesh(new THREE.BoxGeometry(pw, ph, ph), pm);
      m.position.set(dot.x + x, dot.y, dot.z);
      m.name = `drive-${name.replace(' ', '-')}`;
      rig.add(m);
      stage.addPickable(m, { title: name, tag: 'A', rows: [['role', name === 'p contact' ? 'hole injection' : 'electron injection']], note: 'display geometry: the contact layout is not a model input' });
    }
    const a = new THREE.Vector3(dot.x - px + pw / 2, dot.y, dot.z), z = new THREE.Vector3(dot.x + px - pw / 2, dot.y, dot.z);
    const path = fatLine([a.toArray(), z.toArray()], { color: 0xf2c14e, width: 3, opacity: 0.95, dashed: d.temporal === 'pulsed', dashSize: 0.2, gapSize: 0.14 });
    path.material.depthTest = false;
    path.renderOrder = 5;
    path.name = 'drive-current';
    rig.add(path);
    for (const t of [0.22, 0.5, 0.78]) {
      const c = new THREE.Mesh(new THREE.ConeGeometry(Math.max(rb * 2.2, 0.06), Math.max(rb * 6, 0.2), 14), new THREE.MeshBasicMaterial({ color: 0xf2c14e, depthTest: false, toneMapped: false }));
      c.position.copy(a).lerp(z, t);
      orientY(c, new THREE.Vector3(1, 0, 0));
      c.renderOrder = 6;
      c.name = 'drive-current-arrow';
      rig.add(c);
    }
    if (d.temporal === 'pulsed') {
      for (const t of [0.12, 0.36, 0.64, 0.88]) {
        const pip = new THREE.Mesh(new THREE.SphereGeometry(Math.max(rb * 1.4, 0.05), 14, 10), glowMaterial([1.0, 0.82, 0.3], { opacity: 0.9 }));
        pip.position.copy(a).lerp(z, t);
        pip.name = 'drive-current-pip';
        stage.markBloom(pip, 1.2);
        rig.add(pip);
      }
    }
    probes.push(a.clone().lerp(z, 0.36));
    info.hasContacts = true;
    const iTxt = e.I_uA != null ? `I = ${nice(e.I_uA, 3)} µA` : e.mu != null ? `μ = ${nice(e.mu, 3)} carriers per period` : String(e.mechanism);
    labels.add(rig, [dot.x - px, dot.y + ph / 2 + 0.05, dot.z], 'p contact', { tag: 'A', anchor: 'right', priority: 70 });
    labels.add(rig, [dot.x + px, dot.y + ph / 2 + 0.05, dot.z], 'n contact', { tag: 'A', anchor: 'left', priority: 70 });
    labels.add(rig, [dot.x, dot.y - ph / 2 - 0.1, dot.z], `${iTxt} · ${tempoText(d)}`,
      { tag: e.tag || 'A', anchor: 'center', priority: 75, sub: `${e.mechanism}; conventional current p to n` });
  }

  // ------------------------------------------------------------ Stark gate electrode
  if (d.family === 'stark_gated' && d.gate) {
    const g = d.gate;
    // a lateral gate electrode beside the device at the dot's height (display geometry): the Stark shift moves the line, the laser stays fixed
    const gh = Math.max(U * 0.16, 0.6), gt = Math.max(U * 0.02, 0.05);
    const gx = Math.max(Math.min(size.x / 2, U * 0.6) + U * 0.1, 0.9);
    const plate = new THREE.Mesh(new THREE.BoxGeometry(gt, gh, gh),
      new THREE.MeshStandardMaterial({ color: 0xe0b84f, emissive: 0x6a5214, metalness: 0.6, roughness: 0.4, transparent: true, opacity: 0.92 }));
    plate.position.set(dot.x + gx, dot.y, dot.z);
    plate.name = 'drive-gate';
    rig.add(plate);
    const lead = fatLine([[dot.x + gx, dot.y + gh / 2, dot.z], [dot.x + gx, dot.y + gh / 2 + Math.max(U * 0.12, 0.5), dot.z]], { color: 0xc9a24a, width: 2.4 });
    lead.name = 'drive-gate-lead';
    rig.add(lead);
    stage.addPickable(plate, {
      title: 'gate electrode', tag: 'A',
      rows: [['detuning off the gate', `${fmtSigned(g.det_off_meV)} meV`], ['detuning on the gate', `${fmtSigned(g.det_on_meV)} meV`], ['gate width', `${nice(g.t_gate_ps, 4)} ps`]],
      note: 'display geometry; the Stark shift moves the line, the laser stays fixed',
    });
    info.hasGate = true;
    probes.push(new THREE.Vector3(dot.x + gx, dot.y, dot.z));
    labels.add(rig, [dot.x + gx + gt, dot.y + gh / 2, dot.z], 'gate electrode', {
      tag: g.tag || 'A', anchor: 'left', priority: 78,
      sub: `Stark gate: Δ ${fmtSigned(g.det_off_meV)} → ${fmtSigned(g.det_on_meV)} meV for ${nice(g.t_gate_ps, 3)} ps`,
    });
  }

  rig.traverse((o) => { if (o.material?.isLineMaterial) o.material.needsUpdate = true; });
  return { rig, probes, info, dot };
}
