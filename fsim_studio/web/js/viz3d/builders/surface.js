// g2(T, param) surfaces (surface / *): two sheets (lo, hi), never one
// averaged surface; the g2 = 0.5 isoline (the T_c ridge) on each; T_c
// interval points from evaluate_envelope; a floor heat-map shadow of the lo
// sheet for projector legibility.  Contouring here is drawing, not physics.
import * as THREE from 'three';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { fatLine, nice } from './common.js';
import { viridis, viridisGradientCss } from '../core/colormap.js';

const WX = 12, WZ = 9, HY = 5;

function isoSegments(grid, xs, zs, level) {
  const segs = [];
  const ny = grid.length, nx = grid[0].length;
  const ok = (v) => v !== null && Number.isFinite(v);
  for (let j = 0; j < ny - 1; j++) {
    for (let i = 0; i < nx - 1; i++) {
      const c = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]];
      const v = c.map(([a, b]) => grid[b][a]);
      if (!v.every(ok)) continue;
      const pts = [];
      for (let k = 0; k < 4; k++) {
        const a = v[k], b = v[(k + 1) % 4];
        if ((a - level) * (b - level) < 0 || (a === level && b !== level)) {
          const f = (level - a) / (b - a);
          const [ia, ja] = c[k], [ib, jb] = c[(k + 1) % 4];
          pts.push([xs[ia] + (xs[ib] - xs[ia]) * f, zs[ja] + (zs[jb] - zs[ja]) * f]);
        }
      }
      if (pts.length >= 2) segs.push([pts[0], pts[1]]);
      if (pts.length === 4) segs.push([pts[2], pts[3]]);
    }
  }
  return segs;
}

function sheet(grid, xs, zs, opacity) {
  const ny = grid.length, nx = grid[0].length;
  const pos = new Float32Array(nx * ny * 3), col = new Float32Array(nx * ny * 3);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = grid[j][i];
    const k = (j * nx + i) * 3;
    const ok = v !== null && Number.isFinite(v);
    pos.set([xs[i], ok ? v * HY : 0, zs[j]], k);
    col.set(viridis(ok ? v : NaN), k);
  }
  const idx = [];
  const ok = (i, j) => grid[j][i] !== null && Number.isFinite(grid[j][i]);
  for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
    if (!(ok(i, j) && ok(i + 1, j) && ok(i, j + 1) && ok(i + 1, j + 1))) continue;
    const a = j * nx + i, b = a + 1, c = a + nx, d = c + 1;
    idx.push(a, c, b, b, c, d);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  const m = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.45, metalness: 0.05, side: THREE.DoubleSide, transparent: true, opacity, depthWrite: opacity > 0.8 });
  return new THREE.Mesh(g, m);
}

function segLines(segs, y, color, width, dashed = false) {
  if (!segs.length) return null;
  const g = new LineSegmentsGeometry();
  g.setPositions(segs.flatMap(([a, b]) => [a[0], y, a[1], b[0], y, b[1]]));
  const m = new LineMaterial({ color, linewidth: width, dashed, dashSize: 0.18, gapSize: 0.12 });
  const l = new LineSegments2(g, m);
  if (dashed) l.computeLineDistances();
  return l;
}

export function buildSurface(spec, ctx) {
  const { stage, hud, labels } = ctx;
  const o = spec.overlays;
  const ax = spec.axes;
  const root = new THREE.Group();
  const T = o.x, P = o.y;
  const tMin = Math.min(...T), tMax = Math.max(...T);
  const logY = !!ax.y.log;
  const yv = (v) => (logY ? Math.log10(v) : v);
  const pMin = yv(Math.min(...P)), pMax = yv(Math.max(...P));
  const xs = T.map((t) => ((t - tMin) / (tMax - tMin || 1)) * WX);
  const zs = P.map((p) => ((yv(p) - pMin) / (pMax - pMin || 1)) * WZ);
  const xOfT = (t) => ((t - tMin) / (tMax - tMin || 1)) * WX;
  const zOfP = (p) => ((yv(p) - pMin) / (pMax - pMin || 1)) * WZ;

  const lo = sheet(o.lo, xs, zs, 0.92), hi = sheet(o.hi, xs, zs, 0.42);
  stage.addPickable(lo, { title: 'lo sheet', tag: 'A', rows: [['meaning', 'pointwise minimum']], note: o.sheet_note });
  stage.addPickable(hi, { title: 'hi sheet', tag: 'A', rows: [['meaning', 'pointwise maximum']], note: o.sheet_note });
  root.add(lo, hi);
  const wire = new THREE.LineSegments(new THREE.WireframeGeometry(hi.geometry), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.12 }));
  root.add(wire);

  // g2 = 0.5 isolines on each sheet + projected on the floor.
  const isoLo = isoSegments(o.lo, xs, zs, o.iso), isoHi = isoSegments(o.hi, xs, zs, o.iso);
  const Li = segLines(isoLo, o.iso * HY + 0.01, 0xf4f1e6, 3.2);
  const Hi = segLines(isoHi, o.iso * HY + 0.01, 0xf4f1e6, 2.4, true);
  for (const l of [Li, Hi]) if (l) { stage.markBloom(l, 1.6); root.add(l); }
  const Fl = segLines(isoLo, 0.02, 0xf4f1e6, 1.4);
  if (Fl) root.add(Fl);

  // Floor heat-map shadow (lo sheet).
  const nx = T.length, ny = P.length;
  const cv = document.createElement('canvas');
  cv.width = nx; cv.height = ny;
  const c2 = cv.getContext('2d');
  const img = c2.createImageData(nx, ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = o.lo[j][i];
    const [r, g, b] = viridis(v === null ? NaN : v);
    const k = ((ny - 1 - j) * nx + i) * 4;
    img.data.set([r * 255, g * 255, b * 255, v === null ? 60 : 255], k);
  }
  c2.putImageData(img, 0, 0);
  const tex = new THREE.CanvasTexture(cv);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.magFilter = THREE.LinearFilter; tex.minFilter = THREE.LinearFilter;
  // Texel centres sit on the data points; the plane spans half a cell beyond.
  const dx = WX / Math.max(1, nx - 1), dz = WZ / Math.max(1, ny - 1);
  const floor = new THREE.Mesh(new THREE.PlaneGeometry(WX + dx, WZ + dz), new THREE.MeshBasicMaterial({ map: tex, transparent: true, opacity: 0.85 }));
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(WX / 2, -0.005, WZ / 2);
  root.add(floor);

  // T_c interval points (evaluate_envelope scalar band) on the iso height.
  const tcGeo = new THREE.SphereGeometry(0.055, 16, 12);
  const tcMat = new THREE.MeshBasicMaterial({ color: 0xffffff, toneMapped: false });
  let nTc = 0;
  for (const p of o.tc_points || []) {
    for (const t of [p.T_lo, p.T_hi]) {
      if (t === null || !Number.isFinite(t) || t < tMin || t > tMax) continue;
      const s = new THREE.Mesh(tcGeo, tcMat);
      s.position.set(xOfT(t), o.iso * HY, zOfP(p.y));
      stage.markBloom(s, 1.5);
      root.add(s);
      nTc += 1;
    }
  }

  // Axis frame and ticks.
  const frameCol = 0x4a4d55;
  root.add(fatLine([[0, 0, WZ], [WX, 0, WZ]], { color: frameCol, width: 1 }));
  root.add(fatLine([[0, 0, 0], [0, 0, WZ]], { color: frameCol, width: 1 }));
  root.add(fatLine([[0, 0, 0], [0, HY, 0]], { color: frameCol, width: 1 }));
  root.add(fatLine([[0, HY * o.iso, 0], [WX, HY * o.iso, 0]], { color: 0x7f7d77, width: 1, dashed: true }));
  const tStep = (tMax - tMin) > 150 ? 50 : (tMax - tMin) > 40 ? 20 : 10;
  for (let t = Math.ceil(tMin / tStep) * tStep; t <= tMax + 1e-9; t += tStep) {
    labels.add(root, [xOfT(t), 0, WZ + 0.45], `${t}`, { kind: 'tick', anchor: 'center' });
  }
  labels.add(root, [WX / 2, 0, WZ + 1.1], `${ax.x.label} (${ax.x.unit})`, { kind: 'axis', anchor: 'center' });
  const pTicks = P.length <= 6 ? P : P.filter((_, i) => i % Math.ceil(P.length / 6) === 0);
  for (const p of pTicks) labels.add(root, [WX + 0.35, 0, zOfP(p)], `${nice(p, 3)}`, { kind: 'tick', anchor: 'left' });
  labels.add(root, [WX + 0.5, 0, -0.7], `${ax.y.label}${ax.y.unit ? ` (${ax.y.unit})` : ''}${logY ? ', log' : ''}`, { kind: 'axis', anchor: 'left' });
  root.add(fatLine([[WX, 0, 0], [WX, 0, WZ]], { color: frameCol, width: 1 }));
  for (const g of [0, 0.5, 1]) labels.add(root, [-0.15, g * HY, 0], `g2 ${g}`, { kind: 'tick', anchor: 'right' });
  labels.add(root, [WX, o.iso * HY, 0], 'g2 = 0.5', { kind: 'chip', anchor: 'left', sub: isoLo.length || isoHi.length ? 'isoline = T_c ridge (solid lo, dashed hi)' : 'no crossing in this grid' });

  // Legend: colour bar + sheet key.
  const bar = document.createElement('div');
  bar.className = 'v3-colorbar';
  bar.innerHTML = `<div class="v3-colorbar-title">g2(0), colour (viridis)</div><div class="v3-colorbar-bar" style="background:${viridisGradientCss()}"><span class="v3-colorbar-mark" style="left:50%"></span></div><div class="v3-colorbar-ticks"><span>0</span><span>0.5</span><span>1</span></div>`;
  hud.legend.appendChild(bar);
  const key = document.createElement('div');
  key.className = 'v3-colorbar';
  key.innerHTML = '<div class="v3-key"><i style="border-top-style:solid"></i>lo sheet (opaque) + iso</div><div class="v3-key"><i style="border-top-style:dashed"></i>hi sheet (translucent) + iso</div>'
    + `<div class="v3-key">floor: lo sheet heat map${nTc ? ` · ${nTc} T_c points` : ''}</div>`;
  hud.legend.appendChild(key);

  const gS = hud.group('sheets');
  hud.toggle(gS, [['both', 'lo + hi'], ['lo', 'lo'], ['hi', 'hi']], 'both', (v) => {
    lo.visible = v !== 'hi'; hi.visible = wire.visible = v !== 'lo';
    if (Li) Li.visible = v !== 'hi'; if (Hi) Hi.visible = v !== 'lo';
    stage.requestRender();
  });

  const c = [WX / 2, HY * 0.35, WZ / 2];
  return {
    root,
    defaultView: { pos: [WX * 1.75, HY * 3.0, WZ * 2.9], target: c },
    bookmarks: [
      { name: 'Overview', pos: [WX * 1.75, HY * 3.0, WZ * 2.9], target: c },
      { name: 'T axis (front)', pos: [WX / 2, HY * 0.6, WZ * 2.6], target: [WX / 2, HY * 0.5, WZ / 2] },
      { name: 'Floor map (top)', pos: [WX / 2 + 0.2, HY * 4.2, WZ / 2 + 2], target: [WX / 2, 0, WZ / 2] },
      { name: 'Param axis (side)', pos: [WX * 2.1, HY * 0.6, WZ / 2], target: [WX / 2, HY * 0.5, WZ / 2] },
    ],
    noScaleBar: true,
    extraNotes: [o.sheet_note],
  };
}
