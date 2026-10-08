// Shared mesh helpers for the builders (materials by role, break marks,
// axis ticks, magnifier inset).  Geometry is built from SceneSpec numbers.
import * as THREE from 'three';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { ROLE_TINT } from '../core/colormap.js';

export function layerMaterial(role, { clip = null, opacity = null, emissive = 0x000000, emissiveIntensity = 0 } = {}) {
  const transparentRoles = { cladding: 0.2, oxide: 0.35, thermal: 0.45, mirror: 0.75 };
  const op = opacity ?? transparentRoles[role] ?? 1;
  const solid = op >= 1;
  const m = new THREE.MeshPhysicalMaterial({
    color: ROLE_TINT[role] ?? 0x888888,
    roughness: solid ? 0.42 : 0.6, metalness: solid ? 0.12 : 0.0,
    clearcoat: solid ? 0.25 : 0.0, clearcoatRoughness: 0.5,
    transparent: !solid, opacity: op, depthWrite: solid,
    emissive, emissiveIntensity,
    side: solid ? THREE.DoubleSide : THREE.FrontSide,
  });
  if (clip) m.clippingPlanes = [clip];
  return m;
}

// Machined edge outline for a mesh (thin, low-contrast).
export function outline(mesh, { color = 0xaeb4bd, opacity = 0.35, clip = null } = {}) {
  const l = new THREE.LineSegments(new THREE.EdgesGeometry(mesh.geometry, 20),
    new THREE.LineBasicMaterial({ color, transparent: true, opacity, depthWrite: false }));
  if (clip) l.material.clippingPlanes = [clip];
  mesh.add(l);
  return l;
}

export function glowMaterial(rgb, { opacity = 1, additive = true, map = null } = {}) {
  return new THREE.MeshBasicMaterial({
    color: new THREE.Color(rgb[0], rgb[1], rgb[2]),
    transparent: true, opacity, depthWrite: false, map,
    blending: additive ? THREE.AdditiveBlending : THREE.NormalBlending,
    side: THREE.DoubleSide, toneMapped: false,
  });
}

// Light cone that fades from the apex (alpha 1) to the base (alpha 0).
let _fade = null;
export function coneMaterial(rgb, opacity = 0.22) {
  if (!_fade) {
    const c = document.createElement('canvas');
    c.width = 4; c.height = 128;
    const g = c.getContext('2d');
    const gr = g.createLinearGradient(0, 0, 0, 128);
    gr.addColorStop(0, '#fff'); gr.addColorStop(0.35, '#777'); gr.addColorStop(1, '#000');
    g.fillStyle = gr; g.fillRect(0, 0, 4, 128);
    _fade = new THREE.CanvasTexture(c);
  }
  return new THREE.MeshBasicMaterial({
    color: new THREE.Color(rgb[0], rgb[1], rgb[2]), alphaMap: _fade, transparent: true, opacity,
    depthWrite: false, blending: THREE.AdditiveBlending, side: THREE.DoubleSide,
  });
}

export function box(w, h, d, mat, [x, y, z]) {
  const m = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), mat);
  m.position.set(x, y, z);
  return m;
}

// Fat line in world units (LineMaterial width in px).
export function fatLine(points, { color = 0xe9e8e3, width = 1.6, opacity = 1, dashed = false, dashSize = 0.2, gapSize = 0.12 } = {}) {
  const g = new LineGeometry();
  g.setPositions(points.flat());
  const m = new LineMaterial({ color, linewidth: width, transparent: opacity < 1, opacity, dashed, dashSize, gapSize, worldUnits: false });
  const l = new Line2(g, m);
  if (dashed) l.computeLineDistances();
  return l;
}

export function updateLineResolution(root, w, h) {
  root.traverse((o) => { if (o.material?.isLineMaterial) o.material.resolution.set(w, h); });
}

// Zig-zag break mark spanning from a to b (world), amplitude amp, n teeth,
// offset along `normal`.
export function breakMark(a, b, amp, normal, n = 6, opts = {}) {
  const A = new THREE.Vector3(...a), B = new THREE.Vector3(...b), N = new THREE.Vector3(...normal).normalize();
  const pts = [];
  for (let i = 0; i <= n * 2; i++) {
    const t = i / (n * 2);
    const p = A.clone().lerp(B, t).addScaledVector(N, (i % 2 ? 1 : -1) * amp * (i === 0 || i === n * 2 ? 0 : 1));
    pts.push([p.x, p.y, p.z]);
  }
  return fatLine(pts, { color: 0xd8d6cf, width: 1.4, ...opts });
}

// Magnifier frame (DOM) over the stage plus an inset scene with its own
// orthographic camera, magnification computed live against the main view.
export function makeMagnifier(stage, hud, { mag, buildInset, unitLabel = 'nm', mainUnitsPerInset = 1000, title = 'dot magnifier' }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x07080a);
  scene.environment = stage.scene.environment;
  scene.environmentIntensity = 0.8;
  const l1 = new THREE.DirectionalLight(0xffffff, 1.6); l1.position.set(3, 6, 8);
  scene.add(l1, new THREE.AmbientLight(0xffffff, 0.25));
  const cam = new THREE.OrthographicCamera(-1, 1, 1, -1, -1e5, 1e5);
  const info = buildInset(scene, cam);
  const frame = document.createElement('div');
  frame.className = 'v3-inset';
  frame.innerHTML = `<div class="v3-inset-head"><span class="v3-inset-title">${title}</span><span class="v3-badge v3-inset-badge"></span></div>`
    + '<div class="v3-inset-foot"><span class="v3-scalebar"><span class="v3-scalebar-rule"></span><span class="v3-scalebar-text"></span></span><span class="v3-inset-note"></span></div>';
  hud.el.appendChild(frame);
  const badge = frame.querySelector('.v3-inset-badge');
  const rule = frame.querySelector('.v3-scalebar-rule');
  const stext = frame.querySelector('.v3-scalebar-text');
  frame.querySelector('.v3-inset-note').textContent = info.note || '';
  const rect = () => {
    const r = frame.getBoundingClientRect(), p = stage.wrap.getBoundingClientRect();
    return { x: Math.round(r.left - p.left), y: Math.round(r.top - p.top + 26), w: Math.round(r.width), h: Math.round(r.height - 52) };
  };
  const inset = { scene, camera: cam, rect, visible: false, clear: 0x07080a };
  const updateCam = () => {
    const { w, h } = rect();
    if (w < 4 || h < 4) return;
    // Main view: world units per pixel at the focus; inset shows the same
    // pixel footprint `mag` times smaller (in the inset's own units).
    const mainWorldPerPx = 1 / stage.pxPerWorld();
    const insetUnitsPerPx = (mainWorldPerPx * info.insetUnitsPerWorld) / mag;
    const halfW = (w * insetUnitsPerPx) / 2, halfH = (h * insetUnitsPerPx) / 2;
    const [cx, cy] = info.center;
    cam.left = cx - halfW; cam.right = cx + halfW; cam.top = cy + halfH; cam.bottom = cy - halfH;
    cam.position.set(cx, cy, 1000); cam.lookAt(cx, cy, 0);
    cam.updateProjectionMatrix();
    badge.textContent = `×${mag} vs main view`;
    const pxPerUnit = 1 / insetUnitsPerPx;
    let len = Math.pow(10, Math.floor(Math.log10(90 / pxPerUnit)));
    for (const m of [5, 2]) if (m * len * pxPerUnit <= 90) { len *= m; break; }
    rule.style.width = `${(len * pxPerUnit).toFixed(1)}px`;
    stext.textContent = `${+len.toPrecision(3)} ${unitLabel}`;
  };
  stage.onFrame.add(() => { if (inset.visible) updateCam(); });
  const set = (v) => { inset.visible = v; frame.style.display = v ? 'block' : 'none'; stage.setInset(inset); stage.requestRender(); };
  set(false);
  return { set, get visible() { return inset.visible; }, frame };
}

export function nice(v, d = 3) {
  if (!Number.isFinite(v)) return 'n/a';
  return +v.toPrecision(d) + '';
}
