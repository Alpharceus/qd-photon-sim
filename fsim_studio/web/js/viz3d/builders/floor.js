// Enclosure floor: the breadboard the device is posted on.  An anodized
// plate with the M6 hole grid (the world's ground) and one post from the
// floor to the bottom of the device.  Decoration only: the hole pitch is
// set from the device size, never a length scale, and the note says so.
import * as THREE from 'three';

let _holeTex = null;
function holeTexture() {
  if (_holeTex) return _holeTex;
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const g = c.getContext('2d');
  g.fillStyle = '#16191d';
  g.fillRect(0, 0, 64, 64);
  // countersunk M6 hole: rim highlight, dark bore
  g.beginPath(); g.arc(32, 32, 9, 0, Math.PI * 2); g.fillStyle = '#23272d'; g.fill();
  g.beginPath(); g.arc(32, 32, 6.5, 0, Math.PI * 2); g.fillStyle = '#050607'; g.fill();
  _holeTex = new THREE.CanvasTexture(c);
  _holeTex.colorSpace = THREE.SRGBColorSpace;
  _holeTex.wrapS = _holeTex.wrapT = THREE.RepeatWrapping;
  _holeTex.anisotropy = 4;
  return _holeTex;
}

let _fade = null;
function fadeTexture() {
  if (_fade) return _fade;
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d');
  const gr = g.createRadialGradient(64, 64, 10, 64, 64, 64);
  gr.addColorStop(0, '#fff'); gr.addColorStop(0.55, '#bbb'); gr.addColorStop(1, '#000');
  g.fillStyle = gr; g.fillRect(0, 0, 128, 128);
  _fade = new THREE.CanvasTexture(c);
  return _fade;
}

// Solid extent of the device: opaque or depth-writing meshes only (glows,
// light cones, fat lines and CSS labels are left out).
export function solidBox(root) {
  const box = new THREE.Box3();
  root.updateMatrixWorld(true);
  const tmp = new THREE.Box3();
  root.traverse((o) => {
    if (!o.isMesh || !o.visible || o.isInstancedMesh) return;
    const m = Array.isArray(o.material) ? o.material[0] : o.material;
    if (!m || m.isLineMaterial || m.blending === THREE.AdditiveBlending || m.depthWrite === false) return;
    if (!o.geometry.boundingBox) o.geometry.computeBoundingBox();
    tmp.copy(o.geometry.boundingBox).applyMatrix4(o.matrixWorld);
    box.union(tmp);
  });
  return box;
}

export function addFloor(stage, root) {
  const box = solidBox(root);
  if (box.isEmpty()) return null;
  const size = new THREE.Vector3(); box.getSize(size);
  const c = new THREE.Vector3(); box.getCenter(c);
  const span = Math.max(size.x, size.z, 1e-6);
  const gap = Math.max(size.y * 0.55, span * 0.12);
  const y0 = box.min.y - gap;
  const group = new THREE.Group();
  group.name = 'enclosure-floor';

  const W = span * 7;
  const pitch = span / 4;
  const tex = holeTexture().clone();
  tex.needsUpdate = true;
  tex.repeat.set(W / pitch, W / pitch);
  const plate = new THREE.Mesh(new THREE.PlaneGeometry(W, W),
    new THREE.MeshStandardMaterial({ map: tex, alphaMap: fadeTexture(), transparent: true, roughness: 0.82, metalness: 0.25, depthWrite: true }));
  plate.rotation.x = -Math.PI / 2;
  // offset half a pitch so a hole centre sits exactly under the post
  plate.position.set(c.x + pitch / 2, y0, c.z + pitch / 2);
  plate.renderOrder = -1;
  group.add(plate);

  // Post: anodized pillar on the floor, snapped to a hole, carrying the device.
  const r = Math.max(0.08 * span, 0.3 * Math.min(size.x, size.z));
  const post = new THREE.Mesh(new THREE.CylinderGeometry(r, r * 1.08, gap, 40),
    new THREE.MeshPhysicalMaterial({ color: 0x23262b, metalness: 0.7, roughness: 0.45 }));
  post.position.set(c.x, y0 + gap / 2, c.z);
  group.add(post);
  const foot = new THREE.Mesh(new THREE.CylinderGeometry(r * 1.6, r * 1.6, gap * 0.06, 40),
    new THREE.MeshPhysicalMaterial({ color: 0x1b1e22, metalness: 0.6, roughness: 0.5 }));
  foot.position.set(c.x, y0 + gap * 0.03, c.z);
  group.add(foot);

  root.add(group);
  return { group, note: 'enclosure floor (M6 hole grid) and post are scenery, not to scale; the device above is true scale' };
}
