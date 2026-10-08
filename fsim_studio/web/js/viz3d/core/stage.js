// Stage: renderer, camera, controls, selective bloom, CSS2D overlay,
// render-on-demand loop, camera tweens, magnifier inset, picking, export,
// dispose.  One Stage per mounted scene; dispose() frees every GPU object.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer } from 'three/addons/renderers/CSS2DRenderer.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { ENCLOSURE } from './colormap.js';

export const BLOOM_LAYER = 1;
const PIXEL_RATIO_CAP = 1.5;
const TWEEN_MS = 600;
const HDR_GAIN = 4.0;

const easeInOutCubic = (t) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

export class Stage {
  constructor(el, { reducedMotion = false, bloom = true, ortho = false } = {}) {
    this.el = el;
    this.reducedMotion = reducedMotion;
    this.disposed = false;
    this.animators = new Set();
    this.pickables = [];
    this.onFrame = new Set();
    this.inset = null;

    const wrap = document.createElement('div');
    wrap.className = 'v3-canvas-wrap';
    el.appendChild(wrap);
    this.wrap = wrap;

    const r = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
    r.setPixelRatio(Math.min(window.devicePixelRatio || 1, PIXEL_RATIO_CAP));
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.0;
    r.outputColorSpace = THREE.SRGBColorSpace;
    r.localClippingEnabled = true;
    r.setClearColor(ENCLOSURE.bg, 1);
    r.domElement.className = 'v3-canvas';
    wrap.appendChild(r.domElement);
    this.renderer = r;
    r.domElement.addEventListener('webglcontextlost', (e) => { e.preventDefault(); this.contextLost?.(); });

    this.css = new CSS2DRenderer();
    this.css.domElement.className = 'v3-css2d';
    wrap.appendChild(this.css.domElement);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(ENCLOSURE.bg);
    const pm = new THREE.PMREMGenerator(r);
    this.envRT = pm.fromScene(new RoomEnvironment(), 0.04);
    pm.dispose();
    this.scene.environment = this.envRT.texture;
    this.scene.environmentIntensity = 0.35;
    const key = new THREE.DirectionalLight(0xfff4e8, 0.95);
    key.position.set(5, 9, 7);
    const rim = new THREE.DirectionalLight(0x9fb8ff, 0.35);
    rim.position.set(-8, 3, -6);
    this.scene.add(key, rim, new THREE.AmbientLight(0xffffff, 0.08));

    this.camera = ortho ? new THREE.OrthographicCamera(-10, 10, 6, -6, -1000, 1000)
      : new THREE.PerspectiveCamera(30, 1, 0.01, 5000);
    this.camera.position.set(8, 6, 12);

    this.controls = new OrbitControls(this.camera, this.css.domElement);
    this.controls.enableDamping = !reducedMotion;
    this.controls.dampingFactor = 0.08;
    this.controls.minPolarAngle = 0.08;
    this.controls.maxPolarAngle = Math.PI * 0.58;
    this.controls.screenSpacePanning = true;
    this.controls.addEventListener('change', () => this.requestRender());
    this.controls.addEventListener('start', () => { this.tween = null; });

    this.clipPlane = new THREE.Plane(new THREE.Vector3(-1, 0, 0), 1e6);

    // Bloom only on emissive elements: the scene renders to a HalfFloat
    // (HDR) target; markBloom() lifts emissive materials above 1.0 and the
    // bloom threshold sits at 1.0, so lit (non-emissive) surfaces, kept
    // below 1.0 by the low environment intensity, never bloom.
    this.bloomEnabled = bloom;
    this.composer = new EffectComposer(r);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloomPass = new UnrealBloomPass(new THREE.Vector2(256, 256), 0.6, 0.5, 1.35);
    this.bloomPass.enabled = bloom;
    this.composer.addPass(this.bloomPass);
    this.composer.addPass(new OutputPass());

    this.frameTimes = [];
    this._raf = 0;
    this._needs = true;
    this._last = performance.now();

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this._onMove = (e) => this._hover(e);
    this._onLeave = () => this.onHover?.(null);
    this.css.domElement.addEventListener('pointermove', this._onMove);
    this.css.domElement.addEventListener('pointerleave', this._onLeave);

    this.ro = new ResizeObserver(() => this.resize());
    this.ro.observe(wrap);
    this.resize();
  }

  get size() {
    const r = this.wrap.getBoundingClientRect();
    return { w: Math.max(1, Math.round(r.width)), h: Math.max(1, Math.round(r.height)) };
  }

  resize() {
    if (this.disposed) return;
    const { w, h } = this.size;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = `${w}px`;
    this.renderer.domElement.style.height = `${h}px`;
    this.css.setSize(w, h);
    const pr = this.renderer.getPixelRatio();
    this.composer.setPixelRatio(pr);
    this.composer.setSize(w, h);
    // half-resolution bloom (UnrealBloomPass mips start at size/2 already)
    this.bloomPass.setSize(Math.max(1, Math.round(w * pr / 2)), Math.max(1, Math.round(h * pr / 2)));
    if (this.camera.isPerspectiveCamera) {
      this.camera.aspect = w / h;
    } else {
      const hh = (this.camera.top - this.camera.bottom) / 2;
      this.camera.left = -hh * w / h; this.camera.right = hh * w / h;
    }
    this.camera.updateProjectionMatrix();
    for (const f of this.onResize || []) f(w, h);
    this.requestRender();
  }

  // Lift an emissive object's materials into HDR (x HDR_GAIN, once) so it
  // crosses the bloom threshold.
  markBloom(obj, gain = HDR_GAIN) {
    obj.traverse((o) => {
      o.layers.enable(BLOOM_LAYER);
      const mats = Array.isArray(o.material) ? o.material : [o.material];
      for (const m of mats) {
        if (!m || m.userData.hdr) continue;
        m.userData.hdr = gain;
        if (m.emissive && m.emissiveIntensity !== undefined && !m.isMeshBasicMaterial) m.emissiveIntensity *= gain;
        else if (m.color) m.color.multiplyScalar(gain);
      }
    });
  }

  addPickable(obj, info) { obj.userData.info = info; this.pickables.push(obj); }

  addAnimator(fn) { this.animators.add(fn); this.requestRender(); return () => this.animators.delete(fn); }

  requestRender() {
    this._needs = true;
    if (!this._raf && !this.disposed) this._raf = requestAnimationFrame((t) => this._frame(t));
  }

  _frame() {
    this._raf = 0;
    if (this.disposed) return;
    const now = performance.now();
    const dt = Math.min(0.1, (now - this._last) / 1000);
    this._last = now;
    let again = false;
    if (this.tween) again = this._stepTween(now) || again;
    if (this.controls.enableDamping) again = this.controls.update() || again;
    for (const fn of this.animators) again = fn(dt, now) || again;
    this.render();
    const ft = performance.now() - now;
    if (again) {
      this.frameTimes.push(ft);
      if (this.frameTimes.length > 60) {
        const avg = this.frameTimes.reduce((a, b) => a + b, 0) / this.frameTimes.length;
        // Budget: drop bloom when frames cost > 20 ms for ~1 s.
        if (avg > 20 && this.bloomEnabled) { this.bloomEnabled = false; this.onBloomDrop?.(); }
        this.frameTimes = [];
      }
      this.requestRender();
    } else {
      this.frameTimes = [];
    }
  }

  render() {
    if (this.disposed) return;
    for (const f of this.onFrame) f();
    this.bloomPass.enabled = this.bloomEnabled;
    this.composer.render();
    if (this.inset && this.inset.visible) this._renderInset();
    this.css.render(this.scene, this.camera);
    this.afterCss?.();
  }

  // Magnifier inset: {scene, camera, rect: () => {x, y, w, h} in CSS px from
  // the top-left of the stage, visible}
  setInset(inset) { this.inset = inset; this.requestRender(); }

  _renderInset() {
    const r = this.renderer;
    const { h: H } = this.size;
    const { x, y, w, h } = this.inset.rect();
    const yb = H - y - h;
    r.setScissorTest(true);
    r.setScissor(x, yb, w, h);
    r.setViewport(x, yb, w, h);
    r.setClearColor(this.inset.clear ?? 0x07080a, 1);
    r.clear(true, true, false);
    r.render(this.inset.scene, this.inset.camera);
    r.setScissorTest(false);
    const { w: W } = this.size;
    r.setViewport(0, 0, W, H);
    r.setClearColor(ENCLOSURE.bg, 1);
  }

  // Pixels per world unit at the controls target (perspective) or per the
  // ortho frustum.
  pxPerWorld() {
    const { h } = this.size;
    if (this.camera.isOrthographicCamera) {
      return h / ((this.camera.top - this.camera.bottom) / this.camera.zoom);
    }
    const d = this.camera.position.distanceTo(this.controls.target);
    const vis = 2 * d * Math.tan(THREE.MathUtils.degToRad(this.camera.fov) / 2);
    return h / vis;
  }

  setView(pos, target, { instant = false } = {}) {
    const p = new THREE.Vector3(...pos), t = new THREE.Vector3(...target);
    if (instant || this.reducedMotion) {
      this.camera.position.copy(p); this.controls.target.copy(t);
      this.camera.updateProjectionMatrix(); this.controls.update();
      this.tween = null; this.requestRender();
      return;
    }
    this.tween = { t0: performance.now(), p0: this.camera.position.clone(), q0: this.controls.target.clone(), p1: p, q1: t };
    this.requestRender();
  }

  _stepTween(now) {
    const tw = this.tween;
    const k = Math.min(1, (now - tw.t0) / TWEEN_MS);
    const e = easeInOutCubic(k);
    this.camera.position.lerpVectors(tw.p0, tw.p1, e);
    this.controls.target.lerpVectors(tw.q0, tw.q1, e);
    this.camera.lookAt(this.controls.target);
    if (k >= 1) { this.tween = null; this.controls.update(); return false; }
    return true;
  }

  _hover(e) {
    if (!this.onHover || !this.pickables.length) return;
    const rect = this.css.domElement.getBoundingClientRect();
    this.pointer.set(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(this.pointer, this.camera);
    const hits = this.raycaster.intersectObjects(this.pickables, false)
      .filter((h) => this._clipOk(h.point) && h.object.visible);
    this.onHover(hits.length ? { info: hits[0].object.userData.info, instanceId: hits[0].instanceId, x: e.clientX - rect.left, y: e.clientY - rect.top } : null);
  }

  _clipOk(p) { return this.clipPlane.distanceToPoint(p) >= -1e-9; }

  // Render a frame and read back a coarse pixel grid (used by tests to
  // prove the canvas is not blank).
  sample(n = 24) {
    this.render();
    const gl = this.renderer.getContext();
    const W = gl.drawingBufferWidth, H = gl.drawingBufferHeight;
    const px = new Uint8Array(4);
    const seen = new Set();
    let lit = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const x = Math.floor((i + 0.5) * W / n), y = Math.floor((j + 0.5) * H / n);
        gl.readPixels(x, y, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
        const lum = px[0] + px[1] + px[2];
        if (lum > 3 * 22) lit += 1;
        seen.add(`${px[0] >> 3},${px[1] >> 3},${px[2] >> 3}`);
      }
    }
    return { samples: n * n, lit, distinct: seen.size, width: W, height: H, calls: this.renderer.info.render.calls, triangles: this.renderer.info.render.triangles };
  }

  stats() {
    this.renderer.info.autoReset = false;
    this.renderer.info.reset();
    this.composer.render();
    if (this.inset && this.inset.visible) this._renderInset();
    const out = { calls: this.renderer.info.render.calls, triangles: this.renderer.info.render.triangles };
    this.renderer.info.autoReset = true;
    this.render();
    return out;
  }

  async exportPNG(scale = 2) {
    const prev = this.renderer.getPixelRatio();
    const { w, h } = this.size;
    this.renderer.setPixelRatio(scale);
    this.renderer.setSize(w, h, false);
    this.composer.setPixelRatio(scale); this.composer.setSize(w, h);
    this.render();
    const url = this.renderer.domElement.toDataURL('image/png');
    this.renderer.setPixelRatio(prev);
    this.resize();
    return url;
  }

  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    cancelAnimationFrame(this._raf);
    this.ro.disconnect();
    this.controls.dispose();
    this.css.domElement.removeEventListener('pointermove', this._onMove);
    this.css.domElement.removeEventListener('pointerleave', this._onLeave);
    const seen = new Set();
    const free = (o) => {
      if (o.geometry && !seen.has(o.geometry)) { seen.add(o.geometry); o.geometry.dispose(); }
      const mats = Array.isArray(o.material) ? o.material : [o.material];
      for (const m of mats) {
        if (!m || seen.has(m)) continue;
        seen.add(m);
        for (const v of Object.values(m)) if (v && v.isTexture) v.dispose();
        m.dispose();
      }
    };
    this.scene.traverse(free);
    this.inset?.scene?.traverse(free);
    this.envRT.dispose();
    this.composer.dispose();
    this.bloomPass.dispose?.();
    this.renderer.dispose();
    this.renderer.forceContextLoss?.();
    this.wrap.remove();
  }
}
