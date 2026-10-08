// Scale honesty for every 3D view.  An axis drawn at anything other than
// true scale MUST carry a badge; makeAxis throws otherwise, so a builder
// cannot quietly exaggerate.  Scale bars pick a 1-2-5 length that fits.

export class ScaleError extends Error {}

export function makeAxis(name, ax = {}) {
  const scale = ax.scale ?? 1;
  if (!Number.isFinite(scale) || scale <= 0) throw new ScaleError(`axis ${name}: invalid scale ${scale}`);
  if (Math.abs(scale - 1) > 1e-12 && !(typeof ax.badge === 'string' && ax.badge.trim())) {
    throw new ScaleError(`axis ${name}: scale ${scale} without an exaggeration badge`);
  }
  return { name, scale, badge: ax.badge || null, label: ax.label || name, unit: ax.unit || '' };
}

// Validate a spec's axes and return the list of badges to show.
export function axesBadges(spec) {
  const badges = [];
  for (const [k, ax] of Object.entries(spec.axes || {})) {
    const a = makeAxis(k, ax);
    if (a.badge) badges.push(a.badge);
  }
  return badges;
}

// Register an exaggeration used inside a builder (e.g. a gap drawn x20).
export function exaggerate(registry, what, scale, badge) {
  makeAxis(what, { scale, badge });
  if (Math.abs(scale - 1) > 1e-12) registry.add(badge);
  return scale;
}

export function niceLength(maxLen) {
  if (!(maxLen > 0)) return 0;
  const p = Math.pow(10, Math.floor(Math.log10(maxLen)));
  for (const m of [5, 2, 1]) if (m * p <= maxLen) return m * p;
  return p;
}

// Format a length given in um with a sensible unit.
export function fmtLength(um) {
  if (um >= 1000) return `${+(um / 1000).toPrecision(3)} mm`;
  if (um >= 1) return `${+um.toPrecision(3)} µm`;
  return `${+(um * 1000).toPrecision(3)} nm`;
}

// DOM scale bar that tracks the camera: worldPerUnit converts a world
// distance into the spec's length unit (um by default).
export class ScaleBar {
  constructor(parent, { unitToUm = 1, maxPx = 140 } = {}) {
    this.el = document.createElement('div');
    this.el.className = 'v3-scalebar';
    this.el.innerHTML = '<span class="v3-scalebar-rule"></span><span class="v3-scalebar-text"></span>';
    parent.appendChild(this.el);
    this.rule = this.el.firstChild;
    this.text = this.el.lastChild;
    this.unitToUm = unitToUm;
    this.maxPx = maxPx;
  }
  // pxPerWorld: screen pixels per world unit at the focus distance;
  // worldPerSpecUnit: world units per spec length unit.
  update(pxPerWorld, worldPerSpecUnit) {
    if (!(pxPerWorld > 0)) { this.el.style.display = 'none'; return; }
    const pxPerUnit = pxPerWorld * worldPerSpecUnit;
    const len = niceLength(this.maxPx / pxPerUnit);
    this.el.style.display = '';
    this.rule.style.width = `${(len * pxPerUnit).toFixed(1)}px`;
    this.text.textContent = fmtLength(len * this.unitToUm);
  }
  dispose() { this.el.remove(); }
}
