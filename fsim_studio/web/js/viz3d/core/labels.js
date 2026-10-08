// CSS2D annotations (capped at 30 per scene) and provenance tag chips.
// Tags are graded by LINE FORM, never hue: [V] solid, [DR] solid + inset
// rule, [E] outline, [A] dashed outline.
import { CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

export const MAX_LABELS = 30;
const TAGS = { V: 'v', DR: 'dr', E: 'e', A: 'a' };

export function normTag(tag) {
  const t = String(tag ?? 'A').replace(/[[\]]/g, '').toUpperCase().split('/')[0].trim();
  return TAGS[t] ? t : 'A';
}

export function tagChip(tag) {
  const t = normTag(tag);
  const el = document.createElement('span');
  el.className = `v3-tag v3-tag-${TAGS[t]}`;
  el.textContent = t;
  el.title = { V: 'verified against the cited paper', DR: 'derived', E: 'estimate / class range', A: 'assumption' }[t];
  return el;
}

export function esc(s) {
  return String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

export class Labels {
  constructor() { this.items = []; this.dropped = 0; }

  // kind: 'chip' (engraved plate with tag), 'note' (plain engraved text),
  // 'badge' (exaggeration / magnifier badge). anchor: 'left' puts the plate
  // to the right of the point with a leader rule; 'center' centres it.
  add(parent, pos, text, { tag = null, kind = 'chip', anchor = 'left', sub = null, cls = '', priority = null } = {}) {
    if (this.items.length >= MAX_LABELS) { this.dropped += 1; return null; }
    const wrap = document.createElement('div');
    wrap.className = `v3-label v3-label-${kind} v3-anchor-${anchor} ${cls}`;
    const plate = document.createElement('div');
    plate.className = 'v3-label-plate';
    if (kind === 'chip') {
      // Numbered pin shown instead of the callout on compact viewports; the
      // number keys into the HUD's label list.
      const pin = document.createElement('span');
      pin.className = 'v3-label-pin';
      plate.appendChild(pin);
    }
    const main = document.createElement('span');
    main.className = 'v3-label-text';
    main.textContent = text;
    plate.appendChild(main);
    if (tag) plate.appendChild(tagChip(tag));
    if (sub) {
      const s = document.createElement('span');
      s.className = 'v3-label-sub';
      s.textContent = sub;
      plate.appendChild(s);
    }
    wrap.appendChild(plate);
    const obj = new CSS2DObject(wrap);
    obj.position.set(pos[0], pos[1], pos[2]);
    obj.center.set(anchor === 'left' ? 0 : anchor === 'right' ? 1 : 0.5, 0.5);
    obj.userData.kind = kind;
    obj.userData.priority = priority ?? (kind === 'chip' ? 50 : kind === 'tick' || kind === 'axis' ? 20 : 30);
    obj.userData.order = this.items.length;
    obj.userData.text = text;
    obj.userData.tag = tag;
    parent.add(obj);
    this.items.push(obj);
    return obj;
  }

  setVisible(v) { for (const o of this.items) o.visible = v; }

  // Chip labels in priority order, numbered 1..N (the pin numbers).
  pinned() {
    if (!this._pins) {
      this._pins = this.items.filter((o) => o.userData.kind === 'chip')
        .sort((a, b) => (b.userData.priority - a.userData.priority) || (a.userData.order - b.userData.order));
      this._pins.forEach((o, i) => {
        o.userData.pin = i + 1;
        const p = o.element.querySelector('.v3-label-pin');
        if (p) p.textContent = String(i + 1);
      });
    }
    return this._pins;
  }

  // Collision-aware layout, run after each CSS2D render: labels are placed
  // in priority order; a label that overlaps an already placed label or a
  // HUD panel is nudged vertically (up to two rows either way), otherwise
  // hidden (its object stays reachable through the hover tooltip).  At most
  // `budget` chip labels are shown.
  layout(rootEl, obstacleEls, budget = MAX_LABELS) {
    const R = rootEl.getBoundingClientRect();
    const placed = [];
    for (const e of obstacleEls) {
      if (!e || !e.isConnected) continue;
      const r = e.getBoundingClientRect();
      if (r.width > 0 && r.height > 0) placed.push(r);
    }
    const items = this.items.filter((o) => o.visible && o.element.style.display !== 'none')
      .sort((a, b) => (b.userData.priority - a.userData.priority) || (a.userData.order - b.userData.order));
    const hit = (r) => placed.some((p) => r.left < p.right + 2 && r.right > p.left - 2 && r.top < p.bottom + 1 && r.bottom > p.top - 1);
    const inside = (r) => r.left >= R.left + 4 && r.right <= R.right - 4 && r.top >= R.top + 4 && r.bottom <= R.bottom - 4;
    let chips = 0;
    for (const o of items) {
      const plate = o.element.firstChild;
      plate.style.transform = '';
      o.element.style.visibility = '';
      const isChip = o.userData.kind === 'chip';
      if (isChip && chips >= budget) { o.element.style.visibility = 'hidden'; continue; }
      const r0 = plate.getBoundingClientRect();
      const h = r0.height + 3;
      let ok = false;
      for (const dy of [0, h, -h, 2 * h, -2 * h]) {
        const r = { left: r0.left, right: r0.right, top: r0.top + dy, bottom: r0.bottom + dy };
        if (inside(r) && !hit(r)) {
          if (dy) plate.style.transform = `translateY(${dy}px)`;
          placed.push(r);
          ok = true;
          break;
        }
      }
      if (!ok) o.element.style.visibility = 'hidden';
      else if (isChip) chips += 1;
    }
  }

  dispose() {
    for (const o of this.items) { o.element.remove(); o.parent?.remove(o); }
    this.items = [];
  }
}
