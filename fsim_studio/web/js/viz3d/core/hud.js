// DOM heads-up display layered over a 3D stage: title + catalog id, tagged
// readouts, mandatory labels, exaggeration badges, controls, tooltip,
// legend and banner slots.  Every readout carries its provenance chip; a
// value without a tag is not rendered.
import { tagChip, esc } from './labels.js';

export function fmtNum(v, unit = '') {
  if (v === null || v === undefined || (typeof v === 'number' && !Number.isFinite(v))) return 'n/a';
  if (typeof v === 'boolean') return v ? 'yes' : 'no';
  if (typeof v !== 'number') return String(v);
  if (Number.isInteger(v) && Math.abs(v) < 1e6) return String(v);
  const a = Math.abs(v);
  let s;
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) {
    const [m, e] = v.toExponential(2).split('e');
    s = `${m}e${Number(e)}`;
  } else if (a >= 100) s = v.toFixed(unit === 'K' || unit === 'nm' ? 1 : 0);
  else if (a >= 1) s = v.toFixed(3).replace(/0+$/, '').replace(/\.$/, '');
  else s = v.toPrecision(3);
  return s;
}

function el(tag, cls, html) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html !== undefined) e.innerHTML = html;
  return e;
}

export class Hud {
  constructor(root, spec) {
    this.root = root;
    const h = el('div', 'v3-hud');
    root.appendChild(h);
    this.el = h;
    const head = el('div', 'v3-head');
    head.appendChild(el('div', 'v3-title', esc(spec.title || spec.kind)));
    const id = el('div', 'v3-catalog');
    id.textContent = spec.catalog_id || '';
    head.appendChild(id);
    h.appendChild(head);
    this.head = head;
    // Right column: instrument shelf (readouts) with the legend under it.
    this.right = el('div', 'v3-right');
    h.appendChild(this.right);
    this.readouts = el('div', 'v3-readouts');
    this.right.appendChild(this.readouts);
    this.notes = el('div', 'v3-notes');
    h.appendChild(this.notes);
    this.controls = el('div', 'v3-controls');
    h.appendChild(this.controls);
    this.legend = el('div', 'v3-legend');
    this.right.appendChild(this.legend);
    this.tip = el('div', 'v3-tip');
    h.appendChild(this.tip);
    this.banner = null;
  }

  setReadouts(scalars, keys) {
    this.readouts.innerHTML = '';
    const list = keys || Object.keys(scalars || {});
    for (const k of list) {
      const s = scalars?.[k];
      if (!s || !s.tag) continue;
      const row = el('div', 'v3-readout');
      row.dataset.key = k;
      const lab = el('span', 'v3-readout-label');
      lab.textContent = s.label || k;
      const val = el('span', 'v3-readout-value');
      val.textContent = fmtNum(s.value, s.unit);
      if (s.note && (s.value === null || s.value === undefined)) val.title = s.note;
      const unit = el('span', 'v3-readout-unit');
      unit.textContent = s.unit || '';
      row.append(lab, val, unit, tagChip(s.tag));
      this.readouts.appendChild(row);
    }
    // Compact viewports show the first 4 readouts (tiny ones the first 2);
    // the rest open on demand.
    this.moreBtn = null;
    if (this.readouts.querySelectorAll('.v3-readout').length > 2) {
      const b = el('button', 'v3-more-btn');
      b.type = 'button';
      b.addEventListener('click', () => {
        this.root.classList.toggle('v3-readouts-open');
        this._syncMore();
      });
      this.readouts.appendChild(b);
      this.moreBtn = b;
    }
    this.setRowLimit(this.rowLimit ?? 4);
  }

  // Rows past `n` collapse behind the "+N more" button (compact only).
  setRowLimit(n) {
    this.rowLimit = n;
    this.readouts.querySelectorAll('.v3-readout').forEach((r, i) => r.classList.toggle('v3-more', i >= n));
    this._syncMore();
  }

  _syncMore() {
    const b = this.moreBtn;
    if (!b) return;
    const extra = this.readouts.querySelectorAll('.v3-readout.v3-more').length;
    b.hidden = extra === 0;
    b.textContent = this.root.classList.contains('v3-readouts-open') ? 'fewer' : `+${extra} more`;
  }

  updateReadout(key, value) {
    const v = this.readouts.querySelector(`[data-key="${key}"] .v3-readout-value`);
    if (v) v.textContent = fmtNum(value);
  }

  // Badges stay visible; the caveat lines live in one collapsible "notes"
  // chip (open by default on large viewports, closed when compact).
  setNotes(lines, badges = [], { open = true } = {}) {
    this.notes.innerHTML = '';
    for (const b of badges) {
      const e = el('div', 'v3-badge');
      e.textContent = b;
      this.notes.appendChild(e);
    }
    if (!lines || !lines.length) return;
    const d = el('details', 'v3-notes-box');
    d.open = open;
    const sm = el('summary', 'v3-notes-chip');
    sm.textContent = `notes (${lines.length})`;
    d.appendChild(sm);
    for (const line of lines) {
      const e = el('div', 'v3-note');
      e.textContent = line;
      d.appendChild(e);
    }
    this.notes.appendChild(d);
    this.notesBox = d;
  }

  // Pin key: on compact viewports every callout is a numbered pin and the
  // full text lives in this list (bottom left, collapsed by default).
  setKey(pins) {
    this.keyBox?.remove();
    this.keyRows = new Map();
    if (!pins.length) return;
    const d = el('details', 'v3-notes-box v3-key-box');
    const sm = el('summary', 'v3-notes-chip');
    sm.textContent = `labels (${pins.length})`;
    d.appendChild(sm);
    const list = el('ol', 'v3-key-list');
    for (const o of pins) {
      const li = el('li', 'v3-key-row');
      const n = el('span', 'v3-key-num');
      n.textContent = String(o.userData.pin);
      const t = el('span', 'v3-key-text');
      t.textContent = o.userData.text;
      li.append(n, t);
      if (o.userData.tag) li.appendChild(tagChip(o.userData.tag));
      list.appendChild(li);
      this.keyRows.set(o, li);
    }
    d.appendChild(list);
    this.notes.appendChild(d);
    this.keyBox = d;
  }

  // Rows for labels whose object is hidden (e.g. the other stack) drop out.
  syncKey() {
    if (!this.keyRows) return;
    for (const [o, li] of this.keyRows) {
      const off = o.element.style.display === 'none';
      if (li.hidden !== off) li.hidden = off;
    }
  }

  addBadge(text) {
    const e = el('div', 'v3-badge');
    e.textContent = text;
    this.notes.prepend(e);
    return e;
  }

  // `sub` is a neutral second line (e.g. passes of a model that never
  // gates); it never takes the pass style.
  setBanner(text, kind = 'fail', sub = '') {
    if (!this.banner) {
      this.banner = el('div', 'v3-banner');
      this.el.appendChild(this.banner);
    }
    this.banner.className = `v3-banner v3-banner-${kind}`;
    this.banner.innerHTML = `<span class="v3-banner-icon" aria-hidden="true"></span><span class="v3-banner-text"><span class="v3-banner-main">${esc(text)}</span>`
      + (sub ? `<span class="v3-banner-sub">${esc(sub)}</span>` : '') + '</span>';
  }

  // Split verdict: optical PASS with hardware FAIL renders as two halves.
  setVerdict(v) {
    if (!v) return;
    const box = el('div', 'v3-verdict');
    const half = (word, pass, what) => {
      const s = el('span', `v3-verdict-half ${pass ? 'is-pass' : 'is-fail'}`);
      s.innerHTML = `<span class="v3-verdict-icon" aria-hidden="true"></span><b>${word}</b> ${esc(what)}`;
      return s;
    };
    if (v.optical_pass !== null && v.optical_pass !== undefined) box.appendChild(half(v.optical_pass ? 'PASS' : 'FAIL', v.optical_pass, 'optical'));
    if (v.device_pass !== null && v.device_pass !== undefined) box.appendChild(half(v.device_pass ? 'PASS' : 'FAIL', v.device_pass, 'device gate'));
    if (box.children.length) this.head.appendChild(box);
  }

  group(label) {
    const g = el('div', 'v3-ctl-group');
    if (label) g.appendChild(el('span', 'v3-ctl-label', esc(label)));
    this.controls.appendChild(g);
    return g;
  }

  button(group, text, onClick, { title = '', key = '' } = {}) {
    const b = el('button', 'v3-btn');
    b.type = 'button';
    b.innerHTML = key ? `<kbd>${esc(key)}</kbd><span class="v3-btn-text">${esc(text)}</span>` : esc(text);
    if (title) b.title = title;
    b.addEventListener('click', onClick);
    group.appendChild(b);
    return b;
  }

  toggle(group, options, value, onChange) {
    const seg = el('div', 'v3-seg');
    const btns = options.map(([val, text]) => {
      const b = el('button', 'v3-seg-btn');
      b.type = 'button';
      // "optical (waveguide.py)": the parenthetical drops on tiny viewports
      // (full text stays in the tooltip).
      const m = /^(.*?)\s+(\(.*\))$/.exec(text);
      if (m) {
        b.innerHTML = `${esc(m[1])}<span class="v3-seg-sub">&nbsp;${esc(m[2])}</span>`;
        b.title = text;
      } else b.textContent = text;
      b.setAttribute('aria-pressed', String(val === value));
      b.addEventListener('click', () => {
        btns.forEach((x) => x.setAttribute('aria-pressed', 'false'));
        b.setAttribute('aria-pressed', 'true');
        onChange(val);
      });
      seg.appendChild(b);
      return b;
    });
    group.appendChild(seg);
    return seg;
  }

  slider(group, { min = 0, max = 1, step = 0.001, value = 1, label = '', onInput }) {
    const wrap = el('label', 'v3-slider');
    if (label) wrap.appendChild(el('span', 'v3-ctl-label', esc(label)));
    const s = el('input');
    s.type = 'range'; s.min = min; s.max = max; s.step = step; s.value = value;
    s.addEventListener('input', () => onInput(Number(s.value)));
    wrap.appendChild(s);
    group.appendChild(wrap);
    return s;
  }

  showTip(info, x, y) {
    if (!info) { this.tip.style.display = 'none'; return; }
    this.tip.innerHTML = '';
    const t = el('div', 'v3-tip-title');
    t.textContent = info.title || '';
    if (info.tag) t.appendChild(tagChip(info.tag));
    this.tip.appendChild(t);
    for (const [k, v] of info.rows || []) {
      const r = el('div', 'v3-tip-row');
      r.innerHTML = `<span>${esc(k)}</span><span>${esc(v)}</span>`;
      this.tip.appendChild(r);
    }
    if (info.note) this.tip.appendChild(el('div', 'v3-tip-note', esc(info.note)));
    this.tip.style.display = 'block';
    const W = this.root.clientWidth;
    const tw = 300;
    this.tip.style.left = `${Math.min(x + 16, W - tw - 12)}px`;
    this.tip.style.top = `${y + 16}px`;
  }

  dispose() { this.el.remove(); }
}
