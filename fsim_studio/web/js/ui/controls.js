// Softkey-style controls: segmented, toggle, toast, popover, table, and the
// loading / empty / error state plates.
import { h, clear } from "./dom.js";
import { icon } from "./icons.js";

/** Segmented control. options: [{value, label, sub?, disabled?, title?}] */
export function segmented({ label, options, value, onChange }) {
  const root = h("div.segmented", { role: "radiogroup", "aria-label": label });
  const btns = options.map((o) => {
    const b = h("button.seg", {
      type: "button", role: "radio", disabled: !!o.disabled, title: o.title || null,
      dataset: { value: o.value },
      on: { click: () => { if (!o.disabled) { select(o.value); onChange?.(o.value); } } },
    }, h("span.seg-label", o.label), o.sub ? h("span.seg-sub", o.sub) : null);
    return b;
  });
  root.addEventListener("keydown", (e) => {
    if (!["ArrowLeft", "ArrowRight"].includes(e.key)) return;
    const enabled = btns.filter((b) => !b.disabled);
    const i = enabled.findIndex((b) => b.getAttribute("aria-checked") === "true");
    const n = enabled[(i + (e.key === "ArrowRight" ? 1 : enabled.length - 1)) % enabled.length];
    n.click(); n.focus(); e.preventDefault();
  });
  function select(v) {
    for (const b of btns) {
      const on = b.dataset.value === v;
      b.setAttribute("aria-checked", on ? "true" : "false");
      b.tabIndex = on ? 0 : -1;
    }
  }
  select(value);
  root.append(...btns);
  return {
    el: root,
    select,
    setSub(v, text) { const b = btns.find((x) => x.dataset.value === v); if (b) { let s = b.querySelector(".seg-sub"); if (!s) { s = h("span.seg-sub"); b.append(s); } s.textContent = text; } },
    setDisabled(v, dis, title) { const b = btns.find((x) => x.dataset.value === v); if (b) { b.disabled = dis; b.title = title || ""; } },
  };
}

export function switchToggle({ label, checked, onChange, hideLabel = false }) {
  const input = h("input", { type: "checkbox", role: "switch", checked: !!checked,
    on: { change: (e) => onChange?.(e.target.checked) } });
  return h("label.switch", input, h("span.switch-track", { "aria-hidden": "true" }),
    h(hideLabel ? "span.sr-only" : "span.switch-label", label));
}

// ------------------------------------------------------------------ toast
export function toast(message, { kind = "info", timeout = 4200, action = null } = {}) {
  const host = document.getElementById("toasts");
  if (!host) return;
  const ico = { info: "info", error: "warn", ok: "check" }[kind] || "info";
  const el = h("div.toast", { class: `toast-${kind}` }, icon(ico, { size: 18 }), h("span.toast-msg", message),
    action ? h("button.btn.btn-quiet", { type: "button", on: { click: () => { action.fn(); el.remove(); } } }, action.label) : null,
    h("button.icon-btn", { type: "button", "aria-label": "Dismiss", on: { click: () => el.remove() } }, icon("close", { size: 16 })));
  host.appendChild(el);
  if (timeout) setTimeout(() => el.remove(), timeout);
  return el;
}

// ------------------------------------------------------------------ popover
let openPop = null;
export function popover(anchor, content, { label = "Details", onClose } = {}) {
  closePopover();
  const pop = h("div.popover", { role: "dialog", "aria-label": label, tabindex: "-1" },
    h("button.icon-btn.pop-close", { type: "button", "aria-label": "Close", on: { click: () => closePopover() } }, icon("close", { size: 16 })),
    content);
  document.body.appendChild(pop);
  const r = anchor.getBoundingClientRect();
  const w = Math.min(380, window.innerWidth - 24);
  pop.style.width = `${w}px`;
  pop.style.left = `${Math.max(12, Math.min(r.left, window.innerWidth - w - 12))}px`;
  pop.style.top = `${r.bottom + 8}px`;
  const away = (e) => { if (!pop.contains(e.target) && !anchor.contains(e.target)) closePopover(); };
  const esc = (e) => { if (e.key === "Escape") { e.stopPropagation(); closePopover(); anchor.focus(); } };
  setTimeout(() => document.addEventListener("pointerdown", away), 0);
  document.addEventListener("keydown", esc, true);
  openPop = { pop, cleanup: () => { document.removeEventListener("pointerdown", away); document.removeEventListener("keydown", esc, true); onClose?.(); } };
  pop.focus();
  return pop;
}
export function closePopover() {
  if (!openPop) return;
  openPop.cleanup();
  openPop.pop.remove();
  openPop = null;
}

// ------------------------------------------------------------------ table
/** columns: [{key, label, align?, fmt?}], rows: [obj] */
export function dataTable({ columns, rows, caption, className = "" }) {
  return h("div.table-wrap", { class: className },
    h("table.data-table",
      caption ? h("caption", caption) : null,
      h("thead", h("tr", columns.map((c) => h("th", { scope: "col", class: c.align === "right" ? "num" : "" }, c.label)))),
      h("tbody", rows.map((r) => h("tr", columns.map((c) => {
        const v = c.fmt ? c.fmt(r[c.key], r) : r[c.key];
        return h("td", { class: c.align === "right" ? "num" : "" }, v instanceof Node ? v : (v ?? ""));
      }))))));
}

// ------------------------------------------------------------------ state plates
export function skeleton(kind = "block", { lines = 3, height = null } = {}) {
  const el = h("div.skeleton", { class: `sk-${kind}`, "aria-hidden": "true" });
  if (height) el.style.height = `${height}px`;
  for (let i = 0; i < lines; i++) el.appendChild(h("span.sk-line"));
  return el;
}

export function emptyState({ title, body, action = null, iconName = "info" }) {
  return h("div.empty-state",
    icon(iconName, { size: 22 }),
    h("p.empty-title", title),
    body ? h("p.empty-body", body) : null,
    action ? h("button.btn", { type: "button", on: { click: action.fn } }, action.label) : null);
}

export function errorState({ title = "Run failed", message, actions = [] }) {
  return h("div.error-state", { role: "alert" },
    h("div.err-head", icon("warn", { size: 20 }), h("span.err-title", title)),
    h("p.err-msg", message),
    actions.length ? h("div.err-actions", actions.map((a) =>
      h("button.btn", { type: "button", class: a.primary ? "btn-outline" : "", on: { click: a.fn } }, a.label))) : null);
}

export function mountInto(el, ...nodes) { clear(el); el.append(...nodes.filter(Boolean)); return el; }
