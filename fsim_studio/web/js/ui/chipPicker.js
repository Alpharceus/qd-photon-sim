// Chip picker: a popover grid of registry entries for one axis (devlib-p11).
// Each entry shows label, summary, refs and a status badge. An entry that does not fit the rest of
// the chain is greyed with a dashed rim and the reason (from validate) but stays selectable (choosing it makes
// the run bar say INCOMPATIBLE with that reason); an entry whose status is needs_physics reads
// "awaiting physics" and cannot be chosen. Nothing here knows a kind by name:
// the entries and their status come from GET /api/compose/registry.
import { h, clear } from "./dom.js";
import { icon } from "./icons.js";
import { popover, closePopover } from "./controls.js";
import { registryText } from "./math.js";

const STATUS_WORD = { ready: "ready", needs_physics: "awaiting physics", tier: "tier" };

/** A validator message for a greyed entry: numpy noise and the Python exception name removed, one clamp. */
function tidyReason(text) {
  const t = String(text || "").replace(/np\.(?:u?int|float)\d+\(([^)]*)\)/g, "$1").replace(/^(?:[A-Za-z]+Error:\s*)/, "").replace(/\s+/g, " ").trim();
  return t.length > 150 ? `${t.slice(0, 147)}...` : t;
}

/** First sentence-ish reference, trimmed for the card. */
function refLine(refs) {
  if (!refs || !refs.length) return "";
  const r = String(refs[0]).replace(/\s+/g, " ");
  return r.length > 130 ? `${r.slice(0, 127)}...` : r;
}

/**
 * openChipPicker(anchor, {axisLabel, entries, current, check, onPick})
 *   entries: registry entries of one axis
 *   current: the kind now selected
 *   check(kind) -> Promise<{ok, reason}>   async compatibility of swapping that entry in
 *   onPick(kind)
 */
export function openChipPicker(anchor, { axisLabel, entries, current, check, onPick }) {
  const grid = h("div.cp-grid", { role: "listbox", "aria-label": `${axisLabel} kinds` });
  const cards = new Map();
  const note = h("p.cp-note", { "aria-live": "polite" }, "Checking each entry against the rest of the chain");

  const pick = (e) => { closePopover(); onPick(e.kind); };

  for (const e of entries) {
    const awaiting = e.status === "needs_physics";
    const btn = h("button.cp-entry", {
      type: "button", role: "option", dataset: { kind: e.kind, status: e.status },
      "aria-selected": e.kind === current ? "true" : "false",
      "aria-disabled": awaiting ? "true" : null,
      on: { click: () => { if (e.status !== "needs_physics") pick(e); } },
    },
    h("span.cp-top",
      h("span.cp-label", registryText(e.label)),
      h("span.cp-badge", { class: `cp-${e.status}` }, STATUS_WORD[e.status] || e.status)),
    h("span.cp-sum", registryText(e.summary || "")),
    refLine(e.refs) ? h("span.cp-ref", registryText(refLine(e.refs))) : null,
    h("span.cp-why", { hidden: !awaiting }, awaiting ? "Registered, but the solver for it has not landed; it cannot run yet." : ""));
    if (awaiting) btn.title = "awaiting physics";
    cards.set(e.kind, btn);
    grid.append(btn);
  }

  const content = h("div.cp", h("h2.why-title", `${axisLabel}: choose a kind`), note, grid);
  const pop = popover(anchor, content, { label: `${axisLabel} kind picker` });
  pop.classList.add("cp-pop");
  const w = Math.min(760, window.innerWidth - 24);
  pop.style.width = `${w}px`;
  const r = anchor.getBoundingClientRect();
  pop.style.left = `${Math.max(12, Math.min(r.left, window.innerWidth - w - 12))}px`;
  // tall cards: use up to ~620 px, sliding up over the chain when there is not enough room below the chip
  const cap = Math.min(620, window.innerHeight - 72);
  pop.style.maxHeight = `${cap}px`;
  pop.style.top = `${Math.max(56, Math.min(r.bottom + 8, window.innerHeight - cap - 16))}px`;

  // Async compatibility: one validate per candidate (the swap is the only change).
  let open = true;
  const t = setInterval(() => { if (!pop.isConnected) { open = false; clearInterval(t); } }, 400);
  let pending = 0;
  const done = () => {
    pending -= 1;
    if (pending <= 0 && open) {
      const nBad = [...cards.values()].filter((b) => b.dataset.compat === "no").length;
      note.textContent = nBad
        ? `${nBad} greyed: they do not fit the rest of the chain as it stands (hover for the reason). You can still pick one; the run bar then says INCOMPATIBLE.`
        : "Every ready entry fits the rest of the chain as it stands.";
    }
  };
  for (const e of entries) {
    if (e.status === "needs_physics") continue;
    pending += 1;
    Promise.resolve(check(e.kind)).then((res) => {
      if (!open) return;
      const btn = cards.get(e.kind);
      btn.dataset.compat = res.ok ? "yes" : "no";
      if (!res.ok && e.kind !== current) {
        btn.setAttribute("aria-disabled", "true");   // greyed, not blocked: the click handler lets it through
        const badge = btn.querySelector(".cp-badge");
        if (badge) badge.hidden = true;   // a greyed entry shows only its reason, never a READY badge
        const why = btn.querySelector(".cp-why");
        why.hidden = false;
        why.replaceChildren(icon("warn", { size: 13 }), document.createTextNode(` ${tidyReason(res.reason)}`));
        btn.title = res.reason;
      }
    }).catch(() => {}).finally(done);
  }
  if (!pending) done();
  // focus the current entry for keyboard users
  queueMicrotask(() => (cards.get(current) || grid.querySelector("button"))?.focus());
  return pop;
}

/** Clear helper for tests. */
export function closeChipPicker() { closePopover(); }
void clear;
