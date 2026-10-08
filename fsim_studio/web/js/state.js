// Tiny store: one object, shallow-merge updates, per-key subscribers.
// Keys: theme, present, route, card, cardInfo, design, ranged, mode, run (last result),
// runState, slots {A, B, C}, compareRef, selected (inspector target).

const PREF_KEY = "fsim-studio:prefs";

function readPrefs() {
  try { return JSON.parse(localStorage.getItem(PREF_KEY) || "{}") || {}; } catch { return {}; }
}
function writePrefs(p) {
  try { localStorage.setItem(PREF_KEY, JSON.stringify(p)); } catch { /* private window: fine */ }
}

function readSlots() {
  try { return JSON.parse(sessionStorage.getItem("fsim-studio:slots") || "null"); } catch { return null; }
}

const prefs = readPrefs();

const state = {
  theme: prefs.theme === "light" ? "light" : "dark",
  present: false,
  route: { ws: "design", args: [] },
  card: null,
  cardInfo: null,
  design: null,
  ranged: {},
  mode: "point",
  run: null,
  runState: "idle",
  slots: readSlots() || { A: null, B: null, C: null },
  compareRef: "A",
  selected: null,
};

const subs = new Map();

export function get(key) { return key ? state[key] : state; }

export function set(patch) {
  const changed = [];
  for (const [k, v] of Object.entries(patch)) {
    if (state[k] !== v) { state[k] = v; changed.push(k); }
  }
  if (changed.includes("theme")) writePrefs({ ...readPrefs(), theme: state.theme });
  if (changed.includes("slots")) {
    try { sessionStorage.setItem("fsim-studio:slots", JSON.stringify(state.slots)); } catch { /* quota */ }
  }
  for (const k of changed) for (const fn of subs.get(k) || []) fn(state[k], state);
  for (const fn of subs.get("*") || []) if (changed.length) fn(changed, state);
}

export function on(key, fn) {
  if (!subs.has(key)) subs.set(key, new Set());
  subs.get(key).add(fn);
  return () => subs.get(key).delete(fn);
}
