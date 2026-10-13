// Result export (R2 finding / devlib-r2-fix2): the numbers on a shelf leave the app as JSON or CSV with the spec,
// the provenance tag of every number and the code versions in the header. Pure data -> text; one download helper.
// JSON: {schema, kind, name, exported_at, versions, tag_legend, spec, ...payload}. CSV: "# key: value" header lines,
// then a header row and the rows; every value column is followed by its provenance tag column.
import { h } from "./dom.js";
import { getJSON } from "../api.js";
import { toast } from "./controls.js";

export const EXPORT_SCHEMA = "fsim-export/1";
export const TAG_LEGEND = { V: "verified against the cited paper", DR: "derived", E: "estimate / class range", A: "assumption" };

let versionsP = null;
/** {fsim_studio, fsim_core_hash} from /api/health; cached; {} when the server cannot answer. */
export function versions() {
  if (!versionsP) versionsP = getJSON("/api/health").then((r) => ({ fsim_studio: r.version, fsim_core_hash: r.fsim_core_hash })).catch(() => ({}));
  return versionsP;
}

const csvCell = (v) => {
  if (v === null || v === undefined) return "";
  const s = typeof v === "number" ? (Number.isFinite(v) ? String(v) : "") : String(v);
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};

/** columns: [{key, label}]; rows: array of objects. Returns CSV text (header comment lines first). */
export function toCsv({ header = {}, columns, rows }) {
  const lines = Object.entries(header).map(([k, v]) => `# ${k}: ${typeof v === "object" ? JSON.stringify(v) : v}`);
  lines.push(columns.map((c) => csvCell(c.label ?? c.key)).join(","));
  for (const r of rows) lines.push(columns.map((c) => csvCell(r[c.key])).join(","));
  return `${lines.join("\n")}\n`;
}

export function download(filename, text, mime) {
  const url = URL.createObjectURL(new Blob([text], { type: `${mime};charset=utf-8` }));
  const a = h("a", { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 4000);
}

export function jsonHeader({ kind, name, vers, spec }) {
  return { schema: EXPORT_SCHEMA, kind, name, exported_at: new Date().toISOString(), versions: vers, tag_legend: TAG_LEGEND, spec };
}

/**
 * exportResult({format, kind, name, spec, payload, table})
 *   payload: extra JSON members (result, provenance ...); table: {columns, rows} for CSV (tags as columns).
 */
export async function exportResult({ format, kind, name, spec, payload, table }) {
  const vers = await versions();
  const stamp = name || kind;
  if (format === "json") {
    download(`${stamp}.json`, `${JSON.stringify({ ...jsonHeader({ kind, name: stamp, vers, spec }), ...payload }, null, 1)}\n`, "application/json");
  } else {
    if (!table) { toast("Nothing to export as CSV yet."); return; }
    download(`${stamp}.csv`, toCsv({ header: { schema: EXPORT_SCHEMA, kind, name: stamp, versions: vers, tags: "V verified, DR derived, E estimate, A assumption", spec_kind: spec ? ["material", "structure", "drive", "setup"].map((a) => `${a}=${spec[a]?.kind}`).join(" ") : "" }, ...table }), "text/csv");
  }
  toast(`Exported ${format.toUpperCase()}.`, { kind: "ok" });
}

/** Two buttons, JSON and CSV; `get()` returns the export arguments (without format) or null when there is nothing yet. */
export function exportButtons(get, { label = "Export" } = {}) {
  const mk = (format) => h("button.btn.btn-outline.export-btn", {
    type: "button", dataset: { export: format }, title: `Download the result on the shelf as ${format.toUpperCase()} (spec, tags and versions included)`,
    "aria-label": `${label} as ${format.toUpperCase()}`,
    on: { click: async () => { const a = get(); if (!a) { toast("Run first; there is nothing to export yet."); return; } try { await exportResult({ ...a, format }); } catch (e) { toast(`Export failed: ${e.message}`, { kind: "error" }); } } },
  }, format.toUpperCase());
  const root = h("div.export-bar", { role: "group", "aria-label": `${label} the result` }, h("span.engrave.export-k", label), mk("json"), mk("csv"));
  return { el: root, setEnabled: (on) => { for (const b of root.querySelectorAll("button")) b.disabled = !on; } };
}
