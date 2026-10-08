// Verdict badges for parsed VERDICT lines (/api/campaigns). Display mapping only: the word
// and idealized_status come verbatim from out/*.md; nothing is recomputed. Status colour is
// never alone (icon + word), and an optical pass that fails hardware is always a split badge.
import { h } from "./dom.js";
import { icon } from "./icons.js";

export const SPLIT_WORD = "pass_hardware_infeasible";

/** {status: fail|pass|split|info, word} for a verdict line or a bare word. */
export function lineStatus(v) {
  const word = typeof v === "string" ? v : v?.word;
  const f = (typeof v === "object" && v?.fields) || {};
  const w = String(word || "").toLowerCase();
  const ideal = String(f.idealized_status || "").toLowerCase();
  if (w === SPLIT_WORD) return { status: "split", word };
  if (w.includes("fail") || w.startsWith("no_")) return { status: "fail", word, idealSplit: ideal === SPLIT_WORD };
  if (w.includes("pass")) return { status: "pass", word };
  return { status: "info", word };
}

export const humanWord = (w) => String(w || "").replace(/_/g, " ").toUpperCase();

/** Badge for one word/line. size: "" | "sm" | "lg" */
export function verdictBadge(v, { size = "" } = {}) {
  const st = lineStatus(v);
  const cls = `vb ${size ? `vb-${size}` : ""}`;
  if (st.status === "split") {
    return h("span", { class: `${cls} vb-split`, dataset: { status: "split" }, role: "img", "aria-label": "optical pass, hardware fail" },
      h("span.vb-half.is-pass", icon("check", { size: size === "sm" ? 12 : 15 }), "optical PASS"),
      h("span.vb-half.is-fail", icon("cross", { size: size === "sm" ? 12 : 15 }), "hardware FAIL"));
  }
  const ico = st.status === "pass" ? "check" : st.status === "fail" ? "cross" : "info";
  const main = h("span", { class: `${cls} vb-${st.status}`, dataset: { status: st.status } },
    icon(ico, { size: size === "sm" ? 12 : 15 }), h("span.vb-word", humanWord(st.word)));
  if (st.idealSplit) {
    return h("span.vb-pair", main, h("span", { class: `${cls} vb-split`, dataset: { status: "split" }, title: "idealized_status=pass_hardware_infeasible on this line" },
      h("span.vb-half.is-pass", icon("check", { size: 12 }), "optical PASS"),
      h("span.vb-half.is-fail", icon("cross", { size: 12 }), "hardware FAIL")));
  }
  return main;
}

/** Short scope of a verdict line: regime / family / strain bound / screening / model. */
export function lineScope(v) {
  const f = v?.fields || {};
  const bits = [];
  if (f.family) bits.push(f.family.replace(/_/g, " "));
  if (f.regime) bits.push(f.regime === "deterministic_pair" ? "SET (deterministic pair)" : f.regime);
  if (f.strain_bound) bits.push(`${f.strain_bound}${f.bound_role ? ` (${f.bound_role.replace(/_/g, " ")})` : ""}`);
  if (f.rep_rate_hz) bits.push(`${Number(f.rep_rate_hz) / 1e6} MHz`);
  if (f.screening) bits.push(`screening ${f.screening.replace(/_/g, " ")}`);
  if (!bits.length && f.model) bits.push(f.model);
  return bits.join(" · ");
}

/** Count fields of a line, verbatim, each with its own denominator when the line prints one. */
export const COUNT_KEYS = ["eligible", "paired_optical_pass", "quality_pass", "hardware_qualified", "rti_qualified",
  "idealized_pass_count", "hardware_pass_count", "headline_coverage_pulsed", "headline_coverage", "coverage", "eligible_dedup", "invalid"];

/**
 * The line's own denominator for its counts: the N of its coverage field ("128/128" -> 128),
 * else of headline_coverage(_pulsed). Read verbatim from the printed line; null when none printed.
 */
export function lineDenominator(fields) {
  const f = fields || {};
  for (const k of ["coverage", "headline_coverage", "headline_coverage_pulsed"]) {
    const m = /^\s*\d+\s*\/\s*(\d+)\s*$/.exec(String(f[k] ?? ""));
    if (m) return { n: m[1], key: k };
  }
  return null;
}

/** Count text with its denominator: "16/768" stays; "128" becomes "128 / 192" from the line's coverage. */
export function countWithDenom(fields, key) {
  const v = String(fields?.[key] ?? "");
  if (v.includes("/")) return { text: v.replace(/\s*\/\s*/, " / "), denomFrom: null };
  const d = lineDenominator(fields);
  if (!d || key === d.key) return { text: v, denomFrom: null };
  return { text: `${v} / ${d.n}`, denomFrom: d.key };
}

/**
 * Conservative representative of a multi-line campaign: the least favourable line by its own
 * word (fail < fail with an idealized optical pass < info < split < pass), first line on ties.
 */
export function conservativeLine(lines) {
  const rank = (v) => {
    const st = lineStatus(v);
    return st.status === "fail" ? (st.idealSplit ? 1 : 0) : st.status === "info" ? 2 : st.status === "split" ? 3 : 4;
  };
  let best = null;
  for (const v of lines || []) if (!best || rank(v) < rank(best)) best = v;
  return best;
}
