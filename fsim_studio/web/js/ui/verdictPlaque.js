// Verdict plaque: word + reason + n/N, read from a parsed VERDICT line (out/*/results.md or
// verdict.md via /api/campaigns). Never recomputed here. Status colour is never alone: every
// state carries an icon and its word. Optical PASS with hardware FAIL renders as a split badge.
import { h } from "./dom.js";
import { icon } from "./icons.js";
import { lineStatus } from "./verdictBadge.js";

const REASON_KEYS = ["g2_min", "flux_max", "idealized_status", "paired_optical_pass", "hardware_qualified",
  "quality_pass", "regime", "family", "strain_bound"];
const DENOM_KEYS = ["eligible", "headline_coverage_pulsed", "headline_coverage", "coverage"];

// Status comes only from the parsed word and the line's own idealized_status (never recomputed
// from counts): split only for the word pass_hardware_infeasible.
function statusOf(word) {
  return lineStatus(word).status;
}

function humanWord(word) {
  return String(word || "").replace(/_/g, " ").toUpperCase();
}

/**
 * verdict: {word, fields, raw, line} | null; source: "out/rt_edge/verdict.md"; extra: {note, count}
 */
export function verdictPlaque(verdict, { source = "", note = "", count = 0, campaignId = "" } = {}) {
  if (!verdict) {
    return h("section.plaque.plaque-none", { "aria-label": "Campaign verdict" },
      h("div.plaque-word", icon("minus", { size: 18 }), h("span", "NO CAMPAIGN VERDICT")),
      h("p.plaque-reason", note || "No committed sweep in out/ covers this platform. Run results below are single designs, not a verdict."),
    );
  }
  const f = verdict.fields || {};
  const status = statusOf(verdict.word);
  const idealSplit = status === "fail" && String(f.idealized_status || "").toLowerCase() === "pass_hardware_infeasible";
  const denomKey = DENOM_KEYS.find((k) => f[k] != null && String(f[k]).includes("/"))
    || DENOM_KEYS.find((k) => f[k] != null);
  const idealShown = statusOf(verdict.word) === "fail" && String(f.idealized_status || "").toLowerCase() === "pass_hardware_infeasible";
  const reasons = REASON_KEYS.filter((k) => f[k] != null && k !== denomKey && !(idealShown && k === "idealized_status")).slice(0, 3)
    .map((k) => h("span.plaque-kv", h("span.k", k), h("span.v", f[k])));

  let wordEl;
  if (status === "split") {
    wordEl = h("div.plaque-word.split",
      h("span.split-half.is-pass", icon("check", { size: 16 }), "OPTICAL PASS"),
      h("span.split-half.is-fail", icon("cross", { size: 16 }), "HARDWARE FAIL"));
  } else {
    const ico = status === "pass" ? "check" : status === "fail" ? "cross" : "info";
    wordEl = h("div.plaque-word", icon(ico, { size: 18 }), h("span", humanWord(verdict.word)));
  }
  const ideal = idealSplit ? h("div.plaque-ideal", { title: "idealized_status on this line" },
    h("span.k", "idealized"),
    h("span.split-half.is-pass", icon("check", { size: 14 }), "optical PASS"),
    h("span.split-half.is-fail", icon("cross", { size: 14 }), "hardware FAIL")) : null;
  return h("section.plaque", { class: `plaque-${status}`, "aria-label": "Campaign verdict", dataset: { status } },
    h("div.plaque-top", wordEl,
      denomKey ? h("div.plaque-denom", h("span.k", denomKey.replace(/_/g, " ")), h("span.v", f[denomKey])) : null),
    ideal,
    h("div.plaque-reasons", reasons),
    h("p.plaque-source",
      `VERDICT line ${verdict.line ?? ""} · ${source}`,
      count > 1 ? ` · 1 of ${count} lines for this platform` : "",
      note ? ` · ${note}` : ""),
  );
}
