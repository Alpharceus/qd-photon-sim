import { parseRich } from "./math.js";

// ONE place for display units and notation (labels only, no physics). Used by paramForm, readout and viz3d.
// Registry ids stay ASCII ("ueV", "1/ns", "kg/m^3"); the screen shows µeV, ns⁻¹, kg m⁻³, and nothing at all for a dimensionless "1".
const LEGACY = { um2: "µm²",  cm2: "cm⁻²", "1/cm2": "cm⁻²", "/s": "s⁻¹", "um^2": "µm²" };
const ATOM = { ueV: "µeV", um: "µm", uA: "µA", um2: "µm²", ohm: "Ω", pi: "π", deg: "°", us: "µs", uW: "µW" };
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const sup = (n) => String(n).split("").map((c) => SUP[c] ?? c).join("");
const NB = "\u2009"; // thin space between unit factors

function atom(a) {
  const m = String(a).match(/^([A-Za-zµΩ°]+?)(?:\^(-?\d+))?$/);
  if (!m) return { t: a, e: 1 };
  const exp = m[2] != null ? Number(m[2]) : 1;
  return { t: ATOM[m[1]] ?? m[1], e: exp };
}

export function displayUnit(u) {
  if (u == null) return "";
  const s = String(u).trim();
  if (s === "" || s === "1" || s === "enum") return "";
  if (Object.prototype.hasOwnProperty.call(LEGACY, s)) return LEGACY[s];
  if (!/^[A-Za-zµΩ°0-9^/ .*-]+$/.test(s)) return s;
  const parts = s.split("/").map((p) => p.trim()).filter(Boolean);
  const out = [];
  parts.forEach((p, i) => {
    for (const a of p.split(/[\s*.]+/).filter(Boolean)) {
      if (a === "1") continue;
      const { t, e } = atom(a);
      const ex = i === 0 ? e : -e;
      out.push(ex === 1 ? t : `${t}${sup(ex)}`);
    }
  });
  return out.join(NB);
}

/** Notation in prose notes: Z^2 -> Z², "x Z" -> "×Z" (a multiplication sign, not the letter x). `$...$` math is left to KaTeX. */
export function displayNotation(t) {
  const fix = (x) => x
    .replace(/\^2/g, "²").replace(/\^3/g, "³")
    .replace(/\s[x*]\s(?=[A-Za-zΑ-ω(])/g, " ×")
    .replace(/\bueV\b/g, "µeV").replace(/\b(\d+(?:\.\d+)?)\s?um\b/g, "$1 µm");
  return parseRich(t).map((seg) => (seg.t === "text" ? fix(seg.v).replace(/\$/g, "\\$") : seg.t === "display" ? `$$${seg.v}$$` : `$${seg.v}$`)).join("");
}
