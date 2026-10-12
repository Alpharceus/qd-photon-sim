// ONE place for display units and notation (labels only, no physics). Used by paramForm, readout and viz3d.
// Registry ids stay ASCII ("ueV", "um"); the screen shows µeV, µm, and nothing at all for a dimensionless "1".
const MAP = { "1": "", ueV: "µeV", um: "µm", um2: "µm²", "um^2": "µm²", uA: "µA", cm2: "cm⁻²", "1/cm2": "cm⁻²", "/s": "s⁻¹" };

export function displayUnit(u) {
  if (u == null) return "";
  const s = String(u).trim();
  if (Object.prototype.hasOwnProperty.call(MAP, s)) return MAP[s];
  return s.replace(/\bueV\b/g, "µeV").replace(/\bum\b/g, "µm").replace(/\^-1/g, "⁻¹").replace(/\^-2/g, "⁻²").replace(/\^-3/g, "⁻³").replace(/\^2/g, "²").replace(/\^3/g, "³");
}

/** Notation in prose notes: Z^2 -> Z², "x Z" -> "×Z" (a multiplication sign, not the letter x). */
export function displayNotation(t) {
  return String(t == null ? "" : t)
    .replace(/\^2/g, "²").replace(/\^3/g, "³")
    .replace(/\s[x*]\s(?=[A-Za-zΑ-ω(])/g, " ×")
    .replace(/\bueV\b/g, "µeV").replace(/\b(\d+(?:\.\d+)?)\s?um\b/g, "$1 µm");
}
