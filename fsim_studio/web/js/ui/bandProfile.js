// Material band view for the Composer (R2 finish, FIN1): a 2D band-edge / confinement profile of the chosen dot material.
// It is drawn ONLY from the parameters the registry exposes for a material:
//   lambda_nm  -> the e-h transition energy E_X = hc/lambda        (true, the tag of lambda_nm)
//   l_z_nm     -> the well width along the growth axis, to scale   (the tag of l_z_nm)
//   l_xy_nm    -> the in-plane footprint, to scale                 (the tag of l_xy_nm)
//   E_b_meV    -> the biexciton binding energy, shown as a number
//   charge_state, dipole -> shown as words
// The registry carries NO band offsets, so the confinement shift and the barrier heights are drawn as fixed fractions of E_X
// and tagged [A]: the picture says "a well this wide, a transition this high", never "this exact heterostructure".
// (The 3D band builder, viz3d/builders/bandProfile.js, needs a full scene descriptor that exists for two cards only.)
import { h, svgEl } from "./dom.js";
import { tagChip, normTag, widest } from "./tagChip.js";
import { fmt } from "./format.js";

const HC_EV_NM = 1239.84198;          // h c in eV nm [V] CODATA
const CONF_FRAC = 0.1;                // [A] confinement shift of each carrier, as a fraction of E_X
const BARRIER_FRAC = 0.3;             // [A] barrier height of each band above the well edge, as a fraction of E_X
const Z_HALF_NM = 6;                  // common half-axis (nm) so well widths are comparable between materials

const W = 640, H = 330, PADL = 64, PADR = 28, PADT = 26, PADB = 46;

/** Numeric value of a leaf read by paramForm.readLeaf: the value, or the midpoint of a range. */
function num(leaf) {
  if (!leaf) return NaN;
  if (leaf.range) return (Number(leaf.range[0]) + Number(leaf.range[1])) / 2;
  return Number(leaf.value);
}

/**
 * bandPanel(entry, getLeaf) -> { el, sig }
 *   entry   the registry material entry (label, kind, params)
 *   getLeaf(name) -> {value, range, tag} | null  (the spec's current leaf, else the registry default)
 */
export function bandPanel(entry, getLeaf) {
  const L = (n) => getLeaf(n) || null;
  const lam = num(L("lambda_nm")), lz = num(L("l_z_nm")), lxy = num(L("l_xy_nm")), eb = num(L("E_b_meV"));
  const charge = L("charge_state")?.value ?? "X", dipole = L("dipole")?.value ?? "in_plane";
  const EX = Number.isFinite(lam) && lam > 0 ? HC_EV_NM / lam : NaN;
  const root = h("div.band-panel", { dataset: { material: entry.kind } });
  if (!Number.isFinite(EX) || !Number.isFinite(lz)) {
    root.append(h("p.bp-none", "The registry gives this material no wavelength or well width, so no band view can be drawn."));
    return { el: root, sig: "none" };
  }
  // vertical layout in eV above the hole ground level
  const conf = CONF_FRAC * EX, bar = BARRIER_FRAC * EX;
  const eLvl = EX, hLvl = 0;
  const cbWell = eLvl - conf, vbWell = hLvl + conf;
  const cbBar = cbWell + bar, vbBar = vbWell - bar;
  const yMax = cbBar + 0.12 * EX, yMin = vbBar - 0.26 * EX;   // the lower margin holds the well-width dimension
  const X = (z) => PADL + ((z + Z_HALF_NM) / (2 * Z_HALF_NM)) * (W - PADL - PADR);
  const Y = (E) => PADT + (1 - (E - yMin) / (yMax - yMin)) * (H - PADT - PADB);
  const z0 = -lz / 2, z1 = lz / 2;

  const svg = svgEl("svg", { class: "bp-svg", viewBox: `0 0 ${W} ${H}`, role: "img", preserveAspectRatio: "xMidYMid meet",
    "aria-label": `Band-edge profile of ${entry.label}: transition ${fmt(EX)} eV, well ${fmt(lz)} nm` });
  const line = (x1, y1, x2, y2, cls) => svg.append(svgEl("path", { d: `M ${x1} ${y1} L ${x2} ${y2}`, class: cls }));
  const text = (x, y, s, cls = "bp-t", anchor = "start") => { const t = svgEl("text", { x, y, class: cls, "text-anchor": anchor }); t.textContent = s; svg.append(t); return t; };

  // axes: energy ticks on the left (relative to the hole ground level), growth axis in nm along the bottom
  line(PADL, PADT, PADL, H - PADB, "bp-axis");
  line(PADL, H - PADB, W - PADR, H - PADB, "bp-axis");
  const span = yMax - yMin;
  const rawStep = span / 5, mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= rawStep) || 10 * mag;
  for (let e = Math.ceil(yMin / step) * step; e <= yMax + 1e-9; e += step) {
    line(PADL - 4, Y(e), PADL, Y(e), "bp-axis");
    text(PADL - 8, Y(e) + 3.5, fmt(Number(e.toFixed(6))), "bp-tick", "end");
  }
  text(14, PADT + 4, "E (eV)", "bp-t");
  for (const z of [-Z_HALF_NM, -Z_HALF_NM / 2, 0, Z_HALF_NM / 2, Z_HALF_NM]) {
    line(X(z), H - PADB, X(z), H - PADB + 4, "bp-axis");
    text(X(z), H - PADB + 16, fmt(z), "bp-tick", "middle");
  }
  text((PADL + W - PADR) / 2, H - 8, "growth axis z (nm)", "bp-t", "middle");

  // dot slab, band edges (barrier | well | barrier), confined levels
  svg.append(svgEl("rect", { x: X(z0), y: PADT, width: X(z1) - X(z0), height: H - PADT - PADB, class: "bp-slab" }));
  const edge = (eb1, eb2, cls) => {
    svg.append(svgEl("path", { d: `M ${X(-Z_HALF_NM)} ${Y(eb1)} L ${X(z0)} ${Y(eb1)} L ${X(z0)} ${Y(eb2)} L ${X(z1)} ${Y(eb2)} L ${X(z1)} ${Y(eb1)} L ${X(Z_HALF_NM)} ${Y(eb1)}`, class: cls }));
  };
  edge(cbBar, cbWell, "bp-edge bp-cb");
  edge(vbBar, vbWell, "bp-edge bp-vb");
  const wellX0 = X(z0), wellX1 = X(z1);
  line(wellX0 + 3, Y(eLvl), wellX1 - 3, Y(eLvl), "bp-level bp-e");
  line(wellX0 + 3, Y(hLvl), wellX1 - 3, Y(hLvl), "bp-level bp-h");
  text(X(Z_HALF_NM) - 2, Y(cbBar) - 6, "conduction band edge", "bp-t", "end");
  text(X(Z_HALF_NM) - 2, Y(vbBar) + 14, "valence band edge", "bp-t", "end");
  text(wellX1 + 6, Y(eLvl) + 4, "e ground level", "bp-t");
  text(wellX1 + 6, Y(hLvl) + 4, "h ground level", "bp-t");
  // emission: a downward arrow from the electron level to the hole level, labelled with hc/lambda
  const xm = (wellX0 + wellX1) / 2;
  line(xm, Y(eLvl) + 3, xm, Y(hLvl) - 8, "bp-photon");
  svg.append(svgEl("path", { d: `M ${xm - 4.5} ${Y(hLvl) - 13} L ${xm} ${Y(hLvl) - 5} L ${xm + 4.5} ${Y(hLvl) - 13} Z`, class: "bp-photon-head" }));
  text(xm - 10, (Y(eLvl) + Y(hLvl)) / 2 + 4, `hν = ${fmt(EX)} eV  (${fmt(lam)} nm)`, "bp-t bp-hv", "end");
  // well width, to scale
  const yd = Y(vbBar) + 11;
  line(wellX0, yd, wellX1, yd, "bp-dim");
  line(wellX0, yd - 4, wellX0, yd + 4, "bp-dim");
  line(wellX1, yd - 4, wellX1, yd + 4, "bp-dim");
  text(xm, yd + 14, `l_z = ${fmt(lz)} nm`, "bp-t", "middle");
  // the assumptions are on the picture, not hidden in a tooltip
  text(PADL + 6, PADT + 10, `[A] shift ${fmt(CONF_FRAC * 100)}% and barriers ${fmt(BARRIER_FRAC * 100)}% of E_X: the registry has no band offsets`, "bp-note");

  const tagOf = (n) => normTag(L(n)?.tag);
  const row = (k, v, tag) => h("div.bp-kv", h("span.bp-k", k), h("span.bp-v", v), tag ? tagChip(tag, { size: "sm" }) : null);
  const facts = h("div.bp-facts",
    row("Transition hν", `${fmt(EX)} eV`, tagOf("lambda_nm")),
    row("Well width l_z", `${fmt(lz)} nm`, tagOf("l_z_nm")),
    Number.isFinite(lxy) ? row("Footprint l_xy", `${fmt(lxy)} nm`, tagOf("l_xy_nm")) : null,
    Number.isFinite(eb) ? row("Biexciton binding", `${fmt(eb)} meV`, tagOf("E_b_meV")) : null,
    row("Charge state", String(charge), tagOf("charge_state")),
    row("Dipole", String(dipole).replace(/_/g, " "), tagOf("dipole")),
    row("Barrier and shift", "schematic", "A"));
  root.append(svg, facts, h("p.bp-cap", `${entry.label}: drawn from lambda, l_z and l_xy of the registry entry. Barrier heights are not in the registry.`));
  root.dataset.ex = EX.toFixed(4);
  root.dataset.lz = String(lz);
  const sig = [entry.kind, EX.toFixed(4), lz, lxy, eb, charge, dipole].join("|");
  root.dataset.sig = sig;
  root.dataset.tag = widest(["lambda_nm", "l_z_nm", "l_xy_nm"].map(tagOf)) || "A";
  return { el: root, sig };
}
