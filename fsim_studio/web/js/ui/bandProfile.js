// Material band view for the Composer (R2 BAND): a 2D band-edge / confinement profile of the chosen dot material, drawn from the
// registry's literature `band` block (fsim_core/compose/data/band_offsets.yaml, display-only) plus the parameters the material exposes:
//   band.Eg_barrier, band.variants.{unstr,strained,indot}.{Eg_dot,dE_c,dE_v} -> the barrier | dot | barrier edges (each leaf tagged)
//   band.masses            -> the inverse-mass split of the confinement energy between the e and h levels
//   band.flags             -> polar (tilted dot, e and h at opposite interfaces), type2_risk (weak hole confinement), hard_wall (broken barriers)
//   lambda_nm, E_b_meV     -> E_X = hc/lambda and the single-particle gap E_X + E_b that sets the confinement energy
//   l_z_nm, l_xy_nm        -> the well width / footprint, to scale
// Which variant is drawn (caveat 1 of the answer): E_g,dot and dE_c from the in-dot value where the table gives one, else the strained one;
// dE_v from the in-dot value, else the half-way value 0.5 (unstr + str). The barrier gap is E_g,barrier, so the drawn dE_c is the residual
// E_g,barrier - E_g,dot - dE_v; when it differs from the tabulated one by more than 5 meV (mixed variants) the caption says so.
// Level split (the answer's one rule): dE = E_X + E_b - E_g,dot, E_e = dE m_h / (m_e + m_h), E_h = dE m_e / (m_e + m_h), each clamped to 0.85 of its offset.
import { h, svgEl } from "./dom.js";
import { tagChip, normTag, widest } from "./tagChip.js";
import { fmt } from "./format.js";
import { tex, richText } from "./math.js";

const HC_EV_NM = 1239.84198;          // h c in eV nm [V] CODATA
const CLAMP = 0.85;                   // [A] a level sits at most this fraction of the way up its offset
const TILT_FRAC = 0.15;               // [A] polar dots: the field tilt across the well, as a fraction of E_g,dot (schematic, not a computed field)
const HARD_CAP_EV = 0.7;              // [A] hard-wall barriers are drawn this tall and broken: "effectively infinite"
const Z_HALF_NM = 6;                  // common half-axis (nm) so well widths are comparable between materials

const W = 640, H = 340, PADL = 64, PADR = 28, PADT = 30, PADB = 46;

/** Numeric value of a leaf read by paramForm.readLeaf: the value, or the midpoint of a range. */
function num(leaf) {
  if (!leaf) return NaN;
  if (leaf.range) return (Number(leaf.range[0]) + Number(leaf.range[1])) / 2;
  return Number(leaf.value);
}

const val = (leaf) => (leaf && Number.isFinite(Number(leaf.value)) ? Number(leaf.value) : NaN);

/** Which literature variant is drawn for E_g,dot, dE_c and dE_v (see the header). Each result: {value, tag, from, where}. */
export function pickBand(band) {
  const v = band.variants || {};
  const first = (q, order) => { for (const n of order) if (v[n]?.[q]) return { ...v[n][q], from: n }; return null; };
  const Eg = first("Eg_dot", ["indot", "strained", "unstr"]);
  const dEc = first("dE_c", ["indot", "strained", "unstr"]);
  let dEv = first("dE_v", ["indot"]);
  if (!dEv) {
    const u = v.unstr?.dE_v, s = v.strained?.dE_v;
    dEv = u && s
      ? { value: 0.5 * (Number(u.value) + Number(s.value)), tag: "E", from: "half", where: "0.5 (unstrained + strained): caveat 1, the shear strain relaxes toward the apex [E]" }
      : first("dE_v", ["strained", "unstr"]);
  }
  return { Eg, dEc, dEv };
}

/**
 * bandPanel(entry, getLeaf) -> { el, sig }
 *   entry   the registry material entry (label, kind, params, band)
 *   getLeaf(name) -> {value, range, tag} | null  (the spec's current leaf, else the registry default)
 */
export function bandPanel(entry, getLeaf) {
  const L = (n) => getLeaf(n) || null;
  const lam = num(L("lambda_nm")), lz = num(L("l_z_nm")), lxy = num(L("l_xy_nm")), eb = num(L("E_b_meV"));
  const charge = L("charge_state")?.value ?? "X", dipole = L("dipole")?.value ?? "in_plane";
  const EX = Number.isFinite(lam) && lam > 0 ? HC_EV_NM / lam : NaN;
  const band = entry.band;
  const root = h("div.band-panel", { dataset: { material: entry.kind } });
  if (!Number.isFinite(EX) || !Number.isFinite(lz)) {
    root.append(h("p.bp-none", "The registry gives this material no wavelength or well width, so no band view can be drawn."));
    return { el: root, sig: "none" };
  }
  const sel = band ? pickBand(band) : null;
  const Eg = val(sel?.Eg), dEv = val(sel?.dEv) / 1000, Egb = val(band?.Eg_barrier);
  if (!sel || ![Eg, dEv, Egb, val(sel.dEc)].every(Number.isFinite)) {
    root.append(h("p.bp-none", "The registry carries no literature band offsets for this material, so no band view can be drawn."));
    return { el: root, sig: `none|${entry.kind}` };
  }
  const flags = band.flags || {};
  const dEcTab = val(sel.dEc) / 1000;
  const dEc = Egb - Eg - dEv;                         // residual: the drawn barrier gap is E_g,barrier
  const mixed = Math.abs(dEc - dEcTab) * 1000 > 5;
  // the answer's level-split rule
  const me = val(band.masses?.m_e), mh = val(band.masses?.m_h);
  const confTotal = Math.max(0, EX + (Number.isFinite(eb) ? eb / 1000 : 0) - Eg);
  const Ee = Math.min(confTotal * mh / (me + mh), CLAMP * Math.max(dEc, 0));
  const Eh = Math.min(confTotal * me / (me + mh), CLAMP * Math.max(dEv, 0));
  const tilt = flags.polar ? TILT_FRAC * Eg : 0;
  // vertical layout in eV above the dot valence-band edge; hard-wall barriers are capped and broken
  const capC = flags.hard_wall ? Math.min(dEc, HARD_CAP_EV) : dEc, capV = flags.hard_wall ? Math.min(dEv, HARD_CAP_EV) : dEv;
  const cbBar = Eg + capC, vbBar = -capV;
  const span0 = cbBar - vbBar;
  const nNotes = (flags.polar ? 1 : 0) + (flags.type2_risk ? 1 : 0);   // flag notes sit in a band of headroom above the barrier
  const yMax = cbBar + (0.1 + 0.07 * nNotes) * span0 + tilt / 2, yMin = vbBar - 0.2 * span0;   // the lower margin holds the well-width dimension
  const X = (z) => PADL + ((z + Z_HALF_NM) / (2 * Z_HALF_NM)) * (W - PADL - PADR);
  const Y = (E) => PADT + (1 - (E - yMin) / (yMax - yMin)) * (H - PADT - PADB);
  const z0 = -lz / 2, z1 = lz / 2;
  const cbAt = (u) => Eg + tilt * (0.5 - u), vbAt = (u) => tilt * (0.5 - u);   // u = 0 at the left interface, 1 at the right

  const svg = svgEl("svg", { class: "bp-svg", viewBox: `0 0 ${W} ${H}`, role: "img", preserveAspectRatio: "xMidYMid meet",
    "aria-label": `Band-edge profile of ${entry.label}: ${band.dot} in ${band.barrier}, E_g ${fmt(Eg)} eV, offsets dE_c ${fmt(dEc * 1000)} meV and dE_v ${fmt(dEv * 1000)} meV, transition ${fmt(EX)} eV` });
  const line = (x1, y1, x2, y2, cls) => svg.append(svgEl("path", { d: `M ${x1} ${y1} L ${x2} ${y2}`, class: cls }));
  const text = (x, y, s, cls = "bp-t", anchor = "start") => { const t = svgEl("text", { x, y, class: cls, "text-anchor": anchor }); t.textContent = s; svg.append(t); return t; };

  // axes: energy ticks on the left (relative to the dot valence-band edge), growth axis in nm along the bottom
  line(PADL, PADT, PADL, H - PADB, "bp-axis");
  line(PADL, H - PADB, W - PADR, H - PADB, "bp-axis");
  const span = yMax - yMin;
  const rawStep = span / 5, mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= rawStep) || 10 * mag;
  for (let e = Math.ceil(yMin / step) * step; e <= yMax + 1e-9; e += step) {
    line(PADL - 4, Y(e), PADL, Y(e), "bp-axis");
    text(PADL - 8, Y(e) + 3.5, fmt(Number(e.toFixed(6))), "bp-tick", "end");
  }
  text(14, PADT - 10, "E (eV), dot valence-band edge = 0", "bp-t");
  for (const z of [-Z_HALF_NM, -Z_HALF_NM / 2, 0, Z_HALF_NM / 2, Z_HALF_NM]) {
    line(X(z), H - PADB, X(z), H - PADB + 4, "bp-axis");
    text(X(z), H - PADB + 16, fmt(z), "bp-tick", "middle");
  }
  text((PADL + W - PADR) / 2, H - 8, "growth axis z (nm)", "bp-t", "middle");

  // dot slab, band edges (barrier | well | barrier), confined levels
  svg.append(svgEl("rect", { x: X(z0), y: PADT, width: X(z1) - X(z0), height: H - PADT - PADB, class: "bp-slab" }));
  const edge = (bar, at, cls) => {
    const brk = flags.hard_wall ? " bp-break" : "";
    svg.append(svgEl("path", { d: `M ${X(-Z_HALF_NM)} ${Y(bar)} L ${X(z0)} ${Y(bar)} L ${X(z0)} ${Y(at(0))} L ${X(z1)} ${Y(at(1))} L ${X(z1)} ${Y(bar)} L ${X(Z_HALF_NM)} ${Y(bar)}`, class: cls + brk }));
  };
  edge(cbBar, cbAt, "bp-edge bp-cb");
  edge(vbBar, vbAt, "bp-edge bp-vb");
  const wellX0 = X(z0), wellX1 = X(z1), wellW = wellX1 - wellX0;
  // levels: e above the well CB edge, h below the well VB edge; polar dots put them at opposite interfaces
  const eU = flags.polar ? [0.55, 1] : [0, 1], hU = flags.polar ? [0, 0.45] : [0, 1];
  const eY = (flags.polar ? cbAt(0.78) : Eg) + Ee, hY = (flags.polar ? vbAt(0.22) : 0) - Eh;
  line(wellX0 + wellW * eU[0] + 3, Y(eY), wellX0 + wellW * eU[1] - 3, Y(eY), "bp-level bp-e");
  line(wellX0 + wellW * hU[0] + 3, Y(hY), wellX0 + wellW * hU[1] - 3, Y(hY), "bp-level bp-h");
  text(X(Z_HALF_NM) - 2, Y(cbBar) - 6, "conduction band edge", "bp-t", "end");
  text(X(Z_HALF_NM) - 2, Y(vbBar) + 14, "valence band edge", "bp-t", "end");
  text(wellX1 + 6, Y(eY) + 4, "e ground level", "bp-t");
  text(wellX1 + 6, Y(hY) + 4, "h ground level", "bp-t");
  // emission: a downward arrow from the electron level to the hole level, labelled with hc/lambda
  const xe = wellX0 + wellW * (flags.polar ? 0.78 : 0.5), xh = wellX0 + wellW * (flags.polar ? 0.22 : 0.5);
  line(xe, Y(eY) + 3, xh, Y(hY) - 8, "bp-photon");
  svg.append(svgEl("path", { d: `M ${xh - 4.5} ${Y(hY) - 13} L ${xh} ${Y(hY) - 5} L ${xh + 4.5} ${Y(hY) - 13} Z`, class: "bp-photon-head" }));
  text(Math.min(xe, xh) - 10, (Y(eY) + Y(hY)) / 2 + 4, `hν = ${fmt(EX)} eV  (${fmt(lam)} nm)`, "bp-t bp-hv", "end");
  // offsets, each with the tag of its literature number
  const tg = (leaf) => `[${normTag(leaf.tag) || "?"}]`;
  const dim = (xx, y1, y2) => { line(xx, y1, xx, y2, "bp-dim"); line(xx - 4, y1, xx + 4, y1, "bp-dim"); line(xx - 4, y2, xx + 4, y2, "bp-dim"); };
  const xo = X(-Z_HALF_NM) + 12;
  dim(xo, Y(Eg), Y(cbBar));
  const smallC = Math.abs(Y(Eg) - Y(cbBar)) < 26, smallV = Math.abs(Y(0) - Y(vbBar)) < 26;   // a short offset: label beside it, off the edge line
  text(xo + 8, smallC ? Y(cbBar) + 18 : (Y(Eg) + Y(cbBar)) / 2 + 4, flags.hard_wall ? `ΔEc ≈ ${fmt(dEc * 1000)} meV ${tg(sel.dEc)}, effectively infinite` : `ΔEc = ${fmt(dEc * 1000)} meV ${tg(sel.dEc)}`, "bp-t bp-off");
  dim(xo, Y(0), Y(vbBar));
  text(xo + 8, smallV ? Y(vbBar) + 18 : (Y(0) + Y(vbBar)) / 2 + 4, flags.hard_wall ? `ΔEv ≈ ${fmt(dEv * 1000)} meV ${tg(sel.dEv)}, effectively infinite` : `ΔEv = ${fmt(dEv * 1000)} meV ${tg(sel.dEv)}`, "bp-t bp-off");
  text(X(-Z_HALF_NM) + 8, Y(cbBar) - 6, `Eg,barrier = ${fmt(Egb)} eV ${tg(band.Eg_barrier)}`, "bp-t bp-off");
  text(xo + 8, Y(Eg / 2) + 4, `Eg,dot = ${fmt(Eg)} eV ${tg(sel.Eg)}`, "bp-t bp-off");
  if (flags.hard_wall) text(X(Z_HALF_NM) - 2, Y(cbBar) + 14, "barriers broken: effectively infinite (hard wall)", "bp-note", "end");
  // well width, to scale
  const yd = Y(vbBar) + 11;
  line(wellX0, yd, wellX1, yd, "bp-dim");
  line(wellX0, yd - 4, wellX0, yd + 4, "bp-dim");
  line(wellX1, yd - 4, wellX1, yd + 4, "bp-dim");
  text((wellX0 + wellX1) / 2, yd + 14, `l_z = ${fmt(lz)} nm`, "bp-t", "middle");
  // flags are on the picture, not hidden in a tooltip
  const notes = [];
  if (flags.polar) notes.push("polar: tilted dot (schematic tilt [A]), e and h at opposite interfaces");
  if (flags.type2_risk) notes.push("weak hole confinement (type-II risk)");
  notes.forEach((n, i) => text(PADL + 6, PADT + 10 + 13 * i, n, "bp-note bp-flagnote"));

  const tagOf = (n) => normTag(L(n)?.tag);
  const varName = { indot: "in-dot", strained: "strained", unstr: "unstrained", half: "half-way unstrained / strained" };
  const row = (k, v, tag) => h("div.bp-kv", h("span.bp-k", k), h("span.bp-v", v), tag ? tagChip(tag, { size: "sm" }) : null);
  const T = (src) => tex(src);
  const facts = h("div.bp-facts",
    row(T("h\\nu = E_{\\mathrm X}"), `${fmt(EX)} eV`, tagOf("lambda_nm")),
    row(T("E_{g,\\mathrm{dot}}"), `${fmt(Eg)} eV (${varName[sel.Eg.from]})`, sel.Eg.tag),
    row(T("E_{g,\\mathrm{barrier}}"), `${fmt(Egb)} eV`, band.Eg_barrier.tag),
    row(T("\\Delta E_{\\mathrm c}"), `${fmt(dEc * 1000)} meV (${varName[sel.dEc.from]})`, sel.dEc.tag),
    row(T("\\Delta E_{\\mathrm v}"), `${fmt(dEv * 1000)} meV (${varName[sel.dEv.from]})`, sel.dEv.tag),
    row(T("E_{\\mathrm e},\\ E_{\\mathrm h}"), `${fmt(Ee * 1000)}, ${fmt(Eh * 1000)} meV`, "A"),
    row(T("\\ell_z"), `${fmt(lz)} nm`, tagOf("l_z_nm")),
    Number.isFinite(lxy) ? row(T("\\ell_{xy}"), `${fmt(lxy)} nm`, tagOf("l_xy_nm")) : null,
    Number.isFinite(eb) ? row(T("E_{\\mathrm b}"), `${fmt(eb)} meV`, tagOf("E_b_meV")) : null,
    row("Charge state", String(charge), tagOf("charge_state")),
    row("Dipole", String(dipole).replace(/_/g, " "), tagOf("dipole")));
  root.append(svg, facts);
  const flagP = (s) => h("p.bp-flag", s);
  if (flags.polar) root.append(flagP("Polar nitride: spontaneous and piezoelectric fields tilt the bands, so the dot region is drawn tilted with the electron and hole at opposite interfaces. For tall dots E_X can fall below the bulk gap (quantum-confined Stark effect), and the flat-band level split below is then clamped at zero."));
  if (flags.type2_risk) root.append(flagP("Weak hole confinement (type-II risk): the light-hole offset is about 0 and the result hinges on the unmeasured GaP/InP valence offset and on CuPt ordering of the barrier. The heavy-hole offset is drawn."));
  if (flags.hard_wall) root.append(flagP("Hard wall: the shell or glass puts both carriers several eV from the nanocrystal edges. The barriers are drawn broken and are effectively infinite."));
  if (mixed) root.append(flagP(`Mixed strain variants: the drawn conduction offset is E_g,barrier - E_g,dot - ΔE_v = ${fmt(dEc * 1000)} meV; the table gives ${fmt(dEcTab * 1000)} meV (${varName[sel.dEc.from]}).`));
  const allLeaves = [band.Eg_barrier, ...Object.values(band.masses || {}), ...Object.values(band.variants || {}).flatMap((vv) => Object.values(vv))];
  const vLeaves = allLeaves.filter((l) => normTag(l?.tag) === "V");
  const nChecked = vLeaves.filter((l) => l.checked === true).length;
  const nUnchecked = vLeaves.length - nChecked;
  root.append(h("p.bp-cap", richText(`${entry.label}: ${band.dot} in ${band.barrier}, type ${band.type}. Offsets are literature values at 0-10 K (display only; no solver reads them). `
    + `Sources: ${(band.sources || []).join("; ")}. `
    + (nChecked ? `${nChecked} [V] value${nChecked === 1 ? " was" : "s were"} checked against the source tables (2026-10-13). ` : "")
    + (nUnchecked ? `${nUnchecked} [V] value${nUnchecked === 1 ? " is" : "s are"} recalled from the reviews and not yet re-checked against the papers. ` : "")
    + "Level split: inverse-mass rule, clamped to 0.85 of each offset [A].")));
  const tags = [band.Eg_barrier.tag, sel.Eg.tag, sel.dEc.tag, sel.dEv.tag].map((t) => normTag(t));
  root.dataset.ex = EX.toFixed(4);
  root.dataset.lz = String(lz);
  root.dataset.eg = Eg.toFixed(4);
  root.dataset.dec = (dEc * 1000).toFixed(1);
  root.dataset.dev = (dEv * 1000).toFixed(1);
  root.dataset.flags = Object.keys(flags).filter((k) => flags[k]).join(",");
  const sig = [entry.kind, EX.toFixed(4), lz, lxy, eb, charge, dipole, Eg.toFixed(4), (dEc * 1000).toFixed(1), (dEv * 1000).toFixed(1)].join("|");
  root.dataset.sig = sig;
  root.dataset.tag = widest([...tags, ...["lambda_nm", "l_z_nm", "l_xy_nm"].map(tagOf)]) || "A";
  return { el: root, sig };
}
