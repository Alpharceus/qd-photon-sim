// ONE place for typeset mathematics (devlib-r2-math): KaTeX, vendored under vendor/katex (offline).
//   tex(src, {display})   -> element rendered by KaTeX; plain-text fallback when KaTeX did not load
//   richText(str)         -> DocumentFragment: text with inline $...$ and display $$...$$ (\$ = a literal dollar)
//   symbolFor(meta)       -> TeX of a parameter / scalar: meta.tex, else derived from the registry id
//   symbolNode(meta)      -> tex(symbolFor(meta))
//   plainSymbol(texSrc)   -> Unicode approximation (aria labels, tooltips, search text, Plotly axis titles)
// Nothing here ever uses innerHTML on data: text goes in as text nodes, TeX goes through KaTeX with
// trust:false (no \href, \url, \includegraphics, \htmlClass ...), strict ignored.

const KATEX_OPTS = { throwOnError: false, output: "htmlAndMathml", trust: false, strict: "ignore", maxExpand: 200 };

export function katexReady() {
  return typeof window !== "undefined" && typeof window.katex?.render === "function";
}

// ---------------------------------------------------------------- plain-text approximation
const GREEK_TEX = {
  alpha: "α", beta: "β", gamma: "γ", delta: "δ", epsilon: "ε", varepsilon: "ε", zeta: "ζ", eta: "η", theta: "θ",
  vartheta: "ϑ", kappa: "κ", lambda: "λ", mu: "μ", nu: "ν", xi: "ξ", pi: "π", rho: "ρ", sigma: "σ", tau: "τ",
  phi: "φ", varphi: "φ", chi: "χ", psi: "ψ", omega: "ω",
  Gamma: "Γ", Delta: "Δ", Theta: "Θ", Lambda: "Λ", Xi: "Ξ", Pi: "Π", Sigma: "Σ", Phi: "Φ", Psi: "Ψ", Omega: "Ω",
  hbar: "ħ", times: "×", cdot: "·", pm: "±", approx: "≈", sim: "~", le: "≤", leq: "≤", ge: "≥", geq: "≥",
  ne: "≠", infty: "∞", to: "→", rightarrow: "→", ll: "≪", gg: "≫", propto: "∝", partial: "∂", ldots: "…",
  langle: "⟨", rangle: "⟩", perp: "⊥", parallel: "∥", ell: "ℓ", sqrt: "√", int: "∫", sum: "Σ", dot: "", hat: "",
};
const SUPS = { 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹", "-": "⁻", "+": "⁺", "*": "*", "(": "⁽", ")": "⁾", n: "ⁿ", i: "ⁱ" };
const SUBS = { 0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄", 5: "₅", 6: "₆", 7: "₇", 8: "₈", 9: "₉", "-": "₋", "+": "₊", "(": "₍", ")": "₎" };
const mapAll = (s, tbl) => ([...s].every((c) => tbl[c] !== undefined) ? [...s].map((c) => tbl[c]).join("") : null);

/** Best-effort Unicode rendering of a TeX fragment: Greek letters, simple sub/superscripts, \mathrm{...}. */
export function plainSymbol(src) {
  let s = String(src == null ? "" : src);
  s = s.replace(/\\(?:mathrm|mathbf|mathit|text|operatorname|boldsymbol|hat|bar|tilde|dot)\s*\{([^{}]*)\}/g, "$1");
  s = s.replace(/\\(?:mathrm|mathbf)\s*([A-Za-z])/g, "$1");
  s = s.replace(/\\sqrt\s*\{([^{}]*)\}/g, "√($1)");
  s = s.replace(/\\frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}/g, "($1)/($2)");
  s = s.replace(/\\[,;:! ]|\\quad|\\qquad/g, " ");
  s = s.replace(/\\([A-Za-z]+)/g, (m, n) => GREEK_TEX[n] ?? n);
  s = s.replace(/\\([{}$%&#_])/g, "$1");
  s = s.replace(/\^\s*\{([^{}]*)\}|\^\s*([A-Za-z0-9*+-])/g, (m, a, b) => { const t = a ?? b; return mapAll(t, SUPS) ?? `^(${t})`; });
  s = s.replace(/_\s*\{([^{}]*)\}|_\s*([A-Za-z0-9])/g, (m, a, b) => { const t = a ?? b; return mapAll(t, SUBS) ?? `(${t})`; });
  s = s.replace(/[{}]/g, "").replace(/\s+/g, " ").trim();
  return s;
}

/** One element for a multi-part result: in a flex / grid parent the pieces must stay ONE item, as the plain text node was. */
function wrapRt(frag) {
  if (frag.childNodes.length < 2) return frag;
  const span = document.createElement("span");
  span.className = "rt";
  span.appendChild(frag);
  return span;
}

// KaTeX's MathML carries the TeX source in <annotation>; it is not needed for a11y and would put backslashes into textContent.
function stripAnnotation(el) {
  for (const a of el.querySelectorAll("annotation")) a.remove();
}

// ---------------------------------------------------------------- tex()
/** Element rendered by KaTeX. display: block (centred) math, else inline. Never throws. */
const TEX_CACHE = new Map(); // "d|src" -> rendered element; repeated symbols (every animation frame) are a clone, not a KaTeX run
export function tex(src, { display = false, cls = "" } = {}) {
  const s = String(src == null ? "" : src);
  const key = `${display ? "d" : "i"}|${s}`;
  let tpl = TEX_CACHE.get(key);
  if (!tpl) {
    tpl = document.createElement(display ? "div" : "span");
    tpl.className = `math${display ? " math-display" : ""}`;
    tpl.dataset.tex = s;
    let ok = false;
    if (katexReady()) {
      try {
        window.katex.render(s, tpl, { ...KATEX_OPTS, displayMode: !!display });
        stripAnnotation(tpl);
        ok = true;
      } catch { /* fall through to the plain text */ }
    }
    if (!ok) {
      tpl.classList.add("math-plain");
      tpl.textContent = plainSymbol(s);
      return withCls(tpl.cloneNode(true), cls); // not cached: KaTeX may still arrive
    }
    if (TEX_CACHE.size > 3000) TEX_CACHE.clear();
    TEX_CACHE.set(key, tpl);
  }
  return withCls(tpl.cloneNode(true), cls);
}
const withCls = (el, cls) => { if (cls) for (const c of String(cls).split(/\s+/).filter(Boolean)) el.classList.add(c); return el; };

/** Re-render fallbacks once KaTeX has arrived late (the vendored script is deferred and normally runs first). */
export function upgradeMath(root = document) {
  if (!katexReady()) return 0;
  let n = 0;
  for (const el of root.querySelectorAll(".math.math-plain")) {
    try {
      window.katex.render(el.dataset.tex || "", el, { ...KATEX_OPTS, displayMode: el.classList.contains("math-display") });
      stripAnnotation(el);
      el.classList.remove("math-plain");
      n++;
    } catch { /* leave the plain text */ }
  }
  return n;
}

// ---------------------------------------------------------------- richText()
/** Split `str` into [{t:"text"|"math"|"display", v}] ; `\$` is a literal dollar outside math. */
export function parseRich(str) {
  const s = String(str == null ? "" : str);
  const out = [];
  let buf = "";
  let i = 0;
  const flush = () => { if (buf) { out.push({ t: "text", v: buf }); buf = ""; } };
  const findClose = (from, delim) => {
    for (let j = from; j < s.length; j++) {
      if (s[j] === "\\") { j++; continue; }
      if (delim === "$$" ? s.startsWith("$$", j) : s[j] === "$") return j;
    }
    return -1;
  };
  while (i < s.length) {
    const c = s[i];
    if (c === "\\" && s[i + 1] === "$") { buf += "$"; i += 2; continue; }
    if (c === "$") {
      const disp = s[i + 1] === "$";
      const delim = disp ? "$$" : "$";
      const start = i + delim.length;
      const end = findClose(start, delim);
      if (end > start) {
        flush();
        out.push({ t: disp ? "display" : "math", v: s.slice(start, end) });
        i = end + delim.length;
        continue;
      }
    }
    buf += c;
    i++;
  }
  flush();
  return out;
}

/** DocumentFragment: plain segments as text nodes, `$...$` / `$$...$$` through KaTeX. XSS-safe by construction. */
export function richText(str) {
  const frag = document.createDocumentFragment();
  for (const seg of parseRich(str)) {
    if (seg.t === "text") frag.appendChild(autoMath(seg.v) || document.createTextNode(seg.v));
    else frag.appendChild(tex(seg.v, { display: seg.t === "display" }));
  }
  return wrapRt(frag);
}

// Plain-ASCII notation that the app and the backend notes still write (T_hs, g²(0), γ*, Γ_X) is typeset on the way
// into the DOM, so no call site has to be converted by hand and no prose shows an underscore as a subscript.
// A subscript is `_` + 1..6 alphanumerics, NOT followed by another `_word` (so snake_case ids such as a_ac_ns_per_K stay untouched
// unless the registry knows them: then the id is replaced by its symbol).
const AUTO_RE = new RegExp(
  [
    "(?<![A-Za-z0-9_])g(?:²|\\^2|\\(2\\))(?:_(?<gsub>[A-Za-z]{2,6})(?![A-Za-z0-9_]))?(?:\\((?<garg>[0τt])\\))?", // g²(0), g^2(0), g(2)(0), bare g², g²_dot
    "(?<![A-Za-z0-9_.])g2(?:\\((?<g2arg>[0τt])\\))?(?![A-Za-z0-9_])", // g2(0), bare g2 (not g2_vs_T, g2.dat)
    "(?<![A-Za-z0-9_.])(?<base>[A-Za-zΑ-ωℓ])_(?<sub>[A-Za-z0-9]{1,10})(?![A-Za-z0-9_])", // T_hs, T_c, Γ_X, k_nr, V_j, F_cav, t_Lorentz
    "(?<![A-Za-z0-9_.])(?<gbase>alpha|beta|gamma|delta|epsilon|eta|theta|kappa|lambda|mu|nu|rho|sigma|tau|phi|omega|Gamma|Delta|Omega)_(?<gsub2>[A-Za-z0-9]{1,10})(?![A-Za-z0-9_])", // Gamma_zpl
    "(?<![A-Za-z0-9_])(?:γ|gamma)\\*", // γ*
    "(?<mant>\\d(?:\\.\\d+)?)×10(?<exp>[⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)", // 1.2×10⁻⁵ (what fmt() writes)
    "(?<![A-Za-z0-9_.])(?<id>[A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+)(?![A-Za-z0-9_])", // a registry id (eta_out, pol_factor) in prose
  ].join("|"), "g");
const SUPCH = { "⁻": "-", "⁺": "+", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9" };
const subTex = (w) => (/^\d+$/.test(w) ? w : `\\mathrm{${w}}`);

// registry ids (parameters and observables) -> TeX, so a backend note such as "eta_out includes the NA collection" shows the symbol
const ID_TEX = new Map();
/** Learn the symbols of a /api/compose/registry payload: the id's own tex when some entry has one, else the derived symbol. */
export function registerSymbols(reg) {
  if (!reg) return 0;
  const seen = new Map();
  const note = (name, t) => { if (name && (t || !seen.has(name))) seen.set(name, t || null); };
  for (const plural of ["materials", "structures", "drives", "setups"]) {
    for (const e of reg[plural] || []) for (const p of e.params || []) note(p.name, p.tex);
  }
  for (const o of reg.observables || []) note(o.id, o.tex);
  for (const [name, t] of seen) ID_TEX.set(name, t || symbolFor({ name }));
  return seen.size;
}

/** DocumentFragment for a plain string with typeset ASCII notation, or null when nothing in it needs typesetting. */
export function autoMath(str) {
  const s = String(str);
  if (!/[_²^*(×]/.test(s) || s.length > 4000) return null;
  let last = 0;
  let hit = false;
  const frag = document.createDocumentFragment();
  for (const m of s.matchAll(AUTO_RE)) {
    const g = m.groups;
    let t;
    if (g.mant !== undefined) {
      // numbers change every animation frame: plain DOM (m x10 + <sup>), never a KaTeX run
      if (m.index > last) frag.appendChild(document.createTextNode(s.slice(last, m.index)));
      const sup = document.createElement("sup");
      sup.className = "num-exp";
      sup.textContent = [...g.exp].map((c) => SUPCH[c]).join("");
      frag.append(document.createTextNode(`${g.mant} × 10`), sup);
      last = m.index + m[0].length;
      hit = true;
      continue;
    } else if (g.base !== undefined) t = `${g.base}_{${subTex(g.sub)}}`;
    else if (g.gbase !== undefined) t = `\\${g.gbase}_{${subTex(g.gsub2)}}`;
    else if (g.id !== undefined) {
      t = ID_TEX.get(g.id);
      if (!t) continue; // not a registry id: leave the text alone
    } else if (m[0].startsWith("γ") || m[0].startsWith("gamma")) t = "\\gamma^{*}";
    else {
      const arg = g.garg ?? g.g2arg;
      t = `g^{(2)}${g.gsub ? `_{${subTex(g.gsub)}}` : ""}${arg ? `(${arg === "τ" ? "\\tau" : arg})` : ""}`;
    }
    if (m.index > last) frag.appendChild(document.createTextNode(s.slice(last, m.index)));
    const node = tex(t);
    node.dataset.src = m[0]; // the ASCII source, so a verbatim-text check can read the string the backend sent
    frag.appendChild(node);
    last = m.index + m[0].length;
    hit = true;
  }
  if (!hit) return null;
  if (last < s.length) frag.appendChild(document.createTextNode(s.slice(last)));
  return wrapRt(frag);
}

/** Plain-text version of a rich string (search text, aria-labels, tooltips, Plotly titles). */
export function richPlain(str) {
  return parseRich(str).map((seg) => (seg.t === "text" ? seg.v : plainSymbol(seg.v))).join("");
}

// ---------------------------------------------------------------- symbolFor()
const GREEK_NAMES = new Set(["alpha", "beta", "gamma", "delta", "zeta", "eta", "theta", "kappa", "lambda", "mu", "nu", "xi",
  "pi", "rho", "sigma", "tau", "phi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi", "Omega"]);
const GREEK_ALIAS = { eps: "varepsilon", epsilon: "varepsilon", vphi: "varphi" };
// trailing unit tokens that the registry bakes into ids (fwhm_ps, lambda_nm, a_ac_ns_per_K, rho_kg_m3 ...)
const UNIT_TOKENS = new Set(["nm", "um", "ps", "ns", "ueV", "meV", "eV", "MHz", "Hz", "GHz", "K", "nW", "uA", "ohm", "F", "pi", "ps2",
  "s", "cps", "kg", "m", "per", "4K", "W", "V", "mA", "deg", "um2", "cm2", "m3", "pF", "pm"]);
// ids whose derivation is better written by hand (also the common scalars of the results tables)
const OVERRIDES = {
  g2: "g^{(2)}", g2_0: "g^{(2)}(0)", g2_op: "g^{(2)}(0)", g2_raw: "g^{(2)}_{\\mathrm{raw}}(0)", g2_bgc: "g^{(2)}_{\\mathrm{bgc}}(0)",
  g2_irf: "g^{(2)}_{\\mathrm{IRF}}(0)", g2_binned: "g^{(2)}_{\\mathrm{bin}}(0)", g2_XX: "g^{(2)}_{XX}(0)", g2_tau: "g^{(2)}(\\tau)",
  hom_I: "I_{\\mathrm{HOM}}", hom_V: "V_{\\mathrm{HOM}}", hom_V_raw: "V_{\\mathrm{HOM}}^{\\mathrm{raw}}", hom_V_corr: "V_{\\mathrm{HOM}}^{\\mathrm{corr}}",
  hom_I_gated: "I_{\\mathrm{HOM}}^{\\mathrm{gate}}", hom_XX: "I_{\\mathrm{HOM}}^{XX}", hom_zpl_fraction: "q_{\\mathrm{ZPL}}",
  F_P: "F_{\\mathrm P}", F_P_eff: "F_{\\mathrm P}L", F_P_peak: "F_{\\mathrm P}", F_p_rail: "F_{\\mathrm{rail}}", NA: "\\mathrm{NA}", T_NA: "T_{\\mathrm{NA}}",
  Z_phonon: "Z(T)", Z_filt: "Z_{\\mathrm{filt}}", gamma_star_ns: "\\gamma^{*}", Gamma_ns: "\\Gamma", beta_entered: "\\beta",
  eta_out: "\\eta_{\\mathrm{out}}", eta_extract: "\\beta\\,\\eta_{\\mathrm{out}}", eta_rad: "\\eta_{\\mathrm{rad}}", eta_rad_eff: "\\eta_{\\mathrm{rad}}^{\\mathrm{eff}}",
  eta_prep: "\\eta_{\\mathrm{prep}}", eta_pol: "\\eta_{\\mathrm{pol}}", eta_blink: "\\eta_{\\mathrm{blink}}",
  eta_prep_x_eta_rad: "\\eta_{\\mathrm{prep}}\\,\\eta_{\\mathrm{rad}}", n_ph: "n_{\\mathrm{ph}}", n_exc: "n_{\\mathrm{exc}}", n_XX: "n_{XX}",
  P_X_end: "P_X(t_{\\mathrm{end}})", B_first_lens: "B_{\\mathrm{lens}}", B_lens: "B_{\\mathrm{lens}}", B_fibre: "B_{\\mathrm{fibre}}",
  B_first_lens_Hz: "\\dot B_{\\mathrm{lens}}", B_fibre_Hz: "\\dot B_{\\mathrm{fibre}}", rate_detected_Hz: "R_{\\mathrm{det}}",
  tau_ps: "\\tau", tau_1e_ps: "\\tau_{1/e}", tau_XX_ps: "\\tau_{XX}", tau_EL_ps: "\\tau_{\\mathrm{EL}}", tau_purcell_ps: "\\tau_{\\mathrm P}",
  linewidth_ueV: "\\hbar\\Gamma_{\\mathrm{tot}}", detuning_ueV: "\\delta_{\\mathrm{EC}}", spectral_overlap: "\\varepsilon_{XX}",
  inversion: "w", eps_op: "\\varepsilon", rho_op: "\\rho", gamma_op: "\\Gamma", T_c: "T_{\\mathrm c}", T_j_op: "T_{\\mathrm j}", dT_J: "\\Delta T_{\\mathrm J}",
};

const escWord = (w) => String(w).replace(/[\\{}^_$&#%~]/g, " ");
const single = (t) => /^[A-Za-z0-9]$/.test(t);
const WORD_FIX = { xx: "XX", x: "X" };
const subOf = (toks0) => {
  const toks = toks0.map((x) => WORD_FIX[x] ?? x);
  if (!toks.length) return "";
  if (toks.length === 1 && single(toks[0])) return `_{${toks[0]}}`;
  return `_{${/^\d+$/.test(toks.join("")) ? toks.join("") : `\\mathrm{${escWord(toks.join(" ")).trim().replace(/ +/g, "\\,")}}`}}`;
};

/** TeX of a parameter / scalar / observable id: meta.tex when present, else derived. Never raw snake_case. */
export function symbolFor(meta) {
  const m = typeof meta === "string" ? { name: meta } : (meta || {});
  if (m.tex) return String(m.tex);
  const name = String(m.name ?? m.id ?? m.key ?? "");
  if (!name) return "";
  if (OVERRIDES[name]) return OVERRIDES[name];
  const toks = name.split("_").filter(Boolean);
  while (toks.length > 1 && UNIT_TOKENS.has(toks[toks.length - 1])) toks.pop();
  const key = toks.join("_");
  if (OVERRIDES[key]) return OVERRIDES[key];
  let base = toks.shift();
  let num = "";
  const dm = base.match(/^([A-Za-z]+?)(\d+)$/);
  if (dm) { base = dm[1]; num = dm[2]; }
  const star = toks[0] === "star";
  if (star) toks.shift();
  const greek = GREEK_ALIAS[base] || (GREEK_NAMES.has(base) ? base : null);
  const sub = [...(num ? [num] : []), ...toks];
  if (greek) return `\\${greek}${star ? "^{*}" : ""}${subOf(sub)}`;
  if (base.length <= 2) return `${base}${subOf(sub)}`;
  // a plain word: roman text, not an italic product of letters
  return `\\mathrm{${escWord([base + num, ...toks].join(" ")).trim().replace(/ +/g, "\\,")}}`;
}

export function symbolNode(meta, opts) { return tex(symbolFor(meta), opts); }

const AXIS_PLURAL = { material: "materials", structure: "structures", drive: "drives", setup: "setups" };
const AXIS_WORD = { material: "Material", structure: "Structure", drive: "Drive", setup: "Setup", observable: "Observable" };

/** {"drive.fwhm_ps": "\\tau_{\\mathrm{FWHM}}", "observable.g2_0": ...} from a /api/compose/registry payload (first tex per axis + name). */
export function texIndex(reg) {
  const idx = {};
  for (const [axis, plural] of Object.entries(AXIS_PLURAL)) {
    for (const e of reg?.[plural] || []) for (const p of e.params || []) if (p.tex && !idx[`${axis}.${p.name}`]) idx[`${axis}.${p.name}`] = p.tex;
  }
  for (const o of reg?.observables || []) if (o.tex) idx[`observable.${o.id}`] = o.tex;
  return idx;
}

/** TeX of a spec path such as "drive.fwhm_ps" (registry tex when an index is given, else derived). */
export function pathSymbol(path, idx = null) {
  const s = String(path);
  if (s === "T_K") return "T";
  const i = s.indexOf(".");
  const name = i < 0 ? s : s.slice(i + 1);
  return (idx && idx[s]) || symbolFor({ name });
}

/** Rich string "Drive $\\tau$" for a spec path: axis word plus the typeset symbol; the registry id is never shown. */
export function pathRichText(path, idx = null) {
  const s = String(path);
  if (s === "T_K") return "Temperature $T$";
  const axis = s.split(".")[0];
  return `${AXIS_WORD[axis] || axis} $${pathSymbol(s, idx)}$`;
}

/** Unicode symbol text of a parameter (search, aria, Plotly titles). */
export function symbolText(meta) { return plainSymbol(symbolFor(meta)); }

// ---------------------------------------------------------------- Plotly text
// Plotly renders `$...$` only through MathJax (not loaded: the app is offline). Axis titles, annotations,
// legends and hover text instead get Plotly's own HTML subset (<sub>, <sup>, <i>) from the same TeX, and
// plain-ASCII notation (T_hs, g2(0), E_X) is typeset the same way. installPlotlyMath() wraps the Plotly
// entry points once, so every chart gets consistent notation without each caller converting its strings.
const GREEK_UP = new Set(["Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi", "Omega"]);
const escHtml = (c) => (c === "<" ? "&lt;" : c === ">" ? "&gt;" : c === "&" ? "&amp;" : c);

/** TeX -> Plotly-HTML (<i>, <sub>, <sup>, Unicode). Covers the subset the registry symbols use. */
export function texToPlotly(src) {
  return texToPlotlyInner(String(src == null ? "" : src), false);
}

function texToPlotlyInner(s, upright) {
  let out = "";
  let i = 0;
  const readGroup = () => {
    let depth = 1;
    const start = i;
    while (i < s.length && depth) { if (s[i] === "{") depth++; else if (s[i] === "}") depth--; i++; }
    return s.slice(start, depth ? i : i - 1);
  };
  const readArg = () => {
    while (s[i] === " ") i++;
    if (s[i] === "{") { i++; return readGroup(); }
    if (s[i] === "\\") { const m = /^\\[A-Za-z]+/.exec(s.slice(i)); const t = m ? m[0] : s.slice(i, i + 2); i += t.length; return t; }
    return s[i++] ?? "";
  };
  while (i < s.length) {
    const c = s[i];
    if (c === "\\") {
      const m = /^\\([A-Za-z]+|.)/.exec(s.slice(i));
      const name = m[1];
      i += m[0].length;
      if (name === "mathrm" || name === "text" || name === "operatorname" || name === "mathbf") out += texToPlotlyInner(readArg(), true);
      else if (name === "sqrt") out += `√(${texToPlotlyInner(readArg(), upright)})`;
      else if (name === "frac") { const a = readArg(); const b = readArg(); out += `(${texToPlotlyInner(a, upright)})/(${texToPlotlyInner(b, upright)})`; }
      else if (name === "hat" || name === "boldsymbol" || name === "dot" || name === "bar" || name === "left" || name === "right") { if (["hat", "boldsymbol", "dot", "bar"].includes(name)) out += texToPlotlyInner(readArg(), upright); }
      else if (name === "," || name === ";" || name === ":" || name === " " || name === "quad" || name === "!") out += name === "!" ? "" : " ";
      else if (GREEK_TEX[name] !== undefined) { const g = GREEK_TEX[name]; out += upright || GREEK_UP.has(name) || !/[α-ω]/.test(g) ? g : `<i>${g}</i>`; }
      else out += name;
      continue;
    }
    if (c === "{") { i++; out += texToPlotlyInner(readGroup(), upright); continue; }
    if (c === "^" || c === "_") {
      i++;
      const inner = texToPlotlyInner(readArg(), upright);
      out += c === "^" ? `<sup>${inner}</sup>` : `<sub>${inner}</sub>`;
      continue;
    }
    i++;
    if (/[A-Za-z]/.test(c) && !upright) out += `<i>${c}</i>`;
    else out += escHtml(c);
  }
  return out.replace(/<\/i><i>/g, "");
}

/** Plain-ASCII notation in chart strings -> Plotly HTML (T_hs, g2(0), E_X, x^2). Idempotent. */
function plainNotation(s) {
  return s
    .replace(/(^|[^A-Za-z0-9_<])g(?:²|2)\((0|τ|t)\)/g, "$1<i>g</i><sup>(2)</sup>($2)")
    .replace(/(^|[^A-Za-z0-9_<\/])([A-Za-zΓΔΘΛΞΠΣΦΨΩα-ω])_\{?([A-Za-z0-9,]{1,8})\}?(?![A-Za-z0-9_])/g, (m, pre, b, sub) => `${pre}${/[a-zα-ω]/.test(b) ? `<i>${b}</i>` : b}<sub>${sub}</sub>`)
    .replace(/([A-Za-z\)])\^(\d)(?![0-9])/g, "$1<sup>$2</sup>");
}

/** A chart string: `$..$` through texToPlotly, the rest through plainNotation. */
export function plotlyText(str) {
  if (typeof str !== "string" || !str) return str;
  if (!str.includes("$")) return /_|g2|g²|\^/.test(str) ? plainNotation(str) : str;
  return parseRich(str).map((seg) => (seg.t === "text" ? plainNotation(seg.v) : texToPlotly(seg.v))).join("");
}

const TEXT_KEYS = new Set(["text", "name", "title", "hovertemplate", "ticktext", "label", "legendgrouptitle"]);
function walkPlotly(o, depth = 0) {
  if (o == null || depth > 8) return;
  if (Array.isArray(o)) { for (const v of o) if (v && typeof v === "object" && !ArrayBuffer.isView(v)) walkPlotly(v, depth + 1); return; }
  if (typeof o !== "object" || ArrayBuffer.isView(o)) return;
  for (const k of Object.keys(o)) {
    const v = o[k];
    if (TEXT_KEYS.has(k)) {
      if (typeof v === "string") o[k] = plotlyText(v);
      else if (Array.isArray(v) && v.every((x) => typeof x === "string")) o[k] = v.map(plotlyText);
      else if (v && typeof v === "object") walkPlotly(v, depth + 1);
    } else if (v && typeof v === "object" && !["x", "y", "z", "customdata", "marker", "line", "error_x", "error_y", "meta", "transforms"].includes(k)) walkPlotly(v, depth + 1);
    else if (k === "marker" && v && typeof v === "object") walkPlotly(v, depth + 1);
  }
}

export function installPlotlyMath() {
  const P = typeof window !== "undefined" ? window.Plotly : null;
  if (!P || P.__fsimMath) return !!P;
  const wrap = (fn, prep) => {
    const orig = P[fn];
    if (typeof orig !== "function") return;
    P[fn] = function (...args) { try { prep(args); } catch { /* never block a plot over a label */ } return orig.apply(this, args); };
  };
  const dataLayout = (a) => { walkPlotly(a[1]); walkPlotly(a[2]); };
  wrap("newPlot", dataLayout);
  wrap("react", dataLayout);
  wrap("relayout", (a) => { const u = a[1]; if (u && typeof u === "object") for (const k of Object.keys(u)) if (typeof u[k] === "string" && /(^|\.)(text|title)$/.test(k)) u[k] = plotlyText(u[k]); else if (u[k] && typeof u[k] === "object") walkPlotly(u[k]); });
  wrap("update", (a) => { walkPlotly(a[1]); walkPlotly(a[2]); });
  wrap("restyle", (a) => { walkPlotly(a[1]); });
  P.__fsimMath = true;
  return true;
}
if (typeof window !== "undefined") {
  if (!installPlotlyMath()) window.addEventListener("DOMContentLoaded", () => installPlotlyMath(), { once: true });
}

// ---------------------------------------------------------------- typeset what is written to the DOM as plain text
// `el.textContent = "..."` (status lines, notes, captions) bypasses h(); one MutationObserver typesets those text nodes
// (`$...$` and the ASCII notation of autoMath) so no call site is left with a raw subscript or formula.
const OBS_SKIP_TAGS = new Set(["CODE", "PRE", "OPTION", "TEXTAREA", "INPUT", "SCRIPT", "STYLE", "TITLE", "KBD", "TEXT", "TSPAN", "SVG"]);
const OBS_QUICK = /[$_²^*×]/;

function typesetTextNode(n) {
  const p = n.parentNode;
  if (!p || p.nodeType !== 1 || OBS_SKIP_TAGS.has(String(p.tagName).toUpperCase())) return;
  if (p.closest?.(".math, .katex, svg, .js-plotly-plot, [data-no-math]")) return;
  const v = n.nodeValue;
  if (!v || !OBS_QUICK.test(v)) return;
  let frag = null;
  if (v.includes("$") && parseRich(v).some((s) => s.t !== "text")) frag = richText(v);
  else frag = autoMath(v);
  if (frag) p.replaceChild(frag, n);
}

function typesetTree(root) {
  if (root.nodeType === 3) { typesetTextNode(root); return; }
  if (root.nodeType !== 1 || root.closest?.(".math, .katex, svg, .js-plotly-plot")) return;
  const w = document.createTreeWalker(root, 4);
  const nodes = [];
  while (w.nextNode()) nodes.push(w.currentNode);
  for (const n of nodes) typesetTextNode(n);
}

export function installAutoMath(root = document.documentElement) {
  if (typeof MutationObserver === "undefined" || root.__fsimAutoMath) return;
  root.__fsimAutoMath = true;
  // Batched: text that changes every animation frame is typeset once it has been still for 150 ms, so playback never pays for KaTeX.
  const pending = new Set();
  let timer = null;
  const flush = () => {
    timer = null;
    const nodes = [...pending];
    pending.clear();
    for (const n of nodes) if (n.isConnected) typesetTree(n);
  };
  const mo = new MutationObserver((recs) => {
    for (const r of recs) {
      if (r.type === "characterData") pending.add(r.target);
      else for (const n of r.addedNodes) if (n.nodeType === 3 || (n.nodeType === 1 && !n.classList?.contains("math"))) pending.add(n);
    }
    if (pending.size && timer == null) timer = setTimeout(flush, 150);
    else if (pending.size) { clearTimeout(timer); timer = setTimeout(flush, 150); }
  });
  mo.observe(root, { childList: true, subtree: true, characterData: true });
  typesetTree(root);
}
if (typeof document !== "undefined" && document.documentElement) {
  installAutoMath();
  // the registry's ids -> symbols (so notes that name `eta_out` or `pol_factor` show the symbol); fetched once through api.js
  import("../api.js").then((a) => a.composeRegistry?.()).then((reg) => { if (registerSymbols(reg)) typesetTree(document.body); }).catch(() => { /* ids stay as written */ });
}
