#!/usr/bin/env node
/*
 * presentation2/render_pptx.js -- pptxgenjs renderer per the contract in
 * presentation2/SCHEMA.md. Reads every presentation2/sections/NN_*.json
 * (section "00" skipped unless --include-sample), and writes
 * out/presentation2/qd_physics_2h.pptx.
 *
 * Never special-cases section content: every layout branch works from the
 * generic slide object the schema defines. Equation images must already be
 * rendered by presentation2/mathimg.py (build.py's "equations" build step)
 * before this runs -- this script only looks them up by the same content
 * hash mathimg.py uses, it never calls Python.
 *
 * CLI: node presentation2/render_pptx.js [--include-sample] [--sections-dir DIR]
 */
"use strict";

const fs = require("fs");
const path = require("path");
const crypto = require("crypto");

let pptxgen;
try {
  pptxgen = require("pptxgenjs");
} catch (e) {
  pptxgen = require(path.join(__dirname, "node_modules", "pptxgenjs"));
}

const HERE = __dirname;
const ROOT = path.resolve(HERE, "..");
const EQ_DIR = path.join(HERE, "figures", "out", "eq");
const OUT_DIR = path.join(ROOT, "out", "presentation2");
const OUT_PPTX = path.join(OUT_DIR, "qd_physics_2h.pptx");

// ---- CLI args -------------------------------------------------------------

const argv = process.argv.slice(2);
const INCLUDE_SAMPLE = argv.includes("--include-sample");
let SECTIONS_DIR = path.join(HERE, "sections");
const sdIdx = argv.indexOf("--sections-dir");
if (sdIdx !== -1 && argv[sdIdx + 1]) {
  SECTIONS_DIR = path.resolve(argv[sdIdx + 1]);
}

// ---- theme ------------------------------------------------------------

const ACCENT_PALETTE = [
  "1E2761", // navy
  "065A82", // deep blue
  "2C5F2D", // forest
  "6D2E46", // berry
  "028090", // teal
  "B85042", // terracotta
  "36454F", // charcoal
  "990011", // cherry
];
const INK = "212121";
const MUTED = "6B6B6B";
const LIGHT_MUTED = "D6D9E6";
const WHITE = "FFFFFF";

const TITLE_FONT = "Cambria";
const BODY_FONT = "Calibri";

function accentFor(sectionId) {
  const n = parseInt(sectionId, 10);
  return ACCENT_PALETTE[Number.isFinite(n) ? ((n % ACCENT_PALETTE.length) + ACCENT_PALETTE.length) % ACCENT_PALETTE.length : 0];
}

// ---- geometry (LAYOUT_WIDE = 13.333 x 7.5 in) ------------------------------

const MARGIN_X = 0.5;
const CONTENT_X = MARGIN_X;
const CONTENT_W = 13.333 - 2 * MARGIN_X;
const TITLE_Y = 0.4;
const TITLE_H = 0.85;
const BODY_Y = 1.4;
const BODY_BOTTOM = 6.85;
const BODY_H = BODY_BOTTOM - BODY_Y;
const FOOTER_Y = 7.05;
const FOOTER_H = 0.3;

// ---- section loading --------------------------------------------------

function loadSections() {
  const files = fs.readdirSync(SECTIONS_DIR)
    .filter((f) => /^\d{2}_.*\.json$/.test(f))
    .sort();
  const sections = [];
  for (const f of files) {
    const isSample = f.startsWith("00_");
    if (isSample && !INCLUDE_SAMPLE) continue;
    const data = JSON.parse(fs.readFileSync(path.join(SECTIONS_DIR, f), "utf8"));
    sections.push(data);
  }
  return sections;
}

// ---- equation image lookup (must match mathimg.py's hash exactly) ---------

function equationImagePath(latex) {
  const hash = crypto.createHash("sha256").update(latex, "utf8").digest("hex").slice(0, 16);
  return path.join(EQ_DIR, `eq_${hash}.png`);
}

// ---- inline-math-only bullet detection --------------------------------

const ONLY_MATH_RE = /^\$([^$]+)\$$/;

function onlyMathLatex(bullet) {
  const m = ONLY_MATH_RE.exec(bullet.trim());
  return m ? m[1] : null;
}

// ---- small text helpers -------------------------------------------------

function bulletTextRuns(items, opts) {
  return items.map((text, i) => ({
    text,
    options: Object.assign({ bullet: true, breakLine: i < items.length - 1 }, opts),
  }));
}

function provenanceText(repoNumbers) {
  if (!repoNumbers) return null;
  const parts = Object.entries(repoNumbers).map(([key, entry]) => {
    const val = typeof entry.value === "number" ? String(entry.value) : String(entry.value);
    return `${key} = ${val} (${entry.file})`;
  });
  return parts.length ? "Source: " + parts.join("  |  ") : null;
}

// ---- bulleted-column renderer (handles math-only lines as images) ---------

function renderBulletsColumn(slide, bullets, box, opts) {
  opts = opts || {};
  const fontSize = opts.fontSize || 20;
  const color = opts.color || INK;

  const groups = [];
  let current = [];
  for (const b of bullets) {
    const latex = onlyMathLatex(b);
    if (latex !== null) {
      if (current.length) {
        groups.push({ type: "text", items: current });
        current = [];
      }
      groups.push({ type: "eq", latex });
    } else {
      current.push(b);
    }
  }
  if (current.length) groups.push({ type: "text", items: current });

  const EQ_H = 0.85;
  const nText = groups.filter((g) => g.type === "text").reduce((s, g) => s + g.items.length, 0);
  const nEq = groups.filter((g) => g.type === "eq").length;
  const eqTotalH = nEq * EQ_H;
  const textTotalH = Math.max(box.h - eqTotalH, 0.2);

  let y = box.y;
  for (const g of groups) {
    if (g.type === "eq") {
      const imgPath = equationImagePath(g.latex);
      if (fs.existsSync(imgPath)) {
        slide.addImage({
          path: imgPath,
          x: box.x, y, w: box.w, h: EQ_H,
          sizing: { type: "contain", w: box.w, h: EQ_H },
        });
      }
      y += EQ_H;
    } else {
      const h = textTotalH * (g.items.length / Math.max(nText, 1));
      slide.addText(bulletTextRuns(g.items, { fontSize, color, fontFace: BODY_FONT }), {
        x: box.x, y, w: box.w, h,
        valign: "top", margin: 0, paraSpaceAfter: 10, autoFit: true,
      });
      y += h;
    }
  }
}

// ---- footer + notes (applied to every slide) -------------------------

function finishSlide(slide, ctx) {
  const { sectionTitle, pageNum, repoNumbers, dark, notes } = ctx;
  const footerColor = dark ? LIGHT_MUTED : MUTED;

  slide.addText(sectionTitle, {
    x: MARGIN_X, y: FOOTER_Y, w: 4.0, h: FOOTER_H,
    fontSize: 10, color: footerColor, fontFace: BODY_FONT, margin: 0, valign: "middle",
  });

  const prov = provenanceText(repoNumbers);
  if (prov) {
    slide.addText(prov, {
      x: 4.6, y: FOOTER_Y, w: 4.13, h: FOOTER_H,
      fontSize: 8, italic: true, color: footerColor, fontFace: BODY_FONT,
      align: "center", margin: 0, valign: "middle",
    });
  }

  slide.addText(String(pageNum), {
    x: 11.83, y: FOOTER_Y, w: 1.0, h: FOOTER_H,
    fontSize: 10, color: footerColor, fontFace: BODY_FONT,
    align: "right", margin: 0, valign: "middle",
  });

  if (notes) slide.addNotes(notes);
}

// ---- layout renderers ---------------------------------------------------

function addTitleText(slide, text, opts) {
  slide.addText(text, Object.assign({
    x: CONTENT_X, y: TITLE_Y, w: CONTENT_W, h: TITLE_H,
    fontSize: 32, bold: true, fontFace: TITLE_FONT, margin: 0, valign: "top",
  }, opts));
}

function layoutTitle(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  s.addText(slide.title, {
    x: 0.8, y: 2.4, w: 13.333 - 1.6, h: 1.6,
    fontSize: 40, bold: true, color: accent, fontFace: TITLE_FONT,
    align: "center", valign: "middle", margin: 0,
  });
  if (slide.subtitle) {
    s.addText(slide.subtitle, {
      x: 1.3, y: 4.05, w: 13.333 - 2.6, h: 0.8,
      fontSize: 20, italic: true, color: MUTED, fontFace: BODY_FONT,
      align: "center", valign: "top", margin: 0,
    });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutSectionDivider(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: accent };
  s.addText(slide.title, {
    x: 0.8, y: 2.6, w: 13.333 - 1.6, h: 1.6,
    fontSize: 40, bold: true, color: WHITE, fontFace: TITLE_FONT,
    align: "center", valign: "middle", margin: 0,
  });
  if (slide.subtitle) {
    s.addText(slide.subtitle, {
      x: 1.3, y: 4.25, w: 13.333 - 2.6, h: 0.7,
      fontSize: 18, italic: true, color: LIGHT_MUTED, fontFace: BODY_FONT,
      align: "center", valign: "top", margin: 0,
    });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: true, notes: slide.notes });
}

function layoutQuote(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  s.addText(`"${slide.title}"`, {
    x: 1.2, y: 2.0, w: 13.333 - 2.4, h: 2.6,
    fontSize: 30, italic: true, bold: true, color: accent, fontFace: TITLE_FONT,
    align: "center", valign: "middle", margin: 0,
  });
  if (slide.subtitle) {
    s.addText(`- ${slide.subtitle}`, {
      x: 1.5, y: 4.7, w: 13.333 - 3.0, h: 0.6,
      fontSize: 16, color: MUTED, fontFace: BODY_FONT,
      align: "center", valign: "top", margin: 0,
    });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutBullets(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  if (slide.subtitle) {
    s.addText(slide.subtitle, { x: CONTENT_X, y: TITLE_Y + TITLE_H, w: CONTENT_W, h: 0.35, fontSize: 16, italic: true, color: MUTED, fontFace: BODY_FONT, margin: 0 });
  }
  if (slide.bullets && slide.bullets.length) {
    renderBulletsColumn(s, slide.bullets, { x: CONTENT_X, y: BODY_Y, w: CONTENT_W, h: BODY_H });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function addFigureWithCaption(s, figure, box) {
  const imgPath = path.join(ROOT, figure.path);
  if (fs.existsSync(imgPath)) {
    const capH = figure.caption ? 0.4 : 0;
    s.addImage({
      path: imgPath,
      x: box.x, y: box.y, w: box.w, h: box.h - capH,
      sizing: { type: "contain", w: box.w, h: box.h - capH },
    });
    if (figure.caption) {
      s.addText(figure.caption, {
        x: box.x, y: box.y + box.h - capH, w: box.w, h: capH,
        fontSize: 12, italic: true, color: MUTED, fontFace: BODY_FONT,
        align: "center", valign: "top", margin: 0,
      });
    }
  }
}

function layoutBulletsFigure(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  const gap = 0.4;
  const leftW = CONTENT_W * 0.45;
  const rightW = CONTENT_W - leftW - gap;
  if (slide.bullets && slide.bullets.length) {
    renderBulletsColumn(s, slide.bullets, { x: CONTENT_X, y: BODY_Y, w: leftW, h: BODY_H });
  }
  if (slide.figure) {
    addFigureWithCaption(s, slide.figure, { x: CONTENT_X + leftW + gap, y: BODY_Y, w: rightW, h: BODY_H });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutFigure(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  if (slide.figure) {
    addFigureWithCaption(s, slide.figure, { x: CONTENT_X + 1.0, y: BODY_Y, w: CONTENT_W - 2.0, h: BODY_H });
  }
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutEquation(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  const eqs = slide.equations || [];
  const rowH = BODY_H / Math.max(eqs.length, 1);
  eqs.forEach((eq, i) => {
    const y = BODY_Y + i * rowH;
    const capH = 0.4;
    const imgPath = equationImagePath(eq.latex);
    if (fs.existsSync(imgPath)) {
      s.addImage({
        path: imgPath,
        x: CONTENT_X + 1.0, y, w: CONTENT_W - 2.0, h: rowH - capH,
        sizing: { type: "contain", w: CONTENT_W - 2.0, h: rowH - capH },
      });
    }
    if (eq.caption) {
      s.addText(eq.caption, {
        x: CONTENT_X + 1.0, y: y + rowH - capH, w: CONTENT_W - 2.0, h: capH,
        fontSize: 12, italic: true, color: MUTED, fontFace: BODY_FONT,
        align: "center", valign: "top", margin: 0,
      });
    }
  });
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutTwoColumn(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  const cols = slide.columns || [];
  const gap = 0.5;
  const colW = (CONTENT_W - gap * (cols.length - 1)) / Math.max(cols.length, 1);
  cols.forEach((col, i) => {
    const x = CONTENT_X + i * (colW + gap);
    s.addText(col.heading, {
      x, y: BODY_Y, w: colW, h: 0.5,
      fontSize: 20, bold: true, color: accent, fontFace: TITLE_FONT, margin: 0,
    });
    renderBulletsColumn(s, col.bullets || [], { x, y: BODY_Y + 0.55, w: colW, h: BODY_H - 0.55 });
  });
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

function layoutTable(pres, slide, sectionTitle, pageNum, accent) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  addTitleText(s, slide.title, { color: accent });
  const table = slide.table || { header: [], rows: [] };
  const headerRow = table.header.map((h) => ({
    text: h,
    options: { bold: true, color: WHITE, fill: { color: accent }, fontFace: BODY_FONT, fontSize: 14 },
  }));
  const bodyRows = table.rows.map((row, ri) =>
    row.map((cell) => ({
      text: String(cell),
      options: {
        color: INK, fontFace: BODY_FONT, fontSize: 14,
        fill: { color: ri % 2 === 0 ? "F7F7F7" : WHITE },
      },
    }))
  );
  s.addTable([headerRow, ...bodyRows], {
    x: CONTENT_X, y: BODY_Y, w: CONTENT_W, h: BODY_H,
    fontSize: 14, border: { type: "solid", color: "DDDDDD", pt: 0.75 },
    autoPage: false, valign: "middle",
  });
  finishSlide(s, { sectionTitle, pageNum, repoNumbers: slide.repo_numbers, dark: false, notes: slide.notes });
}

const LAYOUTS = {
  title: layoutTitle,
  "section-divider": layoutSectionDivider,
  quote: layoutQuote,
  bullets: layoutBullets,
  "bullets+figure": layoutBulletsFigure,
  figure: layoutFigure,
  equation: layoutEquation,
  "two-column": layoutTwoColumn,
  table: layoutTable,
};

// ---- main ---------------------------------------------------------------

function main() {
  const sections = loadSections();

  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.defineSlideMaster({ title: "QD_MASTER", background: { color: WHITE } });

  let pageNum = 0;
  for (const section of sections) {
    const accent = accentFor(section.section);
    for (const slide of section.slides || []) {
      pageNum += 1;
      const renderer = LAYOUTS[slide.layout];
      if (!renderer) {
        throw new Error(`unknown layout '${slide.layout}' on slide ${slide.id}`);
      }
      renderer(pres, slide, section.title, pageNum, accent);
    }
  }

  fs.mkdirSync(OUT_DIR, { recursive: true });
  pres.writeFile({ fileName: OUT_PPTX }).then(() => {
    console.log(`wrote ${OUT_PPTX} (${pageNum} slides)`);
  }).catch((err) => {
    console.error(err);
    process.exit(1);
  });
}

main();
