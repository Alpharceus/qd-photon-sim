# presentation2: slide data schema (shared by every section worker and both renderers)

The 2-hour deck "How quantum dots emit single photons, from the physics to the electrically driven simulator" is data-driven. Each section is one JSON file `presentation2/sections/NN_<slug>.json`; the renderers (`render_pptx.js`, `render_html.py`) turn the same data into `out/presentation2/qd_physics_2h.pptx` and `out/presentation2/index.html`. Figures are PNG files produced by scripts under `presentation2/figures/` (matplotlib; every figure script must be re-runnable and deterministic) and referenced by relative path. Equations are LaTeX strings; `mathimg.py` renders them to PNG with matplotlib mathtext (no TeX install) at build time.

## Section file

```json
{
  "section": "02",
  "title": "How a quantum dot emits",
  "minutes": 28,
  "slides": [ <slide>, ... ]
}
```

## Slide object

```json
{
  "id": "02-05",
  "layout": "title | bullets | bullets+figure | figure | equation | equation+figure | two-column | table | bullets+table | section-divider | quote",
  "title": "Fermi's golden rule for spontaneous emission",
  "subtitle": "optional one line",
  "bullets": ["plain text; inline math allowed as $...$ (rendered as text in HTML, as an image in pptx if the line is only math)", "..."],
  "equations": [{"latex": "\\Gamma_{rad} = \\frac{n \\omega^3 |d_{cv}|^2}{3 \\pi \\varepsilon_0 \\hbar c^3}", "caption": "spontaneous emission rate in a medium of index n"}],
  "figure": {"path": "presentation2/figures/out/02_05_golden_rule.png", "caption": "one line", "script": "presentation2/figures/fig_02_05_golden_rule.py"},
  "columns": [{"heading": "left", "bullets": ["..."]}, {"heading": "right", "bullets": ["..."]}],
  "table": {"header": ["regime", "ideality n", "physics"], "rows": [["...", "...", "..."]]},
  "notes": "Speaker notes: the full explanation a lecturer would say, 120-250 words, with citations in brackets [Fox, Quantum Optics, ch. 8] [Michler ed., Quantum Dots for Quantum Information Technologies (2017)] [Reischle et al., Opt. Express 16, 12771 (2008)]. Numbers that come from the repository must name the file they come from.",
  "sources": ["Fox 2006 ch. 8", "Reischle 2008 OE 16 12771"],
  "repo_numbers": {"S_300K_gainp": {"value": 0.00461745, "file": "out/rt_edge/verdict.md"}}
}
```

Rules:
- Every slide has `id`, `layout`, `title`, `notes`, `sources` (at least one). `bullets` at most 6 items, each at most ~18 words. Body text large; one idea per slide.
- Equations: at most 3 per slide; each with a caption naming every symbol not already defined on that slide. Equations are rendered on EVERY layout that carries an `equations` array (below the bullets or beside the figure), not only on `equation` slides; likewise `bullets` are rendered on EVERY layout that carries them, including `equation` and `equation+figure` (bullets above the equation block).
- Speaker notes: 120-250 words (the validator enforces 115-255 as hard bounds). Figures: exactly 1600x900 or 1200x1200 px (validator-enforced).
- Figures: generated, never downloaded; the script must write the PNG at 1600x900 px (16:9) or 1200x1200 for square, `dpi=150`, white background, fonts >= 14 pt. Use repository modules for physics curves where they exist (e.g. `fsim_core.linewidth.gamma_anchor`, `fsim_core.dot_levels`, `fsim_core.transport`, `fsim_core.cw_g2`), otherwise closed-form textbook expressions with the source named in the caption.
- `repo_numbers`: any number taken from the repository must be listed here so the renderer can print a provenance footnote and the validator can verify it. Two forms: (a) a number that literally appears in a data/text file: `{"value": 0.0046, "file": "out/rt_edge/verdict.md"}` (the validator searches the file within 1e-6 relative); (b) a computed number: `{"value": 0.4042, "file": "fsim_core/cw_g2.py", "how": "<python expression>"}` where `how` is a single Python expression evaluated with the repository root on sys.path (e.g. `"__import__('fsim_core.cw_g2',fromlist=['x']).g2_cw_zero(r_ns=10, gamma_X_ns=1.0, gamma_XX_ns=1.4, t_X=0.5, t_XX=0.1)"`); the validator evaluates it in a subprocess and compares within 1e-6 relative (1e-3 for values that depend on a sweep CSV). Every `how` must be deterministic and run in under 30 s.
- No invented references. Textbook citations allowed: Fox "Quantum Optics" (OUP 2006); Bimberg, Grundmann, Ledentsov "Quantum Dot Heterostructures" (Wiley 1999); Michler (ed.) "Quantum Dots for Quantum Information Technologies" (Springer 2017); Sze & Ng "Physics of Semiconductor Devices" 3rd ed.; Coldren, Corzine, Masanovic "Diode Lasers and Photonic Integrated Circuits" 2nd ed.; Loudon "The Quantum Theory of Light"; Yu & Cardona "Fundamentals of Semiconductors"; Chuang "Physics of Photonic Devices"; Vurgaftman, Meyer, Ram-Mohan JAP 89, 5815 (2001); plus the six digested papers in `../_goal/paper_digests.md` and the references in `.workers/briefs/digest-qd-quantum-treatment.md` (the owner's slides).
- Pace target: ~70 s per slide; a 28-minute section has 22-26 slides.

## Renderer contract

- `python presentation2/build.py` : validates every section JSON against this schema (`presentation2/validate_sections.py`), runs every figure script whose PNG is missing or older than the script, renders equations, then runs `node presentation2/render_pptx.js` and `python presentation2/render_html.py` with their default output directory `out/presentation2`, then `python ~/.claude/skills/pptx/scripts/office/validate.py out/presentation2/qd_physics_2h.pptx`; `verify/verify_presentation2.py` passes `--out-dir out/presentation2_selftest` to both renderers so its fixture checks never overwrite the deliverables.
- `render_pptx.js`: pptxgenjs, `LAYOUT_WIDE` (13.33 x 7.5 in), one master with a slim footer (section title, slide number, provenance footnote if `repo_numbers`), title 32-40 pt, body 20-24 pt, dark-on-light theme, accent colour per section; equations placed as images centred with captions; speaker notes from `notes`.
- `render_html.py`: single self-contained file (base64 images, inline CSS, no external scripts), one `<section>` per slide, MathJax-free (equations are images; inline `$...$` shown in a monospace span), a sticky table of contents, dark-mode media query, no horizontal scroll.
