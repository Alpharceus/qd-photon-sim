# qd-photon-sim

A Python simulator for electrically driven epitaxial quantum-dot (QD)
single-photon sources. It computes the second-order correlation at zero delay,
g2(0), the brightness (collected photon flux), and how both change with
temperature, and it provides an inverse-design mode that turns a g2(0) target
into required device parameters.

The core is headless numpy/scipy code. Three front ends sit on top of it:
a local web application (FSIM Studio), a Dear PyGui device designer, and a
Streamlit dashboard. None of the front ends contains physics.

## What it models

- **Statistics of the source.** Exciton-biexciton cascade emission with
  Lorentzian lineshapes and spectral filtering, injection background, and the
  resulting g2(0) and brightness as functions of temperature and drive.
- **Three levels of dynamics.**
  - Rate equations (the default, fast path).
  - A hand-written Lindblad master-equation core (`fsim_core/lindblad.py`)
    for density-matrix calculations.
  - Independent-boson-model (polaron) phonon lineshapes with a dot-geometry
    form factor (`fsim_core/qd_gf.py`).
- **Device tiers.**
  - Edge-emitting InP-dot devices at elevated temperature (ridge waveguide,
    p-i-n transport, confinement levels, CW correlation).
  - A planar nitride (InGaN/GaN) cavity tier, including polar and nonpolar
    orientation and the quantum-confined Stark shift.
  - A nitride nanowire tier (vertical and horizontal geometries, injector and
    surface models).
  - A resonantly driven Fock-state source tier (`fsim_core/resonant_source.py`)
    built on the Lindblad core.
- **Thermal and cavity layers.** A spreading-resistance thermal stack with
  junction temperature and runaway detection, and cavity tracking, Purcell
  enhancement and collection gain.
- **Inverse design.** `fsim_core/spec.py` converts a g2(0) target into the
  required linewidth, gain, mode volume and background budget. It produces
  necessary conditions, not a cavity geometry.

Scope limits: no growth modelling, and no indistinguishability metrics outside
the resonant Fock-state tier.

## Provenance discipline

Every parameter in a card carries one of four tags:

| Tag | Meaning |
|---|---|
| `[V]` | Verified against the cited paper |
| `[DR]` | Derived from published data |
| `[E]` | Estimate or class range |
| `[A]` | Assumption |

Every output inherits the **widest** (least certain) tag found anywhere in its
input chain, and the card schema enforces this. Ranges are swept, never
averaged. Inputs tagged `[A]` can be evaluated as envelopes (bands and
sensitivity rankings) instead of point predictions. Literature numbers in the
code are cited by author, journal and year.

## Validation classes

Each check in `verify/` belongs to exactly one of five classes. Keeping them
apart avoids presenting a calibration as if it were a prediction.

| Class | Meaning |
|---|---|
| Numerical verification | One method checked against another (for example a closed form against Monte Carlo) under the same physical assumptions; no external data. |
| Source transcription | A number, formula or claim in the code checked against the paper that states it, or a structural check of the evidence ledger. |
| Calibration | Free parameters fitted to data and evaluated on that same data. Not held-out validation. |
| Unfitted comparison | A previously fixed model or class range compared with published data it was not fitted to, with no parameter adjusted. Still not held-out, because the target was known when the ranges were chosen. |
| Held-out prediction | A model checked against data withheld from fitting and range selection. |

Published datasets used for calibration or comparison are cited by reference
only: Chatzarakis et al., Phys. Rev. Applied 20, 034011 (2023); Laferriere
et al., Nano Lett. 23, 962 (2023); Reischle et al., Appl. Phys. Lett. 97,
143513 (2010) and Opt. Express 16, 12771 (2008); Zhao et al. (laser model);
Deshpande et al. (2014, nitride comparison).

Two points about the data:

- The comparisons run against digitizations of published figures plus values
  printed in the papers, not author-released datasets. The source PDFs are not
  distributed, so reproducing a fit from scratch means obtaining the papers and
  re-extracting the data.
- The rule that an anchor claim rest on two independently verified sources is
  a policy for this repository's evidence ledger, not a general scientific
  requirement.

## Installation

```
pip install -r requirements.txt
```

Requires Python 3.9 or newer. Dependencies: numpy, scipy, matplotlib, PyYAML,
streamlit, plotly, dearpygui. Node.js is needed only for the browser-based
verifiers (`verify/*.cjs`); everything else runs without it.

## Running

### FSIM Studio

A local web application over `fsim_core`. The server binds to 127.0.0.1 only.

```
run-studio.bat                 # Windows; native window
run-studio.bat --browser       # default browser instead
./run-studio.sh                # POSIX equivalent

python -m fsim_studio [--browser] [--no-window] [--port N]
```

Workspaces:

- **Overview**: landing page summarising each stored sweep campaign under
  `out/` (verdict line, counts, generated figures and CSVs).
- **Designer**: edit a device (dot, cavity, drive, injection, thermal stack),
  validate it, and run it, with optional envelope mode for ranged inputs.
- **Compare**: place several runs side by side and export a comparison report.
- **Results**: browse stored sweep campaigns under `out/`, with filtering and
  tables.
- **Explain**: step-by-step breakdowns of how a result was computed, plus model
  explorers for the filtered cascade and drive statistics, phonon lineshapes,
  transport and nitride Stark shift, and the QD laser model.
- **Library**: literature parameter cards with their tags and published-data
  overlays.
- **Story**: scripted, scene-based walkthroughs, including 3D device scenes.

The HTTP API is documented in `fsim_studio/API.md`.

### Legacy front ends

```
python fsim_gui/designer.py        # Dear PyGui device designer
streamlit run fsim_gui/app.py      # Streamlit validation dashboard
```

The designer provides a block-diagram device chain, a fabrication-stack editor
with a live cross-section, envelope mode, A/B/C comparison slots and report
bundles. Launchers are `run-designer.bat|sh` and `run-dashboard.bat|sh`.

### Scripts

`scripts/` holds batch drivers that write figures and CSVs under `out/`, for
example `run_phase0.py` to `run_phase3.py` (validation fits, thermal maps,
requirement envelopes), `run_spec.py` (inverse-design spec sheets),
`run_rt_edge.py` (edge-emitter sweep), `run_nitride_cavity.py`,
`run_nitride_nanowire.py` and `run_fock_source.py`.

## Verification

```
python verify/verify_fsim.py && python verify/audit_physics.py && python verify/verify_materials.py && python verify/verify_cw_g2.py && python verify/verify_dot_levels.py && python verify/verify_waveguide.py && python verify/verify_transport.py
```

Each script exits 0 only if all of its checks pass. `verify/` also contains
per-module suites for the Lindblad core, the nitride and nanowire tiers, the
resonant source, and Studio (`verify_studio_*.py` and the `.cjs` browser
verifiers). Modules with independent Monte Carlo or finite-difference second
methods are checked against them. `tests/oracles/` holds independent
electromagnetic reference formulas used only as test oracles; they are never
imported by the library.

## Repository layout

| Path | Contents |
|---|---|
| `fsim_core/` | Headless physics and assembly: parameter-card schema (`card.py`), spectral filtering, drive statistics, thermal stack, cavity, `device.py` (design blocks, `evaluate()`, `evaluate_envelope()`), `spec.py` (inverse design), `lindblad.py`, `qd_gf.py`, `resonant_source.py`, and the edge-emitter and nitride tier modules |
| `fsim_studio/` | FSIM Studio server, job and cache handling, and web front end (`web/`); API reference in `API.md` |
| `fsim_gui/` | Dear PyGui designer (`designer.py`) and Streamlit dashboard (`app.py`) |
| `fsim_viz/` | Matplotlib figure factory; each figure ships with its CSV |
| `fsim_theme/` | Shared design tokens, plot styles and fonts for all front ends |
| `cards/` | Parameter cards for published datasets and `*-design.yaml` device designs |
| `scripts/` | Batch drivers for sweeps, fits and spec sheets |
| `verify/` | Regression suites, physics audit against published values, and evidence-ledger data |
| `tests/oracles/` | Independent reference implementations used only by tests |
| `skills/` | Citation-verification helpers (see `THIRD_PARTY.md`) |
| `data/` | Small input tables |
| `out/` | Generated report bundles (figures and CSVs), regenerable from `scripts/` |

## License

MIT; see `LICENSE`. The files under `skills/` are adapted from third-party
code under a different license; see `THIRD_PARTY.md`.
