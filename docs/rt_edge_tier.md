# RT edge-emitter tier (branch `rt-edge-emitter`)

This tier asks whether the F-series simulator's InP-dot platform can be
pushed to a room-temperature (300 K heatsink), electrically driven,
edge-emitting single-photon source, and answers with a reproducible
acceptance sweep rather than a hand-picked corner. This document is the
reader's map: what was asked, what came back, why the material platform in
the code differs from the literal wording of the original brief, where the
new physics lives, how to reproduce every number, and what is still missing.

Every number quoted below was read from the file named next to it at the
time this document was written (2026-09-05), by re-running the corresponding
script or `cat`-ing the corresponding output file. None is retyped from
memory.

## Goal and honest outcome

**Goal.** Evaluate, over the full contract-declared parameter grid
(`docs/rt_edge_contract.md`), whether either of the two RT edge-emitter
design cards (`cards/edge-inp-gaasp-design.yaml`,
`cards/edge-inp-gainp-design.yaml`) has a physically eligible operating
corner with pulsed intrinsic `g2(0) < 0.5` at a 300 K heatsink, backed by
paper-verified evidence for every literature anchor the claim depends on.

**Outcome: FAIL**, honestly. From `out/rt_edge/verdict.md`:

```
VERDICT: FAIL g2_min=0.3648 g2_median=0.773 median_pass=false coverage=0 eligible=54/54 evidence=incomplete conditional=false
```

Per-card statistics (same file):

| card | role | g2_pulsed min/median | g2_cw0 min/median | g2_cw0_raw min/median | eligible | favorable rows |
|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.4117 / 0.8555 | 0.4026 / 0.9443 | 0.9784 / 0.9992 | 27/27 | 0 |
| edge-inp-gainp-design | fallback | 0.3648 / 0.7325 | 0.3413 / 0.905 | 0.8907 / 0.9926 | 27/27 | 0 |

The pooled pulsed minimum (0.3648, from the fallback card) does sit under the
0.5 headline threshold, but the median (0.773) does not, and the run reports
`favorable_rows=0` for both cards: no single grid point passes every gate at
once (pulsed intrinsic AND CW intrinsic AND CW raw, same row). The sweep's
own fail reasons (`out/rt_edge/verdict.md`, "Fail reasons"):

- `evidence_incomplete`
- `no_favorable_corner_metrics`

**Literature ceiling.** Quoting `out/rt_edge/verdict.md`, "Literature
ceiling" section: "No published III-V quantum dot has demonstrated g2(0) <
0.5 at 300 K under electrical driving; the best reported room-temperature
values sit around 0.5-0.57." The sweep's 0.3648 is reported *against* that
ceiling, not as a claim of having exceeded it in a real device — it is a
model prediction under `[A]`/`[E]` inputs, not a measurement.

**The two single-source evidence claims.** The acceptance gate requires two
independent verified primary sources for every literature anchor a PASS
would depend on. Two claims currently have zero, per
`python verify/verify_rt_edge_contract.py` (76/78 checks passed; the 2
failures are exactly these) and `python verify/verify_rt_edge_papers.py`
(14/16 checks passed; same 2 failures):

- `reischle2008_g2_80K` — the 80 K electrically driven InP-dot `g2(0)` result
  (Reischle et al., APL 97, 143513 (2010)) has no second independent 80 K
  electrical single-dot measurement. Reischle et al., Opt. Express 16, 12771
  (2008) is the only other electrically driven InP-dot g2 paper on file, but
  it reports 20-40 K (Joule-heated), not 80 K, so it cannot stand in.
- `hkust_inp_gaasp_wavelength` — the HKUST InP/GaAsP 668 nm single-photon
  anchor. A second, DOI-sourced wavelength candidate exists (Luo et al.,
  Optics Express 30, 40750 (2022), 750 nm InP/GaAsP QD laser on Si), but
  promoting it to `verified` would make `verify/verify_rt_edge_cards.py`
  reject `cards/edge-inp-gaasp-design.yaml`'s `emission.lambda_nm` (668 nm),
  which is deliberately pinned to the only tabulated phosphide refractive-
  index point rather than to the unconfirmed HKUST target wavelength. Both
  anchors are kept `missing` (see `docs/rt_edge_contract.md`, "Evidence
  status") and both block PASS.

## The material platform switch, and why

The original brief described the target device as "InP cladding, GaAsP
active region." Taken literally that is inverted: InP (`E_g` = 1.35 eV at
300 K) has a *smaller* gap than any GaAs<sub>1-x</sub>P<sub>x</sub> alloy
(1.42-2.27 eV), so InP cannot act as a cladding/barrier confining a GaAsP
well, and coherent GaAsP-on-InP is tensile-mismatched by 3.7% or more with no
strain-balancing partner and no literature precedent (`../_goal/
materials_research.md`, sections 0 and 1 — conclusions cited here only,
since that file lives outside this repository).

The platform actually implemented, and the one the two design cards use, is
the swapped reading: **InP quantum dots embedded in a GaAs<sub>0.6</sub>
P<sub>0.4</sub> well, with (Al,Ga)InP-class (or AlGaAs-class) barriers and
cladding, on a GaAs substrate.** Three independent reasons, all from
`../_goal/materials_research.md`'s conclusions:

- **Lattice mismatch.** InP on GaAs<sub>0.6</sub>P<sub>0.4</sub> is +5.34%
  mismatched — large, but the same order as InAs/GaAs (7.2%, a system that
  routinely forms dots) and closer than InAs/InP (3.2%) once density and
  height are controlled by growth (V/III ratio). GaAsP-on-InP, the literal
  reading, has no lattice-matched or even close-to-matched composition at
  all: the mismatch is never below 3.7% and there is no strain-compensating
  partner on the InP side.
- **Band alignment and carrier confinement.** From the model-solid valence-
  band offsets, the unstrained InP conduction band sits about 0.5 eV below
  GaAs<sub>0.6</sub>P<sub>0.4</sub>'s and the unstrained InP valence band
  sits only about 0.05 eV above it — deep electron confinement, weak,
  strain-set hole confinement. That is workable (it is the same qualitative
  picture as the already-published InP/GaInP system) and matches published
  type-I InP/GaAs<sub>0.65</sub>P<sub>0.35</sub> dots. The literal
  reading is not workable: GaAsP is the wide-gap partner, so InP (not GaAsP)
  would collect carriers — the device as literally described could not
  confine carriers in the stated active region at all.
- **Transparency.** GaAs<sub>0.6</sub>P<sub>0.4</sub> and (Al,Ga)InP
  claddings are transparent at the dots' ~650-670 nm emission (their band
  gaps sit well above the photon energy at these compositions, keeping P
  content at or below x = 0.4 to stay direct-gap and avoid Γ-X leakage);
  InP itself has a smaller gap than the emission photon energy at these
  wavelengths and would reabsorb rather than clad.

This is exactly the material system the wafers described in the simulator's
own literature base are grown from (Gu et al., Opt. Express 33, 23732
(2025); Luo et al. / HKUST, Opt. Express (2022); both InP-dot-in-GaAsP-well,
AlGaAs-barrier laser structures), so the code follows the growable platform
rather than the literal prose.

## Module map

| Module | Role |
|---|---|
| `fsim_core/materials.py` | III-V material database: Vurgaftman/Meyer/Ram-Mohan (2001) binary parameters, alloy interpolation, strain, band alignment, refractive index. Every RT-edge band-structure and index number traces here. |
| `fsim_core/linewidth.py` | Literature-anchored 300 K linewidth model `Gamma(T) = Gamma_0 + a_ac*T + b_LO/(exp(E_LO/kT)-1)`, opt-in via `dot.linewidth="anchored"`, replacing the class-fit extrapolation. |
| `fsim_core/dot_levels.py` | Module L: single-band effective-mass confinement solver (1D well, 2D disk, BenDaniel-Duke), literature-class reproductions, retention mapping. `verify/verify_dot_levels.py` is its regression + "known deviations" ledger (see below). |
| `fsim_core/transport.py` | Module T: p-i-n diode model of the actual layer stack (carrier statistics, junction voltage/temperature, injection efficiency), replacing free `[A]` drive inputs with a computed `b_e(I, T)`. |
| `fsim_core/waveguide.py` | Scalar effective-index slab solver for the ridge waveguide, guided-mode Purcell factor, edge-emission collection. |
| `fsim_core/cw_g2.py` | Module D-CW: continuous-wave rate-equation `g2(tau)`, Poissonian background, detector-IRF convolution — the CW/DC analogue of the pulsed peak-area `g2(0)`. |
| `fsim_core/device.py` | `DeviceDesign` opt-in blocks and `evaluate()`/`evaluate_envelope()`. New fields: `dot.linewidth` (`"class"` legacy / `"anchored"`), `ret.mode` (`"proxy"` legacy / `"confinement"`), `drive.mode` (`"EL"`/`"PL"` legacy / `"EL-transport"`), `drive.cw` (bool, default `False`), `emission.type` (`"none"` legacy / `"edge"`), `filter.track_material` (`""` legacy hard-coded GaAs / `"dot"`/`"matrix"`), `aperture.compose` (bool, default `False`). Every default reproduces the legacy path (Lemma 1: brightness/collection terms never lower intrinsic `g2`). |
| `cards/edge-inp-gaasp-design.yaml`, `cards/edge-inp-gainp-design.yaml` | The two design cards the sweep evaluates: primary (InP/GaAsP/AlGaAs) and fallback (InP/GaInP/AlGaInP) platforms. |
| `scripts/run_rt_edge.py` | The acceptance sweep: a thin caller of `device.evaluate()` over the contract grid, writing `sweep.csv`, `envelope.png`, `verdict.md`, `manifest.json` to `out/rt_edge/`. Contains no physics of its own. |
| `verify/verify_materials.py`, `verify_dot_levels.py`, `verify_transport.py`, `verify_waveguide.py`, `verify_cw_g2.py` | Per-module regression suites against published numbers the code did not produce. |
| `verify/verify_rt_edge_cards.py`, `verify_rt_edge_contract.py`, `verify_rt_edge_papers.py`, `verify_rt_edge_sweep.py`, `verify_device_rt.py`, `verify_designer_rt.py`, `verify_spec_rt.py` | RT-edge-tier-specific checks: card resolution, contract-heading/anchor validity, evidence-source counting, sweep determinism/schema, `device.py` opt-in-block regression, GUI smoke test, spec-sheet adapter. |

## How to run everything

Legacy + module regression suite (must all pass; from `CLAUDE.md`):

```
python verify/verify_fsim.py && python verify/audit_physics.py && python verify/verify_materials.py && python verify/verify_cw_g2.py && python verify/verify_dot_levels.py && python verify/verify_waveguide.py && python verify/verify_transport.py
```

Current counts observed when this document was written (re-run any of these
to reproduce):

| script | result |
|---|---|
| `python verify/verify_fsim.py` | 51/51 checks passed |
| `python verify/audit_physics.py` | 23/23 quantitative audit items PASS (0 flagged); plus 2 documented notes |
| `python verify/verify_materials.py` | 56/56 materials checks passed |
| `python verify/verify_cw_g2.py` | 17/17 checks passed |
| `python verify/verify_dot_levels.py` | 14/14 checks passed; 8 known deviations (documented, not counted) |
| `python verify/verify_waveguide.py` | 13/13 waveguide checks passed |
| `python verify/verify_transport.py` | 22/22 transport checks passed |

RT-edge-tier-specific checks (not part of the `CLAUDE.md` test command, but
covering the new modules and cards):

| script | result |
|---|---|
| `python verify/verify_rt_edge_cards.py` | 131/131 rt-edge card checks passed |
| `python verify/verify_rt_edge_contract.py` | 76/78 contract checks passed (the 2 failures are the `reischle2008_g2_80K` and `hkust_inp_gaasp_wavelength` evidence gaps above) |
| `python verify/verify_rt_edge_papers.py` | 14/16 rt-edge paper checks passed (same 2 evidence gaps) |
| `python verify/verify_rt_edge_sweep.py` | 43/43 rt-edge sweep checks passed |
| `python verify/verify_device_rt.py` | 206/206 device RT checks passed |
| `python verify/verify_designer_rt.py` | 6/6 designer RT checks passed |
| `python verify/verify_spec_rt.py` | 11/11 spec RT checks passed |

The acceptance sweep (regenerates `out/rt_edge/sweep.csv`, `envelope.png`,
`verdict.md`, `manifest.json`):

```
python scripts/run_rt_edge.py               # full contract grid (54 eligible rows); required for a PASS verdict
python scripts/run_rt_edge.py --quick        # small endpoints-only grid; explicitly incomplete, cannot grant PASS
python scripts/run_rt_edge.py --out-dir out/rt_edge   # (default) where the artifacts land
```

The GUI screenshot self-test (drives `fsim_gui/designer.py` headless and
writes a real PNG; exercised by `verify/verify_designer_rt.py`):

```
python fsim_gui/designer.py --frames 10 --screenshot out/rt_edge/gui-smoke.png
```

The presentation build packages the sweep and figures into a shareable
bundle:

```
python scripts/make_presentation.py
```

## The provenance-tag convention, and the "known deviations" table

Every literature number in the code carries one of four tags in its comment
or docstring, per `CLAUDE.md` and `AGENTS.md`:

- `[V]` — verified directly against the cited paper.
- `[DR]` — derived from published data (e.g. a closed form applied to
  published binary parameters).
- `[E]` — a class estimate or range, used where no measurement exists for
  the exact system.
- `[A]` — an assumption, with no literature anchor at all.

Every derived output inherits the widest (least certain) tag in its input
chain.

`verify/verify_dot_levels.py` adds a second bookkeeping device on top of
this: where the geometry a paper actually specifies (its stated material,
matrix, barrier, height, composition) genuinely misses that paper's own
reported number by more than the model's stated closed-form/finite-
difference accuracy can explain, the comparison is **not** forced to pass by
widening the tolerance window. It is moved to a separate "known deviations"
table, printed by the script itself, that records the paper's number, the
model's number, the size of the miss, and the specific documented model
limitation responsible (see the `fsim_core.dot_levels` module docstring).
Known deviations are counted on their own line — "8 known deviations
(documented, not counted)" in the current run — and never as a pass or a
fail. Example entries from that table (`python verify/verify_dot_levels.py`
output):

- `d-i-a`: In<sub>0.5</sub>Ga<sub>0.5</sub>As/GaAs (Chatzarakis geometry, no
  SSL) — model `E_X` = 1.199 eV vs. published 1.30 eV (miss +101 meV);
  attributed to the omitted piezoelectric field of this (211)B growth class.
- `d-ii-a`/`d-ii-b`: InP/GaInP (Reischle/Pryor geometry) — model `E_X` =
  1.686 eV vs. published 1.815-1.836 eV / ~1.84-1.95 eV (miss +129 to +264
  meV); attributed to combining Pryor's isolated measured VBO with this
  module's own bulk gaps and linear-deformation-potential strain.
- `d-iv-a`: HKUST InP/GaAs<sub>0.65</sub>P<sub>0.35</sub> disk — model `E_X`
  = 1.484 eV (836 nm) vs. published 750-755 nm = 1.643-1.653 eV (Gu et al.,
  2025) (miss +159 to +169 meV); attributed to large-mismatch 3D-island
  relaxation that a coherent single-band disk solver does not capture,
  compounded by the paper's line being an ensemble observable.

## Limitations and next steps

- **Single-band solver misses the HKUST 750 nm line.** The confinement
  solver in `fsim_core/dot_levels.py` is a coherent single-band effective-
  mass model (1D well / 2D disk with BenDaniel-Duke matching); it does not
  capture the 3D-island relaxation or k·p band mixing that a 5%+ lattice-
  mismatched system like InP-on-GaAsP actually undergoes. The HKUST device's
  own reported 750-755 nm line (`d-iv-a`, above) is missed by +159 to +169
  meV as a result.
- **Refractive-index tabulation only covers 650-668 nm.** The design cards'
  emission wavelengths are pinned to the one tabulated phosphide index point
  the materials database carries, not to a continuous dispersion curve —
  this is also why the HKUST wavelength anchor is deliberately kept
  `missing` rather than being reconciled with the card (see "Evidence
  status" above).
- **CW detection is IRF-limited at 300 K.** At 300 K the intrinsic dip is
  fast enough that realistic detector timing jitter (`irf_ps` in the swept
  50-200 ps range) substantially washes out the measured dip: the per-card
  `g2_cw0_raw` minimum/median (0.9784/0.9992 and 0.8907/0.9926) sit far
  above the intrinsic `g2_cw0` values (0.4026/0.9443 and 0.3413/0.905),
  which is why the headline metric is the pulsed intrinsic `g2(0)`, not the
  CW raw value, and why a real measurement of this platform needs a pulsed
  drive rather than DC.
- **Pulsed drive required.** Following from the point above: a CW/DC
  measurement of this platform at 300 K cannot resolve the antibunching
  this model predicts; an actual demonstration needs a pulsed excitation
  scheme, as `scripts/run_rt_edge.py`'s pulsed/CW duty-state pairing already
  assumes.
- **Missing second sources.** The `reischle2008_g2_80K` and
  `hkust_inp_gaasp_wavelength` evidence claims each have zero independent
  verified sources (see above); until a second independent source is found
  and verified for each, `docs/rt_edge_contract.md`'s evidence gate blocks
  PASS regardless of what the numerical sweep finds.
