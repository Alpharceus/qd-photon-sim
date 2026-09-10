# Nitride geometry and Stark contract

This is the piece-5 contract for the nonpolar a-plane and QW-thickness-fluctuation
headline card pairs (`cards/nitride-nonpolar-{pulse,set}-design.yaml`,
`cards/nitride-qw-fluctuation-{pulse,set}-design.yaml`), and the frozen
specification piece 6 (`nitride-sweep-results`) consumes for its geometry
sweep, `--stark` diagnostic mode, and evidence review. It supersedes ONLY
the new opt-in geometry/orientation/QW/bias behavior added by pieces 1-4;
it does not modify or re-litigate `docs/nitride_cavity_contract.md`, the
c-plane cards it governs (`cards/nitride-cavity-{pulse,set}-design.yaml`),
the Deshpande comparison card, or `verify/data/nitride_cavity_anchors.yaml`,
all of which remain unchanged (see "Files out of scope" in
`.workers/specs/nitride-geometry-cards.md`). No production code is edited
by this piece.

## Card pairs (piece 5)

### Nonpolar headline pair

`cards/nitride-nonpolar-pulse-design.yaml` and `cards/nitride-nonpolar-set-design.yaml`
are each a copy of the corresponding round-1 c-plane headline card
(`cards/nitride-cavity-pulse-design.yaml` / `cards/nitride-cavity-set-design.yaml`)
with exactly two intended fields changed: `nitride.dot.orientation`
(`c_plane` -> `a_plane`) and `nitride.dot.polarization_factor` (unset ->
`0.0`, the explicit no-op restatement of `a_plane`'s own zero factor --
see "Accepted field table" below). `H=3 nm` (`nitride.dot.height_nm`),
`R=10 nm` (`nitride.dot.radius_nm`), `x_in=0.25`, `Q=2000`
(`nitride.cavity.Q`), the 100 ps / 80 MHz drive, `nitride.tau_rad0_ns=1.0 ns`
[A], `drive.b_res=0.1`, `thermal` (mesa 20 um, GaN/DBR/Si stack), density,
aperture, and background are held IDENTICAL to the c-plane parent for
controlled comparison. `nitride.dot.orientation=m_plane` is a separately
NAMED core-grid orientation using the SAME reduced material model as
`a_plane` (both resolve to a zero normal polarization factor,
`fsim_core/nitride_materials.py _ORIENTATION_FACTORS`); this contract does
not invent an anisotropic performance difference between them, and no
headline card is issued for `m_plane` (it appears only as a geometry-sweep
axis value in piece 6).

Setting the polarization field's normal component to zero removes this
reduced scalar model's intrinsic QCSE field; it does NOT rotate the
c-plane strain tensor, the c-plane band-edge partition, or the c-plane
effective masses used inside `fsim_core/nitride_levels.py` for the
in-plane (radial) confinement direction -- see "Nonpolar mass/strain/
valence-ordering caveat" below. A scratch copy of either nonpolar headline
card with `polarization_factor=1.0` on `orientation=a_plane` was tried
(`verify/verify_nitride_geometry_cards.py`'s own scratch test) and
confirmed to make `fsim_core.nitride_levels.levels()` return
`valid=False`, `invalid_reasons=("polarization_factor contradicts the
selected orientation",)` -- `fsim_core.nitride_materials.orientation_factor`
raises `ValueError` for exactly this contradiction, caught inside
`nitride_levels._validate` and surfaced as an invalid row, never a
silently-accepted, physically-wrong contradiction and never a hard crash
on card load.

### QW-fluctuation headline pair

`cards/nitride-qw-fluctuation-pulse-design.yaml` and
`cards/nitride-qw-fluctuation-set-design.yaml` are each a copy of the same
c-plane parent with `nitride.dot.geometry_type=qw_fluctuation`,
`nitride.dot.shape=disc` (the only shape this geometry type supports),
`H=3.5 nm`, `R=20 nm`, `w=3 nm` (`nitride.dot.wl_thickness_nm`, the local
thickness-FLUCTUATION width), `x_in=0.25` [A design choices], and the
matching diode fields `drive.diode.d_i_nm=27`, `drive.diode.d_active_nm=3.5`,
`drive.diode.wl_thickness_nm=3.0` (EQUAL to `nitride.dot.wl_thickness_nm`),
`drive.diode.x_in=0.25` (EQUAL to `nitride.dot.x_in`) -- required exactly
by `fsim_core/device.py _evaluate_nitride`'s own QW-card consistency
guard ("QW cards require drive.diode.%s == nitride.dot.%s" for
`wl_thickness_nm`/`x_in`, and `drive.diode.d_i_nm >= dot height`).
`ret.channel=min` is retained (the RetentionBlock default `pair_half` is
still rejected -- see `docs/nitride_cavity_contract.md`). No explicit
`nitride.background.reservoir_energy_eV` override is set: `_evaluate_nitride`
rejects any such override on a `qw_fluctuation` card by construction, so
this card class always reads its reservoir from the confinement solver's
own matched well (see "Reservoir and SRH convention" below). Same-
composition thickness fluctuation confines carriers LATERALLY through the
subband-energy difference between the thin dot column (`height_nm=3.5`)
and the thicker surrounding well (`wl_thickness_nm=3.0`) --
`fsim_core/nitride_levels.py _qw_levels`'s own `de_lat`/`dh_lat` --
never through a composition or barrier-height change.

Within each pair (nonpolar and QW-fluctuation separately), geometry,
current, pulse duration (100 ps), repetition rate (80 MHz), thermal model,
linewidth, cavity, aperture, and background are held IDENTICAL between the
pulse and SET cards; only names, `drive.cycle_loading`, and SET-specific
metadata (`drive.eta_load`, `drive.set_params`) differ. Headline
`thermal.T_hs` default is 300 K with the sweep covering 230/250/273/300 K;
`nitride.dot.screening_fraction` default is 0, paired with 0.5/1 in the
sweep. Each SET card's `drive.set_params.radius_nm` (the classical
charging-energy island) is independent of `nitride.dot.radius_nm` (the
confined optical dot) and is NEVER resized by a geometry sweep row -- a
sweep over `nitride.dot.radius_nm` holds `drive.set_params.radius_nm` fixed
at each card's own declared value.

## Accepted field table (pieces 1-4)

Fields below are accepted by `fsim_core.nitride_levels.NitrideDotSystem`
(pieces 1-2, commit 57ab49d and later) and `fsim_core.device._evaluate_nitride`
(piece 4, commit 9d7ebd3 and later). All are [A] design-choice defaults
unless separately cited on a card's own `provenance.sources` entry.

| Field (`nitride.dot.*`) | Type / accepted values | Default | Notes |
| --- | --- | --- | --- |
| `height_nm`, `radius_nm`, `x_in` | positive float / float in [0,1] | required | core geometry |
| `wl_thickness_nm` | float, `0 <= wl < height_nm` for `qw_fluctuation`; must be `0` for `isolated_dot` | `0.0` | QW fluctuation width `w` |
| `strain_fraction`, `screening_fraction` | float in [0,1] | `1.0`, `0.0` | |
| `external_field_kVcm` | finite float | `0.0` | added to the diode bias field at evaluation time |
| `vbo_InN_GaN_eV`, `strain_c_fraction` | finite float | `1.15`, `0.7` | |
| `orientation` | `c_plane` \| `semipolar_11_22` \| `m_plane` \| `a_plane` | `c_plane` | `fsim_core.nitride_materials._ORIENTATION_FACTORS` |
| `polarization_factor` | `None`, or float exactly matching the orientation's fixed factor (`c_plane`=1, `m_plane`/`a_plane`=0), or any float in [0,1] for `semipolar_11_22` only | `None` (resolves to the orientation's own factor) | a mismatched explicit value raises `ValueError` |
| `shape` | `disc` \| `lens` \| `truncated_cone` | `disc` | `qw_fluctuation` requires `disc` |
| `top_radius_fraction` | float in [0,1], `truncated_cone` only | `1.0` | |
| `shape_height_fraction` | `None`, or float in `[native, 1]` (native = 0.5 lens, `(1+t+t^2)/3` cone, unsupported for `disc`) | `None` | effective-height override |
| `geometry_type` | `isolated_dot` \| `qw_fluctuation` | `isolated_dot` | |

`fsim_core.nitride_transport.NitrideDiode` (`drive.diode.*`, piece 4 wiring)
additionally accepts `d_active_nm` (default 3.0 nm) and `wl_thickness_nm`
(default 3.0 nm) beyond the fields already documented on the c-plane cards'
own `provenance.diode_unset_defaults`; QW-fluctuation cards are the first
to override both explicitly.

`nitride.bias.*` (piece 3/4, `fsim_core.nitride_stark.resolve_bias` via
`_evaluate_nitride`): `mode` (`current` default, or `junction_voltage`),
`field_polarity` (`+1`/`-1`, default `+1`), `V_j_V`/`T_j_K` (required
together, `junction_voltage` mode only), `cavity_reference_V_j_V`
(optional, any mode). None of the four cards in this piece override
`nitride.bias` -- all four run in the default `current` bias mode,
self-consistent `T_j`, `field_polarity=+1`. The `--stark` mode (below)
is the piece-6 consumer of `junction_voltage` mode.

## Fixed-cavity-reference Stark convention

`nitride.cavity.T_track` anchors the SAME dot, at the SAME supplied
current (or, in `--stark` mode, a declared reference bias), evaluated at
`T_track` instead of the row's own `T_j` -- a single fixed design anchor
(`_evaluate_nitride`'s `track0`), not recomputed per operating point. On
every card in this piece, `nitride.cavity.T_track=300.0 K`, matching
`thermal.T_hs`'s own headline default, so the cavity/dot are on-resonance
at the headline operating point by construction; sweeping `T_hs` away from
300 K, or `V_j` in `--stark` mode, moves the dot line away from this fixed
cavity anchor and the reported `detuning_meV`/`Fp_add`/`F_eff_*` track that
real detuning. Source validity (whether the underlying literature number
was retrieved and transcribed -- `evidence_status` in the ledger) and
spectral validity (whether a row's dot/cavity/counting solve converged to
a physically bound, finite state -- `scalars['valid']`/`spectroscopy_valid`)
are SEPARATE axes; a card can be spectrally valid while resting on entirely
[A] inputs, and a literature number can be fully [V] while being [A]/[E]
to TRANSFER onto this model (see the ledger's own schema header). No
current-derived screening law is inferred anywhere in this piece or in
piece 6: `nitride.dot.screening_fraction` is always an explicit scenario
input (0/0.5/1), never fit to a measured current or bias trace.

## `--stark` output and compatibility table (frozen for piece 6)

`scripts/run_nitride_geometry_stark.py` (piece 6, out of scope for this
piece) MUST implement `--stark` as: for each of the c-plane and nonpolar
headline cards, at `nitride.dot.height_nm` in `[1, 2, 3, 5]` nm (other
geometry fields held at the card's own reference), sweep `nitride.bias.mode
=junction_voltage` over a declared `V_j_V` grid at fixed `T_j_K`, and
separately sweep `nitride.bias.mode=current` (self-consistent `T_j`,
labelled `temperature_mode='self_consistent'` per `scalars`) over a current
grid, for `nitride.dot.screening_fraction` in `[0, 0.5, 1]`. Required
per-row outputs (all already exposed by `_evaluate_nitride`'s `scalars`):
`E_X_eV`, `field_kVcm`, `overlap_sq`, `tau_rad_bare_ns`,
`tau_rad_cavity_ns`, `valid`/`invalid_reasons`, `evaluation_kind`
(`stark_diagnostic` for `junction_voltage` mode). `dE_X/dV` in meV/V is the
finite-difference slope of `E_X_eV` against `V_j_V` along a validity-masked
segment of the trace (never differenced across an invalid row). Figures
plot `E_X(V)` and `tau_rad(V)` (bare and cavity) with the three screening
curves overlaid and `zhang2016_stark_slope`'s -10 meV/V marked as a
NON-GATING comparison line (ledger `evidence_kind=non_gating_comparison`
-- never a fit target, never relabelled [V] for this model's transfer).
The compatibility table maps a user-SUPPLIED measured-slope interval
(never a fabricated one, and never the illustrative or Zhang interval
substituted for a user's own) to the subset of `screening_fraction` values
whose predicted `dE_X/dV` interval overlaps it, per geometry and
orientation; a tie across all three screening fractions is reported
`compatible=None` / `screening_unidentifiable`, matching piece 3's own
`nitride_stark` tie-break convention (never silently picking one).

## Frozen geometry grid (piece 6 core sweep)

Heights (nm): `1, 2, 3, 4, 5, 7, 10` (extends the c-plane cards' own
`[1,2,3,4,5]` axis). Radii (nm): `5, 10, 15, 20, 30` (extends `[5,10,15]`).
Orientations: `c_plane`, `semipolar_11_22`, `m_plane`, `a_plane` (all four
names exactly as `fsim_core.nitride_materials._ORIENTATION_FACTORS` keys).
`thermal.T_hs` (K): `230, 250, 273, 300`. `nitride.dot.screening_fraction`:
`0, 0.5, 1`. Both electrical regimes (`rectangular`, `deterministic_pair`).
Core grid size: `7 heights x 5 radii x 4 orientations x 4 T_hs x 3
screening x 2 regimes = 3,360 rows`. Piece 6's total `evaluate()` call
budget across this core grid plus its reduced shape/QW/Stark/sensitivity
sets MUST stay `<= 10,000` calls, and the full run must finish under 30
minutes on a laptop (today's rate is approximately 0.09 s/call, so 10,000
calls is approximately 15 minutes); `--quick` must finish inside a
10-minute test command.

### Auxiliary Wang 2017 nonpolar comparison row

One auxiliary sweep row (NOT part of the 3,360-row core grid, added on top
of it) uses `nitride.dot.height_nm=7.0`, `nitride.dot.radius_nm=17.5`
(half the ledger's `wang2017_geometry` 35 nm AFM diameter),
`orientation=a_plane`, at the nonpolar cards' own `x_in=0.25` -- composition
is NOT transcribed from Wang et al. 2017 (the digest does not report the
dots' In fraction) and this row's `x_in` is therefore an [E]/[A]
composition-transfer disclosure, not a measured value. This row exists to
place this project's own model output alongside the `wang2017_*` ledger
anchors' optical, 220 K, uncapped-AFM measurements for VISUAL comparison
only; it is never treated as a validation of the electrical, 230-300 K,
capped-dot model against those numbers (see `wang2017_geometry`'s own
`transfer_notes`).

## Reservoir and SRH convention (per card)

`fsim_core.device._evaluate_nitride`'s `_nitride_reservoir_energy_eV`
resolves the reservoir spectral edge from `nitride.dot.wl_thickness_nm`
(zero -> GaN barrier Varshni edge; nonzero -> the strained InGaN material
continuum), independently of `drive.diode.wl_thickness_nm` (the SRH
background carrier model's own reservoir thickness) unless
`nitride.background.reservoir_energy_eV` is set explicitly. This
independence is a PRE-EXISTING, documented mismatch for `isolated_dot`
cards; this contract does not resolve it there, only for the new
`qw_fluctuation` class:

| Card | `nitride.dot.wl_thickness_nm` | `drive.diode.wl_thickness_nm` | Reservoir/SRH relationship |
| --- | --- | --- | --- |
| `cards/nitride-cavity-pulse-design.yaml` | 0.0 (GaN barrier) | 0.5 nm | historical mismatch, unchanged, documented in `docs/nitride_cavity_contract.md` |
| `cards/nitride-cavity-set-design.yaml` | 0.0 (GaN barrier) | 0.5 nm | historical mismatch, unchanged |
| `cards/nitride-deshpande2014-comparison-design.yaml` | 0.0 (GaN barrier) | 0.5 nm | historical mismatch, unchanged |
| `cards/nitride-nonpolar-pulse-design.yaml` | 0.0 (GaN barrier) | 0.5 nm | historical mismatch, unchanged (orientation-only card) |
| `cards/nitride-nonpolar-set-design.yaml` | 0.0 (GaN barrier) | 0.5 nm | historical mismatch, unchanged |
| `cards/nitride-qw-fluctuation-pulse-design.yaml` | 3.0 nm (matched InGaN well) | 3.0 nm (EQUAL, enforced by `_evaluate_nitride`'s QW guard) | MATCHED: confinement dot and diode SRH background share the same surrounding-well material and thickness |
| `cards/nitride-qw-fluctuation-set-design.yaml` | 3.0 nm (matched InGaN well) | 3.0 nm (EQUAL, enforced) | MATCHED, same as pulse |

For the QW-fluctuation cards, `fsim_core.nitride_levels._qw_levels` reports
`reservoir_kind` as `'ingan_qw'` when BOTH carriers' own escape channel
selects the surrounding InGaN well below the GaN cladding edge, and
`'gan_barrier'` when either carrier's true escape channel is the wider GaN
barrier plateau instead -- a physically-motivated PER-CARRIER selection
computed by the confinement solver itself (`fsim_core/nitride_levels.py`,
under active revision by a concurrent piece-2 fix at the time this
contract was written; see this piece's own STATUS `decisions` for the
specific `reservoir_kind` this pair's spec-mandated `H=3.5/R=20/w=3 nm`
geometry evaluates to today). Either value is a legitimate, physically
determined QW-branch output; `reservoir_kind` is still always sourced from
the `qw_fluctuation` code path (never inherited from the `isolated_dot`
default) whenever `nitride.dot.geometry_type=qw_fluctuation`, which is
what this piece's verifier checks. QW cards still use the EXISTING
bulk-SRH carrier approximation (`fsim_core.nitride_transport.
evaluate_injection`) for the diode's own background rate, rather than a
newly solved QW carrier distribution -- matched reservoir THICKNESS and
composition does not imply a re-derived carrier transport model.

## Nonpolar mass/strain/valence-ordering caveat

`fsim_core.nitride_materials.orientation_factor`'s own docstring: the
c-plane scalar strain, band-edge partition, and c-plane masses are
retained for every orientation [A/E transfer]. For nonpolar `m_plane`/
`a_plane` growth the confinement axis is perpendicular to c; whether
`fsim_core/nitride_levels.py` swaps the growth-axis and lateral effective
masses for this orientation is owned by piece 2 (`fsim_core/
nitride_levels.py`, out of scope for this piece -- see that module's own
docstring for the current convention and magnitude estimate). Nonpolar
valence-band ordering is not modelled at all [A; Schade et al., phys.
status solidi (b) (2011), see `verify/data/nitride_geometry_stark_anchors.
yaml schade2011_orientation`, which is explicitly NOT a numeric
polarization-factor measurement]. No accuracy is claimed for this c-plane-
parameter transfer to nonpolar growth, and none is claimed for the
volume-mapped shape approximations (`disc`/`lens`/`truncated_cone`,
`fsim_core/nitride_levels.py _geometry`'s own docstring: paraboloid vs.
spherical-cap lens profile choice is part of an UNQUANTIFIED shape
systematic error). Both are reported as unknown systematic error, not
bounded; the planned mitigation is an alternative-height sensitivity sweep
(`nitride.dot.height_nm` one-at-a-time axis, already present on every
card's own `provenance.ranges`), not a corrected mass/strain model.

## Idealized optical vs. hardware outcomes

As on the c-plane cards, the pulse-regime cards report only the idealized
optical counting result (`g2_op`, `collected_flux_pulsed_s`); the SET-regime
cards report the idealized deterministic one-pair-per-cycle counting result
(`fsim_core.pulse_counting.deterministic_cycle_g2`) SEPARATELY from the SET
hardware feasibility screen (`set_feasible`, `fsim_core.drive_mech.
set_feasibility`) -- an infeasible SET screen never overwrites the
idealized counting result, and `scalars['device_pass']` requires both to
pass together for `deterministic_pair` cards. This separation is
unchanged by orientation or geometry_type.

## Evidence ledger schema

`verify/data/nitride_geometry_stark_anchors.yaml` is a SEPARATE ledger from
`verify/data/nitride_cavity_anchors.yaml` (unmodified). Top-level `anchors`
is a mapping keyed by anchor id; each entry carries `citation`,
`doi_or_url`, `location`, `evidence_status` (`full_text` \| `abstract_only`
\| `figure_reading` \| `missing`, same enum as the existing ledger),
`evidence_kind` (`source_transcription` \| `non_gating_comparison` \|
`missing_evidence`, NEW field distinguishing verification role from
retrieval completeness), `tag` (`V`/`DR`/`E`/`A`, the measurement's own
tag, distinct from any card-field TRANSFER tag), `platform`, `excitation`,
`observable`, `value`, `unit`, `tolerance`, `transfer_notes`. `doi_or_url`
uses only identifiers actually printed in `../_goal/nitride_digests.md`
(never a fabricated DOI); the ledger includes the Wang et al., APL 111,
053101 (2017) fine-structure-splitting entry as a `non_gating_comparison`
omitted-physics caution, and Schade et al. (2011) as `missing_evidence`
(no specific number is transcribed from it anywhere in this project).
Schade's paper is never treated as a numeric polarization-factor
measurement.

## Round trip and verifier scope

`verify/verify_nitride_geometry_cards.py` loads every new card through
`fsim_core.device.DeviceDesign.load` (not raw YAML parsing) and evaluates
each at 230 K and 300 K in both regimes with explicit `screening_fraction`
in `{0, 0.5, 1}`, checking finite scalars on valid rows and retained
`invalid_reasons` on invalid ones (a schema error, e.g. a card that fails
to `load()`, is never relabelled as an expected physics invalidity). It
diffs each new card's loaded `design` dict against its c-plane parent,
restricted to the physics-relevant top-level blocks (excluding
`provenance`, which is documentation), against an explicit per-pair
allow-list of intended-to-differ keys, to prove the orientation/geometry
change is isolated and that density, aperture, background, lifetime, and
SET hardware parameters have not been opportunistically improved. It
checks nonpolar screening-degeneracy at fixed bias, QW reservoir wiring,
the polarization-factor contradiction scratch test, ledger transcription
independence (Wang geometry/lifetime/linewidth/g2, Zhang slope, with
excitation/temperature/raw-vs-corrected distinctions enforced), and that
`verify/verify_nitride_cards.py`'s existing 377 checks still pass
unmodified. It round-trips each new card's design dict through an
in-memory YAML dump/reload (never overwriting the on-disk card fixtures)
and confirms the nitride geometry/bias/provenance blocks survive intact.
