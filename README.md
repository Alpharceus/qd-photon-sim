# qd-photon-sim

Simulation and inverse-design stack for electrically driven epitaxial-QD
single-photon devices, built on the F-series mathematics (filtered-cascade
identity — the mu→0 limit of the cap-2 filtered cascade — g²₀ = ε = t_XX/t_X;
background law g² = 1−ρ²(1−ε); master ceiling
ρ(T_c)²[1−ε(T_c)] = ½; electrical-separation theorem). NSF NQVL QCAP SLE.

**Provenance discipline:** every parameter carries a tag — [V] measured,
[DR] derived from published data, [E] estimate, [A] assumed — and every
output inherits the widest tag in its input chain (schema-enforced). Ranges
are swept, never averaged. Unmeasured [A] inputs are evaluated as envelopes,
not point predictions.

## Validation record

Every check in this repository sorts into exactly one of five classes;
conflating them (treating a calibration as if it were a prediction, or a
model-vs-data comparison as if it were an independent transcription check)
is the mistake this section exists to prevent.

- **(N) numerical verification** — one method checked against another under
  the same physical assumptions, with no external data involved: `gate_d1.py`,
  `gate_d2.py`, `gate_d3.py`, `gate_spec.py`, `gate_v11_gui.py` (internal
  wiring, round-trip, and GUI-headless-boot consistency checks);
  `verify_cw_g2.py` (also carries one literature anchor, the (g3) Reischle
  2008 CW-dip check, class T); `verify_dbr.py`; `verify_designer_rt.py`;
  `verify_device_rt.py` (also carries a handful of literature anchors — the
  Bommer et al. 2011 retention activation energy, the Reischle 2008 80 K
  electrical rho anchor, and a Schubert et al. 1995 refractive-index check —
  class T); `verify_dot_levels.py` (also carries one literature-anchored
  acceptance check, the Bommer-class hole-escape window, plus several
  non-scored "known deviation" literature comparisons that are reported but
  never counted pass/fail, class T); `verify_drive.py`; `verify_fsim.py`;
  `verify_gf.py`; `verify_presentation.py`; `verify_presentation2.py`;
  `verify_rt_edge_cards.py` (checks that the design cards round-trip and
  that their citation strings and mode opt-ins match the contract/ledger —
  mechanical, not a re-check of the cited numbers themselves);
  `verify_rt_edge_sweep.py`; `verify_sde.py`; `verify_spec_rt.py`;
  `verify_waveguide.py` (also carries published-class range comparisons for
  the HKUST ridge stack, class T); `verify4.py`; and the `mc_*` Monte-Carlo
  second methods used throughout `verify/`.
- **(T) source transcription** — a number, formula, or claim in the code
  checked against the paper that states it, or a check of the evidence
  ledger's own structure: `audit_physics.py`; `verify_materials.py`;
  `verify_rt_edge_contract.py` (a ledger-structure validator — schema,
  cross-references, required section headings — not a check of the paper
  values themselves); `verify_rt_edge_papers.py`; `verify_transport.py`;
  and the anchors ledger (`verify/data/rt_edge_anchors.yaml`).
- **(C) parameter calibration** — free parameters fit to data and then
  evaluated on that same data. Not held-out validation.
- **(M) model-vs-data comparison, unfitted** — a class-range or previously
  fixed model checked against published data it was not fit to, without
  adjusting any parameter to match it. Still not held-out validation: the
  comparison target was known while the class ranges were chosen: V-b (the
  Laferrière envelope) and V-c below.
- **(P) held-out prediction** — a fitted or class-range model checked
  against data it did not see during fitting or class-range selection.
  **Currently empty**: no result below has been checked against withheld
  data.

The model is calibrated on one published dataset (V-a) and compared, not
fitted, against two others (V-b, V-c); no held-out prediction exists yet
(details and honesty ledgers in `notes/`):

| Arm | Class | Dataset | Result |
|---|---|---|---|
| V-a | (C) | Chatzarakis et al., PRApplied 20, 034011 (2023) + supplement | joint over-determined fit (g² + τ(T) + Γ(T)), max resid 0.028; T_c = 249 K — calibration, not held out |
| V-b | (M) | Laferrière et al., Nano Lett. 23, 962 (2023) | ε→1 limit confirmed; 300 K point at the Theorem-0 edge |
| V-c | (M) | Reischle et al., APL 97, 143513 (2010) + OE 16, 12771 (2008) | ρ-limited with ε small; the 2008 paper's Eq. (1) is the F2 law |

Named open model residuals: Γ(T) high-T shape (Tier-3 independent-boson
candidate); re-excitation/refilling channel (WP-M2′, three-paper convergence).

The requirement that an anchor claim rest on two distinct verified sources
(the "Evidence status" section of `docs/rt_edge_contract.md`) is a project
policy adopted for this repository's evidence ledger, not a general
scientific requirement: a single primary publication can be sufficient
evidence for a claim, and two publications that are both off-platform
(neither on the actual InP-dot device/material system) are not automatically
sufficient just because there are two of them.

## Scope and data provenance

- **The validation data are digitizations.** The V-a/V-b/V-c comparisons run
  against my digitizations of published figures plus values printed in the
  papers, not author-released datasets. The source PDFs are deliberately not
  distributed with this repository, so reproducing the fits from scratch
  means obtaining the papers and re-extracting the data.
- **Inverse design outputs necessary conditions, not geometry.** `spec.py`
  turns a g² target into required linewidth, gain, mode volume, and
  background budget. It does not produce a cavity geometry; the EM-solver
  step that would consume these requirements is planned, not present.
- **Model scope.** Exciton–biexciton cascade with Lorentzian lineshapes.
  No carrier transport, no growth modelling, no indistinguishability
  metrics.

## Layout (three-layer rule: core is headless; GUIs compute no physics)

- `fsim_core/` — physics + assembly:
  - `card.py` — parameter-card schema with [V/DR/E/A] tags
  - `spectral.py` (Module C) — F1/F1a closed forms, transmissions, MC second methods
  - `loading.py` (Module D) — cap-2/F1b drive statistics, injection background, F5 aperture lemma
  - `integrator.py` (Module E) — background law, master-ceiling T_c, sensitivities
  - `thermal.py` (Module A) — spreading-resistance stack, T_j, runaway detection
  - `cavity.py` (Module B) — tracking rule, F_eff, collection gain, SiN β (Lemma 1)
  - `fitting.py` — joint validation-fit driver
  - `device.py` — DeviceDesign blocks + `evaluate()` / `evaluate_envelope()`
  - `spec.py` — **inverse design**: required κ, G, mode volume Ṽ, b_e budget, density, mesa from a g² target
  - `presets.py`, `design_meta.py` — device presets; per-parameter provenance metadata
- `fsim_viz/` — matplotlib figure factory; every figure ships its CSV (Origin-ready)
- `fsim_gui/`
  - `designer.py` — **device designer** (Dear PyGui): block-diagram chain, fab-stack editor with live cross-section, envelope mode ([A] inputs default-ranged → bands + tornado), A/B/C comparison slots, one-click report bundles
  - `app.py` — validation dashboard (Streamlit; frozen at its four panels)
- `cards/` — parameter cards (validation + prediction) and `*-design.yaml` device designs
- `scripts/` — `run_phase0..3.py` (validation fits, thermal maps, V-b/V-c, prediction envelopes + Osinski packet), `run_spec.py` (design-target spec sheets)
- `verify/` — regression suite (50 checks incl. MC/FD second methods), physics audit vs published values (23 items), phase gates (`gate_d1..d3.py`, `gate_spec.py`)
- `notes/` — per-phase results notes, physics audit, Osinski design-review packet
- `out/` — committed report bundles (figures + CSVs), regenerable from the scripts

## Run

```
pip install -r requirements.txt

python verify/verify_fsim.py      # regression suite (50 checks), exit 0 iff green
python verify/audit_physics.py    # known-parameter physics audit (23 items)

python fsim_gui/designer.py       # device designer (configure -> RUN -> graphs/numbers)
streamlit run fsim_gui/app.py     # validation dashboard

python scripts/run_phase0.py      # V-a joint fit          -> out/phase0
python scripts/run_phase1.py      # thermal maps + V-c     -> out/phase1
python scripts/run_phase2.py      # V-b + cavity design    -> out/phase2
python scripts/run_phase3.py      # requirement envelopes  -> out/phase3
python scripts/run_spec.py        # inverse-design spec    -> out/spec
```

Python ≥ 3.9. After editing cards or code, re-run `verify/` — every phase
gate and review finding is encoded as a permanent check.

## RT edge-emitter tier (branch rt-edge-emitter)

Branch `rt-edge-emitter` adds a non-cryogenic (230--300 K heatsink),
electrically driven, edge-emitting InP-dot tier on top of the F-series core:
new `fsim_core` modules for materials, 300 K linewidth, confinement levels,
p-i-n transport, ridge waveguide, and CW correlation; new opt-in
`device.py` blocks; two design cards; and a reproducible acceptance sweep
(`scripts/run_rt_edge.py`) with its own verify suite.

The honest outcome, from `out/rt_edge/verdict.md`, is:

```
VERDICT: FAIL model=finite_pulse:true,tau_cap_density:false g2_min=nan g2_median_eligible=nan diag_g2_min=0.9817 diag_g2_median_diagnostic=0.9993 flux_max=582.8 flux_margin=0.5828 flux_shortfall_deprecated=1.716 median_pass=false coverage_over_eligible=nan eligible_fraction=0 eligible=0/768 flux_floor_excluded=768 evidence=incomplete conditional=false headline_coverage=0/768 headline_coverage_pulsed=0/384 headline_dedup_mismatch_groups=0 eligible_dedup=0/384 eligible_dedup_mismatch_groups=0 rows_scheduled=0/768 cw_raw_coverage=0/0 gamma300_pass_max=nan gamma300_threshold=n/a T_pass_min=none headline_by_T=230:0/192,250:0/192,273:0/192,300:0/192 headline_by_T_pulsed=230:0/96,250:0/96,273:0/96,300:0/96
```

This is the pkg5b-corrected headline model (`drive.finite_pulse=true`, the
real finite-pulse-waveform loading calculation, replacing the legacy static
per-pulse `mu` approximation every earlier package used): under it, **no
sampled corner clears the 1 kHz collected-flux floor at all**
(`eligible=0/768`), so `g2_min`/`T_pass_min` are undefined (`nan`/`none`) and
the previous conditional result no longer holds. The old, uncorrected
`finite_pulse=false, tau_cap_density=false` combination (every package before
this one) now also reports `0/768` eligible on the current physics; only the
opt-in `ret.tau_cap_scales_with_density=true` ("full-cancellation") retention
convention restores any eligible/passing rows when `finite_pulse=false`
(392/768 eligible, 170 headline passes, `g2_min=0.2995`) -- but that
combination never runs with the corrected loading model. With
`finite_pulse=true` (the corrected model) turned on, `tau_cap_density=true`
gives 496/768 eligible rows and still **0 headline passes**
(`g2_min=0.8291`, above the 0.5 gate): **nothing passes g2 < 0.5 once
finite-pulse counting is on**, in either retention convention. The
mechanism is re-excitation inside the 100 ps pulsed-counting gate (a second
capture-emission cycle before the gate closes), which the static
per-pulse loading approximation could not see. See the verdict's "Model
sensitivity" table for all four combinations side by side. See
**`docs/rt_edge_tier.md`** for current counts, assumptions, disclosures, and
council review history.

## Status

All planned phases delivered (0, V, 1, 2, 3 + designer D0–D3 + spec mode).
Headline spec results (tag [A], class-proxy Γ/retention): the staged
77 K / g²≤0.1 target closes on every route; 120 K / g²≤0.1 requires
Δ_XX ≳ 5 meV and either b_e ≤ 0.002 (slit) or a tracked cavity with G ≥ 24,
κ ≈ 1.03–1.06 meV, Ṽ ≤ 14 (λ/n)³; the 300 K route needs G ≥ 5–53,
κ ≤ 3–10 meV, Ṽ ≤ 4–6 (λ/n)³. The decisive in-house measurements, in
order: Δ_XX distribution, injection background b_e(I,T), Γ(T).

Source papers (PDFs) and program planning documents live outside this
repository and are not distributed with it.
