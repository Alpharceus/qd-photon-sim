# qd-photon-sim

Simulation and inverse-design stack for electrically driven epitaxial-QD
single-photon devices, built on the F-series mathematics (filtered-cascade
identity g²₀ = ε = t_XX/t_X, background law g² = 1−ρ²(1−ε), master ceiling
ρ(T_c)²[1−ε(T_c)] = ½, electrical-separation theorem). NSF NQVL QCAP SLE.

**Provenance discipline:** every parameter carries a tag — [V] measured,
[DR] derived from published data, [E] estimate, [A] assumed — and every
output inherits the widest tag in its input chain (schema-enforced). Ranges
are swept, never averaged. Unmeasured [A] inputs are evaluated as envelopes,
not point predictions.

## Validation record

The model is validated against all three published datasets the program
names (details and honesty ledgers in `notes/`):

| Arm | Dataset | Result |
|---|---|---|
| V-a | Chatzarakis et al., PRApplied 20, 034011 (2023) + supplement | joint over-determined fit (g² + τ(T) + Γ(T)), max resid 0.028; T_c = 249 K |
| V-b | Laferrière et al., Nano Lett. 23, 962 (2023) | ε→1 limit confirmed; 300 K point at the Theorem-0 edge |
| V-c | Reischle et al., APL 97, 143513 (2010) + OE 16, 12771 (2008) | ρ-limited with ε small; the 2008 paper's Eq. (1) is the F2 law |

Named open model residuals: Γ(T) high-T shape (Tier-3 independent-boson
candidate); re-excitation/refilling channel (WP-M2′, three-paper convergence).

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

Branch `rt-edge-emitter` adds a room-temperature (300 K heatsink),
electrically driven, edge-emitting InP-dot tier on top of the F-series core:
new `fsim_core` modules for materials, 300 K linewidth, confinement levels,
p-i-n transport, ridge waveguide, and CW correlation; new opt-in
`device.py` blocks; two design cards; and a reproducible acceptance sweep
(`scripts/run_rt_edge.py`) with its own verify suite.

The honest conditional outcome, from `out/rt_edge/verdict.md`, is:

```
VERDICT: FAIL g2_min=0.3986 g2_median_eligible=0.7019 diag_g2_min=0.3986 diag_g2_median_diagnostic=0.7352 flux_max=1544 flux_margin=1.544 flux_shortfall_deprecated=0.6475 median_pass=false coverage_over_eligible=0.25 eligible_fraction=0.1667 eligible=32/192 flux_floor_excluded=160 evidence=incomplete conditional=true headline_coverage=8/192 cw_raw_coverage=0/192 gamma300_pass_max=8 gamma300_threshold=8-10
```

The stop rule is not met as stated. The model's pulsed intrinsic `g2(0)` at
300 K is controlled by an unmeasured InP-dot exciton linewidth; the sampled
threshold is bracketed at 8--10 meV. The explicit 6.5 meV per-card result
fails for the primary InP/GaAsP card and passes for the fallback InP/(Al)GaInP
card. At 6.5 meV, the primary returns `g2_pulsed=0.6523` and 291.6 photons/s,
while the fallback returns `g2_pulsed=0.4563` and 1441 photons/s. The primary card never clears the flux floor, so the conditional result
rests on the fallback card alone. Evidence remains
incomplete: no second verified source supports Reischle 2008's 80 K electrical
`g2` anchor, and the HKUST 750 nm line is not reproduced by the single-band
solver (about 816 nm for the GaAsP card). The CW dip is detector-limited, so
pulsed drive is required. Like-for-like, Reischle's deconvolved 0.25 +/- 0.05
is compared with the model's best intrinsic corner. See **`docs/rt_edge_tier.md`** for current verify
counts, the 80 K/230 K comparisons, assumptions, disclosures, and council
review history.

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
