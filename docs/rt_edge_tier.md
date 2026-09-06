# RT edge-emitter tier

This tier evaluates non-cryogenic heatsink operation (`T_hs >= 230 K`), electrically driven, edge-emitting InP-dot single-photon operation with the contract grid in `docs/rt_edge_contract.md`. It is a model-and-evidence assessment, not a device demonstration. The current results below are read from the generated acceptance artifacts and the current verification scripts.

## Outcome

Under the relaxed stop rule (`T_hs >= 230 K`), the sweep passes its headline g2 and flux gates at 230, 250, and 273 K for the primary InP/GaAsP card, and at 230, 250, 273, and 300 K for the fallback InP/(Al)GaInP card. These headline rows occur at the sampled `gamma300=6.0 meV` endpoint; the primary has no headline pass at 300 K, while the fallback retains 12/192 headline rows there. The pooled linewidth threshold is bracketed at `T_pass_min=230` and `gamma300_threshold=10-12` meV. The 300 K line remains the room-temperature result: it is conditional on the fallback card and on an unmeasured InP-dot linewidth, and evidence is still incomplete, so the verdict remains conditional rather than a PASS.

The current `out/rt_edge/verdict.md` line is quoted verbatim:

```
VERDICT: FAIL g2_min=0.3214 g2_median_eligible=0.6373 diag_g2_min=0.3214 diag_g2_median_diagnostic=0.6706 flux_max=6890 flux_margin=6.89 flux_shortfall_deprecated=0.1451 median_pass=false coverage_over_eligible=0.3066 eligible_fraction=0.5521 eligible=424/768 flux_floor_excluded=344 evidence=incomplete conditional=true headline_coverage=130/768 cw_raw_coverage=0/424 gamma300_pass_max=10 gamma300_threshold=10-12 T_pass_min=230 headline_by_T=230:60/192,250:34/192,273:24/192,300:12/192
```

The fallback card supplies the best eligible corner at `delta_xx=8 meV`, `gamma300=6 meV`, and 50 ps IRF: pulsed `g2(0)=0.3214`, collected flux 6890 photons/s, and 130 headline-passing rows of 768. Per temperature, headline rows are 60/192, 34/192, 24/192, and 12/192 at 230, 250, 273, and 300 K respectively. The CW raw gate remains zero. The primary card contributes passes below 300 K, but the 300 K conditional result rests on the fallback card.

Evidence is still incomplete for two central claims: there is no second verified source for Reischle 2008's 80 K electrical `g2`, and the HKUST 750 nm line is not reproduced by the single-band solver (the GaAsP card computes about 816 nm). The evidence ledger also records the transferred 80 K residual-background claim as single-sourced. These gaps independently block PASS.

## Comparison points and interpretation

The model's relevant comparison values must be interpreted as predictions, not new measurements. Like-for-like against Reischle et al. (2008), the IRF-deconvolved, background-included anchor is `g2(0)=0.25 +/- 0.05`, versus the model's best intrinsic corner `g2_min=0.3214`; the raw and fully corrected values are not the matching convention. Against Chatzarakis et al. (2023), the 230 K comparison is `g2(0)=0.36`; the model's 230 K comparison prediction is 0.36. These numbers and the anchor provenance are checked by `verify/verify_rt_edge_papers.py` and `verify/data/rt_edge_anchors.yaml`.

The key practical result is therefore narrow: non-cryogenic operation first passes at 230 K (`T_pass_min=230`), with the pooled continuous linewidth threshold bracketed between 10 and 12 meV by the sampled grid. The 6.5 meV per-card anchor still passes only on the fallback card; a second independent source for the 80 K electrical anchor is also needed before the evidence gate can pass.

## NA-collection correction

The waveguide review found that the old collection estimate used the wrong separable formula inside a circular NA cone. The corrected model integrates the circular cone and uses the per-axis Gaussian factor `erf(sqrt(2) theta_max/theta_div)`, with `theta_div=lambda/(pi w0)`; the review also standardizes the 1/e² width convention. At the favourable diagnostic corner, the old result was NA collection `0.548244` and 1544.4 photons/s; after the correction these are `0.857853` and 6890.03 photons/s (about 4.46x higher). The updated flux decomposition still identifies confinement retention `S` as the limiting factor. This is a model correction, not an experimental brightness measurement.

## Assumptions carried by the cards

The design cards explicitly identify the inputs carrying `[A]` or `[E]` provenance that control this conclusion:

| input | card status |
|---|---|
| `dot.gamma300` | `[E]` class range; exact InP platform value unmeasured |
| `drive.b_res` | `[E]` 80 K residual-background transfer |
| `emission.NA` | `[A]` collection lever |
| `emission.R_back` | `[A]` HR-back-facet lever |
| `emission.L_um` | `[A]` ridge-length lever |
| `aperture.density_cm2` | `[A]` single-dot-isolation target |
| `emission.lambda_nm` | `[DR]` solver-derived wavelength, with an `[E]` index fallback |

The passing fallback configuration declares `NA=0.75`, `R_back=0.95`, and `L_um=250`; its card also declares a density of `3e8 cm^-2`. The flux floor is 1000 photons/s and the maximum flux margin is 1.544. The brightness factor table in `out/rt_edge/verdict.md` attributes the 1544 photons/s maximum to the reported decomposition. The facet model uses the fused factor `T / (T + (1 - R_back))`, independently matching the back-solved factor at 0.935012.

## Required disclosures

- The filter window automatically matches the exciton linewidth, so `t_X=0.5` for every sweep row. This convention makes reported flux insensitive to `gamma300` by construction; it is not a physical claim that real-device brightness is linewidth-independent.
- `drive.b_res` transfers an approximately 12% residual background from an 80 K Reischle-class data point to a different 300 K platform.
- Pulsed drive is required for a demonstrable result: the predicted CW dip is narrower than the 50--200 ps detector response and convolution raises raw CW `g2(0)` above the threshold.

## Reproduction and current verification counts

Run the acceptance sweep with `python scripts/run_rt_edge.py`; it regenerates `out/rt_edge/sweep.csv` and `out/rt_edge/verdict.md`. The current script outputs are:

| script | current result |
|---|---|
| `verify/verify_fsim.py` | 51/51 checks passed |
| `verify/audit_physics.py` | 23/23 quantitative audit items PASS (0 flagged); 2 documented notes |
| `verify/verify_materials.py` | 76/76 materials checks passed |
| `verify/verify_cw_g2.py` | 17/17 CW checks passed |
| `verify/verify_dot_levels.py` | 14/14 checks passed; 8 known deviations documented |
| `verify/verify_waveguide.py` | 31/31 waveguide checks passed |
| `verify/verify_transport.py` | 34/34 transport checks passed |
| `verify/verify_rt_edge_cards.py` | 215/215 rt-edge card checks passed |
| `verify/verify_rt_edge_contract.py` | 95/97 contract checks passed (two evidence gaps) |
| `verify/verify_rt_edge_papers.py` | 14/17 rt-edge paper checks passed (three single-source evidence gaps) |
| `verify/verify_rt_edge_sweep.py` | 94/94 rt-edge sweep checks passed |
| `verify/verify_device_rt.py` | 206/206 device RT checks passed |
| `verify/verify_designer_rt.py` | 6/6 designer RT checks passed |
| `verify/verify_spec_rt.py` | 11/11 spec RT checks passed |

## Council review history

The review record is in `../_goal/PROGRESS.md`; the associated `.workers/specs/rt-fix-*.md` fix specifications are git-ignored working files.

- Opus pass 1 found stale, overconfident evidence and result wording; fix round 1 introduced explicit evidence gating.
- Opus pass 2 found the linewidth proxy and wavelength treatment insufficient; fix round 2 made the linewidth range and solver deviation explicit.
- Opus pass 3 found hidden brightness/collection assumptions; fix round 3 declared the collection levers and flux floor.
- Opus pass 4 found the decisive eligible corner had to be evaluated with the final cards; fix round 4 regenerated the 192-row sweep and conditional verdict.
- A codex-sol cross-review found that the human-facing documents still described an earlier sweep; this final documentation round refreshes the verdict, counts, disclosures, and conditional conclusion.
- Council review pass 5 refreshed the gamma300 threshold bracket, flux margin, 6.5 meV per-card anchor, like-for-like Reischle comparison, and fused facet-model note.
- Council review pass 6 (Opus, 2026-09-06) confirmed the physics reproduces to float precision (facet factor applied once, 1544 photons/s factor chain, 8-10 meV threshold bracket) and found reporting-attribution errors; fix round 6 made the residual-background citation ledger-driven (Opt. Express 16, 12771), prints `none(flux)` / `none(g2)` instead of a misleading threshold when nothing passes, names the smallest factor of the whole brightness chain (retention S) as the limiter, unifies the coverage/median names with the VERDICT line, adds the IRF-deconvolved Reischle 0.25 +/- 0.05 anchor with [V] tags, and rewrites the cards' signal-fraction sentences (Lemma 1: rho is independent of the collection levers).
- Council review pass 8: lecture-deck physics reviews found the NA-collection bug; the circular NA-cone correction refreshed the flux values and temperature-axis conclusion.

## Limits

The single-band confinement solver is a deliberately limited model of a large-mismatch 3D island; it does not reproduce the HKUST 750 nm line and should not be used to erase that mismatch by substituting a preferred wavelength. The simulator supplies a transparent decision envelope and measurement priorities, not a fabricated device geometry or a measured single-photon claim.
