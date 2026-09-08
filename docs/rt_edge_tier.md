# RT edge-emitter tier

This tier evaluates non-cryogenic heatsink operation (`T_hs >= 230 K`), electrically driven, edge-emitting InP-dot single-photon operation with the contract grid in `docs/rt_edge_contract.md`. It is a model-and-evidence assessment, not a device demonstration. The current results below are read from the generated acceptance artifacts and the current verification scripts.

## Outcome

**pkg5b (2026-09-07/08) rerun -- the headline verdict changed qualitatively, not just numerically.** The sweep now carries two opt-in physics-correction switches as genuine axes: `drive.finite_pulse` (finding 1: a real finite-pulse-waveform loading calculation, replacing the legacy static per-pulse `mu` Poisson-loading approximation every earlier package used) and `ret.tau_cap_scales_with_density` (finding 1b: the opt-in "full-cancellation" retention convention, versus the shipped "no-cancellation" default). The headline combination is `drive.finite_pulse=true, ret.tau_cap_scales_with_density=false` -- the corrected model, left at the shipped retention convention. Under it, **no sampled corner anywhere in the grid clears the 1 kHz collected-flux eligibility floor** (`eligible=0/768`): every one of the 424 rows that used to be reported eligible under the old, uncorrected static-loading approximation is now excluded by the flux floor. `g2_min`, `T_pass_min`, and `gamma300_pass_max` are therefore undefined (`nan`/`none`) rather than favorable numbers, and `conditional=false` (there is no headline-passing row for evidence completion to be conditional on).

The current `out/rt_edge/verdict.md` line is quoted verbatim:

```
VERDICT: FAIL model=finite_pulse:true,tau_cap_density:false g2_min=nan g2_median_eligible=nan diag_g2_min=0.9817 diag_g2_median_diagnostic=0.9993 flux_max=582.8 flux_margin=0.5828 flux_shortfall_deprecated=1.716 median_pass=false coverage_over_eligible=nan eligible_fraction=0 eligible=0/768 flux_floor_excluded=768 evidence=incomplete conditional=false headline_coverage=0/768 headline_coverage_pulsed=0/384 headline_dedup_mismatch_groups=0 eligible_dedup=0/384 eligible_dedup_mismatch_groups=0 rows_scheduled=0/768 cw_raw_coverage=0/0 gamma300_pass_max=nan gamma300_threshold=n/a T_pass_min=none headline_by_T=230:0/192,250:0/192,273:0/192,300:0/192 headline_by_T_pulsed=230:0/96,250:0/96,273:0/96,300:0/96
```

`out/rt_edge/verdict.md`'s "Model sensitivity" table reports all four `(drive.finite_pulse, ret.tau_cap_scales_with_density)` combinations side by side: the old, uncorrected default every package before pkg5b silently used, `finite_pulse=false, tau_cap_density=false`, ALSO now reports `0/768` eligible on the current (already-corrected-elsewhere, e.g. NA-collection and waveguide dispersion) physics at HEAD -- the previous conditional-PASS-adjacent result was already stale before the finite-pulse correction was even applied. Only the opt-in `tau_cap_density=true` ("full-cancellation") retention convention restores any eligible or passing rows at all: `finite_pulse=false, tau_cap_density=true` gives 392/768 eligible, 170 headline passes, `g2_min=0.2995`, `flux_max=4767` photons/s; `finite_pulse=true, tau_cap_density=true` gives 496/768 eligible but 0 headline passes (`g2_min=0.8291`, above the 0.5 gate). Both cards report `eligible=0/384` and `favorable_rows=0` under the headline combination; the maximum collected pulsed flux anywhere in the headline grid is 582.8 photons/s, 1.72x below the 1000 photons/s floor. The measured single-row cost extrapolated to 25.3 minutes for the full four-model-combination grid, well under the 90-minute budget, so the full (non-reduced) grid was run for all four combinations (wall time 1561.1 s = 26.0 minutes).

Evidence is still incomplete for the same three central claims as before (unaffected by this rerun): there is no second verified source for Reischle 2008's 80 K electrical `g2` or for the transferred 80 K residual-background claim, and the HKUST 750 nm line is not reproduced by the single-band solver. Fail reasons are now `evidence_incomplete` AND `no_eligible_rows` (previously just `evidence_incomplete`): even a complete evidence ledger could not turn this run into a PASS.

## Comparison points and interpretation

The model's relevant comparison values must be interpreted as predictions, not new measurements. With no eligible corner under the headline model, there is no finite `g2_min` to compare like-for-like against Reischle et al. (2008)'s IRF-deconvolved, background-included anchor (`g2(0)=0.25 +/- 0.05`); `out/rt_edge/verdict.md`'s Literature ceiling section states this explicitly ("No finite pooled g2_min is available in this run to compare against Reischle's deconvolved value"). The Chatzarakis et al. (2023) 6.5 meV per-card anchor check also now fails on both cards at their own default `delta_xx` (`edge-inp-gaasp-design`: g2_pulsed=0.9999, flux=2.285 photons/s; `edge-inp-gainp-design`: g2_pulsed=0.9979, flux=76.68 photons/s -- neither passes).

The key practical result is therefore stark: under the corrected finite-pulse loading model, non-cryogenic operation does not pass at any sampled heatsink temperature or linewidth in this grid (`T_pass_min=none`, `gamma300_pass_max=nan`, reason `none(flux,g2)` for both cards). The platform's apparent viability in every prior package's reported verdict rested on the uncorrected static-loading approximation; whether the `ret.tau_cap_scales_with_density=true` convention that does restore eligible rows is itself the physically correct retention convention for this platform is exactly the open question finding 1b raised and this sweep does not resolve.

## NA-collection correction

The waveguide review found that the old collection estimate used the wrong separable formula inside a circular NA cone. The corrected model integrates the circular cone and uses the per-axis Gaussian factor `erf(sqrt(2) theta_max/theta_div)`, with `theta_div=lambda/(pi w0)`; the review also standardizes the 1/e² width convention. At the time of that review, the favourable diagnostic corner moved from NA collection `0.548244` (1544.4 photons/s) to `0.857853` (6890.03 photons/s, about 4.46x higher) -- both figures computed under the since-superseded static per-pulse loading approximation. This NA-collection correction itself is unaffected by pkg5b's finite-pulse loading correction (it is a purely geometric/optical factor); what has since changed is the loading term multiplying it: under the current headline model (`drive.finite_pulse=true`) the same corner's collected flux is 582.8 photons/s (NA collection itself is 0.852617, essentially unchanged -- see `out/rt_edge/verdict.md`'s brightness-factor table). The updated flux decomposition still identifies confinement retention `S` as the limiting factor. This is a model correction, not an experimental brightness measurement.

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

Under the pkg5b headline model (`drive.finite_pulse=true, ret.tau_cap_scales_with_density=false`), no configuration passes: the best diagnostic (lowest-g2) corner (NA=0.8, R_back=0.95, L_um=250, both cards' own aperture density declared, `edge-inp-gainp-design` at `delta_xx=8 meV`, `gamma300=6 meV`, `T_hs=230 K`) reaches 578.843 photons/s -- close to, but not the same row as, the grid's actual maximum collected flux (`flux_max=582.8` photons/s, the VERDICT line's value), which occurs at a DIFFERENT corner (`edge-inp-gainp-design`, `delta_xx=4 meV`, `gamma300=20 meV`) with a worse g2_pulsed=0.9997. Re-excitation inside the 100 ps pulsed-counting gate (a second capture-emission cycle before the gate closes) makes higher `gamma300`/lower `delta_xx` corners brighter but NOT more antibunched, so the flux-maximizing and g2-minimizing corners diverge; neither reaches the 1000 photons/s floor. The flux floor is 1000 photons/s, so the maximum flux margin is 0.5828 (`flux_max`/floor, below 1.0 -- the floor is NOT cleared anywhere in this grid). The facet model is a ray-probability series with the dot at mid-ridge: `eta_facet = 0.5*T*exp(-a*L/2)*(1+R_back*exp(-a*L)) / (1-R_back*R_front*exp(-2*a*L))`, with `R_front = 1-T`, `R_back=None` meaning the bare cleaved facet (`R_back=R_front`), and single-pass propagation already folded into `eta_facet` (`[DR]` Coldren, Corzine & Masanovic, *Diode Lasers and Photonic Integrated Circuits*, 2nd ed., ch. 2). At this configuration it evaluates to `eta_facet=0.7840154`, independently matching the sweep's own back-solved `eta_total/(beta*eta_NA)` (`verify/verify_rt_edge_sweep.py`) -- this facet-model self-check still holds; it is the collected-flux magnitude (dominated by the finite-pulse-corrected loading term and confinement retention `S`, not the facet model) that no longer clears the floor.

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
| `verify/verify_rt_edge_sweep.py` | 130/130 rt-edge sweep checks passed |
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
- pkg5b (2026-09-07/08): peer-review packages 1-6 (finite-pulse loading, tau_cap density scaling, and several unrelated fixes) landed but the acceptance sweep had not been rerun against them; `drive.finite_pulse` and `ret.tau_cap_scales_with_density` are added as genuine sweep axes (rather than fixed at their card defaults), the corrected `finite_pulse=true, tau_cap_density=false` combination becomes the headline, and the other three combinations are reported in a new "Model sensitivity" table. The rerun found `eligible=0/768` under the headline model -- no sampled corner clears the collected-flux floor -- superseding every prior package's conditional-near-PASS result, which rested on the uncorrected static-loading approximation.

## Limits

The single-band confinement solver is a deliberately limited model of a large-mismatch 3D island; it does not reproduce the HKUST 750 nm line and should not be used to erase that mismatch by substituting a preferred wavelength. The simulator supplies a transparent decision envelope and measurement priorities, not a fabricated device geometry or a measured single-photon claim.
