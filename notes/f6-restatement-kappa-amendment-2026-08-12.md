# Paper trail: F6 restatement + kappa-ceiling amendment (work order v3, V5)

**Date:** 2026-08-12 · **Trigger:** T-1 slit-held implementation + T-R corner
retarget (commit c65d5f7; baseline tag `tier-runs-2026-08-12`).

## V5.1 — F6 tracking rule, restated

**Old phrasing** (foundation/planning docs): "keep the slit below T_target" /
"external slit filter retained below T_target."

**Problem found by the T-1 implementation:** the phrasing is ambiguous, and
one reading is wrong. A slit *fixed at the absolute design-point energy* (or a
slit that follows the cavity mode — the pre-T-1 model behavior) leaves the XX
line CLOSER to the acceptance center than X below T_target: at 77 K under a
120 K-tracked cavity, X sits +6.8 meV and XX +3.3 meV off the mode — the F6
inversion (eps = 3.25, g2 > 1). Verified in `verify_gf.py` (slit-held device
check) and `out/tier_device/mechanism_runs_slit_held_12pair.csv`.

**Restatement (standing rule):**

> **The filter tracks the *emitter*, never the mode.** The external slit is
> centered on the X line at each operating temperature (the lab monochromator
> convention); the cavity mode walks per its own dn/dT. Below T_target the
> slit therefore rejects XX (at -Delta_XX from center) while the mode sits
> blue of the line; the two acceptance centers coincide only at the design
> point.

**Coincidence caveat:** the design-point coincidence is at **T_junction =
T_target, not T_heatsink = T_target** — junction heating offsets the slit-mode
coincidence by Delta_T_J (~1.6 K at the staged 1 um/10 uA point, worth
-0.31 meV of mode walk). Encoded in `verify_gf.py`.

Model support: `FilterBlock.track = "mode" | "hold"` with
`spectral.transmission2/epsilon2` (independent slit/cavity centers) and the
`ibm_transmission(delta_c_meV=...)` extension; `"mode"` preserves the legacy
shared-center behavior bit-identically.

## V5.2 — kappa-ceiling amendment (written, not silent)

**Amended:** the staged-spec line "kappa <= 1.03-1.06 meV at the 120 K / 0.1
cell" (spec sheets `out/spec/`, work order v1/v2 context).

**Statement:** the staged kappa <= 1.06 meV ceiling was a
**Lorentzian/Poisson-era proxy** for the g2 <= 0.1 target. With the IBM
lineshape and the F8 loading factor in the chain, the binding constraint is
**g2 <= 0.1 per rail at the design point**, which under quiet statistics
(F_p = 0.5 [A]) admits kappa ~ 3 meV — and the brightness floor then prefers
it (t_X falls monotonically with pair count). The T-R solver
(`scripts/run_post_tier_gates.py`) therefore optimizes max t_X subject to
g2 <= 0.1 and Delta_T_J <= 3 K with NO independent kappa ceiling; the
admissible cavity window on today's inputs is 10-12 pairs (baseline and
aggressive-transparency contrast) to 10-14 (conservative contrast), and the
Poisson rail admits no configuration at 120 K at all.

The old ceiling remains correct *as what it was*: the narrowest-acceptance
reading of the Lorentzian 0.1 requirement without loading statistics. Spec
sheets regenerated after V3 (dispersion) land should carry this amendment in
their headers.

## Status of the numbers

Everything above is structural. Specific pair counts, g2, and t_X values
remain gated on V1 (band confirmation), V2 (rail magnitude), V3 (dispersion),
V4 -> staging map, T-4/V6 (Delta_XX) per the v3 guardrail.
