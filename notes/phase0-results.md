# Phase-0 results note — V-a gate: PASS (2026-07-21)

## Claim (framed before the fit, per council rule)

One physically constrained {ε(T), ρ(T)} parameterization reproduces the
Chatzarakis six-point g²(T) series within ±0.03/point, using the published
filter windows, the Γ(≥250 K) = 6–7 meV anchor, and an extracted E_a
consistent with the published activation energy. (Planning doc §1, V-a.)

## Result

**PASS.** Max |residual| = 0.011 across all six points (78 K point treated as
the published upper bound g² ≤ 0.02, one-sided). Γ(250 K) = 6.47 meV (anchor
6.5±0.5). Extracted E_a = 240.4 meV vs the published **240±20 meV for this
exact microcavity dot** (Fig. S5(b)) — note the planning doc's 265±30 meV is
the ⟨Al⟩=65% *macro* sample; the criterion is satisfied under both.
Solved master ceiling: **T_c = 261 K** — consistent with the paper's
observation of resolved X/XX to 260 K and its claim that minor modifications
reach room temperature.

## What the data forced into the model (findings)

1. **Windows move.** The published spectral windows both widen (1.4 → 7.7 meV)
   and slide blue of X (dx −0.4 → −3.4 meV) with temperature — away from the
   red-shifted XX (Δ_XX = 5.9 meV, digitized). Module C gained a window-offset
   parameter `dx`; a centered-window model cannot reproduce the series.
2. **Two-channel retention + coupled background.** A single-Arrhenius ρ(T)
   with constant background cannot satisfy g²(78) ≤ 0.02 and g²(120) = 0.10
   simultaneously. The paper's own p-shell channel (E_b = 35±5 meV, prefactor
   ~100) plus a background that grows with carrier escape,
   B(T) = b0 + β(1−S(T)), resolves it. This is the planning doc's "background
   channels" made concrete.
3. **Decomposition.** At 230 K the model attributes the g² rise ≈ 2/3 to ρ
   (retention/background) and ≈ 1/3 to ε (XX leakage through the widened
   window). The paper's qualitative attribution ("relative increase of the
   background level") is reproduced, with the ε share quantified.

## Honesty ledger (caveats carried forward)

- **Existence, not uniqueness.** 8 free parameters against 5 values + 1 bound
  + 3 anchor/prior residuals: V-a asks whether a physical parameterization
  *exists* (it does); parameter values are not uniquely determined. The
  Γ-decomposition (a_ac vs b_lo) is degenerate given one anchor — a_ac sits at
  its range top (0.02 meV/K). Do not quote fitted parameters as measurements.
- **150 K and 210 K windows are interpolated [E]** (Fig. S4 supplement not
  available); replace when obtained. All other windows digitized from the
  Fig. 4 raster, tick-calibrated (±0.1 meV edges; X-peak ±0.5 meV at high T).
- Tag chain of every output: **[A]** (fit parameters are assumptions
  constrained by this one dataset).
- X⁺ sits 1.5 meV below X, inside the low-T windows: treated as antibunched
  signal (paper's cross-correlation result), not background.

## Gate consequence

Phase-0 exit criterion met → Phase V (GUI v1) and Phase 1 (Modules D/A) are
unblocked per the plan. Laferrière (V-b) and Reischle (V-c) cards are the next
validation targets.
