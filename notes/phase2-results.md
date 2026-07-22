# Phase-2 results note — foundation cross-check, all three validations closed, Module B (2026-07-22)

## 0. F-series foundation doc cross-check

`F-series-spectral-foundation-2026-07-21.md` vs the implementation: **exact
agreement**. F1 matches `f1b_g2` term for term (incl. cap-2 P₂ = 1−e^(−μ)(1+μ));
F1a matches the transmissions (doc in half-width w, code in full width — same
formulas); the doc's F1b expansion ε[1+μ(⅓−ε)] is the first-order Taylor of the
exact closed form (verified algebraically and as a regression check); F2/F3/F5
identical. Newly implemented from the doc: the **F5 aperture lemma**
(g² = 1−Σp²/(Σp)², MC-checked at the doc's verify3 cases) and the **Theorem-0
slice** checks (ε=1: g²→1 at μ→0, →½ at μ→∞). Suite: 43/43.

## 1. V-a upgraded: joint over-determined fit (supplement digitized)

New data: Fig. S4 **printed X decay times** 78–230 K [V] (τ ratio = independent
retention constraint S(T)/S(78)); Fig. S2 **Γ_X and Γ_XX(T)** digitized, 8–260 K
[E: mesa dot used as class proxy for the cavity dot]. Fig. S4 has **no shaded
windows** — the 150/210 K windows remain interpolated [E]. S5(b) confirms
E_a = 240 meV for the cavity dot.

- First joint fit **FAILED** V-a (+0.041 @ 230 K) and localized the failure:
  the Γ(T) model with E_LO fixed at 36.6 meV cannot produce the S2 shape.
  Named Tier-2 refinement: **effective phonon energy freed** → E_lo = 18.3 meV
  (alloy/interface-phonon range), acoustic term → 0 (degeneracy resolved).
- Refit **PASSES**: max g² |resid| = 0.028; r_XX = 0.72 (Γ_XX/Γ_X, matching S2 —
  the Γ_X=Γ_XX ledger assumption is measurably wrong and now fitted);
  E_a = 220.8 meV (1σ of the published 240±20); b_p ≈ 24 ("order 100" prior, 1σ);
  worst τ-ratio error 18.5%, worst Γ error 34%.
- **T_c revises 261 → 249 K.** Known residual: model overshoots Γ(≥240 K)
  (8.4 vs 6–7 meV at 250 K) — the three-term Γ model cannot do both the steep
  90–150 K rise and the high-T flattening. Tier-3 candidate (independent-boson
  lineshape) named, not built (anti-speculation rule).

## 2. V-c(ii) CLOSED (Reischle OE 16, 12771 (2008), 80 K DC)

The paper's own Eq. (1) **is** the F2 background law. F2 inversions reproduce
the published corrected values exactly: dot B (4 K): g2_s = 0.038±0.057 vs
published 0.04; dot C (80 K): 0.032±0.065 vs 0.03. Computed ε bound
(Γ(80 K) = 0.8 meV [V], Δ 4–6 meV [DR], w 1–3 meV swept) ≤ 0.033 — the residual
is consistent with pure cascade leakage. **The 80 K electrical result is
ρ-limited with ε small — V-c fully closed** (V-c(i) 2010 device + V-c(ii) 80 K).

## 3. V-b CLOSED at the ε→1 limit (Laferrière Nano Lett. 23, 962 (2023))

Series 0.11/0.34/0.57 at 77/220/300 K; published windows 0.1 nm → 12 nm
(175 K) → 25 nm (300 K); background fitted out by the authors → ρ ≈ 1; Δ_XX,
Γ(T), per-T windows unpublished → swept [E] (existence check, not a fit).

- **300 K**: ε(sweep) = 0.83–1.00 under the 25 nm window → the ε→1 corridor
  g² ∈ [0.5, 1); measured 0.57 sits at the saturated **Theorem-0 edge** — the
  foundation's unification slice observed in a third independent platform.
- **220 K + 300 K simultaneously covered** for Δ_XX ∈ [4.5, 5.5] meV — a
  falsifiable inference: V-b favors a large-splitting dot (measurable).
- **77 K** exceeds the spectral envelope by ~+0.04: the residual matches the
  paper's own re-excitation channel (0.021 at 4 K saturation, center-dip
  signature, growing with T). Non-cascade, outside cap-2 — the **same WP-M2′
  channel** as the Reischle-2010 refilling residual. Three papers now converge
  on this named refinement tier.

**All three validation arms (V-a, V-b, V-c) are closed at working level.**
The F-series survives contact with every published dataset the plan names,
with two named model residuals (Γ high-T shape; re-excitation/refilling tier).

## 4. Phase 2 — Module B (analytic tier)

`fsim_core/cavity.py`: Varshni emitter walk, cavity-mode drift, tracking rule,
F_eff = F_P·κ/(κ+Γ), collection gain → ρ (F6 ii, wired into the integrator),
SiN β brightness-only (Lemma 1 **enforced by a regression check**: β cannot
appear in any g² path source). Deliverable card `qcap-cavity.yaml`
(T_target = 120 K): mode placed 24.0 meV red of the cryogenic line,
detuning(120 K) = 0 — the tracking rule satisfied. Envelope with tracked
cavity (κ = 1.8 meV): T_c = 120 K (Δ=1.5, G=5) → 185 K (Δ=5, G=15), tag [A].

**New design rules out of Module B:**
1. **Below T_target a red-tracked cavity favors XX** (cavity-only ε = 1.8 at
   20 K — it selects the wrong line). The slit filter must be retained below
   T_target; the cavity is a safe filter only near/above it.
2. **Δ = 1.5 meV exceeds the ceiling at T_target even with the cavity**: the
   splitting requirement, not the cavity, is binding (F7 restated) — the
   in-house Δ_XX measurement stays the program's decisive experiment.
3. F_eff/F_P = 0.11–0.56 over κ = 0.3–3 meV at Γ(120 K) ≈ 2.4 meV: Purcell
   is overlap-throttled; κ ≳ Γ wanted for rate, κ ≪ Δ for filtering — the
   κ window is Γ(T_target) ≲ κ ≪ Δ, quantified per design in the ε map.

MEEP/COMSOL supplies the actual {κ, F_P, G} numbers [A → computed]; the
Phase-2 exit criterion (a design card satisfying the tracking rule at
T_target) is met at the analytic tier.

## Ready for Phase 3

Remaining [E] items that do not block Phase 3: Chatzarakis 150/210 K windows
(only in unpublished form); Γ-model high-T shape (named). The prediction
machinery (envelope maps, sensitivity ranking, aperture statistics via F5) is
in place.
