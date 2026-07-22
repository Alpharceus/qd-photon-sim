# FSIM Phase-3 design-review packet — [FOR OSINSKI]

**Program:** NSF NQVL QCAP SLE · prepared 2026-07-22 · repo `fsim/` (all
numbers regenerable: `python scripts/run_phase3.py`; figure + CSV per panel in
`out/phase3/`)

## What stands behind these predictions (validation record)

The governing F-series mathematics is implemented, MC/FD-verified (43/43
regression checks), and **validated against all three published datasets the
plan names**:

| Arm | Dataset | Result |
|---|---|---|
| V-a | Chatzarakis g²(T) 78–230 K + supplement τ(T) [V] + Γ(T) [E] | joint over-determined fit PASS, max resid 0.028; T_c = 249 K |
| V-b | Laferrière 0.11/0.34/0.57 (77/220/300 K) | ε→1 limit confirmed at 220/300 K; 300 K point at the Theorem-0 edge |
| V-c | Reischle 2010 (~40 K pulsed) + 2008 (80 K DC) | ρ-limited with ε small; the 2008 paper's own Eq. (1) is the F2 law |

Named model residuals (open, none blocking): Γ(T) high-T shape (Tier-3
independent-boson candidate); re-excitation/refilling channel (WP-M2′,
three-paper convergence).

**Every number below carries tag chain [A]** — requirement envelopes over
swept unmeasured inputs, never point predictions (schema-enforced).

## (ii) The 300 K question — piezo-variant requirement map

T_c(Δ, ρ) at the narrow-filter bound, Γ(300 K) anchored to the published
6–7 meV [V]:

- **[FOR OSINSKI] 300 K is a (Δ ≥ ~5 meV) AND (ρ ≥ ~0.81–0.86) problem.**
  At the (211)B yield line (Δ > 5 meV in >50% of dots [V]): ρ ≥ 0.825–0.863.
  At the demonstrated median Δ = 5.4 [V]: ρ ≥ 0.809–0.843. At a selection-tail
  Δ = 8: ρ ≥ 0.755–0.772. (Reproduces the F7 flow-down independently through
  the fitted class Γ(T).)
- The measured devices already sit at or above these ρ values *at their own
  temperatures* (Chatzarakis 230 K: ρ = 0.884 [derived from the V-a fit];
  Reischle 80 K: 0.88 [V]) — **the 300 K problem is holding ρ there at T_j =
  300 K under injection, i.e. background-and-collection engineering, not
  emitter physics.** Direction fact, not a tweak.
- If the InP/GaAsP platform's Δ stays at its class-typical ~1.8–3 meV, 300 K
  is closed at zero background (ρ_req > 0.94 → impossible with any injection
  background); only orientation/strain engineering of Δ [A, unmeasured]
  reopens it. The staged 77–120 K deliverable is untouched by this.

## (i) Staged InP/GaAsP device, 77–120 K (anchored deliverable)

Electrical envelope (class retention proxy [A]; w = Γ operating convention
[A]; ΔT_J = 2.3–2.9 K at 1 µm mesa / 10 µA CW worst epi-k — thermally benign):

- **77 K: g² ≤ 0.5 guaranteed worst-case for Δ ≥ 2.5 meV**; the best-case edge
  reaches g² ≈ 0.05–0.08 at Δ = 5. The g² ≤ 0.1 target is *not guaranteed*
  worst-case anywhere in the swept box — it requires bounding b_e (below).
- **[FOR OSINSKI] 120 K fails structurally if the injection background reaches
  b_e ≈ 0.3**: with the class p-shell quench (S(123 K) ≈ 0.53 [A proxy]),
  ρ² < 0.5 regardless of Δ. The injection-background budget at 120 K is
  **b_e ≲ 0.18** (signal units). Two caveats cut both ways: the InP/GaAsP s-p
  spacing may exceed the arsenide 35 meV (less quench, more headroom) — same
  unmeasured-Γ/level-structure item as the PL queue below.

## (iii) Aperture/density design rules (F5 lemma)

- **The F5 lemma is unforgiving: one 30%-bright in-window competitor costs
  g² = 0.36 by itself.** Penalty ≤ 0.05 at a 1 µm aperture requires density
  ≤ 2×10⁸ cm⁻²; at native 10¹⁰ cm⁻² even 0.3 µm apertures sit near N_w ≈ 1
  at high T (foundation checkpoint reproduced: N_w ≈ 1.6 at 0.1 µm²,
  Γ = 6.5 meV).
- **[FOR OSINSKI] low-density growth (~7×10⁸ cm⁻² class) plus ≤ 0.5 µm
  apertures is a requirement, not an optimization**, for any high-T target.
  (Continuous-N approximation of the Poisson mixture; exact treatment is a
  one-line refinement if a decision ever hinges on it.)

## (iv) Measurement-priority ranking (in-house queue)

Which measurement narrows the staged-device T_c envelope (169 K wide over the
[A] box) the most:

| Rank | Measurement | Envelope narrowing |
|---|---|---|
| 1 | **injection background b_e(I,T)** — EL spectra under injection | 47.6 K |
| 2 | **Δ_XX distribution on InP/GaAsP** — in-house PL (zero published values [V, null]) | 39.2 K |
| 3 | **Γ(T) on InP/GaAsP** — in-house PL | 31.8 K |
| 4 | operating loading μ — power series | 1.4 K |

- **[FOR OSINSKI] b_e outranks Δ_XX** on the *staged* device (its assumed
  range spans two decades). For the *300 K route* the order inverts: Δ is
  binding there (map (ii)). The two queues are consistent: **measure b_e first
  for the staged deliverable; measure Δ_XX first for the 300 K decision.**
  Γ(T) rides along with either PL campaign. μ is not worth dedicated beam time.
- The ranking is itself [A]-conditioned (it depends on the assumed ranges);
  its job is to order the first measurements, after which it should be re-run
  with the narrowed card — the intended loop.

## Standing design rules from Phase 2 (carried into any cavity design)

1. Mode-tracking rule met at T_target = 120 K: mode 24.0 meV red of the
   cryogenic line (`cards/qcap-cavity.yaml`).
2. **A red-tracked cavity favors XX below T_target** (cavity-only ε > 1 at
   20 K): retain the slit filter below T_target.
3. Δ = 1.5 meV exceeds the ceiling at T_target even with the cavity — the
   splitting requirement binds before the cavity does.
4. κ window: Γ(T_target) ≲ κ ≪ Δ (Purcell wants κ large, filtering wants it
   small); F_eff/F_P = 0.11–0.56 over κ = 0.3–3 meV.
5. Thermal (Phase 1): mesa ≥ 1 µm keeps ΔT_J < 3 K at 10 µA CW; substrate
   choice (GaAs vs GaAs/Si) is thermally second-order for etched mesas.

## Reality check (risk ledger §6 restated)

FSIM cannot settle what only measurement supplies. The decisive experiment
remains the in-house InP/GaAsP Δ-distribution + Γ(T) campaign; this packet's
job was to make that experiment maximally informative (ranking above) and
every preceding design decision a computed one (maps and rules above).
