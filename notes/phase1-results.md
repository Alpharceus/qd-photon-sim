# Phase-1 results note — drive realism (2026-07-21)

## Delivered

Module D (loading & drive statistics) and Module A (analytic thermal), wired
into the integrator as the **electrical-separation theorem**: drive enters
g²(T; design) only through ΔT_J (Module A) and Δρ (Module D injection
channels); the ε path is untouched by construction. 34/34 regression checks
(every new closed form has its MC or FD second method).

## F1b (finite-μ operating point) — derived and verified

Cap-2 Poisson loading, peak-area convention:
g²(μ) = 2P₂ε/[P₁+P₂(1+ε)]², with F1 (g²=ε) recovered as μ→0 and saturation
limit 2ε/(1+ε)². The entire finite-μ penalty is a factor f(μ) ∈ [1, 2).
MC-verified, including the exactness of the background law
g² = 1−ρ²(1−g²_dot) for Poissonian background over any dot statistics.

**V-a recheck with μ_op = 0.33** (paper: pump 3× below saturation, [E]):
still PASS, max |resid| = 0.011; T_c moves 261.4 → 260.6 K. The drive
penalty at the published operating point costs ≈ 1 K of ceiling.

## Module A findings (T_j envelope maps, tag chain [A])

Staged-device mesa (etched through a 1.5 µm epi stack, epi-k swept 8–40 W/mK),
CW at 1.9 V, T_hs = 77 K:

- **Mesa diameter is the thermal design lever.** Worst-case ΔT_J at 10 µA:
  ~28 K for a 0.3 µm mesa vs < 0.3 K at 3 µm. At 50 µA a 0.3 µm mesa can sit
  ~190 K above the heatsink; at 200 µA the k(T) feedback predicts **thermal
  runaway** below ~0.4 µm (reported as ∞, plotted as a gap — not extrapolated).
- **Substrate choice is thermally second-order for small mesas.** GaAs vs
  GaAs/Si (with a defective 2 µm buffer) differ by < 2% in ΔT_J across the
  map: the mesa pillar's 1D resistance dominates. The Si advantage only
  appears for large/planar geometries. Consequence: GaAs/Si integration is
  not a thermal penalty for the QCAP mesa geometry at these scales.
- **The dominant unknown is the epi-stack effective conductivity** (alloy
  scattering; anisotropy) — the band spans ~4× in ΔT_J. It heads the
  measurement queue for the thermal path, ahead of any COMSOL refinement.

Caveats: analytic cone/pillar model (FD-checked at the 8% level for the disk
limit; cone approximation unvalidated beyond ~25% class); k(T) exponents for
alloys assumed [A]; heat assumed generated at the junction plane; no surface/
interface (Kapitza) resistance. COMSOL replaces this per plan at Phase-1 exit.

## V-c status: BLOCKED (not failed)

`reischle.yaml` carries placeholders. Needed from Reischle et al., APL 97,
143513 (2010): g²(0) at 80 K under pulsed electrical drive, X/XX splitting,
filter window, drive conditions (I, V, repetition rate). Drop the PDF in the
project root and re-run `scripts/run_phase1.py`. The V-c claim to test: the
80 K electrical result is ρ-limited (injection background + ΔT_J) with ε
subdominant.

## Phase-1 exit criteria

- T_j map for the staged-device mesa on GaAs/Si: **done** (`out/phase1/`).
- V-c consistency: **open, blocked on source data** — the only remaining item.
