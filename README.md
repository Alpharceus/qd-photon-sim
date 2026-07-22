# FSIM — F-series simulator

Phase 0 vertical slice (see `../FSIM-planning-doc-2026-07-21.md`): card schema,
Module C (spectral closed forms), Module E (integrator), the Chatzarakis V-a
fit, and its figure/CSV bundle.

## Layout (three-layer rule)

- `fsim_core/` — headless physics: `card.py` (schema + [V/DR/E/A] tags),
  `spectral.py` (Module C: F1/F1a, top-hat/cavity transmissions, MC second
  methods), `integrator.py` (Module E: background law, master ceiling → T_c,
  OAT sensitivity), `fitting.py` (Phase-0 V-a fit driver).
- `fsim_viz/` — figure factory; every figure ships its underlying CSV.
- `cards/` — parameter cards (YAML). Every entry: value|range, unit, tag, source.
- `verify/verify_fsim.py` — regression suite (Yuki standard: units, limits,
  MC/quadrature second method for every closed form). Run on every change.
- `scripts/run_phase0.py` — Phase-0 driver; writes `out/phase0/` report bundle.
- `out/` — generated bundles (figure + CSVs + fitted params + card snapshot).

## Run

```
python verify/verify_fsim.py      # regression suite, exit 0 iff all pass
python verify/audit_physics.py    # known-parameter physics audit
python scripts/run_phase0.py      # V-a fit + report bundle
python scripts/run_phase1.py      # V-a recheck + T_j maps + V-c analysis
python scripts/run_phase2.py      # V-b + cavity design rules
python scripts/run_phase3.py      # requirement envelopes + Osinski packet
python fsim_gui/designer.py       # DEVICE DESIGNER (Dear PyGui): block-diagram
                                  #   canvas, fab-stack editor + cross-section,
                                  #   RUN -> graphs + numbers; designs saved as
                                  #   cards/<name>-design.yaml
streamlit run fsim_gui/app.py     # validation dashboard (spectral explainer,
                                  #   cascade, V-a fit, card editor)
```

## Status

**Phase-0 gate: CLOSED (V-a PASS, 2026-07-21).** All six published g²(T)
points fit within ±0.03/point with the digitized windows, the Γ(250 K) anchor,
and E_a = 240 meV (published for this dot). Solved ceiling T_c = 261 K. See
`notes/phase0-results.md` for findings and the honesty ledger (existence, not
uniqueness; 150/210 K windows interpolated [E] pending the supplement).

**Phase 1 (drive realism): delivered 2026-07-21.** Module D (cap-2/F1b
operating point, injection background channels, MC-verified) + Module A
(analytic spreading-resistance thermal, FD-checked, runaway detection) +
`g2_electrical` (electrical-separation theorem). V-a re-passes with the F1b
penalty (T_c 261.4 → 260.6 K). ΔT_J envelope maps in `out/phase1/`; findings
in `notes/phase1-results.md`. **V-c(i) PASS** on the digitized Reischle APL
97, 143513 (2010) data: trion → ε = 0, measured g² = 0.37 requires ρ = 0.794,
inside the spectrum-digitized ρ(w) envelope — the electrical result is
ρ-limited, as the F-series claims. V-c(ii) (the 80 K point) awaits the
Opt. Express 16, 12771 (2008) PDF.

**Phase 2 + all validations closed 2026-07-22.** Foundation doc cross-checked
(exact agreement; F5 aperture lemma + Theorem-0 slice added). V-a re-passed as
a joint over-determined fit against the supplement's τ(T) [V] and Γ_X/Γ_XX(T)
[E] data (named Tier-2 refinement: effective phonon energy E_lo = 18.3 meV;
r_XX = 0.72; T_c revised 261 → 249 K). V-c(ii) closed on OE 16, 12771 (2008)
(the paper's Eq. (1) is the F2 law; residuals ≈ ε). V-b closed at the ε→1
limit (300 K point at the Theorem-0 edge; 77 K residual = named re-excitation
channel, three-paper convergence on WP-M2′). Module B (`fsim_core/cavity.py`):
tracking rule met at T_target = 120 K (`cards/qcap-cavity.yaml`); design rules
in `notes/phase2-results.md`. 43/43 checks. Run `python scripts/run_phase2.py`.

**Phase 3 (prediction runs): delivered 2026-07-22.** The four §4 deliverables
— (i) staged 77–120 K envelope, (ii) piezo-variant (Δ,ρ) requirement map vs
300 K with measured devices overplotted, (iii) F5 aperture/density rules,
(iv) measurement-priority ranking — assembled into the design-review packet
`notes/phase3-osinski-packet.md` with every claim tagged. Bundle in
`out/phase3/`. Run `python scripts/run_phase3.py`.

**Phase V (GUI v1): delivered 2026-07-21.** `fsim_gui/app.py` — Streamlit thin
client (three-layer rule enforced by a regression check): spectral explainer
(X/XX Lorentzians, Γ(T) on a temperature slider, shaded window, live ε/g²₀),
cascade diagram, results dashboard (fit-from-card button, T_c readout,
sensitivity tornado), tag-colored card editor that saves the same YAML the CLI
uses. Exit criterion verified headlessly: the GUI regenerates the Phase-0 fit
from `chatzarakis.yaml` (PASS, max |resid| = 0.011).
