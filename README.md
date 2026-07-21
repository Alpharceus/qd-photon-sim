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
python scripts/run_phase0.py      # V-a fit + report bundle
python scripts/run_phase1.py      # V-a recheck + T_j maps + V-c analysis
streamlit run fsim_gui/app.py     # Phase-V GUI (spectral explainer, cascade,
                                  #   dashboard, card editor; zero physics)
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

**Phase V (GUI v1): delivered 2026-07-21.** `fsim_gui/app.py` — Streamlit thin
client (three-layer rule enforced by a regression check): spectral explainer
(X/XX Lorentzians, Γ(T) on a temperature slider, shaded window, live ε/g²₀),
cascade diagram, results dashboard (fit-from-card button, T_c readout,
sensitivity tornado), tag-colored card editor that saves the same YAML the CLI
uses. Exit criterion verified headlessly: the GUI regenerates the Phase-0 fit
from `chatzarakis.yaml` (PASS, max |resid| = 0.011).
