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
```

## Status

**Phase-0 gate: CLOSED (V-a PASS, 2026-07-21).** All six published g²(T)
points fit within ±0.03/point with the digitized windows, the Γ(250 K) anchor,
and E_a = 240 meV (published for this dot). Solved ceiling T_c = 261 K. See
`notes/phase0-results.md` for findings and the honesty ledger (existence, not
uniqueness; 150/210 K windows interpolated [E] pending the supplement).
