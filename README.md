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

## Status / honesty ledger

`cards/chatzarakis.yaml` currently carries **PLACEHOLDER** entries: the four
interior g²(T) points, Δ_XX, and the filter window are not yet digitized from
Phys. Rev. Applied 20, 034011 (paywalled). The pipeline runs end-to-end and the
fit machinery is verified, but **the V-a gate cannot close until the digitized
values replace the placeholders** — the card and every output say so.
