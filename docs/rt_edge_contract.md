# Room-temperature edge-emitter contract

This is the opt-in contract for the room-temperature electrical edge-emitter tier. Legacy cards retain their existing behavior unless they explicitly select a new mode. Literature anchors live in `verify/data/rt_edge_anchors.yaml`; `[V]` means verified published value, `[DR]` derived from published values, `[E]` class estimate, and `[A]` assumption.

## Block schema

The new device fields are `dot.linewidth` (`"class"` is legacy; `"anchored"` opts into the anchor-selected linewidth), `ret.mode` (`"proxy"` is legacy; `"confinement"` opts into confinement retention), and `drive.mode` (`"EL"` or `"PL"` remain valid; `"EL-transport"` opts into transport-aware electrical drive). The new `emission` block has `type`: `"none"` (legacy) or `"edge"`. `aperture.compose` is boolean and defaults to `False`; `filter.track_material` defaults to `""`, meaning legacy hard-coded GaAs filtering; `drive.cw` defaults to `False`. Every stated default must reproduce the legacy path.

## Units and normalization

Energies are meV unless a field ends in `_eV`; temperatures are K; currents are uA; rates are ns^-1; lengths are nm except fields ending `_um` or `_um2`.

`b_e` is background photon counts inside the detection window per collected X photon, evaluated at the operating current and junction temperature. This reconciles the historical meanings: the `device.py` Arrhenius amplitude is first evaluated as its legacy background channel and then normalized to collected X counts; the card phrase “total at op point” is already that ratio; and the `spec.py` raw quantity must be converted to that ratio before use. Thus the legacy Arrhenius channel maps to `b_e = B_arrhenius(I,T_j)/S_X,collect(I,T_j)`; it is not an independently added background term.

## Temperature and rate mapping

`thermal.T_hs` is the heatsink set point. `T_j` is the self-consistent junction temperature and is the temperature consumed by anchored linewidth, retention/escape, electrical transport, and CW escape/refill. `T_hs` remains the legacy input only when the legacy thermal path has no self-heating calculation. Rate inputs and outputs use ns^-1; no implicit ps^-1 conversion is permitted.

## Aperture assumptions

The aperture-derived competitor count is continuous in the new composed model: its excess contribution produces a continuous g2 penalty, so changing area or density cannot create integer discontinuities. Integer rounding is legacy-only, retained solely when `aperture.compose: false`. Composition combines the aperture selection penalty with the normalized in-window background once, rather than treating the competitor count as a second collection loss.

## Lemma 1

**Lemma 1.** Beta factor, out-coupling, and collection efficiency are counted exactly once, in brightness. They must never lower g2, except through the existing spectral-leakage epsilon path. In particular, changing `emission.type` from `none` to `edge` changes collected brightness, not intrinsic multiphoton probability.

## Acceptance gates

The headline metric is pulsed intrinsic `g2(0)`. `g2_cw0` and `g2_cw0_raw` are secondary diagnostics. PASS requires at least one physically eligible corner and every paper gate to pass. Any missing-evidence anchor blocks PASS. A corner relying on `[A]` or `[E]` inputs must be tagged as such in its report, even if its numerical metric is favorable.

## Evidence status

The YAML table is the machine-readable evidence ledger. The temperature trend, GaAs p-i-n I-V anchors, vertical-reference extraction efficiencies, InP/GaInP X-XX proxy, the 300 K linewidth class, and the InP/AlGaInP retention activation energy now have two independent verified sources each. The 300 K linewidth class adds Matsuda et al., Phys. Rev. B 63, 121304(R) (2001) (12 meV, In0.5Ga0.5As/GaAs dot) alongside Chatzarakis et al., Phys. Rev. Applied 20, 034011 (2023) (6-7 meV, InAs/GaAs dot); both remain class proxies, not InP/GaAsP measurements. The retention activation energy adds Schulz et al., physica status solidi (c) 6, 906 (2009) (150 meV, ensemble Arrhenius fit, reference InP/Al0.20Ga0.80InP barrier) alongside Bommer et al., J. Appl. Phys. 110, 063108 (2011) (96 +/- 7 meV, single-dot).

The 80 K electrical g2 still lacks a second independent 80 K electrical single-dot measurement: Reischle et al., Appl. Phys. Lett. 97, 143513 (2010) is the only other electrically driven InP-dot g2 paper in the fact-checked digests, but it reports 20-40 K (Joule-heated), not 80 K, so it cannot honestly stand in for the missing anchor. The HKUST InP/GaAsP wavelength as a single-photon anchor now has a confirmed, DOI-sourced second wavelength candidate -- Luo et al., Optics Express 30, 40750 (2022), 750 nm InP/GaAsP QD laser on Si -- but its anchor's status is kept `missing` rather than `verified`, because promoting the primary anchor (gu25-hkust-wavelength) to `verified` would make verify_rt_edge_cards.py reject cards/edge-inp-gaasp-design.yaml's emission.lambda_nm (668 nm), which is deliberately pinned to the only tabulated phosphide index point rather than to the unconfirmed HKUST target wavelength. Both entries remain `missing` and therefore block PASS. The X-XX entries are explicitly proxies: the required InP/GaAsP single-dot measurement remains absent.
