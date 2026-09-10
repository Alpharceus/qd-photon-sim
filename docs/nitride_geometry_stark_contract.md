# Nitride geometry and Stark contract

These four opt-in cards extend the existing c-plane cards without changing their numerics. Nonpolar cards use `orientation=a_plane` and `polarization_factor=0` with the same H=3 nm, R=10 nm, x_in=0.25, Q=2000, 100 ps, 80 MHz, density, aperture, background, lifetime and cavity assumptions. QW cards use `geometry_type=qw_fluctuation`, `shape=disc`, H=3.5 nm, R=20 nm, w=3 nm, x_in=0.25, `d_i_nm=27`, `d_active_nm=3.5`, matched diode well thickness/composition, and `ret.channel=min`.

The QW reservoir is the GaN spectral reservoir plus the historical bulk-SRH layer mismatch for isolated-dot cards; QW cards share the surrounding-well material and thickness between confinement and diode background, while retaining the bulk-SRH carrier approximation. Same-composition thickness fluctuation confines laterally through the subband-energy difference. The SET island radius is independent of optical radius.

All geometry, orientation, screening and bias fields are [A] design choices unless explicitly marked [V]. Schade et al., phys. status solidi (b) (2011), DOI https://doi.org/10.1002/pssb.201046350, supports only an orientation/valence caveat; polarization factors here are [A], not measurements. Nonpolar strain, mass and valence ordering transfers and volume-mapped shapes have unknown systematic error; alternative-height sensitivity is planned.

The fixed-cavity-reference Stark convention anchors the cavity at one declared voltage (0 V for diagnostic traces), while dot energy, lifetime and overlap vary with bias. Source validity and spectral validity are separate; optical and hardware outcomes are reported separately. No current-derived screening law is inferred.

Evidence is qualified: Wang et al., Sci. Rep. 7, 12089 (2017), DOI https://arxiv.org/abs/1610.00152, reports uncapped AFM ~7 nm/~35 nm, 2.54 eV at 220 K, 19.0+/-0.4 meV linewidth, raw/corrected g2 0.47/0.21, 357+/-20 ps fast and ~4 ns slow decay, optical 76 MHz, 1 ps, 800 nm two-photon excitation. It is not electrical SPS or 300 K validation. Wang et al., APL 111, 053101 (2017), reports FSS 2-12 meV at 200 K (16 dots), omitted by the single-band model. Zhang et al., APL 108, 153102 (2016), Fig. 5, gives an optical-under-bias comparison of about -10 meV/V below 2 V; it is non-gating.

Piece 6 freezes the 230/250/273/300 K grid, screening 0/0.5/1, fixed cavity reference, supplementary budgets, diagnostics, manifest and VERDICT formats. No fabricated measurement interval is used; a user-supplied slope interval remains distinct from illustrative or Zhang intervals.
