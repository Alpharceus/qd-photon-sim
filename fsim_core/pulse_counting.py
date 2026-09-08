"""Module D-FP -- finite-pulse photon counting for the exciton-biexciton
cascade under a rectangular pump waveform.

Why this module exists (peer-review-triage.md finding 1). The pulsed g2(0)
computed by loading.f1b_g2 / f8_g2 assumes a STATIC per-pulse loading
distribution: a mean number of excitons mu = r_dot * tau_pulse_ns is drawn
once per pulse (Poisson or moment-matched), and the cap-2 truncation (n>=2
relaxes to |XX>) stands in for whatever actually happens during the pulse.
That is exact only in the limit where the pulse is instantaneous compared to
every other rate in the problem. At the RT edge-emitter operating point the
escape rates k_X, k_XX (cw_g2.escape_rates_from_retention) and the ~100 ps
pump are FAST relative to the radiative lifetime, so a dot can be re-excited
and re-emit more than once within a single pulse -- repeated X emission, not
merely a fixed cascade photon count. This module instead propagates the
exact state-resolved factorial-moment hierarchy through the pump ("on") and
dark ("off") windows of one pulse period and reads off the per-period photon
statistics directly, retaining nonradiative (escape) transitions in the
generator without counting them as photons. [DR] moment hierarchy built on
the existing cw_g2 three-state generator; [A] rectangular pump waveform (the
diode's actual current pulse shape is not modeled here). Cite Hanschke et
al., npj Quantum Inf. 4, 43 (2018) for re-excitation during a finite pulse,
and Reischle et al., Optics Express 16, 12771 (2008) (DOI
10.1364/OE.16.012771) / Appl. Phys. Lett. 97, 143513 (2010) for the cascade
rate conventions (gamma_XX = 2 gamma_X, cap-2 ladder) this module shares
with cw_g2.py and loading.py. (Hanschke et al., "Origin of antibunching in
resonance fluorescence", is a different paper, PRL 125, 170402 (2020) -- do
not conflate the two.)

See also fsim_core.drive_mech.reexc_g2 for an independent re-excitation
model built on the same counting-moment idea (a closed linear ODE hierarchy
for the state-resolved photon-counting moments, no truncation) but NOT
interchangeable with this module: reexc_g2 has no escape channels (k_X =
k_XX = 0 always), counts the X line only (no XX cascade/filter split), and
has no gate (it integrates to a fixed t_end well past both lifetimes rather
than reading out a periodic steady state at a chosen window). The two are
cross-checked against each other at k=0/t_XX=0 in
verify/verify_pulse_counting.py.

Model. For probability vector p (basis |0>, |X>, |XX>, cw_g2.generator's
ordered basis) and state-resolved first/second factorial count moments m1,
m2 of the DETECTED (filtered) photon stream, propagate:

    dp/dt  = M p
    dm1/dt = M m1 + J p
    dm2/dt = M m2 + 2 J m1

M = cw_g2.generator(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio),
built once with the pump rate on (r_ns = the given rate) and once off
(r_ns = 0). J is the "counted jump" operator, same row/column convention as
M (rows = to-state, columns = from-state): J[0,1] = t_X * gamma_X_ns (an
X->G radiative jump, collected with filter transmission t_X) and
J[1,2] = t_XX * gamma_XX_ns (an XX->X radiative jump, collected with t_XX).
Stacking (p, m1, m2) into one 9-vector, the three equations above are one
linear, block lower-triangular system dv/dt = A v, so each window's
propagator is exact: v(t) = expm(A t) v(0) (scipy.linalg.expm).

One pulse period consists of a pump ("on") window of length tau_on_ns
followed by a dark ("off") window of length tau_dark_ns. The dot's
occupation carries over from one period to the next (a periodic steady
state), but the PHOTON COUNTS do not -- a detector/counter reads out once
per period and resets, so m1 = m2 = 0 at the start of every period, while p
is the fixed point of the one-period map Phi = expm(A_off_p tau_dark_ns)
expm(A_on_p tau_on_ns) restricted to its own p-block (m1, m2 evolve linearly
in p and do not feed back into it, so the periodic p-map is exactly the
3x3 map cw_g2.generator already gives). The fixed point is found by
iterating Phi from the ground state to a tolerance of 1e-12 (max 10000
iterations), matching the review's independent reproduction. The returned
g2 = sum(m2)/sum(m1)**2 and mean_counts = sum(m1) are therefore per-PULSE-
PERIOD quantities (the peak-area convention loading.f1b_g2 also uses, but
now including re-excitation instead of assuming an isolated instantaneous
load).

Counting gate (pr-pkg4-fix, gate-consistent counting). pulse_g2's optional
gate_ns restricts WHEN a jump counts, not the period being propagated: after
the periodic steady state p is found, the same one-period propagation runs
with the counting operator J active only for 0 <= t < gate_ns, t measured
from the PUMP ONSET (t=0 at the start of the on-window) -- outside the gate
the augmented block uses J=0 (a zero counting operator), so p keeps
evolving normally under M but m1/m2 stop accumulating. gate_ns=None (the
default) counts the whole period (gate_ns = tau_on_ns + tau_dark_ns) and is
bit-identical to the pre-gate propagation. A gate_ns at or beyond the period
also counts the whole period (a gate cannot outlast the one period this
call tracks). [A] recommended gate for a real detector: gate_ns =
tau_pulse_ns + 5*tau_rad_ns (catches >99% of a single-exponential decay
tail after the pump turns off).

Caller-side background window (pr-pkg4-fix2, item 1). This module has no
notion of an injection background channel -- that is a device.py-level
addition, folded into rho AFTER pulse_g2 returns mean_counts/mean_counts_x.
device.py integrates that background over min(gate_ns, tau_on_ns), the
overlap of the counting gate with the PUMP window, not the full gate: the
cascade's own photon-counting moments (m1, m2 above) correctly keep
accumulating radiative afterglow after the pump turns off (M's own decay
dynamics), but the transport background rate has no such tail modeled
anywhere in this module or its caller -- background afterglow past the
pulse is neglected [A]. This choice, not a wider or narrower window, is
what makes device.py's rho_pulsed reduce exactly to its CW rho in the
tau_dark_ns -> 0 limit (pr-pkg4-fix3 item 3, corrected by pr-pkg4-fix4
item 1: ONLY when the RESOLVED retention params device.py's evaluate()
actually used -- params["b0"] == params["beta"] == 0, exposed as
scalars["retention_params_used"] -- and the cavity is disabled, true of
every shipped card. This is NOT the same as device.py's raw
design.ret.b0/design.ret.beta fields: in ret.mode="proxy" (the default) an
unset (0.0) ret.b0/ret.beta resolves through the class-proxy Arrhenius fit
to a nonzero b0/beta, and in ret.mode="confinement" an explicit
ret.overrides entry can set a nonzero b0/beta on top of raw fields that
still read 0.0 -- either way the two rho definitions then differ by
construction, since the CW rho's own background/signal never carry a
b0/beta term or a cavity gain factor at all; see device.py's finite_pulse
block for the detail): both quantities then integrate signal and
background over the same (pump-only) window.

Physics honesty (pr-pkg4-fix3 item 4). Neglecting background afterglow
past the pulse end is one-sided and optimistic, never the reverse: it can
only omit background photons a real detector would still see, so
device.py's rho_pulsed is an UPPER BOUND on the true pulsed rho, not a
central estimate. At the FAVOURABLE corner (T_hs=230 K, gamma300=6,
delta_xx=8, NA=0.8, R_back=0.95, L=250 um) n_bg -- the COUNTED injection-window
background, not the neglected afterglow -- is already ~17.42% of the
counted background B_fp (n_bg/B_fp; pr-pkg4-fix4 item 5 correction: an
earlier version of this paragraph mislabeled n_bg/B_fp itself as "the
neglected background", which it is not -- n_bg IS counted, matching
device.py's own phrasing at device.py's DriveBlock.gate_ns docstring); a
~1 ns background carrier lifetime (afterglow decaying on that scale
instead of being cut off at the pulse end) would carry ~11x more
background photons into the counting window and move rho_pulsed from
0.858 to ~0.69, while g2_op moves only from 0.9817 to 0.98826 over the
same change [A] -- g2 depends on the cascade dynamics this module already
tracks exactly, not on the
signal-to-background ratio, so it is far less sensitive to the neglected
background model than rho is. At T_hs=300 K instead (same corner
otherwise) n_bg/B_fp = 0.117152, rho_pulsed = 0.8662227488, g2_op =
0.9976578221 -- quote both temperatures, they are not interchangeable.

Lemma 1 (collection-efficiency invariance): g2 is invariant under scaling
t_X and t_XX by a common factor, because J (and hence m1) scales linearly in
that factor while the 2*J*m1 forcing term in dm2/dt makes m2 scale
quadratically -- sum(m2)/sum(m1)**2 cancels the factor exactly. Checked in
verify/verify_pulse_counting.py.

Units: rates in 1/ns, times in ns (cw_g2.py convention).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from . import cw_g2

_TOL = 1e-12
_MAX_ITER = 10000


def _periodic_steady_state(M_on, M_off, tau_on_ns, tau_dark_ns):
    """Fixed point p of the one-period map Phi = expm(M_off tau_dark) @
    expm(M_on tau_on), found by iterating Phi from the ground state to
    tolerance _TOL (max _MAX_ITER iterations). Returns (p, converged)."""
    Phi = expm(M_off * tau_dark_ns) @ expm(M_on * tau_on_ns)
    p = np.array([1.0, 0.0, 0.0])
    converged = False
    for _ in range(_MAX_ITER):
        p_next = Phi @ p
        if np.max(np.abs(p_next - p)) < _TOL:
            p = p_next
            converged = True
            break
        p = p_next
    return p, converged


def _augmented(M, J):
    """9x9 block lower-triangular generator of the (p, m1, m2) stack for one
    window's constant M, J (see module docstring)."""
    Z = np.zeros((3, 3))
    return np.block([[M, Z, Z],
                     [J, M, Z],
                     [Z, 2.0 * J, M]])


def _propagate_period(M_on, M_off, J, p0, tau_on_ns, tau_dark_ns, gate_ns=None):
    """One period's (p, m1, m2) starting from p=p0, m1=m2=0, pump on for
    tau_on_ns then off for tau_dark_ns, with the counting operator J active
    only for 0 <= t < gate_ns (t from the pump onset, i.e. the start of the
    on-window) -- outside the gate the augmented block runs with J replaced
    by the zero matrix, so p keeps evolving under M while m1/m2 freeze.
    gate_ns=None (default) counts the whole period and is bit-identical to
    the original ungated propagation (two expm calls, on then off, both at
    full window length -- no zero-length windows are introduced). A gate_ns
    at or beyond the period is clamped to the period. Returns
    (p_period, m1, m2)."""
    period = tau_on_ns + tau_dark_ns
    gate = period if gate_ns is None else min(gate_ns, period)
    Z = np.zeros((3, 3))
    v = np.zeros(9)
    v[0:3] = p0
    on_count = min(gate, tau_on_ns)
    on_rest = tau_on_ns - on_count
    if on_count > 0:
        v = expm(_augmented(M_on, J) * on_count) @ v
    if on_rest > 0:
        v = expm(_augmented(M_on, Z) * on_rest) @ v
    dark_count = min(max(gate - tau_on_ns, 0.0), tau_dark_ns)
    dark_rest = tau_dark_ns - dark_count
    if dark_count > 0:
        v = expm(_augmented(M_off, J) * dark_count) @ v
    if dark_rest > 0:
        v = expm(_augmented(M_off, Z) * dark_rest) @ v
    return v[0:3], v[3:6], v[6:9]


def pulse_g2(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, t_X, t_XX,
            tau_on_ns, tau_dark_ns, pump_ratio=1.0, split=False,
            gate_ns=None) -> dict:
    """Per-pulse-period g2 and mean detected counts of the filtered X/XX
    cascade under a rectangular pump of duration tau_on_ns at rate r_ns,
    followed by a dark window of duration tau_dark_ns (see module
    docstring). gate_ns restricts counting to 0 <= t < gate_ns from the
    pump onset (module docstring "Counting gate"); gate_ns=None (default)
    counts the whole period and reproduces the pre-gate numbers bit-
    identically. split=True additionally decomposes the mean counts into
    the X-line-only and XX-line-only shares (finding 4's gated rho needs
    the X-line count alone, undiluted by the XX contribution) -- the
    XX-only share is recovered as mean_counts - mean_counts_x by linearity
    of the m1 ODE in J (p's own trajectory does not depend on J), so only
    one extra gated propagation is needed, not two."""
    M_on = cw_g2.generator(r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    M_off = cw_g2.generator(0.0, gamma_X_ns, gamma_XX_ns, k_X, k_XX, pump_ratio)
    p_ss, converged = _periodic_steady_state(M_on, M_off, tau_on_ns, tau_dark_ns)

    J = np.zeros((3, 3))
    J[0, 1] = t_X * gamma_X_ns
    J[1, 2] = t_XX * gamma_XX_ns
    p_period, m1, m2 = _propagate_period(M_on, M_off, J, p_ss, tau_on_ns, tau_dark_ns, gate_ns)

    mean_counts = float(np.sum(m1))
    if mean_counts < 1e-300:
        result = {"g2": float("nan"), "mean_counts": mean_counts,
                  "p_period": p_period, "converged": False}
    else:
        result = {"g2": float(np.sum(m2) / mean_counts**2), "mean_counts": mean_counts,
                  "p_period": p_period, "converged": bool(converged)}

    if split:
        J_X = np.zeros((3, 3)); J_X[0, 1] = J[0, 1]
        _, m1_x, _ = _propagate_period(M_on, M_off, J_X, p_ss, tau_on_ns, tau_dark_ns, gate_ns)
        result["mean_counts_x"] = float(np.sum(m1_x))
        result["mean_counts_xx"] = mean_counts - result["mean_counts_x"]

    return result
