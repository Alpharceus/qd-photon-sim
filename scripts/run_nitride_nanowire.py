"""Nanowire family/regime/strain-bound sweep (piece 9, formerly piece 8).

Independent predictions for the InGaN/GaN dot-in-nanowire tier
(``fsim_core/nitride_nanowire_device.py``): horizontal_as_built (Deshpande
et al. 2013/2014 as-built geometry) and vertical_photonic (designed HE11
wire), each evaluated under both drive regimes (rectangular 100 ps pulse,
deterministic_pair idealized SET loading) and both strain-bound scenarios
(unrelaxed conservative_lower, relaxed headline_upper), per
docs/nitride_nanowire_contract.md ("Sweep grid, VERDICT format, and output
paths"). Every plotted point and every VERDICT/BEST_PASSING_FLUX number is
an evaluate() row or a declared reduction of one (row_id traceable through
sweep.csv); nothing here invents a standalone physics result. Where this
spec text and the contract disagree, the contract wins (see this file's own
module docstring history and the STATUS decisions this piece reports).

Card selection is a pure filename lookup (house style, matching
scripts/run_nitride_geometry_stark.py's ``_card_name``): the four gated
headline cards (cards/nitride-nanowire-{horizontal,vertical}-{pulse,set}-
design.yaml) are loaded once and deep-copied per row so ``card_hash`` is
always the real sha256 of the file actually evaluated.

Output is a SINGLE sweep.csv (contract "Output paths (piece 9)": no
per-kind CSV split, unlike the geometry_stark piece) -- every row (core
grid, reduced cuts, planar-reference replays, and the Deshpande
comparisons) carries the same ``row_kind``/``row_id`` bookkeeping columns
plus whatever the row's own evaluate() call returned (nanowire and planar
rows return different column sets; the CSV writer takes the union and
blank-fills, exactly like csv.DictWriter's own default).
"""
from __future__ import annotations
import argparse, copy, hashlib, itertools, json, math, os, re, sys, time
import concurrent.futures as cf
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

import fsim_core.device as device  # noqa: E402

# --------------------------------------------------------------- grid axes
FAMILIES = ("horizontal_as_built", "vertical_photonic")
REGIMES = ("rectangular", "deterministic_pair")
STRAIN_BOUNDS = ("unrelaxed", "relaxed")
REP_RATES = (80.0e6, 200.0e6)
HEIGHT_NM = (1.5, 2.0, 3.0, 4.0)
X_IN = (0.25, 0.40)
T_HS = (230.0, 250.0, 273.0, 300.0)
CORE_R_NM = {
    "horizontal_as_built": (10.0, 12.5, 15.0, 20.0, 25.0, 40.0),
    "vertical_photonic": (60.0, 80.0, 100.0, 120.0),
}
# --quick: a declared deterministic subset (acceptance criterion 2), both
# families/regimes/bounds, T in {230,300} and both rates, never a
# full-coverage claim (complete=False always, enforced in main()).
QUICK_CORE_R_NM = {
    "horizontal_as_built": (12.5, 25.0),
    "vertical_photonic": (80.0, 120.0),
}
QUICK_HEIGHT_NM = (2.0, 4.0)
QUICK_X_IN = (0.25, 0.40)
QUICK_T_HS = (230.0, 300.0)

G2, FLUX_FLOOR = 0.5, 1000.0

# Reduced-cut evaluation is always at T in {230,300} and both rates (both
# regimes/bounds), per the contract's "Sweep grid" section -- independent of
# --quick, which only shrinks the MAIN grid (contract: "Reduced cuts ...
# T=230/300 K and both rates").
CUT_T_HS = (230.0, 300.0)

CARD_NAME = {
    ("horizontal_as_built", "rectangular"): "nitride-nanowire-horizontal-pulse-design.yaml",
    ("horizontal_as_built", "deterministic_pair"): "nitride-nanowire-horizontal-set-design.yaml",
    ("vertical_photonic", "rectangular"): "nitride-nanowire-vertical-pulse-design.yaml",
    ("vertical_photonic", "deterministic_pair"): "nitride-nanowire-vertical-set-design.yaml",
}
PLANAR_CARD_NAME = {
    ("c_plane", "rectangular"): "nitride-cavity-pulse-design.yaml",
    ("c_plane", "deterministic_pair"): "nitride-cavity-set-design.yaml",
    ("nonpolar", "rectangular"): "nitride-nonpolar-pulse-design.yaml",
    ("nonpolar", "deterministic_pair"): "nitride-nonpolar-set-design.yaml",
}
DESHPANDE2014_CARD = "nitride-deshpande2014-comparison-design.yaml"

# Reference geometry for reduced cuts and the Deshpande 2014-like comparison
# (contract Card schema defaults: 12.5 nm horizontal disc==core, 80 nm
# vertical core with a 12.5 nm disc, both at height 2 nm, x_in 0.40).
REF_GEOM = {
    "horizontal_as_built": {"core_radius_nm": 12.5, "height_nm": 2.0, "x_in": 0.40},
    "vertical_photonic": {"core_radius_nm": 80.0, "height_nm": 2.0, "x_in": 0.40},
}

# ------------------------------------------------------------------ utils
class _Tee:
    """H1 fix: the runner writes its OWN log file INSIDE the run directory
    (out/nitride_nanowire/<dir>/run.log) so manifest.json's output_hashes
    (which hashes every file actually sitting in the resolved run dir, see
    main()) always has a real, in-directory file to hash -- never a log the
    orchestrator's own shell redirected somewhere OUTSIDE the run dir (the
    prior failure: the artifact-mode verifier's hash check failed 171/172
    because a "full-rerun.log" hash was recorded for a file that did not
    exist inside out/nitride_nanowire/full/)."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
        return len(data)

    def flush(self):
        for s in self.streams:
            s.flush()


def _safe(out):
    p = Path(out).resolve()
    base = (ROOT / "out" / "nitride_nanowire").resolve()
    if p != base and base not in p.parents:
        raise ValueError("out-dir must stay under out/nitride_nanowire")
    return p


def _primitive(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    return json.dumps(v, sort_keys=True, default=str)


def _finite(x):
    try:
        return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)
    except TypeError:
        return False


def _close(a, b, rtol=1e-6, atol=1e-9):
    fa, fb = _finite(a), _finite(b)
    if not fa and not fb:
        return True
    if fa != fb:
        return False
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


# Bookkeeping-only keys on a row's parameter dict (never a physics input):
# excluding them from the cache/dedup identity (fix round: "make the cache
# identity physics-only so duplicate evaluate() inputs deduplicate") lets a
# reduced-cut row whose sampled value happens to equal the main-grid default
# (e.g. the screening_fraction=0.0 cut, which duplicates 32 core-row inputs)
# collapse onto the SAME real evaluate() call as its core-row twin, instead
# of hashing separately just because its sensitivity_axis/sensitivity_value
# labels differ. The row's OWN CSV content still carries its real
# sensitivity_axis/sensitivity_value (row.update(p) uses the untouched p);
# only the identity/cache KEY drops them.
_NONPHYSICS_KEYS = ("sensitivity_axis", "sensitivity_value")


def _physics_only(p):
    return {k: v for k, v in p.items() if k not in _NONPHYSICS_KEYS}


_CARD_CACHE = {}
_CARD_HASH_CACHE = {}


def _card_hash(name):
    if name not in _CARD_HASH_CACHE:
        _CARD_HASH_CACHE[name] = hashlib.sha256((ROOT / "cards" / name).read_bytes()).hexdigest()
    return _CARD_HASH_CACHE[name]


def _load_card(name):
    if name not in _CARD_CACHE:
        _CARD_CACHE[name] = device.DeviceDesign.load(ROOT / "cards" / name)
    return copy.deepcopy(_CARD_CACHE[name])


# ------------------------------------------------------ parallel evaluation
# fsim_core.nitride_nanowire_device.evaluate_nanowire's own injector_
# feasibility resonance search (fsim_core/nitride_nanowire_injector.py's
# _find_resonance/_coarse_scan/T_at, called on EVERY row regardless of
# regime -- the device module always prices the RT screen, never gated on
# cycle_loading) measured at 6-13 s per evaluate() call on this machine,
# independent of anything this sweep script controls (a dependency-module
# performance property, not a bug this piece may patch -- "Dependency
# files may be read but not modified"). At that per-call cost the mandated
# 2560-row core grid cannot finish serially inside 1800 s (measured: a
# single evaluate() call 6.1-9.8 s; a 436-call quick run took 1587 s
# serially). Evaluate calls are independent, side-effect-free (each builds
# its own DeviceDesign from a plain parameter dict/tuple and returns plain
# scalars), so they are dispatched across a process pool -- same evaluate()
# calls, same inputs, same outputs, only wall-clock distributed across
# cores; this changes nothing about which rows are evaluated or what they
# return, only how fast. Benchmarked throughput on this machine: ~2.0
# calls/s for small-core-radius rows, ~1.3 calls/s for large-radius rows,
# at N_WORKERS=os.cpu_count().
N_WORKERS = max(1, os.cpu_count() or 4)


def _evaluate_job(job):
    """Top-level, picklable worker entry point for ProcessPoolExecutor.
    `job` is (kind, payload): kind="nanowire" -> payload is a full_defaults()
    parameter dict, built into a design via _design(); kind="planar" ->
    payload is the (family_tag, regime, T_hs, screening_fraction) tuple for
    _planar_design(); kind="deshpande2014" -> payload is ignored, the
    untouched old planar 2014 card is loaded directly. Runs in a worker
    process; identity/cache/row-assembly bookkeeping stays in the parent
    (see _row/_deshpande2014_replay), so this function does exactly what
    the serial code path did per row, just off the main process."""
    kind, payload = job
    if kind == "nanowire":
        d, _card = _design(payload)
    elif kind == "planar":
        d, _card = _planar_design(*payload)
    else:
        d = device.DeviceDesign.load(ROOT / "cards" / DESHPANDE2014_CARD)
    return device.evaluate(d)["scalars"]


def _job_identity(job):
    """Same identity string _row()/_deshpande2014_replay() used to key their
    in-process cache, computed up front so identical rows (a literal
    duplicate parameter dict or planar tuple) are deduplicated to ONE real
    evaluate() call before dispatch, exactly as the serial cache did."""
    kind, payload = job
    if kind == "nanowire":
        return json.dumps({"planar": None, **_physics_only(payload)}, sort_keys=True, separators=(",", ":"), default=str)
    if kind == "planar":
        return json.dumps({"planar": payload}, sort_keys=True, separators=(",", ":"), default=str)
    return json.dumps({"deshpande2014_replay": True}, sort_keys=True)


def _evaluate_all(jobs, n_workers=N_WORKERS):
    """Deduplicate `jobs` (list of (kind, payload)) by identity, evaluate
    every unique one across a process pool, and return {identity: scalars}.
    Prints periodic progress to stderr (this run takes many minutes; a
    silent foreground process gives no feedback otherwise)."""
    unique = {}
    for job in jobs:
        ident = _job_identity(job)
        if ident not in unique:
            unique[ident] = job
    results = {}
    idents = list(unique.keys())
    total = len(idents)
    if total == 0:
        return results
    t0 = time.time()
    done = 0
    with cf.ProcessPoolExecutor(max_workers=n_workers) as ex:
        fut_to_ident = {ex.submit(_evaluate_job, unique[ident]): ident for ident in idents}
        for fut in cf.as_completed(fut_to_ident):
            ident = fut_to_ident[fut]
            results[ident] = fut.result()
            done += 1
            if done % 100 == 0 or done == total:
                print("progress evaluate_calls=%d/%d elapsed_s=%.1f workers=%d" % (
                    done, total, time.time() - t0, n_workers), file=sys.stderr, flush=True)
    return results


# ------------------------------------------------------- default parameter set
def full_defaults(family):
    """Every field this sweep can vary, at its card-default (main-grid)
    value, for ONE nanowire family. Every row (core or reduced-cut) starts
    from this dict so that (a) two rows with identical physical inputs
    share one cache slot regardless of which builder produced them, and
    (b) sweep.csv carries a stable, complete column set for every nanowire
    row (house style: scripts/run_nitride_geometry_stark.py's `_base`)."""
    ref = REF_GEOM[family]
    return dict(
        family=family, regime="rectangular",
        core_radius_nm=ref["core_radius_nm"], height_nm=ref["height_nm"], x_in=ref["x_in"],
        T_hs=300.0, strain_bound="relaxed", rep_rate_hz=200.0e6,
        screening_fraction=0.0, shell="none",
        S_cm_s=1.0e3, reservoir_access=1.0, occupied_dot_access=0.05,
        NA=0.5, bottom_reflectivity=0.0, dipole_weights=None,
        gamma300=35.0, tau_rad0_ns=1.0, tau_cap_ps=10.0,
        R_s_ohm=(2.38e9 if family == "horizontal_as_built" else 1.0e6),
        Rth_K_W=(1.0e9 if family == "horizontal_as_built" else 1.0e7),
        C_parasitic_F=0.0,
        al_fraction=0.30, electron_barrier_thickness_nm=2.0, hole_barrier_thickness_nm=2.0,
        growth_tolerance_steps=1.0, alignment_uncertainty_meV=15.0,
        occupancy_control_known=False, second_pair_control_known=False,
        I_uA=0.002, tau_pulse_ns=0.1, b_res=0.1,
        sensitivity_axis="", sensitivity_value="",
    )


def _design(p):
    """Build a nanowire DeviceDesign from a full parameter dict `p` (see
    full_defaults). Sweep coupling (DEVICE-PIECE CONSTRAINT 1 / spec
    "Sweep coupled dimensions together"): core/outer/conducting radii
    always follow the unshelled core; horizontal's disc/SET-island radius
    follows the core too; vertical's disc/SET-island radius stays at its
    own fixed value (12.5 nm default) independent of the core sweep."""
    family, regime = p["family"], p["regime"]
    d = _load_card(CARD_NAME[(family, regime)])
    nw, dot = d.nitride["nanowire"], d.nitride["dot"]
    core = float(p["core_radius_nm"])
    shell = p.get("shell", "none")
    nw["family"] = family
    nw["core_radius_nm"] = core
    nw["shell"] = shell
    nw["outer_radius_nm"] = core if shell == "none" else core + 3.0
    nw["strain_bound"] = p["strain_bound"]
    dot["strain_fraction"] = 1.0 if p["strain_bound"] == "unrelaxed" else 0.0
    dot["height_nm"] = float(p["height_nm"])
    dot["x_in"] = float(p["x_in"])
    dot["screening_fraction"] = float(p.get("screening_fraction", 0.0))
    d.drive.diode["conducting_radius_nm"] = core
    if family == "horizontal_as_built":
        disc = core
    else:
        disc = float(dot.get("radius_nm", 12.5))  # vertical: fixed disc, never the core (H6)
    dot["radius_nm"] = disc
    d.drive.set_params["radius_nm"] = disc
    d.drive.set_params.pop("C_sigma_F", None)

    surf = d.nitride["surface"]
    surf["shell"] = shell
    surf["shell_multiplier"] = 1.0 if shell == "none" else 0.1
    surf["S_cm_s"] = float(p.get("S_cm_s", 1.0e3))
    surf["reservoir_access"] = float(p.get("reservoir_access", 1.0))
    surf["occupied_dot_access"] = float(p.get("occupied_dot_access", 0.05))

    ph = d.nitride["photonics"]
    if family == "horizontal_as_built":
        ph["NA"] = float(p.get("NA", 0.5))
    else:
        ph["bottom_reflectivity"] = float(p.get("bottom_reflectivity", 0.0))
    if p.get("dipole_weights") is not None:
        ph["dipole_weights"] = tuple(p["dipole_weights"])

    wt = d.nitride["wire_thermal"]
    wt["R_s_ohm"] = float(p.get("R_s_ohm", wt["R_s_ohm"]))
    wt["Rth_K_W"] = float(p.get("Rth_K_W", wt["Rth_K_W"]))
    wt["C_parasitic_F"] = float(p.get("C_parasitic_F", 0.0))

    d.nitride["tau_rad0_ns"] = float(p.get("tau_rad0_ns", 1.0))
    d.nitride["tau_cap_ps"] = float(p.get("tau_cap_ps", 10.0))
    d.dot.gamma300 = float(p.get("gamma300", 35.0))

    inj = d.nitride["injector"]
    inj["al_fraction"] = float(p.get("al_fraction", 0.30))
    inj["electron_barrier_thickness_nm"] = float(p.get("electron_barrier_thickness_nm", 2.0))
    inj["hole_barrier_thickness_nm"] = float(p.get("hole_barrier_thickness_nm", 2.0))
    inj["growth_tolerance_steps"] = float(p.get("growth_tolerance_steps", 1.0))
    inj["alignment_uncertainty_meV"] = float(p.get("alignment_uncertainty_meV", 15.0))
    inj["occupancy_control_known"] = bool(p.get("occupancy_control_known", False))
    inj["second_pair_control_known"] = bool(p.get("second_pair_control_known", False))

    tau_pulse_ns = float(p.get("tau_pulse_ns", 0.1))
    rep = float(p["rep_rate_hz"])
    d.drive.diode["tau_pulse_ns"] = tau_pulse_ns
    d.drive.rep_rate_hz = rep
    d.drive.duty = tau_pulse_ns * 1e-9 * rep  # spec: "duty=tau*rep"
    d.drive.cycle_loading = regime
    d.drive.I_uA = float(p.get("I_uA", 0.002))
    d.drive.b_res = float(p.get("b_res", 0.1))
    d.thermal.T_hs = float(p["T_hs"])
    return d, CARD_NAME[(family, regime)]


def _planar_design(family_tag, regime, T_hs, screening_fraction):
    name = PLANAR_CARD_NAME[(family_tag, regime)]
    d = _load_card(name)
    d.nitride["dot"]["screening_fraction"] = float(screening_fraction)
    d.thermal.T_hs = float(T_hs)
    return d, name


# ------------------------------------------------------------------ rows
def _row(rid, kind, p, seen, results, planar=None):
    """planar=None -> nanowire row from `p` via _design(); planar=(family_tag,
    regime,T_hs,screening) -> a planar-platform reference replay instead
    (row_kind carries the distinction; the CSV column union covers both
    schemas). `results` is the {identity: scalars} map _evaluate_all()
    already computed (in a process pool, see above); `seen` is a mutable
    set this function uses to reproduce the original in-process cache's
    `cache_hit` semantics (True iff an earlier row in THIS run already used
    the identical identity) without re-calling evaluate() -- the real
    evaluate() call already happened once per unique identity, up front.
    The identity key is PHYSICS-ONLY (see `_physics_only`); the row's own
    CSV content still carries the untouched `p` (sensitivity_axis/
    sensitivity_value included) via `row.update(p)` below."""
    ident = json.dumps({"planar": planar, **_physics_only(p)} if planar is None else {"planar": planar},
                        sort_keys=True, separators=(",", ":"), default=str)
    hit = ident in seen
    seen.add(ident)
    s = results[ident]
    card = CARD_NAME[(p["family"], p["regime"])] if planar is None else PLANAR_CARD_NAME[(planar[0], planar[1])]
    row = {"row_id": rid, "row_kind": kind, "card_file": card, "card_hash": _card_hash(card),
           "cache_hit": hit, "cache_identity": ident}
    if planar is None:
        row.update(p)
    else:
        row.update({"family": planar[0], "regime": planar[1], "T_hs": planar[2], "screening_fraction": planar[3]})
    for k, v in s.items():
        row[k] = _primitive(v)
    row["valid"] = bool(s.get("valid"))
    row["invalid_reasons"] = json.dumps(s.get("invalid_reasons", []))
    row["optical_pass"] = bool(s.get("optical_pass", False))
    row["hardware_qualified"] = bool(s.get("hardware_qualified", False))
    row["rti_qualified"] = bool(s.get("rti_qualified", False))
    row["headline_eligible"] = bool(s.get("headline_eligible", False))
    flux = s.get("collected_flux_pulsed_s")
    row["eligible"] = bool(row["valid"] and _finite(s.get("g2_op")) and _finite(flux) and flux >= FLUX_FLOOR)
    _attach_quality_columns(row, s)
    return row


# H3 gate-anti-monotonicity fix: one_pair_valid (blocked_load_probability<=
# 1e-9) can be satisfied by ADDING loss (faster occupied-dot emptying) as
# readily as by improving device quality, so optical_pass alone is not a
# throughput signal. photons_per_cycle/emission_probability_per_cycle are
# script-computed columns (never re-deriving physics, only dividing/summing
# evaluate()'s own reported outputs) and quality_pass is a script-computed
# gate ORTHOGONAL to the contract's own optical_pass definition (no fsim_core
# or contract change): quality_pass=optical_pass AND photons_per_cycle>=0.01
# [A, orchestrator threshold: one collected photon per hundred cycles].
QUALITY_PHOTONS_PER_CYCLE_FLOOR = 0.01


def _attach_quality_columns(row, s):
    flux = s.get("collected_flux_pulsed_s")
    rep_hz = s.get("rep_rate_hz", row.get("rep_rate_hz"))
    if _finite(flux) and _finite(rep_hz) and float(rep_hz) > 0:
        ppc = float(flux) / float(rep_hz)
    else:
        ppc = float("nan")
    row["photons_per_cycle"] = ppc
    mcx, mcxx = s.get("mean_counts_x"), s.get("mean_counts_xx")
    if _finite(mcx) and _finite(mcxx):
        row["emission_probability_per_cycle"] = min(1.0, float(mcx) + float(mcxx))
    else:
        row["emission_probability_per_cycle"] = float("nan")
    row["quality_pass"] = bool(row.get("optical_pass") and _finite(ppc) and ppc >= QUALITY_PHOTONS_PER_CYCLE_FLOOR)
    return row


def _deshpande2014_replay(rid, seen, results):
    """Untouched OLD planar 2014 comparison card, replayed unmodified
    alongside the new horizontal design rows (spec: "Include the untouched
    old planar 2014 comparison alongside new horizontal design rows")."""
    ident = json.dumps({"deshpande2014_replay": True}, sort_keys=True)
    hit = ident in seen
    seen.add(ident)
    s = results[ident]
    row = {"row_id": rid, "row_kind": "planar_2014_replay", "card_file": DESHPANDE2014_CARD,
           "card_hash": _card_hash(DESHPANDE2014_CARD), "cache_hit": hit, "cache_identity": ident,
           "measured_g2": 0.29, "measured_lifetime_ns": 1.3, "measured_lambda_nm": 630.0}
    for k, v in s.items():
        row[k] = _primitive(v)
    row["valid"] = bool(s.get("valid"))
    row["invalid_reasons"] = json.dumps(s.get("invalid_reasons", []))
    row["optical_pass"] = bool(s.get("optical_pass", False))
    _attach_quality_columns(row, s)
    return row


# ------------------------------------------------------------- combo builders
def build_core(family, quick):
    rs = QUICK_CORE_R_NM[family] if quick else CORE_R_NM[family]
    hs = QUICK_HEIGHT_NM if quick else HEIGHT_NM
    xs = QUICK_X_IN if quick else X_IN
    ts = QUICK_T_HS if quick else T_HS
    out = []
    for r, h, x, t, reg, sb, rate in itertools.product(rs, hs, xs, ts, REGIMES, STRAIN_BOUNDS, REP_RATES):
        p = full_defaults(family)
        p.update(core_radius_nm=r, height_nm=h, x_in=x, T_hs=t, regime=reg,
                  strain_bound=sb, rep_rate_hz=rate)
        out.append(p)
    return out


def build_reduced_cuts(quick):
    """One-at-a-time reduced cuts at each family's reference geometry, both
    regimes/bounds, T in {230,300}, both rates (contract "Sweep grid").
    `quick` trims the per-axis value list (never the two families/regimes/
    bounds/rates/temperatures) so --quick stays well under 600 s."""
    out = []

    def axis(name, key, values):
        for family in FAMILIES:
            for reg, sb, t, rate in itertools.product(REGIMES, STRAIN_BOUNDS, CUT_T_HS, REP_RATES):
                for v in values:
                    p = full_defaults(family)
                    p.update(regime=reg, strain_bound=sb, T_hs=t, rep_rate_hz=rate)
                    p[key] = v
                    p["sensitivity_axis"] = name
                    p["sensitivity_value"] = v
                    out.append(p)

    if quick:
        axis("screening_fraction", "screening_fraction", (1.0,))
        axis("b_res", "b_res", (0.02,))
        axis("occupied_dot_access", "occupied_dot_access", (1.0,))
        axis("R_s_ohm", "R_s_ohm", (1.0e6,))  # horizontal-only meaningful value, evaluated on both (vertical is its own default)
        axis("current_uA", "I_uA", (0.02,))
        return out

    # REDUCED-CUT BUDGET TRIM (decisions: see main()'s `decisions` list and
    # results.md "Reduced-cut coverage" section). Each real evaluate() call
    # measured 6-13 s on this machine (dominated by fsim_core.nitride_
    # nanowire_injector.injector_feasibility's resonance search, run on
    # EVERY row regardless of regime -- a dependency-module cost this piece
    # cannot patch). At that cost the mandated 2560-row core grid alone
    # consumes essentially the whole 1800 s/40-minute operational budget
    # even parallelized across all cores (see N_WORKERS above), so the
    # full-mode reduced-cut axis LIST is trimmed to the orchestrator's
    # explicit reduced-cut deliverable list (screening {0,1};
    # occupied_dot_access's 1.0 conservative partner; S {1e2,1e4}; shell
    # AlGaN; the injector cuts named in the spec text -- barrier thickness,
    # Al fraction, alignment uncertainty, growth tolerance, occupation-
    # control uncertainty; R_s designed 1e6 on the horizontal family only)
    # plus one b_res=0.02 point kept for the contract's own opening-
    # paragraph emphasis on that sensitivity. Each kept axis samples only
    # the ALTERNATIVE value(s), never re-testing a value that already
    # equals the main-grid default (redundant with the core grid).
    # RESTORED this fix round (see the current_pulse_width axis below):
    # current_uA/tau_pulse_ns pulse sensitivity, previously dropped.
    # DROPPED entirely this run (never Cartesian-producted against the
    # kept axes, simply not sampled): reservoir_access, gamma300
    # (linewidth), tau_rad0_ns, tau_cap_ps (capture), C_parasitic_F,
    # Rth_K_W (both families), R_s_ohm on the
    # vertical family, NA (horizontal)/bottom_reflectivity (vertical)
    # collection-envelope cuts. This is a real, reported coverage gap
    # against the contract's fuller "Sweep grid" reduced-cut list, not a
    # silent omission -- see results.md.
    axis("screening_fraction", "screening_fraction", (0.0, 1.0))
    axis("b_res", "b_res", (0.02,))
    axis("S_cm_s", "S_cm_s", (1.0e2, 1.0e4))
    axis("shell", "shell", ("AlGaN",))
    axis("occupied_dot_access", "occupied_dot_access", (1.0,))
    axis("al_fraction", "al_fraction", (0.2,))
    axis("growth_tolerance_steps", "growth_tolerance_steps", (2.0,))
    axis("alignment_uncertainty_meV", "alignment_uncertainty_meV", (30.0,))
    axis("occupation_control_uncertainty", "occupancy_control_known", (True,))

    # injector_barrier_thickness_nm: electron AND hole barrier thickness
    # move together (one "barrier thickness" cut, contract's "injector
    # barrier thickness" reduced axis); one alternative value (thinner).
    for family in FAMILIES:
        for reg, sb, t, rate in itertools.product(REGIMES, STRAIN_BOUNDS, CUT_T_HS, REP_RATES):
            p = full_defaults(family)
            p.update(regime=reg, strain_bound=sb, T_hs=t, rep_rate_hz=rate,
                      electron_barrier_thickness_nm=1.0, hole_barrier_thickness_nm=1.0)
            p["sensitivity_axis"] = "injector_barrier_thickness_nm"
            p["sensitivity_value"] = 1.0
            out.append(p)

    # R_s designed-contact alternative, horizontal family ONLY (spec: "R_s
    # designed 1e6 on the horizontal family"); the vertical family's own
    # default already IS 1e6, so no vertical R_s cut is needed here.
    for reg, sb, t, rate in itertools.product(REGIMES, STRAIN_BOUNDS, CUT_T_HS, REP_RATES):
        p = full_defaults("horizontal_as_built")
        p.update(regime=reg, strain_bound=sb, T_hs=t, rep_rate_hz=rate, R_s_ohm=1.0e6)
        p["sensitivity_axis"] = "R_s_ohm"; p["sensitivity_value"] = 1.0e6
        out.append(p)

    # RESTORED (fix round, H-block "Previous attempt failed because" M9-M10:
    # "restore the current / pulse-width cut"): I in {0.001,0.002,0.02} uA x
    # tau_pulse in {0.01,0.1,1} ns, Cartesian (3x3=9 rows), at the horizontal
    # reference geometry only, rectangular regime (the 100 ps headline lives
    # on this regime), relaxed headline strain bound, T_hs=300 K, the
    # horizontal card's own 200 MHz default rate -- "preserving the 100 ps
    # headline" (one of the 9 combos, I=0.002/tau=0.1, IS the headline
    # point itself, kept in the declared cut list rather than skipped, per
    # spec "as separate cuts"; the physics-only cache identity (see
    # _physics_only) collapses it onto whichever existing row already has
    # identical physics inputs, so it costs zero EXTRA real evaluate() calls
    # even though it is declared/counted as one of the 9). duty=tau*rep and
    # the pair supply move together automatically inside _design() (spec:
    # "Coupled duty and pair supply change consistently").
    for i_ua in (0.001, 0.002, 0.02):
        for tau_ns in (0.01, 0.1, 1.0):
            p = full_defaults("horizontal_as_built")
            p.update(regime="rectangular", strain_bound="relaxed", T_hs=300.0, rep_rate_hz=200.0e6,
                      I_uA=i_ua, tau_pulse_ns=tau_ns)
            p["sensitivity_axis"] = "current_pulse_width"
            p["sensitivity_value"] = f"I={i_ua}uA_tau={tau_ns}ns"
            out.append(p)
    return out


DROPPED_CUT_AXES = (
    "reservoir_access", "gamma300 (linewidth)", "tau_rad0_ns", "tau_cap_ps (capture)",
    "C_parasitic_F", "Rth_K_W (both families)",
    "R_s_ohm (vertical family)", "NA (horizontal)/bottom_reflectivity (vertical) collection envelope",
)


def build_dipole_falsification():
    """A single, cheap demonstration that CPLANE_ONLY_DIPOLE_WEIGHTS
    predicts the wrong-sign degree_of_linear_polarization against the
    deshpande2013_polarization +70% anchor (contract "Composition rules
    for the device piece" bullet 9); non-gating, results-text only."""
    out = []
    for weights, tag in (((1.0 / 3, 1.0 / 3, 1.0 / 3), "isotropic_default"),
                          ((0.0, 0.5, 0.5), "cplane_only_falsified")):
        p = full_defaults("horizontal_as_built")
        p["dipole_weights"] = weights
        p["sensitivity_axis"] = "dipole_weights"
        p["sensitivity_value"] = tag
        out.append(p)
    return out


def build_deshpande2013(quick):
    """2013 diagnostic rows at T_hs=10 K (outside the main T_hs grid):
    x_in=0.25, H=2nm, R in {12.5,15}, I in {1,2} nA, relaxed strain bound
    (2013 relaxed-matched anchor). CW electrical in the paper; this
    evaluator has no CW/HBT path, so these are published as the pulsed
    model's drive_mismatch comparison (contract: "publish the pulsed model
    comparison as drive_mismatch ... never substitute pulse_g2 as
    predicted CW g2")."""
    radii = (12.5,) if quick else (12.5, 15.0)
    currents = (0.001,) if quick else (0.001, 0.002)
    out = []
    for r in radii:
        for i_ua in currents:
            p = full_defaults("horizontal_as_built")
            p.update(core_radius_nm=r, height_nm=2.0, x_in=0.25, T_hs=10.0,
                      regime="rectangular", strain_bound="relaxed", rep_rate_hz=80.0e6,
                      I_uA=i_ua)
            p["sensitivity_axis"] = "deshpande2013_comparison"
            p["sensitivity_value"] = f"R={r}nm_I={i_ua}uA"
            out.append(p)
    return out


def build_planar_reference(quick):
    out = []
    ts = (230.0, 300.0)
    for family_tag in ("c_plane", "nonpolar"):
        for reg in REGIMES:
            for t in ts:
                for s in (0.0, 1.0):
                    out.append((family_tag, reg, t, s))
    return out


def planned_counts(quick):
    core = sum(len(build_core(f, quick)) for f in FAMILIES)
    cuts = len(build_reduced_cuts(quick)) + len(build_dipole_falsification())
    desh2013 = len(build_deshpande2013(quick))
    planar = len(build_planar_reference(quick))
    total = core + cuts + desh2013 + planar + 1  # +1: deshpande2014 replay
    return {"core_rows": core, "reduced_cut_rows": cuts, "deshpande2013_rows": desh2013,
            "planar_reference_rows": planar, "deshpande2014_replay_rows": 1,
            "planned_evaluate_calls_upper_bound": total}


# --------------------------------------------------------------------- csv
def _write_csv(path, rows):
    keys = sorted({k for r in rows for k in r})
    import csv
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)


# --report-only fix (sweep fix 1): the closed set of row columns this
# script's OWN downstream logic (never fsim_core, which only ever sees
# csv-string columns through THIS module's row dicts) compares against the
# literal Python singleton True/False -- `_above_floor`'s `valid is True`,
# results.md's `set_feasible is True`/`rti_feasible is True`, and every
# plain-truthiness `if row.get("optical_pass")` / `if not row.get("valid")`
# check throughout _results_md/verdict_lines/the plot functions. Every
# OTHER csv column, in particular the `bound_reversal_pair`/`bound_reversal`
# STRING sentinels ("True"/"False"/"not_comparable"/"not_computed",
# compared with `== "True"` throughout -- NEVER coerced to bool, or those
# comparisons silently go permanently False), is left as the raw CSV string
# or opportunistically float-converted (every numeric column is already
# consumed via float()/_finite_num() at its use sites, which are tolerant
# of strings too -- this is a belt-and-suspenders exact-type match, not a
# functional requirement for the numeric columns).
_CSV_BOOL_COLUMNS = frozenset((
    "valid", "optical_pass", "hardware_qualified", "rti_qualified", "headline_eligible",
    "eligible", "quality_pass", "cache_hit", "single_mode", "one_pair_valid",
    "pair_supply_possible", "rti_feasible", "rti_transport_feasible", "set_feasible",
    "thermal_converged", "device_pass", "rti_device_pass",
))


def _coerce_csv_row(row):
    out = {}
    for k, v in row.items():
        if v == "":
            out[k] = None
        elif k in _CSV_BOOL_COLUMNS:
            out[k] = (v == "True") if isinstance(v, str) else bool(v)
        else:
            try:
                out[k] = float(v)
            except (TypeError, ValueError):
                out[k] = v
    return out


def _read_csv_rows(path):
    import csv
    with path.open("r", newline="", encoding="utf-8") as f:
        return [_coerce_csv_row(r) for r in csv.DictReader(f)]


# ---------------------------------------------------------------- grouping
def _group(rows, **filters):
    out = []
    for r in rows:
        ok = True
        for k, v in filters.items():
            rv = r.get(k)
            if isinstance(v, float) or isinstance(rv, float):
                try:
                    ok = ok and math.isclose(float(rv), float(v), rel_tol=1e-6, abs_tol=1e-9)
                except (TypeError, ValueError):
                    ok = False
            else:
                ok = ok and (rv == v)
            if not ok:
                break
        if ok:
            out.append(r)
    return out


def _bound_partner(row, core_rows):
    other = "relaxed" if row.get("strain_bound") == "unrelaxed" else "unrelaxed"
    matches = _group(core_rows, family=row.get("family"), regime=row.get("regime"),
                      rep_rate_hz=row.get("rep_rate_hz"), core_radius_nm=row.get("core_radius_nm"),
                      height_nm=row.get("height_nm"), x_in=row.get("x_in"), T_hs=row.get("T_hs"),
                      strain_bound=other)
    return matches[0] if matches else None


# -------------------------------------------------------------------- plots
def _finite_num(v):
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _leg(ax, **kw):
    return ax.legend(fontsize=7, loc="best", **kw)


BOUND_STYLE = {"unrelaxed": ("--s", "unrelaxed (conservative_lower)"), "relaxed": ("-o", "relaxed (headline_upper)")}

# M7 fix: the SAME one-line qualification is now printed on EVERY exported
# figure (previously distributed piecemeal across only some of the six
# PNGs) -- results obligation bullet 11's RC caveat / access-1.0 lifetime
# cap / opposite-endpoint anchor match / E_C/kT charging wall, all in one
# sentence. Any figure-specific caption is APPENDED to it, never replaces
# it. The exact string used is also recorded into each figure's own
# plot_row_mapping contract dict under "__caption__" so the artifact
# verifier can assert every figure actually carries it (never trusting the
# PNG's own pixels for this).
COMMON_FIGURE_QUALIFICATION = (
    "Qualification (applies to every figure this run): RC caveat -- the as-built horizontal contact cannot "
    "deliver a 100 ps step (see 'RC caveat'); access cap -- occupied_dot_access=1.0 caps tau_X/tau_XX (see "
    "'Access-1.0 lifetime cap'); opposite endpoints -- the 2013/2014 strain-bound anchors are matched by "
    "OPPOSITE bounds, never averaged; charging wall -- deterministic loading at 230-300K fails "
    "set_EC_over_kT for every core_radius_nm>=10 nm priced here (see 'E_C/kT wall')."
)


def _figure_footer(caption):
    return f"{caption} {COMMON_FIGURE_QUALIFICATION}" if caption else COMMON_FIGURE_QUALIFICATION


def plot_vs_axis(out, name, core, x_key, y_keys, title, fixed, group_extra=(), ylog=None, guide=None,
                  caption=None, mark_multimode=False):
    """One figure, subplot per (family, y_key); traces are (strain_bound,
    <group_extra combo>), x-axis swept, all other coordinates held at
    `fixed` (dict of family-independent fixed values) -- never a merged
    trace across different fixed coordinates (house style). `caption`
    (L18: "annotate every exported figure with a one-line qualification",
    M7: always combined with COMMON_FIGURE_QUALIFICATION) is printed as a
    figure-level footer; `mark_multimode` (L18: "mark multimode / non-
    resetting rows") DROPS vertical_photonic points with single_mode=False
    from the joined strain-bound trace (M7 fix: they used to be plotted
    TWICE, once joined into the relaxed line and again as an overlay
    cross) and instead overlays them ONLY as distinct gray x-markers,
    never joined into any trace. Points with one_pair_valid=False are
    additionally overlaid once per subplot with a distinct open-circle
    marker and a single legend entry (M7: "mark one_pair_valid False
    points with a distinct marker and legend entry")."""
    import matplotlib.pyplot as plt
    contract = {}
    n_rows_plot = len(FAMILIES)
    n_cols = len(y_keys)
    fig, axes = plt.subplots(n_rows_plot, n_cols, figsize=(6.0 * n_cols, 4.2 * n_rows_plot), squeeze=False)
    for fi, family in enumerate(FAMILIES):
        # Non-swept radius/height/x_in default to THIS family's own
        # reference geometry (REF_GEOM) so a single figure call can serve
        # both families even though their reference core radii differ
        # (12.5 nm horizontal vs 80 nm vertical); `fixed` (regime/T_hs/
        # explicit overrides) is applied on top and never fixes x_key
        # itself, so no trace ever joins two different swept-axis values.
        base_filter = dict(REF_GEOM[family])
        base_filter.pop(x_key, None)
        base_filter.update(fixed)
        base_filter.pop(x_key, None)
        base_filter["family"] = family
        pool = _group(core, **base_filter)
        extra_vals = sorted({tuple(r.get(k) for k in group_extra) for r in pool}) if group_extra else [()]
        for ci, y_key in enumerate(y_keys):
            ax = axes[fi][ci]
            opv_false_pts = []
            for sb in STRAIN_BOUNDS:
                style, label = BOUND_STYLE[sb]
                for ev in extra_vals:
                    ef = dict(zip(group_extra, ev))
                    rr = [r for r in pool if r.get("strain_bound") == sb
                          and all(r.get(k) == v for k, v in ef.items())
                          and _finite_num(r.get(x_key)) and _finite_num(r.get(y_key))]
                    if mark_multimode and family == "vertical_photonic":
                        # M7 fix: multimode (headline_eligible False)
                        # points are DROPPED from the joined line trace --
                        # they are drawn ONLY as the distinct cross
                        # overlay below, never joined into the relaxed/
                        # unrelaxed trace beneath it.
                        rr = [r for r in rr if r.get("single_mode") not in (False, "False")]
                    if not rr:
                        continue
                    rr = sorted(rr, key=lambda z: float(z[x_key]))
                    xv = [float(z[x_key]) for z in rr]; yv = [float(z[y_key]) for z in rr]
                    lab = label if not ev else f"{label} {ev}"
                    ax.plot(xv, yv, style, ms=4, lw=1.3, label=lab)
                    contract[f"{family}|{y_key}|{lab}"] = {"row_ids": [z["row_id"] for z in rr], "x": xv, "y": yv}
                    opv_false_pts += [r for r in rr if r.get("one_pair_valid") in (False, "False")]
            if opv_false_pts:
                oxv = [float(z[x_key]) for z in opv_false_pts]; oyv = [float(z[y_key]) for z in opv_false_pts]
                ax.scatter(oxv, oyv, marker="o", facecolors="none", edgecolors="red", s=55, linewidths=1.1,
                           label="one_pair_valid=False", zorder=5)
                contract[f"{family}|{y_key}|one_pair_valid_False"] = {
                    "row_ids": [z["row_id"] for z in opv_false_pts], "x": oxv, "y": oyv}
            if mark_multimode and family == "vertical_photonic":
                mm = [r for r in pool if r.get("single_mode") in (False, "False")
                      and _finite_num(r.get(x_key)) and _finite_num(r.get(y_key))]
                if mm:
                    mm = sorted(mm, key=lambda z: float(z[x_key]))
                    mxv = [float(z[x_key]) for z in mm]; myv = [float(z[y_key]) for z in mm]
                    ax.scatter(mxv, myv, marker="x", c="0.5", s=22, label="multimode (headline_eligible=False)")
                    contract[f"{family}|{y_key}|multimode"] = {"row_ids": [z["row_id"] for z in mm], "x": mxv, "y": myv}
            if guide is not None and y_key == guide[0]:
                ax.axhline(guide[1], color="k", lw=0.9, ls="--", label=guide[2])
            use_ylog = (y_key in ("collected_flux_pulsed_s", "collected_flux_delivered_s")) if ylog is None else ylog
            if use_ylog:
                ax.set_yscale("log")
            ax.set(xlabel=x_key, ylabel=y_key, title=f"{family} {y_key}")
            _leg(ax)
    fig.suptitle(title, fontsize=9)
    footer_text = _figure_footer(caption)
    fig.text(0.5, 0.01, footer_text, ha="center", fontsize=6.5, wrap=True)
    contract["__caption__"] = footer_text
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(out / name, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return contract


def plot_delivered_vs_commanded(out, core):
    """M6 fix: the prior version plotted ONLY the R_s_ohm=1e6 sensitivity
    rows, so the as-built 2.38 GOhm branch (both families' CORE-grid
    default) never appeared and the figure showed delivered==commanded only.
    Both as-built core rows (per family, at their own default R_s_ohm) AND
    the horizontal R_s_ohm=1e6 designed-contact sensitivity rows are plotted
    here, labelled separately."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    contract = {}
    series = []
    for family in FAMILIES:
        pts = [r for r in core if r.get("row_kind") == "core" and r.get("family") == family
               and r.get("regime") == "deterministic_pair"
               and _finite_num(r.get("collected_flux_pulsed_s")) and _finite_num(r.get("collected_flux_delivered_s"))]
        if pts:
            r_s_val = pts[0].get("R_s_ohm")
            # M7 fix: only horizontal_as_built's core default is genuinely
            # "as-built" (the Deshpande device); vertical_photonic's core
            # default is a DESIGNED contact (R_s_ohm=1e6), never fabricated
            # -- the legend must not call it "as-built" too.
            tag = "as-built" if family == "horizontal_as_built" else "designed"
            series.append((f"{family} {tag} core default (R_s_ohm={float(r_s_val):.3g})",
                            pts, "^" if family == "horizontal_as_built" else "v"))
    sens_pts = [r for r in core if r.get("row_kind") == "sensitivity" and r.get("sensitivity_axis") == "R_s_ohm"
                and r.get("family") == "horizontal_as_built" and r.get("regime") == "deterministic_pair"
                and _finite_num(r.get("collected_flux_pulsed_s")) and _finite_num(r.get("collected_flux_delivered_s"))]
    if sens_pts:
        series.append(("horizontal_as_built designed contact (R_s_ohm=1e6, sensitivity)", sens_pts, "s"))
    all_x = []
    for label, pts, marker in series:
        commanded = [float(z["collected_flux_pulsed_s"]) for z in pts]
        delivered = [float(z["collected_flux_delivered_s"]) for z in pts]
        all_x += commanded + delivered
        ax.scatter(commanded, delivered, marker=marker, alpha=0.6, s=18, label=label)
        contract[label] = {"row_ids": [z["row_id"] for z in pts], "x": commanded, "y": delivered}
    if all_x:
        lo, hi = min(v for v in all_x if v > 0), max(all_x)
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="delivered == commanded")
    ax.set(xlabel="collected_flux_pulsed_s (commanded/idealized)", ylabel="collected_flux_delivered_s (RC-limited)",
           xscale="log", yscale="log", title="delivered vs commanded flux (SET regime, RC diagnostic)")
    _leg(ax)
    footer_text = _figure_footer("RC caveat: the as-built 2.38 GOhm horizontal contact cannot deliver a 100 ps "
                                   "step; only the R_s_ohm=1e6 designed contact and the vertical family's own "
                                   "1e6 default approach delivered==commanded.")
    fig.text(0.5, 0.01, footer_text, ha="center", fontsize=6.5, wrap=True)
    contract["__caption__"] = footer_text
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(out / "delivered_vs_commanded_flux.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract


def _max_over_radius(rr, y_key):
    """Group `rr` by core_radius_nm and take the row that MAXIMIZES y_key at
    each radius; returns (xs, ys, ids) with `ids` the row_id that actually
    PRODUCED each y value (L12 fix: the prior version plotted the max but
    recorded the FIRST row's id at that radius, a traceability mismatch)."""
    seen = {}
    for r in rr:
        seen.setdefault(float(r["core_radius_nm"]), []).append(r)
    xs = sorted(seen)
    best_rows = [max(seen[x], key=lambda z: float(z[y_key])) for x in xs]
    ys = [float(z[y_key]) for z in best_rows]
    ids = [z["row_id"] for z in best_rows]
    return xs, ys, ids


def plot_hardware_screens(out, core):
    """M5 fix: the RTI panel's original `rti_level_margin_kT` is identically
    0 on every priced SET row (a flat, non-informative line -- see results.md
    "Gate anti-monotonicity"/E_C/kT wall sections); it is DROPPED here and
    replaced with two VARYING RTI diagnostics, `rti_bypass_fraction` and
    `rti_alignment_error_e_meV`, so the figure actually shows something."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.4))
    contract = {}
    ax = axes[0]
    for family in FAMILIES:
        rr = [r for r in core if r.get("row_kind") == "core" and r.get("family") == family
              and r.get("regime") == "deterministic_pair" and r.get("strain_bound") == "relaxed"
              and _finite_num(r.get("core_radius_nm")) and _finite_num(r.get("set_EC_over_kT"))]
        xs, ys, ids = _max_over_radius(rr, "set_EC_over_kT")
        if xs:
            ax.plot(xs, ys, "o-", label=family)
            contract[f"set_EC_over_kT|{family}"] = {"row_ids": ids, "x": xs, "y": ys}
    ax.axhline(10.0, color="k", ls="--", lw=0.9, label="ec_margin threshold (10 kT)")
    ax.set(xlabel="core_radius_nm", ylabel="set_EC_over_kT", title="Coulomb-blockade screen", yscale="log")
    _leg(ax)
    for ci, (y_key, title) in enumerate((("rti_bypass_fraction", "RT injector bypass fraction"),
                                          ("rti_alignment_error_e_meV", "RT injector electron alignment error"))):
        ax = axes[ci + 1]
        for family in FAMILIES:
            rr = [r for r in core if r.get("row_kind") == "core" and r.get("family") == family
                  and r.get("regime") == "deterministic_pair" and r.get("strain_bound") == "relaxed"
                  and _finite_num(r.get("core_radius_nm")) and _finite_num(r.get(y_key))]
            xs, ys, ids = _max_over_radius(rr, y_key)
            if xs:
                ax.plot(xs, ys, "o-", label=family)
                contract[f"{y_key}|{family}"] = {"row_ids": ids, "x": xs, "y": ys}
        ax.set(xlabel="core_radius_nm", ylabel=y_key, title=title)
        _leg(ax)
    fig.suptitle("Non-gating hardware screens (conditional engineering screens, not measured devices)", fontsize=9)
    footer_text = _figure_footer("Charging wall: deterministic loading at 230-300K fails the Coulomb-blockade "
                                   "screen at every priced core_radius_nm>=10 nm; the RT injector screen is "
                                   "rti_status=unknown_incomplete on every SET row (conditional-engineering "
                                   "screens, neither is a demonstrated hardware result).")
    fig.text(0.5, 0.01, footer_text, ha="center", fontsize=6.5, wrap=True)
    contract["__caption__"] = footer_text
    fig.tight_layout(rect=(0, 0.10, 1, 0.93))
    fig.savefig(out / "hardware_screens.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract


def plot_reversal_map(out, core):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(FAMILIES), figsize=(6.4 * len(FAMILIES), 4.6), squeeze=False)
    contract = {}
    for fi, family in enumerate(FAMILIES):
        ax = axes[0][fi]
        rr = [r for r in core if r.get("row_kind") == "core" and r.get("family") == family
              and r.get("regime") == "deterministic_pair"]
        radii = sorted({float(r["core_radius_nm"]) for r in rr if _finite_num(r.get("core_radius_nm"))})
        temps = sorted({float(r["T_hs"]) for r in rr if _finite_num(r.get("T_hs"))})
        mat = []
        cellmap = {}
        for t in temps:
            row_vals = []
            for r0 in radii:
                cand = [r for r in rr if _close(float(r["core_radius_nm"]), r0) and _close(float(r["T_hs"]), t)]
                revs = [1.0 if r.get("bound_reversal_pair") == "True" else (0.0 if r.get("bound_reversal_pair") == "False" else float("nan")) for r in cand]
                finite_revs = [v for v in revs if math.isfinite(v)]
                row_vals.append(max(finite_revs) if finite_revs else float("nan"))
                cellmap[f"r{r0:g}_T{t:g}"] = [r["row_id"] for r in cand]
            mat.append(row_vals)
        im = ax.imshow(mat, vmin=0, vmax=1, aspect="auto", origin="lower", cmap="coolwarm")
        fig.colorbar(im, ax=ax, label="bound_reversal (1=reversed)")
        ax.set(xticks=range(len(radii)), xticklabels=[f"{v:g}" for v in radii],
               yticks=range(len(temps)), yticklabels=[f"{v:g}" for v in temps],
               xlabel="core_radius_nm", ylabel="T_hs K", title=f"{family} strain-pair reversal map")
        contract[family] = cellmap
    fig.suptitle("Strain-bound reversal map (M11: legend for not_comparable)", fontsize=9)
    footer_text = _figure_footer("White/blank cells = not_comparable (no pair where BOTH strain-bound rows "
                                   "clear the 1000/s optical floor); colored cells are 0=not reversed, "
                                   "1=reversed on mu/collected_flux_pulsed_s/g2_op (never averaged with "
                                   "unreversed cells).")
    fig.text(0.5, 0.01, footer_text, ha="center", fontsize=6.5, wrap=True)
    contract["__caption__"] = footer_text
    fig.tight_layout(rect=(0, 0.10, 1, 0.92))
    fig.savefig(out / "strain_reversal_map.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    return contract


def _safe_plot(name, fail_log, fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 -- never abort the sweep on a figure failure
        fail_log.append({"figure": name, "reason": f"{type(exc).__name__}: {exc}"})
        return None


# ------------------------------------------------------------------ VERDICT
BOUND_ROLE = {"unrelaxed": "conservative_lower", "relaxed": "headline_upper"}


def verdict_lines(core, quick):
    """One VERDICT line per (family, regime, strain_bound, rep_rate_hz)
    combination -- 2*2*2*2=16 lines (the contract's VERDICT template names
    rep_rate_hz as a per-line field, so the rate is a further split beyond
    the "family x regime x strain_bound" grouping the piece-9 spec text
    itself names). `paired_optical_pass` is defined and printed exactly
    once, matching the contract's own definition."""
    lines = ["Definition: paired_optical_pass counts an optical_pass row "
             "whose strain-bound PARTNER row (same family/regime/rep_rate_hz/"
             "core_radius_nm/height_nm/x_in/T_hs, opposite strain_bound) is "
             "present in this run -- NOT an AND requiring both strain "
             "scenarios to individually pass.", ""]
    all_core = [r for r in core if r.get("row_kind") == "core"]
    verdicts = []
    for family in FAMILIES:
        expected_all = build_core(family, quick)
        for regime in REGIMES:
            for strain_bound in STRAIN_BOUNDS:
                for rate in REP_RATES:
                    group = _group(all_core, family=family, regime=regime,
                                    strain_bound=strain_bound, rep_rate_hz=rate)
                    expected = [p for p in expected_all if p["regime"] == regime
                                and p["strain_bound"] == strain_bound and p["rep_rate_hz"] == rate]
                    n, total = len(group), len(expected)
                    invalid = sum(1 for r in group if not r.get("valid"))
                    eligible = sum(1 for r in group if r.get("eligible"))
                    paired = 0
                    for r in group:
                        if not r.get("optical_pass"):
                            continue
                        if _bound_partner(r, all_core) is not None:
                            paired += 1
                    quality = sum(1 for r in group if r.get("quality_pass"))
                    hw_q = sum(1 for r in group if r.get("hardware_qualified"))
                    rti_q = sum(1 for r in group if r.get("rti_qualified"))
                    any_optical = any(r.get("optical_pass") for r in group)
                    if regime == "deterministic_pair":
                        if hw_q > 0:
                            ideal = "pass_hardware_qualified"
                        elif rti_q > 0:
                            ideal = "pass_rti_qualified"
                        elif any_optical:
                            ideal = "pass_hardware_infeasible"
                        else:
                            ideal = "no_idealized_pass"
                    else:
                        ideal = "pass" if any_optical else "no_idealized_pass"
                    complete = bool((not quick) and n == total)
                    screening = group[0].get("screening_fraction", 0.0) if group else 0.0
                    access = group[0].get("occupied_dot_access", 0.05) if group else 0.05
                    v = dict(idealized_status=ideal, family=family, regime=regime, strain_bound=strain_bound,
                             bound_role=BOUND_ROLE[strain_bound], rep_rate_hz=rate, complete=complete,
                             eligible=eligible, paired_optical_pass=paired, quality_pass=quality,
                             hardware_qualified=hw_q, rti_qualified=rti_q, coverage=f"{n}/{total}", invalid=invalid,
                             screening=screening, access=access)
                    verdicts.append(v)
                    # quality_pass (H3 fix) is printed BESIDE paired_optical_pass:
                    # optical_pass/paired_optical_pass alone can be satisfied by
                    # a row that empties its dot faster only because it LOSES
                    # more photons (surface loss, tight access), so quality_pass
                    # =optical_pass AND photons_per_cycle>=0.01 is a second,
                    # independent throughput screen -- see the "Gate
                    # anti-monotonicity" results.md section. This field is a
                    # script/orchestrator addition; it does NOT change the
                    # contract's own frozen optical_pass definition or its
                    # VERDICT template (docs/nitride_nanowire_contract.md is
                    # out of scope for this piece).
                    lines.append(
                        "VERDICT: idealized_status=%s family=%s regime=%s strain_bound=%s bound_role=%s "
                        "rep_rate_hz=%g complete=%s eligible=%d paired_optical_pass=%d quality_pass=%d "
                        "hardware_qualified=%d rti_qualified=%d coverage=%s invalid=%d flux_floor=1000/s "
                        "screening=%g access=%g" % (
                            ideal, family, regime, strain_bound, BOUND_ROLE[strain_bound], rate, complete,
                            eligible, paired, quality, hw_q, rti_q, f"{n}/{total}", invalid,
                            float(screening), float(access)))
    return lines, verdicts


def _above_floor(r):
    return r.get("valid") is True and _finite_num(r.get("collected_flux_pulsed_s")) and float(r["collected_flux_pulsed_s"]) >= 1000.0


def attach_bound_reversal(core_rows):
    """Declared transform (contract: "aggregate/interpolate only with
    declared transforms"): reproduces nitride_nanowire_device._bound_
    reversal's own ordering rule (mu / collected_flux_pulsed_s / g2_op,
    both partners above the 1000/s optical floor) directly on this run's
    already-evaluated row PAIRS, since this script calls evaluate_nanowire
    once per (family,regime,rate,radius,height,x_in,T_hs,strain_bound)
    combination via fsim_core.device.evaluate() rather than pairing bounds
    through evaluate_strain_pair -- no second evaluate() call and no
    re-derivation of the physics itself, only the same declared comparison
    the device module documents. Written as `bound_reversal_pair` (a new
    column) so it is never confused with the device's own per-row
    `bound_reversal` field, which stays the "not_computed" string sentinel
    on every row here (device.py: that field is only set by
    evaluate_strain_pair, which this script does not call)."""
    only_core = [r for r in core_rows if r.get("row_kind") == "core"]
    for r in only_core:
        r["bound_reversal_pair"] = "not_computed"
    for r in only_core:
        if r.get("strain_bound") != "relaxed":
            continue
        partner = _bound_partner(r, only_core)
        if partner is None:
            r["bound_reversal_pair"] = "not_comparable"
            continue
        if not (_above_floor(r) and _above_floor(partner)):
            reversal = "not_comparable"
        else:
            checks = []
            for key, higher_is_headline in (("mu", True), ("collected_flux_pulsed_s", True), ("g2_op", False)):
                rv, uv = r.get(key), partner.get(key)
                if not (_finite_num(rv) and _finite_num(uv)):
                    continue
                checks.append(float(rv) < float(uv) if higher_is_headline else float(rv) > float(uv))
            reversal = "not_comparable" if not checks else str(any(checks))
        r["bound_reversal_pair"] = reversal
        partner["bound_reversal_pair"] = reversal


def _headline_eligible_of(r):
    return r.get("headline_eligible") in (True, "True")


def best_passing_flux_lines(core):
    """H2 fix: headline selection requires optical_pass AND headline_eligible
    (a vertical_photonic row above the LP11 single-mode cutoff is never
    nominated, even if it is the brightest optical_pass row); optical-only
    counts are kept and printed SEPARATELY so the eligibility filter's
    effect is visible, not hidden. Every BEST line prints the row's full
    geometry/rate/bound, BOTH commanded and RC-delivered flux,
    headline_eligible, hardware_qualified/rti_qualified (False on every
    idealized-loading row), blocked_load_probability (H3: "report
    blocked_load_probability on every SET nomination"), and the row's
    unrelaxed strain-bound partner's row_id/flux/bound_reversal_pair value
    (never just the commanded-flux number in isolation)."""
    lines = []
    all_core = [r for r in core if r.get("row_kind") == "core"]
    for family in FAMILIES:
        for label, pred in (("BEST_PASSING_FLUX", lambda r: True),
                             ("BEST_PASSING_FLUX_300K", lambda r: _close(float(r.get("T_hs", -1)), 300.0))):
            optical_cand = [r for r in all_core if r.get("family") == family and r.get("optical_pass")
                             and pred(r) and _finite_num(r.get("collected_flux_pulsed_s"))]
            elig_cand = [r for r in optical_cand if _headline_eligible_of(r)]
            n_optical, n_elig = len(optical_cand), len(elig_cand)
            if not elig_cand:
                lines.append(f"{label} family={family} value=none row_id=none "
                              f"optical_pass_candidates={n_optical} headline_eligible_candidates={n_elig}")
                continue
            best = max(elig_cand, key=lambda r: float(r["collected_flux_pulsed_s"]))
            partner = _bound_partner(best, all_core)
            delivered = best.get("collected_flux_delivered_s")
            delivered_txt = f"{float(delivered):.6g}" if _finite_num(delivered) else "n/a"
            blocked_txt = (f"{float(best.get('blocked_load_probability')):.6g}"
                            if _finite_num(best.get("blocked_load_probability")) else "n/a")
            # M4 fix: photons_per_cycle/quality_pass use COMMANDED flux
            # (collected_flux_pulsed_s) -- printed here explicitly, beside
            # a REPORT-derived photons_per_cycle_delivered (delivered
            # flux/rep_rate_hz, never written back into sweep.csv) so the
            # BEST line never implies the commanded per-cycle number is
            # what an RC-limited detector would actually see.
            ppc_commanded = best.get("photons_per_cycle")
            ppc_commanded_txt = f"{float(ppc_commanded):.6g}" if _finite_num(ppc_commanded) else "n/a"
            rep_hz_best = best.get("rep_rate_hz")
            if _finite_num(delivered) and _finite_num(rep_hz_best) and float(rep_hz_best) > 0:
                ppc_delivered = float(delivered) / float(rep_hz_best)
                ppc_delivered_txt = f"{ppc_delivered:.6g}"
                cycles_per_photon_txt = f"{1.0 / ppc_delivered:.4g}" if ppc_delivered > 0 else "inf"
            else:
                ppc_delivered_txt = "n/a"; cycles_per_photon_txt = "n/a"
            if partner is not None:
                p_flux = partner.get("collected_flux_pulsed_s")
                p_flux_txt = f"{float(p_flux):.6g}" if _finite_num(p_flux) else "n/a"
                partner_txt = (f"unrelaxed_partner_row_id={partner['row_id']} "
                                f"unrelaxed_partner_flux={p_flux_txt} "
                                f"bound_reversal={best.get('bound_reversal_pair', 'not_computed')}")
            else:
                partner_txt = "unrelaxed_partner_row_id=none unrelaxed_partner_flux=n/a bound_reversal=not_comparable"
            lines.append(
                f"{label} family={family} value={float(best['collected_flux_pulsed_s']):.6g} "
                f"row_id={best['row_id']} core_radius_nm={float(best.get('core_radius_nm', float('nan'))):g} "
                f"height_nm={float(best.get('height_nm', float('nan'))):g} "
                f"x_in={float(best.get('x_in', float('nan'))):g} T_hs={float(best['T_hs']):g} "
                f"rep_rate_hz={float(best['rep_rate_hz']):g} strain_bound={best['strain_bound']} "
                f"screening={float(best.get('screening_fraction', 0.0)):g} regime={best['regime']} "
                f"commanded_flux={float(best['collected_flux_pulsed_s']):.6g} delivered_flux={delivered_txt} "
                f"photons_per_cycle_commanded={ppc_commanded_txt} photons_per_cycle_delivered={ppc_delivered_txt} "
                f"(one delivered photon per {cycles_per_photon_txt} cycles) "
                f"headline_eligible={_headline_eligible_of(best)} "
                f"hardware_qualified={bool(best.get('hardware_qualified'))} "
                f"rti_qualified={bool(best.get('rti_qualified'))} blocked_load_probability={blocked_txt} "
                f"quality_pass={bool(best.get('quality_pass'))} {partner_txt} "
                f"optical_pass_candidates={n_optical} headline_eligible_candidates={n_elig}")
    return lines


def bounds_table(core):
    """Every headline pairing at the reference geometry, side by side --
    conservative unrelaxed vs headline relaxed, both regimes/rates, T in
    {230,300}, per family (contract: "Every headline table/figure pairs
    conservative unrelaxed and headline relaxed strain at identical other
    inputs")."""
    all_core = [r for r in core if r.get("row_kind") == "core"]
    lines = ["| family | regime | rep_rate_hz | T_hs K | unrelaxed g2 | unrelaxed flux/s | relaxed g2 | relaxed flux/s | bound_reversal | row_ids |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for family in FAMILIES:
        ref = REF_GEOM[family]
        for regime in REGIMES:
            for rate in REP_RATES:
                for t in (230.0, 300.0):
                    u = _group(all_core, family=family, regime=regime, rep_rate_hz=rate, T_hs=t,
                                strain_bound="unrelaxed", core_radius_nm=ref["core_radius_nm"],
                                height_nm=ref["height_nm"], x_in=ref["x_in"])
                    rl = _group(all_core, family=family, regime=regime, rep_rate_hz=rate, T_hs=t,
                                 strain_bound="relaxed", core_radius_nm=ref["core_radius_nm"],
                                 height_nm=ref["height_nm"], x_in=ref["x_in"])
                    if not u or not rl:
                        continue
                    u, rl = u[0], rl[0]
                    lines.append(f"| {family} | {regime} | {rate:.3g} | {t:g} | {u.get('g2_op')} | "
                                  f"{u.get('collected_flux_pulsed_s')} | {rl.get('g2_op')} | "
                                  f"{rl.get('collected_flux_pulsed_s')} | {rl.get('bound_reversal_pair')} | "
                                  f"{u.get('row_id')},{rl.get('row_id')} |")
    return lines


_SI_RANGE_LAMBDA_RE = re.compile(r"si_complex_index: lambda_nm=([0-9.eE+-]+) outside")


def _qcse_levels_replay(row):
    """H2 fix: `field_kVcm` stored on an INVALID row is the transport
    DEPLETION field alone -- fsim_core/nitride_nanowire_device.py's
    exception-path fallback sets `_field_kVcm = inj["depletion_field_kVcm"]`
    before the row goes invalid (see its own comment there) -- NEVER the
    built-in piezoelectric/polarization field, so quoting it as "the
    built-in piezoelectric field" is wrong by orders of magnitude. This
    replays ONLY fsim_core.nitride_nanowire_levels (never scripts/
    run_nitride_nanowire.py's own evaluate() path, never a second
    evaluate_nanowire() call -- a read-only replay of an already-committed
    dependency module, not this piece's own production routine) at the
    row's own x_in/height_nm/core_radius_nm/outer_radius_nm/
    screening_fraction, the UNRELAXED strain bound, with
    external_field_kVcm set to the row's OWN recorded field_kVcm --
    reproducing device.py's own post-feedback sysB convention
    (external_field_kVcm = ext_field_kVcm(0.0 [A] card default) +
    depletion_field_kVcm, and depletion_field_kVcm IS this row's own
    field_kVcm) -- so this recovers the row's OWN true total field, not a
    new physics result. Returns the NanowireLevels replay object, or None
    if the row is not horizontal_as_built (the only family with an
    affected invalid row this piece has ever seen; disc_radius_nm for
    vertical_photonic is not carried as its own sweep.csv column) or the
    replay itself fails to converge."""
    if row.get("family") != "horizontal_as_built":
        return None
    depletion = row.get("field_kVcm")
    if not _finite_num(depletion):
        return None
    try:
        from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels
        sys_ = NitrideNanowireSystem(
            height_nm=float(row.get("height_nm")), core_radius_nm=float(row.get("core_radius_nm")),
            outer_radius_nm=float(row.get("outer_radius_nm") or row.get("core_radius_nm")),
            disc_radius_nm=None, x_in=float(row.get("x_in")), strain_bound="unrelaxed",
            screening_fraction=float(row.get("screening_fraction") or 0.0),
            external_field_kVcm=float(depletion))
        t_k = float(row.get("T_j")) if _finite_num(row.get("T_j")) else float(row.get("T_hs", 300.0))
        lv = levels(sys_, T_K=t_k)
    except Exception:
        return None
    return lv if lv.valid else None


def _diagnostic_lambda_from_invalid_reasons(row):
    """When a row is invalid ONLY because the DOWNSTREAM photonics
    Si-substrate complex-index table (fsim_core/nitride_nanowire_photonics.py,
    piece 3) covers only its own tabulated range (read at call sites from
    the module's own `_SI_INDEX_ANCHORS_NM`, 380-750 nm as committed --
    never hardcoded here), the bare-dot emission wavelength that
    fsim_core/nitride_nanowire_levels.py (piece 2) already computed upstream
    is still recoverable from the row's own invalid_reasons diagnostic
    string -- a read of an already-computed, already-reported number, never
    a re-derivation or an assumed n_substrate override. Returns
    (value_or_None, recovered_bool); recovered_bool distinguishes this
    diagnostic-only recovery from a genuinely valid row's own lambda_nm."""
    try:
        reasons = json.loads(row.get("invalid_reasons", "[]"))
    except (TypeError, json.JSONDecodeError):
        return None, False
    for reason in reasons:
        mo = _SI_RANGE_LAMBDA_RE.search(reason)
        if mo:
            return float(mo.group(1)), True
    return None, False


def _rate_pair_rows(all_core, family, T_hs=300.0, strain_bound="relaxed", x_in=0.40):
    ref = REF_GEOM[family]
    out = {}
    for rate in REP_RATES:
        m = _group(all_core, family=family, regime="deterministic_pair", rep_rate_hz=rate, T_hs=T_hs,
                    strain_bound=strain_bound, core_radius_nm=ref["core_radius_nm"],
                    height_nm=ref["height_nm"], x_in=x_in)
        if m:
            out[rate] = m[0]
    return out


def _optical_pass_of(r):
    return bool(r.get("optical_pass") in (True, "True"))


def _pick(rows, **filters):
    m = _group(rows, **filters)
    return m[0] if m else None


# ------------------------------------------- H3/M9-M10 new results.md sections
def _find_axis_flip_pair(sens, family, axis, lo_val, hi_val, strain_bound="relaxed"):
    """M3 fix: search PROGRAMMATICALLY across (T_hs, rep_rate_hz) at the
    family's reference geometry (never hardcode the rate/T that happens to
    demonstrate it) for the `axis` cut pair (lo_val, hi_val) where
    one_pair_valid flips False->True while collected_flux_pulsed_s FALLS
    -- the actual gate anti-monotonicity evidence. A pair where
    one_pair_valid is False (or True) at BOTH endpoints demonstrates
    nothing and is skipped. Returns (T_hs, rep_rate_hz, lo_row, hi_row) or
    (None, None, None, None) if no combination in this run's coverage
    demonstrates the flip."""
    ref = REF_GEOM[family]
    for t in CUT_T_HS:
        for rate in REP_RATES:
            filt = dict(family=family, regime="deterministic_pair", strain_bound=strain_bound,
                        T_hs=t, rep_rate_hz=rate, core_radius_nm=ref["core_radius_nm"],
                        height_nm=ref["height_nm"], x_in=ref["x_in"])
            lo = _pick(sens, sensitivity_axis=axis, sensitivity_value=lo_val, **filt)
            hi = _pick(sens, sensitivity_axis=axis, sensitivity_value=hi_val, **filt)
            if lo is None or hi is None:
                continue
            lo_opv = lo.get("one_pair_valid") in (True, "True")
            hi_opv = hi.get("one_pair_valid") in (True, "True")
            if (not lo_opv) and hi_opv and _finite_num(lo.get("collected_flux_pulsed_s")) \
               and _finite_num(hi.get("collected_flux_pulsed_s")) \
               and float(hi["collected_flux_pulsed_s"]) < float(lo["collected_flux_pulsed_s"]):
                return t, rate, lo, hi
    return None, None, None, None


def _gate_anti_monotonicity_section(core, all_core, quick):
    """H3 fix: one_pair_valid (blocked_load_probability<=1e-9) can be
    cleared by ADDING sidewall loss (faster occupied-dot emptying) as
    readily as by improving device quality -- demonstrated with the S_cm_s
    and occupied_dot_access reduced-cut rows at the horizontal reference
    geometry, plus the temperature/composition pattern of which core rows
    actually pass, all read from the rows themselves, never hardcoded."""
    lines = ["## Gate anti-monotonicity (H3 obligation)", "",
             "one_pair_valid can be satisfied by ADDING sidewall loss (faster occupied-dot emptying) as "
             "readily as by improving device quality -- optical_pass/paired_optical_pass alone therefore do "
             "NOT certify throughput. quality_pass=optical_pass AND photons_per_cycle>=0.01 [A, orchestrator "
             "threshold: one collected photon per hundred cycles] is reported beside optical_pass/"
             "paired_optical_pass in every VERDICT line and the per-temperature tables below. M4 fix: "
             "photons_per_cycle (sweep.csv column) and quality_pass BOTH use COMMANDED flux "
             "(collected_flux_pulsed_s/rep_rate_hz), never the RC-limited delivered flux -- the BEST lines "
             "below additionally print a report-derived photons_per_cycle_delivered "
             "(collected_flux_delivered_s/rep_rate_hz) beside it so the commanded per-cycle number is never "
             "mistaken for what an RC-limited detector would actually see.", ""]
    sens = [r for r in core if r.get("row_kind") == "sensitivity"]
    fam = "horizontal_as_built"
    ref = REF_GEOM[fam]
    base_filter = dict(family=fam, regime="deterministic_pair", strain_bound="relaxed",
                        T_hs=300.0, rep_rate_hz=200.0e6, core_radius_nm=ref["core_radius_nm"],
                        height_nm=ref["height_nm"], x_in=ref["x_in"])
    ref_row = _pick(all_core, **base_filter)
    s_lo = _pick(sens, sensitivity_axis="S_cm_s", sensitivity_value=100.0, **base_filter)
    s_hi = _pick(sens, sensitivity_axis="S_cm_s", sensitivity_value=10000.0, **base_filter)
    acc_hi = _pick(sens, sensitivity_axis="occupied_dot_access", sensitivity_value=1.0, **base_filter)

    def _fmt(r, label):
        if r is None:
            return f"{label}: not in this run's coverage"
        return (f"{label} (row {r.get('row_id')}): one_pair_valid={r.get('one_pair_valid')} "
                f"g2_op={r.get('g2_op')} flux={r.get('collected_flux_pulsed_s')} "
                f"quality_pass={r.get('quality_pass')}")

    # M3 fix: the 300K/200MHz reference-condition S_cm_s triple below does
    # NOT demonstrate the flip (one_pair_valid is False at all three
    # S_cm_s values there) -- the pair that DOES is found programmatically
    # across (T_hs, rep_rate_hz) and quoted FIRST, never the non-
    # demonstrating triple alone.
    t_flip, rate_flip, s_lo_flip, s_hi_flip = _find_axis_flip_pair(sens, fam, "S_cm_s", 100.0, 10000.0)
    lines.append("S_cm_s pair that actually DEMONSTRATES the anti-monotonicity (selected programmatically: "
                 "the (T_hs, rep_rate_hz) combination at the horizontal reference geometry where "
                 "one_pair_valid flips False->True while collected_flux_pulsed_s FALLS; the 300K/200MHz "
                 "reference-condition triple below does NOT demonstrate it -- one_pair_valid is False at "
                 "all three S_cm_s values there):")
    if s_lo_flip is not None:
        lines.append("- " + _fmt(s_lo_flip, f"S_cm_s={s_lo_flip.get('sensitivity_value')} (lower surface loss, "
                                             f"T_hs={t_flip:g}K, {rate_flip/1.0e6:g}MHz)"))
        lines.append("- " + _fmt(s_hi_flip, f"S_cm_s={s_hi_flip.get('sensitivity_value')} (higher surface loss, "
                                             f"T_hs={t_flip:g}K, {rate_flip/1.0e6:g}MHz)"))
        flux_lo = float(s_lo_flip["collected_flux_pulsed_s"]); flux_hi = float(s_hi_flip["collected_flux_pulsed_s"])
        s_lo_val = float(s_lo_flip.get("sensitivity_value")); s_hi_val = float(s_hi_flip.get("sensitivity_value"))
        ratio = flux_lo / flux_hi if flux_hi > 0 else float("nan")
        # M5 fix: quality_pass removes only ABSOLUTELY DIM optical passes;
        # it does NOT repair this loss-induced ordering (the row that
        # newly optical_passes here is the one that lost more flux).
        lines.append(f"quality_pass (photons_per_cycle>=0.01) removes only ABSOLUTELY DIM optical passes; it "
                     f"does NOT repair this loss-induced ordering: row {s_hi_flip.get('row_id')} "
                     f"(S_cm_s={s_hi_flip.get('sensitivity_value')}) keeps quality_pass="
                     f"{s_hi_flip.get('quality_pass')} even though the "
                     f"{(s_hi_val / s_lo_val if s_lo_val else float('nan')):.0f}x surface-recombination "
                     f"increase (S_cm_s {s_lo_flip.get('sensitivity_value')}->{s_hi_flip.get('sensitivity_value')}) "
                     f"cut its commanded flux {ratio:.2g}x ({flux_lo:.4g}->{flux_hi:.4g} /s).")
    else:
        lines.append("- no (T_hs, rep_rate_hz) combination in this run's coverage demonstrates the S_cm_s flip.")
    lines.append("")

    lines.append("S_cm_s and occupied_dot_access triple/pair at the 300K/200MHz reference condition, for "
                 "context only (S_cm_s here does NOT flip; occupied_dot_access DOES):")
    lines.append("- " + _fmt(s_lo, "S_cm_s=100 (lower surface loss)"))
    lines.append("- " + _fmt(ref_row, "S_cm_s=1000 (main-grid default)"))
    lines.append("- " + _fmt(s_hi, "S_cm_s=10000 (higher surface loss)"))
    lines.append("- " + _fmt(ref_row, "occupied_dot_access=0.05 (main-grid default, lower loss)"))
    lines.append("- " + _fmt(acc_hi, "occupied_dot_access=1.0 (higher loss)"))
    lines.append("")

    # L fix: split by strain_bound too (the prior table mixed both bounds
    # into one x_in-only breakdown per T_hs).
    lines.append("Per-T_hs optical_pass counts, horizontal_as_built (both regimes, split by strain_bound), "
                 "with x_in composition (read from the rows, never hardcoded):")
    for regime in REGIMES:
        for sb in STRAIN_BOUNDS:
            parts = []
            for t in T_HS:
                grp = [r for r in all_core if r.get("family") == fam and r.get("regime") == regime
                       and r.get("strain_bound") == sb and _close(float(r.get("T_hs", -1)), t)]
                passing = [r for r in grp if r.get("optical_pass")]
                n25 = sum(1 for r in passing if _close(float(r.get("x_in", -1)), 0.25))
                n40 = sum(1 for r in passing if _close(float(r.get("x_in", -1)), 0.40))
                parts.append(f"T={t:g}K:{len(passing)}(x_in0.25={n25},x_in0.40={n40})")
            lines.append(f"- {fam}/{regime}/{sb}: " + " ".join(parts))
    lines.append("")

    # H1 fix: composition statements are generated PER (family, strain_
    # bound, x_in) COUNT FROM THE ROWS (never hardcoded numbers). The
    # prior text's "these are the ONLY x_in=0.40 optical passes in either
    # family this run" claim was FALSE (the vertical_photonic family
    # passes at x_in=0.40 under the RELAXED bound too) -- this is
    # rewritten per family, never claiming exclusivity across families.
    def _xin_tag(v):
        if not _finite_num(v):
            return None
        fv = float(v)
        if _close(fv, 0.25):
            return 0.25
        if _close(fv, 0.40):
            return 0.40
        return fv

    pass_counts = {}
    for r in all_core:
        if not r.get("optical_pass"):
            continue
        key = (r.get("family"), r.get("strain_bound"), _xin_tag(r.get("x_in")))
        pass_counts.setdefault(key, []).append(r)

    def _n(family, sb, xin):
        return len(pass_counts.get((family, sb, xin), []))

    lines.append("Composition of optical passes, per family (counted from the rows; see the definition of "
                 "'entirely'/'both compositions' inline):")
    for family in FAMILIES:
        n_rel25, n_rel40 = _n(family, "relaxed", 0.25), _n(family, "relaxed", 0.40)
        n_unrel40 = _n(family, "unrelaxed", 0.40)
        n_rel = n_rel25 + n_rel40
        tag = "HORIZONTAL" if family == "horizontal_as_built" else "VERTICAL"
        if n_rel == 0:
            lines.append(f"{tag} ({family}) RELAXED optical passes: 0 total (x_in=0.25: 0, x_in=0.40: 0).")
        elif n_rel40 == 0:
            lines.append(f"{tag} ({family}) RELAXED optical passes: {n_rel} total (x_in=0.25: {n_rel25}, "
                         f"x_in=0.40: {n_rel40}) -- ENTIRELY x_in=0.25, no relaxed x_in=0.40 optical pass "
                         "this run.")
        elif n_rel25 == 0:
            lines.append(f"{tag} ({family}) RELAXED optical passes: {n_rel} total (x_in=0.25: {n_rel25}, "
                         f"x_in=0.40: {n_rel40}) -- ENTIRELY x_in=0.40.")
        else:
            lines.append(f"{tag} ({family}) RELAXED optical passes: {n_rel} total (x_in=0.25: {n_rel25}, "
                         f"x_in=0.40: {n_rel40}) -- this family passes at BOTH compositions this run.")
        rows40 = pass_counts.get((family, "unrelaxed", 0.40), [])
        if n_unrel40:
            fluxes40 = sorted(float(r["collected_flux_pulsed_s"]) for r in rows40
                               if _finite_num(r.get("collected_flux_pulsed_s")))
            ids40 = ",".join(r["row_id"] for r in rows40[:8]) + ("..." if n_unrel40 > 8 else "")
            excl = " (the only x_in=0.40 passes in THIS family)" if n_rel40 == 0 else ""
            lines.append(f"{tag} ({family}) x_in=0.40 UNRELAXED optical passes: {n_unrel40} total{excl}, "
                         f"collected_flux_pulsed_s spans {fluxes40[0]:.4g}-{fluxes40[-1]:.4g} /s, "
                         f"row_ids={ids40}.")
        else:
            lines.append(f"{tag} ({family}) x_in=0.40 UNRELAXED optical passes: 0 total.")
    total_40 = sum(_n(f_, sb, 0.40) for f_ in FAMILIES for sb in STRAIN_BOUNDS)
    total_all = sum(_n(f_, sb, x) for f_ in FAMILIES for sb in STRAIN_BOUNDS for x in (0.25, 0.40))
    lines.append(f"Across BOTH families this run: {total_40}/{total_all} optical passes are x_in=0.40 -- "
                 "the horizontal family's x_in=0.40 passes are UNRELAXED-only, while the vertical family also "
                 "passes at x_in=0.40 under the RELAXED bound (see the per-family lines above); the earlier "
                 "'these are the only x_in=0.40 passes in either family' claim was HORIZONTAL-only and is not "
                 "repeated here.")
    if quick:
        lines.append("Quick-mode coverage note: this run used --quick (a declared deterministic subset); the "
                     "counts above are NOT the full main grid's coverage.")
    lines.append("")
    return lines


def _per_temperature_count_tables(all_core):
    """M9-M10 obligation: per-temperature count tables."""
    lines = ["## Per-temperature counts (M9-M10 obligation)", "",
             "| family | regime | T_hs K | n_core | n_eligible | n_optical_pass | n_quality_pass |",
             "|---|---|---|---|---|---|---|"]
    for family in FAMILIES:
        for regime in REGIMES:
            for t in T_HS:
                grp = [r for r in all_core if r.get("family") == family and r.get("regime") == regime
                       and _close(float(r.get("T_hs", -1)), t)]
                n = len(grp)
                n_elig = sum(1 for r in grp if r.get("eligible"))
                n_opt = sum(1 for r in grp if r.get("optical_pass"))
                n_qual = sum(1 for r in grp if r.get("quality_pass"))
                lines.append(f"| {family} | {regime} | {t:g} | {n} | {n_elig} | {n_opt} | {n_qual} |")
    lines.append("")
    return lines


def _headline_nominations_16(core, all_core):
    """M9-M10 obligation: 'the 16 per-family/regime/bound/rate nominations'
    -- one nominated reference-headline row per (family,regime,strain_bound,
    rep_rate_hz) combination, distinct from the aggregated 4-line
    BEST_PASSING_FLUX(_300K) summary above."""
    lines = ["## 16 per-family/regime/bound/rate reference-headline nominations (M9-M10 obligation)", "",
             "One nominated row per (family,regime,strain_bound,rep_rate_hz) group -- the brightest "
             "optical_pass row in that group (headline_eligible additionally required for vertical_photonic, "
             "H2/bullet 10); 'none' if the group has no such row.", "",
             "| family | regime | strain_bound | rep_rate_hz | row_id | commanded_flux/s | "
             "photons_per_cycle_delivered | blocked_load_probability | g2_op | headline_eligible | "
             "quality_pass |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for family in FAMILIES:
        for regime in REGIMES:
            for sb in STRAIN_BOUNDS:
                for rate in REP_RATES:
                    grp = _group(all_core, family=family, regime=regime, strain_bound=sb, rep_rate_hz=rate)
                    cand = [r for r in grp if r.get("optical_pass")
                            and (family != "vertical_photonic" or _headline_eligible_of(r))
                            and _finite_num(r.get("collected_flux_pulsed_s"))]
                    if not cand:
                        lines.append(f"| {family} | {regime} | {sb} | {rate:.3g} | none | | | | | False | |")
                        continue
                    best = max(cand, key=lambda r: float(r["collected_flux_pulsed_s"]))
                    delivered = best.get("collected_flux_delivered_s")
                    ppc_d = (float(delivered) / rate) if _finite_num(delivered) and rate > 0 else float("nan")
                    ppc_d_txt = f"{ppc_d:.4g}" if math.isfinite(ppc_d) else "n/a"
                    blocked_txt = (f"{float(best.get('blocked_load_probability')):.4g}"
                                    if _finite_num(best.get("blocked_load_probability")) else "n/a")
                    lines.append(f"| {family} | {regime} | {sb} | {rate:.3g} | {best['row_id']} | "
                                 f"{float(best['collected_flux_pulsed_s']):.6g} | {ppc_d_txt} | {blocked_txt} | "
                                 f"{best.get('g2_op')} | {_headline_eligible_of(best)} | {best.get('quality_pass')} |")
    lines.append("")
    return lines


def _ensemble_yield_proxy_section(core):
    """M9-M10 obligation: independent ensemble-yield proxy at 10K/300K via
    the surface module's own yield_ratio helper, non-gating, printed next to
    the measured 0.52 ensemble ratio -- never fit to it."""
    lines = ["## Ensemble-yield proxy at 10K/300K vs measured 0.52 (M9-M10 obligation, non-gating)", ""]
    import fsim_core.nitride_nanowire_surface as _nw_surface
    d13 = [r for r in core if r.get("row_kind") == "sensitivity"
           and r.get("sensitivity_axis") == "deshpande2013_comparison"]
    row10 = _pick(d13, family="horizontal_as_built", core_radius_nm=12.5, I_uA=0.002)
    all_core = [r for r in core if r.get("row_kind") == "core"]
    row300 = _pick(all_core, family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                    T_hs=300.0, rep_rate_hz=80.0e6, core_radius_nm=12.5, height_nm=2.0, x_in=0.25)
    if row10 is None or row300 is None:
        lines.append("Matching 10K/300K row pair (horizontal_as_built, R=12.5nm, h=2nm, x_in=0.25, relaxed, "
                     "2nA pulses/80MHz) not both present in this run's coverage -- proxy not computed.")
        lines.append("")
        return lines

    def _gamma_loss(r):
        # gamma_X_ns (radiative, antenna-corrected) and k_X_ns (nonradiative/
        # surface/escape competition) are ADDITIVE, separate rate channels
        # (levels.rates()'s own convention: total decay = gamma_X_ns +
        # k_X_ns, never k_X_ns already containing gamma_X_ns) -- so
        # loss_ns is k_X_ns directly, never (k_X_ns - gamma_X_ns), which
        # can go negative whenever the antenna-boosted radiative rate
        # exceeds the (small, low-T) nonradiative rate.
        g, k = r.get("gamma_X_ns"), r.get("k_X_ns")
        if not (_finite_num(g) and _finite_num(k)):
            return None, None
        return float(g), float(k)

    g300, l300 = _gamma_loss(row300)
    g10, l10 = _gamma_loss(row10)
    if None in (g300, l300, g10, l10):
        lines.append(f"gamma_X_ns/k_X_ns not finite on row {row300.get('row_id')} (300K) or "
                     f"{row10.get('row_id')} (10K) -- proxy not computed.")
        lines.append("")
        return lines
    lines.append(f"Inputs (this run's own gamma_X_ns radiative / k_X_ns nonradiative additive rate pair, "
                 f"rows {row300.get('row_id')} 300K / {row10.get('row_id')} 10K): gamma_300_ns={g300:.6g}, "
                 f"loss_300_ns={l300:.6g}, gamma_10_ns={g10:.6g}, loss_10_ns={l10:.6g}.")
    try:
        result = _nw_surface.yield_ratio(gamma_300_ns=g300, loss_300_ns=l300, gamma_10_ns=g10, loss_10_ns=l10)
    except _nw_surface.NitrideNanowireSurfaceError as exc:
        lines.append(f"yield_ratio rejected these inputs this run: {exc} -- proxy not computed (never "
                     "silently repaired/clamped).")
        lines.append("")
        return lines
    if result.get("defined"):
        lines.append(f"yield_300K={result['yield_300K']:.6g}, yield_10K={result['yield_10K']:.6g}, MODEL "
                     f"PROXY ratio={result['ratio']:.6g} -- vs measured ensemble PL ratio 0.52 (Deshpande "
                     "2013, ledger deshpande2013_thermal_and_pl). This proxy assumes equal absorption/"
                     "capture/collection between 10K and 300K [A] and is a SINGLE-DOT proxy, never "
                     "verification of the ensemble (many-wire) PL measurement; no channel in this model was "
                     "fit to reproduce 0.52, and none does.")
    else:
        lines.append(f"yield_ratio undefined this run: {result.get('reason')}.")
    lines.append("")
    return lines


def _planar_reference_section(core):
    """M9-M10 obligation: print the 16 planar_reference rows plus a
    reconciliation paragraph against the nanowire tier's own flux range."""
    lines = ["## Planar-reference rows and cross-round reconciliation (M9-M10 obligation)", "",
             "| family_tag | regime | T_hs K | screening | g2_op | collected_flux_pulsed_s | valid | row_id |",
             "|---|---|---|---|---|---|---|---|"]
    planar_rows = [r for r in core if r.get("row_kind") == "planar_reference"]
    for r in sorted(planar_rows, key=lambda z: (str(z.get("family", "")), str(z.get("regime", "")),
                                                   float(z.get("T_hs", 0) or 0), float(z.get("screening_fraction", 0) or 0))):
        lines.append(f"| {r.get('family')} | {r.get('regime')} | {r.get('T_hs')} | "
                     f"{r.get('screening_fraction')} | {r.get('g2_op')} | {r.get('collected_flux_pulsed_s')} | "
                     f"{r.get('valid')} | {r.get('row_id')} |")
    lines.append("")
    cplane_300_screened = _pick(planar_rows, family="c_plane", regime="deterministic_pair",
                                 T_hs=300.0, screening_fraction=1.0)
    all_core = [r for r in core if r.get("row_kind") == "core"]
    nanowire_flux = sorted(float(r["collected_flux_pulsed_s"]) for r in all_core if r.get("optical_pass")
                            and _finite_num(r.get("collected_flux_pulsed_s")))
    lines.append("Reconciliation: planar c-plane SET, screening=1, 300 K (row "
                 f"{cplane_300_screened.get('row_id') if cplane_300_screened else 'n/a'}): "
                 f"{cplane_300_screened.get('collected_flux_pulsed_s') if cplane_300_screened else 'n/a'} /s "
                 "(round 2's own optimized headline was a DIFFERENT, separately optimized design point ~28 "
                 "kHz/s, not this screened-default card row -- the two must not be conflated) vs this run's "
                 "nanowire optical_pass collected_flux_pulsed_s range "
                 + (f"{nanowire_flux[0]:.4g}-{nanowire_flux[-1]:.4g} /s" if nanowire_flux else "n/a")
                 + ". The nanowire tier's higher flux is attributed to: a single-wire supply without the "
                 "planar aperture partition, relaxed (piezoelectric-field-free) strain raising overlap and "
                 "suppressing surface escape, and antenna/waveguide collection geometry differing from the "
                 "planar cavity -- qualitative mechanism attributions, not a controlled one-parameter "
                 "comparison between the two platforms.")
    lines.append("")
    return lines


def _sensitivity_ranking_table(core, all_core):
    """M9-M10/L obligation: a numerical sensitivity table for every reduced
    cut, ranked by MEASURED leverage (never a single narrative-ordered
    'expected ranking' sentence, which mixed commanded flux, delivered
    flux, g2 and an unsampled tau_rad0_ns axis together) -- three SEPARATE
    rankings (commanded flux, delivered flux, g2_op), reported PER FAMILY
    at that family's own reference geometry."""
    lines = ["## Numerical sensitivity table, ranked by MEASURED leverage (M9-M10/L obligation)", "",
             "Ranking is measured from this run's own rows, separately per family and per observable -- "
             "never a single mixed narrative ordering.", ""]
    sens = [r for r in core if r.get("row_kind") == "sensitivity"]
    for fam in FAMILIES:
        ref = REF_GEOM[fam]
        base_filter = dict(family=fam, regime="deterministic_pair", strain_bound="relaxed",
                            T_hs=300.0, rep_rate_hz=200.0e6, core_radius_nm=ref["core_radius_nm"],
                            height_nm=ref["height_nm"], x_in=ref["x_in"])
        ref_row = _pick(all_core, **base_filter)
        ref_flux = (float(ref_row["collected_flux_pulsed_s"])
                    if ref_row and _finite_num(ref_row.get("collected_flux_pulsed_s")) else None)
        ref_flux_d = (float(ref_row["collected_flux_delivered_s"])
                      if ref_row and _finite_num(ref_row.get("collected_flux_delivered_s")) else None)
        ref_g2 = float(ref_row["g2_op"]) if ref_row and _finite_num(ref_row.get("g2_op")) else None

        lines.append(f"### {fam}")
        lines.append("")
        if ref_row:
            lines.append(f"Reference row: {ref_row.get('row_id')} g2_op={ref_row.get('g2_op')} "
                         f"commanded_flux={ref_row.get('collected_flux_pulsed_s')} "
                         f"delivered_flux={ref_row.get('collected_flux_delivered_s')}")
        else:
            lines.append("Reference row: not in this run's coverage.")
        lines.append("")

        axes_present = sorted({r.get("sensitivity_axis") for r in sens
                                if r.get("family") == fam
                                and r.get("sensitivity_axis") not in ("", "deshpande2013_comparison", "current_pulse_width")})
        rows_by_axis = {}
        leverage = {}
        for axis in axes_present:
            rows_axis = [r for r in sens if r.get("sensitivity_axis") == axis and r.get("family") == fam
                         and r.get("regime") == "deterministic_pair" and r.get("strain_bound") == "relaxed"
                         and _close(float(r.get("T_hs", -1)), 300.0) and _close(float(r.get("rep_rate_hz", -1)), 200.0e6)]
            rows_by_axis[axis] = rows_axis
            lev_flux = lev_flux_d = lev_g2 = 0.0
            for r in rows_axis:
                flux, fluxd, g2v = (r.get("collected_flux_pulsed_s"), r.get("collected_flux_delivered_s"),
                                     r.get("g2_op"))
                if ref_flux and _finite_num(flux):
                    lev_flux = max(lev_flux, abs(float(flux) / ref_flux - 1.0))
                if ref_flux_d and _finite_num(fluxd):
                    lev_flux_d = max(lev_flux_d, abs(float(fluxd) / ref_flux_d - 1.0))
                if ref_g2 is not None and _finite_num(g2v):
                    lev_g2 = max(lev_g2, abs(float(g2v) - ref_g2))
            leverage[axis] = {"flux": lev_flux, "flux_delivered": lev_flux_d, "g2": lev_g2}

        for label, key in (("commanded flux |ratio-1| (measured leverage)", "flux"),
                            ("delivered flux |ratio-1| (measured leverage)", "flux_delivered"),
                            ("g2_op |delta| (measured leverage)", "g2")):
            ranked = sorted(axes_present, key=lambda a: -leverage[a][key])
            lines.append(f"Ranked by {label}: " + ", ".join(f"{a}({leverage[a][key]:.3g})" for a in ranked)
                         if ranked else f"Ranked by {label}: (no reduced-cut axes this family this run).")
        lines.append("")
        lines.append("(the contact R_s_ohm row has ratio~1.0 on commanded flux -- it only changes delivered "
                     "flux via the RC time constant, never the idealized/commanded flux.)"
                     if "R_s_ohm" in axes_present else "")
        lines.append("")

        lines.append("| axis | value | g2_op | commanded_flux/s | flux_ratio_to_reference | delivered_flux/s | "
                     "delivered_flux_ratio_to_reference | one_pair_valid | optical_pass | quality_pass | row_id |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        if ref_row:
            lines.append(f"| (reference) | main-grid default | {ref_row.get('g2_op')} | "
                         f"{ref_row.get('collected_flux_pulsed_s')} | 1.0 | "
                         f"{ref_row.get('collected_flux_delivered_s')} | 1.0 | {ref_row.get('one_pair_valid')} | "
                         f"{ref_row.get('optical_pass')} | {ref_row.get('quality_pass')} | {ref_row.get('row_id')} |")
        for axis in sorted(axes_present, key=lambda a: -leverage[a]["flux"]):
            for r in sorted(rows_by_axis[axis], key=lambda z: str(z.get("sensitivity_value"))):
                flux, fluxd = r.get("collected_flux_pulsed_s"), r.get("collected_flux_delivered_s")
                ratio = (f"{float(flux) / ref_flux:.4g}" if ref_flux and _finite_num(flux) else "n/a")
                ratio_d = (f"{float(fluxd) / ref_flux_d:.4g}" if ref_flux_d and _finite_num(fluxd) else "n/a")
                lines.append(f"| {axis} | {r.get('sensitivity_value')} | {r.get('g2_op')} | {flux} | {ratio} | "
                             f"{fluxd} | {ratio_d} | {r.get('one_pair_valid')} | {r.get('optical_pass')} | "
                             f"{r.get('quality_pass')} | {r.get('row_id')} |")
        lines.append("")
    lines.append("tau_rad0_ns/dipole-prior leverage proxy (from the c-plane dipole-prior falsification rows "
                 "above): switching dipole_weights from isotropic to CPLANE_ONLY materially changes "
                 "antenna_rate_factor/tau_rad_photonic_ns/commanded flux (see that section's row values) -- "
                 "a qualitative leverage indicator, not a numeric ranking-table entry, since it sweeps an "
                 "orientation prior rather than the tau_rad0_ns magnitude (not sampled this run).")
    lines.append("")
    lines.append("Restored current/pulse-width cut (I x tau_pulse, 9 rows, horizontal reference geometry, "
                 "rectangular regime, relaxed, 300K, 200MHz -- preserves the 100 ps headline point):")
    lines.append("")
    lines.append("| I_uA | tau_pulse_ns | g2_op | commanded_flux/s | valid | row_id |")
    lines.append("|---|---|---|---|---|---|")
    pulse_rows = [r for r in sens if r.get("sensitivity_axis") == "current_pulse_width"]
    for r in sorted(pulse_rows, key=lambda z: (float(z.get("I_uA", 0)), float(z.get("tau_pulse_ns", 0)))):
        lines.append(f"| {r.get('I_uA')} | {r.get('tau_pulse_ns')} | {r.get('g2_op')} | "
                     f"{r.get('collected_flux_pulsed_s')} | {r.get('valid')} | {r.get('row_id')} |")
    lines.append("")
    return lines


def _definitions_section():
    """L17 obligation: define `eligible` and mention `device_pass`."""
    return ["## Column/gate definitions (L17 obligation)", "",
            "`eligible` = valid AND finite g2_op AND finite collected_flux_pulsed_s AND "
            "collected_flux_pulsed_s>=1000/s (flux_floor) -- NOT the same as `headline_eligible` (the "
            "vertical_photonic single_mode AND approximation_error==0 check, bullet 10) or `optical_pass` "
            "(which additionally requires g2_op<0.5 and, for deterministic_pair rows, one_pair_valid AND "
            "pair_supply_possible). `device_pass` (a device-module output column, sweep.csv) mirrors the "
            "planar convention's Coulomb-screen alias; it is reported as a row column alongside "
            "`set_feasible`/`hardware_qualified` but this results text does not gate any VERDICT/BEST/"
            "quality_pass computation on it beyond the `set_feasible`/`hardware_qualified` values already "
            "discussed above.", ""]


def _fitted_inputs_disclosure():
    """L14 obligation: disclose the transport module's fitted electrical/
    thermal inputs as distinct from the unfitted optical prediction."""
    return ["## Fitted electrical/thermal transport inputs (L14 obligation)", "",
            "Distinct from the unfitted optical prediction (levels/photonics/surface: no parameter there is "
            "chosen to reproduce a held-out optical anchor): fsim_core/nitride_nanowire_transport.py's "
            "NitrideWireDiode.tau_SRH_ns default (0.01 ns, 10 ps) is [A, back-solved] so the GaN-kernel dark "
            "current reproduces the empirically inferred 2.6-3.1 V junction-voltage window at 1 nA/300 K "
            "(Deshpande 2013 terminal range minus the IR drop) -- used in EVERY row's V_j solve this run, "
            "not a measured SRH lifetime. WIRE_RTH_PRESETS_K_W['deshpande_fig4_replay']=2.8e9 K/W is [DR, "
            "re-fit] to both Deshpande 2013 Fig.4 heating rises under the H2 GaN-kernel V_j, at the paper's "
            "10 K bath -- a SEPARATE replay anchor from this sweep's own Rth_K_W card defaults (1.0e9 "
            "horizontal, 1.0e7 vertical, see the Card schema), which are NOT re-fit to Fig.4 and are the "
            "values this sweep's core/sensitivity rows actually use.", ""]


def _results_md(core, quick, complete, runtime_s, evaluate_calls, invalid_by_kind, dipole_rows, decisions):
    all_core = [r for r in core if r.get("row_kind") == "core"]
    lines = ["# Nitride nanowire sweep results", "",
              "Independent predictions (design brief User answer 2): every number below is an "
              "evaluate() output or a declared, traceable reduction of one (row_id in sweep.csv); "
              "no plotted point or headline number is fitted to the Deshpande 2013/2014 held-out "
              "lifetime/g2/wavelength anchors. `family` is horizontal_as_built (Deshpande 2013/2014 "
              "as-built dispersed wire) or vertical_photonic (designed HE11 wire, headline-eligible "
              "only when single_mode AND approximation_error==0); `regime` is rectangular (100 ps "
              "electrical pulse) or deterministic_pair (idealized one-pair-per-cycle SET loading); "
              "`strain_bound` unrelaxed is the conservative_lower scenario, relaxed the "
              "headline_upper scenario -- neither is asserted a rigorous flux bound, and a reversal "
              "is reported (bound_reversal_pair) rather than hidden. Two independent, non-gating "
              "hardware screens (Coulomb-blockade `set_feasible`/`set_EC_over_kT` and resonant-"
              "tunnelling-injector `rti_feasible`) are reported alongside the idealized optical "
              "statistics; neither overwrites `g2_op`/`collected_flux_pulsed_s`.", ""]

    lines.append("## VERDICT lines (16: family x regime x strain_bound x rep_rate_hz)")
    lines.append("")
    vlines, verdicts = verdict_lines(core, quick)
    lines += vlines
    lines.append("")
    lines.append("## Best passing flux per family")
    lines.append("")
    lines += best_passing_flux_lines(core)
    lines.append("")

    lines.append("## Bounds table (conservative unrelaxed vs headline relaxed, reference geometry)")
    lines.append("")
    lines += bounds_table(core)
    lines.append("")
    # M11 fix: print the comparable-pair count/legend and the reversal
    # mechanism across the FULL core grid (not just the 16-row reference
    # table above), independently tallied from each relaxed row's own
    # bound_reversal_pair field (attach_bound_reversal already computed it).
    relaxed_rows = [r for r in all_core if r.get("strain_bound") == "relaxed"]
    n_reversed = sum(1 for r in relaxed_rows if r.get("bound_reversal_pair") == "True")
    n_not_reversed = sum(1 for r in relaxed_rows if r.get("bound_reversal_pair") == "False")
    n_not_comparable = sum(1 for r in relaxed_rows if r.get("bound_reversal_pair") == "not_comparable")
    mech_counts = {"mu": 0, "collected_flux_pulsed_s": 0, "g2_op": 0}
    for r in relaxed_rows:
        if r.get("bound_reversal_pair") != "True":
            continue
        partner = _bound_partner(r, all_core)
        if partner is None:
            continue
        for key, higher_is_headline in (("mu", True), ("collected_flux_pulsed_s", True), ("g2_op", False)):
            rv, uv = r.get(key), partner.get(key)
            if not (_finite_num(rv) and _finite_num(uv)):
                continue
            cond = float(rv) < float(uv) if higher_is_headline else float(rv) > float(uv)
            if cond:
                mech_counts[key] += 1
    lines.append(f"Comparable pairs (relaxed count): {n_reversed + n_not_reversed}/{len(relaxed_rows)} "
                 f"({n_reversed} reversed, {n_not_reversed} not reversed); not_comparable="
                 f"{n_not_comparable}/{len(relaxed_rows)} (no pair where BOTH strain-bound rows clear the "
                 f"1000/s optical floor -- the legend value for every 'not_comparable' cell/entry elsewhere "
                 f"in this report and in strain_reversal_map.png). Reversal-trigger tally among the "
                 f"{n_reversed} reversed pairs: mu={mech_counts['mu']}, "
                 f"collected_flux_pulsed_s={mech_counts['collected_flux_pulsed_s']}, g2_op={mech_counts['g2_op']} -- "
                 + ("every reversed pair triggers on g2_op ONLY, never on mu or collected_flux_pulsed_s."
                    if n_reversed > 0 and mech_counts["g2_op"] == n_reversed and mech_counts["mu"] == 0
                    and mech_counts["collected_flux_pulsed_s"] == 0
                    else "mixed triggers this run; see bound_reversal_pair per row in sweep.csv."))
    lines.append("")

    # ---- 200 MHz vs 80 MHz SET one_pair_valid (results obligation). Both
    # T_hs in {230,300} K are reported (not just 300 K): the two families
    # behave ASYMMETRICALLY at 80 MHz -- vertical_photonic's disc-in-wire
    # geometry clears one_pair_valid at 80 MHz at BOTH temperatures (a real
    # optical_pass). L16 fix: the prior "~60 MHz" claim named an off-grid
    # rate never evaluated this run and is REMOVED (no row backs it); the
    # per-row loop below is also de-templated -- it now reports each row's
    # OWN blocked_load_probability (a real number) instead of repeating an
    # identical boilerplate sentence four times.
    lines.append("## Repetition-rate sensitivity: 80 MHz vs 200 MHz SET one_pair_valid (family-asymmetric)")
    lines.append("")
    lines.append("one_pair_valid requires blocked_load_probability<=1e-9 (the disc must empty between "
                 "cycles); a shorter 200 MHz period leaves less time per cycle for that reset than 80 MHz, "
                 "so a one_pair_valid failure at 200 MHz where it holds at 80 MHz is a loading-window/period "
                 "effect, not a fit. collected_flux_pulsed_s is the idealized/commanded flux, "
                 "collected_flux_delivered_s is the RC-limited delivered flux (bullet 6). No row at any "
                 "OFF-GRID repetition rate was evaluated this run -- only this run's own main-grid values, "
                 "80 MHz and 200 MHz, were ever sampled.")
    lines.append("")
    for family in FAMILIES:
        for t_hs in (230.0, 300.0):
            pair = _rate_pair_rows(all_core, family, T_hs=t_hs)
            if 80.0e6 in pair and 200.0e6 in pair:
                r80, r200 = pair[80.0e6], pair[200.0e6]
                b80 = r80.get("blocked_load_probability"); b200 = r200.get("blocked_load_probability")
                lines.append(f"{family} T_hs={t_hs:g}K: at 80 MHz one_pair_valid={r80.get('one_pair_valid')} "
                              f"optical_pass={_optical_pass_of(r80)} blocked_load_probability={b80} "
                              f"(row {r80.get('row_id')}, g2={r80.get('g2_op')}, "
                              f"flux={r80.get('collected_flux_pulsed_s')}, delivered={r80.get('collected_flux_delivered_s')}); "
                              f"at 200 MHz one_pair_valid={r200.get('one_pair_valid')} optical_pass={_optical_pass_of(r200)} "
                              f"blocked_load_probability={b200} "
                              f"(row {r200.get('row_id')}, g2={r200.get('g2_op')}, flux={r200.get('collected_flux_pulsed_s')}, "
                              f"delivered={r200.get('collected_flux_delivered_s')}).")
            else:
                lines.append(f"{family} T_hs={t_hs:g}K: rate-pair rows not found in this run's coverage (quick={quick}).")
    lines.append("")

    # ---- Deshpande 2014 comparison (300 K)
    lines.append("## Deshpande 2014 comparison (x_in=0.40, T_hs=300 K, 200 MHz, deterministic_pair)")
    lines.append("")
    ref = REF_GEOM["horizontal_as_built"]
    d2014 = {sb: _group(all_core, family="horizontal_as_built", regime="deterministic_pair",
                          rep_rate_hz=200.0e6, T_hs=300.0, strain_bound=sb,
                          core_radius_nm=ref["core_radius_nm"], height_nm=ref["height_nm"], x_in=0.40)
             for sb in STRAIN_BOUNDS}
    lines.append("Measured (abstract-only [V], CONDITIONS INCOMPLETE): lambda~630 nm, lifetime 1.3+/-0.3 ns, g2=0.29.")
    lines.append("")
    lines.append("| strain_bound | lambda_nm (predicted) | tau_rad_bare_ns | tau_rad_photonic_ns | tau_total_X_ns | g2_op | collected_flux_pulsed_s | row_id |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for sb in STRAIN_BOUNDS:
        rr = d2014.get(sb)
        if rr:
            r = rr[0]
            lines.append(f"| {sb} | {r.get('lambda_nm')} | {r.get('tau_rad_bare_ns')} | {r.get('tau_rad_photonic_ns')} | "
                          f"{r.get('tau_total_X_ns')} | {r.get('g2_op')} | {r.get('collected_flux_pulsed_s')} | {r.get('row_id')} |")
        else:
            lines.append(f"| {sb} | not in this run's coverage | | | | | | |")
    old_replay = [r for r in core if r.get("row_kind") == "planar_2014_replay"]
    if old_replay:
        r = old_replay[0]
        lines.append("")
        lines.append(f"Untouched OLD planar 2014 comparison replay (non-gating, unchanged from the prior "
                      f"piece): row {r.get('row_id')}, g2={r.get('g2_op')}, flux={r.get('collected_flux_pulsed_s')} "
                      f"vs measured g2=0.29.")
    lines.append("")

    # ---- Deshpande 2013 comparison (10 K)
    lines.append("## Deshpande 2013 comparison (pulsed replay of a CW measurement, 10 K, drive_mismatch)")
    lines.append("")
    lines.append("Measured (CW electrical, 10 K, [V]): X g2 raw/corrected 0.30/0.16, XX 0.38/0.25 at 1 nA; "
                  "g2-fit lifetimes X 1.1 ns, XX 0.7 ns; TRPL XX 711 ps; emission X=2.84 eV (436.56 nm). This "
                  "evaluator has no CW/HBT drive path (pulsed rectangular / deterministic_pair only); the rows "
                  "below are the pulsed model's own prediction at the paper's geometry/current, published as "
                  "drive_mismatch -- never substituted for a predicted CW g2.")
    lines.append("")
    d2013 = [r for r in core if r.get("row_kind") == "sensitivity" and r.get("sensitivity_axis") == "deshpande2013_comparison"]
    # L13 fix: print bare, photonic, AND total lifetimes together (the prior
    # table printed tau_rad_bare_ns alone, inviting a false "matches the 1.1
    # ns anchor" read); bare's relation to the [A] tau_rad0_ns input is
    # stated as a RATIO computed from the row itself, never hardcoded.
    lines.append("| R nm | I uA | lambda_nm (predicted) | tau_rad_bare_ns | tau_rad_photonic_ns | tau_total_X_ns | "
                 "g2_op (drive_mismatch) | commanded_flux/s | time_avg_current_pA | valid | row_id |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(d2013, key=lambda z: (float(z.get("core_radius_nm", 0)), float(z.get("I_uA", 0)))):
        i_ua, tau_ns, rep = r.get("I_uA"), r.get("tau_pulse_ns"), r.get("rep_rate_hz")
        if _finite_num(i_ua) and _finite_num(tau_ns) and _finite_num(rep):
            duty = float(tau_ns) * 1e-9 * float(rep)
            time_avg_pA = f"{float(i_ua) * 1e6 * duty:.4g}"
        else:
            time_avg_pA = "n/a"
        lines.append(f"| {r.get('core_radius_nm')} | {r.get('I_uA')} | {r.get('lambda_nm')} | {r.get('tau_rad_bare_ns')} | "
                      f"{r.get('tau_rad_photonic_ns')} | {r.get('tau_total_X_ns')} | "
                      f"{r.get('g2_op')} (drive_mismatch) | {r.get('collected_flux_pulsed_s')} | {time_avg_pA} | "
                      f"{r.get('valid')} | {r.get('row_id')} |")
    lines.append("")
    d2013_bare_ratio = [(float(r["tau_rad_bare_ns"]) / float(r["tau_rad0_ns"]))
                         for r in d2013 if _finite_num(r.get("tau_rad_bare_ns")) and _finite_num(r.get("tau_rad0_ns"))
                         and float(r["tau_rad0_ns"]) != 0.0]
    if d2013_bare_ratio:
        lines.append(f"tau_rad_bare_ns / tau_rad0_ns (this run's own [A] input) ratio: "
                     f"{min(d2013_bare_ratio):.4g}-{max(d2013_bare_ratio):.4g} across these rows -- the bare "
                     "radiative lifetime's closeness to the 1.1 ns anchor reflects this fixed [A] input times "
                     "a near-unity geometry factor, NOT an independently verified radiative-rate prediction; "
                     "tau_rad_photonic_ns (after antenna suppression) and tau_total_X_ns (after nonradiative/"
                     "surface competition) are the physically relevant, larger, device-level lifetimes and are "
                     "reported alongside it, never substituted for it.")
    lines.append("Drive AND heating mismatch: these rows run 100 ps rectangular pulses at 80 MHz (this run's "
                 "own I_uA/tau_pulse_ns/rep_rate_hz columns above give the time-average current column), "
                 "roughly 2 orders of magnitude below the paper's CW 1 nA -- both the counting statistics "
                 "(pulsed vs CW g2) and the junction/thermal operating point (heating scales with time-average "
                 "power) differ from the measured device; g2_op above is published as drive_mismatch, never a "
                 "predicted CW g2.")
    lines.append("")
    # Label invalid rows MISSING with their reason, never silently repaired
    # (spec: "Label invalid low-T material/transport rows missing, not
    # silently repaired"). The actual reason(s) are read from each row's own
    # invalid_reasons below -- never assumed/hardcoded here (H4 fix: an
    # earlier draft of this comment wrongly assumed the Si-substrate table
    # cause; the table is actually 380-750 nm, see _SI_INDEX_ANCHORS_NM, and
    # ~446 nm sits INSIDE that range, so if these rows are invalid it is for
    # a different, row-reported reason).
    invalid_d2013 = [r for r in d2013 if not r.get("valid")]
    if invalid_d2013:
        seen_reasons = set()
        for r in invalid_d2013:
            try:
                for reason in json.loads(r.get("invalid_reasons", "[]")):
                    seen_reasons.add(reason)
            except json.JSONDecodeError:
                pass
        lines.append(f"MISSING (not silently repaired): {len(invalid_d2013)}/{len(d2013)} 2013 comparison rows are "
                      f"invalid this run. Reason(s): " + ("; ".join(sorted(seen_reasons)) if seen_reasons else "unknown") + ".")
        lines.append("")

    # ---- opposite-endpoint anchor match (results obligation)
    lines.append("## Opposite-endpoint anchor match (bullet 11 obligation)")
    lines.append("")
    r2013 = next((r for r in d2013 if _close(float(r.get("core_radius_nm", -1)), 12.5)), None)
    r2014_rel = (d2014.get("relaxed") or [None])[0]
    r2014_unrel = (d2014.get("unrelaxed") or [None])[0]
    lam2013 = r2013.get("lambda_nm") if r2013 else None
    lam2013_diag, lam2013_recovered = (_diagnostic_lambda_from_invalid_reasons(r2013) if r2013 else (None, False))
    lam2014_rel = r2014_rel.get("lambda_nm") if r2014_rel else None
    lam2014_unrel = r2014_unrel.get("lambda_nm") if r2014_unrel else None
    if lam2013_recovered:
        lam2013_text = (f"{lam2013_diag:.6g} (DIAGNOSTIC ONLY: the row itself is invalid -- its bare-dot "
                         f"emission wavelength, already computed upstream by the levels module, is recovered "
                         f"from the row's own invalid_reasons message, not a repaired/assumed value; the "
                         f"downstream photonics step fails because this wavelength falls outside the "
                         f"photonics module's own tabulated Si-substrate index range -- see the 'Deshpande "
                         f"2013 comparison' section's MISSING line for the row's actual reason string)")
    elif r2013 is not None and not r2013.get("valid"):
        try:
            _reasons_2013 = json.loads(r2013.get("invalid_reasons", "[]"))
        except (TypeError, json.JSONDecodeError):
            _reasons_2013 = []
        lam2013_text = (f"{lam2013} (row {r2013.get('row_id')} is INVALID this run for a reason OTHER than "
                         f"the Si-substrate index table range: "
                         + ("; ".join(_reasons_2013) if _reasons_2013 else "unknown") + ")")
    else:
        lam2013_text = f"{lam2013}"
    lines.append(f"2013 relaxed predicted lambda_nm={lam2013_text} vs measured ~437 nm (X=2.84 eV, [V] Fig. 3c). "
                 f"2014 relaxed predicted lambda_nm={lam2014_rel} / unrelaxed predicted lambda_nm={lam2014_unrel} "
                 f"vs measured ~630 nm. The two strain-bound anchors are matched by OPPOSITE endpoints (2013 by "
                 f"relaxed, 2014 by whichever endpoint lands closer) and are never averaged: one of x_in transfer, "
                 f"disc thickness, VBO/bowing, or lateral localization carries roughly 0.3 eV of the remaining "
                 f"discrepancy -- stated here, not resolved.")
    lines.append("")

    # ---- RC caveat (bullet 6)
    lines.append("## RC caveat (bullet 6 obligation)")
    lines.append("")
    rc_as_built = [r for r in all_core if r.get("family") == "horizontal_as_built" and r.get("regime") == "deterministic_pair"
                   and _finite_num(r.get("tau_RC_ns"))]
    rc_designed = [r for r in core if r.get("row_kind") == "sensitivity" and r.get("sensitivity_axis") == "R_s_ohm"
                   and r.get("family") == "horizontal_as_built" and _close(float(r.get("R_s_ohm", -1)), 1.0e6)
                   and _finite_num(r.get("tau_RC_ns"))]
    if rc_as_built:
        tau_rc = sorted(float(r["tau_RC_ns"]) for r in rc_as_built)
        dsf = sorted(float(r["delivered_step_fraction"]) for r in rc_as_built if _finite_num(r.get("delivered_step_fraction")))
        dsf_txt = f"{dsf[0]:.4g}-{dsf[-1]:.4g}" if dsf else "n/a"
        lines.append(f"As-built R_s_ohm=2.38e9 (horizontal_as_built, deterministic_pair rows): "
                      f"tau_RC_ns spans {tau_rc[0]:.4g}-{tau_rc[-1]:.4g}, delivered_step_fraction spans "
                      f"{dsf_txt} -- this device CANNOT deliver a 100 ps step.")
    if rc_designed:
        tau_rc_d = sorted(float(r["tau_RC_ns"]) for r in rc_designed)
        dsf_d = sorted(float(r["delivered_step_fraction"]) for r in rc_designed if _finite_num(r.get("delivered_step_fraction")))
        lines.append(f"Designed-contact R_s_ohm=1e6 (sensitivity rows): tau_RC_ns spans "
                      f"{tau_rc_d[0]:.4g}-{tau_rc_d[-1]:.4g}, delivered_step_fraction spans "
                      f"{dsf_d[0]:.4g}-{dsf_d[-1]:.4g} -- only this designed contact can deliver the 100 ps step.")
    if rc_as_built or rc_designed:
        lines.append("These numbers are THIS RUN's own tau_RC_ns/delivered_step_fraction, recomputed from "
                     "sweep.csv rows, not restated from the contract: docs/nitride_nanowire_contract.md "
                     "bullet 6 quotes tau_RC_ns 0.82-1.19 / delivered_step_fraction 0.08-0.12 for the as-built "
                     "device, which this run's own as-built range above SUPERSEDES (the contract text is out "
                     "of scope for this piece and is not edited here).")
    lines.append("")

    # ---- access-1.0 lifetime cap (bullet 8)
    lines.append("## Access-1.0 lifetime cap (bullet 8 obligation)")
    lines.append("")
    access1 = [r for r in core if r.get("row_kind") == "sensitivity" and r.get("sensitivity_axis") == "occupied_dot_access"
               and _close(float(r.get("sensitivity_value", -1)), 1.0) and r.get("family") == "horizontal_as_built"
               and _close(float(r.get("core_radius_nm", -1)), 12.5) and _finite_num(r.get("k_surface_X_ns"))]
    if access1:
        r = access1[0]
        kx, kxx = float(r["k_surface_X_ns"]), float(r["k_surface_XX_ns"])
        tau_x = 1.0 / kx if kx > 0 else float("inf")
        tau_xx = 1.0 / kxx if kxx > 0 else float("inf")
        lines.append(f"occupied_dot_access=1.0 at core_radius_nm=12.5 (row {r.get('row_id')}): "
                      f"1/k_surface_X_ns={tau_x:.4g} ns, 1/k_surface_XX_ns={tau_xx:.4g} ns "
                      "(contract's stated cap: 0.625 ns X, 0.3125 ns XX) -- this caps tau_X/tau_XX regardless of "
                      "any other physics; the headline default occupied_dot_access=0.05 is reported next to it, "
                      "never chosen to reproduce the held-out anchors.")
    else:
        lines.append("occupied_dot_access=1.0 reference row not found in this run's coverage.")
    lines.append("")

    # ---- E_C/kT wall (bullet 11 obligation)
    lines.append("## E_C/kT wall (bullet 11 obligation)")
    lines.append("")
    wall_rows = [r for r in all_core if r.get("regime") == "deterministic_pair" and _finite_num(r.get("core_radius_nm"))
                 and float(r["core_radius_nm"]) >= 10.0 and _finite_num(r.get("set_EC_over_kT"))]
    rti_margin_rows = [r for r in all_core if r.get("regime") == "deterministic_pair" and _finite_num(r.get("core_radius_nm"))
                        and float(r["core_radius_nm"]) >= 10.0 and _finite_num(r.get("rti_level_margin_kT"))]
    if wall_rows:
        ecs = [float(r["set_EC_over_kT"]) for r in wall_rows]
        n_pass_coulomb = sum(1 for r in wall_rows if r.get("set_feasible") is True)
        n_pass_rti = sum(1 for r in rti_margin_rows if r.get("rti_feasible") is True)
        lines.append(f"Across every deterministic_pair core row with core_radius_nm>=10 nm at 230-300 K: "
                      f"set_EC_over_kT spans {min(ecs):.4g}-{max(ecs):.4g} (required ec_margin=10), "
                      f"set_feasible=True count={n_pass_coulomb}/{len(wall_rows)}; rti_feasible=True count="
                      f"{n_pass_rti}/{len(rti_margin_rows)}. Deterministic loading at 230-300 K fails the "
                      "Coulomb-blockade screen for every core_radius_nm>=10 nm priced in this tier.")
        # M5 fix: rti_feasible=False is NOT itself a demonstrated hardware
        # result -- rti_status is unknown_incomplete on essentially every
        # SET row this run (the RT injector screen is a conditional
        # engineering screen, not a measured device; contract "Composition
        # rules for the device piece" bullet 7: "rti_feasible=False must not
        # be reported as a demonstrated physics result"). Counted from the
        # rows, never asserted.
        set_rows = [r for r in all_core if r.get("regime") == "deterministic_pair"]
        n_unknown_incomplete = sum(1 for r in set_rows if r.get("rti_status") == "unknown_incomplete")
        n_transport_infeasible = sum(1 for r in set_rows if r.get("rti_transport_feasible") is False)
        lines.append(f"The RT injector screen is SEPARATELY reported, not folded into the Coulomb wall above: "
                      f"rti_status=unknown_incomplete on {n_unknown_incomplete}/{len(set_rows)} deterministic_pair "
                      f"rows, rti_transport_feasible=False on {n_transport_infeasible}/{len(set_rows)} -- these "
                      "are conditional engineering screens on an unsupported occupation/second-pair control, "
                      "NOT a demonstrated hardware failure by either charging mechanism; deterministic loading "
                      "at 230-300 K fails on the Coulomb-blockade wall (set_feasible) and separately carries an "
                      "incomplete/unresolved RT-injector screen, and neither screen overwrites optical_pass/"
                      "g2_op/collected_flux_pulsed_s computed upstream of it.")
        # M6 fix (contract bullets 7/280 and 11/325): both priced
        # charging-based loading mechanisms share the SAME insufficient
        # disc E_C/kT -- the RT-injector screen's own
        # second_pair_addition_meV input IS set_E_C_meV, the identical
        # Coulomb charging energy the set_EC_over_kT wall above already
        # reports failing at every core_radius_nm>=10 nm priced here.
        # Stated as a CONDITIONAL MODEL LIMITATION shared by both screens'
        # inputs, never a demonstrated hardware failure of the RT-injector
        # mechanism itself.
        lines.append("Both priced charging-based loading mechanisms share the SAME insufficient disc "
                     "E_C/kT: the RT-injector screen's own second_pair_addition_meV input is set_E_C_meV -- "
                     "the identical Coulomb charging energy the set_EC_over_kT wall above already reports as "
                     "an E_C/kT<10 failure at every core_radius_nm>=10 nm priced here. This is a CONDITIONAL "
                     "MODEL LIMITATION shared by both screens' inputs (contract bullet 11's E_C/kT wall "
                     "statement covers 'any charging mechanism priced in this tier'), never a demonstrated "
                     "hardware failure of the RT-injector mechanism specifically (contract bullet 7: "
                     "rti_feasible=False must not be reported as a demonstrated physics result).")
    else:
        lines.append("No deterministic_pair core rows with finite set_EC_over_kT at core_radius_nm>=10 nm.")
    lines.append("")

    # ---- c-plane dipole prior falsification
    lines.append("## c-plane dipole prior falsification (Composition rules bullet 9)")
    lines.append("")
    if dipole_rows:
        by_tag = {r.get("sensitivity_value"): r for r in dipole_rows}
        iso = by_tag.get("isotropic_default"); cpl = by_tag.get("cplane_only_falsified")
        lines.append(f"isotropic default (row {iso.get('row_id') if iso else 'n/a'}): "
                      f"degree_of_linear_polarization={iso.get('degree_of_linear_polarization') if iso else 'n/a'} "
                      f"(anchor: +70% axial, deshpande2013_polarization [V]). "
                      f"CPLANE_ONLY_DIPOLE_WEIGHTS=(0,0.5,0.5) (row {cpl.get('row_id') if cpl else 'n/a'}): "
                      f"degree_of_linear_polarization={cpl.get('degree_of_linear_polarization') if cpl else 'n/a'} "
                      "-- the wrong sign against the +70% anchor, confirming this sensitivity stays a named, "
                      "explicitly falsified alternative and never the headline default.")
    else:
        lines.append("dipole falsification rows not evaluated this run.")
    lines.append("")

    lines += _gate_anti_monotonicity_section(core, all_core, quick)
    lines += _per_temperature_count_tables(all_core)
    lines += _headline_nominations_16(core, all_core)
    lines += _ensemble_yield_proxy_section(core)
    lines += _planar_reference_section(core)
    lines += _sensitivity_ranking_table(core, all_core)
    lines += _definitions_section()
    lines += _fitted_inputs_disclosure()

    lines.append("## Invalid rows (all kinds)")
    lines.append("")
    lines.append(f"invalid_total={sum(invalid_by_kind.values())} across kinds: "
                 + (", ".join(f"{k}={v}" for k, v in sorted(invalid_by_kind.items())) if invalid_by_kind else "none"))
    lines.append("")
    # Breakdown by REASON (not just row_kind): a single dependency-module
    # data-table gap can dominate the invalid count, and that is a finding
    # in its own right, not something to leave implicit in per-kind totals.
    reason_counts = {}
    for r in core:
        if r.get("valid"):
            continue
        try:
            reasons = json.loads(r.get("invalid_reasons", "[]"))
        except (TypeError, json.JSONDecodeError):
            reasons = ["<undecodable invalid_reasons>"]
        for reason in reasons:
            key = reason.split(":", 1)[0] if ":" in reason else reason
            reason_counts[key] = reason_counts.get(key, 0) + 1
    if reason_counts:
        lines.append("Breakdown by reason (leading token before ':', counted per row -- a row may carry "
                     "more than one reason):")
        lines.append("")
        for key, n in sorted(reason_counts.items(), key=lambda kv: -kv[1]):
            lines.append(f"- {key}: {n}")
        top_key, top_n = max(reason_counts.items(), key=lambda kv: kv[1])
        invalid_total = sum(invalid_by_kind.values())
        # H4 fix: this NOTE is now generated FROM THE ROWS (family/bound/
        # x_in/height composition, recovered wavelength range, and the
        # table range read from the photonics module's own
        # _SI_INDEX_ANCHORS_NM), never hardcoded prose -- the prior version
        # hardcoded "covers only 450-630 nm" / "x_in=0.25 ... blue-shifts
        # below 450 nm", both wrong (the real table is 380-750 nm and the
        # affected rows are x_in=0.40 h=3-4nm UNRELAXED rows RED-shifted
        # past 750 nm by QCSE).
        if invalid_total > 0 and top_n / invalid_total > 0.5:
            import fsim_core.nitride_nanowire_photonics as _photonics
            si_lo = float(_photonics._SI_INDEX_ANCHORS_NM[0])
            si_hi = float(_photonics._SI_INDEX_ANCHORS_NM[-1])
            si_rows = []
            for r in core:
                if r.get("row_kind") != "core" or r.get("valid"):
                    continue
                try:
                    reasons = json.loads(r.get("invalid_reasons", "[]"))
                except (TypeError, json.JSONDecodeError):
                    reasons = []
                if not any(reason.startswith(top_key) for reason in reasons):
                    continue
                lam, recovered = _diagnostic_lambda_from_invalid_reasons(r)
                if recovered:
                    si_rows.append((r, lam))
            lines.append("")
            if si_rows:
                lambdas = [lam for _, lam in si_rows]
                families = sorted({r.get("family") for r, _ in si_rows})
                bounds = sorted({r.get("strain_bound") for r, _ in si_rows})
                x_ins = sorted({float(r.get("x_in", -1)) for r, _ in si_rows})
                heights = sorted({float(r.get("height_nm", -1)) for r, _ in si_rows})
                direction = ("red-shifted ABOVE" if min(lambdas) > si_hi
                             else ("blue-shifted BELOW" if max(lambdas) < si_lo else "outside"))
                lines.append(f"NOTE: '{top_key}' alone accounts for {top_n}/{invalid_total} "
                             f"({100.0 * top_n / invalid_total:.0f} percent) of all invalid rows this run -- a "
                             "single dependency-module data-table coverage gap, generated from the rows "
                             "themselves: fsim_core/nitride_nanowire_photonics.py's Si-substrate complex-index "
                             f"table covers {si_lo:g}-{si_hi:g} nm (read from the module's own "
                             "_SI_INDEX_ANCHORS_NM table). The affected core rows (recovered from their own "
                             f"invalid_reasons diagnostic, never a repaired/assumed value) span "
                             f"lambda_nm={min(lambdas):.6g}-{max(lambdas):.6g} nm ({direction} the table) and "
                             f"are composed of family={families}, strain_bound={bounds}, x_in={x_ins}, "
                             f"height_nm={heights} -- not a physics failure this sweep introduces or can fix "
                             "within its own three scoped files. Reported as the interface coordination item "
                             "for the orchestrator/reviewer, not silently repaired (no n_substrate override was "
                             "invented to force these rows valid).")
                lines.append("")
                lines.append("### QCSE excursion physics paragraph (H4 obligation)")
                lines.append("")
                # H2 fix: pick the representative row DETERMINISTICALLY --
                # the smallest row_id among the affected rows (documented
                # here, never the previous `si_rows[0][0]` with no stated
                # rule -- si_rows already iterates in row-id-ascending
                # order since `core` is built that way, so this min() is a
                # no-op on today's data but is now an explicit, checkable
                # rule rather than an accident of iteration order).
                sample = min((r for r, _ in si_rows), key=lambda r: r["row_id"])
                depletion_txt = (f"{float(sample.get('field_kVcm')):.4g} kV/cm"
                                  if _finite_num(sample.get("field_kVcm")) else "n/a")
                # H2 fix: field_kVcm on an INVALID row is the transport
                # DEPLETION field alone (fsim_core/nitride_nanowire_
                # device.py's exception-path fallback), NEVER the built-in
                # piezoelectric/polarization field -- the actual
                # polarization/total field is obtained from a separate
                # levels-only replay (_qcse_levels_replay), never by
                # relabelling field_kVcm.
                lv_replay = _qcse_levels_replay(sample)
                if lv_replay is not None:
                    replay_txt = (f"a levels-only replay (fsim_core.nitride_nanowire_levels, unrelaxed bound, "
                                  f"this row's own x_in/height_nm/core_radius_nm, external_field_kVcm set to "
                                  f"this row's OWN recorded field_kVcm as the resolved depletion field, per "
                                  f"device.py's own post-feedback convention) gives F_pz_kVcm="
                                  f"{lv_replay.F_pz_kVcm:.6g}, total_field_kVcm={lv_replay.field_kVcm:.6g}")
                else:
                    replay_txt = "the levels-only replay could not be computed for this row this run"
                partner_stats = {"valid": 0, "invalid": 0, "missing": 0}
                for r, _ in si_rows:
                    partner = _bound_partner(r, all_core)
                    if partner is None:
                        partner_stats["missing"] += 1
                    elif partner.get("valid"):
                        partner_stats["valid"] += 1
                    else:
                        partner_stats["invalid"] += 1
                lines.append(f"The {len(si_rows)} rows above are the unrelaxed (conservative_lower), "
                             f"x_in=0.40, height_nm in {heights} disc geometries. field_kVcm on an INVALID row "
                             "like these is the transport DEPLETION field alone (fsim_core/nitride_nanowire_"
                             "device.py's invalid-row fallback captures inj['depletion_field_kVcm'] before the "
                             "row goes invalid), NOT the built-in piezoelectric field -- this run's own "
                             f"field_kVcm on the representative row (deterministic rule: smallest row_id among "
                             f"the {len(si_rows)} affected rows), {sample.get('row_id')}: {depletion_txt} "
                             f"depletion field. The actual polarization field driving QCSE is reported "
                             f"separately: {replay_txt}. This field drives the quantum-confined Stark effect "
                             f"(QCSE) far enough to red-shift emission to "
                             f"{min(lambdas)/1000.0:.3g}-{max(lambdas)/1000.0:.3g} um, past the photonics "
                             f"module's own {si_hi:g} nm table ceiling -- an excursion this run reports as "
                             f"MISSING (invalid), never silently repaired. Their relaxed strain-bound partner "
                             f"(strain_fraction=0, no piezoelectric field) is valid for "
                             f"{partner_stats['valid']}/{len(si_rows)} of these rows and invalid/absent for "
                             f"{partner_stats['invalid'] + partner_stats['missing']}/{len(si_rows)} -- i.e. the "
                             "conservative (unrelaxed) bound is ABSENT at exactly these x_in=0.40, thick-disc "
                             "geometries where only the relaxed bound is computable this run, an asymmetry "
                             "the bounds table and headline selection must not paper over by silently "
                             "reporting only the relaxed side.")
            else:
                lines.append(f"NOTE: '{top_key}' accounts for {top_n}/{invalid_total} "
                             f"({100.0 * top_n / invalid_total:.0f} percent) of all invalid rows this run; no "
                             "row's invalid_reasons carried a recoverable diagnostic wavelength for this "
                             "specific reason this run.")
        lines.append("")

    if decisions:
        lines.append("## Decisions (conservative choices under ambiguity)")
        lines.append("")
        for d in decisions:
            lines.append(f"- {d}")
        lines.append("")

    lines.append("## Reduced-cut coverage (budget trim)")
    lines.append("")
    lines.append("evaluate() measured 6.1-9.8 s/call serially on this machine (fsim_core.nitride_nanowire_"
                 "injector.injector_feasibility's resonance search, run on every row regardless of regime, "
                 "dominates); the mandated 2560-row core grid alone consumes most of the 1800 s/40-minute "
                 "operational budget even parallelized across all cores (see manifest n_workers). The "
                 "reduced-cut axis list is therefore trimmed to screening_fraction {0,1}, occupied_dot_"
                 "access's 1.0 conservative partner, S_cm_s {1e2,1e4}, shell AlGaN, b_res 0.02, the "
                 "spec-named injector cuts (barrier thickness, al_fraction, alignment_uncertainty_meV, "
                 "growth_tolerance_steps, occupation-control uncertainty), and R_s_ohm=1e6 on the "
                 "horizontal family only -- each sampling only its non-default alternative value(s). "
                 "DROPPED entirely this run (a real, reported coverage gap against the contract's fuller "
                 "reduced-cut list, never silently omitted): " + ", ".join(DROPPED_CUT_AXES) + ".")
    lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append("Reduced cuts are one-at-a-time and moderately sampled (not exhaustive), always at the "
                 "family's reference geometry (12.5/80 nm core, 2 nm height, x_in=0.40), both regimes/bounds, "
                 "T in {230,300} K, both rates -- never the Cartesian product. bound_reversal_pair is a "
                 "declared post-hoc transform over already-evaluated row pairs (see attach_bound_reversal's "
                 "docstring), not a separate evaluate_strain_pair() call. No fitting to any Deshpande anchor "
                 "occurred anywhere in this sweep. Evaluate calls are dispatched across a process pool "
                 "(N_WORKERS, manifest.json) for wall-clock only -- every row's evaluate() input/output is "
                 "identical to a serial run.")
    lines.append("")
    lines.append(f"runtime_s={runtime_s:.1f} evaluate_calls={evaluate_calls} complete={complete}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------- shared figure/results
def _build_figures_and_results(out, core, quick, complete, runtime_s, evaluate_calls,
                                 invalid_by_kind, dipole_rows, decisions):
    """Shared between the full/quick evaluate path and --report-only mode:
    builds every figure from `core` and writes results.md. NEVER touches
    sweep.csv (the caller owns that). Kept as a single function so the two
    callers cannot silently diverge in which figures/captions get produced."""
    plots = {}
    plot_fail_log = []
    core_only = [r for r in core if r.get("row_kind") == "core"]
    # rep_rate_hz is a nonswept coordinate for all three of these traces and
    # MUST be pinned (spec: "Traces must fix all nonswept coordinates; no
    # joining different radii, temperatures or regimes") -- without it a
    # single (family, strain_bound) trace silently interleaved 80 MHz and
    # 200 MHz rows at the same x value, merging two different fixed
    # coordinates into one line (caught by verify_nitride_nanowire_sweep.py's
    # trace_shares_fixed_coords check). Pinned at 200 MHz, the Deshpande
    # 2014-comparison headline rate.
    fixed_common = {"regime": "deterministic_pair", "x_in": 0.40, "rep_rate_hz": 200.0e6}
    fixed_common_t300 = dict(fixed_common, T_hs=300.0)

    figs = [
        ("g2_flux_vs_core_radius.png", plot_vs_axis,
         (out, "g2_flux_vs_core_radius.png", core_only, "core_radius_nm", ("g2_op", "collected_flux_pulsed_s"),
          "g2 and flux vs core radius, both strain bounds (T_hs=300K, x_in=0.40, rep_rate_hz=200MHz, SET regime)", fixed_common_t300),
         dict(caption="Opposite endpoints: 2013/2014 anchors are matched by opposite strain bounds, never averaged; "
                       "'x' points are vertical_photonic rows with single_mode=False (headline_eligible=False, "
                       "never a nominated headline); open red circles are one_pair_valid=False points.", mark_multimode=True)),
        ("g2_flux_vs_disc_thickness.png", plot_vs_axis,
         (out, "g2_flux_vs_disc_thickness.png", core_only, "height_nm", ("g2_op", "collected_flux_pulsed_s"),
          "g2 and flux vs disc thickness, both strain bounds (T_hs=300K, x_in=0.40, rep_rate_hz=200MHz, SET regime)", fixed_common_t300),
         dict(caption="Access cap: occupied_dot_access=1.0 (not plotted here, default 0.05 shown) caps tau_X/"
                       "tau_XX at 0.625/0.3125 ns at core_radius_nm=12.5 (see 'Access-1.0 lifetime cap' section); "
                       "open red circles are one_pair_valid=False points.")),
        ("g2_flux_vs_ths.png", plot_vs_axis,
         (out, "g2_flux_vs_ths.png", core_only, "T_hs", ("g2_op", "collected_flux_pulsed_s"),
          "g2 and flux vs T_hs, both strain bounds (x_in=0.40, rep_rate_hz=200MHz, SET regime, reference geometry)", fixed_common),
         dict(caption="Charging wall: deterministic loading at 230-300K fails the Coulomb-blockade screen for "
                       "every core_radius_nm>=10 nm priced here; RC caveat applies to the as-built horizontal "
                       "contact (see 'RC caveat' section); open red circles are one_pair_valid=False points.")),
    ]
    for name, fn, args, kwargs in figs:
        result = _safe_plot(name, plot_fail_log, fn, *args, **kwargs)
        if result is not None:
            plots[name] = result

    plots["delivered_vs_commanded_flux.png"] = _safe_plot("delivered_vs_commanded_flux.png", plot_fail_log,
                                                            plot_delivered_vs_commanded, out, core) or {}
    plots["hardware_screens.png"] = _safe_plot("hardware_screens.png", plot_fail_log,
                                                plot_hardware_screens, out, core) or {}
    plots["strain_reversal_map.png"] = _safe_plot("strain_reversal_map.png", plot_fail_log,
                                                   plot_reversal_map, out, core) or {}

    text = _results_md(core, quick, complete, runtime_s, evaluate_calls, invalid_by_kind, dipole_rows, decisions)
    return plots, plot_fail_log, text


def _main_report_only(out, a):
    """Sweep fix 1 (Attempt 2, reporting-only round): regenerate
    results.md, the PNGs and manifest.json from the EXISTING sweep.csv in
    `out` -- NO evaluate() call, NO physics change, sweep.csv itself is
    NEVER rewritten (its sha256 before/after this mode is asserted
    identical below). manifest.json records report_only=True and the
    sweep.csv sha256 it consumed; output_hashes are recomputed for every
    file actually in the run dir; run.log is APPENDED to, never
    truncated, so the original run's own log survives underneath this
    mode's own lines."""
    if not out.is_dir():
        raise SystemExit("--report-only requires an existing --out-dir from a prior run")
    sweep_path = out / "sweep.csv"
    manifest_path = out / "manifest.json"
    if not sweep_path.is_file() or not manifest_path.is_file():
        raise SystemExit("--report-only requires an existing sweep.csv AND manifest.json in --out-dir")
    prior_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    sweep_hash_before = hashlib.sha256(sweep_path.read_bytes()).hexdigest()

    t0 = time.time()
    log_f = open(out / "run.log", "a", encoding="ascii", errors="replace", newline="\n")
    _orig_stdout, _orig_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(_orig_stdout, log_f), _Tee(_orig_stderr, log_f)
    try:
        print("--report-only: regenerating results.md/figures/manifest.json from the existing "
              "sweep.csv (sha256=%s) -- no evaluate() call this run" % sweep_hash_before, flush=True)
        core = _read_csv_rows(sweep_path)
        # attach_bound_reversal is a declared post-hoc transform, a pure
        # function of already-written sweep.csv columns (see its own
        # docstring) -- recomputing it here reproduces the SAME values
        # already on disk; it mutates the in-memory `core` list only,
        # never sweep.csv itself (sweep.csv is not rewritten in this mode).
        attach_bound_reversal(core)
    finally:
        sys.stdout, sys.stderr = _orig_stdout, _orig_stderr
        log_f.close()

    sweep_hash_after = hashlib.sha256(sweep_path.read_bytes()).hexdigest()
    if sweep_hash_after != sweep_hash_before:
        raise SystemExit("--report-only must never modify sweep.csv (hash changed unexpectedly)")

    invalid_by_kind = {}
    for r in core:
        if not r.get("valid"):
            invalid_by_kind[r.get("row_kind", "unknown")] = invalid_by_kind.get(r.get("row_kind", "unknown"), 0) + 1
    dipole_rows = [r for r in core if str(r.get("row_id", "")).startswith("DW")]

    quick = bool(prior_manifest.get("quick", False))
    complete = bool(prior_manifest.get("complete", False))
    runtime_s = float(prior_manifest.get("runtime_s", 0.0))
    evaluate_calls = int(prior_manifest.get("evaluate_calls", 0))
    decisions = list(prior_manifest.get("decisions", []))

    plots, plot_fail_log, text = _build_figures_and_results(
        out, core, quick, complete, runtime_s, evaluate_calls, invalid_by_kind, dipole_rows, decisions)
    (out / "results.md").write_text(text, encoding="utf-8")

    files = [p.name for p in out.iterdir() if p.is_file() and p.name != "manifest.json"]
    hashes = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}
    report_only_runtime_s = time.time() - t0
    manifest = dict(prior_manifest)
    manifest.update({
        "report_only": True,
        "report_only_source_sweep_sha256": sweep_hash_before,
        "report_only_runtime_s": report_only_runtime_s,
        "plot_row_mapping": plots,
        "plot_trace_failures": plot_fail_log,
        "output_hashes": hashes,
    })
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print("report_only=True sweep_csv_sha256=%s report_only_runtime_s=%.2f core_rows=%d invalid_total=%d" % (
        sweep_hash_before, report_only_runtime_s,
        len([r for r in core if r["row_kind"] == "core"]), sum(invalid_by_kind.values())))
    return 0


# -------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out-dir", default=str(ROOT / "out" / "nitride_nanowire"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--max-evaluations", type=int, default=10000)
    a = ap.parse_args(argv)

    if a.report_only:
        return _main_report_only(_safe(a.out_dir), a)

    counts = planned_counts(a.quick)
    horiz_core = len(build_core("horizontal_as_built", a.quick))
    vert_core = len(build_core("vertical_photonic", a.quick))
    if a.dry_run:
        expect_horiz, expect_vert = (1536, 1024) if not a.quick else (horiz_core, vert_core)
        out = dict(counts, horizontal_core_rows=horiz_core, vertical_core_rows=vert_core,
                   expected_horizontal_core_rows=expect_horiz, expected_vertical_core_rows=expect_vert,
                   max_evaluations=a.max_evaluations,
                   within_cap=counts["planned_evaluate_calls_upper_bound"] <= a.max_evaluations
                   and counts["planned_evaluate_calls_upper_bound"] <= 10000,
                   quick=a.quick)
        print(json.dumps(out, sort_keys=True))
        ok = out["within_cap"] and (a.quick or (horiz_core == 1536 and vert_core == 1024))
        return 0 if ok else 1

    if not a.quick and (horiz_core != 1536 or vert_core != 1024):
        raise SystemExit("full-mode core grid does not match the contract's 1536/1024 row counts")
    if counts["planned_evaluate_calls_upper_bound"] > a.max_evaluations or counts["planned_evaluate_calls_upper_bound"] > 10000:
        raise SystemExit("planned evaluate calls (%d) exceed cap" % counts["planned_evaluate_calls_upper_bound"])

    out = _safe(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str((ROOT / "out" / "nitride_nanowire" / "mplconfig"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    t0 = time.time()

    # H1 fix: tee stdout/stderr into run.log INSIDE the run directory for
    # the noisy part of the run (progress lines from _evaluate_all); closed
    # and restored BEFORE output_hashes is computed below so the recorded
    # hash matches the final, static file content. ASCII-only, explicit
    # newline (spec: "ASCII-only ... never print non-ASCII").
    log_path = out / "run.log"
    log_f = open(log_path, "w", encoding="ascii", errors="replace", newline="\n")
    _orig_stdout, _orig_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(_orig_stdout, log_f), _Tee(_orig_stderr, log_f)
    try:
        # Phase A: build every planned row's (kind, payload) job up front
        # (pure parameter-dict construction, no evaluate() calls) so the
        # dedup+dispatch phase below sees the whole run's identity set at
        # once.
        core_jobs = {f: build_core(f, a.quick) for f in FAMILIES}
        cut_jobs = build_reduced_cuts(a.quick)
        dipole_p = build_dipole_falsification()
        d2013_jobs = build_deshpande2013(a.quick)
        planar_jobs = build_planar_reference(a.quick)

        all_jobs = []
        for family in FAMILIES:
            all_jobs += [("nanowire", p) for p in core_jobs[family]]
        all_jobs += [("nanowire", p) for p in cut_jobs]
        all_jobs += [("nanowire", p) for p in dipole_p]
        all_jobs += [("nanowire", p) for p in d2013_jobs]
        all_jobs += [("planar", planar) for planar in planar_jobs]
        all_jobs.append(("deshpande2014", None))

        # Phase B: evaluate every UNIQUE job once, across a process pool
        # (see N_WORKERS/_evaluate_all above -- the per-call cost of this
        # device's injector resonance search makes serial evaluation of the
        # mandated 2560-row core grid alone exceed the operational time
        # budget). The physics-only job identity (_physics_only) means a
        # reduced-cut job whose sampled value equals the main-grid default
        # collapses onto the SAME unique identity as its core-row twin.
        results = _evaluate_all(all_jobs, n_workers=N_WORKERS)
        counter = {"evaluate_calls": len({_job_identity(j) for j in all_jobs})}
        print("evaluate_calls (unique, physics-only identity)=%d total_jobs=%d" % (
            counter["evaluate_calls"], len(all_jobs)), file=sys.stderr, flush=True)
    finally:
        sys.stdout, sys.stderr = _orig_stdout, _orig_stderr
        log_f.close()

    # Phase C: reassemble rows in the ORIGINAL enumeration order, reading
    # each row's scalars back out of `results` (seen/_row reproduces the
    # original in-process cache's cache_hit semantics; see _row's docstring).
    seen = set()
    rid = 0
    core = []

    def add_core(family):
        nonlocal rid
        for p in core_jobs[family]:
            rid += 1
            core.append(_row(f"CO{rid:05d}", "core", p, seen, results))

    for family in FAMILIES:
        add_core(family)

    for p in cut_jobs:
        rid += 1
        core.append(_row(f"SN{rid:05d}", "sensitivity", p, seen, results))

    dipole_rows = []
    for p in dipole_p:
        rid += 1
        r = _row(f"DW{rid:05d}", "sensitivity", p, seen, results)
        core.append(r)
        dipole_rows.append(r)

    for p in d2013_jobs:
        rid += 1
        core.append(_row(f"D13{rid:05d}", "sensitivity", p, seen, results))

    for planar in planar_jobs:
        rid += 1
        core.append(_row(f"PR{rid:05d}", "planar_reference", {"family": "planar", "regime": planar[1]},
                          seen, results, planar=planar))

    rid += 1
    core.append(_deshpande2014_replay(f"LI{rid:05d}", seen, results))

    attach_bound_reversal(core)

    _write_csv(out / "sweep.csv", core)

    runtime_s = time.time() - t0
    complete = bool((not a.quick) and counter["evaluate_calls"] <= a.max_evaluations and runtime_s < 1800.0
                     and len([r for r in core if r["row_kind"] == "core"]) == 2560)

    invalid_by_kind = {}
    for r in core:
        if not r.get("valid"):
            invalid_by_kind[r.get("row_kind", "unknown")] = invalid_by_kind.get(r.get("row_kind", "unknown"), 0) + 1

    decisions = [
        "Measured evaluate() cost on this machine: 6.1-9.8 s/call serially (fsim_core.nitride_nanowire_"
        "injector.injector_feasibility's resonance search dominates, ~85 percent of wall time, run on "
        "EVERY row regardless of regime -- a dependency-module property, not something this piece may "
        "patch). At that cost the mandated 2560-row core grid cannot finish serially inside 1800 s "
        "(a 436-call serial quick run took 1587 s). Rejected alternative: leave evaluation serial and "
        "let the full run overrun the budget. Chosen: dispatch every unique evaluate() call across a "
        "process pool (N_WORKERS=os.cpu_count(), see _evaluate_all) -- identical inputs/outputs per row, "
        "only wall-clock parallelized; benchmarked throughput ~1.0-2.0 calls/s depending on core radius.",
        "Reduced-cut axis LIST trimmed to the orchestrator task's explicit reduced-cut deliverable list "
        "(screening {0,1}; occupied_dot_access's 1.0 partner; S {1e2,1e4}; shell AlGaN; the spec-named "
        "injector cuts -- barrier thickness/Al fraction/alignment/growth tolerance/occupation-control "
        "uncertainty; R_s designed 1e6 on the horizontal family) plus one b_res=0.02 point, each sampling "
        "only its alternative (non-default) value(s), PLUS this fix round's restored current_pulse_width "
        "cut (I in {0.001,0.002,0.02} uA x tau_pulse in {0.01,0.1,1} ns, 9 rows, horizontal reference "
        "geometry, rectangular regime, relaxed bound, 300 K, 200 MHz) -- see DROPPED_CUT_AXES in scripts/"
        "run_nitride_nanowire.py and results.md's 'Reduced-cut coverage' section for the axes still "
        "dropped relative to the contract's fuller reduced-cut list (reservoir_access, gamma300, "
        "tau_rad0_ns, tau_cap_ps, C_parasitic_F, Rth_K_W, vertical R_s_ohm, NA/bottom_reflectivity). "
        "Rejected alternative: the contract's full one-at-a-time axis list (~1350 extra calls), which "
        "single-benchmark projections put close to or beyond the 40-minute stop threshold once added to "
        "the core grid's own runtime; the core/main grid itself is NEVER reduced. The job/cache identity "
        "(_physics_only) was also changed this round to exclude the bookkeeping sensitivity_axis/"
        "sensitivity_value labels, so a reduced-cut row whose sampled value equals the main-grid default "
        "(e.g. the screening_fraction=0.0 cut) now correctly deduplicates onto its core-row twin's real "
        "evaluate() call instead of hashing separately -- rejected alternative: leave the identity as-is "
        "and simply accept the ~32 redundant real evaluate() calls it caused.",
        "bound_reversal is computed as a declared post-hoc transform over paired core rows "
        "(bound_reversal_pair column) rather than via a second evaluate_strain_pair() call, to avoid "
        "doubling the core-grid evaluate() count; the device module's own per-row bound_reversal field "
        "stays its 'not_computed' sentinel on every row here.",
        "The 2013 CW comparison is published as a pulsed-model drive_mismatch (this evaluator has no "
        "CW/HBT path), per the contract's own instruction, never substituted for a predicted CW g2.",
    ]

    plots, plot_fail_log, text = _build_figures_and_results(
        out, core, a.quick, complete, runtime_s, counter["evaluate_calls"], invalid_by_kind, dipole_rows, decisions)
    (out / "results.md").write_text(text, encoding="utf-8")

    files = [p.name for p in out.iterdir() if p.is_file() and p.name != "manifest.json"]
    hashes = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}
    manifest = {
        "quick": a.quick, "complete": complete, "runtime_s": runtime_s,
        "evaluate_calls": counter["evaluate_calls"], "max_evaluations": a.max_evaluations,
        "n_workers": N_WORKERS,
        "dropped_reduced_cut_axes": list(DROPPED_CUT_AXES),
        "row_counts": {"core": len([r for r in core if r["row_kind"] == "core"]),
                        "sensitivity": len([r for r in core if r["row_kind"] == "sensitivity"]),
                        "planar_reference": len([r for r in core if r["row_kind"] == "planar_reference"]),
                        "planar_2014_replay": len([r for r in core if r["row_kind"] == "planar_2014_replay"]),
                        "total": len(core)},
        "invalid_counts_by_kind": invalid_by_kind,
        "invalid_total": sum(invalid_by_kind.values()),
        "grid_specification": {
            "families": list(FAMILIES), "regimes": list(REGIMES), "strain_bounds": list(STRAIN_BOUNDS),
            "rep_rates_hz": list(REP_RATES), "core_radius_nm": {k: list(v) for k, v in CORE_R_NM.items()},
            "height_nm": list(HEIGHT_NM), "x_in": list(X_IN), "T_hs_K": list(T_HS),
            "quick_core_radius_nm": {k: list(v) for k, v in QUICK_CORE_R_NM.items()},
            "quick_height_nm": list(QUICK_HEIGHT_NM), "quick_x_in": list(QUICK_X_IN), "quick_T_hs_K": list(QUICK_T_HS),
        },
        "declared_transforms": {
            "bound_reversal_pair": "attach_bound_reversal(): pairs each relaxed core row with its unrelaxed "
                                     "partner at identical other coordinates and reproduces nitride_nanowire_"
                                     "device._bound_reversal's own mu/collected_flux_pulsed_s/g2_op ordering "
                                     "rule (both partners must clear the 1000/s optical floor); never a new "
                                     "evaluate() call.",
        },
        "plot_row_mapping": plots,
        "plot_trace_failures": plot_fail_log,
        "output_hashes": hashes,
        "card_hashes": {n: _card_hash(n) for n in set(list(CARD_NAME.values()) + list(PLANAR_CARD_NAME.values()) + [DESHPANDE2014_CARD])},
        "source_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in (ROOT / "fsim_core").glob("nitride_nanowire_*.py")},
        "decisions": decisions,
        "versions": {"python": sys.version.split()[0]},
        "resolved_out_dir": str(out),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print("evaluate_calls=%d runtime_s=%.2f complete=%s core_rows=%d invalid_total=%d" % (
        counter["evaluate_calls"], runtime_s, complete,
        len([r for r in core if r["row_kind"] == "core"]), sum(invalid_by_kind.values())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
