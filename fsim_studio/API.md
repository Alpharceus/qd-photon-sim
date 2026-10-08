# FSIM Studio API

Base URL `http://127.0.0.1:<port>` (default 8765; the server binds 127.0.0.1 only). Every `/api/*`
response is JSON except the SSE stream. Launch with `python -m fsim_studio [--browser] [--no-window]
[--port N]`. Checked by `verify/verify_studio_api.py`.

## Conventions

- **Request guard (H13).** The `Host` header must be `127.0.0.1:<port>` or `localhost:<port>` (the
  port `python -m fsim_studio` was started with; the bare names without a port only for
  `create_app()` without a port, i.e. the Flask test client); anything else is **403**. Every
  `POST`, `PUT`, `DELETE` and `PATCH` must carry the header **`X-FSIM-Studio: 1`**, else **403**
  (`fsim_studio/web/js/api.js` sends it on every non-GET request). `T_grid` carries at most 400 points.

- **NaN and inf are `null`.** When a dict has a NaN under key `k` and a non-empty `k_invalid_reason`
  (the fsim_core convention, e.g. `g2_op` / `g2_op_invalid_reason`), the dict also carries
  `k__nan_reason`. Render "n/a + reason", never 0.
- **Numbers are never rounded by the server.** Values from `out/` files are the exact printed text
  (VERDICT fields stay strings, e.g. `"16/768"`, `"0.9817"`); live values are full floats.
- **Errors**: `{"error": "<message>"}` with 400 (bad input), 404 (unknown card/campaign/job/story),
  409 (refused overwrite), 422 (card/design fails `DeviceDesign.load` validation), 501 (scene module
  missing).
- **File paths** returned by campaign endpoints (`files.csv`, `files.figures`, `primary_csv`, `md`)
  are relative to `out/`; fetch them as `/out/<path>`.
- **Static**: `/` serves `fsim_studio/web/index.html`; any other non-`/api` path is served from
  `fsim_studio/web/` (so `/vendor/...`, `/css/...`, `/js/...`).
- **Cache**: results live in `.cache/studio/<sha256>.json` (gitignored; `FSIM_STUDIO_CACHE` env var
  overrides the folder). The key covers the loaded design, mode, ranged, T grid, the sha256 of all
  `fsim_core/*.py` sources (recomputed whenever a source's mtime/size changes, so an edit while the
  server runs takes effect) and the Studio result-schema version (`cache.STUDIO_SCHEMA`), so a physics
  or result-shaping change invalidates every cached result.

## Health and theme

### GET /api/health
```json
{"ok": true, "version": "0.1.0", "fsim_core_hash": "6aa4cb0e95d4bf63c19ac859bb1f2e92f6e705cc50dffbab76b04ccb347bfc4c"}
```

### GET /api/theme/plotly?mode=light|dark
`fsim_theme.plotly_template(mode).to_plotly_json()`: `{"data": {...}, "layout": {...}}`. Use it as
`layout.template`.

## Cards

### GET /api/cards
All `cards/*.yaml`. `kind` is `design` (has a `design:` block, 14 cards) or `param`; `platform` is
`legacy | ingan_gan_planar | ingan_gan_nanowire` for design cards, `null` for parameter cards.
```json
[{"name": "edge-inp-gainp-design", "kind": "design", "platform": "legacy",
  "device": "InP self-assembled QD in Ga0.51In0.49P ...", "title": "edge-inp-gainp-design",
  "role": "device-design"},
 {"name": "chatzarakis", "kind": "param", "platform": null, "device": "(211)B InAs/GaAs ...",
  "title": "chatzarakis", "role": "validation (V-a)"}]
```

### GET /api/cards/<name>
The loaded `DeviceDesign` as a dict (`dataclasses.asdict`, exactly what `PUT`/`/api/run` accept back),
plus `design_meta.META` for every dotted path present (`band` = `[lo, hi]` or null; `choices` when
META has them). 422 for a parameter card.
```json
{"name": "staged-device-design", "platform": "legacy",
 "design": {"name": "staged-device", "dot": {"delta_xx": 3.5, "gamma_scale": 1.0, "...": "..."},
            "ret": {...}, "drive": {...}, "thermal": {"T_hs": 77.0, "layers": [...], "substrate": {...}},
            "cavity": {...}, "filter": {...}, "aperture": {...}, "emission": {...},
            "platform": "legacy", "nitride": {}, "provenance": {}, "quantum": {}},
 "meta": {"dot.delta_xx": {"tag": "A", "unit": "meV", "band": [2.0, 14.0],
                           "source": "(211)B InAs/GaAs 4-13 meV distribution; platform-dependent",
                           "label": "delta_xx"}}}
```

### PUT /api/cards/<name>   body `{"design": {...}, "presets"?: {dot, template, ...}}`
"Save as card" (studio-p2a). `<name>` must be a lower-case slug (`[a-z0-9]+(-[a-z0-9]+)*`, at most 80
characters; `unsaved` and `studio` are reserved), else **400**. A shipped card name (any
`cards/<name>.yaml`) is refused with **409**; a new name is allowed, and a design saved here earlier may
be saved again (overwritten). The file goes to **`cards/studio/<name>.yaml`** (kept apart from the
shipped `cards/*.yaml` that scripts, gates and fsim_gui glob; `FSIM_STUDIO_USER_CARDS` overrides the
folder, the verify scripts use a temp dir): the design is loaded through `DeviceDesign.load`, written
with `DeviceDesign.save`, read back (it must reproduce the design, else 422), then the Studio marker
`meta.saved_by: fsim-studio` (and `meta.presets` when given) is added. A file in that folder without
the marker is never overwritten (409). Saved cards are listed by `GET /api/cards` at once (with
`"saved": true` and `presets`), and `GET /api/cards/<name>`, `/verdict`, `/api/run` and the scene
endpoints resolve them like shipped cards.
```json
{"name": "preset-piezo-variant-micropillar", "path": "cards/studio/preset-piezo-variant-micropillar.yaml", "saved": true}
```

### GET /api/cards/<name>/verdict   (interface I6)
The committed VERDICT line(s) that speak for a design card, matched exactly on what the card carries:
regime (`drive.cycle_loading`), family/orientation (nanowire `family`; planar `dot.orientation`, a
nonpolar/a_plane card -> `family=a_plane` in `nitride_geometry_stark`; c-plane cards ->
`nitride_cavity`), strain bound and rep rate (nanowire), the headline `model=` token (RT edge ->
`rt_edge`). It never picks among several matches: `match` is `"exact"` (one line), `"ambiguous"` (all N
lines returned; render "ambiguous: N lines match") or `"none"` (`reason` says why: legacy non-edge card,
QW-fluctuation supplement, non-gating comparison card). 404 for an unknown card.
```json
{"campaign": "nitride_geometry_stark", "card": "nitride-nonpolar-set-design",
 "md": "nitride_geometry_stark/results.md", "criteria": {"family": "a_plane", "regime": "deterministic_pair"},
 "lines": [{"raw": "VERDICT: idealized_status=pass_hardware_infeasible family=a_plane ...", "line": 11, "...": "..."},
           {"raw": "VERDICT: ...", "line": 12, "...": "..."}],
 "match": "ambiguous", "n_matches": 2,
 "reason": "a_plane orientation: nitride_geometry_stark family a_plane, regime match; ambiguous: 2 lines match (not picked); they differ in ['screening']"}
```

### GET /api/gates   (interface I3)
The two acceptance gates, parsed from `out/rt_edge/verdict.md` (the g2 ceiling cross-checked against
`docs/rt_edge_contract.md`; a disagreement withholds the value). Use only this for the 0.5 line, the
gate box and the flux floor.
```json
{"g2_ceiling": {"value": 0.5, "tag": "A", "unit": "", "definition": "pulsed intrinsic g2(0) < 0.5",
                "source": "out/rt_edge/verdict.md:22", "contract_source": "docs/rt_edge_contract.md:29",
                "tag_source": "declared: acceptance-gate convention (no bracket tag printed next to it)"},
 "flux_floor": {"value": 1000.0, "unit": "photons/s", "tag": "A", "source": "out/rt_edge/verdict.md:90",
                "tag_source": "bracket tag printed in out/rt_edge/verdict.md",
                "definition": "pulsed collected-flux eligibility floor"}}
```

### GET /api/meta
```json
{"META": {"dot.delta_xx": {"unit": "meV", "tag": "A", "lo": 2.0, "hi": 14.0, "source": "..."}},
 "ENV_DEFAULTS": {"dot.delta_xx": [1.5, 5.0], "drive.mu": [0.1, 1.0]},
 "presets": {"dot": {...}, "template": {...}, "cavity": {...}, "drive": {...}, "injection": {...}}}
```

### POST /api/presets/apply   body `{dot, template, cavity, drive, injection?, name?}`
`presets.preset_device(dot, template, cavity, drive)` then `apply_injection_preset` (same order as the
DPG designer). Keys are the preset dict keys from `/api/meta`. `meta` carries the META annotation of
every dotted path in the design (the same shape as `GET /api/cards/<name>`), so an unsaved preset design
renders its tags before it is saved. Picking a preset never upgrades a field's tag.
```json
{"design": {"name": "preset-device", "dot": {"delta_xx": 8.0, "...": "..."}, "...": "..."},
 "meta": {"dot.delta_xx": {"tag": "A", "unit": "meV", "band": [2.0, 14.0], "source": "...", "label": "delta_xx"}}}
```

### POST /api/validate   body `{"design": {...}}`
`violations` = numeric META fields outside `[lo, hi]` (warnings, never clamps). `default_ranged` =
`design_meta.default_ranged(d)` (the envelope a run uses when `ranged` is omitted). `f8_floor` =
`1 - loading.f8b_thin_fano(eta_capture, F_p)` (mu must be >= this).
```json
{"violations": [{"path": "dot.delta_xx", "value": 20.0, "band": [2.0, 14.0]}],
 "default_ranged": {"dot.delta_xx": [1.5, 5.0], "dot.gamma_scale": [0.7, 1.3],
                    "drive.b_e": [0.001, 0.3], "drive.mu": [0.1, 1.0]},
 "f8_floor": 0.0}
```

## Runs and jobs

### POST /api/run
Body: `{card?: name, design?: {...}, mode: "point" | "headline" | "envelope", ranged?: {path: [lo, hi] |
[v1, v2, v3 ...]}, T_grid?: [T, ...] | {lo, hi, n}}`.
- `design` (if given) is the edited design; `card` is then only the label for `run_id`. With only
  `card`, the shipped card is loaded.
- `point`: `evaluate(d, T_grid)`. `headline`: the RT-edge switch selection (copy, `drive.finite_pulse
  = true`, `ret.tau_cap_scales_with_density = false`, `drive.cw = false`) then `evaluate`.
  `envelope`: `evaluate_envelope(d, ranged or default_ranged(d), T_grid)` over the process pool.
- `T_grid` omitted = fsim_core's default 120 points, 4-350 K.
- Legacy non-edge `point` runs with <= 200 T points run inline and return `result` directly. Everything
  else returns a `job_id`; follow it with SSE or polling. A cache hit returns `cached: true` and the
  result immediately (a finished job is still created, so `job_id` is always valid).

Response:
```json
{"job_id": "d1aae84509ad", "cached": false, "eta_s": 0.001, "run_id": "FS-staged-device-9862e2",
 "result": {"...": "only when inline or cached"}}
```
`run_id` = `"FS-" + card name without "-design" + "-" + first 6 hex of the cache key`.

**Point / headline result**
```json
{"curves": {"T_hs": [77.0, 150.0], "Tj": [78.26, 151.79], "g2": [0.0458, 0.3387],
            "eps": [...], "rho2": [...], "gamma": [...]},
 "scalars": {"T_j_op": 78.26, "g2_op": 0.04582197110623831, "brightness_per_pulse": 0.1967,
             "brightness_convention": "pre-retention static (P1+P2)*t_x; ...", "T_c": 175.5,
             "F_eff": null, "tag_chain": "[A]", "provenance": {...}, "...": "~90 keys, schema-driven"},
 "tag_chain": "A", "tag_chain_source": "fsim_core tag_chain",
 "provenance": {"linewidth": {"tag": "A", "note": "arsenide class proxy"}, "...": "..."},
 "labels": ["pre-retention"],
 "run_id": "FS-staged-device-9862e2", "mode": "point", "card": "staged-device-design",
 "platform": "legacy", "T_grid": [77.0, 150.0]}
```
`tag_chain` (interface I1) is never null: one of `"V" | "DR" | "E" | "A"`, the widest tag of the run.
`tag_chain_source` says where it came from: `"fsim_core tag_chain"` (legacy/edge evaluate returns one;
`scalars.tag_chain` keeps fsim_core's raw text, e.g. `"[A]"`), `"parsed from provenance"` (nitride /
nanowire: the widest bracket tag in `scalars.provenance`, e.g. `[A/E/DR]` -> `A`), or `"no tag
returned; assumed widest"` (then `A`).

`labels` (render them, they are mandatory):
- `"static (non-headline)"`: RT edge design with `drive.finite_pulse` false, any mode except
  `headline`. Put it in the g2 and brightness tile titles. For the headline number beside it, issue a
  second run with `mode: "headline", T_grid: [<thermal.T_hs>]`.
- `"headline model (drive.finite_pulse=true, ret.tau_cap_scales_with_density=false)"`: headline runs.
- `"pre-retention"`: `scalars.brightness_convention` starts with "pre-retention"; show the
  convention text verbatim under the brightness number.
- `"flux: commanded (collected_flux_pulsed_s) vs delivered (collected_flux_delivered_s)"`: nanowire
  results; show both numbers, never a bare "flux".
- `"all ranges at midpoint (not a prediction)"`: envelope results (label for the mid curve).
- `"non-headline model"` (interface I2): any RT-edge (`emission.type: edge`) run, outside `mode:
  "headline"`, whose model switches are not the headline set (`drive.finite_pulse` true AND
  `ret.tau_cap_scales_with_density` false; `drive.cw` only adds CW diagnostics and is not part of the
  model). A static edge card carries both this and `"static (non-headline)"`; a re-run of a
  finite_pulse + tau_cap_density row carries only this one.
- `"structural floor set by b_res [A]"` (interface I2): `scalars.g2_op` equals the background floor
  `1 - (1/(1 + drive.b_res))^2` within 1e-9 (`fsim_core.scene_support.is_background_floor`), e.g. the
  nitride SET cards' 0.17355 at b_res 0.1. Put the floor marker on those g2 tiles.

**Envelope result** (adds to the above; `curves` is the all-midpoint curve set, `scalars` come from
an `evaluate()` of the midpoint design, and `scalars_source` says so):
```json
{"curves": {"T_hs": [...], "g2": [...], "eps": [...], "rho2": [...], "Tj": [...], "gamma": [...]},
 "bands": {"g2": {"lo": [...], "hi": [...], "mid": [...]}, "eps": {...}, "rho2": {...},
           "Tj": {...}, "gamma": {...}},
 "scalar_bands": {"T_c": [151.2, 196.0], "g2_op": [0.0213, 0.0939], "eps_op": [...], "rho_op": [...],
                  "dT_J": [...]},
 "tornado": {"dot.delta_xx": 12.3, "dot.gamma_scale": 30.1, "drive.b_e": 0.4, "drive.mu": 2.2},
 "n_samples": 16,
 "ranged": {"dot.delta_xx": [1.5, 5.0], "...": "..."},
 "scalars_source": "mid-design point evaluate(): all ranges at midpoint (not a prediction)",
 "labels": ["pre-retention", "all ranges at midpoint (not a prediction)"], "mode": "envelope"}
```
`scalar_bands` values may be `[null, null]` (no finite sample, e.g. T_c not crossed). A tornado value
of `null` means the collapsed band is not finite: render "band not finite", not 0.

**Progressive results**: for slow platforms (edge emitters, finite-pulse designs, nanowires) with more
than one T point, a `partial` event first carries the single-T (`thermal.T_hs`) result with
`"partial": true`, then `done` carries the full curve.

**Errors** (`state: "error"`, `error` object):
```json
{"kind": "f8_domain",
 "message": "evaluate_envelope: every sample in the ranged box violates the F8 loading domain ...",
 "f8_floor": 0.6,
 "suggestion": "Every sample sits outside the F8 loading domain (mu >= 1 - F_eff). Narrow drive.mu to >= 0.6, or raise drive.F_p / drive.eta_capture."}
```
Other kinds: `value_error`, `error`, `pool` (worker pool broke; retry).

### GET /api/jobs/<id>/events   (Server-Sent Events)
`text/event-stream`; replays all past events, then streams until a terminal event. Comment lines
`: keepalive` every 15 s. Events (`data:` is JSON):
```
event: state     data: {"state": "running", "job_id": "..."}          (queued | running | done | cancelled)
event: progress  data: {"job_id": "...", "k": 17, "n": 49}             (envelope; n = box + midpoint + tornado collapses)
event: partial   data: {"job_id": "...", "result": {...}}              (slow point/headline runs)
event: done      data: {"job_id": "...", "result": {...}}              (terminal)
event: error     data: {"job_id": "...", "error": {...}}               (terminal)
```
A cancel ends the stream with `state: cancelled`. Cached results produce `state: done` + `done` with
`"cached": true`. EventSource reconnects on drop; the replay makes that safe.

### GET /api/jobs/<id>
```json
{"job_id": "cd428f323664", "kind": "envelope", "state": "running",
 "progress": {"k": 17, "n": 49}, "eta_s": 0.03, "run_id": "FS-staged-device-466466"}
```
Adds `partial` (while running, once available), `result` (done) or `error` (error).

### POST /api/animate   (response animation, studio-p2a)
Body: `{card?: name, design?: {...}, param: "<numeric META path>" (default "thermal.T_hs"), lo?, hi?
(default the META band), n?: 2..120 (default 48; 16 on the nanowire tier), mode: "point" | "headline",
grid?: "op" | "full"}`. Computes one frame per value of `numpy.linspace(lo, hi, n)`: the design with
`param` set to that value, run exactly as `POST /api/run` would (same headline switch selection, same
cache key `{design, mode, ranged: null, T_grid}`, same `run_id`), one process-pool task per frame. A
frame is therefore a cache hit for `/api/run` and vice versa, and a repeated animation is served from
the cache at once (`cached: true`, `result` inline). `grid` (default `op` for edge / nanowire classes,
`full` otherwise): `op` evaluates `T_grid = [thermal.T_hs of the frame design]` (T_c is then null),
`full` the grid `/api/run` uses with no grid (fsim_core's default; the nitride display grid). `param`
must be a numeric META field with a `[lo, hi]` band that the design carries; `lo == hi`, `lo`/`hi`
outside that META validity band (T_hs 4..350 K, ...), `n` that is not an integer in 2..120 (`0`, `2.5`,
`"7"`, `true` included: no silent default), `mode: "envelope"` are **400**.

Response: `{job_id, cached, eta_s, n, values, n_cached, result?}`. Follow the job with SSE: one
`frame` event per frame as it completes (out of order; cached frames first), interleaved with
`progress` `{k, n}`, then `done`:
```
event: frame     data: {"job_id": "...", "index": 7, "cached": false, "frame": {...}}
```
Frame (also the entries of `result.frames`; a failed frame is `{index, value, error: {kind, message}}`):
Honesty fields (studio-p2d, physics brief section 7): every frame also carries `grid` (`op` | `full`);
a `grid: "op"` frame always has `T_c: null` with `T_c__nan_reason` (T_c is never read from an op frame), and
`one_pair_valid`, `set_feasible`, `flat_band`, `depletion_regime` ride along verbatim when the result has
them. The final `result` adds `runaway_at` (index of the first frame whose evaluate() reports `runaway`, else
null; that frame and every later one carry `after_runaway: true`), `density` (rule 6, advisory: for
`thermal.T_hs` `{rule: "arrhenius", checked, T_c, window, max_step, required_step: 5, ok, suggest: {lo, hi, n}}`,
for `drive.V` on an InGaN/GaN platform `{rule: "stark", max_step, required_step: 0.2, ok}`) and `slice_note`
("one range at a time, not a prediction" when the animated parameter's tag is `A`, else null). Default `n` for
`drive.V` on an InGaN/GaN platform is one frame per 0.2 V. The browser blends only two adjacent frames that
are `g2_op_valid`, not runaway, with equal label sets and no regime change; everything else snaps to the nearest
computed frame.

```json
{"index": 2, "value": 100.0, "run_id": "FS-staged-device-f97273", "tag_chain": "A",
 "tag_chain_source": "fsim_core tag_chain", "labels": ["pre-retention"], "provenance": {"...": "..."},
 "scalars": {"g2_op": 0.0951512514110624, "collected_flux_pulsed_s": null,
             "brightness_per_pulse": 0.1967346701436833, "T_j_op": 101.44588313687613,
             "eps_op": 0.0540624331889245, "rho_op": 0.9818900186269597, "T_c": 175.5180852366353,
             "brightness_convention": "pre-retention static (P1+P2)*t_x; ...", "g2_op_valid": true,
             "g2_op_invalid_reason": "", "runaway": false, "flux_measurable": false, "rep_rate_hz": null},
 "curves": {"T_hs": [...], "g2": [...]}}
```
`scalars` is the subset `g2_op, collected_flux_pulsed_s, collected_flux_delivered_s,
brightness_per_pulse, T_j_op, eps_op, rho_op, T_c` (each only when the run returns it, with its
`__nan_reason`) plus the verbatim flags `brightness_convention, g2_op_valid, g2_op_invalid_reason,
runaway, flux_measurable, rep_rate_hz`: copied from the run, never interpolated. `curves` appears only
for multi-point grids. Result: `{param, lo, hi, n, values, mode, grid, card, platform, unit, param_tag,
frames: [...], computed: [indices], errors: [indices]}`. Interpolation between computed frames is the
browser's display business and is always marked "interp." there.

### DELETE /api/jobs/<id>
Marks the job cancelled and returns its snapshot (`state: "cancelled"`). Queued futures are cancelled;
an evaluation already running in a worker finishes and its result is discarded (the pool is not
killed).

### POST /api/compare/report   body `{"entries": [{"label": "A", "run_id": "FS-..."}], "title"?: str}`
Loads each run from the cache and calls `fsim_viz.report.designer_report` unchanged (labels
`current`/`A`/`B`/`C` get the fixed slot colours). The bundle is written to
`.cache/studio/reports/<timestamp>-<hash>/` (never `out/`). 404 if a run_id is not cached.
```json
{"path": "C:\\...\\.cache\\studio\\reports\\20261007-110504-3fa2c1",
 "files": ["curves_A.csv", "curves_B.csv", "design_A.yaml", "design_B.yaml", "ranged_B.csv",
           "report.pdf", "report.png", "report.svg", "scalars_comparison.csv"]}
```

## Campaigns (out/)

### GET /api/campaigns
Campaigns = folders with `manifest.json` at `out/*/` or `out/*/*/`, excluding scratch / `quick` /
`c6_scratch_lw_quick` / presentation folders: `rt_edge`, `nitride_cavity`, `nitride_geometry_stark`,
`nitride_nanowire/full`. Then collections without a manifest: `rt_campaign`, `tier_geometry`,
`tier_device` (`kind: "collection"`, "figures + CSV") and `phase0-3`, `spec`, `validation`, `zhao`
(`kind: "legacy"`, "legacy F-series").
```json
[{"id": "rt_edge", "title": "rt edge", "kind": "campaign", "path": "out/rt_edge", "has_manifest": true,
  "verdicts": [{"raw": "VERDICT: FAIL model=finite_pulse:true,tau_cap_density:false g2_min=0.9817 ...",
                "kind": "VERDICT", "word": "FAIL", "line": 6,
                "fields": {"model": "finite_pulse:true,tau_cap_density:false", "g2_min": "0.9817",
                           "flux_max": "1117", "eligible": "16/768", "...": "..."}}],
  "best": [], "status_word": "FAIL", "status_words": {"FAIL": 1},
  "generated": "2026-09-23T19:40:37.059843+00:00", "commit": "5beb620",
  "commit_source": "git log (last commit touching the folder)", "stale": false, "stale_files": null},
 {"id": "tier_geometry", "kind": "collection", "collection_label": "figures + CSV", "stale": true,
  "stale_files": ["corner_frontier.csv"], "stale_note": "pre-audit (2026-08-12), may carry H5 1-D Purcell",
  "verdicts": [], "status_word": null, "...": "..."}]
```
- `word`: the bare token after `VERDICT:` (`FAIL`), else the `idealized_status` field
  (`no_idealized_pass`, `pass_hardware_infeasible`). `status_word` is the single word when all lines
  agree, else `"mixed"`; `status_words` counts each word. Never recompute a verdict in the UI.
- Prose tolerance: a bare token after `key=value` continues that value (`"model": "finite electrical
  pulse"`); a `( ... )` run goes to `note` (BEST lines, e.g. "one delivered photon per 707.3 cycles").
- `commit`: the manifest's generation commit when present (`commit_source: "manifest"`), otherwise the
  last git commit touching the folder. `stale`: pre-audit badge (`tier_geometry/corner_frontier.csv`,
  every `tier_device` file); `stale_files: null` with `stale: true` means the whole folder.
- BEST lines (`BEST_PASSING_FLUX`, `BEST_PASSING_FLUX_300K`) parse the same way, `kind` = the token.

### GET /api/campaigns/<id>      (id may contain "/", e.g. `nitride_nanowire/full`)
```json
{"id": "rt_edge", "manifest": {"schema_version": 1, "grid": {...}, "headline_model": {...}, "...": "..."},
 "verdicts": [...], "best": [...], "status_word": "FAIL", "status_words": {"FAIL": 1},
 "commit": "5beb620", "commit_source": "...", "generated": "...", "md": "rt_edge/verdict.md",
 "files": {"csv": ["rt_edge/sweep.csv"], "figures": ["rt_edge/envelope.png", "rt_edge/gui-smoke.png"],
           "md": ["rt_edge/_verify_fixture_verdict.md", "rt_edge/verdict.md"],
           "other": ["rt_edge/evidence.json", "rt_edge/manifest.json"]},
 "primary_csv": "rt_edge/sweep.csv",
 "columns": [{"name": "card_id", "role": "input", "kind": "string"},
             {"name": "g2_pulsed", "role": "output", "kind": "number"},
             {"name": "eligible_pulsed", "role": "verdict", "kind": "bool"}],
 "stale": false, "stale_files": null, "stale_note": null,
 "headline_filter": {"model_finite_pulse": true, "model_tau_cap_density": false}}
```
`columns` describe `primary_csv` (`sweep.csv`, else the first CSV). `role` is a heuristic
(input = manifest grid/axes names, `*_tag`, known axis names; verdict = bool or gate-like names;
output = the rest); `kind` = `number | bool | string | mixed | empty`. `headline_filter` is non-null
only for `rt_edge`.

### GET /api/campaigns/<id>/sweep?file=&filter=&cols=&limit=5000&offset=0&headline=0|1
Typed in-memory table per CSV (parsed once, reused): `''`, `nan`, `NaN`, `None`, `not_applicable` ->
`null`; `True`/`False` -> booleans; numbers -> floats; else strings.
- `file`: a CSV name inside the campaign folder (default `primary_csv`).
- `filter`: URL-encoded JSON `{col: value | [values] | {"min": a, "max": b}}`, ANDed. Values match by
  type (`true` matches the CSV `True`; `"True"` is also accepted).
- `headline=1` (rt_edge only): pins `model_finite_pulse=True, model_tau_cap_density=False` (768 of
  3072 rows). Default it ON for rt_edge; other model combos belong in a greyed "never gates" view.
- `cols`: comma list; rows are arrays in that column order.
```json
{"file": "rt_edge/sweep.csv", "columns": ["T_hs_K", "g2_pulsed"],
 "rows": [[230.0, 0.9983614825044536], [230.0, 0.9983614825044536]],
 "total": 3072, "filtered": 768, "returned": 2, "offset": 0, "headline": true,
 "filter": {"model_finite_pulse": true, "model_tau_cap_density": false}}
```
The nanowire sweep (3216 x 265) answers in ~0.4 s after its first load.

### GET /out/<path>
Static file from `out/` (figures, CSV, md, json).

## Scenes (3D; builders owned by fsim_studio/scene.py)

`GET /api/scene/device?card=&T=`, `/api/scene/band?card=&T=`, `/api/scene/cascade?card=&T=`,
`/api/scene/surface?card=&param=w&mode=envelope`, `/api/scene/lattice?campaign=&x=&y=&z=&headline=1`
-> SceneSpec JSON (schema in scene.py's docstring). 501 while `fsim_studio/scene.py` is missing; 400
for builder errors (`{"error": "..."}`). `card` is a card NAME, resolved only through
`api_cards.load_design` (H10): a path or anything that is not a valid card name is 400, an unknown card
404, a parameter card 422.

Lattice (interface I4): `overlays.rows` carries per row `model_finite_pulse`, `model_tau_cap_density`
(null for campaigns without model axes) and `gates` (true only for headline-model rows on rt_edge,
true for every row elsewhere). `overlays.banner` is an object built from GATING rows only:
`{text, style: "pass"|"fail", pass, gating_rows, rows, pass_column, non_gating_passes: null |
{count, rows, text}}`; with `headline=0` on rt_edge it reads "0 / 768 headline-model rows pass
headline_pass" and `non_gating_passes` = 232 / 2304 with the text "... never gates (shown greyed)".
Draw non-gating rows greyed; never style the banner from them.

The only job descriptor function the backend runs is `fsim_studio.scene:surface_compute` (exact
allowlist, H16); anything else is 400.

When a builder returns a job descriptor (`{"job": true, "fn": "fsim_studio.scene:...", "kwargs",
"eta_s"}`, e.g. surface on edge / planar-nitride cards), the backend runs it in the pool and answers
**202**:
```json
{"job_id": "969bc07b140c", "cached": false, "eta_s": 51930.0, "run_id": "FS-scene-surface-4b1d0e",
 "scene_job": true, "descriptor": {"job": true, "fn": "fsim_studio.scene:surface_compute", "...": "..."}}
```
Follow it with `/api/jobs/<id>/events`; the `done` result is the SceneSpec. Cached repeats answer 200
with `result`. Show `eta_s` before letting the user wait; offer cancel.

## Stories

### GET /api/stories
```json
[{"name": "rt-single-photons", "title": "Room-temperature single photons: where we stand",
  "n_scenes": 13, "baked_at": "2026-10-07T11:05:04+00:00"}]
```

### GET /api/stories/<name>   (interface I5)
The story JSON plus `current_hash`, `baked` and `stale` (`true` when `baked_hash` differs from the
current hash of fsim_core sources + every `out/**/manifest.json` + every out/ file a number reads +
the card file of every live source + every number's resolved fields). Show a red "stale bake"
indicator when `stale`.

Prose is interpolated: `title`, `claim`, `caveat` and `notes` come back with every `{name}` placeholder
replaced from the baked numbers (story-wide `name` namespace); the templates come back as
`title_template`, `claim_template`, ... Forms: `{name}`, `{name[i]}` (list element), `{name:.3g}`
(format spec), `{name.tag}` (resolved tag), `{name.T}` (temperature a live / csv_cell source was
evaluated at). `prose_errors` lists unresolvable placeholders (empty for the shipped story). Every
scene has non-empty presenter `notes`.

Each number carries `name`, `tag_declared` (the author's tag) and, after bake, `tag` + `tag_source`:
the tag resolved from the source where the source carries one (live: the run's tag_chain; md_regex:
bracket tags inside the matched text; json: the file's `tag_chain`; csv: a `<col>_tag` column), else
the declared tag with `tag_source` saying so.
```json
{"name": "rt-single-photons", "title": "...", "baked_at": "...", "baked_hash": "...", "stale": false,
 "bake_errors": [],
 "scenes": [{"id": "s05", "catalog_id": "rt_edge-verdict", "title": "...and purity fails: VERDICT FAIL",
             "claim": "Even the best eligible row has pulsed g2(0) = 0.9817 (median 0.9879), against the 0.5 gate.",
             "claim_template": "Even the best eligible row has pulsed g2(0) = {g2_min} (median {g2_median}), against the {gate_g2} gate.",
             "numbers": [{"name": "g2_min", "label": "pulsed g2 min (eligible)", "unit": "", "tag": "A",
                          "tag_declared": "A", "tag_source": "declared in the story (the source prints no tag)",
                          "source": {"kind": "verdict", "file": "out/rt_edge/verdict.md", "field": "g2_min"},
                          "value": 0.9817, "raw": "0.9817",
                          "resolved_from": "out/rt_edge/verdict.md:6 VERDICT g2_min"},
                         {"label": "live g2_op at 230 K (headline switches)", "tag": "A",
                          "source": {"kind": "live", "card": "edge-inp-gainp-design", "mode": "headline",
                                     "overrides": {"thermal.T_hs": 230.0, "emission.NA": 0.8, "...": "..."},
                                     "scalar": "g2_op"},
                          "value": 0.9816951150623668, "raw": "0.9816951150623668",
                          "resolved_from": "live fsim_core.device.evaluate(...)"}],
             "visual": {"kind": "chart", "spec": {"preset": "g2_vs_flux_scatter", "campaign": "rt_edge",
                                                  "headline": true, "x": "collected_flux_pulsed_s",
                                                  "y": "g2_pulsed"}},
             "caveat": "FAIL. The gate box is empty for every headline row.", "notes": ""}]}
```
- `value` is a number, string, bool or list (e.g. `"16/768"`, `[0.75, 0.8]`); `raw` is the exact text
  it came from; display `raw`-faithful values (round only in the UI). A number with an `error` key
  failed to resolve: show the error, never a placeholder value.
- `visual.kind` is `chart | scene3d | board | table`; `visual.spec` is a rendering hint (preset name,
  campaign, columns), not data. Chart data comes from the sweep/campaign endpoints.
- `caveat` must stay pinned on screen; `notes` are presenter notes.

### PUT /api/stories/<name>   body = the story JSON (must carry `scenes: [...]`)
Writes `fsim_studio/stories/<name>.json` and returns it as `GET` does. Authoring fields only (H4):
resolved number fields (`value`, `raw`, `resolved_from`, `tag`, `tag_source`, `at_T`, `error`) and the
bake state are stripped from the body and restored from the stored story for numbers whose (scene,
name, source) did not change; a new or changed source stays unresolved until the next bake. When a
scene carries `<field>_template`, that template is stored (the rendered prose is discarded). Set the
author's tag as `tag_declared`.

### POST /api/stories/<name>/bake[?live=0]
Resolves every number's `source` (out/ files, and live fsim_core calls, about 8 s for the default
story; `live=0` skips live sources) and stores `value/raw/resolved_from`, `baked_hash`, `baked_at`,
`bake_errors`. Returns the baked story plus `bake_seconds`. Source kinds are documented at the top of
`fsim_studio/stories.py` (`verdict`, `best`, `md_regex`, `json`, `csv_cell`, `csv_agg`, `commit`,
`live`). File sources must lie under `out/`; `md_regex` patterns are capped at 300 characters and
refused when they carry a nested quantifier or a backreference (H13). Unresolvable prose placeholders
and duplicate number names are reported in `bake_errors`.

## Library: literature parameter cards (fsim_studio/api_library.py, studio-p2b)

Parameter cards are `cards/*.yaml` without a `design:` block. Values come from
`fsim_core.card.load_card`; the two device-tier cards it cannot read (`chatzarakis2023`,
`laferriere2023`: datasets without a per-dataset `source`) are read raw and `loader` says so.
Ranges always come back as `[lo, hi]` with `value: null`; they are never collapsed to a value.

### GET /api/params
`[{name, base, edited, title, device, role, source, paper: {id, label, cite, order}, cls: {code, label,
reason}, tag_counts: {V, DR, E, A}, n_params, n_conflicts, datasets: [{name, tag, n_rows}], dataset_tag_counts,
loader, design, has_device_block}]`. `design` is the design card that derives from this card
("Open as design"), else null. `cls.code` is `C` (calibration), `M` (comparison), `F` (fitted model)
or null, from the physics rules file (`.workers/studio/phase2/physics-brief.md`, "Literature cards rules").

### GET /api/params/<name>
The summary above plus `meta` (whitespace-normalised), `params: [{name, tag, unit, source, value,
range, conflict}]` (`conflict`: the card's own `conflict:` note on that parameter, verbatim, else null; the Library renders it prominently), `datasets: [{name, tag, source, source_note, columns, rows, n_rows, excluded, plot}]`
(`plot` = axis hints: x / y / err / bound / interval columns; labels only), `device_block: [{path,
value | range | text, tag, comment}]` (device-tier cards: the tag is the bracket tag the card writes in
that line's comment, inherited from the nearest tagged parent; `mid*` range centres are listed only as
hidden), `placeholders`, `overlays: {<dataset>: {available, reason?}}`, `editable`, `shipped`.
404 unknown card, 422 design card, 400 bad name.

### GET /api/params/<name>/overlay/<dataset>
The model curve for that dataset, only where the rules file allows one:
`{card, dataset, cls: {code, label}, tag, label, calls: [...], params_source, series: [{name, x, y,
mode?}], bands?: [{name, x, lo, hi}], note?}`. Every array is the return value of the named
fsim_core call (`verify/verify_studio_library.py` checks exact equality against a direct call).
404 `{error, reason}` where no overlay is allowed (and for edited copies). Closed forms and verdicts are
not computed in the Flask layer: the Reischle rho line and the piezo requirement band are read from the
committed `out/phase1/vc_rho_window.csv` and `out/phase3/map_rho_required_300K.csv`, the Zhao verdict word
(`readouts.verdict`, with `verdict_source`) from `out/zhao/zhao_fit_comparison.csv`, and the tau anchor,
F2 inversion and implied dip recovery time from `fsim_core/studio_support.py`.

### PUT /api/params/<name>   body `{params: {<param>: {value} | {lo, hi} | {value, lo, hi}}}`
"Edit copy". `<name>` must end in `-edited` (409 otherwise: shipped cards are never overwritten);
the copy starts from the existing `<name>.yaml` or else the shipped base card, and is written with the
legacy `fsim_gui/app.py` card-editor schema: a non-null `value` replaces the range, otherwise `lo` and
`hi` replace the value with a range. Unknown parameter, `lo > hi` or a non-number -> 400; a card
without a `params` block -> 422; a name that a design saved by Studio already uses (`cards/studio/<name>.yaml`) -> 409. Returns `{name, base, path}`. Requires `X-FSIM-Studio: 1`.

## Model explorers B: phonon, transport + nitride Stark, QD laser (fsim_studio/api_explore_b.py, studio-p2d E2)
Three-layer rule: every array is the return value of the named public `fsim_core` call
(`verify/verify_studio_explore_b.py` checks exact equality against direct calls). Authority:
`.workers/studio/phase2/physics-brief.md` sections 3-5. Each explorer's `/meta` carries the domains,
defaults, tags and the pinned caveat text verbatim. Slow work is a job: the POST returns the usual
`submit_call` shape `{job_id, cached, eta_s, run_id, result?}` (follow it with `/api/jobs/<id>/events`);
the exact job functions are in `jobs.ALLOWED_CALLS` (`fsim_studio.api_explore_b:` `phonon_frames`,
`transport_sweep`, `stark_trace`, `sde_frame`, `sde_li`, `sde_trace`). Out-of-range input -> 400 `{error}`.

### Phonon (`fsim_core.qd_gf`)
- `GET /api/explore/phonon/meta` -> `{domain, defaults, T_frames, alpha_default_ps2, caveat, not_brightness, tags, tag_notes, calls}`.
- `GET /api/explore/phonon/zt?l_xy=&l_z=&alpha=` (live, ~2 ms/point) -> `{T, S_total, Z, alpha_ps2}` on 1..300 K.
- `GET /api/explore/phonon/live?T=&l_xy=&l_z=&alpha=&gamma_zpl=&w=&kappa=|kappa_on=0&F_cav=&delta=` -> a
  `submit_call` envelope `{job_id, cached, eta_s, run_id, result?, T, inputs}` (200 when cached, else 202): a
  ONE-frame `phonon_frames` job on the pool (brief 3(e): T, l_xy, l_z, alpha select precomputed frames, so the
  request thread never integrates the ~1 s cold phi window). `T` must be one of `T_frames` (else 400). The job
  result is `{frames: [<state>], T, inputs}` with the state
  `{T, S_total, Z, sideband_fraction, spectrum: {omega, S}, transmission: {delta, ibm, lorentz},
  point: {t_ibm, t_lorentz, optimism_gap, t_purcell, Z_eff, rate_mult}, gamma_fit, alpha_ps2}`.
- `POST /api/explore/phonon/frames` body = the live arguments except `T` -> job; result
  `{frames: [<state above>, one per T in T_frames], T, inputs}`. One batch task keeps the per-process caches warm.

### Transport and nitride Stark (`fsim_core.transport`, `nitride_transport`, `nitride_stark`, `nitride_levels`)
- `GET /api/explore/transport/meta`, `GET .../point?preset=red|hkust|gaas&I_uA=&T=&n_dot=&aperture=&tau_pulse=&w=&dE_WL=&tau_rad=&E_X_eV=`
  (live, ~18 ms; the `transport.evaluate_injection` row incl. `r_captured`, `r_matrix`, `xi`, `E_U_meV`, `leak_valleys`),
  `GET .../depletion?preset=&T=` (`Diode.depletion` from V_j = 0 to V_bi, plus one flat-band point),
  `GET .../nitride_point?I_uA=&T=&...&polarity=&ext=` (injection row + `bias` from `resolve_bias`).
- `POST /api/explore/transport/sweep` body `{preset, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad}` -> job; result
  `{columns: {I_uA, V_j, V_applied, eta_inj, mu, b_e, P_junction_W, r_dot, r_captured, r_matrix, saturated}, r_dot_max, bg_fit: {max_rel_err, ...}}` (50 points, 1e-3..100 uA).
- `POST /api/explore/transport/stark` body `{orientation: c_plane|a_plane, polarity: +-1, T, ext, screenings?: [0, 0.5, 1]}` ->
  `{traces: [{screening_fraction, job_id, ...}]}`, one job per screening; each result `{rows: [{V_j, E_X_eV, field_kVcm,
  spectroscopy_valid, depletion_regime, flat_band, dE_X_dV_meV_per_V, derivative_valid, ...}], zhang_slope_meV_per_V}`
  at one row per 0.2 V (0..3.2). Slopes come from `stark_derivatives` (never across invalid rows or regime kinks).
  The API never returns a "compatible" verdict: that is `screening_compatibility` in the committed sweep only.

### QD laser, Zhao (`fsim_core.sde`)
- `GET /api/explore/sde/meta` (caveat verbatim, `I_grid`, `F_pumps`), `GET /api/explore/sde/card` (cards/zhao.yaml points [V] and out/zhao CSVs verbatim).
- `POST /api/explore/sde/frames` body `{I_grid?, F_pumps?, n_runs, t_end, seed}` -> `{frames: [{I, F_pump, job_id, ...}]}`; each result
  `{I_over_Ith, F_pump, n_runs, t_end, dt, seed, g2_mean, g2_se, mean_S, R_pump}` (`sde.g2_vs_pump`, seeded, so a frame is reproducible).
- `POST /api/explore/sde/li` -> job `{I, S, rho_ES, rho_GS, residual, I_th}`; `POST /api/explore/sde/trace` body `{I, F_pump, t_end, seed}` -> job
  (decimated `t_ns`, `S`, `n_burn`, `burn_t_ns`, `g2_0`, `g2_0_ext`, `fano`, `fano_ext`).
- `disable_stim`, `beta_sp`, `N_dots`, `E_a_*`, `eps_gain` and `params` in any body -> 400 (tuned Tier-2 parameters and the verify-only estimator control are not exposed).
