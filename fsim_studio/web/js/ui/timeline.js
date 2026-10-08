// Response-animation timeline: a scrubber over n parameter frames, play/pause, speed and loop.
// Frames arrive from the backend one by one (POST /api/animate, `frame` SSE events); the strip
// shows each frame cell as computed (solid), not computed yet (hatched) or failed (dashed rim).
// No numbers are made here. The timeline only says WHERE the playhead is relative to the
// computed frames; the workspace decides what to show there:
//   on a computed frame        -> {kind: "frame", held: i}
//   between two computed frames -> {kind: "interp", held: i, next: i + 1, f}  (display interpolation)
//   between two computed ones the workspace refuses to blend (opts.canInterp returns a reason:
//     invalid g2, runaway, label change ...) -> {kind: "snap", held: the NEARER computed frame, why}
//   touching an uncomputed one  -> {kind: "gap", held: last computed <= pos or null}
// Honesty in playback (spec studio-p2a section 2): the playhead never plays across a gap. On
// reaching one it skips to the next computed frame when one exists, otherwise it waits on the
// last computed frame until the next one arrives (or stops when the job is over).
// Cost: one rAF tick per screen frame while playing, and only while playing (render on
// demand); every tick is timed with the Performance API ("fsim-anim-frame" measures).
import { h, clear, reducedMotion } from "./dom.js";
import { icon } from "./icons.js";

const FPS_1X = 8; // grid frames per second at 1x
export const PERF_MEASURE = "fsim-anim-frame";

function pauseGlyph() {
  return h("span.ico.tl-pause-ico", { "aria-hidden": "true" }, h("span"), h("span"));
}

/**
 * opts: {n, values, label, valueText(v) -> string, onChange(pos, where), onSettle(pos, where),
 *        canInterp?(i, i + 1) -> null | reason string, speed?: 0.5|1|2, loop?: bool}
 */
export function createTimeline(opts) {
  const n = opts.n;
  const state = new Array(n).fill("pending"); // pending | computed | error
  let pos = 0;
  let playing = false;
  let waiting = false;
  let jobOver = false;
  let speed = opts.speed || 1;
  let loop = !!opts.loop;
  let raf = 0;
  let last = 0;
  let nMeasures = 0;
  let limit = null; // last displayable frame (a thermal runaway ends the range); playback stops there
  const end = () => (limit == null ? n - 1 : Math.max(0, Math.min(limit, n - 1)));
  const reduced = reducedMotion();

  const cells = Array.from({ length: n }, (_, i) => h("span.tl-cell", { dataset: { i: String(i), state: "pending" } }));
  const headEl = h("span.tl-head", { "aria-hidden": "true" });
  const track = h("div.tl-track", {
    role: "slider", tabindex: "0",
    "aria-label": `${opts.label || "parameter"} timeline: drag to scrub, space to play`,
    "aria-valuemin": "0", "aria-valuemax": String(n - 1), "aria-valuenow": "0",
  }, h("div.tl-cells", cells), headEl);
  const playBtn = h("button.icon-btn.tl-play", {
    type: "button", "aria-pressed": "false",
    disabled: reduced,
    title: reduced ? "Playback is off while the system asks for reduced motion; drag the scrubber instead" : "Play the computed frames",
    "aria-label": "Play",
    on: { click: () => (playing ? pause() : play()) },
  }, icon("run", { size: 18 }));
  const speedSeg = h("div.segmented.tl-speed", { role: "radiogroup", "aria-label": "Playback speed" },
    [0.5, 1, 2].map((s) => h("button.seg", {
      type: "button", role: "radio", dataset: { speed: String(s) }, "aria-checked": String(s === speed),
      disabled: reduced,
      on: { click: () => setSpeed(s) },
    }, h("span.seg-label", `${s}×`))));
  const loopBtn = h("button.icon-btn.tl-loop", {
    type: "button", "aria-pressed": String(loop), title: "Loop playback", "aria-label": "Loop playback", disabled: reduced,
    on: { click: () => { loop = !loop; loopBtn.setAttribute("aria-pressed", String(loop)); } },
  }, icon("retry", { size: 16 }), h("span.tl-loop-t", "Loop"));
  const posText = h("span.tl-pos", { "aria-live": "off" });
  const countText = h("span.tl-count");
  const legend = h("div.tl-legend", { "aria-hidden": "true" },
    h("span.tl-key", h("span.tl-sw", { dataset: { state: "computed" } }), "computed frame"),
    h("span.tl-key", h("span.tl-sw", { dataset: { state: "pending" } }), "not computed yet"),
    h("span.tl-key", h("span.tl-sw.tl-sw-interp"), "between frames: interp. (display only)"),
    h("span.tl-key", h("span.tl-sw.tl-sw-snap"), "snapped: blending refused"));
  const el = h("div.timeline", { dataset: { playing: "0" } },
    h("div.tl-row", playBtn, track),
    h("div.tl-row.tl-sub", posText, h("span.tl-spacer"), countText, speedSeg, loopBtn),
    legend);

  // ------------------------------------------------------------ where is the playhead
  const isC = (i) => i >= 0 && i < n && state[i] === "computed";
  function where(p) {
    const i = Math.floor(p + 1e-9);
    const f = p - i;
    if (f < 1e-6 || i >= n - 1) {
      const k = Math.min(i, n - 1);
      if (isC(k)) return { kind: "frame", held: k, near: k, f: 0 };
      return { kind: "gap", held: heldBefore(k), at: k };
    }
    if (isC(i) && isC(i + 1)) {
      const why = opts.canInterp?.(i, i + 1) || null;
      if (why) return { kind: "snap", held: f < 0.5 ? i : i + 1, near: f < 0.5 ? i : i + 1, a: i, b: i + 1, f, why };
      return { kind: "interp", held: i, next: i + 1, f, near: f < 0.5 ? i : i + 1 };
    }
    return { kind: "gap", held: heldBefore(i), at: f > 0.5 ? i + 1 : i };
  }
  function heldBefore(i) {
    for (let k = Math.min(i, n - 1); k >= 0; k--) if (isC(k)) return k;
    return null;
  }
  function nextComputed(i) {
    for (let k = Math.max(0, i); k < n; k++) if (isC(k)) return k;
    return null;
  }

  // ------------------------------------------------------------ drawing
  function draw() {
    const pct = n > 1 ? (100 * pos) / (n - 1) : 0;
    headEl.style.left = `${pct}%`;
    track.setAttribute("aria-valuenow", pos.toFixed(2));
    const w = where(pos);
    headEl.dataset.kind = w.kind;
    // a snap shows the value of the frame the readouts come from, not the interpolated parameter value
    const v = w.kind === "snap" ? valueOfFrame(w.near) : valueAt(pos);
    posText.textContent = `${opts.label || ""} ${v}${w.kind === "interp" ? " · interp." : w.kind === "gap" ? " · not computed" : w.kind === "snap" ? ` · snapped to frame ${w.held + 1}/${n}` : ` · frame ${w.held + 1}/${n}`}${waiting ? " · waiting for the next frame" : ""}`;
    track.setAttribute("aria-valuetext", posText.textContent);
    return w;
  }
  function valueAt(p) {
    const i = Math.floor(p), f = p - i;
    const vs = opts.values;
    const v = i >= n - 1 ? vs[n - 1] : vs[i] + (vs[i + 1] - vs[i]) * f;
    return opts.valueText ? opts.valueText(v) : String(v);
  }
  function valueOfFrame(i) {
    const v = opts.values[i];
    return opts.valueText ? opts.valueText(v) : String(v);
  }
  function syncCount() {
    const c = state.filter((s) => s === "computed").length;
    const e = state.filter((s) => s === "error").length;
    countText.textContent = `${c}/${n} frames computed${e ? `, ${e} failed` : ""}`;
  }

  function emit(source) {
    const t0 = performance.now();
    const w = draw();
    opts.onChange?.(pos, w, source);
    const t1 = performance.now();
    try {
      performance.measure(PERF_MEASURE, { start: t0, end: t1 });
      if (++nMeasures > 600) { performance.clearMeasures(PERF_MEASURE); nMeasures = 0; }
    } catch { /* old engines: no measure options */ }
    return w;
  }

  // ------------------------------------------------------------ playback
  function tick(now) {
    raf = 0;
    if (!playing) return;
    const dt = Math.min(0.1, (now - last) / 1000);
    last = now;
    if (!waiting) {
      let p = pos + dt * FPS_1X * speed;
      const i = Math.floor(pos + 1e-9);
      // never play across an uncomputed frame: skip ahead or wait on the last computed one
      const seg = Math.floor(p + 1e-9);
      for (let k = i + 1; k <= Math.min(seg + 1, n - 1); k++) {
        if (k - 1 < p && !isC(k)) {
          const j = nextComputed(k);
          if (j != null) { p = j; break; }
          if (jobOver) { p = k - 1; break; }
          p = k - 1;
          waiting = true;
          break;
        }
      }
      if (p >= end()) {
        p = end();
        if (loop) {
          pos = p;
          emit("play");
          const f = nextComputed(0);
          p = f == null ? 0 : f;
        } else {
          pos = p;
          emit("play");
          pause();
          return;
        }
      }
      pos = p;
      emit("play");
    }
    raf = requestAnimationFrame(tick);
  }
  function play() {
    if (reduced || playing) return;
    if (pos >= end() - 1e-9) pos = nextComputed(0) ?? 0;
    if (!isC(Math.floor(pos))) {
      const j = nextComputed(Math.floor(pos));
      if (j == null) return;
      pos = j;
    }
    playing = true;
    waiting = false;
    el.dataset.playing = "1";
    playBtn.setAttribute("aria-pressed", "true");
    playBtn.setAttribute("aria-label", "Pause");
    clear(playBtn).append(pauseGlyph());
    last = performance.now();
    raf = requestAnimationFrame(tick);
  }
  function pause() {
    const was = playing;
    playing = false;
    waiting = false;
    cancelAnimationFrame(raf);
    raf = 0;
    el.dataset.playing = "0";
    playBtn.setAttribute("aria-pressed", "false");
    playBtn.setAttribute("aria-label", "Play");
    clear(playBtn).append(icon("run", { size: 18 }));
    if (was) { const w = draw(); opts.onSettle?.(pos, w); }
  }
  function setSpeed(s) {
    speed = s;
    for (const b of speedSeg.children) b.setAttribute("aria-checked", String(Number(b.dataset.speed) === s));
  }

  // ------------------------------------------------------------ scrubbing
  function posFromX(x) {
    const r = track.getBoundingClientRect();
    const f = Math.max(0, Math.min(1, (x - r.left) / Math.max(1, r.width)));
    return f * (n - 1);
  }
  let dragging = false;
  track.addEventListener("pointerdown", (e) => {
    if (playing) pause();
    dragging = true;
    track.setPointerCapture?.(e.pointerId);
    seek(posFromX(e.clientX), "scrub");
  });
  track.addEventListener("pointermove", (e) => { if (dragging) seek(posFromX(e.clientX), "scrub"); });
  const endDrag = () => { if (!dragging) return; dragging = false; const w = draw(); opts.onSettle?.(pos, w); };
  track.addEventListener("pointerup", endDrag);
  track.addEventListener("pointercancel", endDrag);
  track.addEventListener("keydown", (e) => {
    let p = null;
    if (e.key === "ArrowRight") p = Math.min(n - 1, Math.floor(pos + 1e-9) + 1);
    else if (e.key === "ArrowLeft") p = Math.max(0, Math.ceil(pos - 1e-9) - 1);
    else if (e.key === "Home") p = 0;
    else if (e.key === "End") p = n - 1;
    else if (e.key === " " || e.key === "Enter") { e.preventDefault(); if (playing) pause(); else play(); return; }
    if (p == null) return;
    e.preventDefault();
    if (playing) pause();
    seek(p, "key");
    opts.onSettle?.(pos, where(pos));
  });

  function seek(p, source = "api") {
    pos = Math.max(0, Math.min(n - 1, p));
    return emit(source);
  }

  function setComputed(i, st = "computed") {
    if (i < 0 || i >= n) return;
    state[i] = st;
    cells[i].dataset.state = st;
    syncCount();
    if (waiting && st === "computed") waiting = false;
    if (!playing) draw();
  }

  syncCount();
  draw();
  return {
    el,
    setComputed,
    seek,
    play,
    pause,
    where: () => where(pos),
    /** Last displayable frame index (null: all). Cells past it are marked, playback stops on it. */
    setLimit(i) {
      limit = i == null ? null : i;
      cells.forEach((c, k) => { if (limit != null && k > limit) c.dataset.beyond = "1"; else delete c.dataset.beyond; });
      if (!playing) draw();
    },
    refresh() { return emit("refresh"); },
    setJobOver(v = true) { jobOver = v; if (v && waiting) { waiting = false; pause(); } },
    get pos() { return pos; },
    get playing() { return playing; },
    get reduced() { return reduced; },
    states: () => state.slice(),
    dispose() { playing = false; cancelAnimationFrame(raf); },
  };
}
