// FSIM Studio API client: JSON fetch helpers, SSE job runner, offline detection.
// Every number the UI shows comes back through here from the Flask backend
// (fsim_studio/server.py); nothing is computed client-side.

export class ApiError extends Error {
  constructor(message, { status = 0, body = null, offline = false } = {}) {
    super(message);
    this.status = status;
    this.body = body;
    this.offline = offline;
  }
}

const listeners = new Set();
let online = true;

/** Subscribe to reachability changes: fn(online:boolean). */
export function onReachability(fn) {
  listeners.add(fn);
  return () => listeners.delete(fn);
}

function setOnline(v) {
  if (v === online) return;
  online = v;
  for (const fn of listeners) fn(v);
}

async function request(method, url, body) {
  let resp;
  try {
    const headers = {};
    if (body !== undefined) headers["Content-Type"] = "application/json";
    // Mutating requests carry the Studio header; the server refuses them without it (H13).
    if (method !== "GET") headers["X-FSIM-Studio"] = "1";
    resp = await fetch(url, {
      method,
      headers,
      body: body !== undefined ? JSON.stringify(body) : undefined,
    });
  } catch (exc) {
    setOnline(false);
    throw new ApiError("The Studio server is not reachable. Start it with run-studio.bat, then retry.",
      { offline: true });
  }
  setOnline(true);
  let data = null;
  const text = await resp.text();
  if (text) {
    try { data = JSON.parse(text); } catch { data = { error: text }; }
  }
  if (!resp.ok) {
    const msg = (data && (data.error || data.message)) || `${method} ${url} failed (${resp.status})`;
    throw new ApiError(msg, { status: resp.status, body: data });
  }
  return data;
}

export const getJSON = (url) => request("GET", url);
export const postJSON = (url, body) => request("POST", url, body ?? {});
export const putJSON = (url, body) => request("PUT", url, body ?? {});
export const del = (url) => request("DELETE", url);

// ------------------------------------------------------------- cached lookups
const memo = new Map();
function once(key, fn) {
  if (!memo.has(key)) {
    const p = fn().catch((e) => { memo.delete(key); throw e; });
    memo.set(key, p);
  }
  return memo.get(key);
}

export const health = () => getJSON("/api/health");
export const listCards = () => once("cards", () => getJSON("/api/cards"));
export const getCard = (name) => once(`card:${name}`, () => getJSON(`/api/cards/${encodeURIComponent(name)}`));
export const getMeta = () => once("meta", () => getJSON("/api/meta"));
export const listCampaigns = () => once("campaigns", () => getJSON("/api/campaigns"));
export const plotlyTemplate = (mode) => once(`tpl:${mode}`, () => getJSON(`/api/theme/plotly?mode=${mode}`));
export const validate = (design) => postJSON("/api/validate", { design });
export const compareReport = (entries, title) => postJSON("/api/compare/report", { entries, title });
/** Verdict line(s) matched to a design card by the backend (I6): {campaign, lines, match, reason}. */
export const cardVerdict = (name) => once(`verdict:${name}`, () => getJSON(`/api/cards/${encodeURIComponent(name)}/verdict`));

// ------------------------------------------------------------- gates (I3)
// The g2 ceiling, the flux floor and the gate box come ONLY from GET /api/gates (parsed by the
// backend from out/rt_edge/verdict.md and the README gate definition). Never typed in the UI.
let gatesValue = null;
let gatesError = null;
/** Resolves {g2_ceiling: {value, tag, source}, flux_floor: {value, unit, tag, source}} or null. */
export function gates() {
  return once("gates", () => getJSON("/api/gates")).then((g) => { gatesValue = g; gatesError = null; return g; })
    .catch((e) => { gatesError = e.message; return null; });
}
/** Synchronous view of the last loaded gates (null until gates() resolved, or when it failed). */
export function gatesNow() { return gatesValue; }
export function gatesProblem() { return gatesError; }

/**
 * Fetch a SceneSpec. Slow scenes come back as a job descriptor
 * ({scene_job: true, job_id}); those are awaited through the job stream.
 */
export async function sceneSpec(kind, params) {
  const qs = new URLSearchParams(params).toString();
  const r = await getJSON(`/api/scene/${kind}?${qs}`);
  if (r && r.scene_job) {
    if (r.result) return r.result;
    return awaitJob(r.job_id).promise;
  }
  return r;
}

/**
 * Submit a run and follow it to completion.
 * body: {card?, design?, mode, ranged?, T_grid?}
 * hooks: {onState(state), onProgress(k, n), onPartial(result), onSubmitted(resp)}
 * Returns {promise, cancel()}; promise resolves {result, cached, eta_s, run_id, elapsed_s}
 * and rejects with ApiError (body.error carries the server's {kind, message, suggestion, f8_floor}).
 */
export function runJob(body, hooks = {}) {
  let cancelled = false;
  let follower = null;
  let jobId = null;
  let deleted = false;
  // Cancel is honoured even before the POST returns: the DELETE goes out once the id arrives.
  const sendDelete = () => {
    if (deleted || !jobId) return;
    deleted = true;
    del(`/api/jobs/${jobId}`).catch(() => {});
  };
  const t0 = performance.now();

  const promise = (async () => {
    const resp = await postJSON("/api/run", body);
    jobId = resp.job_id;
    if (cancelled) {
      if (!resp.result) sendDelete();
      throw new ApiError("Run cancelled", { body: { error: { kind: "cancelled" } } });
    }
    hooks.onSubmitted?.(resp);
    if (resp.error) {
      throw new ApiError(resp.error.message || "Run failed", { body: { error: resp.error } });
    }
    if (resp.result) {
      return { ...resp, elapsed_s: (performance.now() - t0) / 1000 };
    }
    if (cancelled) throw new ApiError("Run cancelled", { body: { error: { kind: "cancelled" } } });
    follower = awaitJob(jobId, hooks);
    const result = await follower.promise;
    return { ...resp, result, elapsed_s: (performance.now() - t0) / 1000 };
  })();

  return {
    promise,
    get jobId() { return jobId; },
    cancel() {
      cancelled = true;
      follower?.close();
      sendDelete();
    },
  };
}

/** Follow /api/jobs/<id>/events (SSE). Falls back to polling if EventSource errors. */
export function awaitJob(jobId, hooks = {}) {
  let es = null;
  let pollTimer = null;
  let settled = false;
  const close = () => {
    settled = true;
    es?.close();
    clearTimeout(pollTimer);
  };
  const promise = new Promise((resolve, reject) => {
    const fail = (error) => {
      if (settled) return;
      close();
      reject(new ApiError(error?.message || "Run failed", { body: { error } }));
    };
    const done = (result) => {
      if (settled) return;
      close();
      resolve(result);
    };
    const poll = async () => {
      if (settled) return;
      try {
        const snap = await getJSON(`/api/jobs/${jobId}`);
        if (snap.progress) hooks.onProgress?.(snap.progress.k, snap.progress.n);
        if (snap.state === "done") return done(snap.result);
        if (snap.state === "error") return fail(snap.error);
        if (snap.state === "cancelled") return fail({ kind: "cancelled", message: "Run cancelled" });
        hooks.onState?.(snap.state);
      } catch (e) {
        if (e.offline) return fail({ kind: "offline", message: e.message });
      }
      pollTimer = setTimeout(poll, 600);
    };
    try {
      es = new EventSource(`/api/jobs/${jobId}/events`);
    } catch {
      poll();
      return;
    }
    const parse = (ev) => { try { return JSON.parse(ev.data); } catch { return {}; } };
    es.addEventListener("state", (ev) => {
      const d = parse(ev);
      hooks.onState?.(d.state);
      if (d.state === "cancelled") fail({ kind: "cancelled", message: "Run cancelled" });
    });
    es.addEventListener("progress", (ev) => { const d = parse(ev); hooks.onProgress?.(d.k, d.n); });
    es.addEventListener("partial", (ev) => { const d = parse(ev); hooks.onPartial?.(d.result); });
    es.addEventListener("done", (ev) => done(parse(ev).result));
    es.addEventListener("error", (ev) => {
      // A server-sent "error" event carries data; a transport error does not.
      if (ev.data) return fail(parse(ev).error);
      if (settled) return;
      es.close();
      poll();
    });
  });
  return { promise, close };
}
