// Shared by the FSIM Studio browser verifiers: start (and stop) a private Studio server whose
// saved-card folder and job cache are temporary, so a verifier never writes into the user's real
// cards/studio/ (review finding 6).
//
//   const { startStudioServer, requireTempUserCards } = require("./studio_server.cjs");
//   const srv = await startStudioServer(ROOT, 8791);   // { base, userCards, cache, stop() }
//
// requireTempUserCards(ROOT) is the guard for a verifier pointed at an already-running server: the
// FSIM_STUDIO_USER_CARDS this process sees (and the server must have been started with) has to be
// set and must not be the repository's cards/studio.
const path = require("path");
const fs = require("fs");
const os = require("os");
const { spawn, spawnSync } = require("child_process");

function requireTempUserCards(root) {
  const uc = process.env.FSIM_STUDIO_USER_CARDS || "";
  const real = path.resolve(root, "cards", "studio");
  if (!uc || path.resolve(uc) === real) {
    console.error("refusing to run: set FSIM_STUDIO_USER_CARDS to a temporary folder, the same one the server was started with\n"
      + "(or run this verifier with no URL and it starts its own server). Saving through a server without it would write\n"
      + `into ${real}.`);
    process.exit(2);
  }
  return path.resolve(uc);
}

async function waitHealth(base, ms) {
  const t0 = Date.now();
  for (;;) {
    try {
      const r = await fetch(`${base}/api/health`);
      if (r.ok) return;
    } catch { /* not up yet */ }
    if (Date.now() - t0 > ms) throw new Error(`Studio server did not answer ${base}/api/health in ${ms} ms`);
    await new Promise((r) => setTimeout(r, 300));
  }
}

async function startStudioServer(root, port = 8791) {
  const userCards = fs.mkdtempSync(path.join(os.tmpdir(), "fsim_studio_user_cards_"));
  const cache = fs.mkdtempSync(path.join(os.tmpdir(), "fsim_studio_cache_"));
  const env = { ...process.env, FSIM_STUDIO_USER_CARDS: userCards, FSIM_STUDIO_CACHE: cache, MSYS_NO_PATHCONV: "1" };
  const child = spawn(process.env.PYTHON || "python", ["-m", "fsim_studio", "--no-window", "--port", String(port)], { cwd: root, env, stdio: "ignore", windowsHide: true });
  const base = `http://127.0.0.1:${port}`;
  const stop = () => {
    try {
      if (process.platform === "win32") spawnSync("taskkill", ["/pid", String(child.pid), "/T", "/F"], { stdio: "ignore" });
      else child.kill("SIGTERM");
    } catch { /* already gone */ }
    for (const d of [userCards, cache]) { try { fs.rmSync(d, { recursive: true, force: true }); } catch { /* temp */ } }
  };
  try { await waitHealth(base, 60000); } catch (e) { stop(); throw e; }
  return { base, userCards, cache, stop };
}

module.exports = { startStudioServer, requireTempUserCards };
