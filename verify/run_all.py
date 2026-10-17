"""The single gate (R2 test audit 1): run EVERY verifier and exit 0 iff all pass.

  python verify/run_all.py [--fast] [--only <glob>] [--port 8840] [--timeout SECONDS] [--list]

The verifier set is the glob verify/verify_*.py + verify/verify_*.cjs (plus verify/audit_physics.py, named in CLAUDE.md), so a new
verifier cannot be left out. Python verifiers run first, one at a time, with OMP_NUM_THREADS=2. Unless --fast, ONE private Studio
server is then started (temporary user-card folder and job cache, FSIM_STUDIO_WORKERS=2, port from --port, default 8840) and the
.cjs verifiers run strictly one after another against it; the server is killed by PID (taskkill /T /F on Windows). A verifier
counts as passed only if its exit code is 0 AND it printed an "N/N ..." line with N == N; a verifier that printed no such line
fails ("did not print its N/N line"). Output: a table  name | N/N | seconds | status.

  --fast      python verifiers only (no server, no browser)
  --only G    only verifiers whose file name matches the glob G (fnmatch, e.g. "verify_compose_*")
  --timeout   per-verifier limit in seconds (default 3600)
"""
from __future__ import annotations

import fnmatch
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "verify"
N_RE = re.compile(r"(\d+)\s*/\s*(\d+)")


def discover(only=None):
    py = sorted(VERIFY.glob("verify_*.py")) + [VERIFY / "audit_physics.py"]
    cjs = sorted(VERIFY.glob("verify_*.cjs"))
    keep = lambda p: only is None or fnmatch.fnmatch(p.name, only) or fnmatch.fnmatch(p.stem, only)   # noqa: E731
    return [p for p in py if p.is_file() and keep(p)], [p for p in cjs if keep(p)]


def parse_counts(text: str):
    """The last 'N/N' pair on a line that reads like a summary ('... checks passed', 'matched', or a bare 'N/N name')."""
    best = None
    for line in text.splitlines():
        if not re.search(r"passed|matched|checks|tests", line, re.I):
            continue
        m = N_RE.search(line)
        if m:
            best = (int(m.group(1)), int(m.group(2)))
    if best is None:
        for line in reversed(text.splitlines()):
            m = N_RE.search(line)
            if m:
                return int(m.group(1)), int(m.group(2))
    return best


def run_one(cmd, env, timeout):
    t0 = time.time()
    try:
        cp = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout)
        out, rc = (cp.stdout or "") + "\n" + (cp.stderr or ""), cp.returncode
    except subprocess.TimeoutExpired as e:
        out, rc = ((e.stdout or b"").decode("utf-8", "replace") if isinstance(e.stdout, bytes) else (e.stdout or "")), -999
    dt = time.time() - t0
    counts = parse_counts(out)
    if rc == -999:
        status = "FAIL (timeout)"
    elif rc != 0:
        status = f"FAIL (exit {rc})"
    elif counts is None:
        status = "FAIL (did not print its N/N line)"
    elif counts[0] != counts[1]:
        status = "FAIL (N != N)"
    else:
        status = "ok"
    return counts, dt, status, out


def kill_tree(pid):
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"], capture_output=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
        else:
            os.kill(pid, 15)
    except Exception:                                                          # noqa: BLE001
        pass


def start_server(port):
    user_cards = tempfile.mkdtemp(prefix="fsim_runall_cards_")
    cache = tempfile.mkdtemp(prefix="fsim_runall_cache_")
    env = {**os.environ, "FSIM_STUDIO_USER_CARDS": user_cards, "FSIM_STUDIO_CACHE": cache, "FSIM_STUDIO_WORKERS": "2",
           "OMP_NUM_THREADS": "2", "MSYS_NO_PATHCONV": "1"}
    flags = 0x08000000 if os.name == "nt" else 0                              # CREATE_NO_WINDOW
    proc = subprocess.Popen([sys.executable, "-m", "fsim_studio", "--no-window", "--port", str(port)], cwd=ROOT, env=env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=flags)
    base = f"http://127.0.0.1:{port}"
    t0 = time.time()
    while time.time() - t0 < 90:
        if proc.poll() is not None:
            raise RuntimeError(f"Studio server exited with {proc.returncode} on port {port} (port busy?)")
        try:
            with urllib.request.urlopen(f"{base}/api/health", timeout=2) as r:
                if r.status == 200:
                    return proc, base, env, (user_cards, cache)
        except Exception:                                                      # noqa: BLE001
            time.sleep(0.4)
    kill_tree(proc.pid)
    raise RuntimeError("Studio server did not answer /api/health in 90 s")


def main(argv):
    def opt(name, default=None):
        return argv[argv.index(name) + 1] if name in argv else default
    fast = "--fast" in argv
    only = opt("--only")
    port = int(opt("--port", "8840"))
    timeout = float(opt("--timeout", "3600"))
    py, cjs = discover(only)
    if "--list" in argv:
        print("\n".join(p.name for p in py + ([] if fast else cjs)))
        return 0
    base_env = {**os.environ, "OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "PYTHONUNBUFFERED": "1"}
    rows = []

    def record(name, counts, dt, status):
        rows.append((name, f"{counts[0]}/{counts[1]}" if counts else "-", dt, status))
        print(f"{name:<44} {rows[-1][1]:>9} {dt:8.1f}s  {status}", flush=True)

    t_all = time.time()
    print(f"{'verifier':<44} {'N/N':>9} {'seconds':>9}  status", flush=True)
    for p in py:
        counts, dt, status, out = run_one([sys.executable, str(p.relative_to(ROOT))], base_env, timeout)
        record(p.name, counts, dt, status)
        if status != "ok":
            print("    last output: " + " | ".join(l.strip() for l in out.strip().splitlines()[-3:])[:400], flush=True)
            for fl in [l.strip() for l in out.splitlines() if l.lstrip().startswith(("FAIL", "[FAIL"))][:12]:
                print("    " + fl[:300], flush=True)
    if not fast and cjs:
        try:
            proc, url, senv, tmp = start_server(port)
        except Exception as e:                                                 # noqa: BLE001
            print(f"cannot start the Studio server: {e}")
            for p in cjs:
                record(p.name, None, 0.0, "FAIL (no server)")
            proc = None
        if proc is not None:
            try:
                benv = {**senv, "PLAYWRIGHT_BROWSERS_PATH": os.environ.get("PLAYWRIGHT_BROWSERS_PATH", "")} if os.environ.get("PLAYWRIGHT_BROWSERS_PATH") else dict(senv)
                for p in cjs:                                                  # strictly one browser verifier at a time
                    counts, dt, status, out = run_one(["node", str(p.relative_to(ROOT)), url], benv, timeout)
                    record(p.name, counts, dt, status)
                    if status != "ok":
                        print("    last output: " + " | ".join(l.strip() for l in out.strip().splitlines()[-3:])[:400], flush=True)
                        for fl in [l.strip() for l in out.splitlines() if l.lstrip().startswith(("FAIL", "[FAIL"))][:12]:
                            print("    " + fl[:300], flush=True)
            finally:
                kill_tree(proc.pid)
                for d in tmp:
                    shutil.rmtree(d, ignore_errors=True)
    bad = [r for r in rows if r[3] != "ok"]
    print("-" * 80)
    print(f"{len(rows) - len(bad)}/{len(rows)} verifiers passed in {time.time() - t_all:.0f} s"
          + (": FAILED " + ", ".join(r[0] for r in bad) if bad else ""))
    return 0 if rows and not bad else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
