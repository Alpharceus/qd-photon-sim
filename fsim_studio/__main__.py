"""python -m fsim_studio [--browser] [--no-window] [--port 8765]

Default: a native pywebview window "FSIM Studio" (1600x1000, min 1280x800).
--browser opens the default browser instead; --no-window only serves.
The server binds 127.0.0.1 only.
"""
from __future__ import annotations

import argparse
import sys
import threading
import time
import webbrowser

from . import ROOT

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _serve(app, port: int) -> None:
    app.run(host="127.0.0.1", port=port, threaded=True, use_reloader=False, debug=False)


def _wait_up(url: str, timeout: float = 30.0) -> bool:
    import urllib.request
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            with urllib.request.urlopen(url + "/api/health", timeout=1.0):
                return True
        except Exception:  # noqa: BLE001 -- polling until the server answers
            time.sleep(0.1)
    return False


def _shutdown(app) -> None:
    jm = app.config.get("JOB_MANAGER")
    if jm is not None:
        jm.terminate()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m fsim_studio", description=__doc__.splitlines()[0])
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--browser", action="store_true", help="open the default browser instead of a window")
    ap.add_argument("--no-window", action="store_true", help="serve only (no window, no browser)")
    args = ap.parse_args(argv)

    from .server import create_app
    app = create_app(port=args.port)
    url = f"http://127.0.0.1:{args.port}"

    try:
        return _run(app, args, url)
    finally:
        _shutdown(app)


def _run(app, args, url) -> int:
    if args.no_window:
        print(f"FSIM Studio serving on {url}", flush=True)
        try:
            _serve(app, args.port)
        except KeyboardInterrupt:
            pass
        return 0

    t = threading.Thread(target=_serve, args=(app, args.port), daemon=True)
    t.start()
    _wait_up(url)
    if args.browser:
        webbrowser.open(url + "/")
        print(f"FSIM Studio serving on {url} (Ctrl+C to stop)", flush=True)
        try:
            while t.is_alive():
                t.join(0.5)
        except KeyboardInterrupt:
            pass
        return 0
    try:
        import webview
    except ImportError:
        print("pywebview not installed; opening the default browser instead", flush=True)
        webbrowser.open(url + "/")
        try:
            while t.is_alive():
                t.join(0.5)
        except KeyboardInterrupt:
            pass
        return 0
    webview.create_window("FSIM Studio", url + "/", width=1600, height=1000,
                          min_size=(1280, 800))
    webview.start()
    return 0


if __name__ == "__main__":
    sys.exit(main())
