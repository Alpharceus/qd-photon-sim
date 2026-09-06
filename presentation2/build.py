"""presentation2/build.py -- orchestrates the presentation2 build per the
renderer contract in presentation2/SCHEMA.md:
validate -> figures -> equations -> pptx -> html -> office validate.

The "figure path exists" check in validate_sections.py can only pass once
the figure scripts have actually run, so this does a fast structural-only
validation pass first (catches schema errors before spending time running
figure scripts), runs the figure scripts, then re-validates in full
(including figure-path existence and repo_numbers) before rendering.

CLI: python presentation2/build.py [--include-sample] [--sections-dir DIR]
Exits non-zero on any failure.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

sys.path.insert(0, str(HERE))
import mathimg  # noqa: E402
import validate_sections  # noqa: E402


def run(cmd: list) -> None:
    print(f"$ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=str(ROOT))
    if result.returncode != 0:
        print(f"FAILED (exit {result.returncode}): {' '.join(str(c) for c in cmd)}", file=sys.stderr)
        sys.exit(result.returncode)


def run_figure_scripts(sections: list) -> None:
    """Run every figure script referenced by `sections` whose PNG is missing
    or older than the script."""
    by_script: dict[str, list[Path]] = {}
    for sec in sections:
        for slide in sec.get("slides", []):
            fig = slide.get("figure")
            if fig and "script" in fig and "path" in fig:
                by_script.setdefault(fig["script"], []).append(ROOT / fig["path"])

    for rel_script in sorted(by_script):
        script_path = ROOT / rel_script
        png_paths = by_script[rel_script]
        stale = not all(p.exists() for p in png_paths)
        if not stale and script_path.exists():
            newest_png_mtime = min(p.stat().st_mtime for p in png_paths)
            stale = script_path.stat().st_mtime > newest_png_mtime
        if stale:
            run([sys.executable, str(script_path)])
        else:
            print(f"skip (up to date): {rel_script}")


def render_equations(sections: list) -> None:
    count = 0
    for sec in sections:
        for slide in sec.get("slides", []):
            for eq in slide.get("equations", []) or []:
                mathimg.render_equation(eq["latex"])
                count += 1
    print(f"equations: {count} rendered/cached under presentation2/figures/out/eq/")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-sample", action="store_true",
                         help="include presentation2/sections/00_sample.json in the build")
    parser.add_argument("--sections-dir", default=None,
                         help="directory of NN_*.json section files (default: presentation2/sections)")
    args = parser.parse_args()

    print("== 1/6 validate (structural) ==")
    ok, sections = validate_sections.validate_all(
        section_dir=args.sections_dir, include_sample=args.include_sample, check_figures=False,
    )
    if not ok:
        sys.exit(1)

    print("== 2/6 figures ==")
    run_figure_scripts(sections)

    print("== 3/6 validate (figures + repo_numbers) ==")
    ok, sections = validate_sections.validate_all(
        section_dir=args.sections_dir, include_sample=args.include_sample, check_figures=True,
    )
    if not ok:
        sys.exit(1)

    print("== 4/6 equations ==")
    render_equations(sections)

    print("== 5/6 pptx + html ==")
    node_cmd = ["node", str(HERE / "render_pptx.js")]
    py_cmd = [sys.executable, str(HERE / "render_html.py")]
    if args.include_sample:
        node_cmd.append("--include-sample")
        py_cmd.append("--include-sample")
    if args.sections_dir:
        node_cmd += ["--sections-dir", str(args.sections_dir)]
        py_cmd += ["--sections-dir", str(args.sections_dir)]
    run(node_cmd)
    run(py_cmd)

    print("== 6/6 office validate ==")
    validate_script = Path.home() / ".claude" / "skills" / "pptx" / "scripts" / "office" / "validate.py"
    pptx_out = ROOT / "out" / "presentation2" / "qd_physics_2h.pptx"
    run([sys.executable, str(validate_script), str(pptx_out)])

    print("BUILD OK")


if __name__ == "__main__":
    main()
