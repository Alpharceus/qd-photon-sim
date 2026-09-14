"""Mechanical card checks for the four nitride nanowire scenarios."""
import os
import sys
import yaml

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, ROOT)
from fsim_core.device import DeviceDesign, evaluate

NAMES = (
    "nitride-nanowire-horizontal-pulse-design.yaml",
    "nitride-nanowire-horizontal-set-design.yaml",
    "nitride-nanowire-vertical-pulse-design.yaml",
    "nitride-nanowire-vertical-set-design.yaml",
)

def main():
    checks = 0
    failures = []
    designs = []
    for name in NAMES:
        path = os.path.join(ROOT, "cards", name)
        raw = yaml.safe_load(open(path, encoding="utf-8"))
        d = DeviceDesign.load(path)
        designs.append(d)
        checks += 1
        if d.platform != "ingan_gan_nanowire": failures.append(name + ": platform")
        checks += 1
        if raw["design"]["drive"]["diode"]["preset"] != "nitride-nanowire": failures.append(name + ": preset")
        for T in (230.0, 300.0):
            for bound in ("unrelaxed", "relaxed"):
                d.nitride["nanowire"]["strain_bound"] = bound
                d.nitride["dot"]["strain_fraction"] = 1.0 if bound == "unrelaxed" else 0.0
                try:
                    result = evaluate(d, T_grid=[T])
                    checks += 1
                    if result["scalars"].get("platform") != "ingan_gan_nanowire": failures.append(name + ": evaluation")
                except Exception as exc:
                    failures.append(name + ": " + str(exc))
    checks += 1
    if designs[0].drive.cycle_loading == designs[1].drive.cycle_loading: failures.append("pulse/set loading")
    checks += 1
    if designs[2].nitride["nanowire"]["family"] == designs[0].nitride["nanowire"]["family"]: failures.append("families")
    print(f"{checks - len(failures)}/{checks} nitride nanowire card checks passed")
    if failures:
        for failure in failures: print("FAIL: " + failure)
        return 1
    return 0

if __name__ == "__main__": sys.exit(main())
