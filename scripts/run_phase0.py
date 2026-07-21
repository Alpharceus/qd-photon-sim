"""Phase 0 driver: load chatzarakis.yaml, run the V-a fit, emit the report bundle.

Usage: python scripts/run_phase0.py [path/to/card.yaml]
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import load_card
from fsim_core.fitting import V_A_TOL, fit_phase0
from fsim_viz.figures import phase0_bundle


def main():
    card_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "cards" / "chatzarakis.yaml"
    card = load_card(card_path)
    print(f"card: {card.name}  ({card_path})")
    ph = card.placeholders()
    if ph:
        print(f"!! PLACEHOLDER entries: {', '.join(ph)}")
        print("!! V-a gate cannot close on this run; replace with digitized values.\n")

    fit = fit_phase0(card)

    print(f"{'T (K)':>7} {'w':>6} {'dx':>6} {'g2 data':>9} {'g2 model':>9} {'resid':>7}  |resid|<= {V_A_TOL}")
    for row, m, r in zip(fit.data, fit.model_g2, fit.residuals):
        ok = "ok" if abs(r) <= V_A_TOL else "FAIL"
        lbl = f"<={row['g2']:.2f}" if row.get("bound") == "upper" else f"{row['g2']:.3f}"
        print(f"{row['T']:>7.0f} {row['w']:>6.2f} {row['dx']:>6.2f} {lbl:>9} {m:>9.3f} {r:>+7.3f}  {ok}")

    gA, gAv, gAt, T_anchor = fit.gamma_anchor
    print(f"\nGamma({T_anchor:.0f} K): model {gA:.2f} meV  vs anchor {gAv}+-{gAt} meV")
    print("fitted parameters:")
    for k, v in fit.params.items():
        print(f"  {k:>10} = {v:.4g}")

    res = phase0_bundle(fit, card_path, ROOT / "out" / "phase0")
    print(f"\nmaster ceiling  rho^2(1-eps)=1/2  ->  Tc = {res['Tc']:.1f} K   tag chain {fit.tag.label}")
    print(f"report bundle -> {res['outdir']}  (figure pdf/svg/png + CSVs + params + card)")
    print("fit status:", "PASS (+-0.03/point)" if fit.passed else "FAIL (+-0.03/point)")
    for n in fit.notes:
        print("NOTE:", n)


if __name__ == "__main__":
    main()
