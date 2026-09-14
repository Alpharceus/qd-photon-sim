"""Independent contract checks for nanowire device dispatch and gates."""
import math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.device import DeviceDesign, DriveBlock, RetentionBlock, CavityBlock, EmissionBlock, ThermalBlock, DotBlock, evaluate


def card(family="horizontal_as_built", regime="rectangular"):
    core = 12.5 if family == "horizontal_as_built" else 80.0
    dotr = core if family == "horizontal_as_built" else 12.5
    return DeviceDesign(platform="ingan_gan_nanowire", dot=DotBlock(linewidth="anchored", lineshape="lorentzian", gamma300=3.0), ret=RetentionBlock(mode="nitride_confinement"),
      drive=DriveBlock(mode="EL-transport",I_uA=.02,duty=.008,rep_rate_hz=80e6,b_res=.1,
        diode={"preset":"nitride-nanowire","tau_pulse_ns":.1}), thermal=ThermalBlock(T_hs=300.), cavity=CavityBlock(enabled=False),
      emission=EmissionBlock(type="nanowire"), nitride={
      "nanowire":{"family":family,"core_radius_nm":core,"outer_radius_nm":core,"strain_bound":"relaxed","barrier_left_nm":15.,"barrier_right_nm":15.},
      "dot":{"radius_nm":dotr,"height_nm":2.,"x_in":.4,"strain_fraction":0.,"tau_rad0_ns":1.,"tau_cap_ps":10.},
      "surface":{"occupied_dot_access":.05},"photonics":{"NA":.5},"wire_thermal":{"Rth_K_W":3.1e9,"eta_total":1.0,"f_Rs_local":1.0,"R_s_ohm":2.38e9,"C_parasitic_F":0.0},
      "injector":{}},)


def main():
    checks=[]
    for family in ("horizontal_as_built","vertical_photonic"):
      d=card(family); r=evaluate(d,T_grid=[230.,300.]); checks += [r["scalars"]["platform"]=="ingan_gan_nanowire",len(r["curves"]["g2_op"])==2,math.isfinite(r["scalars"]["T_j"])]
    s=evaluate(card(regime="deterministic_pair"))["scalars"]
    checks += ["set_E_C_meV" in s, "rti_feasible" in s, s["rti_qualified"] <= s["optical_pass"]]
    try: evaluate(card("horizontal_as_built").__class__(platform="ingan_gan_nanowire"))
    except ValueError: checks.append(True)
    else: checks.append(False)
    assert all(checks), checks
    print(f"{len(checks)}/{len(checks)} nitride nanowire device checks passed")

if __name__ == "__main__": main()
