"""Checks for the drive registry entries of the composable engine (fsim_core/compose/drives.py, P4).

Run: python verify/verify_compose_drives.py   (exit code 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import drive_mech                                              # noqa: E402
from fsim_core import resonant_source as rs                                   # noqa: E402
from fsim_core.compose import drives as D                                     # noqa: E402
from fsim_core.compose import registry as R                                   # noqa: E402
from fsim_core.compose.registry import Material, Registry, Structure          # noqa: E402
from fsim_core.compose.spec import validate_spec                              # noqa: E402
from fsim_core.compose.types import (DrivePlan, EmitterModel, Pulse, Rates,   # noqa: E402
                                     StructureResponse)

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not ok else ""))


def raises(exc, fn, *a, contains=(), **kw):
    try:
        fn(*a, **kw)
    except exc as e:
        msg = str(e)
        return all(c in msg for c in contains), msg
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"
    return False, "no exception"


def em_with(transitions):
    return EmitterModel(E_X_meV=1380.0, lambda_nm=897.4, n_host=3.5, gamma0_ns=1.5, k_nr_ns=0.0,
                        gamma_star_ns=0.1, phonon=None, delta_xx_meV=-2.0, gamma_xx_ratio=2.0,
                        fss_ueV=0.0, transitions=frozenset(transitions), dipole="in-plane",
                        dEdT_meV_K=0.0, sigma_sd_ueV=0.0, tau_sd_ns=1.0)


EM = em_with({"X", "XX"})
EM_NO_XX = em_with({"X"})
RT = Rates(Gamma_ns=1.5, gamma_star_ns=0.1, k_nr_ns=0.1, L=1.0, F_P_eff=0.0, beta=0.5)

NEEDS_PHYSICS = ()                                         # P9 made the three phonon drives ready, FX2c arp
THROWAWAY = D._placeholder("test_needs_physics", "t", "s", "reason x", (D._GAUSS,), ("r",), D.OBS_COHERENT_PULSED)
PHONON_READY = ("phonon_assisted", "super", "dichromatic", "arp")   # drives_phonon.py
READY = ("resonant_pulse", "tpe", "cw_resonant", "above_band", "p_shell", "electrical_pulsed", "electrical_cw")

# ---------------------------------------------------------------- 1. registration / payload
reg = R.REGISTRY
STARK_GATED = ("stark_gated_resonant",)                    # p10e: registered in drives.py, planned in section 3 of the engine checks
check(f"registered drive kinds == the declared set ({len(R.kinds('drive'))} kinds, exact equality with the registry; R2 audit 16)",
      set(reg.kinds("drive")) == set(READY) | set(NEEDS_PHYSICS) | set(PHONON_READY) | set(STARK_GATED)
      and len(reg.kinds("drive")) == len(set(reg.kinds("drive"))),
      f"registry-only: {sorted(set(reg.kinds('drive')) - set(READY) - set(PHONON_READY) - set(STARK_GATED))}; "
      f"declared-only: {sorted((set(READY) | set(PHONON_READY) | set(STARK_GATED)) - set(reg.kinds('drive')))}")
entries = {k: reg.get("drive", k) for k in READY + NEEDS_PHYSICS + PHONON_READY}
check("every drive has full ParamSpecs and a non-empty observables set",
      all(e.params and e.observables for e in entries.values()))

# private registry so the payload check does not depend on other pieces being importable
priv = Registry()
for e in D.DRIVES:
    priv.register(e)
payload = R.registry_payload(priv)
try:
    txt = json.dumps(payload)
    ok_json = True
except TypeError as e:
    ok_json, txt = False, str(e)
check("registry payload (drives) is JSON-serialisable", ok_json, txt[:200])
pd = {d["kind"]: d for d in payload["drives"]}
check("payload carries observables per drive; every shipped drive is ready (FX2c: no needs_physics left)",
      all(pd[k]["observables"] for k in pd) and all(pd[k]["status"] == "ready" for k in READY + PHONON_READY))
check("a needs_physics entry (test-registered throwaway, not shipped) carries a reason and plan() raises "
      "NotImplementedError with it",
      THROWAWAY.status == "needs_physics" and THROWAWAY.reason and raises(
          NotImplementedError, THROWAWAY.plan, {}, EM, RT, contains=("awaiting physics", THROWAWAY.reason))[0])
check("phonon drive ParamSpecs cover the brief ranges (phonon_assisted blue 0.8-1.3 meV incl. Thomas 0.58, 5-20 pi; "
      "SUPER Bracht Table 1 row 2 defaults -5.0/-13.56 meV, Piccinini -3.5/-7.5 meV inside the ranges)",
      (lambda pa, sp: pa.param("detuning_meV").lo <= 0.58 and pa.param("detuning_meV").hi >= 1.3
       and pa.param("area_pi").lo == 5 and pa.param("area_pi").hi == 20
       and sp.param("delta1_meV").default == -5.0 and sp.param("delta2_meV").default == -13.56
       and sp.param("delta1_meV").lo <= -3.5 <= sp.param("delta1_meV").hi and sp.param("delta2_meV").lo <= -7.5 <= sp.param("delta2_meV").hi)(
          entries["phonon_assisted"], entries["super"]))
check("cw drives exclude per-pulse brightness; pulsed drives include it",
      not ({"B_first_lens", "B_fibre"} & (entries["cw_resonant"].observables | entries["electrical_cw"].observables))
      and "B_first_lens" in entries["resonant_pulse"].observables)
check("electrical mechanism choices come from drive_mech.roster() (all roster names covered)",
      set(D.MECH_PULSED) | set(D.MECH_CW) == {i.mechanism for i in drive_mech.roster()}
      and not set(D.MECH_PULSED) & set(D.MECH_CW)
      and tuple(entries["electrical_pulsed"].param("mechanism").choices) == D.MECH_PULSED
      and tuple(entries["electrical_cw"].param("mechanism").choices) == D.MECH_CW)

# ---------------------------------------------------------------- 2. plans well-formed
plans = {}
plans["resonant_pulse"] = entries["resonant_pulse"].plan({}, EM, RT)
plans["tpe"] = entries["tpe"].plan({}, EM, RT)
plans["cw_resonant"] = entries["cw_resonant"].plan({}, EM, RT)
plans["above_band"] = entries["above_band"].plan({}, EM, RT)
plans["p_shell"] = entries["p_shell"].plan({}, EM, RT)
plans["electrical_pulsed"] = entries["electrical_pulsed"].plan({}, EM, RT)
plans["electrical_cw"] = entries["electrical_cw"].plan({}, EM, RT)
check("all default plans are DrivePlan instances", all(isinstance(p, DrivePlan) for p in plans.values()))
rp = plans["resonant_pulse"]
check("resonant_pulse plan: coherent / pulsed / X, one Pulse with positive width, area_rad = area_pi*pi",
      rp.family == "coherent" and rp.timing == "pulsed" and rp.target == "X" and len(rp.pulses) == 1
      and isinstance(rp.pulses[0], Pulse) and rp.pulses[0].fwhm_ps > 0
      and rp.pulses[0].area_rad == 1.0 * math.pi and rp.f_rep_MHz == 80.0)
p25 = entries["resonant_pulse"].plan({"area_pi": 2.5, "shape": "rect", "fwhm_ps": 7.0}, EM, RT)
check("area_rad = area_pi*pi exactly for 2.5 pi, rect 7 ps",
      p25.pulses[0].area_rad == 2.5 * math.pi and p25.pulses[0].shape == "rect"
      and p25.pulses[0].fwhm_ps == 7.0 and "sigma_ps" not in p25.inputs)
tp = plans["tpe"]
check("tpe plan: coherent / pulsed / XX_tpe / needs cascade / relative area flagged",
      tp.target == "XX_tpe" and tp.needs == frozenset({"cascade"}) and tp.family == "coherent"
      and tp.inputs["area_relative_to_tpe_pi"] and tp.pulses[0].area_rad == math.pi
      and tp.inputs["delta_xx_meV"] == -2.0)
cw = plans["cw_resonant"]
check("cw_resonant plan: timing cw, no f_rep, omega_cw_ns = hbar*Omega(1 ueV)/hbar",
      cw.timing == "cw" and cw.f_rep_MHz is None and cw.pulses == ()
      and abs(cw.omega_cw_ns - 1e-3 / 6.582119569e-4) < 1e-12 * cw.omega_cw_ns)
cws = entries["cw_resonant"].plan({"saturation_s": 2.0}, EM, RT)
g1 = RT.Gamma_ns + RT.k_nr_ns
g2c = 0.5 * (g1 + RT.gamma_star_ns)
check("cw_resonant saturation s = Omega^2 T1 T2 round trip (T1=1/(G+k), T2=2/(G+k+g*))",
      abs(cws.omega_cw_ns ** 2 * (1 / g1) * (1 / g2c) - 2.0) < 1e-12)
ab = plans["above_band"]
check("above_band plan: incoherent / pulsed, p*tau_on = P/P_sat, relax = tau_cap, no cascade by default",
      ab.family == "incoherent" and ab.timing == "pulsed" and ab.tau_on_ns > 0
      and abs(ab.pump_rate_ns * ab.tau_on_ns - 1.0) < 1e-12 and abs(ab.relax_ns - 0.03) < 1e-15
      and ab.needs == frozenset() and ab.inputs["reexcitation"] is True)
ab2 = entries["above_band"].plan({"cap2": True, "cw": True, "pump": 2.0}, EM, RT)
check("above_band cap2 needs cascade; cw variant has timing cw, f_rep None, tau_on 0",
      ab2.needs == frozenset({"cascade"}) and ab2.timing == "cw" and ab2.f_rep_MHz is None
      and ab2.tau_on_ns == 0.0 and abs(ab2.pump_rate_ns - 2.0 * g1) < 1e-12)
ps = plans["p_shell"]
check("p_shell plan: no re-excitation flag, b_dark and tau_relax carried",
      ps.inputs["reexcitation"] is False and ps.inputs["b_dark"] == 0.0 and abs(ps.relax_ns - 0.01) < 1e-15)
check("bad inputs rejected (fwhm 0, f_rep -1, window 0, b_dark 2)",
      raises(ValueError, entries["resonant_pulse"].plan, {"fwhm_ps": 0.0}, EM, RT)[0]
      and raises(ValueError, entries["resonant_pulse"].plan, {"f_rep_MHz": -1.0}, EM, RT)[0]
      and raises(ValueError, entries["above_band"].plan, {"window_ps": 0.0}, EM, RT)[0]
      and raises(ValueError, entries["p_shell"].plan, {"b_dark": 2.0}, EM, RT)[0])
check("resonant_pulse target X+ needs 'X+' in transitions",
      raises(ValueError, entries["resonant_pulse"].plan, {"target": "X+"}, EM, RT, contains=("X+",))[0]
      and entries["resonant_pulse"].plan({"target": "X+"}, em_with({"X", "X+"}), RT).target == "X+")
cp = entries["resonant_pulse"].plan({"cross_pol": True, "laser_leak": 0.01}, EM, RT)
check("cross_pol carries laser_leak in plan.inputs (0 when not cross-polarised)",
      cp.inputs["cross_pol"] and cp.inputs["laser_leak"] == 0.01
      and entries["resonant_pulse"].plan({"laser_leak": 0.01}, EM, RT).inputs["laser_leak"] == 0.0)

# ---------------------------------------------------------------- 3. sign convention
blue = entries["resonant_pulse"].plan({"detuning_meV": 0.83}, EM, RT)
red = entries["resonant_pulse"].plan({"detuning_meV": -0.83}, EM, RT)
check("resonant_pulse detuning_meV=+0.83 (laser blue, paper Delta) -> Pulse.detuning_meV = repo det_X = -0.83",
      blue.pulses[0].detuning_meV == -0.83 and blue.inputs["det_X_meV"] == -0.83
      and blue.inputs["Delta_paper_meV"] == 0.83)
check("detuning_meV=-0.83 (laser red) -> det_X = +0.83; default 0 -> det_X == 0",
      red.pulses[0].detuning_meV == 0.83 and rp.pulses[0].detuning_meV == 0.0)
check("cw and tpe use the same sign flip",
      entries["cw_resonant"].plan({"detuning_meV": 0.2}, EM, RT).inputs["det_X_meV"] == -0.2
      and entries["tpe"].plan({"detuning_meV": 0.2}, EM, RT).pulses[0].detuning_meV == -0.2)

# ---------------------------------------------------------------- 4. Gaussian mapping
for fw in (0.7, 3.0, 21.4):
    pl = entries["resonant_pulse"].plan({"fwhm_ps": fw, "shape": "gauss"}, EM, RT)
    sig_ns = rs.gaussian_sigma(fw * 1e-3)            # resonant_source works in ns
    check(f"gauss intensity FWHM {fw} ps -> sigma_ps equals resonant_source.gaussian_sigma (1e-12)",
          abs(pl.inputs["sigma_ps"] * 1e-3 - sig_ns) <= 1e-12 * sig_ns + 1e-15
          and abs(pl.inputs["sigma_ps"] - fw / (2 * math.sqrt(math.log(2)))) < 1e-12)
# intensity (Omega^2) FWHM: Omega^2(t) = exp(-t^2/sigma^2) halves at t = sigma sqrt(ln2) -> FWHM 2 sigma sqrt(ln 2)
s3 = entries["resonant_pulse"].plan({"fwhm_ps": 3.0}, EM, RT).inputs["sigma_ps"]
check("sigma consistent with the Omega^2 half-maximum points (exp(-t^2/sigma^2) = 1/2)",
      abs(math.exp(-((1.5 / s3) ** 2)) - 0.5) < 1e-12)

# ---------------------------------------------------------------- 5. tpe without XX
check("tpe.plan on a material without XX -> ValueError mentioning XX",
      raises(ValueError, entries["tpe"].plan, {}, EM_NO_XX, RT, contains=("XX",))[0])


@dataclass(frozen=True)
class _Mat(Material):
    def build(self, vals, T_K):
        raise NotImplementedError


@dataclass(frozen=True)
class _Str(Structure):
    def response(self, vals, em, T_K, NA):
        raise NotImplementedError


def _spec(drive_kind):
    return {"schema": "fsim-composed/1", "T_K": 4.2,
            "material": {"kind": "m_plain"}, "structure": {"kind": "s_plain"},
            "drive": {"kind": drive_kind, "params": {}}, "setup": {"kind": "lens_fibre_det"},
            "solver": "auto"}


vreg = Registry()
vreg.register(_Mat(kind="m_plain", label="m", summary="no XX", provides=frozenset({"X"})))
vreg.register(_Str(kind="s_plain", label="s", summary="s"))
vreg.register(R.LENS_FIBRE_DET)
vreg.register(entries["tpe"])
vreg.register(entries["resonant_pulse"])
vreg.register(_Mat(kind="m_xx", label="m", summary="with XX", provides=frozenset({"X", "XX"})))
res_tpe = validate_spec(_spec("tpe"), vreg)
check("validate_spec: tpe on a material that does not provide XX -> error issue naming XX",
      not res_tpe["ok"] and any(i["level"] == "error" and "XX" in i["message"] for i in res_tpe["issues"]))
spec_xx = _spec("tpe")
spec_xx["material"] = {"kind": "m_xx"}
check("validate_spec: tpe with an XX-providing material passes the requires check",
      not any("requires" in i["message"] for i in validate_spec(spec_xx, vreg)["issues"]))
vreg.register(THROWAWAY)
rp_ok = validate_spec(_spec("resonant_pulse"), vreg)
check("validate_spec: resonant_pulse resolves defaults cleanly; needs_physics drive is an error issue",
      not any(i["level"] == "error" and i["path"].startswith("drive") for i in rp_ok["issues"])
      and any("awaiting physics" in i["message"] for i in validate_spec(_spec("test_needs_physics"), vreg)["issues"]))

# ---------------------------------------------------------------- 6. electrical == direct drive_mech (bit-for-bit)
ep = entries["electrical_pulsed"]
ec = entries["electrical_cw"]
v = {"mechanism": "M-3 pulsed", "I_uA": 7.5, "tau_on_ns": 0.8, "eta_capture": 2e-5,
     "F_p_rail": 0.9, "R_drive_ohm": 50.0, "C_dep_F": 1e-12}
pl = ep.plan(v, EM, RT)
direct = drive_mech.mech_pulsed(I_uA=7.5, tau_on_ns=0.8, eta_capture=2e-5, F_p_rail=0.9,
                                R_drive_ohm=50.0, C_dep_F=1e-12)
check("electrical_pulsed M-3: interface == direct mech_pulsed (dataclass equality, incl. notes)",
      pl.interface == direct and pl.family == "electrical" and pl.timing == "pulsed"
      and pl.interface.g2_load.__func__ is drive_mech.DriveInterface.g2_load)
check("electrical_pulsed M-3: tau_on_ns = tau_eff after the RC rise time; mu = eta*N_c bit-identical",
      pl.tau_on_ns == direct.notes["tau_eff_ns"] and pl.interface.mu == direct.mu)
pl2 = ep.plan({"mechanism": "M-3 pulsed"}, EM, RT)
check("electrical_pulsed defaults: M-3 with I=10 uA, tau_on 1 ns, eta 1e-5 == direct",
      pl2.interface == drive_mech.mech_pulsed(I_uA=10.0, tau_on_ns=1.0, eta_capture=1e-5))
vq = {"mechanism": "M-2 quiet-rail", "mu": 0.7, "eta_capture": 0.6, "Z_drive_ohm": 3e3,
      "Z_junction_ohm": 7e3, "f_reg_Hz": 1e6, "f_sample_Hz": 1e8}
plq = ec.plan(vq, EM, RT)
dq = drive_mech.mech_quiet_rail(mu=0.7, Z_drive_ohm=3e3, Z_junction_ohm=7e3, eta_capture=0.6,
                                f_reg_Hz=1e6, f_sample_Hz=1e8)
check("electrical_cw M-2 quiet-rail (with regulation bandwidth): interface == direct mech_quiet_rail",
      plq.interface == dq and plq.timing == "cw" and plq.f_rep_MHz is None and plq.tau_on_ns == 0.0)
plr = ec.plan({"mechanism": "M-5a RTI", "mu": 0.4, "eta_capture": 0.05, "eta_boost": 3.0,
               "dg_inj_factor": 0.4, "F_p_rail": 1.0}, EM, RT)
check("electrical_cw M-5a RTI: interface == direct mech_rti",
      plr.interface == drive_mech.mech_rti(mu=0.4, eta_capture_base=0.05, eta_boost=3.0,
                                           dg_inj_factor=0.4, F_p_rail=1.0))
pl1 = ec.plan({"mechanism": "M-1 poisson-rail", "mu": 0.4, "eta_capture": 0.8}, EM, RT)
check("electrical_cw M-1: interface == direct mech_poisson_rail",
      pl1.interface == drive_mech.mech_poisson_rail(mu=0.4, eta_capture=0.8))
pls = ec.plan({"mechanism": "M-4 SET/gated", "T_junction_K": 77.0, "mu": 0.9, "eps_cycle": 0.02,
               "radius_nm": 20.0, "R_T_ohm": 1e6}, EM, RT)
check("electrical_cw M-4 SET/gated: interface == direct mech_set (feasibility notes bit-identical)",
      pls.interface == drive_mech.mech_set("gated", T_K=77.0, eps_cycle=0.02, mu=0.9,
                                           radius_nm=20.0, R_T_ohm=1e6))
plt = ep.plan({"mechanism": "M-4 SET/turnstile", "T_junction_K": 4.0, "R_T_ohm": 2e5}, EM, RT)
check("electrical_pulsed M-4 SET/turnstile (cycled -> plan timing pulsed): interface == direct",
      plt.interface == drive_mech.mech_set("turnstile", T_K=4.0, eps_cycle=0.01, mu=0.5,
                                           radius_nm=20.0, R_T_ohm=2e5) and plt.timing == "pulsed")
check("electrical: a cw mechanism in the pulsed entry (and vice versa) is rejected",
      raises(ValueError, ep.plan, {"mechanism": "M-1 poisson-rail"}, EM, RT)[0]
      and raises(ValueError, ec.plan, {"mechanism": "M-3 pulsed"}, EM, RT)[0])
check("electrical summaries carry the RT-edge passthrough label",
      all("use the RT-edge tier passthrough for the full transport chain" in e.summary for e in (ep, ec)))
check("every plan payload field in inputs is JSON-serialisable",
      all(json.dumps(p.inputs) for p in plans.values()))

n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} compose_drives checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
