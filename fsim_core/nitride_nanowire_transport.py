"""Single-wire axial p-i-n transport and pulse diagnostics for InGaN/GaN.

This is deliberately an opt-in compact model.  It reuses the planar
NitrideDiode voltage/depletion kernel but never calls its aperture-based
loading path.  Every circuit electron supplies one electron-hole pair
equivalent (I/q), rather than two half-pairs.  Deshpande et al., Nat.
Commun. 4, 1675 (2013), reported 2.38 Gohm for their contacted wire [V];
the remaining transport and thermal defaults below are labelled assumptions.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math

from .nitride_materials import KB_EV, EPS0_SI
from .nitride_transport import NitrideDiode
from .transport import Q_SI, qfl_suppression, xi_window
from .thermal import Layer, Stack, t_junction


_PROVENANCE = {
    "R_s_ohm": "[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.4 measured device series resistance",
    "designed_R_s_ohm": "[A] 1e6 ohm designed-contact sensitivity; not a resistivity-derived geometry law",
    "geometry": "[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.2: 2 nm disc and 15 nm GaN spacers; [A] transfer defaults",
    "capture": "[A] capture/matrix competition envelope; tau_cap is an explicit input, not a fitted Deshpande lifetime",
    "thermal": "[A] constant-Rth adapter to thermal.t_junction: 1e9 K/W horizontal, 1e7 K/W vertical design baseline",
    "eta_total": "[A] 0.01 emitted-optical-power approximation; it is not collection efficiency",
    "pulse": "[A] 90 percent step and 10 percent inter-pulse residual diagnostic thresholds",
}


def _finite_positive(name, value, zero=False):
    value = float(value)
    if not math.isfinite(value) or (value < 0.0 if zero else value <= 0.0):
        raise ValueError(name + " must be finite and " + ("nonnegative" if zero else "positive"))
    return value


@dataclass(frozen=True)
class NitrideWireDiode:
    """Axial wire geometry. area_cm2 is derived from conducting radius only."""
    N_A: float = 1e17
    N_D: float = 1e18
    n_ideality: float = 2.0
    tau_SRH_ns: float = 1.0
    eps_r: float = 10.28
    T: float = 300.0
    core_radius_nm: float = 12.5
    conducting_radius_nm: float = 12.5
    barrier_left_nm: float = 15.0
    barrier_right_nm: float = 15.0
    d_active_nm: float = 2.0
    x_in: float = 0.40
    R_s_ohm: float = 2.38e9
    f_Rs_local: float = 1.0
    tau_matrix_ns: float = 1.0
    tau_cap_ps: float = 10.0
    reservoir_surface_ns: float = 0.0

    def __post_init__(self):
        for name in ("N_A", "N_D", "n_ideality", "tau_SRH_ns", "eps_r", "T",
                     "core_radius_nm", "conducting_radius_nm", "barrier_left_nm",
                     "barrier_right_nm", "d_active_nm", "R_s_ohm", "tau_matrix_ns",
                     "tau_cap_ps"):
            _finite_positive(name, getattr(self, name))
        _finite_positive("reservoir_surface_ns", self.reservoir_surface_ns, zero=True)
        if not 0.0 <= self.x_in <= 1.0 or not 0.0 <= self.f_Rs_local <= 1.0:
            raise ValueError("x_in or f_Rs_local outside physical range")
        if self.conducting_radius_nm > self.core_radius_nm:
            raise ValueError("conducting_radius_nm must not exceed core_radius_nm")

    @property
    def d_i_nm(self):
        """[DR] Active disc plus both axial GaN spacers."""
        return self.barrier_left_nm + self.d_active_nm + self.barrier_right_nm

    @property
    def area_cm2(self):
        return math.pi * (self.conducting_radius_nm * 1e-7) ** 2

    def _kernel(self):
        return NitrideDiode(N_A=self.N_A, N_D=self.N_D, d_i_nm=self.d_i_nm,
                            d_active_nm=self.d_active_nm,
                            area_um2=self.area_cm2 * 1e8, R_s_ohm=self.R_s_ohm,
                            n_ideality=self.n_ideality, tau_SRH_ns=self.tau_SRH_ns,
                            eps_r=self.eps_r, T=self.T, x_in=self.x_in,
                            wl_thickness_nm=self.d_active_nm,
                            f_Rs_local=self.f_Rs_local)


def wire_pin(**overrides):
    """Construct the explicit wire preset; no planar-area override exists."""
    if "area_cm2" in overrides or "area_um2" in overrides:
        raise TypeError("wire area is derived from conducting_radius_nm")
    return NitrideWireDiode(**overrides)


def _injection_efficiency(barrier_e_eV, barrier_h_eV, T_K):
    be = _finite_positive("barrier_e_eV", barrier_e_eV, zero=True)
    bh = _finite_positive("barrier_h_eV", barrier_h_eV, zero=True)
    kT = KB_EV * T_K
    re, rh = math.exp(-be / kT), math.exp(-bh / kT)
    return 1.0 / (1.0 + re + rh), re, rh


def evaluate_injection(diode, *, I_uA, T_K, tau_pulse_ns, E_X_eV,
                       reservoir_energy_eV, barrier_e_eV, barrier_h_eV,
                       surface_reservoir_ns, tau_cap_ps, S_dot, w_meV,
                       eta_rad_matrix, eta_total):
    """Resolve serial supply, capture competition, background and accounting.

    ``surface_reservoir_ns`` is a loss *rate* in ns^-1 despite its historic
    suffix; the name is retained as frozen by the transport interface.
    """
    if not isinstance(diode, NitrideWireDiode):
        raise TypeError("diode must be NitrideWireDiode")
    I_uA = _finite_positive("I_uA", I_uA, zero=True)
    T_K = _finite_positive("T_K", T_K)
    tau_pulse_ns = _finite_positive("tau_pulse_ns", tau_pulse_ns, zero=True)
    for name, value in (("E_X_eV", E_X_eV), ("reservoir_energy_eV", reservoir_energy_eV),
                        ("w_meV", w_meV), ("eta_total", eta_total)):
        _finite_positive(name, value, zero=True)
    for name, value in (("surface_reservoir_ns", surface_reservoir_ns),
                        ("tau_cap_ps", tau_cap_ps), ("S_dot", S_dot),
                        ("eta_rad_matrix", eta_rad_matrix)):
        _finite_positive(name, value, zero=True)
    if eta_rad_matrix > 1.0 or eta_total > 1.0 or S_dot > 1.0:
        raise ValueError("efficiencies must be in [0, 1]")
    kernel = diode._kernel()
    I_A = I_uA * 1e-6
    V_terminal, V_j = kernel.v_of_i(I_A, T_K)
    dep = kernel.depletion(V_j, T_K)
    eta_inj, leak_e, leak_h = _injection_efficiency(barrier_e_eV, barrier_h_eV, T_K)
    r_supply = I_A / Q_SI
    r_leakage = r_supply * (1.0 - eta_inj)
    usable = r_supply - r_leakage
    f_dot = qfl_suppression(float(E_X_eV), V_j, KB_EV * T_K)
    f_bg = qfl_suppression(float(reservoir_energy_eV), V_j, KB_EV * T_K)
    k_cap = 1000.0 / float(tau_cap_ps) if tau_cap_ps > 0.0 else 0.0
    k_matrix_rad = eta_rad_matrix / diode.tau_matrix_ns
    k_matrix_nonrad = (1.0 - eta_rad_matrix) / diode.tau_matrix_ns
    k_surface = float(surface_reservoir_ns)
    k_total = k_cap + k_matrix_rad + k_matrix_nonrad + k_surface
    f_capture = k_cap / k_total if k_total else 0.0
    # Dot QFL suppression declines capture and returns that supply to the
    # reservoir, where it competes only with the declared reservoir channels.
    r_captured = usable * f_capture * S_dot * f_dot
    reservoir = usable - r_captured
    # The reservoir QFL factor is evaluated at its own energy and enters the
    # reservoir competition.  Suppressed radiative supply is thereby
    # redistributed to nonradiative/surface reservoir channels, never lost.
    k_matrix_rad_qfl = k_matrix_rad * f_bg
    kres = k_matrix_rad_qfl + k_matrix_nonrad + k_surface
    r_matrix_radiative = reservoir * k_matrix_rad_qfl / kres if kres else 0.0
    r_matrix_nonradiative = reservoir * k_matrix_nonrad / kres if kres else 0.0
    r_surface = reservoir * k_surface / kres if kres else reservoir
    r_other = 0.0
    xi = float(xi_window(float(w_meV), max((reservoir_energy_eV - E_X_eV) * 1e3, 0.0), KB_EV * T_K * 1e3))
    accepted_background = r_matrix_radiative * xi
    residual = r_supply - (r_leakage + r_captured + r_matrix_radiative +
                           r_matrix_nonradiative + r_surface + r_other)
    mu = r_captured * tau_pulse_ns * 1e-9
    hnu = float(E_X_eV)
    power = I_A * V_j + diode.f_Rs_local * I_A * I_A * diode.R_s_ohm - I_A * eta_total * hnu
    return {"valid": math.isfinite(residual), "reasons": [] if math.isfinite(residual) else ["nonfinite accounting"],
            "provenance": dict(_PROVENANCE), "area_cm2": diode.area_cm2, "J_A_cm2": I_A / diode.area_cm2,
            "V_j": V_j, "V_terminal": V_terminal, "depletion_field_kVcm": dep.F_kVcm,
            "C_dep_F": dep.C_dep_pF * 1e-12, "eta_inj": eta_inj, "leakage_e_ratio": leak_e,
            "leakage_h_ratio": leak_h, "f_capture": f_capture, "f_qfl_dot": f_dot,
            "f_qfl_background": f_bg, "r_supply_s": r_supply, "r_captured_s": r_captured,
            "r_matrix_radiative_s": r_matrix_radiative, "r_matrix_nonradiative_s": r_matrix_nonradiative,
            "r_surface_reservoir_s": r_surface, "r_leakage_s": r_leakage, "r_other_declared_loss_s": r_other,
            "raw_background_radiative_s": r_matrix_radiative, "accepted_background_s": accepted_background,
            "background_window_fraction": xi, "mu": mu, "power_on_W": power,
            "accounting_residual_s": residual}


def wire_operating_point(diode, *, I_uA, T_hs_K, duty, Rth_K_W, eta_total, h_nu_eV):
    """Self-consistent constant-Rth pulse-average heating [A] via thermal.py."""
    if not isinstance(diode, NitrideWireDiode):
        raise TypeError("diode must be NitrideWireDiode")
    I_uA = _finite_positive("I_uA", I_uA, zero=True)
    T_hs_K = _finite_positive("T_hs_K", T_hs_K)
    Rth_K_W = _finite_positive("Rth_K_W", Rth_K_W, zero=True)
    if not 0.0 <= float(duty) <= 1.0 or not 0.0 <= float(eta_total) <= 1.0:
        raise ValueError("duty and eta_total must be in [0, 1]")
    _finite_positive("h_nu_eV", h_nu_eV, zero=True)
    if I_uA == 0.0 or Rth_K_W == 0.0 or duty == 0.0:
        return {"valid": True, "reasons": [], "T_j_K": T_hs_K, "P_on_W": 0.0, "P_average_W": 0.0,
                "V_j": 0.0, "V_terminal": 0.0, "iterations": 0, "provenance": dict(_PROVENANCE)}
    radius_m = diode.conducting_radius_nm * 1e-9
    # A zero-temperature-coefficient artificial layer implements exactly a
    # constant Rth while still exercising thermal.t_junction's convergence.
    thickness = Rth_K_W * math.pi * radius_m * radius_m
    stack = Stack([Layer(thickness, 1.0, 0.0, False)], 1e300, 0.0)
    T_j = T_hs_K
    for iteration in range(1, 81):
        I_A = I_uA * 1e-6
        V_terminal, V_j = diode._kernel().v_of_i(I_A, T_j)
        P_on = I_A * V_j + diode.f_Rs_local * I_A * I_A * diode.R_s_ohm - I_A * eta_total * h_nu_eV
        P_average = duty * P_on
        T_new = t_junction(P_average, radius_m, stack, T_hs_K)
        if not math.isfinite(T_new):
            return {"valid": False, "reasons": ["thermal runaway/nonconvergence"], "T_j_K": float("inf"),
                    "P_on_W": P_on, "P_average_W": P_average, "V_j": V_j, "V_terminal": V_terminal,
                    "iterations": iteration, "provenance": dict(_PROVENANCE)}
        if abs(T_new - T_j) < 1e-5:
            return {"valid": True, "reasons": [], "T_j_K": T_new, "P_on_W": P_on, "P_average_W": P_average,
                    "V_j": V_j, "V_terminal": V_terminal, "iterations": iteration, "provenance": dict(_PROVENANCE)}
        T_j = 0.5 * T_j + 0.5 * T_new
    return {"valid": False, "reasons": ["thermal nonconvergence"], "T_j_K": T_j,
            "P_on_W": P_on, "P_average_W": P_average, "V_j": V_j, "V_terminal": V_terminal,
            "iterations": 80, "provenance": dict(_PROVENANCE)}


def pulse_delivery(diode, *, V_j, T_K, tau_pulse_ns, rep_rate_hz, C_parasitic_F=0.0):
    """RC diagnostic only: it never changes ideal pulse loading."""
    if not isinstance(diode, NitrideWireDiode):
        raise TypeError("diode must be NitrideWireDiode")
    V_j = _finite_positive("V_j", V_j, zero=True); T_K = _finite_positive("T_K", T_K)
    tau_pulse_ns = _finite_positive("tau_pulse_ns", tau_pulse_ns, zero=True)
    rep_rate_hz = _finite_positive("rep_rate_hz", rep_rate_hz)
    C_parasitic_F = _finite_positive("C_parasitic_F", C_parasitic_F, zero=True)
    dep = diode._kernel().depletion(V_j, T_K)
    C_total = dep.C_dep_pF * 1e-12 + C_parasitic_F
    tau_rc = diode.R_s_ohm * C_total
    if tau_rc == 0.0:
        step, f3, residual = 1.0, float("inf"), 0.0
    else:
        step = -math.expm1(-(tau_pulse_ns * 1e-9) / tau_rc)
        f3 = 1.0 / (2.0 * math.pi * tau_rc)
        residual = math.exp(-(1.0 / rep_rate_hz - tau_pulse_ns * 1e-9) / tau_rc) if tau_pulse_ns * 1e-9 < 1.0 / rep_rate_hz else 1.0
    feasible = step >= 0.9 and residual <= 0.1
    return {"valid": True, "reasons": [], "provenance": dict(_PROVENANCE), "C_dep_F": dep.C_dep_pF * 1e-12,
            "C_total_F": C_total, "tau_RC_s": tau_rc, "f_3dB_Hz": f3,
            "delivered_step_fraction": step, "inter_pulse_residual_fraction": residual,
            "pulse_delivery_feasible": feasible}
