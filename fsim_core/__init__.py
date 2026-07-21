"""fsim-core: headless F-series physics. No GUI imports, ever (three-layer rule)."""
from .card import Card, DataSet, Param, Tag, load_card, widest
from .fitting import FitResult, fit_phase0, V_A_TOL
from .integrator import g2_curve, g2_from, g2_of_T, oat_sensitivity, rho_of_T, solve_Tc
from .spectral import (
    KB,
    epsilon,
    epsilon_narrow_filter,
    gamma_of_T,
    mc_transmission,
    optimal_width,
    transmission,
)

__version__ = "0.1.0"
