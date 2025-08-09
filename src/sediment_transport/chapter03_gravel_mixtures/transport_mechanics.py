"""
Core Transport Mechanics for Mixed Sediments (Chapter 3, Equations 3-77 to 3-115)

This module implements the major bed-load transport relations for gravel and
sediment mixtures, including hiding functions, abrasion effects, and morphologic
complexity corrections.
"""

import numpy as np
from typing import Union, Tuple, Dict, Literal
from ..utils.validators import validate_positive, validate_array
from ..utils.constants import GRAVITY, WATER_DENSITY


def ashida_michiue_transport(tau_si_star: Union[float, np.ndarray],
                            tau_sci_star: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Ashida-Michiue bed-load transport relation.
    
    Equation 3-77a: q_i* = 17(τ_si* - τ_sci*)√(τ_si* - τ_sci*)
    
    Args:
        tau_si_star: Shields stress for size i (dimensionless)
        tau_sci_star: Critical Shields stress for size i (dimensionless)
        
    Returns:
        Einstein transport number q_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-77a
    """
    validate_positive(tau_si_star, "shields_stress")
    validate_positive(tau_sci_star, "critical_shields_stress")
    
    # Only transport when above threshold
    excess_stress = np.maximum(tau_si_star - tau_sci_star, 0.0)
    
    # Ashida-Michiue relation
    q_i_star = 17.0 * excess_stress * np.sqrt(excess_stress)
    
    return q_i_star


def parker_klingeman_mclean_transport(tau_si_star: Union[float, np.ndarray],
                                     tau_suri_star: Union[float, np.ndarray],
                                     phi_range: Literal['low', 'high'] = 'high') -> Union[float, np.ndarray]:
    """
    Parker-Klingeman-McLean bed-load transport relation.
    
    Equations 3-78c,d,e: W_ui* = G_u(φ), φ = τ_si*/τ_suri*
    
    Args:
        tau_si_star: Shields stress for size i (dimensionless)
        tau_suri_star: Reference Shields stress for size i (dimensionless)
        phi_range: Transport regime ('low' for φ < 1.65, 'high' for φ ≥ 1.65)
        
    Returns:
        Dimensionless transport rate W_ui* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-78c,d,e
    """
    validate_positive(tau_si_star, "shields_stress")
    validate_positive(tau_suri_star, "reference_shields_stress")
    
    # Mobility parameter
    phi = tau_si_star / tau_suri_star
    
    if phi_range == 'low':
        # Equation 3-78e: For 0.95 ≤ φ < 1.65
        phi_clipped = np.clip(phi, 0.95, 1.65)
        G_u = 0.0025 * np.exp(14.2 * (phi_clipped - 1.0) - 9.28 * (phi_clipped - 1.0)**2)
    else:
        # Equation 3-78e: For φ ≥ 1.65
        phi_clipped = np.maximum(phi, 1.65)
        G_u = 11.2 * (1.0 - 0.822 / phi_clipped)**4.5
    
    return G_u


def parker_surface_based_transport(phi: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Parker surface-based bed-load transport relation.
    
    Equations 3-79a,f: W_i* = 0.00218*G(φ) with transport function G(φ)
    
    Args:
        phi: Mobility parameter (dimensionless)
        
    Returns:
        Dimensionless transport rate W_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-79a,f
    """
    validate_positive(phi, "mobility_parameter")
    
    # Transport function G(φ)
    G_phi = np.where(
        phi > 1.59,
        5474.0 * (1.0 - 0.853 / phi)**4.5,
        np.where(
            phi >= 1.0,
            np.exp(14.2 * (phi - 1.0) - 9.28 * (phi - 1.0)**2),
            phi**14.2
        )
    )
    
    # Surface-based Parker relation
    W_i_star = 0.00218 * G_phi
    
    return W_i_star


def tsujimoto_entrainment_transport(tau_si_star: Union[float, np.ndarray],
                                   tau_sci_star: Union[float, np.ndarray],
                                   L_so_star: float = 100.0) -> Union[float, np.ndarray]:
    """
    Tsujimoto entrainment-based transport relation.
    
    Equations 3-80a,e: E_i* = 0.02τ_si*[1 - 0.7(τ_sci*/τ_si*)]³
                       q_i* = E_i*L_si* = 0.02L_so*τ_si*[1 - 0.7(τ_sci*/τ_si*)]³
    
    Args:
        tau_si_star: Shields stress for size i (dimensionless)
        tau_sci_star: Critical Shields stress for size i (dimensionless)
        L_so_star: Dimensionless reference step length, default 100.0
        
    Returns:
        Einstein transport number q_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-80a,e
    """
    validate_positive(tau_si_star, "shields_stress")
    validate_positive(tau_sci_star, "critical_shields_stress")
    validate_positive(L_so_star, "reference_step_length")
    
    # Avoid division by zero and ensure physical validity
    tau_ratio = np.minimum(tau_sci_star / tau_si_star, 1.0)
    
    # Entrainment rate
    E_i_star = 0.02 * tau_si_star * (1.0 - 0.7 * tau_ratio)**3
    
    # Transport rate (step length assumed constant)
    q_i_star = E_i_star * L_so_star
    
    return q_i_star


def hunziker_jaeggi_transport(D_i: Union[float, np.ndarray],
                             D_m: float,
                             tau_sm_star: float,
                             tau_scm_star: float = 0.05) -> Union[float, np.ndarray]:
    """
    Hunziker-Jaeggi bed-load transport relation.
    
    Equation 3-81a: q_i* = 5(D_i/D_m)^(-3/2) * α(D_i/D_m)^1.5 * (τ_sm* - τ_scm*)
    Equation 3-81e: α = 0.011(τ_sm*)^1.5 - 0.3
    
    Args:
        D_i: Grain size for size class i (mm)
        D_m: Arithmetic mean surface grain size (mm)
        tau_sm_star: Mean Shields stress (dimensionless)
        tau_scm_star: Critical mean Shields stress, default 0.05
        
    Returns:
        Einstein transport number q_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-81a,e
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_m, "mean_grain_size")
    validate_positive(tau_sm_star, "mean_shields_stress")
    
    # Only transport when above threshold
    if tau_sm_star <= tau_scm_star:
        return np.zeros_like(D_i)
    
    # Size ratio
    size_ratio = D_i / D_m
    
    # Alpha parameter
    alpha = 0.011 * (tau_sm_star**1.5) - 0.3
    
    # Hunziker-Jaeggi relation
    q_i_star = 5.0 * (size_ratio**(-1.5)) * alpha * (size_ratio**1.5) * (tau_sm_star - tau_scm_star)
    
    return q_i_star


def wilcock_kenworthy_two_fraction(phi: Union[float, np.ndarray],
                                  is_sand: bool = False,
                                  F1: float = 0.5) -> Union[float, np.ndarray]:
    """
    Wilcock-Kenworthy two-fraction model.
    
    Equation 3-83c: G = 0.002/(1-φ)^7.5 for φ < φ', G = A*φ^0.25 for φ ≥ φ'
    
    Args:
        phi: Mobility parameter φ = τ_si*/τ_ssri* (dimensionless)
        is_sand: If True, applies sand enhancement; if False, gravel fraction
        F1: Sand fraction in surface (for sand enhancement), default 0.5
        
    Returns:
        Transport function G (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-83c
    """
    validate_positive(phi, "mobility_parameter")
    
    # Critical mobility parameter (typically around 1.0)
    phi_prime = 1.0
    
    # Transport function
    if np.isscalar(phi):
        if phi < phi_prime:
            G = 0.002 / ((1.0 - phi)**7.5) if phi < 1.0 else np.inf
        else:
            A = 1.0  # Coefficient for high mobility regime
            G = A * (phi**0.25)
    else:
        G = np.where(
            phi < phi_prime,
            np.where(phi < 1.0, 0.002 / ((1.0 - phi)**7.5), 1e6),  # Large value for phi → 1
            phi**0.25
        )
    
    # Sand enhancement (Equation 3-83d)
    if is_sand:
        k = 20.0  # Sand enhancement parameter
        enhancement_factor = np.exp(-k * F1)
        G = G * enhancement_factor
    
    return G


def wilcock_crowe_multisize(phi: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Wilcock-Crowe multi-fraction transport model.
    
    Equation 3-84b: G = 0.002φ^7.5 for φ < 1.35, G = 14[1 - 0.894/φ^0.5]^4.5 for φ ≥ 1.35
    
    Args:
        phi: Mobility parameter (dimensionless)
        
    Returns:
        Transport function G (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-84b
    """
    validate_positive(phi, "mobility_parameter")
    
    # Critical mobility parameter
    phi_crit = 1.35
    
    # Transport function
    G = np.where(
        phi < phi_crit,
        0.002 * (phi**7.5),
        14.0 * (1.0 - 0.894 / np.sqrt(phi))**4.5
    )
    
    return G


def wu_wang_jia_transport(tau_si_star: Union[float, np.ndarray],
                         tau_suri_star: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Wu-Wang-Jia bed-load transport relation.
    
    Equation 3-85a: W_i* = 0.0053 * (1/τ_suri*)^(3/2) * (τ_si*/τ_suri* - 1)^2.2
    
    Args:
        tau_si_star: Shields stress for size i (dimensionless)
        tau_suri_star: Reference Shields stress for size i (dimensionless)
        
    Returns:
        Dimensionless transport rate W_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-85a
    """
    validate_positive(tau_si_star, "shields_stress")
    validate_positive(tau_suri_star, "reference_shields_stress")
    
    # Mobility parameter
    phi = tau_si_star / tau_suri_star
    
    # Only transport when above threshold
    excess_mobility = np.maximum(phi - 1.0, 0.0)
    
    # Wu-Wang-Jia relation
    W_i_star = 0.0053 * (tau_suri_star**(-1.5)) * (excess_mobility**2.2)
    
    return W_i_star


def powell_reid_laronne_transport(tau_si_star: Union[float, np.ndarray],
                                 tau_sci_star: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Powell-Reid-Laronne bed-load transport relation.
    
    Equation 3-86a: W_i* = 11.2 * (1 - 1/φ)^4.5, φ = τ_si*/τ_sci*
    
    Args:
        tau_si_star: Shields stress for size i (dimensionless)
        tau_sci_star: Critical Shields stress for size i (dimensionless)
        
    Returns:
        Dimensionless transport rate W_i* (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-86a
    """
    validate_positive(tau_si_star, "shields_stress")
    validate_positive(tau_sci_star, "critical_shields_stress")
    
    # Mobility parameter
    phi = tau_si_star / tau_sci_star
    
    # Only transport when above threshold
    phi_clipped = np.maximum(phi, 1.0)
    
    # Powell-Reid-Laronne relation
    W_i_star = 11.2 * (1.0 - 1.0 / phi_clipped)**4.5
    
    return W_i_star


def generic_bedload_transport(tau_bs: Union[float, np.ndarray],
                             D: Union[float, np.ndarray],
                             tau_sc_star: Union[float, np.ndarray],
                             n_L: float = 1.5,
                             rho_s: float = 2650.0,
                             rho: float = WATER_DENSITY) -> Union[float, np.ndarray]:
    """
    Generic bed-load transport relation.
    
    Equation 3-89a: q = √(RgDD) * [(τ_bs/ρRgD) - τ*_sc]^n_L
    
    Args:
        tau_bs: Boundary shear stress (Pa)
        D: Grain size (m)
        tau_sc_star: Critical Shields stress (dimensionless)
        n_L: Transport exponent, default 1.5
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Bed-load transport rate (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-89a
    """
    validate_positive(tau_bs, "boundary_shear_stress")
    validate_positive(D, "grain_size")
    validate_positive(tau_sc_star, "critical_shields_stress")
    
    # Submerged specific gravity
    R = (rho_s / rho) - 1.0
    
    # Shields stress
    tau_star = tau_bs / (rho * R * GRAVITY * D)
    
    # Excess stress
    excess_stress = np.maximum(tau_star - tau_sc_star, 0.0)
    
    # Generic transport relation
    q = np.sqrt(R * GRAVITY * D**3) * (excess_stress**n_L)
    
    return q


def complexity_corrected_transport(q_base: Union[float, np.ndarray],
                                  C_comp: Union[float, np.ndarray] = 1.0) -> Union[float, np.ndarray]:
    """
    Apply morphologic complexity correction to transport rate.
    
    Equation 3-89c: q̄ = C_comp * q(τ̄_bs, D̄)
    
    Args:
        q_base: Base transport rate without complexity correction (m²/s)
        C_comp: Complexity coefficient (≥ 1), default 1.0
        
    Returns:
        Complexity-corrected transport rate (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-89c
    """
    validate_positive(C_comp, "complexity_coefficient")
    
    if np.any(C_comp < 1.0):
        raise ValueError("Complexity coefficient must be ≥ 1.0")
    
    return C_comp * q_base


def equal_mobility_condition(F_i: np.ndarray, f_bi: np.ndarray) -> bool:
    """
    Check if equal mobility condition is satisfied.
    
    Equation 3-91: f_bi = F_i (perfect equal mobility)
    
    Args:
        F_i: Surface size fractions (dimensionless)
        f_bi: Bed-load size fractions (dimensionless)
        
    Returns:
        True if equal mobility condition is approximately satisfied
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-91
    """
    validate_array(F_i)
    validate_array(f_bi)
    
    if len(F_i) != len(f_bi):
        raise ValueError("Surface and bed-load fraction arrays must have same length")
    
    # Check if distributions are approximately equal (within 10% tolerance)
    relative_diff = np.abs(f_bi - F_i) / (F_i + 1e-12)  # Avoid division by zero
    
    return np.all(relative_diff < 0.1)


def static_armor_predictor(f_i: np.ndarray, D_i: np.ndarray) -> np.ndarray:
    """
    Predict static armor size distribution.
    
    Equation 3-111: F_ai = (f_i * D_i^1.35) / Σ(f_i * D_i^1.35)
    
    Args:
        f_i: Substrate size fractions (dimensionless)
        D_i: Grain sizes for each fraction (mm)
        
    Returns:
        Static armor size fractions (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-111
    """
    validate_array(f_i)
    validate_array(D_i)
    
    if len(f_i) != len(D_i):
        raise ValueError("Substrate fractions and grain sizes must have same length")
    
    # Weighted by grain size to power 1.35
    weighted_fractions = f_i * (D_i**1.35)
    
    # Normalize to get static armor distribution
    F_ai = weighted_fractions / np.sum(weighted_fractions)
    
    return F_ai


def shields_stress_calculation(H: Union[float, np.ndarray], 
                              S: Union[float, np.ndarray],
                              D50: Union[float, np.ndarray],
                              rho_s: float = 2650.0,
                              rho: float = WATER_DENSITY) -> Union[float, np.ndarray]:
    """
    Calculate Shields stress from flow conditions.
    
    Equation 3-113: τ*50 = HS/(RD50)
    
    Args:
        H: Flow depth (m)
        S: Bed slope (dimensionless)
        D50: Median grain size (m)
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Shields stress τ*50 (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-113
    """
    validate_positive(H, "flow_depth")
    validate_positive(S, "bed_slope")
    validate_positive(D50, "median_grain_size")
    
    # Submerged specific gravity
    R = (rho_s / rho) - 1.0
    
    # Shields stress
    tau_star_50 = (H * S) / (R * D50)
    
    return tau_star_50