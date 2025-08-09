"""
Active Layer Sediment Transport Mechanics (Chapter 3, Equations 3-23 to 3-44)

This module implements the active layer formulation for sediment transport
in gravel-bed rivers, including conservation equations, entrainment-deposition
processes, and dimensionless transport parameters.
"""

import numpy as np
from typing import Union, Tuple, Optional
from ..utils.validators import validate_positive, validate_array
from ..utils.constants import GRAVITY, WATER_DENSITY


def active_layer_conservation(F_i: np.ndarray, 
                             eta_b: float,
                             L_a: float,
                             q_i: np.ndarray,
                             f_Ii: np.ndarray,
                             lambda_p: float = 0.4,
                             dx: float = 1.0,
                             dt: float = 1.0) -> np.ndarray:
    """
    Active layer conservation equation for sediment transport.
    
    Equation 3-23: (1-λ_p)[f_Ii ∂η_b/∂t + ∂/∂t(L_a F_i)] = -∂q_i/∂s
    
    Args:
        F_i: Size fractions in active layer (array)
        eta_b: Elevation of bottom of surface layer (m)
        L_a: Active layer thickness (m)
        q_i: Bed-load transport rate by size class (m²/s)
        f_Ii: Interfacial exchange fractions (array)
        lambda_p: Porosity (dimensionless), default 0.4
        dx: Streamwise grid spacing (m), default 1.0
        dt: Time step (s), default 1.0
        
    Returns:
        Rate of change of active layer composition
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-23
    """
    validate_array(F_i)
    validate_array(q_i) 
    validate_array(f_Ii)
    validate_positive(L_a, "active_layer_thickness")
    
    if len(F_i) != len(q_i) or len(F_i) != len(f_Ii):
        raise ValueError("All size fraction arrays must have same length")
    
    # Spatial gradient of transport rate (∂q_i/∂s)
    dq_dx = np.gradient(q_i, dx)
    
    # Right-hand side: -∂q_i/∂s
    rhs = -dq_dx
    
    # Left-hand side coefficient
    coeff = 1.0 - lambda_p
    
    # Rate of change term
    rate_of_change = rhs / coeff
    
    return rate_of_change


def bed_elevation_change(eta_b: float, L_a: float) -> float:
    """
    Relationship between bed elevation and active layer.
    
    Equation 3-24: η = η_b + L_a
    
    Args:
        eta_b: Bottom of surface layer elevation (m)
        L_a: Active layer thickness (m)
        
    Returns:
        Bed elevation (m)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-24
    """
    return eta_b + L_a


def total_transport_conservation(q_T: float, 
                                lambda_p: float = 0.4,
                                dx: float = 1.0) -> float:
    """
    Total transport mass conservation.
    
    Equation 3-25: (1-λ_p) ∂η/∂t = -∂q_T/∂s
    
    Args:
        q_T: Total bed-load transport rate (m²/s)
        lambda_p: Porosity (dimensionless), default 0.4
        dx: Streamwise grid spacing (m), default 1.0
        
    Returns:
        Rate of bed elevation change (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-25
    """
    validate_positive(q_T, "total_transport_rate")
    
    # Spatial gradient of total transport
    dq_T_dx = np.gradient([q_T], dx)[0] if hasattr(q_T, '__len__') else 0
    
    # Rate of bed elevation change
    deta_dt = -dq_T_dx / (1.0 - lambda_p)
    
    return deta_dt


def total_transport_rate(q_i: np.ndarray) -> float:
    """
    Calculate total transport rate from size class contributions.
    
    Equation 3-26: q_T = Σq_i (i=1 to n)
    
    Args:
        q_i: Bed-load transport rates by size class (m²/s)
        
    Returns:
        Total bed-load transport rate (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-26
    """
    validate_array(q_i)
    return np.sum(q_i)


def bed_load_fractions(q_i: np.ndarray) -> np.ndarray:
    """
    Calculate bed load size fractions from transport rates.
    
    Equation 3-28: f_bi = q_i / Σq_i
    
    Args:
        q_i: Bed-load transport rates by size class (m²/s)
        
    Returns:
        Bed load size fractions (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-28
    """
    validate_array(q_i)
    
    q_total = np.sum(q_i)
    if q_total <= 0:
        raise ValueError("Total transport rate must be positive")
        
    return q_i / q_total


def active_layer_thickness(D90: float, n_a: float = 2.0) -> float:
    """
    Calculate active layer thickness from grain size.
    
    Equation 3-30: L_a = n_a D90
    
    Args:
        D90: 90th percentile grain size (m)
        n_a: Calibration parameter (order of 1), default 2.0
        
    Returns:
        Active layer thickness (m)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-30
    """
    validate_positive(D90, "D90_grain_size")
    validate_positive(n_a, "calibration_parameter")
    
    return n_a * D90


def exchange_fractions_degradation(f_i_substrate: np.ndarray,
                                  eta_b: float) -> np.ndarray:
    """
    Calculate exchange fractions during bed degradation.
    
    Equation 3-31: f_Ii = f_i(z)|_z=η_b for ∂η_b/∂t < 0
    
    Args:
        f_i_substrate: Substrate size fractions at depth (array)
        eta_b: Bottom of surface layer elevation (m)
        
    Returns:
        Exchange fractions for degradation
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-31
    """
    validate_array(f_i_substrate)
    
    # During degradation, exchange fractions equal substrate fractions
    return f_i_substrate.copy()


def entrainment_deposition_rates(E_bi: np.ndarray, 
                                P_si: np.ndarray,
                                step_lengths: np.ndarray,
                                dx: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate entrainment and deposition rates for bed load transport.
    
    Equations 3-34, 3-35: Mean step length and deposition rate
    Equations 3-36, 3-37: Entrainment formulation and transport rate
    
    Args:
        E_bi: Bed-load entrainment rates by size (1/s)
        P_si: Probability density for step lengths (1/m)
        step_lengths: Array of step length values (m)
        dx: Spatial grid spacing (m), default 1.0
        
    Returns:
        Tuple of (deposition_rates, transport_rates)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-34 to 3-37
    """
    validate_array(E_bi)
    validate_array(P_si)
    validate_array(step_lengths)
    
    # Equation 3-34: Mean step length
    L_si = np.sum(step_lengths * P_si * dx)
    
    # Equation 3-35: Deposition rate (simplified)
    D_bi = E_bi * np.exp(-step_lengths / L_si)
    
    # Equation 3-37: Transport rate
    q_i = E_bi * L_si
    
    return D_bi, q_i


def dimensionless_transport_parameters(q_i: np.ndarray,
                                     F_i: np.ndarray, 
                                     D_i: np.ndarray,
                                     rho_s: float = 2650.0,
                                     rho: float = WATER_DENSITY) -> Tuple[np.ndarray, float]:
    """
    Calculate dimensionless transport parameters.
    
    Equations 3-44a,b,c: Submerged specific gravity, Einstein number, particle Reynolds
    
    Args:
        q_i: Transport rates by size class (m²/s)
        F_i: Surface size fractions (dimensionless)
        D_i: Grain sizes for each class (m)
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Tuple of (Einstein_numbers, particle_Reynolds_number)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-44a,b,c
    """
    validate_array(q_i)
    validate_array(F_i)
    validate_array(D_i)
    
    if len(q_i) != len(F_i) or len(q_i) != len(D_i):
        raise ValueError("All arrays must have same length")
    
    # Equation 3-44a: Submerged specific gravity
    R = (rho_s / rho) - 1.0
    
    # Equation 3-44b: Grain-size-specific Einstein number
    q_i_star = np.zeros_like(q_i)
    for i in range(len(q_i)):
        if F_i[i] > 0:
            denominator = F_i[i] * np.sqrt(R * GRAVITY * D_i[i]**3)
            q_i_star[i] = q_i[i] / denominator
        
    # Equation 3-44c: Particle Reynolds number (using geometric mean)
    D_g = np.exp(np.sum(F_i * np.log(D_i + 1e-12)))  # Avoid log(0)
    R_pg = np.sqrt(R * GRAVITY * D_g**3) / 1e-6  # Assuming nu = 1e-6 m²/s
    
    return q_i_star, R_pg


def surface_based_transport_formulation(q_Ui: np.ndarray, 
                                       F_i: np.ndarray) -> np.ndarray:
    """
    Surface-based formulation for bed-load transport.
    
    Equation 3-42a: q_i = q_Ui / F_i
    
    Args:
        q_Ui: Volume transport rate per unit fraction content (m²/s)
        F_i: Surface size fractions (dimensionless)
        
    Returns:
        Transport rates by size class (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-42a
    """
    validate_array(q_Ui)
    validate_array(F_i)
    
    if len(q_Ui) != len(F_i):
        raise ValueError("Arrays must have same length")
    
    # Avoid division by zero
    q_i = np.zeros_like(q_Ui)
    for i in range(len(q_Ui)):
        if F_i[i] > 1e-12:  # Small threshold to avoid division by zero
            q_i[i] = q_Ui[i] / F_i[i]
            
    return q_i


def entrainment_rate_formulation(E_Ubi: np.ndarray, 
                                F_i: np.ndarray) -> np.ndarray:
    """
    Surface-based entrainment rate formulation.
    
    Equation 3-42b: E_bi = E_Ubi / F_i
    
    Args:
        E_Ubi: Entrainment rate per unit fraction content (1/s)
        F_i: Surface size fractions (dimensionless)
        
    Returns:
        Entrainment rates by size class (1/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-42b
    """
    validate_array(E_Ubi)
    validate_array(F_i)
    
    if len(E_Ubi) != len(F_i):
        raise ValueError("Arrays must have same length")
    
    # Avoid division by zero
    E_bi = np.zeros_like(E_Ubi)
    for i in range(len(E_Ubi)):
        if F_i[i] > 1e-12:
            E_bi[i] = E_Ubi[i] / F_i[i]
            
    return E_bi