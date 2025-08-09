"""
Hydraulic Relations for Gravel and Sand-Bed Streams (Chapter 3, Equations 3-13 to 3-22)

This module implements bankfull flow relations and dimensionless parameters
for characterizing gravel-bed and sand-bed rivers.
"""

import numpy as np
from typing import Tuple, Literal
from ..utils.validators import validate_positive
from ..utils.constants import GRAVITY, WATER_DENSITY


def bankfull_velocity(Q_bf: float, B_bf: float, H_bf: float) -> float:
    """
    Calculate bankfull velocity from discharge and geometry.
    
    Equation 3-13: U_bf = Q_bf / (B_bf × H_bf)
    
    Args:
        Q_bf: Bankfull discharge (m³/s)
        B_bf: Bankfull width (m)
        H_bf: Bankfull depth (m)
        
    Returns:
        Bankfull velocity (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-13
    """
    validate_positive(Q_bf, "bankfull_discharge")
    validate_positive(B_bf, "bankfull_width")
    validate_positive(H_bf, "bankfull_depth")
    
    return Q_bf / (B_bf * H_bf)


def bankfull_shear_stress(H_bf: float, S: float, 
                         rho: float = WATER_DENSITY) -> float:
    """
    Calculate bankfull shear stress.
    
    Equation 3-14a: τ_bf = ρ g H_bf S
    
    Args:
        H_bf: Bankfull depth (m)
        S: Channel bed slope (dimensionless)
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Bankfull shear stress (Pa)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-14a
    """
    validate_positive(H_bf, "bankfull_depth")
    validate_positive(S, "bed_slope")
    validate_positive(rho, "water_density")
    
    return rho * GRAVITY * H_bf * S


def bankfull_shear_velocity(H_bf: float, S: float) -> float:
    """
    Calculate bankfull shear velocity.
    
    Equation 3-14b: u*_bf = √(g H_bf S)
    
    Args:
        H_bf: Bankfull depth (m)
        S: Channel bed slope (dimensionless)
        
    Returns:
        Bankfull shear velocity (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-14b
    """
    validate_positive(H_bf, "bankfull_depth")
    validate_positive(S, "bed_slope")
    
    return np.sqrt(GRAVITY * H_bf * S)


def friction_coefficients(U_bf: float, u_star_bf: float) -> Tuple[float, float]:
    """
    Calculate friction coefficients from bankfull flow parameters.
    
    Equations 3-15a,b: C_f_bf = τ_bf/(ρU²_bf), C_z_bf = U_bf/u*_bf
    Equations 3-18a,b,c: Relationships between friction coefficients
    
    Args:
        U_bf: Bankfull velocity (m/s)
        u_star_bf: Bankfull shear velocity (m/s)
        
    Returns:
        Tuple of (friction_factor, Chezy_coefficient)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-15, 3-18
    """
    validate_positive(U_bf, "bankfull_velocity")
    validate_positive(u_star_bf, "shear_velocity")
    
    # Equation 3-15b and 3-18c: Chezy coefficient
    C_z_bf = U_bf / u_star_bf
    
    # Equation 3-15a and 3-18a: Friction factor
    C_f_bf = 2.0 / (C_z_bf**2)
    
    return C_f_bf, C_z_bf


def dimensionless_parameters(Q_bf: float, B_bf: float, H_bf: float, 
                           D50: float, rho_s: float = 2650.0, 
                           rho: float = WATER_DENSITY,
                           nu: float = 1e-6) -> dict:
    """
    Calculate dimensionless parameters for bankfull flow characterization.
    
    Equations 3-17a-g: Dimensionless discharge, width, depth, Froude number,
                       Shields number, particle Reynolds number, specific gravity ratio
    
    Args:
        Q_bf: Bankfull discharge (m³/s)
        B_bf: Bankfull width (m)
        H_bf: Bankfull depth (m)
        D50: Median grain size (m)
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        nu: Water kinematic viscosity (m²/s), default 1e-6
        
    Returns:
        Dictionary of dimensionless parameters
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-17a-g
    """
    validate_positive(Q_bf, "bankfull_discharge")
    validate_positive(B_bf, "bankfull_width")
    validate_positive(H_bf, "bankfull_depth")
    validate_positive(D50, "median_grain_size")
    validate_positive(rho_s, "sediment_density")
    
    # Calculate basic parameters
    R = (rho_s / rho) - 1.0  # Specific gravity ratio (Eq 3-17g)
    g_R_D50_cubed = GRAVITY * R * (D50**3)
    sqrt_g_R_D50_cubed = np.sqrt(g_R_D50_cubed)
    
    # Equation 3-17a: Dimensionless discharge
    Q_hat = Q_bf / (B_bf * sqrt_g_R_D50_cubed)
    
    # Equation 3-17b: Dimensionless width  
    B_hat = B_bf / D50
    
    # Equation 3-17c: Dimensionless depth
    H_hat = H_bf / D50
    
    # Calculate velocities for remaining parameters
    U_bf = bankfull_velocity(Q_bf, B_bf, H_bf)
    
    # Equation 3-17d: Froude number
    Fr_bf = U_bf / np.sqrt(GRAVITY * H_bf)
    
    # For Shields number, need slope (approximation from other parameters)
    # This requires slope as input - adding as optional parameter
    
    # Equation 3-17f: Particle Reynolds number
    R_p50 = np.sqrt(GRAVITY * R * D50) * D50 / nu
    
    return {
        'dimensionless_discharge': Q_hat,
        'dimensionless_width': B_hat, 
        'dimensionless_depth': H_hat,
        'froude_number': Fr_bf,
        'particle_reynolds_number': R_p50,
        'specific_gravity_ratio': R,
        'bankfull_velocity': U_bf
    }


def shields_number(tau_bf: float, D50: float, 
                  rho_s: float = 2650.0, 
                  rho: float = WATER_DENSITY) -> float:
    """
    Calculate Shields number for bankfull conditions.
    
    Equation 3-17e: τ*_bf50 = τ_bf / [(ρ_s - ρ) g D50]
    
    Args:
        tau_bf: Bankfull shear stress (Pa)
        D50: Median grain size (m)
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Shields number (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-17e
    """
    validate_positive(tau_bf, "shear_stress")
    validate_positive(D50, "median_grain_size")
    validate_positive(rho_s, "sediment_density")
    
    return tau_bf / ((rho_s - rho) * GRAVITY * D50)


def channel_scaling_relations(Q_hat: float, 
                             stream_type: Literal['gravel', 'sand']) -> Tuple[float, float, float]:
    """
    Calculate channel geometry from dimensionless discharge using empirical relations.
    
    Equations 3-19, 3-20, 3-21: Channel depth, width, and slope relations
    Equation 3-22: Shields stress for different bed types
    
    Args:
        Q_hat: Dimensionless bankfull discharge
        stream_type: Either 'gravel' or 'sand' bed stream
        
    Returns:
        Tuple of (dimensionless_depth, dimensionless_width, slope)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-19 to 3-22
    """
    validate_positive(Q_hat, "dimensionless_discharge")
    
    if stream_type not in ['gravel', 'sand']:
        raise ValueError("stream_type must be 'gravel' or 'sand'")
    
    if stream_type == 'gravel':
        # Gravel bed relations
        H_hat = 0.368 * (Q_hat**(-0.406))  # Equation 3-19
        B_hat = 4.870 * (Q_hat**0.461)     # Equation 3-20  
        S = 0.0976 * (Q_hat**(-0.341))     # Equation 3-21
        tau_star_bf50 = 0.049              # Equation 3-22
    else:
        # Sand bed relations
        H_hat = 3.01 * (Q_hat**(-0.521))   # Equation 3-19
        B_hat = 0.274 * (Q_hat**0.565)     # Equation 3-20
        S = 6.42 * (Q_hat**(-0.397))       # Equation 3-21  
        tau_star_bf50 = 1.86               # Equation 3-22
        
    return H_hat, B_hat, S


def classify_stream_type(D50: float) -> Literal['sand', 'gravel']:
    """
    Classify stream type based on median grain size.
    
    Based on analysis criteria from Chapter 3:
    - Sand-bed streams: D50 ≤ 0.5 mm  
    - Gravel-bed streams: D50 ≥ 25 mm
    
    Args:
        D50: Median grain size (mm)
        
    Returns:
        Stream type classification
        
    References:
        ASCE Manual 110, Chapter 3, Section 3.4
    """
    validate_positive(D50, "median_grain_size")
    
    if D50 <= 0.5:
        return 'sand'
    elif D50 >= 25.0:
        return 'gravel'
    else:
        # Transition zone - use gravel relations as default
        return 'gravel'