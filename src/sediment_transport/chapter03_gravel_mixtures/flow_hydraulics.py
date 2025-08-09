"""
Flow Hydraulics and 1D Flow Equations (Chapter 3, Equations 3-75 to 3-76)

This module implements flow hydraulic calculations including shear stress,
friction coefficients, and one-dimensional flow equations for sediment transport.
"""

import numpy as np
from typing import Union, Tuple
from ..utils.validators import validate_positive, validate_array
from ..utils.constants import GRAVITY, WATER_DENSITY


def shear_velocity(tau_b: Union[float, np.ndarray], 
                  rho: float = WATER_DENSITY) -> Union[float, np.ndarray]:
    """
    Calculate shear velocity from boundary shear stress.
    
    Equation 3-75a: u* = √(τ_b/ρ)
    
    Args:
        tau_b: Boundary shear stress (Pa)
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Shear velocity (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75a
    """
    validate_positive(tau_b, "boundary_shear_stress")
    validate_positive(rho, "water_density")
    
    return np.sqrt(tau_b / rho)


def chezy_coefficient(U: Union[float, np.ndarray], 
                     u_star: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate Chezy coefficient from velocity and shear velocity.
    
    Equation 3-75b: C_z = U/u*
    
    Args:
        U: Flow velocity (m/s)
        u_star: Shear velocity (m/s)
        
    Returns:
        Chezy coefficient (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75b
    """
    validate_positive(U, "flow_velocity")
    validate_positive(u_star, "shear_velocity")
    
    return U / u_star


def friction_factor(U: Union[float, np.ndarray], 
                   tau_b: Union[float, np.ndarray],
                   rho: float = WATER_DENSITY) -> Union[float, np.ndarray]:
    """
    Calculate friction factor from velocity and shear stress.
    
    Equation 3-75c: C_f = τ_b/(ρU²) = 2/C_z²
    
    Args:
        U: Flow velocity (m/s)
        tau_b: Boundary shear stress (Pa)
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Friction factor (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75c
    """
    validate_positive(U, "flow_velocity")
    validate_positive(tau_b, "boundary_shear_stress")
    validate_positive(rho, "water_density")
    
    return tau_b / (rho * U**2)


def chezy_coefficient_logarithmic(H: Union[float, np.ndarray], 
                                 k_s: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate Chezy coefficient using logarithmic velocity profile.
    
    Equation 3-75d: C_z = 2.5 ln(11H/k_s)
    
    Args:
        H: Flow depth (m)
        k_s: Roughness height (m)
        
    Returns:
        Chezy coefficient (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75d
    """
    validate_positive(H, "flow_depth")
    validate_positive(k_s, "roughness_height")
    
    # Ensure H > k_s for physical validity
    ratio = H / k_s
    if np.any(ratio <= 1):
        raise ValueError("Flow depth must be greater than roughness height")
    
    return 2.5 * np.log(11.0 * ratio)


def chezy_coefficient_power_law(H: Union[float, np.ndarray], 
                               k_s: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate Chezy coefficient using power law relation.
    
    Equation 3-75e: C_z = 8.1(H/k_s)^(1/6)
    
    Args:
        H: Flow depth (m)
        k_s: Roughness height (m)
        
    Returns:
        Chezy coefficient (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75e
    """
    validate_positive(H, "flow_depth")
    validate_positive(k_s, "roughness_height")
    
    ratio = H / k_s
    return 8.1 * (ratio**(1.0/6.0))


def roughness_height(D90: Union[float, np.ndarray], 
                    n_k: float = 3.5) -> Union[float, np.ndarray]:
    """
    Calculate roughness height from grain size.
    
    Equation 3-75f: k_s = n_k D90
    
    Args:
        D90: 90th percentile grain size (m)
        n_k: Roughness coefficient, default 3.5
        
    Returns:
        Roughness height (m)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-75f
    """
    validate_positive(D90, "D90_grain_size")
    validate_positive(n_k, "roughness_coefficient")
    
    return n_k * D90


def boundary_shear_stress(H: Union[float, np.ndarray], 
                         S: Union[float, np.ndarray],
                         rho: float = WATER_DENSITY) -> Union[float, np.ndarray]:
    """
    Calculate boundary shear stress from depth and slope.
    
    Equation 3-76a: τ_b = ρ g H S
    
    Args:
        H: Flow depth (m)
        S: Bed slope (dimensionless)
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Boundary shear stress (Pa)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-76a
    """
    validate_positive(H, "flow_depth")
    validate_positive(S, "bed_slope")
    validate_positive(rho, "water_density")
    
    return rho * GRAVITY * H * S


def shear_velocity_from_geometry(H: Union[float, np.ndarray], 
                                S: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate shear velocity from depth and slope.
    
    Equation 3-76b: u* = √(g H S)
    
    Args:
        H: Flow depth (m)
        S: Bed slope (dimensionless)
        
    Returns:
        Shear velocity (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-76b
    """
    validate_positive(H, "flow_depth")
    validate_positive(S, "bed_slope")
    
    return np.sqrt(GRAVITY * H * S)


def continuity_equation_1d(H: Union[float, np.ndarray],
                          U: Union[float, np.ndarray],
                          dx: float = 1.0,
                          dt: float = 1.0) -> Union[float, np.ndarray]:
    """
    One-dimensional continuity equation for shallow water flow.
    
    Equation 3-76c: ∂H/∂t + ∂(UH)/∂s = 0
    
    Args:
        H: Flow depth (m)
        U: Flow velocity (m/s)
        dx: Spatial grid spacing (m), default 1.0
        dt: Time step (s), default 1.0
        
    Returns:
        Rate of depth change (m/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-76c
    """
    validate_positive(H, "flow_depth")
    
    # Calculate spatial derivative ∂(UH)/∂s
    discharge = U * H
    if hasattr(discharge, '__len__') and len(discharge) > 1:
        dUH_dx = np.gradient(discharge, dx)
    else:
        dUH_dx = 0.0
    
    # ∂H/∂t = -∂(UH)/∂s
    dH_dt = -dUH_dx
    
    return dH_dt


def momentum_equation_1d(U: Union[float, np.ndarray],
                        H: Union[float, np.ndarray],
                        S: Union[float, np.ndarray],
                        C_f: Union[float, np.ndarray],
                        dx: float = 1.0,
                        dt: float = 1.0) -> Union[float, np.ndarray]:
    """
    One-dimensional momentum equation for shallow water flow.
    
    Equation 3-76d: ∂U/∂t + U∂U/∂s = -g∂H/∂s + gS - C_f U²/H
    
    Args:
        U: Flow velocity (m/s)
        H: Flow depth (m) 
        S: Bed slope (dimensionless)
        C_f: Friction factor (dimensionless)
        dx: Spatial grid spacing (m), default 1.0
        dt: Time step (s), default 1.0
        
    Returns:
        Rate of velocity change (m/s²)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-76d
    """
    validate_positive(H, "flow_depth")
    
    # Calculate spatial derivatives
    if hasattr(U, '__len__') and len(U) > 1:
        dU_dx = np.gradient(U, dx)
        dH_dx = np.gradient(H, dx)
    else:
        dU_dx = 0.0
        dH_dx = 0.0
    
    # Momentum equation terms
    advection = U * dU_dx
    pressure_gradient = -GRAVITY * dH_dx
    slope_term = GRAVITY * S
    friction_term = -C_f * U**2 / H
    
    # ∂U/∂t = -U∂U/∂s - g∂H/∂s + gS - C_f U²/H
    dU_dt = -advection + pressure_gradient + slope_term + friction_term
    
    return dU_dt


def flow_hydraulics_complete(H: float, S: float, U: float,
                           D90: float, n_k: float = 3.5,
                           rho: float = WATER_DENSITY) -> dict:
    """
    Complete hydraulic calculations for given flow conditions.
    
    Combines Equations 3-75 and 3-76 for comprehensive flow analysis.
    
    Args:
        H: Flow depth (m)
        S: Bed slope (dimensionless)
        U: Flow velocity (m/s)
        D90: 90th percentile grain size (m)
        n_k: Roughness coefficient, default 3.5
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Dictionary of hydraulic parameters
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-75, 3-76
    """
    validate_positive(H, "flow_depth")
    validate_positive(S, "bed_slope") 
    validate_positive(U, "flow_velocity")
    validate_positive(D90, "D90_grain_size")
    
    # Calculate roughness height
    k_s = roughness_height(D90, n_k)
    
    # Calculate shear parameters
    tau_b = boundary_shear_stress(H, S, rho)
    u_star = shear_velocity(tau_b, rho)
    
    # Calculate friction coefficients
    C_z = chezy_coefficient(U, u_star)
    C_f = friction_factor(U, tau_b, rho)
    
    # Alternative Chezy coefficients
    C_z_log = chezy_coefficient_logarithmic(H, k_s)
    C_z_power = chezy_coefficient_power_law(H, k_s)
    
    return {
        'boundary_shear_stress': tau_b,
        'shear_velocity': u_star,
        'chezy_coefficient': C_z,
        'friction_factor': C_f,
        'chezy_logarithmic': C_z_log,
        'chezy_power_law': C_z_power,
        'roughness_height': k_s,
        'froude_number': U / np.sqrt(GRAVITY * H),
        'reynolds_number': U * H / 1e-6  # Assuming nu = 1e-6 m²/s
    }