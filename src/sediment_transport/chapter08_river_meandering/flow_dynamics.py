"""
Flow and Bed Topography in Meanders

This module implements models for flow distribution, secondary currents,
and bed topography development in meandering channels.

Equations implemented:
    - 8-13, 8-14: Flow distribution in curved channels
    - 8-15: Sediment transport power law
    - 8-16, 8-17: Differential equations for flow development
    - 8-18 to 8-28: Parameters and solutions for bend flow

References:
    ASCE Manual 110, Chapter 8, Section 8.3
    Odgaard (1989), Engelund (1974), Struiksma et al. (1985)
"""

import numpy as np
from typing import Union, Tuple, Dict, Optional
from scipy.integrate import odeint
from scipy.optimize import fsolve
from ..utils.validators import validate_positive, validate_range

def flow_velocity_distribution(transverse_position_n: float,
                             centerline_velocity_uc: float,
                             normalized_velocity_gradient_Utc: float) -> float:
    """
    Calculate velocity distribution across channel section in curved flow.
    
    Assumes linear variation over central portion of cross section.
    
    Parameters:
    -----------
    transverse_position_n : float
        Transverse position from centerline (m), positive toward concave bank
    centerline_velocity_uc : float
        Centerline velocity (m/s)
    normalized_velocity_gradient_Utc : float
        Normalized transverse velocity gradient at centerline
        
    Returns:
    --------
    float
        Velocity at position n (m/s)
        
    Equation: 8-13
        u = uc * (1 + n * Utc)
        
    Reference: Odgaard (1989)
    """
    validate_positive(centerline_velocity_uc, "centerline_velocity_uc")
    
    return centerline_velocity_uc * (1.0 + transverse_position_n * normalized_velocity_gradient_Utc)

def flow_depth_distribution(transverse_position_n: float,
                          centerline_depth_dc: float,
                          transverse_bed_slope_Stc: float) -> float:
    """
    Calculate depth distribution across channel section in curved flow.
    
    Parameters:
    -----------
    transverse_position_n : float
        Transverse position from centerline (m)
    centerline_depth_dc : float
        Centerline depth (m)
    transverse_bed_slope_Stc : float
        Transverse bed slope at centerline (∂d/∂n)c
        
    Returns:
    --------
    float
        Depth at position n (m)
        
    Equation: 8-14
        d = dc * (1 + n * Stc)
        
    Reference: Odgaard (1989)
    """
    validate_positive(centerline_depth_dc, "centerline_depth_dc")
    
    return centerline_depth_dc * (1.0 + transverse_position_n * transverse_bed_slope_Stc)

def sediment_transport_power_law(local_velocity_u: float,
                               centerline_velocity_uc: float,
                               centerline_transport_qsc: float,
                               transport_exponent_M: float) -> float:
    """
    Calculate sediment transport using power law relationship.
    
    Parameters:
    -----------
    local_velocity_u : float
        Local velocity (m/s)
    centerline_velocity_uc : float
        Centerline velocity (m/s)
    centerline_transport_qsc : float
        Centerline sediment transport rate (m²/s)
    transport_exponent_M : float
        Transport exponent (typically 2-4)
        
    Returns:
    --------
    float
        Local sediment transport rate (m²/s)
        
    Equation: 8-15
        qs = qsc * (u/uc)^M
        
    Reference: Simons & Sentürk (1977)
    """
    validate_positive(local_velocity_u, "local_velocity_u")
    validate_positive(centerline_velocity_uc, "centerline_velocity_uc")
    validate_positive(centerline_transport_qsc, "centerline_transport_qsc")
    validate_range(transport_exponent_M, 1.0, 6.0, "transport_exponent_M")
    
    return centerline_transport_qsc * ((local_velocity_u / centerline_velocity_uc) ** transport_exponent_M)

def calculate_flow_parameters(channel_width_b: float,
                            centerline_depth_dc: float,
                            friction_parameter_m: float,
                            transport_exponent_M: float,
                            transverse_bed_slope_factor_B: float,
                            transverse_mass_flux_factor_alpha: float,
                            critical_shields_parameter_theta: float,
                            particle_froude_number_FDc: float) -> Dict[str, float]:
    """
    Calculate dimensionless parameters for bend flow analysis.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    centerline_depth_dc : float
        Centerline depth (m)
    friction_parameter_m : float
        Friction parameter (κ√8/f)
    transport_exponent_M : float
        Sediment transport exponent
    transverse_bed_slope_factor_B : float
        Transverse bed slope factor
    transverse_mass_flux_factor_alpha : float
        Transverse mass flux factor
    critical_shields_parameter_theta : float
        Critical Shields parameter
    particle_froude_number_FDc : float
        Particle densimetric Froude number
        
    Returns:
    --------
    dict
        Dictionary of calculated parameters (a1 through a6)
        
    Equations: 8-18 to 8-23
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(centerline_depth_dc, "centerline_depth_dc")
    validate_positive(friction_parameter_m, "friction_parameter_m")
    validate_range(transport_exponent_M, 1.0, 6.0, "transport_exponent_M")
    validate_positive(transverse_bed_slope_factor_B, "transverse_bed_slope_factor_B")
    validate_positive(transverse_mass_flux_factor_alpha, "transverse_mass_flux_factor_alpha")
    validate_positive(critical_shields_parameter_theta, "critical_shields_parameter_theta")
    validate_positive(particle_froude_number_FDc, "particle_froude_number_FDc")
    
    # von Kármán constant
    kappa = 0.4
    
    # Dimensionless coordinate
    sigma = 0.0  # This would be s/b in actual calculations
    
    # Calculate parameters
    m2 = (1 - transport_exponent_M) / (transport_exponent_M + 1) + 2
    
    # Equation 8-18
    a1 = (kappa**2 * channel_width_b**2) / (2 * friction_parameter_m**2 * centerline_depth_dc**2)
    
    # Equation 8-19
    a2 = 1 - (transport_exponent_M + 1) / (transport_exponent_M + 2)
    
    # Equation 8-20
    term1 = (8 * transverse_bed_slope_factor_B * np.sqrt(critical_shields_parameter_theta) * 
             friction_parameter_m * (transport_exponent_M + 1) * centerline_depth_dc) / \
            (kappa**2 * particle_froude_number_FDc * (transport_exponent_M + 2) * channel_width_b)
    term2 = (2 * kappa**2 * friction_parameter_m) / ((transport_exponent_M + 1) * (transport_exponent_M + 2) * centerline_depth_dc)
    a3 = term1 + term2
    
    # Equation 8-21
    a4 = (2 * kappa**2 * (1 + transverse_mass_flux_factor_alpha * (1 - transport_exponent_M))) / \
         ((transport_exponent_M + 1) * (transport_exponent_M + 2))
    
    # Equation 8-22
    a5 = (8 * critical_shields_parameter_theta * transverse_mass_flux_factor_alpha) / \
         (kappa**2 * (transport_exponent_M + 2) * particle_froude_number_FDc * 
          (1 + 2 * friction_parameter_m**2 / (transport_exponent_M + 1)))
    
    # Equation 8-23
    a6 = (8 * transverse_mass_flux_factor_alpha * friction_parameter_m * (transport_exponent_M + 2) * 
          centerline_depth_dc) / (kappa**2 * channel_width_b)
    
    return {
        'a1': a1,
        'a2': a2,
        'a3': a3,
        'a4': a4,
        'a5': a5,
        'a6': a6,
        'm2': m2
    }

def bend_flow_differential_equations(y: np.ndarray, sigma: float, params: Dict[str, float]) -> np.ndarray:
    """
    System of differential equations for bend flow development.
    
    Parameters:
    -----------
    y : array
        State vector [Utc, Stc] (normalized velocity gradient, transverse bed slope)
    sigma : float
        Dimensionless streamwise coordinate (s/b)
    params : dict
        Flow parameters (a1-a6)
        
    Returns:
    --------
    array
        Derivatives [dUtc/dσ, dStc/dσ]
        
    Equations: 8-16, 8-17
    """
    Utc, Stc = y
    a1, a2, a3, a4, a5, a6 = params['a1'], params['a2'], params['a3'], params['a4'], params['a5'], params['a6']
    
    # Equation 8-16
    dUtc_dsigma = a1 + 0.5 * a2 * Stc
    
    # Equation 8-17
    d2Stc_dsigma2 = a3 + a4 * (Stc + a5 * Utc) + a6
    
    # For first-order system, we need dStc/dσ
    # This requires solving the full second-order equation
    # For simplification, we'll use a quasi-steady approximation
    dStc_dsigma = np.sqrt(abs(d2Stc_dsigma2)) * np.sign(d2Stc_dsigma2)
    
    return np.array([dUtc_dsigma, dStc_dsigma])

def fully_developed_bend_flow(particle_froude_number_FDc: float,
                            radius_curvature_rc: float,
                            friction_parameter_m: float,
                            transverse_bed_slope_factor_B: float,
                            critical_shields_parameter_theta: float) -> Tuple[float, float]:
    """
    Calculate fully developed bend flow characteristics.
    
    For constant-radius bends where d/dσ = 0.
    
    Parameters:
    -----------
    particle_froude_number_FDc : float
        Particle densimetric Froude number
    radius_curvature_rc : float
        Radius of curvature (m)
    friction_parameter_m : float
        Friction parameter
    transverse_bed_slope_factor_B : float
        Transverse bed slope factor
    critical_shields_parameter_theta : float
        Critical Shields parameter
        
    Returns:
    --------
    tuple
        (normalized_velocity_gradient, transverse_bed_slope)
        
    Equations: 8-25, 8-26
    """
    validate_positive(particle_froude_number_FDc, "particle_froude_number_FDc")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_positive(friction_parameter_m, "friction_parameter_m")
    validate_positive(transverse_bed_slope_factor_B, "transverse_bed_slope_factor_B")
    validate_positive(critical_shields_parameter_theta, "critical_shields_parameter_theta")
    
    # Equation 8-25: Fully developed condition
    # u'/uc = d'/dc (continuity requirement)
    velocity_depth_ratio = 1.0  # From continuity
    
    # Equation 8-26: Fully developed transverse bed slope
    H_factor = ((2 * friction_parameter_m + 1) * (friction_parameter_m + 1)) / \
               (transverse_bed_slope_factor_B * 0.4 * np.sqrt(friction_parameter_m * 
                (friction_parameter_m + 1 + 2 * friction_parameter_m**2)))
    
    transverse_bed_slope = particle_froude_number_FDc / radius_curvature_rc * H_factor
    
    return velocity_depth_ratio, transverse_bed_slope

def simulate_bend_flow_development(channel_width_b: float,
                                 centerline_depth_dc: float,
                                 curvature_function: callable,
                                 flow_params: Dict[str, float],
                                 streamwise_distance: np.ndarray,
                                 initial_conditions: Tuple[float, float] = (0.0, 0.0)) -> Dict[str, np.ndarray]:
    """
    Simulate the development of flow and bed topography through a bend.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    centerline_depth_dc : float
        Centerline depth (m)
    curvature_function : callable
        Function returning curvature as function of streamwise distance
    flow_params : dict
        Flow parameters from calculate_flow_parameters()
    streamwise_distance : array
        Streamwise distance coordinates (m)
    initial_conditions : tuple, optional
        Initial (Utc, Stc) values
        
    Returns:
    --------
    dict
        Simulation results with velocity gradients, bed slopes, depths
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(centerline_depth_dc, "centerline_depth_dc")
    
    # Convert to dimensionless coordinate
    sigma_array = streamwise_distance / channel_width_b
    
    # Solve differential equation system
    solution = odeint(bend_flow_differential_equations, 
                     initial_conditions, 
                     sigma_array, 
                     args=(flow_params,))
    
    Utc_array = solution[:, 0]
    Stc_array = solution[:, 1]
    
    # Calculate derived quantities
    results = {
        'distance': streamwise_distance,
        'sigma': sigma_array,
        'velocity_gradient_Utc': Utc_array,
        'bed_slope_Stc': Stc_array,
        'curvature': np.array([curvature_function(s) for s in streamwise_distance])
    }
    
    # Add transverse velocity and depth distributions
    transverse_positions = np.linspace(-channel_width_b/2, channel_width_b/2, 21)
    
    # Calculate distributions at each streamwise location
    velocity_distributions = []
    depth_distributions = []
    
    for i, s in enumerate(streamwise_distance):
        u_centerline = 1.0  # Normalized
        d_centerline = centerline_depth_dc
        
        u_dist = [flow_velocity_distribution(n, u_centerline, Utc_array[i]) 
                  for n in transverse_positions]
        d_dist = [flow_depth_distribution(n, d_centerline, Stc_array[i]) 
                  for n in transverse_positions]
        
        velocity_distributions.append(u_dist)
        depth_distributions.append(d_dist)
    
    results['transverse_positions'] = transverse_positions
    results['velocity_distributions'] = np.array(velocity_distributions)
    results['depth_distributions'] = np.array(depth_distributions)
    
    return results

def secondary_current_strength(radius_curvature_rc: float,
                             channel_width_b: float,
                             centerline_velocity_uc: float,
                             friction_parameter_m: float) -> float:
    """
    Estimate secondary current strength in curved channels.
    
    Parameters:
    -----------
    radius_curvature_rc : float
        Radius of curvature (m)
    channel_width_b : float
        Channel width (m)
    centerline_velocity_uc : float
        Centerline velocity (m/s)
    friction_parameter_m : float
        Friction parameter
        
    Returns:
    --------
    float
        Secondary current velocity scale (m/s)
    """
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(centerline_velocity_uc, "centerline_velocity_uc")
    validate_positive(friction_parameter_m, "friction_parameter_m")
    
    # Characteristic secondary current velocity
    # Based on centrifugal acceleration and friction effects
    secondary_velocity = (centerline_velocity_uc * channel_width_b) / (friction_parameter_m * radius_curvature_rc)
    
    return secondary_velocity

def water_surface_superelevation(centerline_velocity_uc: float,
                                radius_curvature_rc: float,
                                gravity: float = 9.81) -> float:
    """
    Calculate water surface superelevation in curved channels.
    
    Parameters:
    -----------
    centerline_velocity_uc : float
        Centerline velocity (m/s)
    radius_curvature_rc : float
        Radius of curvature (m)
    gravity : float, optional
        Gravitational acceleration (m/s²)
        
    Returns:
    --------
    float
        Water surface superelevation (m)
    """
    validate_positive(centerline_velocity_uc, "centerline_velocity_uc")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_positive(gravity, "gravity")
    
    # Based on centrifugal force balance
    superelevation = (centerline_velocity_uc**2) / (gravity * radius_curvature_rc)
    
    return superelevation