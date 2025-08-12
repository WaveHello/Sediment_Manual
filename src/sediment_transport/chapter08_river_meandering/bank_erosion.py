"""
Bank Erosion and Channel Migration

This module implements models for bank erosion rates and channel migration
based on flow velocity, channel curvature, and hydraulic conditions.

Equations implemented:
    - 8-8: Nanson-Hickin migration rate relationship
    - 8-9: Ikeda velocity-based bank erosion model
    - 8-10: Parker-Andrews migration rate model
    - 8-11: Odgaard depth-based bank erosion model
    - 8-12: Combined velocity-depth bank erosion model

References:
    ASCE Manual 110, Chapter 8, Section 8.2.3
    Hickin & Nanson (1975), Ikeda et al. (1981), Parker (1983), Odgaard (1987)
"""

import numpy as np
from typing import Union, Tuple, Dict, Optional
from ..utils.validators import validate_positive, validate_range

def nanson_hickin_migration_rate(channel_width_b: float,
                               radius_curvature_rc: float,
                               method: str = 'general') -> float:
    """
    Calculate bank erosion rate using Nanson-Hickin relationship.
    
    Based on dendrochronology data from the Beatton River, Canada.
    Shows that channel curvature plays an important role in migration rates.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    radius_curvature_rc : float
        Radius of curvature (m)
    method : str, optional
        Calculation method: 'general' or 'optimized'
        
    Returns:
    --------
    float
        Erosion rate (m/year)
        
    Equation: 8-8
        ve = 2.0 * (b/rc)^0.32  [for b/rc ≤ 0.32]
        ve = 0.2 * (b/rc)^(-0.32) [for b/rc > 0.32]
        
    Reference: Nanson & Hickin (1983)
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    
    width_radius_ratio = channel_width_b / radius_curvature_rc
    
    if method == 'general':
        # General form with transition at b/rc = 0.32
        if width_radius_ratio <= 0.32:
            erosion_rate = 2.0 * (width_radius_ratio ** 0.32)
        else:
            erosion_rate = 0.2 * (width_radius_ratio ** -0.32)
    
    elif method == 'optimized':
        # Optimized form that peaks around b/rc = 0.32
        erosion_rate = 2.0 * np.exp(-((np.log(width_radius_ratio) + 1.139)**2) / 0.5)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return erosion_rate

def ikeda_bank_erosion(near_bank_velocity_ub: float,
                      reach_average_velocity_u: float,
                      erodibility_parameter_E: float) -> float:
    """
    Calculate bank erosion rate using Ikeda velocity-based model.
    
    Assumes bank particles are eroded when near-bank velocity exceeds
    the reach-averaged velocity.
    
    Parameters:
    -----------
    near_bank_velocity_ub : float
        Near-bank depth-averaged velocity (m/s)
    reach_average_velocity_u : float
        Reach-averaged velocity (m/s)
    erodibility_parameter_E : float
        Bank erodibility parameter (1/s)
        
    Returns:
    --------
    float
        Bank erosion rate (m/s)
        
    Equation: 8-9
        ve = E * (ub - u)
        
    Reference: Ikeda et al. (1981)
    """
    validate_positive(near_bank_velocity_ub, "near_bank_velocity_ub")
    validate_positive(reach_average_velocity_u, "reach_average_velocity_u")
    validate_positive(erodibility_parameter_E, "erodibility_parameter_E")
    
    velocity_increment = near_bank_velocity_ub - reach_average_velocity_u
    
    # Erosion only occurs when near-bank velocity exceeds average
    if velocity_increment > 0:
        return erodibility_parameter_E * velocity_increment
    else:
        return 0.0

def parker_andrews_migration_rate(scour_factor_A: float,
                                 reach_average_velocity_u: float,
                                 channel_width_b: float,
                                 radius_curvature_rc: float) -> float:
    """
    Calculate migration rate using Parker-Andrews model for developed bend flow.
    
    For constant-radius curves with developed flow conditions.
    
    Parameters:
    -----------
    scour_factor_A : float
        Scour factor parameterizing secondary currents (order of 1)
    reach_average_velocity_u : float
        Reach-averaged velocity (m/s)
    channel_width_b : float
        Channel width (m)
    radius_curvature_rc : float
        Radius of curvature (m)
        
    Returns:
    --------
    float
        Bank erosion rate (m/s)
        
    Equation: 8-10
        ve = (EA * u * b) / rc
        
    Reference: Parker (1983), Parker & Andrews (1986)
    """
    validate_positive(scour_factor_A, "scour_factor_A")
    validate_positive(reach_average_velocity_u, "reach_average_velocity_u")
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    
    return (scour_factor_A * reach_average_velocity_u * channel_width_b) / radius_curvature_rc

def odgaard_depth_based_erosion(near_bank_depth_db: float,
                               centerline_depth_dc: float,
                               reach_average_velocity_u: float,
                               erosion_parameter_E_prime: float) -> float:
    """
    Calculate bank erosion rate based on near-bank depth increment.
    
    Relates erosion to bank stability as affected by scour depth.
    As bank height increases, stability decreases.
    
    Parameters:
    -----------
    near_bank_depth_db : float
        Near-bank depth (m)
    centerline_depth_dc : float
        Centerline depth (m)
    reach_average_velocity_u : float
        Reach-averaged velocity (m/s)
    erosion_parameter_E_prime : float
        Erosion parameter (1/s)
        
    Returns:
    --------
    float
        Bank erosion rate (m/s)
        
    Equation: 8-11
        ve = E' * u * (db - dc)
        
    Reference: Odgaard (1989a)
    """
    validate_positive(near_bank_depth_db, "near_bank_depth_db")
    validate_positive(centerline_depth_dc, "centerline_depth_dc")
    validate_positive(reach_average_velocity_u, "reach_average_velocity_u")
    validate_positive(erosion_parameter_E_prime, "erosion_parameter_E_prime")
    
    depth_increment = near_bank_depth_db - centerline_depth_dc
    
    # Erosion proportional to depth increment
    return erosion_parameter_E_prime * reach_average_velocity_u * max(0, depth_increment)

def combined_erosion_model(near_bank_velocity_ub: float,
                         reach_average_velocity_u: float,
                         near_bank_depth_db: float,
                         centerline_depth_dc: float,
                         weight_factor_C1: float,
                         weight_factor_C2: float,
                         erodibility_parameter_E: float) -> float:
    """
    Calculate bank erosion using combined velocity-depth model.
    
    Combines the effects of both velocity and depth increments
    with weighting factors.
    
    Parameters:
    -----------
    near_bank_velocity_ub : float
        Near-bank velocity (m/s)
    reach_average_velocity_u : float
        Reach-averaged velocity (m/s)
    near_bank_depth_db : float
        Near-bank depth (m)
    centerline_depth_dc : float
        Centerline depth (m)
    weight_factor_C1 : float
        Velocity weighting factor (positive)
    weight_factor_C2 : float
        Depth weighting factor (can be positive, negative, or zero)
    erodibility_parameter_E : float
        Bank erodibility parameter (1/s)
        
    Returns:
    --------
    float
        Bank erosion rate (m/s)
        
    Equation: 8-12
        ve = E * [C1 * (ub - u) + C2 * (db - dc)]
        
    Reference: Howard (1992)
    """
    validate_positive(near_bank_velocity_ub, "near_bank_velocity_ub")
    validate_positive(reach_average_velocity_u, "reach_average_velocity_u")
    validate_positive(near_bank_depth_db, "near_bank_depth_db")
    validate_positive(centerline_depth_dc, "centerline_depth_dc")
    validate_positive(weight_factor_C1, "weight_factor_C1")
    validate_positive(erodibility_parameter_E, "erodibility_parameter_E")
    
    velocity_increment = near_bank_velocity_ub - reach_average_velocity_u
    depth_increment = near_bank_depth_db - centerline_depth_dc
    
    erosion_term = weight_factor_C1 * velocity_increment + weight_factor_C2 * depth_increment
    
    return erodibility_parameter_E * max(0, erosion_term)

def width_based_migration_rate(channel_width_b: float,
                             empirical_coefficient: float = 0.01) -> float:
    """
    Calculate migration rate based on empirical width relationship (Brice).
    
    Demonstrates that erosion rate increases with channel width.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    empirical_coefficient : float, optional
        Empirical coefficient (default=0.01)
        
    Returns:
    --------
    float
        Mean erosion rate (m/year)
        
    Reference: Brice (1982)
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(empirical_coefficient, "empirical_coefficient")
    
    return empirical_coefficient * channel_width_b

def drainage_area_migration_rate(drainage_area_A: float,
                               empirical_coefficient: float = 0.05) -> float:
    """
    Calculate migration rate based on drainage area relationship (Hooke).
    
    Indirect relationship through correlation with drainage area.
    
    Parameters:
    -----------
    drainage_area_A : float
        Drainage area (km²)
    empirical_coefficient : float, optional
        Empirical coefficient (default=0.05)
        
    Returns:
    --------
    float
        Mean erosion rate (m/year)
        
    Reference: Hooke (1980)
    """
    validate_positive(drainage_area_A, "drainage_area_A")
    validate_positive(empirical_coefficient, "empirical_coefficient")
    
    return empirical_coefficient * np.sqrt(drainage_area_A)

def migration_rate_comparison(channel_width_b: float,
                            radius_curvature_rc: float,
                            velocity_params: Dict[str, float],
                            depth_params: Dict[str, float],
                            erodibility_E: float) -> Dict[str, float]:
    """
    Compare migration rates calculated by different methods.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    radius_curvature_rc : float
        Radius of curvature (m)
    velocity_params : dict
        Parameters for velocity-based models
    depth_params : dict
        Parameters for depth-based models
    erodibility_E : float
        Bank erodibility parameter
        
    Returns:
    --------
    dict
        Comparison of migration rates from different models
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_positive(erodibility_E, "erodibility_E")
    
    results = {}
    
    # Nanson-Hickin method
    results['nanson_hickin'] = nanson_hickin_migration_rate(channel_width_b, radius_curvature_rc)
    
    # Ikeda method
    if all(k in velocity_params for k in ['near_bank_velocity', 'reach_average_velocity']):
        results['ikeda'] = ikeda_bank_erosion(
            velocity_params['near_bank_velocity'],
            velocity_params['reach_average_velocity'],
            erodibility_E
        )
    
    # Parker-Andrews method
    if 'reach_average_velocity' in velocity_params:
        scour_factor = velocity_params.get('scour_factor', 1.0)
        results['parker_andrews'] = parker_andrews_migration_rate(
            scour_factor,
            velocity_params['reach_average_velocity'],
            channel_width_b,
            radius_curvature_rc
        )
    
    # Odgaard depth-based method
    if all(k in depth_params for k in ['near_bank_depth', 'centerline_depth']):
        velocity = velocity_params.get('reach_average_velocity', 1.0)
        results['odgaard_depth'] = odgaard_depth_based_erosion(
            depth_params['near_bank_depth'],
            depth_params['centerline_depth'],
            velocity,
            erodibility_E
        )
    
    # Width-based empirical method
    results['width_based'] = width_based_migration_rate(channel_width_b)
    
    # Add statistics
    if len(results) > 1:
        values = list(results.values())
        results['mean'] = np.mean(values)
        results['std'] = np.std(values)
        results['coefficient_of_variation'] = np.std(values) / np.mean(values) if np.mean(values) > 0 else 0
    
    return results

def bank_failure_mechanisms(bank_height: float,
                          bank_angle: float,
                          cohesion: float,
                          friction_angle: float,
                          unit_weight: float) -> Dict[str, Union[float, bool]]:
    """
    Assess bank failure mechanisms beyond flow-induced erosion.
    
    Includes mass wasting, piping, and slope stability analysis.
    
    Parameters:
    -----------
    bank_height : float
        Bank height (m)
    bank_angle : float
        Bank angle from horizontal (degrees)
    cohesion : float
        Soil cohesion (Pa)
    friction_angle : float
        Soil internal friction angle (degrees)
    unit_weight : float
        Soil unit weight (N/m³)
        
    Returns:
    --------
    dict
        Bank stability assessment
    """
    validate_positive(bank_height, "bank_height")
    validate_range(bank_angle, 0, 90, "bank_angle")
    validate_positive(cohesion, "cohesion")
    validate_range(friction_angle, 0, 45, "friction_angle")
    validate_positive(unit_weight, "unit_weight")
    
    # Convert angles to radians
    bank_angle_rad = np.radians(bank_angle)
    friction_angle_rad = np.radians(friction_angle)
    
    # Critical height for cohesive failure (Coulomb analysis)
    critical_height = (4 * cohesion) / (unit_weight * np.tan(np.pi/4 - friction_angle_rad/2))
    
    # Factor of safety
    if bank_height > 0:
        factor_of_safety = critical_height / bank_height
    else:
        factor_of_safety = float('inf')
    
    # Stability assessment
    stable = factor_of_safety > 1.2  # Standard factor of safety
    failure_risk = 'low' if factor_of_safety > 1.5 else ('moderate' if factor_of_safety > 1.0 else 'high')
    
    results = {
        'critical_height': critical_height,
        'factor_of_safety': factor_of_safety,
        'stable': stable,
        'failure_risk': failure_risk,
        'mechanism': 'cohesive_failure' if not stable else 'stable'
    }
    
    return results