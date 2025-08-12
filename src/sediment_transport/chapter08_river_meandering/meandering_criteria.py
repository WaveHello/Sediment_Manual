"""
Meandering Criteria and Stability Thresholds

This module implements criteria for determining when channels will meander
versus braid, based on relationships between channel slope, discharge, 
sediment size, and bank characteristics.

Equations implemented:
    - 8-1: Leopold-Wolman threshold criterion
    - 8-2: Henderson refined criterion with grain size
    - 8-3: Lane criterion for sand-bed rivers  
    - 8-4: Millar criterion with bank vegetation influence

References:
    ASCE Manual 110, Chapter 8, Section 8.2.1
    Leopold & Wolman (1957), Henderson (1963), Lane (1957), Millar (2000)
"""

import numpy as np
from typing import Union, Tuple
from ..utils.validators import validate_positive, validate_range

def leopold_wolman_threshold(discharge_Q: float) -> float:
    """
    Calculate the Leopold-Wolman threshold slope between meandering and braided channels.
    
    For slopes above this threshold, channels tend to braid.
    For slopes below this threshold, channels tend to meander.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
        
    Returns:
    --------
    float
        Threshold slope (dimensionless)
        
    Equation: 8-1
        S = 0.012 * Q^(-0.44)
        
    Reference: Leopold & Wolman (1957)
    """
    validate_positive(discharge_Q, "discharge_Q")
    
    return 0.012 * (discharge_Q ** -0.44)

def henderson_threshold(discharge_Q: float, median_diameter_D: float) -> float:
    """
    Calculate Henderson's refined threshold criterion accounting for sediment size.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
    median_diameter_D : float
        Median particle diameter (mm)
        
    Returns:
    --------
    float
        Threshold slope (dimensionless)
        
    Equation: 8-2
        S = 0.0002 * D^1.14 * Q^(-0.44)
        
    Reference: Henderson (1963)
    """
    validate_positive(discharge_Q, "discharge_Q")
    validate_positive(median_diameter_D, "median_diameter_D")
    
    return 0.0002 * (median_diameter_D ** 1.14) * (discharge_Q ** -0.44)

def lane_criterion(discharge_Q: float, K_constant: float = 0.0007) -> float:
    """
    Calculate Lane's criterion for sand-bed rivers.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
    K_constant : float, optional
        Lane's constant (default=0.0007 for metric units)
        0.0017 for English units
        
    Returns:
    --------
    float
        Threshold slope (dimensionless)
        
    Equation: 8-3
        S = K * Q^(-0.25)
        
    Reference: Lane (1957)
    """
    validate_positive(discharge_Q, "discharge_Q")
    validate_positive(K_constant, "K_constant")
    
    return K_constant * (discharge_Q ** -0.25)

def millar_criterion(discharge_Q: float, 
                    median_diameter_D: float,
                    bank_friction_angle_phi: float) -> float:
    """
    Calculate Millar's criterion including bank vegetation effects.
    
    This criterion accounts for the stabilizing effect of bank vegetation
    through an enhanced friction angle parameter.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
    median_diameter_D : float
        Median sediment diameter for banks and bed surface (m)
    bank_friction_angle_phi : float
        Bank sediment friction angle (degrees)
        Without vegetation: ~40° for coarse gravel
        With vegetation: can exceed 40°
        
    Returns:
    --------
    float
        Threshold slope (dimensionless)
        
    Equation: 8-4
        S = 0.0002 * (φ^1.75 * D^0.61) / Q^0.25
        
    Reference: Millar (2000)
    """
    validate_positive(discharge_Q, "discharge_Q")
    validate_positive(median_diameter_D, "median_diameter_D")
    validate_range(bank_friction_angle_phi, 20.0, 60.0, "bank_friction_angle_phi")
    
    return 0.0002 * ((bank_friction_angle_phi ** 1.75) * (median_diameter_D ** 0.61)) / (discharge_Q ** -0.25)

def classify_channel_pattern(slope_S: float,
                           discharge_Q: float,
                           criterion: str = 'leopold_wolman',
                           **kwargs) -> Tuple[str, float, str]:
    """
    Classify channel pattern as meandering or braided based on selected criterion.
    
    Parameters:
    -----------
    slope_S : float
        Channel slope (dimensionless)
    discharge_Q : float  
        Bank-full discharge (m³/s)
    criterion : str
        Criterion to use: 'leopold_wolman', 'henderson', 'lane', or 'millar'
    **kwargs : dict
        Additional parameters for specific criteria
        
    Returns:
    --------
    tuple
        (pattern, threshold_slope, confidence)
        pattern: 'meandering' or 'braided'
        threshold_slope: calculated threshold
        confidence: assessment of prediction confidence
        
    Example:
    --------
    >>> pattern, threshold, confidence = classify_channel_pattern(
    ...     slope_S=0.001, discharge_Q=100.0, criterion='leopold_wolman')
    >>> print(f"Pattern: {pattern}, Threshold: {threshold:.6f}")
    """
    validate_positive(slope_S, "slope_S")
    validate_positive(discharge_Q, "discharge_Q")
    
    # Calculate threshold based on selected criterion
    if criterion == 'leopold_wolman':
        threshold = leopold_wolman_threshold(discharge_Q)
        
    elif criterion == 'henderson':
        D = kwargs.get('median_diameter_D', 1.0)  # default 1mm
        threshold = henderson_threshold(discharge_Q, D)
        
    elif criterion == 'lane':
        K = kwargs.get('K_constant', 0.0007)
        threshold = lane_criterion(discharge_Q, K)
        
    elif criterion == 'millar':
        D = kwargs.get('median_diameter_D', 0.001)  # default 1mm in meters
        phi = kwargs.get('bank_friction_angle_phi', 35.0)  # default angle
        threshold = millar_criterion(discharge_Q, D, phi)
        
    else:
        raise ValueError(f"Unknown criterion: {criterion}")
    
    # Classify pattern
    if slope_S > threshold:
        pattern = 'braided'
        ratio = slope_S / threshold
    else:
        pattern = 'meandering'
        ratio = threshold / slope_S
    
    # Assess confidence based on how close to threshold
    if ratio > 2.0:
        confidence = 'high'
    elif ratio > 1.5:
        confidence = 'moderate'
    else:
        confidence = 'low'
    
    return pattern, threshold, confidence

def stability_transition_zone(discharge_Q: float,
                            criterion: str = 'leopold_wolman',
                            factor: float = 1.5,
                            **kwargs) -> Tuple[float, float]:
    """
    Calculate the transition zone boundaries around the stability threshold.
    
    Real rivers often show intermediate behavior near the theoretical threshold.
    This function defines upper and lower bounds for a transition zone.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
    criterion : str
        Stability criterion to use
    factor : float
        Factor defining transition zone width (default=1.5)
    **kwargs : dict
        Additional parameters for specific criteria
        
    Returns:
    --------
    tuple
        (lower_bound, upper_bound) slopes defining transition zone
    """
    validate_positive(discharge_Q, "discharge_Q")
    validate_positive(factor, "factor")
    
    # Get base threshold
    if criterion == 'leopold_wolman':
        threshold = leopold_wolman_threshold(discharge_Q)
    elif criterion == 'henderson':
        D = kwargs.get('median_diameter_D', 1.0)
        threshold = henderson_threshold(discharge_Q, D)
    elif criterion == 'lane':
        K = kwargs.get('K_constant', 0.0007)
        threshold = lane_criterion(discharge_Q, K)
    elif criterion == 'millar':
        D = kwargs.get('median_diameter_D', 0.001)
        phi = kwargs.get('bank_friction_angle_phi', 35.0)
        threshold = millar_criterion(discharge_Q, D, phi)
    else:
        raise ValueError(f"Unknown criterion: {criterion}")
    
    # Define transition zone
    lower_bound = threshold / factor
    upper_bound = threshold * factor
    
    return lower_bound, upper_bound

def compare_criteria(discharge_Q: float,
                    median_diameter_D: float = 1.0,
                    bank_friction_angle_phi: float = 35.0,
                    K_constant: float = 0.0007) -> dict:
    """
    Compare threshold slopes calculated by different criteria.
    
    Parameters:
    -----------
    discharge_Q : float
        Bank-full discharge (m³/s)
    median_diameter_D : float, optional
        Median particle diameter (mm for Henderson, m for Millar)
    bank_friction_angle_phi : float, optional
        Bank friction angle (degrees)
    K_constant : float, optional
        Lane's constant
        
    Returns:
    --------
    dict
        Dictionary with threshold values for each criterion
    """
    validate_positive(discharge_Q, "discharge_Q")
    
    results = {}
    
    # Leopold-Wolman
    results['leopold_wolman'] = leopold_wolman_threshold(discharge_Q)
    
    # Henderson
    results['henderson'] = henderson_threshold(discharge_Q, median_diameter_D)
    
    # Lane
    results['lane'] = lane_criterion(discharge_Q, K_constant)
    
    # Millar (convert diameter to meters)
    D_meters = median_diameter_D / 1000.0
    results['millar'] = millar_criterion(discharge_Q, D_meters, bank_friction_angle_phi)
    
    # Add statistics
    values = list(results.values())
    results['mean'] = np.mean(values)
    results['std'] = np.std(values)
    results['range'] = max(values) - min(values)
    
    return results