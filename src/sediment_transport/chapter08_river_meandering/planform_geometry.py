"""
Meander Planform Geometry

This module implements relationships for meander planform characteristics
including wavelength, amplitude, and radius of curvature based on channel width.

Equations implemented:
    - 8-5: Wavelength-width relationship (Leopold & Wolman)
    - 8-6: Amplitude-width relationship (Leopold & Wolman)
    - 8-7: Wavelength-radius relationship (Leopold & Wolman)

References:
    ASCE Manual 110, Chapter 8, Section 8.2.2
    Leopold & Wolman (1960), Leopold et al. (1964), Zeller (1967)
"""

import numpy as np
from typing import Union, Tuple, Dict
from ..utils.validators import validate_positive, validate_range

def meander_wavelength(channel_width_b: float) -> float:
    """
    Calculate meander wavelength from channel width using Leopold-Wolman relationship.
    
    This relationship is consistent across a wide range of stream sizes,
    from furrow meanders to large rivers and even the Gulf Stream.
    
    Parameters:
    -----------
    channel_width_b : float
        Bank-full channel width (m)
        
    Returns:
    --------
    float
        Meander wavelength (m)
        
    Equation: 8-5
        L = 11.0 * b^1.01
        
    Reference: Leopold & Wolman (1960)
    """
    validate_positive(channel_width_b, "channel_width_b")
    
    return 11.0 * (channel_width_b ** 1.01)

def meander_amplitude(channel_width_b: float) -> float:
    """
    Calculate meander amplitude from channel width using Leopold-Wolman relationship.
    
    Parameters:
    -----------
    channel_width_b : float
        Bank-full channel width (m)
        
    Returns:
    --------
    float
        Meander amplitude (m)
        
    Equation: 8-6
        A = 3.0 * b^1.1
        
    Reference: Leopold & Wolman (1960)
    """
    validate_positive(channel_width_b, "channel_width_b")
    
    return 3.0 * (channel_width_b ** 1.1)

def wavelength_radius_relationship(radius_curvature_rc: float) -> float:
    """
    Calculate meander wavelength from minimum radius of curvature.
    
    Parameters:
    -----------
    radius_curvature_rc : float
        Minimum radius of curvature (m)
        
    Returns:
    --------
    float
        Meander wavelength (m)
        
    Equation: 8-7
        λ = 4.6 * rc^0.98
        
    Reference: Leopold & Wolman (1960)
    """
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    
    return 4.6 * (radius_curvature_rc ** 0.98)

def minimum_radius_curvature(channel_width_b: float) -> float:
    """
    Calculate minimum radius of curvature from channel width.
    
    Derived from the Leopold-Wolman relationships, this gives rc ≈ 2.4*b.
    
    Parameters:
    -----------
    channel_width_b : float
        Bank-full channel width (m)
        
    Returns:
    --------
    float
        Minimum radius of curvature (m)
        
    Note: Derived from combining equations 8-5 and 8-7
    """
    validate_positive(channel_width_b, "channel_width_b")
    
    # From L = 11.0*b^1.01 and λ = 4.6*rc^0.98
    # Solving for rc when L = λ
    wavelength = meander_wavelength(channel_width_b)
    rc = (wavelength / 4.6) ** (1.0 / 0.98)
    
    return rc

def planform_characteristics(channel_width_b: float) -> Dict[str, float]:
    """
    Calculate all meander planform characteristics from channel width.
    
    Parameters:
    -----------
    channel_width_b : float
        Bank-full channel width (m)
        
    Returns:
    --------
    dict
        Dictionary containing:
        - wavelength: meander wavelength (m)
        - amplitude: meander amplitude (m)
        - radius_curvature: minimum radius of curvature (m)
        - width_wavelength_ratio: b/L ratio
        - width_amplitude_ratio: b/A ratio
        - width_radius_ratio: b/rc ratio
        
    Example:
    --------
    >>> props = planform_characteristics(50.0)
    >>> print(f"Wavelength: {props['wavelength']:.1f} m")
    >>> print(f"Amplitude: {props['amplitude']:.1f} m")
    """
    validate_positive(channel_width_b, "channel_width_b")
    
    # Calculate primary characteristics
    wavelength = meander_wavelength(channel_width_b)
    amplitude = meander_amplitude(channel_width_b)
    radius_curvature = minimum_radius_curvature(channel_width_b)
    
    # Calculate dimensionless ratios
    width_wavelength_ratio = channel_width_b / wavelength
    width_amplitude_ratio = channel_width_b / amplitude
    width_radius_ratio = channel_width_b / radius_curvature
    
    return {
        'wavelength': wavelength,
        'amplitude': amplitude,
        'radius_curvature': radius_curvature,
        'width_wavelength_ratio': width_wavelength_ratio,
        'width_amplitude_ratio': width_amplitude_ratio,
        'width_radius_ratio': width_radius_ratio
    }

def sinuosity_calculation(channel_length: float, valley_length: float) -> float:
    """
    Calculate channel sinuosity as the ratio of channel length to valley length.
    
    Parameters:
    -----------
    channel_length : float
        Length along the channel centerline (m)
    valley_length : float
        Straight-line valley length (m)
        
    Returns:
    --------
    float
        Sinuosity ratio (dimensionless, ≥ 1.0)
        
    Note: Sinuosity = 1.0 for straight channels, π ≈ 3.14 for circular loops
    """
    validate_positive(channel_length, "channel_length")
    validate_positive(valley_length, "valley_length")
    
    sinuosity = channel_length / valley_length
    
    if sinuosity < 1.0:
        raise ValueError("Channel length cannot be less than valley length")
    
    return sinuosity

def sinuosity_from_planform(amplitude: float, wavelength: float) -> float:
    """
    Estimate sinuosity from meander amplitude and wavelength.
    
    Assumes sinusoidal planform geometry.
    
    Parameters:
    -----------
    amplitude : float
        Meander amplitude (m)
    wavelength : float
        Meander wavelength (m)
        
    Returns:
    --------
    float
        Estimated sinuosity
    """
    validate_positive(amplitude, "amplitude")
    validate_positive(wavelength, "wavelength")
    
    # For sinusoidal curves: s ≈ 1 + (2π²A²)/(3λ²) for small A/λ
    # More accurate formula using elliptic integral approximation
    k = 2 * np.pi * amplitude / wavelength  # wave steepness parameter
    
    if k < 0.1:  # Small amplitude approximation
        sinuosity = 1.0 + (k**2) / 6.0
    else:  # More accurate approximation
        sinuosity = 1.0 + (k**2) * (1.0/6.0 + (3.0/40.0)*k**2)
    
    return sinuosity

def channel_slope_from_sinuosity(valley_slope: float, sinuosity: float) -> float:
    """
    Calculate channel slope from valley slope and sinuosity.
    
    As sinuosity increases, channel slope decreases below valley slope.
    
    Parameters:
    -----------
    valley_slope : float
        Valley slope (dimensionless)
    sinuosity : float
        Channel sinuosity (≥ 1.0)
        
    Returns:
    --------
    float
        Channel slope (dimensionless)
    """
    validate_positive(valley_slope, "valley_slope")
    validate_range(sinuosity, 1.0, 10.0, "sinuosity")
    
    return valley_slope / sinuosity

def validate_planform_consistency(wavelength: float,
                                amplitude: float,
                                channel_width: float,
                                tolerance: float = 0.2) -> Dict[str, Union[bool, str]]:
    """
    Check if measured planform characteristics are consistent with Leopold-Wolman relationships.
    
    Parameters:
    -----------
    wavelength : float
        Measured wavelength (m)
    amplitude : float
        Measured amplitude (m)
    channel_width : float
        Channel width (m)
    tolerance : float, optional
        Acceptable relative error (default=0.2 for ±20%)
        
    Returns:
    --------
    dict
        Validation results with consistency checks
    """
    validate_positive(wavelength, "wavelength")
    validate_positive(amplitude, "amplitude")
    validate_positive(channel_width, "channel_width")
    validate_range(tolerance, 0.01, 1.0, "tolerance")
    
    # Calculate expected values
    expected_wavelength = meander_wavelength(channel_width)
    expected_amplitude = meander_amplitude(channel_width)
    
    # Calculate relative errors
    wavelength_error = abs(wavelength - expected_wavelength) / expected_wavelength
    amplitude_error = abs(amplitude - expected_amplitude) / expected_amplitude
    
    # Check consistency
    wavelength_consistent = wavelength_error <= tolerance
    amplitude_consistent = amplitude_error <= tolerance
    overall_consistent = wavelength_consistent and amplitude_consistent
    
    results = {
        'overall_consistent': overall_consistent,
        'wavelength_consistent': wavelength_consistent,
        'amplitude_consistent': amplitude_consistent,
        'wavelength_error': wavelength_error,
        'amplitude_error': amplitude_error,
        'expected_wavelength': expected_wavelength,
        'expected_amplitude': expected_amplitude,
        'status': 'PASS' if overall_consistent else 'FAIL'
    }
    
    return results

def scale_planform_characteristics(reference_width: float,
                                 target_width: float,
                                 reference_characteristics: Dict[str, float]) -> Dict[str, float]:
    """
    Scale planform characteristics from a reference width to a target width.
    
    Useful for scaling laboratory results or applying relationships across different scales.
    
    Parameters:
    -----------
    reference_width : float
        Reference channel width (m)
    target_width : float
        Target channel width (m)
    reference_characteristics : dict
        Reference planform characteristics
        
    Returns:
    --------
    dict
        Scaled planform characteristics
    """
    validate_positive(reference_width, "reference_width")
    validate_positive(target_width, "target_width")
    
    scale_factor = target_width / reference_width
    
    # Scale according to Leopold-Wolman exponents
    scaled = {
        'width': target_width,
        'wavelength': reference_characteristics.get('wavelength', 0) * (scale_factor ** 1.01),
        'amplitude': reference_characteristics.get('amplitude', 0) * (scale_factor ** 1.1),
        'radius_curvature': reference_characteristics.get('radius_curvature', 0) * scale_factor,
        'scale_factor': scale_factor
    }
    
    return scaled