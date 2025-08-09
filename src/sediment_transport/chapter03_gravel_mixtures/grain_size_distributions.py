"""
Grain Size Distribution Analysis (Chapter 3, Equations 3-1 to 3-12)

This module implements statistical methods for analyzing sediment grain size
distributions, including ψ scale transformations, probability distributions,
and discretized statistical formulations.
"""

import numpy as np
from typing import Union, Tuple
from ..utils.validators import validate_positive, validate_array


def psi_scale_transform(D: Union[float, np.ndarray], 
                       inverse: bool = False) -> Union[float, np.ndarray]:
    """
    ψ scale transformation for grain size representation.
    
    Equations 3-1a,b: ψ = -ln(D)/ln(2), D = 2^(-ψ)
    
    Args:
        D: Grain size in mm (if inverse=False) or ψ values (if inverse=True)
        inverse: If True, convert from ψ to D; if False, convert D to ψ
        
    Returns:
        ψ values (if inverse=False) or grain sizes in mm (if inverse=True)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-1a,b
    """
    if inverse:
        # D = 2^(-ψ) 
        return 2.0 ** (-D)
    else:
        # ψ = -ln(D)/ln(2)
        validate_positive(D, "grain_size")
        return -np.log(D) / np.log(2.0)


def percentile_size(grain_sizes: np.ndarray, 
                   cumulative_fraction: np.ndarray,
                   percentile: float) -> float:
    """
    Calculate percentile grain size from cumulative distribution.
    
    Equations 3-3, 3-4: Percentile size calculations
    
    Args:
        grain_sizes: Array of grain sizes (mm)
        cumulative_fraction: Cumulative fraction passing each size
        percentile: Percentile to calculate (0-100)
        
    Returns:
        Grain size at specified percentile (mm)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-3, 3-4
    """
    validate_array(grain_sizes)
    validate_array(cumulative_fraction)
    
    if not 0 <= percentile <= 100:
        raise ValueError("Percentile must be between 0 and 100")
        
    target_fraction = percentile / 100.0
    return np.interp(target_fraction, cumulative_fraction, grain_sizes)


def grain_size_statistics(grain_sizes: np.ndarray, 
                         fractions: np.ndarray,
                         use_geometric: bool = True) -> Tuple[float, float]:
    """
    Calculate mean and standard deviation of grain size distribution.
    
    Equations 3-5a,b: Arithmetic mean and standard deviation
    Equations 3-6a,b: Geometric mean and standard deviation  
    Equations 3-9a,b: Geometric statistics for lognormal case
    
    Args:
        grain_sizes: Array of grain sizes (mm)
        fractions: Array of mass fractions for each size class
        use_geometric: If True, calculate geometric statistics; if False, arithmetic
        
    Returns:
        Tuple of (mean, standard_deviation)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-5 to 3-9
    """
    validate_array(grain_sizes)
    validate_array(fractions)
    validate_positive(grain_sizes, "grain_sizes")
    
    if len(grain_sizes) != len(fractions):
        raise ValueError("grain_sizes and fractions must have same length")
        
    if not np.isclose(np.sum(fractions), 1.0, rtol=1e-3):
        raise ValueError("Fractions must sum to 1.0")
    
    if use_geometric:
        # Equations 3-6a,b: Geometric mean and standard deviation
        log_sizes = np.log(grain_sizes)
        mean_log = np.sum(fractions * log_sizes)
        var_log = np.sum(fractions * (log_sizes - mean_log)**2)
        
        geometric_mean = np.exp(mean_log)  # Dg
        geometric_std = np.exp(np.sqrt(var_log))  # σg
        
        return geometric_mean, geometric_std
    else:
        # Equations 3-5a,b: Arithmetic mean and standard deviation
        arithmetic_mean = np.sum(fractions * grain_sizes)
        variance = np.sum(fractions * (grain_sizes - arithmetic_mean)**2)
        arithmetic_std = np.sqrt(variance)
        
        return arithmetic_mean, arithmetic_std


def lognormal_distribution(D: Union[float, np.ndarray],
                          Dg: float, 
                          sigma_g: float,
                          return_density: bool = True) -> Union[float, np.ndarray]:
    """
    Lognormal probability distribution for grain sizes.
    
    Equations 3-7a,b: Lognormal distribution forms
    Equations 3-8a,b: Parameters for lognormal distributions
    
    Args:
        D: Grain size(s) (mm)
        Dg: Geometric mean grain size (mm)
        sigma_g: Geometric standard deviation
        return_density: If True, return probability density; if False, cumulative
        
    Returns:
        Probability density or cumulative distribution
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-7, 3-8
    """
    validate_positive(D, "grain_size")
    validate_positive(Dg, "geometric_mean")
    validate_positive(sigma_g, "geometric_std")
    
    if sigma_g <= 1.0:
        raise ValueError("Geometric standard deviation must be > 1.0")
    
    # Convert to log space
    ln_D = np.log(D)
    ln_Dg = np.log(Dg)
    ln_sigma_g = np.log(sigma_g)
    
    if return_density:
        # Equation 3-7a: Probability density function
        coefficient = 1.0 / (D * ln_sigma_g * np.sqrt(2 * np.pi))
        exponent = -0.5 * ((ln_D - ln_Dg) / ln_sigma_g)**2
        return coefficient * np.exp(exponent)
    else:
        # Equation 3-7b: Cumulative distribution function
        from scipy.stats import norm
        standardized = (ln_D - ln_Dg) / ln_sigma_g
        return norm.cdf(standardized)


def discretized_grain_size_stats(grain_sizes: np.ndarray,
                                fractions: np.ndarray) -> Tuple[float, float]:
    """
    Calculate discretized grain size statistics for size classes.
    
    Equations 3-10: Arithmetic mean grain size using linear scale
    Equations 3-12a,b,c,d: Discretized statistical formulations
    
    Args:
        grain_sizes: Array of grain sizes for size classes (mm)
        fractions: Array of mass fractions for each class
        
    Returns:
        Tuple of (discretized_mean, discretized_std)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-10, 3-12
    """
    validate_array(grain_sizes)
    validate_array(fractions)
    
    if len(grain_sizes) != len(fractions):
        raise ValueError("Arrays must have same length")
    
    # Equation 3-10: Arithmetic mean using linear scale
    mean_size = np.sum(fractions * grain_sizes)
    
    # Equation 3-12: Discretized standard deviation
    variance = np.sum(fractions * (grain_sizes - mean_size)**2)
    std_size = np.sqrt(variance)
    
    return mean_size, std_size


def grain_size_intervals(psi_min: float, 
                        psi_max: float, 
                        n_intervals: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create discretized grain size intervals using ψ scale.
    
    Equations 3-11a,b,c: Discretized grain size intervals
    
    Args:
        psi_min: Minimum ψ value
        psi_max: Maximum ψ value  
        n_intervals: Number of size intervals
        
    Returns:
        Tuple of (psi_values, grain_sizes_mm, interval_widths)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-11
    """
    if n_intervals <= 0:
        raise ValueError("Number of intervals must be positive")
    
    # Create ψ intervals
    psi_values = np.linspace(psi_min, psi_max, n_intervals + 1)
    
    # Convert to grain sizes using Equation 3-1b
    grain_sizes = psi_scale_transform(psi_values, inverse=True)
    
    # Calculate interval widths in ψ space
    delta_psi = (psi_max - psi_min) / n_intervals
    interval_widths = np.full(n_intervals, delta_psi)
    
    return psi_values, grain_sizes, interval_widths


def characteristic_grain_size(Dg: float, sigma_g: float) -> float:
    """
    Calculate characteristic grain size parameter.
    
    Equation 3-29: Dσ = Dg * σg
    
    Args:
        Dg: Geometric mean grain size (mm)
        sigma_g: Geometric standard deviation
        
    Returns:
        Characteristic grain size Dσ (mm)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-29
    """
    validate_positive(Dg, "geometric_mean")
    validate_positive(sigma_g, "geometric_std")
    
    return Dg * sigma_g