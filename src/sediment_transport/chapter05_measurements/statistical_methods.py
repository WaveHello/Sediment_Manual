"""
Statistical methods for sediment transport measurements (Chapter 5, Equations 5-1 to 5-6).

This module implements statistical methods for sampling design, bias corrections,
and confidence interval calculations for sediment transport measurements.
"""

import numpy as np
from typing import Union, Tuple, Optional
from scipy import stats
from ..utils.validators import validate_positive, validate_range, validate_array_like

def areal_to_volumetric_conversion(
    areal_proportion: Union[float, np.ndarray],
    particle_diameter: Union[float, np.ndarray],
    conversion_factor: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Convert areal sample proportions to volumetric proportions.
    
    Equation 5-1: p(VW)i = Cp(A)iDix
    
    Where:
    - p(VW)i is volumetric proportion by weight for size class i
    - p(A)i is areal proportion for size class i  
    - Di is particle diameter for size class i [m]
    - C is conversion factor
    - x is diameter exponent (typically 1.0)
    
    Parameters
    ----------
    areal_proportion : float or array-like
        Areal proportion for each size class [-]
    particle_diameter : float or array-like
        Particle diameter for each size class [m]
    conversion_factor : float, optional
        Conversion factor C, default 1.0
        
    Returns
    -------
    float or ndarray
        Volumetric proportion by weight [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 5, Equation 5-1
    """
    validate_range(areal_proportion, "areal_proportion", 0.0, 1.0)
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(conversion_factor, "conversion_factor")
    
    # Equation 5-1: p(VW)i = Cp(A)iDi
    volumetric_proportion = conversion_factor * areal_proportion * particle_diameter
    
    return volumetric_proportion


def sample_size_binomial(
    expected_proportion: float,
    confidence_level: float = 0.95,
    margin_of_error: float = 0.05,
    finite_population: Optional[int] = None
) -> int:
    """
    Calculate required sample size for binomial sampling.
    
    Equation 5-2: n = (z²pq) / E²
    
    Where:
    - n is required sample size
    - z is critical z-value for confidence level
    - p is expected proportion
    - q = 1 - p
    - E is margin of error
    
    Parameters
    ----------
    expected_proportion : float
        Expected proportion of characteristic [0-1]
    confidence_level : float, optional
        Confidence level, default 0.95 (95%)
    margin_of_error : float, optional
        Acceptable margin of error, default 0.05
    finite_population : int, optional
        Population size for finite population correction
        
    Returns
    -------
    int
        Required sample size
        
    References
    ----------
    ASCE Manual 110, Chapter 5, Equation 5-2
    """
    validate_range(expected_proportion, "expected_proportion", 0.0, 1.0)
    validate_range(confidence_level, "confidence_level", 0.0, 1.0)
    validate_positive(margin_of_error, "margin_of_error")
    
    # Critical z-value for confidence level
    alpha = 1 - confidence_level
    z = stats.norm.ppf(1 - alpha/2)
    
    p = expected_proportion
    q = 1 - p
    E = margin_of_error
    
    # Equation 5-2: Basic sample size calculation
    n_infinite = (z**2 * p * q) / E**2
    
    # Finite population correction if specified
    if finite_population is not None:
        N = finite_population
        n = n_infinite / (1 + (n_infinite - 1) / N)
    else:
        n = n_infinite
        
    return int(np.ceil(n))


def multinomial_confidence_interval(
    sample_proportions: np.ndarray,
    sample_size: int,
    confidence_level: float = 0.95
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate confidence intervals for multinomial proportions.
    
    Equation 5-3: CI = pi ± zα/2 * sqrt(pi(1-pi)/n)
    
    Where:
    - pi is sample proportion for category i
    - zα/2 is critical value
    - n is sample size
    
    Parameters
    ----------
    sample_proportions : array-like
        Sample proportions for each category
    sample_size : int
        Total sample size
    confidence_level : float, optional
        Confidence level, default 0.95
        
    Returns
    -------
    tuple
        (lower_bounds, upper_bounds) for confidence intervals
        
    References
    ----------
    ASCE Manual 110, Chapter 5, Equation 5-3
    """
    proportions = validate_array_like(sample_proportions, "sample_proportions")
    validate_positive(sample_size, "sample_size")
    validate_range(confidence_level, "confidence_level", 0.0, 1.0)
    
    if not np.isclose(np.sum(proportions), 1.0, rtol=1e-3):
        raise ValueError("Sample proportions must sum to 1.0")
    
    alpha = 1 - confidence_level
    z = stats.norm.ppf(1 - alpha/2)
    
    # Equation 5-3: Confidence interval calculation
    standard_errors = np.sqrt(proportions * (1 - proportions) / sample_size)
    margin_of_error = z * standard_errors
    
    lower_bounds = proportions - margin_of_error
    upper_bounds = proportions + margin_of_error
    
    # Ensure bounds are within [0, 1]
    lower_bounds = np.maximum(lower_bounds, 0.0)
    upper_bounds = np.minimum(upper_bounds, 1.0)
    
    return lower_bounds, upper_bounds


def sampling_bias_correction(
    sample_mean: float,
    sample_variance: float,
    sample_size: int,
    bias_type: str = "finite_sample"
) -> Tuple[float, float]:
    """
    Apply bias corrections to sample statistics.
    
    Equation 5-4: Corrected variance = s² * n/(n-1)
    Equation 5-5: Bias-corrected mean for small samples
    
    Parameters
    ----------
    sample_mean : float
        Sample mean
    sample_variance : float
        Sample variance
    sample_size : int
        Sample size
    bias_type : str, optional
        Type of bias correction: 'finite_sample' or 'small_sample'
        
    Returns
    -------
    tuple
        (corrected_mean, corrected_variance)
        
    References
    ----------
    ASCE Manual 110, Chapter 5, Equations 5-4, 5-5
    """
    validate_positive(sample_size, "sample_size")
    validate_positive(sample_variance, "sample_variance")
    
    if sample_size <= 1:
        raise ValueError("Sample size must be > 1 for bias correction")
    
    corrected_mean = sample_mean  # Usually no bias correction needed for mean
    
    if bias_type == "finite_sample":
        # Equation 5-4: Bessel's correction for sample variance
        corrected_variance = sample_variance * sample_size / (sample_size - 1)
    elif bias_type == "small_sample":
        # Equation 5-5: Small sample bias correction
        correction_factor = (sample_size - 1) / sample_size
        corrected_variance = sample_variance / correction_factor
    else:
        raise ValueError("bias_type must be 'finite_sample' or 'small_sample'")
    
    return corrected_mean, corrected_variance


def bootstrap_confidence_interval(
    data: np.ndarray,
    statistic_func: callable,
    confidence_level: float = 0.95,
    n_bootstrap: int = 1000,
    random_state: Optional[int] = None
) -> Tuple[float, float, float]:
    """
    Calculate bootstrap confidence intervals for any statistic.
    
    Equation 5-6: Bootstrap percentile method for confidence intervals
    
    Parameters
    ----------
    data : array-like
        Original data sample
    statistic_func : callable
        Function to compute the statistic (e.g., np.mean, np.std)
    confidence_level : float, optional
        Confidence level, default 0.95
    n_bootstrap : int, optional
        Number of bootstrap samples, default 1000
    random_state : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    tuple
        (statistic_value, lower_bound, upper_bound)
        
    References
    ----------
    ASCE Manual 110, Chapter 5, Equation 5-6
    """
    data = validate_array_like(data, "data")
    validate_range(confidence_level, "confidence_level", 0.0, 1.0)
    validate_positive(n_bootstrap, "n_bootstrap")
    
    if random_state is not None:
        np.random.seed(random_state)
    
    n = len(data)
    bootstrap_statistics = np.zeros(n_bootstrap)
    
    # Generate bootstrap samples and compute statistic
    for i in range(n_bootstrap):
        bootstrap_sample = np.random.choice(data, size=n, replace=True)
        bootstrap_statistics[i] = statistic_func(bootstrap_sample)
    
    # Original statistic
    original_statistic = statistic_func(data)
    
    # Equation 5-6: Percentile method confidence intervals
    alpha = 1 - confidence_level
    lower_percentile = 100 * (alpha / 2)
    upper_percentile = 100 * (1 - alpha / 2)
    
    lower_bound = np.percentile(bootstrap_statistics, lower_percentile)
    upper_bound = np.percentile(bootstrap_statistics, upper_percentile)
    
    return original_statistic, lower_bound, upper_bound


def effective_sample_size(
    nominal_sample_size: int,
    correlation_coefficient: float
) -> float:
    """
    Calculate effective sample size accounting for spatial/temporal correlation.
    
    Used when samples are not independent due to spatial or temporal correlation.
    
    Parameters
    ----------
    nominal_sample_size : int
        Nominal (actual) number of samples collected
    correlation_coefficient : float
        Correlation coefficient between samples [0-1]
        
    Returns
    -------
    float
        Effective sample size
        
    Notes
    -----
    This accounts for reduced effective degrees of freedom when samples
    are spatially or temporally correlated.
    """
    validate_positive(nominal_sample_size, "nominal_sample_size")
    validate_range(correlation_coefficient, "correlation_coefficient", 0.0, 1.0)
    
    # Effective sample size with correlation
    effective_n = nominal_sample_size * (1 - correlation_coefficient) / (1 + correlation_coefficient)
    
    return max(1.0, effective_n)  # Ensure at least 1


def statistical_power_analysis(
    effect_size: float,
    sample_size: int,
    alpha: float = 0.05,
    test_type: str = "two_tailed"
) -> float:
    """
    Calculate statistical power for hypothesis testing.
    
    Parameters
    ----------
    effect_size : float
        Standardized effect size (Cohen's d)
    sample_size : int
        Sample size
    alpha : float, optional
        Type I error rate, default 0.05
    test_type : str, optional
        'one_tailed' or 'two_tailed', default 'two_tailed'
        
    Returns
    -------
    float
        Statistical power (1 - β)
    """
    validate_positive(effect_size, "effect_size")
    validate_positive(sample_size, "sample_size")
    validate_range(alpha, "alpha", 0.0, 1.0)
    
    # Critical z-value
    if test_type == "two_tailed":
        z_critical = stats.norm.ppf(1 - alpha/2)
    elif test_type == "one_tailed":
        z_critical = stats.norm.ppf(1 - alpha)
    else:
        raise ValueError("test_type must be 'one_tailed' or 'two_tailed'")
    
    # Standard error
    standard_error = 1.0 / np.sqrt(sample_size)
    
    # Power calculation
    z_beta = z_critical - effect_size / standard_error
    power = 1 - stats.norm.cdf(z_beta)
    
    return power