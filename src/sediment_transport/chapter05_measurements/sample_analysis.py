"""
Sample analysis methods and procedures (Chapter 5).

This module implements sediment sample analysis techniques including
sieve analysis, particle size distributions, concentration calculations,
and quality control procedures.
"""

import numpy as np
from typing import Union, Tuple, List, Dict, Optional
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

def sieve_analysis_processing(
    sieve_sizes: np.ndarray,
    mass_retained: np.ndarray,
    total_sample_mass: float,
    wash_loss: float = 0.0
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Process sieve analysis data to calculate particle size distribution.
    
    Parameters
    ----------
    sieve_sizes : array-like
        Sieve opening sizes [mm], in descending order
    mass_retained : array-like
        Mass retained on each sieve [g]
    total_sample_mass : float
        Total dry mass of original sample [g]
    wash_loss : float, optional
        Mass lost during washing [g], default 0.0
        
    Returns
    -------
    dict
        Dictionary containing sieve analysis results
        
    References
    ----------
    ASCE Manual 110, Chapter 5: Sediment Sample Analysis
    """
    sieves = validate_array_like(sieve_sizes, "sieve_sizes")
    masses = validate_array_like(mass_retained, "mass_retained")
    validate_positive(total_sample_mass, "total_sample_mass")
    validate_positive(sieves, "sieve_sizes")
    validate_positive(masses, "mass_retained")
    
    if len(sieves) != len(masses):
        raise ValueError("Number of sieve sizes must match number of mass measurements")
    
    # Check that sieves are in descending order
    if not np.all(np.diff(sieves) <= 0):
        # Sort in descending order
        sort_indices = np.argsort(sieves)[::-1]
        sieves = sieves[sort_indices]
        masses = masses[sort_indices]
    
    # Calculate cumulative mass retained
    cumulative_mass = np.cumsum(masses)
    
    # Account for wash loss (fines that passed through finest sieve)
    total_recovered = np.sum(masses) + wash_loss
    
    # Calculate percent retained and percent passing
    percent_retained = (masses / total_sample_mass) * 100
    percent_passing = ((total_sample_mass - cumulative_mass) / total_sample_mass) * 100
    
    # Ensure percent passing starts at 100 and ends at wash_loss percentage
    percent_passing = np.maximum(percent_passing, (wash_loss / total_sample_mass) * 100)
    
    # Calculate size statistics
    # D10, D50, D90 (particle sizes at 10%, 50%, 90% passing)
    d10 = np.interp(10.0, percent_passing[::-1], sieves[::-1])
    d50 = np.interp(50.0, percent_passing[::-1], sieves[::-1])
    d90 = np.interp(90.0, percent_passing[::-1], sieves[::-1])
    
    # Additional percentiles
    d16 = np.interp(16.0, percent_passing[::-1], sieves[::-1])
    d84 = np.interp(84.0, percent_passing[::-1], sieves[::-1])
    
    # Uniformity coefficient and curvature coefficient
    uniformity_coefficient = d60 = np.interp(60.0, percent_passing[::-1], sieves[::-1])
    d30 = np.interp(30.0, percent_passing[::-1], sieves[::-1])
    uniformity_coefficient = d60 / d10 if d10 > 0 else np.inf
    curvature_coefficient = (d30**2) / (d60 * d10) if (d60 > 0 and d10 > 0) else 0.0
    
    # Sorting parameters
    geometric_mean = np.sqrt(d16 * d84)
    geometric_sorting = np.sqrt(d84 / d16)
    
    # Folk & Ward (1957) parameters in phi units
    sieves_phi = -np.log2(sieves)  # Convert mm to phi scale
    d16_phi = -np.log2(d16)
    d50_phi = -np.log2(d50) 
    d84_phi = -np.log2(d84)
    
    folk_mean_phi = (d16_phi + d50_phi + d84_phi) / 3.0
    folk_sorting_phi = (d84_phi - d16_phi) / 4.0
    
    # Mass balance check
    mass_balance_error = abs(total_recovered - total_sample_mass) / total_sample_mass * 100
    
    results = {
        'sieve_sizes_mm': sieves,
        'mass_retained_g': masses,
        'cumulative_mass_g': cumulative_mass,
        'percent_retained': percent_retained,
        'percent_passing': percent_passing,
        'total_sample_mass_g': total_sample_mass,
        'total_recovered_mass_g': total_recovered,
        'wash_loss_g': wash_loss,
        'mass_balance_error_percent': mass_balance_error,
        'd10_mm': d10,
        'd16_mm': d16,
        'd30_mm': d30,
        'd50_mm': d50,
        'd60_mm': d60,
        'd84_mm': d84,
        'd90_mm': d90,
        'uniformity_coefficient': uniformity_coefficient,
        'curvature_coefficient': curvature_coefficient,
        'geometric_mean_mm': geometric_mean,
        'geometric_sorting': geometric_sorting,
        'folk_mean_phi': folk_mean_phi,
        'folk_sorting_phi': folk_sorting_phi
    }
    
    return results


def particle_size_distribution(
    particle_diameters: np.ndarray,
    weights: Optional[np.ndarray] = None,
    bin_method: str = "equal_log",
    n_bins: int = 20
) -> Dict[str, Union[float, np.ndarray, str]]:
    """
    Calculate particle size distribution from individual particle measurements.
    
    Parameters
    ----------
    particle_diameters : array-like
        Individual particle diameters [mm]
    weights : array-like, optional
        Weight of each particle, if None assume equal weights
    bin_method : str, optional
        Binning method: 'equal_log', 'equal_linear', 'standard_sieves'
    n_bins : int, optional
        Number of bins for distribution, default 20
        
    Returns
    -------
    dict
        Dictionary containing particle size distribution
    """
    diameters = validate_array_like(particle_diameters, "particle_diameters")
    validate_positive(diameters, "particle_diameters")
    validate_positive(n_bins, "n_bins")
    
    if weights is not None:
        weights = validate_array_like(weights, "weights")
        validate_positive(weights, "weights")
        if len(weights) != len(diameters):
            raise ValueError("Weights array must have same length as diameters")
    else:
        weights = np.ones(len(diameters))
    
    # Create bins
    if bin_method == "equal_log":
        # Equal intervals in log space
        log_min = np.log10(np.min(diameters))
        log_max = np.log10(np.max(diameters))
        bin_edges = np.logspace(log_min, log_max, n_bins + 1)
    elif bin_method == "equal_linear":
        # Equal intervals in linear space
        bin_edges = np.linspace(np.min(diameters), np.max(diameters), n_bins + 1)
    elif bin_method == "standard_sieves":
        # Standard sieve sizes
        standard_sieves = np.array([
            0.063, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0
        ])
        # Use only sieves within data range
        min_d, max_d = np.min(diameters), np.max(diameters)
        bin_edges = standard_sieves[(standard_sieves >= min_d * 0.5) & 
                                  (standard_sieves <= max_d * 2.0)]
        n_bins = len(bin_edges) - 1
    else:
        raise ValueError("bin_method must be 'equal_log', 'equal_linear', or 'standard_sieves'")
    
    # Calculate histogram
    hist_counts, _ = np.histogram(diameters, bins=bin_edges, weights=weights)
    
    # Bin centers (geometric mean for log bins, arithmetic mean for linear)
    if bin_method == "equal_log" or bin_method == "standard_sieves":
        bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    else:
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
    # Calculate percentages
    total_weight = np.sum(hist_counts)
    percent_in_bin = (hist_counts / total_weight) * 100
    
    # Cumulative distribution
    cumulative_percent = np.cumsum(percent_in_bin)
    
    # Calculate percentiles
    percentiles = [10, 16, 25, 50, 75, 84, 90, 95]
    percentile_sizes = np.interp(percentiles, cumulative_percent, bin_centers)
    
    # Statistical measures
    mean_diameter = np.average(diameters, weights=weights)
    median_diameter = percentile_sizes[3]  # D50
    
    # Weighted standard deviation
    weighted_variance = np.average((diameters - mean_diameter)**2, weights=weights)
    std_diameter = np.sqrt(weighted_variance)
    
    # Skewness and kurtosis (approximated from percentiles)
    d16, d50, d84 = percentile_sizes[1], percentile_sizes[3], percentile_sizes[5]
    skewness_approx = (d16 + d84 - 2*d50) / (d84 - d16)
    
    results = {
        'bin_method': bin_method,
        'n_bins': n_bins,
        'bin_edges_mm': bin_edges,
        'bin_centers_mm': bin_centers,
        'bin_counts': hist_counts,
        'percent_in_bin': percent_in_bin,
        'cumulative_percent': cumulative_percent,
        'total_particles': len(diameters),
        'total_weight': total_weight,
        'mean_diameter_mm': mean_diameter,
        'median_diameter_mm': median_diameter,
        'std_diameter_mm': std_diameter,
        'skewness_approx': skewness_approx,
        'percentile_values': dict(zip([f'd{p}' for p in percentiles], percentile_sizes)),
        'sorting_parameter': d84 / d16,
        'particle_diameters_mm': diameters,
        'particle_weights': weights
    }
    
    return results


def sediment_concentration_calculation(
    sample_volume: float,
    evaporation_dish_mass: float,
    dish_plus_sediment_mass: float,
    sample_temperature: float = 20.0,
    evaporation_temperature: float = 105.0
) -> Dict[str, float]:
    """
    Calculate sediment concentration from evaporation method.
    
    Parameters
    ----------
    sample_volume : float
        Volume of water sample [mL]
    evaporation_dish_mass : float
        Mass of empty evaporation dish [g]
    dish_plus_sediment_mass : float
        Mass of dish plus dried sediment [g]
    sample_temperature : float, optional
        Temperature of sample when collected [°C], default 20.0
    evaporation_temperature : float, optional
        Temperature used for evaporation [°C], default 105.0
        
    Returns
    -------
    dict
        Dictionary containing concentration calculations
    """
    validate_positive(sample_volume, "sample_volume")
    validate_positive(evaporation_dish_mass, "evaporation_dish_mass")
    validate_positive(dish_plus_sediment_mass, "dish_plus_sediment_mass")
    
    if dish_plus_sediment_mass <= evaporation_dish_mass:
        raise ValueError("Dish plus sediment mass must be greater than empty dish mass")
    
    # Calculate sediment mass
    sediment_mass_g = dish_plus_sediment_mass - evaporation_dish_mass
    
    # Convert to concentration units
    concentration_mg_L = (sediment_mass_g * 1000.0) / (sample_volume / 1000.0)  # mg/L
    concentration_g_m3 = concentration_mg_L  # g/m³ (equivalent)
    concentration_kg_m3 = concentration_mg_L / 1000.0  # kg/m³
    
    # Parts per million by weight (approximately equal to mg/L for dilute solutions)
    concentration_ppm = concentration_mg_L
    
    # Temperature corrections (if needed)
    # Water density correction
    water_density_sample = 1000.0 * (1 - 6.8e-5 * (sample_temperature - 20.0))
    water_density_reference = 1000.0
    
    # Concentration corrected to standard temperature
    concentration_corrected = concentration_mg_L * (water_density_sample / water_density_reference)
    
    # Calculate detection limit and precision estimates
    balance_precision = 0.0001  # g (typical analytical balance)
    volume_precision = sample_volume * 0.01  # 1% volume uncertainty
    
    # Uncertainty propagation (simplified)
    mass_uncertainty = balance_precision  # g
    volume_uncertainty = volume_precision / 1000.0  # L
    
    concentration_uncertainty = concentration_mg_L * np.sqrt(
        (mass_uncertainty / sediment_mass_g)**2 + 
        (volume_uncertainty / (sample_volume / 1000.0))**2
    )
    
    # Method detection limit (3x standard deviation of blanks, estimated)
    estimated_detection_limit = (3 * balance_precision * 1000.0) / (sample_volume / 1000.0)
    
    results = {
        'sample_volume_mL': sample_volume,
        'evaporation_dish_mass_g': evaporation_dish_mass,
        'dish_plus_sediment_mass_g': dish_plus_sediment_mass,
        'sediment_mass_g': sediment_mass_g,
        'concentration_mg_L': concentration_mg_L,
        'concentration_g_m3': concentration_g_m3,
        'concentration_kg_m3': concentration_kg_m3,
        'concentration_ppm': concentration_ppm,
        'concentration_corrected_mg_L': concentration_corrected,
        'sample_temperature_C': sample_temperature,
        'evaporation_temperature_C': evaporation_temperature,
        'water_density_correction_factor': water_density_sample / water_density_reference,
        'concentration_uncertainty_mg_L': concentration_uncertainty,
        'relative_uncertainty_percent': (concentration_uncertainty / concentration_mg_L) * 100,
        'estimated_detection_limit_mg_L': estimated_detection_limit,
        'above_detection_limit': concentration_mg_L > estimated_detection_limit
    }
    
    return results


def quality_control_analysis(
    sample_data: Dict[str, np.ndarray],
    control_limits: Optional[Dict[str, Tuple[float, float]]] = None,
    outlier_method: str = "iqr"
) -> Dict[str, Union[np.ndarray, List[int], Dict]]:
    """
    Perform quality control analysis on sediment sample data.
    
    Parameters
    ----------
    sample_data : dict
        Dictionary containing sample measurements with keys as parameter names
        and values as arrays of measurements
    control_limits : dict, optional
        Dictionary of (lower, upper) control limits for each parameter
    outlier_method : str, optional
        Method for outlier detection: 'iqr', 'zscore', 'modified_zscore'
        
    Returns
    -------
    dict
        Dictionary containing quality control results
    """
    if outlier_method not in ['iqr', 'zscore', 'modified_zscore']:
        raise ValueError("outlier_method must be 'iqr', 'zscore', or 'modified_zscore'")
    
    qc_results = {}
    
    for parameter, data in sample_data.items():
        data_array = validate_array_like(data, f"sample_data[{parameter}]")
        
        # Basic statistics
        mean_val = np.mean(data_array)
        median_val = np.median(data_array)
        std_val = np.std(data_array, ddof=1)
        
        # Outlier detection
        outlier_indices = []
        
        if outlier_method == "iqr":
            # Interquartile range method
            q1 = np.percentile(data_array, 25)
            q3 = np.percentile(data_array, 75)
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            outlier_indices = np.where((data_array < lower_bound) | (data_array > upper_bound))[0]
            
        elif outlier_method == "zscore":
            # Z-score method (typically |z| > 3)
            z_scores = np.abs((data_array - mean_val) / std_val)
            outlier_indices = np.where(z_scores > 3.0)[0]
            
        elif outlier_method == "modified_zscore":
            # Modified Z-score using median absolute deviation
            median_abs_dev = np.median(np.abs(data_array - median_val))
            modified_z_scores = 0.6745 * (data_array - median_val) / median_abs_dev
            outlier_indices = np.where(np.abs(modified_z_scores) > 3.5)[0]
        
        # Control limit violations
        control_violations = []
        if control_limits is not None and parameter in control_limits:
            lower_limit, upper_limit = control_limits[parameter]
            control_violations = np.where((data_array < lower_limit) | (data_array > upper_limit))[0]
        
        # Data quality flags
        quality_flags = {
            'outliers_detected': len(outlier_indices) > 0,
            'control_violations': len(control_violations) > 0,
            'high_variability': (std_val / mean_val) > 0.5 if mean_val > 0 else False,
            'small_sample_size': len(data_array) < 10,
            'non_normal_distribution': abs(mean_val - median_val) > 0.5 * std_val
        }
        
        # Overall quality assessment
        quality_issues = sum(quality_flags.values())
        if quality_issues == 0:
            quality_rating = "excellent"
        elif quality_issues <= 2:
            quality_rating = "good"
        elif quality_issues <= 3:
            quality_rating = "acceptable"
        else:
            quality_rating = "poor"
        
        qc_results[parameter] = {
            'data': data_array,
            'n_samples': len(data_array),
            'mean': mean_val,
            'median': median_val,
            'std': std_val,
            'coefficient_variation': (std_val / mean_val) * 100 if mean_val > 0 else np.inf,
            'outlier_indices': outlier_indices.tolist(),
            'n_outliers': len(outlier_indices),
            'outlier_values': data_array[outlier_indices].tolist(),
            'control_violation_indices': control_violations.tolist(),
            'n_control_violations': len(control_violations),
            'quality_flags': quality_flags,
            'quality_rating': quality_rating,
            'outlier_method': outlier_method
        }
    
    # Overall dataset quality
    parameter_ratings = [results['quality_rating'] for results in qc_results.values()]
    rating_scores = {'excellent': 4, 'good': 3, 'acceptable': 2, 'poor': 1}
    avg_score = np.mean([rating_scores[rating] for rating in parameter_ratings])
    
    if avg_score >= 3.5:
        overall_quality = "excellent"
    elif avg_score >= 2.5:
        overall_quality = "good"
    elif avg_score >= 1.5:
        overall_quality = "acceptable"
    else:
        overall_quality = "poor"
    
    qc_results['overall_assessment'] = {
        'quality_rating': overall_quality,
        'average_score': avg_score,
        'n_parameters': len(sample_data),
        'total_outliers': sum(len(results['outlier_indices']) for results in qc_results.values() if 'outlier_indices' in results),
        'total_control_violations': sum(len(results['control_violation_indices']) for results in qc_results.values() if 'control_violation_indices' in results),
        'recommendations': []
    }
    
    # Generate recommendations
    recommendations = qc_results['overall_assessment']['recommendations']
    
    total_outliers = qc_results['overall_assessment']['total_outliers']
    if total_outliers > 0:
        recommendations.append(f"Investigate {total_outliers} outlier(s) - check for measurement errors or unusual conditions")
    
    if qc_results['overall_assessment']['total_control_violations'] > 0:
        recommendations.append("Address control limit violations - check sampling and analysis procedures")
    
    high_variability_params = [param for param, results in qc_results.items() 
                              if isinstance(results, dict) and results.get('quality_flags', {}).get('high_variability', False)]
    if high_variability_params:
        recommendations.append(f"High variability in {high_variability_params} - consider increasing sample size or improving sampling consistency")
    
    if overall_quality == "poor":
        recommendations.append("Overall data quality is poor - conduct comprehensive review of sampling and analysis procedures")
    
    return qc_results


def grain_size_statistics(
    particle_size_data: Dict[str, Union[float, np.ndarray]],
    method: str = "folk_ward"
) -> Dict[str, float]:
    """
    Calculate comprehensive grain size statistics.
    
    Parameters
    ----------
    particle_size_data : dict
        Dictionary containing size percentiles (d10, d16, d25, d50, d75, d84, d90)
    method : str, optional
        Statistical method: 'folk_ward', 'geometric_moments', or 'standard_moments'
        
    Returns
    -------
    dict
        Dictionary containing grain size statistics
    """
    if method == "folk_ward":
        # Folk & Ward (1957) method using phi scale
        d16 = particle_size_data.get('d16', 0)
        d25 = particle_size_data.get('d25', 0)
        d50 = particle_size_data.get('d50', 0)
        d75 = particle_size_data.get('d75', 0)
        d84 = particle_size_data.get('d84', 0)
        d90 = particle_size_data.get('d90', 0)
        d95 = particle_size_data.get('d95', 0)
        
        # Convert to phi scale
        phi16 = -np.log2(d16) if d16 > 0 else 0
        phi25 = -np.log2(d25) if d25 > 0 else 0
        phi50 = -np.log2(d50) if d50 > 0 else 0
        phi75 = -np.log2(d75) if d75 > 0 else 0
        phi84 = -np.log2(d84) if d84 > 0 else 0
        phi90 = -np.log2(d90) if d90 > 0 else 0
        phi95 = -np.log2(d95) if d95 > 0 else 0
        
        # Folk & Ward statistics
        mean_phi = (phi16 + phi50 + phi84) / 3.0
        sorting_phi = (phi84 - phi16) / 4.0 + (phi95 - phi25) / 6.6
        skewness_phi = ((phi16 + phi84 - 2*phi50) / (2*(phi84 - phi16))) + ((phi25 + phi75 - 2*phi50) / (2*(phi75 - phi25)))
        kurtosis_phi = (phi95 - phi25) / (2.44 * (phi75 - phi25))
        
        # Convert back to mm scale
        mean_mm = 2**(-mean_phi)
        
        # Sorting classification
        if sorting_phi < 0.35:
            sorting_class = "very well sorted"
        elif sorting_phi < 0.50:
            sorting_class = "well sorted"
        elif sorting_phi < 0.71:
            sorting_class = "moderately well sorted"
        elif sorting_phi < 1.00:
            sorting_class = "moderately sorted"
        elif sorting_phi < 2.00:
            sorting_class = "poorly sorted"
        elif sorting_phi < 4.00:
            sorting_class = "very poorly sorted"
        else:
            sorting_class = "extremely poorly sorted"
        
        # Skewness classification
        if skewness_phi < -0.3:
            skewness_class = "very coarse skewed"
        elif skewness_phi < -0.1:
            skewness_class = "coarse skewed"
        elif skewness_phi < 0.1:
            skewness_class = "nearly symmetrical"
        elif skewness_phi < 0.3:
            skewness_class = "fine skewed"
        else:
            skewness_class = "very fine skewed"
        
        results = {
            'method': method,
            'mean_phi': mean_phi,
            'mean_mm': mean_mm,
            'sorting_phi': sorting_phi,
            'skewness_phi': skewness_phi,
            'kurtosis_phi': kurtosis_phi,
            'sorting_class': sorting_class,
            'skewness_class': skewness_class,
            'geometric_mean_mm': np.sqrt(d16 * d84),
            'geometric_sorting': np.sqrt(d84 / d16) if d16 > 0 else np.inf,
            'uniformity_coefficient': particle_size_data.get('d60', 0) / particle_size_data.get('d10', 1) if particle_size_data.get('d10', 0) > 0 else np.inf
        }
        
    else:
        raise ValueError("Only 'folk_ward' method is currently implemented")
    
    return results