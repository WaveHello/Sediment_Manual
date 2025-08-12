"""
Bed material sampling methods and protocols (Chapter 5).

This module implements various bed material sampling techniques including
grid sampling, volumetric methods, photographic analysis, and bias corrections.
"""

import numpy as np
from typing import Union, Tuple, List, Dict, Optional
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY

def wolman_walk_sampling(
    grid_spacing: float,
    stream_width: float,
    sample_points_per_transect: int = 11,
    number_of_transects: Optional[int] = None
) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Design Wolman walk (pebble count) sampling grid.
    
    The Wolman walk method involves systematic sampling of surface particles
    along predetermined transects across the stream.
    
    Parameters
    ----------
    grid_spacing : float
        Distance between sampling points [m]
    stream_width : float
        Average stream width [m]
    sample_points_per_transect : int, optional
        Number of sample points per transect, default 11
    number_of_transects : int, optional
        Number of transects, if None calculated based on reach length
        
    Returns
    -------
    tuple
        (total_samples, x_coordinates, y_coordinates)
        
    References
    ----------
    ASCE Manual 110, Chapter 5: Wolman Walk Sampling Protocol
    """
    validate_positive(grid_spacing, "grid_spacing")
    validate_positive(stream_width, "stream_width")
    validate_positive(sample_points_per_transect, "sample_points_per_transect")
    
    if number_of_transects is None:
        # Default: 10 transects minimum or 1 per channel width
        number_of_transects = max(10, int(stream_width / grid_spacing))
    
    # Generate sampling coordinates
    transect_spacing = grid_spacing
    point_spacing = stream_width / (sample_points_per_transect - 1)
    
    x_coords = []
    y_coords = []
    
    for i in range(number_of_transects):
        transect_x = i * transect_spacing
        for j in range(sample_points_per_transect):
            point_y = j * point_spacing
            x_coords.append(transect_x)
            y_coords.append(point_y)
    
    total_samples = len(x_coords)
    x_coordinates = np.array(x_coords)
    y_coordinates = np.array(y_coords)
    
    return total_samples, x_coordinates, y_coordinates


def grid_sampling_design(
    sampling_area_length: float,
    sampling_area_width: float,
    target_sample_size: int,
    sampling_pattern: str = "regular"
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Design regular or random grid sampling pattern for bed material.
    
    Parameters
    ----------
    sampling_area_length : float
        Length of sampling area [m]
    sampling_area_width : float
        Width of sampling area [m]
    target_sample_size : int
        Target number of samples
    sampling_pattern : str, optional
        'regular', 'random', or 'stratified', default 'regular'
        
    Returns
    -------
    tuple
        (x_coordinates, y_coordinates, grid_spacing)
    """
    validate_positive(sampling_area_length, "sampling_area_length")
    validate_positive(sampling_area_width, "sampling_area_width")
    validate_positive(target_sample_size, "target_sample_size")
    
    area = sampling_area_length * sampling_area_width
    grid_spacing = np.sqrt(area / target_sample_size)
    
    if sampling_pattern == "regular":
        # Regular grid pattern
        n_x = int(sampling_area_length / grid_spacing)
        n_y = int(sampling_area_width / grid_spacing)
        
        x_coords = np.linspace(grid_spacing/2, sampling_area_length - grid_spacing/2, n_x)
        y_coords = np.linspace(grid_spacing/2, sampling_area_width - grid_spacing/2, n_y)
        
        x_grid, y_grid = np.meshgrid(x_coords, y_coords)
        x_coordinates = x_grid.flatten()
        y_coordinates = y_grid.flatten()
        
    elif sampling_pattern == "random":
        # Random sampling pattern
        np.random.seed(42)  # For reproducibility
        x_coordinates = np.random.uniform(0, sampling_area_length, target_sample_size)
        y_coordinates = np.random.uniform(0, sampling_area_width, target_sample_size)
        
    elif sampling_pattern == "stratified":
        # Stratified random sampling
        n_strata_x = int(np.sqrt(target_sample_size * sampling_area_length / sampling_area_width))
        n_strata_y = int(target_sample_size / n_strata_x)
        
        strata_length = sampling_area_length / n_strata_x
        strata_width = sampling_area_width / n_strata_y
        
        x_coordinates = []
        y_coordinates = []
        
        np.random.seed(42)
        for i in range(n_strata_x):
            for j in range(n_strata_y):
                x_min = i * strata_length
                x_max = (i + 1) * strata_length
                y_min = j * strata_width
                y_max = (j + 1) * strata_width
                
                x_coordinates.append(np.random.uniform(x_min, x_max))
                y_coordinates.append(np.random.uniform(y_min, y_max))
        
        x_coordinates = np.array(x_coordinates)
        y_coordinates = np.array(y_coordinates)
        
    else:
        raise ValueError("sampling_pattern must be 'regular', 'random', or 'stratified'")
    
    return x_coordinates, y_coordinates, grid_spacing


def volumetric_sampling_protocol(
    sample_volume: float,
    particle_density: float = 2650.0,
    water_content: float = 0.0,
    sieve_sizes: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Calculate volumetric sampling parameters and corrections.
    
    Parameters
    ----------
    sample_volume : float
        Volume of sediment sample [m³]
    particle_density : float, optional
        Particle density [kg/m³], default 2650.0
    water_content : float, optional
        Water content as fraction of dry weight, default 0.0
    sieve_sizes : array-like, optional
        Sieve sizes for mass distribution [mm]
        
    Returns
    -------
    dict
        Dictionary containing sampling parameters
    """
    validate_positive(sample_volume, "sample_volume")
    validate_positive(particle_density, "particle_density")
    validate_range(water_content, "water_content", 0.0, 1.0)
    
    # Calculate sample mass
    dry_mass = sample_volume * particle_density
    total_mass = dry_mass * (1 + water_content)
    
    # Minimum sample mass for various analyses
    min_mass_sieve = 0.1  # kg for sieve analysis
    min_mass_hydrometer = 0.05  # kg for hydrometer analysis
    
    results = {
        'sample_volume_m3': sample_volume,
        'dry_mass_kg': dry_mass,
        'total_mass_kg': total_mass,
        'particle_density_kg_m3': particle_density,
        'water_content': water_content,
        'adequate_for_sieve': dry_mass >= min_mass_sieve,
        'adequate_for_hydrometer': dry_mass >= min_mass_hydrometer,
        'recommended_subsample_mass_kg': min(dry_mass * 0.1, 1.0)
    }
    
    if sieve_sizes is not None:
        sieve_sizes = validate_array_like(sieve_sizes, "sieve_sizes")
        # Estimate mass retained on each sieve (simplified uniform distribution)
        mass_per_sieve = dry_mass / len(sieve_sizes)
        results['estimated_mass_per_sieve_kg'] = mass_per_sieve
        results['sieve_sizes_mm'] = sieve_sizes
    
    return results


def photographic_sampling_analysis(
    image_area: float,
    pixel_resolution: float,
    particle_count: int,
    calibration_length: float
) -> Dict[str, float]:
    """
    Calculate photographic sampling analysis parameters.
    
    Parameters
    ----------
    image_area : float
        Total area covered by photograph [m²]
    pixel_resolution : float
        Resolution in pixels per meter
    particle_count : int
        Number of particles identified in image
    calibration_length : float
        Known length for calibration [m]
        
    Returns
    -------
    dict
        Dictionary containing analysis parameters
    """
    validate_positive(image_area, "image_area")
    validate_positive(pixel_resolution, "pixel_resolution")
    validate_positive(particle_count, "particle_count")
    validate_positive(calibration_length, "calibration_length")
    
    # Calculate sampling statistics
    particle_density_per_m2 = particle_count / image_area
    average_particle_area = image_area / particle_count
    equivalent_diameter = np.sqrt(4 * average_particle_area / np.pi)
    
    # Sampling efficiency metrics
    total_pixels = image_area * pixel_resolution**2
    pixels_per_particle = total_pixels / particle_count
    
    results = {
        'image_area_m2': image_area,
        'particle_count': particle_count,
        'particle_density_per_m2': particle_density_per_m2,
        'average_particle_area_m2': average_particle_area,
        'equivalent_diameter_m': equivalent_diameter,
        'pixel_resolution_per_m': pixel_resolution,
        'pixels_per_particle': pixels_per_particle,
        'calibration_length_m': calibration_length
    }
    
    return results


def surface_oriented_sampling(
    measured_diameters: np.ndarray,
    sampling_method: str = "grid_by_number"
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Process surface-oriented sampling data with bias corrections.
    
    Surface sampling tends to over-represent larger particles due to
    their greater probability of being selected.
    
    Parameters
    ----------
    measured_diameters : array-like
        Measured particle diameters [mm]
    sampling_method : str, optional
        Sampling method: 'grid_by_number' or 'areal_sampling'
        
    Returns
    -------
    dict
        Dictionary containing processed results and bias corrections
    """
    diameters = validate_array_like(measured_diameters, "measured_diameters")
    validate_positive(diameters, "measured_diameters")
    
    n_particles = len(diameters)
    
    # Basic statistics
    d_mean = np.mean(diameters)
    d_median = np.median(diameters)
    d_std = np.std(diameters, ddof=1)
    
    # Percentiles
    d16 = np.percentile(diameters, 16)
    d50 = np.percentile(diameters, 50)  # Same as median
    d84 = np.percentile(diameters, 84)
    d90 = np.percentile(diameters, 90)
    
    # Bias correction for surface sampling
    if sampling_method == "grid_by_number":
        # Grid-by-number samples particles proportional to their number
        # Less bias than areal sampling
        bias_correction_factor = 1.0
    elif sampling_method == "areal_sampling":
        # Areal sampling over-represents larger particles
        # Apply correction based on diameter
        weights = 1.0 / diameters  # Weight inversely proportional to diameter
        weights_normalized = weights / np.sum(weights)
        
        d_mean_corrected = np.sum(diameters * weights_normalized)
        bias_correction_factor = d_mean_corrected / d_mean
    else:
        raise ValueError("sampling_method must be 'grid_by_number' or 'areal_sampling'")
    
    # Geometric statistics
    geometric_mean = np.exp(np.mean(np.log(diameters)))
    geometric_std = np.exp(np.std(np.log(diameters), ddof=1))
    
    results = {
        'n_particles': n_particles,
        'diameter_mean_mm': d_mean,
        'diameter_median_mm': d_median,
        'diameter_std_mm': d_std,
        'diameter_d16_mm': d16,
        'diameter_d50_mm': d50,
        'diameter_d84_mm': d84,
        'diameter_d90_mm': d90,
        'geometric_mean_mm': geometric_mean,
        'geometric_std': geometric_std,
        'sorting_parameter': d84 / d16,  # Geometric sorting
        'skewness': (d16 * d84 - d50**2) / (d50 * (d84 - d16)),
        'bias_correction_factor': bias_correction_factor,
        'sampling_method': sampling_method,
        'measured_diameters_mm': diameters
    }
    
    return results


def sampling_representativeness_analysis(
    sample_statistics: Dict[str, float],
    population_statistics: Optional[Dict[str, float]] = None,
    confidence_level: float = 0.95
) -> Dict[str, float]:
    """
    Analyze representativeness of bed material samples.
    
    Parameters
    ----------
    sample_statistics : dict
        Dictionary containing sample statistics
    population_statistics : dict, optional
        Known population statistics for comparison
    confidence_level : float, optional
        Confidence level for analysis, default 0.95
        
    Returns
    -------
    dict
        Dictionary containing representativeness analysis
    """
    validate_range(confidence_level, "confidence_level", 0.0, 1.0)
    
    results = {
        'confidence_level': confidence_level,
        'sample_size': sample_statistics.get('n_particles', 0)
    }
    
    # Calculate sampling uncertainty
    if 'diameter_std_mm' in sample_statistics and 'n_particles' in sample_statistics:
        n = sample_statistics['n_particles']
        std = sample_statistics['diameter_std_mm']
        mean = sample_statistics['diameter_mean_mm']
        
        # Standard error of the mean
        standard_error = std / np.sqrt(n)
        
        # Confidence interval for mean
        from scipy import stats
        alpha = 1 - confidence_level
        t_critical = stats.t.ppf(1 - alpha/2, df=n-1)
        margin_of_error = t_critical * standard_error
        
        results.update({
            'standard_error_mm': standard_error,
            'margin_of_error_mm': margin_of_error,
            'confidence_interval_lower_mm': mean - margin_of_error,
            'confidence_interval_upper_mm': mean + margin_of_error,
            'relative_error_percent': (margin_of_error / mean) * 100
        })
    
    # Compare with population if available
    if population_statistics is not None:
        pop_mean = population_statistics.get('diameter_mean_mm')
        sample_mean = sample_statistics.get('diameter_mean_mm')
        
        if pop_mean is not None and sample_mean is not None:
            bias = sample_mean - pop_mean
            relative_bias = (bias / pop_mean) * 100
            
            results.update({
                'bias_mm': bias,
                'relative_bias_percent': relative_bias,
                'population_mean_mm': pop_mean,
                'sample_mean_mm': sample_mean
            })
    
    return results


def stratified_sampling_design(
    reach_length: float,
    habitat_types: List[str],
    habitat_proportions: np.ndarray,
    total_sample_size: int
) -> Dict[str, Union[int, List[str]]]:
    """
    Design stratified sampling by habitat type.
    
    Parameters
    ----------
    reach_length : float
        Total reach length [m]
    habitat_types : list
        List of habitat type names
    habitat_proportions : array-like
        Proportion of each habitat type [0-1]
    total_sample_size : int
        Total number of samples desired
        
    Returns
    -------
    dict
        Dictionary containing stratified sampling design
    """
    validate_positive(reach_length, "reach_length")
    validate_positive(total_sample_size, "total_sample_size")
    
    proportions = validate_array_like(habitat_proportions, "habitat_proportions")
    
    if len(habitat_types) != len(proportions):
        raise ValueError("Number of habitat types must match number of proportions")
    
    if not np.isclose(np.sum(proportions), 1.0, rtol=1e-3):
        raise ValueError("Habitat proportions must sum to 1.0")
    
    # Calculate samples per stratum
    samples_per_habitat = (proportions * total_sample_size).astype(int)
    
    # Adjust for rounding errors
    total_allocated = np.sum(samples_per_habitat)
    remaining_samples = total_sample_size - total_allocated
    
    # Allocate remaining samples to largest habitats
    if remaining_samples > 0:
        largest_habitats = np.argsort(proportions)[-remaining_samples:]
        samples_per_habitat[largest_habitats] += 1
    
    results = {
        'reach_length_m': reach_length,
        'habitat_types': habitat_types,
        'habitat_proportions': proportions,
        'total_sample_size': total_sample_size,
        'samples_per_habitat': samples_per_habitat,
        'sampling_density_per_habitat': samples_per_habitat / (proportions * reach_length)
    }
    
    return results