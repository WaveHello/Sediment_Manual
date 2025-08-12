"""
Suspended sediment sampling methods and protocols (Chapter 5).

This module implements suspended sediment sampling techniques including
depth-integrated sampling, point sampling, automatic samplers, and
isokinetic sampling velocity calculations.
"""

import numpy as np
from typing import Union, Tuple, List, Dict, Optional
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

def depth_integrated_sampling(
    water_depth: float,
    flow_velocity: float,
    sampling_rate: float,
    sampler_nozzle_diameter: float,
    transit_rate: Optional[float] = None
) -> Dict[str, float]:
    """
    Calculate depth-integrated sampling parameters.
    
    Depth-integrated sampling collects a representative sample throughout
    the entire water column depth at a given vertical.
    
    Parameters
    ----------
    water_depth : float
        Total water depth [m]
    flow_velocity : float
        Mean flow velocity [m/s]
    sampling_rate : float
        Sampling rate [L/s]
    sampler_nozzle_diameter : float
        Sampler nozzle diameter [mm]
    transit_rate : float, optional
        Sampler transit rate through water column [m/s]
        
    Returns
    -------
    dict
        Dictionary containing sampling parameters
        
    References
    ----------
    ASCE Manual 110, Chapter 5: Suspended Sediment Sampling
    """
    validate_positive(water_depth, "water_depth")
    validate_positive(flow_velocity, "flow_velocity")
    validate_positive(sampling_rate, "sampling_rate")
    validate_positive(sampler_nozzle_diameter, "sampler_nozzle_diameter")
    
    # Convert nozzle diameter to meters
    nozzle_diameter_m = sampler_nozzle_diameter / 1000.0
    nozzle_area = np.pi * (nozzle_diameter_m / 2.0)**2
    
    # Calculate sampling velocity through nozzle
    sampling_velocity = (sampling_rate / 1000.0) / nozzle_area  # m/s
    
    # Isokinetic sampling condition
    isokinetic_ratio = sampling_velocity / flow_velocity
    
    # Calculate transit rate if not provided
    if transit_rate is None:
        # Recommended transit rate for uniform sampling
        transit_rate = water_depth / 60.0  # Complete transit in 60 seconds
    
    # Total sampling time for round trip
    round_trip_time = 2 * water_depth / transit_rate
    
    # Sample volume collected
    sample_volume = sampling_rate * round_trip_time / 1000.0  # m³
    
    # Water volume sampled through nozzle
    water_volume_sampled = nozzle_area * flow_velocity * round_trip_time
    
    results = {
        'water_depth_m': water_depth,
        'flow_velocity_m_s': flow_velocity,
        'sampling_rate_L_s': sampling_rate,
        'nozzle_diameter_mm': sampler_nozzle_diameter,
        'nozzle_area_m2': nozzle_area,
        'sampling_velocity_m_s': sampling_velocity,
        'isokinetic_ratio': isokinetic_ratio,
        'transit_rate_m_s': transit_rate,
        'round_trip_time_s': round_trip_time,
        'sample_volume_L': sample_volume * 1000,
        'water_volume_sampled_m3': water_volume_sampled,
        'is_isokinetic': abs(isokinetic_ratio - 1.0) < 0.1  # Within 10%
    }
    
    return results


def point_sampling_protocol(
    sampling_depths: np.ndarray,
    water_depth: float,
    flow_velocity_profile: np.ndarray,
    sampling_duration: float = 60.0
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Design point sampling protocol at multiple depths.
    
    Point sampling collects sediment samples at specific depths
    rather than integrating over the entire water column.
    
    Parameters
    ----------
    sampling_depths : array-like
        Depths below surface for sampling [m]
    water_depth : float
        Total water depth [m]
    flow_velocity_profile : array-like
        Flow velocity at each sampling depth [m/s]
    sampling_duration : float, optional
        Sampling duration at each point [s], default 60.0
        
    Returns
    -------
    dict
        Dictionary containing point sampling design
    """
    depths = validate_array_like(sampling_depths, "sampling_depths")
    velocities = validate_array_like(flow_velocity_profile, "flow_velocity_profile")
    validate_positive(water_depth, "water_depth")
    validate_positive(sampling_duration, "sampling_duration")
    
    if len(depths) != len(velocities):
        raise ValueError("Number of depths must match number of velocities")
    
    # Check that all depths are within water column
    if np.any(depths >= water_depth):
        raise ValueError("All sampling depths must be less than water depth")
    
    # Calculate relative depths
    relative_depths = depths / water_depth
    
    # Standard sampling depths (if not specified)
    standard_depths = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9]) * water_depth
    
    # Calculate sampling weights based on velocity
    velocity_weights = velocities / np.sum(velocities)
    
    # Estimate concentration variation with depth (simplified)
    # Higher concentrations typically near bed
    depth_factor = np.exp(-2 * relative_depths)  # Exponential decay from bed
    concentration_weights = depth_factor / np.sum(depth_factor)
    
    results = {
        'sampling_depths_m': depths,
        'relative_depths': relative_depths,
        'water_depth_m': water_depth,
        'flow_velocities_m_s': velocities,
        'sampling_duration_s': sampling_duration,
        'velocity_weights': velocity_weights,
        'concentration_weights': concentration_weights,
        'standard_depths_m': standard_depths,
        'number_of_points': len(depths),
        'total_sampling_time_s': len(depths) * sampling_duration
    }
    
    return results


def discharge_weighted_sampling(
    cross_section_areas: np.ndarray,
    flow_velocities: np.ndarray,
    sediment_concentrations: np.ndarray
) -> Dict[str, float]:
    """
    Calculate discharge-weighted sediment concentration.
    
    The discharge-weighted concentration accounts for the varying
    flow velocity across the cross-section.
    
    Parameters
    ----------
    cross_section_areas : array-like
        Area of each subsection [m²]
    flow_velocities : array-like
        Mean velocity in each subsection [m/s]
    sediment_concentrations : array-like
        Sediment concentration in each subsection [mg/L]
        
    Returns
    -------
    dict
        Dictionary containing discharge-weighted results
    """
    areas = validate_array_like(cross_section_areas, "cross_section_areas")
    velocities = validate_array_like(flow_velocities, "flow_velocities")
    concentrations = validate_array_like(sediment_concentrations, "sediment_concentrations")
    
    validate_positive(areas, "cross_section_areas")
    validate_positive(velocities, "flow_velocities")
    validate_positive(concentrations, "sediment_concentrations")
    
    if not (len(areas) == len(velocities) == len(concentrations)):
        raise ValueError("All arrays must have the same length")
    
    # Calculate discharge in each subsection
    discharges = areas * velocities
    total_discharge = np.sum(discharges)
    
    # Calculate sediment load in each subsection
    sediment_loads = concentrations * discharges
    total_sediment_load = np.sum(sediment_loads)
    
    # Discharge-weighted concentration
    discharge_weighted_concentration = total_sediment_load / total_discharge
    
    # Area-weighted concentration (for comparison)
    total_area = np.sum(areas)
    area_weighted_concentration = np.sum(concentrations * areas) / total_area
    
    # Flow-weighted concentration
    velocity_weighted_concentration = np.sum(concentrations * velocities) / np.sum(velocities)
    
    results = {
        'total_discharge_m3_s': total_discharge,
        'total_sediment_load_mg_s': total_sediment_load / 1000.0,  # Convert to kg/s
        'discharge_weighted_concentration_mg_L': discharge_weighted_concentration,
        'area_weighted_concentration_mg_L': area_weighted_concentration,
        'velocity_weighted_concentration_mg_L': velocity_weighted_concentration,
        'concentration_range_mg_L': [np.min(concentrations), np.max(concentrations)],
        'discharge_range_m3_s': [np.min(discharges), np.max(discharges)],
        'number_of_subsections': len(areas),
        'total_cross_section_area_m2': total_area
    }
    
    return results


def automatic_sampler_programming(
    sampling_interval: float,
    sampling_duration: float,
    bottle_volume: float,
    number_of_bottles: int,
    flow_hydrograph: Optional[np.ndarray] = None
) -> Dict[str, Union[float, int, np.ndarray]]:
    """
    Program automatic sampler for suspended sediment collection.
    
    Parameters
    ----------
    sampling_interval : float
        Time interval between samples [hours]
    sampling_duration : float
        Duration of each sample collection [seconds]
    bottle_volume : float
        Volume of each sample bottle [mL]
    number_of_bottles : int
        Total number of sample bottles
    flow_hydrograph : array-like, optional
        Flow values for flow-weighted sampling [m³/s]
        
    Returns
    -------
    dict
        Dictionary containing sampler programming parameters
    """
    validate_positive(sampling_interval, "sampling_interval")
    validate_positive(sampling_duration, "sampling_duration")
    validate_positive(bottle_volume, "bottle_volume")
    validate_positive(number_of_bottles, "number_of_bottles")
    
    # Convert sampling interval to seconds
    sampling_interval_s = sampling_interval * 3600.0
    
    # Total monitoring duration
    total_duration_hours = number_of_bottles * sampling_interval
    total_duration_s = total_duration_hours * 3600.0
    
    # Sample collection rate
    sample_volume_per_s = bottle_volume / sampling_duration  # mL/s
    
    # Generate sampling schedule
    sampling_times_hours = np.arange(0, total_duration_hours, sampling_interval)
    
    # Flow-weighted sampling if hydrograph provided
    if flow_hydrograph is not None:
        flows = validate_array_like(flow_hydrograph, "flow_hydrograph")
        validate_positive(flows, "flow_hydrograph")
        
        if len(flows) != len(sampling_times_hours):
            # Interpolate flows to match sampling times
            time_flow = np.linspace(0, total_duration_hours, len(flows))
            flows_interp = np.interp(sampling_times_hours, time_flow, flows)
        else:
            flows_interp = flows
        
        # Adjust sampling duration based on flow
        mean_flow = np.mean(flows_interp)
        flow_ratios = flows_interp / mean_flow
        adjusted_durations = sampling_duration * flow_ratios
        
        # Ensure durations don't exceed bottle capacity
        max_duration = bottle_volume / sample_volume_per_s
        adjusted_durations = np.minimum(adjusted_durations, max_duration)
        
        flow_weighted = True
    else:
        flows_interp = None
        adjusted_durations = np.full(len(sampling_times_hours), sampling_duration)
        flow_weighted = False
    
    results = {
        'sampling_interval_hours': sampling_interval,
        'sampling_duration_s': sampling_duration,
        'bottle_volume_mL': bottle_volume,
        'number_of_bottles': number_of_bottles,
        'total_duration_hours': total_duration_hours,
        'sample_volume_per_s_mL': sample_volume_per_s,
        'sampling_times_hours': sampling_times_hours,
        'adjusted_durations_s': adjusted_durations,
        'flow_weighted': flow_weighted,
        'flows_m3_s': flows_interp,
        'total_sample_volume_mL': np.sum(adjusted_durations * sample_volume_per_s)
    }
    
    return results


def isokinetic_sampling_velocity(
    stream_velocity: Union[float, np.ndarray],
    nozzle_diameter: float,
    target_sample_rate: float
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate required sampling parameters for isokinetic sampling.
    
    Isokinetic sampling ensures that the velocity of water entering
    the sampler nozzle equals the stream velocity.
    
    Parameters
    ----------
    stream_velocity : float or array-like
        Stream velocity at sampling point [m/s]
    nozzle_diameter : float
        Sampler nozzle diameter [mm]
    target_sample_rate : float
        Desired sample collection rate [L/min]
        
    Returns
    -------
    dict
        Dictionary containing isokinetic sampling parameters
    """
    validate_positive(stream_velocity, "stream_velocity")
    validate_positive(nozzle_diameter, "nozzle_diameter")
    validate_positive(target_sample_rate, "target_sample_rate")
    
    # Convert units
    nozzle_diameter_m = nozzle_diameter / 1000.0
    nozzle_area = np.pi * (nozzle_diameter_m / 2.0)**2
    target_sample_rate_m3_s = target_sample_rate / (1000.0 * 60.0)  # L/min to m³/s
    
    # Calculate required sampling velocity for isokinetic conditions
    required_sampling_velocity = stream_velocity  # For perfect isokinetic sampling
    
    # Calculate actual sample rate that would result from isokinetic sampling
    isokinetic_sample_rate = required_sampling_velocity * nozzle_area  # m³/s
    isokinetic_sample_rate_L_min = isokinetic_sample_rate * 1000.0 * 60.0  # L/min
    
    # Calculate velocity ratio if using target sample rate
    actual_sampling_velocity = target_sample_rate_m3_s / nozzle_area
    velocity_ratio = actual_sampling_velocity / stream_velocity
    
    # Sampling efficiency (isokinetic efficiency)
    # Efficiency decreases as velocity ratio deviates from 1.0
    efficiency = 1.0 - 0.1 * np.abs(velocity_ratio - 1.0)  # Simplified relationship
    efficiency = np.maximum(efficiency, 0.5)  # Minimum 50% efficiency
    
    results = {
        'stream_velocity_m_s': stream_velocity,
        'nozzle_diameter_mm': nozzle_diameter,
        'nozzle_area_m2': nozzle_area,
        'target_sample_rate_L_min': target_sample_rate,
        'required_sampling_velocity_m_s': required_sampling_velocity,
        'isokinetic_sample_rate_L_min': isokinetic_sample_rate_L_min,
        'actual_sampling_velocity_m_s': actual_sampling_velocity,
        'velocity_ratio': velocity_ratio,
        'sampling_efficiency': efficiency,
        'is_isokinetic': np.abs(velocity_ratio - 1.0) < 0.1,  # Within 10%
        'deviation_from_isokinetic_percent': (velocity_ratio - 1.0) * 100
    }
    
    return results


def suspended_sediment_rating_curve(
    discharge_data: np.ndarray,
    concentration_data: np.ndarray,
    rating_curve_type: str = "power_law"
) -> Dict[str, Union[float, np.ndarray, str]]:
    """
    Develop suspended sediment rating curve relating discharge to concentration.
    
    Parameters
    ----------
    discharge_data : array-like
        Discharge measurements [m³/s]
    concentration_data : array-like
        Corresponding sediment concentration measurements [mg/L]
    rating_curve_type : str, optional
        Type of rating curve: 'power_law' or 'linear', default 'power_law'
        
    Returns
    -------
    dict
        Dictionary containing rating curve parameters and statistics
    """
    discharge = validate_array_like(discharge_data, "discharge_data")
    concentration = validate_array_like(concentration_data, "concentration_data")
    
    validate_positive(discharge, "discharge_data")
    validate_positive(concentration, "concentration_data")
    
    if len(discharge) != len(concentration):
        raise ValueError("Discharge and concentration data must have same length")
    
    if rating_curve_type == "power_law":
        # Power law: C = aQ^b
        # Log-transform: log(C) = log(a) + b*log(Q)
        log_Q = np.log(discharge)
        log_C = np.log(concentration)
        
        # Linear regression on log-transformed data
        coefficients = np.polyfit(log_Q, log_C, 1)
        b = coefficients[0]  # Exponent
        log_a = coefficients[1]  # Intercept
        a = np.exp(log_a)  # Coefficient
        
        # Predicted concentrations
        C_predicted = a * discharge**b
        
        # Calculate R² for original data
        C_mean = np.mean(concentration)
        ss_tot = np.sum((concentration - C_mean)**2)
        ss_res = np.sum((concentration - C_predicted)**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        curve_equation = f"C = {a:.3f} * Q^{b:.3f}"
        
    elif rating_curve_type == "linear":
        # Linear: C = aQ + b
        coefficients = np.polyfit(discharge, concentration, 1)
        slope = coefficients[0]
        intercept = coefficients[1]
        
        # Predicted concentrations
        C_predicted = slope * discharge + intercept
        
        # Calculate R²
        C_mean = np.mean(concentration)
        ss_tot = np.sum((concentration - C_mean)**2)
        ss_res = np.sum((concentration - C_predicted)**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        curve_equation = f"C = {slope:.3f} * Q + {intercept:.3f}"
        a, b = slope, intercept  # For consistency
        
    else:
        raise ValueError("rating_curve_type must be 'power_law' or 'linear'")
    
    # Calculate residuals and error statistics
    residuals = concentration - C_predicted
    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))
    
    results = {
        'rating_curve_type': rating_curve_type,
        'coefficient_a': a,
        'exponent_b': b if rating_curve_type == 'power_law' else None,
        'slope': slope if rating_curve_type == 'linear' else None,
        'intercept': intercept if rating_curve_type == 'linear' else None,
        'r_squared': r_squared,
        'rmse_mg_L': rmse,
        'mae_mg_L': mae,
        'n_data_points': len(discharge),

        'concentration_range_mg_L': [np.min(concentration), np.max(concentration)],
        'equation': curve_equation,
        'discharge_data_m3_s': discharge,
        'concentration_data_mg_L': concentration,
        'predicted_concentration_mg_L': C_predicted,
        'residuals_mg_L': residuals
    }
    
    return results


def sampling_frequency_optimization(
    hydrograph_duration: float,
    peak_discharge: float,
    base_discharge: float,
    available_bottles: int,
    priority_periods: Optional[List[Tuple[float, float]]] = None
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Optimize sampling frequency for varying flow conditions.
    
    Parameters
    ----------
    hydrograph_duration : float
        Total duration of hydrograph [hours]
    peak_discharge : float
        Peak discharge value [m³/s]
    base_discharge : float
        Base flow discharge [m³/s]
    available_bottles : int
        Number of sample bottles available
    priority_periods : list of tuples, optional
        List of (start_time, end_time) for high-priority sampling periods [hours]
        
    Returns
    -------
    dict
        Dictionary containing optimized sampling schedule
    """
    validate_positive(hydrograph_duration, "hydrograph_duration")
    validate_positive(peak_discharge, "peak_discharge")
    validate_positive(base_discharge, "base_discharge")
    validate_positive(available_bottles, "available_bottles")
    
    if peak_discharge <= base_discharge:
        raise ValueError("Peak discharge must be greater than base discharge")
    
    # Generate synthetic hydrograph for optimization
    time_hours = np.linspace(0, hydrograph_duration, 1000)
    # Simplified gamma-shaped hydrograph
    peak_time = hydrograph_duration * 0.3  # Peak at 30% of duration
    hydrograph = base_discharge + (peak_discharge - base_discharge) * \
                 np.exp(-(time_hours - peak_time)**2 / (0.1 * hydrograph_duration)**2)
    
    # Calculate sampling priority weights
    # Higher weight for higher discharges and priority periods
    discharge_weights = (hydrograph - base_discharge) / (peak_discharge - base_discharge)
    
    # Add priority period weights
    if priority_periods is not None:
        for start_time, end_time in priority_periods:
            priority_mask = (time_hours >= start_time) & (time_hours <= end_time)
            discharge_weights[priority_mask] *= 2.0  # Double weight in priority periods
    
    # Normalize weights
    total_weight = np.sum(discharge_weights)
    normalized_weights = discharge_weights / total_weight
    
    # Distribute samples based on weights
    samples_per_timestep = normalized_weights * available_bottles
    
    # Select sampling times based on cumulative distribution
    cumulative_samples = np.cumsum(samples_per_timestep)
    sample_indices = np.unique(np.searchsorted(cumulative_samples, 
                                             np.linspace(0.5, available_bottles-0.5, available_bottles)))
    
    sampling_times = time_hours[sample_indices]
    sampling_discharges = hydrograph[sample_indices]
    
    # Calculate sampling intervals
    sampling_intervals = np.diff(np.concatenate([[0], sampling_times, [hydrograph_duration]]))
    
    results = {
        'hydrograph_duration_hours': hydrograph_duration,
        'available_bottles': available_bottles,
        'actual_samples': len(sampling_times),
        'sampling_times_hours': sampling_times,
        'sampling_discharges_m3_s': sampling_discharges,
        'sampling_intervals_hours': sampling_intervals,
        'average_interval_hours': np.mean(sampling_intervals[1:-1]),  # Exclude first/last
        'minimum_interval_hours': np.min(sampling_intervals[1:-1]),
        'maximum_interval_hours': np.max(sampling_intervals[1:-1]),
        'peak_discharge_m3_s': peak_discharge,
        'base_discharge_m3_s': base_discharge,
        'discharge_range_sampled': [np.min(sampling_discharges), np.max(sampling_discharges)],
        'hydrograph_coverage_percent': (len(sampling_times) / len(time_hours)) * 100
    }
    
    return results