"""
Sampling equipment specifications and calibration functions (Chapter 5).

This module provides functions for sediment sampling equipment specifications,
calibration curves, efficiency calculations, and selection criteria.
"""

import numpy as np
from typing import Union, Tuple, List, Dict, Optional
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

def sampler_specifications(
    sampler_type: str,
    nozzle_diameter: float,
    intake_efficiency: float = 1.0
) -> Dict[str, Union[str, float]]:
    """
    Get standard specifications for sediment sampling equipment.
    
    Parameters
    ----------
    sampler_type : str
        Type of sampler: 'US_DH-48', 'US_DH-59', 'US_DH-76', 'US_D-74', 
        'US_P-63', 'US_BM-54', 'Delft_bottle', 'Van_Dorn'
    nozzle_diameter : float
        Nozzle diameter [mm]
    intake_efficiency : float, optional
        Intake efficiency factor [0-1], default 1.0
        
    Returns
    -------
    dict
        Dictionary containing sampler specifications
        
    References
    ----------
    ASCE Manual 110, Chapter 5: Sediment Sampling Equipment
    """
    validate_positive(nozzle_diameter, "nozzle_diameter")
    validate_range(intake_efficiency, "intake_efficiency", 0.0, 1.0)
    
    # Standard sampler specifications database
    sampler_database = {
        'US_DH-48': {
            'description': 'US DH-48 depth-integrating suspended sediment sampler',
            'body_material': 'bronze',
            'weight_kg': 22.0,
            'capacity_L': 0.5,
            'max_velocity_m_s': 4.6,
            'nozzle_options_mm': [3.2, 4.8, 6.4, 9.5, 12.7],
            'applications': ['wadeable streams', 'bridge sampling'],
            'isokinetic_range': [0.7, 1.3]
        },
        'US_DH-59': {
            'description': 'US DH-59 depth-integrating suspended sediment sampler', 
            'body_material': 'aluminum',
            'weight_kg': 30.0,
            'capacity_L': 1.0,
            'max_velocity_m_s': 3.0,
            'nozzle_options_mm': [4.8, 6.4, 9.5, 12.7, 16.0],
            'applications': ['deep water', 'high velocity'],
            'isokinetic_range': [0.8, 1.2]
        },
        'US_DH-76': {
            'description': 'US DH-76 depth-integrating suspended sediment sampler',
            'body_material': 'fiberglass',
            'weight_kg': 45.0,
            'capacity_L': 3.0,
            'max_velocity_m_s': 3.5,
            'nozzle_options_mm': [6.4, 9.5, 12.7, 16.0, 19.0],
            'applications': ['large rivers', 'heavy suspended loads'],
            'isokinetic_range': [0.8, 1.2]
        },
        'US_D-74': {
            'description': 'US D-74 instantaneous suspended sediment sampler',
            'body_material': 'aluminum',
            'weight_kg': 7.0,
            'capacity_L': 0.4,
            'max_velocity_m_s': 5.0,
            'nozzle_options_mm': [4.8, 6.4, 9.5],
            'applications': ['point sampling', 'research'],
            'isokinetic_range': [0.9, 1.1]
        },
        'US_P-63': {
            'description': 'US P-63 point-integrating suspended sediment sampler',
            'body_material': 'bronze',
            'weight_kg': 95.0,
            'capacity_L': 1.0,
            'max_velocity_m_s': 7.6,
            'nozzle_options_mm': [6.4, 9.5, 12.7, 16.0],
            'applications': ['point sampling', 'high velocity'],
            'isokinetic_range': [0.9, 1.1]
        },
        'US_BM-54': {
            'description': 'US BM-54 bed material sampler',
            'body_material': 'bronze',
            'weight_kg': 15.0,
            'capacity_L': 0.2,
            'max_velocity_m_s': 3.0,
            'nozzle_options_mm': [6.4, 9.5],
            'applications': ['bed material', 'coarse sediment'],
            'isokinetic_range': [0.8, 1.2]
        },
        'Delft_bottle': {
            'description': 'Delft bottle sediment sampler',
            'body_material': 'plastic',
            'weight_kg': 2.0,
            'capacity_L': 1.0,
            'max_velocity_m_s': 2.0,
            'nozzle_options_mm': [10.0, 15.0, 20.0],
            'applications': ['low velocity', 'fine sediment'],
            'isokinetic_range': [0.7, 1.3]
        },
        'Van_Dorn': {
            'description': 'Van Dorn water sampler (modified for sediment)',
            'body_material': 'PVC',
            'weight_kg': 3.5,
            'capacity_L': 2.5,
            'max_velocity_m_s': 1.5,
            'nozzle_options_mm': [15.0, 20.0, 25.0],
            'applications': ['lakes', 'reservoirs', 'low velocity'],
            'isokinetic_range': [0.6, 1.4]
        }
    }
    
    if sampler_type not in sampler_database:
        raise ValueError(f"Unknown sampler type: {sampler_type}")
    
    specs = sampler_database[sampler_type].copy()
    
    # Add calculated parameters
    nozzle_area = np.pi * (nozzle_diameter / 2000.0)**2  # m²
    
    # Check if nozzle diameter is available for this sampler
    available_nozzles = specs['nozzle_options_mm']
    nozzle_available = any(abs(d - nozzle_diameter) < 0.1 for d in available_nozzles)
    
    specs.update({
        'sampler_type': sampler_type,
        'nozzle_diameter_mm': nozzle_diameter,
        'nozzle_area_m2': nozzle_area,
        'intake_efficiency': intake_efficiency,
        'nozzle_available': nozzle_available,
        'effective_intake_area_m2': nozzle_area * intake_efficiency
    })
    
    return specs


def isokinetic_sampler_efficiency(
    stream_velocity: Union[float, np.ndarray],
    sampling_velocity: Union[float, np.ndarray],
    particle_size: Union[float, np.ndarray] = 50.0,
    particle_density: float = 2650.0
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate isokinetic sampling efficiency.
    
    Efficiency depends on the ratio of sampling velocity to stream velocity
    and particle characteristics.
    
    Parameters
    ----------
    stream_velocity : float or array-like
        Stream velocity [m/s]
    sampling_velocity : float or array-like
        Velocity of water entering sampler [m/s]
    particle_size : float or array-like, optional
        Particle diameter [μm], default 50.0
    particle_density : float, optional
        Particle density [kg/m³], default 2650.0
        
    Returns
    -------
    dict
        Dictionary containing efficiency calculations
    """
    validate_positive(stream_velocity, "stream_velocity")
    validate_positive(sampling_velocity, "sampling_velocity")
    validate_positive(particle_size, "particle_size")
    validate_positive(particle_density, "particle_density")
    
    # Velocity ratio
    velocity_ratio = sampling_velocity / stream_velocity
    
    # Particle settling velocity (Stokes law for fine particles)
    d_particle = particle_size * 1e-6  # Convert μm to m
    kinematic_viscosity = 1e-6  # m²/s (water at 20°C)
    
    settling_velocity = (particle_density - WATER_DENSITY) * GRAVITY * d_particle**2 / (18 * kinematic_viscosity * WATER_DENSITY)
    
    # Particle Stokes number
    stokes_number = settling_velocity / stream_velocity
    
    # Sampling efficiency model (simplified)
    # Based on theoretical and empirical relationships
    
    if np.isscalar(velocity_ratio):
        velocity_ratios = np.array([velocity_ratio])
    else:
        velocity_ratios = np.array(velocity_ratio)
    
    efficiencies = np.zeros_like(velocity_ratios)
    
    for i, vr in enumerate(velocity_ratios):
        if vr < 0.7:
            # Under-sampling: efficiency decreases
            efficiencies[i] = 0.5 + 0.5 * (vr / 0.7)**2
        elif vr > 1.3:
            # Over-sampling: efficiency decreases
            efficiencies[i] = 1.0 - 0.3 * ((vr - 1.3) / 0.7)**2
        else:
            # Near-isokinetic: high efficiency
            deviation = abs(vr - 1.0)
            efficiencies[i] = 1.0 - 0.2 * (deviation / 0.3)**2
    
    # Adjust for particle size effects
    if np.any(stokes_number > 0.1):
        # Coarse particles more affected by non-isokinetic sampling
        size_correction = 1.0 - 0.5 * np.maximum(0, stokes_number - 0.1)
        efficiencies *= size_correction
    
    # Ensure efficiencies are within reasonable bounds
    efficiencies = np.clip(efficiencies, 0.3, 1.0)
    
    if np.isscalar(velocity_ratio):
        efficiency = float(efficiencies[0])
        stokes_num = float(stokes_number) if np.isscalar(stokes_number) else stokes_number
    else:
        efficiency = efficiencies
        stokes_num = stokes_number
    
    results = {
        'velocity_ratio': velocity_ratio,
        'sampling_efficiency': efficiency,
        'particle_settling_velocity_m_s': settling_velocity,
        'particle_stokes_number': stokes_num,
        'isokinetic_quality': 'excellent' if np.all(np.abs(velocity_ratio - 1.0) < 0.1) 
                             else 'good' if np.all(np.abs(velocity_ratio - 1.0) < 0.2)
                             else 'acceptable' if np.all(np.abs(velocity_ratio - 1.0) < 0.3)
                             else 'poor',
        'stream_velocity_m_s': stream_velocity,
        'sampling_velocity_m_s': sampling_velocity,
        'particle_size_um': particle_size
    }
    
    return results


def equipment_calibration_curve(
    calibration_flows: np.ndarray,
    measured_concentrations: np.ndarray,
    known_concentrations: np.ndarray,
    equipment_type: str = "pump_sampler"
) -> Dict[str, Union[float, np.ndarray, str]]:
    """
    Develop equipment calibration curve for concentration measurements.
    
    Parameters
    ----------
    calibration_flows : array-like
        Flow rates during calibration [L/s]
    measured_concentrations : array-like
        Concentrations measured by equipment [mg/L]
    known_concentrations : array-like
        Known reference concentrations [mg/L]
    equipment_type : str, optional
        Type of equipment being calibrated
        
    Returns
    -------
    dict
        Dictionary containing calibration curve parameters
    """
    flows = validate_array_like(calibration_flows, "calibration_flows")
    measured = validate_array_like(measured_concentrations, "measured_concentrations")
    known = validate_array_like(known_concentrations, "known_concentrations")
    
    validate_positive(flows, "calibration_flows")
    validate_positive(measured, "measured_concentrations")
    validate_positive(known, "known_concentrations")
    
    if not (len(flows) == len(measured) == len(known)):
        raise ValueError("All calibration arrays must have same length")
    
    # Linear calibration: C_true = a * C_measured + b
    calibration_coeffs = np.polyfit(measured, known, 1)
    slope = calibration_coeffs[0]
    intercept = calibration_coeffs[1]
    
    # Predicted values
    predicted = slope * measured + intercept
    
    # Statistical measures
    residuals = known - predicted
    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))
    
    # R-squared
    known_mean = np.mean(known)
    ss_tot = np.sum((known - known_mean)**2)
    ss_res = np.sum(residuals**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Bias and precision
    bias = np.mean(residuals)
    precision = np.std(residuals, ddof=1)
    
    # Flow-dependent analysis
    flow_categories = ['low', 'medium', 'high']
    flow_percentiles = [33, 67]
    flow_thresholds = np.percentile(flows, flow_percentiles)
    
    flow_analysis = {}
    for i, category in enumerate(flow_categories):
        if i == 0:
            mask = flows <= flow_thresholds[0]
        elif i == 1:
            mask = (flows > flow_thresholds[0]) & (flows <= flow_thresholds[1])
        else:
            mask = flows > flow_thresholds[1]
        
        if np.any(mask):
            category_residuals = residuals[mask]
            flow_analysis[f'{category}_flow'] = {
                'count': np.sum(mask),
                'bias_mg_L': np.mean(category_residuals),
                'precision_mg_L': np.std(category_residuals, ddof=1) if np.sum(mask) > 1 else 0.0,
                'flow_range_L_s': [np.min(flows[mask]), np.max(flows[mask])]
            }
    
    results = {
        'equipment_type': equipment_type,
        'calibration_slope': slope,
        'calibration_intercept': intercept,
        'r_squared': r_squared,
        'rmse_mg_L': rmse,
        'mae_mg_L': mae,
        'bias_mg_L': bias,
        'precision_mg_L': precision,
        'n_calibration_points': len(flows),
        'flow_range_L_s': [np.min(flows), np.max(flows)],
        'concentration_range_mg_L': [np.min(known), np.max(known)],
        'calibration_equation': f'C_true = {slope:.3f} * C_measured + {intercept:.3f}',
        'flow_dependent_analysis': flow_analysis,
        'flows_L_s': flows,
        'measured_concentrations_mg_L': measured,
        'known_concentrations_mg_L': known,
        'predicted_concentrations_mg_L': predicted,
        'residuals_mg_L': residuals
    }
    
    return results


def sampler_selection_criteria(
    water_depth: float,
    flow_velocity: float,
    expected_concentration: float,
    particle_size_range: Tuple[float, float],
    sampling_purpose: str = "routine_monitoring"
) -> Dict[str, Union[str, List[str], bool]]:
    """
    Recommend appropriate sampler based on site conditions.
    
    Parameters
    ----------
    water_depth : float
        Water depth at sampling site [m]
    flow_velocity : float
        Flow velocity at sampling site [m/s]
    expected_concentration : float
        Expected sediment concentration [mg/L]
    particle_size_range : tuple
        Expected particle size range (min, max) [μm]
    sampling_purpose : str, optional
        Purpose: 'routine_monitoring', 'research', 'load_calculation', 'calibration'
        
    Returns
    -------
    dict
        Dictionary containing sampler recommendations
    """
    validate_positive(water_depth, "water_depth")
    validate_positive(flow_velocity, "flow_velocity")
    validate_positive(expected_concentration, "expected_concentration")
    validate_positive(particle_size_range[0], "particle_size_range minimum")
    validate_positive(particle_size_range[1], "particle_size_range maximum")
    
    if particle_size_range[1] <= particle_size_range[0]:
        raise ValueError("Maximum particle size must be greater than minimum")
    
    recommendations = []
    suitable_samplers = []
    unsuitable_samplers = []
    
    # Sampler suitability criteria
    sampler_criteria = {
        'US_DH-48': {
            'max_depth': 5.0,
            'max_velocity': 4.6,
            'min_concentration': 10.0,
            'max_particle_size': 2000.0,
            'applications': ['wadeable_streams', 'routine_monitoring']
        },
        'US_DH-59': {
            'max_depth': 10.0,
            'max_velocity': 3.0,
            'min_concentration': 5.0,
            'max_particle_size': 5000.0,
            'applications': ['deep_water', 'load_calculation']
        },
        'US_DH-76': {
            'max_depth': 15.0,
            'max_velocity': 3.5,
            'min_concentration': 20.0,
            'max_particle_size': 10000.0,
            'applications': ['large_rivers', 'high_concentration']
        },
        'US_P-63': {
            'max_depth': 20.0,
            'max_velocity': 7.6,
            'min_concentration': 1.0,
            'max_particle_size': 1000.0,
            'applications': ['research', 'point_sampling']
        },
        'Delft_bottle': {
            'max_depth': 3.0,
            'max_velocity': 2.0,
            'min_concentration': 1.0,
            'max_particle_size': 500.0,
            'applications': ['low_velocity', 'research']
        }
    }
    
    # Evaluate each sampler
    for sampler, criteria in sampler_criteria.items():
        suitability_score = 0
        issues = []
        
        # Check depth
        if water_depth <= criteria['max_depth']:
            suitability_score += 1
        else:
            issues.append(f"depth too great (>{criteria['max_depth']}m)")
        
        # Check velocity
        if flow_velocity <= criteria['max_velocity']:
            suitability_score += 1
        else:
            issues.append(f"velocity too high (>{criteria['max_velocity']}m/s)")
        
        # Check concentration
        if expected_concentration >= criteria['min_concentration']:
            suitability_score += 1
        else:
            issues.append(f"concentration too low (<{criteria['min_concentration']}mg/L)")
        
        # Check particle size
        if particle_size_range[1] <= criteria['max_particle_size']:
            suitability_score += 1
        else:
            issues.append(f"particles too large (>{criteria['max_particle_size']}μm)")
        
        # Check application match
        purpose_mapping = {
            'routine_monitoring': 'routine_monitoring',
            'research': 'research', 
            'load_calculation': 'load_calculation',
            'calibration': 'research'
        }
        
        if purpose_mapping[sampling_purpose] in criteria['applications']:
            suitability_score += 1
        
        # Categorize sampler
        if suitability_score >= 4:
            suitable_samplers.append({
                'sampler': sampler,
                'score': suitability_score,
                'issues': issues
            })
        else:
            unsuitable_samplers.append({
                'sampler': sampler,
                'score': suitability_score,
                'issues': issues
            })
    
    # Sort by suitability score
    suitable_samplers.sort(key=lambda x: x['score'], reverse=True)
    
    # Generate recommendations
    if suitable_samplers:
        best_sampler = suitable_samplers[0]
        recommendations.append(f"Primary recommendation: {best_sampler['sampler']} (score: {best_sampler['score']}/5)")
        
        if len(suitable_samplers) > 1:
            alt_sampler = suitable_samplers[1]
            recommendations.append(f"Alternative: {alt_sampler['sampler']} (score: {alt_sampler['score']}/5)")
    else:
        recommendations.append("No fully suitable sampler found - consider modifying sampling approach")
    
    # Add specific recommendations based on conditions
    if flow_velocity > 5.0:
        recommendations.append("High velocity conditions - consider US P-63 for point sampling")
    
    if expected_concentration < 10.0:
        recommendations.append("Low concentration - ensure adequate sample volume and analysis sensitivity")
    
    if particle_size_range[1] > 2000.0:
        recommendations.append("Coarse particles present - consider bed material sampler or larger nozzle")
    
    results = {
        'water_depth_m': water_depth,
        'flow_velocity_m_s': flow_velocity,
        'expected_concentration_mg_L': expected_concentration,
        'particle_size_range_um': particle_size_range,
        'sampling_purpose': sampling_purpose,
        'suitable_samplers': [s['sampler'] for s in suitable_samplers],
        'unsuitable_samplers': [s['sampler'] for s in unsuitable_samplers],
        'recommendations': recommendations,
        'detailed_analysis': {
            'suitable': suitable_samplers,
            'unsuitable': unsuitable_samplers
        }
    }
    
    return results


def pump_sampler_design(
    intake_depth: float,
    pumping_rate: float,
    tube_diameter: float,
    tube_length: float,
    particle_size: float = 100.0
) -> Dict[str, Union[float, bool, str]]:
    """
    Design parameters for pump-based sediment sampling systems.
    
    Parameters
    ----------
    intake_depth : float
        Depth of intake below surface [m]
    pumping_rate : float
        Pump flow rate [L/min]
    tube_diameter : float
        Internal diameter of intake tube [mm]
    tube_length : float
        Length of intake tube [m]
    particle_size : float, optional
        Representative particle size [μm], default 100.0
        
    Returns
    -------
    dict
        Dictionary containing pump sampler design parameters
    """
    validate_positive(intake_depth, "intake_depth")
    validate_positive(pumping_rate, "pumping_rate")
    validate_positive(tube_diameter, "tube_diameter")
    validate_positive(tube_length, "tube_length")
    validate_positive(particle_size, "particle_size")
    
    # Convert units
    tube_diameter_m = tube_diameter / 1000.0
    tube_area = np.pi * (tube_diameter_m / 2.0)**2
    pumping_rate_m3_s = pumping_rate / (1000.0 * 60.0)
    
    # Flow velocity in tube
    tube_velocity = pumping_rate_m3_s / tube_area
    
    # Reynolds number in tube
    kinematic_viscosity = 1e-6  # m²/s (water at 20°C)
    reynolds_number = tube_velocity * tube_diameter_m / kinematic_viscosity
    
    # Friction factor (Darcy-Weisbach)
    if reynolds_number < 2300:
        friction_factor = 64 / reynolds_number  # Laminar
        flow_regime = "laminar"
    else:
        friction_factor = 0.316 / reynolds_number**0.25  # Turbulent (smooth pipes)
        flow_regime = "turbulent"
    
    # Head loss in tube
    head_loss = friction_factor * (tube_length / tube_diameter_m) * (tube_velocity**2) / (2 * GRAVITY)
    
    # Total head requirement (elevation + friction)
    total_head = intake_depth + head_loss
    
    # Particle settling in tube
    particle_diameter_m = particle_size * 1e-6
    particle_density = 2650.0  # kg/m³
    
    # Settling velocity (Stokes law)
    settling_velocity = ((particle_density - WATER_DENSITY) * GRAVITY * particle_diameter_m**2) / (18 * kinematic_viscosity * WATER_DENSITY)
    
    # Transit time in tube
    transit_time = tube_length / tube_velocity
    
    # Settling distance during transit
    settling_distance = settling_velocity * transit_time
    
    # Particle transport efficiency
    if settling_distance < tube_diameter_m:
        particle_efficiency = 1.0 - (settling_distance / tube_diameter_m) * 0.5
        transport_quality = "excellent"
    elif settling_distance < 2 * tube_diameter_m:
        particle_efficiency = 0.7 - (settling_distance / tube_diameter_m - 1.0) * 0.3
        transport_quality = "good"
    else:
        particle_efficiency = 0.4
        transport_quality = "poor"
    
    # Minimum velocity to prevent settling
    min_velocity_critical = 2.0 * settling_velocity  # Empirical criterion
    velocity_adequate = tube_velocity >= min_velocity_critical
    
    results = {
        'intake_depth_m': intake_depth,
        'pumping_rate_L_min': pumping_rate,
        'tube_diameter_mm': tube_diameter,
        'tube_length_m': tube_length,
        'tube_area_m2': tube_area,
        'tube_velocity_m_s': tube_velocity,
        'reynolds_number': reynolds_number,
        'flow_regime': flow_regime,
        'friction_factor': friction_factor,
        'head_loss_m': head_loss,
        'total_head_m': total_head,
        'particle_settling_velocity_m_s': settling_velocity,
        'transit_time_s': transit_time,
        'settling_distance_m': settling_distance,
        'particle_efficiency': particle_efficiency,
        'transport_quality': transport_quality,
        'minimum_critical_velocity_m_s': min_velocity_critical,
        'velocity_adequate': velocity_adequate,
        'design_recommendations': []
    }
    
    # Generate design recommendations
    recommendations = results['design_recommendations']
    
    if not velocity_adequate:
        recommendations.append(f"Increase pumping rate to >{min_velocity_critical * tube_area * 60000:.1f} L/min")
    
    if head_loss > intake_depth:
        recommendations.append("High friction losses - consider larger tube diameter")
    
    if particle_efficiency < 0.7:
        recommendations.append("Poor particle transport - reduce tube length or increase diameter")
    
    if reynolds_number > 50000:
        recommendations.append("Very high Reynolds number - check for excessive turbulence")
    
    if not recommendations:
        recommendations.append("Design appears adequate for specified conditions")
    
    return results