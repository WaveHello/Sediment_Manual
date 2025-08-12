"""
Field procedures and protocols for sediment transport measurements (Chapter 5).

This module implements field measurement protocols, site selection criteria,
data quality assessment, and equipment maintenance procedures.
"""

import numpy as np
from typing import Union, Tuple, List, Dict, Optional
from datetime import datetime, timedelta
from ..utils.validators import validate_positive, validate_range, validate_array_like

def sampling_site_selection(
    stream_characteristics: Dict[str, float],
    measurement_objectives: List[str],
    site_constraints: Optional[Dict[str, Union[str, float]]] = None
) -> Dict[str, Union[str, List[str], bool, float]]:
    """
    Evaluate and recommend sediment sampling site selection criteria.
    
    Parameters
    ----------
    stream_characteristics : dict
        Dictionary containing stream properties:
        - 'width_m': stream width [m]
        - 'depth_m': mean depth [m]
        - 'velocity_m_s': mean velocity [m/s]
        - 'slope': channel slope [m/m]
        - 'discharge_m3_s': discharge [m³/s]
    measurement_objectives : list
        List of measurement objectives:
        - 'suspended_load'
        - 'bed_load'
        - 'total_load'
        - 'size_distribution'
        - 'load_calculation'
        - 'rating_curve'
    site_constraints : dict, optional
        Site-specific constraints and requirements
        
    Returns
    -------
    dict
        Dictionary containing site selection analysis and recommendations
        
    References
    ----------
    ASCE Manual 110, Chapter 5: Field Sampling Procedures
    """
    # Validate required stream characteristics
    required_params = ['width_m', 'depth_m', 'velocity_m_s', 'discharge_m3_s']
    for param in required_params:
        if param not in stream_characteristics:
            raise ValueError(f"Missing required stream characteristic: {param}")
        validate_positive(stream_characteristics[param], param)
    
    width = stream_characteristics['width_m']
    depth = stream_characteristics['depth_m']
    velocity = stream_characteristics['velocity_m_s']
    discharge = stream_characteristics['discharge_m3_s']
    slope = stream_characteristics.get('slope', 0.001)
    
    # Initialize site evaluation
    site_quality_score = 0
    recommendations = []
    site_issues = []
    
    # Evaluate hydraulic conditions
    froude_number = velocity / np.sqrt(9.81 * depth)
    reynolds_number = velocity * depth / 1e-6  # Assuming water at 20°C
    
    # Channel geometry evaluation
    width_depth_ratio = width / depth
    
    # Hydraulic evaluation criteria
    if 0.1 <= froude_number <= 0.8:
        site_quality_score += 2
        recommendations.append("Good flow conditions - subcritical flow suitable for sampling")
    elif froude_number > 0.8:
        site_quality_score += 1
        site_issues.append("High Froude number - potential flow instability")
    else:
        site_issues.append("Very low Froude number - potential stagnant conditions")
    
    if reynolds_number > 2000:
        site_quality_score += 1
        recommendations.append("Turbulent flow conditions - good for mixing")
    else:
        site_issues.append("Low Reynolds number - potential laminar flow issues")
    
    # Channel geometry evaluation
    if 5 <= width_depth_ratio <= 50:
        site_quality_score += 2
        recommendations.append("Appropriate width-to-depth ratio")
    else:
        site_issues.append(f"Unusual width-to-depth ratio ({width_depth_ratio:.1f})")
    
    # Depth evaluation for sampling equipment
    if depth >= 0.5:
        site_quality_score += 2
        recommendations.append("Adequate depth for standard sampling equipment")
    elif depth >= 0.2:
        site_quality_score += 1
        recommendations.append("Shallow depth - use appropriate sampling equipment")
    else:
        site_issues.append("Very shallow depth - limited sampling options")
    
    # Velocity evaluation
    if 0.1 <= velocity <= 3.0:
        site_quality_score += 2
        recommendations.append("Good velocity range for most sampling equipment")
    elif velocity > 3.0:
        site_quality_score += 1
        site_issues.append("High velocity - use heavy-duty sampling equipment")
    else:
        site_issues.append("Very low velocity - potential settling in samplers")
    
    # Objective-specific evaluations
    objective_compatibility = {}
    
    for objective in measurement_objectives:
        compatible = True
        objective_notes = []
        
        if objective == "suspended_load":
            if velocity < 0.05:
                compatible = False
                objective_notes.append("Velocity too low for significant suspension")
            if depth < 0.3:
                objective_notes.append("Shallow depth limits vertical sampling")
            
        elif objective == "bed_load":
            if velocity < 0.1:
                compatible = False
                objective_notes.append("Velocity too low for bed load transport")
            if depth > 10:
                objective_notes.append("Great depth complicates bed load sampling")
                
        elif objective == "total_load":
            if velocity < 0.1:
                compatible = False
                objective_notes.append("Velocity too low for significant transport")
                
        elif objective == "rating_curve":
            if discharge < 0.1:
                objective_notes.append("Very low discharge may limit rating curve range")
            
        objective_compatibility[objective] = {
            'compatible': compatible,
            'notes': objective_notes
        }
    
    # Site accessibility and practical considerations
    if site_constraints:
        access_score = 0
        
        if site_constraints.get('road_access', False):
            access_score += 1
            recommendations.append("Good road access facilitates equipment transport")
        
        if site_constraints.get('bridge_available', False):
            access_score += 2
            recommendations.append("Bridge access enables safer high-flow sampling")
        
        if site_constraints.get('wading_safe', False):
            access_score += 1
            recommendations.append("Safe wading conditions")
        
        power_available = site_constraints.get('power_available', False)
        if power_available:
            access_score += 1
            recommendations.append("Power availability supports automatic sampling")
        
        site_quality_score += access_score
    
    # Overall site rating
    max_score = 10  # Adjust based on criteria
    site_rating_percent = (site_quality_score / max_score) * 100
    
    if site_rating_percent >= 80:
        overall_rating = "excellent"
    elif site_rating_percent >= 60:
        overall_rating = "good"
    elif site_rating_percent >= 40:
        overall_rating = "acceptable"
    else:
        overall_rating = "poor"
    
    # Generate specific recommendations
    if width > 20:
        recommendations.append("Wide channel - use multiple verticals for cross-sectional sampling")
    
    if velocity > 2.0:
        recommendations.append("High velocity conditions - use streamlined samplers")
    
    if froude_number > 0.7:
        recommendations.append("Near-critical flow - monitor for flow variations")
    
    results = {
        'stream_width_m': width,
        'stream_depth_m': depth,
        'flow_velocity_m_s': velocity,
        'discharge_m3_s': discharge,
        'froude_number': froude_number,
        'reynolds_number': reynolds_number,
        'width_depth_ratio': width_depth_ratio,
        'site_quality_score': site_quality_score,
        'site_rating_percent': site_rating_percent,
        'overall_rating': overall_rating,
        'measurement_objectives': measurement_objectives,
        'objective_compatibility': objective_compatibility,
        'recommendations': recommendations,
        'site_issues': site_issues,
        'suitable_for_objectives': all(obj['compatible'] for obj in objective_compatibility.values())
    }
    
    return results


def data_quality_assessment(
    field_measurements: Dict[str, np.ndarray],
    measurement_metadata: Dict[str, Union[str, float]],
    quality_thresholds: Optional[Dict[str, float]] = None
) -> Dict[str, Union[str, float, List[str], Dict]]:
    """
    Assess quality of field measurement data.
    
    Parameters
    ----------
    field_measurements : dict
        Dictionary containing field measurement arrays
    measurement_metadata : dict
        Metadata including conditions, equipment, procedures
    quality_thresholds : dict, optional
        Custom quality control thresholds
        
    Returns
    -------
    dict
        Dictionary containing data quality assessment
    """
    # Default quality thresholds
    default_thresholds = {
        'concentration_cv_max': 30.0,  # Maximum coefficient of variation (%)
        'temperature_range_max': 10.0,  # Maximum temperature range (°C)
        'velocity_cv_max': 20.0,  # Maximum velocity coefficient of variation (%)
        'min_sample_size': 5,  # Minimum number of samples
        'duplicate_agreement_max': 15.0  # Maximum disagreement in duplicates (%)
    }
    
    if quality_thresholds is None:
        quality_thresholds = default_thresholds
    else:
        # Merge with defaults
        thresholds = default_thresholds.copy()
        thresholds.update(quality_thresholds)
        quality_thresholds = thresholds
    
    quality_flags = []
    quality_metrics = {}
    data_completeness = {}
    
    # Evaluate each measurement parameter
    for parameter, data in field_measurements.items():
        data_array = validate_array_like(data, f"field_measurements[{parameter}]")
        
        n_samples = len(data_array)
        if n_samples == 0:
            quality_flags.append(f"No data for {parameter}")
            continue
            
        # Basic statistics
        mean_val = np.mean(data_array)
        std_val = np.std(data_array, ddof=1) if n_samples > 1 else 0
        cv = (std_val / mean_val * 100) if mean_val != 0 else np.inf
        
        # Data completeness
        non_zero_count = np.count_nonzero(data_array)
        completeness_percent = (non_zero_count / n_samples) * 100
        
        data_completeness[parameter] = {
            'total_samples': n_samples,
            'non_zero_samples': non_zero_count,
            'completeness_percent': completeness_percent
        }
        
        # Quality checks
        parameter_flags = []
        
        if n_samples < quality_thresholds['min_sample_size']:
            parameter_flags.append(f"Small sample size ({n_samples})")
        
        if 'concentration' in parameter.lower():
            if cv > quality_thresholds['concentration_cv_max']:
                parameter_flags.append(f"High variability (CV={cv:.1f}%)")
        
        if 'velocity' in parameter.lower():
            if cv > quality_thresholds['velocity_cv_max']:
                parameter_flags.append(f"High velocity variability (CV={cv:.1f}%)")
        
        if 'temperature' in parameter.lower():
            temp_range = np.max(data_array) - np.min(data_array)
            if temp_range > quality_thresholds['temperature_range_max']:
                parameter_flags.append(f"Large temperature variation ({temp_range:.1f}°C)")
        
        # Check for outliers (simple z-score method)
        if n_samples > 3:
            z_scores = np.abs((data_array - mean_val) / std_val) if std_val > 0 else np.zeros_like(data_array)
            outliers = np.sum(z_scores > 3)
            if outliers > 0:
                parameter_flags.append(f"{outliers} potential outliers")
        
        quality_metrics[parameter] = {
            'n_samples': n_samples,
            'mean': mean_val,
            'std': std_val,
            'cv_percent': cv,
            'flags': parameter_flags
        }
        
        quality_flags.extend([f"{parameter}: {flag}" for flag in parameter_flags])
    
    # Environmental condition evaluation
    environmental_quality = "good"
    environmental_notes = []
    
    if 'weather_conditions' in measurement_metadata:
        weather = measurement_metadata['weather_conditions']
        if 'rain' in weather.lower() or 'storm' in weather.lower():
            environmental_quality = "poor"
            environmental_notes.append("Adverse weather conditions during sampling")
        elif 'wind' in weather.lower():
            environmental_quality = "fair"
            environmental_notes.append("Windy conditions may affect sampling")
    
    if 'water_temperature_c' in measurement_metadata:
        temp = measurement_metadata['water_temperature_c']
        if temp < 1 or temp > 35:
            environmental_notes.append(f"Extreme water temperature ({temp}°C)")
    
    # Equipment and procedural evaluation
    equipment_quality = "good"
    equipment_notes = []
    
    if 'calibration_date' in measurement_metadata:
        # Check if calibration is recent (within 6 months)
        cal_date = measurement_metadata['calibration_date']
        if isinstance(cal_date, str):
            try:
                cal_datetime = datetime.strptime(cal_date, '%Y-%m-%d')
                days_since_cal = (datetime.now() - cal_datetime).days
                if days_since_cal > 180:
                    equipment_quality = "fair"
                    equipment_notes.append("Equipment calibration >6 months old")
            except ValueError:
                equipment_notes.append("Invalid calibration date format")
    
    if 'sampler_type' in measurement_metadata:
        sampler = measurement_metadata['sampler_type']
        if 'unknown' in sampler.lower():
            equipment_quality = "poor"
            equipment_notes.append("Unknown sampler type")
    
    # Overall quality assessment
    total_flags = len(quality_flags)
    total_parameters = len(field_measurements)
    
    if total_flags == 0:
        overall_quality = "excellent"
    elif total_flags <= total_parameters * 0.2:
        overall_quality = "good"
    elif total_flags <= total_parameters * 0.5:
        overall_quality = "fair"
    else:
        overall_quality = "poor"
    
    # Generate recommendations
    recommendations = []
    
    if total_flags > 0:
        recommendations.append(f"Address {total_flags} quality issues identified")
    
    if environmental_quality == "poor":
        recommendations.append("Consider resampling under better environmental conditions")
    
    if equipment_quality == "fair":
        recommendations.append("Update equipment calibration")
    
    high_cv_params = [param for param, metrics in quality_metrics.items() 
                     if metrics['cv_percent'] > 25]
    if high_cv_params:
        recommendations.append(f"Investigate high variability in {high_cv_params}")
    
    if not recommendations:
        recommendations.append("Data quality appears acceptable for analysis")
    
    results = {
        'overall_quality': overall_quality,
        'environmental_quality': environmental_quality,
        'equipment_quality': equipment_quality,
        'total_quality_flags': total_flags,
        'quality_flags': quality_flags,
        'quality_metrics': quality_metrics,
        'data_completeness': data_completeness,
        'environmental_notes': environmental_notes,
        'equipment_notes': equipment_notes,
        'recommendations': recommendations,
        'measurement_metadata': measurement_metadata
    }
    
    return results


def field_measurement_protocol(
    measurement_type: str,
    site_conditions: Dict[str, float],
    equipment_list: List[str],
    safety_requirements: Optional[List[str]] = None
) -> Dict[str, Union[List[str], Dict, str]]:
    """
    Generate field measurement protocol based on conditions and objectives.
    
    Parameters
    ----------
    measurement_type : str
        Type of measurement: 'suspended_sediment', 'bed_material', 'bed_load'
    site_conditions : dict
        Site-specific conditions (depth, velocity, temperature, etc.)
    equipment_list : list
        Available equipment
    safety_requirements : list, optional
        Specific safety requirements
        
    Returns
    -------
    dict
        Dictionary containing detailed measurement protocol
    """
    if measurement_type not in ['suspended_sediment', 'bed_material', 'bed_load']:
        raise ValueError("measurement_type must be 'suspended_sediment', 'bed_material', or 'bed_load'")
    
    # Basic site condition validation
    depth = site_conditions.get('depth_m', 1.0)
    velocity = site_conditions.get('velocity_m_s', 0.5)
    temperature = site_conditions.get('temperature_c', 15.0)
    
    validate_positive(depth, "depth_m")
    validate_positive(velocity, "velocity_m_s")
    
    protocol = {
        'measurement_type': measurement_type,
        'pre_sampling_checklist': [],
        'sampling_procedure': [],
        'post_sampling_tasks': [],
        'safety_considerations': [],
        'equipment_requirements': [],
        'data_recording': [],
        'quality_control': []
    }
    
    # Pre-sampling checklist
    protocol['pre_sampling_checklist'].extend([
        "Check weather conditions and safety",
        "Verify equipment calibration dates", 
        "Prepare sample containers and labels",
        "Record site coordinates and time",
        "Measure water temperature and pH",
        "Document flow conditions and stage"
    ])
    
    # Measurement-specific procedures
    if measurement_type == 'suspended_sediment':
        protocol['sampling_procedure'].extend([
            "Select representative sampling verticals across channel",
            "Use depth-integrating sampler for each vertical",
            "Maintain constant transit rate (30-60 seconds per round trip)",
            "Collect minimum 500 mL sample volume",
            "Avoid disturbing bed during sampling",
            "Record exact sampling locations and times"
        ])
        
        protocol['equipment_requirements'].extend([
            "Depth-integrating suspended sediment sampler",
            "Sample bottles (500-1000 mL)",
            "Wading rod or cable-and-reel system",
            "Current meter for velocity measurement",
            "Field notebook and data sheets"
        ])
        
        if depth > 1.5:
            protocol['equipment_requirements'].append("Bridge board or boat for deep water sampling")
            protocol['safety_considerations'].append("Use safety rope for deep water operations")
        
        if velocity > 2.0:
            protocol['sampling_procedure'].append("Use heavier sampler (DH-76 or P-63) for high velocity")
            protocol['safety_considerations'].append("Exercise caution in high velocity conditions")
    
    elif measurement_type == 'bed_material':
        protocol['sampling_procedure'].extend([
            "Use systematic grid sampling pattern",
            "Sample surface layer (0-2 cm depth) for armor analysis", 
            "Collect subsurface samples (2-10 cm) for substrate analysis",
            "Minimum 100 particles for pebble counts",
            "Avoid sampling in flow separation zones",
            "Photograph sampling locations"
        ])
        
        protocol['equipment_requirements'].extend([
            "Shovel or scoop for bulk samples",
            "Sieve set for field screening",
            "Sample bags and permanent markers",
            "Measuring tape for grid layout",
            "Camera for site documentation"
        ])
        
    elif measurement_type == 'bed_load':
        protocol['sampling_procedure'].extend([
            "Position sampler on channel bed",
            "Sample for consistent time periods (2-10 minutes)",
            "Multiple samples across channel width",
            "Account for sampler efficiency (typically 0.5-1.5)",
            "Measure exact sampling area and time",
            "Handle samples carefully to avoid loss"
        ])
        
        protocol['equipment_requirements'].extend([
            "Bed load sampler (Helley-Smith or similar)",
            "Strong mesh sample bags",
            "Stopwatch for timing",
            "Scale for field weighing",
            "Protective containers for transport"
        ])
    
    # Site-specific modifications
    if depth < 0.5:
        protocol['safety_considerations'].append("Shallow water - watch for obstacles and uneven bottom")
        protocol['sampling_procedure'].append("Use lightweight sampling equipment for shallow conditions")
    
    if velocity < 0.1:
        protocol['quality_control'].append("Low velocity conditions - check for sample settling in equipment")
    
    if temperature < 5:
        protocol['safety_considerations'].extend([
            "Cold water conditions - wear appropriate protection",
            "Monitor for hypothermia symptoms",
            "Keep equipment from freezing"
        ])
    
    # Standard safety requirements
    default_safety = [
        "Wear personal flotation device (PFD)",
        "Use buddy system - never sample alone",
        "Carry communication device (radio/cell phone)",
        "Inform supervisor of sampling schedule",
        "Check local hazards (currents, debris, wildlife)"
    ]
    
    if safety_requirements:
        protocol['safety_considerations'].extend(safety_requirements)
    protocol['safety_considerations'].extend(default_safety)
    
    # Data recording requirements
    protocol['data_recording'].extend([
        "Site identification and coordinates",
        "Date, time, and weather conditions",
        "Water stage and flow conditions",
        "Sample identification numbers",
        "Equipment type and settings",
        "Sampling locations and methods",
        "Field observations and notes"
    ])
    
    # Quality control procedures
    protocol['quality_control'].extend([
        "Collect duplicate samples (10% of total)",
        "Include field blanks for contamination check",
        "Cross-check measurements between team members",
        "Document any deviations from protocol",
        "Perform immediate visual inspection of samples"
    ])
    
    # Post-sampling tasks
    protocol['post_sampling_tasks'].extend([
        "Label all samples with waterproof markers",
        "Complete chain of custody forms",
        "Store samples at appropriate temperature",
        "Clean and inspect equipment",
        "Upload field notes and photos",
        "Schedule laboratory analysis"
    ])
    
    return protocol


def equipment_maintenance_schedule(
    equipment_list: List[Dict[str, Union[str, float]]],
    usage_frequency: str = "monthly",
    environmental_conditions: str = "normal"
) -> Dict[str, Union[str, List[Dict], Dict]]:
    """
    Generate equipment maintenance schedule for sediment sampling equipment.
    
    Parameters
    ----------
    equipment_list : list of dict
        List of equipment with type, age, and usage information
    usage_frequency : str, optional
        Usage frequency: 'weekly', 'monthly', 'quarterly', 'annually'
    environmental_conditions : str, optional
        Environmental conditions: 'normal', 'harsh', 'marine'
        
    Returns
    -------
    dict
        Dictionary containing maintenance schedule and procedures
    """
    frequency_multipliers = {
        'weekly': 0.25,
        'monthly': 1.0,
        'quarterly': 3.0,
        'annually': 12.0
    }
    
    condition_multipliers = {
        'normal': 1.0,
        'harsh': 0.7,  # More frequent maintenance
        'marine': 0.5   # Much more frequent due to salt corrosion
    }
    
    if usage_frequency not in frequency_multipliers:
        raise ValueError("usage_frequency must be 'weekly', 'monthly', 'quarterly', or 'annually'")
    
    if environmental_conditions not in condition_multipliers:
        raise ValueError("environmental_conditions must be 'normal', 'harsh', or 'marine'")
    
    base_interval = frequency_multipliers[usage_frequency] * condition_multipliers[environmental_conditions]
    
    # Standard maintenance procedures by equipment type
    maintenance_procedures = {
        'suspended_sediment_sampler': {
            'daily': ['Rinse with clean water', 'Check nozzle for clogs', 'Inspect gaskets'],
            'weekly': ['Disassemble and clean thoroughly', 'Check moving parts', 'Lubricate hinges'],
            'monthly': ['Calibrate intake efficiency', 'Replace worn gaskets', 'Check weight and balance'],
            'quarterly': ['Professional inspection', 'Replace worn components', 'Efficiency testing'],
            'annually': ['Complete overhaul', 'Replace all gaskets and seals', 'Recertification']
        },
        'current_meter': {
            'daily': ['Clean propeller', 'Check battery charge', 'Verify readings'],
            'weekly': ['Calibrate against standard', 'Clean electrical contacts', 'Check cables'],
            'monthly': ['Replace batteries', 'Lubricate bearings', 'Check mounting hardware'],
            'quarterly': ['Professional calibration', 'Replace worn bearings', 'Update firmware'],
            'annually': ['Factory service', 'Complete electrical check', 'Accuracy verification']
        },
        'bed_load_sampler': {
            'daily': ['Empty and clean bag', 'Check frame integrity', 'Inspect mesh'],
            'weekly': ['Check mounting hardware', 'Inspect welds', 'Clean thoroughly'],
            'monthly': ['Replace mesh if worn', 'Check frame alignment', 'Lubricate moving parts'],
            'quarterly': ['Efficiency calibration', 'Structural inspection', 'Replace worn parts'],
            'annually': ['Complete rebuild if needed', 'Welding inspection', 'Recertification']
        },
        'pump_sampler': {
            'daily': ['Check pump operation', 'Clean intake lines', 'Test flow rates'],
            'weekly': ['Service pump motor', 'Clean filters', 'Check electrical connections'],
            'monthly': ['Replace pump seals', 'Calibrate flow rates', 'Service controls'],
            'quarterly': ['Professional pump service', 'Replace wear parts', 'Electrical inspection'],
            'annually': ['Complete system overhaul', 'Replace major components', 'Recertification']
        }
    }
    
    # Generate equipment-specific schedules
    equipment_schedules = []
    
    for equipment in equipment_list:
        eq_type = equipment.get('type', 'unknown')
        eq_age = equipment.get('age_years', 0)
        eq_id = equipment.get('id', f"equipment_{len(equipment_schedules)}")
        
        if eq_type in maintenance_procedures:
            procedures = maintenance_procedures[eq_type]
            
            # Adjust intervals based on age
            age_factor = 1.0 + (eq_age * 0.1)  # 10% more frequent per year of age
            
            schedule = {
                'equipment_id': eq_id,
                'equipment_type': eq_type,
                'age_years': eq_age,
                'next_daily': 'Every use',
                'next_weekly': f"Every {int(base_interval * 7 / age_factor)} days",
                'next_monthly': f"Every {int(base_interval * 30 / age_factor)} days", 
                'next_quarterly': f"Every {int(base_interval * 90 / age_factor)} days",
                'next_annually': f"Every {int(base_interval * 365 / age_factor)} days",
                'procedures': procedures,
                'priority': 'high' if eq_age > 5 else 'normal'
            }
            
            equipment_schedules.append(schedule)
    
    # Generate overall maintenance calendar
    maintenance_calendar = {
        'usage_frequency': usage_frequency,
        'environmental_conditions': environmental_conditions,
        'total_equipment_count': len(equipment_list),
        'high_priority_count': sum(1 for eq in equipment_schedules if eq.get('priority') == 'high'),
        'recommended_actions': []
    }
    
    # Generate recommendations
    recommendations = maintenance_calendar['recommended_actions']
    
    if environmental_conditions == 'marine':
        recommendations.extend([
            "Increase cleaning frequency due to salt exposure",
            "Use marine-grade lubricants and protective coatings",
            "Monitor for corrosion more frequently"
        ])
    
    if environmental_conditions == 'harsh':
        recommendations.extend([
            "Increase inspection frequency due to harsh conditions",
            "Stock extra replacement parts",
            "Consider equipment protection measures"
        ])
    
    old_equipment = [eq for eq in equipment_schedules if eq['age_years'] > 10]
    if old_equipment:
        recommendations.append(f"Consider replacement for {len(old_equipment)} equipment items >10 years old")
    
    if usage_frequency == 'weekly':
        recommendations.append("High usage frequency - consider backup equipment")
    
    results = {
        'maintenance_calendar': maintenance_calendar,
        'equipment_schedules': equipment_schedules,
        'general_recommendations': recommendations,
        'emergency_procedures': [
            "Always carry backup equipment for critical measurements",
            "Maintain emergency repair kit",
            "Know equipment suppliers for rapid replacement",
            "Document all maintenance activities"
        ]
    }
    
    return results


def measurement_uncertainty_analysis(
    measurement_data: Dict[str, np.ndarray],
    uncertainty_sources: Dict[str, float],
    confidence_level: float = 0.95
) -> Dict[str, Union[float, Dict, str]]:
    """
    Analyze measurement uncertainty for sediment transport data.
    
    Parameters
    ----------
    measurement_data : dict
        Dictionary containing measurement arrays
    uncertainty_sources : dict
        Dictionary of uncertainty source magnitudes (as percentages)
    confidence_level : float, optional
        Confidence level for uncertainty intervals, default 0.95
        
    Returns
    -------
    dict
        Dictionary containing uncertainty analysis results
    """
    validate_range(confidence_level, "confidence_level", 0.0, 1.0)
    
    # Standard uncertainty sources and typical values
    default_uncertainties = {
        'instrument_precision': 2.0,  # %
        'calibration_drift': 1.5,     # %
        'sampling_variability': 10.0, # %
        'operator_bias': 5.0,         # %
        'environmental_effects': 3.0, # %
        'sample_handling': 2.0        # %
    }
    
    # Merge with provided uncertainties
    uncertainties = default_uncertainties.copy()
    uncertainties.update(uncertainty_sources)
    
    uncertainty_results = {}
    
    for parameter, data in measurement_data.items():
        data_array = validate_array_like(data, f"measurement_data[{parameter}]")
        
        if len(data_array) == 0:
            continue
        
        mean_value = np.mean(data_array)
        std_value = np.std(data_array, ddof=1) if len(data_array) > 1 else 0
        
        # Calculate combined standard uncertainty
        # Assume uncertainties are independent and combine in quadrature
        combined_uncertainty_pct = np.sqrt(sum(u**2 for u in uncertainties.values()))
        combined_uncertainty_abs = (combined_uncertainty_pct / 100.0) * mean_value
        
        # Statistical uncertainty from replicate measurements
        statistical_uncertainty = std_value / np.sqrt(len(data_array))
        
        # Total uncertainty (systematic + statistical)
        total_uncertainty = np.sqrt(combined_uncertainty_abs**2 + statistical_uncertainty**2)
        total_uncertainty_pct = (total_uncertainty / mean_value) * 100 if mean_value != 0 else 0
        
        # Expanded uncertainty (coverage factor for desired confidence level)
        from scipy import stats
        alpha = 1 - confidence_level
        if len(data_array) > 1:
            # Use t-distribution for small samples
            coverage_factor = stats.t.ppf(1 - alpha/2, df=len(data_array)-1)
        else:
            # Use normal distribution for large samples or single measurement
            coverage_factor = stats.norm.ppf(1 - alpha/2)
        
        expanded_uncertainty = coverage_factor * total_uncertainty
        expanded_uncertainty_pct = (expanded_uncertainty / mean_value) * 100 if mean_value != 0 else 0
        
        # Uncertainty budget breakdown
        uncertainty_budget = {}
        total_variance = sum(((u/100.0 * mean_value)**2) for u in uncertainties.values())
        
        for source, uncertainty_pct in uncertainties.items():
            uncertainty_abs = (uncertainty_pct / 100.0) * mean_value
            contribution_pct = (uncertainty_abs**2 / total_variance) * 100 if total_variance > 0 else 0
            uncertainty_budget[source] = {
                'uncertainty_percent': uncertainty_pct,
                'uncertainty_absolute': uncertainty_abs,
                'contribution_to_total_percent': contribution_pct
            }
        
        # Quality indicators
        uncertainty_class = "A" if expanded_uncertainty_pct < 10 else \
                          "B" if expanded_uncertainty_pct < 20 else \
                          "C" if expanded_uncertainty_pct < 30 else "D"
        
        measurement_quality = "excellent" if uncertainty_class == "A" else \
                            "good" if uncertainty_class == "B" else \
                            "acceptable" if uncertainty_class == "C" else "poor"
        
        uncertainty_results[parameter] = {
            'mean_value': mean_value,
            'standard_deviation': std_value,
            'sample_size': len(data_array),
            'statistical_uncertainty': statistical_uncertainty,
            'combined_systematic_uncertainty': combined_uncertainty_abs,
            'combined_systematic_uncertainty_pct': combined_uncertainty_pct,
            'total_uncertainty': total_uncertainty,
            'total_uncertainty_pct': total_uncertainty_pct,
            'expanded_uncertainty': expanded_uncertainty,
            'expanded_uncertainty_pct': expanded_uncertainty_pct,
            'coverage_factor': coverage_factor,
            'confidence_level': confidence_level,
            'uncertainty_class': uncertainty_class,
            'measurement_quality': measurement_quality,
            'uncertainty_budget': uncertainty_budget,
            'confidence_interval_lower': mean_value - expanded_uncertainty,
            'confidence_interval_upper': mean_value + expanded_uncertainty
        }
    
    # Overall assessment
    avg_uncertainty = np.mean([result['expanded_uncertainty_pct'] 
                              for result in uncertainty_results.values() 
                              if result['expanded_uncertainty_pct'] < np.inf])
    
    overall_quality = "excellent" if avg_uncertainty < 10 else \
                     "good" if avg_uncertainty < 20 else \
                     "acceptable" if avg_uncertainty < 30 else "poor"
    
    # Identify dominant uncertainty sources
    all_contributions = {}
    for parameter_results in uncertainty_results.values():
        for source, budget in parameter_results['uncertainty_budget'].items():
            if source not in all_contributions:
                all_contributions[source] = []
            all_contributions[source].append(budget['contribution_to_total_percent'])
    
    dominant_sources = []
    for source, contributions in all_contributions.items():
        avg_contribution = np.mean(contributions)
        if avg_contribution > 20:  # More than 20% average contribution
            dominant_sources.append((source, avg_contribution))
    
    dominant_sources.sort(key=lambda x: x[1], reverse=True)
    
    results = {
        'uncertainty_analysis': uncertainty_results,
        'overall_quality': overall_quality,
        'average_uncertainty_percent': avg_uncertainty,
        'confidence_level': confidence_level,
        'dominant_uncertainty_sources': dominant_sources[:3],  # Top 3
        'uncertainty_sources_used': uncertainties,
        'recommendations': []
    }
    
    # Generate recommendations
    recommendations = results['recommendations']
    
    if avg_uncertainty > 25:
        recommendations.append("High measurement uncertainty - review sampling and analysis procedures")
    
    for source, contribution in dominant_sources[:2]:  # Top 2 sources
        if contribution > 30:
            recommendations.append(f"Address {source} - major contributor to uncertainty ({contribution:.1f}%)")
    
    if any(result['sample_size'] < 5 for result in uncertainty_results.values()):
        recommendations.append("Increase sample size to reduce statistical uncertainty")
    
    if not recommendations:
        recommendations.append("Uncertainty levels appear acceptable for intended use")
    
    return results