"""
Chapter 9: Stream Restoration - Stability Assessment Module

This module implements channel stability assessment methods including
qualitative and quantitative approaches, incipient motion analysis, 
stream power analysis, and bank stability assessment for stream restoration
projects based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.4
- Thorne (1999), Simon et al. (2000, 2003)
- Lane (1955), Shields (1936), Meyer-Peter and Muller (1948)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass
import warnings

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class AssessmentType(Enum):
    """Types of stability assessment approaches"""
    QUALITATIVE = "qualitative"           # Visual, reconnaissance
    SEMI_QUANTITATIVE = "semi_quantitative"  # Empirical relationships
    QUANTITATIVE = "quantitative"         # Detailed hydraulic analysis


class ChannelEvolutionStage(Enum):
    """Channel evolution model stages (Simon 1989)"""
    STAGE_I = "pre_modification"         # Stable, pre-disturbance
    STAGE_II = "construction"            # Active construction/disturbance  
    STAGE_III = "degradation"            # Active incision
    STAGE_IV = "degradation_widening"    # Incision with bank retreat
    STAGE_V = "aggradation_widening"     # Bed aggradation, continued widening
    STAGE_VI = "quasi_equilibrium"       # New equilibrium condition


class BedMaterial(Enum):
    """Classification of bed material types"""
    COBBLE_BOULDER = "cobble_boulder"    # >64 mm
    GRAVEL_COBBLE = "gravel_cobble"      # 2-256 mm  
    SAND = "sand"                        # 0.062-2 mm
    SILT_CLAY = "silt_clay"             # <0.062 mm


@dataclass
class StabilityMetrics:
    """Collection of stability assessment metrics"""
    stream_power: float                  # W/m²
    shear_stress: float                 # N/m²
    shields_parameter: float            # dimensionless
    velocity: float                     # m/s
    froude_number: float               # dimensionless
    critical_values: Dict[str, float]   # threshold values


@dataclass
class IncipientMotionResult:
    """Results of incipient motion analysis"""
    critical_grain_size: float         # mm
    method_used: str
    safety_factor: float               # actual/critical
    stability_assessment: str          # stable/marginal/unstable
    limitations: List[str]


@dataclass
class BankStabilityResult:
    """Results of bank stability analysis"""
    factor_of_safety: float           # >1.0 = stable
    critical_height: float            # m
    failure_mode: str                 # rotational, planar, etc.
    stability_class: str              # stable/marginal/unstable
    recommendations: List[str]


def calculate_stream_power(discharge: float,
                         channel_width: float,
                         channel_slope: float,
                         water_density: float = 1000.0) -> float:
    """
    Calculate unit stream power for stability assessment.
    
    Args:
        discharge: Flow discharge (m³/s)
        channel_width: Channel width (m)
        channel_slope: Channel slope (dimensionless)
        water_density: Water density (kg/m³)
        
    Returns:
        Unit stream power (W/m²)
        
    Reference: ASCE Manual 110, Table 9-8
    """
    validate_positive(discharge, "discharge")
    validate_positive(channel_width, "channel_width")
    validate_positive(channel_slope, "channel_slope")
    
    # Unit stream power = γ * Q * S / W
    # where γ = specific weight of water (N/m³)
    gamma_w = water_density * GRAVITY  # N/m³
    
    unit_stream_power = (gamma_w * discharge * channel_slope) / channel_width
    
    return unit_stream_power


def calculate_shear_stress(discharge: float,
                         channel_width: float,
                         channel_depth: float,
                         channel_slope: float,
                         water_density: float = 1000.0) -> float:
    """
    Calculate average bed shear stress.
    
    Args:
        discharge: Flow discharge (m³/s)
        channel_width: Channel width (m)
        channel_depth: Channel depth (m)
        channel_slope: Channel slope (dimensionless)
        water_density: Water density (kg/m³)
        
    Returns:
        Average bed shear stress (N/m²)
        
    Reference: ASCE Manual 110, Section 9.4.3.4
    """
    validate_positive(discharge, "discharge")
    validate_positive(channel_width, "channel_width")
    validate_positive(channel_depth, "channel_depth")
    validate_positive(channel_slope, "channel_slope")
    
    # Calculate hydraulic radius
    area = channel_width * channel_depth
    wetted_perimeter = channel_width + 2 * channel_depth
    hydraulic_radius = area / wetted_perimeter
    
    # Average shear stress = γ * R * S
    gamma_w = water_density * GRAVITY
    shear_stress = gamma_w * hydraulic_radius * channel_slope
    
    return shear_stress


def calculate_shields_parameter(shear_stress: float,
                              grain_size: float,
                              sediment_density: float = 2650.0,
                              water_density: float = 1000.0) -> float:
    """
    Calculate Shields dimensionless parameter.
    
    Args:
        shear_stress: Bed shear stress (N/m²)
        grain_size: Representative grain size D50 (m)
        sediment_density: Sediment particle density (kg/m³)
        water_density: Water density (kg/m³)
        
    Returns:
        Shields parameter (dimensionless)
        
    Reference: ASCE Manual 110, Table 9-6
    """
    validate_positive(shear_stress, "shear_stress")
    validate_positive(grain_size, "grain_size")
    
    gamma_s = sediment_density * GRAVITY
    gamma_w = water_density * GRAVITY
    
    # Shields parameter = τ / [(γs - γw) * D]
    shields_param = shear_stress / ((gamma_s - gamma_w) * grain_size)
    
    return shields_param


def assess_incipient_motion_lane(discharge: float,
                               channel_width: float,
                               channel_depth: float,
                               channel_slope: float,
                               concentration: float = 0.0) -> IncipientMotionResult:
    """
    Assess incipient motion using Lane (1955) method for fine sediments.
    
    Args:
        discharge: Flow discharge (m³/s)
        channel_width: Channel width (m)
        channel_depth: Channel depth (m)  
        channel_slope: Channel slope (dimensionless)
        concentration: Suspended sediment concentration (mg/L)
        
    Returns:
        IncipientMotionResult with critical grain size and assessment
        
    Reference: ASCE Manual 110, Table 9-6, Figure 9-8
    """
    # Calculate hydraulic parameters
    shear_stress = calculate_shear_stress(discharge, channel_width, 
                                        channel_depth, channel_slope)
    
    # Lane relationship: τc = γw * H * S = 13,000 * H * S * Dcritical (for Dcrit > 6mm)
    # Rearranging: Dcritical = τc / (13,000 * H * S)
    H = channel_depth
    S = channel_slope
    
    # Apply Lane's empirical relationship
    if concentration < 1000:  # mg/L, clear water
        critical_grain_size_mm = shear_stress / (13000 * H * S)
    elif concentration < 2000:  # moderate sediment load
        critical_grain_size_mm = shear_stress / (10000 * H * S)  
    else:  # heavy sediment load
        critical_grain_size_mm = shear_stress / (8000 * H * S)
    
    # Adjust for channel sinuosity (assume moderately sinuous)
    critical_grain_size_mm *= 0.75  # 25% reduction per USDA guidance
    
    # Safety assessment
    # Note: This method gives the critical size that would be stable
    # For assessment, we need to compare with actual bed material
    safety_factor = 1.0  # Default - would need actual bed material size
    
    if critical_grain_size_mm < 0.1:
        stability = "unstable"
        assessment_note = "Very fine material, high erosion potential"
    elif critical_grain_size_mm < 1.0:
        stability = "marginal"  
        assessment_note = "Fine material, moderate erosion risk"
    else:
        stability = "stable"
        assessment_note = "Coarse material can be retained"
    
    limitations = [
        "Based on observations from straight canals",
        "Not applicable where much sand is carried",
        "Reduced values needed for curved channels",
        "Limited to non-cohesive materials <5mm"
    ]
    
    return IncipientMotionResult(
        critical_grain_size=critical_grain_size_mm,
        method_used="Lane (1955)",
        safety_factor=safety_factor,
        stability_assessment=stability,
        limitations=limitations
    )


def assess_incipient_motion_shields(discharge: float,
                                  channel_width: float,
                                  channel_depth: float,
                                  channel_slope: float,
                                  bed_d50: float) -> IncipientMotionResult:
    """
    Assess incipient motion using Shields (1936) method.
    
    Args:
        discharge: Flow discharge (m³/s)
        channel_width: Channel width (m)
        channel_depth: Channel depth (m)
        channel_slope: Channel slope (dimensionless)
        bed_d50: Median grain size of bed material (mm)
        
    Returns:
        IncipientMotionResult with critical assessment
        
    Reference: ASCE Manual 110, Table 9-6
    """
    # Calculate hydraulic parameters
    shear_stress = calculate_shear_stress(discharge, channel_width,
                                        channel_depth, channel_slope)
    
    # Convert bed_d50 from mm to m
    d50_m = bed_d50 / 1000.0
    
    # Calculate Shields parameter
    shields_param = calculate_shields_parameter(shear_stress, d50_m)
    
    # Critical Shields parameter (typically 0.06 for coarse sediments)
    shields_critical = 0.06
    
    # Calculate critical grain size that would be stable
    gamma_s = 2650 * GRAVITY  # N/m³
    gamma_w = 1000 * GRAVITY  # N/m³
    
    critical_grain_size_m = shear_stress / (shields_critical * (gamma_s - gamma_w))
    critical_grain_size_mm = critical_grain_size_m * 1000
    
    # Safety factor
    safety_factor = shields_critical / shields_param if shields_param > 0 else float('inf')
    
    # Stability assessment
    if shields_param > shields_critical * 1.5:
        stability = "unstable"
    elif shields_param > shields_critical:
        stability = "marginal"
    else:
        stability = "stable"
    
    limitations = [
        "Assumes Shields constant = 0.06",
        "Most applicable to coarse, non-cohesive sediments",
        "Does not account for hiding effects in mixed sizes"
    ]
    
    return IncipientMotionResult(
        critical_grain_size=critical_grain_size_mm,
        method_used="Shields (1936)",
        safety_factor=safety_factor,
        stability_assessment=stability,
        limitations=limitations
    )


def assess_channel_evolution_stage(channel_width_ratio: float,
                                 channel_depth_ratio: float,
                                 bank_height: float,
                                 bank_angle: float,
                                 incision_rate: float) -> ChannelEvolutionStage:
    """
    Assess channel evolution stage using Simon (1989) model.
    
    Args:
        channel_width_ratio: Current/historic width ratio
        channel_depth_ratio: Current/historic depth ratio  
        bank_height: Bank height (m)
        bank_angle: Bank angle (degrees from horizontal)
        incision_rate: Recent incision rate (m/year)
        
    Returns:
        ChannelEvolutionStage classification
        
    Reference: ASCE Manual 110, Table 9-8
    """
    # Stage determination logic based on channel characteristics
    if abs(incision_rate) < 0.01 and channel_width_ratio < 1.2:
        # Minimal change, stable condition
        if channel_depth_ratio < 1.1:
            return ChannelEvolutionStage.STAGE_I  # Pre-modification
        else:
            return ChannelEvolutionStage.STAGE_VI  # New equilibrium
    
    elif incision_rate < -0.05:  # Significant downcutting
        if channel_width_ratio < 1.5:
            return ChannelEvolutionStage.STAGE_III  # Active degradation
        else:
            return ChannelEvolutionStage.STAGE_IV  # Degradation + widening
    
    elif incision_rate > 0.02:  # Aggradation
        return ChannelEvolutionStage.STAGE_V  # Aggradation + widening
    
    elif channel_width_ratio > 2.0:  # Significant widening
        if bank_angle > 70:  # Near vertical banks
            return ChannelEvolutionStage.STAGE_IV  # Active widening
        else:
            return ChannelEvolutionStage.STAGE_V  # Transitioning to equilibrium
    
    else:
        return ChannelEvolutionStage.STAGE_II  # Construction/transitional


def assess_bank_stability_simple(bank_height: float,
                               bank_angle: float,
                               cohesion: float,
                               friction_angle: float,
                               unit_weight: float = 18.0) -> BankStabilityResult:
    """
    Simple bank stability analysis using infinite slope method.
    
    Args:
        bank_height: Height of bank (m)
        bank_angle: Bank angle (degrees from horizontal)
        cohesion: Soil cohesion (kPa)
        friction_angle: Internal friction angle (degrees)
        unit_weight: Soil unit weight (kN/m³)
        
    Returns:
        BankStabilityResult with factor of safety
        
    Reference: ASCE Manual 110, Section 9.4.3.5
    """
    validate_positive(bank_height, "bank_height")
    validate_range(bank_angle, "bank_angle", 0.0, 90.0)
    validate_positive(friction_angle, "friction_angle")
    
    # Convert angles to radians
    beta = np.radians(bank_angle)
    phi = np.radians(friction_angle)
    
    # Infinite slope stability analysis
    # Factor of safety = (c + σ * tan(φ)) / (σ * tan(β))
    # where σ = γ * h * cos²(β) for normal stress
    
    normal_stress = unit_weight * bank_height * np.cos(beta)**2
    shear_stress = unit_weight * bank_height * np.sin(beta) * np.cos(beta)
    
    # Shear strength
    shear_strength = cohesion + normal_stress * np.tan(phi)
    
    # Factor of safety
    if shear_stress > 0:
        factor_of_safety = shear_strength / shear_stress
    else:
        factor_of_safety = float('inf')
    
    # Critical height for vertical banks (c / γ)
    if cohesion > 0:
        critical_height = (4 * cohesion) / unit_weight  # Approximate for vertical cut
    else:
        critical_height = 0
    
    # Stability classification
    if factor_of_safety > 1.5:
        stability_class = "stable"
        failure_mode = "none expected"
    elif factor_of_safety > 1.2:
        stability_class = "marginal"
        failure_mode = "possible shallow failures"
    else:
        stability_class = "unstable"
        if bank_angle > 60:
            failure_mode = "rotational failure likely"
        else:
            failure_mode = "shallow slope failures"
    
    # Recommendations
    recommendations = []
    if factor_of_safety < 1.3:
        recommendations.append("Consider bank protection measures")
        if bank_height > 3.0:
            recommendations.append("Detailed geotechnical analysis recommended")
    
    if bank_angle > 45 and cohesion < 10:
        recommendations.append("Consider bank slope reduction")
    
    if factor_of_safety < 1.0:
        recommendations.append("Immediate stabilization required")
    
    return BankStabilityResult(
        factor_of_safety=factor_of_safety,
        critical_height=critical_height,
        failure_mode=failure_mode,
        stability_class=stability_class,
        recommendations=recommendations
    )


def calculate_stability_metrics(discharge: float,
                              channel_width: float,
                              channel_depth: float,
                              channel_slope: float,
                              bed_d50: float) -> StabilityMetrics:
    """
    Calculate comprehensive stability metrics for assessment.
    
    Args:
        discharge: Flow discharge (m³/s)
        channel_width: Channel width (m)
        channel_depth: Channel depth (m)
        channel_slope: Channel slope (dimensionless)
        bed_d50: Median grain size (mm)
        
    Returns:
        StabilityMetrics with all calculated parameters
        
    Reference: ASCE Manual 110, Section 9.4
    """
    # Calculate basic hydraulic parameters
    area = channel_width * channel_depth
    velocity = discharge / area
    
    # Calculate stability metrics
    stream_power = calculate_stream_power(discharge, channel_width, channel_slope)
    shear_stress = calculate_shear_stress(discharge, channel_width, 
                                        channel_depth, channel_slope)
    shields_param = calculate_shields_parameter(shear_stress, bed_d50/1000.0)
    froude_number = velocity / np.sqrt(GRAVITY * channel_depth)
    
    # Critical values for assessment
    critical_values = {
        "stream_power_threshold": 35.0,      # W/m² (Brookes 1990)
        "shear_stress_sand": 2.0,            # N/m² for sand
        "shear_stress_gravel": 20.0,         # N/m² for gravel
        "shields_critical": 0.06,            # Dimensionless
        "froude_critical": 1.0,              # Dimensionless
        "velocity_erosion": 1.0              # m/s for fine sediments
    }
    
    return StabilityMetrics(
        stream_power=stream_power,
        shear_stress=shear_stress,
        shields_parameter=shields_param,
        velocity=velocity,
        froude_number=froude_number,
        critical_values=critical_values
    )


def comprehensive_stability_assessment(discharge: float,
                                     channel_width: float,
                                     channel_depth: float,
                                     channel_slope: float,
                                     bed_d50: float,
                                     bed_material_type: BedMaterial,
                                     bank_height: Optional[float] = None,
                                     bank_properties: Optional[Dict[str, float]] = None) -> Dict[str, Union[str, float, List[str]]]:
    """
    Comprehensive stability assessment using multiple methods.
    
    Args:
        discharge: Design discharge (m³/s)
        channel_width: Channel width (m)
        channel_depth: Channel depth (m)
        channel_slope: Channel slope (dimensionless)
        bed_d50: Median grain size (mm)
        bed_material_type: Classification of bed material
        bank_height: Bank height (m), optional
        bank_properties: Dict with cohesion, friction_angle, unit_weight
        
    Returns:
        Dictionary with comprehensive assessment results
        
    Reference: ASCE Manual 110, Table 9-8
    """
    results = {
        "overall_assessment": "",
        "stability_score": 0.0,
        "critical_issues": [],
        "recommendations": [],
        "confidence_level": "moderate"
    }
    
    # Calculate stability metrics
    metrics = calculate_stability_metrics(discharge, channel_width, 
                                        channel_depth, channel_slope, bed_d50)
    
    stability_score = 0.0
    issues = []
    recommendations = []
    
    # Stream power assessment
    if metrics.stream_power < metrics.critical_values["stream_power_threshold"]:
        stability_score += 2
    elif metrics.stream_power < metrics.critical_values["stream_power_threshold"] * 1.5:
        stability_score += 1
        issues.append("Elevated stream power may cause erosion")
    else:
        issues.append("High stream power - significant erosion potential")
        recommendations.append("Consider energy dissipation measures")
    
    # Shear stress assessment
    if bed_material_type == BedMaterial.SAND:
        threshold = metrics.critical_values["shear_stress_sand"]
    elif bed_material_type in [BedMaterial.GRAVEL_COBBLE, BedMaterial.COBBLE_BOULDER]:
        threshold = metrics.critical_values["shear_stress_gravel"]
    else:  # Silt/clay
        threshold = 15.0  # N/m² for cohesive materials
    
    if metrics.shear_stress < threshold:
        stability_score += 2
    elif metrics.shear_stress < threshold * 1.2:
        stability_score += 1
        issues.append("Shear stress approaches critical values")
    else:
        issues.append("Shear stress exceeds critical values for bed material")
        recommendations.append("Consider bed protection or coarser material")
    
    # Shields parameter assessment (for non-cohesive materials)
    if bed_material_type != BedMaterial.SILT_CLAY:
        if metrics.shields_parameter < metrics.critical_values["shields_critical"]:
            stability_score += 2
        elif metrics.shields_parameter < metrics.critical_values["shields_critical"] * 1.2:
            stability_score += 1
        else:
            issues.append("Shields parameter indicates bed material movement")
    
    # Froude number assessment
    if metrics.froude_number < 0.5:
        stability_score += 1  # Subcritical, stable
    elif metrics.froude_number < 1.0:
        # Transitional flow regime
        issues.append("Froude number indicates transitional flow conditions")
    else:
        issues.append("Supercritical flow - high erosion potential")
        recommendations.append("Consider flow control structures")
    
    # Bank stability assessment (if data provided)
    if bank_height is not None and bank_properties is not None:
        bank_result = assess_bank_stability_simple(
            bank_height=bank_height,
            bank_angle=bank_properties.get("bank_angle", 45),
            cohesion=bank_properties.get("cohesion", 10),
            friction_angle=bank_properties.get("friction_angle", 30)
        )
        
        if bank_result.factor_of_safety > 1.3:
            stability_score += 1
        elif bank_result.factor_of_safety < 1.0:
            issues.append(f"Bank instability: {bank_result.failure_mode}")
            recommendations.extend(bank_result.recommendations)
    
    # Overall assessment
    max_score = 6 if bank_height is not None else 5
    stability_score_normalized = stability_score / max_score
    
    if stability_score_normalized > 0.8:
        overall_assessment = "STABLE - Channel appears stable under design conditions"
        confidence_level = "high"
    elif stability_score_normalized > 0.6:
        overall_assessment = "MARGINALLY STABLE - Some stability concerns identified"
        confidence_level = "moderate"
    else:
        overall_assessment = "UNSTABLE - Significant stability issues identified"
        confidence_level = "low"
        recommendations.append("Detailed stability analysis recommended")
        recommendations.append("Consider alternative design approaches")
    
    # Add general recommendations
    if len(issues) > 2:
        recommendations.append("Multiple stability issues - comprehensive redesign needed")
    
    if bed_material_type == BedMaterial.SILT_CLAY:
        recommendations.append("Cohesive bed analysis may require specialized methods")
        confidence_level = "low"
    
    results.update({
        "overall_assessment": overall_assessment,
        "stability_score": stability_score_normalized,
        "critical_issues": issues,
        "recommendations": recommendations,
        "confidence_level": confidence_level,
        "metrics": {
            "stream_power": metrics.stream_power,
            "shear_stress": metrics.shear_stress,
            "shields_parameter": metrics.shields_parameter,
            "froude_number": metrics.froude_number,
            "velocity": metrics.velocity
        }
    })
    
    return results


def select_assessment_tools(project_complexity: str,
                          data_availability: Dict[str, bool],
                          bed_material_type: BedMaterial,
                          channel_conditions: str) -> List[str]:
    """
    Select appropriate stability assessment tools based on project characteristics.
    
    Args:
        project_complexity: "simple", "moderate", or "complex"
        data_availability: Dict indicating available data types
        bed_material_type: Classification of bed material
        channel_conditions: "stable", "unstable", or "unknown"
        
    Returns:
        List of recommended assessment tools
        
    Reference: ASCE Manual 110, Section 9.4.4
    """
    tools = []
    
    # Always start with basic assessments
    if project_complexity in ["simple", "moderate", "complex"]:
        tools.extend(["lane_relations", "stream_power_analysis"])
    
    # Lane-type relations for preliminary assessment
    if channel_conditions == "stable" or data_availability.get("basic_geometry", False):
        tools.append("hydraulic_geometry_relations")
    
    # Process-based classification for disturbed systems  
    if channel_conditions == "unstable":
        tools.append("channel_evolution_model")
    
    # Incipient motion analysis for coarse beds
    if bed_material_type in [BedMaterial.GRAVEL_COBBLE, BedMaterial.COBBLE_BOULDER]:
        if data_availability.get("bed_material_gradation", False):
            tools.append("shields_parameter_analysis")
            tools.append("meyer_peter_muller_analysis")
    
    # Sediment budgets for sand-bed streams
    if (bed_material_type == BedMaterial.SAND and 
        data_availability.get("sediment_data", False)):
        tools.append("sediment_budget_analysis")
    
    # Bank stability for high banks
    if data_availability.get("bank_properties", False):
        tools.append("bank_stability_analysis")
    
    # Complex assessments for detailed studies
    if project_complexity == "complex":
        tools.extend([
            "multiple_assessment_methods",
            "uncertainty_analysis",
            "sensitivity_analysis"
        ])
        
        if data_availability.get("detailed_hydraulics", False):
            tools.append("numerical_modeling")
    
    # Empirical tools for cohesive boundaries
    if bed_material_type == BedMaterial.SILT_CLAY:
        tools.extend([
            "empirical_velocity_criteria",
            "regional_stream_power_indices",
            "slope_area_relations"
        ])
    
    return tools