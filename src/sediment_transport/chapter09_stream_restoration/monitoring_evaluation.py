"""
Chapter 9: Stream Restoration - Monitoring and Evaluation Module

This module implements monitoring and evaluation frameworks for stream
restoration projects including performance metrics, adaptive management,
and long-term assessment protocols based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.8
- FISRWG (1998), Clary and Webster (1989)
- Seal et al. (1999)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass
from datetime import datetime, timedelta

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class MonitoringType(Enum):
    """Types of monitoring parameters"""
    PHYSICAL = "physical"              # Geomorphic and hydraulic
    BIOLOGICAL = "biological"          # Ecosystem response
    WATER_QUALITY = "water_quality"    # Chemical parameters
    SOCIAL = "social"                  # Stakeholder response


class PerformanceLevel(Enum):
    """Project performance assessment levels"""
    EXCELLENT = "excellent"           # Exceeds objectives
    GOOD = "good"                    # Meets objectives
    MARGINAL = "marginal"            # Partially meets objectives
    POOR = "poor"                    # Fails to meet objectives


class AdaptiveAction(Enum):
    """Types of adaptive management actions"""
    CONTINUE_MONITORING = "continue_monitoring"
    MINOR_MAINTENANCE = "minor_maintenance"
    DESIGN_MODIFICATION = "design_modification"
    MAJOR_INTERVENTION = "major_intervention"


@dataclass
class MonitoringParameter:
    """Definition of a monitoring parameter"""
    parameter_name: str
    measurement_units: str
    measurement_frequency: str
    target_value: Optional[float]
    acceptable_range: Optional[Tuple[float, float]]
    monitoring_method: str
    equipment_required: List[str]


@dataclass
class PerformanceMetric:
    """Performance assessment metric"""
    metric_name: str
    baseline_value: float
    target_value: float
    current_value: float
    trend: str                       # improving, stable, declining
    performance_level: PerformanceLevel
    years_to_target: Optional[int]


@dataclass
class MonitoringEvent:
    """Record of a monitoring event"""
    event_date: datetime
    parameters_measured: Dict[str, float]
    field_conditions: str
    equipment_used: List[str]
    data_quality: str
    notes: str


def develop_monitoring_plan(project_objectives: List[str],
                          project_scale: str,
                          budget_available: float,
                          monitoring_duration: int) -> Dict[str, List[MonitoringParameter]]:
    """
    Develop comprehensive monitoring plan based on project characteristics.
    
    Args:
        project_objectives: List of project objectives
        project_scale: "local", "reach", or "watershed"
        budget_available: Available monitoring budget ($)
        monitoring_duration: Planned monitoring duration (years)
        
    Returns:
        Dictionary of monitoring parameters by type
        
    Reference: ASCE Manual 110, Section 9.8
    """
    validate_positive(budget_available, "budget_available")
    validate_positive(monitoring_duration, "monitoring_duration")
    
    monitoring_plan = {
        "physical": [],
        "biological": [],
        "water_quality": [],
        "social": []
    }
    
    # Standard physical monitoring (always included)
    monitoring_plan["physical"].extend([
        MonitoringParameter(
            parameter_name="discharge",
            measurement_units="m³/s",
            measurement_frequency="continuous",
            target_value=None,
            acceptable_range=None,
            monitoring_method="stream_gauge_or_flow_measurement",
            equipment_required=["flow_meter", "stage_recorder"]
        ),
        MonitoringParameter(
            parameter_name="channel_width",
            measurement_units="m",
            measurement_frequency="annual",
            target_value=None,
            acceptable_range=None,
            monitoring_method="cross_sectional_survey",
            equipment_required=["survey_equipment", "GPS"]
        ),
        MonitoringParameter(
            parameter_name="channel_depth",
            measurement_units="m", 
            measurement_frequency="annual",
            target_value=None,
            acceptable_range=None,
            monitoring_method="cross_sectional_survey",
            equipment_required=["survey_equipment", "GPS"]
        ),
        MonitoringParameter(
            parameter_name="bed_material_d50",
            measurement_units="mm",
            measurement_frequency="annual",
            target_value=None,
            acceptable_range=None,
            monitoring_method="pebble_count_or_bulk_sampling",
            equipment_required=["sampling_grid", "sieves"]
        )
    ])
    
    # Add suspended sediment monitoring
    monitoring_plan["physical"].append(
        MonitoringParameter(
            parameter_name="suspended_sediment_concentration",
            measurement_units="mg/L",
            measurement_frequency="monthly",
            target_value=None,
            acceptable_range=None,
            monitoring_method="water_sampling_and_analysis",
            equipment_required=["water_samplers", "laboratory_analysis"]
        )
    )
    
    # Water quality monitoring (if objectives include habitat improvement)
    if any("habitat" in obj.lower() or "water quality" in obj.lower() for obj in project_objectives):
        monitoring_plan["water_quality"].extend([
            MonitoringParameter(
                parameter_name="temperature",
                measurement_units="°C",
                measurement_frequency="continuous",
                target_value=17.0,  # For coldwater species
                acceptable_range=(5.0, 20.0),
                monitoring_method="continuous_temperature_logger",
                equipment_required=["temperature_loggers"]
            ),
            MonitoringParameter(
                parameter_name="dissolved_oxygen",
                measurement_units="mg/L",
                measurement_frequency="weekly",
                target_value=8.0,
                acceptable_range=(6.0, 12.0),
                monitoring_method="handheld_DO_meter",
                equipment_required=["DO_meter", "calibration_solutions"]
            ),
            MonitoringParameter(
                parameter_name="turbidity",
                measurement_units="NTU",
                measurement_frequency="monthly",
                target_value=10.0,
                acceptable_range=(0.0, 25.0),
                monitoring_method="turbidity_meter",
                equipment_required=["turbidity_meter"]
            )
        ])
    
    # Biological monitoring (if objectives include ecosystem restoration)
    if any("habitat" in obj.lower() or "fish" in obj.lower() or "ecosystem" in obj.lower() 
           for obj in project_objectives):
        monitoring_plan["biological"].extend([
            MonitoringParameter(
                parameter_name="fish_species_richness",
                measurement_units="number_of_species",
                measurement_frequency="annual",
                target_value=None,
                acceptable_range=None,
                monitoring_method="electrofishing_survey",
                equipment_required=["electrofishing_equipment", "nets", "ID_guides"]
            ),
            MonitoringParameter(
                parameter_name="macroinvertebrate_EPT_index",
                measurement_units="EPT_taxa_count",
                measurement_frequency="annual",
                target_value=15.0,
                acceptable_range=(10.0, 25.0),
                monitoring_method="kick_net_sampling",
                equipment_required=["kick_nets", "sorting_trays", "microscopes"]
            ),
            MonitoringParameter(
                parameter_name="riparian_vegetation_cover",
                measurement_units="percent",
                measurement_frequency="annual",
                target_value=80.0,
                acceptable_range=(60.0, 100.0),
                monitoring_method="transect_surveys",
                equipment_required=["measuring_tapes", "plant_ID_guides"]
            )
        ])
    
    # Social monitoring (if community involvement is significant)
    if budget_available > 100000:  # Only for well-funded projects
        monitoring_plan["social"].append(
            MonitoringParameter(
                parameter_name="stakeholder_satisfaction",
                measurement_units="satisfaction_score",
                measurement_frequency="annual",
                target_value=4.0,  # Out of 5
                acceptable_range=(3.0, 5.0),
                monitoring_method="surveys_and_interviews",
                equipment_required=["survey_forms", "interview_protocols"]
            )
        )
    
    # Adjust monitoring intensity based on budget and scale
    if budget_available < 50000 or project_scale == "local":
        # Reduce monitoring frequency for budget-constrained projects
        for param_list in monitoring_plan.values():
            for param in param_list:
                if param.measurement_frequency == "continuous":
                    param.measurement_frequency = "monthly"
                elif param.measurement_frequency == "monthly":
                    param.measurement_frequency = "quarterly"
    
    return monitoring_plan


def calculate_performance_metrics(baseline_data: Dict[str, float],
                                target_values: Dict[str, float],
                                current_data: Dict[str, float],
                                years_since_construction: int) -> List[PerformanceMetric]:
    """
    Calculate performance metrics for restoration project assessment.
    
    Args:
        baseline_data: Pre-project baseline measurements
        target_values: Target values for each parameter
        current_data: Current measurements
        years_since_construction: Years since project completion
        
    Returns:
        List of PerformanceMetric objects
        
    Reference: ASCE Manual 110, Section 9.8
    """
    validate_positive(years_since_construction, "years_since_construction")
    
    performance_metrics = []
    
    for parameter in baseline_data.keys():
        if parameter in current_data and parameter in target_values:
            baseline = baseline_data[parameter]
            target = target_values[parameter]
            current = current_data[parameter]
            
            # Calculate progress toward target
            total_change_needed = target - baseline
            actual_change = current - baseline
            
            if total_change_needed != 0:
                progress_fraction = actual_change / total_change_needed
            else:
                progress_fraction = 1.0 if abs(current - target) < 0.1 * abs(target) else 0.0
            
            # Determine trend
            if years_since_construction > 1:
                # Would need historical data to determine trend
                # For now, assume based on progress
                if progress_fraction > 0.8:
                    trend = "improving"
                elif progress_fraction > 0.2:
                    trend = "stable" 
                else:
                    trend = "declining"
            else:
                trend = "insufficient_data"
            
            # Assess performance level
            if progress_fraction >= 1.0:
                performance_level = PerformanceLevel.EXCELLENT
            elif progress_fraction >= 0.8:
                performance_level = PerformanceLevel.GOOD
            elif progress_fraction >= 0.5:
                performance_level = PerformanceLevel.MARGINAL
            else:
                performance_level = PerformanceLevel.POOR
            
            # Estimate years to target (if not already achieved)
            if progress_fraction < 1.0 and progress_fraction > 0:
                years_to_target = int((1.0 - progress_fraction) / (progress_fraction / years_since_construction))
                years_to_target = min(years_to_target, 20)  # Cap at 20 years
            elif progress_fraction >= 1.0:
                years_to_target = 0
            else:
                years_to_target = None  # Unable to estimate
            
            performance_metrics.append(PerformanceMetric(
                metric_name=parameter,
                baseline_value=baseline,
                target_value=target,
                current_value=current,
                trend=trend,
                performance_level=performance_level,
                years_to_target=years_to_target
            ))
    
    return performance_metrics


def assess_adaptive_management_needs(performance_metrics: List[PerformanceMetric],
                                   project_objectives: List[str],
                                   years_since_construction: int,
                                   maintenance_history: List[str]) -> Dict[str, Union[AdaptiveAction, List[str]]]:
    """
    Assess need for adaptive management actions based on performance.
    
    Args:
        performance_metrics: List of performance metrics
        project_objectives: Original project objectives
        years_since_construction: Years since construction
        maintenance_history: History of maintenance actions
        
    Returns:
        Dictionary with recommended adaptive actions
        
    Reference: ASCE Manual 110, Section 9.8
    """
    validate_positive(years_since_construction, "years_since_construction")
    
    # Count performance levels
    performance_counts = {level: 0 for level in PerformanceLevel}
    for metric in performance_metrics:
        performance_counts[metric.performance_level] += 1
    
    total_metrics = len(performance_metrics)
    
    # Calculate overall performance score
    if total_metrics > 0:
        score = (performance_counts[PerformanceLevel.EXCELLENT] * 4 + 
                performance_counts[PerformanceLevel.GOOD] * 3 +
                performance_counts[PerformanceLevel.MARGINAL] * 2 +
                performance_counts[PerformanceLevel.POOR] * 1) / total_metrics
    else:
        score = 2.0  # Default marginal
    
    # Determine primary adaptive action
    if score >= 3.5:
        primary_action = AdaptiveAction.CONTINUE_MONITORING
    elif score >= 2.5:
        if years_since_construction < 3:
            primary_action = AdaptiveAction.CONTINUE_MONITORING  # Allow time for establishment
        else:
            primary_action = AdaptiveAction.MINOR_MAINTENANCE
    elif score >= 1.5:
        primary_action = AdaptiveAction.DESIGN_MODIFICATION
    else:
        primary_action = AdaptiveAction.MAJOR_INTERVENTION
    
    # Specific recommendations
    recommendations = []
    
    # Physical performance issues
    physical_issues = [m for m in performance_metrics if "channel" in m.metric_name.lower() 
                      and m.performance_level in [PerformanceLevel.POOR, PerformanceLevel.MARGINAL]]
    if physical_issues:
        recommendations.append("Assess channel stability and sediment transport")
        if years_since_construction > 2:
            recommendations.append("Consider structural modifications or maintenance")
    
    # Biological performance issues
    bio_issues = [m for m in performance_metrics if any(term in m.metric_name.lower() 
                  for term in ["fish", "invertebrate", "vegetation"])
                  and m.performance_level == PerformanceLevel.POOR]
    if bio_issues:
        recommendations.append("Evaluate habitat limiting factors")
        recommendations.append("Consider supplemental habitat enhancements")
    
    # Water quality issues
    wq_issues = [m for m in performance_metrics if any(term in m.metric_name.lower()
                 for term in ["temperature", "oxygen", "turbidity"])
                 and m.performance_level == PerformanceLevel.POOR]
    if wq_issues:
        recommendations.append("Investigate water quality sources")
        recommendations.append("Consider watershed-scale interventions")
    
    # Time-based recommendations
    if years_since_construction < 2:
        recommendations.append("Allow additional time for system establishment")
    elif years_since_construction > 10:
        recommendations.append("Evaluate long-term sustainability")
        recommendations.append("Plan for major maintenance or redesign")
    
    # Maintenance history considerations
    if "sediment_removal" in maintenance_history:
        recommendations.append("Evaluate upstream sediment sources")
    if "structure_repair" in maintenance_history:
        recommendations.append("Assess structure design adequacy")
    if len(maintenance_history) > 5:
        recommendations.append("High maintenance frequency indicates design issues")
    
    return {
        "primary_action": primary_action,
        "overall_performance_score": score,
        "specific_recommendations": recommendations,
        "priority_metrics": [m.metric_name for m in performance_metrics 
                           if m.performance_level == PerformanceLevel.POOR],
        "monitoring_adjustments": _recommend_monitoring_adjustments(performance_metrics, years_since_construction)
    }


def _recommend_monitoring_adjustments(performance_metrics: List[PerformanceMetric],
                                    years_since_construction: int) -> List[str]:
    """
    Recommend adjustments to monitoring program based on performance.
    
    Args:
        performance_metrics: Performance metrics
        years_since_construction: Years since construction
        
    Returns:
        List of monitoring adjustment recommendations
    """
    adjustments = []
    
    # Reduce monitoring for well-performing parameters
    excellent_metrics = [m for m in performance_metrics if m.performance_level == PerformanceLevel.EXCELLENT]
    if len(excellent_metrics) > 0 and years_since_construction > 3:
        adjustments.append("Reduce monitoring frequency for well-performing parameters")
    
    # Increase monitoring for poor performers
    poor_metrics = [m for m in performance_metrics if m.performance_level == PerformanceLevel.POOR]
    if poor_metrics:
        adjustments.append("Increase monitoring frequency for underperforming parameters")
        adjustments.append("Add process-level monitoring to understand failure mechanisms")
    
    # Long-term considerations
    if years_since_construction > 5:
        adjustments.append("Consider transition to long-term monitoring protocol")
        adjustments.append("Focus on sustainability indicators")
    
    return adjustments


def design_monitoring_schedule(monitoring_parameters: Dict[str, List[MonitoringParameter]],
                             project_duration: int,
                             seasonal_constraints: Optional[Dict[str, List[str]]] = None) -> Dict[str, List[Tuple[str, str]]]:
    """
    Design detailed monitoring schedule with timing considerations.
    
    Args:
        monitoring_parameters: Monitoring parameters by category
        project_duration: Project monitoring duration (years)
        seasonal_constraints: Optional seasonal timing constraints
        
    Returns:
        Dictionary with monitoring schedule by year
        
    Reference: ASCE Manual 110, Section 9.8
    """
    validate_positive(project_duration, "project_duration")
    
    schedule = {}
    
    for year in range(1, project_duration + 1):
        year_schedule = []
        
        # Physical monitoring
        for param in monitoring_parameters.get("physical", []):
            if param.measurement_frequency == "continuous":
                year_schedule.append((f"{param.parameter_name}_continuous", "year_round"))
            elif param.measurement_frequency == "monthly":
                year_schedule.extend([(f"{param.parameter_name}_monthly", f"month_{month}") 
                                     for month in range(1, 13)])
            elif param.measurement_frequency == "quarterly":
                year_schedule.extend([(f"{param.parameter_name}_quarterly", f"quarter_{q}") 
                                     for q in range(1, 5)])
            elif param.measurement_frequency == "annual":
                # Schedule annual monitoring based on parameter type
                if "bed_material" in param.parameter_name:
                    timing = "late_summer"  # Low flow period
                elif "vegetation" in param.parameter_name:
                    timing = "late_spring"  # After growing season start
                else:
                    timing = "summer"  # General field season
                year_schedule.append((f"{param.parameter_name}_annual", timing))
        
        # Biological monitoring (seasonal considerations)
        for param in monitoring_parameters.get("biological", []):
            if "fish" in param.parameter_name:
                timing = "late_summer"  # Optimal for fish sampling
            elif "invertebrate" in param.parameter_name:
                timing = "spring_fall"  # Avoid extreme temperatures
            elif "vegetation" in param.parameter_name:
                timing = "growing_season"
            else:
                timing = "field_season"
            year_schedule.append((f"{param.parameter_name}_biological", timing))
        
        # Water quality monitoring
        for param in monitoring_parameters.get("water_quality", []):
            if param.measurement_frequency == "continuous":
                year_schedule.append((f"{param.parameter_name}_continuous", "year_round"))
            elif param.measurement_frequency == "weekly":
                year_schedule.append((f"{param.parameter_name}_weekly", "year_round"))
            elif param.measurement_frequency == "monthly":
                year_schedule.extend([(f"{param.parameter_name}_monthly", f"month_{month}") 
                                     for month in range(1, 13)])
        
        # Apply seasonal constraints if provided
        if seasonal_constraints:
            filtered_schedule = []
            for task, timing in year_schedule:
                parameter_type = task.split("_")[0]
                if parameter_type in seasonal_constraints:
                    constraints = seasonal_constraints[parameter_type]
                    if timing in constraints or any(constraint in timing for constraint in constraints):
                        filtered_schedule.append((task, timing))
                else:
                    filtered_schedule.append((task, timing))
            year_schedule = filtered_schedule
        
        schedule[f"year_{year}"] = year_schedule
    
    return schedule


def estimate_monitoring_costs(monitoring_parameters: Dict[str, List[MonitoringParameter]],
                            project_duration: int,
                            site_characteristics: Dict[str, str]) -> Dict[str, float]:
    """
    Estimate monitoring program costs.
    
    Args:
        monitoring_parameters: Monitoring parameters by category
        project_duration: Monitoring duration (years)
        site_characteristics: Site access, location, etc.
        
    Returns:
        Dictionary with cost estimates by category
        
    Reference: Based on typical monitoring program costs
    """
    validate_positive(project_duration, "project_duration")
    
    costs = {
        "equipment": 0.0,
        "personnel": 0.0,
        "laboratory": 0.0,
        "data_management": 0.0,
        "reporting": 0.0,
        "annual_total": 0.0,
        "total_program": 0.0
    }
    
    # Equipment costs (one-time)
    equipment_costs = {
        "flow_meter": 3000,
        "stage_recorder": 5000,
        "temperature_loggers": 200,  # per logger
        "DO_meter": 1500,
        "turbidity_meter": 2000,
        "electrofishing_equipment": 15000,
        "survey_equipment": 10000,
        "GPS": 2000
    }
    
    all_equipment = []
    for param_list in monitoring_parameters.values():
        for param in param_list:
            all_equipment.extend(param.equipment_required)
    
    unique_equipment = list(set(all_equipment))
    for equipment in unique_equipment:
        if equipment in equipment_costs:
            costs["equipment"] += equipment_costs[equipment]
        else:
            costs["equipment"] += 1000  # Default cost
    
    # Annual personnel costs
    field_days_per_year = 0
    
    for param_list in monitoring_parameters.values():
        for param in param_list:
            if param.measurement_frequency == "continuous":
                field_days_per_year += 4  # Quarterly maintenance
            elif param.measurement_frequency == "monthly":
                field_days_per_year += 12
            elif param.measurement_frequency == "quarterly":
                field_days_per_year += 4
            elif param.measurement_frequency == "weekly":
                field_days_per_year += 52
            elif param.measurement_frequency == "annual":
                field_days_per_year += 2
    
    # Personnel cost per field day
    daily_rate = 800  # Including overhead, travel
    if site_characteristics.get("access", "moderate") == "difficult":
        daily_rate *= 1.5
    if site_characteristics.get("location", "local") == "remote":
        daily_rate *= 2.0
    
    costs["personnel"] = field_days_per_year * daily_rate
    
    # Laboratory costs
    lab_samples_per_year = 0
    for param_list in monitoring_parameters.values():
        for param in param_list:
            if "analysis" in param.monitoring_method.lower():
                if param.measurement_frequency == "monthly":
                    lab_samples_per_year += 12
                elif param.measurement_frequency == "quarterly":
                    lab_samples_per_year += 4
                elif param.measurement_frequency == "annual":
                    lab_samples_per_year += 1
    
    costs["laboratory"] = lab_samples_per_year * 50  # $50 per sample average
    
    # Data management and reporting
    costs["data_management"] = 5000  # Annual database and QA/QC
    costs["reporting"] = 10000  # Annual report preparation
    
    # Calculate totals
    costs["annual_total"] = (costs["personnel"] + costs["laboratory"] + 
                           costs["data_management"] + costs["reporting"])
    costs["total_program"] = costs["equipment"] + costs["annual_total"] * project_duration
    
    return costs