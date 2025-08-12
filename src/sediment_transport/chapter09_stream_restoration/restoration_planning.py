"""
Chapter 9: Stream Restoration - Restoration Planning Module

This module implements restoration project planning frameworks including
objective setting, team coordination, habitat assessment, and project scoping
based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration
- Federal Interagency Stream Restoration Working Group (FISRWG 1998)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class RestorationType(Enum):
    """Classification of restoration project types based on ASCE Manual 110 Table 9-1"""
    RESTORATION = "restoration"  # Return to predisturbance conditions
    REHABILITATION = "rehabilitation"  # Partial recovery of functions
    PRESERVATION = "preservation"  # Maintain current functions
    MITIGATION = "mitigation"  # Compensate for environmental damage
    NATURALIZATION = "naturalization"  # Establish dynamically stable systems
    CREATION = "creation"  # Form new system where none existed
    ENHANCEMENT = "enhancement"  # Improve existing environmental quality
    RECLAMATION = "reclamation"  # Change biophysical capacity


class ProjectScale(Enum):
    """Project scale classification affecting restoration objectives"""
    LOCAL = "local"  # 20-100 channel widths, single bends
    REACH = "reach"  # Multiple local measures over long reaches
    WATERSHED = "watershed"  # Systemic restoration affecting land use


class DegradationCause(Enum):
    """Typical causes of river corridor degradation from Table 9-2"""
    NATURAL_EVENTS = "natural_events"  # Floods, landslides, earthquakes
    LAND_USE_CHANGES = "land_use_changes"  # Urbanization, logging, grazing
    FLOW_REGULATION = "flow_regulation"  # Dams, withdrawals, diversions
    CHANNEL_MODIFICATIONS = "channel_modifications"  # Channelization, protection


@dataclass
class ProjectObjective:
    """Individual project objective with measurable criteria"""
    description: str
    parameter: str  # Physical parameter being targeted
    target_value: float
    units: str
    measurement_method: str
    success_criteria: str


@dataclass
class HabitatAssessment:
    """Habitat quality assessment and degradation factors"""
    current_quality_score: float  # 0-100 scale
    degradation_factors: List[DegradationCause]
    limiting_factors: List[str]
    target_species: List[str]
    habitat_requirements: Dict[str, Dict[str, float]]


def define_restoration_terminology() -> Dict[str, str]:
    """
    Define restoration terminology for stakeholder communication.
    
    Returns:
        Dictionary of restoration terms and definitions from Table 9-1
        
    Reference: ASCE Manual 110, Table 9-1
    """
    definitions = {
        "restoration": (
            "Reestablishment of the structure and function of ecosystems. "
            "Ecological restoration returns an ecosystem as closely as possible "
            "to predisturbance conditions and functions."
        ),
        "rehabilitation": (
            "Partial recovery of ecosystem functions and processes. "
            "Establishes geologically and hydrologically stable landscapes "
            "that support the natural ecosystem mosaic."
        ),
        "preservation": (
            "Activities to maintain current functions and characteristics "
            "of an ecosystem or to protect it from future damage."
        ),
        "mitigation": (
            "Activity to compensate for or alleviate environmental damage. "
            "May restore a site to socially acceptable condition, but not "
            "necessarily to a natural condition."
        ),
        "naturalization": (
            "Management aimed at establishing hydraulically and morphologically "
            "varied, yet dynamically stable fluvial systems capable of supporting "
            "healthy, biologically diverse aquatic ecosystems."
        ),
        "creation": (
            "Forming a new system where one did not formerly exist "
            "(e.g., constructing a wetland)."
        ),
        "enhancement": (
            "Subjective term for activities undertaken to improve existing "
            "environmental quality."
        ),
        "reclamation": (
            "Activities intended to change the biophysical capacity of an "
            "ecosystem. The resulting ecosystem differs from the pre-recovery state."
        )
    }
    return definitions


def assess_project_scale(channel_width: float, project_length: float, 
                        watershed_area: float) -> ProjectScale:
    """
    Determine project scale based on geometric characteristics.
    
    Args:
        channel_width: Representative channel width (m)
        project_length: Length of project reach (m)
        watershed_area: Contributing watershed area (km²)
        
    Returns:
        ProjectScale classification
        
    Reference: ASCE Manual 110, Section 9.1.2.1
    """
    validate_positive(channel_width, "channel_width")
    validate_positive(project_length, "project_length")
    validate_positive(watershed_area, "watershed_area")
    
    # Calculate project length in channel widths
    length_in_widths = project_length / channel_width
    
    if length_in_widths < 20:
        return ProjectScale.LOCAL
    elif length_in_widths <= 100 and watershed_area < 500:
        return ProjectScale.REACH
    else:
        return ProjectScale.WATERSHED


def evaluate_river_dynamism(flow_coefficient_variation: float,
                           flood_frequency_variance: float,
                           geomorphic_context: str) -> str:
    """
    Evaluate river dynamism characteristics for restoration planning.
    
    Args:
        flow_coefficient_variation: Coefficient of variation of annual flows
        flood_frequency_variance: Variance in flood frequency distribution
        geomorphic_context: Description of geomorphic setting
        
    Returns:
        Classification of system dynamics
        
    Reference: ASCE Manual 110, Section 9.1.2.2
    """
    validate_positive(flow_coefficient_variation, "flow_coefficient_variation")
    validate_positive(flood_frequency_variance, "flood_frequency_variance")
    
    if flow_coefficient_variation < 0.3 and flood_frequency_variance < 1.0:
        if "humid" in geomorphic_context.lower() or "mild relief" in geomorphic_context.lower():
            return "dynamic_equilibrium"
    
    if flow_coefficient_variation > 0.8 or flood_frequency_variance > 2.0:
        if any(term in geomorphic_context.lower() for term in ["arid", "semiarid", "proglacial"]):
            return "flood_dominated_transient"
    
    return "intermediate_dynamics"


def identify_team_disciplines() -> Dict[str, List[str]]:
    """
    Identify required disciplines for multidisciplinary restoration team.
    
    Returns:
        Dictionary of disciplines and their expertise areas
        
    Reference: ASCE Manual 110, Section 9.1.3.1
    """
    team_disciplines = {
        "sedimentation_engineer": [
            "sediment transport analysis",
            "hydraulic modeling",
            "channel stability assessment",
            "erosion and deposition prediction"
        ],
        "hydrologist": [
            "hydrologic modeling",
            "flow duration analysis",
            "flood frequency analysis",
            "water balance studies"
        ],
        "geomorphologist": [
            "channel evolution assessment",
            "landform analysis",
            "historical channel changes",
            "regional geomorphic context"
        ],
        "aquatic_ecologist": [
            "habitat assessment",
            "species requirements",
            "ecosystem function evaluation",
            "monitoring protocol design"
        ],
        "fisheries_biologist": [
            "fish population assessment",
            "spawning habitat requirements",
            "migration corridor needs",
            "water quality criteria"
        ],
        "riparian_ecologist": [
            "vegetation assessment",
            "revegetation planning",
            "wildlife habitat evaluation",
            "succession planning"
        ]
    }
    return team_disciplines


def set_project_objectives(general_goals: List[str], 
                          habitat_assessment: HabitatAssessment,
                          project_scale: ProjectScale) -> List[ProjectObjective]:
    """
    Convert general goals into specific, measurable objectives.
    
    Args:
        general_goals: List of broad project goals
        habitat_assessment: Current habitat conditions
        project_scale: Project scale classification
        
    Returns:
        List of specific project objectives with measurable criteria
        
    Reference: ASCE Manual 110, Section 9.1.3.2
    """
    objectives = []
    
    for goal in general_goals:
        if "water quality" in goal.lower():
            objectives.append(ProjectObjective(
                description="Maintain temperature suitable for aquatic organisms",
                parameter="daily_maximum_temperature",
                target_value=17.0,  # °C as specified in SRSRT 1994
                units="°C",
                measurement_method="continuous temperature monitoring",
                success_criteria="90% of summer days < 17°C"
            ))
            
        elif "habitat" in goal.lower():
            objectives.append(ProjectObjective(
                description="Increase pool volume for fish refuge",
                parameter="pool_volume_per_km",
                target_value=500.0,  # m³/km
                units="m³/km",
                measurement_method="topographic survey",
                success_criteria="25% increase from baseline"
            ))
            
        elif "flow" in goal.lower():
            objectives.append(ProjectObjective(
                description="Maintain adequate streamflow for aquatic life",
                parameter="minimum_flow",
                target_value=0.1,  # fraction of mean annual flow
                units="fraction of MAF",
                measurement_method="continuous flow monitoring",
                success_criteria="Flow > 0.1 MAF 95% of time"
            ))
    
    # Add scale-specific objectives
    if project_scale == ProjectScale.WATERSHED:
        objectives.append(ProjectObjective(
            description="Restore dynamic ecosystem processes",
            parameter="floodplain_connectivity",
            target_value=0.8,  # fraction of historic connectivity
            units="dimensionless",
            measurement_method="hydraulic modeling",
            success_criteria="80% of historic floodplain accessible"
        ))
    
    return objectives


def assess_habitat_degradation(physical_parameters: Dict[str, float],
                             biological_indicators: Dict[str, float],
                             degradation_causes: List[DegradationCause]) -> HabitatAssessment:
    """
    Assess current habitat quality and degradation factors.
    
    Args:
        physical_parameters: Dict of physical habitat parameters
        biological_indicators: Dict of biological condition indicators  
        degradation_causes: List of identified degradation causes
        
    Returns:
        HabitatAssessment object with quality score and limiting factors
        
    Reference: ASCE Manual 110, Section 9.1.3.2.1, Table 9-2
    """
    # Calculate habitat quality score (0-100)
    physical_score = 0.0
    bio_score = 0.0
    
    # Physical habitat scoring
    if "pool_frequency" in physical_parameters:
        pool_freq = physical_parameters["pool_frequency"]
        physical_score += min(pool_freq / 5.0, 1.0) * 25  # Target: 5 pools per km
    
    if "riffle_frequency" in physical_parameters:
        riffle_freq = physical_parameters["riffle_frequency"]
        physical_score += min(riffle_freq / 10.0, 1.0) * 25  # Target: 10 riffles per km
    
    if "bank_stability" in physical_parameters:
        stability = physical_parameters["bank_stability"]
        physical_score += stability * 25  # 0-1 scale
    
    if "substrate_quality" in physical_parameters:
        substrate = physical_parameters["substrate_quality"]
        physical_score += substrate * 25  # 0-1 scale
    
    # Biological indicator scoring
    if "fish_diversity" in biological_indicators:
        diversity = biological_indicators["fish_diversity"]
        bio_score += min(diversity / 10.0, 1.0) * 50  # Target: 10 species
    
    if "macroinvertebrate_score" in biological_indicators:
        macro_score = biological_indicators["macroinvertebrate_score"]
        bio_score += min(macro_score / 100.0, 1.0) * 50  # Target: 100 EPT score
    
    # Overall quality score
    quality_score = (physical_score + bio_score) / 2.0
    
    # Identify limiting factors based on degradation causes
    limiting_factors = []
    for cause in degradation_causes:
        if cause == DegradationCause.FLOW_REGULATION:
            limiting_factors.extend(["altered flow regime", "thermal stratification"])
        elif cause == DegradationCause.CHANNEL_MODIFICATIONS:
            limiting_factors.extend(["reduced complexity", "bank instability"])
        elif cause == DegradationCause.LAND_USE_CHANGES:
            limiting_factors.extend(["excess sedimentation", "nutrient loading"])
        elif cause == DegradationCause.NATURAL_EVENTS:
            limiting_factors.extend(["physical disturbance", "habitat fragmentation"])
    
    return HabitatAssessment(
        current_quality_score=quality_score,
        degradation_factors=degradation_causes,
        limiting_factors=limiting_factors,
        target_species=["salmonids", "native cyprinids", "EPT taxa"],
        habitat_requirements={
            "salmonids": {"temperature_max": 17.0, "dissolved_oxygen_min": 7.0, "velocity_range": [0.3, 1.5]},
            "native_cyprinids": {"temperature_max": 25.0, "dissolved_oxygen_min": 5.0, "depth_min": 0.3},
            "EPT_taxa": {"substrate_d50_mm": [2.0, 64.0], "velocity_max": 0.8, "embeddedness_max": 0.3}
        }
    )


def evaluate_project_constraints(project_scale: ProjectScale,
                               urban_proximity: bool,
                               budget_available: float,
                               regulatory_requirements: List[str]) -> Dict[str, Union[str, float, bool]]:
    """
    Evaluate project constraints affecting restoration options.
    
    Args:
        project_scale: Scale of restoration project
        urban_proximity: Whether project is near urban development
        budget_available: Available project budget ($)
        regulatory_requirements: List of applicable regulations
        
    Returns:
        Dictionary of constraint assessments
        
    Reference: ASCE Manual 110, Section 9.1.3.2.2
    """
    validate_positive(budget_available, "budget_available")
    
    constraints = {
        "scale_limitations": project_scale == ProjectScale.LOCAL,
        "urban_constraints": urban_proximity,
        "budget_adequate": budget_available > 100000,  # Minimum for meaningful restoration
        "permitting_complexity": len(regulatory_requirements) > 3,
        "design_flexibility": not urban_proximity and project_scale != ProjectScale.LOCAL
    }
    
    # Assess restoration strategy implications
    if urban_proximity:
        constraints["channel_boundary_fixed"] = True
        constraints["floodplain_access_limited"] = True
        constraints["natural_processes_constrained"] = True
    else:
        constraints["channel_boundary_fixed"] = False
        constraints["floodplain_access_limited"] = False
        constraints["natural_processes_constrained"] = False
    
    # Budget implications
    cost_per_meter = {
        ProjectScale.LOCAL: 1000,     # $1,000/m for local measures
        ProjectScale.REACH: 5000,     # $5,000/m for reach-scale
        ProjectScale.WATERSHED: 10000  # $10,000/m for watershed-scale
    }
    
    constraints["estimated_cost_per_meter"] = cost_per_meter[project_scale]
    constraints["maximum_project_length"] = budget_available / cost_per_meter[project_scale]
    
    return constraints


def develop_restoration_strategy(objectives: List[ProjectObjective],
                               habitat_assessment: HabitatAssessment,
                               constraints: Dict[str, Union[str, float, bool]],
                               project_scale: ProjectScale) -> Dict[str, List[str]]:
    """
    Develop restoration strategy based on objectives and constraints.
    
    Args:
        objectives: List of project objectives
        habitat_assessment: Current habitat conditions
        constraints: Project constraints
        project_scale: Scale of restoration project
        
    Returns:
        Dictionary of restoration strategies by category
        
    Reference: ASCE Manual 110, Section 9.1.3.2.3
    """
    strategy = {
        "channel_modifications": [],
        "flow_management": [],
        "habitat_enhancement": [],
        "monitoring_requirements": [],
        "risk_mitigation": []
    }
    
    # Channel modifications based on scale and constraints
    if project_scale == ProjectScale.WATERSHED and not constraints["urban_constraints"]:
        strategy["channel_modifications"].extend([
            "meandering_planform_restoration",
            "floodplain_reconnection",
            "natural_channel_evolution"
        ])
    elif project_scale == ProjectScale.REACH:
        strategy["channel_modifications"].extend([
            "riffle_pool_enhancement",
            "bank_stabilization",
            "in_stream_structures"
        ])
    else:  # LOCAL scale
        strategy["channel_modifications"].extend([
            "erosion_control_structures",
            "habitat_structures",
            "local_channel_modifications"
        ])
    
    # Flow management strategies
    if any("flow" in obj.description.lower() for obj in objectives):
        strategy["flow_management"].extend([
            "minimum_flow_maintenance",
            "flow_variability_restoration",
            "flood_pulse_management"
        ])
    
    # Habitat enhancement based on assessment
    if habitat_assessment.current_quality_score < 50:
        strategy["habitat_enhancement"].extend([
            "riparian_vegetation_restoration",
            "spawning_habitat_creation",
            "refuge_habitat_enhancement"
        ])
    
    # Monitoring requirements
    strategy["monitoring_requirements"] = [
        "pre_project_baseline_monitoring",
        "construction_phase_monitoring", 
        "post_project_performance_monitoring",
        "adaptive_management_protocols"
    ]
    
    # Risk mitigation
    if constraints["budget_adequate"] == False:
        strategy["risk_mitigation"].append("phased_implementation")
    if constraints["permitting_complexity"]:
        strategy["risk_mitigation"].append("regulatory_coordination")
    
    return strategy


def estimate_implementation_timeline(strategy: Dict[str, List[str]],
                                   project_scale: ProjectScale,
                                   regulatory_complexity: int) -> Dict[str, Tuple[int, int]]:
    """
    Estimate implementation timeline for restoration project phases.
    
    Args:
        strategy: Restoration strategy components
        project_scale: Scale of restoration project
        regulatory_complexity: Number of regulatory requirements (1-10)
        
    Returns:
        Dictionary of project phases with (min_months, max_months) estimates
        
    Reference: ASCE Manual 110, Section 9.2.7
    """
    validate_range(regulatory_complexity, 1, 10, "regulatory_complexity")
    
    # Base timeline factors
    scale_factors = {
        ProjectScale.LOCAL: 1.0,
        ProjectScale.REACH: 2.0, 
        ProjectScale.WATERSHED: 3.0
    }
    
    regulatory_delay = max(1.0, regulatory_complexity / 5.0)
    scale_factor = scale_factors[project_scale]
    
    timeline = {
        "planning_and_design": (
            int(6 * scale_factor * regulatory_delay),
            int(12 * scale_factor * regulatory_delay)
        ),
        "permitting_and_approvals": (
            int(3 * regulatory_delay),
            int(18 * regulatory_delay)
        ),
        "construction_implementation": (
            int(3 * scale_factor),
            int(12 * scale_factor)
        ),
        "establishment_monitoring": (
            int(12 * scale_factor),
            int(60 * scale_factor)
        )
    }
    
    return timeline