"""
Chapter 9: Stream Restoration - Sediment Studies Plan Module

This module implements sediment studies plan methodology for stream restoration
projects including study area boundary determination, stability assessment planning,
data inventory, and resource estimation based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.2
- Federal Interagency Sedimentation Project (2005)
- Edwards and Glysson (1988)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class StudyComplexity(Enum):
    """Complexity level of sediment studies plan"""
    SIMPLE = "simple"      # Qualitative assessment, minimal data
    MODERATE = "moderate"  # Semi-quantitative, standard approaches
    COMPLEX = "complex"    # Quantitative, multiple approaches


class GeomorphicBoundary(Enum):
    """Types of major geomorphic boundaries for study area definition"""
    WATERSHED_DIVIDE = "watershed_divide"
    RESERVOIR = "reservoir"
    MAJOR_CONFLUENCE = "major_confluence"
    BEDROCK_CONTROL = "bedrock_control"
    GEOLOGIC_CONTACT = "geologic_contact"


class ProblemArea(Enum):
    """Potential problem areas for sediment issues in restoration projects"""
    EXPANSIONS = "expansions"
    BRIDGE_CROSSINGS = "bridge_crossings"
    SLOPE_CHANGES = "slope_changes"
    CUTOFFS = "cutoffs"
    UPSTREAM_APPROACH = "upstream_approach"
    DOWNSTREAM_TRANSITION = "downstream_transition"
    CHANNEL_STRUCTURES = "channel_structures"
    TRIBUTARY_JUNCTIONS = "tributary_junctions"
    WATER_DIVERSIONS = "water_diversions"
    UPSTREAM_RESERVOIRS = "upstream_reservoirs"
    DOWNSTREAM_DAMS = "downstream_dams"


@dataclass
class StudyBoundary:
    """Definition of study area boundaries"""
    upstream_limit: float      # Distance upstream of project (m)
    downstream_limit: float    # Distance downstream of project (m)
    lateral_extent: float      # Lateral extent from channel (m)
    boundary_type: GeomorphicBoundary
    justification: str


@dataclass
class DataInventory:
    """Catalog of available data for sediment studies"""
    geometric_data: Dict[str, List[str]]      # Survey data, cross sections, etc.
    hydrologic_data: Dict[str, Tuple[int, int]]  # Data type: (start_year, end_year)
    hydraulic_data: Dict[str, str]            # Hydraulic model data and reports
    sediment_data: Dict[str, List[str]]       # Bed material, suspended sediment
    land_use_data: Dict[str, List[str]]       # Historic and projected land use
    data_quality: Dict[str, str]              # Quality assessment for each dataset


@dataclass
class StudyResource:
    """Resource requirements for study components"""
    personnel_months: float
    equipment_costs: float
    laboratory_costs: float
    data_acquisition_costs: float
    total_cost: float
    timeline_months: int


def define_study_area_boundary(project_reach_length: float,
                              channel_width: float,
                              watershed_area: float,
                              major_controls: List[str],
                              resource_limitations: bool = False) -> StudyBoundary:
    """
    Define boundaries of the sediment study area.
    
    Args:
        project_reach_length: Length of project reach (m)
        channel_width: Representative channel width (m)  
        watershed_area: Contributing watershed area (km²)
        major_controls: List of major geomorphic controls
        resource_limitations: Whether resource constraints limit study area
        
    Returns:
        StudyBoundary object with recommended limits
        
    Reference: ASCE Manual 110, Section 9.2.1
    """
    validate_positive(project_reach_length, "project_reach_length")
    validate_positive(channel_width, "channel_width")
    validate_positive(watershed_area, "watershed_area")
    
    # Ideal study area extends to major geomorphic boundaries
    if not resource_limitations:
        # Upstream limit: 10-20 times project length or to major control
        upstream_limit = max(10 * project_reach_length, 5000.0)  # Min 5 km
        
        # Downstream limit: 5-10 times project length
        downstream_limit = max(5 * project_reach_length, 2000.0)  # Min 2 km
        
        # Lateral extent: 5-10 times channel width minimum
        lateral_extent = max(10 * channel_width, 100.0)  # Min 100 m
        
        # Adjust based on major controls
        if "reservoir" in [control.lower() for control in major_controls]:
            boundary_type = GeomorphicBoundary.RESERVOIR
        elif "confluence" in [control.lower() for control in major_controls]:
            boundary_type = GeomorphicBoundary.MAJOR_CONFLUENCE
        elif "bedrock" in [control.lower() for control in major_controls]:
            boundary_type = GeomorphicBoundary.BEDROCK_CONTROL
        else:
            boundary_type = GeomorphicBoundary.WATERSHED_DIVIDE
            
    else:
        # Resource-limited study area
        upstream_limit = max(5 * project_reach_length, 1000.0)
        downstream_limit = max(3 * project_reach_length, 500.0)
        lateral_extent = max(5 * channel_width, 50.0)
        boundary_type = GeomorphicBoundary.MAJOR_CONFLUENCE
    
    justification = f"""
    Study area boundaries based on:
    - Project reach length: {project_reach_length:.0f} m
    - Geomorphic controls: {', '.join(major_controls)}
    - Watershed area: {watershed_area:.1f} km²
    - Resource limitations: {'Yes' if resource_limitations else 'No'}
    """
    
    return StudyBoundary(
        upstream_limit=upstream_limit,
        downstream_limit=downstream_limit,
        lateral_extent=lateral_extent,
        boundary_type=boundary_type,
        justification=justification.strip()
    )


def identify_potential_problem_areas(project_features: List[str],
                                   channel_characteristics: Dict[str, float],
                                   infrastructure_present: List[str]) -> Dict[ProblemArea, float]:
    """
    Identify potential sediment problem areas in the study reach.
    
    Args:
        project_features: List of planned project features
        channel_characteristics: Dict of channel geometry and properties
        infrastructure_present: List of existing infrastructure
        
    Returns:
        Dictionary mapping problem areas to risk levels (0-1)
        
    Reference: ASCE Manual 110, Section 9.2.3
    """
    problem_areas = {}
    
    # Channel expansions and contractions
    if "width_ratio" in channel_characteristics:
        width_ratio = channel_characteristics["width_ratio"]
        if width_ratio > 1.5 or width_ratio < 0.7:
            problem_areas[ProblemArea.EXPANSIONS] = min(abs(width_ratio - 1.0), 1.0)
    
    # Bridge crossings and constrictions
    if any("bridge" in feature.lower() for feature in infrastructure_present):
        problem_areas[ProblemArea.BRIDGE_CROSSINGS] = 0.8
    
    # Abrupt slope changes
    if "slope_change" in channel_characteristics:
        slope_change = channel_characteristics["slope_change"]
        if slope_change > 0.001:  # >0.1% change
            problem_areas[ProblemArea.SLOPE_CHANGES] = min(slope_change / 0.005, 1.0)
    
    # Cutoffs and alignment changes
    if any("cutoff" in feature.lower() or "realign" in feature.lower() 
           for feature in project_features):
        problem_areas[ProblemArea.CUTOFFS] = 0.9
    
    # Project transitions
    problem_areas[ProblemArea.UPSTREAM_APPROACH] = 0.6  # Always moderate risk
    problem_areas[ProblemArea.DOWNSTREAM_TRANSITION] = 0.6
    
    # Channel structures
    if any("weir" in feature.lower() or "grade control" in feature.lower() or 
           "structure" in feature.lower() for feature in project_features):
        problem_areas[ProblemArea.CHANNEL_STRUCTURES] = 0.7
    
    # Tributary junctions
    if "tributary_count" in channel_characteristics:
        trib_count = channel_characteristics["tributary_count"]
        if trib_count > 0:
            problem_areas[ProblemArea.TRIBUTARY_JUNCTIONS] = min(trib_count / 5.0, 1.0)
    
    # Water diversions
    if any("diversion" in feature.lower() for feature in infrastructure_present):
        problem_areas[ProblemArea.WATER_DIVERSIONS] = 0.8
    
    # Reservoirs and dams
    if any("dam" in feature.lower() or "reservoir" in feature.lower() 
           for feature in infrastructure_present):
        problem_areas[ProblemArea.UPSTREAM_RESERVOIRS] = 0.9
        problem_areas[ProblemArea.DOWNSTREAM_DAMS] = 0.9
    
    return problem_areas


def develop_data_inventory(study_area: StudyBoundary,
                         available_data_sources: Dict[str, List[str]],
                         project_timeline: int) -> DataInventory:
    """
    Catalog available geometric, hydrologic, hydraulic, and sedimentary data.
    
    Args:
        study_area: Study area boundaries
        available_data_sources: Dict of available data by category
        project_timeline: Project timeline in years
        
    Returns:
        DataInventory object with cataloged data and quality assessment
        
    Reference: ASCE Manual 110, Section 9.2.4
    """
    # Initialize data inventory categories
    geometric_data = available_data_sources.get("geometric", [])
    hydrologic_sources = available_data_sources.get("hydrologic", [])
    hydraulic_data = available_data_sources.get("hydraulic", [])
    sediment_sources = available_data_sources.get("sediment", [])
    land_use_sources = available_data_sources.get("land_use", [])
    
    # Assess hydrologic data period of record
    hydrologic_data = {}
    current_year = 2024  # Would use actual current year in practice
    
    for source in hydrologic_sources:
        if "USGS" in source:
            # Assume USGS data has good period of record
            start_year = max(1950, current_year - 50)
            hydrologic_data[source] = (start_year, current_year)
        elif "local" in source.lower():
            # Shorter record for local data
            start_year = max(1990, current_year - 20)
            hydrologic_data[source] = (start_year, current_year)
        else:
            # Default shorter record
            start_year = current_year - 10
            hydrologic_data[source] = (start_year, current_year)
    
    # Assess data quality
    data_quality = {}
    
    # Geometric data quality
    if geometric_data:
        if any("lidar" in source.lower() for source in geometric_data):
            data_quality["geometric"] = "high"
        elif any("survey" in source.lower() for source in geometric_data):
            data_quality["geometric"] = "moderate"
        else:
            data_quality["geometric"] = "low"
    else:
        data_quality["geometric"] = "none"
    
    # Hydrologic data quality
    min_record_length = min([end - start for start, end in hydrologic_data.values()] + [0])
    if min_record_length >= 30:
        data_quality["hydrologic"] = "high"
    elif min_record_length >= 10:
        data_quality["hydrologic"] = "moderate" 
    elif min_record_length > 0:
        data_quality["hydrologic"] = "low"
    else:
        data_quality["hydrologic"] = "none"
    
    # Sediment data quality (typically limited)
    if sediment_sources:
        data_quality["sediment"] = "moderate" if len(sediment_sources) > 2 else "low"
    else:
        data_quality["sediment"] = "none"
    
    # Land use data quality
    if any("GIS" in source for source in land_use_sources):
        data_quality["land_use"] = "high"
    elif land_use_sources:
        data_quality["land_use"] = "moderate"
    else:
        data_quality["land_use"] = "low"
    
    return DataInventory(
        geometric_data={"surveys": geometric_data},
        hydrologic_data=hydrologic_data,
        hydraulic_data={"models": hydraulic_data},
        sediment_data={"samples": sediment_sources},
        land_use_data={"sources": land_use_sources},
        data_quality=data_quality
    )


def determine_study_approach(project_complexity: StudyComplexity,
                           problem_areas: Dict[ProblemArea, float],
                           data_inventory: DataInventory,
                           risk_tolerance: float) -> Dict[str, List[str]]:
    """
    Determine appropriate study methodology based on project characteristics.
    
    Args:
        project_complexity: Complexity level of restoration project
        problem_areas: Identified problem areas with risk levels
        data_inventory: Available data inventory
        risk_tolerance: Project risk tolerance (0-1, higher = more tolerant)
        
    Returns:
        Dictionary of recommended study approaches by category
        
    Reference: ASCE Manual 110, Section 9.2.5
    """
    validate_range(risk_tolerance, 0.0, 1.0, "risk_tolerance")
    
    study_approaches = {
        "stability_assessment": [],
        "sediment_budget": [],
        "hydraulic_modeling": [],
        "field_investigation": [],
        "monitoring": []
    }
    
    # Maximum risk level from problem areas
    max_risk = max(problem_areas.values()) if problem_areas else 0.0
    
    # Stability assessment approach
    if project_complexity == StudyComplexity.SIMPLE and max_risk < 0.5:
        study_approaches["stability_assessment"] = ["qualitative_assessment"]
    elif project_complexity == StudyComplexity.MODERATE or max_risk < 0.7:
        study_approaches["stability_assessment"] = [
            "semi_quantitative_assessment",
            "threshold_analysis",
            "historical_trend_analysis"
        ]
    else:  # Complex or high risk
        study_approaches["stability_assessment"] = [
            "quantitative_stability_analysis",
            "multiple_assessment_methods",
            "uncertainty_analysis",
            "sensitivity_analysis"
        ]
    
    # Sediment budget requirements
    if max_risk > 0.6 or project_complexity != StudyComplexity.SIMPLE:
        if data_inventory.data_quality.get("sediment", "none") != "none":
            study_approaches["sediment_budget"] = [
                "reach_scale_sediment_budget",
                "with_without_project_scenarios"
            ]
        else:
            study_approaches["sediment_budget"] = [
                "simplified_sediment_balance",
                "order_of_magnitude_analysis"
            ]
    
    # Hydraulic modeling requirements
    data_available = (data_inventory.data_quality.get("geometric", "none") != "none" and
                     data_inventory.data_quality.get("hydrologic", "none") != "none")
    
    if data_available and (max_risk > 0.5 or project_complexity != StudyComplexity.SIMPLE):
        study_approaches["hydraulic_modeling"] = [
            "1D_hydraulic_modeling",
            "design_discharge_analysis",
            "water_surface_profile_analysis"
        ]
        
        if project_complexity == StudyComplexity.COMPLEX:
            study_approaches["hydraulic_modeling"].extend([
                "2D_hydraulic_modeling",
                "unsteady_flow_analysis"
            ])
    
    # Field investigation needs
    if data_inventory.data_quality.get("geometric", "none") == "none":
        study_approaches["field_investigation"].append("topographic_survey")
    
    if data_inventory.data_quality.get("sediment", "none") == "none":
        study_approaches["field_investigation"].extend([
            "bed_material_sampling",
            "suspended_sediment_monitoring"
        ])
    
    if max_risk > 0.7:
        study_approaches["field_investigation"].extend([
            "geomorphic_assessment",
            "bank_stability_evaluation"
        ])
    
    # Monitoring requirements
    study_approaches["monitoring"] = [
        "baseline_monitoring",
        "construction_monitoring"
    ]
    
    if risk_tolerance < 0.5:  # Low risk tolerance = more monitoring
        study_approaches["monitoring"].extend([
            "intensive_post_project_monitoring",
            "adaptive_management_triggers"
        ])
    else:
        study_approaches["monitoring"].append("standard_post_project_monitoring")
    
    return study_approaches


def estimate_study_resources(study_approaches: Dict[str, List[str]],
                           study_area: StudyBoundary,
                           data_collection_required: List[str]) -> StudyResource:
    """
    Estimate personnel, equipment, and cost requirements for sediment studies.
    
    Args:
        study_approaches: Dictionary of study approaches by category
        study_area: Study area boundaries
        data_collection_required: List of required data collection activities
        
    Returns:
        StudyResource object with detailed resource estimates
        
    Reference: ASCE Manual 110, Section 9.2.7
    """
    # Base resource requirements
    personnel_months = 0.0
    equipment_costs = 0.0
    laboratory_costs = 0.0
    data_acquisition_costs = 0.0
    
    # Study area size factor
    total_reach_length = study_area.upstream_limit + study_area.downstream_limit
    size_factor = max(1.0, total_reach_length / 5000.0)  # Normalize to 5 km baseline
    
    # Personnel requirements by study component
    component_months = {
        "qualitative_assessment": 0.5,
        "semi_quantitative_assessment": 1.5,
        "quantitative_stability_analysis": 3.0,
        "multiple_assessment_methods": 4.0,
        "reach_scale_sediment_budget": 2.0,
        "1D_hydraulic_modeling": 1.5,
        "2D_hydraulic_modeling": 3.0,
        "topographic_survey": 1.0,
        "bed_material_sampling": 0.5,
        "suspended_sediment_monitoring": 1.0,
        "geomorphic_assessment": 1.5,
        "baseline_monitoring": 0.5
    }
    
    # Calculate personnel requirements
    for category, approaches in study_approaches.items():
        for approach in approaches:
            if approach in component_months:
                personnel_months += component_months[approach] * size_factor
    
    # Equipment costs
    for data_type in data_collection_required:
        if "topographic" in data_type.lower():
            equipment_costs += 15000.0  # Survey equipment and crew
        elif "sediment" in data_type.lower():
            equipment_costs += 8000.0   # Sampling equipment
        elif "monitoring" in data_type.lower():
            equipment_costs += 25000.0  # Monitoring instruments
    
    # Laboratory costs
    if "bed_material_sampling" in data_collection_required:
        laboratory_costs += 2000.0 * size_factor  # Sieve analysis, etc.
    if "suspended_sediment_monitoring" in data_collection_required:
        laboratory_costs += 5000.0 * size_factor  # Concentration analysis
    
    # Data acquisition costs
    if any("hydraulic_modeling" in approach for approaches in study_approaches.values()
           for approach in approaches):
        data_acquisition_costs += 10000.0  # Software licenses, data sources
    
    # Calculate total cost (assume $12,000/month loaded personnel cost)
    personnel_cost = personnel_months * 12000.0
    total_cost = personnel_cost + equipment_costs + laboratory_costs + data_acquisition_costs
    
    # Timeline estimate (assume 75% efficiency due to scheduling, weather, etc.)
    timeline_months = max(int(personnel_months / 0.75), 6)
    
    return StudyResource(
        personnel_months=personnel_months,
        equipment_costs=equipment_costs,
        laboratory_costs=laboratory_costs,
        data_acquisition_costs=data_acquisition_costs,
        total_cost=total_cost,
        timeline_months=timeline_months
    )


def develop_data_collection_plan(data_inventory: DataInventory,
                               study_approaches: Dict[str, List[str]],
                               project_timeline: int) -> Dict[str, Dict[str, Union[str, int, List[str]]]]:
    """
    Develop detailed data collection plan for required data.
    
    Args:
        data_inventory: Current data inventory
        study_approaches: Selected study approaches
        project_timeline: Project timeline in months
        
    Returns:
        Dictionary of data collection plans by data type
        
    Reference: ASCE Manual 110, Section 9.2.6
    """
    collection_plan = {}
    
    # Geometric data collection
    if data_inventory.data_quality.get("geometric", "none") in ["none", "low"]:
        collection_plan["geometric"] = {
            "method": "Lidar and field survey",
            "standard": "USGS topographic mapping standards",
            "timeline_months": 3,
            "deliverables": [
                "Digital terrain model",
                "Channel cross sections",
                "Longitudinal profile",
                "Planform mapping"
            ],
            "quality_control": [
                "Survey grade GPS control",
                "Multiple survey crews for verification",
                "Independent check surveys"
            ]
        }
    
    # Hydrologic data
    if data_inventory.data_quality.get("hydrologic", "none") in ["none", "low"]:
        collection_plan["hydrologic"] = {
            "method": "Stream gauging and flow measurement",
            "standard": "USGS Water Supply Paper 2175",
            "timeline_months": min(12, project_timeline - 6),
            "deliverables": [
                "Flow duration curves",
                "Flood frequency analysis", 
                "Stage-discharge relationships",
                "Peak flow records"
            ],
            "quality_control": [
                "Calibrated equipment",
                "Regular gauge maintenance",
                "Independent flow measurements"
            ]
        }
    
    # Sediment data collection
    if ("sediment_budget" in [item for sublist in study_approaches.values() for item in sublist] or
        data_inventory.data_quality.get("sediment", "none") == "none"):
        collection_plan["sediment"] = {
            "method": "Standardized sediment sampling",
            "standard": "Federal Interagency Sedimentation Project methods",
            "timeline_months": 6,
            "deliverables": [
                "Bed material gradations",
                "Suspended sediment concentrations",
                "Sediment rating curves",
                "Bedload transport measurements"
            ],
            "quality_control": [
                "Replicate samples",
                "Split sample analysis",
                "Certified laboratory procedures",
                "Chain of custody protocols"
            ]
        }
    
    # Water quality monitoring (if required for habitat objectives)
    if any("monitoring" in approach for approaches in study_approaches.values()
           for approach in approaches):
        collection_plan["water_quality"] = {
            "method": "Continuous and discrete monitoring",
            "standard": "EPA monitoring protocols",
            "timeline_months": 12,
            "deliverables": [
                "Temperature records",
                "Dissolved oxygen profiles",
                "Turbidity measurements",
                "Nutrient concentrations"
            ],
            "quality_control": [
                "Calibrated instruments",
                "Quality assurance samples",
                "Data validation procedures"
            ]
        }
    
    return collection_plan


def prepare_study_report_outline() -> Dict[str, List[str]]:
    """
    Prepare outline for sediment studies report.
    
    Returns:
        Dictionary with report sections and topics
        
    Reference: ASCE Manual 110, Table 9-4
    """
    report_outline = {
        "executive_summary": [
            "Project overview",
            "Key findings",
            "Recommendations",
            "Risk assessment"
        ],
        "geography": [
            "Project and study area boundaries",
            "Current watershed land use",
            "Projected future watershed land use",
            "Regional geomorphic setting"
        ],
        "data": [
            "Available data sources",
            "Data quality assessment",
            "Recommendations for additional data collection",
            "Data gaps and limitations"
        ],
        "history": [
            "Historic land use in contributing watershed",
            "Hydrologic record analysis",
            "Stream behavior in study reach",
            "Aggrading and degrading trends",
            "Flood event impacts",
            "Historical changes to river system"
        ],
        "bed_and_banks": [
            "Bed controls and material characteristics",
            "Bank heights, angles, and vegetation",
            "Bank stability assessment",
            "Substrate conditions"
        ],
        "channel_stability": [
            "Existing channel conditions",
            "Problems upstream and downstream",
            "Knickpoints and knickzones",
            "Stability trends and projections"
        ],
        "physical_habitat": [
            "Current habitat conditions",
            "Features to preserve or modify",
            "Habitat limiting factors",
            "Restoration opportunities"
        ],
        "project_effects": [
            "Water surface elevation changes",
            "Sediment transport capacity impacts",
            "Upstream effects",
            "Downstream effects", 
            "Tributary impacts"
        ],
        "recommendations": [
            "Project alternatives analysis",
            "Preferred restoration approach",
            "Future data collection needs",
            "Monitoring requirements",
            "Adaptive management strategies"
        ]
    }
    
    return report_outline