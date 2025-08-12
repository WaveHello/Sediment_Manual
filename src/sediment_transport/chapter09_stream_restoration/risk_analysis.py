"""
Chapter 9: Stream Restoration - Risk Analysis Module

This module implements risk analysis methods for stream restoration projects
including failure mode analysis, uncertainty quantification, and risk mitigation
strategies based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.1.3.5
- Johnson and Rinaldi (1998), Brookes and Shields (1996)
- FISRWG (1998)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class RiskCategory(Enum):
    """Categories of restoration project risks"""
    TECHNICAL = "technical"           # Design and engineering risks
    ENVIRONMENTAL = "environmental"   # Ecosystem and climate risks
    ECONOMIC = "economic"            # Cost and funding risks
    SOCIAL = "social"                # Stakeholder and regulatory risks
    OPERATIONAL = "operational"      # Construction and maintenance risks


class RiskLevel(Enum):
    """Risk severity levels"""
    LOW = "low"                      # Minor impact, low probability
    MODERATE = "moderate"            # Moderate impact or probability
    HIGH = "high"                    # High impact or probability
    CRITICAL = "critical"            # Severe impact and/or high probability


class FailureMode(Enum):
    """Types of restoration project failures"""
    CHANNEL_INSTABILITY = "channel_instability"     # Erosion/sedimentation
    STRUCTURE_FAILURE = "structure_failure"         # Habitat structure failure
    BIOLOGICAL_FAILURE = "biological_failure"       # No ecological response
    MAINTENANCE_FAILURE = "maintenance_failure"     # Unsustainable maintenance
    STAKEHOLDER_REJECTION = "stakeholder_rejection" # Social/political failure


@dataclass
class RiskFactor:
    """Individual risk factor definition"""
    risk_name: str
    category: RiskCategory
    probability: float               # 0-1 scale
    impact_severity: float          # 0-1 scale
    risk_level: RiskLevel
    description: str
    potential_causes: List[str]
    mitigation_strategies: List[str]
    monitoring_indicators: List[str]


@dataclass
class FailureModeAnalysis:
    """Failure mode and effects analysis"""
    failure_mode: FailureMode
    probability_of_occurrence: float  # 0-1 scale
    severity_of_consequences: float   # 0-1 scale
    detectability: float             # 0-1 scale (1 = easily detected)
    risk_priority_number: float      # Combined risk score
    failure_mechanisms: List[str]
    warning_signs: List[str]
    preventive_measures: List[str]


def assess_technical_risks(design_parameters: Dict[str, float],
                         site_conditions: Dict[str, str],
                         data_quality: Dict[str, str]) -> List[RiskFactor]:
    """
    Assess technical and design-related risks.
    
    Args:
        design_parameters: Design discharge, geometry, etc.
        site_conditions: Geological, hydrological conditions
        data_quality: Quality of design data
        
    Returns:
        List of technical risk factors
        
    Reference: ASCE Manual 110, Section 9.1.3.5
    """
    risks = []
    
    # Channel design risks
    design_discharge = design_parameters.get("design_discharge", 0)
    channel_slope = design_parameters.get("channel_slope", 0)
    
    if data_quality.get("hydrologic", "poor") in ["poor", "limited"]:
        risks.append(RiskFactor(
            risk_name="inadequate_hydrologic_data",
            category=RiskCategory.TECHNICAL,
            probability=0.7,
            impact_severity=0.8,
            risk_level=RiskLevel.HIGH,
            description="Insufficient hydrologic data for reliable design",
            potential_causes=[
                "Short period of record",
                "Ungaged watershed",
                "Climate change uncertainty"
            ],
            mitigation_strategies=[
                "Use multiple estimation methods",
                "Include safety factors in design",
                "Plan for adaptive management"
            ],
            monitoring_indicators=[
                "Actual vs predicted flows",
                "Channel stability during events"
            ]
        ))
    
    # Sediment transport risks
    if data_quality.get("sediment", "poor") == "poor":
        risks.append(RiskFactor(
            risk_name="sediment_transport_uncertainty",
            category=RiskCategory.TECHNICAL,
            probability=0.8,
            impact_severity=0.7,
            risk_level=RiskLevel.HIGH,
            description="High uncertainty in sediment transport predictions",
            potential_causes=[
                "No local sediment data",
                "Complex transport relationships",
                "Upstream land use changes"
            ],
            mitigation_strategies=[
                "Conservative channel design",
                "Include sediment management provisions",
                "Monitor sediment transport closely"
            ],
            monitoring_indicators=[
                "Bed elevation changes",
                "Sediment deposition patterns",
                "Channel capacity changes"
            ]
        ))
    
    # Geotechnical risks
    if site_conditions.get("foundation", "unknown") == "weak":
        risks.append(RiskFactor(
            risk_name="geotechnical_instability",
            category=RiskCategory.TECHNICAL,
            probability=0.6,
            impact_severity=0.9,
            risk_level=RiskLevel.HIGH,
            description="Poor foundation conditions for structures",
            potential_causes=[
                "Weak soil layers",
                "High groundwater",
                "Expansive clays"
            ],
            mitigation_strategies=[
                "Deep foundations",
                "Soil improvement",
                "Alternative structure types"
            ],
            monitoring_indicators=[
                "Structure settlement",
                "Cracking or deformation",
                "Groundwater levels"
            ]
        ))
    
    # Channel stability risks
    if channel_slope > 0.01:  # Steep channels
        risks.append(RiskFactor(
            risk_name="channel_incision_risk",
            category=RiskCategory.TECHNICAL,
            probability=0.5,
            impact_severity=0.8,
            risk_level=RiskLevel.MODERATE,
            description="Risk of channel downcutting in steep reaches",
            potential_causes=[
                "Excess stream power",
                "Inadequate grade control",
                "Upstream flow changes"
            ],
            mitigation_strategies=[
                "Install grade control structures",
                "Use energy dissipation features",
                "Armor critical sections"
            ],
            monitoring_indicators=[
                "Bed elevation changes",
                "Headcut migration",
                "Scour depth measurements"
            ]
        ))
    
    return risks


def assess_environmental_risks(watershed_characteristics: Dict[str, float],
                             climate_projections: Dict[str, float],
                             ecological_context: Dict[str, str]) -> List[RiskFactor]:
    """
    Assess environmental and ecological risks.
    
    Args:
        watershed_characteristics: Size, land use, etc.
        climate_projections: Projected changes in climate
        ecological_context: Species, habitat conditions
        
    Returns:
        List of environmental risk factors
        
    Reference: ASCE Manual 110, Section 9.8
    """
    risks = []
    
    # Climate change risks
    temperature_change = climate_projections.get("temperature_change", 0)
    precipitation_change = climate_projections.get("precipitation_change", 0)
    
    if abs(temperature_change) > 2.0:  # >2°C change
        risks.append(RiskFactor(
            risk_name="climate_change_impacts",
            category=RiskCategory.ENVIRONMENTAL,
            probability=0.8,
            impact_severity=0.7,
            risk_level=RiskLevel.HIGH,
            description="Significant climate change impacts on project performance",
            potential_causes=[
                "Temperature increases",
                "Altered precipitation patterns",
                "Extreme event frequency changes"
            ],
            mitigation_strategies=[
                "Design for climate uncertainty",
                "Include adaptation measures",
                "Monitor climate indicators"
            ],
            monitoring_indicators=[
                "Water temperature trends",
                "Flow regime changes",
                "Species composition shifts"
            ]
        ))
    
    # Extreme event risks
    if precipitation_change > 20:  # >20% increase
        risks.append(RiskFactor(
            risk_name="extreme_flood_risk",
            category=RiskCategory.ENVIRONMENTAL,
            probability=0.6,
            impact_severity=0.9,
            risk_level=RiskLevel.HIGH,
            description="Increased risk of extreme flood events",
            potential_causes=[
                "Climate change",
                "Urbanization",
                "Land use changes"
            ],
            mitigation_strategies=[
                "Design for larger floods",
                "Include emergency spillways",
                "Plan for post-flood recovery"
            ],
            monitoring_indicators=[
                "Peak flow magnitudes",
                "Flood damage assessment",
                "Structure performance"
            ]
        ))
    
    # Ecological risks
    if ecological_context.get("invasive_species", "none") != "none":
        risks.append(RiskFactor(
            risk_name="invasive_species_colonization",
            category=RiskCategory.ENVIRONMENTAL,
            probability=0.7,
            impact_severity=0.6,
            risk_level=RiskLevel.MODERATE,
            description="Risk of invasive species establishment",
            potential_causes=[
                "Disturbed substrate",
                "Altered flow conditions",
                "Nearby invasive populations"
            ],
            mitigation_strategies=[
                "Use native species only",
                "Rapid revegetation",
                "Invasive species monitoring"
            ],
            monitoring_indicators=[
                "Species composition surveys",
                "Vegetation cover assessment",
                "Growth rate monitoring"
            ]
        ))
    
    # Watershed-scale risks
    urban_fraction = watershed_characteristics.get("urban_fraction", 0)
    if urban_fraction > 0.3:
        risks.append(RiskFactor(
            risk_name="continued_watershed_degradation",
            category=RiskCategory.ENVIRONMENTAL,
            probability=0.8,
            impact_severity=0.8,
            risk_level=RiskLevel.HIGH,
            description="Ongoing watershed development impacts",
            potential_causes=[
                "Continued urbanization",
                "Stormwater runoff increases",
                "Water quality degradation"
            ],
            mitigation_strategies=[
                "Watershed-scale planning",
                "Stormwater management",
                "Land use controls"
            ],
            monitoring_indicators=[
                "Runoff patterns",
                "Water quality trends",
                "Land use change rates"
            ]
        ))
    
    return risks


def perform_failure_mode_analysis(project_components: List[str],
                                design_complexity: str,
                                maintenance_capability: str) -> List[FailureModeAnalysis]:
    """
    Perform failure mode and effects analysis (FMEA).
    
    Args:
        project_components: List of major project components
        design_complexity: "simple", "moderate", or "complex"
        maintenance_capability: "high", "moderate", or "low"
        
    Returns:
        List of failure mode analyses
        
    Reference: ASCE Manual 110, Section 9.1.3.5
    """
    analyses = []
    
    # Channel instability failure mode
    channel_prob = 0.3 if design_complexity == "simple" else 0.5 if design_complexity == "moderate" else 0.7
    channel_severity = 0.8  # High consequences
    channel_detect = 0.7    # Moderately detectable
    
    analyses.append(FailureModeAnalysis(
        failure_mode=FailureMode.CHANNEL_INSTABILITY,
        probability_of_occurrence=channel_prob,
        severity_of_consequences=channel_severity,
        detectability=1.0 - channel_detect,  # Lower is better for risk
        risk_priority_number=channel_prob * channel_severity * (1.0 - channel_detect) * 1000,
        failure_mechanisms=[
            "Sediment transport imbalance",
            "Extreme flood events",
            "Upstream flow changes",
            "Bank erosion and widening"
        ],
        warning_signs=[
            "Bed elevation changes",
            "Bank scour development",
            "Increased turbidity",
            "Structure undermining"
        ],
        preventive_measures=[
            "Regular sediment monitoring",
            "Adaptive management protocols",
            "Upstream flow management",
            "Bank protection maintenance"
        ]
    ))
    
    # Structure failure mode (if structures present)
    if any("structure" in comp.lower() for comp in project_components):
        struct_prob = 0.2 if maintenance_capability == "high" else 0.4 if maintenance_capability == "moderate" else 0.6
        struct_severity = 0.6  # Moderate to high consequences
        struct_detect = 0.8    # Usually detectable
        
        analyses.append(FailureModeAnalysis(
            failure_mode=FailureMode.STRUCTURE_FAILURE,
            probability_of_occurrence=struct_prob,
            severity_of_consequences=struct_severity,
            detectability=1.0 - struct_detect,
            risk_priority_number=struct_prob * struct_severity * (1.0 - struct_detect) * 1000,
            failure_mechanisms=[
                "Material degradation",
                "Flood damage",
                "Poor installation",
                "Inadequate maintenance"
            ],
            warning_signs=[
                "Structural cracking",
                "Settlement or tilting",
                "Scour around foundations",
                "Material loss"
            ],
            preventive_measures=[
                "Regular structural inspections",
                "Preventive maintenance",
                "Flood damage assessment",
                "Material selection review"
            ]
        ))
    
    # Biological failure mode
    bio_prob = 0.4 if design_complexity == "complex" else 0.6  # Higher for simple projects
    bio_severity = 0.7  # High ecological consequences
    bio_detect = 0.5    # Often hard to detect early
    
    analyses.append(FailureModeAnalysis(
        failure_mode=FailureMode.BIOLOGICAL_FAILURE,
        probability_of_occurrence=bio_prob,
        severity_of_consequences=bio_severity,
        detectability=1.0 - bio_detect,
        risk_priority_number=bio_prob * bio_severity * (1.0 - bio_detect) * 1000,
        failure_mechanisms=[
            "Inadequate habitat creation",
            "Water quality limitations",
            "Species colonization barriers",
            "Competition from invasives"
        ],
        warning_signs=[
            "Low species diversity",
            "Poor survival rates",
            "Invasive species dominance",
            "Continued habitat degradation"
        ],
        preventive_measures=[
            "Comprehensive habitat design",
            "Water quality management",
            "Species introduction programs",
            "Invasive species control"
        ]
    ))
    
    # Maintenance failure mode
    maint_prob = 0.1 if maintenance_capability == "high" else 0.3 if maintenance_capability == "moderate" else 0.7
    maint_severity = 0.8  # Can lead to other failures
    maint_detect = 0.6    # Moderately detectable
    
    analyses.append(FailureModeAnalysis(
        failure_mode=FailureMode.MAINTENANCE_FAILURE,
        probability_of_occurrence=maint_prob,
        severity_of_consequences=maint_severity,
        detectability=1.0 - maint_detect,
        risk_priority_number=maint_prob * maint_severity * (1.0 - maint_detect) * 1000,
        failure_mechanisms=[
            "Insufficient funding",
            "Lack of technical expertise",
            "Access difficulties",
            "Institutional changes"
        ],
        warning_signs=[
            "Deferred maintenance",
            "Budget reductions",
            "Staff turnover",
            "Performance degradation"
        ],
        preventive_measures=[
            "Establish maintenance endowment",
            "Train local personnel",
            "Simplify maintenance requirements",
            "Create maintenance partnerships"
        ]
    ))
    
    return analyses


def develop_risk_mitigation_plan(risk_factors: List[RiskFactor],
                               failure_modes: List[FailureModeAnalysis],
                               project_budget: float,
                               risk_tolerance: str) -> Dict[str, List[str]]:
    """
    Develop comprehensive risk mitigation plan.
    
    Args:
        risk_factors: Identified risk factors
        failure_modes: Failure mode analyses
        project_budget: Available budget for mitigation ($)
        risk_tolerance: "low", "moderate", or "high"
        
    Returns:
        Dictionary with mitigation strategies by category
        
    Reference: ASCE Manual 110, Section 9.1.3.5
    """
    validate_positive(project_budget, "project_budget")
    
    mitigation_plan = {
        "design_modifications": [],
        "monitoring_enhancements": [],
        "contingency_planning": [],
        "stakeholder_engagement": [],
        "adaptive_management": []
    }
    
    # Priority ranking based on risk level and RPN
    high_priority_risks = [r for r in risk_factors if r.risk_level in [RiskLevel.HIGH, RiskLevel.CRITICAL]]
    high_priority_failures = sorted(failure_modes, key=lambda x: x.risk_priority_number, reverse=True)[:2]
    
    # Design modifications for high-priority risks
    for risk in high_priority_risks:
        if risk.category == RiskCategory.TECHNICAL:
            mitigation_plan["design_modifications"].extend([
                "Include safety factors in critical design parameters",
                "Add redundancy to critical systems",
                "Use proven design approaches where possible"
            ])
        elif risk.category == RiskCategory.ENVIRONMENTAL:
            mitigation_plan["design_modifications"].extend([
                "Design for climate change scenarios",
                "Include flexible/adaptable features",
                "Plan for extreme event resilience"
            ])
    
    # Enhanced monitoring for critical failure modes
    for failure in high_priority_failures:
        mitigation_plan["monitoring_enhancements"].extend(failure.warning_signs)
        mitigation_plan["monitoring_enhancements"].extend([
            f"Implement early warning system for {failure.failure_mode.value}",
            f"Increase monitoring frequency for {failure.failure_mode.value} indicators"
        ])
    
    # Contingency planning
    if risk_tolerance == "low":
        mitigation_plan["contingency_planning"].extend([
            "Develop detailed emergency response procedures",
            "Establish material and equipment stockpiles",
            "Create rapid response team protocols",
            "Plan for temporary project modifications"
        ])
    
    # Budget-based mitigation prioritization
    if project_budget > 1000000:  # Large project
        mitigation_plan["design_modifications"].append("Include backup systems for critical components")
        mitigation_plan["monitoring_enhancements"].append("Implement real-time monitoring systems")
    elif project_budget > 500000:  # Medium project
        mitigation_plan["monitoring_enhancements"].append("Include automated data collection where feasible")
    else:  # Small project
        mitigation_plan["adaptive_management"].append("Focus on simple, low-cost mitigation measures")
    
    # Stakeholder engagement for social risks
    social_risks = [r for r in risk_factors if r.category == RiskCategory.SOCIAL]
    if social_risks:
        mitigation_plan["stakeholder_engagement"].extend([
            "Maintain regular communication with stakeholders",
            "Provide project performance updates",
            "Address concerns promptly and transparently",
            "Include stakeholders in monitoring activities"
        ])
    
    # Adaptive management strategies
    mitigation_plan["adaptive_management"].extend([
        "Establish clear trigger points for management actions",
        "Develop decision trees for common scenarios",
        "Create flexible project modification procedures",
        "Plan for iterative design improvements"
    ])
    
    # Remove duplicates and prioritize
    for category in mitigation_plan:
        mitigation_plan[category] = list(set(mitigation_plan[category]))
    
    return mitigation_plan


def calculate_project_risk_score(risk_factors: List[RiskFactor],
                               failure_modes: List[FailureModeAnalysis],
                               mitigation_effectiveness: float = 0.5) -> Dict[str, Union[float, str]]:
    """
    Calculate overall project risk score and classification.
    
    Args:
        risk_factors: List of identified risk factors
        failure_modes: List of failure mode analyses
        mitigation_effectiveness: Effectiveness of mitigation (0-1)
        
    Returns:
        Dictionary with risk scores and classifications
        
    Reference: Based on risk assessment methodologies
    """
    validate_range(mitigation_effectiveness, 0.0, 1.0, "mitigation_effectiveness")
    
    # Calculate weighted risk score from risk factors
    risk_weights = {
        RiskLevel.LOW: 1,
        RiskLevel.MODERATE: 2, 
        RiskLevel.HIGH: 4,
        RiskLevel.CRITICAL: 8
    }
    
    total_risk_score = 0.0
    risk_count = len(risk_factors)
    
    for risk in risk_factors:
        individual_score = risk.probability * risk.impact_severity * risk_weights[risk.risk_level]
        total_risk_score += individual_score
    
    average_risk_score = total_risk_score / max(risk_count, 1)
    
    # Calculate failure mode risk contribution
    total_rpn = sum([fm.risk_priority_number for fm in failure_modes])
    average_rpn = total_rpn / max(len(failure_modes), 1)
    
    # Combined risk score (normalized to 0-100 scale)
    combined_score = (average_risk_score * 10) + (average_rpn / 10)
    combined_score = min(combined_score, 100)
    
    # Apply mitigation effectiveness
    mitigated_score = combined_score * (1.0 - mitigation_effectiveness)
    
    # Risk classification
    if mitigated_score < 20:
        risk_class = "LOW"
        recommendation = "Proceed with standard monitoring"
    elif mitigated_score < 40:
        risk_class = "MODERATE"  
        recommendation = "Implement enhanced monitoring and mitigation measures"
    elif mitigated_score < 70:
        risk_class = "HIGH"
        recommendation = "Require comprehensive risk mitigation plan"
    else:
        risk_class = "CRITICAL"
        recommendation = "Consider project redesign or cancellation"
    
    return {
        "raw_risk_score": combined_score,
        "mitigated_risk_score": mitigated_score,
        "risk_classification": risk_class,
        "recommendation": recommendation,
        "dominant_risk_categories": [r.category.value for r in risk_factors if r.risk_level == RiskLevel.CRITICAL],
        "critical_failure_modes": [fm.failure_mode.value for fm in failure_modes if fm.risk_priority_number > 200]
    }