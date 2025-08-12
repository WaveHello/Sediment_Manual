"""
Chapter 9: Stream Restoration - Restoration Design Module

This module implements restoration design methods including channel design,
habitat structures, floodplain connectivity, and construction considerations
for stream restoration projects based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.5
- FISRWG (1998), Thompson (2002b)
- Shields (1983), ASCE Task Committee (1992)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from enum import Enum
from dataclasses import dataclass

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class ChannelDesignMethod(Enum):
    """Channel design methods for restoration"""
    STABLE_CHANNEL = "stable_channel"        # Traditional stable design
    THRESHOLD_APPROACH = "threshold_approach"  # Incipient motion based
    ACTIVE_BED_APPROACH = "active_bed_approach"  # Mobile bed design
    NATURAL_ANALOGUE = "natural_analogue"    # Reference reach based


class HabitatStructureType(Enum):
    """Types of in-channel habitat structures"""
    SILLS = "sills"                   # Cross-channel structures
    DEFLECTORS = "deflectors"         # Bank-attached flow deflection
    RANDOM_ROCKS = "random_rocks"     # Isolated flow obstructions
    COVERS = "covers"                 # Shade and hiding structures
    LARGE_WOODY_DEBRIS = "large_woody_debris"  # Natural wood structures


class RestorationApproach(Enum):
    """Overall restoration approach philosophy"""
    ASSISTED_RECOVERY = "assisted_recovery"    # Remove constraints, let nature work
    HYBRID_APPROACH = "hybrid_approach"        # Combination of natural and engineered
    ENGINEERED_DESIGN = "engineered_design"    # Full reconstruction


@dataclass
class ChannelGeometry:
    """Channel geometric parameters for design"""
    width: float                    # m
    depth: float                   # m
    slope: float                   # dimensionless
    sinuosity: float              # dimensionless
    bed_d50: float                # mm
    bank_angle: float             # degrees
    floodplain_width: float       # m


@dataclass
class HabitatStructure:
    """Definition of habitat structure design"""
    structure_type: HabitatStructureType
    location: str                 # Description of placement
    dimensions: Dict[str, float]  # Key dimensions
    materials: List[str]          # Construction materials
    design_discharge: float       # m³/s
    expected_effects: List[str]   # Intended hydraulic/habitat effects
    maintenance_requirements: List[str]


def calculate_stable_channel_geometry(design_discharge: float,
                                    bed_d50: float,
                                    valley_slope: float,
                                    design_method: str = "threshold") -> ChannelGeometry:
    """
    Calculate stable channel geometry using threshold or regime approaches.
    
    Args:
        design_discharge: Design discharge (m³/s)
        bed_d50: Median bed material size (mm)
        valley_slope: Valley slope (dimensionless)
        design_method: "threshold" or "regime" approach
        
    Returns:
        ChannelGeometry object with calculated dimensions
        
    Reference: ASCE Manual 110, Section 9.5.1.1
    """
    validate_positive(design_discharge, "design_discharge")
    validate_positive(bed_d50, "bed_d50")
    validate_positive(valley_slope, "valley_slope")
    
    if design_method == "threshold":
        # Threshold approach - design for incipient motion of bed material
        
        # Convert D50 to meters
        d50_m = bed_d50 / 1000.0
        
        # Use Shields criterion for critical shear stress
        shields_critical = 0.06  # For coarse sediments
        gamma_s = 2650 * GRAVITY  # N/m³
        gamma_w = 1000 * GRAVITY  # N/m³
        tau_critical = shields_critical * (gamma_s - gamma_w) * d50_m
        
        # Estimate Manning's n from bed material size
        # Strickler relation: n = 0.047 * D50^(1/6) where D50 in meters
        manning_n = 0.047 * (d50_m**(1/6))
        
        # Initial estimate for width using regime relationships
        # W = a * Q^b (typical values: a=4-8, b=0.5)
        width = 6.0 * (design_discharge**0.5)
        
        # Iteratively solve for depth using Manning equation and shear stress
        depth = 1.0  # Initial guess
        for iteration in range(10):
            # Hydraulic radius for rectangular channel
            hydraulic_radius = (width * depth) / (width + 2 * depth)
            
            # Required slope for critical shear stress
            required_slope = tau_critical / (gamma_w * hydraulic_radius)
            
            # Check if required slope exceeds valley slope
            if required_slope > valley_slope:
                # Need larger channel to reduce slope requirement
                width *= 1.1
                depth *= 1.05
            else:
                # Calculate depth from Manning equation
                area = width * depth
                velocity = design_discharge / area
                calculated_depth = ((design_discharge * manning_n) / 
                                  (width * required_slope**0.5))**(3/5)
                
                # Update depth estimate
                depth = 0.7 * depth + 0.3 * calculated_depth
        
        channel_slope = min(required_slope, valley_slope * 0.8)
        
    elif design_method == "regime":
        # Regime approach using hydraulic geometry relationships
        
        # Typical regime relationships for gravel-bed streams:
        # W = 4.33 * Q^0.5, D = 0.17 * Q^0.4, S = 0.0014 * Q^(-0.44)
        width = 4.33 * (design_discharge**0.5)
        depth = 0.17 * (design_discharge**0.4)
        channel_slope = min(0.0014 * (design_discharge**(-0.44)), valley_slope * 0.8)
        
    else:
        raise ValueError(f"Unknown design method: {design_method}")
    
    # Calculate sinuosity based on valley slope and channel slope
    sinuosity = min(valley_slope / channel_slope, 2.5) if channel_slope > 0 else 1.2
    
    # Estimate bank angle (typically 30-45 degrees for stable channels)
    if bed_d50 > 20:  # Coarse material
        bank_angle = 35.0
    elif bed_d50 > 2:  # Gravel
        bank_angle = 30.0  
    else:  # Sand
        bank_angle = 25.0
    
    # Estimate floodplain width (typically 5-10 times channel width)
    floodplain_width = 7.0 * width
    
    return ChannelGeometry(
        width=width,
        depth=depth,
        slope=channel_slope,
        sinuosity=sinuosity,
        bed_d50=bed_d50,
        bank_angle=bank_angle,
        floodplain_width=floodplain_width
    )


def design_meandering_planform(channel_geometry: ChannelGeometry,
                             reach_length: float,
                             constraints: Optional[Dict[str, float]] = None) -> Dict[str, Union[float, List[Tuple[float, float]]]]:
    """
    Design meandering planform geometry for channel reconstruction.
    
    Args:
        channel_geometry: Basic channel geometry
        reach_length: Length of restoration reach (m)
        constraints: Optional constraints (right_of_way_width, etc.)
        
    Returns:
        Dictionary with meander geometry and centerline coordinates
        
    Reference: ASCE Manual 110, Section 9.5.1.2
    """
    validate_positive(reach_length, "reach_length")
    
    # Calculate meander wavelength using Leopold et al. (1964) relationship
    # λ = 11 * W^1.01 (where W is channel width)
    meander_wavelength = 11.0 * (channel_geometry.width**1.01)
    
    # Alternative using Ackers and Charlton (1970) for varying flows
    # λ = 61.21 * Q^0.467 (where Q is formative discharge)
    # For validation, we'll use the width-based approach
    
    # Calculate meander amplitude (typically wavelength/4 to wavelength/6)
    meander_amplitude = meander_wavelength / 5.0
    
    # Apply constraints if provided
    if constraints is not None:
        max_amplitude = constraints.get("right_of_way_width", float('inf')) / 2.0
        meander_amplitude = min(meander_amplitude, max_amplitude)
        
        if "minimum_radius" in constraints:
            min_radius = constraints["minimum_radius"]
            # Radius of curvature ≈ wavelength²/(8 × amplitude)
            required_amplitude = (meander_wavelength**2) / (8 * min_radius)
            meander_amplitude = max(meander_amplitude, required_amplitude)
    
    # Calculate number of complete meanders in reach
    straight_line_length = reach_length / channel_geometry.sinuosity
    num_meanders = straight_line_length / meander_wavelength
    
    # Generate meander centerline coordinates
    centerline_coords = []
    x = 0.0
    y = 0.0
    
    for i in range(int(num_meanders * 20)):  # 20 points per meander
        phase = (i / 20.0) * 2 * np.pi
        dx = meander_wavelength / 20.0
        dy = meander_amplitude * np.sin(phase) - meander_amplitude * np.sin(phase - dx/meander_wavelength * 2 * np.pi)
        
        x += dx
        y += dy
        
        if x >= straight_line_length:
            break
            
        centerline_coords.append((x, y))
    
    # Add variability to meander dimensions (±20% variation)
    np.random.seed(42)  # For reproducible results
    wavelength_variation = np.random.normal(1.0, 0.1, max(1, int(num_meanders)))
    amplitude_variation = np.random.normal(1.0, 0.1, max(1, int(num_meanders)))
    
    return {
        "meander_wavelength": meander_wavelength,
        "meander_amplitude": meander_amplitude,
        "number_of_meanders": num_meanders,
        "centerline_coordinates": centerline_coords,
        "wavelength_variation": wavelength_variation.tolist(),
        "amplitude_variation": amplitude_variation.tolist(),
        "design_sinuosity": channel_geometry.sinuosity,
        "actual_sinuosity": len(centerline_coords) * (meander_wavelength/20.0) / straight_line_length
    }


def design_habitat_structure(structure_type: HabitatStructureType,
                           channel_geometry: ChannelGeometry,
                           design_discharge: float,
                           target_effects: List[str]) -> HabitatStructure:
    """
    Design habitat structure based on type and channel conditions.
    
    Args:
        structure_type: Type of habitat structure
        channel_geometry: Channel geometry at structure location
        design_discharge: Design discharge (m³/s)
        target_effects: List of desired effects
        
    Returns:
        HabitatStructure object with design specifications
        
    Reference: ASCE Manual 110, Table 9-11
    """
    validate_positive(design_discharge, "design_discharge")
    
    if structure_type == HabitatStructureType.SILLS:
        # Cross-channel sills/weirs
        dimensions = {
            "length": channel_geometry.width,
            "height": min(0.3 * channel_geometry.depth, 1.0),  # Max 1m high
            "crest_width": 0.5,  # m
            "notch_width": channel_geometry.width * 0.3  # For fish passage
        }
        
        materials = ["stone", "gabion", "log_weirs"]
        location = "Extending across channel from bank to bank"
        expected_effects = ["increased_scour_downstream", "pool_formation", "grade_control"]
        maintenance_requirements = ["check_for_flanking", "repair_erosion_damage", "debris_removal"]
        
    elif structure_type == HabitatStructureType.DEFLECTORS:
        # Bank-attached flow deflectors
        dimensions = {
            "length": 2.0 * channel_geometry.width,
            "height": 0.5 * channel_geometry.depth,
            "projection": channel_geometry.width * 0.25,  # Into channel
            "angle": 30.0  # Degrees upstream from perpendicular
        }
        
        materials = ["stone", "root_wads", "boulder_clusters"]
        location = "Along banks, extending into channel"
        expected_effects = ["flow_deflection", "bank_scour_reduction", "opposite_bank_deposition"]
        maintenance_requirements = ["check_structural_integrity", "monitor_opposite_bank"]
        
    elif structure_type == HabitatStructureType.RANDOM_ROCKS:
        # Isolated mid-channel obstructions
        dimensions = {
            "boulder_diameter": max(0.5, channel_geometry.bed_d50 * 0.05),  # m
            "cluster_spacing": 3.0 * channel_geometry.width,
            "cluster_size": 3,  # Number of boulders per cluster
            "protrusion": 0.3 * channel_geometry.depth
        }
        
        materials = ["natural_boulders", "concrete_blocks", "root_wads"]
        location = "Isolated mid-channel flow obstructions"
        expected_effects = ["localized_scour", "velocity_refuge_zones", "habitat_diversity"]
        maintenance_requirements = ["reset_displaced_rocks", "clear_debris_accumulation"]
        
    elif structure_type == HabitatStructureType.COVERS:
        # Overhead cover structures
        dimensions = {
            "length": 2.0 * channel_geometry.width,
            "width": channel_geometry.width * 0.5,
            "clearance": channel_geometry.depth + 0.5,
            "support_spacing": 3.0  # m
        }
        
        materials = ["lumber", "natural_logs", "brush_bundles"]
        location = "Along undercut banks or over pools"
        expected_effects = ["shade_provision", "predator_protection", "temperature_control"]
        maintenance_requirements = ["structural_inspection", "vegetation_management"]
        
    elif structure_type == HabitatStructureType.LARGE_WOODY_DEBRIS:
        # Natural wood structures
        dimensions = {
            "log_length": min(15.0, 2.0 * channel_geometry.width),
            "log_diameter": 0.3,  # m minimum
            "rootwad_diameter": 2.0,  # m
            "burial_depth": 1.0,  # m for stability
        }
        
        materials = ["native_hardwood_trees", "rootwad_complexes"]
        location = "Anchored to banks or spanning channel"
        expected_effects = ["pool_formation", "debris_trapping", "bank_protection"]
        maintenance_requirements = ["anchor_inspection", "flood_damage_assessment"]
        
    else:
        raise ValueError(f"Unknown structure type: {structure_type}")
    
    return HabitatStructure(
        structure_type=structure_type,
        location=location,
        dimensions=dimensions,
        materials=materials,
        design_discharge=design_discharge,
        expected_effects=expected_effects,
        maintenance_requirements=maintenance_requirements
    )


def assess_floodplain_connectivity(channel_geometry: ChannelGeometry,
                                 bankfull_discharge: float,
                                 flood_discharges: Dict[int, float]) -> Dict[str, Union[float, List[str]]]:
    """
    Assess channel-floodplain connectivity for restoration design.
    
    Args:
        channel_geometry: Current channel geometry
        bankfull_discharge: Bank-full discharge (m³/s)
        flood_discharges: Dictionary of return period: discharge (m³/s)
        
    Returns:
        Dictionary with connectivity assessment and recommendations
        
    Reference: ASCE Manual 110, Section 9.5.3
    """
    validate_positive(bankfull_discharge, "bankfull_discharge")
    
    # Calculate channel capacity
    area = channel_geometry.width * channel_geometry.depth
    velocity = 1.0  # m/s, rough estimate for capacity calculation
    channel_capacity = area * velocity
    
    # Assess overbank frequency
    overbank_assessment = {}
    recommendations = []
    
    for return_period, discharge in flood_discharges.items():
        if discharge > channel_capacity:
            overbank_assessment[f"{return_period}_year"] = "overbank_flow"
        else:
            overbank_assessment[f"{return_period}_year"] = "contained_in_channel"
    
    # Calculate connectivity ratio
    floodplain_area = channel_geometry.floodplain_width * 1000  # Assume 1 km reach
    channel_area = channel_geometry.width * 1000
    connectivity_ratio = channel_area / floodplain_area
    
    # Recommendations based on connectivity
    if connectivity_ratio > 0.5:
        recommendations.append("Channel too large relative to floodplain")
        recommendations.append("Consider channel narrowing or floodplain expansion")
    elif connectivity_ratio < 0.1:
        recommendations.append("Poor channel-floodplain connectivity")
        recommendations.append("Consider floodplain lowering or berm creation")
    
    # Check for incised conditions
    if channel_geometry.depth > 3.0:
        recommendations.extend([
            "Potentially incised channel",
            "Assess historical channel dimensions",
            "Consider fill or artificial floodplain creation"
        ])
    
    # Flood frequency recommendations
    two_year_flood = flood_discharges.get(2, bankfull_discharge * 2)
    if two_year_flood > channel_capacity * 2:
        recommendations.append("Frequent overbank flooding expected")
    elif two_year_flood < channel_capacity * 1.2:
        recommendations.append("Infrequent overbank flooding - consider channel modifications")
    
    return {
        "connectivity_ratio": connectivity_ratio,
        "channel_capacity": channel_capacity,
        "overbank_assessment": overbank_assessment,
        "recommendations": recommendations,
        "floodplain_access_frequency": "2-10 years" if two_year_flood > channel_capacity else ">10 years"
    }


def calculate_sediment_budget(upstream_yield: float,
                            channel_capacity: float,
                            project_length: float,
                            bed_material_d50: float) -> Dict[str, float]:
    """
    Calculate sediment budget for restoration design.
    
    Args:
        upstream_yield: Annual bed material yield from upstream (m³/year)
        channel_capacity: Channel sediment transport capacity (m³/year)
        project_length: Length of project reach (m)
        bed_material_d50: Median bed material size (mm)
        
    Returns:
        Dictionary with sediment budget analysis
        
    Reference: ASCE Manual 110, Equation 9-15
    """
    validate_positive(upstream_yield, "upstream_yield")
    validate_positive(channel_capacity, "channel_capacity")
    validate_positive(project_length, "project_length")
    
    # Calculate trap efficiency using Equation 9-15
    # E = 100 * (Ys_in - Ys_out) / Ys_in
    outflow_yield = min(channel_capacity, upstream_yield)
    trap_efficiency = 100 * (upstream_yield - outflow_yield) / upstream_yield if upstream_yield > 0 else 0
    
    # Estimate annual deposition/erosion
    annual_deposition = upstream_yield - outflow_yield  # m³/year
    
    # Convert to average bed elevation change
    channel_width = 10.0  # Assume typical width, would be provided in practice
    elevation_change = annual_deposition / (project_length * channel_width)  # m/year
    
    # Sustainability assessment
    if abs(trap_efficiency) < 5:
        stability_forecast = "stable"
    elif trap_efficiency > 20:
        stability_forecast = "aggradational"
    elif trap_efficiency < -20:
        stability_forecast = "degradational"
    else:
        stability_forecast = "marginally_stable"
    
    # Maintenance implications
    if annual_deposition > 100:  # m³/year
        maintenance_frequency = "annual"
    elif annual_deposition > 50:
        maintenance_frequency = "every_2-3_years"
    elif annual_deposition > 0:
        maintenance_frequency = "every_5-10_years"
    else:
        maintenance_frequency = "monitoring_only"
    
    return {
        "upstream_yield": upstream_yield,
        "channel_capacity": channel_capacity,
        "trap_efficiency_percent": trap_efficiency,
        "annual_deposition": annual_deposition,
        "elevation_change_per_year": elevation_change,
        "stability_forecast": stability_forecast,
        "maintenance_frequency": maintenance_frequency,
        "outflow_yield": outflow_yield
    }


def develop_construction_specifications(channel_geometry: ChannelGeometry,
                                      habitat_structures: List[HabitatStructure],
                                      construction_constraints: Dict[str, Union[str, float]]) -> Dict[str, List[str]]:
    """
    Develop construction specifications and sequencing.
    
    Args:
        channel_geometry: Designed channel geometry
        habitat_structures: List of habitat structures to install
        construction_constraints: Constraints like access, timing windows
        
    Returns:
        Dictionary with construction specifications and sequences
        
    Reference: ASCE Manual 110, Section 9.7
    """
    specifications = {
        "earthwork": [],
        "structures": [],
        "revegetation": [],
        "sequencing": [],
        "quality_control": []
    }
    
    # Earthwork specifications
    specifications["earthwork"].extend([
        f"Channel width: {channel_geometry.width:.1f}m ± 0.5m",
        f"Channel depth: {channel_geometry.depth:.1f}m ± 0.2m",
        f"Bank slopes: {channel_geometry.bank_angle}° ± 5°",
        f"Channel gradient: {channel_geometry.slope:.4f} ± 0.0005",
        "Excavated material to be stockpiled for backfill",
        "Minimize disturbance to existing vegetation"
    ])
    
    # Structure specifications
    for structure in habitat_structures:
        specs = [f"{structure.structure_type.value} installation:"]
        for dim_name, dim_value in structure.dimensions.items():
            specs.append(f"  {dim_name}: {dim_value}")
        specs.extend([f"  Materials: {', '.join(structure.materials)}"])
        specifications["structures"].extend(specs)
    
    # Revegetation specifications
    specifications["revegetation"].extend([
        "Native species only - obtain from local seed sources",
        "Plant during appropriate season (spring/fall)",
        "Erosion control fabric during establishment period",
        "Maintenance for minimum 2 growing seasons",
        "Achieve 80% survival rate by end of establishment period"
    ])
    
    # Construction sequencing
    work_window = construction_constraints.get("work_window", "all_year")
    if work_window == "low_flow_only":
        specifications["sequencing"].extend([
            "Schedule work during low-flow season",
            "Install upstream cofferdam if needed",
            "Dewater work areas progressively"
        ])
    
    specifications["sequencing"].extend([
        "1. Install erosion and sediment controls",
        "2. Excavate new channel alignment",
        "3. Install grade control/habitat structures", 
        "4. Connect old and new channels",
        "5. Install bank protection/revegetation",
        "6. Remove temporary diversions",
        "7. Monitor initial performance"
    ])
    
    # Quality control
    specifications["quality_control"].extend([
        "Survey channel geometry at 50m intervals",
        "Document as-built conditions with GPS",
        "Photo documentation of all structures",
        "Materials testing for imported materials",
        "Weekly progress inspections",
        "Post-construction performance monitoring"
    ])
    
    # Site-specific adjustments
    if construction_constraints.get("urban_proximity", False):
        specifications["quality_control"].extend([
            "Noise monitoring during construction",
            "Dust control measures",
            "Traffic management plan"
        ])
    
    return specifications


def estimate_project_costs(channel_geometry: ChannelGeometry,
                         project_length: float,
                         habitat_structures: List[HabitatStructure],
                         site_conditions: Dict[str, str]) -> Dict[str, float]:
    """
    Estimate project costs for restoration design.
    
    Args:
        channel_geometry: Designed channel geometry
        project_length: Length of restoration reach (m)
        habitat_structures: List of habitat structures
        site_conditions: Site access, materials availability, etc.
        
    Returns:
        Dictionary with cost estimates by category
        
    Reference: Based on typical restoration project costs
    """
    validate_positive(project_length, "project_length")
    
    costs = {
        "earthwork": 0.0,
        "structures": 0.0,
        "revegetation": 0.0,
        "design_engineering": 0.0,
        "construction_management": 0.0,
        "monitoring": 0.0,
        "contingency": 0.0,
        "total": 0.0
    }
    
    # Earthwork costs ($/m³)
    channel_volume = (channel_geometry.width * channel_geometry.depth * project_length)
    
    earthwork_unit_cost = 15.0  # $/m³ base cost
    if site_conditions.get("access", "good") == "difficult":
        earthwork_unit_cost *= 1.5
    if site_conditions.get("rock_excavation", "none") == "some":
        earthwork_unit_cost *= 2.0
    
    costs["earthwork"] = channel_volume * earthwork_unit_cost
    
    # Structure costs
    structure_costs = {
        HabitatStructureType.SILLS: 5000,        # $ per structure
        HabitatStructureType.DEFLECTORS: 3000,
        HabitatStructureType.RANDOM_ROCKS: 1000,
        HabitatStructureType.COVERS: 2000,
        HabitatStructureType.LARGE_WOODY_DEBRIS: 1500
    }
    
    for structure in habitat_structures:
        costs["structures"] += structure_costs.get(structure.structure_type, 2500)
    
    # Revegetation costs ($/m²)
    revegetation_area = 2 * project_length * channel_geometry.width  # Both banks
    revegetation_unit_cost = 5.0  # $/m²
    costs["revegetation"] = revegetation_area * revegetation_unit_cost
    
    # Professional services (% of construction)
    construction_total = costs["earthwork"] + costs["structures"] + costs["revegetation"]
    costs["design_engineering"] = construction_total * 0.15
    costs["construction_management"] = construction_total * 0.10
    
    # Monitoring (5 years)
    costs["monitoring"] = 50000.0  # Base monitoring program
    
    # Contingency
    subtotal = sum(costs.values())
    costs["contingency"] = subtotal * 0.20
    
    # Total project cost
    costs["total"] = sum(costs.values())
    
    return costs