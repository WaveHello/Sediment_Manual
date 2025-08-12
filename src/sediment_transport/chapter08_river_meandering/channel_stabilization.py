"""
Channel Stabilization and Engineering Applications

This module implements engineering methods for channel stabilization,
including submerged vanes, revetments, and flow training structures.

Based on ASCE Manual 110, Chapter 8, Section 8.7
Applications of flow and stability relations for river engineering practice.

References:
    ASCE Manual 110, Chapter 8, Section 8.7
    Odgaard & Kennedy (1983), Pokrefke (1993), Biedenharn et al. (1997)
"""

import numpy as np
from typing import Union, Tuple, Dict, List, Optional
from ..utils.validators import validate_positive, validate_range

def submerged_vane_design(channel_width_b: float,
                        design_depth_d: float,
                        design_velocity_u: float,
                        radius_curvature_rc: float,
                        vane_angle: float = 20.0,
                        vane_height_ratio: float = 0.3) -> Dict[str, Union[float, int]]:
    """
    Design submerged vanes for channel stabilization in bends.
    
    Submerged vanes generate secondary circulation to counteract
    centrifugally induced secondary currents that cause bank erosion.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    design_depth_d : float
        Design water depth (m)
    design_velocity_u : float
        Design velocity (m/s)
    radius_curvature_rc : float
        Radius of curvature (m)
    vane_angle : float, optional
        Vane angle of attack (degrees, default=20°)
    vane_height_ratio : float, optional
        Vane height to depth ratio (default=0.3)
        
    Returns:
    --------
    dict
        Vane design parameters
        
    Reference: Odgaard & Kennedy (1983)
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(design_depth_d, "design_depth_d")
    validate_positive(design_velocity_u, "design_velocity_u")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_range(vane_angle, 10.0, 30.0, "vane_angle")
    validate_range(vane_height_ratio, 0.2, 0.5, "vane_height_ratio")
    
    # Vane dimensions
    vane_height = vane_height_ratio * design_depth_d
    vane_length = min(channel_width_b / 10.0, 2.0 * design_depth_d)  # Typical length
    
    # Vane spacing
    # Based on secondary circulation development length
    vane_spacing = 3.0 * channel_width_b  # Typical spacing
    
    # Number of vanes for bend protection
    bend_arc_length = np.pi * radius_curvature_rc  # Quarter circle approximation
    number_of_vanes = max(1, int(bend_arc_length / vane_spacing))
    
    # Vane placement - distance from outer bank
    bank_offset = 0.3 * channel_width_b  # Place vanes away from outer bank
    
    # Secondary circulation strength (estimated)
    # Based on centrifugal acceleration
    centrifugal_acceleration = design_velocity_u**2 / radius_curvature_rc
    circulation_strength = centrifugal_acceleration * vane_length * np.sin(np.radians(vane_angle))
    
    # Vane effectiveness factor
    width_radius_ratio = channel_width_b / radius_curvature_rc
    effectiveness_factor = min(1.0, 2.0 * width_radius_ratio)  # Empirical relationship
    
    return {
        'vane_height': vane_height,
        'vane_length': vane_length,
        'vane_angle': vane_angle,
        'vane_spacing': vane_spacing,
        'number_of_vanes': number_of_vanes,
        'bank_offset': bank_offset,
        'circulation_strength': circulation_strength,
        'effectiveness_factor': effectiveness_factor,
        'total_bend_length': bend_arc_length
    }

def revetment_design(bank_height_H: float,
                   design_velocity_u: float,
                   stone_specific_gravity: float = 2.65,
                   safety_factor: float = 1.5,
                   side_slope: float = 1.5) -> Dict[str, float]:
    """
    Design riprap revetment for bank protection.
    
    Parameters:
    -----------
    bank_height_H : float
        Bank height (m)
    design_velocity_u : float
        Design velocity near bank (m/s)
    stone_specific_gravity : float, optional
        Specific gravity of riprap stone (default=2.65)
    safety_factor : float, optional
        Safety factor for stone size (default=1.5)
    side_slope : float, optional
        Side slope ratio (horizontal:vertical, default=1.5:1)
        
    Returns:
    --------
    dict
        Revetment design parameters
        
    Reference: Petersen (1986), ASCE (1998)
    """
    validate_positive(bank_height_H, "bank_height_H")
    validate_positive(design_velocity_u, "design_velocity_u")
    validate_positive(stone_specific_gravity, "stone_specific_gravity")
    validate_positive(safety_factor, "safety_factor")
    validate_positive(side_slope, "side_slope")
    
    # Stone size calculation (Isbash formula modified)
    gravity = 9.81
    stone_density = stone_specific_gravity * 1000.0  # kg/m³
    water_density = 1000.0  # kg/m³
    
    # Critical velocity for stone stability
    relative_density = (stone_density - water_density) / water_density
    
    # Minimum stone diameter
    stone_diameter = (design_velocity_u**2) / (2 * gravity * relative_density) * safety_factor
    
    # Stone weight (spherical approximation)
    stone_volume = (np.pi / 6.0) * stone_diameter**3
    stone_weight = stone_volume * stone_density
    
    # Revetment thickness
    revetment_thickness = 1.5 * stone_diameter  # Typically 1.5 to 2 times stone diameter
    
    # Filter layer design
    filter_diameter = stone_diameter / 10.0  # Filter stone size
    filter_thickness = 0.3  # Standard filter thickness (m)
    
    # Toe protection
    toe_depth = max(0.5, stone_diameter)  # Toe protection depth
    toe_width = 2.0 * stone_diameter  # Toe protection width
    
    # Area calculations per unit length
    slope_length = np.sqrt(bank_height_H**2 + (side_slope * bank_height_H)**2)
    revetment_area_per_meter = slope_length * revetment_thickness
    
    # Volume calculations per meter of bank
    stone_volume_per_meter = revetment_area_per_meter * 0.65  # Account for voids
    filter_volume_per_meter = slope_length * filter_thickness * 0.65
    
    return {
        'stone_diameter': stone_diameter,
        'stone_weight': stone_weight,
        'revetment_thickness': revetment_thickness,
        'filter_diameter': filter_diameter,
        'filter_thickness': filter_thickness,
        'toe_depth': toe_depth,
        'toe_width': toe_width,
        'slope_length': slope_length,
        'stone_volume_per_meter': stone_volume_per_meter,
        'filter_volume_per_meter': filter_volume_per_meter,
        'total_volume_per_meter': stone_volume_per_meter + filter_volume_per_meter
    }

def bendway_weir_design(channel_width_b: float,
                       bank_full_depth_d: float,
                       design_discharge_Q: float,
                       radius_curvature_rc: float,
                       weir_angle: float = 70.0) -> Dict[str, Union[float, int]]:
    """
    Design bendway weirs for channel stabilization.
    
    Bendway weirs are rock structures oriented upstream at an angle
    to redirect flow away from the outer bank.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    bank_full_depth_d : float
        Bank-full depth (m)
    design_discharge_Q : float
        Design discharge (m³/s)
    radius_curvature_rc : float
        Radius of curvature (m)
    weir_angle : float, optional
        Weir angle with bank (degrees, default=70°)
        
    Returns:
    --------
    dict
        Bendway weir design parameters
        
    Reference: Pokrefke (1993), USACE (2002)
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(bank_full_depth_d, "bank_full_depth_d")
    validate_positive(design_discharge_Q, "design_discharge_Q")
    validate_positive(radius_curvature_rc, "radius_curvature_rc")
    validate_range(weir_angle, 60.0, 80.0, "weir_angle")
    
    # Weir dimensions
    weir_length = 0.25 * channel_width_b  # Typical length: 25% of channel width
    weir_height = 0.5 * bank_full_depth_d  # Submerged weir
    weir_crest_width = 3.0  # Typical crest width (m)
    
    # Side slopes (typical for rock weirs)
    upstream_slope = 2.0  # 2H:1V
    downstream_slope = 1.5  # 1.5H:1V
    
    # Weir spacing
    # Based on flow redirection length
    redirection_length = 2.0 * channel_width_b
    weir_spacing = redirection_length
    
    # Number of weirs for bend
    bend_arc_length = np.pi * radius_curvature_rc / 2.0  # Quarter circle
    number_of_weirs = max(1, int(bend_arc_length / weir_spacing))
    
    # Rock size for weir construction
    design_velocity = design_discharge_Q / (channel_width_b * bank_full_depth_d)
    
    # Use revetment design function for rock sizing
    rock_design = revetment_design(weir_height, design_velocity * 1.2)  # 20% higher velocity
    rock_diameter = rock_design['stone_diameter']
    
    # Volume calculations per weir
    weir_volume = (weir_length * weir_crest_width * weir_height + 
                  0.5 * weir_length * weir_height * 
                  (upstream_slope + downstream_slope) * weir_height)
    
    # Total rock volume (accounting for voids)
    total_rock_volume = weir_volume * 1.4  # Add 40% for voids and filter
    
    return {
        'weir_length': weir_length,
        'weir_height': weir_height,
        'weir_crest_width': weir_crest_width,
        'weir_angle': weir_angle,
        'upstream_slope': upstream_slope,
        'downstream_slope': downstream_slope,
        'weir_spacing': weir_spacing,
        'number_of_weirs': number_of_weirs,
        'rock_diameter': rock_diameter,
        'weir_volume': weir_volume,
        'total_rock_volume': total_rock_volume,
        'total_project_volume': total_rock_volume * number_of_weirs
    }

def dike_design(channel_width_b: float,
               design_depth_d: float,
               design_velocity_u: float,
               deflection_angle: float = 30.0,
               dike_permeability: float = 0.4) -> Dict[str, float]:
    """
    Design permeable dikes for flow redirection.
    
    Parameters:
    -----------
    channel_width_b : float
        Channel width (m)
    design_depth_d : float
        Design depth (m)
    design_velocity_u : float
        Design velocity (m/s)
    deflection_angle : float, optional
        Flow deflection angle (degrees, default=30°)
    dike_permeability : float, optional
        Dike permeability factor (0-1, default=0.4)
        
    Returns:
    --------
    dict
        Dike design parameters
    """
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(design_depth_d, "design_depth_d")
    validate_positive(design_velocity_u, "design_velocity_u")
    validate_range(deflection_angle, 10.0, 45.0, "deflection_angle")
    validate_range(dike_permeability, 0.1, 0.8, "dike_permeability")
    
    # Dike length based on desired flow deflection
    dike_length = 0.3 * channel_width_b * np.tan(np.radians(deflection_angle))
    
    # Dike height (typically extends above water surface)
    dike_height = design_depth_d * 1.2
    
    # Dike crest width
    dike_crest_width = max(2.0, dike_height / 3.0)
    
    # Side slopes
    water_side_slope = 1.5  # 1.5H:1V
    land_side_slope = 2.0   # 2H:1V
    
    # Rock size calculation
    rock_design = revetment_design(dike_height, design_velocity_u)
    rock_diameter = rock_design['stone_diameter']
    
    # Adjust for permeability
    if dike_permeability > 0.3:
        # Permeable dike - can use smaller stones
        rock_diameter *= 0.8
    
    # Volume calculations
    cross_sectional_area = (dike_crest_width * dike_height + 
                           0.5 * dike_height**2 * (water_side_slope + land_side_slope))
    
    dike_volume = cross_sectional_area * dike_length
    
    # Rock volume (accounting for voids)
    rock_volume = dike_volume * (1.0 - dike_permeability) * 1.3  # Add for filter
    
    return {
        'dike_length': dike_length,
        'dike_height': dike_height,
        'dike_crest_width': dike_crest_width,
        'water_side_slope': water_side_slope,
        'land_side_slope': land_side_slope,
        'rock_diameter': rock_diameter,
        'cross_sectional_area': cross_sectional_area,
        'dike_volume': dike_volume,
        'rock_volume': rock_volume,
        'permeability': dike_permeability
    }

def channel_alignment_design(valley_slope_Sv: float,
                           design_discharge_Q: float,
                           channel_width_b: float,
                           available_corridor_width: float,
                           stability_criterion: str = 'leopold_wolman') -> Dict[str, float]:
    """
    Design stable channel alignment using stability theory.
    
    Parameters:
    -----------
    valley_slope_Sv : float
        Valley slope (dimensionless)
    design_discharge_Q : float
        Design discharge (m³/s)
    channel_width_b : float
        Proposed channel width (m)
    available_corridor_width : float
        Available corridor width for meandering (m)
    stability_criterion : str, optional
        Stability criterion to use (default='leopold_wolman')
        
    Returns:
    --------
    dict
        Recommended channel alignment parameters
    """
    validate_positive(valley_slope_Sv, "valley_slope_Sv")
    validate_positive(design_discharge_Q, "design_discharge_Q")
    validate_positive(channel_width_b, "channel_width_b")
    validate_positive(available_corridor_width, "available_corridor_width")
    
    # Import stability analysis functions
    from .meandering_criteria import classify_channel_pattern
    from .planform_geometry import planform_characteristics
    from .stability_analysis import dominant_wavelength_calculation
    
    # Check channel pattern tendency
    pattern, threshold_slope, confidence = classify_channel_pattern(
        valley_slope_Sv, design_discharge_Q, stability_criterion
    )
    
    # Calculate natural planform characteristics
    planform = planform_characteristics(channel_width_b)
    
    # Check if natural meandering fits available corridor
    natural_amplitude = planform['amplitude']
    corridor_constraint = available_corridor_width / 2.0  # Half-width for amplitude
    
    if natural_amplitude > corridor_constraint:
        # Need to constrain meandering
        constrained_amplitude = 0.8 * corridor_constraint  # Safety factor
        
        # Calculate required wavelength for constrained amplitude
        # Using A = 3.0 * b^1.1 relationship inversely
        amplitude_ratio = constrained_amplitude / natural_amplitude
        wavelength_adjustment = amplitude_ratio**(1.01/1.1)  # Scale wavelength proportionally
        
        recommended_wavelength = planform['wavelength'] * wavelength_adjustment
        recommended_amplitude = constrained_amplitude
        
        # Calculate required sinuosity
        recommended_sinuosity = 1.0 + (2 * np.pi * recommended_amplitude / recommended_wavelength)**2 / 6.0
        
        # Channel slope with sinuosity
        channel_slope = valley_slope_Sv / recommended_sinuosity
        
        alignment_type = 'constrained_meandering'
        
    else:
        # Natural meandering can be accommodated
        recommended_wavelength = planform['wavelength']
        recommended_amplitude = planform['amplitude']
        recommended_sinuosity = 1.0 + (2 * np.pi * recommended_amplitude / recommended_wavelength)**2 / 6.0
        
        channel_slope = valley_slope_Sv / recommended_sinuosity
        
        alignment_type = 'natural_meandering'
    
    # Check if resulting slope is stable
    revised_pattern, _, _ = classify_channel_pattern(
        channel_slope, design_discharge_Q, stability_criterion
    )
    
    # Calculate minimum radius of curvature
    min_radius_curvature = channel_width_b * 2.4  # From Leopold-Wolman relationships
    
    return {
        'alignment_type': alignment_type,
        'recommended_wavelength': recommended_wavelength,
        'recommended_amplitude': recommended_amplitude,
        'recommended_sinuosity': recommended_sinuosity,
        'channel_slope': channel_slope,
        'valley_slope': valley_slope_Sv,
        'min_radius_curvature': min_radius_curvature,
        'natural_pattern': pattern,
        'revised_pattern': revised_pattern,
        'stability_confidence': confidence,
        'corridor_utilization': (2 * recommended_amplitude) / available_corridor_width
    }

def biotechnical_stabilization(bank_height_H: float,
                             bank_angle: float,
                             soil_type: str,
                             design_velocity_u: float,
                             climate_zone: str = 'temperate') -> Dict[str, Union[str, float, List[str]]]:
    """
    Design biotechnical bank stabilization system.
    
    Combines vegetation with structural elements for sustainable stabilization.
    
    Parameters:
    -----------
    bank_height_H : float
        Bank height (m)
    bank_angle : float
        Bank angle from horizontal (degrees)
    soil_type : str
        Soil type: 'cohesive', 'granular', 'mixed'
    design_velocity_u : float
        Near-bank velocity (m/s)
    climate_zone : str, optional
        Climate zone for vegetation selection
        
    Returns:
    --------
    dict
        Biotechnical stabilization design
    """
    validate_positive(bank_height_H, "bank_height_H")
    validate_range(bank_angle, 0.0, 60.0, "bank_angle")
    
    # Determine if bank regrading is needed
    max_stable_angle = {'cohesive': 45.0, 'granular': 35.0, 'mixed': 40.0}
    stable_angle = max_stable_angle.get(soil_type, 35.0)
    
    if bank_angle > stable_angle:
        regrading_required = True
        recommended_angle = stable_angle
        cut_volume = bank_height_H**2 * (np.tan(np.radians(bank_angle)) - 
                                        np.tan(np.radians(stable_angle))) / 2.0
    else:
        regrading_required = False
        recommended_angle = bank_angle
        cut_volume = 0.0
    
    # Vegetation selection based on climate and velocity
    if climate_zone == 'temperate':
        if design_velocity_u < 1.0:
            vegetation_types = ['willow', 'cottonwood', 'native_grasses']
        elif design_velocity_u < 2.0:
            vegetation_types = ['willow_stakes', 'brush_bundles', 'native_shrubs']
        else:
            vegetation_types = ['live_stakes', 'brush_layers', 'erosion_control_fabric']
    elif climate_zone == 'arid':
        vegetation_types = ['drought_tolerant_shrubs', 'native_grasses', 'erosion_control_matting']
    else:
        vegetation_types = ['adapted_local_species', 'erosion_control_fabric']
    
    # Structural elements needed
    structural_elements = []
    
    if design_velocity_u > 1.5:
        structural_elements.append('toe_protection')
    
    if bank_height_H > 3.0:
        structural_elements.append('terracing')
    
    if soil_type == 'granular':
        structural_elements.append('geotextile_reinforcement')
    
    # Maintenance requirements
    maintenance_frequency = 'annually' if design_velocity_u > 1.0 else 'every_2_years'
    
    establishment_time = {
        'grasses': '1-2_growing_seasons',
        'shrubs': '2-3_growing_seasons',
        'trees': '3-5_growing_seasons'
    }
    
    return {
        'regrading_required': regrading_required,
        'recommended_angle': recommended_angle,
        'cut_volume_per_meter': cut_volume,
        'vegetation_types': vegetation_types,
        'structural_elements': structural_elements,
        'maintenance_frequency': maintenance_frequency,
        'establishment_time': establishment_time,
        'sustainability_rating': 'high' if design_velocity_u < 1.0 else 'moderate'
    }