"""
Bed forms and flow resistance equations.

Implements equations for bed form geometry and flow resistance
from Chapter 2, Sections 2.7 and 2.8.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range

def dune_height_relation(
    flow_depth: Union[float, np.ndarray],
    empirical_coefficient: float = 0.3
) -> Union[float, np.ndarray]:
    """
    Calculate dune height using empirical relation.
    
    Equation 2.40: Hd = C * H
    
    Where:
    - Hd is dune height [m]
    - H is flow depth [m]  
    - C is empirical coefficient (typically 0.1-0.3)
    
    Parameters
    ----------
    flow_depth : float or array-like
        Flow depth [m]
    empirical_coefficient : float, optional
        Empirical coefficient [-], default 0.3
        
    Returns
    -------
    float or ndarray
        Dune height [m]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.7.2
    """
    validate_positive(flow_depth, "flow_depth")
    validate_range(empirical_coefficient, "empirical_coefficient", 0.05, 0.5)
    
    return empirical_coefficient * flow_depth


def dune_length_relation(
    flow_depth: Union[float, np.ndarray],
    dune_height: Union[float, np.ndarray] = None,
    length_to_height_ratio: float = 6.0
) -> Union[float, np.ndarray]:
    """
    Calculate dune length using empirical relations.
    
    Equation 2.41: Ld = C * H  or  Ld = α * Hd
    
    Where:
    - Ld is dune length [m]
    - H is flow depth [m]
    - Hd is dune height [m]
    - α is length-to-height ratio (typically 5-10)
    
    Parameters
    ----------
    flow_depth : float or array-like
        Flow depth [m]
    dune_height : float or array-like, optional
        Dune height [m]. If provided, uses length-height relation.
    length_to_height_ratio : float, optional
        Length-to-height ratio [-], default 6.0
        
    Returns
    -------
    float or ndarray
        Dune length [m]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.7.2
    """
    validate_positive(flow_depth, "flow_depth")
    validate_positive(length_to_height_ratio, "length_to_height_ratio")
    
    if dune_height is not None:
        validate_positive(dune_height, "dune_height")
        return length_to_height_ratio * dune_height
    else:
        # Use depth-based relation (Ld ≈ 6H)
        return length_to_height_ratio * flow_depth


def ripple_geometry(
    particle_diameter: Union[float, np.ndarray],
    shields_parameter: Union[float, np.ndarray]
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Calculate ripple height and length using empirical relations.
    
    For ripples (typically D < 0.6 mm and θ < 0.6):
    Hr ≈ 1000 * D
    Lr ≈ 1000 * D
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    shields_parameter : float or array-like
        Shields parameter [-]
        
    Returns
    -------
    tuple
        (ripple_height [m], ripple_length [m])
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.7.1
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(shields_parameter, "shields_parameter")
    
    # Check if conditions are appropriate for ripples
    D_mm = particle_diameter * 1000  # Convert to mm
    
    # Ripple height (limited by grain size)
    ripple_height = np.minimum(1000 * particle_diameter, 0.05)  # Max 5 cm
    
    # Ripple length
    ripple_length = np.maximum(500 * particle_diameter, 0.06)  # Min 6 cm
    
    # Modify based on Shields parameter
    ripple_height *= np.minimum(shields_parameter / 0.1, 1.0)
    ripple_length *= np.minimum(shields_parameter / 0.1, 1.0)
    
    return ripple_height, ripple_length


def bed_roughness_height(
    bed_form_height: Union[float, np.ndarray],
    bed_form_type: str = "dune",
    form_factor: float = 1.1
) -> Union[float, np.ndarray]:
    """
    Calculate bed roughness height from bed form geometry.
    
    Equation 2.42: ks = α * Hf
    
    Where:
    - ks is roughness height [m]
    - Hf is bed form height [m]
    - α is form factor (depends on bed form type)
    
    Parameters
    ----------
    bed_form_height : float or array-like
        Bed form height [m]
    bed_form_type : str, optional
        Type of bed form ("dune", "ripple", "antidune"), default "dune"
    form_factor : float, optional
        Form factor [-], default 1.1
        
    Returns
    -------
    float or ndarray
        Bed roughness height [m]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.8.2
    """
    validate_positive(bed_form_height, "bed_form_height")
    validate_positive(form_factor, "form_factor")
    
    # Adjust form factor based on bed form type
    if bed_form_type.lower() == "ripple":
        alpha = 1.0
    elif bed_form_type.lower() == "dune":
        alpha = 1.1
    elif bed_form_type.lower() == "antidune":
        alpha = 0.8
    else:
        alpha = form_factor
    
    return alpha * bed_form_height


def manning_roughness_coefficient(
    roughness_height: Union[float, np.ndarray],
    hydraulic_radius: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate Manning's roughness coefficient from bed roughness.
    
    Equation 2.43: n = (ks^(1/6)) / (8.5 * R^(1/6))
    
    Where:
    - n is Manning's coefficient [s/m^(1/3)]
    - ks is roughness height [m]
    - R is hydraulic radius [m]
    
    Parameters
    ----------
    roughness_height : float or array-like
        Bed roughness height [m]
    hydraulic_radius : float or array-like
        Hydraulic radius [m]
        
    Returns
    -------
    float or ndarray
        Manning's roughness coefficient [s/m^(1/3)]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.8.3
    """
    validate_positive(roughness_height, "roughness_height")
    validate_positive(hydraulic_radius, "hydraulic_radius")
    
    return (roughness_height**(1/6)) / (8.5 * hydraulic_radius**(1/6))


def darcy_weisbach_friction_factor(
    roughness_height: Union[float, np.ndarray],
    hydraulic_radius: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate Darcy-Weisbach friction factor.
    
    Equation 2.44: f = 8g / [5.75 * log10(12.2 * R/ks)]²
    
    Where:
    - f is Darcy-Weisbach friction factor [-]
    - R is hydraulic radius [m]
    - ks is roughness height [m]
    
    Parameters
    ----------
    roughness_height : float or array-like
        Bed roughness height [m]
    hydraulic_radius : float or array-like
        Hydraulic radius [m]
        
    Returns
    -------
    float or ndarray
        Darcy-Weisbach friction factor [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.8.1
    """
    validate_positive(roughness_height, "roughness_height")
    validate_positive(hydraulic_radius, "hydraulic_radius")
    
    # Relative roughness
    relative_roughness = hydraulic_radius / roughness_height
    
    # Avoid log of negative numbers
    relative_roughness = np.maximum(relative_roughness, 1.1)
    
    denominator = 5.75 * np.log10(12.2 * relative_roughness)
    
    return 8 * 9.81 / denominator**2


def bed_form_migration_rate(
    bed_load_discharge: Union[float, np.ndarray],
    bed_form_height: Union[float, np.ndarray],
    porosity: Union[float, np.ndarray] = 0.4
) -> Union[float, np.ndarray]:
    """
    Estimate bed form migration rate.
    
    Equation 2.45: Vf = qb / [(1-n) * Hf]
    
    Where:
    - Vf is bed form migration velocity [m/s]
    - qb is bed load discharge [m²/s]
    - n is bed porosity [-]
    - Hf is bed form height [m]
    
    Parameters
    ----------
    bed_load_discharge : float or array-like
        Volumetric bed load discharge per unit width [m²/s]
    bed_form_height : float or array-like
        Bed form height [m]
    porosity : float or array-like, optional
        Bed porosity [-], default 0.4
        
    Returns
    -------
    float or ndarray
        Bed form migration velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.7.3
    """
    validate_positive(bed_load_discharge, "bed_load_discharge")
    validate_positive(bed_form_height, "bed_form_height")
    validate_range(porosity, "porosity", 0.0, 1.0)
    
    return bed_load_discharge / ((1 - porosity) * bed_form_height)