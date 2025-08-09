"""
Sediment Characterization for Fine-Grained Sediment (Chapter 4)

This module implements characterization methods for fine-grained cohesive sediment,
including size classification, cohesion assessment, and key properties from 
Tables 4-1, 4-2, and 4-3 from ASCE Manual 110, Chapter 4.
"""

import numpy as np
from typing import Union, Optional, Dict, Tuple, Any
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

# Clay mineral properties from Table 4-2
CLAY_MINERAL_PROPERTIES = {
    'kaolinite': {
        'nominal_diameter': 0.36e-6,  # m
        'cec_range': (3, 15),  # meq/100g
        'critical_salinity': 0.6  # ppt
    },
    'illite': {
        'nominal_diameter': 0.062e-6,  # m
        'cec_range': (10, 40),  # meq/100g
        'critical_salinity': 1.1  # ppt
    },
    'chlorite': {
        'nominal_diameter': 0.062e-6,  # m
        'cec_range': (24, 35),  # meq/100g
        'critical_salinity': None  # Not reported
    },
    'smectite': {
        'nominal_diameter': 0.011e-6,  # m
        'cec_range': (80, 150),  # meq/100g
        'critical_salinity': 2.4  # ppt
    },
    'montmorillonite': {  # Alias for smectite
        'nominal_diameter': 0.011e-6,  # m
        'cec_range': (80, 150),  # meq/100g
        'critical_salinity': 2.4  # ppt
    }
}

def classify_sediment_by_size(
    particle_diameter: Union[float, np.ndarray]
) -> Union[str, np.ndarray]:
    """
    Classify sediment particles by size according to Table 4-1.
    
    Size Classification (Table 4-1):
    - >62 μm: Coarse-grained (cohesionless)
    - 40-62 μm: Fine-grained coarse silt (practically cohesionless)
    - 20-40 μm: Fine-grained coarse silt (cohesion increasingly important)
    - 2-20 μm: Fine-grained medium and fine silt (cohesion important)
    - <2 μm: Fine-grained clay (cohesion very important)
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
        
    Returns
    -------
    str or ndarray
        Size classification string(s)
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-1
    """
    validate_positive(particle_diameter, "particle_diameter")
    
    # Convert to micrometers for classification
    d_um = np.atleast_1d(particle_diameter) * 1e6
    
    classification = np.empty(d_um.shape, dtype=object)
    
    # Apply size classifications
    classification[d_um > 62] = "coarse_grained"
    classification[(d_um >= 40) & (d_um <= 62)] = "fine_grained_coarse_silt_low_cohesion"
    classification[(d_um >= 20) & (d_um < 40)] = "fine_grained_coarse_silt_moderate_cohesion"
    classification[(d_um >= 2) & (d_um < 20)] = "fine_grained_silt_high_cohesion"
    classification[d_um < 2] = "fine_grained_clay_very_high_cohesion"
    
    if np.isscalar(particle_diameter):
        return classification.item()
    else:
        return classification


def assess_cohesion_degree(
    particle_diameter: Union[float, np.ndarray]
) -> Union[str, np.ndarray]:
    """
    Assess degree of cohesion based on particle size (Table 4-1).
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
        
    Returns
    -------
    str or ndarray
        Cohesion degree classification
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-1
    """
    validate_positive(particle_diameter, "particle_diameter")
    
    d_um = np.atleast_1d(particle_diameter) * 1e6
    
    cohesion_degree = np.empty(d_um.shape, dtype=object)
    
    cohesion_degree[d_um > 62] = "cohesionless"
    cohesion_degree[(d_um >= 40) & (d_um <= 62)] = "practically_cohesionless"
    cohesion_degree[(d_um >= 20) & (d_um < 40)] = "cohesion_increasingly_important"
    cohesion_degree[(d_um >= 2) & (d_um < 20)] = "cohesion_important"
    cohesion_degree[d_um < 2] = "cohesion_very_important"
    
    if np.isscalar(particle_diameter):
        return cohesion_degree.item()
    else:
        return cohesion_degree


def get_clay_mineral_properties(
    mineral_type: str
) -> Dict[str, Any]:
    """
    Get properties of clay minerals from Table 4-2.
    
    Parameters
    ----------
    mineral_type : str
        Type of clay mineral ('kaolinite', 'illite', 'chlorite', 'smectite', 'montmorillonite')
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing mineral properties
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-2
    """
    mineral_type = mineral_type.lower()
    
    if mineral_type not in CLAY_MINERAL_PROPERTIES:
        available_types = list(CLAY_MINERAL_PROPERTIES.keys())
        raise ValueError(f"Unknown mineral type '{mineral_type}'. Available types: {available_types}")
    
    return CLAY_MINERAL_PROPERTIES[mineral_type].copy()


def estimate_cation_exchange_capacity(
    mineral_type: str,
    use_midpoint: bool = True
) -> float:
    """
    Estimate cation exchange capacity for clay minerals.
    
    Parameters
    ----------
    mineral_type : str
        Type of clay mineral
    use_midpoint : bool, optional
        Use midpoint of CEC range, default True
        
    Returns
    -------
    float
        Estimated CEC [meq/100g]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-2
    """
    properties = get_clay_mineral_properties(mineral_type)
    cec_range = properties['cec_range']
    
    if use_midpoint:
        return (cec_range[0] + cec_range[1]) / 2.0
    else:
        return cec_range


def critical_salinity_for_flocculation(
    mineral_type: str
) -> Optional[float]:
    """
    Get critical salinity for flocculation of clay minerals.
    
    Parameters
    ----------
    mineral_type : str
        Type of clay mineral
        
    Returns
    -------
    float or None
        Critical salinity [ppt] or None if not reported
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-2
    """
    properties = get_clay_mineral_properties(mineral_type)
    return properties['critical_salinity']


def floc_density_estimate(
    floc_diameter: Union[float, np.ndarray],
    primary_particle_density: float = 2650.0,
    water_density: float = WATER_DENSITY,
    porosity: Optional[float] = None
) -> Union[float, np.ndarray]:
    """
    Estimate floc density considering water content.
    
    Example from text: A floc with density 1,090 kg/m³ composed of clay particles
    with density 2,650 kg/m³ contains nearly 95% water by volume.
    
    Parameters
    ----------
    floc_diameter : float or array-like
        Floc diameter [m]
    primary_particle_density : float, optional
        Density of primary particles [kg/m³], default 2650
    water_density : float, optional
        Water density [kg/m³], default 1000
    porosity : float, optional
        Floc porosity (water fraction by volume), default calculated
        
    Returns
    -------
    float or ndarray
        Estimated floc density [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    """
    validate_positive(floc_diameter, "floc_diameter")
    validate_positive(primary_particle_density, "primary_particle_density")
    
    df = np.atleast_1d(floc_diameter)
    
    if porosity is None:
        # Empirical relationship: larger flocs have higher porosity
        # Based on fractal concepts and observations
        df_um = df * 1e6  # Convert to micrometers
        porosity = 0.5 + 0.4 * np.tanh((df_um - 50) / 100)  # Ranges from ~0.5 to ~0.9
        porosity = np.clip(porosity, 0.1, 0.98)
    
    # Floc density from porosity
    rho_floc = porosity * water_density + (1 - porosity) * primary_particle_density
    
    return rho_floc.item() if np.isscalar(floc_diameter) else rho_floc


def migniot_floc_settling_ratio(
    particle_diameter: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate floc to particle settling velocity ratio based on Migniot (1968).
    
    Empirical relationship showing ratio increases from ~1 at 40 μm 
    to 300-400 at 1 μm due to flocculation effects.
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Primary particle diameter [m]
        
    Returns
    -------
    float or ndarray
        Ratio of floc settling velocity to particle settling velocity
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    Migniot (1968)
    """
    validate_positive(particle_diameter, "particle_diameter")
    
    d_um = np.atleast_1d(particle_diameter) * 1e6  # Convert to micrometers
    
    # Empirical fit to Migniot's data
    # Ratio increases significantly as particle size decreases
    ratio = 1.0 + 350 * np.exp(-(d_um - 1) / 10) * (40 / d_um)**2
    ratio = np.clip(ratio, 1.0, 400.0)  # Physical bounds
    
    return ratio.item() if np.isscalar(particle_diameter) else ratio


def characterization_test_requirements(
    organic_content: float
) -> Dict[str, str]:
    """
    Provide test requirements based on organic content (Table 4-3).
    
    Parameters
    ----------
    organic_content : float
        Organic content [%]
        
    Returns
    -------
    Dict[str, str]
        Dictionary of test requirements and recommendations
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-3
    """
    validate_range(organic_content, "organic_content", 0.0, 100.0)
    
    recommendations = {
        'particle_size_test': 'Use hydrometer test or settling column bottom withdrawal',
        'sample_preparation': 'Use naturally wet samples, do not air-dry initially',
        'settling_velocity': 'Measure in multiport settling column with untreated sample',
        'mineral_composition': 'Use X-ray diffraction for clay and non-clay minerals',
        'organic_content': 'Measure loss on ignition or total organic carbon'
    }
    
    if organic_content > 10:
        recommendations['particle_size_test'] = 'DO NOT perform standard particle size test'
        recommendations['alternative_test'] = 'Measure settling velocity of untreated sample instead'
        recommendations['cec_warning'] = 'CEC may be excessively high, not representative of clay'
    
    recommendations['salinity'] = 'Report if less than 10 ppt (higher salinities have minor effect)'
    
    return recommendations


def salinity_effect_significance(
    salinity: float
) -> Tuple[bool, str]:
    """
    Determine significance of salinity effect on flocculation.
    
    Based on Table 4-3: At salinities >10 ppt, influence on floc structure 
    is comparatively minor.
    
    Parameters
    ----------
    salinity : float
        Water salinity [ppt]
        
    Returns
    -------
    Tuple[bool, str]
        (is_significant, description)
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-3
    Krone (1962)
    """
    validate_positive(salinity, "salinity")
    
    if salinity < 10.0:
        return True, f"Salinity ({salinity:.1f} ppt) significantly affects flocculation"
    else:
        return False, f"Salinity ({salinity:.1f} ppt) has minor effect on flocculation"


def mud_classification_rheological(
    particle_size: float,
    shows_yield_stress: bool,
    viscosity_behavior: str
) -> Dict[str, Any]:
    """
    Classify mud based on rheological definition (Section 4.2.2).
    
    Mehta (2002) definition: Mud is a sediment-water mixture of grains 
    predominantly <63 μm, exhibiting poroelastic or viscoelastic behavior
    when particle-supported, and highly viscous non-Newtonian when fluid-like.
    
    Parameters
    ----------
    particle_size : float
        Characteristic particle size [m]
    shows_yield_stress : bool
        Whether material shows yield stress
    viscosity_behavior : str
        'newtonian', 'shear_thinning', 'shear_thickening'
        
    Returns
    -------
    Dict[str, Any]
        Classification results
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.2.2
    Mehta (2002)
    """
    validate_positive(particle_size, "particle_size")
    
    d_um = particle_size * 1e6
    
    classification = {
        'is_mud': False,
        'size_criterion': d_um < 63,
        'particle_size_um': d_um,
        'rheological_state': 'unknown'
    }
    
    if d_um < 63:
        classification['is_mud'] = True
        
        if shows_yield_stress:
            if viscosity_behavior == 'newtonian':
                classification['rheological_state'] = 'bingham_plastic'
            else:
                classification['rheological_state'] = 'viscoplastic'
        else:
            if viscosity_behavior == 'newtonian':
                classification['rheological_state'] = 'newtonian_fluid'
            elif viscosity_behavior == 'shear_thinning':
                classification['rheological_state'] = 'pseudoplastic'
            elif viscosity_behavior == 'shear_thickening':
                classification['rheological_state'] = 'dilatant'
    
    return classification


def specific_surface_area_estimate(
    particle_diameter: Union[float, np.ndarray],
    particle_density: float = 2650.0,
    shape_factor: float = 6.0
) -> Union[float, np.ndarray]:
    """
    Estimate specific surface area for particles.
    
    Specific surface area = surface area / particle weight
    For spherical particles: SSA = 6/(ρ*d)
    For plate-like clay particles, shape factor adjusts for non-spherical geometry.
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    particle_density : float, optional
        Particle density [kg/m³], default 2650
    shape_factor : float, optional
        Shape factor (6.0 for spheres, higher for plates), default 6.0
        
    Returns
    -------
    float or ndarray
        Specific surface area [m²/kg]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(particle_density, "particle_density")
    validate_positive(shape_factor, "shape_factor")
    
    d = np.atleast_1d(particle_diameter)
    ssa = shape_factor / (particle_density * d)
    
    return ssa.item() if np.isscalar(particle_diameter) else ssa


def cohesion_surface_area_relationship(
    specific_surface_area: Union[float, np.ndarray]
) -> Union[str, np.ndarray]:
    """
    Relate cohesion to specific surface area.
    
    Higher specific surface area (smaller, more platy particles) leads to 
    higher cohesion due to electrochemical forces.
    
    Parameters
    ----------
    specific_surface_area : float or array-like
        Specific surface area [m²/kg]
        
    Returns
    -------
    str or ndarray
        Qualitative cohesion assessment
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    """
    validate_positive(specific_surface_area, "specific_surface_area")
    
    ssa = np.atleast_1d(specific_surface_area)
    
    cohesion_level = np.empty(ssa.shape, dtype=object)
    
    # Typical ranges for different materials
    cohesion_level[ssa < 1] = "very_low_cohesion"  # Coarse particles
    cohesion_level[(ssa >= 1) & (ssa < 10)] = "low_cohesion"  # Fine sand/coarse silt
    cohesion_level[(ssa >= 10) & (ssa < 100)] = "moderate_cohesion"  # Silt
    cohesion_level[(ssa >= 100) & (ssa < 1000)] = "high_cohesion"  # Clay-silt mix
    cohesion_level[ssa >= 1000] = "very_high_cohesion"  # Pure clay
    
    if np.isscalar(specific_surface_area):
        return cohesion_level.item()
    else:
        return cohesion_level