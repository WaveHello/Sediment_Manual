"""
Environmental Effects on Fine-Grained Sediment Transport (Chapter 4)

This module implements environmental factors affecting cohesive sediment behavior
including temperature, salinity, organic content, and pH effects. Based on 
various relationships and observations from ASCE Manual 110, Chapter 4.
"""

import numpy as np
from typing import Union, Optional, Tuple, Dict, Any
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import WATER_TEMPERATURE_20C, WATER_VISCOSITY_20C

def temperature_effect_on_viscosity(
    temperature: Union[float, np.ndarray],
    reference_temp: float = 20.0,
    reference_viscosity: float = WATER_VISCOSITY_20C
) -> Union[float, np.ndarray]:
    """
    Calculate temperature effect on water viscosity.
    
    Uses empirical relationship for water viscosity with temperature.
    
    Parameters
    ----------
    temperature : float or array-like
        Water temperature [°C]
    reference_temp : float, optional
        Reference temperature [°C], default 20.0
    reference_viscosity : float, optional
        Viscosity at reference temperature [Pa·s], default 1.002e-3
        
    Returns
    -------
    float or ndarray
        Water viscosity at given temperature [Pa·s]
        
    References
    ----------
    Standard water properties relationships
    """
    validate_range(temperature, "temperature", 0.0, 50.0)
    
    T = np.atleast_1d(temperature)
    T_ref = reference_temp
    
    # Empirical relationship for water viscosity
    # μ(T) = μ_ref * exp(A * (1/(T+273) - 1/(T_ref+273)))
    A = 1800  # Empirical constant
    
    viscosity = reference_viscosity * np.exp(A * (1/(T + 273.15) - 1/(T_ref + 273.15)))
    
    return viscosity.item() if np.isscalar(temperature) else viscosity


def temperature_effect_on_flocculation(
    temperature: Union[float, np.ndarray],
    reference_settling_velocity: Union[float, np.ndarray],
    reference_temp: float = 15.0
) -> Union[float, np.ndarray]:
    """
    Calculate temperature correction for settling velocity.
    
    Equations 4-22a,b: Temperature-dependent flocculation
    Based on Jiang (1999) and Lau (1994) findings
    
    Parameters
    ----------
    temperature : float or array-like
        Water temperature [°C]
    reference_settling_velocity : float or array-like
        Settling velocity at reference temperature [m/s]
    reference_temp : float, optional
        Reference temperature [°C], default 15.0
        
    Returns
    -------
    float or ndarray
        Temperature-corrected settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equations 4-22a,b
    Jiang (1999), Lau (1994)
    """
    validate_range(temperature, "temperature", 0.0, 50.0)
    validate_positive(reference_settling_velocity, "reference_settling_velocity")
    
    T = np.atleast_1d(temperature)
    ws_ref = np.atleast_1d(reference_settling_velocity)
    
    # Temperature correction factor
    T_prime = T / reference_temp
    phi = 1.776 * (1 + 0.875 * T_prime)
    
    # Higher temperature reduces floc size and settling velocity
    ws_corrected = ws_ref / phi
    
    return ws_corrected.item() if (np.isscalar(temperature) and 
                                  np.isscalar(reference_settling_velocity)) else ws_corrected


def salinity_flocculation_threshold(
    clay_type: str
) -> float:
    """
    Get critical salinity for flocculation of different clay minerals.
    
    From Table 4-2: Critical salinities for clay flocculation
    
    Parameters
    ----------
    clay_type : str
        Type of clay mineral ('kaolinite', 'illite', 'smectite', 'montmorillonite')
        
    Returns
    -------
    float
        Critical salinity [ppt]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-2
    """
    critical_salinities = {
        'kaolinite': 0.6,
        'illite': 1.1,
        'chlorite': None,  # Not reported
        'smectite': 2.4,
        'montmorillonite': 2.4
    }
    
    clay_type = clay_type.lower()
    if clay_type not in critical_salinities:
        available_types = [k for k in critical_salinities.keys() if critical_salinities[k] is not None]
        raise ValueError(f"Unknown clay type: {clay_type}. Available: {available_types}")
    
    salinity = critical_salinities[clay_type]
    if salinity is None:
        raise ValueError(f"Critical salinity not reported for {clay_type}")
    
    return salinity


def salinity_effect_significance(
    salinity: Union[float, np.ndarray],
    critical_threshold: float = 10.0
) -> Union[bool, np.ndarray]:
    """
    Determine if salinity significantly affects floc structure.
    
    Based on observation that above ~10 ppt, salinity effects are minor.
    
    Parameters
    ----------
    salinity : float or array-like
        Water salinity [ppt]
    critical_threshold : float, optional
        Threshold salinity [ppt], default 10.0
        
    Returns
    -------
    bool or ndarray
        True if salinity effects are significant
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Table 4-3
    Krone (1962, 1986)
    """
    validate_positive(salinity, "salinity")
    
    s = np.atleast_1d(salinity)
    significant = s < critical_threshold
    
    return significant.item() if np.isscalar(salinity) else significant


def organic_content_effect_on_settling(
    organic_content: Union[float, np.ndarray],
    base_velocity_coeff: float = 0.2
) -> Union[float, np.ndarray]:
    """
    Calculate effect of organic content on settling velocity coefficient.
    
    Based on Ganju (2001) findings for Florida sediments.
    aw = aw0 + a1*OC + a2*OC² + a3*OC³ + a4*OC⁴
    
    Parameters
    ----------
    organic_content : float or array-like
        Organic content as fraction (0-1)
    base_velocity_coeff : float, optional
        Base velocity coefficient aw0, default 0.2
        
    Returns
    -------
    float or ndarray
        Modified velocity coefficient
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.5.3
    Ganju (2001)
    """
    validate_range(organic_content, "organic_content", 0.0, 1.0)
    
    OC = np.atleast_1d(organic_content)
    
    # Coefficients from Ganju (2001)
    a0 = base_velocity_coeff
    a1 = -6.67e-4
    a2 = 1.7e-4
    a3 = -7.1e-6
    a4 = 1.3e-7
    
    aw = a0 + a1*OC + a2*OC**2 + a3*OC**3 + a4*OC**4
    
    return aw.item() if np.isscalar(organic_content) else aw


def organic_content_effect_on_aggregation(
    organic_content: Union[float, np.ndarray]
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate effects of organic content on aggregation processes.
    
    Organic matter can enhance or reduce aggregation depending on type.
    
    Parameters
    ----------
    organic_content : float or array-like
        Organic content as fraction (0-1)
        
    Returns
    -------
    Dict[str, Union[float, np.ndarray]]
        Dictionary of effects on aggregation properties
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    Dennett et al. (1998), Wells and Goldberg (1993)
    """
    validate_range(organic_content, "organic_content", 0.0, 1.0)
    
    OC = np.atleast_1d(organic_content)
    
    effects = {}
    
    # Biopolymers can enhance flocculation through bridging
    effects['bridge_formation_factor'] = 1.0 + 2.0 * OC
    
    # But high organic content can reduce density and strength
    effects['density_reduction_factor'] = 1.0 - 0.5 * OC
    effects['strength_reduction_factor'] = 1.0 - 0.7 * OC
    
    # Organic coatings can affect surface properties
    effects['surface_modification'] = OC
    
    return effects


def mixed_organic_inorganic_settling(
    organic_fraction: Union[float, np.ndarray],
    inorganic_settling_velocity: Union[float, np.ndarray],
    pure_organic_settling_velocity: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate settling velocity for organic-inorganic mixtures.
    
    Based on Kranck and Milligan (1980) observation that 50% organic/50% inorganic
    mixture settled an order of magnitude faster than pure inorganic.
    
    Parameters
    ----------
    organic_fraction : float or array-like
        Organic fraction by weight (0-1)
    inorganic_settling_velocity : float or array-like
        Settling velocity of pure inorganic sediment [m/s]
    pure_organic_settling_velocity : float or array-like
        Settling velocity of pure organic material [m/s]
        
    Returns
    -------
    float or ndarray
        Mixed material settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    Kranck and Milligan (1980)
    """
    validate_range(organic_fraction, "organic_fraction", 0.0, 1.0)
    validate_positive(inorganic_settling_velocity, "inorganic_settling_velocity")
    validate_positive(pure_organic_settling_velocity, "pure_organic_settling_velocity")
    
    f_org = np.atleast_1d(organic_fraction)
    ws_inorg = np.atleast_1d(inorganic_settling_velocity)
    ws_org = np.atleast_1d(pure_organic_settling_velocity)
    
    # Enhancement factor peaks around 50% organic content
    enhancement = 1.0 + 9.0 * f_org * (1 - f_org)  # Peaks at f_org = 0.5
    
    # Weighted mixture with enhancement
    ws_mixed = ((1 - f_org) * ws_inorg + f_org * ws_org) * enhancement
    
    return ws_mixed.item() if (np.isscalar(organic_fraction) and
                              np.isscalar(inorganic_settling_velocity)) else ws_mixed


def ph_effect_on_aggregation(
    ph: Union[float, np.ndarray],
    optimal_ph: float = 7.0
) -> Union[float, np.ndarray]:
    """
    Calculate pH effects on particle aggregation.
    
    Slightly acidic conditions can enhance aggregation.
    
    Parameters
    ----------
    ph : float or array-like
        Water pH
    optimal_ph : float, optional
        pH for optimal aggregation, default 7.0
        
    Returns
    -------
    float or ndarray
        Aggregation enhancement factor (1.0 = neutral effect)
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.2.2
    Tsai and Hu (1997)
    """
    validate_range(ph, "ph", 4.0, 12.0)
    
    pH = np.atleast_1d(ph)
    
    # Enhancement factor (empirical relationship)
    # Slightly acidic pH enhances aggregation
    deviation = np.abs(pH - optimal_ph)
    enhancement = 1.0 + 0.2 * np.exp(-deviation**2 / 2.0)
    
    return enhancement.item() if np.isscalar(ph) else enhancement


def biogenic_effects_on_flocculation(
    biological_activity_level: Union[float, np.ndarray],
    chlorophyll_a: Optional[Union[float, np.ndarray]] = None,
    carbohydrate_content: Optional[Union[float, np.ndarray]] = None
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate biogenic effects on flocculation and erodibility.
    
    Biological activity affects sediment properties through:
    - Mucous filament production
    - Biopolymer secretion
    - Particle filtering by zooplankton
    
    Parameters
    ----------
    biological_activity_level : float or array-like
        Relative biological activity (0-1)
    chlorophyll_a : float or array-like, optional
        Chlorophyll-a concentration [μg/L]
    carbohydrate_content : float or array-like, optional
        Colloidal carbohydrate content [mg/L]
        
    Returns
    -------
    Dict[str, Union[float, np.ndarray]]
        Dictionary of biological effects
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.1
    Montague et al. (1993), Amos et al. (1998)
    """
    validate_range(biological_activity_level, "biological_activity_level", 0.0, 1.0)
    
    bio_level = np.atleast_1d(biological_activity_level)
    
    effects = {}
    
    # Mucous filament reinforcement
    effects['floc_reinforcement'] = 1.0 + 2.0 * bio_level
    
    # Enhanced aggregation through filtering
    effects['aggregation_enhancement'] = 1.0 + bio_level
    
    # Effects on critical shear stress for erosion
    if chlorophyll_a is not None:
        chl_a = np.atleast_1d(chlorophyll_a)
        # Higher chlorophyll-a typically increases erosion resistance
        effects['erosion_resistance'] = 1.0 + 0.1 * np.log(1 + chl_a)
    
    if carbohydrate_content is not None:
        carb = np.atleast_1d(carbohydrate_content)
        # Carbohydrates can increase cohesion
        effects['cohesion_enhancement'] = 1.0 + 0.2 * carb / (10 + carb)
    
    return effects


def seasonal_variation_effects(
    day_of_year: Union[int, np.ndarray],
    latitude: float = 40.0
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate seasonal effects on sediment transport processes.
    
    Parameters
    ----------
    day_of_year : int or array-like
        Day of year (1-365)
    latitude : float, optional
        Latitude in degrees, default 40.0
        
    Returns
    -------
    Dict[str, Union[float, np.ndarray]]
        Dictionary of seasonal factors
        
    References
    ----------
    Based on general seasonal patterns in aquatic systems
    """
    validate_range(day_of_year, "day_of_year", 1, 365)
    validate_range(latitude, "latitude", -90.0, 90.0)
    
    day = np.atleast_1d(day_of_year)
    
    # Temperature variation (simplified sinusoidal)
    temp_variation = np.sin(2 * np.pi * (day - 80) / 365)  # Peak around day 170 (June)
    
    # Biological activity (follows temperature with lag)
    bio_activity = np.sin(2 * np.pi * (day - 120) / 365)  # Peak around day 210 (July)
    bio_activity = np.maximum(0, bio_activity)  # Only positive values
    
    effects = {
        'temperature_factor': 1.0 + 0.3 * temp_variation,
        'biological_activity': bio_activity,
        'settling_velocity_modifier': 1.0 - 0.2 * temp_variation,  # Faster settling in winter
        'flocculation_modifier': 1.0 + 0.1 * bio_activity  # Enhanced in summer
    }
    
    return effects


def ionic_strength_effects(
    ionic_strength: Union[float, np.ndarray],
    dominant_cation: str = 'sodium'
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Calculate effects of ionic strength on clay behavior.
    
    Parameters
    ----------
    ionic_strength : float or array-like
        Ionic strength [mol/L]
    dominant_cation : str, optional
        Dominant cation type ('sodium', 'calcium', 'potassium'), default 'sodium'
        
    Returns
    -------
    Dict[str, Union[float, np.ndarray]]
        Dictionary of ionic effects
        
    References
    ----------
    Based on clay chemistry principles
    """
    validate_positive(ionic_strength, "ionic_strength")
    
    I = np.atleast_1d(ionic_strength)
    
    # Double layer thickness decreases with ionic strength
    double_layer_thickness = 1.0 / np.sqrt(I + 1e-6)
    
    # Flocculation efficiency increases with ionic strength
    flocculation_efficiency = 1.0 - np.exp(-I / 0.1)
    
    effects = {
        'double_layer_thickness': double_layer_thickness,
        'flocculation_efficiency': flocculation_efficiency,
        'aggregation_rate': 1.0 + 2.0 * flocculation_efficiency
    }
    
    # Cation-specific effects
    if dominant_cation.lower() == 'calcium':
        effects['bridging_enhancement'] = 1.5  # Ca²⁺ enhances bridging
    elif dominant_cation.lower() == 'sodium':
        effects['bridging_enhancement'] = 1.0  # Na⁺ baseline
    elif dominant_cation.lower() == 'potassium':
        effects['bridging_enhancement'] = 1.2  # K⁺ moderate enhancement
    
    return effects


def consolidation_environmental_effects(
    temperature: float,
    salinity: float,
    organic_content: float
) -> Dict[str, float]:
    """
    Calculate environmental effects on consolidation processes.
    
    Parameters
    ----------
    temperature : float
        Water temperature [°C]
    salinity : float
        Water salinity [ppt]
    organic_content : float
        Organic content fraction (0-1)
        
    Returns
    -------
    Dict[str, float]
        Dictionary of consolidation modifiers
        
    References
    ----------
    Based on consolidation theory and observations
    """
    validate_range(temperature, "temperature", 0.0, 50.0)
    validate_positive(salinity, "salinity")
    validate_range(organic_content, "organic_content", 0.0, 1.0)
    
    effects = {}
    
    # Temperature affects consolidation rate through viscosity
    temp_factor = np.exp(0.02 * (temperature - 20))  # Higher T = faster consolidation
    effects['consolidation_rate_modifier'] = temp_factor
    
    # Salinity affects floc structure and consolidation
    if salinity > 10:
        effects['salinity_modifier'] = 1.0  # Minor effect above 10 ppt
    else:
        effects['salinity_modifier'] = 0.5 + 0.5 * salinity / 10  # Stronger effect below 10 ppt
    
    # Organic content reduces consolidation efficiency
    effects['organic_modifier'] = 1.0 - 0.8 * organic_content
    
    # Combined effect
    effects['overall_consolidation_modifier'] = (effects['consolidation_rate_modifier'] * 
                                               effects['salinity_modifier'] * 
                                               effects['organic_modifier'])
    
    return effects