"""
Bed load transport equations.

Implements equations for calculating bed load sediment transport
from Chapter 2, Section 2.6.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range
from ..utils.constants import GRAVITY

def meyer_peter_muller(
    shields_parameter: Union[float, np.ndarray],
    critical_shields_parameter: Union[float, np.ndarray] = 0.047,
    hiding_factor: Union[float, np.ndarray] = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate dimensionless bed load transport rate using Meyer-Peter Müller equation.
    
    q*b = 8(θ - θcr)^1.5
    
    Where:
    - q*b is dimensionless bed load transport rate
    - θ is Shields parameter
    - θcr is critical Shields parameter
    
    Parameters
    ----------
    shields_parameter : float or array-like
        Shields parameter [-]
    critical_shields_parameter : float or array-like, optional
        Critical Shields parameter [-], default 0.047
    hiding_factor : float or array-like, optional
        Hiding factor for mixed-size sediments [-], default 1.0
        
    Returns
    -------
    float or ndarray
        Dimensionless bed load transport rate [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.35
    Meyer-Peter, E., & Müller, R. (1948). Formulas for bed-load transport.
    """
    validate_positive(shields_parameter, "shields_parameter")
    validate_positive(critical_shields_parameter, "critical_shields_parameter")
    validate_positive(hiding_factor, "hiding_factor")
    
    excess_shear = shields_parameter - critical_shields_parameter * hiding_factor
    
    # Only transport when threshold is exceeded
    excess_shear = np.maximum(excess_shear, 0)
    
    return 8 * excess_shear**1.5


def einstein_brown(
    shields_parameter: Union[float, np.ndarray],
    critical_shields_parameter: Union[float, np.ndarray] = 0.047
) -> Union[float, np.ndarray]:
    """
    Calculate dimensionless bed load transport rate using Einstein-Brown equation.
    
    q*b = 40 * θ^3 / (1 + (θcr/θ)^3)
    
    Parameters
    ----------
    shields_parameter : float or array-like
        Shields parameter [-]
    critical_shields_parameter : float or array-like, optional
        Critical Shields parameter [-], default 0.047
        
    Returns
    -------
    float or ndarray
        Dimensionless bed load transport rate [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.36
    Einstein, H. A., & Brown, C. B. (1950). Fluid resistance of composite roughness.
    """
    validate_positive(shields_parameter, "shields_parameter")
    validate_positive(critical_shields_parameter, "critical_shields_parameter")
    
    denominator = 1 + (critical_shields_parameter / shields_parameter)**3
    
    return 40 * shields_parameter**3 / denominator


def bed_load_discharge_volumetric(
    dimensionless_transport_rate: Union[float, np.ndarray],
    particle_diameter: Union[float, np.ndarray],
    sediment_density: Union[float, np.ndarray] = 2650.0,
    fluid_density: Union[float, np.ndarray] = 1000.0,
    gravity: float = GRAVITY
) -> Union[float, np.ndarray]:
    """
    Calculate volumetric bed load discharge per unit width.
    
    qb = q*b * sqrt[(ρs/ρf - 1)gD³]
    
    Where:
    - qb is volumetric bed load discharge per unit width [m²/s]
    - q*b is dimensionless transport rate [-]
    - ρs is sediment density [kg/m³]
    - ρf is fluid density [kg/m³]
    - g is gravitational acceleration [m/s²]
    - D is particle diameter [m]
    
    Parameters
    ----------
    dimensionless_transport_rate : float or array-like
        Dimensionless bed load transport rate [-]
    particle_diameter : float or array-like
        Particle diameter [m]
    sediment_density : float or array-like, optional
        Sediment density [kg/m³], default 2650.0
    fluid_density : float or array-like, optional
        Fluid density [kg/m³], default 1000.0
    gravity : float, optional
        Gravitational acceleration [m/s²], default 9.81
        
    Returns
    -------
    float or ndarray
        Volumetric bed load discharge per unit width [m²/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.34
    """
    validate_positive(dimensionless_transport_rate, "dimensionless_transport_rate")
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(sediment_density, "sediment_density") 
    validate_positive(fluid_density, "fluid_density")
    
    relative_density = sediment_density / fluid_density - 1
    scaling_factor = np.sqrt(relative_density * gravity * particle_diameter**3)
    
    return dimensionless_transport_rate * scaling_factor


def engelund_hansen(
    shields_parameter: Union[float, np.ndarray],
    friction_factor: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate total sediment transport using Engelund-Hansen equation.
    
    qt* = 0.05 * θ^2.5 / f
    
    Where:
    - qt* is dimensionless total transport rate [-]
    - θ is Shields parameter [-]
    - f is friction factor [-]
    
    Parameters
    ----------
    shields_parameter : float or array-like
        Shields parameter [-]
    friction_factor : float or array-like
        Friction factor [-]
        
    Returns
    -------
    float or ndarray
        Dimensionless total transport rate [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.37
    Engelund, F., & Hansen, E. (1967). A monograph on sediment transport in alluvial streams.
    """
    validate_positive(shields_parameter, "shields_parameter")
    validate_positive(friction_factor, "friction_factor")
    
    return 0.05 * shields_parameter**2.5 / friction_factor


def bed_load_particle_velocity(
    dimensionless_transport_rate: Union[float, np.ndarray],
    porosity: Union[float, np.ndarray] = 0.4,
    reference_concentration: Union[float, np.ndarray] = 0.65
) -> Union[float, np.ndarray]:
    """
    Estimate bed load particle velocity from transport rate.
    
    Vp ≈ qb / (ca * δ)
    
    Where:
    - Vp is particle velocity [m/s]
    - qb is bed load discharge [m²/s]
    - ca is reference concentration [-]
    - δ is bed load layer thickness [m]
    
    Parameters
    ----------
    dimensionless_transport_rate : float or array-like
        Dimensionless transport rate [-]
    porosity : float or array-like, optional
        Bed porosity [-], default 0.4
    reference_concentration : float or array-like, optional
        Reference concentration [-], default 0.65
        
    Returns
    -------
    float or ndarray
        Estimated particle velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.6.3
    """
    validate_range(porosity, "porosity", 0.0, 1.0)
    validate_range(reference_concentration, "reference_concentration", 0.0, 1.0)
    
    # Simplified estimate assuming bed load layer thickness ~ 2D
    # This is a rough approximation for illustration
    effective_concentration = reference_concentration * (1 - porosity)
    
    return dimensionless_transport_rate / (2 * effective_concentration)