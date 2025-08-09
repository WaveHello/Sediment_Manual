"""
Suspended load transport equations.

Implements equations for calculating suspended sediment transport
from Chapter 2, Section 2.9.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range
from ..utils.constants import VON_KARMAN_CONSTANT

def rouse_parameter(
    settling_velocity: Union[float, np.ndarray],
    shear_velocity: Union[float, np.ndarray],
    kappa: float = VON_KARMAN_CONSTANT
) -> Union[float, np.ndarray]:
    """
    Calculate Rouse parameter for suspended sediment.
    
    Equation 2.62: P = ws/(κu*)
    
    Where:
    - ws is settling velocity [m/s]
    - κ is von Kármán constant [-]
    - u* is shear velocity [m/s]
    
    Parameters
    ----------
    settling_velocity : float or array-like
        Particle settling velocity [m/s]
    shear_velocity : float or array-like
        Shear velocity [m/s]
    kappa : float, optional
        von Kármán constant [-], default 0.41
        
    Returns
    -------
    float or ndarray
        Rouse parameter [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.62
    """
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(shear_velocity, "shear_velocity")
    validate_positive(kappa, "kappa")
    
    return settling_velocity / (kappa * shear_velocity)


def concentration_profile(
    reference_concentration: Union[float, np.ndarray],
    height: Union[float, np.ndarray],
    reference_height: Union[float, np.ndarray],
    flow_depth: Union[float, np.ndarray],
    rouse_param: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate suspended sediment concentration using Rouse profile.
    
    Equation 2.61: C(z)/Ca = [(H-z)/z * a/(H-a)]^P
    
    Where:
    - C(z) is concentration at height z [kg/m³]
    - Ca is reference concentration [kg/m³]  
    - H is flow depth [m]
    - z is height above bed [m]
    - a is reference height [m]
    - P is Rouse parameter [-]
    
    Parameters
    ----------
    reference_concentration : float or array-like
        Reference concentration [kg/m³]
    height : float or array-like
        Height above bed [m]
    reference_height : float or array-like
        Reference height [m]
    flow_depth : float or array-like
        Flow depth [m]
    rouse_param : float or array-like
        Rouse parameter [-]
        
    Returns
    -------
    float or ndarray
        Concentration at specified height [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.61
    """
    validate_positive(reference_concentration, "reference_concentration")
    validate_positive(height, "height")
    validate_positive(reference_height, "reference_height")
    validate_positive(flow_depth, "flow_depth")
    validate_positive(rouse_param, "rouse_param")
    
    # Ensure height is within flow depth
    height = np.minimum(height, flow_depth * 0.999)  # Avoid singularity at surface
    
    ratio1 = (flow_depth - height) / height
    ratio2 = reference_height / (flow_depth - reference_height)
    
    return reference_concentration * (ratio1 * ratio2)**rouse_param


def suspended_load_discharge(
    velocity_profile: Union[float, np.ndarray],
    concentration_profile_func: callable,
    flow_depth: Union[float, np.ndarray],
    reference_height: Union[float, np.ndarray],
    num_points: int = 50
) -> Union[float, np.ndarray]:
    """
    Calculate suspended load discharge by integrating velocity and concentration profiles.
    
    Equation 2.60: qs = ∫[a to H] u(z)C(z) dz
    
    Where:
    - qs is suspended load discharge per unit width [kg/(m·s)]
    - u(z) is velocity at height z [m/s]
    - C(z) is concentration at height z [kg/m³]
    - integration from reference height a to flow depth H
    
    Parameters
    ----------
    velocity_profile : callable or float/array-like
        Velocity profile function u(z) or constant velocity
    concentration_profile_func : callable
        Function returning concentration as function of height
    flow_depth : float or array-like
        Flow depth [m]
    reference_height : float or array-like
        Reference height (bottom of integration) [m]
    num_points : int, optional
        Number of integration points, default 50
        
    Returns
    -------
    float or ndarray
        Suspended load discharge per unit width [kg/(m·s)]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.60
    """
    validate_positive(flow_depth, "flow_depth")
    validate_positive(reference_height, "reference_height")
    
    # Create integration points
    z = np.linspace(reference_height, flow_depth * 0.999, num_points)
    dz = z[1] - z[0]
    
    # Get velocity at each point
    if callable(velocity_profile):
        u = velocity_profile(z)
    else:
        u = velocity_profile  # Constant velocity
    
    # Get concentration at each point
    C = concentration_profile_func(z)
    
    # Integrate using trapezoidal rule
    integrand = u * C
    
    return np.trapz(integrand, z)


def logarithmic_velocity_profile(
    height: Union[float, np.ndarray],
    shear_velocity: Union[float, np.ndarray],
    roughness_height: Union[float, np.ndarray],
    kappa: float = VON_KARMAN_CONSTANT
) -> Union[float, np.ndarray]:
    """
    Calculate velocity using logarithmic profile.
    
    Equation 2.5: u(z) = (u*/κ) * ln(z/z0)
    
    Where:
    - u* is shear velocity [m/s]
    - κ is von Kármán constant [-]
    - z is height above bed [m]
    - z0 is roughness height [m]
    
    Parameters
    ----------
    height : float or array-like
        Height above bed [m]
    shear_velocity : float or array-like
        Shear velocity [m/s]
    roughness_height : float or array-like
        Roughness height [m]
    kappa : float, optional
        von Kármán constant [-], default 0.41
        
    Returns
    -------
    float or ndarray
        Velocity at specified height [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.5
    """
    validate_positive(height, "height")
    validate_positive(shear_velocity, "shear_velocity")
    validate_positive(roughness_height, "roughness_height")
    validate_positive(kappa, "kappa")
    
    # Ensure height is greater than roughness height
    height = np.maximum(height, roughness_height * 1.1)
    
    return (shear_velocity / kappa) * np.log(height / roughness_height)


def van_rijn_suspended_load(
    flow_velocity: Union[float, np.ndarray],
    flow_depth: Union[float, np.ndarray],
    particle_diameter: Union[float, np.ndarray],
    settling_velocity: Union[float, np.ndarray],
    shields_parameter: Union[float, np.ndarray],
    critical_shields_parameter: Union[float, np.ndarray] = 0.047
) -> Union[float, np.ndarray]:
    """
    Calculate suspended load using van Rijn (1984) method.
    
    This is a simplified implementation of the van Rijn approach.
    
    Parameters
    ----------
    flow_velocity : float or array-like
        Mean flow velocity [m/s]
    flow_depth : float or array-like
        Flow depth [m]
    particle_diameter : float or array-like
        Particle diameter [m]
    settling_velocity : float or array-like
        Settling velocity [m/s]
    shields_parameter : float or array-like
        Shields parameter [-]
    critical_shields_parameter : float or array-like, optional
        Critical Shields parameter [-], default 0.047
        
    Returns
    -------
    float or ndarray
        Suspended load discharge [kg/(m·s)]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.9.4
    van Rijn, L.C. (1984). Sediment transport, part II: suspended load transport.
    """
    validate_positive(flow_velocity, "flow_velocity")
    validate_positive(flow_depth, "flow_depth")
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(shields_parameter, "shields_parameter")
    
    # Transport parameter
    T = (shields_parameter - critical_shields_parameter) / critical_shields_parameter
    T = np.maximum(T, 0)  # No transport below threshold
    
    # Dimensionless particle size
    D_star = particle_diameter * (1650 * 9.81 / (1e-6)**2)**(1/3)  # Simplified
    
    # Reference concentration (simplified)
    Ca = 0.015 * (particle_diameter / flow_depth) * T**1.5 / D_star**0.3
    
    # Suspended load coefficient (simplified)
    Z = settling_velocity / (0.41 * 0.1 * flow_velocity)  # Simplified shear velocity
    
    # Integration factor
    if np.any(Z > 2.5):
        # For large Z, most sediment near bed
        integration_factor = flow_depth / (Z + 1)
    else:
        # For small Z, more uniform distribution
        integration_factor = flow_depth / (2 * Z + 1)
    
    return Ca * flow_velocity * integration_factor * 2650  # Convert to mass