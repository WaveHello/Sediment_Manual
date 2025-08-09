"""
Fluid mechanics equations for sediment transport.

Implements fundamental fluid mechanics relationships needed for 
sediment transport calculations from Chapter 2, Section 2.2.
"""

import numpy as np
from typing import Union, Optional
from ..utils.validators import validate_positive, validate_range

def reynolds_number(
    velocity: Union[float, np.ndarray], 
    length_scale: Union[float, np.ndarray],
    kinematic_viscosity: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate Reynolds number.
    
    Equation 2.1: Re = VL/ν
    
    Where:
    - V is characteristic velocity [m/s]
    - L is characteristic length [m]  
    - ν is kinematic viscosity [m²/s]
    
    Parameters
    ----------
    velocity : float or array-like
        Characteristic velocity [m/s]
    length_scale : float or array-like
        Characteristic length scale [m]
    kinematic_viscosity : float or array-like
        Kinematic viscosity [m²/s]
        
    Returns
    -------
    float or ndarray
        Reynolds number [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.1
    """
    validate_positive(velocity, "velocity")
    validate_positive(length_scale, "length_scale")  
    validate_positive(kinematic_viscosity, "kinematic_viscosity")
    
    return velocity * length_scale / kinematic_viscosity


def froude_number(
    velocity: Union[float, np.ndarray],
    depth: Union[float, np.ndarray], 
    gravity: float = 9.81
) -> Union[float, np.ndarray]:
    """
    Calculate Froude number.
    
    Equation 2.2: Fr = V/√(gH)
    
    Where:
    - V is flow velocity [m/s]
    - g is gravitational acceleration [m/s²]
    - H is flow depth [m]
    
    Parameters
    ----------
    velocity : float or array-like
        Flow velocity [m/s]
    depth : float or array-like
        Flow depth [m]
    gravity : float, optional
        Gravitational acceleration [m/s²], default 9.81
        
    Returns
    -------
    float or ndarray
        Froude number [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.2
    """
    validate_positive(velocity, "velocity")
    validate_positive(depth, "depth")
    validate_positive(gravity, "gravity")
    
    return velocity / np.sqrt(gravity * depth)


def shear_velocity(
    boundary_shear_stress: Union[float, np.ndarray],
    fluid_density: Union[float, np.ndarray] = 1000.0
) -> Union[float, np.ndarray]:
    """
    Calculate shear velocity from boundary shear stress.
    
    Equation 2.8: u* = √(τ₀/ρ)
    
    Where:
    - τ₀ is boundary shear stress [N/m²]
    - ρ is fluid density [kg/m³]
    
    Parameters
    ----------
    boundary_shear_stress : float or array-like
        Boundary shear stress [N/m²]
    fluid_density : float or array-like, optional
        Fluid density [kg/m³], default 1000.0 (water)
        
    Returns
    -------
    float or ndarray
        Shear velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.8
    """
    validate_positive(boundary_shear_stress, "boundary_shear_stress")
    validate_positive(fluid_density, "fluid_density")
    
    return np.sqrt(boundary_shear_stress / fluid_density)


def boundary_shear_stress(
    fluid_density: Union[float, np.ndarray],
    gravity: float,
    hydraulic_radius: Union[float, np.ndarray],
    slope: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate boundary shear stress.
    
    Equation 2.7: τ₀ = ρgRS
    
    Where:
    - ρ is fluid density [kg/m³]
    - g is gravitational acceleration [m/s²]
    - R is hydraulic radius [m]
    - S is energy slope [-]
    
    Parameters
    ----------
    fluid_density : float or array-like
        Fluid density [kg/m³]
    gravity : float
        Gravitational acceleration [m/s²]
    hydraulic_radius : float or array-like
        Hydraulic radius [m]
    slope : float or array-like
        Energy slope [-]
        
    Returns
    -------
    float or ndarray
        Boundary shear stress [N/m²]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.7
    """
    validate_positive(fluid_density, "fluid_density")
    validate_positive(gravity, "gravity")
    validate_positive(hydraulic_radius, "hydraulic_radius")
    validate_positive(slope, "slope")
    
    return fluid_density * gravity * hydraulic_radius * slope