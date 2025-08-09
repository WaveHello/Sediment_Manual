"""
Sediment properties and particle settling equations.

Implements equations for sediment characterization and settling velocity
calculations from Chapter 2, Section 2.3.
"""

import numpy as np
from typing import Union, Optional
from ..utils.validators import validate_positive, validate_range

def settling_velocity_stokes(
    particle_diameter: Union[float, np.ndarray],
    sediment_density: Union[float, np.ndarray] = 2650.0,
    fluid_density: Union[float, np.ndarray] = 1000.0,
    kinematic_viscosity: Union[float, np.ndarray] = 1.0e-6,
    gravity: float = 9.81
) -> Union[float, np.ndarray]:
    """
    Calculate particle settling velocity using Stokes' law.
    
    Valid for particles with Re_p < 0.5 (typically D < 0.1 mm)
    
    Equation 2.17: ws = (ρs - ρf)gD²/(18μ)
    
    Where:
    - ρs is sediment density [kg/m³]
    - ρf is fluid density [kg/m³]
    - g is gravitational acceleration [m/s²]
    - D is particle diameter [m]
    - μ is dynamic viscosity [kg/(m·s)]
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    sediment_density : float or array-like, optional
        Sediment density [kg/m³], default 2650.0 (quartz)
    fluid_density : float or array-like, optional
        Fluid density [kg/m³], default 1000.0 (water)
    kinematic_viscosity : float or array-like, optional
        Kinematic viscosity [m²/s], default 1.0e-6 (water at 20°C)
    gravity : float, optional
        Gravitational acceleration [m/s²], default 9.81
        
    Returns
    -------
    float or ndarray
        Settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.17
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(sediment_density, "sediment_density")
    validate_positive(fluid_density, "fluid_density")
    validate_positive(kinematic_viscosity, "kinematic_viscosity")
    
    dynamic_viscosity = kinematic_viscosity * fluid_density
    
    return ((sediment_density - fluid_density) * gravity * particle_diameter**2) / (18 * dynamic_viscosity)


def settling_velocity_dietrich(
    particle_diameter: Union[float, np.ndarray],
    sediment_density: Union[float, np.ndarray] = 2650.0,
    fluid_density: Union[float, np.ndarray] = 1000.0,
    kinematic_viscosity: Union[float, np.ndarray] = 1.0e-6,
    gravity: float = 9.81
) -> Union[float, np.ndarray]:
    """
    Calculate particle settling velocity using Dietrich (1982) method.
    
    Valid for natural sediment particles across full size range.
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    sediment_density : float or array-like, optional
        Sediment density [kg/m³], default 2650.0
    fluid_density : float or array-like, optional
        Fluid density [kg/m³], default 1000.0
    kinematic_viscosity : float or array-like, optional
        Kinematic viscosity [m²/s], default 1.0e-6
    gravity : float, optional
        Gravitational acceleration [m/s²], default 9.81
        
    Returns
    -------
    float or ndarray
        Settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Section 2.3.4
    Dietrich, W.E. (1982). Settling velocity of natural particles. 
    Water Resources Research, 18(6), 1615-1626.
    """
    validate_positive(particle_diameter, "particle_diameter")
    
    # Calculate dimensionless diameter
    D_star = dimensionless_diameter(
        particle_diameter, sediment_density, fluid_density, kinematic_viscosity, gravity
    )
    
    # Dietrich coefficients for natural particles
    a1 = -3.76715
    a2 = 1.92944  
    a3 = -0.09815
    a4 = -0.00575
    a5 = 0.00056
    
    # Calculate dimensionless settling velocity
    log_D_star = np.log10(D_star)
    
    log_w_star = (a1 + a2 * log_D_star + a3 * log_D_star**2 + 
                  a4 * log_D_star**3 + a5 * log_D_star**4)
    
    w_star = 10**log_w_star
    
    # Convert to dimensional settling velocity
    settling_velocity = w_star * np.sqrt(
        (sediment_density / fluid_density - 1) * gravity * kinematic_viscosity
    )
    
    return settling_velocity


def particle_reynolds_number(
    settling_velocity: Union[float, np.ndarray],
    particle_diameter: Union[float, np.ndarray], 
    kinematic_viscosity: Union[float, np.ndarray] = 1.0e-6
) -> Union[float, np.ndarray]:
    """
    Calculate particle Reynolds number.
    
    Equation 2.18: Re_p = wsD/ν
    
    Where:
    - ws is settling velocity [m/s]
    - D is particle diameter [m]
    - ν is kinematic viscosity [m²/s]
    
    Parameters
    ----------
    settling_velocity : float or array-like
        Particle settling velocity [m/s]
    particle_diameter : float or array-like
        Particle diameter [m]
    kinematic_viscosity : float or array-like, optional
        Kinematic viscosity [m²/s], default 1.0e-6
        
    Returns
    -------
    float or ndarray
        Particle Reynolds number [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.18
    """
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(kinematic_viscosity, "kinematic_viscosity")
    
    return settling_velocity * particle_diameter / kinematic_viscosity


def dimensionless_diameter(
    particle_diameter: Union[float, np.ndarray],
    sediment_density: Union[float, np.ndarray] = 2650.0,
    fluid_density: Union[float, np.ndarray] = 1000.0,
    kinematic_viscosity: Union[float, np.ndarray] = 1.0e-6,
    gravity: float = 9.81
) -> Union[float, np.ndarray]:
    """
    Calculate dimensionless particle diameter.
    
    Equation 2.19: D* = D[(ρs/ρf - 1)g/ν²]^(1/3)
    
    Where:
    - D is particle diameter [m]
    - ρs is sediment density [kg/m³]
    - ρf is fluid density [kg/m³]
    - g is gravitational acceleration [m/s²]
    - ν is kinematic viscosity [m²/s]
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    sediment_density : float or array-like, optional
        Sediment density [kg/m³], default 2650.0
    fluid_density : float or array-like, optional
        Fluid density [kg/m³], default 1000.0
    kinematic_viscosity : float or array-like, optional
        Kinematic viscosity [m²/s], default 1.0e-6
    gravity : float, optional
        Gravitational acceleration [m/s²], default 9.81
        
    Returns
    -------
    float or ndarray
        Dimensionless diameter [-]
        
    References
    ----------
    ASCE Manual 110, Chapter 2, Equation 2.19
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(sediment_density, "sediment_density")
    validate_positive(fluid_density, "fluid_density")
    validate_positive(kinematic_viscosity, "kinematic_viscosity")
    
    relative_density = sediment_density / fluid_density - 1
    
    return particle_diameter * (relative_density * gravity / kinematic_viscosity**2)**(1/3)