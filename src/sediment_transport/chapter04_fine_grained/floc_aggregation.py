"""
Floc Aggregation and Fractal Models for Fine-Grained Sediment (Chapter 4)

This module implements floc aggregation processes, fractal models, and collision
mechanisms for cohesive sediment transport. Based on equations and models from
ASCE Manual 110, Chapter 4, Sections 4.4 and 4.5.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY, WATER_VISCOSITY_20C

# Physical constants
BOLTZMANN_CONSTANT = 1.38e-23  # J/K
WATER_TEMPERATURE_20C = 293.15  # K

def fractal_floc_density(
    floc_diameter: Union[float, np.ndarray],
    fractal_dimension: float = 2.1,
    proportionality_constant: float = 0.0008
) -> Union[float, np.ndarray]:
    """
    Calculate floc density using fractal relationship.
    
    Equation 4-9: ρf = σf * df^(nf-3)
    
    For three-dimensional space: 1 < nf < 3
    - nf ≈ 1.78 for Brownian motion aggregation
    - nf ≈ 1.9-2.1 for shear-induced aggregation
    - nf ≈ 2.1-2.6 for second order aggregates
    
    Parameters
    ----------
    floc_diameter : float or array-like
        Floc diameter [m]
    fractal_dimension : float, optional
        Fractal dimension nf, default 2.1
    proportionality_constant : float, optional
        Proportionality constant σf, default 0.0008
        
    Returns
    -------
    float or ndarray
        Floc density [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-9
    Meakin (1988), Kranenburg (1994a), Winterwerp (1998)
    """
    validate_positive(floc_diameter, "floc_diameter")
    validate_range(fractal_dimension, "fractal_dimension", 1.0, 3.0)
    validate_positive(proportionality_constant, "proportionality_constant")
    
    df = np.atleast_1d(floc_diameter)
    
    # Convert to cm and g/cm³ if using typical literature values
    df_cm = df * 100  # m to cm
    rho_f_g_cm3 = proportionality_constant * df_cm**(fractal_dimension - 3)
    
    # Convert back to kg/m³
    rho_f = rho_f_g_cm3 * 1000
    
    return rho_f.item() if np.isscalar(floc_diameter) else rho_f


def floc_strength(
    floc_density: Union[float, np.ndarray],
    water_density: float = WATER_DENSITY,
    strength_coeff: float = 1.5243e-7,
    strength_exponent: float = 3.0
) -> Union[float, np.ndarray]:
    """
    Calculate floc strength based on density.
    
    Equation 4-10: τf = αc * Δρf^βc
    
    Where Δρf = ρf - ρw (excess density)
    
    Parameters
    ----------
    floc_density : float or array-like
        Floc density [kg/m³]
    water_density : float, optional
        Water density [kg/m³], default 1000
    strength_coeff : float, optional
        Strength coefficient αc, default 1.5243e-7
    strength_exponent : float, optional
        Strength exponent βc, default 3.0
        
    Returns
    -------
    float or ndarray
        Floc strength [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-10
    Partheniades (1993), Krone (1963)
    """
    validate_positive(floc_density, "floc_density")
    validate_positive(water_density, "water_density")
    validate_positive(strength_coeff, "strength_coeff")
    
    delta_rho = floc_density - water_density
    # Ensure positive excess density
    delta_rho = np.maximum(delta_rho, 1.0)
    
    return strength_coeff * delta_rho**strength_exponent


def collision_frequency_brownian(
    diameter_i: Union[float, np.ndarray],
    diameter_j: Union[float, np.ndarray],
    temperature: float = WATER_TEMPERATURE_20C,
    viscosity: float = WATER_VISCOSITY_20C,
    collision_factor: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate Brownian motion collision frequency.
    
    Equation 4-13: βc = (2κT/3μ) * Fc * (di + dj)²
    
    Parameters
    ----------
    diameter_i : float or array-like
        Diameter of particle type i [m]
    diameter_j : float or array-like
        Diameter of particle type j [m]
    temperature : float, optional
        Absolute temperature [K], default 293.15 K (20°C)
    viscosity : float, optional
        Dynamic viscosity [Pa·s], default 1.002e-3
    collision_factor : float, optional
        Collision diameter correction factor Fc, default 1.0
        
    Returns
    -------
    float or ndarray
        Brownian collision frequency [m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-13
    McAnally (1999), Burban et al. (1989)
    """
    validate_positive(diameter_i, "diameter_i")
    validate_positive(diameter_j, "diameter_j")
    validate_positive(temperature, "temperature")
    validate_positive(viscosity, "viscosity")
    validate_range(collision_factor, "collision_factor", 0.0, 1.0)
    
    di = np.atleast_1d(diameter_i)
    dj = np.atleast_1d(diameter_j)
    
    coefficient = 2 * BOLTZMANN_CONSTANT * temperature / (3 * viscosity)
    beta_brownian = coefficient * collision_factor * (di + dj)**2
    
    return beta_brownian.item() if np.isscalar(diameter_i) and np.isscalar(diameter_j) else beta_brownian


def collision_frequency_shear(
    diameter_i: Union[float, np.ndarray],
    diameter_j: Union[float, np.ndarray],
    shear_rate: Union[float, np.ndarray],
    collision_factor: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate fluid shear collision frequency.
    
    Equation 4-13: βc = (4/15) * π * Gs * Fc * (di + dj)³
    
    Where Gs is the flow shear rate measure
    
    Parameters
    ----------
    diameter_i : float or array-like
        Diameter of particle type i [m]
    diameter_j : float or array-like
        Diameter of particle type j [m]
    shear_rate : float or array-like
        Flow shear rate Gs [s⁻¹]
    collision_factor : float, optional
        Collision diameter correction factor Fc, default 1.0
        
    Returns
    -------
    float or ndarray
        Shear collision frequency [m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-13
    Saffman and Turner (1956)
    """
    validate_positive(diameter_i, "diameter_i")
    validate_positive(diameter_j, "diameter_j")
    validate_positive(shear_rate, "shear_rate")
    validate_range(collision_factor, "collision_factor", 0.0, 1.0)
    
    di = np.atleast_1d(diameter_i)
    dj = np.atleast_1d(diameter_j)
    G = np.atleast_1d(shear_rate)
    
    coefficient = 4 * np.pi / 15
    beta_shear = coefficient * G * collision_factor * (di + dj)**3
    
    return beta_shear.item() if (np.isscalar(diameter_i) and 
                               np.isscalar(diameter_j) and 
                               np.isscalar(shear_rate)) else beta_shear


def collision_frequency_settling(
    diameter_i: Union[float, np.ndarray],
    diameter_j: Union[float, np.ndarray],
    settling_velocity_i: Union[float, np.ndarray],
    settling_velocity_j: Union[float, np.ndarray],
    collision_factor: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate differential settling collision frequency.
    
    Equation 4-13: βc = (π/4) * Fc * (di + dj)² * |wsi - wsj|
    
    Parameters
    ----------
    diameter_i : float or array-like
        Diameter of particle type i [m]
    diameter_j : float or array-like
        Diameter of particle type j [m]
    settling_velocity_i : float or array-like
        Settling velocity of particle type i [m/s]
    settling_velocity_j : float or array-like
        Settling velocity of particle type j [m/s]
    collision_factor : float, optional
        Collision diameter correction factor Fc, default 1.0
        
    Returns
    -------
    float or ndarray
        Differential settling collision frequency [m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-13
    Delichatsios and Probestein (1975)
    """
    validate_positive(diameter_i, "diameter_i")
    validate_positive(diameter_j, "diameter_j")
    validate_positive(settling_velocity_i, "settling_velocity_i")
    validate_positive(settling_velocity_j, "settling_velocity_j")
    validate_range(collision_factor, "collision_factor", 0.0, 1.0)
    
    di = np.atleast_1d(diameter_i)
    dj = np.atleast_1d(diameter_j)
    wsi = np.atleast_1d(settling_velocity_i)
    wsj = np.atleast_1d(settling_velocity_j)
    
    coefficient = np.pi / 4
    beta_settling = coefficient * collision_factor * (di + dj)**2 * np.abs(wsi - wsj)
    
    return beta_settling.item() if (np.isscalar(diameter_i) and 
                                   np.isscalar(diameter_j) and
                                   np.isscalar(settling_velocity_i) and
                                   np.isscalar(settling_velocity_j)) else beta_settling


def total_collision_frequency(
    diameter_i: Union[float, np.ndarray],
    diameter_j: Union[float, np.ndarray],
    settling_velocity_i: Union[float, np.ndarray],
    settling_velocity_j: Union[float, np.ndarray],
    shear_rate: Union[float, np.ndarray] = 0.0,
    temperature: float = WATER_TEMPERATURE_20C,
    viscosity: float = WATER_VISCOSITY_20C,
    collision_factor: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate total collision frequency from all mechanisms.
    
    Total βc = βc_Brownian + βc_shear + βc_settling
    
    Parameters
    ----------
    diameter_i : float or array-like
        Diameter of particle type i [m]
    diameter_j : float or array-like
        Diameter of particle type j [m]
    settling_velocity_i : float or array-like
        Settling velocity of particle type i [m/s]
    settling_velocity_j : float or array-like
        Settling velocity of particle type j [m/s]
    shear_rate : float or array-like, optional
        Flow shear rate [s⁻¹], default 0.0
    temperature : float, optional
        Absolute temperature [K], default 293.15 K
    viscosity : float, optional
        Dynamic viscosity [Pa·s], default 1.002e-3
    collision_factor : float, optional
        Collision diameter correction factor, default 1.0
        
    Returns
    -------
    float or ndarray
        Total collision frequency [m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-13
    """
    # Brownian motion contribution
    beta_brownian = collision_frequency_brownian(
        diameter_i, diameter_j, temperature, viscosity, collision_factor
    )
    
    # Fluid shear contribution
    beta_shear = collision_frequency_shear(
        diameter_i, diameter_j, shear_rate, collision_factor
    )
    
    # Differential settling contribution
    beta_settling = collision_frequency_settling(
        diameter_i, diameter_j, settling_velocity_i, settling_velocity_j, collision_factor
    )
    
    return beta_brownian + beta_shear + beta_settling


def limiting_floc_size(
    floc_strength: Union[float, np.ndarray],
    water_viscosity: float = WATER_VISCOSITY_20C,
    shear_rate: Union[float, np.ndarray] = 10.0,
    interpenetration_distance: float = 2e-6
) -> Union[float, np.ndarray]:
    """
    Calculate limiting floc size under shear.
    
    Equation 4-12: d_flim = 2τf ΔR / (μw γ̇)
    
    Parameters
    ----------
    floc_strength : float or array-like
        Floc strength [Pa]
    water_viscosity : float, optional
        Water viscosity [Pa·s], default 1.002e-3
    shear_rate : float or array-like, optional
        Local flow shear rate [s⁻¹], default 10.0
    interpenetration_distance : float, optional
        Interpenetration distance ΔR for colliding flocs [m], default 2e-6
        
    Returns
    -------
    float or ndarray
        Limiting floc diameter [m]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-12
    Krone (1963)
    """
    validate_positive(floc_strength, "floc_strength")
    validate_positive(water_viscosity, "water_viscosity")
    validate_positive(shear_rate, "shear_rate")
    validate_positive(interpenetration_distance, "interpenetration_distance")
    
    return (2 * floc_strength * interpenetration_distance) / (water_viscosity * shear_rate)


def floc_size_distribution(
    diameter: Union[float, np.ndarray],
    modal_diameter: float,
    cv0: float = 3.123e-6,
    xi: float = 0.608,
    kappa: float = 0.0055,
    floc_density: float = 1200.0,
    water_density: float = WATER_DENSITY,
    kinematic_viscosity: float = 1e-6
) -> Union[float, np.ndarray]:
    """
    Calculate floc size distribution using modal relationship.
    
    Equation 4-11: Cv = Cv0 * df^ξ * exp(-κf * df * ζ)
    Where ζ = g(ρf-ρw)/(8νρw)
    
    Parameters
    ----------
    diameter : float or array-like
        Floc diameter [m]
    modal_diameter : float
        Modal floc diameter [m]
    cv0 : float, optional
        Volume concentration coefficient [ppm], default 3.123e-6
    xi : float, optional
        Size distribution exponent, default 0.608
    kappa : float, optional
        Distribution shape parameter [s/μm], default 0.0055
    floc_density : float, optional
        Floc density [kg/m³], default 1200
    water_density : float, optional
        Water density [kg/m³], default 1000
    kinematic_viscosity : float, optional
        Kinematic viscosity [m²/s], default 1e-6
        
    Returns
    -------
    float or ndarray
        Volume concentration of sediment in size class
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-11
    Kranck (1986), Kranck and Milligan (1992)
    """
    validate_positive(diameter, "diameter")
    validate_positive(modal_diameter, "modal_diameter")
    validate_positive(cv0, "cv0")
    validate_positive(floc_density, "floc_density")
    
    df = np.atleast_1d(diameter) * 1e6  # Convert to μm for equation units
    
    # Calculate ζ parameter
    zeta = (GRAVITY * (floc_density - water_density)) / (8 * kinematic_viscosity * water_density)
    
    # Volume concentration
    Cv = cv0 * df**xi * np.exp(-kappa * df * zeta)
    
    return Cv.item() if np.isscalar(diameter) else Cv


def aggregate_orders_properties(
    order: int,
    base_density: float = 1269.0,
    base_strength: float = 2.2
) -> Tuple[float, float]:
    """
    Calculate floc properties for different orders of aggregation.
    
    Based on Krone's conceptual model of aggregate orders:
    - Order 0: Primary grains with strong bonds
    - Order 1: Aggregates of order 0 aggregates
    - Higher orders: Progressively weaker bonds
    
    Parameters
    ----------
    order : int
        Order of aggregation (0, 1, 2, ...)
    base_density : float, optional
        Base density for order 0 [kg/m³], default 1269
    base_strength : float, optional
        Base strength for order 0 [Pa], default 2.2
        
    Returns
    -------
    Tuple[float, float]
        (Aggregate density [kg/m³], Aggregate strength [Pa])
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.4.2
    Krone (1963), Table 4-7
    """
    if order < 0:
        raise ValueError("Aggregate order must be non-negative")
    
    # Empirical relationships from San Francisco Bay sediment (Krone 1963)
    # Density decreases with increasing order
    density_factors = [1.0, 0.928, 0.896, 0.877, 0.865, 0.857, 0.850]
    strength_factors = [1.0, 0.177, 0.064, 0.064, 0.037, 0.016, 0.009]
    
    if order >= len(density_factors):
        # Extrapolate for higher orders
        density_factor = density_factors[-1] * (0.99)**max(0, order - len(density_factors) + 1)
        strength_factor = strength_factors[-1] * (0.5)**max(0, order - len(strength_factors) + 1)
    else:
        density_factor = density_factors[order]
        strength_factor = strength_factors[order]
    
    aggregate_density = base_density * density_factor
    aggregate_strength = base_strength * strength_factor
    
    return aggregate_density, aggregate_strength


def floc_equilibrium_size(
    concentration: Union[float, np.ndarray],
    shear_rate: Union[float, np.ndarray],
    formation_coeff: float = 1.0,
    breakup_coeff: float = 1.0,
    formation_exponent: float = 1.0,
    breakup_exponent: float = 0.5
) -> Union[float, np.ndarray]:
    """
    Calculate equilibrium floc size from balance of formation and breakup.
    
    Simplified model: d_eq = (formation_rate/breakup_rate)^(1/n)
    Formation rate ∝ C^a, Breakup rate ∝ G^b
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    shear_rate : float or array-like
        Flow shear rate [s⁻¹]
    formation_coeff : float, optional
        Formation coefficient, default 1.0
    breakup_coeff : float, optional
        Breakup coefficient, default 1.0
    formation_exponent : float, optional
        Formation exponent a, default 1.0
    breakup_exponent : float, optional
        Breakup exponent b, default 0.5
        
    Returns
    -------
    float or ndarray
        Equilibrium floc diameter [m]
        
    References
    ----------
    Based on concepts from Winterwerp (1998), Krone (1963)
    """
    validate_positive(concentration, "concentration")
    validate_positive(shear_rate, "shear_rate")
    
    C = np.atleast_1d(concentration)
    G = np.atleast_1d(shear_rate)
    
    formation_rate = formation_coeff * C**formation_exponent
    breakup_rate = breakup_coeff * G**breakup_exponent
    
    # Equilibrium when formation = breakup
    # Size scales as (formation/breakup)^(1/n)
    equilibrium_size = 1e-4 * (formation_rate / breakup_rate)**(1.0/3.0)  # Empirical scaling
    
    return equilibrium_size.item() if (np.isscalar(concentration) and 
                                      np.isscalar(shear_rate)) else equilibrium_size