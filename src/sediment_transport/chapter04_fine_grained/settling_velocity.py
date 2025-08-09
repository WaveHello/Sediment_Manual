"""
Settling Velocity Models for Fine-Grained Sediment (Chapter 4)

This module implements settling velocity equations for fine-grained cohesive sediment,
including concentration-dependent settling, flocculation effects, and hindered settling.
Based on equations 4-14 to 4-23 from ASCE Manual 110, Chapter 4.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY, QUARTZ_DENSITY

def concentration_dependent_settling_velocity(
    concentration: Union[float, np.ndarray],
    free_settling_velocity: float,
    velocity_scaling_coeff: float = 0.048,
    flocculation_exponent: float = 1.5,
    hindered_coeff: float = 10.0,
    hindered_exponent: float = 1.0,
    c1: float = 0.2,
    c2: float = 5.0,
    c3: float = 75.0
) -> Union[float, np.ndarray]:
    """
    Calculate concentration-dependent settling velocity using four-zone model.
    
    Equation 4-14: Four-zone settling velocity model
    - Free settling: ws = wsf (C < C1)
    - Flocculation settling: ws = aw * C^nw / (C^2 + bw^2)^mw (C1 < C < C2)  
    - Hindered settling: ws = aw * C^(nw-2mw) (C2 < C < C3)
    - Negligible settling: ws ≈ 0 (C > C3)
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    free_settling_velocity : float
        Free settling velocity for individual particles [m/s]
    velocity_scaling_coeff : float, optional
        Velocity scaling coefficient aw [m/s], default 0.048
    flocculation_exponent : float, optional
        Flocculation settling exponent nw, default 1.5
    hindered_coeff : float, optional
        Hindered settling coefficient bw [kg/m³], default 10.0
    hindered_exponent : float, optional
        Hindered settling exponent mw, default 1.0
    c1 : float, optional
        Concentration limit for free settling [kg/m³], default 0.2
    c2 : float, optional
        Concentration for peak settling velocity [kg/m³], default 5.0
    c3 : float, optional
        Concentration limit for hindered settling [kg/m³], default 75.0
        
    Returns
    -------
    float or ndarray
        Settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-14
    Hwang (1989), Krone (1962)
    """
    validate_positive(concentration, "concentration")
    validate_positive(free_settling_velocity, "free_settling_velocity")
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    validate_positive(hindered_coeff, "hindered_coeff")
    
    C = np.atleast_1d(concentration)
    ws = np.zeros_like(C, dtype=float)
    
    # Free settling zone (C < C1)
    free_mask = C < c1
    ws[free_mask] = free_settling_velocity
    
    # Flocculation settling zone (C1 <= C <= C2)
    floc_mask = (C >= c1) & (C <= c2)
    ws[floc_mask] = velocity_scaling_coeff * (
        C[floc_mask]**flocculation_exponent / 
        (C[floc_mask]**2 + hindered_coeff**2)**(hindered_exponent/2)
    )
    
    # Hindered settling zone (C2 < C < C3)
    hindered_mask = (C > c2) & (C < c3)
    ws[hindered_mask] = velocity_scaling_coeff * C[hindered_mask]**(
        flocculation_exponent - 2*hindered_exponent
    )
    
    # Negligible settling zone (C >= C3)
    negligible_mask = C >= c3
    ws[negligible_mask] = 1e-8  # Essentially zero but avoid division issues
    
    return ws.item() if np.isscalar(concentration) else ws


def krone_flocculation_settling(
    concentration: Union[float, np.ndarray],
    velocity_scaling_coeff: float = 0.048,
    hindered_coeff: float = 25.0
) -> Union[float, np.ndarray]:
    """
    Calculate flocculation settling velocity using Krone's formula.
    
    Equation 4-15: ws = aw * bw^(-2mw) * C^nw
    Derived for low concentration condition (C << bw) from Equation 4-14.
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    velocity_scaling_coeff : float, optional
        Velocity scaling coefficient aw [m/s], default 0.048
    hindered_coeff : float, optional
        Hindered settling coefficient bw [kg/m³], default 25.0
        
    Returns
    -------
    float or ndarray
        Settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-15
    Krone (1962)
    """
    validate_positive(concentration, "concentration")
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    validate_positive(hindered_coeff, "hindered_coeff")
    
    # Using nw = 4/3 and mw = 1 from Krone's analysis
    nw = 4.0/3.0
    mw = 1.0
    
    return velocity_scaling_coeff * hindered_coeff**(-2*mw) * concentration**nw


def hindered_settling_velocity(
    concentration: Union[float, np.ndarray],
    velocity_scaling_coeff: float = 0.048,
    flocculation_exponent: float = 1.5,
    hindered_exponent: float = 1.0
) -> Union[float, np.ndarray]:
    """
    Calculate hindered settling velocity for high concentration conditions.
    
    Equation 4-16: ws = aw * C^(nw-2mw)
    Derived for high concentration condition (C >> bw) from Equation 4-14.
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    velocity_scaling_coeff : float, optional
        Velocity scaling coefficient aw [m/s], default 0.048
    flocculation_exponent : float, optional
        Flocculation settling exponent nw, default 1.5
    hindered_exponent : float, optional
        Hindered settling exponent mw, default 1.0
        
    Returns
    -------
    float or ndarray
        Settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-16
    Richardson and Zaki (1954)
    """
    validate_positive(concentration, "concentration")
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    
    # Note: nw - 2mw must be < 0 for hindered settling
    exponent = flocculation_exponent - 2*hindered_exponent
    if exponent >= 0:
        raise ValueError("For hindered settling, nw - 2mw must be negative")
    
    return velocity_scaling_coeff * concentration**exponent


def peak_settling_velocity(
    velocity_scaling_coeff: float,
    flocculation_exponent: float,
    hindered_coeff: float,
    hindered_exponent: float
) -> float:
    """
    Calculate peak settling velocity and corresponding concentration.
    
    Equation 4-17: wsm = aw * bw^(nw-2mw) * [(2mw)/(nw) - 1]^(nw/2) * [(2mw)/(nw)]^(-mw)
    
    Parameters
    ----------
    velocity_scaling_coeff : float
        Velocity scaling coefficient aw [m/s]
    flocculation_exponent : float
        Flocculation settling exponent nw
    hindered_coeff : float
        Hindered settling coefficient bw [kg/m³]
    hindered_exponent : float
        Hindered settling exponent mw
        
    Returns
    -------
    float
        Peak settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-17
    """
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    validate_positive(hindered_coeff, "hindered_coeff")
    
    nw = flocculation_exponent
    mw = hindered_exponent
    bw = hindered_coeff
    
    if nw <= 2*mw:
        raise ValueError("For peak velocity calculation, nw must be > 2mw")
    
    term1 = velocity_scaling_coeff * bw**(nw - 2*mw)
    term2 = ((2*mw/nw) - 1)**(nw/2)
    term3 = (2*mw/nw)**(-mw)
    
    return term1 * term2 * term3


def peak_concentration(
    hindered_coeff: float,
    flocculation_exponent: float,
    hindered_exponent: float
) -> float:
    """
    Calculate concentration at peak settling velocity.
    
    Equation 4-18: C2 = bw * [(2mw)/(nw) - 1]^(1/2)
    
    Parameters
    ----------
    hindered_coeff : float
        Hindered settling coefficient bw [kg/m³]
    flocculation_exponent : float
        Flocculation settling exponent nw
    hindered_exponent : float
        Hindered settling exponent mw
        
    Returns
    -------
    float
        Concentration at peak settling velocity [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-18
    """
    validate_positive(hindered_coeff, "hindered_coeff")
    
    nw = flocculation_exponent
    mw = hindered_exponent
    
    if nw <= 2*mw:
        raise ValueError("For peak concentration calculation, nw must be > 2mw")
    
    return hindered_coeff * ((2*mw/nw) - 1)**(0.5)


def settling_flux(
    concentration: Union[float, np.ndarray],
    velocity_scaling_coeff: float,
    flocculation_exponent: float,
    hindered_coeff: float,
    hindered_exponent: float
) -> Union[float, np.ndarray]:
    """
    Calculate settling flux Fs = ws * C.
    
    Equation 4-19: Fs = aw * C^(nw+1) / (C^2 + bw^2)^mw
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    velocity_scaling_coeff : float
        Velocity scaling coefficient aw [m/s]
    flocculation_exponent : float
        Flocculation settling exponent nw
    hindered_coeff : float
        Hindered settling coefficient bw [kg/m³]
    hindered_exponent : float
        Hindered settling exponent mw
        
    Returns
    -------
    float or ndarray
        Settling flux [kg/m²/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-19
    """
    validate_positive(concentration, "concentration")
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    validate_positive(hindered_coeff, "hindered_coeff")
    
    C = np.atleast_1d(concentration)
    numerator = velocity_scaling_coeff * C**(flocculation_exponent + 1)
    denominator = (C**2 + hindered_coeff**2)**hindered_exponent
    
    flux = numerator / denominator
    return flux.item() if np.isscalar(concentration) else flux


def maximum_settling_flux(
    velocity_scaling_coeff: float,
    flocculation_exponent: float,
    hindered_coeff: float,
    hindered_exponent: float
) -> Tuple[float, float]:
    """
    Calculate maximum settling flux and corresponding concentration.
    
    Equation 4-20: Maximum settling flux
    Equation 4-21: Concentration at maximum flux C2'
    
    Parameters
    ----------
    velocity_scaling_coeff : float
        Velocity scaling coefficient aw [m/s]
    flocculation_exponent : float
        Flocculation settling exponent nw
    hindered_coeff : float
        Hindered settling coefficient bw [kg/m³]
    hindered_exponent : float
        Hindered settling exponent mw
        
    Returns
    -------
    Tuple[float, float]
        (Maximum settling flux [kg/m²/s], Concentration at max flux [kg/m³])
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equations 4-20, 4-21
    """
    validate_positive(velocity_scaling_coeff, "velocity_scaling_coeff")
    validate_positive(hindered_coeff, "hindered_coeff")
    
    nw = flocculation_exponent
    mw = hindered_exponent
    bw = hindered_coeff
    
    # Concentration at maximum flux (Equation 4-21)
    c2_prime = bw * ((2*mw)/(nw + 1) - 1)**(0.5)
    
    # Maximum settling flux (Equation 4-20)
    term1 = velocity_scaling_coeff * bw**(nw + 1 - 2*mw)
    term2 = ((2*mw)/(nw + 1) - 1)**((nw + 1)/2)
    term3 = ((2*mw)/(nw + 1))**(-mw)
    
    fsm = term1 * term2 * term3
    
    return fsm, c2_prime


def temperature_corrected_settling_velocity(
    concentration: Union[float, np.ndarray],
    temperature: float,
    reference_temp: float = 15.0,
    reference_velocity_func: callable = None,
    **velocity_params
) -> Union[float, np.ndarray]:
    """
    Calculate temperature-corrected settling velocity.
    
    Equations 4-22a,b: Temperature correction for settling velocity
    ws50(C,Tc) = Φ * ws50(C,15)
    Φ = 1.776(1 + 0.875T')
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    temperature : float
        Water temperature [°C]
    reference_temp : float, optional
        Reference temperature [°C], default 15.0
    reference_velocity_func : callable, optional
        Function to calculate reference settling velocity
    **velocity_params
        Additional parameters for velocity function
        
    Returns
    -------
    float or ndarray
        Temperature-corrected settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equations 4-22a,b
    Jiang (1999), Lau (1994)
    """
    validate_positive(concentration, "concentration")
    validate_range(temperature, "temperature", -5.0, 50.0)
    
    # Default to concentration-dependent settling if no function provided
    if reference_velocity_func is None:
        reference_velocity_func = concentration_dependent_settling_velocity
    
    # Calculate reference settling velocity at 15°C
    ws_ref = reference_velocity_func(concentration, **velocity_params)
    
    # Temperature correction factor
    T_prime = temperature / reference_temp
    phi = 1.776 * (1 + 0.875 * T_prime)
    
    return phi * ws_ref


def shear_corrected_settling_velocity(
    concentration: Union[float, np.ndarray],
    shear_rate: Union[float, np.ndarray],
    reference_velocity_func: callable = None,
    lambda1: float = 266.0,
    lambda2: float = 9.0,
    **velocity_params
) -> Union[float, np.ndarray]:
    """
    Calculate shear-corrected settling velocity.
    
    Equation 4-23: Shear rate correction for settling velocity
    ws|γ = ws|γ=0 * [1 + λ1*γ] / [1 + λ2*γ^2]
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    shear_rate : float or array-like
        Flow shear rate [Hz or s⁻¹]
    reference_velocity_func : callable, optional
        Function to calculate reference settling velocity
    lambda1 : float, optional
        First shear correction coefficient, default 266.0
    lambda2 : float, optional
        Second shear correction coefficient, default 9.0
    **velocity_params
        Additional parameters for velocity function
        
    Returns
    -------
    float or ndarray
        Shear-corrected settling velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-23
    van Leussen (1994), Teeter (2001a)
    """
    validate_positive(concentration, "concentration")
    validate_positive(shear_rate, "shear_rate")
    
    # Default to concentration-dependent settling if no function provided
    if reference_velocity_func is None:
        reference_velocity_func = concentration_dependent_settling_velocity
    
    # Calculate reference settling velocity at zero shear
    ws_ref = reference_velocity_func(concentration, **velocity_params)
    
    # Shear correction factor
    gamma = np.atleast_1d(shear_rate)
    correction_factor = (1 + lambda1 * gamma) / (1 + lambda2 * gamma**2)
    
    ws_corrected = ws_ref * correction_factor
    
    if np.isscalar(shear_rate):
        return ws_corrected.item() if np.isscalar(concentration) else ws_corrected
    else:
        return ws_corrected