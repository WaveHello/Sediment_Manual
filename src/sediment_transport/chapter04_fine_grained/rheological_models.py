"""
Rheological Models for Fine-Grained Sediment (Chapter 4)

This module implements rheological models for mud behavior including viscosity
models, viscoelastic properties, and constitutive relationships. Based on 
equations 4-1 to 4-8 from ASCE Manual 110, Chapter 4.
"""

import numpy as np
from typing import Union, Optional, Tuple
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY, WATER_VISCOSITY_20C

def sisko_viscosity_model(
    shear_rate: Union[float, np.ndarray],
    infinite_viscosity: float,
    consistency: float,
    flow_behavior_index: float
) -> Union[float, np.ndarray]:
    """
    Calculate dynamic viscosity using Sisko power-law model.
    
    Equation 4-1: μ = μ∞ + cr * γ̇^(nr-1)
    
    Where:
    - nr < 1: shear-thinning behavior
    - nr > 1: shear-thickening behavior  
    - nr = 1: Newtonian behavior (cr = 0)
    
    Parameters
    ----------
    shear_rate : float or array-like
        Shear rate γ̇ [s⁻¹]
    infinite_viscosity : float
        Viscosity at infinite shear rate μ∞ [Pa·s]
    consistency : float
        Consistency of non-Newtonian fluid cr [Pa·s^n]
    flow_behavior_index : float
        Flow behavior coefficient nr [-]
        
    Returns
    -------
    float or ndarray
        Dynamic viscosity [Pa·s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-1
    Sisko (1958)
    """
    validate_positive(shear_rate, "shear_rate")
    validate_positive(infinite_viscosity, "infinite_viscosity")
    validate_positive(consistency, "consistency")
    
    gamma_dot = np.atleast_1d(shear_rate)
    
    # Handle special case of nr = 1 (Newtonian)
    if np.abs(flow_behavior_index - 1.0) < 1e-10:
        viscosity = infinite_viscosity
    else:
        viscosity = infinite_viscosity + consistency * gamma_dot**(flow_behavior_index - 1)
    
    return viscosity.item() if np.isscalar(shear_rate) else viscosity


def concentration_dependent_viscosity(
    concentration: Union[float, np.ndarray],
    water_viscosity: float = WATER_VISCOSITY_20C,
    alpha_coeff: float = 1.68,
    beta_coeff: float = 0.346
) -> Union[float, np.ndarray]:
    """
    Calculate viscosity as function of suspension concentration.
    
    Equation 4-2: μ∞ = μw * (1 + αr * C^βr)
    
    Where C is dry sediment mass per unit volume of suspension
    
    Parameters
    ----------
    concentration : float or array-like
        Suspension concentration [kg/m³]
    water_viscosity : float, optional
        Water viscosity μw [Pa·s], default 1.002e-3
    alpha_coeff : float, optional
        Concentration coefficient αr, default 1.68
    beta_coeff : float, optional
        Concentration exponent βr, default 0.346
        
    Returns
    -------
    float or ndarray
        Viscosity at infinite shear rate [Pa·s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-2
    Engelund and Zhaohui (1984) - for kaolinite
    """
    validate_positive(concentration, "concentration")
    validate_positive(water_viscosity, "water_viscosity")
    
    C = np.atleast_1d(concentration)
    viscosity = water_viscosity * (1 + alpha_coeff * C**beta_coeff)
    
    return viscosity.item() if np.isscalar(concentration) else viscosity


def amazon_kinematic_viscosity(
    concentration: Union[float, np.ndarray],
    shear_rate: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate kinematic viscosity for Amazon sediment.
    
    Equations 4-3a,b: 
    ν = C * exp(-0.78γ̇ + 10.24) for γ̇ < 3.9 Hz
    ν = C * exp(-0.017γ̇ + 12.95) for γ̇ ≥ 3.9 Hz
    
    Parameters
    ----------
    concentration : float or array-like
        Sediment concentration [kg/m³]
    shear_rate : float or array-like
        Shear rate [Hz]
        
    Returns
    -------
    float or ndarray
        Kinematic viscosity [m²/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equations 4-3a,b
    Vinzon (1998)
    """
    validate_positive(concentration, "concentration")
    validate_positive(shear_rate, "shear_rate")
    
    C = np.atleast_1d(concentration)
    gamma = np.atleast_1d(shear_rate)
    
    # Conditions for low and high shear rates
    low_shear = gamma < 3.9
    high_shear = gamma >= 3.9
    
    viscosity = np.zeros_like(gamma * C)
    
    # Low shear rate regime
    if np.any(low_shear):
        viscosity[low_shear] = C * np.exp(-0.78 * gamma[low_shear] + 10.24)
    
    # High shear rate regime
    if np.any(high_shear):
        viscosity[high_shear] = C * np.exp(-0.017 * gamma[high_shear] + 12.95)
    
    return viscosity.item() if (np.isscalar(concentration) and 
                               np.isscalar(shear_rate)) else viscosity


def mixture_cec_viscosity_parameters(
    kaolinite_fraction: float,
    attapulgite_fraction: float,
    bentonite_fraction: float,
    water_fraction: float,
    kaolinite_cec: float = 6.0,
    attapulgite_cec: float = 28.0,
    bentonite_cec: float = 105.0
) -> Tuple[float, float, float]:
    """
    Calculate viscosity parameters for clay mixtures based on CEC.
    
    Equation 4-4: CECs = fK*CECK + fA*CECA + fB*CECB
    
    Jinchai (1998) relationships:
    μ∞ = 0.05*CECs + 0.001
    nr = -0.033*CECs + 0.28
    log cr = 0.13*CECs - 0.22
    
    Parameters
    ----------
    kaolinite_fraction : float
        Weight fraction of kaolinite fK
    attapulgite_fraction : float
        Weight fraction of attapulgite fA
    bentonite_fraction : float
        Weight fraction of bentonite fB
    water_fraction : float
        Weight fraction of water fw
    kaolinite_cec : float, optional
        CEC of kaolinite [meq/100g], default 6.0
    attapulgite_cec : float, optional
        CEC of attapulgite [meq/100g], default 28.0
    bentonite_cec : float, optional
        CEC of bentonite [meq/100g], default 105.0
        
    Returns
    -------
    Tuple[float, float, float]
        (μ∞ [Pa·s], nr [-], cr [Pa·s^n])
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-4
    Jinchai (1998)
    """
    validate_range(kaolinite_fraction, "kaolinite_fraction", 0.0, 1.0)
    validate_range(attapulgite_fraction, "attapulgite_fraction", 0.0, 1.0)
    validate_range(bentonite_fraction, "bentonite_fraction", 0.0, 1.0)
    validate_range(water_fraction, "water_fraction", 0.0, 1.0)
    
    # Check that fractions sum to 1
    total_fraction = kaolinite_fraction + attapulgite_fraction + bentonite_fraction + water_fraction
    if not np.isclose(total_fraction, 1.0, rtol=1e-3):
        raise ValueError("Weight fractions must sum to 1.0")
    
    # Calculate mixture CEC
    cec_s = (kaolinite_fraction * kaolinite_cec + 
             attapulgite_fraction * attapulgite_cec + 
             bentonite_fraction * bentonite_cec)
    
    # Calculate viscosity parameters
    mu_infinity = 0.05 * cec_s + 0.001
    nr = -0.033 * cec_s + 0.28
    log_cr = 0.13 * cec_s - 0.22
    cr = 10**log_cr
    
    return mu_infinity, nr, cr


def bulk_density_from_concentration(
    concentration: Union[float, np.ndarray],
    water_density: float = WATER_DENSITY,
    grain_density: float = 2650.0
) -> Union[float, np.ndarray]:
    """
    Calculate bulk density from dry mass concentration.
    
    Equation 4-5: ρ = (ρw * ρs) / (ρs - C * (ρs - ρw) / ρs)
    
    Parameters
    ----------
    concentration : float or array-like
        Dry mass concentration C [kg/m³]
    water_density : float, optional
        Water density ρw [kg/m³], default 1000
    grain_density : float, optional
        Grain density ρs [kg/m³], default 2650
        
    Returns
    -------
    float or ndarray
        Bulk density ρ [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-5
    """
    validate_positive(concentration, "concentration")
    validate_positive(water_density, "water_density")
    validate_positive(grain_density, "grain_density")
    
    C = np.atleast_1d(concentration)
    
    # Solids volume fraction
    phi = C / grain_density
    
    # Bulk density
    rho = water_density + C * (grain_density - water_density) / grain_density
    
    return rho.item() if np.isscalar(concentration) else rho


def standard_solid_model(
    strain: Union[float, np.ndarray],
    strain_rate: Union[float, np.ndarray],
    stress: Union[float, np.ndarray],
    stress_rate: Union[float, np.ndarray],
    viscosity: float,
    shear_modulus_1: float,
    shear_modulus_2: float
) -> Union[float, np.ndarray]:
    """
    Standard solid viscoelastic constitutive equation.
    
    Equation 4-6: τ + (μ/(G1+G2)) * τ̇ = (G1*G2/(G1+G2)) * γ + (μ*G1/(G1+G2)) * γ̇
    
    Parameters
    ---------- 
    \u00b5
    strain : float or array-like
        Shear strain γ [-]
    strain_rate : float or array-like
        Strain rate γ̇ [s⁻¹]
    stress : float or array-like
        Shear stress τ [Pa]
    stress_rate : float or array-like
        Stress rate τ̇ [Pa/s]
    viscosity : float
        Viscosity μ [Pa·s]
    shear_modulus_1 : float
        First shear modulus G1 [Pa]
    shear_modulus_2 : float
        Second shear modulus G2 [Pa]
        
    Returns
    -------
    float or ndarray
        Constitutive relationship residual
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-6
    Keedwell (1984)
    """
    validate_positive(viscosity, "viscosity")
    validate_positive(shear_modulus_1, "shear_modulus_1")
    validate_positive(shear_modulus_2, "shear_modulus_2")
    
    tau = np.atleast_1d(stress)
    tau_dot = np.atleast_1d(stress_rate)
    gamma = np.atleast_1d(strain)
    gamma_dot = np.atleast_1d(strain_rate)
    
    G1 = shear_modulus_1
    G2 = shear_modulus_2
    mu = viscosity
    
    # Left side of equation
    lhs = tau + (mu / (G1 + G2)) * tau_dot
    
    # Right side of equation
    rhs = (G1 * G2 / (G1 + G2)) * gamma + (mu * G1 / (G1 + G2)) * gamma_dot
    
    residual = lhs - rhs
    
    return residual.item() if (np.isscalar(stress) and np.isscalar(strain)) else residual


def voigt_model(
    strain: Union[float, np.ndarray],
    strain_rate: Union[float, np.ndarray],
    viscosity: float,
    shear_modulus: float
) -> Union[float, np.ndarray]:
    """
    Voigt (Kelvin) viscoelastic constitutive equation.
    
    Equation 4-7: τ = G * γ + μ * γ̇
    
    Parameters
    ----------
    strain : float or array-like
        Shear strain γ [-]
    strain_rate : float or array-like
        Strain rate γ̇ [s⁻¹]
    viscosity : float
        Viscosity μ [Pa·s]
    shear_modulus : float
        Shear modulus G [Pa]
        
    Returns
    -------
    float or ndarray
        Shear stress τ [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-7
    """
    validate_positive(viscosity, "viscosity")
    validate_positive(shear_modulus, "shear_modulus")
    
    gamma = np.atleast_1d(strain)
    gamma_dot = np.atleast_1d(strain_rate)
    
    stress = shear_modulus * gamma + viscosity * gamma_dot
    
    return stress.item() if (np.isscalar(strain) and np.isscalar(strain_rate)) else stress


def frequency_dependent_parameters(
    frequency: Union[float, np.ndarray],
    alpha_rh: float,
    beta_rh: float
) -> Union[float, np.ndarray]:
    """
    Calculate frequency-dependent viscoelastic parameters.
    
    Equation 4-8: G1, G2, μ = exp(αrh) * f^βrh
    
    Parameters
    ----------
    frequency : float or array-like
        Forcing frequency f [Hz]
    alpha_rh : float
        Frequency coefficient αrh
    beta_rh : float
        Frequency exponent βrh
        
    Returns
    -------
    float or ndarray
        Parameter value at given frequency
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-8
    Jiang and Mehta (1995)
    """
    validate_positive(frequency, "frequency")
    validate_range(frequency, "frequency", 0.02, 40.0)
    
    f = np.atleast_1d(frequency)
    parameter = np.exp(alpha_rh) * f**beta_rh
    
    return parameter.item() if np.isscalar(frequency) else parameter


def shear_wave_velocity_modulus(
    density: Union[float, np.ndarray],
    shear_wave_velocity: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate shear modulus from shear wave velocity.
    
    For Voigt solid: G ≈ ρ * V²
    
    Parameters
    ----------
    density : float or array-like
        Material density [kg/m³]
    shear_wave_velocity : float or array-like
        Shear wave velocity [m/s]
        
    Returns
    -------
    float or ndarray
        Shear modulus [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.2.2
    Mehta et al. (1995)
    """
    validate_positive(density, "density")
    validate_positive(shear_wave_velocity, "shear_wave_velocity")
    
    rho = np.atleast_1d(density)
    V = np.atleast_1d(shear_wave_velocity)
    
    shear_modulus = rho * V**2
    
    return shear_modulus.item() if (np.isscalar(density) and 
                                   np.isscalar(shear_wave_velocity)) else shear_modulus


def temperature_viscosity_relationship(
    temperature: Union[float, np.ndarray],
    reference_temp: float = 293.15,
    reference_viscosity: float = WATER_VISCOSITY_20C
) -> Union[float, np.ndarray]:
    """
    Calculate temperature-dependent viscosity using Arrhenius relationship.
    
    Based on momentum exchange theory: ln(μ) varies linearly with 1/T
    
    Parameters
    ----------
    temperature : float or array-like
        Absolute temperature [K]
    reference_temp : float, optional
        Reference temperature [K], default 293.15 (20°C)
    reference_viscosity : float, optional
        Viscosity at reference temperature [Pa·s], default 1.002e-3
        
    Returns
    -------
    float or ndarray
        Temperature-corrected viscosity [Pa·s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.2.2
    Krone (1983)
    """
    validate_positive(temperature, "temperature")
    validate_range(temperature, "temperature", 273.0, 323.0)  # 0°C to 50°C
    validate_positive(reference_temp, "reference_temp")
    validate_positive(reference_viscosity, "reference_viscosity")
    
    T = np.atleast_1d(temperature)
    
    # Arrhenius-type relationship for water viscosity
    # μ(T) = μ_ref * exp(A * (1/T - 1/T_ref))
    # A ≈ 1800 K for water
    A = 1800.0  # Activation energy parameter
    
    viscosity = reference_viscosity * np.exp(A * (1/T - 1/reference_temp))
    
    return viscosity.item() if np.isscalar(temperature) else viscosity


def yield_stress_bingham_model(
    shear_rate: Union[float, np.ndarray],
    yield_stress: float,
    plastic_viscosity: float
) -> Union[float, np.ndarray]:
    """
    Calculate shear stress using Bingham plastic model.
    
    τ = τB + μp * γ̇ for τ > τB
    τ = 0 for τ ≤ τB (solid behavior)
    
    Parameters
    ----------
    shear_rate : float or array-like
        Shear rate [s⁻¹]
    yield_stress : float
        Bingham yield stress τB [Pa]
    plastic_viscosity : float
        Plastic viscosity μp [Pa·s]
        
    Returns
    -------
    float or ndarray
        Shear stress [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Figure 4-1
    """
    validate_positive(shear_rate, "shear_rate")
    validate_positive(yield_stress, "yield_stress")
    validate_positive(plastic_viscosity, "plastic_viscosity")
    
    gamma_dot = np.atleast_1d(shear_rate)
    
    # Bingham model: linear relationship above yield stress
    stress = yield_stress + plastic_viscosity * gamma_dot
    
    return stress.item() if np.isscalar(shear_rate) else stress