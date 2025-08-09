"""
Measurement conversions and standardization functions for sediment transport data.

This module provides unit conversion functions and standardization methods
for sediment transport measurements, including concentration units, discharge,
particle sizes, and temperature corrections.
"""

import numpy as np
from typing import Union, Optional
from ..utils.validators import validate_positive, validate_range
from ..utils.constants import WATER_DENSITY, GRAVITY

def concentration_unit_conversion(
    concentration: Union[float, np.ndarray],
    from_units: str,
    to_units: str,
    water_density: float = WATER_DENSITY,
    sediment_density: float = 2650.0
) -> Union[float, np.ndarray]:
    """
    Convert between different sediment concentration units.
    
    Supported units:
    - mg/L (milligrams per liter)
    - ppm (parts per million by weight)
    - g/m³ (grams per cubic meter)
    - kg/m³ (kilograms per cubic meter)
    - percent_by_weight (%)
    - percent_by_volume (%)
    
    Parameters
    ----------
    concentration : float or array-like
        Concentration value(s) to convert
    from_units : str
        Source units
    to_units : str
        Target units
    water_density : float, optional
        Water density [kg/m³], default 1000.0
    sediment_density : float, optional
        Sediment density [kg/m³], default 2650.0
        
    Returns
    -------
    float or ndarray
        Converted concentration values
    """
    validate_positive(concentration, "concentration")
    validate_positive(water_density, "water_density")
    validate_positive(sediment_density, "sediment_density")
    
    # Conversion factors to mg/L (reference unit)
    to_mg_per_L = {
        "mg/L": 1.0,
        "ppm": 1.0,  # For dilute solutions, ppm ≈ mg/L
        "g/m³": 1.0,  # 1 g/m³ = 1 mg/L
        "kg/m³": 1000.0,
        "percent_by_weight": 10000.0 * water_density / 1000.0,  # % to mg/L
        "percent_by_volume": 10000.0 * sediment_density / 1000.0  # vol% to mg/L
    }
    
    # Convert from source units to mg/L
    if from_units not in to_mg_per_L:
        raise ValueError(f"Unsupported source units: {from_units}")
    
    concentration_mg_L = concentration * to_mg_per_L[from_units]
    
    # Convert from mg/L to target units
    if to_units not in to_mg_per_L:
        raise ValueError(f"Unsupported target units: {to_units}")
    
    result = concentration_mg_L / to_mg_per_L[to_units]
    
    return result


def discharge_unit_conversion(
    discharge: Union[float, np.ndarray],
    from_units: str,
    to_units: str
) -> Union[float, np.ndarray]:
    """
    Convert between different discharge units.
    
    Supported units:
    - m³/s (cubic meters per second)
    - L/s (liters per second)
    - ft³/s or cfs (cubic feet per second)
    - gpm (gallons per minute)
    - MGD (million gallons per day)
    
    Parameters
    ----------
    discharge : float or array-like
        Discharge value(s) to convert
    from_units : str
        Source units
    to_units : str
        Target units
        
    Returns
    -------
    float or ndarray
        Converted discharge values
    """
    validate_positive(discharge, "discharge")
    
    # Conversion factors to m³/s (reference unit)
    to_m3_per_s = {
        "m³/s": 1.0,
        "m3/s": 1.0,
        "L/s": 0.001,
        "ft³/s": 0.0283168,
        "cfs": 0.0283168,
        "gpm": 6.30902e-5,
        "MGD": 0.0437764
    }
    
    # Convert from source units to m³/s
    if from_units not in to_m3_per_s:
        raise ValueError(f"Unsupported source units: {from_units}")
    
    discharge_m3_s = discharge * to_m3_per_s[from_units]
    
    # Convert from m³/s to target units
    if to_units not in to_m3_per_s:
        raise ValueError(f"Unsupported target units: {to_units}")
    
    result = discharge_m3_s / to_m3_per_s[to_units]
    
    return result


def particle_size_unit_conversion(
    particle_size: Union[float, np.ndarray],
    from_units: str,
    to_units: str
) -> Union[float, np.ndarray]:
    """
    Convert between different particle size units.
    
    Supported units:
    - m (meters)
    - mm (millimeters)
    - μm or um (micrometers)
    - nm (nanometers)
    - in (inches)
    - phi (phi scale)
    
    Parameters
    ----------
    particle_size : float or array-like
        Particle size value(s) to convert
    from_units : str
        Source units
    to_units : str
        Target units
        
    Returns
    -------
    float or ndarray
        Converted particle size values
    """
    validate_positive(particle_size, "particle_size")
    
    # Special handling for phi scale
    if from_units == "phi" and to_units == "mm":
        return 2.0 ** (-particle_size)  # D = 2^(-φ)
    elif from_units == "mm" and to_units == "phi":
        return -np.log2(particle_size)  # φ = -log₂(D)
    
    # Conversion factors to meters (reference unit)
    to_meters = {
        "m": 1.0,
        "mm": 0.001,
        "μm": 1e-6,
        "um": 1e-6,
        "nm": 1e-9,
        "in": 0.0254
    }
    
    # Convert from source units to meters
    if from_units not in to_meters:
        raise ValueError(f"Unsupported source units: {from_units}")
    
    size_m = particle_size * to_meters[from_units]
    
    # Convert from meters to target units
    if to_units not in to_meters:
        raise ValueError(f"Unsupported target units: {to_units}")
    
    result = size_m / to_meters[to_units]
    
    return result


def temperature_correction_factor(
    temperature: Union[float, np.ndarray],
    reference_temperature: float = 20.0,
    property_type: str = "viscosity"
) -> Union[float, np.ndarray]:
    """
    Calculate temperature correction factors for water properties.
    
    Parameters
    ----------
    temperature : float or array-like
        Water temperature [°C]
    reference_temperature : float, optional
        Reference temperature [°C], default 20.0
    property_type : str, optional
        Property to correct: 'viscosity' or 'density'
        
    Returns
    -------
    float or ndarray
        Temperature correction factor
    """
    validate_range(temperature, "temperature", -5.0, 40.0)
    
    T = np.atleast_1d(temperature)
    T_ref = reference_temperature
    
    if property_type == "viscosity":
        # Viscosity temperature correction (simplified)
        # μ(T) = μ(20°C) * correction_factor
        correction_factor = np.exp(1.3272 * (20 - T) / (T + 105))
    elif property_type == "density":
        # Density temperature correction (simplified)
        # ρ(T) = ρ(20°C) * correction_factor  
        correction_factor = 1 - 6.8e-5 * (T - T_ref) - 8.5e-6 * (T - T_ref)**2
    else:
        raise ValueError("property_type must be 'viscosity' or 'density'")
    
    if np.isscalar(temperature):
        return float(correction_factor[0])
    else:
        return correction_factor


def standard_reference_conditions(
    value: Union[float, np.ndarray],
    actual_temperature: float,
    actual_pressure: float = 101325.0,
    reference_temperature: float = 20.0,
    reference_pressure: float = 101325.0,
    correction_type: str = "gas_volume"
) -> Union[float, np.ndarray]:
    """
    Correct measurements to standard reference conditions.
    
    Parameters
    ----------
    value : float or array-like
        Measured value to correct
    actual_temperature : float
        Actual temperature during measurement [°C]
    actual_pressure : float, optional
        Actual pressure [Pa], default 101325.0 (1 atm)
    reference_temperature : float, optional
        Reference temperature [°C], default 20.0
    reference_pressure : float, optional
        Reference pressure [Pa], default 101325.0
    correction_type : str, optional
        Type of correction: 'gas_volume' or 'liquid_density'
        
    Returns
    -------
    float or ndarray
        Value corrected to standard conditions
    """
    validate_positive(actual_pressure, "actual_pressure")
    validate_positive(reference_pressure, "reference_pressure")
    
    T_actual = actual_temperature + 273.15  # Convert to Kelvin
    T_ref = reference_temperature + 273.15
    
    if correction_type == "gas_volume":
        # Ideal gas law correction: V₂ = V₁ * (P₁/P₂) * (T₂/T₁)
        correction_factor = (actual_pressure / reference_pressure) * (T_ref / T_actual)
    elif correction_type == "liquid_density":
        # Simplified liquid density correction
        temp_correction = temperature_correction_factor(
            actual_temperature, reference_temperature, "density"
        )
        pressure_correction = 1.0  # Assume incompressible for small pressure changes
        correction_factor = temp_correction * pressure_correction
    else:
        raise ValueError("correction_type must be 'gas_volume' or 'liquid_density'")
    
    return value * correction_factor


def sediment_transport_rate_conversion(
    transport_rate: Union[float, np.ndarray],
    from_units: str,
    to_units: str,
    sediment_density: float = 2650.0,
    water_density: float = WATER_DENSITY
) -> Union[float, np.ndarray]:
    """
    Convert between sediment transport rate units.
    
    Supported units:
    - kg/s (mass per time)
    - tons/day
    - m³/s (volume per time)
    - kg/s/m (mass per time per width)
    - m³/s/m (volume per time per width)
    
    Parameters
    ----------
    transport_rate : float or array-like
        Transport rate value(s) to convert
    from_units : str
        Source units
    to_units : str
        Target units
    sediment_density : float, optional
        Sediment density [kg/m³], default 2650.0
    water_density : float, optional
        Water density [kg/m³], default 1000.0
        
    Returns
    -------
    float or ndarray
        Converted transport rate values
    """
    validate_positive(transport_rate, "transport_rate")
    validate_positive(sediment_density, "sediment_density")
    
    # Conversion factors to kg/s (reference unit)
    to_kg_per_s = {
        "kg/s": 1.0,
        "tons/day": 1000.0 / (24 * 3600),  # metric tons per day to kg/s
        "m³/s": sediment_density,  # volume to mass using sediment density
        "kg/s/m": 1.0,  # per unit width, same mass rate
        "m³/s/m": sediment_density  # volume per width to mass per width
    }
    
    # Convert from source units to kg/s
    if from_units not in to_kg_per_s:
        raise ValueError(f"Unsupported source units: {from_units}")
    
    rate_kg_s = transport_rate * to_kg_per_s[from_units]
    
    # Convert from kg/s to target units
    if to_units not in to_kg_per_s:
        raise ValueError(f"Unsupported target units: {to_units}")
    
    result = rate_kg_s / to_kg_per_s[to_units]
    
    return result


def velocity_unit_conversion(
    velocity: Union[float, np.ndarray],
    from_units: str,
    to_units: str
) -> Union[float, np.ndarray]:
    """
    Convert between velocity units.
    
    Supported units:
    - m/s (meters per second)
    - cm/s (centimeters per second)
    - ft/s (feet per second)
    - km/h (kilometers per hour)
    - mph (miles per hour)
    
    Parameters
    ----------
    velocity : float or array-like
        Velocity value(s) to convert
    from_units : str
        Source units
    to_units : str
        Target units
        
    Returns
    -------
    float or ndarray
        Converted velocity values
    """
    validate_positive(velocity, "velocity")
    
    # Conversion factors to m/s (reference unit)
    to_m_per_s = {
        "m/s": 1.0,
        "cm/s": 0.01,
        "ft/s": 0.3048,
        "km/h": 1.0 / 3.6,
        "mph": 0.44704
    }
    
    # Convert from source units to m/s
    if from_units not in to_m_per_s:
        raise ValueError(f"Unsupported source units: {from_units}")
    
    velocity_m_s = velocity * to_m_per_s[from_units]
    
    # Convert from m/s to target units
    if to_units not in to_m_per_s:
        raise ValueError(f"Unsupported target units: {to_units}")
    
    result = velocity_m_s / to_m_per_s[to_units]
    
    return result


def normalize_concentration(
    concentration: Union[float, np.ndarray],
    reference_condition: str = "standard_temperature",
    actual_temperature: Optional[float] = None,
    reference_temperature: float = 20.0
) -> Union[float, np.ndarray]:
    """
    Normalize concentration measurements to standard conditions.
    
    Parameters
    ----------
    concentration : float or array-like
        Measured concentration values
    reference_condition : str, optional
        Reference condition: 'standard_temperature' or 'field_temperature'
    actual_temperature : float, optional
        Actual measurement temperature [°C]
    reference_temperature : float, optional
        Reference temperature [°C], default 20.0
        
    Returns
    -------
    float or ndarray
        Normalized concentration values
    """
    if reference_condition == "standard_temperature" and actual_temperature is not None:
        # Apply temperature correction
        temp_factor = temperature_correction_factor(
            actual_temperature, reference_temperature, "density"
        )
        return concentration * temp_factor
    else:
        # No correction needed
        return concentration