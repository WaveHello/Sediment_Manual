"""
Unit conversion utilities for sediment transport calculations.
"""

import numpy as np
from typing import Union

# Length conversions
def meters_to_feet(meters: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert meters to feet."""
    return meters * 3.28084

def feet_to_meters(feet: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert feet to meters."""
    return feet / 3.28084

def mm_to_meters(millimeters: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert millimeters to meters."""
    return millimeters * 1e-3

def meters_to_mm(meters: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert meters to millimeters."""
    return meters * 1e3

# Flow rate conversions  
def cms_to_cfs(cms: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert cubic meters per second to cubic feet per second."""
    return cms * 35.3147

def cfs_to_cms(cfs: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert cubic feet per second to cubic meters per second."""
    return cfs / 35.3147

# Grain size conversions
def phi_to_mm(phi: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert phi scale to millimeters.
    
    D_mm = 2^(-φ)
    """
    return 2**(-phi)

def mm_to_phi(mm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert millimeters to phi scale.
    
    φ = -log₂(D_mm)
    """
    return -np.log2(mm)

# Density conversions
def specific_gravity_to_density(
    specific_gravity: Union[float, np.ndarray],
    reference_density: float = 1000.0
) -> Union[float, np.ndarray]:
    """Convert specific gravity to density."""
    return specific_gravity * reference_density

def density_to_specific_gravity(
    density: Union[float, np.ndarray], 
    reference_density: float = 1000.0
) -> Union[float, np.ndarray]:
    """Convert density to specific gravity."""
    return density / reference_density

# Temperature conversions
def celsius_to_fahrenheit(celsius: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert Celsius to Fahrenheit."""
    return celsius * 9/5 + 32

def fahrenheit_to_celsius(fahrenheit: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert Fahrenheit to Celsius."""
    return (fahrenheit - 32) * 5/9

# Sediment concentration conversions
def mg_per_l_to_kg_per_m3(mg_l: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert mg/L to kg/m³."""
    return mg_l * 1e-3

def kg_per_m3_to_mg_per_l(kg_m3: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert kg/m³ to mg/L."""
    return kg_m3 * 1e3

def ppm_to_kg_per_m3(ppm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert parts per million to kg/m³ (assuming water density = 1000 kg/m³)."""
    return ppm * 1e-3