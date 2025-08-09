"""
Chapter 4: Fine-Grained Sediment Transport

This package implements comprehensive models for fine-grained cohesive sediment
transport processes based on ASCE Manual 110, Chapter 4.

Modules:
--------
- sediment_characterization: Size classification, cohesion assessment, clay properties
- rheological_models: Viscosity models, viscoelastic behavior, constitutive equations  
- floc_aggregation: Fractal models, collision mechanisms, aggregate properties
- settling_velocity: Concentration-dependent settling, flocculation effects
- deposition_transport: Deposition rates, transport formulations, bed exchange
- consolidation_gelling: Consolidation mechanics, gelling processes, time effects
- environmental_effects: Temperature, salinity, organic content, pH effects

Key Features:
------------
- Complete implementation of equations 4-1 through 4-37
- Support for both single-class and multi-class sediment
- Fractal-based floc models with realistic density relationships
- Temperature and salinity corrections for all processes
- Comprehensive rheological models from Newtonian to viscoelastic
- Environmental factor effects on all transport processes

References:
----------
ASCE Manual 110: Sedimentation Engineering - Processes, Measurements, 
Modeling, and Practice. Chapter 4: Fine-Grained Sediment Transport.
Authored by Ashish J. Mehta and William H. McAnally.

Example Usage:
-------------
>>> from sediment_transport.chapter04_fine_grained import settling_velocity as sv
>>> from sediment_transport.chapter04_fine_grained import floc_aggregation as fa
>>> 
>>> # Calculate concentration-dependent settling velocity
>>> concentration = 5.0  # kg/m³
>>> ws = sv.concentration_dependent_settling_velocity(concentration)
>>> 
>>> # Calculate floc density using fractal relationship
>>> floc_diameter = 100e-6  # m
>>> rho_floc = fa.fractal_floc_density(floc_diameter)
>>> 
>>> # Calculate collision frequency
>>> beta = fa.total_collision_frequency(10e-6, 50e-6, 1e-4, 2e-4, shear_rate=5.0)
"""

from . import sediment_characterization
from . import rheological_models  
from . import floc_aggregation
from . import settling_velocity
from . import deposition_transport
from . import consolidation_gelling
from . import environmental_effects

__all__ = [
    'sediment_characterization',
    'rheological_models',
    'floc_aggregation', 
    'settling_velocity',
    'deposition_transport',
    'consolidation_gelling',
    'environmental_effects'
]

# Version information
__version__ = '1.0.0'
__author__ = 'Based on ASCE Manual 110, Chapter 4 by Ashish J. Mehta and William H. McAnally'

# Key equations implemented
IMPLEMENTED_EQUATIONS = [
    '4-1: Sisko power-law viscosity model',
    '4-2: Concentration-dependent viscosity', 
    '4-3: Amazon sediment kinematic viscosity',
    '4-4: Mixture CEC calculation',
    '4-5: Bulk density from concentration',
    '4-6: Standard solid viscoelastic model',
    '4-7: Voigt viscoelastic model',
    '4-8: Frequency-dependent parameters',
    '4-9: Fractal floc density relationship',
    '4-10: Floc strength model',
    '4-11: Floc size distribution',
    '4-12: Limiting floc size under shear',
    '4-13: Collision frequency mechanisms',
    '4-14: Four-zone settling velocity model',
    '4-15: Krone flocculation settling',
    '4-16: Hindered settling velocity',
    '4-17: Peak settling velocity',
    '4-18: Peak concentration',
    '4-19: Settling flux',
    '4-20: Maximum settling flux',
    '4-21: Concentration at maximum flux',
    '4-22: Temperature correction factors',
    '4-23: Shear rate correction',
    '4-24: Single-class deposition rate',
    '4-25: Exponential concentration decay',
    '4-26: Multi-class deposition rates',
    '4-27: Multi-class concentration evolution',
    '4-28: Consolidation continuity equation',
    '4-29: Two-mode consolidation model'
]

# Key tables implemented
IMPLEMENTED_TABLES = [
    'Table 4-1: Sediment size and cohesion classification',
    'Table 4-2: Clay mineral properties and critical salinities',
    'Table 4-3: Characterization test requirements',
    'Table 4-4: Mud properties for rheometry',
    'Table 4-5: Viscoelastic model coefficients', 
    'Table 4-6: Rheological parameters at reference frequency',
    'Table 4-7: Aggregate order properties',
    'Table 4-8: Settling velocity model coefficients',
    'Table 4-9: Consolidation model parameters',
    'Table 4-10: Shear strength vs. solid volume fraction',
    'Table 4-11: Soil consistency classification',
    'Table 4-12: Organic-rich sediment parameters',
    'Table 4-13: Erosion rate constant coefficients'
]

# Quick access to commonly used functions
from .settling_velocity import concentration_dependent_settling_velocity
from .floc_aggregation import fractal_floc_density, total_collision_frequency
from .rheological_models import sisko_viscosity_model, concentration_dependent_viscosity
from .sediment_characterization import classify_sediment_by_size, get_clay_mineral_properties
from .deposition_transport import single_class_deposition_rate, multiclass_deposition_rate
from .consolidation_gelling import two_mode_consolidation_velocity
from .environmental_effects import temperature_effect_on_flocculation, organic_content_effect_on_settling

# Expose commonly used functions at package level
__all__.extend([
    'concentration_dependent_settling_velocity',
    'fractal_floc_density',
    'total_collision_frequency', 
    'sisko_viscosity_model',
    'concentration_dependent_viscosity',
    'classify_sediment_by_size',
    'get_clay_mineral_properties',
    'single_class_deposition_rate',
    'multiclass_deposition_rate',
    'two_mode_consolidation_velocity',
    'temperature_effect_on_flocculation',
    'organic_content_effect_on_settling'
])