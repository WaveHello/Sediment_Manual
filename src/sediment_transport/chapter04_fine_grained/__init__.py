"""
Chapter 4: Fine-Grained Sediment Transport

This package implements mathematical models and equations for fine-grained sediment
transport, including cohesive sediment behavior, flocculation, rheological properties,
and transport processes.

Modules:
    sediment_characterization: Clay minerals, cohesion, particle size relationships
    rheological_models: Viscosity models, viscoelastic properties, mud rheology
    floc_aggregation: Fractal models, aggregation mechanics, collision functions
    settling_velocity: Concentration-dependent settling, flocculation effects
    deposition_transport: Erosion, deposition, consolidation, mass exchange
    consolidation_gelling: Consolidation mechanics and gelling processes
    environmental_effects: Temperature, salinity, organic content effects

References:
    ASCE Manual 110, Chapter 4: Fine-Grained Sediment Transport
"""

from .sediment_characterization import *
from .rheological_models import *
from .floc_aggregation import *
from .settling_velocity import *
from .deposition_transport import *
from .consolidation_gelling import *
from .environmental_effects import *

__all__ = [
    # Sediment characterization
    'clay_cation_exchange_capacity',
    'bulk_density_concentration',
    'cohesion_classification',
    
    # Rheological models
    'sisko_viscosity_model',
    'concentration_viscosity_relation',
    'viscoelastic_standard_solid',
    'voigt_constitutive_model',
    
    # Floc aggregation
    'fractal_floc_density',
    'floc_strength_relation',
    'collision_frequency_function',
    'limiting_floc_size',
    
    # Settling velocity
    'multizone_settling_velocity',
    'concentration_dependent_settling',
    'temperature_settling_correction',
    'shear_rate_settling_effect',
    
    # Deposition and transport
    'deposition_rate_model',
    'multiclass_deposition',
    'critical_shear_deposition',
    'transport_flux_calculation',
    
    # Consolidation and gelling
    'consolidation_settling',
    'hindered_settling_dense',
    'effective_stress_consolidation',
    
    # Environmental effects
    'salinity_flocculation_effect',
    'temperature_viscosity_correction',
    'organic_content_influence',
]