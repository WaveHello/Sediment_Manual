"""
Chapter 3: Transport of Gravel and Sediment Mixtures

This module implements equations from ASCE Manual 110, Chapter 3, covering:
- Grain size distributions and statistical analysis
- Hydraulic relations for gravel-bed and sand-bed streams  
- Active layer sediment transport mechanics
- Hiding functions and threshold conditions
- Flow hydraulics and dimensionless parameters

References:
    ASCE Manual 110: Sedimentation Engineering: Processes, Measurements, 
    Modeling, and Practice (2008)
"""

from .grain_size_distributions import *
from .hydraulic_relations import *
from .active_layer_mechanics import *
from .hiding_functions import *
from .flow_hydraulics import *
from .transport_mechanics import *

__all__ = [
    # Grain size distributions
    'psi_scale_transform', 'grain_size_statistics', 'lognormal_distribution',
    'discretized_grain_size_stats',
    
    # Hydraulic relations  
    'bankfull_velocity', 'bankfull_shear_stress', 'dimensionless_parameters',
    'channel_scaling_relations',
    
    # Active layer mechanics
    'active_layer_conservation', 'bed_elevation_change', 'exchange_fractions',
    'entrainment_deposition_rates',
    
    # Hiding functions
    'egiazaroff_hiding_function', 'power_law_hiding_function', 
    'shields_threshold_relations', 'substrate_hiding_functions',
    
    # Flow hydraulics
    'shear_velocity', 'friction_coefficients', 'one_dimensional_flow',
    
    # Transport mechanics
    'ashida_michiue_transport', 'parker_surface_based_transport', 
    'wilcock_crowe_multisize', 'generic_bedload_transport', 
    'complexity_corrected_transport', 'equal_mobility_condition',
    'static_armor_predictor', 'shields_stress_calculation'
]