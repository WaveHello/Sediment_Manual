"""
Chapter 5: Sediment Transport Measurements

This package implements statistical methods, sampling procedures, and measurement
techniques for sediment transport data collection and analysis.

Modules:
    statistical_methods: Sample size calculations, bias corrections, confidence intervals
    sampling_equipment: Equipment specifications, calibration functions
    bed_material_sampling: Grid sampling, volumetric methods, bias corrections
    suspended_sediment_sampling: Isokinetic samplers, collection procedures
    sample_analysis: Sieve analysis, concentration calculations
    field_procedures: Protocols, QA/QC methods, data validation
    measurement_conversions: Unit conversions, standardization functions

References:
    ASCE Manual 110, Chapter 5: Sediment Transport Measurements
"""

from .statistical_methods import *
from .sampling_equipment import *
from .bed_material_sampling import *
from .suspended_sediment_sampling import *
from .sample_analysis import *
from .field_procedures import *
from .measurement_conversions import *

__all__ = [
    # Statistical methods
    'areal_to_volumetric_conversion',
    'sample_size_binomial',
    'multinomial_confidence_interval',
    'sampling_bias_correction',
    'bootstrap_confidence_interval',
    
    # Sampling equipment
    'sampler_specifications',
    'isokinetic_sampler_efficiency',
    'equipment_calibration_curve',
    'sampler_selection_criteria',
    
    # Bed material sampling
    'wolman_walk_sampling',
    'grid_sampling_design',
    'volumetric_sampling_protocol',
    'photographic_sampling_analysis',
    'surface_oriented_sampling',
    
    # Suspended sediment sampling
    'depth_integrated_sampling',
    'point_sampling_protocol',
    'discharge_weighted_sampling',
    'automatic_sampler_programming',
    'isokinetic_sampling_velocity',
    
    # Sample analysis
    'sieve_analysis_processing',
    'particle_size_distribution',
    'sediment_concentration_calculation',
    'quality_control_analysis',
    'grain_size_statistics',
    
    # Field procedures
    'sampling_site_selection',
    'data_quality_assessment',
    'field_measurement_protocol',
    'equipment_maintenance_schedule',
    'measurement_uncertainty_analysis',
    
    # Measurement conversions
    'concentration_unit_conversion',
    'discharge_unit_conversion',
    'particle_size_unit_conversion',
    'temperature_correction_factor',
    'standard_reference_conditions',
]