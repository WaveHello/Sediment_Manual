"""
Chapter 8: River Meandering and Channel Stability

This module implements mathematical models for river meandering processes,
channel stability analysis, and flow dynamics in curved channels based on
ASCE Manual 110, Chapter 8.

Modules:
    meandering_criteria: Stability thresholds and meandering criteria (Equations 8-1 to 8-4)
    planform_geometry: Meander planform relationships (Equations 8-5 to 8-7) 
    bank_erosion: Bank erosion and migration rates (Equations 8-8 to 8-12)
    flow_dynamics: Flow and bed topography in meanders (Equations 8-13 to 8-28)
    stability_analysis: Perturbation stability analysis (Equations 8-29 to 8-44)
    channel_migration: Channel migration simulation (Equations 8-45 to 8-48)
    channel_stabilization: Engineering applications and stabilization techniques

References:
    ASCE Manual 110: Sedimentation Engineering - Chapter 8
    "River Meandering and Channel Stability" by A. Jacob Odgaard and Jorge D. Abad
"""

from .meandering_criteria import *
from .planform_geometry import *
from .bank_erosion import *
from .flow_dynamics import *
from .stability_analysis import *
from .channel_migration import *
from .channel_stabilization import *

__version__ = "1.0.0"
__author__ = "ASCE Manual 110 - Chapter 8"

# Module metadata
CHAPTER_TITLE = "River Meandering and Channel Stability"
CHAPTER_NUMBER = 8
EQUATION_COUNT = 48
MODULE_COUNT = 7

# Key physical constants for meandering analysis
VON_KARMAN_CONSTANT = 0.4  # von Kármán constant
GRAVITY = 9.81  # m/s²
WATER_DENSITY = 1000.0  # kg/m³
SEDIMENT_DENSITY = 2650.0  # kg/m³ (quartz)
DENSIMETRIC_RATIO = 1.65  # (ρs - ρ)/ρ

# Typical parameter ranges
TYPICAL_RANGES = {
    'friction_parameter_m': (3.0, 5.0),
    'transverse_bed_slope_factor_B': (3.0, 6.0),
    'transverse_mass_flux_factor_alpha': (0.3, 0.5),
    'critical_shields_parameter': (0.03, 0.06),
    'sediment_transport_exponent_M': (2.0, 4.0),
    'bank_friction_angle': (30.0, 45.0),  # degrees
    'width_depth_ratio': (10.0, 60.0),
    'particle_froude_number': (5.0, 15.0)
}