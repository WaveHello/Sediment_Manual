"""
SandWave: A Python library for sediment transport calculations.

This library provides implementations of equations from the ASCE Manual 110:
"Sedimentation Engineering: Processes, Measurements, Modeling, and Practice"
"""

__version__ = "0.1.0"
__author__ = "Generated from ASCE Manual 110"

# Import main modules for easy access
from . import chapter02_transport_morphodynamics
from . import utils

__all__ = [
    "chapter02_transport_morphodynamics",
    "utils",
]