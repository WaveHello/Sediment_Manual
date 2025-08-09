"""
Physical constants used in sediment transport calculations.

All constants are in SI units unless otherwise specified.
"""

# Fluid properties
WATER_DENSITY = 1000.0  # kg/m³ at 20°C
WATER_VISCOSITY_20C = 1.002e-3  # dynamic viscosity, kg/(m·s) at 20°C
WATER_KINEMATIC_VISCOSITY_20C = 1.002e-6  # kinematic viscosity, m²/s at 20°C

# Sediment properties  
QUARTZ_DENSITY = 2650.0  # kg/m³
LIMESTONE_DENSITY = 2710.0  # kg/m³
COAL_DENSITY = 1300.0  # kg/m³

# Physical constants
GRAVITY = 9.81  # m/s²
VON_KARMAN_CONSTANT = 0.41  # von Kármán constant

# Common particle sizes (in mm, converted to m)
CLAY_MAX_SIZE = 0.002e-3  # m
SILT_MAX_SIZE = 0.063e-3  # m  
SAND_MAX_SIZE = 2.0e-3    # m
GRAVEL_MAX_SIZE = 64.0e-3 # m

# Standard grain size classes (phi scale)
PHI_SCALE = {
    "clay": 8,      # φ > 8
    "silt": 4,      # 4 < φ ≤ 8  
    "sand": -1,     # -1 < φ ≤ 4
    "gravel": -6,   # -6 < φ ≤ -1
    "cobble": -8,   # -8 < φ ≤ -6
    "boulder": -10  # φ ≤ -8
}

# Temperature-dependent water properties (at 1 atm)
WATER_PROPERTIES = {
    # Temperature (°C): (density kg/m³, kinematic_viscosity m²/s)
    0: (999.8, 1.787e-6),
    10: (999.7, 1.307e-6), 
    20: (998.2, 1.004e-6),
    30: (995.7, 0.801e-6),
    40: (992.2, 0.658e-6)
}