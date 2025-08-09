"""
Consolidation and Gelling Processes for Fine-Grained Sediment (Chapter 4)

This module implements consolidation mechanics, settling behavior, and gelling
processes for cohesive sediment. Based on equations 4-28, 4-29 and consolidation
theory from ASCE Manual 110, Chapter 4, Section 4.7.
"""

import numpy as np
from typing import Union, Optional, Tuple, Dict, Any
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

def sediment_continuity_consolidation(
    layer_thickness: float,
    consolidation_velocity: Union[float, np.ndarray],
    concentration: Union[float, np.ndarray],
    settling_flux: Union[float, np.ndarray] = 0.0
) -> Union[float, np.ndarray]:
    """
    Calculate sediment continuity during consolidation.
    
    Equation 4-28: ∂h'C/∂t = ∂wsc C/∂z' + q
    
    Where:
    - h' = thickness of consolidating layer
    - wsc = consolidation velocity  
    - z' = z/h' (normalized coordinate)
    - q = net settling flux at top of layer
    
    Parameters
    ----------
    layer_thickness : float
        Thickness of consolidating layer h' [m]
    consolidation_velocity : float or array-like
        Consolidation velocity wsc [m/s]
    concentration : float or array-like
        Sediment concentration C [kg/m³]
    settling_flux : float or array-like, optional
        Net settling flux q [kg/m²/s], default 0.0
        
    Returns
    -------
    float or ndarray
        Rate of change of sediment content [kg/m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-28
    Jiang (1999)
    """
    validate_positive(layer_thickness, "layer_thickness")
    validate_positive(consolidation_velocity, "consolidation_velocity")
    validate_positive(concentration, "concentration")
    
    h_prime = layer_thickness
    wsc = np.atleast_1d(consolidation_velocity)
    C = np.atleast_1d(concentration)
    q = np.atleast_1d(settling_flux)
    
    # Simplified form assuming uniform properties
    # Full solution requires solving the continuity equation
    rate_change = (wsc * C) / h_prime + q / h_prime
    
    return rate_change.item() if (np.isscalar(consolidation_velocity) and 
                                 np.isscalar(concentration)) else rate_change


def two_mode_consolidation_velocity(
    concentration: Union[float, np.ndarray],
    loose_soil_velocity: float,
    loose_soil_concentration: float,
    compact_soil_velocity: float,
    compact_soil_concentration: float,
    transition_concentration: float,
    transition_parameter: float = 5.0,
    transition_exponent: float = 10.0
) -> Union[float, np.ndarray]:
    """
    Calculate consolidation velocity using two-mode model.
    
    Equation 4-29: Two-mode consolidation model
    wsc = wsc1 * F1 + wsc2 * F2 * [1 - Ft]
    
    Where Ft is the mode transition function
    
    Parameters
    ----------
    concentration : float or array-like
        Sediment concentration [kg/m³]
    loose_soil_velocity : float
        Consolidation velocity for loose soil mode [m/s]
    loose_soil_concentration : float
        Reference concentration for loose soil [kg/m³]
    compact_soil_velocity : float
        Consolidation velocity for compact soil mode [m/s]
    compact_soil_concentration : float
        Saturation concentration for compact soil [kg/m³]
    transition_concentration : float
        Mode transition concentration [kg/m³]
    transition_parameter : float, optional
        Transition sharpness parameter mt, default 5.0
    transition_exponent : float, optional
        Transition exponent nt, default 10.0
        
    Returns
    -------
    float or ndarray
        Consolidation velocity [m/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-29
    Toorman and Berlamont (1993)
    """
    validate_positive(concentration, "concentration")
    validate_positive(loose_soil_velocity, "loose_soil_velocity")
    validate_positive(compact_soil_velocity, "compact_soil_velocity")
    validate_positive(loose_soil_concentration, "loose_soil_concentration")
    validate_positive(compact_soil_concentration, "compact_soil_concentration")
    validate_positive(transition_concentration, "transition_concentration")
    
    C = np.atleast_1d(concentration)
    wsc1 = loose_soil_velocity
    wsc2 = compact_soil_velocity
    Cs1 = loose_soil_concentration
    Cs2 = compact_soil_concentration
    Ct = transition_concentration
    mt = transition_parameter
    nt = transition_exponent
    
    # Mode transition function
    Ft = np.exp(-mt * (C / Ct)**nt)
    
    # Mode 1 function (loose soil)
    F1 = C / Cs1
    
    # Mode 2 function (compact soil)  
    F2 = (1 - C / Cs2)**(-1)
    
    # Combined consolidation velocity
    wsc = wsc1 * F1 * Ft + wsc2 * F2 * (1 - Ft)
    
    return wsc.item() if np.isscalar(concentration) else wsc


def consolidation_time_scales(
    bed_thickness: Union[float, np.ndarray],
    initial_phase: bool = True
) -> Union[float, np.ndarray]:
    """
    Calculate characteristic time scales for consolidation.
    
    Initial consolidation: t ∝ 1/h (inversely proportional to thickness)
    Later consolidation: t ∝ 1/h² (inversely proportional to thickness squared)
    
    Parameters
    ----------
    bed_thickness : float or array-like
        Bed thickness [m]
    initial_phase : bool, optional
        True for initial rapid consolidation, False for later slow phase
        
    Returns
    -------
    float or ndarray
        Characteristic consolidation time [s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.7.1
    Toorman (1996), Winterwerp (1999)
    """
    validate_positive(bed_thickness, "bed_thickness")
    
    h = np.atleast_1d(bed_thickness)
    
    # Characteristic parameters (order of magnitude estimates)
    if initial_phase:
        # Initial consolidation: t ~ h/wsc
        characteristic_velocity = 1e-6  # m/s
        time_scale = h / characteristic_velocity
    else:
        # Later consolidation: t ~ h²/D where D is consolidation coefficient
        consolidation_coefficient = 1e-8  # m²/s
        time_scale = h**2 / consolidation_coefficient
    
    return time_scale.item() if np.isscalar(bed_thickness) else time_scale


def flocculation_settling_consolidation_stages(
    time: Union[float, np.ndarray],
    slurry_type: str = 'kaolinite'
) -> Union[str, np.ndarray]:
    """
    Identify consolidation stages based on time progression.
    
    Three stages:
    1. Flocculation stage (seconds to minutes)
    2. Primary consolidation with settling (minutes to hours)  
    3. Secondary consolidation (hours to days/months)
    
    Parameters
    ----------
    time : float or array-like
        Time since start of process [s]
    slurry_type : str, optional
        Type of slurry ('kaolinite', 'bentonite'), default 'kaolinite'
        
    Returns
    -------
    str or ndarray
        Consolidation stage identifier
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Figure 4-16
    Schiffman et al. (1985)
    """
    validate_positive(time, "time")
    
    t = np.atleast_1d(time)
    stages = np.empty(t.shape, dtype=object)
    
    if slurry_type.lower() == 'kaolinite':
        # Fast-flocculating material
        flocculation_time = 60  # 1 minute
        primary_time = 3600  # 1 hour
    elif slurry_type.lower() == 'bentonite':
        # Slow-flocculating material
        flocculation_time = 1800  # 30 minutes  
        primary_time = 24 * 3600  # 24 hours
    else:
        # Default values
        flocculation_time = 300  # 5 minutes
        primary_time = 7200  # 2 hours
    
    stages[t <= flocculation_time] = "flocculation"
    stages[(t > flocculation_time) & (t <= primary_time)] = "primary_consolidation"
    stages[t > primary_time] = "secondary_consolidation"
    
    return stages.item() if np.isscalar(time) else stages


def effective_stress_profile(
    depth: Union[float, np.ndarray],
    bulk_density: Union[float, np.ndarray],
    water_density: float = WATER_DENSITY,
    lutocline_depth: float = 0.0
) -> Union[float, np.ndarray]:
    """
    Calculate effective stress profile in consolidating bed.
    
    Effective stress = Total stress - Pore pressure
    Above lutocline: σ'e ≈ 0 (no particle contact)
    Below lutocline: σ'e increases with depth due to overburden
    
    Parameters
    ----------
    depth : float or array-like
        Depth below lutocline [m]
    bulk_density : float or array-like
        Bulk density of sediment-water mixture [kg/m³]
    water_density : float, optional
        Water density [kg/m³], default 1000
    lutocline_depth : float, optional
        Depth of lutocline from surface [m], default 0.0
        
    Returns
    -------
    float or ndarray
        Effective stress [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Figure 4-17
    Sills and Elder (1986)
    """
    validate_positive(depth, "depth")
    validate_positive(bulk_density, "bulk_density")
    
    z = np.atleast_1d(depth)
    rho = np.atleast_1d(bulk_density)
    
    # Effective stress due to particle-supported load
    # σ'e = (ρ - ρw) * g * z for z > lutocline
    effective_stress = np.maximum(0, (rho - water_density) * GRAVITY * z)
    
    return effective_stress.item() if (np.isscalar(depth) and 
                                     np.isscalar(bulk_density)) else effective_stress


def gelling_process_parameters(
    clay_type: str,
    temperature: float = 20.0
) -> Dict[str, float]:
    """
    Get parameters for thixotropic gelling process.
    
    Gelling characteristics vary by clay type and temperature.
    
    Parameters
    ----------
    clay_type : str
        Type of clay ('kaolinite', 'bentonite', 'illite')
    temperature : float, optional
        Water temperature [°C], default 20.0
        
    Returns
    -------
    Dict[str, float]
        Dictionary of gelling parameters
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.7.2
    Day and Ripple (1966), Mitchell (1993)
    """
    validate_range(temperature, "temperature", 0.0, 50.0)
    
    # Base parameters for different clay types
    base_params = {
        'kaolinite': {
            'peak_strength_ratio': 5.0,  # SP/SR
            'gelling_time': 86400,  # 24 hours in seconds
            'remold_strength_ratio': 0.2
        },
        'bentonite': {
            'peak_strength_ratio': 15.0,
            'gelling_time': 259200,  # 72 hours
            'remold_strength_ratio': 0.05
        },
        'illite': {
            'peak_strength_ratio': 8.0,
            'gelling_time': 172800,  # 48 hours  
            'remold_strength_ratio': 0.1
        }
    }
    
    clay_type = clay_type.lower()
    if clay_type not in base_params:
        raise ValueError(f"Unknown clay type: {clay_type}")
    
    params = base_params[clay_type].copy()
    
    # Temperature correction (higher temperature reduces gelling time)
    temp_factor = np.exp(-0.02 * (temperature - 20))  # Empirical
    params['gelling_time'] *= temp_factor
    
    # Sensitivity
    params['sensitivity'] = params['peak_strength_ratio'] / params['remold_strength_ratio']
    
    return params


def strength_recovery_with_time(
    time: Union[float, np.ndarray],
    remolded_strength: float,
    peak_strength: float,
    recovery_time_constant: float = 86400.0
) -> Union[float, np.ndarray]:
    """
    Calculate strength recovery during gelling process.
    
    Exponential recovery model: S(t) = SR + (SP - SR) * [1 - exp(-t/tc)]
    
    Parameters
    ----------
    time : float or array-like
        Time since remolding [s]
    remolded_strength : float
        Initial remolded strength SR [Pa]
    peak_strength : float
        Peak undisturbed strength SP [Pa]
    recovery_time_constant : float, optional
        Time constant for recovery [s], default 86400 (24 hours)
        
    Returns
    -------
    float or ndarray
        Strength at time t [Pa]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Figure 4-19
    Mitchell (1976)
    """
    validate_positive(time, "time")
    validate_positive(remolded_strength, "remolded_strength")
    validate_positive(peak_strength, "peak_strength")
    validate_positive(recovery_time_constant, "recovery_time_constant")
    
    if peak_strength <= remolded_strength:
        raise ValueError("Peak strength must be greater than remolded strength")
    
    t = np.atleast_1d(time)
    tc = recovery_time_constant
    
    # Exponential recovery
    strength = remolded_strength + (peak_strength - remolded_strength) * (1 - np.exp(-t / tc))
    
    return strength.item() if np.isscalar(time) else strength


def pore_water_pressure_evolution(
    time: Union[float, np.ndarray],
    initial_pressure: float,
    consolidation_coefficient: float = 1e-7,
    drainage_path: float = 0.1
) -> Union[float, np.ndarray]:
    """
    Calculate pore water pressure during consolidation.
    
    Simplified consolidation theory: exponential decay of excess pore pressure
    
    Parameters
    ----------
    time : float or array-like
        Time [s]
    initial_pressure : float
        Initial excess pore pressure [Pa]
    consolidation_coefficient : float, optional
        Consolidation coefficient [m²/s], default 1e-7
    drainage_path : float, optional
        Drainage path length [m], default 0.1
        
    Returns
    -------
    float or ndarray
        Pore water pressure [Pa]
        
    References
    ----------
    Based on consolidation theory, Terzaghi (1943)
    """
    validate_positive(time, "time")
    validate_positive(initial_pressure, "initial_pressure")
    validate_positive(consolidation_coefficient, "consolidation_coefficient")
    validate_positive(drainage_path, "drainage_path")
    
    t = np.atleast_1d(time)
    cv = consolidation_coefficient
    H = drainage_path
    
    # Time factor
    Tv = cv * t / H**2
    
    # Simplified solution (first term of Fourier series)
    pressure = initial_pressure * np.exp(-np.pi**2 * Tv / 4)
    
    return pressure.item() if np.isscalar(time) else pressure


def aggregate_order_crushing(
    overburden_pressure: Union[float, np.ndarray],
    initial_order: int = 0,
    crushing_threshold: float = 1000.0  # Pa
) -> Union[int, np.ndarray]:
    """
    Determine aggregate order changes due to overburden crushing.
    
    Based on Krone's concept: higher order aggregates form initially,
    then are crushed to lower orders by overburden pressure.
    
    Parameters
    ----------
    overburden_pressure : float or array-like
        Overburden pressure [Pa]
    initial_order : int, optional
        Initial aggregate order, default 0
    crushing_threshold : float, optional
        Pressure threshold for crushing [Pa], default 1000
        
    Returns
    -------
    int or ndarray
        Final aggregate order after crushing
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Section 4.7.1
    Krone (1963)
    """
    validate_positive(overburden_pressure, "overburden_pressure")
    
    P = np.atleast_1d(overburden_pressure)
    
    # Number of crushing events based on pressure level
    crushing_events = np.floor(P / crushing_threshold).astype(int)
    
    # Final order (cannot go below 0)
    final_order = np.maximum(0, initial_order - crushing_events)
    
    return final_order.item() if np.isscalar(overburden_pressure) else final_order


def consolidation_settlement(
    time: Union[float, np.ndarray],
    initial_thickness: float,
    final_thickness: float,
    consolidation_time: float
) -> Union[float, np.ndarray]:
    """
    Calculate settlement during consolidation process.
    
    Exponential approach to final settlement
    
    Parameters
    ----------
    time : float or array-like
        Time [s]
    initial_thickness : float
        Initial layer thickness [m]
    final_thickness : float
        Final consolidated thickness [m]  
    consolidation_time : float
        Characteristic consolidation time [s]
        
    Returns
    -------
    float or ndarray
        Layer thickness at time t [m]
        
    References
    ----------
    Based on consolidation settlement theory
    """
    validate_positive(time, "time")
    validate_positive(initial_thickness, "initial_thickness")
    validate_positive(final_thickness, "final_thickness")
    validate_positive(consolidation_time, "consolidation_time")
    
    if final_thickness >= initial_thickness:
        raise ValueError("Final thickness must be less than initial thickness")
    
    t = np.atleast_1d(time)
    h0 = initial_thickness
    hf = final_thickness
    tc = consolidation_time
    
    # Exponential settlement
    thickness = hf + (h0 - hf) * np.exp(-t / tc)
    
    return thickness.item() if np.isscalar(time) else thickness