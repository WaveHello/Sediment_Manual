"""
Deposition and Transport Processes for Fine-Grained Sediment (Chapter 4)

This module implements deposition rate equations and transport formulations for
cohesive sediment, including single-class and multi-class deposition models.
Based on equations 4-24 to 4-27 from ASCE Manual 110, Chapter 4.
"""

import numpy as np
from typing import Union, Optional, Tuple, List, Dict, Any
from ..utils.validators import validate_positive, validate_range, validate_array_like
from ..utils.constants import GRAVITY, WATER_DENSITY

def single_class_deposition_rate(
    concentration: Union[float, np.ndarray],
    settling_velocity: Union[float, np.ndarray],
    water_depth: float,
    bed_shear_stress: Union[float, np.ndarray],
    critical_stress_deposition: float
) -> Union[float, np.ndarray]:
    """
    Calculate deposition rate for single size class.
    
    Equation 4-24: dC/dt = -(ws*C/h) * (1 - τb/τd) for τb < τd
                    dC/dt = 0 for τb ≥ τd
    
    Parameters
    ----------
    concentration : float or array-like
        Depth-averaged suspended sediment concentration [kg/m³]
    settling_velocity : float or array-like
        Settling velocity [m/s]
    water_depth : float
        Water depth [m]
    bed_shear_stress : float or array-like
        Bed shear stress [Pa]
    critical_stress_deposition : float
        Critical stress for deposition [Pa]
        
    Returns
    -------
    float or ndarray
        Rate of concentration change due to deposition [kg/m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-24
    """
    validate_positive(concentration, "concentration")
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(water_depth, "water_depth")
    validate_positive(bed_shear_stress, "bed_shear_stress")
    validate_positive(critical_stress_deposition, "critical_stress_deposition")
    
    C = np.atleast_1d(concentration)
    ws = np.atleast_1d(settling_velocity)
    tau_b = np.atleast_1d(bed_shear_stress)
    tau_d = critical_stress_deposition
    
    # Deposition rate
    deposition_rate = np.zeros_like(C)
    
    # Deposition occurs only when τb < τd
    deposition_mask = tau_b < tau_d
    if np.any(deposition_mask):
        deposition_rate[deposition_mask] = -(ws * C * (1 - tau_b / tau_d) / water_depth)[deposition_mask]
    
    return deposition_rate.item() if (np.isscalar(concentration) and 
                                     np.isscalar(bed_shear_stress)) else deposition_rate


def exponential_concentration_decay(
    time: Union[float, np.ndarray],
    initial_concentration: float,
    settling_velocity: float,
    water_depth: float,
    bed_shear_stress: float,
    critical_stress_deposition: float
) -> Union[float, np.ndarray]:
    """
    Calculate exponential concentration decay during deposition.
    
    Equation 4-25: C/C0 = exp[-(ws/h) * (1 - τb/τd) * t] for τb < τd
    
    Parameters
    ----------
    time : float or array-like
        Time [s]
    initial_concentration : float
        Initial concentration C0 [kg/m³]
    settling_velocity : float
        Settling velocity [m/s]
    water_depth : float
        Water depth [m]
    bed_shear_stress : float
        Bed shear stress [Pa]
    critical_stress_deposition : float
        Critical stress for deposition [Pa]
        
    Returns
    -------
    float or ndarray
        Concentration at time t [kg/m³]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-25
    """
    validate_positive(time, "time")
    validate_positive(initial_concentration, "initial_concentration")
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(water_depth, "water_depth")
    validate_positive(bed_shear_stress, "bed_shear_stress")
    validate_positive(critical_stress_deposition, "critical_stress_deposition")
    
    t = np.atleast_1d(time)
    
    if bed_shear_stress >= critical_stress_deposition:
        # No deposition
        concentration = np.full_like(t, initial_concentration)
    else:
        # Exponential decay
        decay_rate = (settling_velocity / water_depth) * (1 - bed_shear_stress / critical_stress_deposition)
        concentration = initial_concentration * np.exp(-decay_rate * t)
    
    return concentration.item() if np.isscalar(time) else concentration


def multiclass_deposition_rate(
    concentrations: np.ndarray,
    settling_velocities: np.ndarray,
    water_depth: float,
    bed_shear_stress: Union[float, np.ndarray],
    critical_stresses: np.ndarray
) -> np.ndarray:
    """
    Calculate deposition rates for multiple size classes.
    
    Equation 4-26: dCi/dt = -(wsi*Ci/h) * (1 - τb/τdi) for τb < τdi
    
    Parameters
    ----------
    concentrations : array-like
        Concentrations for each size class [kg/m³]
    settling_velocities : array-like
        Settling velocities for each size class [m/s]
    water_depth : float
        Water depth [m]
    bed_shear_stress : float or array-like
        Bed shear stress [Pa]
    critical_stresses : array-like
        Critical stresses for deposition for each class [Pa]
        
    Returns
    -------
    ndarray
        Deposition rates for each size class [kg/m³/s]
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-26
    """
    C = validate_array_like(concentrations, "concentrations")
    ws = validate_array_like(settling_velocities, "settling_velocities")
    tau_d = validate_array_like(critical_stresses, "critical_stresses")
    
    if len(C) != len(ws) or len(C) != len(tau_d):
        raise ValueError("Arrays must have same length")
    
    validate_positive(water_depth, "water_depth")
    
    tau_b = np.atleast_1d(bed_shear_stress)
    n_classes = len(C)
    deposition_rates = np.zeros(n_classes)
    
    for i in range(n_classes):
        if tau_b < tau_d[i]:
            deposition_rates[i] = -(ws[i] * C[i] / water_depth) * (1 - tau_b / tau_d[i])
    
    return deposition_rates


def multiclass_concentration_evolution(
    time: Union[float, np.ndarray],
    initial_concentrations: np.ndarray,
    settling_velocities: np.ndarray,
    water_depth: float,
    bed_shear_stress: float,
    critical_stresses: np.ndarray,
    size_distribution: Optional[np.ndarray] = None
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Calculate concentration evolution for multiple size classes.
    
    Equation 4-27: Complex multiclass evolution with size distribution effects
    
    Parameters
    ----------
    time : float or array-like
        Time values [s]
    initial_concentrations : array-like
        Initial concentrations for each size class [kg/m³]
    settling_velocities : array-like
        Settling velocities for each size class [m/s]
    water_depth : float
        Water depth [m]
    bed_shear_stress : float
        Bed shear stress [Pa]
    critical_stresses : array-like
        Critical stresses for deposition [Pa]
    size_distribution : array-like, optional
        Frequency distribution for each size class
        
    Returns
    -------
    ndarray or Tuple[ndarray, ndarray]
        Concentration evolution for each class, optionally with total
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Equation 4-27
    Mehta and Lott (1987)
    """
    C0 = validate_array_like(initial_concentrations, "initial_concentrations")
    ws = validate_array_like(settling_velocities, "settling_velocities")
    tau_d = validate_array_like(critical_stresses, "critical_stresses")
    
    t = np.atleast_1d(time)
    n_classes = len(C0)
    n_times = len(t)
    
    if size_distribution is None:
        phi = np.ones(n_classes) / n_classes  # Uniform distribution
    else:
        phi = validate_array_like(size_distribution, "size_distribution")
        if len(phi) != n_classes:
            raise ValueError("Size distribution must match number of classes")
        phi = phi / np.sum(phi)  # Normalize
    
    # Calculate concentration evolution for each class
    concentrations = np.zeros((n_times, n_classes))
    
    for i in range(n_classes):
        if bed_shear_stress < critical_stresses[i]:
            # Class deposits
            decay_rate = (ws[i] / water_depth) * (1 - bed_shear_stress / critical_stresses[i])
            concentrations[:, i] = C0[i] * np.exp(-decay_rate * t)
        else:
            # Class remains in suspension
            concentrations[:, i] = np.full_like(t, C0[i])
    
    # Total concentration
    total_concentration = np.sum(concentrations, axis=1)
    
    if np.isscalar(time):
        return concentrations[0, :], total_concentration[0]
    else:
        return concentrations, total_concentration


def deposition_probability(
    particle_diameter: Union[float, np.ndarray],
    floc_strength: Union[float, np.ndarray],
    near_bed_shear_rate: Union[float, np.ndarray],
    collision_frequency: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Calculate probability of particle deposition based on floc properties.
    
    Deposition probability depends on:
    - Floc size and strength
    - Near-bed shearing conditions
    - Collision frequency with bed
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle/floc diameter [m]
    floc_strength : float or array-like
        Floc strength [Pa]
    near_bed_shear_rate : float or array-like
        Near-bed shear rate [s⁻¹]
    collision_frequency : float or array-like
        Collision frequency with bed [s⁻¹]
        
    Returns
    -------
    float or ndarray
        Deposition probability [0-1]
        
    References
    ----------
    Based on concepts from Stolzenbach et al. (1992)
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(floc_strength, "floc_strength")
    validate_positive(near_bed_shear_rate, "near_bed_shear_rate")
    validate_positive(collision_frequency, "collision_frequency")
    
    d = np.atleast_1d(particle_diameter)
    tau_f = np.atleast_1d(floc_strength)
    gamma = np.atleast_1d(near_bed_shear_rate)
    beta = np.atleast_1d(collision_frequency)
    
    # Simplified model: probability increases with floc strength and collision frequency
    # but decreases with shear rate
    p_stick = tau_f / (tau_f + 0.1 * gamma)  # Sticking probability
    p_collision = 1 - np.exp(-beta * 1e-3)  # Collision probability
    
    p_deposition = p_stick * p_collision
    p_deposition = np.clip(p_deposition, 0.0, 1.0)
    
    return p_deposition.item() if (np.isscalar(particle_diameter) and 
                                  np.isscalar(floc_strength)) else p_deposition


def critical_stress_for_deposition(
    particle_diameter: Union[float, np.ndarray],
    particle_density: float = 2650.0,
    water_density: float = WATER_DENSITY,
    empirical_coeff: float = 0.1
) -> Union[float, np.ndarray]:
    """
    Estimate critical shear stress for deposition.
    
    Simplified relationship based on particle properties.
    
    Parameters
    ----------
    particle_diameter : float or array-like
        Particle diameter [m]
    particle_density : float, optional
        Particle density [kg/m³], default 2650
    water_density : float, optional
        Water density [kg/m³], default 1000
    empirical_coeff : float, optional
        Empirical coefficient, default 0.1
        
    Returns
    -------
    float or ndarray
        Critical stress for deposition [Pa]
        
    References
    ----------
    Based on empirical observations from literature
    """
    validate_positive(particle_diameter, "particle_diameter")
    validate_positive(particle_density, "particle_density")
    
    d = np.atleast_1d(particle_diameter)
    
    # Smaller, lighter particles deposit more easily
    relative_density = (particle_density - water_density) / water_density
    tau_d = empirical_coeff * relative_density * GRAVITY * d
    
    return tau_d.item() if np.isscalar(particle_diameter) else tau_d


def sorting_during_deposition(
    initial_size_distribution: np.ndarray,
    transport_distances: np.ndarray,
    settling_velocities: np.ndarray,
    flow_velocity: float
) -> np.ndarray:
    """
    Calculate size sorting effects during transport and deposition.
    
    Parameters
    ----------
    initial_size_distribution : array-like
        Initial size distribution (fractions)
    transport_distances : array-like
        Transport distances [m]
    settling_velocities : array-like
        Settling velocities for each size class [m/s]
    flow_velocity : float
        Mean flow velocity [m/s]
        
    Returns
    -------
    ndarray
        Size distribution at each distance
        
    References
    ----------
    Based on concepts from Lin (1986), Lin and Mehta (1997)
    """
    f0 = validate_array_like(initial_size_distribution, "initial_size_distribution")
    x = validate_array_like(transport_distances, "transport_distances")
    ws = validate_array_like(settling_velocities, "settling_velocities")
    
    if len(f0) != len(ws):
        raise ValueError("Size distribution and settling velocities must have same length")
    
    validate_positive(flow_velocity, "flow_velocity")
    
    n_distances = len(x)
    n_classes = len(f0)
    size_distributions = np.zeros((n_distances, n_classes))
    
    # Travel time to each location
    travel_times = x / flow_velocity
    
    for i, t in enumerate(travel_times):
        # Fraction remaining in suspension for each size class
        remaining_fractions = np.exp(-ws * t / 1.0)  # Simplified model
        
        # Renormalize to get distribution
        total_remaining = np.sum(f0 * remaining_fractions)
        if total_remaining > 1e-10:
            size_distributions[i, :] = f0 * remaining_fractions / total_remaining
        else:
            size_distributions[i, :] = np.zeros(n_classes)
    
    return size_distributions


def entrainment_from_lutocline(
    richardson_number: Union[float, np.ndarray],
    entrainment_coeff: float = 2.5e-4
) -> Union[float, np.ndarray]:
    """
    Calculate entrainment rate from fluid mud layer (lutocline).
    
    Equation 4-39 (referenced): E = A * Ri^(-n)
    where E is nondimensional entrainment rate
    
    Parameters
    ----------
    richardson_number : float or array-like
        Richardson number characterizing density stratification
    entrainment_coeff : float, optional
        Entrainment coefficient, default 2.5e-4
        
    Returns
    -------
    float or ndarray
        Nondimensional entrainment rate
        
    References
    ----------
    ASCE Manual 110, Chapter 4, Figure 4-28
    Mehta and Srinivas (1993)
    """
    validate_positive(richardson_number, "richardson_number")
    validate_positive(entrainment_coeff, "entrainment_coeff")
    
    Ri = np.atleast_1d(richardson_number)
    
    # Power law relationship
    entrainment_rate = entrainment_coeff * Ri**(-1.0)
    
    return entrainment_rate.item() if np.isscalar(richardson_number) else entrainment_rate


def bed_exchange_continuous_model(
    concentration: Union[float, np.ndarray],
    settling_velocity: Union[float, np.ndarray],
    erosion_rate: Union[float, np.ndarray],
    water_depth: float
) -> Union[float, np.ndarray]:
    """
    Calculate net bed exchange assuming continuous erosion and deposition.
    
    Net exchange = Deposition flux - Erosion flux
    
    Parameters
    ----------
    concentration : float or array-like
        Suspended sediment concentration [kg/m³]
    settling_velocity : float or array-like
        Settling velocity [m/s]
    erosion_rate : float or array-like
        Erosion rate [kg/m²/s]
    water_depth : float
        Water depth [m]
        
    Returns
    -------
    float or ndarray
        Net exchange rate (positive = net deposition) [kg/m³/s]
        
    References
    ----------
    Alternative to mutually exclusive model
    """
    validate_positive(concentration, "concentration")
    validate_positive(settling_velocity, "settling_velocity")
    validate_positive(water_depth, "water_depth")
    
    C = np.atleast_1d(concentration)
    ws = np.atleast_1d(settling_velocity)
    E = np.atleast_1d(erosion_rate)
    
    # Deposition flux
    deposition_flux = ws * C
    
    # Net exchange rate per unit volume
    net_exchange = (deposition_flux - E) / water_depth
    
    return net_exchange.item() if (np.isscalar(concentration) and 
                                  np.isscalar(erosion_rate)) else net_exchange