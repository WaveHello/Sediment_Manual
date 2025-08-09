"""
Hiding Functions and Threshold Conditions (Chapter 3, Equations 3-60 to 3-74)

This module implements hiding functions that account for the effects of grain size
distributions on sediment entrainment thresholds in mixed-size sediments.
"""

import numpy as np
from typing import Union, Literal
from ..utils.validators import validate_positive, validate_array


def egiazaroff_hiding_function(D_i: Union[float, np.ndarray], 
                              D_g: float,
                              return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Egiazaroff hiding function for threshold conditions.
    
    Equations 3-72a,b: Original Egiazaroff formulation
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_g: Geometric mean grain size of surface (mm)
        return_reduced: If True, return reduced form; if False, standard form
        
    Returns:
        Hiding function value(s)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-72a,b
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_g, "geometric_mean_size")
    
    # Size ratio
    size_ratio = D_i / D_g
    
    # Egiazaroff hiding function
    log_19 = np.log(19.0)
    log_term = np.log(19.0 * size_ratio)
    
    # Avoid numerical issues with log
    log_term = np.maximum(log_term, 1e-6)
    
    if return_reduced:
        # Equation 3-72b: τ*_bsci/τ*_bscg = (D_i/D_g)[log(19)/log(19D_i/D_g)]²
        hiding_function = size_ratio * (log_19 / log_term)**2
    else:
        # Equation 3-72a: τ*_sci/τ*_scg = [log(19)/log(19D_i/D_g)]²
        hiding_function = (log_19 / log_term)**2
    
    return hiding_function


def modified_egiazaroff_hiding_function(D_i: Union[float, np.ndarray], 
                                       D_g: float,
                                       return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Modified Egiazaroff hiding function (Ashida and Michiue).
    
    Equations 3-73a,b: Modified form for fine particles
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_g: Geometric mean grain size of surface (mm)
        return_reduced: If True, return reduced form; if False, standard form
        
    Returns:
        Modified hiding function value(s)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-73a,b
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_g, "geometric_mean_size")
    
    size_ratio = D_i / D_g
    
    # For fine particles (D_i/D_g ≤ 0.4), use constant value
    if np.isscalar(size_ratio):
        if size_ratio <= 0.4:
            return 0.843 if not return_reduced else 0.843 * size_ratio
        else:
            # Use regular Egiazaroff for coarser particles
            return egiazaroff_hiding_function(D_i, D_g, return_reduced)
    else:
        # Array case
        hiding_function = np.where(
            size_ratio <= 0.4,
            0.843 if not return_reduced else 0.843 * size_ratio,
            egiazaroff_hiding_function(D_i, D_g, return_reduced)
        )
        return hiding_function


def power_law_hiding_function(D_i: Union[float, np.ndarray], 
                             D50: float, 
                             gamma: float = -0.982,
                             return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Power law hiding function.
    
    Equations 3-74a,b: Power relation hiding functions
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D50: Median grain size (mm)
        gamma: Power law exponent, default -0.982
        return_reduced: If True, return reduced form; if False, standard form
        
    Returns:
        Power law hiding function value(s)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-74a,b
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D50, "median_grain_size")
    
    size_ratio = D_i / D50
    
    if return_reduced:
        # Equation 3-74b: τ*_bssri/τ*_bssr50 = (D_i/D50)^(γ-1)
        hiding_function = size_ratio**(gamma - 1.0)
    else:
        # Equation 3-74a: τ*_ssri/τ*_ssr50 = (D_i/D50)^γ
        hiding_function = size_ratio**gamma
    
    return hiding_function


def reduced_hiding_function(D_i: Union[float, np.ndarray],
                           D_g: float,
                           F_hc: Union[float, np.ndarray],
                           function_type: Literal['critical', 'reference'] = 'critical') -> Union[float, np.ndarray]:
    """
    General reduced hiding function formulation.
    
    Equations 3-60a,b: τ*_bsci/τ*_bscg = F_hc(D_i/D_g) × (D_i/D_g)
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_g: Geometric mean grain size (mm)
        F_hc: Hiding function values
        function_type: Either 'critical' or 'reference'
        
    Returns:
        Reduced hiding function values
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-60a,b
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_g, "geometric_mean_size")
    
    size_ratio = D_i / D_g
    
    # Equation 3-60: Reduced hiding function
    reduced_function = F_hc * size_ratio
    
    return reduced_function


def shields_threshold_relations(R_pg: float) -> float:
    """
    Calculate Shields threshold stress using modified Brownlie fit.
    
    Equation 3-71: τ*_scg = (1/2)[0.22R_pg^-0.6 + 0.06×10^(-7.7R_pg^-0.6)]
    
    Args:
        R_pg: Particle Reynolds number
        
    Returns:
        Shields threshold stress (dimensionless)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-71
    """
    validate_positive(R_pg, "particle_reynolds_number")
    
    # Modified Brownlie fit coefficients
    term1 = 0.22 * (R_pg**(-0.6))
    term2 = 0.06 * (10**(-7.7 * R_pg**(-0.6)))
    
    tau_star_scg = 0.5 * (term1 + term2)
    
    return tau_star_scg


def size_independence_limiting_case(D_i: Union[float, np.ndarray],
                                   D_g: float,
                                   R_pg: float,
                                   return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Size-independence limiting case (F_hc = F_hr = 1).
    
    Equations 3-61a-d: Size-independence case
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_g: Geometric mean grain size (mm)
        R_pg: Particle Reynolds number
        return_reduced: If True, return reduced form
        
    Returns:
        Size-independent hiding function values
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-61a-d
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_g, "geometric_mean_size")
    validate_positive(R_pg, "particle_reynolds_number")
    
    size_ratio = D_i / D_g
    
    if return_reduced:
        # Equations 3-61c,d: τ*_bsci/τ*_bscg = D_i/D_g
        return size_ratio
    else:
        # Equations 3-61a,b: τ*_sci/τ*_scg = (R_pg)^-1
        return np.full_like(size_ratio, R_pg**(-1))


def equal_threshold_limiting_case(D_i: Union[float, np.ndarray],
                                 D_g: float,
                                 R_pg: float,
                                 return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Equal-threshold limiting case (F_hc = F_hr = (D_i/D_g)^-1).
    
    Equations 3-62a-d: Equal-threshold case
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_g: Geometric mean grain size (mm)
        R_pg: Particle Reynolds number
        return_reduced: If True, return reduced form
        
    Returns:
        Equal-threshold hiding function values
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-62a-d
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_g, "geometric_mean_size")
    validate_positive(R_pg, "particle_reynolds_number")
    
    size_ratio = D_i / D_g
    
    if return_reduced:
        # Equations 3-62c,d: τ*_bsci/τ*_bscg = 1
        return np.ones_like(size_ratio)
    else:
        # Equations 3-62a,b: τ*_sci/τ*_scg = (D_i/D_g) × (R_pg)^-1
        return size_ratio * (R_pg**(-1))


def substrate_hiding_functions(D_i: Union[float, np.ndarray],
                              D_ug: float,
                              R_pug: float,
                              F_uhc: Union[float, np.ndarray],
                              return_reduced: bool = False) -> Union[float, np.ndarray]:
    """
    Substrate-based hiding functions.
    
    Equations 3-67a,b, 3-68a,b: Substrate hiding functions
    
    Args:
        D_i: Grain size(s) for size class i (mm)
        D_ug: Geometric mean grain size of substrate (mm)
        R_pug: Particle Reynolds number for substrate
        F_uhc: Substrate hiding function values
        return_reduced: If True, return reduced form
        
    Returns:
        Substrate hiding function values
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-67, 3-68
    """
    validate_positive(D_i, "grain_size")
    validate_positive(D_ug, "substrate_geometric_mean")
    validate_positive(R_pug, "substrate_reynolds_number")
    
    size_ratio = D_i / D_ug
    
    if return_reduced:
        # Equations 3-68a,b: τ*_bsuci/τ*_bsucg = F_uhc(D_i/D_ug) × (D_i/D_ug)
        return F_uhc * size_ratio
    else:
        # Equations 3-67a,b: τ*_suci/τ*_sucg = F_uhc(D_i/D_ug) × (R_pug)^-1
        return F_uhc * (R_pug**(-1))