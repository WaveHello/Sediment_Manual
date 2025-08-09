"""
Engineering Applications for Mixed Sediment Transport (Chapter 3, Equations 3-116 to 3-161)

This module implements advanced engineering applications including 2D transport,
subsidence effects, suspended sediment interactions, and overbank deposition
for practical sediment management problems.
"""

import numpy as np
from typing import Union, Tuple, Dict, Optional
from ..utils.validators import validate_positive, validate_array
from ..utils.constants import GRAVITY, WATER_DENSITY, VON_KARMAN


def subsidence_modified_exner(q_i: np.ndarray, 
                             f_i: np.ndarray,
                             F_i: np.ndarray,
                             eta_b: float,
                             L_a: float,
                             sigma_sub: float,
                             A_i: Optional[np.ndarray] = None,
                             lambda_p: float = 0.4,
                             dx: float = 1.0,
                             dt: float = 1.0) -> Tuple[float, np.ndarray]:
    """
    Modified Exner equation accounting for subsidence effects.
    
    Equation 3-116: (1-λ_p)[f_i(∂η_b/∂t + σ_sub) + ∂/∂t(L_a F_i)] = -∂q_i/∂s - A_i
    
    Args:
        q_i: Bed-load transport rates by size class (m²/s)
        f_i: Substrate size fractions (dimensionless)
        F_i: Active layer size fractions (dimensionless)
        eta_b: Bed elevation (m)
        L_a: Active layer thickness (m)
        sigma_sub: Subsidence rate (m/year), negative for uplift
        A_i: Abrasion loss rates by size class (m/s), optional
        lambda_p: Porosity (dimensionless), default 0.4
        dx: Spatial grid spacing (m), default 1.0
        dt: Time step (s), default 1.0
        
    Returns:
        Tuple of (bed_elevation_change_rate, layer_fraction_change_rates)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-116
    """
    validate_array(q_i)
    validate_array(f_i)
    validate_array(F_i)
    validate_positive(L_a, "active_layer_thickness")
    
    if A_i is None:
        A_i = np.zeros_like(q_i)
    
    # Convert subsidence rate from m/year to m/s
    sigma_sub_per_sec = sigma_sub / (365.25 * 24 * 3600)
    
    # Spatial gradient of transport rates
    dq_dx = np.gradient(q_i, dx)
    
    # Bed elevation change rate (including subsidence)
    deta_b_dt = (-np.sum(dq_dx) - np.sum(A_i)) / ((1.0 - lambda_p) * np.sum(f_i)) - sigma_sub_per_sec
    
    # Active layer fraction change rates
    dF_i_dt = np.zeros_like(F_i)
    for i in range(len(F_i)):
        # Rate of change for each size fraction
        numerator = -dq_dx[i] - A_i[i] - f_i[i] * (deta_b_dt + sigma_sub_per_sec)
        denominator = (1.0 - lambda_p) * L_a
        dF_i_dt[i] = numerator / denominator
    
    return deta_b_dt, dF_i_dt


def vectorial_bedload_transport_2d(tau_bs_s: Union[float, np.ndarray],
                                  tau_bs_n: Union[float, np.ndarray],
                                  tau_sci_star: Union[float, np.ndarray],
                                  D_i: Union[float, np.ndarray],
                                  grad_eta_s: Union[float, np.ndarray],
                                  grad_eta_n: Union[float, np.ndarray],
                                  beta: float = 2.0,
                                  rho_s: float = 2650.0,
                                  rho: float = WATER_DENSITY) -> Tuple[np.ndarray, np.ndarray]:
    """
    2D vectorial bed-load transport for river bends and planform sorting.
    
    Equations 3-117 to 3-122: 2D bed-load transport with bed slope effects
    
    Args:
        tau_bs_s: Streamwise boundary shear stress (Pa)
        tau_bs_n: Normal boundary shear stress (Pa)
        tau_sci_star: Critical Shields stress for size i (dimensionless)
        D_i: Grain sizes (m)
        grad_eta_s: Streamwise bed slope (dimensionless)
        grad_eta_n: Normal bed slope (dimensionless)
        beta: Slope effect parameter, default 2.0
        rho_s: Sediment density (kg/m³), default 2650
        rho: Water density (kg/m³), default 1000
        
    Returns:
        Tuple of (q_i_streamwise, q_i_normal) transport rates (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-117 to 3-122
    """
    validate_positive(D_i, "grain_sizes")
    
    # Submerged specific gravity
    R = (rho_s / rho) - 1.0
    
    # Shields stresses
    tau_si_star = tau_bs_s / (rho * R * GRAVITY * D_i)
    tau_sn_star = tau_bs_n / (rho * R * GRAVITY * D_i)
    
    # Streamwise transport (standard relation)
    excess_stress_s = np.maximum(tau_si_star - tau_sci_star, 0.0)
    q_i_s = 17.0 * np.sqrt(R * GRAVITY * D_i**3) * excess_stress_s**(3/2)
    
    # Normal transport with bed slope effects (Ikeda-Parker linearization)
    # Equation 3-120b: q_i,n = q_i,n - β(τ*_si/τ*_sci)^(1/2) * (τ_bs,n/τ_bs,s) * (∂η/∂n)
    mobility_ratio = np.maximum(tau_si_star / tau_sci_star, 1.0)
    shear_ratio = tau_bs_n / (tau_bs_s + 1e-12)  # Avoid division by zero
    
    # Normal transport component
    slope_correction = beta * np.sqrt(mobility_ratio) * shear_ratio * grad_eta_n
    q_i_n = -slope_correction * q_i_s  # Negative because transport opposes slope
    
    return q_i_s, q_i_n


def suspended_sediment_advection_diffusion(c_i: np.ndarray,
                                          u: np.ndarray,
                                          w: np.ndarray,
                                          v_si: Union[float, np.ndarray],
                                          D_d: np.ndarray,
                                          z: np.ndarray,
                                          dx: float = 1.0,
                                          dz: float = 0.1,
                                          dt: float = 1.0) -> np.ndarray:
    """
    Advection-diffusion equation for suspended sediment transport.
    
    Equation 3-128: ∂c_i/∂t + u∂c_i/∂s + (w - v_si)∂c_i/∂z = ∂/∂z(D_d ∂c_i/∂z)
    
    Args:
        c_i: Sediment concentration profile for size i (kg/m³)
        u: Streamwise velocity profile (m/s)
        w: Vertical velocity profile (m/s)
        v_si: Settling velocity for size i (m/s)
        D_d: Eddy diffusivity profile (m²/s)
        z: Vertical coordinate array (m)
        dx: Streamwise grid spacing (m), default 1.0
        dz: Vertical grid spacing (m), default 0.1
        dt: Time step (s), default 1.0
        
    Returns:
        Rate of concentration change ∂c_i/∂t (kg/m³/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-128
    """
    validate_array(c_i)
    validate_array(u)
    validate_array(w)
    validate_array(D_d)
    validate_array(z)
    
    # Spatial derivatives
    dc_dx = np.gradient(c_i, dx)
    dc_dz = np.gradient(c_i, dz)
    
    # Diffusion term: ∂/∂z(D_d ∂c_i/∂z)
    diffusion_flux = D_d * dc_dz
    diffusion_term = np.gradient(diffusion_flux, dz)
    
    # Advection terms
    streamwise_advection = u * dc_dx
    vertical_advection = (w - v_si) * dc_dz
    
    # Full advection-diffusion equation
    dc_dt = -streamwise_advection - vertical_advection + diffusion_term
    
    return dc_dt


def garcia_parker_entrainment(u_star: Union[float, np.ndarray],
                             v_si: Union[float, np.ndarray],
                             D_i: Union[float, np.ndarray],
                             D50: float,
                             R_pi: Union[float, np.ndarray],
                             sigma: float,
                             modified: bool = False,
                             S: Optional[float] = None) -> Union[float, np.ndarray]:
    """
    Garcia-Parker entrainment relation for suspension.
    
    Equations 3-143a-e: E*_si = A*Z_ui^5/(1 + A*Z_ui^5/0.3)
    
    Args:
        u_star: Shear velocity (m/s)
        v_si: Settling velocity for size i (m/s)
        D_i: Grain size for size i (m)
        D50: Median grain size (m)
        R_pi: Particle Reynolds number for size i
        sigma: Geometric standard deviation of bed material
        modified: If True, use modified version with slope effect
        S: Bed slope (for modified version only)
        
    Returns:
        Dimensionless entrainment rate E*_si
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-143a-e
    """
    validate_positive(u_star, "shear_velocity")
    validate_positive(v_si, "settling_velocity")
    validate_positive(D_i, "grain_size")
    validate_positive(D50, "median_grain_size")
    
    # Shear velocity ratio
    u_star_vs_ratio = u_star / v_si
    
    # Size ratio
    size_ratio = D_i / D50
    
    # Correction factor
    lambda_m = 1.0 - 0.298 * sigma
    
    if modified and S is not None:
        # Modified version with slope effect
        A = 7.8e-7
        Z_ui = lambda_m * (u_star_vs_ratio**0.6) * (S**0.08) * (R_pi**0.2) * (size_ratio**0.02)
    else:
        # Original version
        A = 1.3e-7
        Z_ui = lambda_m * (u_star_vs_ratio**0.6) * (R_pi**0.2) * (size_ratio**0.02)
    
    # Garcia-Parker entrainment relation
    E_star_si = A * (Z_ui**5) / (1.0 + A * (Z_ui**5) / 0.3)
    
    return E_star_si


def mclean_near_bed_concentration(tau_bs: Union[float, np.ndarray],
                                 tau_bsc: Union[float, np.ndarray],
                                 f_bi: np.ndarray,
                                 u_star_s: Union[float, np.ndarray],
                                 v_si: Union[float, np.ndarray],
                                 D84: float) -> Tuple[np.ndarray, float]:
    """
    McLean near-bed concentration model for suspended sediment.
    
    Equations 3-144a-n: Complex near-bed concentration relations
    
    Args:
        tau_bs: Boundary shear stress due to skin friction (Pa)
        tau_bsc: Critical boundary shear stress (Pa)
        f_bi: Bed-load size fractions (dimensionless)
        u_star_s: Skin friction shear velocity (m/s)
        v_si: Settling velocities by size class (m/s)
        D84: 84th percentile grain size (m)
        
    Returns:
        Tuple of (near_bed_concentrations, reference_elevation)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-144a-n
    """
    validate_positive(tau_bs, "boundary_shear_stress")
    validate_positive(tau_bsc, "critical_shear_stress")
    validate_array(f_bi)
    validate_positive(u_star_s, "shear_velocity")
    validate_array(v_si)
    validate_positive(D84, "D84_grain_size")
    
    # Total near-bed concentration (Equation 3-144b)
    gamma_o = 0.004
    stress_ratio = tau_bs / tau_bsc
    c_bT = gamma_o * (stress_ratio - 1.0) / (1.0 + gamma_o * (stress_ratio - 1.0))
    
    # Availability factors (Equations 3-144d,e)
    phi_i = np.zeros_like(f_bi)
    u_star_sc = np.sqrt(tau_bsc / WATER_DENSITY)
    
    for i in range(len(f_bi)):
        ratio = u_star_s / v_si[i]
        if ratio > 1.0:
            phi_i[i] = 1.0
        else:
            phi_i[i] = (u_star_s - u_star_sc) / (v_si[i] - u_star_sc)
    
    # Size-specific fractions (Equation 3-144c)
    denominator = np.sum(phi_i * f_bi)
    if denominator > 0:
        fs_bi = (phi_i * f_bi) / denominator
    else:
        fs_bi = f_bi  # Fallback to bed-load fractions
    
    # Near-bed concentrations by size class
    c_bi = fs_bi * c_bT
    
    # Reference elevation (Equations 3-144h,i,j)
    a_D = 0.12
    a_o = 0.056
    A1 = 0.68
    
    # Simplified reference elevation using D84
    ln_D84 = np.log(D84 * 1000)  # Convert to mm for empirical relation
    A2 = 0.0204 * (ln_D84**2) + 0.022 * ln_D84 + 0.0709
    
    delta_B = D84 * A1 * (stress_ratio - 1.0) / (1.0 + A2 * (stress_ratio - 1.0))
    z_b = max(a_D * D84, a_o * delta_B)
    
    return c_bi, z_b


def overbank_deposition_rate(C_ucfi: Union[float, np.ndarray],
                            Q_f: float,
                            B_f: float,
                            v_si: Union[float, np.ndarray],
                            alpha_f: float = 1.0,
                            F_l: float = 1.0) -> Union[float, np.ndarray]:
    """
    Calculate overbank deposition rate on floodplains.
    
    Equations 3-152, 3-153: Volume deposition rate with floodplain scaling
    
    Args:
        C_ucfi: Concentration in channel flow above bank-full (dimensionless)
        Q_f: Floodplain discharge (m³/s)
        B_f: Floodplain width (m)
        v_si: Settling velocities by size class (m/s)
        alpha_f: Attenuation coefficient, default 1.0
        F_l: Floodplain number, default 1.0
        
    Returns:
        Volume deposition rate q_obi (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equations 3-152, 3-153
    """
    validate_positive(Q_f, "floodplain_discharge")
    validate_positive(B_f, "floodplain_width")
    validate_positive(v_si, "settling_velocities")
    
    # Deposition efficiency
    settling_parameter = alpha_f * v_si * (B_f**2) / Q_f
    deposition_efficiency = 1.0 - np.exp(-settling_parameter)
    
    # Volume deposition rate per unit floodplain area
    D_fpi = (F_l * C_ucfi * Q_f / (B_f**2)) * deposition_efficiency
    
    # Total deposition rate per unit channel length
    q_obi = D_fpi * B_f
    
    return q_obi


def tracer_based_bedload_transport(v_bi: Union[float, np.ndarray],
                                  L_a: float,
                                  F_i: np.ndarray,
                                  lambda_p: float = 0.4) -> np.ndarray:
    """
    Tracer-based bed-load transport calculation.
    
    Equation 3-155: q_i = (1-λ_p) v_bi L_a F_i
    
    Args:
        v_bi: Mean virtual velocity of grain size i (m/s)
        L_a: Active layer thickness (m)
        F_i: Grain-size fractions in active layer (dimensionless)
        lambda_p: Porosity, default 0.4
        
    Returns:
        Bed-load transport rates by size class (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-155
    """
    validate_positive(L_a, "active_layer_thickness")
    validate_array(F_i)
    
    # Tracer-based transport relation
    q_i = (1.0 - lambda_p) * v_bi * L_a * F_i
    
    return q_i


def stratification_correction(D_do: Union[float, np.ndarray],
                             c_i: np.ndarray,
                             u_profile: np.ndarray,
                             z: np.ndarray,
                             R_T: float = 1.0) -> Union[float, np.ndarray]:
    """
    Apply density stratification correction to eddy diffusivity.
    
    Equation 3-133: D_d = D_do * F_strat(Ri_g), Ri_g = -R_T g Σ(∂c_i/∂z) / (∂u/∂z)²
    
    Args:
        D_do: Eddy diffusivity without stratification (m²/s)
        c_i: Concentration profiles by size class (kg/m³)
        u_profile: Velocity profile (m/s)
        z: Vertical coordinate (m)
        R_T: Turbulent Schmidt number, default 1.0
        
    Returns:
        Stratification-corrected eddy diffusivity (m²/s)
        
    References:
        ASCE Manual 110, Chapter 3, Equation 3-133
    """
    validate_array(c_i)
    validate_array(u_profile)
    validate_array(z)
    
    # Vertical gradients
    du_dz = np.gradient(u_profile, z)
    
    # Total concentration gradient
    total_c = np.sum(c_i, axis=0) if c_i.ndim > 1 else c_i
    dc_dz = np.gradient(total_c, z)
    
    # Gradient Richardson number
    Ri_g = -R_T * GRAVITY * dc_dz / (du_dz**2 + 1e-12)  # Avoid division by zero
    
    # Stratification function
    F_strat = 1.0 / (1.0 + 4.7 * np.maximum(Ri_g, 0.0))
    
    # Corrected diffusivity
    D_d = D_do * F_strat
    
    return D_d