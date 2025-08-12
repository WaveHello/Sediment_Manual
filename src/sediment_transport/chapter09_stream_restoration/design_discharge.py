"""
Chapter 9: Stream Restoration - Design Discharge Module

This module implements channel-forming discharge calculation methods including
effective discharge (Qeff), bank-full discharge (Qbf), and return interval 
discharge (Qri) for stream restoration design based on ASCE Manual 110 Chapter 9.

References:
- ASCE Manual 110, Chapter 9: Stream Restoration, Section 9.3
- Andrews (1980), Wolman and Miller (1960)
- Soar and Thorne (2001), Biedenharn et al. (2000)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Callable
from enum import Enum
from dataclasses import dataclass
from scipy import stats, optimize
import warnings

from ..utils.validators import validate_positive, validate_range
from ..utils.constants import *


class DischargeMethod(Enum):
    """Methods for determining channel-forming discharge"""
    EFFECTIVE_DISCHARGE = "effective_discharge"    # Qeff - most sediment transport
    BANKFULL_DISCHARGE = "bankfull_discharge"     # Qbf - field indicators
    RETURN_INTERVAL = "return_interval"           # Qri - statistical frequency


class StreamType(Enum):
    """Stream type classification affecting discharge relationships"""
    PERENNIAL_SNOWMELT = "perennial_snowmelt"    # Montane west streams
    PERENNIAL_RAINFALL = "perennial_rainfall"    # Eastern humid streams  
    EPHEMERAL_FLASHY = "ephemeral_flashy"       # Arid/semiarid streams
    REGULATED = "regulated"                      # Dam-affected streams


@dataclass
class DischargeFrequency:
    """Discharge frequency distribution data"""
    discharges: np.ndarray          # Discharge values (m³/s)
    frequencies: np.ndarray         # Frequency of occurrence (fraction)
    flow_duration_curve: Optional[np.ndarray] = None
    period_of_record: int = 0


@dataclass
class SedimentRating:
    """Sediment transport rating curve"""
    discharges: np.ndarray          # Discharge values (m³/s)  
    sediment_loads: np.ndarray      # Sediment transport rates (kg/s)
    rating_equation: Optional[Callable] = None
    r_squared: float = 0.0


@dataclass
class DesignDischarge:
    """Design discharge calculation results"""
    method: DischargeMethod
    discharge_value: float          # m³/s
    confidence_interval: Tuple[float, float]
    data_quality: str              # high, moderate, low
    limitations: List[str]
    recommended_use: str


def create_discharge_histogram(daily_flows: np.ndarray, 
                             num_bins: int = 50,
                             min_discharge: Optional[float] = None) -> DischargeFrequency:
    """
    Create discharge frequency histogram from daily flow data.
    
    Args:
        daily_flows: Array of daily discharge values (m³/s)
        num_bins: Number of histogram bins
        min_discharge: Minimum discharge threshold
        
    Returns:
        DischargeFrequency object with histogram data
        
    Reference: ASCE Manual 110, Section 9.3.1.1.2, Step 1
    """
    if len(daily_flows) == 0:
        raise ValueError("No flow data provided")
    
    # Remove zero or negative flows if min_discharge specified
    if min_discharge is not None:
        daily_flows = daily_flows[daily_flows >= min_discharge]
    
    # Create histogram
    hist_counts, bin_edges = np.histogram(daily_flows, bins=num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Convert counts to frequencies
    total_days = len(daily_flows)
    frequencies = hist_counts / total_days
    
    # Create flow duration curve
    sorted_flows = np.sort(daily_flows)[::-1]  # Descending order
    exceedance_probs = np.arange(1, len(sorted_flows) + 1) / len(sorted_flows)
    flow_duration_curve = np.column_stack([exceedance_probs, sorted_flows])
    
    return DischargeFrequency(
        discharges=bin_centers,
        frequencies=frequencies,
        flow_duration_curve=flow_duration_curve,
        period_of_record=len(daily_flows) // 365
    )


def develop_sediment_rating_curve(discharge_data: np.ndarray,
                                sediment_data: np.ndarray,
                                rating_form: str = "power_law") -> SedimentRating:
    """
    Develop sediment transport rating curve from measured data.
    
    Args:
        discharge_data: Measured discharge values (m³/s)
        sediment_data: Measured sediment loads (kg/s)
        rating_form: Form of rating equation ('power_law', 'exponential')
        
    Returns:
        SedimentRating object with fitted curve
        
    Reference: ASCE Manual 110, Section 9.3.1.1.2, Step 2
    """
    validate_positive(discharge_data.min(), "discharge_data")
    validate_positive(sediment_data.min(), "sediment_data")
    
    if len(discharge_data) != len(sediment_data):
        raise ValueError("Discharge and sediment data must have same length")
    
    # Remove zero values for log transformation
    valid_mask = (discharge_data > 0) & (sediment_data > 0)
    Q_valid = discharge_data[valid_mask]
    Qs_valid = sediment_data[valid_mask]
    
    if len(Q_valid) < 3:
        raise ValueError("Insufficient valid data points for rating curve")
    
    if rating_form == "power_law":
        # Fit: Qs = a * Q^b
        log_Q = np.log10(Q_valid)
        log_Qs = np.log10(Qs_valid)
        
        # Linear regression in log space
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_Q, log_Qs)
        
        # Convert back to power law parameters
        a = 10**intercept
        b = slope
        r_squared = r_value**2
        
        # Create rating equation function
        def rating_equation(Q):
            return a * (Q**b)
            
    elif rating_form == "exponential":
        # Fit: Qs = a * exp(b * Q)
        def exp_func(Q, a, b):
            return a * np.exp(b * Q)
        
        try:
            from scipy.optimize import curve_fit
            popt, pcov = curve_fit(exp_func, Q_valid, Qs_valid, 
                                 p0=[1.0, 0.01], maxfev=1000)
            a, b = popt
            
            # Calculate R²
            y_pred = exp_func(Q_valid, a, b)
            ss_res = np.sum((Qs_valid - y_pred)**2)
            ss_tot = np.sum((Qs_valid - np.mean(Qs_valid))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            def rating_equation(Q):
                return exp_func(Q, a, b)
                
        except:
            # Fall back to power law if exponential fit fails
            return develop_sediment_rating_curve(discharge_data, sediment_data, "power_law")
    else:
        raise ValueError(f"Unknown rating form: {rating_form}")
    
    # Generate rating curve points for full discharge range
    Q_range = np.logspace(np.log10(Q_valid.min()), np.log10(Q_valid.max()), 100)
    Qs_range = rating_equation(Q_range)
    
    return SedimentRating(
        discharges=Q_range,
        sediment_loads=Qs_range,
        rating_equation=rating_equation,
        r_squared=r_squared
    )


def calculate_effective_discharge(discharge_freq: DischargeFrequency,
                                sediment_rating: SedimentRating,
                                method: str = "histogram") -> DesignDischarge:
    """
    Calculate effective discharge that transports the most sediment.
    
    Args:
        discharge_freq: Discharge frequency distribution
        sediment_rating: Sediment transport rating curve
        method: Calculation method ('histogram' or 'analytical')
        
    Returns:
        DesignDischarge object with effective discharge
        
    Reference: ASCE Manual 110, Section 9.3.1.1, Figure 9-2
    """
    if method == "histogram":
        # Method using discharge histogram (Figure 9-2 approach)
        Q_bins = discharge_freq.discharges
        freq = discharge_freq.frequencies
        
        # Calculate sediment load for each discharge bin
        if sediment_rating.rating_equation is not None:
            Qs_bins = sediment_rating.rating_equation(Q_bins)
        else:
            # Interpolate from rating curve data
            Qs_bins = np.interp(Q_bins, sediment_rating.discharges, 
                               sediment_rating.sediment_loads)
        
        # Calculate collective sediment discharge (frequency × load)
        collective_sediment = freq * Qs_bins
        
        # Find discharge bin with maximum collective sediment
        max_idx = np.argmax(collective_sediment)
        Qeff = Q_bins[max_idx]
        
        # Estimate confidence interval based on neighboring bins
        if max_idx > 0 and max_idx < len(Q_bins) - 1:
            ci_lower = Q_bins[max_idx - 1]
            ci_upper = Q_bins[max_idx + 1]
        else:
            ci_lower = Qeff * 0.8
            ci_upper = Qeff * 1.2
            
    elif method == "analytical":
        # Analytical method using flow duration curve
        if discharge_freq.flow_duration_curve is None:
            raise ValueError("Flow duration curve required for analytical method")
        
        exceedance = discharge_freq.flow_duration_curve[:, 0]
        flows = discharge_freq.flow_duration_curve[:, 1]
        
        # Calculate sediment loads
        if sediment_rating.rating_equation is not None:
            sediment_loads = sediment_rating.rating_equation(flows)
        else:
            sediment_loads = np.interp(flows, sediment_rating.discharges,
                                     sediment_rating.sediment_loads)
        
        # Calculate incremental sediment transport
        # dQ/dt for each flow increment
        flow_increments = np.diff(flows)
        prob_increments = np.diff(exceedance)
        
        # Avoid division by zero
        prob_increments = np.where(prob_increments == 0, 1e-10, prob_increments)
        
        # Sediment transport per unit time for each increment  
        incremental_transport = sediment_loads[:-1] * np.abs(prob_increments)
        
        # Find flow with maximum incremental transport
        max_idx = np.argmax(incremental_transport)
        Qeff = flows[max_idx]
        
        # Confidence interval
        ci_lower = Qeff * 0.85
        ci_upper = Qeff * 1.15
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Assess data quality
    if discharge_freq.period_of_record >= 30 and sediment_rating.r_squared >= 0.8:
        data_quality = "high"
    elif discharge_freq.period_of_record >= 10 and sediment_rating.r_squared >= 0.6:
        data_quality = "moderate"
    else:
        data_quality = "low"
    
    # Limitations
    limitations = []
    if discharge_freq.period_of_record < 10:
        limitations.append("Short period of record (<10 years)")
    if sediment_rating.r_squared < 0.7:
        limitations.append("Poor sediment rating curve fit")
    limitations.append("Sensitive to number of discharge bins")
    limitations.append("Assumes steady sediment transport relations")
    
    return DesignDischarge(
        method=DischargeMethod.EFFECTIVE_DISCHARGE,
        discharge_value=Qeff,
        confidence_interval=(ci_lower, ci_upper),
        data_quality=data_quality,
        limitations=limitations,
        recommended_use="Channel design and restoration planning"
    )


def estimate_bankfull_discharge_field(channel_width: float,
                                    channel_depth: float,
                                    channel_slope: float,
                                    manning_n: float,
                                    field_indicators: Dict[str, bool]) -> DesignDischarge:
    """
    Estimate bank-full discharge using field indicators and hydraulics.
    
    Args:
        channel_width: Bank-full channel width (m)
        channel_depth: Bank-full channel depth (m) 
        channel_slope: Channel slope (dimensionless)
        manning_n: Manning's roughness coefficient
        field_indicators: Dict of field indicator observations
        
    Returns:
        DesignDischarge object with bank-full discharge estimate
        
    Reference: ASCE Manual 110, Table 9-5
    """
    validate_positive(channel_width, "channel_width")
    validate_positive(channel_depth, "channel_depth")
    validate_positive(channel_slope, "channel_slope")
    validate_positive(manning_n, "manning_n")
    
    # Calculate hydraulic properties
    area = channel_width * channel_depth
    wetted_perimeter = channel_width + 2 * channel_depth
    hydraulic_radius = area / wetted_perimeter
    
    # Calculate bank-full discharge using Manning's equation
    Qbf = (1 / manning_n) * area * (hydraulic_radius**(2/3)) * (channel_slope**0.5)
    
    # Adjust based on field indicators
    confidence_multiplier = 1.0
    limitations = []
    
    if field_indicators.get("active_floodplain", False):
        confidence_multiplier *= 1.1
    else:
        limitations.append("No clear active floodplain identified")
        confidence_multiplier *= 0.9
    
    if field_indicators.get("channel_scour_line", False):
        confidence_multiplier *= 1.05
    else:
        limitations.append("Channel scour line not clearly defined")
    
    if field_indicators.get("vegetation_break", False):
        confidence_multiplier *= 1.05
    else:
        limitations.append("Vegetation break not distinct")
        
    if field_indicators.get("terrace_edge", False):
        confidence_multiplier *= 1.1
    else:
        limitations.append("Terrace edge not well-defined")
    
    # Assess channel stability for data quality
    if len([ind for ind in field_indicators.values() if ind]) >= 3:
        data_quality = "high"
        ci_factor = 0.15  # ±15%
    elif len([ind for ind in field_indicators.values() if ind]) >= 2:
        data_quality = "moderate" 
        ci_factor = 0.25  # ±25%
    else:
        data_quality = "low"
        ci_factor = 0.40  # ±40%
        
    # Apply confidence adjustment
    Qbf_adjusted = Qbf * confidence_multiplier
    
    # Calculate confidence interval
    ci_lower = Qbf_adjusted * (1 - ci_factor)
    ci_upper = Qbf_adjusted * (1 + ci_factor)
    
    # Additional limitations
    limitations.extend([
        "Field indicators can be misleading in unstable channels",
        "Requires stable, alluvial channel conditions",
        "Manning's n estimation introduces uncertainty"
    ])
    
    return DesignDischarge(
        method=DischargeMethod.BANKFULL_DISCHARGE,
        discharge_value=Qbf_adjusted,
        confidence_interval=(ci_lower, ci_upper),
        data_quality=data_quality,
        limitations=limitations,
        recommended_use="Stability assessment and Qeff estimation"
    )


def calculate_return_interval_discharge(annual_peak_flows: np.ndarray,
                                      return_periods: List[int] = [1.5, 2, 5, 10, 25],
                                      distribution: str = "log_pearson_3") -> Dict[int, DesignDischarge]:
    """
    Calculate return interval discharges using flood frequency analysis.
    
    Args:
        annual_peak_flows: Array of annual maximum flows (m³/s)
        return_periods: List of return periods in years
        distribution: Statistical distribution ('log_pearson_3', 'gumbel', 'lognormal')
        
    Returns:
        Dictionary mapping return periods to DesignDischarge objects
        
    Reference: ASCE Manual 110, Table 9-5
    """
    if len(annual_peak_flows) < 10:
        warnings.warn("Less than 10 years of peak flow data - results unreliable")
    
    results = {}
    
    # Remove zero or negative values
    valid_flows = annual_peak_flows[annual_peak_flows > 0]
    
    for return_period in return_periods:
        if distribution == "log_pearson_3":
            # Log-Pearson Type III distribution (standard in US)
            log_flows = np.log10(valid_flows)
            mean_log = np.mean(log_flows)
            std_log = np.std(log_flows, ddof=1)
            
            # Calculate skewness
            skew_log = stats.skew(log_flows, bias=False)
            
            # Frequency factor for Log-Pearson III
            exceedance_prob = 1.0 / return_period
            
            if abs(skew_log) < 0.01:  # Essentially normal
                k_factor = stats.norm.ppf(1 - exceedance_prob)
            else:
                # Use approximate formula for skewed distribution
                k_factor = stats.norm.ppf(1 - exceedance_prob)
                k_factor += (skew_log / 6.0) * (k_factor**2 - 1)
                k_factor += (skew_log**2 / 36.0) * (k_factor**3 - 6*k_factor)
            
            log_Q = mean_log + k_factor * std_log
            Q_return = 10**log_Q
            
        elif distribution == "gumbel":
            # Gumbel (Extreme Value Type I) distribution
            params = stats.gumbel_r.fit(valid_flows)
            Q_return = stats.gumbel_r.ppf(1 - 1/return_period, *params)
            
        elif distribution == "lognormal":
            # 2-parameter lognormal distribution
            log_flows = np.log(valid_flows)
            mean_log = np.mean(log_flows)
            std_log = np.std(log_flows, ddof=1)
            
            # Standard normal quantile
            z = stats.norm.ppf(1 - 1/return_period)
            log_Q = mean_log + z * std_log
            Q_return = np.exp(log_Q)
            
        else:
            raise ValueError(f"Unknown distribution: {distribution}")
        
        # Estimate confidence interval (approximate)
        # Based on standard error of estimate
        n_years = len(valid_flows)
        if n_years >= 30:
            ci_factor = 0.10  # ±10% for long records
            data_quality = "high"
        elif n_years >= 15:
            ci_factor = 0.20  # ±20% for moderate records
            data_quality = "moderate"
        else:
            ci_factor = 0.35  # ±35% for short records  
            data_quality = "low"
        
        ci_lower = Q_return * (1 - ci_factor)
        ci_upper = Q_return * (1 + ci_factor)
        
        # Limitations
        limitations = [
            "No direct physical basis for channel formation",
            f"Based on {n_years} years of peak flow data",
            "Assumes stationary climate conditions",
            "Relations to Qeff and Qbf inconsistent in literature"
        ]
        
        if n_years < 20:
            limitations.append("Short record length reduces reliability")
        
        results[return_period] = DesignDischarge(
            method=DischargeMethod.RETURN_INTERVAL,
            discharge_value=Q_return,
            confidence_interval=(ci_lower, ci_upper),
            data_quality=data_quality,
            limitations=limitations,
            recommended_use=f"First approximation for {return_period}-year design events"
        )
    
    return results


def validate_design_discharge(discharge_value: float,
                            drainage_area: float,
                            mean_annual_flow: float,
                            regional_relationships: Optional[Dict[str, float]] = None) -> Dict[str, Union[bool, str]]:
    """
    Validate computed design discharge against regional relationships and checks.
    
    Args:
        discharge_value: Computed design discharge (m³/s)
        drainage_area: Contributing drainage area (km²)
        mean_annual_flow: Mean annual flow (m³/s)
        regional_relationships: Dict of regional regression parameters
        
    Returns:
        Dictionary of validation checks and results
        
    Reference: ASCE Manual 110, Section 9.3.1.4
    """
    validate_positive(discharge_value, "discharge_value")
    validate_positive(drainage_area, "drainage_area")
    validate_positive(mean_annual_flow, "mean_annual_flow")
    
    validation_results = {}
    
    # Check against mean annual flow
    q_to_mean_ratio = discharge_value / mean_annual_flow
    if 0.5 <= q_to_mean_ratio <= 3.0:
        validation_results["mean_annual_check"] = True
        validation_results["mean_annual_comment"] = f"Ratio Q/Qmean = {q_to_mean_ratio:.2f} is reasonable"
    else:
        validation_results["mean_annual_check"] = False
        validation_results["mean_annual_comment"] = f"Ratio Q/Qmean = {q_to_mean_ratio:.2f} is unusual"
    
    # Check specific discharge (discharge per unit area)
    specific_discharge = discharge_value / drainage_area  # m³/s/km²
    
    if 0.001 <= specific_discharge <= 0.1:  # Typical range for most climates
        validation_results["specific_discharge_check"] = True
        validation_results["specific_discharge_comment"] = f"Specific discharge {specific_discharge:.4f} m³/s/km² is reasonable"
    else:
        validation_results["specific_discharge_check"] = False
        validation_results["specific_discharge_comment"] = f"Specific discharge {specific_discharge:.4f} m³/s/km² is unusual"
    
    # Regional regression check (if provided)
    if regional_relationships is not None:
        # Typical form: Q = a * A^b * P^c (where A=area, P=precipitation)
        a = regional_relationships.get("coefficient", 1.0)
        b = regional_relationships.get("area_exponent", 0.7)
        
        predicted_Q = a * (drainage_area ** b)
        ratio = discharge_value / predicted_Q
        
        if 0.5 <= ratio <= 2.0:
            validation_results["regional_check"] = True
            validation_results["regional_comment"] = f"Ratio to regional estimate = {ratio:.2f}"
        else:
            validation_results["regional_check"] = False
            validation_results["regional_comment"] = f"Ratio to regional estimate = {ratio:.2f} is outside expected range"
    
    # Overall validation
    checks_passed = sum([
        validation_results.get("mean_annual_check", False),
        validation_results.get("specific_discharge_check", False),
        validation_results.get("regional_check", True)  # Default to True if not checked
    ])
    
    total_checks = len([k for k in validation_results.keys() if k.endswith("_check")])
    
    if checks_passed == total_checks:
        validation_results["overall_assessment"] = "PASS - Design discharge appears reasonable"
    elif checks_passed >= total_checks * 0.7:
        validation_results["overall_assessment"] = "CAUTION - Some validation checks failed"
    else:
        validation_results["overall_assessment"] = "FAIL - Multiple validation checks failed"
    
    return validation_results


def select_design_discharge(discharge_estimates: Dict[str, DesignDischarge],
                          project_objectives: List[str],
                          data_availability: Dict[str, str],
                          risk_tolerance: float) -> DesignDischarge:
    """
    Select appropriate design discharge based on project requirements.
    
    Args:
        discharge_estimates: Dict of discharge estimates by method
        project_objectives: List of project objectives
        data_availability: Dict of data quality by type  
        risk_tolerance: Risk tolerance (0-1, higher = more tolerant)
        
    Returns:
        Selected DesignDischarge with justification
        
    Reference: ASCE Manual 110, Table 9-5
    """
    validate_range(risk_tolerance, 0.0, 1.0, "risk_tolerance")
    
    # Scoring system for method selection
    method_scores = {}
    
    for method_name, discharge_est in discharge_estimates.items():
        score = 0
        
        # Data quality scoring
        quality_scores = {"high": 3, "moderate": 2, "low": 1}
        score += quality_scores.get(discharge_est.data_quality, 0)
        
        # Method-specific considerations
        if discharge_est.method == DischargeMethod.EFFECTIVE_DISCHARGE:
            if "channel_design" in [obj.lower() for obj in project_objectives]:
                score += 3  # Preferred for design
            if data_availability.get("sediment", "none") != "none":
                score += 2  # Requires sediment data
            else:
                score -= 3  # Penalize if no sediment data
                
        elif discharge_est.method == DischargeMethod.BANKFULL_DISCHARGE:
            if "stability_assessment" in [obj.lower() for obj in project_objectives]:
                score += 2  # Good for stability assessment
            if data_availability.get("field_survey", "none") != "none":
                score += 1
                
        elif discharge_est.method == DischargeMethod.RETURN_INTERVAL:
            if "preliminary_design" in [obj.lower() for obj in project_objectives]:
                score += 1  # OK for preliminary estimates
            score -= 1  # Generally less preferred
        
        # Risk tolerance adjustment
        num_limitations = len(discharge_est.limitations)
        if risk_tolerance > 0.7:  # High risk tolerance
            score -= num_limitations * 0.5
        else:  # Low risk tolerance
            score -= num_limitations * 1.0
            
        method_scores[method_name] = score
    
    # Select method with highest score
    best_method = max(method_scores.items(), key=lambda x: x[1])
    selected_discharge = discharge_estimates[best_method[0]]
    
    # Add selection justification
    justification = f"Selected based on: data quality ({selected_discharge.data_quality}), "
    justification += f"project objectives, risk tolerance ({risk_tolerance:.1f}). "
    justification += f"Score: {best_method[1]:.1f}"
    
    # Create modified discharge result with selection justification
    return DesignDischarge(
        method=selected_discharge.method,
        discharge_value=selected_discharge.discharge_value,
        confidence_interval=selected_discharge.confidence_interval,
        data_quality=selected_discharge.data_quality,
        limitations=selected_discharge.limitations + [f"Selection: {justification}"],
        recommended_use=selected_discharge.recommended_use
    )