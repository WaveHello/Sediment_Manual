"""
Channel Migration and Evolution Simulation

This module implements models for simulating long-term channel migration
and meander evolution, including Hickin mapping and process-based models.

Equations implemented:
    - 8-45, 8-46: Transverse bed slope variations in laboratory model
    - 8-47, 8-48: Hickin mapping for channel migration simulation

References:
    ASCE Manual 110, Chapter 8, Sections 8.5-8.6
    Johannesson & Parker (1985), Garcia et al. (1994), Howard (1992)
"""

import numpy as np
from typing import Union, Tuple, Dict, List, Optional, Callable
from scipy.integrate import solve_ivp
from ..utils.validators import validate_positive, validate_range

def laboratory_bed_slope_odg(streamwise_position_s: float,
                           wavelength_L: float,
                           amplitude_factor: float = 0.071,
                           phase_offset: float = 0.8) -> float:
    """
    Calculate transverse bed slope in laboratory model using ODG predictions.
    
    Based on IIHR laboratory bend experiments analyzed by Odgaard.
    
    Parameters:
    -----------
    streamwise_position_s : float
        Streamwise position along centerline (m)
    wavelength_L : float
        Meander wavelength (m)
    amplitude_factor : float, optional
        Amplitude factor (default=0.071)
    phase_offset : float, optional
        Phase offset (default=0.8)
        
    Returns:
    --------
    float
        Transverse bed slope at position s
        
    Equation: 8-45
        Stc(ODG) = 0.071 * sin(2πs/L - 0.8)
        
    Reference: Odgaard & Bergs (1988)
    """
    validate_positive(wavelength_L, "wavelength_L")
    validate_positive(amplitude_factor, "amplitude_factor")
    
    phase_argument = 2 * np.pi * streamwise_position_s / wavelength_L - phase_offset
    return amplitude_factor * np.sin(phase_argument)

def laboratory_bed_slope_ikd(streamwise_position_s: float,
                           wavelength_L: float,
                           amplitude_factor: float = 0.039,
                           phase_offset: float = 0.8) -> float:
    """
    Calculate transverse bed slope in laboratory model using IKD predictions.
    
    Parameters:
    -----------
    streamwise_position_s : float
        Streamwise position along centerline (m)
    wavelength_L : float
        Meander wavelength (m)
    amplitude_factor : float, optional
        Amplitude factor for IKD model (default=0.039)
    phase_offset : float, optional
        Phase offset (default=0.8)
        
    Returns:
    --------
    float
        Transverse bed slope at position s
        
    Equation: 8-46
        Stc(IKD) = 0.039 * sin(2πs/L - 0.8)
        
    Reference: Ikeda et al. (1981)
    """
    validate_positive(wavelength_L, "wavelength_L")
    validate_positive(amplitude_factor, "amplitude_factor")
    
    phase_argument = 2 * np.pi * streamwise_position_s / wavelength_L - phase_offset
    return amplitude_factor * np.sin(phase_argument)

def hickin_mapping_x_component(centerline_point_xp: float,
                              erosion_velocity_ve: float,
                              local_angle_theta: float) -> float:
    """
    Calculate x-component of centerline migration using Hickin mapping.
    
    Parameters:
    -----------
    centerline_point_xp : float
        Current x-coordinate of centerline point (m)
    erosion_velocity_ve : float
        Bank erosion velocity (m/s or m/year)
    local_angle_theta : float
        Local angle of centerline with x-axis (radians)
        
    Returns:
    --------
    float
        Rate of change of x-coordinate (m/s or m/year)
        
    Equation: 8-47
        dxp/dt = ve * sin(θ)
        
    Reference: Hickin (1974), Johannesson & Parker (1985)
    """
    validate_positive(erosion_velocity_ve, "erosion_velocity_ve")
    
    return erosion_velocity_ve * np.sin(local_angle_theta)

def hickin_mapping_y_component(centerline_point_yp: float,
                              erosion_velocity_ve: float,
                              local_angle_theta: float) -> float:
    """
    Calculate y-component of centerline migration using Hickin mapping.
    
    Parameters:
    -----------
    centerline_point_yp : float
        Current y-coordinate of centerline point (m)
    erosion_velocity_ve : float
        Bank erosion velocity (m/s or m/year)
    local_angle_theta : float
        Local angle of centerline with x-axis (radians)
        
    Returns:
    --------
    float
        Rate of change of y-coordinate (m/s or m/year)
        
    Equation: 8-48
        dyp/dt = -ve * cos(θ)
        
    Reference: Hickin (1974), Johannesson & Parker (1985)
    """
    validate_positive(erosion_velocity_ve, "erosion_velocity_ve")
    
    return -erosion_velocity_ve * np.cos(local_angle_theta)

class MeanderEvolutionSimulator:
    """
    Class for simulating long-term meander evolution using various models.
    """
    
    def __init__(self, 
                 erosion_model: str = 'hickin',
                 erosion_parameters: Optional[Dict[str, float]] = None,
                 time_step_years: float = 1.0):
        """
        Initialize meander evolution simulator.
        
        Parameters:
        -----------
        erosion_model : str
            Erosion model type: 'hickin', 'ikeda', 'odgaard'
        erosion_parameters : dict, optional
            Parameters for erosion model
        time_step_years : float
            Time step for simulation (years)
        """
        self.erosion_model = erosion_model
        self.erosion_parameters = erosion_parameters or {}
        self.time_step_years = time_step_years
        
        # Default parameters
        self.default_params = {
            'erodibility_constant': 1e-7,  # m/s
            'width_exponent': 1.0,
            'curvature_exponent': 1.0,
            'discharge_exponent': 0.5
        }
        
        # Update with user parameters
        self.default_params.update(self.erosion_parameters)
    
    def calculate_erosion_rate(self, 
                             channel_width: float,
                             radius_curvature: float,
                             discharge: Optional[float] = None) -> float:
        """
        Calculate bank erosion rate based on selected model.
        
        Parameters:
        -----------
        channel_width : float
            Channel width (m)
        radius_curvature : float
            Local radius of curvature (m)
        discharge : float, optional
            Discharge for discharge-dependent models (m³/s)
            
        Returns:
        --------
        float
            Erosion rate (m/year)
        """
        validate_positive(channel_width, "channel_width")
        validate_positive(radius_curvature, "radius_curvature")
        
        if self.erosion_model == 'hickin':
            # Hickin-Nanson type relationship
            from .bank_erosion import nanson_hickin_migration_rate
            return nanson_hickin_migration_rate(channel_width, radius_curvature)
            
        elif self.erosion_model == 'width_based':
            # Simple width-based erosion
            return self.default_params['erodibility_constant'] * \
                   (channel_width ** self.default_params['width_exponent'])
                   
        elif self.erosion_model == 'curvature_based':
            # Curvature-dependent erosion
            curvature = 1.0 / radius_curvature
            return self.default_params['erodibility_constant'] * \
                   (curvature ** self.default_params['curvature_exponent']) * \
                   (channel_width ** self.default_params['width_exponent'])
                   
        elif self.erosion_model == 'discharge_based' and discharge is not None:
            # Discharge-dependent erosion (Hooke-type)
            from .bank_erosion import drainage_area_migration_rate
            # Convert discharge to approximate drainage area
            drainage_area = discharge * 2.0  # Rough approximation
            return drainage_area_migration_rate(drainage_area)
            
        else:
            raise ValueError(f"Unknown erosion model: {self.erosion_model}")
    
    def simulate_centerline_evolution(self,
                                    initial_centerline: np.ndarray,
                                    channel_width: float,
                                    simulation_time_years: float,
                                    output_interval_years: float = 10.0,
                                    cutoff_threshold: Optional[float] = None) -> Dict[str, np.ndarray]:
        """
        Simulate centerline evolution over time.
        
        Parameters:
        -----------
        initial_centerline : array
            Initial centerline coordinates [[x1,y1], [x2,y2], ...]
        channel_width : float
            Channel width (m)
        simulation_time_years : float
            Total simulation time (years)
        output_interval_years : float
            Time interval for output (years)
        cutoff_threshold : float, optional
            Distance threshold for automatic cutoffs (m)
            
        Returns:
        --------
        dict
            Simulation results with centerlines at different times
        """
        validate_positive(channel_width, "channel_width")
        validate_positive(simulation_time_years, "simulation_time_years")
        validate_positive(output_interval_years, "output_interval_years")
        
        if initial_centerline.shape[1] != 2:
            raise ValueError("Initial centerline must have shape (n_points, 2)")
        
        # Initialize arrays
        n_points = len(initial_centerline)
        current_centerline = initial_centerline.copy()
        
        # Time arrays
        time_years = np.arange(0, simulation_time_years + output_interval_years, output_interval_years)
        n_outputs = len(time_years)
        
        # Storage arrays
        centerlines = np.zeros((n_outputs, n_points, 2))
        centerlines[0] = current_centerline
        
        # Migration tracking
        migration_rates = []
        sinuosities = []
        
        # Main simulation loop
        for i, t in enumerate(time_years[1:], 1):
            # Calculate local properties along centerline
            for j in range(1, n_points - 1):  # Skip endpoints
                # Local angle calculation
                dx = current_centerline[j+1, 0] - current_centerline[j-1, 0]
                dy = current_centerline[j+1, 1] - current_centerline[j-1, 1]
                local_angle = np.arctan2(dy, dx)
                
                # Local radius of curvature estimation
                # Using three-point circle fitting
                x1, y1 = current_centerline[j-1]
                x2, y2 = current_centerline[j]
                x3, y3 = current_centerline[j+1]
                
                # Calculate radius of curvature
                try:
                    radius_curvature = self._calculate_radius_of_curvature(x1, y1, x2, y2, x3, y3)
                    radius_curvature = max(radius_curvature, channel_width)  # Minimum radius
                except:
                    radius_curvature = 10 * channel_width  # Default large radius
                
                # Calculate erosion rate
                erosion_rate = self.calculate_erosion_rate(channel_width, radius_curvature)
                
                # Convert to velocity (m/year to m/time_step)
                erosion_velocity = erosion_rate * self.time_step_years
                
                # Apply Hickin mapping
                dx_dt = hickin_mapping_x_component(0, erosion_velocity, local_angle)
                dy_dt = hickin_mapping_y_component(0, erosion_velocity, local_angle)
                
                # Update position
                current_centerline[j, 0] += dx_dt
                current_centerline[j, 1] += dy_dt
            
            # Check for cutoffs if threshold specified
            if cutoff_threshold is not None:
                current_centerline = self._check_and_apply_cutoffs(current_centerline, cutoff_threshold)
            
            # Store results
            centerlines[i] = current_centerline
            
            # Calculate metrics
            sinuosity = self._calculate_sinuosity(current_centerline)
            sinuosities.append(sinuosity)
            
            # Average migration rate
            if i > 1:
                migration_distance = np.mean(np.sqrt(np.sum((centerlines[i] - centerlines[i-1])**2, axis=1)))
                migration_rate = migration_distance / output_interval_years
                migration_rates.append(migration_rate)
        
        return {
            'time_years': time_years,
            'centerlines': centerlines,
            'sinuosities': np.array(sinuosities),
            'migration_rates': np.array(migration_rates),
            'channel_width': channel_width,
            'simulation_parameters': self.default_params
        }
    
    def _calculate_radius_of_curvature(self, x1: float, y1: float, 
                                     x2: float, y2: float, 
                                     x3: float, y3: float) -> float:
        """
        Calculate radius of curvature from three points using circle fitting.
        """
        # Calculate circle through three points
        ax = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)
        if abs(ax) < 1e-10:
            return 1e6  # Nearly straight - large radius
        
        bx = -(x1**2 + y1**2) * (y2 - y3) - (x2**2 + y2**2) * (y3 - y1) - (x3**2 + y3**2) * (y1 - y2)
        by = -(x1**2 + y1**2) * (x3 - x2) - (x2**2 + y2**2) * (x1 - x3) - (x3**2 + y3**2) * (x2 - x1)
        
        # Center of circle
        cx = bx / (2 * ax)
        cy = by / (2 * ax)
        
        # Radius
        radius = np.sqrt((x1 - cx)**2 + (y1 - cy)**2)
        
        return radius
    
    def _calculate_sinuosity(self, centerline: np.ndarray) -> float:
        """
        Calculate sinuosity of centerline.
        """
        # Channel length
        distances = np.sqrt(np.sum(np.diff(centerline, axis=0)**2, axis=1))
        channel_length = np.sum(distances)
        
        # Valley length (straight-line distance)
        valley_length = np.sqrt((centerline[-1, 0] - centerline[0, 0])**2 + 
                               (centerline[-1, 1] - centerline[0, 1])**2)
        
        if valley_length < 1e-10:
            return 1.0
        
        return channel_length / valley_length
    
    def _check_and_apply_cutoffs(self, centerline: np.ndarray, threshold: float) -> np.ndarray:
        """
        Check for and apply cutoffs when channel loops come too close.
        """
        n_points = len(centerline)
        cutoff_applied = False
        
        for i in range(n_points - 10):  # Avoid checking very close points
            for j in range(i + 10, n_points):
                distance = np.sqrt((centerline[i, 0] - centerline[j, 0])**2 + 
                                 (centerline[i, 1] - centerline[j, 1])**2)
                
                if distance < threshold:
                    # Apply cutoff - connect points i and j directly
                    new_centerline = np.vstack([
                        centerline[:i+1],
                        centerline[j:]
                    ])
                    cutoff_applied = True
                    break
            
            if cutoff_applied:
                break
        
        return new_centerline if cutoff_applied else centerline

def stochastic_migration_model(initial_centerline: np.ndarray,
                             channel_width: float,
                             simulation_time_years: float,
                             random_factor: float = 0.1,
                             correlation_length: float = 100.0,
                             random_seed: Optional[int] = None) -> Dict[str, np.ndarray]:
    """
    Simulate meander evolution with stochastic components.
    
    Combines deterministic migration with random perturbations
    representing natural variability.
    
    Parameters:
    -----------
    initial_centerline : array
        Initial centerline coordinates
    channel_width : float
        Channel width (m)
    simulation_time_years : float
        Simulation time (years)
    random_factor : float
        Magnitude of random perturbations
    correlation_length : float
        Spatial correlation length for random perturbations (m)
    random_seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    dict
        Stochastic simulation results
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    validate_positive(channel_width, "channel_width")
    validate_positive(simulation_time_years, "simulation_time_years")
    validate_range(random_factor, 0.0, 1.0, "random_factor")
    validate_positive(correlation_length, "correlation_length")
    
    # Initialize deterministic simulator
    simulator = MeanderEvolutionSimulator(erosion_model='hickin')
    
    # Run deterministic simulation
    deterministic_result = simulator.simulate_centerline_evolution(
        initial_centerline, channel_width, simulation_time_years
    )
    
    # Add stochastic perturbations
    stochastic_centerlines = deterministic_result['centerlines'].copy()
    n_times, n_points, _ = stochastic_centerlines.shape
    
    # Generate spatially correlated random perturbations
    for t in range(1, n_times):
        # Generate random perturbations
        perturbations_x = np.random.normal(0, random_factor * channel_width, n_points)
        perturbations_y = np.random.normal(0, random_factor * channel_width, n_points)
        
        # Apply spatial correlation (simple moving average)
        window_size = max(3, int(correlation_length / (channel_width * 5)))
        
        # Smooth perturbations to create spatial correlation
        perturbations_x = np.convolve(perturbations_x, 
                                    np.ones(window_size)/window_size, 
                                    mode='same')
        perturbations_y = np.convolve(perturbations_y, 
                                    np.ones(window_size)/window_size, 
                                    mode='same')
        
        # Add perturbations to centerline
        stochastic_centerlines[t, :, 0] += perturbations_x
        stochastic_centerlines[t, :, 1] += perturbations_y
    
    # Update results
    result = deterministic_result.copy()
    result['centerlines'] = stochastic_centerlines
    result['random_factor'] = random_factor
    result['correlation_length'] = correlation_length
    result['model_type'] = 'stochastic'
    
    return result