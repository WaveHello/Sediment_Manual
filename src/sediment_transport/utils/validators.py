"""
Input validation functions for sediment transport calculations.
"""

import numpy as np
from typing import Union, Optional, Any

def validate_positive(
    value: Union[float, np.ndarray], 
    name: str
) -> None:
    """
    Validate that a value or array is positive.
    
    Parameters
    ----------
    value : float or array-like
        Value(s) to validate
    name : str
        Parameter name for error message
        
    Raises
    ------
    ValueError
        If any value is not positive
    """
    value_array = np.atleast_1d(value)
    
    if np.any(value_array <= 0):
        raise ValueError(f"{name} must be positive, got {value}")


def validate_range(
    value: Union[float, np.ndarray],
    name: str,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    min_inclusive: bool = True,
    max_inclusive: bool = True
) -> None:
    """
    Validate that a value or array is within a specified range.
    
    Parameters
    ----------
    value : float or array-like
        Value(s) to validate
    name : str
        Parameter name for error message
    min_val : float, optional
        Minimum allowed value
    max_val : float, optional
        Maximum allowed value
    min_inclusive : bool, optional
        Whether minimum is inclusive, default True
    max_inclusive : bool, optional
        Whether maximum is inclusive, default True
        
    Raises
    ------
    ValueError
        If any value is outside the specified range
    """
    value_array = np.atleast_1d(value)
    
    if min_val is not None:
        if min_inclusive:
            if np.any(value_array < min_val):
                raise ValueError(f"{name} must be >= {min_val}, got {value}")
        else:
            if np.any(value_array <= min_val):
                raise ValueError(f"{name} must be > {min_val}, got {value}")
    
    if max_val is not None:
        if max_inclusive:
            if np.any(value_array > max_val):
                raise ValueError(f"{name} must be <= {max_val}, got {value}")
        else:
            if np.any(value_array >= max_val):
                raise ValueError(f"{name} must be < {max_val}, got {value}")


def validate_array_like(
    value: Any,
    name: str,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None
) -> np.ndarray:
    """
    Validate and convert input to numpy array.
    
    Parameters
    ----------
    value : array-like
        Input to convert
    name : str
        Parameter name for error message
    min_length : int, optional
        Minimum array length
    max_length : int, optional
        Maximum array length
        
    Returns
    -------
    ndarray
        Validated numpy array
        
    Raises
    ------
    ValueError
        If input cannot be converted to array or length is invalid
    """
    try:
        arr = np.atleast_1d(value)
    except (ValueError, TypeError):
        raise ValueError(f"{name} must be array-like, got {type(value)}")
    
    if min_length is not None and len(arr) < min_length:
        raise ValueError(f"{name} must have at least {min_length} elements, got {len(arr)}")
    
    if max_length is not None and len(arr) > max_length:
        raise ValueError(f"{name} must have at most {max_length} elements, got {len(arr)}")
        
    return arr


def validate_units(
    value: Union[float, np.ndarray],
    name: str,
    expected_units: str,
    typical_range: Optional[tuple] = None
) -> None:
    """
    Validate that a value is in reasonable range for expected units.
    
    Parameters
    ----------
    value : float or array-like
        Value(s) to validate
    name : str
        Parameter name
    expected_units : str
        Expected units (for documentation)
    typical_range : tuple, optional
        (min, max) typical range for the parameter
        
    Raises
    ------
    ValueError
        If values are outside typical range
    """
    if typical_range is None:
        return
        
    value_array = np.atleast_1d(value)
    min_val, max_val = typical_range
    
    if np.any(value_array < min_val) or np.any(value_array > max_val):
        import warnings
        warnings.warn(
            f"{name} values {value} are outside typical range "
            f"{typical_range} for units {expected_units}. "
            f"Check units and values.",
            UserWarning
        )