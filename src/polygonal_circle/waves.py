"""Pluggable periodic waveform generators for recursive radial perturbation."""

from __future__ import annotations

from enum import Enum
import math
from typing import Callable, Union
import numpy as np


class WaveformType(str, Enum):
    """Available microscopic waveform generator types."""
    TRIANGLE = "triangle"
    SINE = "sine"
    PARABOLIC = "parabolic"
    MULTI_HARMONIC = "multi_harmonic"
    SMOOTHED_SQUARE = "smoothed_square"


def triangle_wave(u: np.ndarray | float, phase: float = 0.0) -> np.ndarray:
    """Evaluate periodic symmetric triangle wave with period 1 and range [0, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    u_shifted = np.mod(u_arr + phase, 1.0)
    return 1.0 - 2.0 * np.abs(u_shifted - 0.5)


def zero_mean_triangle_wave(u: np.ndarray | float) -> np.ndarray:
    """Zero-mean periodic triangle wave taking values in [-1, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    return 1.0 - 4.0 * np.abs(np.mod(u_arr, 1.0) - 0.5)


def zero_mean_sine_wave(u: np.ndarray | float) -> np.ndarray:
    """Zero-mean pure sinusoidal wave taking values in [-1, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    return np.cos(2.0 * np.pi * u_arr)


def zero_mean_parabolic_wave(u: np.ndarray | float) -> np.ndarray:
    """C^1 smooth piecewise quadratic wave with zero mean and range [-1, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    x = np.mod(u_arr, 1.0) - 0.5  # in [-0.5, 0.5]
    abs_x = np.abs(x)

    # Piecewise quadratic construction:
    # Inner region |x| <= 0.25: 1 - 16 * x^2
    # Outer region |x| > 0.25: -1 + 16 * (|x| - 0.5)^2
    inner = abs_x <= 0.25
    result = np.empty_like(x)
    result[inner] = 1.0 - 16.0 * (x[inner] ** 2)
    result[~inner] = -1.0 + 16.0 * ((abs_x[~inner] - 0.5) ** 2)
    return result


def zero_mean_multi_harmonic_wave(u: np.ndarray | float) -> np.ndarray:
    """Asymmetric multi-frequency harmonic wave with zero mean and range [-1, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    raw = (
        np.cos(2.0 * np.pi * u_arr)
        + 0.40 * np.cos(4.0 * np.pi * u_arr + 0.5)
        - 0.25 * np.sin(6.0 * np.pi * u_arr)
    )
    # Peak-to-peak normalization factor computed over fine grid
    norm_factor = 1.55012
    return raw / norm_factor


def zero_mean_smoothed_square_wave(u: np.ndarray | float, steepness: float = 4.0) -> np.ndarray:
    """Smooth tanh-based pseudo-square wave with zero mean and range [-1, 1]."""
    u_arr = np.asarray(u, dtype=np.float64)
    raw = np.tanh(steepness * np.sin(2.0 * np.pi * u_arr))
    return raw / np.tanh(steepness)


_ZERO_MEAN_REGISTRY: dict[WaveformType, Callable[[np.ndarray | float], np.ndarray]] = {
    WaveformType.TRIANGLE: zero_mean_triangle_wave,
    WaveformType.SINE: zero_mean_sine_wave,
    WaveformType.PARABOLIC: zero_mean_parabolic_wave,
    WaveformType.MULTI_HARMONIC: zero_mean_multi_harmonic_wave,
    WaveformType.SMOOTHED_SQUARE: zero_mean_smoothed_square_wave,
}


def get_waveform(
    waveform: Union[WaveformType, str, Callable[[np.ndarray | float], np.ndarray]] = WaveformType.TRIANGLE,
    zero_mean: bool = False,
) -> Callable[[np.ndarray | float], np.ndarray]:
    """Retrieve or construct a normalized waveform generator.
    
    If zero_mean is True: returns function with range [-1, 1] and integral 0.
    If zero_mean is False: returns function with range [0, 1].
    """
    if callable(waveform):
        base_fn = waveform
        w_type = None
    elif isinstance(waveform, WaveformType):
        w_type = waveform
        base_fn = _ZERO_MEAN_REGISTRY[w_type]
    else:
        # String case: e.g. "triangle", "sine"
        val = str(waveform).lower().split(".")[-1]
        w_type = WaveformType(val)
        base_fn = _ZERO_MEAN_REGISTRY[w_type]

    if zero_mean:
        return base_fn

    # For triangle wave with zero_mean=False, preserve exact original formula
    if w_type == WaveformType.TRIANGLE:
        return triangle_wave

    # Transform zero-mean [-1, 1] to non-zero-mean [0, 1]
    def non_zero_mean_fn(u: np.ndarray | float) -> np.ndarray:
        return 0.5 + 0.5 * base_fn(u)

    return non_zero_mean_fn
