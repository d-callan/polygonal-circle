"""Waveform generators for recursive radial perturbation."""

from __future__ import annotations

import numpy as np


def triangle_wave(u: np.ndarray | float, phase: float = 0.0) -> np.ndarray:
    """Evaluate periodic symmetric triangle wave with period 1 and range [0, 1].
    
    Formula:
        u_shifted = (u + phase) mod 1
        phi(u) = 1.0 - 2.0 * |u_shifted - 0.5|
    
    Properties:
        phi(0) = 0.0 (trough)
        phi(0.5) = 1.0 (crest)
        phi(1.0) = 0.0
        Slope is +2 on (0, 0.5) and -2 on (0.5, 1.0).
        Breakpoints occur at half-integers: u + phase in {0, 0.5, 1.0, 1.5, ...}.
    """
    u_arr = np.asarray(u, dtype=np.float64)
    u_shifted = np.mod(u_arr + phase, 1.0)
    return 1.0 - 2.0 * np.abs(u_shifted - 0.5)
