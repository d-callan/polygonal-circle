"""Tests for ScaleSpaceObserver: proving Euclidean geometry as the Infrared observer limit."""

import math
import numpy as np
import pytest

from polygonal_circle.model import PolygonalCircle
from polygonal_circle.observer import ScaleSpaceObserver
from polygonal_circle.waves import WaveformType


def test_observer_ir_convergence_to_pi():
    """Verify that a macro-observer (sigma >= 1.0) measures Euclidean pi, Q=1, and D=1."""
    pc = PolygonalCircle(radius=1.0, eta=0.10, a=0.50, b=3)
    observer = ScaleSpaceObserver(pc, level=7, num_points=2048)

    obs = observer.observe(sigma=1.2)

    # 1. Observed pi converges to mathematical pi within 1e-4
    assert abs(obs.pi_obs - math.pi) < 1e-4, f"pi_obs {obs.pi_obs} deviated from pi"

    # 2. Isoperimetric quotient converges to 1.0000 (perfect circle)
    assert abs(obs.isoperimetric_quotient - 1.0) < 1e-4, f"Q {obs.isoperimetric_quotient} deviated from 1.0"

    # 3. Apparent dimension is exactly 1.000 (smooth curve)
    assert abs(obs.apparent_dimension - 1.0) < 0.02, f"D {obs.apparent_dimension} deviated from 1.0"


def test_observer_invariance_across_all_waveforms():
    """Verify that all 5 microscopic waveforms project onto the exact same Euclidean circle for macro-observers."""
    waveforms = [
        WaveformType.TRIANGLE,
        WaveformType.SINE,
        WaveformType.PARABOLIC,
        WaveformType.MULTI_HARMONIC,
        WaveformType.SMOOTHED_SQUARE,
    ]

    for w in waveforms:
        pc = PolygonalCircle(radius=1.0, eta=0.10, a=0.50, b=3, waveform=w)
        observer = ScaleSpaceObserver(pc, level=7, num_points=2048)
        obs = observer.observe(sigma=1.0)

        assert abs(obs.pi_obs - math.pi) < 5e-4, f"{w.value} failed pi convergence: {obs.pi_obs}"
        assert abs(obs.isoperimetric_quotient - 1.0) < 5e-4, f"{w.value} failed Q=1: {obs.isoperimetric_quotient}"


def test_observer_micro_fractal_crossover():
    """Verify that as observer resolution increases (sigma -> 0), fractal roughness emerges."""
    pc = PolygonalCircle(radius=1.0, eta=0.10, a=0.50, b=3)
    observer = ScaleSpaceObserver(pc, level=7, num_points=2048)

    obs_macro = observer.observe(sigma=1.0)
    obs_mid = observer.observe(sigma=0.05)
    obs_micro = observer.observe(sigma=0.005)

    # Perimeter and pi_obs strictly increase as resolution uncovers fractal facets
    assert obs_macro.perimeter < obs_mid.perimeter < obs_micro.perimeter
    assert obs_macro.pi_obs < obs_mid.pi_obs < obs_micro.pi_obs

    # Isoperimetric quotient strictly decreases away from Euclidean circle (Q=1)
    assert obs_macro.isoperimetric_quotient > obs_mid.isoperimetric_quotient > obs_micro.isoperimetric_quotient
    assert obs_micro.isoperimetric_quotient < 0.95

    # Apparent dimension increases from 1.000 towards fractal dimension
    assert obs_macro.apparent_dimension < obs_micro.apparent_dimension
    assert obs_micro.apparent_dimension > 1.02


def test_observer_filter_curve_output():
    """Verify filter_curve returns valid coordinate arrays."""
    pc = PolygonalCircle(radius=1.0, eta=0.10, a=0.50, b=3)
    observer = ScaleSpaceObserver(pc, level=5, num_points=1024)
    xs, ys = observer.filter_curve(sigma=0.1)

    assert len(xs) == 1024
    assert len(ys) == 1024
    assert not np.isnan(xs).any()
    assert not np.isnan(ys).any()
    assert not np.isinf(xs).any()
    assert not np.isinf(ys).any()
