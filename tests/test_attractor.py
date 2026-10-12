"""Tests for skew-product dynamical systems, Lyapunov spectrum, and strange attractors."""

import math
import numpy as np
import pytest

from polygonal_circle.attractor import (
    SkewProductCircleMap,
    compute_lyapunov_spectrum,
)


def test_lyapunov_spectrum_subcritical():
    """For a < 1/b: contraction dominates expansion, area-contracting, D_KY = 1."""
    spec = compute_lyapunov_spectrum(a=0.2, b=3)
    assert spec.is_chaotic is True
    assert spec.lambda_theta > 0
    assert np.isclose(spec.lambda_theta, math.log(3.0))

    assert spec.lambda_y < 0
    assert np.isclose(spec.lambda_y, math.log(0.2))

    assert spec.lyapunov_sum < 0
    assert spec.is_area_contracting is True
    assert spec.kaplan_yorke_dimension == 1.0


def test_lyapunov_spectrum_critical():
    """At a = 1/b: exact balance between chaotic expansion and contraction."""
    spec = compute_lyapunov_spectrum(a=1.0 / 3.0, b=3)
    assert np.isclose(spec.lyapunov_sum, 0.0, atol=1e-12)
    assert spec.kaplan_yorke_dimension == 1.0


def test_lyapunov_spectrum_supercritical():
    """For a > 1/b: chaotic expansion dominates, strange attractor with D_KY > 1."""
    spec = compute_lyapunov_spectrum(a=0.5, b=3)
    assert spec.lyapunov_sum > 0
    assert spec.is_area_contracting is False

    expected_d = 1.0 + math.log(1.5) / math.log(3.0)
    assert np.isclose(spec.kaplan_yorke_dimension, expected_d)
    assert spec.kaplan_yorke_dimension > 1.3


def test_skew_product_forward_step():
    sp = SkewProductCircleMap(a=0.5, b=3)
    t_next, y_next = sp.forward_step(theta=0.0, y=0.0)
    assert t_next == 0.0
    # At u = 0, zero-mean triangle wave is 0 - 4|0 - 0.5|?
    # Let's check zero_mean_triangle_wave(0.0): 1 - 4|0 - 0.5| = 1 - 4(0.5) = -1.0
    assert np.isclose(y_next, -1.0)


def test_invariant_graph_attractor():
    sp = SkewProductCircleMap(a=0.5, b=3)
    thetas, ys = sp.sample_invariant_attractor(n_points=50, depth=10)
    assert len(thetas) == 50
    assert len(ys) == 50
    # Because |phi| <= 1, |y*| <= sum_{k=1}^10 a^k = a(1 - a^10)/(1 - a) < 1.0 for a=0.5
    assert np.all(np.abs(ys) < 1.0)


def test_invalid_parameters():
    with pytest.raises(ValueError):
        compute_lyapunov_spectrum(a=1.5, b=3)
    with pytest.raises(ValueError):
        compute_lyapunov_spectrum(a=0.5, b=1)
