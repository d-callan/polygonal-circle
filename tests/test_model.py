"""Tests for the PolygonalCircle mathematical model."""

import numpy as np
import pytest

from polygonal_circle.model import PolygonalCircle, ScalingRegime


def test_invalid_parameters():
    with pytest.raises(ValueError, match="Radius"):
        PolygonalCircle(radius=-1.0)
    with pytest.raises(ValueError, match="eta"):
        PolygonalCircle(eta=1.5)
    with pytest.raises(ValueError, match="Shrinkage"):
        PolygonalCircle(a=1.2)
    with pytest.raises(ValueError, match="violate positivity"):
        # eta * a / (1 - a) = 0.8 * 0.6 / 0.4 = 1.2 >= 1.0
        PolygonalCircle(eta=0.8, a=0.6)


def test_circular_envelope_containment():
    pc = PolygonalCircle(radius=2.5, eta=0.15, a=0.5, b=3)
    thetas = np.linspace(0, 2 * np.pi, 500)

    for level in [0, 1, 2, 3, 4]:
        radii = pc.radius_at(thetas, level)
        # 1. Radius never exceeds envelope R
        assert np.all(radii <= 2.5 + 1e-12)
        # 2. Radius never drops below theoretical minimum R * (1 - eta * a / (1 - a))
        min_bound = 2.5 * (1.0 - 0.15 * (0.5 / 0.5))  # 2.5 * 0.85 = 2.125
        assert np.all(radii >= min_bound - 1e-12)


def test_polygon_convex_envelope():
    # Because all vertices have ||v|| <= R and the disk of radius R is convex,
    # all line segments between vertices must lie entirely within ||v|| <= R.
    pc = PolygonalCircle(radius=1.0, eta=0.1, a=0.5, b=3)
    poly = pc.polygon(level=3)

    assert poly.max_radius() <= 1.0 + 1e-12

    # Sample interior points along each edge
    v1 = poly.vertices
    v2 = np.roll(poly.vertices, -1, axis=0)
    for t in [0.25, 0.5, 0.75]:
        midpoints = (1.0 - t) * v1 + t * v2
        norms = np.linalg.norm(midpoints, axis=1)
        assert np.all(norms <= 1.0 + 1e-12)


def test_scaling_properties():
    # Convergent regime: a < 1/3
    pc_conv = PolygonalCircle(a=0.25, b=3)
    assert pc_conv.scaling.regime == ScalingRegime.CONVERGENT
    assert pc_conv.scaling.theoretical_dimension == 1.0

    # Critical regime: a = 1/3
    pc_crit = PolygonalCircle(a=1.0 / 3.0, b=3)
    assert pc_crit.scaling.regime == ScalingRegime.CRITICAL
    assert pc_crit.scaling.theoretical_dimension == 1.0

    # Divergent regime: a = 0.5 > 1/3
    pc_div = PolygonalCircle(a=0.5, b=3)
    assert pc_div.scaling.regime == ScalingRegime.DIVERGENT
    expected_dim = 1.0 + np.log(1.5) / np.log(3.0)
    assert np.isclose(pc_div.scaling.theoretical_dimension, expected_dim)
