"""Tests for exact polygon geometry computations."""

import numpy as np
import pytest

from polygonal_circle.geometry import Polygon


def test_square_geometry():
    # Square with side 2 centered at origin, area = 4, perimeter = 8
    verts = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    poly = Polygon(verts)

    assert poly.num_vertices == 4
    assert np.isclose(poly.perimeter(), 8.0)
    assert np.isclose(poly.area(), 4.0)
    # Isoperimetric quotient 4 * pi * A / P^2 = 4 * pi * 4 / 64 = pi / 4
    assert np.isclose(poly.isoperimetric_quotient(), np.pi / 4.0)
    np.testing.assert_allclose(poly.centroid(), [0.0, 0.0], atol=1e-12)


def test_regular_hexagon_geometry():
    # Regular hexagon inscribed in unit circle: side = 1, perimeter = 6, area = 3*sqrt(3)/2
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    verts = np.column_stack([np.cos(angles), np.sin(angles)])
    poly = Polygon(verts)

    assert poly.num_vertices == 6
    assert np.isclose(poly.perimeter(), 6.0)
    expected_area = 1.5 * np.sqrt(3.0)
    assert np.isclose(poly.area(), expected_area)
    assert poly.max_radius() <= 1.0 + 1e-12
    assert poly.min_radius() >= 1.0 - 1e-12


def test_invalid_polygon():
    with pytest.raises(ValueError):
        Polygon([[0, 0], [1, 1]])  # < 3 vertices
