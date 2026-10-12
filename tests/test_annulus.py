"""Tests for AnnularPolygonalCircle and annular boundary confinement."""

import numpy as np
import pytest

from polygonal_circle.annulus import (
    AnnularPolygonalCircle,
    zero_mean_triangle_wave,
)


def test_zero_mean_triangle_wave():
    u = np.linspace(-3.0, 3.0, 600)
    phi = zero_mean_triangle_wave(u)

    assert np.all(phi >= -1.0 - 1e-12)
    assert np.all(phi <= 1.0 + 1e-12)

    # Check key values
    test_u = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    expected = np.array([-1.0, 0.0, 1.0, 0.0, -1.0])
    np.testing.assert_allclose(zero_mean_triangle_wave(test_u), expected, atol=1e-12)

    # Check mean over one period
    fine_u = np.linspace(0.0, 1.0, 10000, endpoint=False)
    assert np.isclose(np.mean(zero_mean_triangle_wave(fine_u)), 0.0, atol=1e-3)


def test_annulus_parameter_validation():
    with pytest.raises(ValueError, match="Inner radius"):
        AnnularPolygonalCircle(r_inner=-0.5, r_outer=1.0)
    with pytest.raises(ValueError, match="Outer radius"):
        AnnularPolygonalCircle(r_inner=1.2, r_outer=1.0)
    with pytest.raises(ValueError, match="too large"):
        AnnularPolygonalCircle.from_center_and_width(center_radius=1.0, corridor_width=2.5)


def test_vertex_confinement():
    # Test across various corridor widths
    for width in [0.40, 0.20, 0.05, 0.01]:
        ac = AnnularPolygonalCircle.from_center_and_width(
            center_radius=1.0, corridor_width=width, a=0.5, b=3, base_sides=12
        )
        thetas = np.linspace(0, 2 * np.pi, 500)

        for lvl in [0, 1, 2, 3, 4]:
            radii = ac.radius_at(thetas, lvl)
            assert np.all(radii >= ac.r_inner - 1e-12)
            assert np.all(radii <= ac.r_outer + 1e-12)


def test_segment_confinement():
    # With base_sides=12, at level >= 1, chords are short enough that segments stay confined
    ac = AnnularPolygonalCircle.from_center_and_width(
        center_radius=1.0, corridor_width=0.20, a=0.5, b=3, base_sides=12
    )

    for lvl in [1, 2, 3, 4]:
        report = ac.confinement_report(lvl)
        assert report.vertices_confined
        assert report.segments_confined
        assert report.strictly_confined


def test_wiggle_room_knob_limit():
    # As corridor width W -> 0, boundary converges to smooth circle of radius R0
    center_r = 1.0
    ac_narrow = AnnularPolygonalCircle.from_center_and_width(
        center_radius=center_r, corridor_width=1e-5, a=0.5, b=3, base_sides=12
    )
    poly = ac_narrow.polygon(level=3)

    # Perimeter should match regular polygon with 12 * 3^3 = 324 vertices on unit circle
    n_v = poly.num_vertices
    expected_p = 2.0 * n_v * center_r * np.sin(np.pi / n_v)
    assert np.isclose(poly.perimeter(), expected_p, rtol=1e-4)
    assert np.isclose(poly.mean_radius(), center_r, atol=1e-4)


def test_annulus_divergent_perimeter():
    # Fixed corridor width W=0.10 with a=0.5 > 1/3 yields growing perimeter
    ac = AnnularPolygonalCircle.from_center_and_width(
        center_radius=1.0, corridor_width=0.10, a=0.5, b=3, base_sides=12
    )
    perims = [ac.polygon(n).perimeter() for n in range(5)]

    # Perimeter grows as level increases
    assert perims[-1] > perims[0]
