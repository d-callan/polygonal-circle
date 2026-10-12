"""Tests for ParametricCarriers and ScaleSpaceObserver projection onto arbitrary Euclidean curves."""

import math
import numpy as np
import pytest

from polygonal_circle.carriers import (
    CircleCarrier,
    EllipseCarrier,
    RegularPolygonCarrier,
    SquareCarrier,
    SuperellipseCarrier,
    TrefoilCarrier,
    TriangleCarrier,
)
from polygonal_circle.observer import ScaleSpaceObserver


def test_circle_carrier():
    c = CircleCarrier(radius=1.0)
    assert c.name == "circle"
    assert math.isclose(c.exact_perimeter(), 2.0 * math.pi, rel_tol=1e-7)
    assert math.isclose(c.exact_area(), math.pi, rel_tol=1e-7)
    assert math.isclose(c.isoperimetric_quotient(), 1.0, rel_tol=1e-7)


def test_ellipse_carrier():
    e = EllipseCarrier(semi_major=1.25, semi_minor=0.80)
    assert e.name == "ellipse"
    assert math.isclose(e.exact_area(), math.pi * 1.25 * 0.80, rel_tol=1e-7)

    # Ramanujan perimeter matches high-precision numerical quadrature within 1e-5
    thetas = np.linspace(0, 2 * math.pi, 16384, endpoint=False)
    r = e.radius_at(thetas)
    xs = r * np.cos(thetas)
    ys = r * np.sin(thetas)
    p_num = float(np.sum(np.sqrt((np.roll(xs, -1) - xs) ** 2 + (np.roll(ys, -1) - ys) ** 2)))

    assert math.isclose(e.exact_perimeter(), p_num, rel_tol=1e-5)
    assert e.isoperimetric_quotient() < 1.0


def test_superellipse_carrier():
    s = SuperellipseCarrier(semi_major=1.0, semi_minor=1.0, exponent=4.0)
    assert s.name == "superellipse"

    # Gamma-function squircle area
    expected_area = 4.0 * (math.gamma(1.25) ** 2) / math.gamma(1.5)
    assert math.isclose(s.exact_area(), expected_area, rel_tol=1e-5)
    assert s.isoperimetric_quotient() < 1.0


def test_trefoil_carrier():
    t = TrefoilCarrier(base_radius=1.0, modulation=0.15)
    assert t.name == "trefoil"
    expected_area = math.pi * (1.0**2) * (1.0 + 0.5 * (0.15**2))
    assert math.isclose(t.exact_area(), expected_area, rel_tol=1e-6)


def test_square_carrier_exact():
    sq = SquareCarrier(inradius=1.0)
    assert sq.name == "square"
    assert math.isclose(sq.exact_perimeter(), 8.0, rel_tol=1e-7)
    assert math.isclose(sq.exact_area(), 4.0, rel_tol=1e-7)
    assert math.isclose(sq.isoperimetric_quotient(), math.pi / 4.0, rel_tol=1e-7)


def test_triangle_carrier_exact():
    tri = TriangleCarrier(inradius=1.0)
    assert tri.name == "triangle"
    expected_p = 6.0 * math.sqrt(3.0)
    expected_a = 3.0 * math.sqrt(3.0)
    expected_q = math.pi * math.sqrt(3.0) / 9.0
    assert math.isclose(tri.exact_perimeter(), expected_p, rel_tol=1e-7)
    assert math.isclose(tri.exact_area(), expected_a, rel_tol=1e-7)
    assert math.isclose(tri.isoperimetric_quotient(), expected_q, rel_tol=1e-7)


def test_polygon_corner_rounding_and_optimal_window():
    """Verify corner rounding at large sigma and fractal divergence at small sigma for square & triangle."""
    for poly in [SquareCarrier(inradius=1.0), TriangleCarrier(inradius=1.0)]:
        obs = ScaleSpaceObserver.from_carrier(
            carrier=poly,
            eta=0.10,
            a=0.50,
            b=3,
            level=7,
            num_points=4096,
        )

        p_exact = poly.exact_perimeter()

        # Macro: corners rounded -> P < P_exact
        macro = obs.observe(sigma=0.20)
        assert macro.perimeter < p_exact, f"{poly.name} expected corner deficit at sigma=0.20"

        # Optimal observer window: sigma ~ 0.05 rad recovers Euclidean perimeter to < 2.5%
        meso = obs.observe(sigma=0.05)
        rel_err = abs(meso.perimeter - p_exact) / p_exact
        assert rel_err < 0.025, f"{poly.name} failed recovery in optimal window: err={rel_err:.4e}"

        # Micro: fractal teeth -> P > P_exact
        micro = obs.observe(sigma=0.005)
        assert micro.perimeter > p_exact, f"{poly.name} expected fractal divergence at sigma=0.005"
        assert micro.perimeter > meso.perimeter


@pytest.mark.parametrize(
    "carrier",
    [
        CircleCarrier(radius=1.0),
        EllipseCarrier(semi_major=1.25, semi_minor=0.80),
        SuperellipseCarrier(semi_major=1.0, semi_minor=1.0, exponent=4.0),
        TrefoilCarrier(base_radius=1.0, modulation=0.15),
    ],
)
def test_scale_space_observer_recovers_arbitrary_carrier(carrier):
    """Verify that any arbitrary macroscopic Euclidean carrier is recovered in the observer window."""
    obs = ScaleSpaceObserver.from_carrier(
        carrier=carrier,
        eta=0.10,
        a=0.50,
        b=3,
        level=7,
        waveform="triangle",
        num_points=2048,
    )

    # In the intermediate scale-space observer window (lambda_ripple << sigma << lambda_carrier)
    # sigma ~ 0.2 rad filters out fractal ripples while preserving the carrier shape
    meso = obs.observe(sigma=0.20)
    micro = obs.observe(sigma=0.005)

    # 1. Observer in the scale-space window recovers carrier's theoretical perimeter within 1.5%
    assert meso.perimeter_relative_error < 0.015, (
        f"{carrier.name} failed perimeter recovery: "
        f"measured={meso.perimeter:.4f}, target={carrier.exact_perimeter():.4f}, "
        f"err={meso.perimeter_relative_error:.4e}"
    )

    # 2. In the observer window, apparent dimension is close to Euclidean D = 1.000
    assert abs(meso.apparent_dimension - 1.0) < 0.05, (
        f"{carrier.name} failed Euclidean D=1: {meso.apparent_dimension}"
    )

    # 3. Observer recovers carrier's isoperimetric quotient within 2.5%
    target_q = carrier.isoperimetric_quotient()
    assert abs(meso.isoperimetric_quotient - target_q) / target_q < 0.025

    # 4. Micro-observer perceives fractal roughening (perimeter strictly greater than in observer window)
    assert micro.perimeter > meso.perimeter
    assert micro.apparent_dimension > 1.02
