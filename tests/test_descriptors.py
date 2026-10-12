"""Tests for Fourier descriptor ingestion and non-star-shaped carriers."""

import math
import numpy as np
import pytest

from polygonal_circle.carriers import EllipseCarrier, SquareCarrier, TrefoilCarrier
from polygonal_circle.descriptors import (
    CrescentCarrier,
    FourierDescriptorCarrier,
    KidneyBeanCarrier,
)
from polygonal_circle.planar import (
    PlanarScaleSpaceObserver,
    RadialCarrierAdapter,
)


def _carrier_points(carrier, n=2048):
    """Sample raw (x, y) boundary points from a radial carrier."""
    u = np.linspace(0.0, 1.0, n, endpoint=False)
    adapted = RadialCarrierAdapter(carrier)
    return adapted.point_at(u)


def test_descriptor_closed_form_area_circle():
    """Descriptor of a circle: only k=1 survives; A = pi(a1 d1 - b1 c1) = pi."""
    n = 4096
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs, ys = np.cos(t), np.sin(t)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=8)
    assert math.isclose(fd.exact_area(), math.pi, rel_tol=1e-6)
    assert math.isclose(fd.exact_perimeter(), 2 * math.pi, rel_tol=1e-3)


def test_descriptor_roundtrip_square():
    """Point cloud of a square -> descriptor recovers P and A within 2%."""
    sq = SquareCarrier(inradius=1.0)
    xs, ys = _carrier_points(sq, 4096)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=48)
    assert math.isclose(fd.exact_area(), sq.exact_area(), rel_tol=0.02)
    assert math.isclose(fd.exact_perimeter(), sq.exact_perimeter(), rel_tol=0.02)


def test_descriptor_roundtrip_trefoil():
    tr = TrefoilCarrier(base_radius=1.0, modulation=0.15)
    xs, ys = _carrier_points(tr, 4096)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=24)
    assert math.isclose(fd.exact_area(), tr.exact_area(), rel_tol=0.02)


def test_descriptor_closed_form_area_matches_shoelace():
    """Coefficient-identity area agrees with shoelace on the reconstructed curve."""
    xs, ys = _carrier_points(EllipseCarrier(semi_major=1.25, semi_minor=0.80), 4096)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=24)
    curve = fd.sample(8192)
    assert math.isclose(fd.exact_area(), curve.area, rel_tol=1e-3)


def test_crescent_exact_area_invariant_under_dent():
    """The crescent's area is exactly pi * scale^2 for every dent amplitude."""
    for d in (0.0, 0.5, 1.2):
        c = CrescentCarrier(scale=1.0, dent=d)
        assert math.isclose(c.exact_area(), math.pi, rel_tol=1e-12)
        # numerical check
        assert math.isclose(c.numerical_area(16384), math.pi, rel_tol=1e-4)


def test_kidney_bean_closed_form_area():
    kb = KidneyBeanCarrier()
    assert math.isclose(kb.exact_area(), kb.numerical_area(16384), rel_tol=1e-4)


def test_descriptor_carrier_observation_window():
    """Observer recovers a descriptor-ingested carrier's P and A."""
    sq = SquareCarrier(inradius=1.0)
    xs, ys = _carrier_points(sq, 4096)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=48)
    obs = PlanarScaleSpaceObserver(fd, eta=0.05, a=0.5, b=3, level=6)
    sweep = obs.sweep_resolutions(np.logspace(-3.0, -0.2, 16))
    best_p = min(o.perimeter_relative_error for o in sweep)
    best_a = min(o.area_relative_error for o in sweep)
    assert best_p < 0.025
    assert best_a < 0.025


def test_point_at_ccw():
    """Descriptor reconstruction traverses counter-clockwise (positive signed area)."""
    c = CrescentCarrier(scale=1.0, dent=1.2)
    u = np.linspace(0, 1, 1024, endpoint=False)
    xs, ys = c.point_at(u)
    fd = FourierDescriptorCarrier.from_points(xs, ys, num_modes=16)
    xs2, ys2 = fd.point_at(u)
    from polygonal_circle.planar import signed_area

    assert signed_area(xs2, ys2) > 0
