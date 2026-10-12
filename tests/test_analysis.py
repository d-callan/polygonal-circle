"""Tests for multi-level geometric and fractal analysis."""

import numpy as np
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.analysis import PolygonalCircleAnalysis


def test_area_boundedness():
    pc = PolygonalCircle(radius=1.0, eta=0.1, a=0.5, b=3)
    analysis = PolygonalCircleAnalysis(pc)
    df = analysis.metrics_dataframe(max_level=5)

    # Area must be strictly less than pi * R^2 = pi
    assert np.all(df["area"] < np.pi)
    # Area must be strictly positive
    assert np.all(df["area"] > 2.0)


def test_scaling_regimes_behavior():
    # Convergent regime (a = 0.2 < 1/3): perimeter ratios approach 1.0
    pc_conv = PolygonalCircle(a=0.2, eta=0.1, b=3)
    df_conv = PolygonalCircleAnalysis(pc_conv).metrics_dataframe(max_level=5)
    last_ratio_conv = df_conv["perimeter_ratio"].iloc[-1]
    assert np.isclose(last_ratio_conv, 1.0, atol=1e-4)

    # Divergent regime (a = 0.5 > 1/3): perimeters increase steadily
    pc_div = PolygonalCircle(a=0.5, eta=0.1, b=3)
    df_div = PolygonalCircleAnalysis(pc_div).metrics_dataframe(max_level=7)
    # At level 7, ratio has grown above 1.14 towards 1.5
    assert df_div["perimeter_ratio"].iloc[-1] > 1.14
    assert df_div["perimeter"].iloc[-1] > df_div["perimeter"].iloc[0]


def test_richardson_dimension():
    # At coarse scales, Euclidean circle perimeter dominates (dim ~ 1.0).
    # As resolution refines (levels 6-10), fractal deviations dominate and
    # empirical dimension approaches theoretical D = 1 + log(1.5)/log(3) ~ 1.369.
    pc = PolygonalCircle(a=0.5, eta=0.1, b=3)
    analysis = PolygonalCircleAnalysis(pc)
    dim, r2, slope = analysis.richardson_dimension(max_level=10, min_level=6)

    assert r2 > 0.95
    # The empirical dimension approaches theoretical ~1.369
    assert 1.25 <= dim <= 1.40


def test_box_counting_dimension():
    pc = PolygonalCircle(a=0.5, eta=0.1, b=3)
    analysis = PolygonalCircleAnalysis(pc)
    dim, r2, deltas, counts = analysis.box_counting_dimension(level=4, num_scales=6)

    assert r2 > 0.95
    assert 1.0 <= dim <= 1.8
    assert len(counts) == 6
