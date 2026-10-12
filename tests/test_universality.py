"""Tests for universality and waveform invariance across microscopic generators."""

import math
import numpy as np
import pytest

from polygonal_circle.model import PolygonalCircle
from polygonal_circle.universality import (
    analyze_waveform_universality,
    compute_multi_waveform_collapse,
    evaluate_waveform_integrals,
    get_all_waveform_types,
)
from polygonal_circle.waves import WaveformType, get_waveform


def test_waveform_periodicity_and_zero_mean():
    """Verify that all waveforms have period 1 and mean 0 (in zero-mean mode)."""
    u_samples = np.linspace(0.0, 1.0, 50000, endpoint=False)

    for w_type in get_all_waveform_types():
        fn = get_waveform(w_type, zero_mean=True)
        vals = fn(u_samples)

        # Mean should be very close to 0 (integral = 0)
        mean_val = float(np.mean(vals))
        assert abs(mean_val) < 1e-2, f"{w_type} failed zero-mean with mean={mean_val}"

        # Range should be within [-1.1, 1.1]
        assert np.min(vals) >= -1.05, f"{w_type} min too low: {np.min(vals)}"
        assert np.max(vals) <= 1.05, f"{w_type} max too high: {np.max(vals)}"

        # Periodicity check: phi(u + 1) == phi(u)
        test_u = np.array([0.1, 0.35, 0.72])
        np.testing.assert_allclose(fn(test_u), fn(test_u + 1.0), atol=1e-10)
        np.testing.assert_allclose(fn(test_u), fn(test_u + 2.0), atol=1e-10)


def test_universality_invariance_across_regimes():
    """Verify that all 5 waveforms share identical critical a_c, Euclidean limit, and fractal divergence."""
    summaries = analyze_waveform_universality(b=3, eta=0.1, max_level=7)
    assert len(summaries) == 5

    theoretical_d = 1.0 + math.log(1.5) / math.log(3.0)

    for s in summaries:
        # 1. Critical threshold is universally 1/3
        assert np.isclose(s.critical_a, 1.0 / 3.0)

        # 2. In convergent regime (a = 0.20), dimension is Euclidean D ~ 1.0
        assert abs(s.convergent_dimension_a02 - 1.0) < 0.05, (
            f"{s.waveform_name} convergent dimension {s.convergent_dimension_a02} deviates from 1.0"
        )

        # 3. In divergent regime (a = 0.50), dimension is strictly fractal (> 1.05) and <= theoretical ~ 1.369
        assert 1.05 <= s.divergent_dimension_a05 <= theoretical_d + 0.05, (
            f"{s.waveform_name} divergent dimension {s.divergent_dimension_a05} not in valid fractal range"
        )

        # 4. Growth ratio P_{n}/P_{n-1} is strictly expanding (> 1.05)
        assert s.perimeter_ratio_a05 > 1.05, (
            f"{s.waveform_name} ratio {s.perimeter_ratio_a05} not expanding"
        )


def test_multi_waveform_collapse():
    """Verify that multi-waveform collapse produces positive, non-empty scaling data."""
    data = compute_multi_waveform_collapse(a=0.50, b=3, eta=0.1, levels=[3, 4, 5, 6])
    assert len(data.waveforms) == 5
    for w in data.waveforms:
        vals = data.rescaled_excess[w]
        assert len(vals) == 4
        assert all(v > 0 for v in vals)
