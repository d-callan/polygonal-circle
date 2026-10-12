"""Universality and waveform invariance analysis for recursive polygonal circles.

This module demonstrates that the microscopic tooth geometry is an irrelevant operator:
Whether using:
1. Piecewise linear triangle waves (C^0)
2. Pure sinusoidal harmonics (C^infinity)
3. Piecewise quadratic parabolic waves (C^1)
4. Asymmetric multi-harmonic waves (C^infinity)
5. Smoothed hyperbolic tangent square waves (C^infinity)

The macroscopic geometric and statistical-physics properties are strictly universal:
- Critical threshold: a_c = 1 / b (identical across all waveforms)
- Order parameter exponent: beta = 1 (linear onset Psi ~ (a - a_c)^1 across all waveforms)
- Asymptotic fractal dimension: D = 1 + ln(ba) / ln(b) (identical across all waveforms)
- Infrared fixed point: Euclidean circle D = 1 for all a < a_c
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from polygonal_circle.analysis import PolygonalCircleAnalysis
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.renormalization import critical_a, order_parameter
from polygonal_circle.waves import WaveformType, get_waveform


@dataclass(frozen=True)
class WaveformUniversalSummary:
    """Summary of scaling and fractal properties for a specific waveform."""
    waveform_name: str
    continuity_class: str
    critical_a: float
    convergent_dimension_a02: float
    divergent_dimension_a05: float
    theoretical_dimension_a05: float
    perimeter_ratio_a05: float
    expected_perimeter_ratio: float


def get_all_waveform_types() -> list[WaveformType]:
    """Return all available standard microscopic waveforms."""
    return [
        WaveformType.TRIANGLE,
        WaveformType.SINE,
        WaveformType.PARABOLIC,
        WaveformType.MULTI_HARMONIC,
        WaveformType.SMOOTHED_SQUARE,
    ]


def evaluate_waveform_integrals(n_samples: int = 10000) -> dict[str, dict[str, float]]:
    """Compute empirical mean and bounds for all zero-mean waveforms to verify calibration."""
    u = np.linspace(0.0, 1.0, n_samples, endpoint=False)
    results = {}
    for w_type in get_all_waveform_types():
        fn = get_waveform(w_type, zero_mean=True)
        vals = fn(u)
        results[w_type.value] = {
            "mean": float(np.mean(vals)),
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
        }
    return results


def analyze_waveform_universality(
    b: int = 3,
    eta: float = 0.1,
    max_level: int = 8,
) -> list[WaveformUniversalSummary]:
    """Compute empirical dimensions and perimeter growth ratios across all waveforms."""
    a_c = critical_a(b)
    expected_ratio = float(b) * 0.50
    theoretical_d = 1.0 + float(np.log(b * 0.50) / np.log(b))

    continuity_map = {
        WaveformType.TRIANGLE: "C^0 (Piecewise Linear)",
        WaveformType.SINE: "C^inf (Harmonic)",
        WaveformType.PARABOLIC: "C^1 (Piecewise Quadratic)",
        WaveformType.MULTI_HARMONIC: "C^inf (Asymmetric Sum)",
        WaveformType.SMOOTHED_SQUARE: "C^inf (Steep Tanh)",
    }

    summaries = []
    for w_type in get_all_waveform_types():
        # 1. Convergent regime (a = 0.20 < a_c = 1/3)
        pc_conv = PolygonalCircle(radius=1.0, eta=eta, a=0.20, b=b, waveform=w_type)
        an_conv = PolygonalCircleAnalysis(pc_conv)
        dim_conv, _, _ = an_conv.richardson_dimension(max_level=max_level, min_level=max_level - 3)

        # 2. Divergent regime (a = 0.50 > a_c = 1/3)
        pc_div = PolygonalCircle(radius=1.0, eta=eta, a=0.50, b=b, waveform=w_type)
        an_div = PolygonalCircleAnalysis(pc_div)
        dim_div, _, _ = an_div.richardson_dimension(max_level=max_level, min_level=max_level - 3)

        # Growth ratio P_n / P_{n-1} at highest level
        metrics_div = an_div.evaluate_levels(max_level)
        p_prev = metrics_div[-2].perimeter
        p_curr = metrics_div[-1].perimeter
        ratio = p_curr / p_prev

        summaries.append(
            WaveformUniversalSummary(
                waveform_name=w_type.value,
                continuity_class=continuity_map[w_type],
                critical_a=a_c,
                convergent_dimension_a02=dim_conv,
                divergent_dimension_a05=dim_div,
                theoretical_dimension_a05=theoretical_d,
                perimeter_ratio_a05=ratio,
                expected_perimeter_ratio=expected_ratio,
            )
        )

    return summaries


@dataclass(frozen=True)
class MultiWaveformCollapseData:
    """Universal scaling collapse data aggregated across all 5 waveforms."""
    levels: list[int]
    waveforms: list[str]
    rescaled_excess: dict[str, list[float]]


def compute_multi_waveform_collapse(
    a: float = 0.50,
    b: int = 3,
    eta: float = 0.1,
    levels: list[int] = [3, 4, 5, 6, 7],
) -> MultiWaveformCollapseData:
    """Evaluate rescaled excess perimeter (P_n - P_0) * b^{-n * Psi(a)} for all waveforms."""
    psi = order_parameter(a, b)
    rescaled = {}

    for w_type in get_all_waveform_types():
        pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=b, waveform=w_type)
        p0 = pc.polygon(0).perimeter()
        y_vals = []
        for lvl in levels:
            pn = pc.polygon(lvl).perimeter()
            excess = pn - p0
            # Universal scaling coordinate: excess * b^{-n * Psi}
            y_scaled = excess / (float(b) ** (lvl * psi))
            y_vals.append(float(y_scaled))
        rescaled[w_type.value] = y_vals

    return MultiWaveformCollapseData(
        levels=levels,
        waveforms=[w.value for w in get_all_waveform_types()],
        rescaled_excess=rescaled,
    )
