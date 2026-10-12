"""polygonal-circle: recursive polygonal models of circles and fractal boundaries."""

from polygonal_circle.geometry import Polygon
from polygonal_circle.carriers import (
    CircleCarrier,
    EllipseCarrier,
    ParametricCarrier,
    SuperellipseCarrier,
    TrefoilCarrier,
)
from polygonal_circle.waves import (
    WaveformType,
    get_waveform,
    triangle_wave,
    zero_mean_triangle_wave,
)
from polygonal_circle.model import PolygonalCircle, ScalingRegime, ScalingProperties
from polygonal_circle.analysis import PolygonalCircleAnalysis, LevelMetric
from polygonal_circle.annulus import (
    AnnularPolygonalCircle,
    AnnularConfinementReport,
)
from polygonal_circle.renormalization import (
    OperatorRelevance,
    PhaseState,
    classify_phase,
    critical_a,
    order_parameter,
    order_parameter_derivative,
    rg_flow_trajectory,
)
from polygonal_circle.attractor import (
    LyapunovSpectrum,
    SkewProductCircleMap,
    compute_lyapunov_spectrum,
)
from polygonal_circle.observer import (
    ObserverObservation,
    ScaleSpaceObserver,
)
from polygonal_circle.universality import (
    WaveformUniversalSummary,
    analyze_waveform_universality,
    compute_multi_waveform_collapse,
    get_all_waveform_types,
)

__version__ = "0.1.0"
__all__ = [
    "Polygon",
    "ParametricCarrier",
    "CircleCarrier",
    "EllipseCarrier",
    "SuperellipseCarrier",
    "TrefoilCarrier",
    "WaveformType",
    "get_waveform",
    "triangle_wave",
    "zero_mean_triangle_wave",
    "PolygonalCircle",
    "ScalingRegime",
    "ScalingProperties",
    "PolygonalCircleAnalysis",
    "LevelMetric",
    "AnnularPolygonalCircle",
    "AnnularConfinementReport",
    "OperatorRelevance",
    "PhaseState",
    "classify_phase",
    "critical_a",
    "order_parameter",
    "order_parameter_derivative",
    "rg_flow_trajectory",
    "LyapunovSpectrum",
    "SkewProductCircleMap",
    "compute_lyapunov_spectrum",
    "ObserverObservation",
    "ScaleSpaceObserver",
    "WaveformUniversalSummary",
    "analyze_waveform_universality",
    "compute_multi_waveform_collapse",
    "get_all_waveform_types",
]
