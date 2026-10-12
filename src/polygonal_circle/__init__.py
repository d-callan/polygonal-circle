"""polygonal-circle: recursive polygonal models of circles and fractal boundaries."""

from polygonal_circle.geometry import Polygon
from polygonal_circle.waves import triangle_wave
from polygonal_circle.model import PolygonalCircle, ScalingRegime, ScalingProperties
from polygonal_circle.analysis import PolygonalCircleAnalysis, LevelMetric
from polygonal_circle.annulus import (
    AnnularPolygonalCircle,
    AnnularConfinementReport,
    zero_mean_triangle_wave,
)

__version__ = "0.1.0"
__all__ = [
    "Polygon",
    "triangle_wave",
    "zero_mean_triangle_wave",
    "PolygonalCircle",
    "ScalingRegime",
    "ScalingProperties",
    "PolygonalCircleAnalysis",
    "LevelMetric",
    "AnnularPolygonalCircle",
    "AnnularConfinementReport",
]
