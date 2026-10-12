"""Recursive polygonal circle mathematical model."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Union
import numpy as np

from polygonal_circle.geometry import Polygon
from polygonal_circle.waves import WaveformType, get_waveform, triangle_wave


class ScalingRegime(str, Enum):
    """Regime of perimeter scaling as recursion level n -> infinity."""
    CONVERGENT = "convergent"      # a < 1/b: perimeter converges to finite limit
    CRITICAL = "critical"          # a = 1/b: perimeter increments stay roughly constant
    DIVERGENT = "divergent"        # a > 1/b: perimeter diverges geometrically


@dataclass(frozen=True)
class ScalingProperties:
    """Theoretical scaling properties for given parameters (a, b)."""
    a: float
    b: int
    product_ab: float
    regime: ScalingRegime
    critical_a: float
    theoretical_dimension: float

    @classmethod
    def from_parameters(cls, a: float, b: int = 3) -> ScalingProperties:
        crit = 1.0 / b
        prod = a * b
        tol = 1e-9
        if abs(a - crit) < tol:
            regime = ScalingRegime.CRITICAL
            dim = 1.0
        elif a < crit:
            regime = ScalingRegime.CONVERGENT
            dim = 1.0
        else:
            regime = ScalingRegime.DIVERGENT
            # Richardson / box-counting dimension: D = 1 + log(a*b)/log(b) = 2 - log(1/a)/log(b)
            dim = 1.0 + float(np.log(prod) / np.log(b))

        return cls(
            a=a,
            b=b,
            product_ab=prod,
            regime=regime,
            critical_a=crit,
            theoretical_dimension=dim,
        )


class PolygonalCircle:
    """Model of a recursive polygonal circle.
    
    The radial boundary is defined by:
        r_n(theta) = R * [ 1 - eta * sum_{k=1}^n a^k * phi(b^k * theta / (2*pi)) ]
    
    where phi is a repeating triangular wave with range [0, 1].
    
    Attributes:
        radius: Reference circle envelope radius R (default 1.0).
        eta: Irregularity amplitude (default 0.1). Must satisfy eta * a / (1 - a) < 1.
        a: Amplitude shrinkage per level (default 0.5).
        b: Spatial frequency multiplier per level (default 3).
        base_sides: Number of sides for level 0 inscribed polygon (default 6).
    """

    def __init__(
        self,
        radius: float = 1.0,
        eta: float = 0.1,
        a: float = 0.5,
        b: int = 3,
        base_sides: int = 6,
        waveform: Union[WaveformType, str, Callable[[np.ndarray | float], np.ndarray]] = WaveformType.TRIANGLE,
    ):
        if radius <= 0:
            raise ValueError(f"Radius R must be positive, got {radius}")
        if not (0.0 <= eta < 1.0):
            raise ValueError(f"Amplitude eta must be in [0, 1), got {eta}")
        if not (0.0 < a < 1.0):
            raise ValueError(f"Shrinkage a must be in (0, 1), got {a}")
        if b < 2:
            raise ValueError(f"Frequency multiplier b must be integer >= 2, got {b}")
        if base_sides < 3:
            raise ValueError(f"Base polygon must have >= 3 sides, got {base_sides}")

        # Geometric series bound check for positive radius
        max_radial_reduction = eta * a / (1.0 - a)
        if max_radial_reduction >= 1.0:
            raise ValueError(
                f"Parameters violate positivity: eta * a / (1 - a) = {max_radial_reduction:.4f} >= 1.0. "
                "Radius would become non-positive for high levels."
            )

        self.radius = float(radius)
        self.eta = float(eta)
        self.a = float(a)
        self.b = int(b)
        self.base_sides = int(base_sides)
        self.waveform = waveform
        self.phi_fn = get_waveform(waveform, zero_mean=False)

    @property
    def scaling(self) -> ScalingProperties:
        return ScalingProperties.from_parameters(self.a, self.b)

    def radial_perturbation_sum(self, theta: np.ndarray | float, level: int) -> np.ndarray:
        """Compute the perturbation sum sum_{k=1}^n a^k * phi(b^k * theta / (2*pi))."""
        theta_arr = np.asarray(theta, dtype=np.float64)
        if level <= 0:
            return np.zeros_like(theta_arr)

        total = np.zeros_like(theta_arr)
        for k in range(1, level + 1):
            ak = self.a ** k
            freq = self.b ** k
            u = freq * theta_arr / (2.0 * np.pi)
            phi = self.phi_fn(u)
            total += ak * phi
        return total

    def radius_at(self, theta: np.ndarray | float, level: int) -> np.ndarray:
        """Compute r_n(theta) = R * [ 1 - eta * sum_{k=1}^n a^k * phi(b^k * theta / (2*pi)) ]."""
        theta_arr = np.asarray(theta, dtype=np.float64)
        pert = self.radial_perturbation_sum(theta_arr, level)
        return self.radius * (1.0 - self.eta * pert)

    def breakpoint_angles(self, level: int) -> np.ndarray:
        """Return the exact breakpoint angles for level n.
        
        At level n, the finest frequency wave b^n has breakpoints spaced at
        theta = 2*pi*j / (base_sides * b^n).
        Since b^k divides base_sides * b^n for all k <= n, this grid contains
        all kinks of all constituent waves k = 1, ..., n.
        """
        if level < 0:
            raise ValueError(f"Level must be non-negative, got {level}")
        num_vertices = self.base_sides * (self.b ** level)
        return np.linspace(0.0, 2.0 * np.pi, num_vertices, endpoint=False)

    def polygon(self, level: int) -> Polygon:
        """Construct the exact Polygon at recursion level n.
        
        Vertices are computed at all exact breakpoints:
            (x_j, y_j) = (r_n(theta_j) * cos(theta_j), r_n(theta_j) * sin(theta_j))
        """
        thetas = self.breakpoint_angles(level)
        radii = self.radius_at(thetas, level)
        x = radii * np.cos(thetas)
        y = radii * np.sin(thetas)
        verts = np.column_stack([x, y])
        return Polygon(verts)

    def sample_curve(self, level: int, num_points: int = 2000) -> tuple[np.ndarray, np.ndarray]:
        """Dense continuous sampling in polar coordinates for visualization."""
        thetas = np.linspace(0.0, 2.0 * np.pi, num_points, endpoint=False)
        radii = self.radius_at(thetas, level)
        x = radii * np.cos(thetas)
        y = radii * np.sin(thetas)
        return x, y
