"""Annular polygonal circle: polygons confined between two nested circular perimeters."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from polygonal_circle.geometry import Polygon
from polygonal_circle.model import ScalingProperties


def zero_mean_triangle_wave(u: np.ndarray | float, phase: float = 0.0) -> np.ndarray:
    """Evaluate periodic symmetric triangle wave with range [-1, 1] and zero mean.
    
    Formula:
        u_shifted = (u + phase) mod 1
        phi(u) = 1.0 - 4.0 * |u_shifted - 0.5|
    
    Properties:
        phi(0) = -1.0 (trough)
        phi(0.25) = 0.0 (zero-crossing)
        phi(0.5) = +1.0 (crest)
        phi(0.75) = 0.0 (zero-crossing)
        phi(1.0) = -1.0 (trough)
    """
    u_arr = np.asarray(u, dtype=np.float64)
    u_shifted = np.mod(u_arr + phase, 1.0)
    return 1.0 - 4.0 * np.abs(u_shifted - 0.5)


@dataclass(frozen=True)
class AnnularConfinementReport:
    """Report on geometric confinement of vertices and segments in the annulus."""
    level: int
    r_inner: float
    r_outer: float
    vertex_min_radius: float
    vertex_max_radius: float
    segment_min_distance: float
    segment_max_distance: float
    vertices_confined: bool
    segments_confined: bool

    @property
    def strictly_confined(self) -> bool:
        return self.vertices_confined and self.segments_confined


class AnnularPolygonalCircle:
    """Polygonal circle confined within an annular corridor [r_inner, r_outer].
    
    The radial boundary is centered at R0 = (r_inner + r_outer) / 2 with
    half-width delta = (r_outer - r_inner) / 2 ("wiggle room knob"):
    
        r_n(theta) = R0 + delta * psi_n(theta)
        
    where:
        psi_n(theta) = ((1 - a) / a) * sum_{k=1}^n a^k * phi_zero_mean(b^k * theta / (2*pi))
    
    Because |phi_zero_mean| <= 1, the normalized sum satisfies:
        |psi_n(theta)| <= 1 - a^n < 1
    
    Therefore, for ALL levels n >= 1 and ALL angles theta:
        r_inner < r_n(theta) < r_outer
        
    Attributes:
        r_inner: Inner circular boundary radius (> 0).
        r_outer: Outer circular boundary radius (> r_inner).
        a: Amplitude shrinkage per level in (0, 1).
        b: Spatial frequency multiplier per level (default 3).
        base_sides: Number of sides at level 0 (default 12 for tight annular chords).
    """

    def __init__(
        self,
        r_inner: float = 0.85,
        r_outer: float = 1.15,
        a: float = 0.5,
        b: int = 3,
        base_sides: int = 12,
    ):
        if r_inner <= 0:
            raise ValueError(f"Inner radius r_inner must be positive, got {r_inner}")
        if r_outer <= r_inner:
            raise ValueError(f"Outer radius r_outer ({r_outer}) must be > r_inner ({r_inner})")
        if not (0.0 < a < 1.0):
            raise ValueError(f"Shrinkage a must be in (0, 1), got {a}")
        if b < 2:
            raise ValueError(f"Frequency multiplier b must be >= 2, got {b}")
        if base_sides < 3:
            raise ValueError(f"Base polygon must have >= 3 sides, got {base_sides}")

        self.r_inner = float(r_inner)
        self.r_outer = float(r_outer)
        self.center_radius = 0.5 * (self.r_inner + self.r_outer)
        self.corridor_width = self.r_outer - self.r_inner
        self.half_width = 0.5 * self.corridor_width
        self.a = float(a)
        self.b = int(b)
        self.base_sides = int(base_sides)

    @classmethod
    def from_center_and_width(
        cls,
        center_radius: float = 1.0,
        corridor_width: float = 0.2,
        a: float = 0.5,
        b: int = 3,
        base_sides: int = 12,
    ) -> AnnularPolygonalCircle:
        """Create an AnnularPolygonalCircle from center radius R0 and corridor width W."""
        half_w = 0.5 * corridor_width
        if half_w >= center_radius:
            raise ValueError(
                f"Corridor width {corridor_width} too large for center radius {center_radius}; "
                "inner radius would be non-positive."
            )
        return cls(
            r_inner=center_radius - half_w,
            r_outer=center_radius + half_w,
            a=a,
            b=b,
            base_sides=base_sides,
        )

    @classmethod
    def from_relative_wiggle(
        cls,
        center_radius: float = 1.0,
        wiggle_fraction: float = 0.1,
        a: float = 0.5,
        b: int = 3,
        base_sides: int = 12,
    ) -> AnnularPolygonalCircle:
        """Create from center radius R0 and relative wiggle fraction w in (0, 1).
        
        r_inner = R0 * (1 - w), r_outer = R0 * (1 + w).
        """
        if not (0.0 < wiggle_fraction < 1.0):
            raise ValueError(f"Wiggle fraction must be in (0, 1), got {wiggle_fraction}")
        return cls(
            r_inner=center_radius * (1.0 - wiggle_fraction),
            r_outer=center_radius * (1.0 + wiggle_fraction),
            a=a,
            b=b,
            base_sides=base_sides,
        )

    @property
    def scaling(self) -> ScalingProperties:
        return ScalingProperties.from_parameters(self.a, self.b)

    def normalized_perturbation(self, theta: np.ndarray | float, level: int) -> np.ndarray:
        """Evaluate normalized recursive perturbation psi_n(theta) in (-1, 1)."""
        theta_arr = np.asarray(theta, dtype=np.float64)
        if level <= 0:
            return np.zeros_like(theta_arr)

        scale_factor = (1.0 - self.a) / self.a
        total = np.zeros_like(theta_arr)
        for k in range(1, level + 1):
            ak = self.a ** k
            freq = self.b ** k
            u = freq * theta_arr / (2.0 * np.pi)
            phi = zero_mean_triangle_wave(u)
            total += ak * phi

        return scale_factor * total

    def radius_at(self, theta: np.ndarray | float, level: int) -> np.ndarray:
        """Compute r_n(theta) = R0 + delta * psi_n(theta)."""
        theta_arr = np.asarray(theta, dtype=np.float64)
        psi = self.normalized_perturbation(theta_arr, level)
        return self.center_radius + self.half_width * psi

    def breakpoint_angles(self, level: int) -> np.ndarray:
        """Sample exact breakpoint angles for level n."""
        if level < 0:
            raise ValueError(f"Level must be >= 0, got {level}")
        num_vertices = self.base_sides * (self.b ** level)
        return np.linspace(0.0, 2.0 * np.pi, num_vertices, endpoint=False)

    def polygon(self, level: int) -> Polygon:
        """Construct the exact Polygon at recursion level n."""
        thetas = self.breakpoint_angles(level)
        radii = self.radius_at(thetas, level)
        x = radii * np.cos(thetas)
        y = radii * np.sin(thetas)
        return Polygon(np.column_stack([x, y]))

    def confinement_report(self, level: int) -> AnnularConfinementReport:
        """Examine whether vertices and all connecting straight segments stay inside [r_inner, r_outer]."""
        poly = self.polygon(level)
        v = poly.vertices
        n_v = len(v)

        # 1. Vertices check
        v_radii = poly.radii()
        v_min = float(np.min(v_radii))
        v_max = float(np.max(v_radii))
        tol = 1e-12
        v_confined = (v_min >= self.r_inner - tol) and (v_max <= self.r_outer + tol)

        # 2. Segments check
        # Segments maximum distance from origin is at endpoints (by convexity of norm)
        seg_max = v_max

        # Segments minimum distance from origin
        v_next = np.roll(v, -1, axis=0)
        edge = v_next - v
        edge_sq = np.sum(edge ** 2, axis=1)

        # Projection parameter t* of origin onto line v + t*edge
        t_star = -np.sum(v * edge, axis=1) / edge_sq

        # Clamped to segment [0, 1]
        t_clamped = np.clip(t_star, 0.0, 1.0)
        closest_pts = v + t_clamped[:, np.newaxis] * edge
        seg_dists = np.linalg.norm(closest_pts, axis=1)
        seg_min = float(np.min(seg_dists))

        seg_confined = (seg_min >= self.r_inner - tol) and (seg_max <= self.r_outer + tol)

        return AnnularConfinementReport(
            level=level,
            r_inner=self.r_inner,
            r_outer=self.r_outer,
            vertex_min_radius=v_min,
            vertex_max_radius=v_max,
            segment_min_distance=seg_min,
            segment_max_distance=seg_max,
            vertices_confined=v_confined,
            segments_confined=seg_confined,
        )
