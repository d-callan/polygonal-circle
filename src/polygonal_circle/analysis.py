"""Analysis routines: perimeters, areas, scaling laws, Richardson and box-counting dimensions."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy import stats

from polygonal_circle.geometry import Polygon
from polygonal_circle.model import PolygonalCircle


@dataclass
class LevelMetric:
    """Metrics evaluated for a specific recursion level n."""
    level: int
    num_vertices: int
    perimeter: float
    perimeter_ratio: float | None
    area: float
    isoperimetric_quotient: float
    max_radius: float
    min_radius: float
    mean_radius: float
    radial_amplitude: float
    mean_segment_length: float
    total_turn_angle: float


class PolygonalCircleAnalysis:
    """Multi-level geometric and fractal analysis of a PolygonalCircle."""

    def __init__(self, circle: PolygonalCircle):
        self.circle = circle

    def evaluate_levels(self, max_level: int = 5) -> list[LevelMetric]:
        """Compute comprehensive geometric metrics from level 0 up to max_level."""
        metrics: list[LevelMetric] = []
        prev_perimeter: float | None = None

        for n in range(max_level + 1):
            poly = self.circle.polygon(n)
            p = poly.perimeter()
            ratio = (p / prev_perimeter) if prev_perimeter is not None else None
            prev_perimeter = p

            seg_lens = poly.segment_lengths()
            mean_seg = float(np.mean(seg_lens))

            # Compute turning angles between consecutive edges
            edge_vecs = poly.edge_vectors()
            # Normalize vectors to unit directions
            unit_edges = edge_vecs / seg_lens[:, np.newaxis]
            next_unit_edges = np.roll(unit_edges, -1, axis=0)
            # Dot products clamped to [-1, 1] for acos
            dots = np.clip(np.sum(unit_edges * next_unit_edges, axis=1), -1.0, 1.0)
            turn_angles = np.arccos(dots)
            total_turn = float(np.sum(np.abs(turn_angles)))

            r_max = poly.max_radius()
            r_min = poly.min_radius()

            metrics.append(
                LevelMetric(
                    level=n,
                    num_vertices=poly.num_vertices,
                    perimeter=p,
                    perimeter_ratio=ratio,
                    area=poly.area(),
                    isoperimetric_quotient=poly.isoperimetric_quotient(),
                    max_radius=r_max,
                    min_radius=r_min,
                    mean_radius=poly.mean_radius(),
                    radial_amplitude=r_max - r_min,
                    mean_segment_length=mean_seg,
                    total_turn_angle=total_turn,
                )
            )

        return metrics

    def metrics_dataframe(self, max_level: int = 5) -> pd.DataFrame:
        """Return metrics across levels as a clean pandas DataFrame."""
        records = [m.__dict__ for m in self.evaluate_levels(max_level)]
        return pd.DataFrame(records)

    def richardson_dimension(
        self, max_level: int = 6, min_level: int = 1
    ) -> tuple[float, float, float]:
        """Estimate fractal dimension via Richardson coastline scaling.
        
        Fits log(Perimeter) = (1 - D) * log(mean_segment_length) + const.
        Returns:
            dimension: D = 1 - slope
            r_squared: coefficient of determination R^2
            slope: fitted slope
        """
        metrics = self.evaluate_levels(max_level)
        sub = [m for m in metrics if min_level <= m.level <= max_level]
        if len(sub) < 2:
            return 1.0, 1.0, 0.0

        log_eps = np.log([m.mean_segment_length for m in sub])
        log_p = np.log([m.perimeter for m in sub])

        res = stats.linregress(log_eps, log_p)
        slope = float(res.slope)
        r2 = float(res.rvalue ** 2)
        dim = 1.0 - slope
        return dim, r2, slope

    def box_counting_dimension(
        self, level: int = 5, num_scales: int = 10
    ) -> tuple[float, float, np.ndarray, np.ndarray]:
        """Compute 2D box-counting dimension of polygon boundary at given level.
        
        Subdivides the bounding box into grid cells of size delta and counts
        the number of boxes N(delta) traversed by the polygon's boundary segments.
        
        Returns:
            dimension: D_box = -slope of log(N) vs log(delta)
            r_squared: R^2 of linear fit
            deltas: array of box sizes
            box_counts: array of occupied box counts N(delta)
        """
        poly = self.circle.polygon(level)
        verts = poly.vertices
        n_v = len(verts)

        # Scale range: from large fraction of radius down to sub-segment resolution
        seg_lens = poly.segment_lengths()
        delta_max = self.circle.radius * 0.5
        delta_min = max(float(np.min(seg_lens)) * 0.25, 1e-4)

        deltas = np.geomspace(delta_max, delta_min, num_scales)
        counts = []

        # Prepare segment endpoints
        p1 = verts
        p2 = np.roll(verts, -1, axis=0)

        for delta in deltas:
            occupied_boxes: set[tuple[int, int]] = set()

            # For each segment, sample points along it at spacing <= delta / 2
            for i in range(n_v):
                v1 = p1[i]
                v2 = p2[i]
                length = seg_lens[i]
                num_steps = max(int(np.ceil(2.0 * length / delta)), 1)
                t = np.linspace(0.0, 1.0, num_steps + 1)[:, np.newaxis]
                pts = (1.0 - t) * v1 + t * v2

                # Quantize coordinates to box indices
                bx = np.floor(pts[:, 0] / delta).astype(np.int64)
                by = np.floor(pts[:, 1] / delta).astype(np.int64)
                for k in range(len(bx)):
                    occupied_boxes.add((int(bx[k]), int(by[k])))

            counts.append(len(occupied_boxes))

        counts_arr = np.array(counts, dtype=np.float64)
        log_inv_delta = np.log(1.0 / deltas)
        log_n = np.log(counts_arr)

        res = stats.linregress(log_inv_delta, log_n)
        dim = float(res.slope)
        r2 = float(res.rvalue ** 2)

        return dim, r2, deltas, counts_arr
