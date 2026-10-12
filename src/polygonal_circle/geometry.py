"""Geometric representations and exact polygon computations."""

from __future__ import annotations

import numpy as np


class Polygon:
    """A 2D closed polygon defined by an ordered array of vertices.
    
    Vertices are ordered cyclically counterclockwise.
    Edges connect vertex i to vertex (i + 1) % V.
    """

    def __init__(self, vertices: np.ndarray | list[list[float]]):
        verts = np.asarray(vertices, dtype=np.float64)
        if verts.ndim != 2 or verts.shape[1] != 2:
            raise ValueError(f"Vertices must have shape (N, 2), got {verts.shape}")
        if len(verts) < 3:
            raise ValueError(f"A polygon must have at least 3 vertices, got {len(verts)}")
        self.vertices = verts

    @property
    def num_vertices(self) -> int:
        return len(self.vertices)

    @property
    def x(self) -> np.ndarray:
        return self.vertices[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.vertices[:, 1]

    def closed_vertices(self) -> np.ndarray:
        """Return vertices with the first vertex appended to the end."""
        return np.vstack([self.vertices, self.vertices[0:1]])

    def edge_vectors(self) -> np.ndarray:
        """Return vectors v_{i+1} - v_i for each edge."""
        next_verts = np.roll(self.vertices, -1, axis=0)
        return next_verts - self.vertices

    def segment_lengths(self) -> np.ndarray:
        """Return Euclidean length of each edge."""
        return np.linalg.norm(self.edge_vectors(), axis=1)

    def perimeter(self) -> float:
        """Exact perimeter of the polygon: sum of segment lengths."""
        return float(np.sum(self.segment_lengths()))

    def area(self) -> float:
        """Exact enclosed area using the Shoelace formula."""
        x = self.x
        y = self.y
        x_next = np.roll(x, -1)
        y_next = np.roll(y, -1)
        signed_area = 0.5 * np.sum(x * y_next - x_next * y)
        return float(abs(signed_area))

    def isoperimetric_quotient(self) -> float:
        """Compute circularity metric 4 * pi * Area / Perimeter^2.
        
        Equals 1 for a perfect Euclidean circle; strictly < 1 for all other shapes.
        For a fractal boundary with divergent perimeter and bounded area, Q -> 0.
        """
        p = self.perimeter()
        if p == 0:
            return 0.0
        return float(4.0 * np.pi * self.area() / (p * p))

    def radii(self) -> np.ndarray:
        """Euclidean distance of each vertex from origin (0, 0)."""
        return np.linalg.norm(self.vertices, axis=1)

    def max_radius(self) -> float:
        return float(np.max(self.radii()))

    def min_radius(self) -> float:
        return float(np.min(self.radii()))

    def mean_radius(self) -> float:
        return float(np.mean(self.radii()))

    def centroid(self) -> np.ndarray:
        """Centroid of polygon using standard moment formulas."""
        x = self.x
        y = self.y
        x_next = np.roll(x, -1)
        y_next = np.roll(y, -1)
        cross = x * y_next - x_next * y
        a6 = 3.0 * np.sum(cross)
        if abs(a6) < 1e-15:
            return np.mean(self.vertices, axis=0)
        cx = np.sum((x + x_next) * cross) / a6
        cy = np.sum((y + y_next) * cross) / a6
        return np.array([cx, cy], dtype=np.float64)

    def bounding_box(self) -> tuple[float, float, float, float]:
        """Return (xmin, xmax, ymin, ymax)."""
        return (
            float(np.min(self.x)),
            float(np.max(self.x)),
            float(np.min(self.y)),
            float(np.max(self.y)),
        )
