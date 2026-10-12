"""Fourier descriptor ingestion and non-star-shaped planar carriers.

FourierDescriptorCarrier turns an arbitrary ordered boundary point cloud into
a smooth PlanarCarrier by truncating the discrete Fourier series of the
uniformly arc-length-resampled contour. Signed area is available in closed
form from the descriptor coefficients:

    A = pi * sum_{k>=1} k * (a_k d_k - b_k c_k)

for x(u) = a0 + sum a_k cos(2 pi k u) + b_k sin(2 pi k u),
    y(u) = c0 + sum c_k cos(2 pi k u) + d_k sin(2 pi k u).

CrescentCarrier and KidneyBeanCarrier provide genuinely non-star-shaped
analytic carriers: no single-valued radial function r(theta) about any
interior origin can represent a crescent, so they are impossible inputs for
the radial ScaleSpaceObserver.
"""

from __future__ import annotations

import math

import numpy as np

from polygonal_circle.planar import (
    PlanarCarrier,
    PlanarCurve,
    resample_arclength,
    signed_area,
)


class FourierDescriptorCarrier(PlanarCarrier):
    """Smooth carrier reconstructed from an arbitrary boundary point cloud.

    The truncated Fourier series both denoises the input and yields a carrier
    with analytic geometry: exact_area() uses the closed-form coefficient
    identity, exact_perimeter() uses high-resolution quadrature.
    """

    def __init__(
        self,
        a_cos: np.ndarray,
        b_sin: np.ndarray,
        c_cos: np.ndarray,
        d_sin: np.ndarray,
        a0: float,
        c0: float,
        name: str = "fourier_descriptor",
    ):
        self.a_cos = np.asarray(a_cos, dtype=float)
        self.b_sin = np.asarray(b_sin, dtype=float)
        self.c_cos = np.asarray(c_cos, dtype=float)
        self.d_sin = np.asarray(d_sin, dtype=float)
        self.a0 = float(a0)
        self.c0 = float(c0)
        self.name = name

    @classmethod
    def from_points(
        cls,
        xs: np.ndarray,
        ys: np.ndarray,
        num_modes: int = 32,
        num_samples: int = 4096,
        name: str = "fourier_descriptor",
    ) -> "FourierDescriptorCarrier":
        """Fit a K-mode Fourier descriptor to an ordered closed point cloud.

        Parameters
        ----------
        xs, ys : ordered boundary vertices (curve need not be star-shaped)
        num_modes : number of Fourier harmonics K retained in the descriptor
        num_samples : arc-length-uniform resample size used for the DFT
        """
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        if xs.shape != ys.shape or xs.shape[0] < 3:
            raise ValueError("Point cloud must contain >= 3 vertices.")
        if signed_area(xs, ys) < 0:
            xs, ys = xs[::-1].copy(), ys[::-1].copy()

        xs_r, ys_r = resample_arclength(xs, ys, num_samples)
        n = num_samples
        k_max = min(num_modes, n // 2 - 1)

        fx = np.fft.fft(xs_r) / n
        fy = np.fft.fft(ys_r) / n

        a_cos = 2.0 * np.real(fx[1 : k_max + 1])
        b_sin = -2.0 * np.imag(fx[1 : k_max + 1])
        c_cos = 2.0 * np.real(fy[1 : k_max + 1])
        d_sin = -2.0 * np.imag(fy[1 : k_max + 1])

        return cls(a_cos, b_sin, c_cos, d_sin, np.real(fx[0]), np.real(fy[0]), name=name)

    @property
    def num_modes(self) -> int:
        return self.a_cos.shape[0]

    def point_at(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        u = np.asarray(u, dtype=float)
        k = np.arange(1, self.num_modes + 1)
        ang = 2.0 * np.pi * np.outer(u, k)
        # x = a0 + sum a_k cos(2 pi k u) + b_k sin(2 pi k u); same for y
        xs = self.a0 + np.cos(ang) @ self.a_cos + np.sin(ang) @ self.b_sin
        ys = self.c0 + np.cos(ang) @ self.c_cos + np.sin(ang) @ self.d_sin
        return xs, ys

    def exact_area(self) -> float:
        """Closed-form signed area from descriptor coefficients."""
        k = np.arange(1, self.num_modes + 1)
        return float(abs(math.pi * np.sum(k * (self.a_cos * self.d_sin - self.b_sin * self.c_cos))))

    def exact_perimeter(self, num_points: int = 16384) -> float:
        return self.numerical_perimeter(num_points)


class CrescentCarrier(PlanarCarrier):
    """Analytic crescent moon: x = cos t + d cos 2t, y = sin t, d > 1.

    For dent > 1 the parametric curve develops a deep concave bite: rays from
    any interior point along the dented direction intersect the boundary
    twice, so no radial function r(theta) can represent it. The curve is
    nevertheless simple (non-self-intersecting) for all dent, since points
    with equal y occur only at mirrored parameters where x differs.

    Signed area is exactly pi for every dent (the k=2 cosine mode in x is
    orthogonal to the k=1 sine mode in y in the Green's-theorem integral).
    """

    def __init__(self, scale: float = 1.0, dent: float = 1.2):
        if scale <= 0:
            raise ValueError("Scale must be positive.")
        if dent < 0:
            raise ValueError("Dent must be non-negative.")
        self.scale = float(scale)
        self.dent = float(dent)
        self.name = "crescent"

    def point_at(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        t = 2.0 * np.pi * np.asarray(u, dtype=float)
        xs = self.scale * (np.cos(t) + self.dent * np.cos(2.0 * t))
        ys = self.scale * np.sin(t)
        return xs, ys

    def exact_area(self) -> float:
        return float(math.pi * self.scale**2)

    def exact_perimeter(self, num_points: int = 16384) -> float:
        return self.numerical_perimeter(num_points)

    def is_star_shaped(self, num_points: int = 4096, num_angles: int = 720) -> bool:
        """Check whether the contour is a radial graph about its centroid."""
        u = np.linspace(0.0, 1.0, num_points, endpoint=False)
        xs, ys = self.point_at(u)
        cx, cy = float(np.mean(xs)), float(np.mean(ys))
        angs = np.arctan2(ys - cy, xs - cx)
        dists = np.hypot(xs - cx, ys - cy)
        # A radial graph must have a unique distance per ray direction.
        order = np.argsort(angs)
        sampled = np.interp(
            np.linspace(-math.pi, math.pi, num_angles, endpoint=False),
            angs[order],
            dists[order],
            period=2.0 * math.pi,
        )
        # Count sign changes of local distance deviations would not detect
        # multi-valuedness directly; instead verify that the map t -> angle is
        # monotone, which is necessary and sufficient for a radial graph.
        d_ang = np.diff(np.unwrap(angs))
        return bool(np.all(d_ang > 0) or np.all(d_ang < 0))


class KidneyBeanCarrier(PlanarCarrier):
    """Asymmetric bean: x = cos t + c1 cos 2t + c2 cos 3t, y = sin t (1 - e cos t).

    Non-convex and (for the default coefficients) non-star-shaped, but
    simple: as with the crescent, equal-y points occur only at mirrored
    parameters where x differs.
    """

    def __init__(self, scale: float = 1.0, c1: float = 0.55, c2: float = 0.18, e: float = 0.45):
        if scale <= 0:
            raise ValueError("Scale must be positive.")
        self.scale = float(scale)
        self.c1 = float(c1)
        self.c2 = float(c2)
        self.e = float(e)
        self.name = "kidney_bean"

    def point_at(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        t = 2.0 * np.pi * np.asarray(u, dtype=float)
        xs = self.scale * (np.cos(t) + self.c1 * np.cos(2.0 * t) + self.c2 * np.cos(3.0 * t))
        ys = self.scale * np.sin(t) * (1.0 - self.e * np.cos(t))
        return xs, ys

    def exact_area(self) -> float:
        """Closed form: expand y = sin t - e sin t cos t = sin t - (e/2) sin 2t."""
        # x harmonics: a1=1, a2=c1, a3=c2; y harmonics: d1=1, d2=-e/2
        a = math.pi * self.scale**2 * (1.0 * 1.0 * 1.0)
        # Only the k=1 term survives: a_2 d_2 - b_2 c_2 = c1*(-e/2) - 0 = -c1 e/2 at k=2
        a += math.pi * self.scale**2 * (2.0 * (self.c1 * (-self.e / 2.0)))
        return float(abs(a))

    def exact_perimeter(self, num_points: int = 16384) -> float:
        return self.numerical_perimeter(num_points)
