"""Skew-product dynamical systems, Lyapunov spectrum, and strange attractors.

This module formalizes the dynamical systems interpretation:
The recursive polygonal circle is the invariant graph attractor of a 2D skew-product map
on the cylinder M = S^1 x R:
    theta_{n+1} = (b * theta_n) mod (2pi)
    y_{n+1}     = a * y_n + phi(theta_n / (2pi))

Properties:
1. Lyapunov exponents:
   - lambda_theta = ln(b) > 0 (chaotic angular expansion / Bernoulli shift)
   - lambda_y     = ln(a) < 0 (radial contraction)
2. Lyapunov sum (phase volume expansion rate):
   - Lambda = lambda_theta + lambda_y = ln(b * a)
   - Lambda < 0 (a < 1/b): Area-contracting, smooth 1D attractor (Euclidean circle).
   - Lambda = 0 (a = 1/b): Area-preserving critical boundary.
   - Lambda > 0 (a > 1/b): Strange fractal attractor (chaotic stretching dominates).
3. Ledrappier-Young / Kaplan-Yorke dimension:
   - D_KY = 1 + max(0, (lambda_theta + lambda_y) / lambda_theta) = 1 + max(0, ln(ba) / ln(b)).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from polygonal_circle.annulus import zero_mean_triangle_wave


@dataclass(frozen=True)
class LyapunovSpectrum:
    """Lyapunov spectrum and dynamical properties of the skew-product map."""
    a: float
    b: int
    lambda_theta: float
    lambda_y: float
    lyapunov_sum: float
    is_chaotic: bool
    is_area_contracting: bool
    kaplan_yorke_dimension: float
    attractor_type: str


def compute_lyapunov_spectrum(a: float, b: int = 3) -> LyapunovSpectrum:
    """Compute exact analytical Lyapunov exponents and Kaplan-Yorke dimension."""
    if not (0.0 < a < 1.0):
        raise ValueError(f"a must be in (0, 1), got {a}")
    if b < 2:
        raise ValueError(f"b must be >= 2, got {b}")

    lam_theta = math.log(float(b))
    lam_y = math.log(a)
    lam_sum = lam_theta + lam_y

    is_chaotic = lam_theta > 0.0
    is_contracting = lam_sum < 0.0

    if lam_sum <= 0.0:
        d_ky = 1.0
        attractor_type = "Smooth 1D Invariant Manifold (Euclidean Circle)"
    else:
        d_ky = 1.0 + lam_sum / lam_theta
        attractor_type = "Strange Fractal Attractor (Nowhere Differentiable)"

    return LyapunovSpectrum(
        a=a,
        b=b,
        lambda_theta=lam_theta,
        lambda_y=lam_y,
        lyapunov_sum=lam_sum,
        is_chaotic=is_chaotic,
        is_area_contracting=is_contracting,
        kaplan_yorke_dimension=d_ky,
        attractor_type=attractor_type,
    )


class SkewProductCircleMap:
    """2D discrete skew-product cylinder map generating the circular fractal boundary.
    
    Coordinates: (theta, y) where theta in [0, 2pi), y in R.
    Map T:
        theta_{n+1} = (b * theta_n) mod (2pi)
        y_{n+1}     = a * y_n + phi(theta_n / (2pi))
    """

    def __init__(self, a: float = 0.5, b: int = 3):
        if not (0.0 < a < 1.0):
            raise ValueError(f"a must be in (0, 1), got {a}")
        if b < 2:
            raise ValueError(f"b must be >= 2, got {b}")
        self.a = float(a)
        self.b = int(b)

    def forward_step(self, theta: float, y: float) -> tuple[float, float]:
        """Single forward iteration of the skew product map."""
        theta_next = (self.b * theta) % (2.0 * math.pi)
        u = (theta / (2.0 * math.pi)) % 1.0
        y_next = self.a * y + zero_mean_triangle_wave(u)
        return theta_next, y_next

    def iterate_trajectory(
        self,
        theta0: float,
        y0: float,
        n_steps: int = 1000,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute a forward trajectory of length n_steps."""
        thetas = np.zeros(n_steps + 1)
        ys = np.zeros(n_steps + 1)

        thetas[0] = theta0 % (2.0 * math.pi)
        ys[0] = y0

        curr_t = thetas[0]
        curr_y = ys[0]

        for i in range(1, n_steps + 1):
            curr_t, curr_y = self.forward_step(curr_t, curr_y)
            thetas[i] = curr_t
            ys[i] = curr_y

        return thetas, ys

    def invariant_graph_point(self, theta: float, depth: int = 12) -> float:
        """Compute point on the invariant graph attractor via backward pullback.
        
        y*(theta) = sum_{k=1}^depth a^k * phi(b^k * theta / (2pi))
        """
        val = 0.0
        for k in range(1, depth + 1):
            u = (self.b**k * theta / (2.0 * math.pi)) % 1.0
            val += (self.a**k) * zero_mean_triangle_wave(u)
        return val

    def sample_invariant_attractor(
        self,
        n_points: int = 1000,
        depth: int = 12,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Sample points (theta, y*(theta)) lying on the invariant attractor."""
        thetas = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
        ys = np.array([self.invariant_graph_point(t, depth=depth) for t in thetas])
        return thetas, ys
