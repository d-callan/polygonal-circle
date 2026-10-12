"""Parametric Euclidean carrier curves for general multi-scale fractal perturbations."""

from abc import ABC, abstractmethod
import math
import numpy as np


class ParametricCarrier(ABC):
    """Abstract base class for macroscopic Euclidean carrier curves."""

    name: str

    @abstractmethod
    def radius_at(self, theta: np.ndarray) -> np.ndarray:
        """Evaluate carrier radius r_0(theta) for theta in [0, 2pi)."""
        pass

    @abstractmethod
    def exact_perimeter(self) -> float:
        """Theoretical Euclidean perimeter L_0."""
        pass

    @abstractmethod
    def exact_area(self) -> float:
        """Theoretical enclosed Euclidean area A_0."""
        pass

    def isoperimetric_quotient(self) -> float:
        """Classical isoperimetric quotient Q = 4 * pi * A / L^2."""
        l0 = self.exact_perimeter()
        a0 = self.exact_area()
        return float(4.0 * math.pi * a0 / (l0**2)) if l0 > 0 else 0.0


class CircleCarrier(ParametricCarrier):
    """Euclidean circle with constant radius R."""

    name = "circle"

    def __init__(self, radius: float = 1.0):
        if radius <= 0:
            raise ValueError(f"Radius must be positive, got {radius}")
        self.radius = float(radius)

    def radius_at(self, theta: np.ndarray) -> np.ndarray:
        return np.full_like(theta, self.radius, dtype=float)

    def exact_perimeter(self) -> float:
        return float(2.0 * math.pi * self.radius)

    def exact_area(self) -> float:
        return float(math.pi * (self.radius**2))


class EllipseCarrier(ParametricCarrier):
    """Euclidean ellipse with semi-axes a (major) and b (minor)."""

    name = "ellipse"

    def __init__(self, semi_major: float = 1.25, semi_minor: float = 0.80):
        if semi_major <= 0 or semi_minor <= 0:
            raise ValueError(f"Semi-axes must be positive, got a={semi_major}, b={semi_minor}")
        self.a = float(semi_major)
        self.b = float(semi_minor)

    def radius_at(self, theta: np.ndarray) -> np.ndarray:
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        denom = np.sqrt((self.b * cos_t) ** 2 + (self.a * sin_t) ** 2)
        return (self.a * self.b) / denom

    def exact_perimeter(self) -> float:
        """Perimeter computed via Ramanujan's second approximation (accurate to O(h^5))."""
        a, b = self.a, self.b
        h = ((a - b) / (a + b)) ** 2
        return float(math.pi * (a + b) * (1.0 + (3.0 * h) / (10.0 + math.sqrt(4.0 - 3.0 * h))))

    def exact_area(self) -> float:
        return float(math.pi * self.a * self.b)


class SuperellipseCarrier(ParametricCarrier):
    """Lamé superellipse / squircle (|x/a|^p + |y/b|^p = 1)."""

    name = "superellipse"

    def __init__(self, semi_major: float = 1.0, semi_minor: float = 1.0, exponent: float = 4.0):
        if semi_major <= 0 or semi_minor <= 0 or exponent <= 0:
            raise ValueError("Dimensions and exponent must be positive.")
        self.a = float(semi_major)
        self.b = float(semi_minor)
        self.p = float(exponent)

    def radius_at(self, theta: np.ndarray) -> np.ndarray:
        cos_t = np.abs(np.cos(theta))
        sin_t = np.abs(np.sin(theta))
        denom = ((cos_t / self.a) ** self.p + (sin_t / self.b) ** self.p) ** (1.0 / self.p)
        return 1.0 / denom

    def exact_area(self) -> float:
        """Area via Dirichlet-Euler beta/gamma integral."""
        gamma_1_p = math.gamma(1.0 + 1.0 / self.p)
        gamma_2_p = math.gamma(1.0 + 2.0 / self.p)
        return float(4.0 * self.a * self.b * (gamma_1_p**2) / gamma_2_p)

    def exact_perimeter(self, n_points: int = 16384) -> float:
        """High-precision numerical quadrature for superellipse perimeter."""
        thetas = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
        r = self.radius_at(thetas)
        xs = r * np.cos(thetas)
        ys = r * np.sin(thetas)
        dx = np.roll(xs, -1) - xs
        dy = np.roll(ys, -1) - ys
        return float(np.sum(np.sqrt(dx**2 + dy**2)))


class TrefoilCarrier(ParametricCarrier):
    """Harmonic 3-lobed closed curve r_0(theta) = R_0 * (1 + alpha * cos(3 * theta))."""

    name = "trefoil"

    def __init__(self, base_radius: float = 1.0, modulation: float = 0.15):
        if base_radius <= 0:
            raise ValueError("Base radius must be positive.")
        if not (0.0 <= modulation < 0.33):
            raise ValueError("Modulation must be in [0, 0.33) to preserve convexity and avoid self-intersection.")
        self.r0 = float(base_radius)
        self.alpha = float(modulation)

    def radius_at(self, theta: np.ndarray) -> np.ndarray:
        return self.r0 * (1.0 + self.alpha * np.cos(3.0 * theta))

    def exact_area(self) -> float:
        """Area by polar integral: (1/2) * integral r^2 dtheta = pi * R0^2 * (1 + alpha^2 / 2)."""
        return float(math.pi * (self.r0**2) * (1.0 + 0.5 * (self.alpha**2)))

    def exact_perimeter(self, n_points: int = 16384) -> float:
        """High-precision numerical integration for trefoil perimeter."""
        thetas = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
        r = self.radius_at(thetas)
        xs = r * np.cos(thetas)
        ys = r * np.sin(thetas)
        dx = np.roll(xs, -1) - xs
        dy = np.roll(ys, -1) - ys
        return float(np.sum(np.sqrt(dx**2 + dy**2)))
