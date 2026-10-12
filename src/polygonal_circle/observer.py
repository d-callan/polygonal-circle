"""Scale-space observer model demonstrating Euclidean geometry as the Infrared fixed point of fractal geometry."""

from dataclasses import dataclass
from typing import Iterable, Sequence
import math
import numpy as np

from polygonal_circle.carriers import CircleCarrier, ParametricCarrier
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.waves import WaveformType, get_waveform


@dataclass(frozen=True)
class ObserverObservation:
    """Measurement outcome of an independent observer with aperture sigma."""

    sigma: float
    mean_radius: float
    perimeter: float
    area: float
    pi_obs: float
    isoperimetric_quotient: float
    apparent_dimension: float
    target_perimeter: float = 0.0
    perimeter_relative_error: float = 0.0


class ScaleSpaceObserver:
    """Simulates an observer with spatial resolution aperture sigma viewing a fractal curve.

    The observer's measurement apparatus applies a Gaussian spatial filter G_sigma
    on the angular boundary coordinate theta in S^1.
    """

    def __init__(
        self,
        pc: PolygonalCircle | None = None,
        carrier: ParametricCarrier | None = None,
        eta: float = 0.10,
        a: float = 0.50,
        b: int = 3,
        level: int = 7,
        waveform: str | WaveformType = "triangle",
        num_points: int = 4096,
    ):
        self.level = level
        self.num_points = num_points
        self.thetas = np.linspace(0.0, 2.0 * np.pi, num_points, endpoint=False)
        self.d_theta = 2.0 * np.pi / num_points

        if pc is not None:
            self.pc = pc
            self.carrier: ParametricCarrier = CircleCarrier(radius=pc.radius)
            self.b = pc.b
            self.r_raw = pc.radius_at(self.thetas, level=level)
        elif carrier is not None:
            self.pc = None
            self.carrier = carrier
            self.b = b
            wf_func = get_waveform(waveform, zero_mean=True)
            # Evaluate multi-scale fractal perturbations on carrier radius
            r0 = carrier.radius_at(self.thetas)
            perturbation = np.zeros(num_points, dtype=float)
            for k in range(1, level + 1):
                u = (self.thetas * (b**k)) / (2.0 * np.pi)
                perturbation += (a**k) * wf_func(u)
            self.r_raw = r0 * (1.0 - eta * perturbation)
        else:
            raise ValueError("Either pc or carrier must be provided.")

        self.target_perimeter = self.carrier.exact_perimeter()
        self.target_area = self.carrier.exact_area()

    @classmethod
    def from_carrier(
        cls,
        carrier: ParametricCarrier,
        eta: float = 0.10,
        a: float = 0.50,
        b: int = 3,
        level: int = 7,
        waveform: str | WaveformType = "triangle",
        num_points: int = 4096,
    ) -> "ScaleSpaceObserver":
        """Instantiate observer for an arbitrary macroscopic Euclidean carrier curve."""
        return cls(
            carrier=carrier,
            eta=eta,
            a=a,
            b=b,
            level=level,
            waveform=waveform,
            num_points=num_points,
        )

    def filter_radius(self, sigma: float) -> np.ndarray:
        """Apply circular Gaussian convolution of standard deviation sigma to radial profile."""
        if sigma <= 1e-6:
            return self.r_raw.copy()

        # Gaussian kernel on circular domain
        k_half = int(min(self.num_points // 2, int(4.0 * sigma / self.d_theta) + 1))
        offsets = np.arange(-k_half, k_half + 1)
        kernel = np.exp(-0.5 * (offsets * self.d_theta / sigma) ** 2)
        kernel /= np.sum(kernel)

        # 3-period tiling for exact periodic circular convolution
        r_tiled = np.tile(self.r_raw, 3)
        r_filtered = np.convolve(r_tiled, kernel, mode="same")
        return r_filtered[self.num_points : 2 * self.num_points]

    def filter_curve(self, sigma: float) -> tuple[np.ndarray, np.ndarray]:
        """Return (x, y) coordinates of the boundary as measured by observer sigma."""
        r_filtered = self.filter_radius(sigma)
        xs = r_filtered * np.cos(self.thetas)
        ys = r_filtered * np.sin(self.thetas)
        return xs, ys

    def observe(self, sigma: float, delta_sigma_ratio: float = 1.05) -> ObserverObservation:
        """Compute the full geometric observation for resolution aperture sigma."""
        r_filtered = self.filter_radius(sigma)
        xs = r_filtered * np.cos(self.thetas)
        ys = r_filtered * np.sin(self.thetas)

        # Perimeter via segment lengths
        dx = np.roll(xs, -1) - xs
        dy = np.roll(ys, -1) - ys
        perimeter = float(np.sum(np.sqrt(dx**2 + dy**2)))

        # Shoelace area
        area = float(0.5 * np.abs(np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys)))

        mean_radius = float(np.mean(r_filtered))
        pi_obs = float(perimeter / (2.0 * mean_radius)) if mean_radius > 0 else 0.0
        q = float(4.0 * np.pi * area / (perimeter**2)) if perimeter > 0 else 0.0

        # Estimate apparent dimension via local numerical derivative d(ln P)/d(ln sigma)
        sigma_hi = sigma * delta_sigma_ratio
        sigma_lo = sigma / delta_sigma_ratio
        r_hi = self.filter_radius(sigma_hi)
        r_lo = self.filter_radius(sigma_lo)

        xs_hi = r_hi * np.cos(self.thetas)
        ys_hi = r_hi * np.sin(self.thetas)
        p_hi = float(np.sum(np.sqrt((np.roll(xs_hi, -1) - xs_hi) ** 2 + (np.roll(ys_hi, -1) - ys_hi) ** 2)))

        xs_lo = r_lo * np.cos(self.thetas)
        ys_lo = r_lo * np.sin(self.thetas)
        p_lo = float(np.sum(np.sqrt((np.roll(xs_lo, -1) - xs_lo) ** 2 + (np.roll(ys_lo, -1) - ys_lo) ** 2)))

        d_log_p = (math.log(p_hi) - math.log(p_lo)) / (math.log(sigma_hi) - math.log(sigma_lo))
        apparent_dim = float(1.0 - d_log_p)

        rel_err = abs(perimeter - self.target_perimeter) / self.target_perimeter if self.target_perimeter > 0 else 0.0

        return ObserverObservation(
            sigma=float(sigma),
            mean_radius=mean_radius,
            perimeter=perimeter,
            area=area,
            pi_obs=pi_obs,
            isoperimetric_quotient=q,
            apparent_dimension=apparent_dim,
            target_perimeter=self.target_perimeter,
            perimeter_relative_error=rel_err,
        )

    def sweep_resolutions(self, sigmas: Iterable[float]) -> list[ObserverObservation]:
        """Perform observations across a spectrum of observer resolutions."""
        return [self.observe(s) for s in sigmas]

    def theoretical_crossover_scale(self) -> float:
        """Estimate the angular resolution sigma* where microscopic perturbations become visible."""
        # The base frequency is b^1. Microscopic fluctuations become detectable
        # when the observer's aperture is smaller than the first harmonic wavelength 2pi / b.
        return float(2.0 * np.pi / (self.b * 3.0))
