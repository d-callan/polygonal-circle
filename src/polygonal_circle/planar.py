"""Planar scale-space machinery for arbitrary closed curves.

Generalizes the radial ScaleSpaceObserver (which filters the scalar profile
r(theta) and therefore can only represent star-shaped curves) to fully
arbitrary closed planar curves x(t) = (x(t), y(t)) via per-coordinate
periodic heat-kernel filtering.

Design notes:
- Curves are always stored uniformly sampled in arc length, so convolution in
  the curve parameter approximates geometric (arc-length) diffusion rather
  than parameterization-dependent smoothing (Sapiro-Tannenbaum).
- The heat kernel is applied in Fourier space: mode m with wavenumber
  k_m = 2*pi*m/L is attenuated by exp(-0.5 * sigma^2 * k_m^2), which gives an
  EXACT semigroup G_s1 * G_s2 = G_sqrt(s1^2 + s2^2).
- Sigma is measured in arc-length units (physical), not radians of theta.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Iterable, Sequence
import math
import warnings

import numpy as np

from polygonal_circle.waves import WaveformType, get_waveform


# ---------------------------------------------------------------------------
# Resampling & basic curve container
# ---------------------------------------------------------------------------

def resample_arclength(
    xs: np.ndarray, ys: np.ndarray, num_points: int
) -> tuple[np.ndarray, np.ndarray]:
    """Resample a closed polygonal curve uniformly in arc length.

    Parameters
    ----------
    xs, ys : vertex coordinates of a closed curve (first point need not repeat)
    num_points : number of output samples

    Returns
    -------
    (xs_new, ys_new) sampled at uniform arc-length spacing, same closure sense.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    n = xs.shape[0]
    if n < 3:
        raise ValueError("Closed curve requires at least 3 points.")

    # cumulative chord lengths with periodic closure
    dxs = np.roll(xs, -1) - xs
    dys = np.roll(ys, -1) - ys
    seg = np.sqrt(dxs**2 + dys**2)
    s_nodes = np.concatenate([[0.0], np.cumsum(seg)])  # length n+1, s_nodes[-1]=L
    total = s_nodes[-1]
    if total <= 0:
        raise ValueError("Degenerate curve with zero perimeter.")

    x_ext = np.concatenate([xs, xs[:1]])
    y_ext = np.concatenate([ys, ys[:1]])
    s_targets = np.linspace(0.0, total, num_points, endpoint=False)
    xs_new = np.interp(s_targets, s_nodes, x_ext)
    ys_new = np.interp(s_targets, s_nodes, y_ext)
    return xs_new, ys_new


def signed_area(xs: np.ndarray, ys: np.ndarray) -> float:
    """Signed shoelace area; positive iff the curve is counter-clockwise."""
    return float(0.5 * np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys))


def ensure_ccw(xs: np.ndarray, ys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (xs, ys) re-ordered counter-clockwise if necessary."""
    if signed_area(xs, ys) < 0:
        return xs[::-1].copy(), ys[::-1].copy()
    return xs, ys


class PlanarCurve:
    """A closed planar curve uniformly sampled in arc length on s in [0, L)."""

    def __init__(self, xs: np.ndarray, ys: np.ndarray):
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        if xs.shape != ys.shape or xs.ndim != 1 or xs.shape[0] < 3:
            raise ValueError("PlanarCurve requires equal-length 1-D arrays (>=3 pts).")
        self.xs, self.ys = ensure_ccw(xs, ys)

    @classmethod
    def from_points(
        cls, xs: np.ndarray, ys: np.ndarray, num_points: int = 4096
    ) -> "PlanarCurve":
        """Build a uniform arc-length sampled curve from arbitrary vertices."""
        xs_r, ys_r = resample_arclength(xs, ys, num_points)
        return cls(xs_r, ys_r)

    # -- geometry -----------------------------------------------------------

    @property
    def num_points(self) -> int:
        return self.xs.shape[0]

    @property
    def perimeter(self) -> float:
        dx = np.roll(self.xs, -1) - self.xs
        dy = np.roll(self.ys, -1) - self.ys
        return float(np.sum(np.sqrt(dx**2 + dy**2)))

    @property
    def area(self) -> float:
        return float(abs(signed_area(self.xs, self.ys)))

    @property
    def isoperimetric_quotient(self) -> float:
        p = self.perimeter
        return float(4.0 * math.pi * self.area / (p**2)) if p > 0 else 0.0

    @property
    def arc_length_step(self) -> float:
        """Uniform arc-length spacing ds between consecutive samples."""
        return self.perimeter / self.num_points

    # -- spectral calculus (exact for band-limited, uniform-in-s samples) ----

    def _spectral_derivative(self, values: np.ndarray) -> np.ndarray:
        """d(values)/ds via FFT for uniformly arc-length-sampled values."""
        n = values.shape[0]
        period = self.perimeter
        freqs = np.fft.fftfreq(n, d=period / n)  # cycles per unit s
        k = 2.0 * np.pi * freqs
        return np.fft.ifft(1j * k * np.fft.fft(values)).real

    def tangents(self) -> tuple[np.ndarray, np.ndarray]:
        """Unit tangent vectors (dx/ds, dy/ds)."""
        tx = self._spectral_derivative(self.xs)
        ty = self._spectral_derivative(self.ys)
        norm = np.sqrt(tx**2 + ty**2)
        norm[norm == 0] = 1.0
        return tx / norm, ty / norm

    def outward_normals(self) -> tuple[np.ndarray, np.ndarray]:
        """Outward unit normals for a CCW curve: rotate tangent by -90 deg."""
        tx, ty = self.tangents()
        return ty, -tx

    def curvature(self) -> np.ndarray:
        """Signed curvature kappa = x'y'' - y'x'' (arc-length param)."""
        xp = self._spectral_derivative(self.xs)
        yp = self._spectral_derivative(self.ys)
        xpp = self._spectral_derivative(xp)
        ypp = self._spectral_derivative(yp)
        return xp * ypp - yp * xpp

    def min_curvature_radius(self) -> float:
        """Minimum radius of curvature min(1/|kappa|) along the curve."""
        kappa = self.curvature()
        k_max = float(np.max(np.abs(kappa)))
        return float(1.0 / k_max) if k_max > 0 else math.inf

    def filtered(self, sigma: float) -> "PlanarCurve":
        """Apply the periodic heat kernel to both coordinates."""
        return PlanarCurve(
            heat_kernel_filter(self.xs, sigma, self.perimeter),
            heat_kernel_filter(self.ys, sigma, self.perimeter),
        )


def heat_kernel_filter(values: np.ndarray, sigma: float, period: float) -> np.ndarray:
    """Convolve a periodic 1-D signal with the wrapped heat kernel via FFT.

    Fourier mode m (wavenumber k_m = 2*pi*m/period) is attenuated by
    exp(-0.5 * sigma^2 * k_m^2). This realizes the exact heat semigroup:
    filter(s1) o filter(s2) = filter(sqrt(s1^2 + s2^2)).
    """
    if sigma <= 0:
        return np.asarray(values, dtype=float).copy()
    vals = np.asarray(values, dtype=float)
    n = vals.shape[0]
    freqs = np.fft.fftfreq(n, d=period / n)
    k = 2.0 * np.pi * freqs
    attenuation = np.exp(-0.5 * (sigma * k) ** 2)
    return np.fft.ifft(np.fft.fft(vals) * attenuation).real


# ---------------------------------------------------------------------------
# Carrier protocol
# ---------------------------------------------------------------------------

class PlanarCarrier(ABC):
    """Macroscopic Euclidean carrier given as a parametric closed curve.

    Subclasses implement point_at(u) for u in [0, 1) traversing the curve
    once counter-clockwise (any smooth parameterization is acceptable;
    sampling always re-parameterizes to uniform arc length).
    """

    name: str

    @abstractmethod
    def point_at(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate (x(u), y(u)) for u in [0, 1)."""
        ...

    @abstractmethod
    def exact_perimeter(self) -> float:
        """Theoretical Euclidean perimeter L_0."""
        ...

    @abstractmethod
    def exact_area(self) -> float:
        """Theoretical enclosed Euclidean area A_0."""
        ...

    def isoperimetric_quotient(self) -> float:
        l0 = self.exact_perimeter()
        a0 = self.exact_area()
        return float(4.0 * math.pi * a0 / (l0**2)) if l0 > 0 else 0.0

    def sample(self, num_points: int = 4096) -> PlanarCurve:
        """Sample the carrier and re-parameterize to uniform arc length."""
        u = np.linspace(0.0, 1.0, num_points, endpoint=False)
        xs, ys = self.point_at(u)
        return PlanarCurve.from_points(xs, ys, num_points)

    def numerical_perimeter(self, num_points: int = 16384) -> float:
        """High-precision chord-sum perimeter on the raw parametric samples."""
        u = np.linspace(0.0, 1.0, num_points, endpoint=False)
        xs, ys = self.point_at(u)
        dx = np.roll(xs, -1) - xs
        dy = np.roll(ys, -1) - ys
        return float(np.sum(np.sqrt(dx**2 + dy**2)))

    def numerical_area(self, num_points: int = 16384) -> float:
        """High-precision shoelace area on the raw parametric samples."""
        u = np.linspace(0.0, 1.0, num_points, endpoint=False)
        xs, ys = self.point_at(u)
        return float(abs(signed_area(xs, ys)))


class RadialCarrierAdapter(PlanarCarrier):
    """Adapt any existing star-shaped ParametricCarrier r(theta) to planar form."""

    def __init__(self, radial_carrier):
        self._rc = radial_carrier
        self.name = radial_carrier.name

    def point_at(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        theta = 2.0 * np.pi * np.asarray(u, dtype=float)
        r = self._rc.radius_at(theta)
        return r * np.cos(theta), r * np.sin(theta)

    def exact_perimeter(self) -> float:
        return self._rc.exact_perimeter()

    def exact_area(self) -> float:
        return self._rc.exact_area()


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PlanarObservation:
    """Geometric measurement of an arbitrary closed curve at aperture sigma."""

    sigma: float
    perimeter: float
    area: float
    isoperimetric_quotient: float
    apparent_dimension: float
    target_perimeter: float = 0.0
    perimeter_relative_error: float = 0.0
    target_area: float = 0.0
    area_relative_error: float = 0.0


class PlanarScaleSpaceObserver:
    """Scale-space observer on arbitrary closed planar curves.

    The microscopic boundary is a PlanarCarrier perturbed along its outward
    normal field by the multi-scale zero-mean wave sum
        d(s) = eta * sum_k a^k * phi(b^k * s / L),
    where s is arc length and L the carrier perimeter. The observer's
    apparatus is the per-coordinate heat kernel with aperture sigma in
    arc-length units.
    """

    def __init__(
        self,
        carrier: PlanarCarrier,
        eta: float = 0.10,
        a: float = 0.50,
        b: int = 3,
        level: int = 7,
        waveform: str | WaveformType = "triangle",
        num_points: int = 4096,
    ):
        self.carrier = carrier
        self.b = b
        self.level = level
        self.num_points = num_points

        base = carrier.sample(num_points)
        self.base_length = base.perimeter
        self.curve = self._perturb(base, eta, a, b, level, waveform)

        self.target_perimeter = carrier.exact_perimeter()
        self.target_area = carrier.exact_area()

    def _perturb(
        self,
        base: PlanarCurve,
        eta: float,
        a: float,
        b: int,
        level: int,
        waveform: str | WaveformType,
    ) -> PlanarCurve:
        """Displace the base curve along its outward normal field."""
        wf = get_waveform(waveform, zero_mean=True)
        u = np.arange(base.num_points) / base.num_points  # s / L in [0, 1)
        disp = np.zeros(base.num_points)
        for k in range(1, level + 1):
            disp += (a**k) * wf(b**k * u)
        # Physical amplitude ~ eta * (L / 2 pi): matches the radial convention
        # eta * R, since L = 2 pi R for a circle.
        disp *= eta * base.perimeter / (2.0 * math.pi)

        nx, ny = base.outward_normals()

        # Embeddedness guard (generalization of the annular corridor):
        # a normal displacement d(s) folds the offset curve locally whenever
        # kappa(s) * d(s) >= 1 (the offset crosses its own evolute). Clamping
        # to 90% of the signed radius of curvature guarantees no local fold
        # and, as a side effect, suppresses ripple exactly at sharp carrier
        # corners where kappa diverges.
        kappa = base.curvature()
        rho_signed = np.divide(
            1.0, kappa, out=np.full_like(kappa, np.inf), where=np.abs(kappa) > 1e-12
        )
        margin = 0.9
        upper = np.where(rho_signed > 0, margin * rho_signed, np.inf)
        lower = np.where(rho_signed < 0, margin * rho_signed, -np.inf)
        n_clamped = int(np.sum((disp > upper) | (disp < lower)))
        if n_clamped > 0:
            warnings.warn(
                f"Perturbation clipped at {n_clamped}/{base.num_points} samples "
                "to keep the boundary embedded (amplitude exceeds the local "
                "radius of curvature).",
                RuntimeWarning,
            )
        disp = np.clip(disp, lower, upper)

        xs = base.xs + disp * nx
        ys = base.ys + disp * ny

        # Perturbed samples are no longer uniform in arc length: resample.
        return PlanarCurve.from_points(xs, ys, base.num_points)

    def observe(self, sigma: float, delta_sigma_ratio: float = 1.05) -> PlanarObservation:
        """Measure perimeter, area, Q, and apparent dimension at aperture sigma."""
        c = self.curve.filtered(sigma)
        perimeter = c.perimeter
        area = c.area
        q = c.isoperimetric_quotient

        sigma_hi = sigma * delta_sigma_ratio
        sigma_lo = sigma / delta_sigma_ratio
        p_hi = self.curve.filtered(sigma_hi).perimeter
        p_lo = self.curve.filtered(sigma_lo).perimeter
        if p_hi > 0 and p_lo > 0:
            d_log_p = (math.log(p_hi) - math.log(p_lo)) / (
                math.log(sigma_hi) - math.log(sigma_lo)
            )
            apparent_dim = float(1.0 - d_log_p)
        else:
            apparent_dim = float("nan")

        rel_p = (
            abs(perimeter - self.target_perimeter) / self.target_perimeter
            if self.target_perimeter > 0
            else 0.0
        )
        rel_a = (
            abs(area - self.target_area) / self.target_area
            if self.target_area > 0
            else 0.0
        )
        return PlanarObservation(
            sigma=float(sigma),
            perimeter=perimeter,
            area=area,
            isoperimetric_quotient=q,
            apparent_dimension=apparent_dim,
            target_perimeter=self.target_perimeter,
            perimeter_relative_error=rel_p,
            target_area=self.target_area,
            area_relative_error=rel_a,
        )

    def sweep_resolutions(self, sigmas: Iterable[float]) -> list[PlanarObservation]:
        return [self.observe(s) for s in sigmas]

    def theoretical_crossover_scale(self) -> float:
        """Aperture where the base ripple harmonic (wavelength L/b) washes out."""
        return float(self.base_length / (3.0 * self.b))


# ---------------------------------------------------------------------------
# Optimal observer window estimation
# ---------------------------------------------------------------------------

def estimate_sigma_star_oracle(
    observations: Sequence[PlanarObservation],
) -> float:
    """Oracle sigma*: argmin of perimeter relative error (validation only)."""
    valid = [o for o in observations if o.target_perimeter > 0]
    if not valid:
        raise ValueError("Oracle sigma* requires a known target perimeter.")
    return min(valid, key=lambda o: o.perimeter_relative_error).sigma


def estimate_sigma_star(
    observations: Sequence[PlanarObservation],
) -> float:
    """Data-driven sigma*: the perimeter-plateau aperture.

    The observed perimeter P(sigma) is monotone decreasing under heat
    filtering, but its logarithmic slope |d ln P / d ln sigma| attains a local
    minimum inside the observer window: microscopic fractal content and
    macroscopic corner/carrier content each erode faster on their own sides
    of the plateau. No ground truth is required.
    """
    obs = sorted(observations, key=lambda o: o.sigma)
    if len(obs) < 5:
        raise ValueError("Need at least 5 sweep points to locate the plateau.")
    log_s = np.log(np.array([o.sigma for o in obs]))
    log_p = np.log(np.array([o.perimeter for o in obs]))
    slope = np.gradient(log_p, log_s)
    # Exclude the two endpoints on each side where the gradient estimate is
    # one-sided and the physics is dominated by pure UV / pure collapse.
    interior = np.argsort(np.abs(slope[2:-2])) + 2
    return float(obs[interior[0]].sigma)
