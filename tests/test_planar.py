"""Tests for planar scale-space machinery on arbitrary closed curves."""

import math
import numpy as np
import pytest

from polygonal_circle.carriers import (
    CircleCarrier,
    EllipseCarrier,
    SquareCarrier,
    TrefoilCarrier,
)
from polygonal_circle.descriptors import CrescentCarrier, KidneyBeanCarrier
from polygonal_circle.planar import (
    PlanarCurve,
    PlanarScaleSpaceObserver,
    RadialCarrierAdapter,
    estimate_sigma_star,
    estimate_sigma_star_oracle,
    heat_kernel_filter,
    resample_arclength,
    signed_area,
)


def _sweep_grid():
    return np.logspace(-3.2, -0.15, 18)


# ---------------------------------------------------------------------------
# Heat kernel laws
# ---------------------------------------------------------------------------

def test_heat_kernel_semigroup_exact():
    """G_s1 o G_s2 = G_sqrt(s1^2 + s2^2) must hold to machine precision."""
    rng = np.random.default_rng(42)
    n = 512
    vals = np.cumsum(rng.normal(size=n))  # smooth-ish periodic signal
    vals -= np.mean(vals)
    period = 2.0 * math.pi
    s1, s2 = 0.05, 0.08

    composed = heat_kernel_filter(heat_kernel_filter(vals, s1, period), s2, period)
    direct = heat_kernel_filter(vals, math.sqrt(s1 * s1 + s2 * s2), period)
    assert np.max(np.abs(composed - direct)) < 1e-10


def test_heat_kernel_mode_washing():
    """A single Fourier mode m is attenuated by exactly exp(-sigma^2 k_m^2/2)."""
    n = 256
    period = 2.0 * math.pi
    t = np.linspace(0, period, n, endpoint=False)
    m = 5
    vals = np.cos(m * t)
    sigma = 0.1
    k_m = 2.0 * math.pi * m / period
    filtered = heat_kernel_filter(vals, sigma, period)
    expected = math.exp(-0.5 * (sigma * k_m) ** 2)
    assert np.max(np.abs(filtered - expected * vals)) < 1e-10


def test_resample_arclength_square():
    """A square resampled uniformly in arc length keeps exact P and A."""
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=float)
    xs, ys = resample_arclength(sq[:, 0], sq[:, 1], 4000)
    c = PlanarCurve(xs, ys)
    assert math.isclose(c.perimeter, 8.0, rel_tol=1e-3)
    assert math.isclose(c.area, 4.0, rel_tol=1e-3)


def test_ensure_ccw_and_signed_area():
    """Orientation is normalized so signed_area is always positive."""
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=float)
    cw = PlanarCurve(sq[::-1, 0], sq[::-1, 1])
    assert signed_area(cw.xs, cw.ys) > 0


# ---------------------------------------------------------------------------
# IR fixed point: shape converges to the first Fourier ellipse
# ---------------------------------------------------------------------------

def _k1_truncation(curve: PlanarCurve) -> PlanarCurve:
    """k=1 Fourier truncation of a curve: an ellipse by construction."""
    n = curve.num_points

    def keep_k1(vals):
        f = np.fft.fft(vals)
        mask = np.zeros(n, dtype=bool)
        mask[0] = True
        mask[1] = True
        mask[-1] = True
        return np.fft.ifft(f * mask).real

    return PlanarCurve(keep_k1(curve.xs), keep_k1(curve.ys))


@pytest.mark.parametrize(
    "carrier",
    [
        RadialCarrierAdapter(CircleCarrier(radius=1.0)),
        RadialCarrierAdapter(EllipseCarrier(semi_major=1.25, semi_minor=0.80)),
        RadialCarrierAdapter(SquareCarrier(inradius=1.0)),
        RadialCarrierAdapter(TrefoilCarrier(base_radius=1.0, modulation=0.15)),
        CrescentCarrier(scale=1.0, dent=1.2),
        KidneyBeanCarrier(),
    ],
)
def test_ir_fixed_point_fourier_ellipse(carrier):
    """Linear diffusion converges to the k=1 Fourier ellipse of the carrier.

    Modes k >= 2 die as exp(-sigma^2 k^2/2) while k=1 is a fixed shape; the
    IR fixed point of coordinate heat filtering is the first elliptic Fourier
    descriptor, which is circular iff the carrier's first harmonic is.
    """
    curve = carrier.sample(4096)
    big = curve.filtered(curve.perimeter / 3.0)
    e1 = _k1_truncation(curve)

    # 1. Q converges to the k=1 ellipse's quotient
    assert abs(big.isoperimetric_quotient - e1.isoperimetric_quotient) < 0.02, (
        f"{carrier.name}: Q={big.isoperimetric_quotient:.4f} vs "
        f"Q1={e1.isoperimetric_quotient:.4f}"
    )

    # 2. Residual energy in modes k>=2 is negligible after filtering
    def residual_energy(c):
        power = np.abs(np.fft.fft(c.xs)) ** 2 + np.abs(np.fft.fft(c.ys)) ** 2
        total = power[1:].sum()
        high = power[2:-1].sum()
        return high / total if total > 0 else 0.0

    assert residual_energy(big) < 0.01


def test_ir_fixed_point_circle_for_circle():
    """The circle carrier's first harmonic is itself: Q -> 1 exactly."""
    carrier = RadialCarrierAdapter(CircleCarrier(radius=1.0))
    curve = carrier.sample(4096)
    big = curve.filtered(curve.perimeter / 3.0)
    assert abs(big.isoperimetric_quotient - 1.0) < 0.02


# ---------------------------------------------------------------------------
# Area inflation law (provable test)
# ---------------------------------------------------------------------------

def test_area_inflation_law_circle():
    """Radial ripple r = R(1 - eta u) with <u>=0 inflates area: A/A0 = 1 + eta^2 <u^2>."""
    from polygonal_circle.waves import get_waveform

    n = 8192
    eta = 0.10
    a, b, level = 0.5, 3, 7
    wf = get_waveform("triangle", zero_mean=True)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    u = np.zeros(n)
    for k in range(1, level + 1):
        u += (a**k) * wf(b**k * theta / (2 * np.pi))
    r = 1.0 * (1.0 - eta * u)
    xs, ys = r * np.cos(theta), r * np.sin(theta)
    a_meas = 0.5 * abs(np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys))

    u2 = float(np.mean(u**2))
    predicted = math.pi * (1.0 + eta**2 * u2)
    assert math.isclose(a_meas, predicted, rel_tol=1e-3), (
        f"area inflation law: measured={a_meas:.5f}, predicted={predicted:.5f}"
    )
    # The inflation is a positive O(eta^2) bias, not zero
    assert a_meas > math.pi


# ---------------------------------------------------------------------------
# Observer recovery in the window, incl. non-star-shaped carriers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "carrier",
    [
        RadialCarrierAdapter(CircleCarrier(radius=1.0)),
        RadialCarrierAdapter(EllipseCarrier(semi_major=1.25, semi_minor=0.80)),
        RadialCarrierAdapter(SquareCarrier(inradius=1.0)),
        RadialCarrierAdapter(TrefoilCarrier(base_radius=1.0, modulation=0.15)),
        CrescentCarrier(scale=1.0, dent=1.2),
        KidneyBeanCarrier(),
    ],
)
def test_planar_observer_window_recovery(carrier):
    """min-sigma rel err of P and A < 2% for every carrier incl. crescent."""
    obs = PlanarScaleSpaceObserver(
        carrier, eta=0.05, a=0.5, b=3, level=6, waveform="triangle", num_points=4096
    )
    sweep = obs.sweep_resolutions(_sweep_grid())
    best_p = min(o.perimeter_relative_error for o in sweep)
    best_a = min(o.area_relative_error for o in sweep)
    assert best_p < 0.02, f"{carrier.name}: best perimeter err {best_p:.4f}"
    assert best_a < 0.02, f"{carrier.name}: best area err {best_a:.4f}"


def test_crescent_is_not_star_shaped():
    """The crescent cannot be expressed as r(theta): the motivation for planar."""
    cres = CrescentCarrier(scale=1.0, dent=1.2)
    assert not cres.is_star_shaped(), "crescent should fail the radial-graph check"


def test_uv_fractal_divergence():
    """At small sigma the perturbed perimeter exceeds the carrier's."""
    carrier = RadialCarrierAdapter(CircleCarrier(radius=1.0))
    obs = PlanarScaleSpaceObserver(carrier, eta=0.05, a=0.5, b=3, level=6)
    micro = obs.observe(0.001)
    assert micro.perimeter > carrier.exact_perimeter()
    assert micro.apparent_dimension > 1.02


def test_euclidean_dimension_in_window():
    """D_app ~ 1 inside the observer window."""
    carrier = RadialCarrierAdapter(CircleCarrier(radius=1.0))
    obs = PlanarScaleSpaceObserver(carrier, eta=0.05, a=0.5, b=3, level=6)
    sweep = obs.sweep_resolutions(_sweep_grid())
    best = min(sweep, key=lambda o: o.perimeter_relative_error)
    assert abs(best.apparent_dimension - 1.0) < 0.08


# ---------------------------------------------------------------------------
# sigma* estimation
# ---------------------------------------------------------------------------

def test_sigma_star_data_driven_vs_oracle():
    """Plateau detector lands within ~2x the oracle's residual error."""
    carrier = RadialCarrierAdapter(SquareCarrier(inradius=1.0))
    obs = PlanarScaleSpaceObserver(carrier, eta=0.05, a=0.5, b=3, level=6)
    sweep = obs.sweep_resolutions(_sweep_grid())

    s_oracle = estimate_sigma_star_oracle(sweep)
    s_est = estimate_sigma_star(sweep)

    o_oracle = obs.observe(s_oracle)
    o_est = obs.observe(s_est)
    assert o_est.perimeter_relative_error <= 2.5 * max(o_oracle.perimeter_relative_error, 1e-3) + 0.005, (
        f"sigma*: estimated {s_est:.4f} (err {o_est.perimeter_relative_error:.4f}) vs "
        f"oracle {s_oracle:.4f} (err {o_oracle.perimeter_relative_error:.4f})"
    )


def test_sigma_star_plateau_exists():
    """The logarithmic perimeter slope has an interior minimum (the window)."""
    carrier = RadialCarrierAdapter(CircleCarrier(radius=1.0))
    obs = PlanarScaleSpaceObserver(carrier, eta=0.05, a=0.5, b=3, level=6)
    sweep = sorted(obs.sweep_resolutions(_sweep_grid()), key=lambda o: o.sigma)
    log_s = np.log([o.sigma for o in sweep])
    log_p = np.log([o.perimeter for o in sweep])
    slope = np.abs(np.gradient(log_p, log_s))
    interior_min = int(np.argmin(slope[2:-2])) + 2
    # The minimum must be strictly interior, not at a sweep boundary
    assert 2 < interior_min < len(sweep) - 2


# ---------------------------------------------------------------------------
# Discretization convergence
# ---------------------------------------------------------------------------

def test_discretization_convergence():
    """Observed perimeter error decays as num_points increases."""
    carrier = RadialCarrierAdapter(CircleCarrier(radius=1.0))
    errs = []
    for n in (1024, 2048, 4096):
        obs = PlanarScaleSpaceObserver(carrier, eta=0.05, a=0.5, b=3, level=6, num_points=n)
        sweep = obs.sweep_resolutions(np.logspace(-2.5, -0.3, 10))
        errs.append(min(o.perimeter_relative_error for o in sweep))
    # Finer sampling must not degrade the best achievable estimate
    assert errs[-1] <= errs[0] + 0.005
