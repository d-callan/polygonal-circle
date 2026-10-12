"""Tests for triangle waveform generation."""

import numpy as np
from polygonal_circle.waves import triangle_wave


def test_triangle_wave_bounds():
    u = np.linspace(-5.0, 5.0, 1000)
    phi = triangle_wave(u)
    assert np.all(phi >= 0.0)
    assert np.all(phi <= 1.0)


def test_triangle_wave_periodicity():
    u = np.linspace(0.0, 1.0, 100)
    phi1 = triangle_wave(u)
    phi2 = triangle_wave(u + 1.0)
    phi3 = triangle_wave(u + 5.0)
    np.testing.assert_allclose(phi1, phi2, atol=1e-14)
    np.testing.assert_allclose(phi1, phi3, atol=1e-14)


def test_triangle_wave_key_values():
    u = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    phi = triangle_wave(u)
    expected = np.array([0.0, 0.5, 1.0, 0.5, 0.0])
    np.testing.assert_allclose(phi, expected, atol=1e-14)


def test_triangle_wave_linearity():
    # Segment between 0 and 0.5 should have slope +2
    u_rise = np.linspace(0.05, 0.45, 20)
    phi_rise = triangle_wave(u_rise)
    slopes_rise = np.diff(phi_rise) / np.diff(u_rise)
    np.testing.assert_allclose(slopes_rise, 2.0, atol=1e-12)

    # Segment between 0.5 and 1.0 should have slope -2
    u_fall = np.linspace(0.55, 0.95, 20)
    phi_fall = triangle_wave(u_fall)
    slopes_fall = np.diff(phi_fall) / np.diff(u_fall)
    np.testing.assert_allclose(slopes_fall, -2.0, atol=1e-12)
