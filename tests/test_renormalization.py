"""Tests for renormalization group formulations and phase transition properties."""

import math
import numpy as np
import pytest

from polygonal_circle.renormalization import (
    OperatorRelevance,
    classify_phase,
    critical_a,
    order_parameter,
    order_parameter_derivative,
    rg_flow_trajectory,
    compute_scaling_collapse,
)


def test_critical_value():
    assert np.isclose(critical_a(3), 1.0 / 3.0)
    assert np.isclose(critical_a(4), 0.25)
    with pytest.raises(ValueError):
        critical_a(1)


def test_order_parameter_phases():
    b = 3
    a_c = 1.0 / 3.0

    # Sub-critical / Euclidean phase: Psi == 0
    assert order_parameter(0.1, b) == 0.0
    assert order_parameter(0.2, b) == 0.0
    assert order_parameter(a_c, b) == 0.0

    # Super-critical / Fractal phase: Psi > 0
    psi_half = order_parameter(0.5, b)
    expected_psi = math.log(1.5) / math.log(3.0)
    assert np.isclose(psi_half, expected_psi)
    assert psi_half > 0.0


def test_critical_exponent_beta():
    """Verify that Psi(a_c + eps) ~ eps^beta with beta = 1."""
    b = 3
    a_c = 1.0 / 3.0

    epsilons = np.logspace(-7, -4, 20)
    psi_vals = np.array([order_parameter(a_c + eps, b) for eps in epsilons])

    # Log-log slope ln(Psi) vs ln(eps) should be 1.0
    log_eps = np.log(epsilons)
    log_psi = np.log(psi_vals)
    slope, intercept = np.polyfit(log_eps, log_psi, deg=1)

    assert np.isclose(slope, 1.0, atol=1e-3)

    # Intercept should be ln(b / ln(b)) because Psi(a_c + eps) ~ (b / ln(b)) * eps
    expected_c = b / math.log(b)
    assert np.isclose(np.exp(intercept), expected_c, rtol=1e-3)


def test_derivative_discontinuity():
    """Verify second-order continuous phase transition via finite jump in dPsi/da."""
    b = 3
    a_c = 1.0 / 3.0

    # Left derivative is 0
    d_left = order_parameter_derivative(a_c - 1e-5, b)
    assert d_left == 0.0

    # Right derivative is b / ln(b)
    d_right = order_parameter_derivative(a_c + 1e-5, b)
    expected_d_right = b / math.log(b)
    assert np.isclose(d_right, expected_d_right, rtol=1e-4)

    # Value at a_c represents the right-hand limit
    assert np.isclose(order_parameter_derivative(a_c, b), expected_d_right)


def test_rg_flow_trajectories():
    b = 3
    # Sub-critical: flows to 0 (IR fixed point)
    traj_ir = rg_flow_trajectory(a=0.2, b=b, initial_perturbation=0.1, steps=6)
    assert traj_ir[-1] < traj_ir[0]
    assert np.isclose(traj_ir[-1], 0.1 * (0.6 ** 6))

    # Critical: constant trajectory
    traj_crit = rg_flow_trajectory(a=1.0 / 3.0, b=b, initial_perturbation=0.1, steps=6)
    np.testing.assert_allclose(traj_crit, 0.1, atol=1e-12)

    # Super-critical: flows to infinity (UV fixed point)
    traj_uv = rg_flow_trajectory(a=0.5, b=b, initial_perturbation=0.1, steps=6)
    assert traj_uv[-1] > traj_uv[0]
    assert np.isclose(traj_uv[-1], 0.1 * (1.5 ** 6))


def test_phase_classification():
    p_euclid = classify_phase(0.2, 3)
    assert p_euclid.relevance == OperatorRelevance.IRRELEVANT
    assert p_euclid.order_parameter == 0.0

    p_crit = classify_phase(1.0 / 3.0, 3)
    assert p_crit.relevance == OperatorRelevance.MARGINAL

    p_fractal = classify_phase(0.5, 3)
    assert p_fractal.relevance == OperatorRelevance.RELEVANT
    assert p_fractal.order_parameter > 0.0


def test_scaling_collapse_data():
    data = compute_scaling_collapse(a_values=[0.45, 0.50, 0.55], levels=[3, 4, 5], b=3, eta=0.1)
    assert len(data.scaling_variable_x) == 9
    assert len(data.rescaled_excess_perimeter_y) == 9
    assert np.all(data.rescaled_excess_perimeter_y > 0)
