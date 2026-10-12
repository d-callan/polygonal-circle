"""Renormalization group analysis and phase transitions for recursive polygonal circles.

This module formalizes the exact statistical physics and RG relationships:
1. Critical threshold: a_c = 1 / b.
2. Order parameter: Psi(a) = D - 1 = max(0, ln(ba) / ln(b)).
3. Critical exponent: Psi(a) ~ (a - a_c)^beta with beta = 1.
4. Continuous phase transition: Jump discontinuity in dPsi/da at a_c.
5. RG Flow: Linearized eigenvalue lambda = b * a, classifying perturbations into:
   - Irrelevant (lambda < 1, a < a_c) -> IR stable fixed point (Euclidean circle).
   - Marginal (lambda = 1, a = a_c) -> Logarithmic scaling.
   - Relevant (lambda > 1, a > a_c) -> UV unstable fixed point (Fractal boundary).
6. Universal Scaling Collapse: Master curve for perimeter growth across levels and widths.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
import numpy as np

from polygonal_circle.annulus import AnnularPolygonalCircle
from polygonal_circle.model import PolygonalCircle


class OperatorRelevance(str, Enum):
    """RG classification of the radial perturbation operator."""
    IRRELEVANT = "irrelevant (flows to smooth circle)"
    MARGINAL = "marginal (constant increments)"
    RELEVANT = "relevant (flows to fractal boundary)"


@dataclass(frozen=True)
class PhaseState:
    """Exact phase state and RG properties for parameter a."""
    a: float
    b: int
    a_c: float
    rg_eigenvalue: float
    relevance: OperatorRelevance
    order_parameter: float
    order_parameter_derivative: float
    phase_name: str


def critical_a(b: int = 3) -> float:
    """Exact critical shrinkage parameter a_c = 1 / b."""
    if b < 2:
        raise ValueError(f"Spatial multiplier b must be >= 2, got {b}")
    return 1.0 / float(b)


def order_parameter(a: float, b: int = 3) -> float:
    """Compute the excess fractal dimension order parameter Psi(a) = D - 1."""
    if not (0.0 < a < 1.0):
        raise ValueError(f"a must be in (0, 1), got {a}")
    a_c = critical_a(b)
    if a <= a_c:
        return 0.0
    return math.log(b * a) / math.log(b)


def order_parameter_derivative(a: float, b: int = 3) -> float:
    """Compute the analytical derivative dPsi/da.
    
    For a < a_c: dPsi/da = 0.
    For a > a_c: dPsi/da = 1 / (a * ln(b)).
    At a_c^+: limit is b / ln(b).
    """
    a_c = critical_a(b)
    if a < a_c:
        return 0.0
    if a == a_c:
        # Return right-sided derivative limit at a_c^+
        return float(b) / math.log(b)
    return 1.0 / (a * math.log(b))


def classify_phase(a: float, b: int = 3) -> PhaseState:
    """Compute exact phase classification and RG properties."""
    a_c = critical_a(b)
    lam = float(b) * a
    psi = order_parameter(a, b)
    d_psi = order_parameter_derivative(a, b)

    tol = 1e-12
    if a < a_c - tol:
        relevance = OperatorRelevance.IRRELEVANT
        phase = "Euclidean (Rectifiable Curve)"
    elif abs(a - a_c) <= tol:
        relevance = OperatorRelevance.MARGINAL
        phase = "Critical Boundary"
    else:
        relevance = OperatorRelevance.RELEVANT
        phase = "Fractal (Non-Rectifiable Boundary)"

    return PhaseState(
        a=a,
        b=b,
        a_c=a_c,
        rg_eigenvalue=lam,
        relevance=relevance,
        order_parameter=psi,
        order_parameter_derivative=d_psi,
        phase_name=phase,
    )


def rg_flow_trajectory(
    a: float,
    b: int = 3,
    initial_perturbation: float = 0.1,
    steps: int = 8,
) -> np.ndarray:
    """Compute RG discrete flow of perturbation amplitude across scales.
    
    Delta P_{n+1} = (b * a) * Delta P_n
    """
    lam = float(b) * a
    n = np.arange(steps + 1)
    return initial_perturbation * (lam ** n)


@dataclass(frozen=True)
class ScalingCollapseData:
    """Data structures for universal scaling collapse."""
    scaling_variable_x: np.ndarray
    rescaled_excess_perimeter_y: np.ndarray
    labels: list[str]


def compute_scaling_collapse(
    a_values: list[float] = [0.40, 0.45, 0.50, 0.55, 0.60],
    levels: list[int] = [3, 4, 5, 6, 7],
    b: int = 3,
    eta: float = 0.1,
) -> ScalingCollapseData:
    """Verify universal data collapse for the divergent fractal phase.
    
    The asymptotic excess perimeter scales as:
        (P_n - P_0) ~ b^{n * Psi(a)}
    
    Rescaling Y = (P_n - P_0) * b^{-n * Psi(a)} against X = n (or (a - a_c)*b^n)
    demonstrates that the asymptotic growth rate collapses across all a > a_c.
    """
    x_list = []
    y_list = []
    labels = []

    a_c = critical_a(b)

    for a in a_values:
        if a <= a_c:
            continue
        pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=b)
        psi = order_parameter(a, b)

        p0 = pc.polygon(0).perimeter()
        for lvl in levels:
            p_n = pc.polygon(lvl).perimeter()
            excess = p_n - p0
            # Rescaled quantity
            y_val = excess / (b ** (lvl * psi))
            x_val = lvl

            x_list.append(x_val)
            y_list.append(y_val)
            labels.append(f"a={a:.2f}")

    return ScalingCollapseData(
        scaling_variable_x=np.array(x_list),
        rescaled_excess_perimeter_y=np.array(y_list),
        labels=labels,
    )
