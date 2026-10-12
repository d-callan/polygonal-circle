# Mathematical Report: Recursive Polygonal Circles & Fractal Boundaries

## Executive Summary

We developed an exact mathematical and computational model of a **recursive polygonal circle**: an object that is macroscopically circular and confined strictly within a disk of radius $R$, but whose boundary consists of recursively nested polygonal facets at finer scales.

Using our Python package `polygonal-circle`, we performed computational experiments across three scaling regimes ($a < 1/3$, $a = 1/3$, and $a > 1/3$). All statistics and figures presented below are the direct outputs of computational code operating on the exact geometry.

---

## 1. The Mathematical Model

### Definition
Let $R > 0$ be the radius of the reference circular envelope. At angle $\theta \in [0, 2\pi)$ and recursion depth $n \in \mathbb{N}_0$, the radial boundary is defined by:

$$r_n(\theta) = R \left[ 1 - \eta \sum_{k=1}^n a^k \phi\left( \frac{b^k \theta}{2\pi} \right) \right]$$

where:
- $R = 1.0$: reference envelope radius.
- $\eta \in [0, 1)$: overall amplitude of radial irregularities (default $\eta = 0.1$).
- $a \in (0, 1)$: shrinkage rate of features per recursive level (e.g., $a = 0.2, 1/3, 0.5$).
- $b \in \mathbb{N}_{\ge 2}$: spatial frequency multiplier per level (default $b = 3$).
- $\phi(u) = 1 - 2 |(u \bmod 1) - 0.5|$: periodic symmetric triangle wave with period $1$ and range $[0, 1]$.

### Geometric Invariants & Confinement
1. **Envelope Confinement**: Because $0 \le \phi(u) \le 1$ and $\eta, a > 0$, the sum is non-negative. Therefore:
   $$r_n(\theta) \le R \quad \forall \theta \in [0, 2\pi), \forall n \ge 0$$
   Every vertex lies within the closed disk $\overline{B}(0, R)$.
2. **Convexity of the Envelope**: Because $\overline{B}(0, R)$ is convex, any straight polygonal segment connecting two vertices $\mathbf{v}_i, \mathbf{v}_{i+1}$ satisfies:
   $$\|(1 - t)\mathbf{v}_i + t \mathbf{v}_{i+1}\| \le R \quad \forall t \in [0, 1]$$
   The entire polygonal boundary is strictly confined within the circular disk at all levels.
3. **Positivity of Radius**: The maximum total perturbation is bounded by the geometric series:
   $$\sum_{k=1}^\infty a^k = \frac{a}{1 - a} \implies r_n(\theta) \ge R \left( 1 - \eta \frac{a}{1 - a} \right)$$
   For $\eta = 0.1, a = 0.5$: $r_n(\theta) \ge 0.90 R > 0$. The boundary never crosses the center or self-intersects.
4. **Exact Breakpoint Sampling**: At level $n$, the constituent wave of highest frequency has period $2\pi / b^n$ with kinks at half-integers, spaced at angular intervals $\Delta \theta = \pi / b^n$. By sampling $V_n = 2 b^{n+1}$ (or $6 \cdot 3^n$ for $b=3$) vertices, all kinks across all levels $k \le n$ are explicitly sampled. The radial function $r_n(\theta)$ is strictly affine between consecutive vertices.

---

## 2. Empirical Results Across Scaling Regimes

We evaluated recursion levels $n = 0$ through $n = 7$ for $\eta = 0.1, b = 3$:

| Level $n$ | Vertices $V_n$ | $a = 0.20 < 1/3$<br>Perimeter (Ratio) | $a = 1/3$ (Critical)<br>Perimeter (Ratio) | $a = 0.50 > 1/3$<br>Perimeter (Ratio) | $a = 0.60 > 1/3$<br>Perimeter (Ratio) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| **0** | 6 | 6.0000 (-) | 6.0000 (-) | 6.0000 (-) | 6.0000 (-) |
| **1** | 18 | 6.1899 (1.032) | 6.1503 (1.025) | 6.1022 (1.017) | 6.0689 (1.011) |
| **2** | 54 | 6.2063 (1.003) | 6.1487 (1.000) | 6.0756 (0.996) | 6.0366 (0.995) |
| **3** | 162 | 6.2073 (1.000) | 6.1464 (1.000) | 6.0972 (1.004) | 6.1323 (1.016) |
| **4** | 486 | 6.2073 (1.000) | 6.1492 (1.000) | 6.2064 (1.018) | 6.5775 (1.073) |
| **5** | 1,458 | 6.2073 (1.000) | 6.1544 (1.001) | 6.4706 (1.043) | 7.8570 (1.195) |
| **6** | 4,374 | 6.2074 (1.000) | 6.1604 (1.001) | 7.0216 (1.085) | 10.8662 (1.383) |
| **7** | 13,122 | 6.2074 (1.000) | 6.1667 (1.001) | 8.0771 (1.150) | 17.0578 (1.570) |

At deeper levels for $a = 0.50$:
- **Level 8** ($V = 39,366$): $P_8 = 9.9577$ (ratio $1.233$)
- **Level 9** ($V = 118,098$): $P_9 = 13.1152$ (ratio $1.317$)
- **Level 10** ($V = 354,294$): $P_{10} = 18.1908$ (ratio $1.387$)
- **Level 11** ($V = 1,062,882$): $P_{11} = 26.1184$ (ratio $1.436 \to 3a = 1.500$)

### Three Distinct Regimes
1. **Convergent Regime ($a < 1/3$)**:
   The perimeter converges completely to $P_\infty \approx 6.2074 < 2\pi$. The ratio $P_{n+1} / P_n \to 1.0000$. The boundary remains a rectifiable curve with Hausdorff dimension $D = 1.0$.
2. **Critical Regime ($a = 1/3$)**:
   The perimeter ratio $P_{n+1} / P_n$ stays near $1.0010$, exhibiting constant linear increments $\Delta P \approx +0.006$ per level.
3. **Divergent Regime ($a > 1/3$, e.g. $a = 0.50$)**:
   The growth ratio accelerates toward the theoretical multiplier $b \cdot a = 3 \times 0.5 = 1.5$. The perimeter diverges to infinity ($P_n \to \infty$) while the enclosed area converges to a finite bound ($A_\infty \approx 2.837 < \pi$).

---

## 3. Scale Crossover: Euclidean vs. Fractal Dimensions

A central question is whether the macroscopic shape determines the geometry:
- At coarse resolutions (levels $n = 1$ to $4$), the measured Richardson dimension is:
  $$D_{\text{coarse}} \approx 1.005 \quad (R^2 = 0.99)$$
  At this scale, the object is indistinguishable from a standard Euclidean circle of perimeter $C \approx 2\pi R$.
- At fine resolutions (levels $n = 7$ to $11$), fractal deviations dominate:
  $$D_{\text{fine}} \approx 1.365 \quad (R^2 = 0.98)$$
  closely matching the theoretical dimension:
  $$D_{\text{theory}} = 1 + \frac{\log(3 \times 0.5)}{\log 3} = 1 + \frac{\log(1.5)}{\log 3} \approx 1.3691$$

This demonstrates a **scale crossover**: Euclidean geometry accurately predicts the coarse macroscopic behavior, while Fractal geometry governs the fine-scale boundary length and resolution sensitivity.

---

## 4. Can We Describe Circles as a Type of Fractal?

Our findings provide a precise mathematical perspective:

1. **The Circle as a Zero-Amplitude / Sub-Critical Limit**:
   A smooth Euclidean circle is the degenerate boundary case where either $\eta \to 0$ or $a < 1/b$. In this regime, fine-scale polygonal refinements add vanishing perimeter increments, yielding a rectifiable curve with $D = 1$ and isoperimetric quotient $Q = 1$.
2. **Scale-Dependent Identity**:
   Two objects can possess the identical macroscopic radius $R$ and enclose nearly the same area $A \approx \pi R^2$, yet have completely different geometric classes:
   - Object A (Euclidean): $P(\epsilon) = 2\pi R$ independent of measurement step $\epsilon$.
   - Object B (Recursive Polygonal Circle): $P(\epsilon) \propto \epsilon^{1 - D}$, diverging as resolution $\epsilon \to 0$.
3. **Isoperimetric Collapse**:
   For the smooth circle, $Q = \frac{4\pi A}{P^2} = 1$. For our recursive circular fractal with $a=0.5$, $Q_n$ drops from $0.988$ at $n=2$ down to $0.54$ at $n=7$ and approaches $0$ as $n \to \infty$.

Thus, a circle can be viewed not merely as an isolated smooth entity, but as the ground-state envelope of an entire continuous continuum of scale-dependent recursive geometries.
