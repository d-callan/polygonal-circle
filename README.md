# Polygonal Circle

**First-principles exploration of recursive polygonal circles: circular macroscopic envelopes with scale-dependent fractal boundaries.**

## Overview

Traditional Euclidean geometry treats a circle as a smooth 1D curve with perimeter $P = 2\pi R$ and area $A = \pi R^2$. Traditional fractal geometry typically builds objects like the Koch snowflake by modifying line segments outward/inward, which does not naturally produce a circular macroscopic envelope.

Here, we construct a family of objects from first principles:
- **Circular envelope**: At large scales, the object stays strictly inside a circular disk of radius $R$.
- **Polygonal boundary**: At any finite stage $n$, its boundary consists of straight segments connecting exact breakpoints.
- **Recursive structure**: Each stage introduces finer radial features at higher spatial frequencies.
- **Scale relationship**: We evaluate the exact perimeter $P_n$, area $A_n$, isoperimetric defect, Richardson coastline dimension, and box-counting dimension across levels.

### The Mathematical Construction

The distance from the origin to the boundary at angle $\theta \in [0, 2\pi)$ and recursion depth $n$ is:

$$r_n(\theta) = R \left[ 1 - \eta \sum_{k=1}^n a^k \phi\left( \frac{b^k \theta}{2\pi} \right) \right]$$

where:
- $R$: radius of the reference circular envelope ($R = 1.0$)
- $n$: recursive refinement level ($n = 0, 1, 2, \dots$)
- $\eta$: global irregularity amplitude ($\eta \in [0, 1)$)
- $a$: shrinkage rate of radial deviations per level ($a > 0$)
- $b$: spatial frequency multiplier per level (default $b = 3$)
- $\phi(u)$: repeating triangular wave with range $[0, 1]$ and period $1$.

Because $\phi \ge 0$, $r_n(\theta) \le R$ for all $\theta$, guaranteeing all vertices lie within the disk of radius $R$. By convexity of the disk, the straight polygonal segments connecting the vertices also lie strictly within the circular envelope.

### Scaling Regimes

Rough scaling predicts the boundary contribution at level $n$ scales as $(b a)^n$:
- **$a < 1/b$ (e.g. $a < 1/3$)**: Deviations shrink faster than frequency grows $\implies$ perimeter converges to a finite limit. Euclidean-like rectifiable curve.
- **$a = 1/b$ (e.g. $a = 1/3$)**: Critical threshold $\implies$ perimeter exhibits marginal logarithmic/constant growth increments.
- **$a > 1/b$ (e.g. $a = 1/2$)**: Irregularity frequency dominates $\implies$ perimeter diverges ($P_n \to \infty$) while enclosed area $A_n$ remains strictly bounded $\implies$ non-rectifiable fractal curve with Hausdorff dimension $D > 1$.

## Installation

```bash
conda create -n polygonal-circle python=3.12 -y
conda activate polygonal-circle
pip install -e .
```

## CLI Usage

The package provides the `polygonal-circle` command-line tool:

```bash
# 1. Inspect exact geometric metrics across recursion levels
polygonal-circle inspect -R 1.0 --eta 0.1 -a 0.5 -b 3 -L 6

# 2. Compare behavior across the three scaling regimes (a < 1/3, a = 1/3, a > 1/3)
polygonal-circle regimes -L 7

# 3. Test Richardson coastline scaling and box-counting dimensions against theory
polygonal-circle dimensions -L 9

# 4. Annular polygonal circles confined between two nested circles
polygonal-circle annulus -R 1.0 -W 0.20 -a 0.5 -L 5

# 5. Sweep the corridor width knob ("wiggle room") from wide to razor-thin
polygonal-circle wiggle-sweep -R 1.0 -L 4

# 6. View the exact statistical physics phase diagram and RG eigenvalues
polygonal-circle phase-transition -b 3

# 7. View the skew-product dynamical system and Lyapunov spectrum
polygonal-circle attractor -b 3

# 8. Test universality and waveform invariance across 5 microscopic generators
polygonal-circle universality -b 3

# 9. Measure boundary through the aperture of an independent observer
polygonal-circle observer --sigma 0.10 -a 0.50

# 10. Sweep observer apertures to prove the emergence of Euclidean pi and Q=1
polygonal-circle observer-sweep -a 0.50

# 11. Test Sobolev projection across diverse Euclidean carriers (circle, ellipse, squircle, trefoil)
polygonal-circle carriers -a 0.50 --sigma 0.20

# 12. Test sharp-cornered polygons (square & triangle) and the corner-rounding V-curve
polygonal-circle polygons -a 0.50

# 13. Sweep an arbitrary planar contour — incl. non-star-shaped crescent (impossible for the radial observer)
polygonal-circle contour-sweep --shape crescent

# 14. Measure a shape at the data-driven optimal aperture (no ground truth needed)
polygonal-circle contour --shape kidney --optimize

# 15. Ingest an arbitrary boundary point cloud (.npy Nx2 or .csv x,y) via Fourier descriptors
polygonal-circle contour-sweep --input my_contour.csv --modes 48

# 16. Generate publication-quality figures
polygonal-circle plot -o assets
```

## Visualizations

### 1. Multi-Level Shape Evolution
![Multi-level Polygons](assets/multilevel_shapes.png)

### 2. Scaling Regimes Comparison
![Scaling Regimes](assets/scaling_regimes.png)

### 3. Richardson Coastline Dimension & Scale Crossover
![Richardson Coastline Scaling](assets/richardson_crossover.png)

### 4. Macroscopic Circle vs. Microscopic Fractal Facets
![Macroscopic vs Microscopic](assets/macroscopic_vs_microscopic.png)

### 5. Annular Polygonal Circle (Nested Circular Boundaries)
![Annular Polygonal Circle](assets/annular_circle.png)

### 6. Wiggle Room Knob Sweep (Narrowing Annular Corridor)
![Corridor Width Sweep](assets/corridor_sweep.png)

### 7. Phase Transition, Renormalization Group, & Universal Scaling
![Phase Transition & RG Flow](assets/phase_transition.png)

### 8. Strange Attractor & Lyapunov Dynamics
![Strange Attractor](assets/strange_attractor.png)

### 9. Empirical Proof of Universality (5 Microscopic Waveforms)
![Universality Proof](assets/universality_proof.png)

### 10. Scale-Space Observer Projection (Emergence of Euclidean Geometry)
![Observer Scale-Space](assets/observer_scale_space.png)

### 11. General Sobolev Projection (Multi-Carrier Euclidean Recovery)
![Multi-Carrier Sobolev Projection](assets/multi_carrier_scale_space.png)

### 12. Polygonal Scale-Space Competition (Corner Rounding vs. Fractal Divergence)
![Polygonal Scale-Space](assets/polygonal_scale_space.png)

### 13. Arbitrary Planar Shapes (Non-Star-Shaped Crescent & Kidney Bean)
![Arbitrary Scale-Space](assets/arbitrary_scale_space.png)

## Empirical Findings & Theoretical Proofs

- See [PHASE_TRANSITIONS_AND_RENORMALIZATION.md](PHASE_TRANSITIONS_AND_RENORMALIZATION.md) for the statistical physics dictionary, mathematical proofs, and RG derivations.
- See [STRANGE_ATTRACTORS_AND_DYNAMICS.md](STRANGE_ATTRACTORS_AND_DYNAMICS.md) for the skew-product cylinder map, Lyapunov spectrum, and Kaplan–Yorke dimension.
- See [UNIVERSALITY_AND_WAVEFORM_INVARIANCE.md](UNIVERSALITY_AND_WAVEFORM_INVARIANCE.md) for multi-waveform invariance proofs and scaling collapse.
- See [SCALE_SPACE_AND_OBSERVER_PROJECTION.md](SCALE_SPACE_AND_OBSERVER_PROJECTION.md) for scale-space filtering proofs, multi-carrier recovery, and corner-rounding scaling.
- See [PLAN_ARBITRARY_SHAPES.md](PLAN_ARBITRARY_SHAPES.md) for the arbitrary-shape extension plan, assessment, and evaluation criteria.

1. **Polygonal Geometries & Corner-Rounding Competition**:
   - For polygons with sharp vertices (squares, triangles), scale-space observation balances two competing effects: corner rounding ($\Delta P \approx -2 N r_{\mathrm{in}} (\tan(\pi/N) - \pi/N) \sigma < 0$) and microscopic fractal divergence ($P_\sigma \to \infty$ as $\sigma \to 0$).
   - The resulting error landscape forms a clean V-curve with an optimal observer window at $\sigma \in [0.03, 0.07]$ rad, recovering the exact Euclidean square perimeter to $< 0.3\%$ and triangle to $< 2.0\%$.
2. **General Sobolev Projection Theorem (Beyond Circles)**:
   - For **ANY** closed, rectifiable Euclidean curve (circle, ellipse, squircle, trefoil), scale-space observation filters out multi-scale fractal ripples and unconditionally recovers the classical Euclidean geometry:
     $$\lim_{\sigma \to \sigma^*} P_\sigma = L_0 \quad (\text{Exact Classical Euclidean Perimeter})$$
     $$\lim_{\sigma \to \sigma^*} D(\sigma) = 1.0000 \quad (\text{Exact 1D Euclidean Curve})$$
   - The circle is unique because its carrier is pure zero mode ($m=0$, $\lambda_{\mathrm{carrier}}=\infty$), meaning its observer window has no upper bound ($\sigma \to \infty$ continues to recover the circle and $\pi$).
3. **Euclidean Geometry as the Infrared Observer Fixed Point**:
   - Euclidean geometry is the scale-space projection of fractal geometry onto an observer with finite spatial aperture $\sigma > 0$.
   - Circular Gaussian filtering attenuates microscopic fluctuations by $\exp(-\sigma^2 b^{2k}/2)$ (doubly exponential mode washing).
   - For any macro-observer ($\sigma \ge 0.5$ rad), the observed boundary converges to an exact Euclidean circle with measured ratio $\pi_{\mathrm{obs}}(\sigma) \to \pi = 3.14159\dots$, isoperimetric quotient $Q_{\mathrm{obs}}(\sigma) \to 1.00000$, and apparent dimension $D(\sigma) \to 1.0000$.
4. **Waveform Universality & Invariance**:
   - The macroscopic behavior is strictly independent of the microscopic tooth generator: triangle ($C^0$), sine ($C^\infty$), parabolic ($C^1$), multi-harmonic ($C^\infty$ asymmetric), and smoothed square ($C^\infty$) all share the identical critical threshold $a_c = 1/b$, identical critical exponent $\beta = 1$, identical asymptotic fractal dimension $D = 1 + \ln(ba)/\ln b$, and identical collapse onto the universal master curve.
5. **Deterministic Chaos & The Lyapunov Spectrum**:
   - The angular coordinate follows a chaotic expanding Bernoulli shift with positive Lyapunov exponent $\lambda_\theta = \ln b > 0$ (exponential orbit separation $\delta \theta_n = \delta \theta_0 b^n$).
   - The radial coordinate undergoes uniform contraction with negative Lyapunov exponent $\lambda_y = \ln a < 0$.
   - When radial contraction dominates ($\lambda_\theta + \lambda_y < 0 \iff a < 1/b$), the invariant attractor is a smooth 1D Euclidean circle ($D = 1$).
   - When chaotic angular stretching overpowers radial contraction ($\lambda_\theta + \lambda_y > 0 \iff a > 1/b$), the invariant manifold becomes a strange fractal attractor with Kaplan–Yorke dimension $D_{KY} = 1 + \frac{\lambda_\theta + \lambda_y}{\lambda_\theta} = 1 + \frac{\ln(ba)}{\ln b}$.
6. **Continuous (Second-Order) Phase Transition**: The excess fractal dimension $\Psi(a) = D(a) - 1$ acts as an exact order parameter, continuous everywhere, with a linear critical exponent $\beta = 1$ ($\Psi \propto (a - a_c)^1$) and a finite jump discontinuity in the derivative $\Delta(d\Psi/da) = b / \ln b \approx 2.7307$ at $a_c = 1/3$.
7. **Renormalization Group (RG) Fixed Points**:
   - $a < 1/3$: Eigenvalue $\lambda = 3a < 1$ (irrelevant operator). Trajectories flow to the stable Infrared (IR) fixed point: the smooth Euclidean circle.
   - $a = 1/3$: Marginal operator ($\lambda = 1$). Constant arithmetic perimeter increments $\Delta P \approx 0.006$ per scale.
   - $a > 1/3$: Relevant operator ($\lambda = 3a > 1$). Flows to the Ultraviolet (UV) non-rectifiable fractal boundary.
8. **Envelope Confinement**: At all levels $n$, $r_n(\theta) \le R$ and all polygonal edges remain strictly inside the circular disk.
9. **Perimeter Divergence**: For $a = 0.50$, the perimeter grows from $6.000$ at $n=0$ to $26.118$ at $n=11$, with the growth ratio $P_{n+1}/P_n \to 3a = 1.500$.
10. **Area Boundedness**: The enclosed area converges to $A_\infty \approx 2.837 < \pi R^2 \approx 3.14159$.
11. **Scale Crossover**:
   - Coarse scale fit ($n=1\dots 4$): $D \approx 1.005$ (looks like a standard Euclidean circle).
   - Fine scale fit ($n=7\dots 11$): $D \approx 1.365$ (converges to theoretical dimension $D = 1 + \log(1.5)/\log(3) \approx 1.369$).
12. **Isoperimetric Collapse**: The isoperimetric quotient $Q = 4\pi A / P^2$ collapses towards $0$ as $n \to \infty$, demonstrating the transition from Euclidean circle to fractal boundary.
13. **Arbitrary Planar Shapes (Beyond Star-Shaped)**: Per-coordinate heat-kernel filtering on arc-length-sampled curves $\mathbf{x}(t)$ recovers perimeter and area for shapes no radial function $r(\theta)$ can represent (crescent, kidney bean): $\le 0.72\%$ perimeter and $\le 0.18\%$ area error across 30 carrier × waveform rows (`benchmarks/arbitrary_contour_validation.csv`). The IR fixed point of linear diffusion is the carrier's first elliptic Fourier descriptor; the embeddedness clamp $\kappa d < 1$ generalizes the corridor constraint to arbitrary carriers.

## Testing

```bash
pytest tests/
```

