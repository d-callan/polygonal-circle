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
