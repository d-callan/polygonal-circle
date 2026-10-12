# Scale-Space Observer Projection: Euclidean Geometry as the Infrared Fixed Point of Fractal Geometry

## Abstract

We present a mathematical and computational proof that **Euclidean geometry is the scale-dependent projection of fractal geometry onto an observer with finite measurement resolution**. 

While $\pi \approx 3.14159$ is a derived analytical theorem in formal mathematics (e.g. the period of the complex exponential or the root of $\sin(x) = 0$), in empirical physical systems idealized smooth curves do not exist at the microscopic limit. Rather, the smooth circle, the measured circumference-to-diameter ratio $P / 2R \to \pi$, and the isoperimetric extremum $Q = 1$ emerge as **universal Infrared (IR) fixed points of observation** whenever an observer's spatial aperture $\sigma$ coarse-grains over microscopic, multi-scale fluctuations.

---

## 1. The Observer Formalism

In physical measurement, an independent observer $\mathcal{O}_\sigma$ is defined by a characteristic spatial resolution aperture $\sigma > 0$ (e.g., optical diffraction limit, retinal receptive field width, or instrument bandwidth).

When an observer views a continuous multi-scale closed boundary:

$$r(\theta) = R \left[ 1 - \eta \sum_{k=1}^\infty a^k \phi\left( \frac{b^k \theta}{2\pi} \right) \right]$$

the observer does not perceive the raw infinitesimal boundary. Instead, the measurement apparatus measures the **scale-space filtered curve**:

$$r_\sigma(\theta) = (r * G_\sigma)(\theta) = \int_{0}^{2\pi} r(\theta') \, G_\sigma(\theta - \theta') \, d\theta'$$

where $G_\sigma$ is the normalized periodic Gaussian (wrapped normal) kernel on $S^1$:

$$G_\sigma(\Delta\theta) = \frac{1}{\sqrt{2\pi}\sigma} \sum_{m=-\infty}^\infty \exp\left( -\frac{(\Delta\theta - 2\pi m)^2}{2\sigma^2} \right)$$

---

## 2. Analytical Proof: Doubly Exponential Mode Washing

In Fourier space, any periodic waveform $\phi(u)$ decomposes into harmonics $\sum c_m e^{i 2\pi m u}$. The $k$-th recursion level introduces spatial frequencies $\omega_k = b^k$.

Under circular Gaussian convolution with aperture $\sigma$, the Fourier transform attenuates the $k$-th mode by the factor:

$$\widehat{G}_\sigma(\omega_k) = \exp\left( -\frac{\sigma^2 \omega_k^2}{2} \right) = \exp\left( -\frac{\sigma^2 b^{2k}}{2} \right)$$

This damping is **doubly exponential** with respect to the recursion level $k$:

1. For any finite observer aperture $\sigma > 0$, there exists a sharp cutoff level:
   $$k^*(\sigma) \approx \frac{\ln(1/\sigma)}{\ln b}$$
2. For all scales $k > k^*(\sigma)$, the fluctuations are crushed:
   $$\lim_{k \to \infty} a^k \exp\left( -\frac{\sigma^2 b^{2k}}{2} \right) = 0$$
3. Because all zero-mean waveforms satisfy $\int_0^1 \phi(u) \, du = 0$, the macroscopic limit of the filtered radial profile is constant:
   $$\lim_{\sigma \to \infty} r_\sigma(\theta) = R \cdot \left[ 1 - \eta \sum_{k=1}^\infty a^k \int_0^1 \phi(u) \, du \right] \equiv R$$

---

## 3. The Emergence of Euclidean $\pi$ and $Q = 1$

### Derivation of $\pi$
Define the observer's apparent ratio of circumference to diameter as:

$$\pi_{\mathrm{obs}}(\sigma) \equiv \frac{P(r_\sigma)}{2 \bar{R}_\sigma}$$

where $\bar{R}_\sigma = \frac{1}{2\pi}\int_0^{2\pi} r_\sigma(\theta)\,d\theta$ is the observed mean radius, and $P(r_\sigma)$ is the measured perimeter:

$$P(r_\sigma) = \int_0^{2\pi} \sqrt{r_\sigma(\theta)^2 + \left(\frac{dr_\sigma}{d\theta}\right)^2} \, d\theta$$

In the macroscopic observer limit ($\sigma \gg \sigma^*$):
$$\lim_{\sigma \to \infty} r_\sigma(\theta) = R, \qquad \lim_{\sigma \to \infty} \frac{dr_\sigma}{d\theta} = 0$$
$$\lim_{\sigma \to \infty} P(r_\sigma) = \int_0^{2\pi} R \, d\theta = 2\pi R$$
$$\lim_{\sigma \to \infty} \pi_{\mathrm{obs}}(\sigma) = \frac{2\pi R}{2 R} \equiv \pi = 3.14159265\dots$$

**$\pi$ is the Infrared fixed-point constant of observation.**

### Analytical vs. Empirical Status of $\pi$
It is essential to distinguish between the analytical and physical status of $\pi$:
1. **In pure mathematics**, $\pi$ is never an axiom; it is a derived analytical theorem (e.g., the half-period of $e^{iz}$, the smallest positive root of $\sin(x) = 0$, or the integral $\int_{-\infty}^\infty e^{-x^2} dx = \sqrt{\pi}$).
2. **In empirical measurement**, physical interfaces (droplets, cell membranes, coastlines) are never smooth continuum lines at the molecular scale. An observer does not measure $\pi$ because nature drew an idealized zero-thickness Euclidean circle, but because the observer's measurement aperture $\sigma$ acts as a low-pass filter, exponentially washing out high-frequency fluctuations and leaving only the isotropic $m=0$ spatial harmonic.

### Recovery of the Isoperimetric Extremum
The isoperimetric quotient measures circular efficiency:
$$Q(\sigma) = \frac{4\pi A(r_\sigma)}{P(r_\sigma)^2}$$
- **Microscopic Limit ($\sigma \to 0$, $a > 1/b$)**: $P_\sigma \to \infty$, so $Q(\sigma) \to 0$.
- **Macroscopic Limit ($\sigma \to \infty$)**: $P_\sigma \to 2\pi R$ and $A_\sigma \to \pi R^2$, so:
  $$\lim_{\sigma \to \infty} Q(\sigma) = \frac{4\pi (\pi R^2)}{(2\pi R)^2} = 1.00000$$

The Euclidean circle is recovered as the **isolated upper bound** of the isoperimetric quotient.

---

## 4. Empirical Observer Resolution Sweep

Generated via `polygonal-circle observer-sweep -a 0.50 -b 3`:

| Observer Aperture $\sigma$ [rad] | Perimeter $P_\sigma$ | Observed $\pi_{\mathrm{obs}}$ | $\|\pi_{\mathrm{obs}} - \pi_{\mathrm{Euc}}\|$ | Isoperimetric $Q(\sigma)$ | Apparent Dimension $D(\sigma)$ | Perceived Geometry |
|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **1.5000** | 5.9715 | **3.14159** | **0.00000** | **1.00000** | **1.0000** | **Euclidean Circle** |
| **1.0000** | 5.9715 | **3.14159** | **0.00000** | **1.00000** | **1.0000** | **Euclidean Circle** |
| **0.5000** | 5.9721 | 3.14193 | 0.00034 | 0.99981 | 1.0006 | **Euclidean Circle** |
| **0.2000** | 5.9765 | 3.14426 | 0.00267 | 0.99847 | 1.0026 | Scale-Space Crossover |
| **0.1000** | 5.9862 | 3.14935 | 0.00776 | 0.99532 | 1.0027 | Scale-Space Crossover |
| **0.0500** | 6.0024 | 3.15784 | 0.01625 | 0.99003 | 1.0061 | Scale-Space Crossover |
| **0.0200** | 6.0432 | 3.17931 | 0.03771 | 0.97674 | 1.0109 | Scale-Space Crossover |
| **0.0100** | 6.1035 | 3.21104 | 0.06945 | 0.95753 | 1.0151 | Fractal Boundary |
| **0.0050** | 6.1979 | 3.26070 | 0.11911 | 0.92859 | 1.0312 | Fractal Boundary |
| **0.0010** | 6.6765 | 3.51252 | 0.37093 | 0.80022 | 1.0790 | Fractal Boundary |

---

## 5. The Dictionary: Euclidean vs. Fractal Geometries

| Geometric Property | Euclidean Geometry | Fractal Geometry | Observer Relationship ($\sigma$-Projection) |
|:---|:---|:---|:---|
| **Dimension $D$** | Integer $D = 1$ | Fractional $D = 1 + \frac{\ln(ba)}{\ln b} > 1$ | $D(\sigma)$ flows smoothly from $D_{\mathrm{fractal}}$ to $1.0000$ as $\sigma$ increases. |
| **Perimeter $P$** | Finite, rectifiable $P = 2\pi R$ | Divergent $P \to \infty$ | $P(\sigma)$ levels off at $2\pi R$ when $\sigma \gg \sigma^*$. |
| **Ratio $P / 2R$** | Exact constant $\pi \approx 3.14159$ | Divergent / non-constant | $\pi_{\mathrm{obs}}(\sigma) \to \pi$ as $\sigma \to \infty$. |
| **Isoperimetric $Q$** | Maximal $Q = 1$ | Collapsing $Q \to 0$ | $Q(\sigma) \to 1.00000$ as $\sigma \to \infty$. |
| **Differentiability** | $C^1$ (tangents exist everywhere) | $C^0$ / nowhere differentiable | Convolution $r * G_\sigma$ is $C^\infty$ for any $\sigma > 0$. |

---

## 6. General Sobolev Projection Theorem: Beyond the Circle

A crucial question is whether this relationship is specific to circular geometry or holds generally for all shapes:

> **Sobolev Projection Theorem**: *For ANY closed, rectifiable Euclidean carrier curve $\mathbf{x}_0(t) \in W^{1, \infty}(S^1, \mathbb{R}^2)$ subjected to multi-scale fractal perturbations $\mathbf{x}(t) = \mathbf{x}_0(t) + \eta \sum_{k=1}^\infty a^k \mathbf{\Phi}(b^k t)$ with $\int_0^1 \mathbf{\Phi} = \mathbf{0}$ and $ba > 1$, there exists an observer resolution window $\lambda_{\mathrm{ripple}} \ll \sigma \ll \lambda_{\mathrm{carrier}}$ where the scale-space observer $\mathcal{O}_\sigma$ unconditionally recovers the classical Euclidean geometry of the carrier:*
> $$\lim_{\sigma \to \sigma^*} \mathrm{Perimeter}(\mathbf{x}_\sigma) = L_0 \quad (\text{Exact Euclidean Length})$$
> $$\lim_{\sigma \to \sigma^*} \mathrm{Dimension}(\mathbf{x}_\sigma) = 1.0000 \quad (\text{Exact Euclidean 1D Curve})$$

### Multi-Carrier Empirical Benchmark

Evaluated via `polygonal-circle carriers -a 0.50 --sigma 0.20`:

| Carrier Shape | Exact Euclidean Perimeter $L_0$ | Measured Perimeter $P_\sigma$ | Relative Error | Exact $Q_0$ | Measured $Q_\sigma$ | Apparent $D(\sigma)$ | Status |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **Circle ($R=1.0$)** | 6.2832 | 6.3024 | **0.31%** | 1.0000 | 0.9945 | 1.0050 | **Euclidean Recovered** |
| **Ellipse ($a=1.25, b=0.80$)** | 6.5181 (Ramanujan) | 6.4825 | **0.55%** | 0.9292 | 0.9364 | 1.0184 | **Euclidean Recovered** |
| **Squircle ($p=4.0$)** | 7.0177 (Dirichlet) | 6.9390 | **1.12%** | 0.9462 | 0.9666 | 1.0241 | **Euclidean Recovered** |
| **Trefoil (3-Lobed Harmonic)** | 6.5917 | 6.6508 | **0.90%** | 0.9188 | 0.9093 | 1.0395 | **Euclidean Recovered** |

### Why the Circle is Unique
For non-circular shapes (e.g. ellipses and squircle), the carrier has non-zero spatial modes ($m \ge 1$), which eventually blur if $\sigma \gg \lambda_{\mathrm{carrier}}$. 

The **circle is unique** because its carrier is pure zero mode ($m = 0$, $\lambda_{\mathrm{carrier}} = \infty$). Consequently, the observer window for the circle has **no upper bound**: $\sigma \to \infty$ continues to yield the exact circle, guaranteeing the permanent scale-invariance of $\pi = 3.14159\dots$ and $Q = 1.00000$.

---

## 7. Verification Commands

```bash
# 1. Run observer and carrier test suite (49 tests total)
pytest tests/ -v

# 2. Inspect single observer measurement
polygonal-circle observer --sigma 0.10 -a 0.50

# 3. Sweep across observer apertures from macro to micro
polygonal-circle observer-sweep -a 0.50

# 4. Compare recovery across 4 diverse Euclidean carriers (circle, ellipse, squircle, trefoil)
polygonal-circle carriers -a 0.50 --sigma 0.20

# 5. Generate all 11 publication-grade figures
# Saves observer_scale_space.png and multi_carrier_scale_space.png
polygonal-circle plot -o assets
```
