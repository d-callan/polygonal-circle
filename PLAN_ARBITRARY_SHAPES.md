# Plan: Multi-Resolution Area & Perimeter Estimation for Arbitrary Shapes

Extends the scale-space observer framework from **star-shaped radial graphs** $r(\theta)$ to **fully arbitrary closed planar curves** $\mathbf{x}(t) = (x(t), y(t))$, including non-star-shaped contours that no radial function can represent.

## 1. Assessment of the External Proposal

The motivating proposal (external agent) is ~70% accurate. Recorded here for provenance:

**Correct:**

- Per-coordinate Gaussian convolution of $\mathbf{x}(t) = (x(t), y(t))$ is the right generalization — the current observer filters the scalar profile $r(\theta)$, which only represents star-shaped curves.
- Green's theorem area $A = \tfrac{1}{2} \oint (x\,dy - y\,dx)$ is correct (already implemented as the shoelace formula in `observer.py`).
- Fourier descriptors are the right tool for ingesting arbitrary boundary point clouds.

**Stale or wrong:**

- The "circle-centric model" framing is outdated: `carriers.py` already supports arbitrary star-shaped carriers. The true limitation is **star-shapedness** — any contour where a ray from the origin intersects the boundary more than once cannot be written as $r(\theta)$.
- The proposed $\sigma^*$ optimizer is circular: minimizing error against "true underlying metrics" requires ground truth, which is unavailable on real data. It is retained as an *oracle* for validation only; a data-driven criterion is required for deployment.
- It ignores the parameterization subtlety: convolution in a fixed parameter $t$ is not geometric diffusion unless $t \propto$ arc length (Sapiro–Tannenbaum). Arc-length resampling is mandatory.
- It misses two provable laws that become regression tests (§4).

## 2. Mathematical Foundation

### Planar heat-kernel scale space

For a closed curve sampled uniformly in arc length $t \in [0, L)$, apply the periodic heat kernel independently per coordinate:

$$\mathbf{x}_\sigma(t) = \bigl(x * G_\sigma\bigr)(t),\qquad \widehat{G_\sigma}(k) = e^{-\sigma^2 k^2 / 2}$$

Properties exploited as tests:

- **Exact semigroup**: $G_{\sigma_1} * G_{\sigma_2} = G_{\sqrt{\sigma_1^2+\sigma_2^2}}$.
- **IR fixed point**: Fourier modes die as $e^{-\sigma^2 k^2/2}$, so only $k{=}1$ survives at large $\sigma$ — the normalized shape converges to the **first elliptic Fourier descriptor** (any pure $k{=}1$ truncation is an ellipse), with $Q \to Q_1 \le 1$. The circle is the fixed point iff the carrier's first harmonic is circular. This is the linear analog — not the full strength — of Gage–Hamilton–Grayson, whose nonlinear curvature flow does converge to a round point.
- **Mode washing**: a perturbation at harmonic $k$ is attenuated by $e^{-\sigma^2 k^2/2}$ — identical mode-washing law as the radial observer, now coordinatewise.

### Area under zero-mean ripple (provable test law)

For a radial perturbation $r = r_0(1 - \eta u)$ with $\langle u \rangle = 0$:

$$A(0) = \tfrac{1}{2}\int r^2\,d\theta = A_0\left(1 + \eta^2 \langle u^2 \rangle\right)$$

Zero-mean fractal ripples **inflate** apparent area by an $O(\eta^2)$ bias even before filtering — a falsifiable prediction the harness must reproduce.

### Fourier descriptor ingestion

Arbitrary boundary point cloud $\to$ uniform arc-length resample $\to$ DFT $\to$ truncate to $K$ modes $\to$ smooth base carrier $\mathbf{x}_0(t)$. Signed area has a closed form in the descriptor coefficients:

$$A = \pi \sum_{k=1}^{K} k\left(a_k d_k - b_k c_k\right)$$

for $x(t) = \sum a_k\cos kt + b_k\sin kt$, $y(t) = \sum c_k\cos kt + d_k\sin kt$. Perimeter is evaluated by high-resolution quadrature on the reconstructed curve.

### Perturbation model

Default **normal-field** perturbation generalizing the corridor model:

$$\mathbf{x}(t) = \mathbf{x}_0(t) + \eta \sum_{k=1}^{n} a^k \phi\!\left(\frac{b^k t}{2\pi}\right)\,\hat{\mathbf{n}}(t)$$

with $\hat{\mathbf{n}}$ the outward unit normal of the carrier. Amplitude guard: $\eta \sum_k a^k |\phi|_\infty < \min_t \rho(t)$ where $\rho$ is the carrier's radius of curvature, preventing self-intersection. Radial mode remains available for star-shaped carriers (backward compatibility).

### Optimal observer window

- **Oracle** $\sigma^*_{\text{oracle}} = \arg\min_\sigma |P(\sigma) - P_0|/P_0$ — validation only.
- **Data-driven** $\sigma^*$: perimeter-plateau detector $\arg\min_\sigma |dP/d\ln\sigma|$ over the mid-range, cross-checked against the $D_{\text{app}}$ crossover point. Acceptance: data-driven pick yields error within ~2x of oracle minimum.

## 3. Milestones

| # | Deliverable | Files |
|---|-------------|-------|
| M0 | `area_relative_error` + `target_area` in `ObserverObservation`, surfaced in CLI sweeps | `observer.py`, `cli.py` |
| M1 | `PlanarCurve` (arc-length resample, P/A/Q), FFT heat-kernel `filter_planar`, `ParametricCarrier2D` ABC + radial adapter | `planar.py` (new) |
| M2 | `FourierDescriptorCarrier` (closed-form area), `KidneyBeanCarrier`, `CrescentCarrier` | `descriptors.py` (new) |
| M3 | Normal-field perturbation + `PlanarScaleSpaceObserver` | `planar.py` |
| M4 | Oracle + data-driven $\sigma^*$ estimators | `planar.py` |
| M5 | CLI `contour-sweep` / `contour`, `plot_arbitrary_scale_space`, docs | `cli.py`, `plotting.py`, `*.md` |
| M6 | Test suite + `benchmarks/arbitrary_contour_validation.csv` | `tests/test_planar.py`, `tests/test_descriptors.py` |

## 4. Proof & Evaluation Criteria

| Test | Law checked | Tolerance |
|------|-------------|-----------|
| Semigroup | $\|G_{\sigma_1}\circ G_{\sigma_2} - G_{\sqrt{\sigma_1^2+\sigma_2^2}}\|$ | $< 10^{-10}$ |
| IR fixed point | $Q(\sigma) \to Q_1$ (k=1 Fourier ellipse) as $\sigma$ grows; $Q \to 1$ for circle carrier | $\|Q - Q_1\| < 0.02$ |
| Area inflation | $A(0)/A_0 = 1 + \eta^2\langle u^2\rangle$, radial ripple | rel. err $< 1\%$ |
| Window recovery | $\min_\sigma$ rel-err of $P$ **and** $A$, carriers x waveforms | $< 2\%$ |
| Non-star-shaped | crescent/kidney recovery — impossible under $r(\theta)$ | $< 2\%$ |
| Descriptor roundtrip | point cloud $\to$ carrier: P, A | $< 2\%$ |
| $\sigma^*$ estimator | data-driven vs oracle error | $\le 2\times$ oracle min |
| Discretization | error vs $N \in \{1024, 2048, 4096, 8192\}$ | monotone decay |
| Regression | all pre-existing tests | 52 green |

All numbers reported in docs must be outputs of `benchmarks/run_arbitrary_validation.py` (repo convention: empirical claims are script outputs, never hand-written).
