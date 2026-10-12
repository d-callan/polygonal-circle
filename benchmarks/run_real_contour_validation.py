"""Cross-resolution reproducibility experiment on real coastlines.

Pre-registered hypotheses (and what the data actually showed):
  H1 (REJECTED): the perimeter-plateau sigma* of the synthetic pipeline is
     reproducible on real contours. Real coastlines have no microscopic ->
     macroscopic crossover; they are self-similar at all scales above the
     digitization granularity, so |d ln P / d ln sigma| has no interior
     plateau minimum and argmin falls to the sweep boundary.
  H2 (HELD): D_app(sigma) IS reproducible -- but only at apertures in the
     shared convergence regime. Each product's D_app curve is a sigmoid
     shifted by its digitization granularity (sigma_onset: Britain 110m ~
     19 km vs 10m ~ 1.7 km). At sigma_ref = sqrt(A)/10 -- a tenth of the
     island's characteristic diameter -- all three products agree to a few
     percent (Madagascar spread 0.2%), and Britain lands on the classic
     Richardson value D ~ 1.25.
  H3 (REJECTED, informative): P and Q become reproducible at large sigma.
     Natural Earth resolutions are independent cartographic generalizations
     (omitted bays are absent, not blurred), so perimeter and circularity
     stay product-dependent at every aperture (16% and 43% spread persists
     at sigma = 100 km for Britain).

Output: benchmarks/real_contour_validation.csv,
        assets/real_contour_scale_space.png

Data: Natural Earth land rings (public domain), data/coastlines/.
Transform: per-island equirectangular tangent projection to km
(src/polygonal_circle/contours.py), identical at all resolutions.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from polygonal_circle.contours import (  # noqa: E402
    island_curve,
    measure_sweep,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "coastlines"
RESOLUTIONS = ["110m", "50m", "10m"]
NUM_POINTS = 8192
SIGMAS = np.logspace(-1.0, 2.48, 70)  # 0.1 .. 300 km, absolute arc-length units
ONSET_THRESHOLD_D = 1.05  # D_app departure from the digitization floor

# Query points (lon, lat) inside each island.
ISLANDS = {
    "britain": (-3.0, 54.0),
    "ireland": (-8.0, 53.4),
    "iceland": (-19.0, 64.9),
    "madagascar": (46.7, -19.5),
    "tasmania": (146.6, -42.0),
    "cuba": (-79.5, 21.9),
}


def cv(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    m = np.mean(arr)
    return float(np.std(arr) / m) if m != 0 else float("nan")


def fractal_onset(sweep) -> float:
    """Smallest sigma where D_app >= ONSET_THRESHOLD_D -- end of the
    digitization floor (below it, polygon chords always measure D ~ 1).
    Returns nan if the product never leaves the floor in the sweep."""
    above = [o.sigma for o in sweep if o.apparent_dimension >= ONSET_THRESHOLD_D]
    return float(above[0]) if above else float("nan")


def at_sigma(sweep, sigma: float):
    return min(sweep, key=lambda o: abs(o.sigma - sigma))


# ---------------- macroscopic limit (IR fixed point) ------------------------

def k1_ellipse(curve):
    """The IR fixed point: reconstruct the contour's k=0,+-1 Fourier modes
    (an ellipse) and return (perimeter, area) of the UN-decayed ellipse.

    This is the only well-defined 'macroscopic' perimeter/area: as
    sigma -> infinity all higher modes vanish and the curve approaches this
    ellipse, but at any finite sigma where fine structure is already gone
    the ellipse itself has partially decayed -- the regimes overlap, so
    macro observables are only defined by extrapolation."""
    fx = np.fft.rfft(curve.xs)
    fy = np.fft.rfft(curve.ys)
    kx = np.zeros_like(fx)
    ky = np.zeros_like(fy)
    kx[:2], ky[:2] = fx[:2], fy[:2]
    from polygonal_circle.planar import PlanarCurve

    e = PlanarCurve(np.fft.irfft(kx, len(curve.xs)),
                    np.fft.irfft(ky, len(curve.ys)))
    return e.perimeter, e.area


def residual_fraction(curve, sigma: float) -> tuple[float, float]:
    """At aperture sigma: fraction of boundary energy in modes k>=2, and the
    k=1 amplitude ratio vs the unfiltered curve (the macro-decay factor)."""
    f0x, f0y = np.fft.rfft(curve.xs), np.fft.rfft(curve.ys)
    c = curve.filtered(sigma)
    fx, fy = np.fft.rfft(c.xs), np.fft.rfft(c.ys)
    resid = float((np.abs(fx[2:]) ** 2 + np.abs(fy[2:]) ** 2).sum())
    tot = float((np.abs(fx[1:]) ** 2 + np.abs(fy[1:]) ** 2).sum())
    k1_decay = float(
        np.hypot(np.abs(fx[1]), np.abs(fy[1]))
        / np.hypot(np.abs(f0x[1]), np.abs(f0y[1]))
    )
    return resid / tot, k1_decay


def macro_window(curve, tol: float = 0.01) -> dict:
    """Smallest sigma where residual (non-k1) energy < tol, plus how much the
    ellipse has itself decayed there -- quantifies the overlap problem."""
    for sigma in np.logspace(0.5, 3.5, 100):
        frac, decay = residual_fraction(curve, sigma)
        if frac < tol:
            return {"sigma_1pct": float(sigma), "k1_decay": decay}
    return {"sigma_1pct": float("nan"), "k1_decay": float("nan")}


def main() -> None:
    rows: list[dict] = []
    curves: dict[str, dict[str, object]] = {}
    sweeps: dict[str, dict[str, list]] = {}

    # First pass: curves, sweeps, raw areas.
    for name, (lon, lat) in ISLANDS.items():
        curves[name], sweeps[name] = {}, {}
        for res in RESOLUTIONS:
            curve = island_curve(DATA_DIR / f"ne_{res}_land.geojson", lon, lat, NUM_POINTS)
            curve = curve.filtered(0.0)  # no-op; keeps code uniform
            sweep = measure_sweep(curve, SIGMAS)
            curves[name][res], sweeps[name][res] = curve, sweep

        # Scale-free reference aperture: a tenth of the characteristic
        # island diameter sqrt(A), using the mean raw area across products.
        a_mean = float(np.mean([s[0].area for s in sweeps[name].values()]))
        sigma_ref = math.sqrt(a_mean) / 10.0

        for res in RESOLUTIONS:
            curve = curves[name][res]
            sweep = sweeps[name][res]
            s_on = fractal_onset(sweep)
            obs_ref = at_sigma(sweep, sigma_ref)
            obs_1km = at_sigma(sweep, 1.0)
            p_e, a_e = k1_ellipse(curve)
            mw = macro_window(curve)
            rows.append(
                {
                    "island": name,
                    "resolution": res,
                    "sigma_ref_km": sigma_ref,
                    "sigma_onset_km": s_on,
                    "D_app_at_ref": obs_ref.apparent_dimension,
                    "D_app_at_1km": obs_1km.apparent_dimension,
                    "P_km_at_ref": obs_ref.perimeter,
                    "Q_at_ref": obs_ref.isoperimetric_quotient,
                    "A_km2_raw": sweep[0].area,
                    "P_km_ellipse": p_e,
                    "A_km2_ellipse": a_e,
                    "sigma_1pct_km": mw["sigma_1pct"],
                    "k1_decay_at_1pct": mw["k1_decay"],
                }
            )
            print(
                f"{name:>11s} {res:>4s} sigma_onset={s_on:7.2f}km "
                f"D(ref={sigma_ref:5.1f}km)={obs_ref.apparent_dimension:.3f} "
                f"P(ref)={obs_ref.perimeter:8.0f}km A={sweep[0].area:10.0f}km^2 "
                f"| P_ellipse={p_e:6.0f} sigma_1%={mw['sigma_1pct']:6.0f}km "
                f"k1_decay={mw['k1_decay']:.2f}"
            )

    df = pd.DataFrame(rows)
    out_csv = Path(__file__).resolve().parent / "real_contour_validation.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")

    print("\nCross-resolution CV (lower = reproducible observable):")
    print(f"{'island':>11s} {'CV[D@ref]':>10s} {'CV[area]':>9s} "
          f"{'CV[P@ref]':>10s} {'CV[Q@ref]':>10s} "
          f"{'CV[P_ell]':>10s} {'CV[A_ell]':>10s}")
    for name in ISLANDS:
        sub = df[df.island == name]
        print(
            f"{name:>11s} {cv(sub.D_app_at_ref.tolist())*100:9.2f}% "
            f"{cv(sub.A_km2_raw.tolist())*100:8.2f}% "
            f"{cv(sub.P_km_at_ref.tolist())*100:9.2f}% "
            f"{cv(sub.Q_at_ref.tolist())*100:9.2f}% "
            f"{cv(sub.P_km_ellipse.tolist())*100:9.2f}% "
            f"{cv(sub.A_km2_ellipse.tolist())*100:9.2f}%"
        )

    # ---------------- figure ----------------
    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    res_colors = {"110m": "#dd6b20", "50m": "#805ad5", "10m": "#38a169"}
    show = ["britain", "madagascar", "iceland"]

    for row, name in enumerate(show):
        sub = df[df.island == name]
        s_ref = float(sub.sigma_ref_km.iloc[0])

        ax = axes[row, 0]
        for res in RESOLUTIONS:
            c = curves[name][res]
            ax.plot(c.xs, c.ys, color=res_colors[res], linewidth=0.8,
                    alpha=0.8, label=res)
        ax.set_aspect("equal")
        ax.set_title(f"(A) {name.title()}: three cartographic products")
        ax.legend(fontsize=7)
        ax.set_xlabel("km")
        ax.set_ylabel("km")

        ax = axes[row, 1]
        for res in RESOLUTIONS:
            sweep = sweeps[name][res]
            ax.plot([o.sigma for o in sweep], [o.perimeter for o in sweep],
                    color=res_colors[res], linewidth=1.6, label=res)
            p_e = float(sub[sub.resolution == res].P_km_ellipse.iloc[0])
            ax.axhline(p_e, color=res_colors[res], linestyle="-.",
                       linewidth=0.9, alpha=0.6)
        ax.axvline(s_ref, color="k", linestyle="--", linewidth=1.0, alpha=0.6,
                   label=r"$\sigma_{ref}=\sqrt{A}/10$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"aperture $\sigma$ [km]")
        ax.set_ylabel(r"$P(\sigma)$ [km]")
        ax.set_title("(B) Power law $+$ extrapolated ellipse asymptotes")
        ax.legend(fontsize=7)

        ax = axes[row, 2]
        for res in RESOLUTIONS:
            sweep = sweeps[name][res]
            ax.plot([o.sigma for o in sweep],
                    [o.apparent_dimension for o in sweep],
                    color=res_colors[res], linewidth=1.6, label=res)
            s_on = float(sub[sub.resolution == res].sigma_onset_km.iloc[0])
            if np.isfinite(s_on):
                ax.axvline(s_on, color=res_colors[res], linestyle=":",
                           linewidth=0.9, alpha=0.7)
        ax.axvline(s_ref, color="k", linestyle="--", linewidth=1.0, alpha=0.6)
        ax.set_xscale("log")
        ax.set_ylim(0.98, 1.45)
        ax.set_xlabel(r"aperture $\sigma$ [km]")
        ax.set_ylabel(r"$D_{app}(\sigma)$")
        ax.set_title("(C) Sigmoid shifted by granularity; converges at ref")
        ax.legend(fontsize=7)

    fig.suptitle(
        "Real Coastlines Under the Heat-Kernel Observer: "
        "Dimension Converges Across Products, Length Does Not",
        y=0.995,
    )
    fig.tight_layout()
    out_png = Path(__file__).resolve().parents[1] / "assets" / "real_contour_scale_space.png"
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
