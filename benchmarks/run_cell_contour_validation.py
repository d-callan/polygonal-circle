"""Reproducibility experiment on cell nuclei (BBBC038) at nested resolutions.

Why this is the CLEAN version of the coastline test:
Natural Earth resolutions are independent cartographic generalizations, so
coarse structure simply absent at 110m can never be recovered by filtering.
Here, a binary nucleus mask is downsampled by block-mean + threshold -- a
genuine low-pass + re-threshold of the SAME boundary. All products are
nested discretizations of one object, so filtering CAN in principle
reconcile them: the falsifiable question is whether the heat-kernel
observer finds an aperture where P, Q and D_app agree.

Coarse contours are rescaled x factor back into fine-pixel units, so sigma
is always in units of the finest pixel pitch.

Output: benchmarks/cell_contour_validation.csv,
        assets/cell_contour_scale_space.png

Data: BBBC038 stage1_train masks (Broad Bioimage Benchmark Collection,
      public), data/cells/*//masks/*.png
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
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from polygonal_circle.contours import (  # noqa: E402
    downsample_mask,
    mask_boundary_ring,
    measure_sweep,
)
from polygonal_circle.planar import PlanarCurve  # noqa: E402

MASK_GLOB = "*/masks/*.png"
FACTORS = [1, 2, 4, 8]  # nested pixel-pitch factors
NUM_POINTS = 4096
SIGMAS = np.logspace(-1.0, 1.9, 60)  # 0.1 .. 80 fine-pixel units
N_CELLS = 6
AREA_RANGE = (900, 4000)  # px -- large enough to survive 8x downsampling
ONSET_THRESHOLD_D = 1.02  # nuclei are near-convex; gentler threshold


def cv(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    m = np.mean(arr)
    return float(np.std(arr) / m) if m != 0 else float("nan")


def pick_cells(data_dir: Path, n: int) -> list[Path]:
    """Deterministic selection of n masks with areas in AREA_RANGE."""
    candidates = []
    for path in sorted(data_dir.rglob(MASK_GLOB)):
        m = np.asarray(Image.open(path)) > 0
        area = int(m.sum())
        if AREA_RANGE[0] <= area <= AREA_RANGE[1]:
            candidates.append((path, area))
        if len(candidates) >= n:
            break
    return candidates


def fractal_onset(sweep) -> float:
    above = [o.sigma for o in sweep if o.apparent_dimension >= ONSET_THRESHOLD_D]
    return float(above[0]) if above else float("nan")


def at_sigma(sweep, sigma: float):
    return min(sweep, key=lambda o: abs(o.sigma - sigma))


def main() -> None:
    data_dir = Path(__file__).resolve().parents[1] / "data" / "cells"
    cells = pick_cells(data_dir, N_CELLS)
    print(f"Selected {len(cells)} nuclei: " +
          ", ".join(f"{p.parent.parent.name[:8]}({a}px)" for p, a in cells))

    rows: list[dict] = []
    cell_curves: dict[str, dict[int, object]] = {}
    cell_sweeps: dict[str, dict[int, list]] = {}

    for idx, (path, area) in enumerate(cells):
        name = f"nucleus_{idx}"
        mask_full = np.asarray(Image.open(path)) > 0
        cell_curves[name], cell_sweeps[name] = {}, {}

        for f in FACTORS:
            mask = downsample_mask(mask_full, f)
            xs, ys = mask_boundary_ring(mask)
            # Rescale coarse contour into fine-pixel units.
            curve = PlanarCurve.from_points(xs * f, ys * f, NUM_POINTS)
            sweep = measure_sweep(curve, SIGMAS)
            cell_curves[name][f], cell_sweeps[name][f] = curve, sweep

        a_mean = float(np.mean([s[0].area for s in cell_sweeps[name].values()]))
        sigma_ref = math.sqrt(a_mean) / 10.0

        for f in FACTORS:
            sweep = cell_sweeps[name][f]
            s_on = fractal_onset(sweep)
            obs_ref = at_sigma(sweep, sigma_ref)
            obs_1px = at_sigma(sweep, 1.0)
            rows.append(
                {
                    "cell": name,
                    "factor": f,
                    "sigma_ref_px": sigma_ref,
                    "sigma_onset_px": s_on,
                    "D_app_at_ref": obs_ref.apparent_dimension,
                    "D_app_at_1px": obs_1px.apparent_dimension,
                    "P_px_at_ref": obs_ref.perimeter,
                    "Q_at_ref": obs_ref.isoperimetric_quotient,
                    "A_px_raw": sweep[0].area,
                }
            )
            print(
                f"  {name} x{f}: onset={s_on:6.2f}px "
                f"D(ref={sigma_ref:4.1f}px)={obs_ref.apparent_dimension:.3f} "
                f"P(ref)={obs_ref.perimeter:6.1f}px Q(ref)={obs_ref.isoperimetric_quotient:.4f}"
            )

    df = pd.DataFrame(rows)
    out_csv = Path(__file__).resolve().parent / "cell_contour_validation.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")

    print("\nCross-resolution CV across nested factors (lower = reproducible):")
    print(f"{'cell':>11s} {'CV[D@ref]':>10s} {'CV[area]':>9s} "
          f"{'CV[P@ref]':>10s} {'CV[Q@ref]':>10s}")
    for name in [f"nucleus_{i}" for i in range(len(cells))]:
        sub = df[df.cell == name]
        print(
            f"{name:>11s} {cv(sub.D_app_at_ref.tolist())*100:9.2f}% "
            f"{cv(sub.A_px_raw.tolist())*100:8.2f}% "
            f"{cv(sub.P_px_at_ref.tolist())*100:9.2f}% "
            f"{cv(sub.Q_at_ref.tolist())*100:9.2f}%"
        )

    # ---------------- figure: 2 cells x (contours, P(s), D(s), Q(s)) --------
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    f_colors = {1: "#38a169", 2: "#2b6cb0", 4: "#805ad5", 8: "#dd6b20"}
    show = [f"nucleus_{i}" for i in range(min(2, len(cells)))]

    for row, name in enumerate(show):
        sub = df[df.cell == name]
        s_ref = float(sub.sigma_ref_px.iloc[0])

        ax = axes[row, 0]
        for f in FACTORS:
            c = cell_curves[name][f]
            ax.plot(c.xs, c.ys, color=f_colors[f], linewidth=0.9, alpha=0.85,
                    label=f"x{f}")
        ax.set_aspect("equal")
        ax.set_title(f"(A) {name}: nested pixel pitches")
        ax.legend(fontsize=7)
        ax.set_xlabel("px")
        ax.set_ylabel("px")

        for col, (key, ylab, title) in enumerate(
            [
                ("perimeter", r"$P(\sigma)$ [px]", "(B) Perimeter: nested products reconcile"),
                ("apparent_dimension", r"$D_{app}(\sigma)$", "(C) Dimension sigmoids"),
                ("isoperimetric_quotient", r"$Q(\sigma)$", "(D) Circularity"),
            ],
            start=1,
        ):
            ax = axes[row, col]
            for f in FACTORS:
                sweep = cell_sweeps[name][f]
                ax.plot([o.sigma for o in sweep], [getattr(o, key) for o in sweep],
                        color=f_colors[f], linewidth=1.6, label=f"x{f}")
            ax.axvline(s_ref, color="k", linestyle="--", linewidth=1.0, alpha=0.6)
            ax.set_xscale("log")
            if key == "apparent_dimension":
                ax.set_ylim(0.98, 1.35)
            ax.set_xlabel(r"aperture $\sigma$ [px]")
            ax.set_ylabel(ylab)
            ax.set_title(title)
            ax.legend(fontsize=7)

    fig.suptitle(
        "Cell Nuclei at Nested Pixel Pitches: "
        "Filtering CAN Reconcile What Cartography Cannot",
        y=1.0,
    )
    fig.tight_layout()
    out_png = Path(__file__).resolve().parents[1] / "assets" / "cell_contour_scale_space.png"
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
