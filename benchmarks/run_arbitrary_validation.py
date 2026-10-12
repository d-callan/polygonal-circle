"""Empirical validation harness for arbitrary-shape scale-space observation.

Sweeps every PlanarCarrier x waveform combination through the planar observer
and records best-case perimeter/area recovery plus sigma* estimates. Every
number quoted in the documentation is produced by this script.

Usage:
    python benchmarks/run_arbitrary_validation.py
    # -> benchmarks/arbitrary_contour_validation.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from polygonal_circle.carriers import (  # noqa: E402
    CircleCarrier,
    EllipseCarrier,
    SquareCarrier,
    TrefoilCarrier,
)
from polygonal_circle.descriptors import (  # noqa: E402
    CrescentCarrier,
    KidneyBeanCarrier,
)
from polygonal_circle.planar import (  # noqa: E402
    PlanarScaleSpaceObserver,
    RadialCarrierAdapter,
    estimate_sigma_star,
    estimate_sigma_star_oracle,
)

ETA = 0.05
A_SHRINK = 0.50
B_FREQ = 3
LEVEL = 6
NUM_POINTS = 4096
SIGMAS = np.logspace(-3.0, -0.15, 20)
WAVEFORMS = ["triangle", "sine", "parabolic", "multi_harmonic", "smoothed_square"]

CARRIERS = [
    RadialCarrierAdapter(CircleCarrier(radius=1.0)),
    RadialCarrierAdapter(EllipseCarrier(semi_major=1.25, semi_minor=0.80)),
    RadialCarrierAdapter(SquareCarrier(inradius=1.0)),
    RadialCarrierAdapter(TrefoilCarrier(base_radius=1.0, modulation=0.15)),
    CrescentCarrier(scale=1.0, dent=1.2),
    KidneyBeanCarrier(),
]


def main() -> None:
    rows = []
    for carrier in CARRIERS:
        for wf in WAVEFORMS:
            obs = PlanarScaleSpaceObserver(
                carrier, eta=ETA, a=A_SHRINK, b=B_FREQ, level=LEVEL,
                waveform=wf, num_points=NUM_POINTS,
            )
            sweep = obs.sweep_resolutions(SIGMAS)
            best_p = min(o.perimeter_relative_error for o in sweep)
            best_a = min(o.area_relative_error for o in sweep)
            s_star = estimate_sigma_star(sweep)
            s_oracle = estimate_sigma_star_oracle(sweep)
            o_star = obs.observe(s_star)
            rows.append(
                {
                    "carrier": carrier.name,
                    "waveform": wf,
                    "target_perimeter": carrier.exact_perimeter(),
                    "target_area": carrier.exact_area(),
                    "best_perimeter_rel_err": best_p,
                    "best_area_rel_err": best_a,
                    "sigma_star_data": s_star,
                    "sigma_star_oracle": s_oracle,
                    "perimeter_rel_err_at_sigma_star": o_star.perimeter_relative_error,
                    "area_rel_err_at_sigma_star": o_star.area_relative_error,
                    "d_app_at_oracle": obs.observe(s_oracle).apparent_dimension,
                }
            )
            print(
                f"{carrier.name:>14s} {wf:>15s} "
                f"P_err={best_p*100:5.2f}% A_err={best_a*100:5.2f}% "
                f"sigma*={s_star:.4f} (oracle {s_oracle:.4f})"
            )

    out = Path(__file__).resolve().parent / "arbitrary_contour_validation.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nWrote {out} ({len(rows)} rows)")

    df = pd.DataFrame(rows)
    print("\nAcceptance summary:")
    print(f"  max best_perimeter_rel_err = {df.best_perimeter_rel_err.max()*100:.2f}%")
    print(f"  max best_area_rel_err      = {df.best_area_rel_err.max()*100:.2f}%")


if __name__ == "__main__":
    main()
