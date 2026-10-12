"""Ingestion and measurement of real-world boundary contours.

Data provenance (see AGENTS.md: record source, format, transformations):
  Source  : Natural Earth land polygons, public domain
            (https://github.com/nvkelso/natural-earth-vector geojson branch)
  Format  : GeoJSON FeatureCollection of Polygon/MultiPolygon in CRS84
            (lon/lat degrees); we take the OUTER ring of each polygon only
            (interior lake holes are ignored consistently at all resolutions).
  Transform: each island ring is projected to a local equirectangular
            tangent plane about its own vertex centroid and expressed in km.
            The identical transform is applied at every digitization
            resolution, so projection distortion is held fixed in the
            cross-resolution reproducibility test.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np

from polygonal_circle.planar import PlanarCurve

EARTH_RADIUS_KM = 6371.0088


def load_land_rings(path: str | Path) -> list[np.ndarray]:
    """Load outer boundary rings (lon/lat, degrees) from a Natural Earth geojson.

    Returns one (N, 2) array per polygon outer ring; interior holes are
    skipped (documented above).
    """
    with open(path) as fh:
        data = json.load(fh)
    rings: list[np.ndarray] = []
    for feature in data["features"]:
        geom = feature["geometry"]
        polys = (
            geom["coordinates"]
            if geom["type"] == "MultiPolygon"
            else [geom["coordinates"]]
        )
        for poly in polys:
            rings.append(np.asarray(poly[0], dtype=float))
    return rings


def select_ring_near(
    rings: Sequence[np.ndarray], lon: float, lat: float
) -> np.ndarray:
    """Pick the ring containing the vertex nearest to a query point (deg)."""
    def min_dist(ring: np.ndarray) -> float:
        return float(np.min(np.hypot(ring[:, 0] - lon, ring[:, 1] - lat)))

    return min(rings, key=min_dist)


def project_km(ring_lonlat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Equirectangular tangent-plane projection about the ring centroid (km)."""
    lon0 = math.radians(float(np.mean(ring_lonlat[:, 0])))
    lat0 = math.radians(float(np.mean(ring_lonlat[:, 1])))
    xs = EARTH_RADIUS_KM * np.cos(lat0) * (np.radians(ring_lonlat[:, 0]) - lon0)
    ys = EARTH_RADIUS_KM * (np.radians(ring_lonlat[:, 1]) - lat0)
    return xs, ys


def island_curve(
    path: str | Path,
    lon: float,
    lat: float,
    num_points: int = 8192,
) -> PlanarCurve:
    """Load an island ring and return a uniform-arc-length PlanarCurve in km."""
    ring = select_ring_near(load_land_rings(path), lon, lat)
    xs, ys = project_km(ring)
    return PlanarCurve.from_points(xs, ys, num_points)


@dataclass(frozen=True)
class ContourMeasurement:
    """Observable record for one aperture on a real (target-free) contour."""

    sigma: float
    perimeter: float
    area: float
    isoperimetric_quotient: float
    apparent_dimension: float


def measure_sweep(
    curve: PlanarCurve,
    sigmas: Iterable[float],
    delta_sigma_ratio: float = 1.05,
) -> list[ContourMeasurement]:
    """Sweep the heat-kernel aperture on a raw contour, no perturbation.

    D_app is computed exactly as in PlanarScaleSpaceObserver.observe: a local
    Richardson slope between sigma / 1.05 and sigma * 1.05, reported as
    1 - d ln P / d ln sigma.
    """
    out: list[ContourMeasurement] = []
    for sigma in sigmas:
        c = curve.filtered(sigma)
        p_hi = curve.filtered(sigma * delta_sigma_ratio).perimeter
        p_lo = curve.filtered(sigma / delta_sigma_ratio).perimeter
        if p_hi > 0 and p_lo > 0:
            d_log_p = (math.log(p_hi) - math.log(p_lo)) / (
                2.0 * math.log(delta_sigma_ratio)
            )
            d_app = float(1.0 - d_log_p)
        else:
            d_app = float("nan")
        out.append(
            ContourMeasurement(
                sigma=float(sigma),
                perimeter=c.perimeter,
                area=c.area,
                isoperimetric_quotient=c.isoperimetric_quotient,
                apparent_dimension=d_app,
            )
        )
    return out


def mask_boundary_ring(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Outer boundary ring of a binary mask, ordered into a closed polygon.

    Pure numpy: boundary = mask minus its 4-neighbour erosion, then boundary
    pixels are sorted by angle about their centroid. This is exact for
    star-shaped sets; nuclei are near-convex in practice, and the shoelace
    area of the ring is checked against the pixel count by callers (a ~5%
    deficit is expected from pixel-corner clipping).
    """
    m = mask.astype(bool)
    up = np.roll(m, 1, 0)
    up[0, :] = False
    dn = np.roll(m, -1, 0)
    dn[-1, :] = False
    lf = np.roll(m, 1, 1)
    lf[:, 0] = False
    rt = np.roll(m, -1, 1)
    rt[:, -1] = False
    edge = m & ~(up & dn & lf & rt)
    ys, xs = np.nonzero(edge)
    order = np.argsort(np.arctan2(ys - ys.mean(), xs - xs.mean()))
    return xs[order].astype(float), ys[order].astype(float)


def downsample_mask(mask: np.ndarray, factor: int) -> np.ndarray:
    """Coarsen a binary mask by integer block-mean + 0.5 threshold.

    Unlike cartographic generalization, this IS a nested discretization:
    the coarse mask is a low-passed version of the fine one, which is what
    makes the cell-boundary reproducibility test cleaner than coastlines.
    """
    h, w = mask.shape
    h2, w2 = h // factor * factor, w // factor * factor
    block = mask[:h2, :w2].reshape(h2 // factor, factor, w2 // factor, factor)
    return block.mean(axis=(1, 3)) > 0.5


def cell_curve(mask: np.ndarray, num_points: int = 4096) -> PlanarCurve:
    """Binary mask -> uniform-arc-length PlanarCurve in pixel units."""
    xs, ys = mask_boundary_ring(mask)
    return PlanarCurve.from_points(xs, ys, num_points)


def plateau_sigma_star(
    sweep: Sequence[ContourMeasurement], smooth_window: int = 5
) -> float:
    """Data-driven sigma* = argmin |d ln P / d ln sigma| on the interior.

    Real-coastline slopes are noisier than synthetic carrier sweeps, so the
    log-slope is smoothed with a small uniform window before the argmin
    (same estimator as estimate_sigma_star, plus smoothing).
    """
    obs = sorted(sweep, key=lambda o: o.sigma)
    log_s = np.log(np.array([o.sigma for o in obs]))
    log_p = np.log(np.array([o.perimeter for o in obs]))
    slope = np.abs(np.gradient(log_p, log_s))
    if smooth_window > 1:
        kernel = np.ones(smooth_window) / smooth_window
        slope = np.convolve(slope, kernel, mode="same")
    interior = np.argsort(slope[3:-3]) + 3
    return float(obs[interior[0]].sigma)
