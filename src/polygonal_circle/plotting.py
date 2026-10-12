"""Publication-quality plotting for recursive polygonal circles and fractal analysis."""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from polygonal_circle.analysis import PolygonalCircleAnalysis
from polygonal_circle.model import PolygonalCircle


def set_plot_style():
    """Configure clean scientific aesthetics."""
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 13,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": ":",
        }
    )


def plot_multilevel_shapes(
    circle: PolygonalCircle,
    levels: list[int] = [0, 1, 2, 3, 4],
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Side-by-side visualization of polygon boundary evolution across levels."""
    set_plot_style()
    n_plots = len(levels)
    fig, axes = plt.subplots(1, n_plots, figsize=(3.2 * n_plots, 3.4), sharex=True, sharey=True)
    if n_plots == 1:
        axes = [axes]

    # Reference circle
    phi = np.linspace(0, 2 * np.pi, 500)
    circ_x = circle.radius * np.cos(phi)
    circ_y = circle.radius * np.sin(phi)

    colors = ["#2b6cb0", "#2c7a7b", "#d69e2e", "#dd6b20", "#e53e3e", "#805ad5"]

    for idx, (ax, lvl) in enumerate(zip(axes, levels)):
        poly = circle.polygon(lvl)
        closed_v = poly.closed_vertices()
        color = colors[idx % len(colors)]

        # Background circular envelope
        ax.plot(circ_x, circ_y, color="#cbd5e0", linestyle="--", linewidth=1.0, label="Envelope $R$")

        # Polygon boundary
        ax.plot(closed_v[:, 0], closed_v[:, 1], color=color, linewidth=1.2, label=f"Level {lvl}")
        if poly.num_vertices <= 60:
            ax.scatter(poly.x, poly.y, s=12, color=color, zorder=4)

        ax.set_aspect("equal", "box")
        ax.set_title(f"Level $n={lvl}$\n$V={poly.num_vertices}, P={poly.perimeter():.3f}$")
        ax.set_xlim(-1.15 * circle.radius, 1.15 * circle.radius)
        ax.set_ylim(-1.15 * circle.radius, 1.15 * circle.radius)
        if idx == 0:
            ax.set_ylabel("$y$")
        ax.set_xlabel("$x$")

    fig.suptitle(
        f"Recursive Polygonal Circle ($a={circle.a}, \\eta={circle.eta}, b={circle.b}$)",
        y=1.03,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_scaling_regimes(
    a_values: list[float] = [0.2, 1.0 / 3.0, 0.5, 0.6],
    eta: float = 0.1,
    b: int = 3,
    max_level: int = 7,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Four-panel comparison of geometric behavior across scaling regimes."""
    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    palette = {
        0.2: ("#3182ce", "$a=0.20 < 1/3$ (Convergent)"),
        1.0 / 3.0: ("#38a169", "$a=1/3$ (Critical)"),
        0.5: ("#dd6b20", "$a=0.50 > 1/3$ (Divergent, $3a=1.5$)"),
        0.6: ("#e53e3e", "$a=0.60 > 1/3$ (Divergent, $3a=1.8$)"),
    }

    for a in a_values:
        pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=b)
        analysis = PolygonalCircleAnalysis(pc)
        df = analysis.metrics_dataframe(max_level=max_level)

        color, label = palette.get(a, ("#718096", f"$a={a:.2f}$"))

        # Panel 1: Perimeter vs level
        axes[0, 0].plot(df["level"], df["perimeter"], marker="o", markersize=4, color=color, label=label)

        # Panel 2: Perimeter ratio P_{n+1} / P_n vs level
        sub_ratio = df.dropna(subset=["perimeter_ratio"])
        axes[0, 1].plot(sub_ratio["level"], sub_ratio["perimeter_ratio"], marker="s", markersize=4, color=color, label=label)

        # Panel 3: Enclosed Area vs level
        axes[1, 0].plot(df["level"], df["area"], marker="^", markersize=4, color=color, label=label)

        # Panel 4: Isoperimetric Quotient 4*pi*A / P^2 vs level
        axes[1, 1].plot(df["level"], df["isoperimetric_quotient"], marker="d", markersize=4, color=color, label=label)

    # Reference lines
    axes[0, 0].axhline(2.0 * np.pi, color="#718096", linestyle="--", alpha=0.7, label="Smooth Circle $2\\pi R$")
    axes[0, 0].set_xlabel("Recursion Level $n$")
    axes[0, 0].set_ylabel("Perimeter $P_n$")
    axes[0, 0].set_title("(A) Boundary Perimeter Evolution")
    axes[0, 0].legend(loc="upper left", framealpha=0.9)

    axes[0, 1].axhline(1.0, color="#718096", linestyle=":", alpha=0.7)
    axes[0, 1].set_xlabel("Recursion Level $n$")
    axes[0, 1].set_ylabel("Ratio $P_n / P_{n-1}$")
    axes[0, 1].set_title("(B) Perimeter Growth Ratio (Asymptote $\\to b \\cdot a$)")
    axes[0, 1].legend(loc="upper left", framealpha=0.9)

    axes[1, 0].axhline(np.pi, color="#718096", linestyle="--", alpha=0.7, label="Circle Area $\\pi R^2$")
    axes[1, 0].set_xlabel("Recursion Level $n$")
    axes[1, 0].set_ylabel("Enclosed Area $A_n$")
    axes[1, 0].set_title("(C) Area Boundedness ($A_n < \\pi R^2$)")
    axes[1, 0].legend(loc="lower right", framealpha=0.9)

    axes[1, 1].axhline(1.0, color="#718096", linestyle="--", alpha=0.7, label="Euclidean Circle ($Q=1$)")
    axes[1, 1].set_xlabel("Recursion Level $n$")
    axes[1, 1].set_ylabel("Isoperimetric Quotient $4\\pi A_n / P_n^2$")
    axes[1, 1].set_title("(D) Isoperimetric Collapse ($Q_n \\to 0$ for Fractals)")
    axes[1, 1].legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_richardson_crossover(
    circle: PolygonalCircle,
    max_level: int = 10,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Richardson coastline scaling plot showing crossover from Euclidean to fractal regimes."""
    set_plot_style()
    analysis = PolygonalCircleAnalysis(circle)
    df = analysis.metrics_dataframe(max_level=max_level)
    sub = df[df["level"] >= 1].copy()

    eps = sub["mean_segment_length"].values
    perim = sub["perimeter"].values
    log_inv_eps = np.log(1.0 / eps)
    log_p = np.log(perim)

    # Fits
    # Coarse regime: levels 1 to 4
    mask_coarse = (sub["level"] >= 1) & (sub["level"] <= 4)
    res_coarse = np.polyfit(log_inv_eps[mask_coarse], log_p[mask_coarse], deg=1)
    d_coarse = 1.0 + res_coarse[0]

    # Asymptotic fractal regime: levels 6 to 10
    mask_fine = (sub["level"] >= 6) & (sub["level"] <= max_level)
    res_fine = np.polyfit(log_inv_eps[mask_fine], log_p[mask_fine], deg=1)
    d_fine = 1.0 + res_fine[0]

    theo_d = circle.scaling.theoretical_dimension

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(log_inv_eps, log_p, color="#2b6cb0", s=35, zorder=5, label="Empirical Measurements")

    for i, row in sub.iterrows():
        ax.annotate(
            f"n={int(row['level'])}",
            (log_inv_eps[i - 1], log_p[i - 1]),
            textcoords="offset points",
            xytext=(6, -2),
            fontsize=8,
            color="#4a5568",
        )

    # Plot regression lines
    grid_coarse = np.linspace(log_inv_eps[mask_coarse].min() - 0.2, log_inv_eps[mask_coarse].max() + 0.5, 50)
    ax.plot(
        grid_coarse,
        np.polyval(res_coarse, grid_coarse),
        color="#718096",
        linestyle="--",
        label=f"Coarse scale fit: $D \\approx {d_coarse:.3f}$ (Euclidean circle)",
    )

    grid_fine = np.linspace(log_inv_eps[mask_fine].min() - 0.5, log_inv_eps[mask_fine].max() + 0.3, 50)
    ax.plot(
        grid_fine,
        np.polyval(res_fine, grid_fine),
        color="#e53e3e",
        linestyle="-",
        linewidth=1.8,
        label=f"Fine scale fit: $D \\approx {d_fine:.3f}$ (Theory: $D={theo_d:.3f}$)",
    )

    ax.set_xlabel("Resolution $\\log(1 / \\bar{\\ell}_n)$")
    ax.set_ylabel("Measured Perimeter $\\log(P_n)$")
    ax.set_title(
        f"Richardson Coastline Scaling & Scale Crossover ($a={circle.a}, \\eta={circle.eta}, b={circle.b}$)"
    )
    ax.legend(loc="upper left", framealpha=0.9)

    plt.tight_layout()
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_macroscopic_vs_microscopic(
    circle: PolygonalCircle,
    level: int = 5,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Dual-panel figure contrasting macroscopic circular appearance with microscopic polygonal facets."""
    set_plot_style()
    fig, (ax_macro, ax_micro) = plt.subplots(1, 2, figsize=(10, 4.8))

    poly = circle.polygon(level)
    closed = poly.closed_vertices()

    # Macro panel
    phi = np.linspace(0, 2 * np.pi, 500)
    ax_macro.plot(circle.radius * np.cos(phi), circle.radius * np.sin(phi), color="#cbd5e0", linestyle="--", label="Envelope $R=1$")
    ax_macro.plot(closed[:, 0], closed[:, 1], color="#2b6cb0", linewidth=1.0, label=f"Boundary (n={level}, V={poly.num_vertices})")
    ax_macro.set_aspect("equal", "box")
    ax_macro.set_title(f"Macroscopic View\nAppears smooth and circular ($R=1$)")
    ax_macro.set_xlabel("$x$")
    ax_macro.set_ylabel("$y$")
    ax_macro.legend(loc="lower left")

    # Zoom window at theta = 0: x in [0.85, 1.05], y in [-0.15, 0.15]
    zoom_box = [0.85, 1.02, -0.12, 0.12]
    rect = plt.Rectangle(
        (zoom_box[0], zoom_box[2]),
        zoom_box[1] - zoom_box[0],
        zoom_box[3] - zoom_box[2],
        linewidth=1.2,
        edgecolor="#e53e3e",
        facecolor="none",
        linestyle="--",
    )
    ax_macro.add_patch(rect)

    # Micro panel
    ax_micro.plot(circle.radius * np.cos(phi), circle.radius * np.sin(phi), color="#cbd5e0", linestyle="--", label="Envelope $R=1$")
    ax_micro.plot(closed[:, 0], closed[:, 1], color="#2b6cb0", linewidth=1.5, marker="o", markersize=3, label="Polygonal Segments")
    ax_micro.set_xlim(zoom_box[0], zoom_box[1])
    ax_micro.set_ylim(zoom_box[2], zoom_box[3])
    ax_micro.set_aspect("equal", "box")
    ax_micro.set_title(f"Microscopic Magnification (Zoom 10x)\nReveals nested polygonal facets")
    ax_micro.set_xlabel("$x$")
    ax_micro.set_ylabel("$y$")
    ax_micro.legend(loc="upper left")

    fig.suptitle(
        f"Scale-Dependent Geometry: Macroscopic Circle vs. Microscopic Fractal Boundary",
        y=1.02,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig
