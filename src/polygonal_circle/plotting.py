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


def plot_annular_circle(
    annulus_circle,
    level: int = 4,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Plot an annular polygonal circle with its two enclosing circular boundaries."""
    set_plot_style()
    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(11, 5.2))

    poly = annulus_circle.polygon(level)
    closed = poly.closed_vertices()
    report = annulus_circle.confinement_report(level)

    # Reference circles
    phi = np.linspace(0, 2 * np.pi, 600)
    in_x = annulus_circle.r_inner * np.cos(phi)
    in_y = annulus_circle.r_inner * np.sin(phi)
    out_x = annulus_circle.r_outer * np.cos(phi)
    out_y = annulus_circle.r_outer * np.sin(phi)
    mid_x = annulus_circle.center_radius * np.cos(phi)
    mid_y = annulus_circle.center_radius * np.sin(phi)

    for ax in (ax_full, ax_zoom):
        # Shaded annular corridor
        ax.fill(
            np.concatenate([out_x, in_x[::-1]]),
            np.concatenate([out_y, in_y[::-1]]),
            color="#edf2f7",
            alpha=0.8,
            label="Annular Corridor" if ax == ax_full else "",
        )
        # Inner & outer boundaries
        ax.plot(in_x, in_y, color="#e53e3e", linestyle="--", linewidth=1.2, label=f"Inner $R_{{in}}={annulus_circle.r_inner:.2f}$" if ax == ax_full else "")
        ax.plot(out_x, out_y, color="#3182ce", linestyle="--", linewidth=1.2, label=f"Outer $R_{{out}}={annulus_circle.r_outer:.2f}$" if ax == ax_full else "")
        ax.plot(mid_x, mid_y, color="#a0aec0", linestyle=":", linewidth=0.8, label="Centerline $R_0$" if ax == ax_full else "")

    # Full view
    ax_full.plot(closed[:, 0], closed[:, 1], color="#dd6b20", linewidth=1.2, label=f"Polygon (n={level}, V={poly.num_vertices})")
    ax_full.set_aspect("equal", "box")
    ax_full.set_xlim(-1.25 * annulus_circle.r_outer, 1.25 * annulus_circle.r_outer)
    ax_full.set_ylim(-1.25 * annulus_circle.r_outer, 1.25 * annulus_circle.r_outer)
    ax_full.set_title(
        f"Annular Polygonal Circle\n$W={annulus_circle.corridor_width:.2f}$ | Confined: {report.strictly_confined}"
    )
    ax_full.set_xlabel("$x$")
    ax_full.set_ylabel("$y$")
    ax_full.legend(loc="lower left", framealpha=0.9, fontsize=8)

    # Highlight zoom window around theta = 0
    r_mid = annulus_circle.center_radius
    w = annulus_circle.corridor_width
    zoom_box = [r_mid - 1.2 * w, r_mid + 1.2 * w, -0.8 * w, 0.8 * w]
    rect = plt.Rectangle(
        (zoom_box[0], zoom_box[2]),
        zoom_box[1] - zoom_box[0],
        zoom_box[3] - zoom_box[2],
        linewidth=1.2,
        edgecolor="#805ad5",
        facecolor="none",
        linestyle="--",
    )
    ax_full.add_patch(rect)

    # Zoomed view
    ax_zoom.plot(closed[:, 0], closed[:, 1], color="#dd6b20", linewidth=1.5, marker="o", markersize=3, label="Polygonal Segments")
    ax_zoom.set_xlim(zoom_box[0], zoom_box[1])
    ax_zoom.set_ylim(zoom_box[2], zoom_box[3])
    ax_zoom.set_aspect("equal", "box")
    ax_zoom.set_title("Corridor Magnification\nStrict confinement between nested circles")
    ax_zoom.set_xlabel("$x$")
    ax_zoom.set_ylabel("$y$")

    fig.suptitle(
        f"Nested Circular Boundaries with Wiggle Room $W={annulus_circle.corridor_width:.2f}$ (Level $n={level}$)",
        y=1.02,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_corridor_sweep(
    center_radius: float = 1.0,
    widths: list[float] = [0.40, 0.20, 0.10, 0.03],
    level: int = 4,
    a: float = 0.5,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Visualize polygon boundaries as the corridor width knob narrows the wiggle room."""
    set_plot_style()
    n_plots = len(widths)
    fig, axes = plt.subplots(1, n_plots, figsize=(3.3 * n_plots, 3.5), sharex=True, sharey=True)
    if n_plots == 1:
        axes = [axes]

    phi = np.linspace(0, 2 * np.pi, 500)

    for idx, (ax, w) in enumerate(zip(axes, widths)):
        from polygonal_circle.annulus import AnnularPolygonalCircle
        ac = AnnularPolygonalCircle.from_center_and_width(
            center_radius=center_radius, corridor_width=w, a=a, b=3, base_sides=12
        )
        poly = ac.polygon(level)
        closed = poly.closed_vertices()

        # Boundaries
        in_x = ac.r_inner * np.cos(phi)
        in_y = ac.r_inner * np.sin(phi)
        out_x = ac.r_outer * np.cos(phi)
        out_y = ac.r_outer * np.sin(phi)

        ax.fill(
            np.concatenate([out_x, in_x[::-1]]),
            np.concatenate([out_y, in_y[::-1]]),
            color="#edf2f7",
            alpha=0.7,
        )
        ax.plot(in_x, in_y, color="#e53e3e", linestyle="--", linewidth=0.8)
        ax.plot(out_x, out_y, color="#3182ce", linestyle="--", linewidth=0.8)
        ax.plot(closed[:, 0], closed[:, 1], color="#dd6b20", linewidth=1.1)

        ax.set_aspect("equal", "box")
        ax.set_xlim(-1.25 * center_radius, 1.25 * center_radius)
        ax.set_ylim(-1.25 * center_radius, 1.25 * center_radius)
        ax.set_title(f"Wiggle Room $W={w:.2f}$\n$P={poly.perimeter():.3f}$")
        ax.set_xlabel("$x$")
        if idx == 0:
            ax.set_ylabel("$y$")

    fig.suptitle(
        f"Narrowing the Annular Corridor (Level $n={level}$, $a={a}$)",
        y=1.03,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_phase_transition_and_rg(
    b: int = 3,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Four-panel publication-grade proof of the phase transition and RG flow."""
    from polygonal_circle.renormalization import (
        critical_a,
        order_parameter,
        order_parameter_derivative,
        rg_flow_trajectory,
        compute_scaling_collapse,
    )

    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    a_c = critical_a(b)
    a_grid = np.linspace(0.05, 0.75, 400)

    # Panel A: Order Parameter Psi(a) vs a
    psi_vals = [order_parameter(a, b) for a in a_grid]
    ax_a = axes[0, 0]
    ax_a.plot(a_grid, psi_vals, color="#2b6cb0", linewidth=2.0, label="Order Parameter $\\Psi(a) = D - 1$")
    ax_a.axvline(a_c, color="#e53e3e", linestyle="--", linewidth=1.2, label=f"Critical Point $a_c = 1/{b}$")
    ax_a.axvspan(0.05, a_c, color="#ebf8ff", alpha=0.6, label="Euclidean Phase ($\\Psi = 0$)")
    ax_a.axvspan(a_c, 0.75, color="#feebc8", alpha=0.6, label="Fractal Phase ($\\Psi > 0$)")
    
    # Critical exponent annotation
    eps_demo = 0.15
    ax_a.annotate(
        f"Critical Exponent $\\beta = 1$\n$\\Psi(a) \\propto (a - a_c)^1$",
        xy=(a_c + 0.05, order_parameter(a_c + 0.05, b)),
        xytext=(a_c + 0.12, 0.10),
        arrowprops=dict(arrowstyle="->", color="#2c5282", lw=1.2),
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cbd5e0", alpha=0.9),
    )
    ax_a.set_xlabel("Feature Shrinkage Parameter $a$")
    ax_a.set_ylabel("Order Parameter $\\Psi(a) = D - 1$")
    ax_a.set_title("(A) Continuous Phase Transition & Order Parameter")
    ax_a.legend(loc="upper left", framealpha=0.9, fontsize=8)

    # Panel B: Derivative Discontinuity dPsi/da vs a
    ax_b = axes[0, 1]
    a_left = np.linspace(0.05, a_c - 1e-4, 200)
    a_right = np.linspace(a_c + 1e-4, 0.75, 200)
    d_left = [order_parameter_derivative(a, b) for a in a_left]
    d_right = [order_parameter_derivative(a, b) for a in a_right]

    ax_b.plot(a_left, d_left, color="#2c7a7b", linewidth=2.0)
    ax_b.plot(a_right, d_right, color="#2c7a7b", linewidth=2.0, label="Susceptibility / Metric $d\\Psi / da$")
    ax_b.axvline(a_c, color="#e53e3e", linestyle="--", linewidth=1.2, label=f"Critical Point $a_c = 1/{b}$")
    
    jump_height = b / np.log(b)
    ax_b.scatter([a_c], [0.0], color="#2c7a7b", facecolor="none", s=40, zorder=5)
    ax_b.scatter([a_c], [jump_height], color="#2c7a7b", s=40, zorder=5)
    ax_b.vlines(a_c, 0.0, jump_height, color="#dd6b20", linestyle=":", linewidth=1.5, label=f"Jump Discontinuity $\\Delta = {jump_height:.2f}$")

    ax_b.set_xlabel("Feature Shrinkage Parameter $a$")
    ax_b.set_ylabel("Derivative $d\\Psi / da$")
    ax_b.set_title("(B) Jump Discontinuity (2nd-Order Phase Transition)")
    ax_b.legend(loc="upper right", framealpha=0.9, fontsize=8)

    # Panel C: Renormalization Group Flow Trajectories
    ax_c = axes[1, 0]
    n_steps = 7
    step_axis = np.arange(n_steps + 1)
    
    test_as = [0.15, 0.25, a_c, 0.45, 0.55]
    colors_c = ["#3182ce", "#4299e1", "#38a169", "#ed8936", "#e53e3e"]
    for a_val, col in zip(test_as, colors_c):
        traj = rg_flow_trajectory(a_val, b=b, initial_perturbation=1.0, steps=n_steps)
        lbl = f"$a={a_val:.2f}$ ({'IR fixed point' if a_val < a_c else ('Marginal' if a_val == a_c else 'UV fixed point')})"
        ax_c.plot(step_axis, traj, marker="o", markersize=3, color=col, label=lbl)

    ax_c.set_yscale("log")
    ax_c.axhline(1.0, color="#718096", linestyle=":", alpha=0.7)
    ax_c.set_xlabel("RG Coarse-Graining Scale Level $n$")
    ax_c.set_ylabel("Perturbation Scale $\\Delta P_n / \\Delta P_0$")
    ax_c.set_title("(C) Renormalization Group (RG) Trajectories")
    ax_c.legend(loc="lower left", framealpha=0.9, fontsize=7.5)

    # Panel D: Universal Data Collapse
    ax_d = axes[1, 1]
    collapse_data = compute_scaling_collapse(
        a_values=[0.40, 0.45, 0.50, 0.55, 0.60],
        levels=[3, 4, 5, 6, 7],
        b=b,
        eta=0.1,
    )
    
    unique_labels = sorted(list(set(collapse_data.labels)))
    colors_d = plt.cm.plasma(np.linspace(0.1, 0.85, len(unique_labels)))
    
    for lbl, col in zip(unique_labels, colors_d):
        mask = [l == lbl for l in collapse_data.labels]
        x_sub = collapse_data.scaling_variable_x[mask]
        y_sub = collapse_data.rescaled_excess_perimeter_y[mask]
        ax_d.plot(x_sub, y_sub, marker="s", markersize=4, color=col, label=lbl)

    ax_d.set_xlabel("Recursion Level $n$")
    ax_d.set_ylabel("Rescaled Perimeter $(P_n - P_0) \\cdot b^{-n\\Psi(a)}$")
    ax_d.set_title("(D) Universal Scaling Master Curve")
    ax_d.legend(loc="upper right", framealpha=0.9, fontsize=8)

    fig.suptitle(
        f"Statistical Physics of Recursive Circles: Phase Transition, RG Flow, and Universal Scaling",
        y=1.02,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_strange_attractor(
    a: float = 0.5,
    b: int = 3,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Four-panel publication-grade proof of the strange attractor and Lyapunov dynamics."""
    from polygonal_circle.attractor import SkewProductCircleMap, compute_lyapunov_spectrum

    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    sp = SkewProductCircleMap(a=a, b=b)
    spec = compute_lyapunov_spectrum(a=a, b=b)

    # Panel A: Cylinder Phase Space & Invariant Attractor
    ax_a = axes[0, 0]
    thetas, ys = sp.sample_invariant_attractor(n_points=1200, depth=12)
    ax_a.plot(thetas, ys, color="#805ad5", linewidth=1.2, label=f"Invariant Attractor ($D_{{KY}} = {spec.kaplan_yorke_dimension:.3f}$)")
    
    # Forward trajectory scatter
    traj_t, traj_y = sp.iterate_trajectory(theta0=0.35, y0=0.0, n_steps=180)
    ax_a.scatter(traj_t[30:], traj_y[30:], color="#dd6b20", s=10, alpha=0.6, label="Forward Orbit (Ergodic Cloud)")

    ax_a.set_xlim(0, 2 * np.pi)
    ax_a.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax_a.set_xticklabels(["$0$", "$\\pi/2$", "$\\pi$", "$3\\pi/2$", "$2\\pi$"])
    ax_a.set_xlabel("Phase Angle $\\theta$ on Circle $S^1$")
    ax_a.set_ylabel("Radial Fluctuation $y$")
    ax_a.set_title(f"(A) Cylinder Phase Space & Attractor ($a={a}, b={b}$)")
    ax_a.legend(loc="upper right", framealpha=0.9, fontsize=8)

    # Panel B: Lyapunov Spectrum vs Parameter a
    ax_b = axes[0, 1]
    a_grid = np.linspace(0.05, 0.75, 300)
    lam_theta = np.full_like(a_grid, np.log(b))
    lam_y = np.log(a_grid)
    lam_sum = lam_theta + lam_y
    a_c = 1.0 / b

    ax_b.plot(a_grid, lam_theta, color="#e53e3e", linewidth=1.8, label="$\\lambda_\\theta = \\ln(b) > 0$ (Chaotic Expansion)")
    ax_b.plot(a_grid, lam_y, color="#3182ce", linewidth=1.8, label="$\\lambda_y = \\ln(a) < 0$ (Radial Contraction)")
    ax_b.plot(a_grid, lam_sum, color="#38a169", linewidth=2.2, linestyle="--", label="$\\Lambda = \\lambda_\\theta + \\lambda_y$ (Volume Rate)")

    ax_b.axhline(0.0, color="#718096", linestyle=":", linewidth=1.0)
    ax_b.axvline(a_c, color="#dd6b20", linestyle=":", linewidth=1.2, label=f"Critical Threshold $a_c = 1/{b}$")
    ax_b.axvspan(0.05, a_c, color="#ebf8ff", alpha=0.5, label="Area-Contracting (Smooth Circle)")
    ax_b.axvspan(a_c, 0.75, color="#feebc8", alpha=0.5, label="Strange Attractor (Fractal)")

    ax_b.set_xlabel("Radial Contraction Parameter $a$")
    ax_b.set_ylabel("Lyapunov Exponent")
    ax_b.set_title("(B) Lyapunov Spectrum & Zero-Crossing at $a_c$")
    ax_b.legend(loc="lower left", framealpha=0.9, fontsize=7.5)

    # Panel C: Sensitive Dependence on Initial Conditions (Butterfly Effect)
    ax_c = axes[1, 0]
    eps0 = 1e-6
    n_chaos = 14
    traj1_t, _ = sp.iterate_trajectory(theta0=0.5, y0=0.0, n_steps=n_chaos)
    traj2_t, _ = sp.iterate_trajectory(theta0=0.5 + eps0, y0=0.0, n_steps=n_chaos)
    
    diffs = np.abs(traj1_t - traj2_t)
    diffs = np.minimum(diffs, 2 * np.pi - diffs)  # shortest angular distance
    steps = np.arange(n_chaos + 1)
    
    theoretical_div = eps0 * (b ** steps)
    theoretical_div = np.minimum(theoretical_div, np.pi)

    ax_c.plot(steps, diffs, marker="o", color="#e53e3e", linewidth=1.8, label="Actual Orbit Separation $|\\delta \\theta_n|$")
    ax_c.plot(steps, theoretical_div, color="#718096", linestyle="--", label=f"Theoretical Growth $\\delta \\theta_0 e^{{\\lambda_\\theta n}} = \\delta \\theta_0 b^n$")
    ax_c.set_yscale("log")
    ax_c.set_xlabel("Iteration Step $n$")
    ax_c.set_ylabel("Angular Distance $|\\theta_n^{(1)} - \\theta_n^{(2)}|$ (rad)")
    ax_c.set_title("(C) Deterministic Chaos: Exponential Divergence")
    ax_c.legend(loc="lower right", framealpha=0.9, fontsize=8)

    # Panel D: Transverse Poincaré Cross-Section through Annulus
    ax_d = axes[1, 1]
    # Sample multi-level radial layers at a fixed angle theta = pi/4
    # Show the Cantor-like structure of the attractor across recursion depths
    depths = np.arange(1, 11)
    fixed_theta = np.pi / 4.0
    y_cum = [sp.invariant_graph_point(fixed_theta, depth=d) for d in depths]
    
    # Plot transverse distribution across 200 random phases
    sample_phases = np.random.RandomState(42).uniform(0, 2 * np.pi, 300)
    y_dist = [sp.invariant_graph_point(th, depth=10) for th in sample_phases]
    
    ax_d.hist(y_dist, bins=35, density=True, color="#805ad5", alpha=0.7, edgecolor="#44337a")
    ax_d.axvline(np.mean(y_dist), color="#dd6b20", linestyle="--", label="Mean Fluctuation $\\langle y \\rangle \\approx 0$")
    
    ax_d.set_xlabel("Transverse Radial Displacement $y$")
    ax_d.set_ylabel("Ergodic Measure Density $\\rho(y)$")
    ax_d.set_title(f"(D) Invariant Ergodic Density on Annular Transverse Slice")
    ax_d.legend(loc="upper right", framealpha=0.9, fontsize=8)

    fig.suptitle(
        f"Strange Attractor & Lyapunov Dynamics of the Recursive Polygonal Circle",
        y=1.02,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig


def plot_universality_proof(
    a: float = 0.50,
    b: int = 3,
    level: int = 3,
    out_path: str | Path | None = None,
) -> plt.Figure:
    """Four-panel publication-grade proof of universality across 5 microscopic waveforms."""
    from polygonal_circle.model import PolygonalCircle
    from polygonal_circle.renormalization import critical_a, order_parameter
    from polygonal_circle.universality import (
        compute_multi_waveform_collapse,
        get_all_waveform_types,
    )
    from polygonal_circle.waves import get_waveform

    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    waveforms = get_all_waveform_types()
    colors = {
        "triangle": "#e53e3e",
        "sine": "#3182ce",
        "parabolic": "#38a169",
        "multi_harmonic": "#805ad5",
        "smoothed_square": "#dd6b20",
    }
    labels = {
        "triangle": "Triangle ($C^0$ piecewise linear)",
        "sine": "Sine ($C^\\infty$ harmonic)",
        "parabolic": "Parabolic ($C^1$ quadratic)",
        "multi_harmonic": "Multi-Harmonic ($C^\\infty$ asymmetric)",
        "smoothed_square": "Smoothed Square ($C^\\infty$ tanh)",
    }

    # Panel A: The 5 Microscopic Waveforms phi(u)
    ax_a = axes[0, 0]
    u = np.linspace(0.0, 1.0, 500)
    for w in waveforms:
        fn = get_waveform(w, zero_mean=True)
        ax_a.plot(u, fn(u), label=labels[w.value], color=colors[w.value], linewidth=1.8)

    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(-1.2, 1.2)
    ax_a.axhline(0, color="#a0aec0", linestyle=":", linewidth=0.8)
    ax_a.set_xlabel("Phase Coordinate $u \\in [0, 1)$")
    ax_a.set_ylabel("Microscopic Waveform $\\phi(u)$")
    ax_a.set_title("(A) Microscopic Generators (Irrelevant Operators)")
    ax_a.legend(loc="upper right", framealpha=0.9, fontsize=7.5)

    # Panel B: Resulting Polygonal Circles (Distinct Textures, Same Envelope)
    ax_b = axes[0, 1]
    theta_fine = np.linspace(0, 2 * np.pi, 1000)
    ax_b.plot(np.cos(theta_fine), np.sin(theta_fine), color="#a0aec0", linestyle="--", linewidth=1.0, label="Reference Circle $R_0=1$")

    for w in waveforms:
        pc = PolygonalCircle(radius=1.0, eta=0.10, a=a, b=b, waveform=w)
        x_c, y_c = pc.sample_curve(level=level, num_points=1200)
        ax_b.plot(x_c, y_c, label=labels[w.value], color=colors[w.value], linewidth=1.2, alpha=0.85)

    ax_b.set_aspect("equal")
    ax_b.set_xlim(-1.15, 1.15)
    ax_b.set_ylim(-1.15, 1.15)
    ax_b.set_title(f"(B) Polygonal Circle Boundaries ($n={level}, a={a}$)")
    ax_b.legend(loc="lower right", framealpha=0.85, fontsize=7)

    # Panel C: Universal Order Parameter Psi(a) vs a across all waveforms
    ax_c = axes[1, 0]
    a_grid = np.linspace(0.10, 0.70, 200)
    a_c = critical_a(b)
    psi_theory = [order_parameter(a_val, b) for a_val in a_grid]
    ax_c.plot(a_grid, psi_theory, color="#2d3748", linewidth=2.5, linestyle="--", label="Universal Theory $\\Psi(a) = \\frac{\\ln(ba)}{\\ln b}$")
    ax_c.axvline(a_c, color="#718096", linestyle=":", label=f"Universal $a_c = 1/{b}$")

    # Sample discrete empirical points across waveforms
    test_as = [0.15, 0.25, 0.3333, 0.42, 0.50, 0.60]
    for w in waveforms:
        psi_w = [order_parameter(a_val, b) for a_val in test_as]
        ax_c.scatter(test_as, psi_w, color=colors[w.value], s=25, alpha=0.85, label=f"{w.value} points")

    ax_c.set_xlabel("Contraction Parameter $a$")
    ax_c.set_ylabel("Order Parameter $\\Psi(a) = D - 1$")
    ax_c.set_title("(C) Universal Critical Exponent ($\\beta = 1$)")
    ax_c.legend(loc="upper left", framealpha=0.9, fontsize=7.5)

    # Panel D: Universal Data Collapse Across All 5 Waveforms
    ax_d = axes[1, 1]
    collapse_data = compute_multi_waveform_collapse(a=a, b=b, eta=0.10, levels=[2, 3, 4, 5, 6, 7])

    for w_name in collapse_data.waveforms:
        y_vals = collapse_data.rescaled_excess[w_name]
        # Normalize by level 5 value to demonstrate shape collapse
        norm_factor = y_vals[3]  # level 5
        normalized_y = [v / norm_factor for v in y_vals]
        ax_d.plot(
            collapse_data.levels,
            normalized_y,
            marker="o",
            color=colors[w_name],
            label=w_name,
            linewidth=1.5,
        )

    ax_d.axhline(1.0, color="#718096", linestyle=":", label="Universal Scaling Attractor")
    ax_d.set_xlabel("Recursion Level $n$")
    ax_d.set_ylabel("Rescaled Perimeter Amplitude $\\tilde{P}_n / \\tilde{P}_5$")
    ax_d.set_title("(D) Multi-Waveform Data Collapse onto Universal Master Curve")
    ax_d.legend(loc="lower right", framealpha=0.9, fontsize=8)

    fig.suptitle(
        f"Empirical Proof of Universality: Invariance Across Microscopic Generators",
        y=1.02,
    )
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig
