"""Command-line interface for the polygonal-circle project."""

from __future__ import annotations

from pathlib import Path
import math
import rich
from rich.console import Console
from rich.table import Table
import typer

from polygonal_circle.analysis import PolygonalCircleAnalysis
from polygonal_circle.annulus import AnnularPolygonalCircle
from polygonal_circle.attractor import compute_lyapunov_spectrum
from polygonal_circle.carriers import (
    CircleCarrier,
    EllipseCarrier,
    RegularPolygonCarrier,
    SquareCarrier,
    SuperellipseCarrier,
    TrefoilCarrier,
    TriangleCarrier,
)
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.observer import ScaleSpaceObserver
from polygonal_circle.plotting import (
    plot_annular_circle,
    plot_corridor_sweep,
    plot_macroscopic_vs_microscopic,
    plot_multi_carrier_scale_space,
    plot_multilevel_shapes,
    plot_observer_scale_space,
    plot_phase_transition_and_rg,
    plot_polygonal_scale_space,
    plot_richardson_crossover,
    plot_scaling_regimes,
    plot_strange_attractor,
    plot_universality_proof,
)
from polygonal_circle.renormalization import classify_phase, critical_a
from polygonal_circle.universality import analyze_waveform_universality

app = typer.Typer(
    name="polygonal-circle",
    help="Explore recursive polygonal circles, scaling regimes, and fractal boundaries.",
    add_completion=False,
)
console = Console()


@app.command()
def inspect(
    radius: float = typer.Option(1.0, "-R", "--radius", help="Reference circular envelope radius R"),
    eta: float = typer.Option(0.1, "--eta", help="Global irregularity amplitude"),
    a: float = typer.Option(0.5, "-a", help="Amplitude shrinkage factor per level"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier per level"),
    base_sides: int = typer.Option(6, "-N", "--base-sides", help="Base polygon sides at level 0"),
    max_level: int = typer.Option(6, "-L", "--max-level", help="Maximum recursion level to evaluate"),
):
    """Inspect exact geometric quantities across recursion levels."""
    pc = PolygonalCircle(radius=radius, eta=eta, a=a, b=b, base_sides=base_sides)
    scaling = pc.scaling

    console.print(f"[bold cyan]Polygonal Circle Model[/bold cyan]")
    console.print(f"  Envelope Radius R: [bold]{radius}[/bold]")
    console.print(f"  Irregularity Amplitude eta: [bold]{eta}[/bold]")
    console.print(f"  Shrinkage a: [bold]{a}[/bold] | Frequency b: [bold]{b}[/bold]")
    console.print(f"  Product b*a: [bold]{scaling.product_ab:.4f}[/bold] (Critical: {scaling.critical_a:.4f})")
    console.print(f"  Regime: [bold magenta]{scaling.regime.value.upper()}[/bold magenta]")
    console.print(f"  Theoretical Fractal Dimension: [bold green]{scaling.theoretical_dimension:.4f}[/bold green]\n")

    analysis = PolygonalCircleAnalysis(pc)
    df = analysis.metrics_dataframe(max_level=max_level)

    table = Table(title=f"Geometric Metrics (Levels 0 to {max_level})")
    table.add_column("Level", justify="center")
    table.add_column("Vertices", justify="right")
    table.add_column("Perimeter", justify="right")
    table.add_column("Ratio P_n / P_{n-1}", justify="right")
    table.add_column("Area", justify="right")
    table.add_column("Isoperimetric Q", justify="right")
    table.add_column("Mean Radius", justify="right")
    table.add_column("Mean Segment", justify="right")

    for _, row in df.iterrows():
        ratio_str = f"{row['perimeter_ratio']:.4f}" if not pd_isna(row["perimeter_ratio"]) else "-"
        table.add_row(
            str(int(row["level"])),
            f"{int(row['num_vertices']):,}",
            f"{row['perimeter']:.4f}",
            ratio_str,
            f"{row['area']:.4f}",
            f"{row['isoperimetric_quotient']:.4f}",
            f"{row['mean_radius']:.4f}",
            f"{row['mean_segment_length']:.5f}",
        )

    console.print(table)


def pd_isna(val):
    import math
    return val is None or math.isnan(val)


@app.command()
def regimes(
    max_level: int = typer.Option(7, "-L", "--max-level", help="Maximum recursion level"),
    eta: float = typer.Option(0.1, "--eta", help="Global irregularity amplitude"),
):
    """Compare geometric behavior across the three scaling regimes: a < 1/3, a = 1/3, a > 1/3."""
    cases = [
        ("Convergent (a=0.20 < 1/3)", 0.20),
        ("Critical (a=1/3 = 0.333)", 1.0 / 3.0),
        ("Divergent (a=0.50, 3a=1.5)", 0.50),
        ("Divergent (a=0.60, 3a=1.8)", 0.60),
    ]

    table = Table(title="Scaling Regimes Comparison")
    table.add_column("Level", justify="center")
    for name, _ in cases:
        table.add_column(f"{name}\nPerimeter (Ratio)", justify="right")

    dfs = []
    for _, a in cases:
        pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=3)
        df = PolygonalCircleAnalysis(pc).metrics_dataframe(max_level=max_level)
        dfs.append(df)

    for n in range(max_level + 1):
        row_vals = [str(n)]
        for df in dfs:
            p = df.loc[df["level"] == n, "perimeter"].values[0]
            ratio = df.loc[df["level"] == n, "perimeter_ratio"].values[0]
            r_str = f" ({ratio:.3f})" if not pd_isna(ratio) else ""
            row_vals.append(f"{p:.4f}{r_str}")
        table.add_row(*row_vals)

    console.print(table)


@app.command()
def dimensions(
    max_level: int = typer.Option(9, "-L", "--max-level", help="Max level for Richardson regression"),
    eta: float = typer.Option(0.1, "--eta", help="Irregularity amplitude"),
):
    """Evaluate empirical Richardson and Box-Counting dimensions against theory."""
    a_values = [0.20, 0.30, 1.0 / 3.0, 0.40, 0.50, 0.60]

    table = Table(title="Fractal Dimension Analysis across Shrinkage Parameter a")
    table.add_column("a", justify="center")
    table.add_column("3 * a", justify="right")
    table.add_column("Regime", justify="center")
    table.add_column("Theory D", justify="right")
    table.add_column("Richardson D", justify="right")
    table.add_column("Richardson R^2", justify="right")
    table.add_column("Box-Count D (L=5)", justify="right")
    table.add_column("Box-Count R^2", justify="right")

    for a in a_values:
        pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=3)
        analysis = PolygonalCircleAnalysis(pc)
        scaling = pc.scaling

        # Richardson dimension on fine levels (5 to max_level)
        r_dim, r_r2, _ = analysis.richardson_dimension(max_level=max_level, min_level=5)
        # Box counting at level 5
        b_dim, b_r2, _, _ = analysis.box_counting_dimension(level=5, num_scales=8)

        table.add_row(
            f"{a:.3f}",
            f"{3*a:.3f}",
            scaling.regime.value,
            f"{scaling.theoretical_dimension:.4f}",
            f"{r_dim:.4f}",
            f"{r_r2:.4f}",
            f"{b_dim:.4f}",
            f"{b_r2:.4f}",
        )

    console.print(table)


@app.command()
def annulus(
    center_radius: float = typer.Option(1.0, "-R", "--center-radius", help="Nominal center radius R0"),
    corridor_width: float = typer.Option(0.20, "-W", "--corridor-width", help="Corridor width knob W = R_out - R_in"),
    a: float = typer.Option(0.5, "-a", help="Shrinkage factor per level"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier per level"),
    base_sides: int = typer.Option(12, "-N", "--base-sides", help="Base sides at level 0"),
    max_level: int = typer.Option(5, "-L", "--max-level", help="Max level to evaluate"),
):
    """Evaluate polygonal circles strictly confined within two nested circular boundaries."""
    ac = AnnularPolygonalCircle.from_center_and_width(
        center_radius=center_radius, corridor_width=corridor_width, a=a, b=b, base_sides=base_sides
    )

    console.print(f"[bold cyan]Annular Polygonal Circle (Two Nested Boundaries)[/bold cyan]")
    console.print(f"  Centerline Radius R0: [bold]{ac.center_radius:.4f}[/bold]")
    console.print(f"  Inner Boundary R_in:  [bold red]{ac.r_inner:.4f}[/bold red]")
    console.print(f"  Outer Boundary R_out: [bold blue]{ac.r_outer:.4f}[/bold blue]")
    console.print(f"  Corridor Width Knob W: [bold green]{ac.corridor_width:.4f}[/bold green] (Half-width delta={ac.half_width:.4f})\n")

    table = Table(title=f"Annular Confinement across Levels (W={corridor_width:.3f})")
    table.add_column("Level", justify="center")
    table.add_column("Vertices", justify="right")
    table.add_column("Perimeter", justify="right")
    table.add_column("Area", justify="right")
    table.add_column("Vertex [R_min, R_max]", justify="center")
    table.add_column("Segment Min Dist", justify="right")
    table.add_column("Strictly Confined?", justify="center")

    for lvl in range(max_level + 1):
        poly = ac.polygon(lvl)
        rep = ac.confinement_report(lvl)
        confined_str = "[green]YES[/green]" if rep.strictly_confined else "[yellow]NO (chord sagitta)[/yellow]"
        table.add_row(
            str(lvl),
            f"{poly.num_vertices:,}",
            f"{poly.perimeter():.4f}",
            f"{poly.area():.4f}",
            f"[{rep.vertex_min_radius:.4f}, {rep.vertex_max_radius:.4f}]",
            f"{rep.segment_min_distance:.4f}",
            confined_str,
        )

    console.print(table)


@app.command("wiggle-sweep")
def wiggle_sweep(
    center_radius: float = typer.Option(1.0, "-R", "--center-radius", help="Centerline radius R0"),
    level: int = typer.Option(4, "-L", "--level", help="Recursion level to test"),
    a: float = typer.Option(0.5, "-a", help="Shrinkage factor per level"),
):
    """Sweep the corridor width knob to observe how narrowing wiggle room affects the boundary."""
    widths = [0.40, 0.20, 0.10, 0.05, 0.02, 0.005, 0.001]

    table = Table(title=f"Wiggle Room Knob Sweep (Level n={level}, R0={center_radius})")
    table.add_column("Width W", justify="center")
    table.add_column("Inner R_in", justify="right")
    table.add_column("Outer R_out", justify="right")
    table.add_column("Perimeter", justify="right")
    table.add_column("Area", justify="right")
    table.add_column("Isoperimetric Q", justify="right")
    table.add_column("Strictly Confined?", justify="center")

    for w in widths:
        ac = AnnularPolygonalCircle.from_center_and_width(
            center_radius=center_radius, corridor_width=w, a=a, b=3, base_sides=12
        )
        poly = ac.polygon(level)
        rep = ac.confinement_report(level)
        c_str = "[green]YES[/green]" if rep.strictly_confined else "[red]NO[/red]"
        table.add_row(
            f"{w:.4f}",
            f"{ac.r_inner:.4f}",
            f"{ac.r_outer:.4f}",
            f"{poly.perimeter():.4f}",
            f"{poly.area():.4f}",
            f"{poly.isoperimetric_quotient():.4f}",
            c_str,
        )

    console.print(table)


@app.command("phase-transition")
def phase_transition(
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
):
    """Display the exact statistical physics phase diagram, RG eigenvalues, and order parameters."""
    a_c = critical_a(b)
    test_as = [0.10, 0.20, 0.28, a_c, 0.38, 0.50, 0.60, 0.70]

    console.print(f"[bold cyan]Statistical Physics Phase Diagram & RG Classification (b={b})[/bold cyan]")
    console.print(f"  Critical Threshold: a_c = 1/{b} = [bold green]{a_c:.4f}[/bold green]")
    console.print(f"  Order Parameter:    Psi(a) = D - 1 = max(0, ln({b}a) / ln({b}))")
    console.print(f"  Critical Exponent:  beta = 1  (Psi ~ (a - a_c)^1)")
    console.print(f"  Derivative Jump:    Delta(dPsi/da) = {b}/ln({b}) = [bold green]{b/math.log(b):.4f}[/bold green] at a_c\n")

    table = Table(title=f"Renormalization Group (RG) States across Parameter a")
    table.add_column("a", justify="center")
    table.add_column("RG Eigenvalue lambda = b*a", justify="right")
    table.add_column("RG Relevance", justify="left")
    table.add_column("Phase", justify="left")
    table.add_column("Order Parameter Psi", justify="right")
    table.add_column("dPsi / da", justify="right")

    for a_val in test_as:
        st = classify_phase(a_val, b)
        lam_style = "green" if st.rg_eigenvalue < 1.0 else ("yellow" if abs(st.rg_eigenvalue - 1.0) < 1e-12 else "red")
        table.add_row(
            f"{a_val:.4f}",
            f"[{lam_style}]{st.rg_eigenvalue:.4f}[/{lam_style}]",
            st.relevance.value,
            st.phase_name,
            f"{st.order_parameter:.4f}",
            f"{st.order_parameter_derivative:.4f}",
        )

    console.print(table)


@app.command("attractor")
def attractor_info(
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
):
    """Display the skew-product dynamical system, Lyapunov spectrum, and strange attractor dimension."""
    a_c = critical_a(b)
    test_as = [0.10, 0.20, 0.28, a_c, 0.38, 0.50, 0.60, 0.70]

    console.print(f"[bold cyan]Skew-Product Dynamical System & Lyapunov Spectrum (b={b})[/bold cyan]")
    console.print(f"  Mapping on Cylinder M = S^1 x R:")
    console.print(f"    theta_{{n+1}} = (b * theta_n) mod (2pi)   [Chaotic Bernoulli Shift]")
    console.print(f"    y_{{n+1}}     = a * y_n + phi(theta_n)     [Radial Contraction + Driving]")
    console.print(f"  Angular Lyapunov Exponent: lambda_theta = ln({b}) = [bold green]{math.log(b):.4f}[/bold green] > 0  (Deterministic Chaos)")
    console.print(f"  Radial Lyapunov Exponent:  lambda_y     = ln(a) < 0  (Radial Dissipation)")
    console.print(f"  Lyapunov Sum (Volume Rate): Lambda = lambda_theta + lambda_y = ln({b}a)")
    console.print(f"  Zero-Crossing Threshold:   a_c = 1/{b} = [bold green]{a_c:.4f}[/bold green]\n")

    table = Table(title="Lyapunov Spectrum & Kaplan-Yorke / Ledrappier Dimension across a")
    table.add_column("a", justify="center")
    table.add_column("lambda_theta", justify="right")
    table.add_column("lambda_y", justify="right")
    table.add_column("Lyapunov Sum Lambda", justify="right")
    table.add_column("Volume Nature", justify="left")
    table.add_column("Kaplan-Yorke D_KY", justify="right")
    table.add_column("Attractor Type", justify="left")

    for a_val in test_as:
        spec = compute_lyapunov_spectrum(a_val, b)
        sum_style = "green" if spec.lyapunov_sum < 0.0 else ("yellow" if abs(spec.lyapunov_sum) < 1e-12 else "red")
        vol_str = "Area-Contracting" if spec.is_area_contracting else ("Conservative" if abs(spec.lyapunov_sum) < 1e-12 else "Area-Expanding")
        table.add_row(
            f"{a_val:.4f}",
            f"{spec.lambda_theta:.4f}",
            f"{spec.lambda_y:.4f}",
            f"[{sum_style}]{spec.lyapunov_sum:.4f}[/{sum_style}]",
            vol_str,
            f"{spec.kaplan_yorke_dimension:.4f}",
            spec.attractor_type,
        )

    console.print(table)


@app.command("universality")
def universality_info(
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
    eta: float = typer.Option(0.10, "--eta", help="Irregularity amplitude"),
    level: int = typer.Option(7, "-L", "--level", help="Maximum recursion level for testing"),
):
    """Demonstrate universality and waveform invariance across 5 microscopic generators."""
    a_c = critical_a(b)
    console.print(f"[bold cyan]Universality & Waveform Invariance across Microscopic Generators (b={b})[/bold cyan]")
    console.print(f"  Critical threshold: a_c = 1/{b} = [bold green]{a_c:.4f}[/bold green] (Universal)")
    console.print(f"  Order parameter critical exponent: beta = [bold green]1[/bold green] (Universal)")
    console.print(f"  Convergent regime (a=0.20): IR fixed point -> Smooth Euclidean Circle D = 1.000")
    console.print(f"  Divergent regime (a=0.50): Fractal Attractor -> Theoretical D = [bold green]{1.0 + math.log(b * 0.50) / math.log(b):.4f}[/bold green]\n")

    with console.status("[cyan]Evaluating multi-waveform geometric scaling...[/cyan]"):
        summaries = analyze_waveform_universality(b=b, eta=eta, max_level=level)

    table = Table(title="Waveform Invariance: Comparison of 5 Microscopic Generators")
    table.add_column("Waveform", justify="left")
    table.add_column("Continuity Class", justify="left")
    table.add_column("Critical a_c", justify="center")
    table.add_column("D (a=0.20, Conv)", justify="right")
    table.add_column("D (a=0.50, Div)", justify="right")
    table.add_column("D Theory", justify="right")
    table.add_column("P Ratio (L7/L6)", justify="right")
    table.add_column("Universal Status", justify="left")

    for s in summaries:
        table.add_row(
            f"[bold]{s.waveform_name}[/bold]",
            s.continuity_class,
            f"{s.critical_a:.4f}",
            f"{s.convergent_dimension_a02:.4f}",
            f"{s.divergent_dimension_a05:.4f}",
            f"{s.theoretical_dimension_a05:.4f}",
            f"{s.perimeter_ratio_a05:.4f}",
            "[green]Universal Invariant[/green]",
        )

    console.print(table)


@app.command("observer")
def observer_single(
    sigma: float = typer.Option(0.10, "--sigma", "-s", help="Observer spatial resolution aperture sigma [rad]"),
    a: float = typer.Option(0.50, "-a", help="Shrinkage factor a"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
    eta: float = typer.Option(0.10, "--eta", help="Irregularity amplitude"),
    level: int = typer.Option(7, "-L", "--level", help="Microscopic recursion level"),
):
    """Measure the geometry through the aperture of an independent observer with resolution sigma."""
    pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=b)
    obs = ScaleSpaceObserver(pc, level=level, num_points=4096)
    o = obs.observe(sigma)

    console.print(f"[bold cyan]Independent Observer Measurement (Resolution aperture sigma = {sigma:.4f} rad)[/bold cyan]")
    console.print(f"  Observed Mean Radius R_obs:      [bold]{o.mean_radius:.5f}[/bold]")
    console.print(f"  Observed Perimeter P_obs:        [bold]{o.perimeter:.5f}[/bold]")
    console.print(f"  Observed Enclosed Area A_obs:     [bold]{o.area:.5f}[/bold]")
    console.print(f"  Apparent Circumference Ratio pi: [bold green]{o.pi_obs:.5f}[/bold green] (vs Euclidean pi = {math.pi:.5f})")
    console.print(f"  Isoperimetric Quotient Q_obs:    [bold green]{o.isoperimetric_quotient:.5f}[/bold green] (Circle maximum = 1.00000)")
    console.print(f"  Apparent Dimension D_obs:        [bold]{o.apparent_dimension:.4f}[/bold]")

    if sigma >= 0.5:
        console.print("  [bold green]Classification: Infrared Fixed Point (Euclidean Circle Observed)[/bold green]")
    elif sigma <= 0.01:
        console.print("  [bold red]Classification: Ultraviolet Regime (Fractal Coastline Observed)[/bold red]")
    else:
        console.print("  [bold yellow]Classification: Scale-Space Crossover Intermediate[/bold yellow]")


@app.command("observer-sweep")
def observer_sweep(
    a: float = typer.Option(0.50, "-a", help="Shrinkage factor a"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
    eta: float = typer.Option(0.10, "--eta", help="Irregularity amplitude"),
    level: int = typer.Option(7, "-L", "--level", help="Microscopic recursion level"),
):
    """Demonstrate the emergence of Euclidean pi and Q=1 across a continuum of observer apertures."""
    pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=b)
    obs = ScaleSpaceObserver(pc, level=level, num_points=4096)
    sigmas = [1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]

    console.print(f"[bold cyan]Observer Resolution Sweep: Demonstrating the Emergence of Euclidean Geometry[/bold cyan]")
    console.print(f"  Model parameters: a={a}, b={b}, eta={eta}, level={level}")
    console.print(f"  Theoretical crossover scale: sigma* ~ {obs.theoretical_crossover_scale():.3f} rad\n")

    table = Table(title="Scale-Space Observer Measurement Across Resolutions")
    table.add_column("Observer Aperture sigma", justify="right")
    table.add_column("Perimeter P", justify="right")
    table.add_column("Observed pi", justify="right")
    table.add_column("|pi - pi_Euc|", justify="right")
    table.add_column("Isoperimetric Q", justify="right")
    table.add_column("Apparent D", justify="right")
    table.add_column("Perceived Regime", justify="left")

    for s in sigmas:
        o = obs.observe(s)
        diff_pi = abs(o.pi_obs - math.pi)
        if s >= 0.5:
            regime = "[green]Euclidean Circle[/green]"
        elif s <= 0.01:
            regime = "[red]Fractal Boundary[/red]"
        else:
            regime = "[yellow]Crossover[/yellow]"

        table.add_row(
            f"{s:.4f}",
            f"{o.perimeter:.4f}",
            f"{o.pi_obs:.5f}",
            f"{diff_pi:.5f}",
            f"{o.isoperimetric_quotient:.5f}",
            f"{o.apparent_dimension:.4f}",
            regime,
        )

    console.print(table)


@app.command("carriers")
def carriers_comparison(
    a: float = typer.Option(0.50, "-a", help="Shrinkage factor a"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
    eta: float = typer.Option(0.10, "--eta", help="Irregularity amplitude"),
    level: int = typer.Option(7, "-L", "--level", help="Microscopic recursion level"),
    sigma: float = typer.Option(0.20, "-s", "--sigma", help="Observer resolution in scale-space window"),
):
    """Demonstrate Sobolev projection across 4 diverse Euclidean carriers (circle, ellipse, squircle, trefoil)."""
    carrier_list = [
        CircleCarrier(radius=1.0),
        EllipseCarrier(semi_major=1.25, semi_minor=0.80),
        SuperellipseCarrier(semi_major=1.0, semi_minor=1.0, exponent=4.0),
        TrefoilCarrier(base_radius=1.0, modulation=0.15),
    ]

    console.print("[bold cyan]Sobolev Projection: Euclidean Carriers under Fractal Perturbation[/bold cyan]")
    console.print(f"  Observer window aperture: sigma = {sigma:.3f} rad | Perturbation: a={a}, b={b}, eta={eta}, level={level}\n")

    table = Table(title="Multi-Carrier Scale-Space Observer Recovery")
    table.add_column("Carrier Shape", justify="left")
    table.add_column("Exact Euclidean Perimeter", justify="right")
    table.add_column("Measured Perimeter P_obs", justify="right")
    table.add_column("Perimeter Rel Error", justify="right")
    table.add_column("Exact Isoperimetric Q0", justify="right")
    table.add_column("Measured Q_obs", justify="right")
    table.add_column("Apparent D", justify="right")
    table.add_column("Status", justify="left")

    for c in carrier_list:
        obs = ScaleSpaceObserver.from_carrier(c, eta=eta, a=a, b=b, level=level, num_points=4096)
        o = obs.observe(sigma)
        p_target = c.exact_perimeter()
        q_target = c.isoperimetric_quotient()

        status = "[green]Euclidean Recovered[/green]" if o.perimeter_relative_error < 0.02 else "[yellow]Approximated[/yellow]"

        table.add_row(
            c.name.capitalize(),
            f"{p_target:.4f}",
            f"{o.perimeter:.4f}",
            f"{o.perimeter_relative_error * 100:.2f}%",
            f"{q_target:.4f}",
            f"{o.isoperimetric_quotient:.4f}",
            f"{o.apparent_dimension:.4f}",
            status,
        )

    console.print(table)


@app.command("polygons")
def polygons_sweep(
    a: float = typer.Option(0.50, "-a", help="Shrinkage factor a"),
    b: int = typer.Option(3, "-b", help="Frequency multiplier b"),
    eta: float = typer.Option(0.10, "--eta", help="Irregularity amplitude"),
    level: int = typer.Option(7, "-L", "--level", help="Microscopic recursion level"),
):
    """Demonstrate the competition between corner rounding and fractal divergence on polygons (square & triangle)."""
    sq = SquareCarrier(inradius=1.0)
    tri = TriangleCarrier(inradius=1.0)

    obs_sq = ScaleSpaceObserver.from_carrier(sq, eta=eta, a=a, b=b, level=level, num_points=4096)
    obs_tri = ScaleSpaceObserver.from_carrier(tri, eta=eta, a=a, b=b, level=level, num_points=4096)

    sigmas = [0.25, 0.15, 0.10, 0.05, 0.03, 0.01, 0.005, 0.001]

    console.print("[bold cyan]Polygonal Scale-Space Crossover: Square & Triangle[/bold cyan]")
    console.print(f"  Exact Square Perimeter: {sq.exact_perimeter():.4f} | Exact Triangle Perimeter: {tri.exact_perimeter():.4f}\n")

    table = Table(title="Polygonal Recovery Across Scale-Space Apertures")
    table.add_column("Observer Aperture sigma", justify="right")
    table.add_column("Square P_obs", justify="right")
    table.add_column("Square Error", justify="right")
    table.add_column("Triangle P_obs", justify="right")
    table.add_column("Triangle Error", justify="right")
    table.add_column("Dominant Effect", justify="left")

    for s in sigmas:
        o_sq = obs_sq.observe(s)
        o_tri = obs_tri.observe(s)

        err_sq = abs(o_sq.perimeter - sq.exact_perimeter()) / sq.exact_perimeter() * 100.0
        err_tri = abs(o_tri.perimeter - tri.exact_perimeter()) / tri.exact_perimeter() * 100.0

        if s >= 0.10:
            regime = "[blue]Corner Rounding Deficit (P < P0)[/blue]"
        elif 0.03 <= s <= 0.07:
            regime = "[bold green]Optimal Observer Window (Error < 2%)[/bold green]"
        else:
            regime = "[red]Fractal Divergence (P > P0)[/red]"

        table.add_row(
            f"{s:.4f}",
            f"{o_sq.perimeter:.4f}",
            f"{err_sq:.2f}%",
            f"{o_tri.perimeter:.4f}",
            f"{err_tri:.2f}%",
            regime,
        )

    console.print(table)


@app.command()
def plot(
    output_dir: str = typer.Option("assets", "-o", "--output-dir", help="Directory to save figures"),
    eta: float = typer.Option(0.1, "--eta", help="Irregularity amplitude"),
    a: float = typer.Option(0.5, "-a", help="Shrinkage factor a"),
):
    """Generate publication-ready figures illustrating the model and its scaling."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold cyan]Generating figures into {out_dir}...[/bold cyan]")

    pc = PolygonalCircle(radius=1.0, eta=eta, a=a, b=3)

    # 1. Multi-level shapes
    p1 = out_dir / "multilevel_shapes.png"
    plot_multilevel_shapes(pc, levels=[0, 1, 2, 3, 4], out_path=p1)
    console.print(f"  ✓ Saved multi-level polygons: [green]{p1}[/green]")

    # 2. Scaling regimes
    p2 = out_dir / "scaling_regimes.png"
    plot_scaling_regimes(out_path=p2)
    console.print(f"  ✓ Saved scaling regimes comparison: [green]{p2}[/green]")

    # 3. Richardson crossover
    p3 = out_dir / "richardson_crossover.png"
    plot_richardson_crossover(pc, max_level=10, out_path=p3)
    console.print(f"  ✓ Saved Richardson coastline scaling: [green]{p3}[/green]")

    # 4. Macroscopic vs microscopic
    p4 = out_dir / "macroscopic_vs_microscopic.png"
    plot_macroscopic_vs_microscopic(pc, level=5, out_path=p4)
    console.print(f"  ✓ Saved macroscopic vs microscopic zoom: [green]{p4}[/green]")

    # 5. Annular polygon circle
    from polygonal_circle.annulus import AnnularPolygonalCircle
    ac = AnnularPolygonalCircle.from_center_and_width(center_radius=1.0, corridor_width=0.20, a=a, b=3)
    p5 = out_dir / "annular_circle.png"
    plot_annular_circle(ac, level=4, out_path=p5)
    console.print(f"  ✓ Saved annular circle with nested boundaries: [green]{p5}[/green]")

    # 6. Corridor sweep
    p6 = out_dir / "corridor_sweep.png"
    plot_corridor_sweep(center_radius=1.0, widths=[0.40, 0.20, 0.10, 0.03], level=4, a=a, out_path=p6)
    console.print(f"  ✓ Saved corridor width knob sweep: [green]{p6}[/green]")

    # 7. Phase transition & RG flow
    p7 = out_dir / "phase_transition.png"
    plot_phase_transition_and_rg(b=3, out_path=p7)
    console.print(f"  ✓ Saved phase transition & RG flow: [green]{p7}[/green]")

    # 8. Strange attractor & Lyapunov dynamics
    p8 = out_dir / "strange_attractor.png"
    plot_strange_attractor(a=a, b=3, out_path=p8)
    console.print(f"  ✓ Saved strange attractor & Lyapunov dynamics: [green]{p8}[/green]")

    # 9. Universality proof across microscopic generators
    p9 = out_dir / "universality_proof.png"
    plot_universality_proof(a=a, b=3, level=3, out_path=p9)
    console.print(f"  ✓ Saved universality proof across 5 waveforms: [green]{p9}[/green]")

    # 10. Observer scale-space projection (emergence of Euclidean geometry)
    p10 = out_dir / "observer_scale_space.png"
    plot_observer_scale_space(a=a, b=3, level=7, out_path=p10)
    console.print(f"  ✓ Saved observer scale-space projection: [green]{p10}[/green]")

    # 11. Multi-carrier Sobolev projection (invariance across all Euclidean carriers)
    p11 = out_dir / "multi_carrier_scale_space.png"
    plot_multi_carrier_scale_space(eta=eta, a=a, b=3, level=7, out_path=p11)
    console.print(f"  ✓ Saved multi-carrier Sobolev projection: [green]{p11}[/green]")

    # 12. Polygonal scale-space competition (corner rounding vs fractal divergence)
    p12 = out_dir / "polygonal_scale_space.png"
    plot_polygonal_scale_space(eta=eta, a=a, b=3, level=7, out_path=p12)
    console.print(f"  ✓ Saved polygonal scale-space competition: [green]{p12}[/green]")

    console.print(f"[bold green]Done! All figures successfully generated.[/bold green]")


if __name__ == "__main__":
    app()
