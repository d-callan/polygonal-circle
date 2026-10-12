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
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.plotting import (
    plot_annular_circle,
    plot_corridor_sweep,
    plot_macroscopic_vs_microscopic,
    plot_multilevel_shapes,
    plot_phase_transition_and_rg,
    plot_richardson_crossover,
    plot_scaling_regimes,
)
from polygonal_circle.renormalization import classify_phase, critical_a

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

    console.print(f"[bold green]Done! All figures successfully generated.[/bold green]")


if __name__ == "__main__":
    app()
