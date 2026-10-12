"""Command-line interface for the polygonal-circle project."""

from __future__ import annotations

from pathlib import Path
import rich
from rich.console import Console
from rich.table import Table
import typer

from polygonal_circle.analysis import PolygonalCircleAnalysis
from polygonal_circle.model import PolygonalCircle
from polygonal_circle.plotting import (
    plot_macroscopic_vs_microscopic,
    plot_multilevel_shapes,
    plot_richardson_crossover,
    plot_scaling_regimes,
)

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

    console.print(f"[bold green]Done! All figures successfully generated.[/bold green]")


if __name__ == "__main__":
    app()
