"""
cli.py — Command-line interface for the ideoplot package.

Usage examples
--------------
# Genome-wide ideoplot with gene markers
$ python -m ideoplot genome --genes BRCA1 BRCA2 TP53 --save out.png

# Local region view
$ python -m ideoplot region --gene MSH2 --margin 500000 --exons --save msh2.png

# List built-in themes
$ python -m ideoplot themes
"""
from __future__ import annotations

import argparse
import sys


def _cmd_genome(args: argparse.Namespace) -> None:
    from .core import IdeoplotPlotter

    chroms = None
    if args.chromosomes:
        chroms = [c if c.startswith("chr") else f"chr{c}" for c in args.chromosomes]

    plotter = IdeoplotPlotter(
        assembly=args.assembly,
        theme=args.theme,
        chromosomes=chroms,
        figsize=(args.width, args.height),
        chrom_width=args.chrom_width,
        font_size=args.font_size,
        layout_rows=args.rows,
        karyotype=args.karyotype,
    )
    plotter.plot(
        genes=args.genes or None,
        scale_bar_mb=args.scale_bar,
        title=args.title or "",
        save=args.save or None,
        dpi=args.dpi,
        show=not args.no_show,
    )


def _cmd_region(args: argparse.Namespace) -> None:
    from .core import GenomeViewer

    viewer = GenomeViewer(
        assembly=args.assembly,
        theme=args.theme,
        figsize_width=args.width,
        font_size=args.font_size,
        include_pseudogenes=not args.no_pseudo,
    )
    viewer.plot_region(
        target_gene=args.gene,
        margin_bp=args.margin,
        show_exons=args.exons,
        highlight_genes=args.highlight or None,
        title=args.title or None,
        save=args.save or None,
        dpi=args.dpi,
        show=not args.no_show,
    )


def _cmd_themes(_args: argparse.Namespace) -> None:
    from .config import THEMES

    print("Built-in colour themes:")
    for name in THEMES:
        print(f"  {name}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        prog="ideoplot",
        description="Human genome visualisation toolkit",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── genome ──────────────────────────────────────────────────────────
    p_genome = sub.add_parser("genome", help="Genome-wide ideoplot")
    p_genome.add_argument("--genes", nargs="*", help="Gene symbols to mark")
    p_genome.add_argument("--chromosomes", nargs="*", help="Specific chromosomes to plot (e.g., 1 2 X)")
    p_genome.add_argument("--assembly", default="hg38")
    p_genome.add_argument("--theme", default="classic")
    p_genome.add_argument("--width", type=float, default=18)
    p_genome.add_argument("--height", type=float, default=9)
    p_genome.add_argument("--chrom-width", dest="chrom_width", type=float, default=0.3)
    p_genome.add_argument("--font-size", dest="font_size", type=int, default=10)
    p_genome.add_argument("--scale-bar", dest="scale_bar", type=float, default=50)
    p_genome.add_argument("--title", default="")
    p_genome.add_argument("--save", default="")
    p_genome.add_argument("--dpi", type=int, default=300)
    p_genome.add_argument("--rows", type=int, default=1, help="Number of layout rows")
    p_genome.add_argument("--karyotype", default="XY", help="Karyotype to display (build default chromosomes): XY or XX")
    p_genome.add_argument("--no-show", action="store_true")
    p_genome.set_defaults(func=_cmd_genome)

    # ── region ──────────────────────────────────────────────────────────
    p_region = sub.add_parser("region", help="Local region view around a gene")
    p_region.add_argument("--gene", required=True, help="Target gene symbol")
    p_region.add_argument("--margin", type=int, default=500_000)
    p_region.add_argument("--exons", action="store_true")
    p_region.add_argument("--highlight", nargs="*", help="Extra genes to highlight")
    p_region.add_argument("--assembly", default="hg38")
    p_region.add_argument("--theme", default="classic")
    p_region.add_argument("--width", type=float, default=16)
    p_region.add_argument("--font-size", dest="font_size", type=int, default=10)
    p_region.add_argument("--no-pseudo", action="store_true")
    p_region.add_argument("--title", default="")
    p_region.add_argument("--save", default="")
    p_region.add_argument("--dpi", type=int, default=300)
    p_region.add_argument("--no-show", action="store_true")
    p_region.set_defaults(func=_cmd_region)

    # ── themes ──────────────────────────────────────────────────────────
    p_themes = sub.add_parser("themes", help="List built-in colour themes")
    p_themes.set_defaults(func=_cmd_themes)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
