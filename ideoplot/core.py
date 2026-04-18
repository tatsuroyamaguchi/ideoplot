"""
core.py — High-level plotting classes.

IdeoplotPlotter  – genome-wide view (all chromosomes) with gene markers
GenomeViewer     – local region view with gene bodies and optional exon structure
"""
from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

from .config import CHROMS_ALL, CHROMS_AUTOSOMES, ColorTheme, THEMES
from .fetch import CytobandFetcher, EnsemblClient, fetch_gene_dataframe
from .utils import (
    add_scale_bar,
    assign_tracks,
    compute_aspect_correction,
    draw_chromosome_body,
    draw_arrow_gene,
    draw_exon_gene,
)


# ---------------------------------------------------------------------------
# IdeoplotPlotter — genome-wide ideoplot
# ---------------------------------------------------------------------------

class IdeoplotPlotter:
    """Plot a genome-wide ideoplot (all chromosomes) with optional gene markers.

    Parameters
    ----------
    assembly : str
        Genome assembly: ``"hg38"`` (default) or ``"hg19"``.
    theme : str or ColorTheme
        Visual theme. Built-in names: ``"classic"``, ``"dark"``,
        ``"colorblind"``, ``"pastel"``.
    chromosomes : list of str, optional
        Which chromosomes to show.  Default: chr1–22, X, Y.
    figsize : tuple of float, optional
        Figure size ``(width_inches, height_inches)``.  Default ``(18, 9)``.
    chrom_width : float
        Width of each chromosome bar in x-data-units.  Default ``0.3``.
    font_size : int
        Label font size.  Default ``10``.
    cache : bool
        Cache network data to ``~/.cache/ideoplot/``.  Default ``True``.
    """

    def __init__(
        self,
        assembly: str = "hg38",
        theme: Union[str, ColorTheme] = "classic",
        chromosomes: Optional[List[str]] = None,
        figsize: Tuple[float, float] = (18, 9),
        chrom_width: float = 0.3,
        font_size: int = 10,
        cache: bool = True,
        layout_rows: int = 1,
        karyotype: str = "XY",
    ) -> None:
        self.assembly = assembly
        self.theme: ColorTheme = THEMES[theme] if isinstance(theme, str) else theme
        if chromosomes is None:
            base_chroms = list(CHROMS_AUTOSOMES) + (["chrX", "chrX"] if karyotype.upper() == "XX" else ["chrX", "chrY"])
        else:
            base_chroms = list(chromosomes)

        self.chromosomes = base_chroms
        self.layout_rows = layout_rows
        self.figsize = (figsize[0], figsize[1] * layout_rows) if layout_rows > 1 else figsize
        self.chrom_width = chrom_width
        self.font_size = font_size
        self._cyto = CytobandFetcher(assembly=assembly, cache=cache)
        self._ensembl = EnsemblClient(cache=cache)

    # ------------------------------------------------------------------
    def plot(
        self,
        genes: Optional[List[str]] = None,
        gene_df: Optional[pd.DataFrame] = None,
        marker_color: Optional[str] = None,
        marker_size: float = 5,
        label_offset: float = 0.25,
        scale_bar_mb: float = 50,
        title: str = "",
        save: Optional[str] = None,
        dpi: int = 300,
        show: bool = True,
        ax: Optional[plt.Axes] = None,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Draw the genome-wide ideoplot.

        Parameters
        ----------
        genes : list of str, optional
            Gene symbols to mark on the ideoplot.  The positions are fetched
            from Ensembl automatically.
        gene_df : pd.DataFrame, optional
            Pre-computed gene positions (columns: ``gene``, ``chrom``,
            ``pos`` in Mb).  Overrides *genes*.
        marker_color : str, optional
            Dot colour for gene markers.  Defaults to ``theme.target_gene``.
        marker_size : float
            Marker dot size.
        label_offset : float
            Horizontal label offset from the marker dot (x-units).
        scale_bar_mb : float
            Scale bar length in Mb.  Pass ``0`` to hide.
        title : str
            Figure title.
        save : str, optional
            File path to save the figure.  Format is inferred from the
            extension (e.g. ``.png``, ``.svg``, ``.pdf``).
        dpi : int
            DPI for raster export.
        show : bool
            Call ``plt.show()`` when done.
        ax : matplotlib Axes, optional
            Draw into an existing Axes (useful for subplots).

        Returns
        -------
        (fig, ax)
        """
        color_map = self.theme.cytoband
        marker_color = marker_color or self.theme.target_gene

        # ---- Data --------------------------------------------------------
        cyto_df = self._cyto.load()
        chrom_sizes = self._cyto.chromosome_sizes()
        chroms = [c for c in self.chromosomes if c in cyto_df["chrom"].unique()]

        # Gene positions
        if gene_df is None and genes:
            gene_df = fetch_gene_dataframe(genes, client=self._ensembl, assembly=self.assembly)
        if gene_df is None:
            gene_df = pd.DataFrame(columns=["gene", "chrom", "pos"])

        # ---- Layout -------------------------------------------------------
        max_bp_mb = max(
            (chrom_sizes[c][1] / 1_000_000 for c in chroms if c in chrom_sizes),
            default=250,
        )

        import math
        cols = math.ceil(len(chroms) / self.layout_rows)
        row_height = max_bp_mb * 1.25

        # ---- Figure -------------------------------------------------------
        if ax is None:
            standard_cols = math.ceil(24 / self.layout_rows)
            # Adjust width proportionally to maintain horizontal aspect ratio.
            # Using +1.0 to account for the -0.8 layout offset.
            target_width = self.figsize[0] * ((cols + 1.0) / (standard_cols + 1.0))
            fig, ax = plt.subplots(figsize=(target_width, self.figsize[1]))
        else:
            fig = ax.get_figure()

        ax.set_xlim(-0.8, cols - 0.2)
        ax.set_ylim((self.layout_rows - 1) * row_height + max_bp_mb * 1.02, -max_bp_mb * 0.05)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        if title:
            ax.set_title(title, fontsize=self.font_size + 4, pad=12)

        plt.tight_layout()
        fig.canvas.draw()
        aspect = compute_aspect_correction(ax, fig)

        # ---- Draw chromosomes --------------------------------------------
        offsets_list: List[float] = [0.0] * len(chroms)
        for i, chrom in enumerate(chroms):
            cd = cyto_df[cyto_df["chrom"] == chrom].copy()
            cd_mb = cd.assign(
                chromStart=cd["chromStart"] / 1_000_000,
                chromEnd=cd["chromEnd"] / 1_000_000,
            )
            if cd_mb.empty:
                continue
            
            xi = i % cols
            ri = i // cols
            base_y = ri * row_height

            y0 = cd_mb["chromStart"].min()
            y1 = cd_mb["chromEnd"].max()
            
            offset = max_bp_mb - (y1 - y0)
            offsets_list[i] = offset + base_y
            cd_mb["chromStart"] += (offset + base_y)
            cd_mb["chromEnd"] += (offset + base_y)
            draw_chromosome_body(
                ax, xi, y0 + offset + base_y, y1 + offset + base_y, cd_mb,
                color_map, width=self.chrom_width, aspect_correction=aspect,
            )
            
            # Label placed below the chromosome
            bottom_y = base_y + max_bp_mb
            ax.text(xi, bottom_y + max_bp_mb * 0.05, chrom.replace("chr", ""),
                    ha="center", va="top", fontsize=self.font_size + 2)

        # ---- Gene markers ------------------------------------------------
        for _, row in gene_df.iterrows():
            chrom = row["chrom"]
            
            indices = [idx for idx, c in enumerate(chroms) if c == chrom]
            for i in indices:
                xi = i % cols
                yi = row["pos"] + offsets_list[i]
                ax.plot(
                    xi + self.chrom_width / 2 + 0.05, yi,
                    "o", color=marker_color, markersize=marker_size, zorder=10,
                )
                ax.text(
                    xi + self.chrom_width / 2 + label_offset, yi,
                    row["gene"],
                    fontsize=self.font_size, va="center", fontstyle="italic", zorder=11,
                )

        # ---- Scale bar ---------------------------------------------------
        if scale_bar_mb > 0:
            x_bar = cols - 0.3
            y_bar = max_bp_mb * 0.1
            add_scale_bar(ax, x_bar, y_bar, scale_bar_mb,
                          color=self.theme.scalebar,
                          fontsize=self.font_size)

        # ---- Save / show -------------------------------------------------
        if save:
            fmt = save.rsplit(".", 1)[-1] if "." in save else "png"
            fig.savefig(save, dpi=dpi, bbox_inches="tight", format=fmt)
            print(f"✅ Saved → {save}")

        if show:
            plt.show()

        return fig, ax

    # ------------------------------------------------------------------
    def add_gene_highlights(
        self,
        ax: plt.Axes,
        chromosomes: List[str],
        gene_df: pd.DataFrame,
        offsets_list: List[float],
        colors: Optional[Dict[str, str]] = None,
    ) -> None:
        """Add colour-coded gene markers to an existing axes.

        Useful for layering multiple gene groups with different colours
        after calling :meth:`plot`.
        """
        import math
        cols = math.ceil(len(chromosomes) / self.layout_rows)
        for _, row in gene_df.iterrows():
            chrom = row["chrom"]
                
            indices = [idx for idx, c in enumerate(chromosomes) if c == chrom]
            for i in indices:
                color = (colors or {}).get(row["gene"], self.theme.target_gene)
                
                xi = i % cols
                yi = row["pos"] + offsets_list[i]
                ax.plot(xi + self.chrom_width / 2 + 0.05, yi, "o",
                        color=color, markersize=5, zorder=10)
                ax.text(xi + self.chrom_width / 2 + 0.25, yi, row["gene"],
                        fontsize=self.font_size, va="center",
                        fontstyle="italic", color=color, zorder=11)


# ---------------------------------------------------------------------------
# GenomeViewer — local region view
# ---------------------------------------------------------------------------

class GenomeViewer:
    """Visualise a genomic region around a target gene.

    The view shows:
    - A cytoband strip at the top for chromsomal context.
    - All protein-coding genes and pseudogenes in the region, drawn
      with arrow bodies (showing strand direction) or exon/intron
      structure.

    Parameters
    ----------
    assembly : str
        Genome assembly: ``"hg38"`` (default) or ``"hg19"``.
    theme : str or ColorTheme
        Visual theme.
    figsize_width : float
        Figure width in inches.  Height is auto-scaled.  Default ``16``.
    gene_height : float
        Height of gene bodies in y-units.  Default ``0.25``.
    font_size : int
        Label font size.  Default ``10``.
    include_pseudogenes : bool
        Whether to draw pseudogenes.  Default ``True``.
    cache : bool
        Cache network data.  Default ``True``.
    """

    def __init__(
        self,
        assembly: str = "hg38",
        theme: Union[str, ColorTheme] = "classic",
        figsize_width: float = 16,
        gene_height: float = 0.25,
        font_size: int = 10,
        include_pseudogenes: bool = True,
        cache: bool = True,
    ) -> None:
        self.assembly = assembly
        self.theme: ColorTheme = THEMES[theme] if isinstance(theme, str) else theme
        self.figsize_width = figsize_width
        self.gene_height = gene_height
        self.font_size = font_size
        self.include_pseudogenes = include_pseudogenes
        self._cyto = CytobandFetcher(assembly=assembly, cache=cache)
        self._ensembl = EnsemblClient(cache=cache)

    # ------------------------------------------------------------------
    def plot_region(
        self,
        target_gene: str,
        margin_bp: int = 500_000,
        show_exons: bool = True,
        highlight_genes: Optional[List[str]] = None,
        title: Optional[str] = None,
        save: Optional[str] = None,
        dpi: int = 300,
        show: bool = True,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Draw a local genomic region centred on *target_gene*.

        Parameters
        ----------
        target_gene : str
            Gene symbol.  The region ``[gene_start - margin_bp,
            gene_end + margin_bp]`` is displayed.
        margin_bp : int
            Flanking margin on each side (default 500 kb).
        show_exons : bool
            When ``True``, fetch exon coordinates and draw exon–intron
            structure.  When ``False``, use filled arrow polygons.
        highlight_genes : list of str, optional
            Additional gene names to highlight in the target colour.
        title : str, optional
            Figure title.  Defaults to auto-generated.
        save : str, optional
            Output file path.
        dpi : int
            DPI for raster output.
        show : bool
            Whether to call ``plt.show()``.

        Returns
        -------
        (fig, ax)
        """
        highlight_set = set(highlight_genes or []) | {target_gene}

        # ---- Fetch target gene position ---------------------------------
        print(f"🔍 Fetching position for {target_gene}…")
        try:
            info = self._ensembl.lookup_gene(target_gene)
            chrom = str(info["seq_region_name"])
            tgt_start_bp = info["start"]
            tgt_end_bp = info["end"]
        except Exception as exc:
            from .fetch import UcscClient
            ucsc = UcscClient(assembly=self.assembly, cache=self._ensembl.cache, cache_dir=self._ensembl.cache_dir)
            res = ucsc.search_position(target_gene)
            if not res:
                raise ValueError(f"Target '{target_gene}' not found in Ensembl or UCSC search.") from exc
            chrom = res["chrom"].replace("chr", "")
            tgt_start_bp = res["start"]
            tgt_end_bp = res["end"]
            
        region_start_bp = max(1, tgt_start_bp - margin_bp)
        region_end_bp = tgt_end_bp + margin_bp

        print(f"  📍 {target_gene}: chr{chrom}:{tgt_start_bp:,}–{tgt_end_bp:,}")
        print(f"  📏 Region: chr{chrom}:{region_start_bp:,}–{region_end_bp:,}")

        # ---- Genes in region ---------------------------------------------
        print("🧬 Fetching genes in region…")
        genes_raw = self._ensembl.genes_in_region(chrom, region_start_bp, region_end_bp)
        
        found_in_ensembl = any(g.get("external_name") == target_gene for g in genes_raw)
        if not found_in_ensembl:
            genes_raw.append({
                "id": target_gene,
                "external_name": target_gene,
                "biotype": "marker",
                "start": tgt_start_bp,
                "end": tgt_end_bp,
                "strand": 1,
            })
            
        genes_raw = sorted(genes_raw, key=lambda g: g["start"])

        def _keep(g: dict) -> bool:
            bt = g.get("biotype", "")
            if bt == "marker":
                return True
            if bt == "protein_coding":
                return True
            if self.include_pseudogenes and "pseudogene" in bt:
                return True
            return False

        genes = [g for g in genes_raw if _keep(g)]

        # ---- Exon data (optional) ----------------------------------------
        gene_details: Dict[str, dict] = {}
        if show_exons:
            gene_ids = [g["id"] for g in genes]
            print(f"🔬 Fetching exon structure for {len(gene_ids)} genes…")
            gene_details = self._ensembl.lookup_ids_with_exons(gene_ids)

        # ---- Cytoband strip ----------------------------------------------
        bands_df = self._cyto.for_region(chrom, region_start_bp, region_end_bp)
        color_map = self.theme.cytoband

        # ---- Track assignment --------------------------------------------
        region_start_mb = region_start_bp / 1_000_000
        region_end_mb = region_end_bp / 1_000_000
        track_indices = assign_tracks(genes, region_start_mb, region_end_mb)
        num_tracks = max(track_indices, default=-1) + 1

        # ---- Figure layout -----------------------------------------------
        fig_height = max(6, 2 + num_tracks * 0.6)
        fig, ax = plt.subplots(figsize=(self.figsize_width, fig_height))

        cyto_y = num_tracks + 1
        cyto_height = 0.5

        # ---- Draw cytoband -----------------------------------------------
        for _, row in bands_df.iterrows():
            bs = max(region_start_bp, row["chromStart"]) / 1_000_000
            be = min(region_end_bp, row["chromEnd"]) / 1_000_000
            color = color_map.get(row["gieStain"], "#FFFFFF")
            rect = mpatches.Rectangle(
                (bs, cyto_y), be - bs, cyto_height,
                facecolor=color, edgecolor="black", linewidth=0.8,
            )
            ax.add_patch(rect)
            cx = (bs + be) / 2
            ax.text(cx, cyto_y + cyto_height + 0.1, row["name"],
                    ha="center", va="bottom",
                    fontsize=self.font_size - 2, rotation=45)

        # ---- Draw genes --------------------------------------------------
        arrow_head_mb = (region_end_mb - region_start_mb) * 0.005
        min_exon_mb = (region_end_mb - region_start_mb) * 0.001

        for gene, track_idx in zip(genes, track_indices):
            g_start_mb = max(region_start_bp, gene["start"]) / 1_000_000
            g_end_mb = min(region_end_bp, gene["end"]) / 1_000_000
            g_name = gene.get("external_name", "Unknown")
            biotype = gene.get("biotype", "")
            strand = gene.get("strand", 1)
            gene_id = gene.get("id", "")

            is_target = g_name in highlight_set
            is_pseudo = "pseudogene" in biotype

            if is_target:
                color = self.theme.target_gene
                fw = "bold"
                zorder = 10
            elif is_pseudo:
                color = self.theme.pseudogene
                fw = "normal"
                zorder = 4
            else:
                color = self.theme.gene_default
                fw = "normal"
                zorder = 5

            y_pos = num_tracks - track_idx - 1

            if show_exons:
                details = gene_details.get(gene_id, {}) or {}
                transcripts = details.get("Transcript", [])
                canonical = next(
                    (tx for tx in transcripts if tx.get("is_canonical") == 1), None
                )
                if not canonical and transcripts:
                    canonical = transcripts[0]
                exons = canonical.get("Exon", []) if canonical else []

                draw_exon_gene(
                    ax, g_start_mb, g_end_mb, y_pos, strand, exons,
                    region_start_mb, region_end_mb, color,
                    height=self.gene_height,
                    min_exon_width_mb=min_exon_mb,
                    arrow_head_mb=(region_end_mb - region_start_mb) * 0.003,
                    zorder=zorder,
                )
            else:
                draw_arrow_gene(
                    ax, g_start_mb, g_end_mb, y_pos, strand, color,
                    height=self.gene_height,
                    arrow_head_mb=arrow_head_mb,
                    zorder=zorder,
                )

            ax.text(
                g_start_mb, y_pos + self.gene_height + 0.05, g_name,
                ha="left", va="bottom", color=color,
                fontsize=self.font_size, fontweight=fw, fontstyle="italic",
                zorder=zorder,
            )

            if is_target:
                ax.axvline(g_start_mb, color=color, linestyle="--",
                           linewidth=0.8, alpha=0.35, zorder=1)
                ax.axvline(g_end_mb, color=color, linestyle="--",
                           linewidth=0.8, alpha=0.35, zorder=1)

        # ---- Axis styling ------------------------------------------------
        ax.set_xlim(region_start_mb, region_end_mb)
        ax.set_ylim(-0.5, cyto_y + 1.5)
        ax.set_yticks([])
        ax.set_xlabel(f"Chromosome {chrom} (Mb)", fontsize=self.font_size + 2)

        mode = "exon structure" if show_exons else "strand direction"
        ax.set_title(
            title or f"{target_gene} — {mode} ({margin_bp // 1_000} kb margin)",
            fontsize=self.font_size + 4, pad=20,
        )

        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)

        plt.tight_layout()

        if save:
            fmt = save.rsplit(".", 1)[-1] if "." in save else "png"
            fig.savefig(save, dpi=dpi, bbox_inches="tight", format=fmt)
            print(f"✅ Saved → {save}")

        if show:
            plt.show()

        return fig, ax

    # ------------------------------------------------------------------
    def compare_regions(
        self,
        genes: List[str],
        margin_bp: int = 500_000,
        show_exons: bool = False,
        save: Optional[str] = None,
        dpi: int = 300,
        show: bool = True,
    ) -> Tuple[plt.Figure, list]:
        """Plot multiple gene regions side-by-side.

        Parameters
        ----------
        genes : list of str
            Gene symbols; each gets its own subplot row.
        margin_bp : int
            Flanking margin.
        show_exons : bool
            Use exon structure (slower, requires extra API calls).
        save : str, optional
            Output file path.
        dpi : int
            DPI for raster output.
        show : bool
            Whether to call ``plt.show()``.

        Returns
        -------
        (fig, list_of_axes)
        """
        fig, axes = plt.subplots(
            len(genes), 1,
            figsize=(self.figsize_width, 5 * len(genes)),
        )
        if len(genes) == 1:
            axes = [axes]

        for ax, gene in zip(axes, genes):
            # Temporarily redirect show/save
            self.plot_region(
                gene, margin_bp=margin_bp, show_exons=show_exons,
                show=False, save=None,
            )
            plt.close("all")
            self.plot_region(
                gene, margin_bp=margin_bp, show_exons=show_exons,
                show=False, save=None,
            )
            # Re-use the sub-axes directly
            sub_fig, sub_ax = self.plot_region(
                gene, margin_bp=margin_bp, show_exons=show_exons,
                show=False, save=None,
            )
            # Transfer content to the target axes
            for item in sub_ax.get_children():
                try:
                    item_copy = item.__class__.__new__(item.__class__)
                    item_copy.__dict__.update(item.__dict__)
                except Exception:
                    pass
            plt.close(sub_fig)

        plt.tight_layout()
        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
            print(f"✅ Saved → {save}")
        if show:
            plt.show()
        return fig, axes

    # ------------------------------------------------------------------
    def plot_multiple_genes(
        self,
        target_gene: str,
        highlight_genes: List[str],
        margin_bp: int = 500_000,
        show_exons: bool = True,
        save: Optional[str] = None,
        dpi: int = 300,
        show: bool = True,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Like :meth:`plot_region` but highlight multiple genes.

        All names in *highlight_genes* (plus *target_gene*) are coloured
        using ``theme.target_gene``.
        """
        return self.plot_region(
            target_gene,
            margin_bp=margin_bp,
            show_exons=show_exons,
            highlight_genes=highlight_genes,
            save=save,
            dpi=dpi,
            show=show,
        )
