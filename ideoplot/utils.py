"""
utils.py — Low-level drawing helpers and layout algorithms.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib.patches import PathPatch


# ---------------------------------------------------------------------------
# Chromosome pill shape (rounded rectangle)
# ---------------------------------------------------------------------------

def _bezier_pill_path(
    xc: float,
    y0: float,
    y1: float,
    rx: float,
    ry: float,
) -> Path:
    """Build a Bezier-approximated pill (vertical rounded rectangle).

    Parameters
    ----------
    xc : float   Centre x of the pill.
    y0 : float   Top y  (smaller value ≡ top of the plot when y increases downward).
    y1 : float   Bottom y.
    rx : float   Half-width in x-units.
    ry : float   Corner radius in y-units.
    """
    k = 0.5523  # Bezier approximation constant for a circle
    xl, xr = xc - rx, xc + rx
    verts = [
        (xl, y0 + ry), (xl, y1 - ry),
        (xl, y1 - ry + k * ry), (xc - k * rx, y1), (xc, y1),
        (xc + k * rx, y1), (xr, y1 - ry + k * ry), (xr, y1 - ry),
        (xr, y0 + ry),
        (xr, y0 + ry - k * ry), (xc + k * rx, y0), (xc, y0),
        (xc - k * rx, y0), (xl, y0 + ry - k * ry), (xl, y0 + ry),
        (xl, y0 + ry),
    ]
    codes = (
        [Path.MOVETO, Path.LINETO]
        + [Path.CURVE4] * 6
        + [Path.LINETO]
        + [Path.CURVE4] * 6
        + [Path.CLOSEPOLY]
    )
    return Path(verts, codes)


def _bezier_pinched_pill_path(
    xc: float,
    y0: float,
    y1: float,
    cy_top: float,
    cy_bot: float,
    rx: float,
    ry: float,
) -> Path:
    """Build a vertical pill with an inward pinch at the centromere."""
    k = 0.5523
    xl, xr = xc - rx, xc + rx
    xl_pinch, xr_pinch = xc - rx * 0.25, xc + rx * 0.25
    
    # If the centromere is too close to the end caps, don't overlap radii
    cy_top = max(cy_top, y0 + ry + 0.01)
    cy_bot = min(cy_bot, y1 - ry - 0.01)
    if cy_bot <= cy_top:
        return _bezier_pill_path(xc, y0, y1, rx, ry)
        
    cy_mid = (cy_top + cy_bot) / 2
    c_off_t = (cy_mid - cy_top) * 0.5
    c_off_b = (cy_bot - cy_mid) * 0.5

    verts = [
        # Start top-left
        (xl, y0 + ry),
        # Line to left centromere top
        (xl, cy_top),
        # Curve into left pinch
        (xl, cy_top + c_off_t), (xl_pinch, cy_mid - c_off_t), (xl_pinch, cy_mid),
        # Curve out of left pinch
        (xl_pinch, cy_mid + c_off_b), (xl, cy_bot - c_off_b), (xl, cy_bot),
        # Line to bottom-left
        (xl, y1 - ry),
        # Bottom curve
        (xl, y1 - ry + k * ry), (xc - k * rx, y1), (xc, y1),
        (xc + k * rx, y1), (xr, y1 - ry + k * ry), (xr, y1 - ry),
        # Line to right centromere bottom
        (xr, cy_bot),
        # Curve into right pinch (going UP! so y decreases)
        (xr, cy_bot - c_off_b), (xr_pinch, cy_mid + c_off_b), (xr_pinch, cy_mid),
        # Curve out of right pinch
        (xr_pinch, cy_mid - c_off_t), (xr, cy_top + c_off_t), (xr, cy_top),
        # Line to top-right
        (xr, y0 + ry),
        # Top curve
        (xr, y0 + ry - k * ry), (xc + k * rx, y0), (xc, y0),
        (xc - k * rx, y0), (xl, y0 + ry - k * ry), (xl, y0 + ry),
        
        # Close
        (xl, y0 + ry),
    ]
    
    codes = (
        [Path.MOVETO, Path.LINETO]
        + [Path.CURVE4] * 3
        + [Path.CURVE4] * 3
        + [Path.LINETO]
        + [Path.CURVE4] * 6
        + [Path.LINETO]
        + [Path.CURVE4] * 3
        + [Path.CURVE4] * 3
        + [Path.LINETO]
        + [Path.CURVE4] * 6
        + [Path.CLOSEPOLY]
    )
    
    return Path(verts, codes)


def draw_chromosome_body(
    ax: Axes,
    x: float,
    y0: float,
    y1: float,
    bands_df,
    color_map: dict,
    width: float = 0.3,
    aspect_correction: float = 1.0,
) -> PathPatch:
    """Draw a single chromosome body with Giemsa bands.

    Parameters
    ----------
    ax
        Target axes.
    x
        Horizontal centre position.
    y0, y1
        Top and bottom of the chromosome in data coordinates.
    bands_df
        DataFrame with columns ``chromStart``, ``chromEnd``, ``gieStain``
        (already in Mb).
    color_map
        Mapping ``{gieStain: hex_color}``.
    width
        Full width of the chromosome bar.
    aspect_correction
        Ratio ``px_per_x_unit / px_per_y_unit`` — used to make the rounded
        cap look circular despite non-square axes.

    Returns
    -------
    PathPatch
        The outline patch (used as clip-path for band rectangles).
    """
    rx = width / 2
    ry = rx * aspect_correction

    cen_bands = bands_df[bands_df["gieStain"] == "acen"]
    has_cen = not cen_bands.empty
    if has_cen:
        cy_top = cen_bands["chromStart"].min()
        cy_bot = cen_bands["chromEnd"].max()
        path = _bezier_pinched_pill_path(x, y0, y1, cy_top, cy_bot, rx, ry)
    else:
        path = _bezier_pill_path(x, y0, y1, rx, ry)

    outline = PathPatch(path, fill=False, edgecolor="black", linewidth=0.6, zorder=5)
    ax.add_patch(outline)

    for _, row in bands_df.iterrows():
        color = color_map.get(row["gieStain"])
        if color is None:
            continue
        rect = mpatches.Rectangle(
            (x - rx, row["chromStart"]),
            rx * 2,
            row["chromEnd"] - row["chromStart"],
            facecolor=color,
            edgecolor="none",
            zorder=4,
        )
        ax.add_patch(rect)
        rect.set_clip_path(outline)

    return outline


# ---------------------------------------------------------------------------
# Aspect correction
# ---------------------------------------------------------------------------

def compute_aspect_correction(ax: Axes, fig: plt.Figure) -> float:
    """Return ``px_per_x_unit / px_per_y_unit`` for the given axes.

    Must be called **after** ``fig.canvas.draw()`` so that the transform
    is finalised.
    """
    trans = ax.transData
    p0 = trans.transform((0, 0))
    px_per_x = abs(trans.transform((1, 0))[0] - p0[0])
    px_per_y = abs(trans.transform((0, 1))[1] - p0[1])
    if px_per_y == 0:
        return 1.0
    return px_per_x / px_per_y


# ---------------------------------------------------------------------------
# Gene track layout
# ---------------------------------------------------------------------------

def assign_tracks(
    genes: List[dict],
    region_start_mb: float,
    region_end_mb: float,
    char_width_fraction: float = 0.012,
) -> List[int]:
    """Greedily assign genes to non-overlapping horizontal tracks.

    Parameters
    ----------
    genes
        List of gene dicts (each must have ``start``, ``end``,
        ``external_name`` keys, all in bp).
    region_start_mb, region_end_mb
        Displayed region in Mb (used to scale label width estimate).
    char_width_fraction
        Width (as a fraction of the region span) attributed to each
        character of the gene name label.

    Returns
    -------
    list of int
        Track index for each gene (same order as *genes*).
    """
    char_width_mb = (region_end_mb - region_start_mb) * char_width_fraction
    tracks: List[float] = []   # end-position of last gene on each track
    result: List[int] = []

    for gene in genes:
        g_start_mb = max(region_start_mb, gene["start"] / 1_000_000)
        g_end_mb = min(region_end_mb, gene["end"] / 1_000_000)
        g_name = gene.get("external_name", "Unknown")
        effective_end = max(g_end_mb, g_start_mb + len(g_name) * char_width_mb)

        placed = False
        for i, track_end in enumerate(tracks):
            if g_start_mb > track_end + 0.01:
                tracks[i] = effective_end
                result.append(i)
                placed = True
                break

        if not placed:
            tracks.append(effective_end)
            result.append(len(tracks) - 1)

    return result


# ---------------------------------------------------------------------------
# Gene body drawing
# ---------------------------------------------------------------------------

def draw_arrow_gene(
    ax: Axes,
    g_start_mb: float,
    g_end_mb: float,
    y_pos: float,
    strand: int,
    color: str,
    height: float = 0.25,
    arrow_head_mb: float = 0.005,
    zorder: int = 5,
) -> None:
    """Draw a gene as a filled arrow polygon.

    Parameters
    ----------
    ax
        Target axes.
    g_start_mb, g_end_mb
        Gene bounds in Mb.
    y_pos
        Bottom y-coordinate for the gene body.
    strand
        +1 for forward, -1 for reverse.
    color
        Fill and edge colour.
    height
        Arrow body height in data y-units.
    arrow_head_mb
        Arrow head length in Mb.
    zorder
        Drawing order.
    """
    g_len = g_end_mb - g_start_mb
    actual_head = min(arrow_head_mb, g_len) if g_len > 0 else 0.001

    if strand == 1:
        verts = [
            (g_start_mb, y_pos),
            (g_end_mb - actual_head, y_pos),
            (g_end_mb, y_pos + height / 2),
            (g_end_mb - actual_head, y_pos + height),
            (g_start_mb, y_pos + height),
        ]
    else:
        verts = [
            (g_end_mb, y_pos),
            (g_start_mb + actual_head, y_pos),
            (g_start_mb, y_pos + height / 2),
            (g_start_mb + actual_head, y_pos + height),
            (g_end_mb, y_pos + height),
        ]

    poly = mpatches.Polygon(
        verts, facecolor=color, edgecolor=color, linewidth=0.5, zorder=zorder
    )
    ax.add_patch(poly)


def draw_exon_gene(
    ax: Axes,
    g_start_mb: float,
    g_end_mb: float,
    y_pos: float,
    strand: int,
    exons: List[dict],
    region_start_mb: float,
    region_end_mb: float,
    color: str,
    height: float = 0.25,
    min_exon_width_mb: float = 0.001,
    arrow_head_mb: float = 0.003,
    zorder: int = 5,
) -> None:
    """Draw a gene with intron line, strand arrow, and exon blocks.

    Parameters
    ----------
    exons
        List of exon dicts with ``start`` and ``end`` keys (in bp).
    min_exon_width_mb
        Minimum rendered exon width so thin exons remain visible.
    """
    line_y = y_pos + height / 2
    g_len = g_end_mb - g_start_mb

    # Intron line
    ax.plot(
        [g_start_mb, g_end_mb], [line_y, line_y],
        color=color, linewidth=1.0, zorder=zorder
    )

    # Strand arrow
    head_len = min(arrow_head_mb, g_len / 2) if g_len > 0 else 0.001
    head_h = height * 0.35
    if strand == 1:
        ax.plot(
            [g_end_mb - head_len, g_end_mb, g_end_mb - head_len],
            [line_y + head_h, line_y, line_y - head_h],
            color=color, linewidth=1.5, zorder=zorder
        )
    else:
        ax.plot(
            [g_start_mb + head_len, g_start_mb, g_start_mb + head_len],
            [line_y + head_h, line_y, line_y - head_h],
            color=color, linewidth=1.5, zorder=zorder
        )

    # Exon blocks
    for ex in exons:
        ex_s = max(region_start_mb, ex["start"] / 1_000_000)
        ex_e = min(region_end_mb, ex["end"] / 1_000_000)
        if ex_e < region_start_mb or ex_s > region_end_mb:
            continue
        ex_w = max(ex_e - ex_s, min_exon_width_mb)
        rect = mpatches.Rectangle(
            (ex_s, line_y - height * 0.3),
            ex_w,
            height * 0.6,
            facecolor=color, edgecolor="none", zorder=zorder + 1
        )
        ax.add_patch(rect)


# ---------------------------------------------------------------------------
# Scale bar
# ---------------------------------------------------------------------------

def add_scale_bar(
    ax: Axes,
    x: float,
    y: float,
    length_mb: float,
    color: str = "black",
    fontsize: int = 10,
) -> None:
    """Draw a vertical scale bar with a length label."""
    ax.plot([x, x], [y, y + length_mb], color=color, linewidth=1.5)
    ax.text(
        x + 0.15, (y + y + length_mb) / 2,
        f"{length_mb:.0f} Mb",
        va="center", fontsize=fontsize, color=color,
    )
