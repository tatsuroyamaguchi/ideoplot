"""
Configuration: color themes, default parameters, constants.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional


# ---------------------------------------------------------------------------
# Cytoband colour maps
# ---------------------------------------------------------------------------
DEFAULT_CYTOBAND_COLORS: Dict[str, str] = {
    "gneg": "#FFFFFF",
    "gpos25": "#D9D9D9",
    "gpos50": "#A6A6A6",
    "gpos75": "#737373",
    "gpos100": "#262626",
    "acen": "#B22222",
    "gvar": "#E0E0E0",
    "stalk": "#708090",
}


@dataclass
class ColorTheme:
    """A set of colors used for ideoplot plotting.

    Parameters
    ----------
    target_gene : str
        Color for the highlighted (target) gene.
    gene_default : str
        Color for regular protein-coding genes.
    pseudogene : str
        Color for pseudogenes.
    chromosome_outline : str
        Color for chromosome outlines.
    scalebar : str
        Color for the scale bar.
    cytoband : dict, optional
        Mapping from gieStain labels to hex colors.
        Defaults to the standard Giemsa stain palette.
    """
    target_gene: str = "#B22222"
    gene_default: str = "#333333"
    pseudogene: str = "#A0A0A0"
    chromosome_outline: str = "#000000"
    scalebar: str = "#000000"
    cytoband: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_CYTOBAND_COLORS))

    def get_cytoband_color(self, stain: str) -> Optional[str]:
        return self.cytoband.get(stain)


# ---------------------------------------------------------------------------
# Built-in themes
# ---------------------------------------------------------------------------
_THEME_CLASSIC = ColorTheme()

_THEME_DARK = ColorTheme(
    target_gene="#FF6B6B",
    gene_default="#C0C0C0",
    pseudogene="#666666",
    chromosome_outline="#AAAAAA",
    scalebar="#AAAAAA",
    cytoband={
        "gneg": "#1E1E2E",
        "gpos25": "#3A3A5C",
        "gpos50": "#555580",
        "gpos75": "#7070A0",
        "gpos100": "#9090C0",
        "acen": "#FF6B6B",
        "gvar": "#2E2E4E",
        "stalk": "#4A4A7A",
    },
)

_THEME_COLORBLIND = ColorTheme(
    target_gene="#D55E00",   # vermillion
    gene_default="#0072B2",  # blue
    pseudogene="#999999",
    chromosome_outline="#000000",
    scalebar="#000000",
    cytoband={
        "gneg": "#FFFFFF",
        "gpos25": "#E8E8E8",
        "gpos50": "#BBBBBB",
        "gpos75": "#777777",
        "gpos100": "#333333",
        "acen": "#D55E00",
        "gvar": "#EEEEEE",
        "stalk": "#888888",
    },
)

_THEME_PASTEL = ColorTheme(
    target_gene="#E05C5C",
    gene_default="#5B8DB8",
    pseudogene="#C0B9AA",
    chromosome_outline="#999999",
    scalebar="#555555",
    cytoband={
        "gneg": "#FAF9F6",
        "gpos25": "#EDE8DD",
        "gpos50": "#CFC4B0",
        "gpos75": "#AFA090",
        "gpos100": "#7A6E60",
        "acen": "#E8937A",
        "gvar": "#F0EBE0",
        "stalk": "#B0A898",
    },
)

THEMES: Dict[str, ColorTheme] = {
    "classic": _THEME_CLASSIC,
    "dark": _THEME_DARK,
    "colorblind": _THEME_COLORBLIND,
    "pastel": _THEME_PASTEL,
}

# ---------------------------------------------------------------------------
# Chromosome order (hg38 / GRCh38)
# ---------------------------------------------------------------------------
CHROMS_AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMS_ALL = CHROMS_AUTOSOMES + ["chrX", "chrY"]

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------
ENSEMBL_SERVER = "https://rest.ensembl.org"
UCSC_CYTOBAND_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
)
UCSC_CYTOBAND_HG19_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
)
