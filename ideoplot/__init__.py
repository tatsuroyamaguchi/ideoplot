"""
ideoplot: Human genome visualization library.

Quick start
-----------
>>> from ideoplot import IdeoplotPlotter, GenomeViewer
>>> # Genome-wide ideoplot with gene markers
>>> plotter = IdeoplotPlotter()
>>> plotter.plot(genes=["BRCA1", "BRCA2", "TP53"], save="my_ideoplot.png")
>>>
>>> # Local gene view with exon structure
>>> viewer = GenomeViewer()
>>> viewer.plot_region("MSH2", margin_bp=500_000, show_exons=True, save="msh2_region.png")
"""

from .core import IdeoplotPlotter, GenomeViewer
from .fetch import EnsemblClient, CytobandFetcher
from .config import ColorTheme, THEMES

__version__ = "1.0.0"
__author__ = "ideoplot package"

__all__ = [
    "IdeoplotPlotter",
    "GenomeViewer",
    "EnsemblClient",
    "CytobandFetcher",
    "ColorTheme",
    "THEMES",
]
