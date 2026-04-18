"""
Basic unit tests (no network required).
Run with: pytest tests/
"""
import pytest
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # headless

from ideoplot.config import ColorTheme, THEMES
from ideoplot.utils import assign_tracks


# ---------------------------------------------------------------------------
# config tests
# ---------------------------------------------------------------------------

class TestColorTheme:
    def test_defaults(self):
        t = ColorTheme()
        assert t.target_gene == "#B22222"
        assert t.get_cytoband_color("gneg") == "#FFFFFF"

    def test_unknown_stain(self):
        t = ColorTheme()
        assert t.get_cytoband_color("xyz") is None

    def test_themes_dict(self):
        assert "classic" in THEMES
        assert "dark" in THEMES
        assert "colorblind" in THEMES
        assert "pastel" in THEMES
        for name, theme in THEMES.items():
            assert isinstance(theme, ColorTheme), name


# ---------------------------------------------------------------------------
# utils tests
# ---------------------------------------------------------------------------

class TestAssignTracks:
    def _make_gene(self, start_bp, end_bp, name="GENE"):
        return {"start": start_bp, "end": end_bp, "external_name": name}

    def test_non_overlapping_go_to_same_track(self):
        genes = [
            self._make_gene(1_000_000, 2_000_000, "A"),
            self._make_gene(5_000_000, 6_000_000, "B"),
        ]
        tracks = assign_tracks(genes, 0.0, 10.0)
        assert tracks[0] == tracks[1] == 0

    def test_overlapping_go_to_different_tracks(self):
        genes = [
            self._make_gene(1_000_000, 4_000_000, "A"),
            self._make_gene(2_000_000, 5_000_000, "B"),
        ]
        tracks = assign_tracks(genes, 0.0, 10.0)
        assert tracks[0] != tracks[1]

    def test_three_genes_two_tracks(self):
        genes = [
            self._make_gene(0, 2_000_000, "A"),
            self._make_gene(3_000_000, 5_000_000, "B"),  # goes to track 0
            self._make_gene(1_000_000, 4_000_000, "C"),  # overlaps both → track 1 or 2
        ]
        tracks = assign_tracks(genes, 0.0, 10.0)
        assert len(tracks) == 3
        assert max(tracks) >= 1

    def test_single_gene(self):
        genes = [self._make_gene(1_000_000, 2_000_000)]
        tracks = assign_tracks(genes, 0.0, 5.0)
        assert tracks == [0]

    def test_empty(self):
        assert assign_tracks([], 0.0, 10.0) == []
