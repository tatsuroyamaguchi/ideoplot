"""
fetch.py — Data access layer for Ensembl REST API and UCSC cytoband data.

All network-heavy results are cached in memory (and optionally on disk)
so repeated calls don't hit the network every time.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import pickle
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import requests

from .config import (
    ENSEMBL_SERVER,
    UCSC_CYTOBAND_URL,
    UCSC_CYTOBAND_HG19_URL,
)
import re

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_SESSION = requests.Session()
_SESSION.headers.update({
    "Content-Type": "application/json",
    "Accept": "application/json",
})

_DEFAULT_CACHE_DIR = Path.home() / ".cache" / "ideoplot"
_MEMORY_CACHE: Dict[str, object] = {}


def _get_cache_path(key: str, cache_dir: Path) -> Path:
    h = hashlib.md5(key.encode()).hexdigest()
    return cache_dir / f"{h}.pkl"


def _load_disk_cache(key: str, cache_dir: Path) -> Optional[object]:
    p = _get_cache_path(key, cache_dir)
    if p.exists():
        with open(p, "rb") as f:
            return pickle.load(f)
    return None


def _save_disk_cache(key: str, data: object, cache_dir: Path) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = _get_cache_path(key, cache_dir)
    with open(p, "wb") as f:
        pickle.dump(data, f)


# ---------------------------------------------------------------------------
# Cytoband fetcher
# ---------------------------------------------------------------------------
class CytobandFetcher:
    """Fetch and cache Giemsa cytoband data from UCSC.

    Parameters
    ----------
    assembly : str
        Genome assembly. ``"hg38"`` (default) or ``"hg19"``.
    cache : bool
        Cache data to disk in ``~/.cache/ideoplot/``. Default ``True``.
    cache_dir : Path-like, optional
        Override the default cache directory.
    """

    def __init__(
        self,
        assembly: str = "hg38",
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
    ) -> None:
        if assembly not in ("hg38", "hg19"):
            raise ValueError(f"assembly must be 'hg38' or 'hg19', got {assembly!r}")
        self.assembly = assembly
        self.cache = cache
        self.cache_dir = Path(cache_dir) if cache_dir else _DEFAULT_CACHE_DIR
        self._url = UCSC_CYTOBAND_URL if assembly == "hg38" else UCSC_CYTOBAND_HG19_URL
        self._df: Optional[pd.DataFrame] = None

    def load(self) -> pd.DataFrame:
        """Return the full cytoband DataFrame (cached after first call).

        Columns: ``chrom``, ``chromStart``, ``chromEnd``, ``name``, ``gieStain``
        (start/end in base-pairs).
        """
        if self._df is not None:
            return self._df

        key = f"cytoband_{self.assembly}"
        if key in _MEMORY_CACHE:
            self._df = _MEMORY_CACHE[key]
            return self._df

        if self.cache:
            cached = _load_disk_cache(key, self.cache_dir)
            if cached is not None:
                self._df = cached
                _MEMORY_CACHE[key] = cached
                return self._df

        cols = ["chrom", "chromStart", "chromEnd", "name", "gieStain"]
        resp = _SESSION.get(self._url, stream=True, timeout=60)
        resp.raise_for_status()

        with gzip.open(io.BytesIO(resp.content)) as gz:
            df = pd.read_csv(gz, sep="\t", names=cols)

        self._df = df
        _MEMORY_CACHE[key] = df
        if self.cache:
            _save_disk_cache(key, df, self.cache_dir)
        return df

    def for_region(
        self,
        chrom: str,
        start_bp: int,
        end_bp: int,
    ) -> pd.DataFrame:
        """Return cytoband rows that overlap a genomic region.

        Parameters
        ----------
        chrom : str
            Chromosome name, e.g. ``"chr2"`` or ``"2"``.
        start_bp, end_bp : int
            Region bounds in base-pairs.
        """
        df = self.load()
        chrom_norm = chrom if chrom.startswith("chr") else f"chr{chrom}"
        mask = (
            (df["chrom"] == chrom_norm)
            & (df["chromEnd"] >= start_bp)
            & (df["chromStart"] <= end_bp)
        )
        return df[mask].copy()

    def chromosome_sizes(self) -> Dict[str, Tuple[int, int]]:
        """Return ``{chrom: (min_start_bp, max_end_bp)}`` for every chromosome."""
        df = self.load()
        grp = df.groupby("chrom")
        return {
            ch: (int(sub["chromStart"].min()), int(sub["chromEnd"].max()))
            for ch, sub in grp
        }


# ---------------------------------------------------------------------------
# Ensembl client
# ---------------------------------------------------------------------------
class EnsemblClient:
    """Light wrapper around the Ensembl REST API.

    Parameters
    ----------
    species : str
        Ensembl species string, default ``"homo_sapiens"``.
    cache : bool
        Cache responses to disk. Default ``True``.
    cache_dir : Path-like, optional
        Override the default cache directory.
    max_retries : int
        Number of retries on HTTP 429 / 503. Default 3.
    retry_delay : float
        Seconds to wait between retries. Default 2.0.
    """

    def __init__(
        self,
        species: str = "homo_sapiens",
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
        max_retries: int = 3,
        retry_delay: float = 2.0,
    ) -> None:
        self.species = species
        self.cache = cache
        self.cache_dir = Path(cache_dir) if cache_dir else _DEFAULT_CACHE_DIR
        self.max_retries = max_retries
        self.retry_delay = retry_delay

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _get(self, path: str, params: Optional[dict] = None) -> dict:
        url = f"{ENSEMBL_SERVER}/{path.lstrip('/')}"
        key = f"GET|{url}|{json.dumps(params, sort_keys=True)}"
        return self._cached_request(key, lambda: _SESSION.get(url, params=params, timeout=30))

    def _post(self, path: str, payload: dict, params: Optional[dict] = None) -> dict:
        url = f"{ENSEMBL_SERVER}/{path.lstrip('/')}"
        key = f"POST|{url}|{json.dumps(payload, sort_keys=True)}"
        return self._cached_request(key, lambda: _SESSION.post(url, json=payload, params=params, timeout=60))

    def _cached_request(self, key: str, requester) -> dict:
        if key in _MEMORY_CACHE:
            return _MEMORY_CACHE[key]

        if self.cache:
            cached = _load_disk_cache(key, self.cache_dir)
            if cached is not None:
                _MEMORY_CACHE[key] = cached
                return cached

        resp = None
        for attempt in range(self.max_retries):
            try:
                resp = requester()
                if resp.status_code in (429, 503):
                    wait = self.retry_delay * (attempt + 1)
                    warnings.warn(f"Ensembl rate-limited, retrying in {wait}s…")
                    time.sleep(wait)
                    continue
                resp.raise_for_status()
                data = resp.json()
                break
            except requests.RequestException as exc:
                if attempt == self.max_retries - 1:
                    raise RuntimeError(f"Ensembl request failed: {exc}") from exc
                time.sleep(self.retry_delay)
        else:
            raise RuntimeError(f"Ensembl returned {resp.status_code} after retries")

        _MEMORY_CACHE[key] = data
        if self.cache:
            _save_disk_cache(key, data, self.cache_dir)
        return data

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def lookup_gene(self, symbol: str) -> dict:
        """Look up a single gene by symbol.

        Returns the Ensembl lookup dict (contains ``seq_region_name``,
        ``start``, ``end``, ``strand``, ``id``, ``biotype``, etc.).
        """
        data = self._get(
            f"lookup/symbol/{self.species}/{symbol}",
            params={"expand": 0},
        )
        if not data or "error" in data:
            raise ValueError(
                f"Gene '{symbol}' not found in Ensembl for species '{self.species}'."
            )
        return data

    def lookup_genes_batch(self, symbols: List[str]) -> Dict[str, dict]:
        """Look up multiple genes in a single POST request.

        Returns ``{symbol: info_dict}``; missing symbols are omitted.
        """
        # Ensembl batch endpoint has a ~1000-symbol limit
        result: Dict[str, dict] = {}
        chunk_size = 500
        for i in range(0, len(symbols), chunk_size):
            chunk = symbols[i : i + chunk_size]
            data = self._post(
                f"lookup/symbol/{self.species}",
                {"symbols": chunk},
            )
            for sym, info in data.items():
                if info:
                    result[sym] = info
        return result

    def genes_in_region(
        self,
        chrom: str,
        start_bp: int,
        end_bp: int,
        features: List[str] = ("gene",),
    ) -> List[dict]:
        """Return all features overlapping a genomic region.

        Parameters
        ----------
        chrom : str
            Bare chromosome number or name (e.g. ``"2"`` or ``"chr2"``).
        start_bp, end_bp : int
            Region bounds in base-pairs.
        features : list of str
            Ensembl feature types to retrieve. Default ``["gene"]``.
        """
        chrom_norm = chrom.lstrip("chr")
        feature_str = ";".join(f"feature={f}" for f in features)
        path = f"overlap/region/{self.species}/{chrom_norm}:{start_bp}-{end_bp}"
        # The overlap endpoint accepts repeated query params; we rebuild
        url = f"{ENSEMBL_SERVER}/{path}?{'&'.join(f'feature={f}' for f in features)}"
        key = f"GET|{url}"
        if key in _MEMORY_CACHE:
            return _MEMORY_CACHE[key]
        resp = _SESSION.get(url, timeout=60)
        resp.raise_for_status()
        result = resp.json()
        _MEMORY_CACHE[key] = result
        if self.cache:
            _save_disk_cache(key, result, self.cache_dir)
        return result

    def lookup_ids_with_exons(self, gene_ids: List[str]) -> Dict[str, dict]:
        """Batch-fetch gene details including transcripts and exons.

        Returns ``{gene_id: detail_dict}`` where each detail_dict has
        a ``"Transcript"`` key containing exon positions.
        """
        result: Dict[str, dict] = {}
        chunk_size = 200
        for i in range(0, len(gene_ids), chunk_size):
            chunk = gene_ids[i : i + chunk_size]
            try:
                data = self._post("lookup/id", {"ids": chunk, "expand": 1})
                result.update(data)
            except Exception as exc:
                warnings.warn(f"Failed to fetch exon data for a batch: {exc}")
        return result

    @staticmethod
    def clear_memory_cache() -> None:
        """Clear the in-memory request cache (disk cache is unchanged)."""
        _MEMORY_CACHE.clear()


# ---------------------------------------------------------------------------
# UCSC Client
# ---------------------------------------------------------------------------
class UcscClient:
    """Fallback client for looking up markers (e.g. STS markers) via UCSC search.

    Parameters
    ----------
    assembly : str
        Genome assembly. ``"hg38"`` (default) or ``"hg19"``.
    cache : bool
        Cache data to disk. Default ``True``.
    cache_dir : Path-like, optional
        Override the default cache directory.
    """

    def __init__(
        self,
        assembly: str = "hg38",
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
    ) -> None:
        self.assembly = assembly
        self.cache = cache
        self.cache_dir = Path(cache_dir) if cache_dir else _DEFAULT_CACHE_DIR

    def search_position(self, term: str) -> Optional[dict]:
        url = f"https://api.genome.ucsc.edu/search?genome={self.assembly}&search={term}"
        key = f"GET|{url}"

        if key in _MEMORY_CACHE:
            data = _MEMORY_CACHE[key]
        else:
            if self.cache:
                cached = _load_disk_cache(key, self.cache_dir)
                if cached is not None:
                    _MEMORY_CACHE[key] = cached
                    data = cached
                else:
                    try:
                        resp = _SESSION.get(url, timeout=30)
                        resp.raise_for_status()
                        data = resp.json()
                        _MEMORY_CACHE[key] = data
                        _save_disk_cache(key, data, self.cache_dir)
                    except Exception as exc:
                        warnings.warn(f"UCSC API search failed for '{term}': {exc}")
                        return None
            else:
                try:
                    resp = _SESSION.get(url, timeout=30)
                    resp.raise_for_status()
                    data = resp.json()
                except Exception as exc:
                    warnings.warn(f"UCSC API search failed for '{term}': {exc}")
                    return None

        matches = data.get("positionMatches", [])
        if not matches:
            return None

        for track_match in matches:
            if track_match.get("matches"):
                first_match = track_match["matches"][0]
                highlight = first_match.get("highlight")
                pos_str = first_match.get("position")
                
                # Highlight contains unpadded exact coords
                if highlight:
                    m = re.search(r'(chr[\w\d]+):(\d+)-(\d+)', highlight)
                    if m:
                        return {"chrom": m.group(1), "start": int(m.group(2)), "end": int(m.group(3))}
                if pos_str:
                    m = re.search(r'(chr[\w\d]+):(\d+)-(\d+)', pos_str)
                    if m:
                        return {"chrom": m.group(1), "start": int(m.group(2)), "end": int(m.group(3))}
        return None


# ---------------------------------------------------------------------------
# High-level helper
# ---------------------------------------------------------------------------
def fetch_gene_dataframe(
    genes: List[str],
    client: Optional[EnsemblClient] = None,
    assembly: str = "hg38",
) -> pd.DataFrame:
    """Return a DataFrame with genomic positions for the given gene symbols.

    Columns: ``gene``, ``chrom``, ``start`` (Mb), ``end`` (Mb), ``pos`` (Mb),
    ``strand``, ``ensembl_id``, ``biotype``.

    ``pos`` is ``start`` for forward-strand genes and ``end`` for reverse-strand
    genes (i.e. the TSS).
    """
    if client is None:
        client = EnsemblClient()

    raw = client.lookup_genes_batch(genes)
    rows = []
    
    missing_genes = [sym for sym in genes if sym not in raw or raw[sym] is None]
    ucsc_results = {}
    if missing_genes:
        ucsc = UcscClient(assembly=assembly, cache=client.cache, cache_dir=client.cache_dir)
        for sym in missing_genes:
            res = ucsc.search_position(sym)
            if res:
                ucsc_results[sym] = res

    for sym in genes:
        if sym in ucsc_results:
            res = ucsc_results[sym]
            start_bp = res["start"]
            end_bp = res["end"]
            rows.append({
                "gene": sym,
                "chrom": res["chrom"],
                "start": start_bp / 1_000_000,
                "end": end_bp / 1_000_000,
                "pos": start_bp / 1_000_000,
                "strand": 1,
                "ensembl_id": "",
                "biotype": "marker",
            })
            continue

        info = raw.get(sym)
        if not info:
            warnings.warn(f"Gene '{sym}' not found — skipped.")
            continue
            
        strand = info.get("strand", 1)
        start_bp = info["start"]
        end_bp = info["end"]
        pos_bp = start_bp if strand == 1 else end_bp
        rows.append(
            {
                "gene": sym,
                "chrom": "chr" + str(info["seq_region_name"]),
                "start": start_bp / 1_000_000,
                "end": end_bp / 1_000_000,
                "pos": pos_bp / 1_000_000,
                "strand": strand,
                "ensembl_id": info.get("id", ""),
                "biotype": info.get("biotype", ""),
            }
        )
    return pd.DataFrame(rows)

