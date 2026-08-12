"""Swappable interval-overlap backend.

Nothing outside this module and annotation.py should import pyranges or
pyranges1 directly, so the backend can change without touching classification,
scoring, or reporting code.

Note: on this project's dev machine (Python 3.14, no python3.14-dev headers
installed), neither `pyranges` (classic, Cython-based) nor `pyranges1` (its
current successor on PyPI; `pyranges2` does not exist as a package) can build
from source, so the pure pandas/numpy `_CustomIndex` fallback is what actually
runs here. Installing `python3.14-dev` and re-running `pip install .[pyranges]`
/ `.[pyranges1]` may enable the faster backend later without any code changes.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class IntervalFeature:
    chrom: str
    start: int
    end: int
    name: str | None = None


class IntervalIndex(Protocol):
    def overlaps_any(self, queries: pd.DataFrame) -> pd.Series:
        """queries has columns [chrom, start, end]; returns a bool Series
        aligned to queries.index."""

    def first_name(self, queries: pd.DataFrame) -> pd.Series:
        """Returns the `name` of an overlapping feature per query row (first
        match, arbitrary among ties), or None if no overlap."""


def select_backend() -> str:
    forced = os.environ.get("PROBEK_INTERVAL_BACKEND")
    if forced:
        return forced
    try:
        import pyranges  # noqa: F401

        return "pyranges"
    except ImportError:
        pass
    try:
        import pyranges1  # noqa: F401

        return "pyranges1"
    except ImportError:
        pass
    return "custom"


def build_index(features: list[IntervalFeature]) -> IntervalIndex:
    backend = select_backend()
    if backend == "pyranges":
        return _PyRangesIndex(features)
    if backend == "pyranges1":
        return _PyRanges1Index(features)
    return _CustomIndex(features)


def _features_to_df(features: list[IntervalFeature]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Chromosome": [f.chrom for f in features],
            "Start": [f.start for f in features],
            "End": [f.end for f in features],
            "Name": [f.name for f in features],
        }
    )


class _PyRangesIndex:
    def __init__(self, features: list[IntervalFeature]):
        import pyranges as pr

        self._pr = pr.PyRanges(_features_to_df(features))

    def _query_pr(self, queries: pd.DataFrame):
        import pyranges as pr

        qdf = queries.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
        qdf = qdf.reset_index().rename(columns={"index": "_qidx"})
        return pr.PyRanges(qdf)

    def overlaps_any(self, queries: pd.DataFrame) -> pd.Series:
        joined = self._query_pr(queries).join(self._pr)
        hit_idx = set(joined.df["_qidx"]) if len(joined) > 0 else set()
        return pd.Series([i in hit_idx for i in queries.index], index=queries.index)

    def first_name(self, queries: pd.DataFrame) -> pd.Series:
        joined = self._query_pr(queries).join(self._pr)
        result = pd.Series([None] * len(queries), index=queries.index, dtype=object)
        if len(joined) == 0:
            return result
        jdf = joined.df.drop_duplicates(subset="_qidx")
        for qidx, name in zip(jdf["_qidx"], jdf["Name"]):
            result.loc[qidx] = name
        return result


class _PyRanges1Index:
    def __init__(self, features: list[IntervalFeature]):
        import pyranges1 as pr1

        self._pr = pr1.PyRanges(_features_to_df(features))

    def _query_pr(self, queries: pd.DataFrame):
        import pyranges1 as pr1

        qdf = queries.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
        qdf = qdf.reset_index().rename(columns={"index": "_qidx"})
        return pr1.PyRanges(qdf)

    def overlaps_any(self, queries: pd.DataFrame) -> pd.Series:
        joined = self._query_pr(queries).join(self._pr)
        jdf = joined.df if hasattr(joined, "df") else joined
        hit_idx = set(jdf["_qidx"]) if len(jdf) > 0 else set()
        return pd.Series([i in hit_idx for i in queries.index], index=queries.index)

    def first_name(self, queries: pd.DataFrame) -> pd.Series:
        joined = self._query_pr(queries).join(self._pr)
        jdf = joined.df if hasattr(joined, "df") else joined
        result = pd.Series([None] * len(queries), index=queries.index, dtype=object)
        if len(jdf) == 0:
            return result
        jdf = jdf.drop_duplicates(subset="_qidx")
        for qidx, name in zip(jdf["_qidx"], jdf["Name"]):
            result.loc[qidx] = name
        return result


class _CustomIndex:
    """Pure pandas/numpy fallback: per-chromosome sorted-start + running
    max-end array, queried via binary search. No third-party interval-tree
    dependency required."""

    def __init__(self, features: list[IntervalFeature]):
        self._by_chrom: dict[str, dict[str, np.ndarray]] = {}
        by_chrom: dict[str, list[IntervalFeature]] = {}
        for f in features:
            by_chrom.setdefault(f.chrom, []).append(f)

        for chrom, feats in by_chrom.items():
            feats = sorted(feats, key=lambda f: f.start)
            starts = np.array([f.start for f in feats], dtype=np.int64)
            ends = np.array([f.end for f in feats], dtype=np.int64)
            names = np.array([f.name for f in feats], dtype=object)
            max_end = np.maximum.accumulate(ends)
            self._by_chrom[chrom] = {
                "starts": starts,
                "ends": ends,
                "names": names,
                "max_end": max_end,
            }

    def _find_overlaps(self, chrom: str, start: int, end: int) -> list[int]:
        """Return indices (within this chrom's arrays) of features overlapping
        [start, end], via the augmented max-end binary-search trick."""
        idx = self._by_chrom.get(chrom)
        if idx is None:
            return []
        starts, ends, max_end = idx["starts"], idx["ends"], idx["max_end"]
        # candidate window: any feature with start <= end could overlap
        hi = int(np.searchsorted(starts, end, side="right"))
        if hi == 0:
            return []
        # Among features[0:hi], only those whose max_end-so-far >= start are
        # possibly relevant; find the earliest position where that holds.
        lo = int(np.searchsorted(max_end[:hi], start, side="left"))
        matches = []
        for i in range(lo, hi):
            if ends[i] >= start and starts[i] <= end:
                matches.append(i)
        return matches

    def overlaps_any(self, queries: pd.DataFrame) -> pd.Series:
        result = []
        for chrom, start, end in zip(queries["chrom"], queries["start"], queries["end"]):
            result.append(len(self._find_overlaps(chrom, start, end)) > 0)
        return pd.Series(result, index=queries.index)

    def first_name(self, queries: pd.DataFrame) -> pd.Series:
        result = []
        for chrom, start, end in zip(queries["chrom"], queries["start"], queries["end"]):
            matches = self._find_overlaps(chrom, start, end)
            if matches:
                names = self._by_chrom[chrom]["names"]
                result.append(names[matches[0]])
            else:
                result.append(None)
        return pd.Series(result, index=queries.index, dtype=object)
