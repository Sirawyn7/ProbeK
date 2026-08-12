"""Tier classification for BLAST off-target hits.

Precedence (checked in this order per hit):
  1. Overlaps a RepeatMasker ERVK interval -> Tier A (HERV-K family; desired).
  2. Else overlaps an exon of an unrelated gene -> Tier C (real off-target risk).
  3. Else (intron only, or no gene feature at all) -> Tier B (low risk).

An ERVK-overlapping hit is Tier A even if it also happens to overlap an exon
(HERV-K loci can sit inside another gene's intron/exon) — the desired-signal
classification wins, per the spec.
"""

from __future__ import annotations

import pandas as pd

from .interval_backend import IntervalIndex
from .models import BlastHit, ClassifiedHit


def _hits_to_query_df(hits: list[BlastHit]) -> pd.DataFrame:
    rows = []
    for h in hits:
        start, end = h.span
        rows.append({"chrom": h.sseqid, "start": start, "end": end})
    return pd.DataFrame(rows)


def classify_hits(
    hits: list[BlastHit],
    ervk_index: IntervalIndex,
    exon_index: IntervalIndex,
    gene_index: IntervalIndex,
) -> list[ClassifiedHit]:
    if not hits:
        return []

    queries = _hits_to_query_df(hits)

    in_ervk = ervk_index.overlaps_any(queries)
    in_exon = exon_index.overlaps_any(queries)
    ervk_name = ervk_index.first_name(queries)
    exon_name = exon_index.first_name(queries)
    gene_name = gene_index.first_name(queries)

    results = []
    for i, hit in enumerate(hits):
        idx = queries.index[i]
        if bool(in_ervk.loc[idx]):
            results.append(ClassifiedHit(hit=hit, tier="A", locus_label=ervk_name.loc[idx]))
        elif bool(in_exon.loc[idx]):
            results.append(ClassifiedHit(hit=hit, tier="C", locus_label=exon_name.loc[idx]))
        else:
            results.append(ClassifiedHit(hit=hit, tier="B", locus_label=gene_name.loc[idx]))
    return results


def tier_counts(classified: list[ClassifiedHit]) -> dict[str, int]:
    counts = {"A": 0, "B": 0, "C": 0}
    for c in classified:
        counts[c.tier] += 1
    return counts
