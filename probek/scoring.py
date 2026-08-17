"""Per-target ranking and top-N selection.

Rules (in priority order):
  1. Any probe with a real off-target (an exon hit in an unrelated gene) is
     excluded from normal ranking -- it's a disqualifying risk, not something
     to trade off against other factors.
  2. Among the rest ("eligible" probes), rank by eFISHent `quality` descending
     -- higher design/thermodynamic quality wins.
  3. Only on an exact `quality` tie: fewer HERV-K-family hits wins (a probe
     cross-reacting with fewer other HERV-K/HML-2 loci is more specific to
     its own locus).
  4. Only if `quality` AND HERV-K-family hits are both exactly tied: fewer
     intron hits wins (intron hits are spliced out and low-risk, but not
     entirely free).
  5. Hits landing entirely outside any gene never affect ranking.

Excluded (exon-hit) probes are still ranked among themselves the same way
(with fewest exon hits breaking ties first), so that if a target doesn't have
enough eligible candidates to fill `--top-n`, the least-bad excluded probes
backfill the remainder rather than leaving the target short. Ties within
each group fall back to original row order (Python's sort is stable), but
this is a per-group guarantee only -- exclusion always wins regardless of
either probe's original position.
"""

from __future__ import annotations

import pandas as pd

from .classify import feature_counts
from .models import ClassifiedHit, TargetResult


def rank_and_select(
    target: str,
    df: pd.DataFrame,
    classified_by_sequence: dict[str, list[ClassifiedHit]],
    top_n: int,
) -> TargetResult:
    """`df` is this target's surviving (non-FAIL) eFISHent rows."""
    rows = []
    for _, row in df.iterrows():
        seq = str(row["sequence"]).upper()
        counts = feature_counts(classified_by_sequence.get(seq, []))
        enriched = row.to_dict()
        enriched["_hervk_count"] = counts["hervk"]
        enriched["_exon_count"] = counts["exon"]
        enriched["_intron_count"] = counts["intron"]
        rows.append(enriched)

    eligible = [r for r in rows if r["_exon_count"] == 0]
    excluded = [r for r in rows if r["_exon_count"] > 0]

    eligible.sort(key=lambda r: (-float(r["quality"]), r["_hervk_count"], r["_intron_count"]))
    excluded.sort(
        key=lambda r: (r["_exon_count"], -float(r["quality"]), r["_hervk_count"], r["_intron_count"])
    )

    ranked = eligible + excluded
    for i, r in enumerate(ranked, start=1):
        r["rank"] = i

    short = len(ranked) < top_n
    selected = ranked[:top_n]
    selected_names = {r["name"] for r in selected}

    return TargetResult(
        target=target,
        ranked_rows=ranked,
        selected_names=selected_names,
        short=short,
    )
