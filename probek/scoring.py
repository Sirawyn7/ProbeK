"""Per-target ranking and top-N selection.

Sort order (spec, in priority order):
  1. Tier C count ascending  (fewest real off-targets first)
  2. Tier A count descending (more HERV-K-family coverage is better)
  3. Tier B count ascending
  4. eFISHent `quality` descending (final tiebreaker)
"""

from __future__ import annotations

import pandas as pd

from .classify import tier_counts
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
        counts = tier_counts(classified_by_sequence.get(seq, []))
        enriched = row.to_dict()
        enriched["tier_a_count"] = counts["A"]
        enriched["tier_b_count"] = counts["B"]
        enriched["tier_c_count"] = counts["C"]
        rows.append(enriched)

    ranked = sorted(
        rows,
        key=lambda r: (
            r["tier_c_count"],
            -r["tier_a_count"],
            r["tier_b_count"],
            -float(r["quality"]),
        ),
    )

    short = len(ranked) < top_n
    for i, r in enumerate(ranked, start=1):
        r["rank"] = i

    selected = ranked if short else ranked[:top_n]
    selected_names = {r["name"] for r in selected}

    return TargetResult(
        target=target,
        ranked_rows=ranked,
        selected_names=selected_names,
        short=short,
    )
