from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .classify import feature_counts
from .io_utils import write_fasta
from .models import ClassifiedHit, TargetResult

logger = logging.getLogger(__name__)


def _plain_feature_label(hit: ClassifiedHit) -> str:
    """Plain-English description of where one off-target hit landed --
    matches the four categories a reader actually cares about: HERV-K
    family, exon of a gene, intron of a gene, or outside any gene."""
    if hit.tier == "A":
        return f"HERV-K family: {hit.locus_label}" if hit.locus_label else "HERV-K family"
    if hit.tier == "C":
        return f"exon of {hit.locus_label}" if hit.locus_label else "exon"
    return f"intron of {hit.locus_label}" if hit.locus_label else "outside any gene"


def format_off_target_loci(classified: list[ClassifiedHit]) -> str:
    """One entry per off-target hit: '<accession>:<start>-<end> (<feature>)',
    e.g. 'NC_000005.10:151268868-151268889 (exon of GM2A)'."""
    parts = []
    for c in classified:
        start, end = c.hit.span
        parts.append(f"{c.hit.sseqid}:{start}-{end} ({_plain_feature_label(c)})")
    return "; ".join(parts)


def flagged_genes(classified: list[ClassifiedHit]) -> str:
    """Comma-separated, de-duplicated gene symbols with a real (exon)
    off-target hit -- the one detail worth seeing without reading the full
    off_target_loci breakdown. Empty string if there are none."""
    genes: list[str] = []
    seen: set[str] = set()
    for c in classified:
        if c.tier == "C" and c.locus_label and c.locus_label not in seen:
            seen.add(c.locus_label)
            genes.append(c.locus_label)
    return ", ".join(genes)


def off_target_risk_label(exon_hits: int) -> str:
    """Coarse at-a-glance risk category, driven purely by exon hits -- the
    only off-target category that represents real cross-hybridization risk.
    0 = Low, 1-2 = Moderate, 3+ = High."""
    if exon_hits == 0:
        return "Low"
    if exon_hits <= 2:
        return "Moderate"
    return "High"


def enrich_target_result(
    result: TargetResult, classified_by_sequence: dict[str, list[ClassifiedHit]]
) -> None:
    """Replaces each row's internal sort-key counts (used only for ranking,
    in probek.scoring's vocabulary) with plain-English off-target detail: a
    per-locus feature+location breakdown, separate HERV-K family / exon /
    intron / outside-gene hit counts, an at-a-glance risk category, and the
    specific genes flagged for real (exon) off-target risk. Mutates `result`
    in place so both the per-target audit CSV and the combined
    final-selection CSV pick up the same enriched columns."""
    for row in result.ranked_rows:
        seq = str(row["sequence"]).upper()
        classified = classified_by_sequence.get(seq, [])
        row["off_target_loci"] = format_off_target_loci(classified)
        counts = feature_counts(classified)
        row["hervk_family_hits"] = counts["hervk"]
        row["exon_hits"] = counts["exon"]
        row["intron_hits"] = counts["intron"]
        row["outside_gene_hits"] = counts["outside_gene"]
        row["off_target_risk"] = off_target_risk_label(counts["exon"])
        row["flagged_genes"] = flagged_genes(classified)
        for internal_key in ("_hervk_count", "_exon_count", "_intron_count"):
            row.pop(internal_key, None)
        row["selected"] = row["name"] in result.selected_names


# Columns worth seeing without scrolling, in the order they should appear.
# "target" and "target_short" only exist in final_selection.csv; harmless to
# list them here since reordering silently skips columns that aren't present.
_PRIORITY_COLUMNS = [
    "target",
    "name",
    "sequence",
    "selected",
    "rank",
    "target_short",
    "off_target_risk",
    "flagged_genes",
    "exon_hits",
    "intron_hits",
    "hervk_family_hits",
    "outside_gene_hits",
    "quality",
    "recommendation",
]


def _reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Puts the glanceable summary columns first and the dense, detailed
    off_target_loci breakdown last, with everything else (the original
    eFISHent columns) in between in their existing order."""
    priority = [c for c in _PRIORITY_COLUMNS if c in df.columns]
    tail = [c for c in ("off_target_loci",) if c in df.columns]
    middle = [c for c in df.columns if c not in priority and c not in tail]
    return df[priority + middle + tail]


def write_target_audit_csv(result: TargetResult, output_dir: Path) -> Path:
    df = _reorder_columns(pd.DataFrame(result.ranked_rows))
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{result.target}_audit.csv"
    df.to_csv(path, index=False)
    return path


def write_final_selection_csv(results: list[TargetResult], output_dir: Path) -> Path:
    rows = []
    for result in results:
        target_rows = []
        for r in result.ranked_rows:
            if r["name"] in result.selected_names:
                row = dict(r)
                row["target"] = result.target
                row["target_short"] = result.short
                target_rows.append(row)
        rows.extend(target_rows)

        if result.short:
            logger.warning(
                "Target '%s' has only %d non-FAIL candidate(s) (< requested top-N); "
                "all were selected.",
                result.target,
                len(result.selected_names),
            )

        # `short` (not enough candidates at all) and this (enough candidates,
        # but not enough off-target-clean ones) are independent conditions --
        # a target can trigger either, both, or neither.
        backfilled = sum(1 for row in target_rows if row.get("exon_hits", 0) > 0)
        if backfilled:
            logger.warning(
                "Target '%s' didn't have enough off-target-clean candidates; "
                "%d selection(s) include a real (exon) off-target hit.",
                result.target,
                backfilled,
            )

    df = _reorder_columns(pd.DataFrame(rows))
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "final_selection.csv"
    df.to_csv(path, index=False)
    return path


def write_fasta_outputs(results: list[TargetResult], output_dir: Path) -> list[Path]:
    fasta_dir = output_dir / "fasta"
    written = []
    combined: list[tuple[str, str]] = []

    for result in results:
        records = [
            (r["name"], str(r["sequence"]).upper())
            for r in result.ranked_rows
            if r["name"] in result.selected_names
        ]
        path = fasta_dir / f"{result.target}.fasta"
        write_fasta(path, records)
        written.append(path)
        combined.extend(records)

    combined_path = fasta_dir / "all_targets.fasta"
    write_fasta(combined_path, combined)
    written.append(combined_path)
    return written
