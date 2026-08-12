from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

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


def feature_counts(classified: list[ClassifiedHit]) -> dict[str, int]:
    counts = {"hervk": 0, "exon": 0, "intron": 0, "outside_gene": 0}
    for c in classified:
        if c.tier == "A":
            counts["hervk"] += 1
        elif c.tier == "C":
            counts["exon"] += 1
        elif c.locus_label:
            counts["intron"] += 1
        else:
            counts["outside_gene"] += 1
    return counts


def enrich_target_result(
    result: TargetResult, classified_by_sequence: dict[str, list[ClassifiedHit]]
) -> None:
    """Replaces each row's internal Tier A/B/C sort-key counts (used only for
    ranking, in probek.scoring's vocabulary) with plain-English off-target
    detail: a per-locus feature+location breakdown, and separate HERV-K
    family / exon / intron / outside-gene hit counts. Mutates `result` in
    place so both the per-target audit CSV and the combined final-selection
    CSV pick up the same enriched columns."""
    for row in result.ranked_rows:
        seq = str(row["sequence"]).upper()
        classified = classified_by_sequence.get(seq, [])
        row["off_target_loci"] = format_off_target_loci(classified)
        counts = feature_counts(classified)
        row["hervk_family_hits"] = counts["hervk"]
        row["exon_hits"] = counts["exon"]
        row["intron_hits"] = counts["intron"]
        row["outside_gene_hits"] = counts["outside_gene"]
        for internal_key in ("tier_a_count", "tier_b_count", "tier_c_count"):
            row.pop(internal_key, None)
        row["selected"] = row["name"] in result.selected_names


def write_target_audit_csv(result: TargetResult, output_dir: Path) -> Path:
    df = pd.DataFrame(result.ranked_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{result.target}_audit.csv"
    df.to_csv(path, index=False)
    return path


def write_final_selection_csv(results: list[TargetResult], output_dir: Path) -> Path:
    rows = []
    for result in results:
        for r in result.ranked_rows:
            if r["name"] in result.selected_names:
                row = dict(r)
                row["target"] = result.target
                row["target_short"] = result.short
                rows.append(row)
        if result.short:
            logger.warning(
                "Target '%s' has only %d non-FAIL candidate(s) (< requested top-N); "
                "all were selected.",
                result.target,
                len(result.selected_names),
            )

    df = pd.DataFrame(rows)
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
