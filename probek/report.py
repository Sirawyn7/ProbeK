from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .io_utils import write_fasta
from .models import ClassifiedHit, TargetResult

logger = logging.getLogger(__name__)


def format_off_target_loci(classified: list[ClassifiedHit]) -> str:
    parts = []
    for c in classified:
        start, end = c.hit.span
        label = f":{c.locus_label}" if c.locus_label else ""
        parts.append(f"{c.hit.sseqid}:{start}-{end}(Tier{c.tier}{label})")
    return "; ".join(parts)


def write_target_audit_csv(
    result: TargetResult,
    classified_by_sequence: dict[str, list[ClassifiedHit]],
    output_dir: Path,
) -> Path:
    rows = []
    for r in result.ranked_rows:
        row = dict(r)
        seq = str(r["sequence"]).upper()
        row["off_target_loci"] = format_off_target_loci(classified_by_sequence.get(seq, []))
        row["selected"] = r["name"] in result.selected_names
        rows.append(row)

    df = pd.DataFrame(rows)
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
