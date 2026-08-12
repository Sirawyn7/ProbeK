from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .config import EFISHENT_COLUMNS
from .exceptions import InputDataError
from .models import ProbeOccurrence, SequenceGroup

_TRAILING_INDEX_RE = re.compile(r"-\d+$")


def derive_target_label(names: list[str], filename_stem: str) -> str:
    """Derive a target/gene label for a file.

    Primary rule: strip a trailing "-<number>" suffix from each `name` value
    (e.g. pNRV101_gag-1 -> pNRV101_gag) and use it if every row in the file
    agrees on the same stripped value. Falls back to the filename stem if
    names disagree or the pattern doesn't match.
    """
    if names:
        stripped = {_TRAILING_INDEX_RE.sub("", n) for n in names}
        if len(stripped) == 1:
            candidate = next(iter(stripped))
            if candidate:
                return candidate
    return filename_stem


def load_csv_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [c for c in EFISHENT_COLUMNS if c not in df.columns]
    if missing:
        raise InputDataError(
            f"{path}: missing expected eFISHent column(s): {', '.join(missing)}"
        )
    return df


def resolve_input_files(input_path: Path) -> list[Path]:
    if input_path.is_dir():
        files = sorted(input_path.glob("*.csv"))
    elif input_path.is_file():
        files = [input_path]
    else:
        raise InputDataError(f"Input path does not exist: {input_path}")
    if not files:
        raise InputDataError(f"No CSV files found in: {input_path}")
    return files


def load_input_csvs(input_path: Path) -> dict[str, pd.DataFrame]:
    """Load all eFISHent CSVs under input_path, keyed by derived target label.

    If two files derive the same target label, their rows are concatenated
    under that one target (this is legitimate — nothing requires one file per
    target).
    """
    frames: dict[str, list[pd.DataFrame]] = {}
    for path in resolve_input_files(input_path):
        df = load_csv_file(path)
        names = df["name"].astype(str).tolist()
        target = derive_target_label(names, path.stem)
        frames.setdefault(target, []).append(df)

    return {target: pd.concat(dfs, ignore_index=True) for target, dfs in frames.items()}


def drop_failed(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["recommendation"].astype(str).str.strip().str.upper() != "FAIL"].reset_index(
        drop=True
    )


def dedup_sequences(frames: dict[str, pd.DataFrame]) -> dict[str, SequenceGroup]:
    """Group rows across all targets by exact sequence string.

    The same sequence can legitimately appear under multiple targets/names
    (e.g. overlapping HERV-K reading frames), so a SequenceGroup tracks every
    occurrence rather than assuming one.
    """
    groups: dict[str, SequenceGroup] = {}
    for target, df in frames.items():
        for _, row in df.iterrows():
            seq = str(row["sequence"]).upper()
            occurrence = ProbeOccurrence(target=target, name=str(row["name"]), row=row.to_dict())
            if seq not in groups:
                groups[seq] = SequenceGroup(sequence=seq)
            groups[seq].occurrences.append(occurrence)
    return groups


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for name, sequence in records:
            f.write(f">{name}\n{sequence}\n")
