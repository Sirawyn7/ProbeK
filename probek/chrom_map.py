"""Reconciles chromosome naming between RefSeq accessions (NC_000001.11) and
UCSC-style names (chr1), using NCBI's own assembly report file rather than a
hand-maintained mapping table.

BLAST hits (sseqid) and the RefSeq GFF3 annotation already share the RefSeq
accession namespace, so no translation is needed between them. The only
source that needs translating is the UCSC RepeatMasker (rmsk) track, whose
genoName column uses chr1-style names.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import GenomeBuild


def parse_assembly_report(path: Path) -> pd.DataFrame:
    """Return columns [refseq_accn, ucsc_name] from an NCBI assembly report.

    Format: comment lines start with '#'; the LAST comment line is the real
    tab-delimited header. Rows where UCSC-style-name is 'na' (unplaced/unmapped
    scaffolds) are dropped since they have no UCSC counterpart to reconcile.
    """
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("#"):
            header_idx = i
        else:
            break
    if header_idx is None:
        raise ValueError(f"No header line found in assembly report: {path}")

    header = lines[header_idx].lstrip("#").strip().split("\t")
    data_lines = lines[header_idx + 1 :]

    from io import StringIO

    df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", names=header, dtype=str)
    df = df.rename(columns={"RefSeq-Accn": "refseq_accn", "UCSC-style-name": "ucsc_name"})
    df = df[["refseq_accn", "ucsc_name"]]
    df = df[df["ucsc_name"].notna() & (df["ucsc_name"] != "na")]
    return df.reset_index(drop=True)


class ChromMap:
    def __init__(self, df: pd.DataFrame):
        self._to_refseq = dict(zip(df["ucsc_name"], df["refseq_accn"]))
        self._to_ucsc = dict(zip(df["refseq_accn"], df["ucsc_name"]))

    def to_refseq(self, name: str) -> str | None:
        return self._to_refseq.get(name)

    def to_ucsc(self, name: str) -> str | None:
        return self._to_ucsc.get(name)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "ChromMap":
        return cls(df)

    @classmethod
    def from_reference_dir(cls, reference_dir: Path, build: GenomeBuild) -> "ChromMap":
        path = (
            Path(reference_dir)
            / "assembly_report"
            / f"{build.assembly_dir_name}_assembly_report.txt"
        )
        return cls(parse_assembly_report(path))
