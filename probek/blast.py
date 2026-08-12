from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .exceptions import MissingToolError
from .models import BlastHit

BLAST_TASK = "blastn-short"
WORD_SIZE = 7
EVALUE = 1000
OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "qlen",
]

BLAST_PARAMS = {
    "task": BLAST_TASK,
    "word_size": WORD_SIZE,
    "evalue": EVALUE,
    "outfmt_fields": OUTFMT_FIELDS,
    "dust": "no",
}

REQUIRED_TOOLS = ["blastn", "makeblastdb", "update_blastdb.pl"]


def check_blast_tools() -> None:
    missing = [tool for tool in REQUIRED_TOOLS if shutil.which(tool) is None]
    if missing:
        raise MissingToolError(
            "Missing required BLAST+ tool(s): "
            + ", ".join(missing)
            + "\n\nInstall with either:\n"
            "  sudo apt install ncbi-blast+\n"
            "or (conda/mamba):\n"
            "  conda install -c bioconda blast"
        )


def blastn_version() -> str:
    out = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
    return out.stdout.strip().splitlines()[0]


def write_batch_fasta(sequences: list[str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for seq in sequences:
            f.write(f">{seq}\n{seq}\n")  # sequence itself is the cache key -> use as FASTA id
    return path


def run_blastn(fasta_path: Path, db_path: Path, out_path: Path) -> Path:
    cmd = [
        "blastn",
        "-task",
        BLAST_TASK,
        "-word_size",
        str(WORD_SIZE),
        "-evalue",
        str(EVALUE),
        "-dust",
        "no",
        "-query",
        str(fasta_path),
        "-db",
        str(db_path),
        "-outfmt",
        "6 " + " ".join(OUTFMT_FIELDS),
        "-out",
        str(out_path),
    ]
    subprocess.run(cmd, check=True)
    return out_path


def parse_outfmt6(path: Path) -> dict[str, list[BlastHit]]:
    """Parse -outfmt 6 tabular output into hits grouped by query sequence (qseqid)."""
    hits_by_seq: dict[str, list[BlastHit]] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            row = dict(zip(OUTFMT_FIELDS, fields))
            hit = BlastHit(
                sseqid=row["sseqid"],
                pident=float(row["pident"]),
                length=int(row["length"]),
                qstart=int(row["qstart"]),
                qend=int(row["qend"]),
                sstart=int(row["sstart"]),
                send=int(row["send"]),
                evalue=float(row["evalue"]),
                qlen=int(row["qlen"]),
            )
            hits_by_seq.setdefault(row["qseqid"], []).append(hit)
    return hits_by_seq


def filter_hits(
    hits: list[BlastHit], min_identity: float, min_coverage: float
) -> list[BlastHit]:
    return [
        h
        for h in hits
        if h.pident >= min_identity and h.coverage * 100 >= min_coverage
    ]
