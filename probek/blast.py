from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

import requests

from . import download_utils
from .exceptions import MissingToolError
from .models import BlastHit
from .prompts import prompt_yes_no

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

BLAST_RELEASE_INDEX_URL = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"

_MANUAL_INSTALL_HINT = (
    "Install manually:\n"
    "  sudo apt install ncbi-blast+\n"
    "or (conda/mamba):\n"
    "  conda install -c bioconda blast"
)


def check_blast_tools() -> list[str]:
    """Returns the subset of REQUIRED_TOOLS not currently found on PATH."""
    return [tool for tool in REQUIRED_TOOLS if shutil.which(tool) is None]


def _local_bin_dir(reference_dir: Path) -> Path:
    return Path(reference_dir) / "tools" / "blast+" / "bin"


def _tools_present_in(bin_dir: Path) -> bool:
    return bin_dir.is_dir() and all((bin_dir / tool).exists() for tool in REQUIRED_TOOLS)


def _prepend_to_path(bin_dir: Path) -> None:
    # Must be absolute: some subprocess calls elsewhere (e.g. update_blastdb.pl,
    # invoked with cwd=reference_data/blastdb/) run with a different working
    # directory, and a relative PATH entry would resolve against *that* cwd
    # instead of the project root, making the tool unfindable.
    os.environ["PATH"] = str(Path(bin_dir).resolve()) + os.pathsep + os.environ.get("PATH", "")


def _find_latest_linux_tarball_url() -> tuple[str, str]:
    """Returns (tarball_url, md5_url) for NCBI's latest prebuilt Linux x64 BLAST+ build."""
    resp = requests.get(BLAST_RELEASE_INDEX_URL, timeout=30)
    resp.raise_for_status()
    match = re.search(r'href="(ncbi-blast-[\d.]+\+-x64-linux\.tar\.gz)"', resp.text)
    if not match:
        raise MissingToolError(
            "Could not automatically determine the latest BLAST+ release for "
            "this platform.\n\n" + _MANUAL_INSTALL_HINT
        )
    filename = match.group(1)
    return BLAST_RELEASE_INDEX_URL + filename, BLAST_RELEASE_INDEX_URL + filename + ".md5"


def install_blast_plus_locally(reference_dir: Path) -> Path:
    """Downloads and extracts NCBI's official prebuilt BLAST+ build for Linux
    x64 into reference_dir/tools/blast+/ — no admin/sudo privileges required.
    Verifies the download against NCBI's published MD5 checksum. Returns the
    resulting bin/ directory."""
    tarball_url, md5_url = _find_latest_linux_tarball_url()

    tools_dir = Path(reference_dir) / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    tarball_path = tools_dir / "blast_plus_download.tar.gz"

    download_utils.download_with_progress(tarball_url, tarball_path, desc="BLAST+ tools")
    download_utils.verify_md5(tarball_path, md5_url)

    with tarfile.open(tarball_path) as tar:
        try:
            tar.extractall(tools_dir, filter="data")
        except TypeError:
            tar.extractall(tools_dir)  # Python <3.12 without the `filter` argument
    tarball_path.unlink()

    extracted_dirs = sorted(
        d for d in tools_dir.iterdir() if d.is_dir() and d.name.startswith("ncbi-blast-")
    )
    if not extracted_dirs:
        raise MissingToolError(
            "BLAST+ download completed but the expected directory was not found "
            f"under {tools_dir}.\n\n" + _MANUAL_INSTALL_HINT
        )

    target_dir = _local_bin_dir(reference_dir).parent
    if target_dir.exists():
        shutil.rmtree(target_dir)
    extracted_dirs[-1].rename(target_dir)

    bin_dir = target_dir / "bin"
    if not _tools_present_in(bin_dir):
        raise MissingToolError(
            f"BLAST+ was installed to {target_dir} but the expected tools were "
            f"not found in its bin/ directory.\n\n" + _MANUAL_INSTALL_HINT
        )
    return bin_dir


def ensure_blast_tools(reference_dir: Path, non_interactive: bool = False) -> None:
    """Ensures blastn/makeblastdb/update_blastdb.pl are available on PATH for
    this process, with minimal manual setup:

    1. If already on PATH (system install via apt/conda/etc.), use that.
    2. Else, if previously auto-installed under reference_dir/tools/blast+/,
       reuse it (just prepend to PATH for this process).
    3. Else, offer to download NCBI's official prebuilt Linux build locally —
       no admin/sudo privileges required — and use that.

    In --non-interactive mode, a missing install is a hard failure with both
    the local-install pointer and the manual system-install commands.
    """
    if not check_blast_tools():
        return

    local_bin = _local_bin_dir(reference_dir)
    if _tools_present_in(local_bin):
        _prepend_to_path(local_bin)
        return

    if non_interactive:
        raise MissingToolError(
            "Missing required BLAST+ tool(s): "
            + ", ".join(check_blast_tools())
            + "\n\nRun without --non-interactive to have ProbeK install these "
            "automatically (no admin/sudo required), or " + _MANUAL_INSTALL_HINT
        )

    print(
        "\nBLAST+ command-line tools (blastn, makeblastdb, update_blastdb.pl) were not found."
    )
    if prompt_yes_no(
        f"Download and install them automatically now to {local_bin.parent}? "
        "(official NCBI build, no admin/sudo required)"
    ):
        bin_dir = install_blast_plus_locally(reference_dir)
        _prepend_to_path(bin_dir)
        return

    raise MissingToolError("BLAST+ is required to continue.\n\n" + _MANUAL_INSTALL_HINT)


def blastn_version() -> str:
    out = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
    return out.stdout.strip().splitlines()[0]


def _fasta_id_for_sequence(seq: str) -> str:
    return "seq_" + hashlib.md5(seq.encode("utf-8")).hexdigest()[:16]


def write_batch_fasta(sequences: list[str], path: Path) -> tuple[Path, dict[str, str]]:
    """Writes a multi-FASTA of `sequences`, each under a short synthetic ID
    rather than the raw sequence itself — using the sequence as its own FASTA
    header trips blastn's "sequence in title line?" sanity-check warning on
    every single query. Returns (fasta_path, id_to_sequence) so callers can
    map blastn's qseqid output back to the original sequence."""
    path.parent.mkdir(parents=True, exist_ok=True)
    id_to_sequence: dict[str, str] = {}
    with open(path, "w", encoding="utf-8") as f:
        for seq in sequences:
            seq_id = _fasta_id_for_sequence(seq)
            id_to_sequence[seq_id] = seq
            f.write(f">{seq_id}\n{seq}\n")
    return path, id_to_sequence


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


_SSEQID_REF_RE = re.compile(r"\|ref\|([^|]+)\|?")


def _normalize_sseqid(raw: str) -> str:
    """BLAST's outfmt 6 `sseqid` is either a bare accession (e.g.
    "NC_000001.11") or, for databases built with NCBI's older-style deflines
    (as the `human_genome` preformatted db currently is), the full pipe-
    delimited form "gi|1234|ref|NC_000001.11|". Everything downstream
    (GFF3 seqid, rmsk genoName translated via chrom_map, assembly report)
    is keyed on the bare accession, so it must be extracted here — otherwise
    every hit silently fails to match any annotation/ERVK interval."""
    match = _SSEQID_REF_RE.search(raw)
    return match.group(1) if match else raw


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
                sseqid=_normalize_sseqid(row["sseqid"]),
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
