"""First-run interactive reference data setup: NCBI genome BLAST db, RefSeq
gene annotation, UCSC RepeatMasker ERVK track, and the NCBI assembly report
used for chromosome-name reconciliation.

Downloads are never silent — each missing piece is confirmed with the user
(y/n) unless --non-interactive, in which case a missing piece is a hard
failure with the exact manual command to run.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import requests

from .config import GenomeBuild
from .exceptions import ReferenceDataMissingError
from .prompts import prompt_yes_no

MANIFEST_NAME = "manifest.json"


def _manifest_path(reference_dir: Path) -> Path:
    return Path(reference_dir) / MANIFEST_NAME


def load_manifest(reference_dir: Path) -> dict:
    path = _manifest_path(reference_dir)
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def save_manifest(reference_dir: Path, manifest: dict) -> None:
    path = _manifest_path(reference_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def blastdb_dir(reference_dir: Path) -> Path:
    return Path(reference_dir) / "blastdb"


def _discover_blastdb_basename(d: Path) -> str | None:
    """update_blastdb.pl's target name (e.g. "human_genome") is just the
    identifier you ask it to fetch — the files it actually writes are named
    after the underlying assembly (e.g. "GCF_000001405.39_top_level.*"), so
    the real queryable -db basename has to be discovered on disk rather than
    assumed to match the identifier we asked for."""
    if not d.exists():
        return None
    nal_files = sorted(d.glob("*.nal"))
    if nal_files:
        return nal_files[0].stem
    # Single-volume db has no .nal alias; a multi-volume one's .NN.nin files
    # would falsely match here too, so exclude those.
    nin_files = sorted(f for f in d.glob("*.nin") if not re.search(r"\.\d+\.nin$", f.name))
    if nin_files:
        return nin_files[0].stem
    return None


def blastdb_path(reference_dir: Path, build: GenomeBuild) -> Path:
    d = blastdb_dir(reference_dir)
    basename = _discover_blastdb_basename(d)
    if basename is None:
        raise ReferenceDataMissingError(f"No usable BLAST database found in {d}.")
    return d / basename


def gff_path(reference_dir: Path, build: GenomeBuild) -> Path:
    return Path(reference_dir) / "annotation" / f"{build.assembly_dir_name}_genomic.gff.gz"


def rmsk_path(reference_dir: Path) -> Path:
    return Path(reference_dir) / "rmsk" / "rmsk.txt.gz"


def assembly_report_path(reference_dir: Path, build: GenomeBuild) -> Path:
    return (
        Path(reference_dir)
        / "assembly_report"
        / f"{build.assembly_dir_name}_assembly_report.txt"
    )


def blastdb_present(reference_dir: Path, build: GenomeBuild) -> bool:
    return _discover_blastdb_basename(blastdb_dir(reference_dir)) is not None


def annotation_present(reference_dir: Path, build: GenomeBuild) -> bool:
    return gff_path(reference_dir, build).exists()


def rmsk_present(reference_dir: Path, build: GenomeBuild) -> bool:
    return rmsk_path(reference_dir).exists()


def assembly_report_present(reference_dir: Path, build: GenomeBuild) -> bool:
    return assembly_report_path(reference_dir, build).exists()


def _download_file(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)


def download_blastdb(reference_dir: Path, build: GenomeBuild) -> None:
    d = blastdb_dir(reference_dir)
    d.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["update_blastdb.pl", build.blast_db_name, "--decompress"],
        cwd=d,
        check=True,
    )


def download_annotation(reference_dir: Path, build: GenomeBuild) -> None:
    _download_file(build.gff_url, gff_path(reference_dir, build))


def download_rmsk(reference_dir: Path, build: GenomeBuild) -> None:
    _download_file(build.rmsk_url, rmsk_path(reference_dir))


def download_assembly_report(reference_dir: Path, build: GenomeBuild) -> None:
    _download_file(build.assembly_report_url, assembly_report_path(reference_dir, build))


def check_reference_data(
    reference_dir: Path, build: GenomeBuild, non_interactive: bool = False
) -> None:
    """Ensures all four reference-data pieces are present, prompting for any
    that are missing (or failing clearly in --non-interactive mode)."""
    reference_dir = Path(reference_dir)
    checks = [
        (
            "genome BLAST database",
            blastdb_present,
            download_blastdb,
            f"update_blastdb.pl {build.blast_db_name} --decompress   (run inside {blastdb_dir(reference_dir)})",
        ),
        (
            "RefSeq gene annotation (GFF3)",
            annotation_present,
            download_annotation,
            f"curl -o {gff_path(reference_dir, build)} {build.gff_url}",
        ),
        (
            "UCSC RepeatMasker (rmsk) track",
            rmsk_present,
            lambda rd, b: download_rmsk(rd, b),
            f"curl -o {rmsk_path(reference_dir)} {build.rmsk_url}",
        ),
        (
            "NCBI assembly report",
            assembly_report_present,
            download_assembly_report,
            f"curl -o {assembly_report_path(reference_dir, build)} {build.assembly_report_url}",
        ),
    ]

    for label, present_fn, download_fn, manual_cmd in checks:
        if present_fn(reference_dir, build):
            continue
        if non_interactive:
            raise ReferenceDataMissingError(
                f"Missing reference data: {label}.\nRun manually:\n  {manual_cmd}"
            )
        print(f"\nReference data missing: {label}")
        if not prompt_yes_no(f"Download it now to {reference_dir}?"):
            raise ReferenceDataMissingError(
                f"Reference data required: {label}.\nRun manually:\n  {manual_cmd}"
            )
        download_fn(reference_dir, build)

    manifest = load_manifest(reference_dir)
    manifest.setdefault("build", build.name)
    manifest.setdefault("blastdb_downloaded_at", datetime.now(timezone.utc).isoformat())
    save_manifest(reference_dir, manifest)


def current_version(reference_dir: Path) -> str:
    """Fingerprint used to detect stale BLAST cache entries when reference
    data is later re-downloaded/updated."""
    manifest = load_manifest(Path(reference_dir))
    key = f"{manifest.get('build', '')}:{manifest.get('blastdb_downloaded_at', '')}"
    return hashlib.sha256(key.encode("utf-8")).hexdigest()
