"""Static configuration: expected eFISHent columns and genome-build registry.

Genome build is kept as a lookup table (not hardcoded inline in reference_data.py
or chrom_map.py) so a future non-human build only requires a new entry here.
"""

from __future__ import annotations

from dataclasses import dataclass

# Columns eFISHent's per-gene CSV output is expected to contain.
EFISHENT_COLUMNS = [
    "name",
    "sequence",
    "length",
    "start",
    "end",
    "GC",
    "TM",
    "deltaG",
    "kmers",
    "count",
    "cpg_fraction",
    "low_complexity",
    "accessibility",
    "on_target_dg",
    "txome_off_targets",
    "off_target_genes",
    "worst_match",
    "expression_risk",
    "quality",
    "recommendation",
]


@dataclass(frozen=True)
class GenomeBuild:
    name: str
    blast_db_name: str  # name passed to update_blastdb.pl
    assembly_accession: str  # e.g. GCF_000001405.40
    assembly_dir_name: str  # e.g. GCF_000001405.40_GRCh38.p14
    gff_url: str
    assembly_report_url: str
    rmsk_url: str


GENOME_BUILDS: dict[str, GenomeBuild] = {
    "GRCh38": GenomeBuild(
        name="GRCh38",
        blast_db_name="human_genome",
        assembly_accession="GCF_000001405.40",
        assembly_dir_name="GCF_000001405.40_GRCh38.p14",
        gff_url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
        ),
        assembly_report_url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
        ),
        rmsk_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
    ),
}

DEFAULT_BUILD = "GRCh38"


@dataclass
class AppConfig:
    input_dir: str
    output_dir: str
    reference_dir: str
    build: str
    top_n: int
    min_identity: float
    min_coverage: float
    non_interactive: bool
    force_reblast: bool
    input_was_explicit: bool

    @property
    def genome_build(self) -> GenomeBuild:
        return GENOME_BUILDS[self.build]
