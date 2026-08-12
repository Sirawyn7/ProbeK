"""Loads gene/exon intervals (RefSeq GFF3) and HERV-K/HML-2 loci (UCSC rmsk,
filtered to repFamily == 'ERVK') into IntervalIndex objects.

This is the only module besides interval_backend.py that constructs indexes;
classify.py only ever calls IntervalIndex methods on what this module returns.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .chrom_map import ChromMap
from .config import GenomeBuild
from .interval_backend import IntervalFeature, IntervalIndex, build_index

GFF_COLUMNS = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]

RMSK_COLUMNS = [
    "bin",
    "swScore",
    "milliDiv",
    "milliDel",
    "milliIns",
    "genoName",
    "genoStart",
    "genoEnd",
    "genoLeft",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repStart",
    "repEnd",
    "repLeft",
    "id",
]

_GENE_ATTR_RE = re.compile(r"(?:^|;)gene=([^;]+)")


def _extract_gene(attributes: str) -> str | None:
    m = _GENE_ATTR_RE.search(attributes)
    return m.group(1) if m else None


def _gff_cache_path(reference_dir: Path, build: GenomeBuild) -> Path:
    return Path(reference_dir) / "annotation" / "exon_gene_cache.parquet"


def load_gff_features(gff_path: Path, reference_dir: Path, build: GenomeBuild) -> pd.DataFrame:
    """Parse exon/gene rows out of the RefSeq GFF3, caching the (much smaller)
    parsed result as parquet so repeat runs skip re-parsing the full genome GFF."""
    cache_path = _gff_cache_path(reference_dir, build)
    if cache_path.exists() and cache_path.stat().st_mtime >= gff_path.stat().st_mtime:
        return pd.read_parquet(cache_path)

    chunks = []
    for chunk in pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=GFF_COLUMNS,
        dtype={"start": "int64", "end": "int64"},
        chunksize=200_000,
        compression="infer",
    ):
        chunk = chunk[chunk["type"].isin(("exon", "gene"))]
        if len(chunk):
            chunk = chunk.assign(gene=chunk["attributes"].map(_extract_gene))
            chunks.append(chunk[["seqid", "type", "start", "end", "gene"]])

    features = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(
        columns=["seqid", "type", "start", "end", "gene"]
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    features.to_parquet(cache_path)
    return features


def load_rmsk_ervk(rmsk_path: Path, chrom_map: ChromMap) -> list[IntervalFeature]:
    """Load the UCSC rmsk track, filtered to repFamily == 'ERVK', translated
    from UCSC chrN naming to RefSeq accessions via chrom_map."""
    df = pd.read_csv(
        rmsk_path,
        sep="\t",
        header=None,
        names=RMSK_COLUMNS,
        compression="infer",
        usecols=["genoName", "genoStart", "genoEnd", "repName", "repFamily"],
    )
    df = df[df["repFamily"] == "ERVK"]

    features = []
    for row in df.itertuples(index=False):
        refseq = chrom_map.to_refseq(row.genoName)
        if refseq is None:
            continue
        # UCSC genoStart is 0-based; BLAST/GFF3 coordinates used elsewhere are
        # 1-based, so shift by one to keep all sources in the same frame.
        features.append(
            IntervalFeature(chrom=refseq, start=row.genoStart + 1, end=row.genoEnd, name=row.repName)
        )
    return features


def build_exon_index(features_df: pd.DataFrame) -> IntervalIndex:
    exons = features_df[features_df["type"] == "exon"]
    feats = [
        IntervalFeature(chrom=r.seqid, start=r.start, end=r.end, name=r.gene)
        for r in exons.itertuples(index=False)
    ]
    return build_index(feats)


def build_gene_index(features_df: pd.DataFrame) -> IntervalIndex:
    genes = features_df[features_df["type"] == "gene"]
    feats = [
        IntervalFeature(chrom=r.seqid, start=r.start, end=r.end, name=r.gene)
        for r in genes.itertuples(index=False)
    ]
    return build_index(feats)


def load_indexes(
    reference_dir: Path, build: GenomeBuild, chrom_map: ChromMap
) -> tuple[IntervalIndex, IntervalIndex, IntervalIndex]:
    """Returns (ervk_index, exon_index, gene_index)."""
    gff_path = Path(reference_dir) / "annotation" / f"{build.assembly_dir_name}_genomic.gff.gz"
    rmsk_path = Path(reference_dir) / "rmsk" / "rmsk.txt.gz"

    features_df = load_gff_features(gff_path, reference_dir, build)
    exon_index = build_exon_index(features_df)
    gene_index = build_gene_index(features_df)
    ervk_index = build_index(load_rmsk_ervk(rmsk_path, chrom_map))

    return ervk_index, exon_index, gene_index
