from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import annotation, blast, cache, chrom_map, classify, io_utils, reference_data, report, scoring
from .config import DEFAULT_BUILD, GENOME_BUILDS, AppConfig
from .exceptions import ProbekError
from .logging_setup import setup_logging
from .startup_checks import DEFAULT_INPUT_DIR, check_input_dir

logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="probek",
        description=(
            "Rank HERV-K FISH probe candidates by genomic off-target risk "
            "(exon / intron / HERV-K-family) beyond eFISHent's transcriptome-only screen."
        ),
    )
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT_DIR,
        help=f"Directory or file of eFISHent-format CSVs (default: ./{DEFAULT_INPUT_DIR}/)",
    )
    parser.add_argument("--output", default="results", help="Output directory (default: ./results/)")
    parser.add_argument("--top-n", type=int, default=12, help="Top candidates to select per target (default: 12)")
    parser.add_argument("--min-identity", type=float, default=90.0, help="Minimum %% identity for a genuine hit (default: 90)")
    parser.add_argument("--min-coverage", type=float, default=90.0, help="Minimum %% query coverage for a genuine hit (default: 90)")
    parser.add_argument("--reference-dir", default="reference_data", help="Reference data directory (default: ./reference_data/)")
    parser.add_argument("--build", default=DEFAULT_BUILD, choices=list(GENOME_BUILDS), help=f"Genome build (default: {DEFAULT_BUILD})")
    parser.add_argument("--non-interactive", action="store_true", help="Never prompt; fail clearly if reference data is missing")
    parser.add_argument("--force-reblast", action="store_true", help="Bypass the sequence cache and re-BLAST everything")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args(argv)
    args.input_was_explicit = any(a == "--input" or a.startswith("--input=") for a in (argv or sys.argv[1:]))
    return args


def build_config(args: argparse.Namespace) -> AppConfig:
    return AppConfig(
        input_dir=args.input,
        output_dir=args.output,
        reference_dir=args.reference_dir,
        build=args.build,
        top_n=args.top_n,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        non_interactive=args.non_interactive,
        force_reblast=args.force_reblast,
        input_was_explicit=args.input_was_explicit,
    )


def run(config: AppConfig) -> int:
    reference_dir = Path(config.reference_dir)
    output_dir = Path(config.output_dir)
    input_path = Path(config.input_dir)
    build = config.genome_build

    blast.check_blast_tools()
    reference_data.check_reference_data(reference_dir, build, non_interactive=config.non_interactive)
    if not check_input_dir(input_path, explicit=config.input_was_explicit):
        return 0

    frames = io_utils.load_input_csvs(input_path)
    frames = {target: io_utils.drop_failed(df) for target, df in frames.items()}
    groups = io_utils.dedup_sequences(frames)

    conn = cache.open_cache(reference_dir / "probek_cache.sqlite3")
    ref_version = reference_data.current_version(reference_dir)
    params_hash = cache.blast_params_hash(blast.BLAST_PARAMS)

    to_blast, cached_hits = cache.partition_by_cache(
        list(groups.keys()), conn, ref_version, params_hash, force=config.force_reblast
    )

    if to_blast:
        logger.info("BLASTing %d unique new sequence(s)...", len(to_blast))
        work_dir = output_dir / "_work"
        work_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = blast.write_batch_fasta(to_blast, work_dir / "batch_query.fasta")
        raw_out = blast.run_blastn(
            fasta_path, reference_data.blastdb_path(reference_dir, build), work_dir / "batch_hits.tsv"
        )
        new_hits = blast.parse_outfmt6(raw_out)
        blastn_version = blast.blastn_version()
        for seq in to_blast:
            hits = new_hits.get(seq, [])
            cache.store_hits(conn, seq, len(seq), hits, ref_version, params_hash, blastn_version)
        cached_hits.update({seq: new_hits.get(seq, []) for seq in to_blast})
    else:
        logger.info("All %d unique sequence(s) served from cache.", len(groups))

    cmap = chrom_map.ChromMap.from_reference_dir(reference_dir, build)
    ervk_index, exon_index, gene_index = annotation.load_indexes(reference_dir, build, cmap)

    classified_by_sequence = {}
    for seq, hits in cached_hits.items():
        genuine = blast.filter_hits(hits, config.min_identity, config.min_coverage)
        classified_by_sequence[seq] = classify.classify_hits(
            genuine, ervk_index, exon_index, gene_index
        )

    results = []
    for target, df in frames.items():
        result = scoring.rank_and_select(target, df, classified_by_sequence, config.top_n)
        report.write_target_audit_csv(result, classified_by_sequence, output_dir)
        results.append(result)

    report.write_final_selection_csv(results, output_dir)
    report.write_fasta_outputs(results, output_dir)

    logger.info("Done. Results written to %s", output_dir)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    setup_logging(args.verbose)
    config = build_config(args)
    try:
        return run(config)
    except ProbekError as e:
        logger.error(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
