# ProbeK

Genomic off-target ranking for HERV-K (HML-2) FISH probe candidates.

[eFISHent](https://github.com/BBQuercus/eFISHent) screens candidate smFISH
probes against the human *transcriptome*, but says nothing about *genomic*
off-targets. HERV-K/HML-2 is a ~90-copy repetitive family, so almost every
probe hits dozens of other HERV-K loci genome-wide — that's expected and
desirable. What eFISHent's screen misses is the distinction between:

- a hit landing in the **exon** of an unrelated gene (real cross-hybridization risk),
- a hit landing in an unrelated gene's **intron** or non-genic space (low risk), and
- a hit landing in another **HERV-K/HML-2 locus** (desired signal — more of this is better).

ProbeK automates the local-BLAST-plus-annotation workflow that previously
meant clicking through the NCBI BLAST web UI one probe at a time, and ranks
candidates accordingly.

## What it does

1. **Load & pre-filter** — reads all eFISHent CSVs in the input directory, drops rows already flagged `FAIL`, and deduplicates identical probe sequences across files/targets.
2. **BLAST** — runs `blastn -task blastn-short` for every unique sequence against a local copy of NCBI's `human_genome` database, batched into one call.
3. **Filter** — discards alignments below a configurable identity/coverage floor (default ≥90%/≥90%).
4. **Classify** — for each surviving hit, checks interval overlap against RepeatMasker ERVK loci, then exons, then falls back to intron/non-genic:
   - **Tier A** — HERV-K family locus (counts *in favor* of the probe)
   - **Tier B** — intron or non-genic (low risk, tolerated)
   - **Tier C** — exon of an unrelated gene (real off-target risk)
5. **Rank & select** — sorts by fewest Tier C, then most Tier A, then fewest Tier B, then eFISHent's `quality` score as a tiebreaker, and selects the top N (default 12) per target.
6. **Output** — per-target audit CSVs, a combined final-selection table, and vendor-ready FASTA files.

## Requirements

- WSL2 (Ubuntu) or native Linux.
- Python 3.10+.
- No NCBI account, API key, or Entrez email needed — BLAST runs entirely locally.

You do **not** need to install BLAST+ yourself — see below.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
# optional, for faster interval-overlap classification — ProbeK falls back to
# a built-in pandas/numpy implementation automatically if this isn't available
# on your system:
pip install -e .[pyranges]    # or: pip install -e .[pyranges1]
```

On first run, ProbeK checks for everything else it needs and offers to set it
up for you, prompting (y/n) before downloading anything — nothing happens
silently:

1. **BLAST+ command-line tools** (`blastn`, `makeblastdb`, `update_blastdb.pl`).
   If these aren't already on your system, ProbeK offers to download NCBI's
   official prebuilt Linux build straight into the project's
   `reference_data/tools/` folder — no `sudo`/admin access or package manager
   needed. (If you'd rather install it system-wide yourself:
   `sudo apt install ncbi-blast+`, or `conda install -c bioconda blast`.)
2. **NCBI's preformatted `human_genome` BLAST database** (~1 GB, GRCh38). Note:
   this is NCBI's own rolling alias for "the current human genome build," so
   its exact patch level can lag behind the GRCh38.p14 (GCF_000001405.40)
   pinned for the annotation below — as of writing it resolves to GRCh38.p13
   (GCF_000001405.39). This doesn't affect correctness: patch releases only
   add alternate/patch scaffolds and never move primary chromosome (chr1-22,
   X, Y, MT) coordinates, which is what classification actually relies on.
3. **The matching RefSeq gene annotation** (GFF3, pinned to GRCh38.p14 / GCF_000001405.40).
4. **The UCSC RepeatMasker track**, filtered to the ERVK family, plus the NCBI
   assembly report used to reconcile RefSeq/UCSC chromosome naming.

Use `--non-interactive` in scripted/CI contexts — it fails with the exact
manual setup command instead of prompting.

## Usage

Drop your eFISHent-format CSVs into [`input_csvs/`](input_csvs/) (see that
directory's README for the expected format and how target labels are
derived), then:

```bash
probek --input input_csvs/
```

If `input_csvs/` is empty, ProbeK prints a reminder of where to put your
files and exits cleanly.

### Flags

| Flag | Default | Description |
|---|---|---|
| `--input` | `./input_csvs/` | Directory or file of eFISHent CSVs |
| `--output` | `./results/` | Output directory |
| `--top-n` | `12` | Candidates to select per target |
| `--min-identity` | `90` | Minimum % identity for a genuine BLAST hit |
| `--min-coverage` | `90` | Minimum % query coverage for a genuine BLAST hit |
| `--reference-dir` | `./reference_data/` | Reference data location |
| `--build` | `GRCh38` | Genome build |
| `--non-interactive` | off | Never prompt; fail clearly if reference data is missing |
| `--force-reblast` | off | Bypass the sequence cache |
| `--verbose` | off | Verbose logging |

## Output

Written to `results/` (configurable via `--output`):

```
results/
├── <target>_audit.csv       # every surviving candidate, full tier breakdown + rank
├── final_selection.csv      # combined top-N selection across all targets
└── fasta/
    ├── <target>.fasta       # selected sequences for one target
    └── all_targets.fasta    # selected sequences, combined
```

## Development

```bash
pytest
```

Tests use synthetic fixtures only — no real BLAST execution or network
downloads. To force a specific interval-overlap backend when testing:

```bash
PROBEK_INTERVAL_BACKEND=custom pytest
```
