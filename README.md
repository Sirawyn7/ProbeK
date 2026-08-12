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
- Python 3.10+ installed and on `PATH` (everything else — BLAST+, Python packages, reference data — ProbeK sets up for you).
- No NCBI account, API key, or Entrez email needed — BLAST runs entirely locally.

## Quick start

1. Extract or clone this folder anywhere.
2. Drop your eFISHent-format CSVs into [`input_csvs/`](input_csvs/) (see that directory's README for the expected format and how target labels are derived).
3. Run:
   ```bash
   bash run.sh
   ```
   (`bash run.sh` always works. `./run.sh` also works, but only if the file's
   executable bit survived getting onto your machine — GitHub's "Download
   ZIP" button, for one, does not preserve it, which shows up as
   `Permission denied`. If you'd rather use `./run.sh`, fix that once with
   `chmod +x run.sh`.)

That's it. The first time you run it, ProbeK sets itself up step by step,
asking before it does anything (nothing happens silently):

1. **Local Python environment** — `run.sh` creates a `.venv/` folder here and installs ProbeK's Python dependencies into it. Nothing is installed system-wide.
2. **BLAST+ command-line tools** — if not already on your system, ProbeK offers to download NCBI's official prebuilt Linux build straight into this project's `reference_data/tools/` folder. No `sudo`/admin access needed.
3. **Reference data** (~1 GB+) — NCBI's genome BLAST database (GRCh38; note
   this tracks whatever patch level NCBI currently ships under the
   `human_genome` name, which as of writing is GRCh38.p13/GCF_000001405.39 —
   this doesn't affect correctness, since patch releases only add
   alternate/patch scaffolds and never move primary chromosome coordinates),
   the matching RefSeq gene annotation (pinned to GRCh38.p14/GCF_000001405.40),
   and the UCSC RepeatMasker track used to identify HERV-K loci.

Every run after the first skips straight to processing your CSVs — only
missing pieces trigger a prompt. If `input_csvs/` is empty, ProbeK just
prints a reminder of where to put your files and exits cleanly.

Pass flags the same way: `bash run.sh --top-n 15`.

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
| `--non-interactive` | off | Never prompt; fail clearly if reference data (or the local Python environment) is missing |
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

<details>
<summary><strong>Manual setup (for development, or if you'd rather manage the Python environment yourself)</strong></summary>

`run.sh` is just a convenience wrapper. You can set up and run ProbeK by hand instead:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
# optional, for faster interval-overlap classification — ProbeK falls back to
# a built-in pandas/numpy implementation automatically if this isn't available
# on your system:
pip install -e .[pyranges]    # or: pip install -e .[pyranges1]
```

Then run `probek` directly (it takes the same flags as `./run.sh`):

```bash
probek --input input_csvs/
```

Use `--non-interactive` in scripted/CI contexts — it fails with the exact
manual setup command instead of prompting, for both this Python environment
step and ProbeK's own BLAST+/reference-data setup.

### Development

```bash
pytest
```

Tests use synthetic fixtures only — no real BLAST execution or network
downloads. To force a specific interval-overlap backend when testing:

```bash
PROBEK_INTERVAL_BACKEND=custom pytest
```

</details>
