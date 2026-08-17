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
2. **BLAST** — runs `blastn -task blastn-short` for every unique sequence against a local copy of NCBI's `human_genome` database, batched into one call. Shows a live elapsed-time indicator while it runs — blastn itself reports no percentage, so this is "still working, Xs elapsed" rather than a true progress bar.
3. **Filter** — discards alignments below a configurable identity/coverage floor (default ≥90%/≥90%).
4. **Classify** — for each surviving hit, checks interval overlap against RepeatMasker ERVK loci, then exons, then falls back to intron/non-genic:
   - **HERV-K family locus** — desired signal, more of this is fine
   - **Intron or non-genic** — low risk, tolerated
   - **Exon of an unrelated gene** — real off-target risk; any probe with at least one of these is excluded from normal ranking (see below)
5. **Rank & select** — any probe with one or more off-target hits landing in an exon of an unrelated gene is excluded from normal selection (a real cross-hybridization risk, not something to rank around). Among the rest, probes are ranked by eFISHent's `quality` score (higher is better); ties are broken by HERV-K-family hit count (fewer is better — a probe cross-reacting with fewer other HERV-K/HML-2 loci is more specific to its own locus), and remaining ties by intron-hit count (fewer is better). Hits landing entirely outside any gene are informational only and never affect ranking. The top N (default 12) per target are selected this way; if a target doesn't have enough exon-hit-free candidates to fill that quota, the shortfall is backfilled from the excluded pool (least-bad — fewest exon hits — first, then the same tiebreak chain), and this is logged so it's never silent. See [Output](#output) for how to spot a backfilled selection.
6. **Output** — per-target audit CSVs, a combined final-selection table, and vendor-ready FASTA files.

## Requirements

- WSL2 (Ubuntu) or native Linux.
- Python 3.10+ installed and on `PATH` (everything else — BLAST+, Python packages, reference data — ProbeK sets up for you).
- No NCBI account, API key, or Entrez email needed — BLAST runs entirely locally.

## Quick start

1. Clone this folder anywhere:
   ```bash
   git clone https://github.com/Sirawyn7/ProbeK.git
   ```
   (You can also grab it via GitHub's "Download ZIP" button instead, but a
   `git clone` is what lets [`bash update.sh`](#updating) pull in future
   updates later. See [Updating](#updating) below if you start from a ZIP
   and want to switch.)
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
2. **BLAST+ command-line tools** — if not already on your system, ProbeK offers to download NCBI's official prebuilt Linux build straight into this project's `reference_data/tools/` folder, with a live progress bar. No `sudo`/admin access needed.
3. **Reference data** (~1 GB+) — NCBI's genome BLAST database (GRCh38; note
   this tracks whatever patch level NCBI currently ships under the
   `human_genome` name, which as of writing is GRCh38.p13/GCF_000001405.39 —
   this doesn't affect correctness, since patch releases only add
   alternate/patch scaffolds and never move primary chromosome coordinates),
   the matching RefSeq gene annotation (pinned to GRCh38.p14/GCF_000001405.40),
   and the UCSC RepeatMasker track used to identify HERV-K loci.

Every download — including the ~1 GB genome database — shows a live progress
bar tracking bytes downloaded against the exact expected total, so you can
always tell it's actually working rather than hung. Every run after the
first skips straight to processing your CSVs — only missing pieces trigger a
prompt. If `input_csvs/` is empty, ProbeK just prints a reminder of where to
put your files and exits cleanly.

Pass flags the same way: `bash run.sh --top-n 15`.

## Updating

```bash
bash update.sh
```

Checks GitHub for a newer version, and if one exists, asks
`Update now? [y/N]:` before doing anything. On yes, it pulls the latest code
and reinstalls dependencies if they changed. It never touches your
`input_csvs/`, `reference_data/`, or `results/` — those aren't tracked by
git, so they're untouched by a pull either way.

This needs the folder to have been set up via `git clone` (see
[Quick start](#quick-start)) rather than a ZIP download. If you started from
a ZIP, `bash update.sh` will detect that and print steps to switch to a git
clone — in short: `git clone` a fresh copy elsewhere, then copy your
`input_csvs/*.csv` and `reference_data/` folder into it so you don't have to
re-download the ~1 GB+ reference data.

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
├── <target>_audit.csv       # every surviving candidate, full off-target breakdown + rank
├── final_selection.csv      # combined top-N selection across all targets (self-contained, same columns)
└── fasta/
    ├── <target>.fasta       # selected sequences for one target
    └── all_targets.fasta    # selected sequences, combined
```

Both CSVs share the same columns and column order — `final_selection.csv` is
self-contained, so you don't need to cross-reference the per-target audit
files to see why a probe was chosen. The glanceable columns come first; the
original eFISHent columns and the full off-target detail trail behind. See
[CSV column reference](#csv-column-reference) below for what every column means.

### CSV column reference

#### Columns ProbeK adds

| Column | Meaning |
|---|---|
| `target` | Which gene target this probe belongs to (e.g. `pNRV101_gag`). Only in `final_selection.csv` — each per-target audit CSV is already just one target. |
| `selected` | Whether this probe is one of the chosen top-N for its target. |
| `rank` | This probe's rank within its target, 1 = best. See sort order below. |
| `target_short` | `True` if this target had fewer than `--top-n` non-FAIL candidates *in total* (counting exon-hit candidates too) — meaning every available candidate was selected. This is independent of whether any selected probe was backfilled from the excluded (exon-hit) pool — check `exon_hits > 0` on a `selected` row for that; both conditions are logged separately when they happen. Only in `final_selection.csv`. |
| `off_target_risk` | At-a-glance category — `Low` / `Moderate` / `High` — based purely on `exon_hits`, the only off-target category that represents real cross-hybridization risk: 0 → Low, 1–2 → Moderate, 3+ → High. Any row above `Low` is normally excluded from selection; it only appears as `selected` if it was backfilled (see [Rank & select](#what-it-does)). |
| `flagged_genes` | Comma-separated gene symbols that had a real (exon) off-target hit, e.g. `GM2A`. Empty if none. The short version of `off_target_loci` for exon hits specifically. |
| `hervk_family_hits` | Off-target hits landing in another annotated HERV-K/HML-2 (RepeatMasker ERVK) locus. This is the *desired* signal for a pan-HERV-K-family probe. It's also this ranking's tiebreaker when two probes have the exact same `quality` — fewer wins, since a probe cross-reacting with fewer other HERV-K/HML-2 loci is more specific to its own locus. |
| `exon_hits` | Off-target hits landing in an exon of an unrelated gene. Real cross-hybridization risk — any probe with `exon_hits > 0` is excluded from normal ranking/selection. It's only selected as a logged fallback if a target doesn't have enough exon-hit-free candidates to fill `--top-n` otherwise. |
| `intron_hits` | Off-target hits landing in an intron of an unrelated gene. Low risk, tolerated. Also this ranking's final tiebreaker, engaged only when `quality` and `hervk_family_hits` are both exactly tied — fewer wins. |
| `outside_gene_hits` | Off-target hits landing outside any annotated gene entirely (intergenic). Low risk, tolerated, and ignored entirely for ranking purposes. |
| `off_target_loci` | The full audit trail: every individual off-target hit, one entry per hit, as `<accession>:<start>-<end> (<feature>)` — e.g. `NC_000005.10:151268868-151268889 (exon of GM2A)`. Feature is always one of `HERV-K family: <element>`, `exon of <gene>`, `intron of <gene>`, or `outside any gene`. |

Probes are ranked (`rank`) by eFISHent's `quality` score (descending), tied
by `hervk_family_hits` (ascending — fewer is better), tied by `intron_hits`
(ascending). `outside_gene_hits` never affects rank. Any probe with
`exon_hits > 0` is excluded from this ranking and only appears in the top-N
selection as a logged fallback if the target has too few exon-hit-free
candidates to fill `--top-n` otherwise.

Two different things can make a target's selection less than ideal, and
they're logged separately: `target_short` means there weren't enough
candidates *in total*; a backfill warning means there were enough
candidates, just not enough *exon-hit-free* ones. A target can hit either,
both, or neither.

#### Columns from eFISHent (passed through unchanged)

These come straight from your input CSVs — ProbeK doesn't recompute or
alter them. Brief descriptions below; see
[eFISHent's own documentation](https://github.com/BBQuercus/eFISHent) for
authoritative definitions, particularly for `kmers`/`count`, which are
inputs to eFISHent's own transcriptome-level off-target scoring that ProbeK
doesn't otherwise use.

| Column | Meaning |
|---|---|
| `name` | Probe identifier, e.g. `pNRV101_gag-17`. |
| `sequence` | The probe's nucleotide sequence — what ProbeK actually BLASTs. |
| `length` | Probe length in bases. |
| `start` / `end` | Probe position within eFISHent's target transcript. |
| `GC` | GC content (%). |
| `TM` | Predicted melting temperature. |
| `deltaG` | Predicted free energy (secondary structure/duplex stability). |
| `kmers` / `count` | Inputs to eFISHent's own k-mer-based off-target scoring. |
| `cpg_fraction` | Fraction of the probe that is CpG dinucleotides. |
| `low_complexity` | Sequence-complexity score (low-complexity/repetitive content). |
| `accessibility` | Predicted accessibility of the target region (secondary structure). |
| `on_target_dg` | Predicted on-target binding free energy. |
| `txome_off_targets` | Count of *transcriptome*-level off-targets eFISHent found (separate from, and not recomputed by, ProbeK's genomic classification). |
| `off_target_genes` | Ensembl transcript IDs eFISHent flagged from the transcriptome screen — informational context only. |
| `worst_match` | eFISHent's worst observed off-target alignment, formatted `<identity%>/<length>bp/<mismatches>mm`. |
| `expression_risk` | eFISHent's expression-level risk assessment. |
| `quality` | eFISHent's composite quality score — higher is better. Used as ProbeK's final ranking tiebreaker. |
| `recommendation` | eFISHent's own categorical flag: blank/OK, `FLAG(...)`, or `FAIL`. Rows already marked `FAIL` are dropped before ProbeK does anything (see [What it does](#what-it-does)). |

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
