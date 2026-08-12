# Drop your eFISHent CSVs here

This is where ProbeK looks for input by default (`probek` with no `--input`
flag reads this directory).

Drop one CSV per gene target here (e.g. `pNRV101_gag.csv`, `pNRV102_pro.csv`),
or a combined file — ProbeK doesn't care how many files you give it. Each CSV
must be eFISHent's standard output format with these columns:

```
name,sequence,length,start,end,GC,TM,deltaG,kmers,count,cpg_fraction,low_complexity,
accessibility,on_target_dg,txome_off_targets,off_target_genes,worst_match,
expression_risk,quality,recommendation
```

**How target labels are derived:** ProbeK strips a trailing `-<number>` from
each row's `name` column (e.g. `pNRV101_gag-1` → `pNRV101_gag`) and uses that
as the target label, provided every row in the file agrees. If they don't
agree, it falls back to the filename (without `.csv`). Two files that happen
to derive the same target label are simply combined.

Rows with `recommendation == FAIL` are dropped automatically before any
BLAST search runs — no need to pre-filter them yourself.

This README is the only file in this directory that's tracked in git — the
CSVs you drop here are treated as your local working data and are gitignored.
