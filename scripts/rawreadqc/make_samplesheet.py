#!/usr/bin/env python3
"""
make_samplesheet.py — build RAWREAD_QC_cleanifier samplesheets from a run manifest.

Manifest (TSV, header) must contain: run_accession, sample_accession, cohort_accession,
library_layout; base_count is used for chunk balancing when present.
FASTQ files are expected at <fastq-dir>/<cohort>/<run>_1.fastq.gz and _2.fastq.gz
(override the pattern with --fastq-pattern).

Examples
  # one cohort, one samplesheet
  make_samplesheet.py --manifest run_manifest.tsv --fastq-dir ~/metaCRC/fastq \
      --cohort PRJNA531273 --out samplesheet.csv

  # split a large cohort into 3 chunks balanced by base_count -> samplesheet.1.csv .. samplesheet.3.csv
  make_samplesheet.py --manifest run_manifest.tsv --fastq-dir ~/metaCRC/fastq \
      --cohort PRJEB6070 --chunks 3 --out samplesheet.csv

  # cap each chunk at 3e12 bases instead of a fixed chunk count
  make_samplesheet.py ... --cohort PRJEB6070 --max-bases 3e12 --out samplesheet.csv

  # every PAIRED cohort: <out-dir>/<cohort>/samplesheet.csv (+ .k.csv chunks if --max-bases)
  make_samplesheet.py --manifest run_manifest.tsv --fastq-dir ~/metaCRC/fastq --all --out-dir ~/metaCRC/runs/qc

Rows whose FASTQ pair does not exist are skipped with a warning unless --allow-missing is given.
The optional 6th column base_count (from the manifest) lets the workflow detect truncated FASTQ
input (fastp raw bases must reach 95 % of it).
Standard library only.
"""
import argparse
import csv
import math
import os
import sys

COLS = ["run_accession", "sample_accession", "cohort_accession", "fastq_1", "fastq_2", "base_count"]


def read_manifest(path):
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
    need = {"run_accession", "sample_accession", "cohort_accession", "library_layout"}
    missing = need - set(rdr.fieldnames or [])
    if missing:
        sys.exit(f"manifest lacks columns: {sorted(missing)}")
    return rows


def balance(rows, nchunks):
    """Greedy longest-processing-time assignment by base_count (largest runs first)."""
    def bases(r):
        try:
            return float(r.get("base_count") or 0)
        except ValueError:
            return 0.0
    chunks = [[] for _ in range(nchunks)]
    load = [0.0] * nchunks
    for r in sorted(rows, key=bases, reverse=True):
        i = load.index(min(load))
        chunks[i].append(r)
        load[i] += bases(r)
    return [c for c in chunks if c], load


def write_sheet(path, rows, fastq_dir, pattern, allow_missing, excluded):
    n_written = n_missing = 0
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for r in rows:
            run, sample, cohort = r["run_accession"], r["sample_accession"], r["cohort_accession"]
            if run in excluded:
                continue
            fq1 = os.path.join(fastq_dir, pattern.format(cohort=cohort, run=run, mate=1))
            fq2 = os.path.join(fastq_dir, pattern.format(cohort=cohort, run=run, mate=2))
            if not allow_missing and not (os.path.isfile(fq1) and os.path.isfile(fq2)):
                n_missing += 1
                continue
            w.writerow([run, sample, cohort, fq1, fq2, r.get("base_count") or ""])
            n_written += 1
    return n_written, n_missing


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--fastq-pattern", default="{cohort}/{run}_{mate}.fastq.gz",
                    help="relative to --fastq-dir; placeholders {cohort} {run} {mate}")
    ap.add_argument("--layout", default="PAIRED", help="library_layout to keep (default PAIRED)")
    ap.add_argument("--cohort", action="append", default=[], help="cohort_accession (repeatable)")
    ap.add_argument("--all", action="store_true", help="one samplesheet per cohort into --out-dir")
    ap.add_argument("--out", default="samplesheet.csv", help="output CSV (chunks: <stem>.<k>.csv)")
    ap.add_argument("--out-dir", default=".", help="with --all: <out-dir>/<cohort>/samplesheet[.k].csv")
    ap.add_argument("--chunks", type=int, default=1, help="split into this many chunks (balanced by base_count)")
    ap.add_argument("--max-bases", type=float, default=0,
                    help="instead of --chunks: smallest chunk count so that each chunk has <= this many bases")
    ap.add_argument("--exclude", help="file with run_accession per line to leave out")
    ap.add_argument("--allow-missing", action="store_true", help="keep rows even if the FASTQ pair is absent")
    args = ap.parse_args()

    rows = [r for r in read_manifest(args.manifest) if r.get("library_layout", "").upper() == args.layout.upper()]
    excluded = set()
    if args.exclude:
        with open(args.exclude) as fh:
            excluded = {l.strip().split()[0] for l in fh if l.strip() and not l.startswith("#")}

    fastq_dir = os.path.abspath(os.path.expanduser(args.fastq_dir))
    cohorts = sorted({r["cohort_accession"] for r in rows}) if args.all else args.cohort
    if not cohorts:
        sys.exit("give --cohort <id> (repeatable) or --all")

    for cohort in cohorts:
        crows = [r for r in rows if r["cohort_accession"] == cohort]
        if not crows:
            print(f"[{cohort}] no {args.layout} runs in manifest", file=sys.stderr)
            continue
        total_bases = sum(float(r.get("base_count") or 0) for r in crows)
        nchunks = args.chunks
        if args.max_bases > 0:
            nchunks = max(1, math.ceil(total_bases / args.max_bases))
        chunks, load = balance(crows, nchunks) if nchunks > 1 else ([crows], [total_bases])

        if args.all:
            out_dir = os.path.join(os.path.abspath(os.path.expanduser(args.out_dir)), cohort)
            os.makedirs(out_dir, exist_ok=True)
            stem, ext = os.path.join(out_dir, "samplesheet"), ".csv"
        else:
            base = os.path.abspath(os.path.expanduser(args.out))
            stem, ext = os.path.splitext(base)
            ext = ext or ".csv"

        for k, chunk in enumerate(chunks, start=1):
            path = f"{stem}{ext}" if len(chunks) == 1 else f"{stem}.{k}{ext}"
            n, miss = write_sheet(path, chunk, fastq_dir, args.fastq_pattern, args.allow_missing, excluded)
            tag = f"{cohort}" if len(chunks) == 1 else f"{cohort}.{k}"
            print(f"[{tag}] {n} runs, {load[k-1]/1e9:.1f} Gbases -> {path}"
                  + (f"  ({miss} run(s) skipped: FASTQ pair missing)" if miss else ""))


if __name__ == "__main__":
    main()
