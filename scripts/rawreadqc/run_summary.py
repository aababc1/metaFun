#!/usr/bin/env python3
"""
run_summary.py — one TSV row per run for RAWREAD_QC_cleanifier (fastp -> Cleanifier).

Reads the fastp JSON report and the Cleanifier filter log, checks that the two tools
agree on the number of read pairs (fastp pairs after filtering == Cleanifier keep + filter),
checks that both output FASTQ files end with the BGZF EOF block (i.e. bgzip finished),
and prints a header + one data row to stdout. Exits non-zero on any inconsistency,
which makes the Nextflow task fail *before* the outputs are moved into place.

Only the Python standard library is used.
"""
import argparse
import datetime as dt
import json
import os
import re
import sys

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

COLUMNS = [
    "run", "sample", "cohort",
    "raw_pairs", "raw_bases", "pct_bases_vs_manifest",
    "fastp_pairs", "fastp_bases", "pct_fastp_pass",
    "host_pairs", "kept_pairs", "pct_host",
    "kept_bases_est", "pct_kept_of_raw",
    "q30_after", "mean_len_r1_after", "mean_len_r2_after",
    "adapter_trimmed_reads", "low_quality_reads", "too_short_reads", "too_many_N_reads",
    "wall_sec", "cpus", "attempt", "fastp_version", "cleanifier_version", "finished",
]


def die(msg, code=2):
    print(f"run_summary.py: ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def parse_cleanifier_log(path):
    """Return (keep, filter) read-pair counts from `cleanifier filter` stdout."""
    keep = flt = None
    with open(path, errors="replace") as fh:
        lines = fh.read().splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "prefix\tkeep\tfilter" and i + 1 < len(lines):
            parts = lines[i + 1].rstrip("\n").split("\t")
            if len(parts) >= 3:
                try:
                    keep, flt = int(parts[1]), int(parts[2])
                except ValueError:
                    pass
            break
    if keep is None:
        # fallback: the pretty table "| keep      | 1032 | 51.60% |"
        for line in lines:
            m = re.match(r"\|\s*(keep|filter)\s*\|\s*(\d+)\s*\|", line)
            if m:
                if m.group(1) == "keep":
                    keep = int(m.group(2))
                else:
                    flt = int(m.group(2))
    if keep is None or flt is None:
        die(f"could not find 'Classification Statistics' (keep/filter counts) in {path}")
    return keep, flt


def check_bgzf(path, expect_empty_ok):
    if not os.path.isfile(path):
        die(f"output file missing: {path}")
    size = os.path.getsize(path)
    if size < len(BGZF_EOF):
        die(f"output file too small ({size} bytes): {path}")
    with open(path, "rb") as fh:
        fh.seek(-len(BGZF_EOF), os.SEEK_END)
        tail = fh.read(len(BGZF_EOF))
    if tail != BGZF_EOF:
        die(f"BGZF EOF block missing (bgzip did not finish cleanly): {path}")
    if size == len(BGZF_EOF) and not expect_empty_ok:
        die(f"output file contains no reads although kept_pairs > 0: {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--fastp-json", required=True)
    ap.add_argument("--cleanifier-log", required=True)
    ap.add_argument("--fastq", nargs=2, required=True, metavar=("R1_GZ", "R2_GZ"),
                    help="the two (temporary) output FASTQ files written by Cleanifier/bgzip")
    ap.add_argument("--wall-sec", type=int, default=-1)
    ap.add_argument("--cpus", type=int, default=0)
    ap.add_argument("--attempt", type=int, default=1)
    ap.add_argument("--cleanifier-version", default="")
    ap.add_argument("--expected-bases", type=float, default=0,
                    help="base_count from the SRA run manifest; raw bases seen by fastp must be >= (1 - tolerance) x this, "
                         "otherwise the input FASTQ is considered truncated (fastp silently accepts a truncated gzip stream)")
    ap.add_argument("--bases-tolerance", type=float, default=0.05)
    ap.add_argument("--no-check", action="store_true", help="skip consistency checks (not recommended)")
    args = ap.parse_args()

    with open(args.fastp_json) as fh:
        fp = json.load(fh)
    s = fp.get("summary", {})
    before = s.get("before_filtering", {})
    after = s.get("after_filtering", {})
    fr = fp.get("filtering_result", {})
    ad = fp.get("adapter_cutting", {})

    raw_reads = int(before.get("total_reads", 0))
    raw_bases = int(before.get("total_bases", 0))
    fp_reads = int(after.get("total_reads", 0))
    fp_bases = int(after.get("total_bases", 0))
    if raw_reads % 2 or fp_reads % 2:
        die(f"fastp read counts are odd (not paired?): before={raw_reads} after={fp_reads}")
    raw_pairs, fp_pairs = raw_reads // 2, fp_reads // 2

    keep, flt = parse_cleanifier_log(args.cleanifier_log)

    pct_bases_vs_manifest = "NA"
    if args.expected_bases > 0:
        pct_bases_vs_manifest = f"{100.0 * raw_bases / args.expected_bases:.2f}"
    if not args.no_check:
        if args.expected_bases > 0 and raw_bases < args.expected_bases * (1 - args.bases_tolerance):
            die(f"input looks truncated: fastp saw {raw_bases} raw bases but the manifest expects {args.expected_bases:.0f} "
                f"({100.0 * raw_bases / args.expected_bases:.1f}% < {100 * (1 - args.bases_tolerance):.0f}%)")
        if keep + flt != fp_pairs:
            die(f"pair count mismatch: fastp after_filtering = {fp_pairs} pairs, "
                f"Cleanifier keep+filter = {keep}+{flt} = {keep + flt} (stream truncated?)")
        for f in args.fastq:
            check_bgzf(f, expect_empty_ok=(keep == 0))

    def pct(a, b):
        return f"{(100.0 * a / b):.3f}" if b else "NA"

    kept_bases_est = int(round(fp_bases * keep / fp_pairs)) if fp_pairs else 0
    row = {
        "run": args.run, "sample": args.sample, "cohort": args.cohort,
        "raw_pairs": raw_pairs, "raw_bases": raw_bases, "pct_bases_vs_manifest": pct_bases_vs_manifest,
        "fastp_pairs": fp_pairs, "fastp_bases": fp_bases, "pct_fastp_pass": pct(fp_pairs, raw_pairs),
        "host_pairs": flt, "kept_pairs": keep, "pct_host": pct(flt, keep + flt),
        "kept_bases_est": kept_bases_est, "pct_kept_of_raw": pct(keep, raw_pairs),
        "q30_after": after.get("q30_rate", "NA"),
        "mean_len_r1_after": after.get("read1_mean_length", "NA"),
        "mean_len_r2_after": after.get("read2_mean_length", "NA"),
        "adapter_trimmed_reads": ad.get("adapter_trimmed_reads", "NA"),
        "low_quality_reads": fr.get("low_quality_reads", "NA"),
        "too_short_reads": fr.get("too_short_reads", "NA"),
        "too_many_N_reads": fr.get("too_many_N_reads", "NA"),
        "wall_sec": args.wall_sec, "cpus": args.cpus, "attempt": args.attempt,
        "fastp_version": s.get("fastp_version", "NA"),
        "cleanifier_version": args.cleanifier_version or "NA",
        "finished": dt.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
    }
    print("\t".join(COLUMNS))
    print("\t".join(str(row[c]) for c in COLUMNS))


if __name__ == "__main__":
    main()
