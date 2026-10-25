#!/usr/bin/env python3
"""
qc_summary.py — merge per-run TSV rows (<qc_logs>/*.qc.tsv, written by run_summary.py)
into one cohort table <cohort>.qc_summary.tsv, sorted by run, one row per run
(the most recent row wins if a run was processed more than once, e.g. after a retry
or a re-submission). The file is written atomically (temp file + rename).

Also prints a short distribution summary (host fraction, wall time) to stdout.
Standard library only.
"""
import argparse
import glob
import os
import statistics
import sys
import tempfile


def read_rows(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return None, []
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:] if not l.startswith("run\t")]
    return header, rows


def quantiles(values):
    if not values:
        return "n=0"
    v = sorted(values)
    n = len(v)
    q = lambda p: v[min(n - 1, int(round(p * (n - 1))))]
    return (f"n={n} min={v[0]:.3g} p25={q(0.25):.3g} median={q(0.5):.3g} "
            f"p75={q(0.75):.3g} max={v[-1]:.3g} mean={statistics.fmean(v):.3g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-logs", required=True, help="directory with <run>.qc.tsv files")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--out", required=True, help="output path (<dest>/<cohort>.qc_summary.tsv)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.qc_logs, "*.qc.tsv")))
    header = None
    by_run = {}
    for f in files:
        h, rows = read_rows(f)
        if h is None:
            continue
        header = header or h
        for r in rows:
            if r.get("cohort", args.cohort) != args.cohort:
                continue
            run = r.get("run")
            if not run:
                continue
            prev = by_run.get(run)
            if prev is None or r.get("finished", "") >= prev.get("finished", ""):
                by_run[run] = r
    if header is None:
        print(f"qc_summary.py: no *.qc.tsv found in {args.qc_logs}", file=sys.stderr)
        header = ["run"]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{args.cohort}.qc_summary.", suffix=".tmp",
                               dir=os.path.dirname(os.path.abspath(args.out)))
    with os.fdopen(fd, "w") as out:
        out.write("\t".join(header) + "\n")
        for run in sorted(by_run):
            out.write("\t".join(str(by_run[run].get(c, "NA")) for c in header) + "\n")
    umask = os.umask(0)
    os.umask(umask)
    os.chmod(tmp, 0o666 & ~umask)   # mkstemp creates 0600; use the normal umask instead
    os.replace(tmp, args.out)

    if not args.quiet:
        def col(name, cast=float):
            vals = []
            for r in by_run.values():
                try:
                    vals.append(cast(r[name]))
                except (KeyError, ValueError, TypeError):
                    pass
            return vals
        print(f"{args.cohort}: {len(by_run)} run(s) -> {args.out}")
        print(f"  pct_host        : {quantiles(col('pct_host'))}")
        print(f"  pct_kept_of_raw : {quantiles(col('pct_kept_of_raw'))}")
        print(f"  wall_sec        : {quantiles(col('wall_sec'))}")
        print(f"  raw_bases (Gb)  : {quantiles([b / 1e9 for b in col('raw_bases')])}")


if __name__ == "__main__":
    main()
