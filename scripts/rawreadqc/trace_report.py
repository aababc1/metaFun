#!/usr/bin/env python3
"""
trace_report.py — summarise a Nextflow trace.txt (raw=true, fields as set in the kbds profile)
plus the cohort qc_summary.tsv into the numbers needed for a test report:
wall time, per-run realtime, %cpu, memory, node-hours and a budget extrapolation.

  trace_report.py --trace ~/metaCRC/runs/qc/PRJNA531273/trace.txt \
      --qc-summary ~/metaCRC/qc/PRJNA531273/PRJNA531273.qc_summary.tsv \
      --node-rate 0.5 --total-runs 8420 [--total-bases 2.2e14]

Standard library only. Prints Markdown.
"""
import argparse
import csv
import statistics
from collections import Counter, defaultdict


def num(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def q(values, p):
    if not values:
        return float("nan")
    v = sorted(values)
    return v[min(len(v) - 1, int(round(p * (len(v) - 1))))]


def fmt_dur(sec):
    if sec is None or sec != sec:
        return "NA"
    if sec < 90:
        return f"{sec:.0f} s"
    if sec < 5400:
        return f"{sec/60:.1f} min"
    return f"{sec/3600:.2f} h"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", required=True)
    ap.add_argument("--qc-summary")
    ap.add_argument("--node-rate", type=float, default=0.5, help="node-hours charged per wall hour (cpu32-only 0.5, cpu64 1)")
    ap.add_argument("--slots", type=int, default=0, help="concurrent tasks per node (nproc/cpus); inferred from the trace if 0")
    ap.add_argument("--total-runs", type=int, default=0, help="extrapolate node-hours to this many runs")
    ap.add_argument("--total-bases", type=float, default=0, help="extrapolate by bases instead of by run count")
    ap.add_argument("--process", default="fastp_cleanifier", help="process name of the per-run task")
    args = ap.parse_args()

    with open(args.trace) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        print("empty trace"); return

    procs = sorted({r["process"] for r in rows})
    print(f"## trace: {args.trace}")
    print(f"- process names (`cut -f4 trace.txt | sort -u`): " + ", ".join(f"`{p}`" for p in procs))
    print()
    print("| process | status | n | realtime p50 | realtime mean | realtime max | %cpu p50 | %cpu mean | peak_rss max | read (rchar) | written (wchar) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for p in procs:
        by_status = defaultdict(list)
        for r in rows:
            if r["process"] == p:
                by_status[r["status"]].append(r)
        for st, rs in sorted(by_status.items()):
            rt = [num(r["realtime"]) / 1000 for r in rs if num(r["realtime"]) is not None]
            cpu = [num(r["%cpu"]) for r in rs if num(r["%cpu"]) is not None]
            rss = [num(r["peak_rss"]) for r in rs if num(r["peak_rss"]) is not None]
            rch = sum(num(r["rchar"], 0) for r in rs)
            wch = sum(num(r["wchar"], 0) for r in rs)
            print(f"| {p} | {st} | {len(rs)} | {fmt_dur(q(rt, .5))} | {fmt_dur(statistics.fmean(rt) if rt else None)} | "
                  f"{fmt_dur(max(rt) if rt else None)} | {q(cpu, .5):.0f}% | {statistics.fmean(cpu) if cpu else float('nan'):.0f}% | "
                  f"{(max(rss)/2**30) if rss else float('nan'):.1f} GiB | {rch/2**30:.1f} GiB | {wch/2**30:.1f} GiB |")
    print()

    task = [r for r in rows if r["process"] == args.process]
    done = [r for r in task if r["status"] in ("COMPLETED", "CACHED")]
    failed = [r for r in task if r["status"] == "FAILED"]
    starts = [num(r["start"]) for r in rows if num(r["start"])]
    ends = [num(r["complete"]) for r in rows if num(r["complete"])]
    submits = [num(r["submit"]) for r in rows if num(r["submit"])]
    wall_s = (max(ends) - min(submits)) / 1000 if ends and submits else float("nan")
    rt_done = [num(r["realtime"]) / 1000 for r in done if num(r["realtime"]) is not None]
    cpus = Counter(r["cpus"] for r in task).most_common(1)[0][0] if task else "?"

    # concurrency actually observed: max number of per-run tasks overlapping at any instant
    events = []
    for r in task:
        s, e = num(r["start"]), num(r["complete"])
        if s and e:
            events += [(s, 1), (e, -1)]
    events.sort()
    cur = peak = 0
    for _, d in events:
        cur += d
        peak = max(peak, cur)
    slots = args.slots or peak

    node_h = wall_s / 3600 * args.node_rate if wall_s == wall_s else float("nan")
    per_run_node_h = node_h / len(done) if done else float("nan")
    print(f"## per-run task `{args.process}`")
    print(f"- runs completed: {len(done)}  failed rows: {len(failed)}  cpus/task: {cpus}  peak concurrency: {peak}")
    print(f"- realtime per run: p50 {fmt_dur(q(rt_done, .5))}, mean {fmt_dur(statistics.fmean(rt_done) if rt_done else None)}, "
          f"max {fmt_dur(max(rt_done) if rt_done else None)}, sum {fmt_dur(sum(rt_done) if rt_done else None)}")
    print(f"- wall (first submit -> last complete): {fmt_dur(wall_s)}  ->  node-hours = {node_h:.2f} at rate {args.node_rate}/h")
    if done and wall_s == wall_s:
        print(f"- throughput: {len(done) / (wall_s / 3600):.1f} runs / node-hour-of-wall;  "
              f"{wall_s / len(done) / 60:.2f} wall-min per run;  {per_run_node_h * 60:.2f} node-min per run")
        slot_min = sum(rt_done) / len(done) / 60
        print(f"- slot-minutes per run (realtime mean) = {slot_min:.1f}  (budget target: <= 17 min at 4 threads)")

    if args.qc_summary:
        with open(args.qc_summary) as fh:
            qs = list(csv.DictReader(fh, delimiter="\t"))
        host = [num(r["pct_host"]) for r in qs if num(r["pct_host"]) is not None]
        kept = [num(r["pct_kept_of_raw"]) for r in qs if num(r["pct_kept_of_raw"]) is not None]
        gb = [num(r["raw_bases"], 0) / 1e9 for r in qs]
        wall = [num(r["wall_sec"]) for r in qs if num(r["wall_sec"]) is not None]
        print()
        print(f"## qc_summary: {args.qc_summary} ({len(qs)} runs)")
        if host:
            print(f"- host (human) read-pair fraction %: min {min(host):.2f}, p25 {q(host,.25):.2f}, median {q(host,.5):.2f}, "
                  f"p75 {q(host,.75):.2f}, max {max(host):.2f}, mean {statistics.fmean(host):.2f}")
        if kept:
            print(f"- pairs kept after fastp+Cleanifier (% of raw): median {q(kept,.5):.2f}, min {min(kept):.2f}, max {max(kept):.2f}")
        if gb:
            print(f"- raw Gbases per run: median {q(gb,.5):.2f}, total {sum(gb):.1f}")
            if wall and sum(gb) > 0:
                print(f"- task wall seconds per raw Gbase (median over runs): "
                      f"{q([w / g for w, g in zip(wall, gb) if g > 0], .5):.0f} s/Gb")
            if args.total_bases > 0 and node_h == node_h:
                print(f"- extrapolation by bases: {args.total_bases:.3g} bases -> "
                      f"{node_h * args.total_bases / (sum(gb) * 1e9):.0f} node-hours")
    if args.total_runs > 0 and per_run_node_h == per_run_node_h:
        print(f"- extrapolation by run count: {args.total_runs} runs -> {per_run_node_h * args.total_runs:.0f} node-hours "
              f"(cohorts differ in run size; prefer the per-base estimate)")


if __name__ == "__main__":
    main()
