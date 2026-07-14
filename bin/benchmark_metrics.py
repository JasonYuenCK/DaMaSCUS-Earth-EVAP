#!/usr/bin/env python3
"""Extract and compare DaMaSCUS-SUN performance/physics regression metrics."""

import argparse
import csv
import math
import os
import re
import statistics
import sys
from collections import defaultdict


FIELDS = [
    "variant",
    "case",
    "ranks",
    "repeat",
    "commit",
    "host",
    "run_status",
    "wall_seconds",
    "simulation_seconds",
    "max_rss_kb",
    "trajectories",
    "valid_trajectories",
    "captured",
    "capture_rate_valid",
    "capture_ci_lower",
    "capture_ci_upper",
    "numerical_failures",
    "sync_interval",
    "sync_rounds",
    "final_round_trajectories",
    "capture_overshoot",
    "total_scatterings",
    "trajectory_rate_per_second",
]


def read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def scalar(text, label, number_type=int):
    match = re.search(r"^" + re.escape(label) + r"\s*([^\s]+)", text, re.MULTILINE)
    if not match:
        return None
    try:
        return number_type(match.group(1))
    except ValueError:
        return None


def elapsed_seconds(value):
    parts = value.strip().split(":")
    try:
        if len(parts) == 2:
            return 60.0 * float(parts[0]) + float(parts[1])
        if len(parts) == 3:
            return 3600.0 * float(parts[0]) + 60.0 * float(parts[1]) + float(parts[2])
        return float(value)
    except ValueError:
        return None


def time_metrics(text):
    elapsed = None
    max_rss = None
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            elapsed = elapsed_seconds(line.rsplit("): ", 1)[-1])
        elif "Maximum resident set size (kbytes):" in line:
            try:
                max_rss = int(line.rsplit(":", 1)[-1].strip())
            except ValueError:
                pass
    return elapsed, max_rss


def optional(value):
    return "" if value is None else value


def extract(args):
    log_text = read_text(args.log)
    time_text = read_text(args.time)
    wall_seconds, max_rss_kb = time_metrics(time_text)

    trajectories = scalar(log_text, "Simulated trajectories:")
    valid = scalar(log_text, "Capture-classified trajectories:")
    captured = scalar(log_text, "Captured count:")
    simulation_seconds = scalar(log_text, "Simulation time [s]:", float)
    numerical_failures = scalar(log_text, "Numerical failure count:")
    sync_interval = scalar(log_text, "Normal-mode MPI sync interval:")
    sync_rounds = scalar(log_text, "MPI sync rounds:")
    final_round = scalar(log_text, "Final MPI round trajectories:")
    capture_overshoot = scalar(log_text, "Capture target overshoot:")
    total_scatterings = scalar(log_text, "Total # of scatterings:")

    capture_rate = None
    if valid not in (None, 0) and captured is not None:
        capture_rate = float(captured) / float(valid)
    if capture_rate is None:
        capture_rate = scalar(log_text, "Capture rate valid:", float)

    ci_lower = None
    ci_upper = None
    ci_match = re.search(
        r"^Capture rate valid 95% CI:\s*\[\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\]",
        log_text,
        re.MULTILINE,
    )
    if ci_match:
        ci_lower = float(ci_match.group(1))
        ci_upper = float(ci_match.group(2))

    trajectory_rate = None
    if trajectories is not None and simulation_seconds not in (None, 0.0):
        trajectory_rate = float(trajectories) / simulation_seconds

    row = {
        "variant": args.variant,
        "case": args.case,
        "ranks": args.ranks,
        "repeat": args.repeat,
        "commit": args.commit,
        "host": args.host,
        "run_status": args.run_status,
        "wall_seconds": optional(wall_seconds),
        "simulation_seconds": optional(simulation_seconds),
        "max_rss_kb": optional(max_rss_kb),
        "trajectories": optional(trajectories),
        "valid_trajectories": optional(valid),
        "captured": optional(captured),
        "capture_rate_valid": optional(capture_rate),
        "capture_ci_lower": optional(ci_lower),
        "capture_ci_upper": optional(ci_upper),
        "numerical_failures": optional(numerical_failures),
        "sync_interval": optional(sync_interval),
        "sync_rounds": optional(sync_rounds),
        "final_round_trajectories": optional(final_round),
        "capture_overshoot": optional(capture_overshoot),
        "total_scatterings": optional(total_scatterings),
        "trajectory_rate_per_second": optional(trajectory_rate),
    }

    needs_header = not os.path.exists(args.output) or os.path.getsize(args.output) == 0
    with open(args.output, "a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        if needs_header:
            writer.writeheader()
        writer.writerow(row)


def numeric(rows, field, number_type=float):
    values = []
    for row in rows:
        value = row.get(field, "")
        if value == "":
            continue
        try:
            values.append(number_type(value))
        except ValueError:
            continue
    return values


def median_value(rows, field, number_type=float):
    values = numeric(rows, field, number_type)
    return statistics.median(values) if values else None


def format_value(value, digits=6):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        return f"{value:.{digits}g}"
    return str(value)


def load_rows(paths):
    rows = []
    for path in paths:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def compare(args):
    grouped = defaultdict(list)
    for row in load_rows(args.metrics):
        grouped[(row["variant"], row["case"], int(row["ranks"]))].append(row)

    baseline_keys = {(case, ranks) for variant, case, ranks in grouped if variant == args.baseline}
    candidate_keys = {(case, ranks) for variant, case, ranks in grouped if variant == args.candidate}
    common_keys = sorted(baseline_keys & candidate_keys, key=lambda item: (item[0], item[1]))
    if not common_keys:
        print("No matching case/rank groups were found for the requested variants.", file=sys.stderr)
        return 2

    print(
        "case\tranks\tbaseline_wall_s\tcandidate_wall_s\tspeedup"
        "\tbaseline_capture_rate\tcandidate_capture_rate\tz_score\tfailure_rate\tstatus"
    )
    failed = False
    for case, ranks in common_keys:
        baseline_rows = grouped[(args.baseline, case, ranks)]
        candidate_rows = grouped[(args.candidate, case, ranks)]

        baseline_wall = median_value(baseline_rows, "wall_seconds")
        candidate_wall = median_value(candidate_rows, "wall_seconds")
        speedup = None
        if baseline_wall is not None and candidate_wall not in (None, 0.0):
            speedup = baseline_wall / candidate_wall

        # Repeats intentionally use the same seed for timing stability, so they
        # are not pooled as independent Monte Carlo samples. Median counters
        # represent one typical repeat for the statistical comparison.
        baseline_valid = median_value(baseline_rows, "valid_trajectories")
        candidate_valid = median_value(candidate_rows, "valid_trajectories")
        baseline_captured = median_value(baseline_rows, "captured")
        candidate_captured = median_value(candidate_rows, "captured")
        baseline_rate = None
        candidate_rate = None
        z_score = None
        if baseline_valid not in (None, 0.0) and baseline_captured is not None:
            baseline_rate = baseline_captured / baseline_valid
        if candidate_valid not in (None, 0.0) and candidate_captured is not None:
            candidate_rate = candidate_captured / candidate_valid
        if baseline_rate is not None and candidate_rate is not None:
            variance = baseline_rate * (1.0 - baseline_rate) / baseline_valid
            variance += candidate_rate * (1.0 - candidate_rate) / candidate_valid
            if variance > 0.0:
                z_score = abs(candidate_rate - baseline_rate) / math.sqrt(variance)
            else:
                z_score = 0.0 if candidate_rate == baseline_rate else float("inf")

        candidate_failures = median_value(candidate_rows, "numerical_failures")
        candidate_trajectories = median_value(candidate_rows, "trajectories")
        failure_rate = None
        if candidate_failures is not None and candidate_trajectories not in (None, 0.0):
            failure_rate = candidate_failures / candidate_trajectories

        reasons = []
        baseline_statuses = numeric(baseline_rows, "run_status", int)
        candidate_statuses = numeric(candidate_rows, "run_status", int)
        baseline_hosts = {row.get("host", "") for row in baseline_rows if row.get("host", "")}
        candidate_hosts = {row.get("host", "") for row in candidate_rows if row.get("host", "")}
        if not baseline_statuses or any(status != 0 for status in baseline_statuses):
            reasons.append("baseline-run")
        if not candidate_statuses or any(status != 0 for status in candidate_statuses):
            reasons.append("candidate-run")
        if not args.allow_host_mismatch and (not baseline_hosts or baseline_hosts != candidate_hosts):
            reasons.append("host-mismatch")
        if baseline_wall is None or candidate_wall is None:
            reasons.append("timing-metrics")
        if z_score is None:
            reasons.append("physics-metrics")
        if failure_rate is None:
            reasons.append("failure-metrics")
        if z_score is not None and z_score > args.max_z_score:
            reasons.append("physics")
        if failure_rate is not None and failure_rate > args.max_failure_rate:
            reasons.append("failures")
        if args.min_speedup is not None:
            if speedup is None:
                reasons.append("performance-metrics")
            elif speedup < args.min_speedup:
                reasons.append("performance")
        status = "FAIL:" + ",".join(reasons) if reasons else "PASS"
        failed = failed or bool(reasons)

        print(
            "\t".join(
                [
                    case,
                    str(ranks),
                    format_value(baseline_wall),
                    format_value(candidate_wall),
                    format_value(speedup),
                    format_value(baseline_rate),
                    format_value(candidate_rate),
                    format_value(z_score),
                    format_value(failure_rate),
                    status,
                ]
            )
        )
    return 1 if failed else 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    # CHPC currently provides Python 3.6, where add_subparsers(required=True)
    # is not available yet. main() performs the equivalent validation.
    subparsers = parser.add_subparsers(dest="command")

    extract_parser = subparsers.add_parser("extract", help="append metrics from one benchmark repeat")
    extract_parser.add_argument("--log", required=True)
    extract_parser.add_argument("--time", required=True)
    extract_parser.add_argument("--output", required=True)
    extract_parser.add_argument("--variant", required=True)
    extract_parser.add_argument("--case", required=True)
    extract_parser.add_argument("--ranks", required=True, type=int)
    extract_parser.add_argument("--repeat", required=True, type=int)
    extract_parser.add_argument("--commit", required=True)
    extract_parser.add_argument("--host", required=True)
    extract_parser.add_argument("--run-status", required=True, type=int)
    extract_parser.set_defaults(func=extract)

    compare_parser = subparsers.add_parser("compare", help="compare matched baseline/candidate result groups")
    compare_parser.add_argument("--baseline", required=True)
    compare_parser.add_argument("--candidate", required=True)
    compare_parser.add_argument("--max-z-score", type=float, default=3.0)
    compare_parser.add_argument("--max-failure-rate", type=float, default=1.0e-4)
    compare_parser.add_argument("--min-speedup", type=float)
    compare_parser.add_argument("--allow-host-mismatch", action="store_true")
    compare_parser.add_argument("metrics", nargs="+")
    compare_parser.set_defaults(func=compare)
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    result = args.func(args)
    return 0 if result is None else result


if __name__ == "__main__":
    sys.exit(main())
