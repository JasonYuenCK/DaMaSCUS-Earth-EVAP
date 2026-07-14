#!/usr/bin/env python3

import csv
import importlib.util
import io
import pathlib
import tempfile
import types
import unittest
from contextlib import redirect_stdout


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
MODULE_PATH = PROJECT_ROOT / "bin" / "benchmark_metrics.py"
SPEC = importlib.util.spec_from_file_location("benchmark_metrics", MODULE_PATH)
benchmark_metrics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark_metrics)


class BenchmarkMetricsTests(unittest.TestCase):
    def test_extracts_runtime_physics_and_batch_metrics(self):
        log = """\
Simulated trajectories:\t\t32000000
Capture-classified trajectories:\t31999990
Captured count:\t\t\t10003
Capture rate valid 95% CI:\t[0.00030650, 0.00031870]
Normal-mode MPI sync interval:\t1048576
MPI sync rounds:\t\t2
Final MPI round trajectories:\t1500000
Capture target overshoot:\t3
Total # of scatterings:\t64000000
Numerical failure count:\t10
Simulation time [s]:\t\t120.500000
"""
        timing = """\
\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:03.50
\tMaximum resident set size (kbytes): 456789
"""
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            log_path = root / "run.log"
            time_path = root / "time.txt"
            output_path = root / "metrics.tsv"
            log_path.write_text(log, encoding="utf-8")
            time_path.write_text(timing, encoding="utf-8")
            args = types.SimpleNamespace(
                log=str(log_path),
                time=str(time_path),
                output=str(output_path),
                variant="perf-initial-validation",
                case="m0p5_sigmam40",
                ranks=32,
                repeat=1,
                commit="abc123",
                host="node-a",
                run_status=0,
            )

            benchmark_metrics.extract(args)

            with output_path.open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["wall_seconds"], "123.5")
            self.assertEqual(row["max_rss_kb"], "456789")
            self.assertEqual(row["sync_rounds"], "2")
            self.assertEqual(row["final_round_trajectories"], "1500000")
            self.assertEqual(row["capture_overshoot"], "3")
            self.assertAlmostEqual(float(row["capture_rate_valid"]), 10003.0 / 31999990.0)
            self.assertAlmostEqual(float(row["trajectory_rate_per_second"]), 32000000.0 / 120.5)

    def test_compare_reports_speedup_without_pooling_repeats(self):
        with tempfile.TemporaryDirectory() as directory:
            metrics_path = pathlib.Path(directory) / "metrics.tsv"
            rows = []
            for variant, wall_seconds in (("main", 20.0), ("perf", 10.0)):
                for repeat in (1, 2):
                    row = {field: "" for field in benchmark_metrics.FIELDS}
                    row.update(
                        variant=variant,
                        case="case-a",
                        ranks="4",
                        repeat=str(repeat),
                        host="node-a",
                        run_status="0",
                        wall_seconds=str(wall_seconds),
                        trajectories="1000",
                        valid_trajectories="1000",
                        captured="100",
                        numerical_failures="0",
                    )
                    rows.append(row)
            with metrics_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=benchmark_metrics.FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)

            args = types.SimpleNamespace(
                metrics=[str(metrics_path)],
                baseline="main",
                candidate="perf",
                max_z_score=3.0,
                max_failure_rate=1.0e-4,
                min_speedup=1.5,
                allow_host_mismatch=False,
            )
            output = io.StringIO()
            with redirect_stdout(output):
                result = benchmark_metrics.compare(args)

            self.assertEqual(result, 0)
            self.assertIn("\t2\t", output.getvalue())
            self.assertIn("\tPASS", output.getvalue())


if __name__ == "__main__":
    unittest.main()
