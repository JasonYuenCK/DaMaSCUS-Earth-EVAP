import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from validation.physics_validation import (
    compare_runs,
    load_run_metrics,
    rewrite_config,
)


BINCOUNT_TEMPLATE = """# valid_trajectories = 100
# captured_particles = 50
# capture_rate_valid = {capture_rate}
# capture_rate_valid_err = 0.04
# capture_rate_valid_CI_95_lower = {ci_lower}
# capture_rate_valid_CI_95_upper = {ci_upper}
# numerical_failure_rate = {failure_rate}
# average_scatterings = {average_scatterings}
# bin_index cap_dt other
0 {bin0} 0
1 {bin1} 0
"""


def write_result(directory, capture_rate=0.5, ci_lower=0.4, ci_upper=0.6,
                 failure_rate=0.0, average_scatterings=10.0,
                 bins=(6.0, 4.0), lifetimes=range(1, 41)):
    directory.mkdir(parents=True)
    (directory / "bincount.txt").write_text(
        BINCOUNT_TEMPLATE.format(
            capture_rate=capture_rate,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            failure_rate=failure_rate,
            average_scatterings=average_scatterings,
            bin0=bins[0],
            bin1=bins[1],
        ),
        encoding="utf-8",
    )
    with (directory / "evaporation_times.txt").open("w", encoding="utf-8") as handle:
        handle.write("# rank trajectory_id lifetime_unbinding_sec\n")
        for index, lifetime in enumerate(lifetimes):
            handle.write(f"0 {index} {lifetime}\n")


class PhysicsValidationWorkflowTests(unittest.TestCase):
    def test_rewrite_config_replaces_and_adds_settings(self):
        source = 'run_mode = "Capture";\nsample_size = 10; // keep\n'
        rewritten = rewrite_config(
            source,
            {
                "run_mode": "Parameter point",
                "sample_size": 100,
                "fixed_seed": 42,
                "snapshot_enabled": False,
            },
        )
        self.assertIn('run_mode = "Parameter point";', rewritten)
        self.assertIn("sample_size = 100; // keep", rewritten)
        self.assertIn("fixed_seed = 42;", rewritten)
        self.assertIn("snapshot_enabled = false;", rewritten)

    def test_load_metrics_normalizes_radial_histogram(self):
        with tempfile.TemporaryDirectory() as temp:
            result = Path(temp) / "result"
            write_result(result)
            metrics = load_run_metrics(result, grid=0, seed=7)
        self.assertEqual(metrics.valid_trajectories, 100)
        self.assertEqual(metrics.captured_particles, 50)
        self.assertEqual(metrics.evaporation_events, 40)
        self.assertAlmostEqual(metrics.complete_evaporation_fraction, 0.8)
        self.assertAlmostEqual(metrics.evaporation_median_sec, 20.5)
        self.assertEqual(metrics.captured_radial_distribution, [0.6, 0.4])

    def test_compatible_runs_pass_physics_gates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline_dir = root / "baseline"
            candidate_dir = root / "candidate"
            write_result(baseline_dir)
            write_result(
                candidate_dir,
                capture_rate=0.52,
                ci_lower=0.42,
                ci_upper=0.62,
                average_scatterings=10.5,
                bins=(5.8, 4.2),
                lifetimes=range(2, 42),
            )
            baseline = load_run_metrics(baseline_dir, 0, 7)
            candidate = load_run_metrics(candidate_dir, 1000, 7)
            comparison = compare_runs(baseline, candidate)
        self.assertEqual(comparison["status"], "pass")

    def test_disjoint_capture_intervals_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline_dir = root / "baseline"
            candidate_dir = root / "candidate"
            write_result(baseline_dir, ci_lower=0.4, ci_upper=0.6)
            write_result(candidate_dir, capture_rate=0.8, ci_lower=0.72, ci_upper=0.87)
            baseline = load_run_metrics(baseline_dir, 0, 7)
            candidate = load_run_metrics(candidate_dir, 1000, 7)
            comparison = compare_runs(baseline, candidate)
        self.assertEqual(comparison["status"], "fail")
        failed = [item["name"] for item in comparison["checks"] if item["status"] == "fail"]
        self.assertIn("capture_95pct_intervals_overlap", failed)


if __name__ == "__main__":
    unittest.main()
