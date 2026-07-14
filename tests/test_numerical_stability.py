import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import evaporation
import evaporation_theory as theory


class NumericalStabilityTests(unittest.TestCase):
    def test_solar_model_loads_from_repository(self):
        data = theory.load_solar_model()
        self.assertGreater(data.shape[0], 100)
        self.assertGreaterEqual(data.shape[1], 20)
        self.assertTrue(np.all(np.isfinite(data)))

    def test_capture_saturation_uses_stable_exponential_difference(self):
        with mock.patch.object(theory, "capture_rate_weak_SD", return_value=1.0e300), \
             mock.patch.object(theory, "capture_rate_geom", return_value=1.0):
            rate = theory.capture_rate_SD(1.0, 1.0, {})
        self.assertTrue(np.isfinite(rate))
        self.assertAlmostEqual(rate, 1.0, places=14)

    def test_population_solution_is_stable_in_limiting_regimes(self):
        self.assertAlmostEqual(theory.N_total(4.0, 1.0, 0.0), 2.0, places=14)
        self.assertAlmostEqual(theory.N_total(4.0, 1.0, 3.0), 1.0, places=14)
        with self.assertRaises(ValueError):
            theory.N_total(1.0, 1.0, np.nan)

    def test_bound_time_analysis_returns_the_computed_average(self):
        with tempfile.TemporaryDirectory() as directory:
            for index, duration in enumerate((3.0, 5.0)):
                rows = np.zeros((2, evaporation.N_COLS), dtype=np.float64)
                rows[:, 0] = (0.0, duration)
                rows[:, 7] = (-1.0, -0.5)
                rows.tofile(os.path.join(directory, f"trajectory_{index}.dat"))

            average = evaporation.analyze_trajectory_folder(
                directory, max_samples=10, num_processes=1
            )
        self.assertAlmostEqual(average, 4.0, places=14)


if __name__ == "__main__":
    unittest.main()
