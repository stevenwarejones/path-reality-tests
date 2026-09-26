"""Adversarial mapping and numerical checks; fixtures are not experimental data."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import bisect
from decimal import Decimal, localcontext
import json
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import numpy as np

STUDY = Path(__file__).resolve().parents[1] / "studies" / "wen-2026-propagator"
sys.path.insert(0, str(STUDY))
import replicate as rep


class ReplicationTests(unittest.TestCase):
    def test_sparse_blocks_use_their_own_coordinates(self):
        rows = [("dense_x", "dense_y", "sparse_x", "sparse_y"),
                (0, 10, -1, 20), (1, 11, 1, 21), (2, 12, None, None)]
        np.testing.assert_array_equal(rep.block(rows, "CD", 2), [[-1, 20], [1, 21]])
        with self.assertRaises(ValueError):
            rep.block(rows, "CD", 3)
        with self.assertRaises(ValueError):
            rep.block(rows, "AB", 2)

    def test_alignment_is_coordinate_based_and_never_extrapolates(self):
        np.testing.assert_allclose(rep.align_theory([0, 2, 7], [1, 5, 15], [1, 6]), [3, 13])
        for x, target in [([0, 0, 7], [1]), ([0, 7, 2], [1]), ([0, 2, 7], [8])]:
            with self.assertRaises(ValueError):
                rep.align_theory(x, [1, 5, 15], target)
        with self.assertRaises(ValueError):
            rep.align_theory([0, 1], [1, float("nan")], [0])

    def test_local_cubic_recovers_known_polynomial(self):
        x = np.array([-2., -1., 0., 1., 2.])
        target = np.array([-1.5, .25, 1.5])
        np.testing.assert_allclose(rep.align_theory(x, x ** 3 - x + 4, target, "local_cubic"),
                                   target ** 3 - target + 4, atol=1e-12)

    def test_phase_wrap_does_not_confuse_branch_cut_with_large_error(self):
        np.testing.assert_allclose(rep.wrap_phase_pi([1.98, -1.98, .02, -.02]), [-.02, .02, .02, -.02])

    def test_metric_has_hand_calculated_scale_and_separates_comparisons(self):
        self.assertAlmostEqual(rep.metrics([1., 2.], [1., 1.])["smape_percent"], 100 / 3)
        self.assertAlmostEqual(rep.metrics([1., 2.], [1., 1.])["rmse"], np.sqrt(.5))
        self.assertNotIn("smape_percent", rep.metrics([-1.], [1.], False))
        with self.assertRaises(ValueError):
            rep.metrics([1., 2.], [1.])

    def test_mean_then_metric_is_not_metric_then_mean(self):
        # Counterexample only: no claim these are the authors' repeated measurements.
        trials = np.array([.8, 1.2])
        metric_after_mean = rep.metrics([trials.mean()], [1.])["smape_percent"]
        mean_of_metrics = np.mean([rep.metrics([t], [1.])["smape_percent"] for t in trials])
        self.assertEqual(metric_after_mean, 0)
        self.assertGreater(mean_of_metrics, 20)

    def test_normalization_removes_gain_but_not_shape(self):
        a = np.array([1., 3.]); b = 2 * a
        self.assertGreater(rep.metrics(a, b)["smape_percent"], 0)
        self.assertEqual(rep.metrics(a / a.sum(), b / b.sum())["smape_percent"], 0)

    def test_replication_from_another_directory_matches_reviewed_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run([sys.executable, str(STUDY / "replicate.py"), "--check"],
                                    cwd=tmp, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_input_directory_cannot_be_used_for_output(self):
        result = subprocess.run([sys.executable, str(STUDY / "replicate.py"),
                                 "--output-dir", str(STUDY / "raw" / "dataset")],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)

    def test_figure3_independent_ooxml_decimal_calculation(self):
        # Independent of openpyxl, numpy interpolation and the implementation's SMAPE.
        with ZipFile(STUDY / "raw" / "dataset" / "Fig3.xlsx") as z:
            root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        cells = {}
        for cell in root.findall(".//m:c", ns):
            value = cell.find("m:v", ns)
            if value is not None and cell.get("t") != "s":
                cells[cell.attrib["r"]] = Decimal(value.text)
        def col(c, n):
            return [cells[f"{c}{i}"] for i in range(2, n + 2)]
        x, y, target, e, q, c = col("A", 850), col("B", 850), col("D", 17), col("E", 17), col("F", 17), col("H", 17)
        with localcontext() as context:
            context.prec = 40
            theory = []
            for t in target:
                i = min(max(bisect.bisect_right(x, t) - 1, 0), len(x) - 2)
                theory.append(y[i] + (y[i + 1] - y[i]) * (t - x[i]) / (x[i + 1] - x[i]))
            reviewed = json.loads((STUDY / "results" / "replication" / "summary.json").read_text())["figure3"]["comparisons"]
            for name, a, b in [("E_vs_Q", e, q), ("E_vs_C", e, c), ("Q_vs_deposited_theory", q, theory)]:
                error = sum(2 * abs(u - v) / (abs(u) + abs(v)) for u, v in zip(a, b)) * 100 / Decimal(17)
                self.assertAlmostEqual(float(error), reviewed[name]["smape_percent"], places=11)


if __name__ == "__main__":
    unittest.main()
