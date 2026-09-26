"""Behavioral tests; tiny workbooks below are synthetic test fixtures only."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import urllib.error
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import openpyxl

STUDY = Path(__file__).resolve().parents[1] / "studies" / "wen-2026-propagator"
sys.path.insert(0, str(STUDY))
import fetch_data as fetch
import extract_workbooks as extract
import synthetic_checks as synthetic


def entry(name, data):
    return {"name": name, "size": len(data), "sha256": hashlib.sha256(data).hexdigest(),
            "url": "https://example.invalid/" + name}


class VerificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = self.root / "dataset"
        self.data.mkdir()
        self.original = b"fixture-only"
        self.entries = [entry("fixture.bin", self.original)]

    def run_cli(self, verify_only=True, response=None):
        (self.root / "manifest.json").write_text(json.dumps({"files": self.entries}))
        argv = ["fetch_data.py", "--directory", str(self.data)]
        if verify_only:
            argv.append("--verify-only")
        with mock.patch.object(fetch, "ROOT", self.root), mock.patch.object(sys, "argv", argv), \
                mock.patch.object(fetch.urllib.request, "urlopen", side_effect=response) as network, \
                contextlib.redirect_stdout(io.StringIO()):
            code = fetch.main()
        return code, network

    def test_verified_copy_is_unchanged_and_offline(self):
        file = self.data / "fixture.bin"
        file.write_bytes(self.original)
        before = file.stat().st_mtime_ns
        code, network = self.run_cli()
        self.assertEqual(code, 0)
        network.assert_not_called()
        self.assertEqual(file.read_bytes(), self.original)
        self.assertEqual(file.stat().st_mtime_ns, before)

    def test_missing_file_fails_without_download_or_repair(self):
        code, network = self.run_cli()
        self.assertEqual(code, 2)
        network.assert_not_called()
        self.assertFalse((self.data / "fixture.bin").exists())

    def test_same_size_corruption_fails_without_repair(self):
        damaged = b"X" + self.original[1:]
        file = self.data / "fixture.bin"
        file.write_bytes(damaged)
        code, network = self.run_cli()
        self.assertEqual(code, 2)
        network.assert_not_called()
        self.assertEqual(file.read_bytes(), damaged)

    def test_wrong_size_fails(self):
        (self.data / "fixture.bin").write_bytes(self.original + b"extra")
        self.assertEqual(self.run_cli()[0], 2)

    def test_unlisted_file_fails(self):
        (self.data / "fixture.bin").write_bytes(self.original)
        (self.data / "unexpected.pdf").write_bytes(b"fixture")
        self.assertEqual(self.run_cli()[0], 2)

    def test_access_control_stops_further_requests(self):
        self.entries.append(entry("second.bin", b"second"))
        error = urllib.error.HTTPError("https://example.invalid/fixture", 403, "Forbidden", {}, None)
        code, network = self.run_cli(verify_only=False, response=error)
        self.assertEqual(code, 2)
        self.assertEqual(network.call_count, 1)
        self.assertEqual(list(self.data.iterdir()), [])

    def test_html_response_is_not_saved_as_original(self):
        self.entries.append(entry("second.bin", b"second"))
        code, network = self.run_cli(verify_only=False,
                                    response=lambda *a, **k: io.BytesIO(b"<html>Denied</html>"))
        self.assertEqual(code, 2)
        self.assertEqual(network.call_count, 1)
        self.assertEqual(list(self.data.iterdir()), [])


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = self.root / "dataset"
        self.data.mkdir()
        self.file = self.data / "fixture.xlsx"
        book = openpyxl.Workbook()
        book.active["A1"] = "fixture"
        book.active["B1"] = "=1+1"
        hidden = book.create_sheet("hidden fixture")
        hidden.sheet_state = "hidden"
        hidden["C3"] = 42
        book.save(self.file)
        book.close()
        # Supply a known cached formula result without relying on an Excel engine.
        with zipfile.ZipFile(self.file) as z:
            parts = {name: z.read(name) for name in z.namelist()}
        ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        xml = ET.fromstring(parts["xl/worksheets/sheet1.xml"])
        cell = xml.find(".//s:c[@r='B1']/s:v", ns)
        self.assertIsNotNone(cell)
        cell.text = "2"
        parts["xl/worksheets/sheet1.xml"] = ET.tostring(xml)
        with zipfile.ZipFile(self.file, "w", zipfile.ZIP_DEFLATED) as z:
            for name, data in parts.items():
                z.writestr(name, data)
        self.original = self.file.read_bytes()
        self.entries = [entry(self.file.name, self.original)]
        self.out = self.root / "output"

    def run_extraction(self):
        (self.root / "manifest.json").write_text(json.dumps({"files": self.entries}))
        argv = ["extract_workbooks.py", "--directory", str(self.data), "--output", str(self.out)]
        with mock.patch.object(extract, "ROOT", self.root), mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stdout(io.StringIO()):
            extract.main()

    def test_hidden_sheet_formula_cache_and_original_preserved(self):
        self.run_extraction()
        report = json.loads((self.out / "fixture.json").read_text())
        self.assertEqual(report["source"]["sha256"], hashlib.sha256(self.original).hexdigest())
        self.assertEqual([s["name"] for s in report["sheets"]], ["Sheet", "hidden fixture"])
        self.assertEqual(report["sheets"][1]["state"], "hidden")
        self.assertEqual(report["sheets"][1]["cells"][0]["coordinate"], "C3")
        formula = next(c for c in report["sheets"][0]["cells"] if c["coordinate"] == "B1")
        self.assertEqual(formula["value_or_formula"], "=1+1")
        self.assertEqual(formula["cached"], 2)
        self.assertEqual(self.file.read_bytes(), self.original)

    def test_entire_set_is_verified_before_any_extract(self):
        self.entries.append(entry("missing.xlsx", b"fixture-only"))
        with self.assertRaises(FileNotFoundError):
            self.run_extraction()
        self.assertFalse(self.out.exists())


class NumericalTests(unittest.TestCase):
    def test_interference_and_incoherent_sum_are_distinct(self):
        h = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        paths, _, amplitudes = synthetic.enumerate_amplitudes(np.stack([h, h]), 0)
        coherent = [abs(amplitudes[paths[:, -1] == i].sum()) ** 2 for i in range(2)]
        incoherent = [sum(abs(amplitudes[paths[:, -1] == i]) ** 2) for i in range(2)]
        np.testing.assert_allclose(coherent, [1, 0], atol=1e-12)
        np.testing.assert_allclose(incoherent, [.5, .5], atol=1e-12)

    def test_complex_matrix_composition(self):
        self.assertLess(synthetic.check_composition(), 1e-12)

    def test_smape_zero_and_one_sided_zero(self):
        self.assertEqual(synthetic.smape(np.array([0., 0.]), np.array([0., 0.])), 0)
        self.assertEqual(synthetic.smape(np.array([0., 2.]), np.array([2., 0.])), 2)

    def test_snapshot_rejects_changed_value_count_and_nonfinite(self):
        expected = {"n": 17, "values": [.5, .5]}
        synthetic.compare_snapshot({"n": 17, "values": [.5 + 1e-14, .5]}, expected)
        for actual in ({"n": 18, "values": [.5, .5]}, {"n": 17, "values": [.5]},
                       {"n": 17, "values": [.6, .4]}, {"n": 17, "values": [float('nan'), .5]}):
            with self.subTest(actual=actual), self.assertRaises(ValueError):
                synthetic.compare_snapshot(actual, expected)


class BundledDatasetTests(unittest.TestCase):
    def test_all_original_hashes_and_cli_from_another_directory(self):
        # A changed manifest must not silently redefine the audited deposit.
        self.assertEqual(hashlib.sha256((STUDY / "manifest.json").read_bytes()).hexdigest(),
                         "8d6e3bdf538bbbd124935d71c2f68096e7d9f98ecefefed96c575e109ea3b135")
        entries = json.loads((STUDY / "manifest.json").read_text())["files"]
        self.assertEqual(len(entries), 19)
        for e in entries:
            fetch.verify(STUDY / "raw" / "dataset" / e["name"], e)
        with tempfile.TemporaryDirectory() as cwd:
            result = subprocess.run([sys.executable, str(STUDY / "fetch_data.py"), "--verify-only"],
                                    cwd=cwd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.count(" verified"), 19)


if __name__ == "__main__":
    unittest.main()
