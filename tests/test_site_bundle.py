import json
import shutil
from pathlib import Path
import tempfile
import unittest
from scripts.build_site_bundle import HANDOFF, verify_results

class SiteBundleTests(unittest.TestCase):
    def test_empty_or_missing_artifacts_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):verify_results(d)
    def test_false_benchmark_status_not_promoted(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'data'
            shutil.copytree(HANDOFF / 'data', path)
            catalog = path / 'catalog.json'
            record = json.loads(catalog.read_text())
            record['entries'][0]['evidence_status'] = 'measured'
            catalog.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, 'invalid catalog identity'):
                verify_results(path)

    def test_unlinked_result_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'data'
            shutil.copytree(HANDOFF / 'data', path)
            result = path / 'anatase_phonons.json'
            result.write_bytes(result.read_bytes() + b' ')
            with self.assertRaisesRegex(ValueError, 'inconsistent computational evidence'):
                verify_results(path)


if __name__ == '__main__':
    unittest.main()
