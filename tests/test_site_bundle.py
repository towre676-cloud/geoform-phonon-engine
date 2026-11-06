import json
from pathlib import Path
import tempfile
import unittest
from scripts.build_site_bundle import verify_results

class SiteBundleTests(unittest.TestCase):
    def test_empty_or_missing_artifacts_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):verify_results(d)
    def test_false_benchmark_status_not_promoted(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)
            (path/'catalog.json').write_text(json.dumps({'schema_version':'geoform.phonons.catalog/1','entries':[
                {'slug':'halite','evidence_status':'measured','result_file':'halite_phonons.json'},
                {'slug':'periclase','evidence_status':'calculated','result_file':'periclase_phonons.json'}]}))
            (path/'benchmark_receipts.json').write_text(json.dumps({'schema_version':'geoform.phonons.benchmark_receipts/1','results':[
                {'material':'halite'},{'material':'periclase'}]}))
            with self.assertRaises(ValueError):verify_results(d)

if __name__=='__main__':unittest.main()
