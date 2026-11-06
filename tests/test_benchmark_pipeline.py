import unittest
from scripts.benchmark_minerals import Q_POINTS, MAX_DMAT_ATOL, MAX_FREQ_ATOL_THZ
from scripts.fetch_reference_sources import MANIFEST
import json

class PipelineTests(unittest.TestCase):
    def test_q_coverage_contains_gamma_edge_and_general(self):
        self.assertGreaterEqual(len(Q_POINTS),12)
        self.assertIn((0,0,0),Q_POINTS)
        self.assertIn((0.5,0,0),Q_POINTS)
        self.assertIn((0.17,0.31,0.09),Q_POINTS)
        self.assertLessEqual(MAX_DMAT_ATOL,1e-5)
        self.assertLessEqual(MAX_FREQ_ATOL_THZ,1e-3)
    def test_pinned_dataset_provenance(self):
        manifest=json.loads(MANIFEST.read_text())
        self.assertEqual(manifest['materials']['halite']['formula'],'NaCl')
        self.assertEqual(manifest['materials']['periclase']['formula'],'MgO')
        self.assertIn('500 eV',manifest['materials']['periclase']['origin'])
if __name__=='__main__':unittest.main()
