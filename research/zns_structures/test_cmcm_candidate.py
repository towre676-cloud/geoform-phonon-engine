import copy
import json
import unittest

import numpy as np

from geoform_phonons.core import Crystal
from validate_cmcm_candidate import PATH, build_candidate, sensitivity
from validate_structures import replay


class CmcmCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = replay()

    def test_replay_and_engine_crystal(self):
        stored = json.loads(PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, build_candidate(self.manifest))
        cell = Crystal(**stored["crystal"])
        self.assertEqual(cell.symbols, ("Zn", "Zn", "S", "S"))
        self.assertEqual(cell.n_atoms, 4)
        self.assertAlmostEqual(cell.volume_angstrom3, 43.918962808, places=9)
        self.assertEqual(stored["hall_number"], 298)
        np.testing.assert_allclose(
            stored["first_seven_unlike_neighbors_angstrom"],
            [2.219, 2.219, 2.243, 2.243, 2.25115036, 2.25115036, 3.60697496],
            rtol=0, atol=1e-8,
        )
        self.assertTrue(stored["research_only"])
        self.assertIsNone(stored["coordinate_model_temperature_K"])
        self.assertIsNone(self.manifest["phases"]["cmcm"]["sites"])
        self.assertFalse(self.manifest["four_phase_study_ready"])

    def test_half_digit_grid_includes_higher_symmetry_limit(self):
        lengths = self.manifest["phases"]["cmcm"]["cell_lengths_angstrom"]
        rows = sensitivity(lengths, (.7, .2))
        self.assertEqual(len(rows), 25)
        self.assertTrue(all(row[4] < row[5] for row in rows))
        self.assertEqual([(zn, s, group) for zn, s, group, *_ in rows if group != 63],
                         [(.75, .25, 69)])
        self.assertAlmostEqual(min(row[3] for row in rows), 1.7944, places=6)
        self.assertAlmostEqual(max(row[4] for row in rows), 2.6916, places=6)

    def test_candidate_rejects_partial_occupancy_and_promotion(self):
        for field, value in (("occupancies", {"Zn": .9, "S": 1.0}), ("research_only", False)):
            with self.subTest(field=field):
                manifest = copy.deepcopy(self.manifest)
                manifest["research_candidates"]["cmcm_jyoti_96gpa"][field] = value
                with self.assertRaisesRegex(ValueError, "ordered research-only"):
                    build_candidate(manifest)
