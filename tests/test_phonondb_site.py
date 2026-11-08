import importlib.util
from types import SimpleNamespace
import unittest

import numpy as np

if importlib.util.find_spec('phonopy') and importlib.util.find_spec('spglib'):
    from scripts import phonondb_site as site
else:
    site = None


@unittest.skipUnless(site, 'requires the optional PhononDB dependencies')
class PhononDBSiteTests(unittest.TestCase):
    def imaginary(self):
        q = np.concatenate([np.linspace(a, b, 12, endpoint=i == len(site.VERTICES)-2)
                            for i, (a, b) in enumerate(zip(site.VERTICES[:-1], site.VERTICES[1:]))])
        x = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0)*2*np.pi, axis=1))]
        comparison = dict(matrix_q_points=16, frequency_q_points=280,
                          max_abs_matrix_error_ev_a2_amu=0, max_abs_frequency_error_thz=0,
                          signed_minimum_thz=-1, reference_signed_minimum_thz=-1,
                          imaginary_sampled_modes=1, gamma_acoustic_roundoff_modes_excluded=0,
                          minimum_q_fractional=[0, 0, 0], stability='sampled_imaginary_modes')
        result = dict(crystal=dict(lattice_angstrom=np.eye(3).tolist()),
                      dispersion=dict(q_fractional=q.tolist(), distance_inverse_angstrom=x.tolist(),
                                      frequency_thz_signed=[[-1, 0, 2]]+[[0, 1, 2]]*47),
                      physics=dict(q_mesh=[6]*3, sampled_stability=comparison,
                                   nonanalytic_correction='omitted on both sides'),
                      density_of_states={'status': 'withheld: sampled_imaginary_modes'},
                      thermodynamics={'status': 'withheld: sampled_imaginary_modes'})
        return result, comparison

    def test_signed_payload_is_accepted_with_thermal_withheld(self):
        result, comparison = self.imaginary()
        self.assertFalse(site.verify_outputs(result, comparison, 1, 6))
        result['thermodynamic_curve'] = [{'temperature_k': 0}]
        with self.assertRaisesRegex(ValueError, 'Unexpected DOS or thermal'):
            site.verify_outputs(result, comparison, 1, 6)

    def test_invalid_numeric_receipt_is_rejected(self):
        result, comparison = self.imaginary()
        comparison['max_abs_matrix_error_ev_a2_amu'] = float('nan')
        with self.assertRaisesRegex(ValueError, 'independent comparison'):
            site.verify_outputs(result, comparison, 1, 6)

    def test_distance_must_match_archived_lattice(self):
        result, comparison = self.imaginary()
        result['crystal']['lattice_angstrom'][0][0] = 2
        with self.assertRaisesRegex(ValueError, 'signed dispersion'):
            site.verify_outputs(result, comparison, 1, 6)

    def test_license_conflict_never_falls_back_to_collection(self):
        row = {'mp_id': 'mp-123'}
        self.assertIsNone(site.record_rights({'status': 'record_license_review'}, row, True))
        self.assertIsNone(site.record_rights({'status': 'collection_license_only'}, row, False))
        self.assertEqual(site.record_rights({'status': 'collection_license_only'}, row, True),
                         'collection_license_only')
        self.assertIsNone(site.record_rights({'status': 'record_license_unconfirmed',
                                             'error': 'ValueError: wrong source record'}, row, True))
        with self.assertRaisesRegex(ValueError, 'rights receipt'):
            site.record_rights({'status': 'record_cc_by_4_0', 'rights_ids': ['cc-by-4.0'],
                                'page_sha256': '1'*64, 'record_title': 'materials id 124'}, row, True)

    def test_mapping_requires_correct_ordered_motif_and_reviewed_phase_id(self):
        row = {'record_id': 'gf06g7088'}
        structure = {'primitive_atoms': 2, 'reduced_composition': {'Na': 1, 'Cl': 1}}
        fcc = np.array([(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)])
        symmetry = SimpleNamespace(number=225, std_lattice=np.eye(3),
                                   std_types=np.array([11]*4+[17]*4),
                                   std_positions=np.vstack((fcc, (fcc+.5)%1)))
        self.assertEqual(site.mineral_mapping(row, structure, symmetry)['name'], 'Halite')
        self.assertIsNone(site.mineral_mapping({'record_id': 'other'}, structure, symmetry))
        symmetry.std_positions[4] = [.25, .25, .25]
        with self.assertRaisesRegex(ValueError, 'motif mismatch'):
            site.mineral_mapping(row, structure, symmetry)
