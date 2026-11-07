"""Static integration contract for the GeoForm ChatGPT Site handoff."""
from pathlib import Path
import hashlib
import json
import unittest

import numpy as np

from scripts.build_site_bundle import verify_results

BASE = Path(__file__).resolve().parents[1] / 'site' / 'geoform-phonons'

class SiteHandoffTests(unittest.TestCase):
    def test_integrity(self):
        manifest=json.loads((BASE/'site-handoff.json').read_text())
        self.assertEqual(manifest['schema_version'],'geoform.phonons.site_handoff/1')
        for name, meta in manifest['assets'].items():
            b=(BASE/name).read_bytes()
            self.assertEqual(len(b),meta['bytes'])
            self.assertEqual(hashlib.sha256(b).hexdigest(),meta['sha256'])
    def test_supported_cases_only(self):
        entries = verify_results(BASE / 'data')
        self.assertEqual({x['slug'] for x in entries}, {
            'anatase', 'cassiterite', 'corundum', 'cubic-catio3', 'halite', 'periclase', 'zincite',
            'stishovite', 'hexagonal-mgb2'})
        manifest = json.loads((BASE / 'site-handoff.json').read_text())
        self.assertNotIn('cubic-catio3', manifest['supported_minerals'])
        self.assertEqual(manifest['computational_phases_without_mineral_link'], ['cubic-catio3', 'hexagonal-mgb2'])

    def test_source_audit_and_complex_comparison(self):
        receipts = json.loads((BASE / 'data/benchmark_receipts.json').read_text())['results']
        audits = [receipt for receipt in receipts if receipt['source_audit']]
        self.assertEqual(len(audits), 7)
        self.assertEqual(sum(len(r['source_audit']['vasp_force_runs']) for r in audits), 27)
        for receipt in audits:
            with self.subTest(material=receipt['material']):
                audit = receipt['source_audit']
                if receipt['material'] not in {'cubic-catio3', 'hexagonal-mgb2'}:
                    self.assertGreater(max(q['max_abs_imaginary_matrix_entry']
                                           for q in receipt['qpoint_comparisons']), 1e-4)
                self.assertEqual(len(receipt['full_ifc_qpoint_comparisons']), 16)
                self.assertLess(receipt['max_abs_frequency_error_thz'], 1e-5)
                self.assertTrue(all(run['max_force_error_ev_angstrom'] == 0
                                    for run in audit['vasp_force_runs']))
        anatase = next(r for r in audits if r['material'] == 'anatase')['source_audit']
        self.assertEqual(anatase['supercell_matrix'], [[4, 0, 0], [0, 4, 0], [0, 0, 1]])
        self.assertEqual((anatase['source_yaml_basis_atoms'], anatase['primitive_atoms']), (12, 6))

    def test_unstable_case_withholds_thermal_and_dos(self):
        record = json.loads((BASE / 'data/cubic-catio3_phonons.json').read_text())
        self.assertIsNone(record['crystal']['mineral_name'])
        self.assertLess(np.min(record['dispersion']['frequency_thz_signed']), -5.0)
        self.assertIn('withheld', record['density_of_states']['status'])
        self.assertIn('withheld', record['thermodynamics']['status'])
        self.assertNotIn('thermodynamic_curve', record)

    def test_zincite_path_instability_is_not_hidden_by_stable_mesh(self):
        record = json.loads((BASE / 'data/zincite_phonons.json').read_text())
        self.assertEqual(record['physics']['unstable_modes_on_shifted_mesh'], 0)
        self.assertEqual(record['physics']['sampled_stability']['dispersion']['imaginary_modes'], 1)
        self.assertIn('withheld', record['density_of_states']['status'])
        self.assertIn('withheld', record['thermodynamics']['status'])
        receipts = json.loads((BASE / 'data/benchmark_receipts.json').read_text())['results']
        reference = next(r for r in receipts if r['material'] == 'zincite')['dispersion_reference']
        self.assertLess(reference['min_frequency_thz'], -0.14)
        self.assertLess(reference['max_abs_frequency_error_thz'], 1e-5)

    def test_new_phase_provenance_and_stishovite_instability(self):
        receipts = json.loads((BASE / 'data/benchmark_receipts.json').read_text())['results']
        records = {r['material']: r for r in receipts}
        silica = records['stishovite']
        self.assertFalse(silica['dos_thermal_available'])
        self.assertEqual(silica['sampled_stability']['shifted_mesh']['imaginary_modes'], 0)
        self.assertLess(silica['dispersion_reference']['min_frequency_thz'], -0.49)
        self.assertTrue(all(run['external_pressure_kbar'] is None
                            for run in silica['source_audit']['vasp_force_runs']))
        mgb2 = records['hexagonal-mgb2']
        self.assertIsNone(mgb2['name'])
        self.assertTrue(mgb2['dos_thermal_available'])
        self.assertTrue(all(float(run['smearing_width_ev']) == 0.2
                            for run in mgb2['source_audit']['vasp_force_runs']))

    def test_stable_dos_and_thermal_contract(self):
        for item in verify_results(BASE / 'data'):
            if item['dos_thermal_available']:
                with self.subTest(material=item['slug']):
                    record = json.loads((BASE / 'data' / item['result_file']).read_text())
                    dos = record['density_of_states']
                    projected = np.asarray(dos['atom_projected_states_per_thz_cell'])
                    np.testing.assert_allclose(projected.sum(axis=0), dos['states_per_thz_cell'])
                    np.testing.assert_allclose(
                        np.trapezoid(projected, dos['frequency_thz'], axis=1), 3, atol=1e-4)
                    self.assertEqual(record['thermodynamics']['temperature_k'], 300)
                    self.assertGreater(record['thermodynamics']['cv_ev_k_cell'], 0)

    def test_viewer_deep_links(self):
        html=(BASE/'index.html').read_text()
        self.assertIn("get('mineral')",html)
        self.assertIn("entry.slug === requestedMineral",html)
        self.assertIn("new Event('change')",html)
        self.assertIn('data/catalog.json',html)
    def test_embed_contract(self):
        js=(BASE/'geoform-phonon-card.js').read_text()
        self.assertIn('mountCard',js)
        self.assertIn('hrefForMineral',js)
        self.assertIn('UNRESOLVED',js)
        self.assertIn('CALCULATED',js)

if __name__ == '__main__': unittest.main()
