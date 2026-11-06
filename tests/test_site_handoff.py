"""Static integration contract for the GeoForm ChatGPT Site handoff."""
from pathlib import Path
import hashlib
import json
import unittest

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
        catalog=json.loads((BASE/'data/catalog.json').read_text())
        receipts=json.loads((BASE/'data/benchmark_receipts.json').read_text())
        self.assertEqual({x['slug'] for x in catalog['entries']},{'halite','periclase'})
        self.assertEqual({x['material'] for x in receipts['results']},{'halite','periclase'})
        for r in receipts['results']:
            self.assertEqual(len(r['qpoint_comparisons']),16)
            self.assertLess(r['max_abs_frequency_error_thz'],3e-6)
            self.assertEqual(r['evidence_status'],'calculated')
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
