"""Offline reproducibility and invalid-input handling for pinned sources."""
import json
from pathlib import Path
import tempfile
import unittest
from scripts.fetch_reference_sources import git_blob_sha1, validate_manifest, fetch_benchmarks, MANIFEST


class SourceTests(unittest.TestCase):
    def test_git_blob_framing(self):
        self.assertEqual(git_blob_sha1(b'hello world\n'), '3b18e512dba79e4c8300dd08aeb37f8e728b8dad')

    def test_reference_manifest_is_pinned(self):
        data=json.loads(MANIFEST.read_text())
        validate_manifest(data)
        self.assertEqual(set(data['materials']), {
            'halite', 'periclase', 'corundum', 'anatase', 'zincite', 'cubic-catio3', 'cassiterite',
            'stishovite', 'hexagonal-mgb2'})
        self.assertEqual(len({x['git_blob_sha1'] for m in data['materials'].values() for x in m['files'].values()}),18)

    def test_unexpected_source_rejected(self):
        r=json.loads(MANIFEST.read_text())
        r['upstream_repo']='random-unknown/repo'
        with self.assertRaises(ValueError):validate_manifest(r)
        r['upstream_repo']='phonopy/phonopy'
        r['materials']['halite']['files']['FORCE_SETS']['path']='../../elsewhere'
        with self.assertRaises(ValueError):validate_manifest(r)

    def test_hash_failure_does_not_write_unverified_file(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):fetch_benchmarks(d,downloader=lambda _:b'wrong')
            self.assertFalse(list(Path(d).rglob('*')))

if __name__=='__main__':unittest.main()
