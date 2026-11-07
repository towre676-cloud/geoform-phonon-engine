"""Acquire *only pinned and hash-verified* public Phonopy reference inputs.

No user-supplied URL is accepted; access is restricted to the immutable
Phonopy commit declared in the manifest. Files are small and source-licensed;
keep the upstream BSD-3 notice with redistribution. A local file is always
verified, not trusted merely because its name exists.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.request import urlopen, Request

MANIFEST = Path(__file__).resolve().parents[1] / 'benchmarks' / 'reference_sources.json'
MAX_BYTES = 200_000


def git_blob_sha1(data: bytes):
    """The exact Git object hash, including `blob {n}\0` framing."""
    return hashlib.sha1(f'blob {len(data)}\0'.encode('ascii') + data).hexdigest()


REQUIRED_MATERIALS = {'halite', 'periclase'}
PROVENANCE_FILES = {'POSCAR-unitcell', 'vaspruns.tar.lzma', 'README.md', 'BORN', 'INCAR'}


def material_files(material):
    return {**material['files'], **material.get('provenance_files', {})}


def validate_manifest(record):
    if record.get('schema_version') != 'geoform.phonons.reference_sources/1':
        raise ValueError('unsupported benchmark source manifest')
    if (record.get('upstream_repo') != 'phonopy/phonopy'
            or not re.fullmatch(r'[0-9a-f]{40}', record.get('upstream_commit', ''))):
        raise ValueError('benchmark commit/repo not pinned')

    def validate_file(info, expected_path):
        if (info.get('path') != expected_path
                or not re.fullmatch(r'[0-9a-f]{40}', info.get('git_blob_sha1', ''))
                or not re.fullmatch(r'[0-9a-f]{64}', info.get('sha256', ''))
                or not 0 < info.get('bytes', 0) <= MAX_BYTES):
            raise ValueError(f'untrusted benchmark source: {expected_path}')

    materials = record.get('materials', {})
    if not REQUIRED_MATERIALS <= materials.keys():
        raise ValueError('missing original benchmark materials')
    validate_file(record['license_file'], 'LICENSE')
    for slug, material in materials.items():
        directory = material['source_directory']
        provenance_files = set(material.get('provenance_files', {}))
        if (not re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', slug)
                or not re.fullmatch(r'example/[A-Za-z0-9_-]+', directory)
                or material.get('evidence_status') != 'calculated'
                or set(material.get('files', {})) != {'phonopy_disp.yaml', 'FORCE_SETS'}
                or not provenance_files <= PROVENANCE_FILES
                or ('structure' in material and
                    not {'POSCAR-unitcell', 'vaspruns.tar.lzma'} <= provenance_files)):
            raise ValueError(f'invalid {slug} sources')
        for name, info in material_files(material).items():
            validate_file(info, f'{directory}/{name}')


def matches_pin(data, info):
    return (isinstance(data, bytes) and len(data) == info['bytes']
            and git_blob_sha1(data) == info['git_blob_sha1']
            and hashlib.sha256(data).hexdigest() == info['sha256'])


def fetch_benchmarks(destination, *, downloader=None):
    """Return a verified digest ledger. Redownload modified files."""
    record = json.loads(MANIFEST.read_text(encoding='utf-8'))
    validate_manifest(record)
    dest = Path(destination)
    if downloader is None:
        def downloader(url):
            request = Request(url, headers={'User-Agent':'GeoFormPhononBenchmark/0.3'})
            with urlopen(request, timeout=30) as response:
                data = response.read(MAX_BYTES + 1)
            if len(data) > MAX_BYTES:
                raise ValueError('upstream reference exceeds size cap')
            return data
    def fetch_file(out, info):
        data = out.read_bytes() if out.is_file() and out.stat().st_size <= MAX_BYTES else b''
        if not matches_pin(data, info):
            url = ('https://raw.githubusercontent.com/phonopy/phonopy/'
                   + record['upstream_commit'] + '/' + info['path'])
            data = downloader(url)
            if not matches_pin(data, info):
                raise ValueError(f'upstream hash mismatch for {info["path"]}')
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(data)
        return {key: info[key] for key in ('git_blob_sha1', 'sha256', 'bytes')}

    ledger = {'LICENSE': fetch_file(dest / 'LICENSE', record['license_file'])}
    for slug, material in sorted(record['materials'].items()):
        ledger[slug] = {
            name: fetch_file(dest / slug / name, info)
            for name, info in sorted(material_files(material).items())
        }
    return ledger


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination', default='benchmarks/upstream')
    args=p.parse_args(argv)
    hashes=fetch_benchmarks(args.destination)
    print(json.dumps(hashes,indent=2))

if __name__ == '__main__':main()
