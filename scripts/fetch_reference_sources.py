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
from pathlib import Path
from urllib.request import urlopen, Request

MANIFEST = Path(__file__).resolve().parents[1] / 'benchmarks' / 'reference_sources.json'
MAX_BYTES = 200_000


def git_blob_sha1(data: bytes):
    """The exact Git object hash, including `blob {n}\0` framing."""
    return hashlib.sha1(f'blob {len(data)}\0'.encode('ascii') + data).hexdigest()


def validate_manifest(record):
    if record.get('schema_version') != 'geoform.phonons.reference_sources/1':
        raise ValueError('unsupported benchmark source manifest')
    commit = record.get('upstream_commit', '')
    if (record.get('upstream_repo') != 'phonopy/phonopy' or len(commit) != 40
            or not all(c in '0123456789abcdef' for c in commit)):
        raise ValueError('benchmark commit/repo not pinned')
    if set(record.get('materials', {})) != {'halite', 'periclase'}:
        raise ValueError('unexpected benchmark material inventory')
    for name, mat in record['materials'].items():
        if mat.get('evidence_status') != 'calculated' or set(mat.get('files', {})) != {'phonopy_disp.yaml', 'FORCE_SETS'}:
            raise ValueError(f'invalid {name} sources')
        for info in mat['files'].values():
            if (not info['path'].startswith('example/') or '..' in info['path']
                    or len(info.get('git_blob_sha1','')) != 40):
                raise ValueError('untrusted benchmark source location')


def fetch_benchmarks(destination, *, downloader=None):
    """Return verified per-file SHA256 ledger. Redownload modified files."""
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
    ledger = {}
    for slug, material in sorted(record['materials'].items()):
        ledger[slug] = {}
        for basename, info in sorted(material['files'].items()):
            out = dest / slug / basename
            data = out.read_bytes() if out.is_file() else b''
            if git_blob_sha1(data) != info['git_blob_sha1']:
                url = ('https://raw.githubusercontent.com/phonopy/phonopy/'
                       + record['upstream_commit'] + '/' + info['path'])
                data = downloader(url)
                if not isinstance(data, bytes) or len(data) > MAX_BYTES or git_blob_sha1(data) != info['git_blob_sha1']:
                    raise ValueError(f'upstream hash mismatch for {slug}/{basename}')
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(data)
            ledger[slug][basename] = {'git_blob_sha1': info['git_blob_sha1'],
                                      'sha256': hashlib.sha256(data).hexdigest(),
                                      'bytes': len(data)}
    return ledger


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination', default='benchmarks/upstream')
    args=p.parse_args(argv)
    hashes=fetch_benchmarks(args.destination)
    print(json.dumps(hashes,indent=2))

if __name__ == '__main__':main()
