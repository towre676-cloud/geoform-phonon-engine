"""Source index and checksum-locked acquisition for the NIMS MDR PhononDB."""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from http.client import HTTPException
from pathlib import Path
import re
from time import monotonic
from urllib.request import urlopen
from uuid import uuid4


INDEX_PATH = 'mdr/phonondb/README.md'
INDEX_ROOT = 'https://raw.githubusercontent.com/atztogo/phonondb'
MDR = 'https://mdr.nims.go.jp'
COLLECTION = 'https://doi.org/10.48505/nims.4197'
LICENSE = 'https://creativecommons.org/licenses/by/4.0/'
RECORD_ID = re.compile(r'[a-z0-9]{9}')
MAX_ARCHIVE_BYTES = 32 * 1024**2
ROW = re.compile(
    r'\|\s*(\d+)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\((\d+)\)\s*\|\s*'
    + re.escape(MDR)
    + r'/concern/datasets/([a-z0-9]{9})\s*\|\s*\[\5\.zip\]\('
    + re.escape(MDR)
    + r'/download_all/\5\.zip\)\s*\|\s*'
)
SHA256 = re.compile(r'[0-9a-f]{64}')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{uuid4().hex}.tmp')
    try:
        tmp.write_bytes(data)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def write_json(path, value):
    atomic_write(path, (json.dumps(value, sort_keys=True, allow_nan=False) + '\n').encode())


def parse_index(data, revision):
    if not re.fullmatch(r'[0-9a-f]{40}', revision):
        raise ValueError('Use the full Git commit of the upstream index')
    rows = []
    for line in data.decode().splitlines():
        if not line.startswith('|'):
            continue
        if line.startswith(('| Materials Project ID', '| ---')):
            continue
        match = ROW.fullmatch(line)
        if not match:
            raise ValueError(f'Unrecognized upstream index row: {line[:180]}')
        mp, formula, symbol, number, record = match.groups()
        if not 1 <= int(number) <= 230:
            raise ValueError(f'Invalid space group for {record}')
        rows.append({
            'mp_id': f'mp-{mp}', 'record_id': record, 'formula': formula,
            'space_group_symbol': symbol, 'space_group_number': int(number),
            'record_url': f'{MDR}/concern/datasets/{record}',
            'archive_url': f'{MDR}/download_all/{record}.zip',
        })
    if not rows or len({row['record_id'] for row in rows}) != len(rows):
        raise ValueError('Empty index or repeated MDR record ID')
    mp_counts = Counter(row['mp_id'] for row in rows)
    return {
        'schema': 'geoform.phonondb.index/1',
        'index_url': f'{INDEX_ROOT}/{revision}/{INDEX_PATH}',
        'index_sha256': digest(data),
        'collection_doi': COLLECTION, 'collection_license_url': LICENSE,
        'attribution': 'Atsushi Togo, National Institute for Materials Science, MDR phonon calculation database',
        'repeated_mp_ids': sorted(mp for mp, count in mp_counts.items() if count > 1),
        'rows': rows,
    }


def locked_archive(row, cache, *, acquire=False, expected_pin=None, timeout=45):
    """Verify previously pinned bytes or explicitly acquire and pin new bytes."""
    record = row['record_id']
    archive = cache / f'{record}.zip'
    lock = cache / f'{record}.json'
    if lock.exists():
        pin = json.loads(lock.read_text())
        if expected_pin is not None and pin != expected_pin:
            raise ValueError('Local pin differs from supplied pin')
    elif acquire:
        pin = expected_pin
    else:
        raise ValueError('No source checksum pin. Run fetch first.')
    if pin is not None and (pin['url'] != row['archive_url'] or not SHA256.fullmatch(pin['sha256'])
                            or not 0 < pin['bytes'] <= MAX_ARCHIVE_BYTES):
        raise ValueError('Invalid archive pin or source URL mismatch')
    if archive.exists() and pin is not None:
        if archive.stat().st_size > MAX_ARCHIVE_BYTES:
            raise ValueError('Cached archive exceeds size limit')
        data = archive.read_bytes()
    elif acquire:
        with urlopen(row['archive_url'], timeout=timeout) as response:
            data = response.read(MAX_ARCHIVE_BYTES + 1)
        if len(data) > MAX_ARCHIVE_BYTES or not data.startswith(b'PK\x03\x04'):
            raise ValueError('Oversized download or response is not ZIP')
    else:
        raise ValueError('Missing pinned archive. Run fetch first.')
    observed = digest(data)
    if pin is not None:
        if observed != pin['sha256'] or len(data) != pin['bytes']:
            raise ValueError('Archive differs from pinned bytes')
    else:
        pin = {'url': row['archive_url'], 'sha256': observed, 'bytes': len(data),
               'checksum_authority': 'first successful HTTPS retrieval; no publisher checksum assumed'}
    if acquire and not (archive.exists() and lock.exists()):
        atomic_write(archive, data)
        write_json(lock, pin)
    return data, pin


def fetch_batch(rows, cache, workers, timeout, *, pins=None):
    def one(row):
        start = monotonic()
        try:
            expected = pins[row['record_id']] if pins is not None else None
            _, pin = locked_archive(row, cache, acquire=True, expected_pin=expected, timeout=timeout)
            result = {'status': 'checksum_pinned', 'archive': pin}
        except (OSError, HTTPException, ValueError, KeyError, TypeError) as error:
            result = {'status': 'fetch_failed', 'error': f'{type(error).__name__}: {error}'}
        return {**row, **result, 'seconds': round(monotonic() - start, 3)}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        yield from pool.map(one, rows)


def inspect_record_license(row, *, timeout=45):
    """Read the record's own dataset JSON-LD, separate from collection rights."""
    with urlopen(row['record_url'], timeout=timeout) as response:
        page = response.read(2 * 1024**2 + 1)
    if len(page) > 2 * 1024**2:
        raise ValueError('MDR record page exceeds size budget')
    scripts = re.findall(r'<script\b[^>]*type="application/ld\+json"[^>]*>(.*?)</script>',
                         page.decode(), flags=re.S)
    for script in scripts:
        cleaned = script.strip().removeprefix('//<![CDATA[').removesuffix('//]]>').strip()
        for item in json.loads(cleaned).get('@graph', []):
            if item.get('@type') != 'Dataset':
                continue
            expected_suffix = f"materials id {row['mp_id'].removeprefix('mp-')}"
            if not item.get('name', '').endswith(expected_suffix):
                raise ValueError('Record page does not identify the indexed Materials Project phase')
            ids = [entry.get('@id') for entry in item.get('license', [])]
            return {'status': 'record_cc_by_4_0' if ids == ['cc-by-4.0'] else 'record_license_review',
                    'rights_ids': ids, 'record_title': item['name'], 'page_sha256': digest(page)}
    raise ValueError('No MDR Dataset JSON-LD found')
