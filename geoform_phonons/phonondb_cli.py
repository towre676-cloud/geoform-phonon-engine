"""Build a local NIMS PhononDB calculation catalog without changing the site."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import importlib.metadata
import json
from pathlib import Path
from time import monotonic

from . import phonondb as db


PACKAGES = ('numpy', 'phonopy', 'spglib', 'symfc', 'scipy', 'PyYAML')
MODULES = ('core.py', 'phonopy_native.py', 'phonondb.py', 'phonondb_verify.py', 'phonondb_cli.py')


def select(manifest, ids, start, limit):
    rows = manifest['rows']
    if ids is not None:
        wanted = set(json.loads(ids.read_text()))
        rows = [row for row in rows if row['record_id'] in wanted]
        if len(rows) != len(wanted):
            raise ValueError('Selection contains unknown MDR IDs')
    selected = rows[start:start + limit] if limit is not None else rows[start:]
    if not selected:
        raise ValueError('No selected index rows')
    return selected


def license_one(args):
    row, timeout = args
    try:
        result = db.inspect_record_license(row, timeout=timeout)
    except (OSError, ValueError, KeyError, TypeError) as error:
        result = {'status': 'record_license_unconfirmed', 'error': f'{type(error).__name__}: {error}'}
    return row['record_id'], result


def verify_one(args):
    from .phonondb_verify import load_harmonic, read_package, sample_harmonic

    row, cache, output, settings = args
    record = row['record_id']
    start = monotonic()
    receipt_path = output / 'records' / f'{record}.json'
    product = output / 'results' / f'{record}.json'
    stage = 'archive'
    fields = {}
    try:
        data, pin = db.locked_archive(row, cache)
        fields['archive'] = pin
        if receipt_path.exists() and product.exists():
            previous = json.loads(receipt_path.read_text())
            if (previous['status'] == 'verified_calculation' and previous['settings'] == settings
                    and previous['source'] == row and previous['archive'] == pin
                    and previous['result_sha256'] == db.digest(product.read_bytes())):
                return {**previous, 'reused': True, 'seconds': round(monotonic() - start, 3)}
        stage = 'package'
        document, packed, unpacked = read_package(data)
        fields.update(yaml_xz_sha256=packed, yaml_sha256=unpacked)
        stage = 'structure_and_ifc'
        phonon, model, structure = load_harmonic(document, row)
        fields['structure'] = structure
        stage = 'numerical_comparison'
        comparison, result = sample_harmonic(phonon, model, settings['mesh'])
        crystal = model.crystal
        result.update(
            schema_version='geoform.phonons.result/1',
            evidence={
                'declared_status': 'calculated', 'source': row['record_url'],
                'source_url': row['record_url'], 'attribution': settings['attribution'],
                'collection_doi': settings['collection_doi'], 'collection_license_url': settings['license_url'],
                'archive': pin, 'yaml_sha256': unpacked, 'ifc_origin': structure['ifc_origin'],
                'verification': 'Independent Fourier/eigenvalue comparison with Phonopy on shared calculated inputs',
                'experimental_validation': False, 'changes': 'Derived GeoForm non-NAC harmonic calculation',
            },
            crystal={
                'phase': f"{row['formula']} · {structure['space_group_symbol']} · {row['mp_id']}",
                'phase_id': f'nims-mdr:{record}', 'mineral_identity_status': 'unassigned',
                'mineral_name': None, 'formula': row['formula'],
                'lattice_angstrom': crystal.lattice_angstrom.tolist(),
                'fractional_positions': crystal.fractional_positions.tolist(),
                'masses_amu': crystal.masses_amu.tolist(), 'symbols': list(crystal.symbols),
                'volume_angstrom3': crystal.volume_angstrom3, 'n_atoms_primitive': crystal.n_atoms,
            },
        )
        result['physics']['acoustic_sum_rule_max_abs_ev_a2'] = structure['asr_max_abs_ev_a2']
        stage = 'export'
        db.write_json(product, result)
        fields.update(status='verified_calculation', comparison=comparison,
                      result_file=f'results/{record}.json', result_sha256=db.digest(product.read_bytes()),
                      dos_thermal_available='thermodynamic_curve' in result)
    except Exception as error:
        # A malformed source must not abort a batch. Preserve the exception type
        # and stage so implementation faults remain visible in the receipt.
        product.unlink(missing_ok=True)
        fields.update(status='verification_failed', stage=stage,
                      error=f'{type(error).__name__}: {error}')
    receipt = {'source': row, 'settings': settings, **fields,
               'seconds': round(monotonic() - start, 3), 'reused': False}
    db.write_json(receipt_path, receipt)
    return receipt


def verified_batch(args, manifest, rows):
    settings = {
        'mesh': args.mesh, 'software': {name: importlib.metadata.version(name) for name in PACKAGES},
        'source_index_sha256': manifest['index_sha256'],
        'attribution': manifest['attribution'], 'collection_doi': manifest['collection_doi'],
        'license_url': manifest['collection_license_url'],
        'pipeline_sha256': {name: db.digest(Path(__file__).with_name(name).read_bytes()) for name in MODULES},
    }
    start = monotonic()
    jobs = [(row, args.cache, args.output, settings) for row in rows]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        receipts = list(pool.map(verify_one, jobs))
    license_checks = json.loads(args.rights.read_text()) if args.rights else {}
    catalogs, fingerprints, payloads = [], defaultdict(list), defaultdict(list)
    for receipt in receipts:
        row = receipt['source']
        record = row['record_id']
        receipt['record_rights'] = license_checks.get(record, {'status': 'collection_license_only'})
        if receipt['status'] != 'verified_calculation':
            continue
        structure = receipt['structure']
        fingerprints[structure['structure_fingerprint']].append(record)
        payloads[receipt['yaml_sha256']].append(record)
        catalogs.append({
            **row, 'phase_id': f'nims-mdr:{record}', 'mineral_name': None,
            'mineral_identity_status': 'unassigned', 'evidence_status': 'calculated',
            'comparison': receipt['comparison'], 'result_file': receipt['result_file'],
            'result_sha256': receipt['result_sha256'],
            'dos_thermal_available': receipt['dos_thermal_available'],
            'record_rights': receipt['record_rights'],
        })
    exact = {key: ids for key, ids in payloads.items() if len(ids) > 1}
    candidates = {key: ids for key, ids in fingerprints.items() if len(ids) > 1}
    # Only byte-identical calculation payloads are aliased. Structure similarity
    # alone never aliases or overwrites separately sourced force data.
    aliases = {alias: ids[0] for ids in exact.values() for alias in ids[1:]}
    db.write_json(args.output / 'catalog.json', {
        'schema': 'geoform.phonondb.catalog/1', 'index_url': manifest['index_url'],
        'source_index_sha256': manifest['index_sha256'], 'attribution': manifest['attribution'],
        'collection_doi': manifest['collection_doi'], 'collection_license_url': manifest['collection_license_url'],
        'entries': [entry for entry in catalogs if entry['record_id'] not in aliases],
        'byte_identical_payload_aliases': aliases, 'structure_duplicate_candidates': candidates,
    })
    elapsed = monotonic() - start
    summary = {
        'schema': 'geoform.phonondb.receipts/1', 'index_url': manifest['index_url'], 'settings': settings,
        'index_rows': len(manifest['rows']), 'attempted': len(rows),
        'verified_calculations': len(catalogs), 'unique_verified_payloads': len(catalogs) - len(aliases),
        'reused': sum(row['reused'] for row in receipts), 'seconds': round(elapsed, 3),
        'stability_counts': dict(Counter(row['comparison']['stability'] for row in catalogs)),
        'failure_stages': dict(Counter(row['stage'] for row in receipts if 'stage' in row)),
        'record_rights_counts': dict(Counter(row['record_rights']['status'] for row in catalogs)),
        'with_dos_thermal': sum(row['dos_thermal_available'] for row in catalogs),
        'exact_duplicate_payload_groups': exact, 'structure_duplicate_candidates': candidates,
        'records': receipts,
    }
    db.write_json(args.output / 'receipts.json', summary)
    return {key: value for key, value in summary.items() if key not in ('records', 'settings')}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    index = commands.add_parser('index', help='Parse a pinned upstream index snapshot')
    index.add_argument('markdown', type=Path)
    index.add_argument('revision')
    index.add_argument('output', type=Path)
    index.add_argument('--sha256', help='Require the SHA-256 of the input index snapshot')
    for name in ('fetch', 'rights', 'verify'):
        command = commands.add_parser(name)
        command.add_argument('index', type=Path)
        command.add_argument('output', type=Path)
        command.add_argument('--ids', type=Path, help='JSON array of MDR record IDs')
        command.add_argument('--start', type=int, default=0)
        command.add_argument('--limit', type=int)
        command.add_argument('--workers', type=int, default=4)
        if name != 'verify':
            command.add_argument('--timeout', type=float, default=45)
        if name != 'rights':
            command.add_argument('--cache', type=Path, required=True)
        if name == 'fetch':
            command.add_argument('--pins', type=Path, help='Replay an archive checksum manifest')
        if name == 'verify':
            command.add_argument('--mesh', type=int, default=6)
            command.add_argument('--rights', type=Path)
    args = parser.parse_args(argv)
    if args.command == 'index':
        manifest = db.parse_index(args.markdown.read_bytes(), args.revision)
        if args.sha256 is not None and args.sha256 != manifest['index_sha256']:
            parser.error('Index bytes differ from the expected SHA-256')
        db.write_json(args.output, manifest)
        print(f"{len(manifest['rows'])} index rows")
        return 0
    if not 1 <= args.workers <= 32 or args.start < 0 or (args.limit is not None and args.limit < 1):
        parser.error('Require workers 1..32, start >= 0 and positive limit')
    manifest = json.loads(args.index.read_text())
    rows = select(manifest, args.ids, args.start, args.limit)
    # Revalidate externally supplied manifest URLs before any network request or file access.
    for row in rows:
        record = row['record_id']
        if (not db.RECORD_ID.fullmatch(record)
                or row['archive_url'] != f'{db.MDR}/download_all/{record}.zip'
                or row['record_url'] != f'{db.MDR}/concern/datasets/{record}'):
            parser.error('Unrecognized source row')
    if len({row['record_id'] for row in rows}) != len(rows):
        parser.error('Repeated record IDs')
    if args.command == 'fetch':
        started = monotonic()
        pins = json.loads(args.pins.read_text()) if args.pins else None
        receipts = list(db.fetch_batch(rows, args.cache, args.workers, args.timeout, pins=pins))
        counts = dict(Counter(row['status'] for row in receipts))
        db.write_json(args.output, {'counts': counts, 'seconds': round(monotonic() - started, 3),
                                    'records': receipts})
        print(json.dumps(counts))
        return int(bool(counts.get('fetch_failed')))
    if args.command == 'rights':
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            records = dict(pool.map(license_one, [(row, args.timeout) for row in rows]))
        db.write_json(args.output, records)
        print(json.dumps(dict(Counter(row['status'] for row in records.values()))))
        return int(any(row['status'] != 'record_cc_by_4_0' for row in records.values()))
    if not 2 <= args.mesh <= 12:
        parser.error('Require mesh 2..12')
    summary = verified_batch(args, manifest, rows)
    print(json.dumps(summary))
    return int(bool(summary['failure_stages']))


if __name__ == '__main__':
    raise SystemExit(main())
