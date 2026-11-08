"""Validate a pinned PhononDB replay and emit the hosted GeoForm wire format."""
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path
import re

import numpy as np
import spglib
from phonopy.structure.atoms import PhonopyAtoms

from geoform_phonons import phonondb as db
from geoform_phonons.phonondb_cli import MODULES, PACKAGES
from geoform_phonons.phonondb_verify import (
    FREQUENCY_ATOL, GAMMA_ACOUSTIC_ATOL_THZ, MATRIX_ATOL, MAX_PRIMITIVE_ATOMS,
    Q_CHECKS, VERTICES, indexed_formula, read_package, reduced,
)

ROOT = Path(__file__).resolve().parents[1]
REVISION = '8cabbe0aee3dde48a19f1235864b0d3ab5cb9330'
INDEX_SHA = '15e71b6511f5b3b071b7e28cda9b5b987575bd47956bb96c5ccd8fb58a6f475f'
THERMAL_LIMITS = {
    'free_energy_ev_cell': 2e-6, 'internal_energy_ev_cell': 2e-6,
    'entropy_ev_k_cell': 5e-9, 'cv_ev_k_cell': 5e-9,
}
TEMPERATURES = [0, 50, 100, 200, 300, 500, 800, 1200]
ROCK_SALT = {
    'gf06g7088': ('Halite', {'Cl', 'Na'}, 'https://aflow.dev.materials.duke.edu/p/AB_cF8_225_a_b-001/'),
    'bk128j14x': ('Periclase', {'Mg', 'O'}, 'https://pubs.usgs.gov/bul/1248/report.pdf'),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_bytes())


def encode(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False, separators=(',', ':')) + '\n').encode()


def finite_array(values, shape=None):
    array = np.asarray(values)
    require(array.dtype.kind in 'fiu' and (shape is None or array.shape == shape)
            and np.isfinite(array).all(), 'Invalid numerical array')
    return array.astype(float)


def verify_structure(document, crystal, structure, row):
    n = structure['primitive_atoms']
    require(type(n) is int and 1 <= n <= MAX_PRIMITIVE_ATOMS
            and crystal['n_atoms_primitive'] == n, 'Invalid primitive atom count')
    cell = document['primitive_cell']
    sites = cell['points']
    symbols = [site['symbol'] for site in sites]
    composition = indexed_formula(row['formula'])
    require(len(sites) == n and symbols == crystal['symbols']
            and reduced(Counter(symbols)) == composition == structure['reduced_composition'],
            'Source/payload site identity or composition mismatch')
    lattice = finite_array(crystal['lattice_angstrom'], (3, 3))
    positions = finite_array(crystal['fractional_positions'], (n, 3))
    offsets = positions - finite_array([site['coordinates'] for site in sites], (n, 3))
    masses = finite_array(crystal['masses_amu'], (n,))
    volume = abs(float(np.linalg.det(lattice)))
    require(np.allclose(lattice, cell['lattice'], atol=2e-7, rtol=0)
            and np.allclose(offsets, np.rint(offsets), atol=2e-6, rtol=0)
            and np.all(masses > 0)
            and np.allclose(masses, [site['mass'] for site in sites], atol=2e-6, rtol=0)
            and volume > 1e-8
            and math.isclose(volume, crystal['volume_angstrom3'], rel_tol=1e-10)
            and crystal['volume_angstrom3'] == structure['primitive_volume_angstrom3']
            and structure['supercell_atoms'] == len(document['supercell']['points'])
            and 0 <= structure['asr_max_abs_ev_a2'] <= 1e-6
            and db.SHA256.fullmatch(structure['structure_fingerprint']), 'Source/payload primitive geometry mismatch')
    tolerance = float(document['phonopy'].get('symmetry_tolerance', 1e-5))
    require(0 < tolerance <= 1e-3 and tolerance == structure['source_symmetry_tolerance_angstrom'],
            'Unexpected symmetry tolerance')
    # Use the original unit cell, as the verifier does, to avoid changing the
    # tolerance behavior for rounded primitive cells close to a symmetry boundary.
    unit = document['unit_cell']
    cell = PhonopyAtoms(cell=unit['lattice'], symbols=[site['symbol'] for site in unit['points']],
                        scaled_positions=[site['coordinates'] for site in unit['points']])
    symmetry = spglib.get_symmetry_dataset((cell.cell, cell.scaled_positions, cell.numbers),
                                          symprec=tolerance)
    require(symmetry is not None and symmetry.international == structure['space_group_symbol']
            and symmetry.number == row['space_group_number'] == structure['space_group_number']
            == document['space_group']['number'], 'Source symmetry or receipt mismatch')
    return symmetry


def mineral_mapping(row, structure, symmetry):
    reference = ROCK_SALT.get(row['record_id'])
    if reference is None:
        return None
    name, elements, source = reference
    require(set(structure['reduced_composition']) == elements
            and list(structure['reduced_composition'].values()) == [1, 1]
            and structure['primitive_atoms'] == 2 and symmetry.number == 225,
            'Curated mineral candidate changed structure')
    # Standardize to the conventional cell, then check the full ordered motif.
    # An origin shift may exchange the labels of the two fcc sublattices.
    labels = sorted(set(symmetry.std_types))
    fcc = np.array([(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)])
    sites = symmetry.std_positions
    types = symmetry.std_types
    expected = np.vstack((fcc, (fcc + .5) % 1))
    require(len(labels) == 2 and len(sites) == 8 and np.allclose(
        np.linalg.norm(symmetry.std_lattice, axis=1), np.linalg.norm(symmetry.std_lattice[0]),
        rtol=1e-6), 'Invalid rock-salt conventional cell')
    origin = sites[types == labels[0]][0]
    observed = (np.round((sites - origin) % 1, 5) % 1).tolist()
    require(all(sorted(site for site, label in zip(observed, types) if label == element)
                == sorted(expected[i*4:(i+1)*4].tolist()) for i, element in enumerate(labels)),
            'Rock-salt ordered 4a/4b motif mismatch')
    return {'name': name, 'reference': source, 'basis': 'Ordered rock-salt cF8, Fm-3m, 4a/4b motif'}


def verify_outputs(result, comparison, n, mesh):
    physics = result['physics']
    dispersion = result['dispersion']
    lattice = finite_array(result['crystal']['lattice_angstrom'], (3, 3))
    path = np.concatenate([np.linspace(a, b, 12, endpoint=i == len(VERTICES) - 2)
                           for i, (a, b) in enumerate(zip(VERTICES[:-1], VERTICES[1:]))])
    q = finite_array(dispersion['q_fractional'], path.shape)
    x = finite_array(dispersion['distance_inverse_angstrom'], (len(q),))
    bands = finite_array(dispersion['frequency_thz_signed'], (len(q), 3*n))
    steps = np.linalg.norm(np.diff(q, axis=0) @ (2*np.pi*np.linalg.inv(lattice).T), axis=1)
    require(np.array_equal(q, path) and x[0] == 0 and np.all(np.diff(x) > 0)
            and np.allclose(x[1:], np.cumsum(steps), rtol=1e-10, atol=1e-12)
            and np.all(np.diff(bands, axis=1) >= 0), 'Invalid signed dispersion')
    mesh_points = len(Q_CHECKS) + len(q) + mesh**3
    stability = comparison['stability']
    minimum = comparison['signed_minimum_thz']
    reference_minimum = comparison['reference_signed_minimum_thz']
    imaginary = comparison['imaginary_sampled_modes']
    ignored = comparison['gamma_acoustic_roundoff_modes_excluded']
    require(physics['q_mesh'] == [mesh]*3 and physics['sampled_stability'] == comparison
            and physics['nonanalytic_correction'] == 'omitted on both sides'
            and comparison['matrix_q_points'] == len(Q_CHECKS)
            and comparison['frequency_q_points'] == mesh_points
            and 0 <= comparison['max_abs_matrix_error_ev_a2_amu'] <= MATRIX_ATOL
            and 0 <= comparison['max_abs_frequency_error_thz'] <= FREQUENCY_ATOL
            and math.isfinite(minimum) and minimum <= bands.min()
            and abs(minimum - reference_minimum) <= comparison['max_abs_frequency_error_thz'] + 1e-12
            and type(imaginary) is int and type(ignored) is int
            and 0 <= imaginary <= mesh_points*3*n and 0 <= ignored <= 9
            and stability in ('sampled_positive', 'sampled_imaginary_modes', 'sampled_zero_mesh_modes')
            and (stability == 'sampled_imaginary_modes') == (imaginary > 0),
            'Invalid independent comparison or stability receipt')
    finite_array(comparison['minimum_q_fractional'], (3,))
    if not imaginary:
        gamma = np.all(q == np.rint(q), axis=1)
        require(np.all(bands[~gamma] >= 0) and minimum >= -GAMMA_ACOUSTIC_ATOL_THZ
                and reference_minimum >= -GAMMA_ACOUSTIC_ATOL_THZ
                and (min(minimum, reference_minimum) >= 0 or ignored > 0),
                'Unsupported classification of negative modes')
    dos, thermal = result['density_of_states'], result['thermodynamics']
    if stability != 'sampled_positive':
        withheld = {'status': f'withheld: {stability}'}
        require(dos == thermal == withheld and 'thermodynamic_curve' not in result,
                'Unexpected DOS or thermal output for a non-positive classification')
        return False
    axis = finite_array(dos['frequency_thz'])
    require(axis.ndim == 1 and 400 <= len(axis) <= 8000 and np.all(np.diff(axis) > 0)
            and dos['gaussian_sigma_thz'] == .1, 'Invalid DOS axis or smearing')
    density = finite_array(dos['states_per_thz_cell'], axis.shape)
    integral = float(np.trapezoid(density, axis))
    require(np.all(density >= 0) and abs(integral - 3*n) <= 3*n*1e-5
            and math.isclose(integral, comparison['dos_integral_states_per_cell'], rel_tol=1e-10)
            and 0 <= comparison['max_abs_dos_error_states_per_thz_cell'] <= 1e-4*max(1, density.max()),
            'Invalid DOS normalization or comparison')
    curve = result['thermodynamic_curve']
    errors = comparison['max_abs_thermal_errors']
    require([p['temperature_k'] for p in curve] == TEMPERATURES and thermal == curve[4]
            and errors.keys() == THERMAL_LIMITS.keys()
            and all(0 <= errors[k] <= n*limit for k, limit in THERMAL_LIMITS.items()),
            'Invalid thermal grid or comparison')
    for p in curve:
        require(p.keys() == {*THERMAL_LIMITS, 'temperature_k', 'zero_point_energy_ev_cell'},
                'Unexpected thermal property')
        values = finite_array([p[k] for k in (*THERMAL_LIMITS, 'zero_point_energy_ev_cell')])
        f, u, s, cv, zpe = values
        require(min(s, cv, zpe) >= 0 and zpe == curve[0]['zero_point_energy_ev_cell']
                and math.isclose(f, u-p['temperature_k']*s, rel_tol=1e-10, abs_tol=1e-12),
                'Inconsistent constant-volume thermal values')
    return True


def source_entry(row):
    return dict(zip(('id', 'mp', 'formula', 'group', 'number'),
                    (row[k] for k in ('record_id', 'mp_id', 'formula', 'space_group_symbol', 'space_group_number'))))


def record_rights(rights, row, allow_collection):
    status = rights['status']
    if status == 'record_cc_by_4_0':
        require(rights['rights_ids'] == ['cc-by-4.0'] and db.SHA256.fullmatch(rights['page_sha256'])
                and rights['record_title'].endswith('materials id ' + row['mp_id'][3:]),
                'Malformed record-level rights receipt')
        return status
    require(status in ('collection_license_only', 'record_license_review', 'record_license_unconfirmed'),
            'Unknown record license status')
    unavailable = (status == 'record_license_unconfirmed' and
                   rights.get('error', '').startswith(('TimeoutError:', 'HTTPError: HTTP Error 5')))
    if allow_collection and (status == 'collection_license_only' or unavailable):
        return 'collection_license_only'
    return None


def verify_record(row, receipt, entry, pin, replay, cache):
    record = row['record_id']
    require(receipt['archive'] == pin, 'Receipt archive differs from the release pin')
    archive, _ = db.locked_archive(row, cache, expected_pin=pin)
    document, packed_sha, yaml_sha = read_package(archive)
    require((packed_sha, yaml_sha) == (receipt['yaml_xz_sha256'], receipt['yaml_sha256']),
            'Archived YAML hashes differ from receipt')
    if receipt['status'] == 'verification_failed':
        count = len(document['primitive_cell']['points'])
        require(entry is None and count > MAX_PRIMITIVE_ATOMS and receipt['stage'] == 'structure_and_ifc'
                and receipt['error'] == f'PackageRejected: primitive_cell has {count} sites; supported range is 1..64',
                'Unexpected verification failure. Replay and inspect this record.')
        return None
    require(receipt['status'] == 'verified_calculation', 'Unexpected verification status')
    relative = f'results/{record}.json'
    require(receipt['result_file'] == relative, 'Invalid result path')
    blob = (replay / relative).read_bytes()
    require(db.digest(blob) == receipt['result_sha256'], 'Result digest differs from receipt')
    result = json.loads(blob)
    evidence, crystal = result['evidence'], result['crystal']
    comparison, structure = receipt['comparison'], receipt['structure']
    expected_entry = {**row, 'phase_id': f'nims-mdr:{record}', 'mineral_name': None,
                      'mineral_identity_status': 'unassigned', 'evidence_status': 'calculated',
                      **{k: receipt[k] for k in ('comparison', 'result_file', 'result_sha256',
                                                  'dos_thermal_available', 'record_rights')}}
    require(entry is None or entry == expected_entry, 'Catalog entry differs from receipt')
    require(result['schema_version'] == 'geoform.phonons.result/1'
            and evidence['declared_status'] == 'calculated' and evidence['experimental_validation'] is False
            and evidence['source'] == evidence['source_url'] == row['record_url']
            and evidence['archive'] == pin and evidence['yaml_sha256'] == yaml_sha
            and evidence['attribution'] == receipt['settings']['attribution']
            and evidence['collection_doi'] == db.COLLECTION and evidence['collection_license_url'] == db.LICENSE
            and evidence['ifc_origin'] == structure['ifc_origin']
            and result['physics']['acoustic_sum_rule_max_abs_ev_a2'] == structure['asr_max_abs_ev_a2']
            and crystal['phase_id'] == expected_entry['phase_id'] and crystal['formula'] == row['formula']
            and crystal['mineral_name'] is None and crystal['mineral_identity_status'] == 'unassigned',
            'Invalid result provenance or phase identity')
    symmetry = verify_structure(document, crystal, structure, row)
    dos = verify_outputs(result, comparison, structure['primitive_atoms'], receipt['settings']['mesh'])
    require(receipt['dos_thermal_available'] is dos, 'Receipt output flag differs from payload')
    return blob, dos, mineral_mapping(row, structure, symmetry)


def export_phases(replay, cache, index_path, collection_page, *, allow_collection=False, allow_partial=False):
    manifest = db.parse_index(index_path.read_bytes(), REVISION)
    require(manifest['index_sha256'] == INDEX_SHA, 'Source index pin mismatch')
    collection = collection_page.read_bytes()
    require(len(collection) < 2 * 1024**2 and all(marker in collection for marker in
            (db.COLLECTION.encode(), b'Atsushi Togo', f'href="{db.LICENSE}"'.encode())),
            'Collection landing page is missing DOI, author or CC BY 4.0 evidence')
    catalog, receipts = (load(replay / f'{name}.json') for name in ('catalog', 'receipts'))
    pins = {}
    for batch in ('first_batch', 'breadth'):
        root = ROOT / 'benchmarks'
        batch_pins = load(root / f'phonondb_{batch}_pins.json')
        ids = load(root / f'phonondb_{batch}_ids.json')
        require(set(ids) == batch_pins.keys() and not (pins.keys() & batch_pins.keys()), 'Invalid release pins')
        pins.update(batch_pins)
    settings = receipts['settings']
    require(settings['pipeline_sha256'] == {
        module: db.digest((ROOT / 'geoform_phonons' / module).read_bytes()) for module in MODULES
    } and settings['software'].keys() == set(PACKAGES)
            and all(re.fullmatch(r'[0-9][0-9A-Za-z.+-]*', v) for v in settings['software'].values())
            and settings['mesh'] == 6 and settings['license_url'] == db.LICENSE
            and settings['source_index_sha256'] == catalog['source_index_sha256'] == INDEX_SHA
            and settings['attribution'] == catalog['attribution'] == manifest['attribution']
            and settings['collection_doi'] == catalog['collection_doi'] == db.COLLECTION
            and catalog['collection_license_url'] == db.LICENSE
            and catalog['index_url'] == receipts['index_url'] == manifest['index_url']
            and catalog['schema'] == 'geoform.phonondb.catalog/1'
            and receipts['schema'] == 'geoform.phonondb.receipts/1'
            and receipts['index_rows'] == len(manifest['rows']), 'Catalog or receipt provenance mismatch')
    rows = {row['record_id']: row for row in manifest['rows']}
    by_id = {r['source']['record_id']: r for r in receipts['records']}
    entries = {r['record_id']: r for r in catalog['entries']}
    require(len(by_id) == len(receipts['records']) == receipts['attempted'] and by_id
            and by_id.keys() <= pins.keys() and len(entries) == len(catalog['entries'])
            and entries.keys() <= by_id.keys(), 'Invalid catalog or receipt inventory')
    missing = sorted(pins.keys()-by_id.keys())
    require(allow_partial or not missing, f'Missing {len(missing)} pinned attempts; replay both batches')
    payloads = defaultdict(list)
    for r in receipts['records']:
        if r['status'] == 'verified_calculation':
            payloads[r['yaml_sha256']].append(r['source']['record_id'])
    aliases = {alias: group[0] for group in payloads.values() for alias in group[1:]}
    verified = {record for group in payloads.values() for record in group}
    require(aliases == catalog['byte_identical_payload_aliases']
            and entries.keys() == verified - aliases.keys()
            and len(verified) == receipts['verified_calculations']
            and len(entries) == receipts['unique_verified_payloads'], 'Invalid alias or verified inventory')
    exported, links, assets, withheld = [], {}, {}, {}
    stability = Counter()
    prefix = 'phonons/nims/'
    for record in sorted(by_id):
        row, receipt = rows[record], by_id[record]
        try:
            require(receipt['source'] == row and receipt['settings'] == settings, 'Source row or settings mismatch')
            product = verify_record(row, receipt, entries.get(record), pins[record], replay, cache)
            if product is None:
                withheld[record] = 'primitive_cell_budget'
                continue
            rights = record_rights(receipt['record_rights'], row, allow_collection)
            if record in aliases or rights is None:
                withheld[record] = 'identical_source_payload' if record in aliases else receipt['record_rights']['status']
                continue
            blob, dos, mapping = product
            packed = gzip.compress(blob, compresslevel=9, mtime=0)
            assets[f'{prefix}{record}.bin'] = packed
            check = receipt['comparison']['stability']
            require(check != 'sampled_zero_mesh_modes', 'Hosted schema 1 does not display zero-mesh-mode cases')
            stability[check] += 1
            exported.append({**source_entry(row), 'dos': dos, 'stability': check, 'rights': rights,
                             'sha256': db.digest(blob), 'bytes': len(blob), 'packedBytes': len(packed),
                             **({'mineralStructure': mapping} if mapping else {})})
            if mapping:
                links[mapping['name']] = f'#phonons?mineral=nims-{record}'
        except (ValueError, KeyError, TypeError, OSError) as exc:
            raise ValueError(f'{record}: {exc}') from exc
    require(exported, 'No redistributable verified phases; inspect rights receipts')
    require(allow_partial or len(exported) == len(entries),
            f'{len(entries) - len(exported)} verified phases lack acceptable rights evidence')
    compositions = {tuple(sorted(by_id[r['id']]['structure']['reduced_composition'].items()))
                    for r in exported}
    groups = {r['number'] for r in exported}
    if not allow_partial:
        require(len(verified) == 542 and len(compositions) == 372 and len(groups) == 105
                and stability == {
            'sampled_imaginary_modes': 277, 'sampled_positive': 265},
            'Full replay disagrees with the frozen release coverage')
    exported.sort(key=lambda row: (int(row['mp'][3:]), row['id']))
    source = {'schema': 1, 'sourceIndexSha256': INDEX_SHA}
    assets[prefix+'index.bin'] = gzip.compress(encode({
        **source, 'attempted': len(by_id), 'indexRecords': len(rows),
        'coverage': stability, 'entries': exported,
    }), mtime=0)
    assets[prefix+'source-index.bin'] = gzip.compress(encode({
        **source, 'entries': [source_entry(row) for row in rows.values()],
    }), mtime=0)
    assets['phonon-nims-links.mjs'] = b'export const nimsMineralLinks = new Map(Object.entries(' + encode(links).strip() + b'));\n'
    assets['PHONONDB_EXPORT_REPORT.json'] = encode({
        'attempted': len(by_id), 'verified': len(verified), 'exported': len(exported),
        'partial_review': allow_partial, 'collection_license_fallback_enabled': allow_collection,
        'coverage': stability, 'rights': Counter(row['rights'] for row in exported),
        'reduced_compositions': len(compositions), 'space_group_numbers': len(groups),
        'mineral_links': links, 'withheld': withheld, 'missing_pinned_attempts': missing,
        'settings': settings, 'collection_page_sha256': db.digest(collection),
        'inputs': {name: db.digest((replay / name).read_bytes()) for name in ('catalog.json', 'receipts.json')},
        'assets': {path: {'bytes': len(data), 'sha256': db.digest(data)} for path, data in assets.items()},
    })
    return assets
