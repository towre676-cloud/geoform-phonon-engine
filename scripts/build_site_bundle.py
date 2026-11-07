"""Package only *verified computational* mineral results for a static GeoForm viewer.

Refuses incomplete result schemas and missing independent comparison receipts;
never silently promotes an illustrative model to a mineral calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geoform_phonons.io import save_json
from scripts.benchmark_minerals import MAX_DMAT_ATOL, MAX_FREQ_ATOL_THZ, Q_POINTS
from scripts.fetch_reference_sources import MANIFEST, material_files

HANDOFF = ROOT / 'site' / 'geoform-phonons'


def verify_results(results):
    """Return the catalog entries after validating evidence and cross-checks."""
    results = Path(results)
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    catalog = json.loads((results / 'catalog.json').read_text(encoding='utf-8'))
    receipts = json.loads((results / 'benchmark_receipts.json').read_text(encoding='utf-8'))
    pin = manifest['upstream_commit']
    if (catalog.get('schema_version') != 'geoform.phonons.catalog/1'
            or receipts.get('schema_version') != 'geoform.phonons.benchmark_receipts/1'
            or catalog.get('upstream_commit') != pin or receipts.get('upstream_commit') != pin):
        raise ValueError('invalid catalog/receipt schema or upstream pin')
    materials = manifest['materials']
    indexed = {item['material']: item for item in receipts['results']}
    if (len(receipts['results']) != len(materials) or indexed.keys() != materials.keys()
            or len(catalog['entries']) != len(materials)
            or {item['slug'] for item in catalog['entries']} != materials.keys()):
        raise ValueError('catalog and receipt inventory must match the pinned benchmark manifest')

    def valid_qpoints(checks):
        return (isinstance(checks, list) and len(checks) == len(Q_POINTS)
                and all(check['q_fractional'] == list(q)
                        and 0 <= check['max_abs_dynamical_matrix_error'] <= MAX_DMAT_ATOL
                        and 0 <= check['max_abs_frequency_error_thz'] <= MAX_FREQ_ATOL_THZ
                        for check, q in zip(checks, Q_POINTS, strict=True)))

    for item in catalog['entries']:
        slug = item['slug']
        material = materials[slug]
        receipt = indexed[slug]
        if (item['evidence_status'] != 'calculated' or receipt['evidence_status'] != 'calculated'
                or item['result_file'] != f'{slug}_phonons.json'
                or receipt['reference_file'] != item['result_file']
                or item['mineral'] != material['name']
                or item['display_name'] != material['display_name']
                or item['formula'] != material['formula']
                or not valid_qpoints(receipt['qpoint_comparisons'])
                or ('structure' in material and (
                    not receipt['source_audit'] or not valid_qpoints(receipt['full_ifc_qpoint_comparisons'])))):
            raise ValueError(f'{slug}: invalid catalog identity or reference checks')
        blob = (results / item['result_file']).read_bytes()
        record = json.loads(blob)
        evidence = record.get('evidence', {})
        check = evidence.get('benchmark', {})
        hashes = {name: info['sha256'] for name, info in material_files(material).items()}
        available = 'temperature_k' in record['thermodynamics']
        dispersion = receipt['dispersion_reference']
        stability = receipt['sampled_stability']
        if (hashlib.sha256(blob).hexdigest() != receipt['result_sha256']
                or record.get('schema_version') != 'geoform.phonons.result/1'
                or evidence.get('declared_status') != 'calculated'
                or check.get('independent_comparison') != 'passed; computational, not experimental'
                or check.get('comparison_points') != len(Q_POINTS)
                or check.get('source_sha256') != hashes or receipt['source_sha256'] != hashes
                or record['crystal'].get('mineral_name') != material['name']
                or record['physics'].get('sampled_stability') != stability
                or dispersion['q_points'] != len(record['dispersion']['q_fractional'])
                or not 0 <= dispersion['max_abs_frequency_error_thz'] <= MAX_FREQ_ATOL_THZ
                or record['physics'].get('nonanalytic_correction') != 'omitted; polar LO-TO splitting not reproduced'
                or available != item['dos_thermal_available'] or available != receipt['dos_thermal_available']
                or available != ('frequency_thz' in record['density_of_states'])
                or available != ('thermodynamic_curve' in record)
                or (available and any(p['imaginary_modes'] for p in stability.values()))):
            raise ValueError(f'{slug}: missing or inconsistent computational evidence')
    return catalog['entries']


def result_names(items):
    return ['catalog.json', 'benchmark_receipts.json', *(item['result_file'] for item in items)]


def stage_handoff(results):
    items = verify_results(results)
    data_dir = HANDOFF / 'data'
    for name in result_names(items):
        source = Path(results) / name
        target = data_dir / name
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
    path = HANDOFF / 'site-handoff.json'
    handoff = json.loads(path.read_text(encoding='utf-8'))
    previous_minerals = handoff['supported_minerals']
    handoff['supported_minerals'] = {
        item['slug']: {**previous_minerals.get(item['slug'], {}), 'formula': item['formula'],
                       'viewer_query': f'?mineral={item["slug"]}'}
        for item in items if item['mineral'] is not None
    }
    handoff['computational_phases_without_mineral_link'] = [
        item['slug'] for item in items if item['mineral'] is None
    ]
    paths = ['index.html', 'geoform-phonon-card.js', 'README.md',
             *('data/' + name for name in result_names(items))]
    handoff['assets'] = {}
    for name in sorted(paths):
        data = (HANDOFF / name).read_bytes()
        handoff['assets'][name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    save_json(handoff, path)


def create_bundle(results, output):
    results=Path(results)
    items=verify_results(results)
    out=Path(output)
    out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        z.write(ROOT/'web'/'phonon_explorer.html',arcname='phonon_explorer.html')
        for name in result_names(items):
            z.write(results / name, arcname=f'data/{name}')
        z.write(MANIFEST, arcname='REFERENCE_SOURCES.json')
        z.write(ROOT/'docs'/'REAL_MINERAL_BENCHMARKS.md',arcname='SCIENCE_AND_PROVENANCE.md')
    return out


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=ROOT/'benchmarks'/'generated')
    p.add_argument('--output',type=Path,default=ROOT/'benchmarks'/'geoform_site_phonons.zip')
    p.add_argument('--stage-handoff', action='store_true',
                   help='refresh checked-in Site handoff assets and their local digests')
    args=p.parse_args(argv)
    if args.stage_handoff:
        stage_handoff(args.results)
    output=create_bundle(args.results,args.output)
    print(f'Verified site bundle: {output}, {output.stat().st_size:,} bytes')

if __name__=='__main__':main()
