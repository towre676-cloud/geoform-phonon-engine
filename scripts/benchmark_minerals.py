"""Reproduce source-pinned *calculated* phonons against independent Phonopy.

Requires: `pip install -e . phonopy` and pinned reference source downloads.
This script fails loudly on any matrix/frequency discrepancy; no substitution
of synthetic force constants, mocked Phonopy calls, or inferred IFC labels.

The calculation omits polar non-analytic LO-TO correction on BOTH sides.
The exported result is not an experimental spectrum or validation of DFT.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geoform_phonons.core import shifted_q_mesh, THZ_FACTOR
from geoform_phonons.io import export_result, save_json
from geoform_phonons.phonopy_native import from_phonopy
from scripts.audit_vasp_source import check_source
from scripts.fetch_reference_sources import MANIFEST, matches_pin, material_files, validate_manifest

# Fixed rational, off-symmetry points in primitive reciprocal coordinates.
# Include Gamma, zone edges and general complex dynamical matrices.
Q_POINTS = (
    (0.,0.,0.), (0.5,0.,0.), (0.5,0.5,0.), (0.5,0.5,0.5),
    (0.125,0.25,0.375), (0.17,0.31,0.09), (0.41,0.02,0.28),
    (0.21,0.14,0.07), (0.36,0.23,0.16), (0.49,0.48,0.11),
    (0.13,0.44,0.27), (0.33,0.37,0.19), (0.07,0.26,0.45),
    (-0.22,0.14,0.31), (1.17,0.31,0.09), (0.25,0.,0.5),
)
MAX_DMAT_ATOL = 5e-6  # unit (eV/Ang^2)/amu (after shortest image reconstruction)
MAX_FREQ_ATOL_THZ = 1e-3  # includes independent CODATA conversion rounding and near Gamma


def check_input_hashes(slug, data_dir, entry):
    digests = {}
    for filename, info in material_files(entry).items():
        data = (data_dir / slug / filename).read_bytes()
        if not matches_pin(data, info):
            raise ValueError(f'{slug}/{filename}: pinned source digest mismatch')
        digests[filename] = info['sha256']
    return digests


def compare_qpoints(slug, model, reference):
    checks = []
    for q, external_mat, external_freq in zip(
            Q_POINTS, reference.dynamical_matrices, reference.frequencies, strict=True):
        our_mat = model.dynamical_matrix(q)
        computed_freq = model.modes(q)[0]
        matrix_error = float(np.max(np.abs(our_mat - external_mat)))
        freq_error = float(np.max(np.abs(computed_freq - np.sort(external_freq))))
        if not np.isfinite(matrix_error) or matrix_error > MAX_DMAT_ATOL:
            raise AssertionError(f'{slug} at {q}: dynamical matrix mismatch {matrix_error:.8g}')
        if not np.isfinite(freq_error) or freq_error > MAX_FREQ_ATOL_THZ:
            raise AssertionError(f'{slug} at {q}: frequency mismatch {freq_error:.8g} THz')
        checks.append({
            'q_fractional': list(q), 'max_abs_dynamical_matrix_error': matrix_error,
            'max_abs_frequency_error_thz': freq_error,
            'max_abs_imaginary_matrix_entry': float(np.max(np.abs(external_mat.imag))),
            'min_frequency_thz': float(computed_freq.min()),
        })
    return checks


def signed_mode_summary(qpoints, frequencies):
    frequencies = np.asarray(frequencies)
    index = int(np.argmin(frequencies) // frequencies.shape[1])
    return {'q_points': len(qpoints), 'imaginary_modes': int(np.count_nonzero(frequencies < 0)),
            'min_frequency_thz': float(frequencies.min()), 'min_q_fractional': list(qpoints[index])}


def compare_one(slug, entry, data_dir, out_dir, *, mesh=8, path_points=24):
    import phonopy
    from phonopy.harmonic.force_constants import compact_fc_to_full_fc

    path = data_dir / slug
    digests = check_input_hashes(slug, data_dir, entry)
    phonon = phonopy.load(
        phonopy_yaml=str(path / 'phonopy_disp.yaml'),
        force_sets_filename=str(path / 'FORCE_SETS'),
        primitive_matrix='auto',  # stored YAML matrix takes precedence where present
        fc_calculator=None,  # Type I finite differences with Phonopy's symfc projector
        symmetrize_fc=True,
        is_compact_fc=True,
        is_nac=False,
        lang='C',
    )
    if phonon.force_constants is None or phonon.nac_params is not None:
        raise RuntimeError('expected bare harmonic force constants from the pinned FORCE_SETS')
    source_audit = check_source(path, entry, phonon) if 'structure' in entry else None
    model = from_phonopy(phonon, tolerance=4e-6)
    asr = model.acoustic_sum_rule_error()
    if asr > 1e-6:
        raise AssertionError(f'{slug}: ASR residual too large ({asr:.9g} eV/A^2)')

    reference = phonon.run_qpoints(Q_POINTS, with_dynamical_matrices=True)
    checks = compare_qpoints(slug, model, reference)
    max_matrix = max(check['max_abs_dynamical_matrix_error'] for check in checks)
    max_frequency = max(check['max_abs_frequency_error_thz'] for check in checks)
    compact_shape = list(phonon.force_constants.shape)
    full_shape = None
    full_checks = None
    if source_audit is not None:
        phonon.force_constants = compact_fc_to_full_fc(phonon.primitive, phonon.force_constants)
        full_shape = list(phonon.force_constants.shape)
        full_model = from_phonopy(phonon, tolerance=4e-6)
        full_checks = compare_qpoints(slug + ' full IFC', full_model, reference)

    upstream_commit = json.loads(MANIFEST.read_text())['upstream_commit']
    source_url = ('https://github.com/phonopy/phonopy/tree/'
                  + upstream_commit + '/' + entry['source_directory'])
    source = f'Phonopy {entry["source_directory"]} at {upstream_commit}'
    provenance = {'evidence_status': 'calculated', 'source': source}
    vertices = [[0,0,0], [0.5,0,0], [0.5,0.5,0], [0,0,0], [0.5,0.5,0.5]]
    mesh_shape = (mesh,) * 3
    result = export_result(model, provenance, vertices,
                           points_per_segment=path_points, mesh=mesh_shape, temperature_k=300)
    band = result['dispersion']
    mesh_q = shifted_q_mesh(mesh_shape)
    mesh_freqs = [model.modes(q)[0] for q in mesh_q]
    stability = {
        'dispersion': signed_mode_summary(band['q_fractional'], band['frequency_thz_signed']),
        'shifted_mesh': signed_mode_summary(mesh_q, mesh_freqs),
        'comparison_points': signed_mode_summary(Q_POINTS, [model.modes(q)[0] for q in Q_POINTS]),
    }
    # Check the whole exported path independently, including soft modes that
    # may be missed by both the uniform mesh and the fixed matrix comparisons.
    reference_bands = phonon.run_qpoints(band['q_fractional']).frequencies
    path_error = float(np.max(np.abs(np.asarray(band['frequency_thz_signed'])
                                     - np.sort(reference_bands, axis=1))))
    if not np.isfinite(path_error) or path_error > MAX_FREQ_ATOL_THZ:
        raise AssertionError(f'{slug}: dispersion frequency mismatch {path_error:.8g} THz')
    dispersion_reference = {
        'q_points': len(reference_bands),
        'max_abs_frequency_error_thz': path_error,
        'min_frequency_thz': float(reference_bands.min()),
    }
    if any(probe['imaginary_modes'] for probe in stability.values()):
        result['density_of_states'] = {'status': 'withheld: imaginary modes on sampled reciprocal points'}
        result['thermodynamics'] = {'status': 'withheld: imaginary modes on sampled reciprocal points'}
        result.pop('thermodynamic_curve', None)
    result['physics']['sampled_stability'] = stability
    result['physics']['nonanalytic_correction'] = 'omitted; polar LO-TO splitting not reproduced'
    result['physics']['force_constants_preprocessing'] = 'Phonopy Type I finite differences with symfc projection; symmetrize_fc=True'
    band['path_description'] = 'Deterministic primitive reciprocal-coordinate path, without symmetry labels'
    result['crystal'].update(mineral_name=entry['name'], display_name=entry['display_name'],
                            formula=entry['formula'], phase=entry.get('phase', 'Rock-salt structure'),
                            lattice_angstrom=model.crystal.lattice_angstrom.tolist(),
                            fractional_positions=model.crystal.fractional_positions.tolist(),
                            masses_amu=model.crystal.masses_amu.tolist())
    result['evidence'].update(
        source_url=source_url, source_note=entry.get('source_note'),
        verification='Compared with Phonopy using the same calculated forces. Experimental validation has not been performed.',
        benchmark={'reference': 'Phonopy ' + phonopy.__version__, 'source_sha256': digests,
                   'comparison_points': len(Q_POINTS), 'max_abs_frequency_error_thz': max_frequency,
                   'max_abs_dynamical_matrix_error': max_matrix,
                   'independent_comparison': 'passed; computational, not experimental'})
    filename = f'{slug}_phonons.json'
    out_dir.mkdir(parents=True, exist_ok=True)
    save_json(result, out_dir / filename)
    return {
        'material': slug, 'name': entry['name'], 'display_name': entry['display_name'],
        'formula': entry['formula'], 'phase': result['crystal']['phase'],
        'evidence_status': 'calculated', 'independent_reference': 'Phonopy ' + phonopy.__version__,
        'from_supplied_forces': True, 'source_url': source_url, 'source_note': entry.get('source_note'),
        'source_sha256': digests, 'source_audit': source_audit,
        'supercell_ifc_shape': compact_shape, 'full_ifc_shape': full_shape,
        'primitive_atoms': model.crystal.n_atoms, 'periodic_ifc_blocks': len(model.blocks),
        'asr_error_ev_a2': asr, 'max_abs_dynamical_matrix_error': max_matrix,
        'max_abs_frequency_error_thz': max_frequency, 'qpoint_comparisons': checks,
        'full_ifc_qpoint_comparisons': full_checks, 'sampled_stability': stability,
        'dispersion_reference': dispersion_reference,
        'dos_thermal_available': 'temperature_k' in result['thermodynamics'],
        'reference_file': filename,
        'result_sha256': hashlib.sha256((out_dir / filename).read_bytes()).hexdigest(),
        'lo_to_limitation': result['physics']['nonanalytic_correction'],
        'frequency_conversion_thz': {'geoform': THZ_FACTOR, 'phonopy': phonon.unit_conversion_factor},
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,default=ROOT/'benchmarks'/'upstream')
    p.add_argument('--output-dir',type=Path,default=ROOT/'benchmarks'/'generated')
    p.add_argument('--mesh',type=int,default=8)
    p.add_argument('--points',type=int,default=24)
    args=p.parse_args(argv)
    if not 2<=args.mesh<=40 or not 2<=args.points<=200:
        raise ValueError('mesh 2..40 and path points 2..200 required')
    manifest=json.loads(MANIFEST.read_text())
    validate_manifest(manifest)
    if not matches_pin((args.data_dir / 'LICENSE').read_bytes(), manifest['license_file']):
        raise ValueError('pinned license digest mismatch')
    receipts=[]
    for slug, entry in sorted(manifest['materials'].items()):
        receipt=compare_one(slug,entry,args.data_dir,args.output_dir,
                            mesh=args.mesh,path_points=args.points)
        receipts.append(receipt)
        print(f'{slug}: PASS {len(receipt["qpoint_comparisons"])} reciprocal q points; '
              f'max frequency discrepancy {receipt["max_abs_frequency_error_thz"]:.8g} THz; '
              f'D(q) error {receipt["max_abs_dynamical_matrix_error"]:.8g}; '
              f'ASR {receipt["asr_error_ev_a2"]:.4g}')
    catalog={'schema_version':'geoform.phonons.catalog/1',
             'catalog_kind':'validated_calculated_benchmarks',
             'upstream_repository':manifest['upstream_repo'],
             'upstream_commit':manifest['upstream_commit'],
             'benchmark_scope':'non-NAC harmonic finite-displacement DFT; computational comparison only',
             'entries':[{'slug':r['material'],'mineral':r['name'],'formula':r['formula'],
                         'display_name':r['display_name'],'phase':r['phase'],
                         'dos_thermal_available':r['dos_thermal_available'],
                         'evidence_status':r['evidence_status'],
                         'result_file':r['reference_file'],
                         'comparison_points':len(r['qpoint_comparisons']),
                         'max_frequency_error_thz':r['max_abs_frequency_error_thz']} for r in receipts]}
    save_json(catalog,args.output_dir/'catalog.json')
    save_json({'schema_version':'geoform.phonons.benchmark_receipts/1',
               'upstream_commit':manifest['upstream_commit'],
               'software_versions':{name: importlib.metadata.version(name) for name in ('numpy', 'phonopy', 'spglib', 'symfc', 'scipy')},
               'results':receipts},
              args.output_dir/'benchmark_receipts.json')
    print('PASS: all Phonopy real-input computational cross-checks')

if __name__=='__main__':main()
