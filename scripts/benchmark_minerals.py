"""Reproduce halite/periclase *calculated* phonons against independent Phonopy.

Requires: `pip install -e . phonopy` and pinned reference source downloads.
This script fails loudly on any matrix/frequency discrepancy; no substitution
of synthetic force constants, mocked Phonopy calls, or inferred IFC labels.

The calculation omits polar non-analytic LO-TO correction on BOTH sides.
The exported result is not an experimental spectrum or validation of DFT.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geoform_phonons.core import PhononInputError
from geoform_phonons.io import export_result, save_json
from geoform_phonons.phonopy_native import from_phonopy
from scripts.fetch_reference_sources import git_blob_sha1, MANIFEST

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
    file_digests={}
    for filename, info in entry['files'].items():
        path=data_dir/slug/filename
        data=path.read_bytes()
        sha=git_blob_sha1(data)
        if sha != info['git_blob_sha1']:
            raise ValueError(f'{slug}/{filename}: expected pinned Git blob {info["git_blob_sha1"]}, got {sha}')
        file_digests[filename]=hashlib.sha256(data).hexdigest()
    return file_digests


def compare_one(slug, entry, data_dir, out_dir, *, mesh=4, path_points=12):
    path=data_dir/slug
    digests=check_input_hashes(slug, data_dir, entry)
    from phonopy import load as phonopy_load
    phonon=phonopy_load(
        phonopy_yaml=str(path/'phonopy_disp.yaml'),
        force_sets_filename=str(path/'FORCE_SETS'),
        symmetrize_fc=True,
        is_compact_fc=True,
        is_nac=False,
        lang='C',
    )
    if phonon.force_constants is None:
        raise RuntimeError('Phonopy failed to construct IFCs from the reference FORCE_SETS')
    model=from_phonopy(phonon, tolerance=4e-6)
    asr=model.acoustic_sum_rule_error()
    if asr > 1e-6:
        raise AssertionError(f'{slug}: ASR residual too large ({asr:.9g} eV/A^2)')
    ref_dm=phonon.dynamical_matrix
    if ref_dm is None:
        raise RuntimeError('Phonopy did not expose an independent dynamical-matrix calculator')
    max_matrix=0.0
    max_frequency=0.0
    checks=[]
    for q in Q_POINTS:
        external_freq=np.sort(np.asarray(phonon.get_frequencies(q),float))
        ref_dm.run(q)
        external_mat=np.asarray(ref_dm.dynamical_matrix)
        our_mat=model.dynamical_matrix(q)
        matrix_error=float(np.max(np.abs(our_mat - external_mat)))
        computed_freq=model.modes(q)[0]
        freq_error=float(np.max(np.abs(computed_freq-external_freq)))
        if matrix_error > MAX_DMAT_ATOL:
            raise AssertionError(f'{slug} at {q}: dynamical matrix mismatch {matrix_error:.8g}')
        if freq_error > MAX_FREQ_ATOL_THZ:
            raise AssertionError(f'{slug} at {q}: frequency mismatch {freq_error:.8g} THz')
        max_matrix=max(max_matrix,matrix_error)
        max_frequency=max(max_frequency,freq_error)
        checks.append({'q_fractional': list(q), 'max_abs_dynamical_matrix_error':matrix_error,
                       'max_abs_frequency_error_thz':freq_error})
    import phonopy
    provenance={
        'evidence_status':'calculated',
        'source':f'Phonopy {entry["origin"]}, locked to {MANIFEST.as_posix()} and pinned Git commit',
        'source_kind':'external finite-displacement DFT example',
        'source_url':'https://github.com/phonopy/phonopy/tree/' +
                     json.loads(MANIFEST.read_text())['upstream_commit'] + '/example/' +
                     ('NaCl' if slug=='halite' else 'MgO'),
        'data_sha256':digests,
        'reference_calculator':'Phonopy '+phonopy.__version__,
        'comparison': 'non-NAC external dynamical matrices and frequencies; no experimental validation',
        'nonanalytic_correction':'intentionally omitted on both implementations',
    }
    result=export_result(model,provenance,
        [[0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.5]],
        points_per_segment=path_points, mesh=(mesh,)*3, temperature_k=300)
    # The result schema preserves source and declared status. A separate
    # certificate carries independent q comparisons and exact hashes.
    result['evidence']['benchmark']={'reference':'phonopy ' + phonopy.__version__,
         'source_sha256':digests,'comparison_points':len(Q_POINTS),
         'max_abs_frequency_error_thz':max_frequency,
         'max_abs_dynamical_matrix_error':max_matrix,
         'independent_comparison':'passed; computational, not experimental'}
    result['crystal']['mineral_name']=entry['name']
    result['crystal']['formula']=entry['formula']
    result['physics']['nonanalytic_correction']='omitted; polar LO-TO splitting not reproduced'
    result['physics']['force_constants_preprocessing']='Phonopy symmetrize_fc=True'
    out_dir.mkdir(parents=True,exist_ok=True)
    filename=f'{slug}_phonons.json'
    save_json(result,out_dir/filename)
    cert={'material':slug, 'name':entry['name'], 'formula':entry['formula'],
          'evidence_status':'calculated', 'independent_reference':'Phonopy '+phonopy.__version__,
          'from_supplied_forces':True, 'supercell_ifc_shape':list(phonon.force_constants.shape),
          'primitive_atoms':model.crystal.n_atoms,'periodic_ifc_blocks':len(model.blocks),
          'asr_error_ev_a2':asr,'max_abs_dynamical_matrix_error':max_matrix,
          'max_abs_frequency_error_thz':max_frequency,'qpoint_comparisons':checks,
          'source_sha256':digests, 'reference_file':filename,
          'lo_to_limitation':'nonanalytic correction off in both codes; not a full polar crystal spectrum'}
    return cert


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,default=ROOT/'benchmarks'/'upstream')
    p.add_argument('--output-dir',type=Path,default=ROOT/'benchmarks'/'generated')
    p.add_argument('--mesh',type=int,default=4)
    p.add_argument('--points',type=int,default=12)
    args=p.parse_args(argv)
    if not 2<=args.mesh<=40 or not 2<=args.points<=200:
        raise ValueError('mesh 2..40 and path points 2..200 required')
    manifest=json.loads(MANIFEST.read_text())
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
                         'evidence_status':r['evidence_status'],
                         'result_file':r['reference_file'],
                         'comparison_points':len(r['qpoint_comparisons']),
                         'max_frequency_error_thz':r['max_abs_frequency_error_thz']} for r in receipts]}
    save_json(catalog,args.output_dir/'catalog.json')
    save_json({'schema_version':'geoform.phonons.benchmark_receipts/1','results':receipts},
              args.output_dir/'benchmark_receipts.json')
    print('PASS: all Phonopy real-input computational cross-checks')

if __name__=='__main__':main()
