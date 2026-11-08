"""Validate MDR calculation packages and produce signed harmonic evidence."""
from __future__ import annotations

from collections import Counter
from io import BytesIO
import lzma
import math
import re
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np
import phonopy
import spglib
import yaml
from phonopy.interface.phonopy_yaml import load_phonopy_yaml
from phonopy.phonon.dos import TotalDos
from phonopy.phonon.thermal_properties import ThermalProperties
from phonopy.physical_units import get_physical_units

from .core import THZ_FACTOR, density_of_states, shifted_q_mesh, signed_frequencies, thermodynamic_curve
from .phonondb import digest
from .phonopy_native import from_phonopy


MAX_YAML_BYTES = 24 * 1024**2
MAX_PRIMITIVE_ATOMS = 64
MAX_SUPERCELL_ATOMS = 650
GAMMA_ACOUSTIC_ATOL_THZ = 1e-3
Q_CHECKS = [(0, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5),
            (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5), (0.5, 0.5, 0.5),
            (0.17, 0.31, 0.09), (0.41, 0.02, 0.28), (0.13, 0.44, 0.27),
            (0.07, 0.26, 0.45), (-0.22, 0.14, 0.31), (1.17, 0.31, 0.09),
            (0.25, 0, 0.5), (0.125, 0.25, 0.375)]
VERTICES = np.asarray([(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0),
                       (0.5, 0.5, 0.5), (0, 0, 0)])
MATRIX_ATOL = 5e-6
FREQUENCY_ATOL = 1e-3


class PackageRejected(ValueError):
    """A record is outside the verified input envelope."""


def require(condition, reason):
    if not condition:
        raise PackageRejected(reason)


def read_package(data):
    with ZipFile(BytesIO(data)) as archive:
        members = [info for info in archive.infolist()
                   if not info.is_dir() and info.filename.rsplit('/', 1)[-1] == 'phonopy_params.yaml.xz']
        require(len(members) == 1, 'Expected one phonopy_params.yaml.xz member')
        require(members[0].file_size <= MAX_YAML_BYTES, 'Compressed YAML exceeds size budget')
        packed = archive.read(members[0])
    stream = lzma.LZMADecompressor(format=lzma.FORMAT_XZ, memlimit=128 * 1024**2)
    decoded = stream.decompress(packed, max_length=MAX_YAML_BYTES + 1)
    require(len(decoded) <= MAX_YAML_BYTES and stream.eof and not stream.unused_data,
            'XZ input is incomplete, trailing or exceeds size budget')
    document = yaml.load(decoded, Loader=yaml.CSafeLoader)
    require(isinstance(document, dict), 'Expected YAML mapping')
    return document, digest(packed), digest(decoded)


def reduced(counts):
    common = math.gcd(*counts.values())
    return {symbol: count // common for symbol, count in sorted(counts.items())}


def indexed_formula(formula):
    tokens = re.findall(r'[A-Z][a-z]?|[0-9]+|[()]', formula)
    require(tokens and ''.join(tokens) == formula, 'Unsupported index formula syntax')
    groups = [Counter()]
    for i, token in enumerate(tokens):
        if token == '(':
            groups.append(Counter())
            continue
        if token.isdecimal():
            require(i > 0 and tokens[i - 1] != '(' and not tokens[i - 1].isdecimal(),
                    'Unexpected formula multiplier')
            continue
        multiplier = int(tokens[i + 1]) if i + 1 < len(tokens) and tokens[i + 1].isdecimal() else 1
        require(multiplier > 0, 'Nonpositive formula count')
        if token == ')':
            require(len(groups) > 1 and groups[-1], 'Unmatched or empty formula group')
            completed = groups.pop()
            groups[-1].update({symbol: count * multiplier for symbol, count in completed.items()})
        else:
            groups[-1][token] += multiplier
    require(len(groups) == 1 and groups[0], 'Unmatched formula group')
    return reduced(groups[0])


def verify_cells(stored, rebuilt):
    require(stored is not None and len(stored) == len(rebuilt), 'Stored/rebuilt cell size mismatch')
    require(list(stored.symbols) == list(rebuilt.symbols), 'Stored/rebuilt site order or species mismatch')
    require(np.allclose(stored.cell, rebuilt.cell, rtol=0, atol=2e-7), 'Stored/rebuilt lattice mismatch')
    require(np.allclose(stored.masses, rebuilt.masses, rtol=0, atol=2e-6), 'Stored/rebuilt masses mismatch')
    offset = stored.scaled_positions - rebuilt.scaled_positions
    require(np.allclose(offset, np.rint(offset), rtol=0, atol=2e-6),
            'Stored/rebuilt fractional site order or positions mismatch')


def load_harmonic(document, row):
    unit = document.get('physical_unit', {})
    require(unit.get('length') == 'angstrom' and unit.get('atomic_mass', '').lower() == 'amu'
            and unit.get('force_constants') == 'eV/angstrom^2', 'Unknown source physical units')
    for key, limit in [('unit_cell', MAX_SUPERCELL_ATOMS), ('primitive_cell', MAX_PRIMITIVE_ATOMS),
                       ('supercell', MAX_SUPERCELL_ATOMS)]:
        cell = document.get(key, {})
        count = len(cell.get('points', []))
        require(0 < count <= limit, f'{key} has {count} sites; supported range is 1..{limit}')
        lattice = np.asarray(cell.get('lattice'), dtype=float)
        require(lattice.shape == (3, 3) and np.isfinite(lattice).all()
                and abs(np.linalg.det(lattice)) > 1e-8, f'{key} invalid lattice')
        for point in cell['points']:
            coordinates = np.asarray(point.get('coordinates'), dtype=float)
            require(coordinates.shape == (3,) and np.isfinite(coordinates).all()
                    and np.isfinite(point.get('mass', np.nan)) and point['mass'] > 0,
                    f'{key} invalid site coordinates or mass')

    smat = np.asarray(document.get('supercell_matrix'), dtype=float)
    require(smat.shape == (3, 3) and np.isfinite(smat).all() and np.all(smat == np.rint(smat))
            and np.max(np.abs(smat)) <= MAX_SUPERCELL_ATOMS, 'Invalid supercell matrix')
    ns = len(document['supercell']['points'])
    require(round(abs(np.linalg.det(smat))) * len(document['unit_cell']['points']) == ns,
            'Supercell matrix does not match recorded site count')

    metadata = document.get('phonopy', {})
    symprec = float(metadata.get('symmetry_tolerance', 1e-5))
    factor = float(metadata['frequency_unit_conversion_factor'])
    require(0 < symprec <= 1e-3 and abs(factor - THZ_FACTOR) < 1e-4,
            'Unexpected source symmetry tolerance or THz conversion')
    require(document['space_group']['number'] == row['space_group_number'],
            'Index and YAML space groups disagree')
    parsed = load_phonopy_yaml(document)
    primitive_matrix = parsed.primitive_matrix
    if primitive_matrix is None:
        primitive_matrix = (parsed.primitive.cell @ np.linalg.inv(parsed.unitcell.cell)).T
    phonon = phonopy.Phonopy(parsed.unitcell, smat.astype(int), primitive_matrix=primitive_matrix,
                             calculator='vasp', symprec=symprec, lang='C')
    verify_cells(parsed.unitcell, phonon.unitcell)
    verify_cells(parsed.primitive, phonon.primitive)
    verify_cells(parsed.supercell, phonon.supercell)
    composition = reduced(Counter(phonon.unitcell.symbols))
    require(composition == indexed_formula(row['formula']), 'Index and cell stoichiometry disagree')
    for cell in [phonon.primitive, phonon.supercell]:
        require(reduced(Counter(cell.symbols)) == composition, 'Cells disagree on stoichiometry')
    symmetry = spglib.get_symmetry_dataset(
        (phonon.unitcell.cell, phonon.unitcell.scaled_positions, phonon.unitcell.numbers),
        symprec=symprec)
    require(symmetry is not None and symmetry.number == row['space_group_number'],
            'Recomputed space group disagrees with index')

    # A matching fingerprint flags a structure for review. It is not a
    # crystallographic equivalence proof and never merges calculations.
    types = symmetry.std_types
    order = sorted(set(types), key=lambda t: (np.count_nonzero(types == t), t))
    motifs = []
    for origin in symmetry.std_positions[types == order[0]]:
        wrapped = np.round((symmetry.std_positions - origin) % 1, 5) % 1
        motifs.append(sorted((int(t), *xyz.tolist()) for t, xyz in zip(types, wrapped)))
    metric = symmetry.std_lattice @ symmetry.std_lattice.T
    fingerprint = digest(repr((int(symmetry.number), np.round(metric, 4).tolist(), min(motifs))).encode())

    source = {}
    if parsed.force_constants is not None:
        phonon.force_constants = parsed.force_constants
        source['ifc_origin'] = 'embedded force constants; no projection'
    else:
        dataset = parsed.dataset
        require(dataset is not None and dataset.get('natom') == ns and dataset.get('first_atoms'),
                'Missing supported Type-I displacement/force arrays')
        for snapshot in dataset['first_atoms']:
            forces = np.asarray(snapshot.get('forces'), dtype=float)
            displacement = np.asarray(snapshot.get('displacement'), dtype=float)
            require(forces.shape == (ns, 3) and np.isfinite(forces).all()
                    and displacement.shape == (3,) and np.isfinite(displacement).all()
                    and np.linalg.norm(displacement) > 0 and 0 <= snapshot['number'] < ns,
                    'Incomplete or invalid displacement/force snapshot')
        phonon.dataset = dataset
        phonon.produce_force_constants(calculate_full_force_constants=False,
                                       fc_calculator='traditional', show_drift=False)
        phonon.symmetrize_force_constants(use_symfc_projector=True, show_drift=False)
        source.update(ifc_origin='embedded Type-I forces; traditional finite differences with symfc projection',
                      force_snapshots=len(dataset['first_atoms']))
    require(phonon.nac_params is None, 'Unexpected NAC on short-range comparison')
    require(np.max(np.abs(phonon.force_constants)) > 0, 'All-zero IFCs')
    model = from_phonopy(phonon, tolerance=4e-6)
    asr = model.acoustic_sum_rule_error()
    require(asr <= 1e-6, 'IFC acoustic sum rule exceeds 1e-6 eV/angstrom^2')
    return phonon, model, {
        **source, 'source_phonopy_version': str(metadata['version']),
        'source_nac_data_present': parsed.nac_params is not None,
        'source_frequency_factor': factor, 'source_symmetry_tolerance_angstrom': symprec,
        'ifc_shape': list(phonon.force_constants.shape), 'asr_max_abs_ev_a2': asr,
        'reduced_composition': composition, 'space_group_number': int(symmetry.number),
        'space_group_symbol': symmetry.international,
        'structure_fingerprint': fingerprint, 'primitive_atoms': len(phonon.primitive), 'supercell_atoms': ns,
        'primitive_volume_angstrom3': model.crystal.volume_angstrom3,
    }


def sample_harmonic(phonon, model, mesh_size):
    mesh = (mesh_size,) * 3
    path = np.concatenate([np.linspace(a, b, 12, endpoint=i == len(VERTICES) - 2)
                           for i, (a, b) in enumerate(zip(VERTICES[:-1], VERTICES[1:]))])
    q = np.vstack([Q_CHECKS, path, shifted_q_mesh(mesh)])
    # Bound dense matrix storage even when a user requests a larger mesh.
    check_matrices = model.dynamical_matrices(Q_CHECKS)
    eigenvalues = np.concatenate([
        np.linalg.eigvalsh(check_matrices),
        *[np.linalg.eigvalsh(model.dynamical_matrices(q[i:i+32]))
          for i in range(len(Q_CHECKS), len(q), 32)],
    ])
    frequencies = signed_frequencies(eigenvalues)

    reference = phonon.run_qpoints(Q_CHECKS, with_dynamical_matrices=True)
    matrix_error = float(np.max(np.abs(check_matrices - reference.dynamical_matrices)))
    require(matrix_error <= MATRIX_ATOL, f'Phonopy dynamical matrix discrepancy {matrix_error:.8g}')
    external = np.sort(phonon.run_qpoints(q).frequencies, axis=1)
    frequency_error = float(np.max(np.abs(frequencies - external)))
    require(frequency_error <= FREQUENCY_ATOL, f'Phonopy signed-frequency discrepancy {frequency_error:.8g} THz')

    # Ignore a small signed Gamma artifact only if both eigenspaces identify
    # it as a mass-weighted uniform translation. Keep all raw signed values.
    negative = (frequencies < 0) | (external < 0)
    ignored = 0
    translations = np.kron(np.sqrt(model.crystal.masses_amu)[:, None], np.eye(3))
    translations /= np.sqrt(model.crystal.masses_amu.sum())
    for index in np.flatnonzero(np.all(q == np.rint(q), axis=1)):
        _, ours = np.linalg.eigh(model.dynamical_matrix(q[index]))
        # With NAC off, the independent Gamma eigenspace comes from Phonopy.
        ext = phonon.run_qpoints([q[index]], with_eigenvectors=True).eigenvectors[0]
        ours_weight = (np.abs(translations.T @ ours)**2).sum(axis=0)
        ext_weight = (np.abs(translations.T @ ext)**2).sum(axis=0)
        roundoff = ((np.abs(frequencies[index]) <= GAMMA_ACOUSTIC_ATOL_THZ)
                    & (np.abs(external[index]) <= GAMMA_ACOUSTIC_ATOL_THZ)
                    & (ours_weight > 0.999) & (ext_weight > 0.999))
        ignored += int(np.count_nonzero(negative[index] & roundoff))
        negative[index] &= ~roundoff
    offset = len(Q_CHECKS) + len(path)
    grid = frequencies[offset:]
    available = not negative.any() and bool(np.all(grid > 0) and np.all(external[offset:] > 0))
    stability = ('sampled_imaginary_modes' if negative.any() else
                 'sampled_positive' if available else 'sampled_zero_mesh_modes')
    minimum = np.unravel_index(np.argmin(frequencies), frequencies.shape)[0]
    comparison = {
        'matrix_q_points': len(Q_CHECKS), 'frequency_q_points': len(q),
        'max_abs_matrix_error_ev_a2_amu': matrix_error, 'max_abs_frequency_error_thz': frequency_error,
        'stability': stability, 'signed_minimum_thz': float(frequencies.min()),
        'minimum_q_fractional': q[minimum].tolist(), 'imaginary_sampled_modes': int(negative.sum()),
        'reference_signed_minimum_thz': float(external.min()),
        'gamma_acoustic_roundoff_modes_excluded': ignored,
    }
    distances = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path, axis=0)
                                                @ model.crystal.reciprocal_vectors, axis=1))]
    dispersion = {'q_fractional': path.tolist(), 'distance_inverse_angstrom': distances.tolist(),
                  'frequency_thz_signed': frequencies[len(Q_CHECKS):offset].tolist(),
                  'path_description': 'Generic primitive reciprocal-coordinate path; no symmetry labels',
                  'band_order': 'ascending eigenvalue (not continuous mode tracking)'}
    result = {
        'physics': {'method': 'harmonic real-space IFC diagonalization',
                    'q_mesh': list(mesh), 'nonanalytic_correction': 'omitted on both sides',
                    'sampled_stability': comparison},
        'dispersion': dispersion,
        'density_of_states': {'status': f'withheld: {stability}'},
        'thermodynamics': {'status': f'withheld: {stability}'},
    }
    if available:
        points = max(400, math.ceil((float(grid.max()) + 1.0) / 0.04) + 1)
        require(points <= 8000, 'DOS frequency range exceeds output budget')
        x, y = density_of_states(grid, points=points)
        temperatures = [0, 50, 100, 200, 300, 500, 800, 1200]
        curve = thermodynamic_curve(grid, temperatures)
        comparison.update(compare_derived(phonon, q[offset:], external[offset:], x, y, curve))
        result.update(density_of_states={'frequency_thz': x.tolist(), 'states_per_thz_cell': y.tolist(),
                                         'gaussian_sigma_thz': 0.1},
                      thermodynamics=curve[4], thermodynamic_curve=curve)
    return comparison, result


def compare_derived(phonon, q, frequencies, axis, dos, curve):
    """Use Phonopy's DOS and thermal kernels on its own unreduced q-point frequencies."""
    full_mesh = SimpleNamespace(frequencies=frequencies, qpoints=q, gamma_index=None,
                                weights=np.ones(len(q), dtype=np.int64),
                                dynamical_matrix=phonon.dynamical_matrix)
    external_dos = TotalDos(full_mesh, sigma=0.1, lang='C')
    external_dos.set_draw_area(float(axis[0]), float(axis[-1]), float(axis[1] - axis[0]))
    external_dos.run()
    require(np.allclose(axis, external_dos.frequency_points, rtol=0, atol=1e-10),
            'DOS frequency axes disagree')
    dos_error = float(np.max(np.abs(dos - external_dos.dos)))
    require(dos_error <= 1e-4 * max(1, float(dos.max())), 'Phonopy Gaussian DOS discrepancy')
    dos_integral = float(np.trapezoid(dos, axis))
    require(abs(dos_integral - frequencies.shape[1]) <= 1e-5 * frequencies.shape[1],
            'DOS mode-count normalization failed')

    thermal = ThermalProperties(full_mesh, lang='C')
    thermal.temperatures = [point['temperature_k'] for point in curve]
    thermal.run(lang='Python')
    _, free_energy, entropy, cv = thermal.thermal_properties
    ev_to_kj_mol = get_physical_units().EvTokJmol
    free_energy /= ev_to_kj_mol
    entropy /= 1000 * ev_to_kj_mol
    cv /= 1000 * ev_to_kj_mol
    expected = {
        'free_energy_ev_cell': free_energy,
        'internal_energy_ev_cell': free_energy + thermal.temperatures * entropy,
        'entropy_ev_k_cell': entropy,
        'cv_ev_k_cell': cv,
    }
    errors = {key: float(np.max(np.abs(reference - [point[key] for point in curve])))
              for key, reference in expected.items()}
    for key, error in errors.items():
        # Phonopy uses older CODATA constants. Bound the discrepancy per
        # primitive atom so the comparison has the same scale across phases.
        per_atom = error / (frequencies.shape[1] // 3)
        require(per_atom <= (5e-9 if '_k_cell' in key else 2e-6),
                f'Phonopy thermal discrepancy for {key}: {error:.8g}')
    return {'max_abs_dos_error_states_per_thz_cell': dos_error, 'dos_integral_states_per_cell': dos_integral,
            'max_abs_thermal_errors': errors}
