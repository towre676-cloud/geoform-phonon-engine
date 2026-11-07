"""Cross-check a pinned Phonopy example against its original VASP force runs."""
from collections import Counter
import re
import tarfile
from xml.etree import ElementTree

import numpy as np


def fractional_error(actual, expected):
    return float(np.max(np.abs((actual - expected + 0.5) % 1 - 0.5)))


def check_source(path, entry, phonon):
    import spglib
    from phonopy.file_IO import parse_FORCE_SETS
    from phonopy.interface.phonopy_yaml import PhonopyYaml
    from phonopy.interface.vasp import read_vasp

    declared = entry['structure']
    source = PhonopyYaml().read(path / 'phonopy_disp.yaml')
    cell = read_vasp(path / 'POSCAR-unitcell')
    sym = spglib.get_symmetry_dataset(
        (cell.cell, cell.scaled_positions, cell.numbers), symprec=1e-5)
    expected_wyckoff = set(declared['wyckoff_by_element'].items())
    if (sym is None or sym.number != declared['number'] or
            sym.international != declared['space_group'] or
            set(zip(cell.symbols, sym.wyckoffs)) != expected_wyckoff or
            cell.symbols != source.unitcell.symbols or
            not np.allclose(cell.cell, source.unitcell.cell, rtol=0, atol=1e-8) or
            fractional_error(cell.scaled_positions, source.unitcell.scaled_positions) > 1e-8):
        raise ValueError(f'{path.name}: POSCAR/YAML structure or phase mismatch')

    factors = {s: int(n or '1') for s, n in re.findall(r'([A-Z][a-z]?)(\d*)', entry['formula'])}
    counts = Counter(phonon.primitive.symbols)
    if set(counts) != set(factors) or any(
            count * sum(factors.values()) != factors[s] * len(phonon.primitive)
            for s, count in counts.items()):
        raise ValueError(f'{path.name}: primitive stoichiometry mismatch')
    supercell = source.supercell
    force_sets = parse_FORCE_SETS(filename=path / 'FORCE_SETS')['first_atoms']
    displacements = source.dataset['first_atoms']
    if (not np.array_equal(phonon.supercell_matrix, np.diag(declared['supercell_diagonal'])) or
            len(phonon.primitive) != declared['primitive_atoms'] or
            len(supercell) != declared['supercell_atoms'] or
            len(force_sets) != len(displacements) or len(force_sets) != declared['displacements'] or
            phonon.supercell.symbols != supercell.symbols or
            not np.allclose(phonon.supercell.cell, supercell.cell, rtol=0, atol=1e-8) or
            fractional_error(phonon.supercell.scaled_positions, supercell.scaled_positions) > 1e-8):
        raise ValueError(f'{path.name}: supercell, primitive or displacement mismatch')

    runs = []
    excluded = []
    with tarfile.open(path / 'vaspruns.tar.lzma') as archive:
        for member in sorted(archive.getmembers(), key=lambda member: member.name):
            if not member.isfile():
                continue
            match = re.fullmatch(r'vaspruns/vasprun(?:\.xml)?-(\d{3})(?:\.xml)?', member.name)
            if not match or int(match[1]) == 0:
                excluded.append(member.name)
                continue
            index = int(match[1])
            if index != len(runs) + 1 or index > len(force_sets):
                raise ValueError(f'{path.name}: unexpected VASP displacement run {index}')
            with archive.extractfile(member) as stream:
                xml = ElementTree.parse(stream)
            run = check_force_run(xml, member.name, force_sets[index - 1],
                                  displacements[index - 1], supercell, entry['vasp'])
            runs.append(run)
    if len(runs) != len(force_sets):
        raise ValueError(f'{path.name}: missing VASP force runs')
    return {
        'space_group': sym.international, 'space_group_number': int(sym.number),
        'wyckoff_by_element': declared['wyckoff_by_element'],
        'conventional_atoms': len(cell), 'source_yaml_basis_atoms': len(source.primitive),
        'primitive_atoms': len(phonon.primitive), 'supercell_atoms': len(supercell),
        'supercell_matrix': phonon.supercell_matrix.tolist(),
        'primitive_matrix': phonon.primitive_matrix.tolist(),
        'born_file_present': entry['vasp']['born_file_present'],
        'excluded_archive_members': excluded, 'vasp_force_runs': runs,
    }


def check_force_run(xml, name, force_set, displacement, supercell, settings):
    def vectors(xpath):
        return np.array([node.text.split() for node in xml.findall(xpath)], dtype=float)

    if (force_set['number'] != displacement['number'] or
            not np.allclose(force_set['displacement'], displacement['displacement'], rtol=0, atol=1e-12)):
        raise ValueError(f'{name}: YAML/FORCE_SETS displacement mismatch')

    expected = supercell.scaled_positions.copy()
    expected[displacement['number']] += displacement['displacement'] @ np.linalg.inv(supercell.cell)
    initial = './/structure[@name="initialpos"]'
    lattice_error = float(np.max(np.abs(vectors(initial + '/crystal/varray[@name="basis"]/v') - supercell.cell)))
    positions_error = fractional_error(vectors(initial + '/varray[@name="positions"]/v'), expected)
    symbols = [node.text.strip() for node in xml.findall('./atominfo/array[@name="atoms"]/set/rc/c[1]')]
    force_blocks = xml.findall('./calculation/varray[@name="forces"]')
    forces = np.array([node.text.split() for node in force_blocks[-1]], dtype=float)
    error = float(np.max(np.abs(forces - force_set['forces'])))
    encut = float(xml.findtext('./incar/i[@name="ENCUT"]'))
    kmesh = list(map(int, xml.findtext('./kpoints/generation/v[@name="divisions"]').split()))
    if (symbols != supercell.symbols or lattice_error > 1e-7 or positions_error > 1e-7
            or not np.array_equal(forces, force_set['forces'])
            or encut != settings['encut_ev'] or kmesh != settings['force_k_mesh']):
        raise ValueError(f'{name}: archived VASP geometry, forces or settings mismatch')
    return {
        'file': name, 'version': xml.findtext('./generator/i[@name="version"]').strip(),
        'encut_ev': encut, 'force_k_mesh': kmesh,
        'smearing_method': (xml.findtext('./incar/i[@name="ISMEAR"]') or '').strip() or None,
        'smearing_width_ev': (xml.findtext('./incar/i[@name="SIGMA"]') or '').strip() or None,
        'external_pressure_kbar': (xml.findtext('./incar/i[@name="PSTRESS"]') or '').strip() or None,
        'gga_tag': (xml.findtext('.//parameters//i[@name="GGA"]') or '').strip(),
        'max_lattice_error_angstrom': lattice_error,
        'max_fractional_position_error': positions_error,
        'max_force_error_ev_angstrom': error,
    }
