"""Replay the research-only 96 GPa ZnS Cmcm geometry and sensitivity checks."""

import argparse
from itertools import product
import json

from ase import Atoms
from ase.data import atomic_numbers
from ase.neighborlist import neighbor_list
from ase.spacegroup import crystal
import numpy as np
import spglib

from validate_structures import MASS_AMU, ROOT, inspect_geometry, replay, require

PATH = ROOT / "candidates" / "cmcm_jyoti_96gpa.json"
SYMBOLS = ("Zn", "Zn", "S", "S")
C_PRIMITIVE = np.array([[.5, -.5, 0], [.5, .5, 0], [0, 0, 1]])
HALL_NUMBER = 298


def neighbor_shells(lattice, positions, symbols):
    atoms = Atoms(symbols, cell=lattice, scaled_positions=positions, pbc=True)
    centers, partners, distances = neighbor_list("ijd", atoms, 4.1)
    shells = {}
    for index, symbol in enumerate(symbols):
        mask = centers == index
        unlike = np.sort(distances[mask & (atoms.numbers[partners] != atoms.numbers[index])])
        like = distances[mask & (atoms.numbers[partners] == atoms.numbers[index])]
        require(len(unlike) >= 7 and len(like), "Incomplete periodic neighbor search")
        shell = np.r_[unlike[:7], min(like)]
        if symbol in shells:
            np.testing.assert_allclose(shell, shells[symbol], atol=1e-10, rtol=0)
        shells[symbol] = shell
    np.testing.assert_allclose(shells["Zn"][:7], shells["S"][:7], atol=1e-10, rtol=0)
    return shells


def analyze_geometry(lengths, y_zn, y_s):
    conventional_lattice = np.diag(lengths)
    sites = [(0, y_zn, .25), (0, y_s, .25)]
    conventional = crystal(("Zn", "S"), sites, spacegroup=63, cell=conventional_lattice)
    conv_symbols = conventional.get_chemical_symbols()
    conv_positions = conventional.get_scaled_positions()

    # Two representatives of each 4c orbit survive after removing C centering.
    representatives = np.array([[0, y_zn, .25], [0, -y_zn, .75],
                                [0, y_s, .25], [0, -y_s, .75]])
    positions = np.mod(representatives @ np.linalg.inv(C_PRIMITIVE), 1)
    lattice = C_PRIMITIVE @ conventional_lattice
    expanded = np.vstack([representatives, representatives + [.5, .5, 0]])
    for symbol in ("Zn", "S"):
        expected = {tuple(row) for row in np.mod(expanded, 1).round(10)[
            np.tile(np.array(SYMBOLS) == symbol, 2)]}
        actual = {tuple(row) for row in conv_positions.round(10)[np.array(conv_symbols) == symbol]}
        require(expected == actual and len(actual) == 4, f"Incorrect {symbol} orbit")

    cell = (lattice, positions, [atomic_numbers[symbol] for symbol in SYMBOLS])
    detected = spglib.get_symmetry_dataset(cell, symprec=1e-5)
    require(detected is not None, "Symmetry detection failed")
    group = detected.number
    conventional_data = inspect_geometry(conventional_lattice, conv_positions, conv_symbols, group, 4)
    primitive_data = inspect_geometry(lattice, positions, SYMBOLS, group, 2)
    if group == 63:
        for dataset, symbols in ((conventional_data, conv_symbols), (primitive_data, SYMBOLS)):
            require(dataset.hall_number == HALL_NUMBER and set(dataset.wyckoffs) == {"c"},
                    "Expected standard Cmcm with 4c sites")
            for symbol in ("Zn", "S"):
                orbit = {i for i, element in zip(dataset.equivalent_atoms, symbols) if element == symbol}
                require(len(orbit) == 1, f"Expected one {symbol} orbit")
    primitive = spglib.find_primitive(cell, symprec=1e-5)
    np.testing.assert_allclose(np.linalg.det(lattice) * 2, np.linalg.det(conventional_lattice), rtol=1e-12)
    expected_atoms = 2 if group == 69 else 4
    require(primitive is not None and len(primitive[1]) == expected_atoms,
            "Primitive normalization failed")
    np.testing.assert_allclose(abs(np.linalg.det(primitive[0])) / expected_atoms,
                               np.linalg.det(lattice) / 4, rtol=1e-12)

    shells = neighbor_shells(lattice, positions, SYMBOLS)
    conv_shells = neighbor_shells(conventional_lattice, conv_positions, conv_symbols)
    for symbol in shells:
        np.testing.assert_allclose(shells[symbol], conv_shells[symbol], atol=1e-10, rtol=0)
    return lattice, positions, group, shells


def sensitivity(lengths, coordinates):
    # If the displayed values were rounded to nearest 0.1, these are the
    # corresponding half-digit offsets. The source gives no uncertainty model.
    offsets = [-.05, -.025, 0, .025, .05]
    rows = []
    for dz, ds in product(offsets, repeat=2):
        y_zn, y_s = np.round(np.asarray(coordinates) + [dz, ds], 10)
        _, _, group, shells = analyze_geometry(lengths, y_zn, y_s)
        rows.append((y_zn, y_s, group, shells["Zn"][0], shells["Zn"][5], shells["Zn"][6]))
    return rows


def build_candidate(manifest):
    phase = manifest["phases"]["cmcm"]
    model = manifest["research_candidates"]["cmcm_jyoti_96gpa"]
    require(model["research_only"] is True and model["occupancies"] == {"Zn": 1.0, "S": 1.0},
            "Expected ordered research-only candidate")
    require(model["pressure_GPa"] == phase["conditions"]["pressure_GPa"] == 96,
            "Inconsistent pressure tags")
    lengths = phase["cell_lengths_angstrom"]
    coordinates = (model["fractional_y"]["Zn"], model["fractional_y"]["S"])
    lattice, positions, group, shells = analyze_geometry(lengths, *coordinates)
    require(group == 63, "Nominal candidate must be Cmcm")
    require(shells["Zn"][5] < shells["Zn"][6], "Expected distinct six-neighbor shell")
    return {
        "schema": "geoform.zns-research-crystal-candidate/1",
        "candidate": "cmcm_jyoti_96gpa",
        "evidence": "theoretical_internal_coordinates_on_measured_lattice",
        "research_only": True,
        "pressure_GPa": model["pressure_GPa"],
        "coordinate_model_temperature_K": model["coordinate_model_temperature_K"],
        "source_manifest": "../manifest.json",
        "crystal": {
            "lattice_angstrom": lattice.round(12).tolist(),
            "fractional_positions": positions.round(12).tolist(),
            "masses_amu": [MASS_AMU[symbol] for symbol in SYMBOLS],
            "symbols": list(SYMBOLS),
        },
        "space_group_number": group,
        "hall_number": HALL_NUMBER,
        "formula": "ZnS",
        "formula_units_in_crystal": 2,
        "volume_angstrom3_per_formula_unit": round(float(np.linalg.det(lattice) / 2), 10),
        "first_seven_unlike_neighbors_angstrom": shells["Zn"][:7].round(8).tolist(),
        "nearest_like_neighbor_angstrom": {symbol: round(float(shell[-1]), 8) for symbol, shell in shells.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Regenerate the research candidate JSON")
    args = parser.parse_args()
    manifest = replay()
    candidate = build_candidate(manifest)
    if args.write:
        PATH.write_text(json.dumps(candidate, indent=2) + "\n", encoding="utf-8")
    else:
        require(json.loads(PATH.read_text(encoding="utf-8")) == candidate, "Stale Cmcm candidate")

    model = manifest["research_candidates"]["cmcm_jyoti_96gpa"]
    lengths = manifest["phases"]["cmcm"]["cell_lengths_angstrom"]
    coordinates = (model["fractional_y"]["Zn"], model["fractional_y"]["S"])
    print("y_Zn  y_S   space_group  d1(Zn-S)  d6(Zn-S)  d7(Zn-S) [Å]")
    for y_zn, y_s, group, d1, d6, d7 in sensitivity(lengths, coordinates):
        print(f"{y_zn:.3f} {y_s:.3f} {group:11d}  {d1:8.5f}  {d6:8.5f}  {d7:8.5f}")


if __name__ == "__main__":
    main()
