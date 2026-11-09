"""Replay the ZnS structure package without starting any simulation."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import warnings

from ase.data import atomic_numbers
from ase.io.cif import parse_cif
from ase.spacegroup import crystal
import numpy as np
import spglib

ROOT = Path(__file__).resolve().parent
MASS_AMU = {"Zn": 65.38, "S": 32.06}
F_PRIMITIVE = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
REFERENCE_GROUPS = {"zincblende": 216, "wurtzite": 186, "rocksalt": 225, "cmcm": 63}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def birch_murnaghan_pressure(volume, eos):
    eta = (eos["V0_angstrom3_per_formula_unit"] / volume) ** (1 / 3)
    return (1.5 * eos["B0_GPa"] * (eta**7 - eta**5)
            * (1 + .75 * (eos["B0_prime"] - 4) * (eta**2 - 1)))


def rocksalt_cell(phase):
    eos = phase["eos"]
    pressure = phase["conditions"]["pressure_GPa"]
    p_min, p_max = eos["source_fit_single_phase_interval_GPa"]
    require(p_min < pressure < p_max, "B1 reconstruction outside the allowed pressure interval")
    require(eos["B0_prime"] == 4, "Only the constrained BM fit is approved for replay")
    hi = eos["V0_angstrom3_per_formula_unit"]
    lo = hi / 2
    require(birch_murnaghan_pressure(lo, eos) > pressure > birch_murnaghan_pressure(hi, eos),
            "Cannot bracket B1 volume")
    for _ in range(70):
        midpoint = (lo + hi) / 2
        if birch_murnaghan_pressure(midpoint, eos) > pressure:
            lo = midpoint
        else:
            hi = midpoint
    volume = (lo + hi) / 2
    return [(volume * phase["Z_conventional"]) ** (1 / 3)] * 3


def inspect_geometry(lattice, positions, symbols, group, z, symprec=1e-5):
    require(np.isfinite(lattice).all() and np.linalg.det(lattice) > 0, "Invalid lattice")
    require(np.isfinite(positions).all(), "Invalid positions")
    require(Counter(symbols) == {"Zn": z, "S": z}, f"Expected {z} ZnS formula units")
    delta = positions[:, None, :] - positions[None, :, :]
    distances = np.max(np.abs(delta - np.rint(delta)), axis=2)
    require(np.count_nonzero(distances < 1e-7) == len(positions), "Duplicate periodic site")
    cell = (lattice, positions, [atomic_numbers[symbol] for symbol in symbols])
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
    require(dataset is not None and dataset.number == group, f"Expected space group {group}")
    return dataset


def positions_match(left, right, tolerance):
    delta = np.asarray(left) - np.asarray(right)
    return np.max(np.abs(delta - np.rint(delta))) <= tolerance


def cif_positions(block):
    return np.column_stack([block[f"_atom_site_fract_{axis}"] for axis in "xyz"])


def load_original(root, source):
    path = root / source["file"]
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == source["sha256"], f"{path.name}: changed source bytes")
    require(re.search(rb"\$Revision: " + str(source["cod_revision"]).encode() + rb" \$", raw),
            f"{path.name}: wrong source revision")
    with path.open("rb") as handle:
        blocks = list(parse_cif(handle))
    require(len(blocks) == 1, f"{path.name}: expected one CIF block")
    block = blocks[0]
    expected_cod = path.stem.removeprefix("avogadro-")
    require(block.name == expected_cod and str(block["_cod_database_code"]) == expected_cod,
            f"{path.name}: COD identity mismatch")
    tag = "_atom_site_occupancy"
    if source["occupancy_field"] == "explicit":
        require(block.get(tag) == [1.0, 1.0], f"{path.name}: occupancy mismatch")
    else:
        require(source["occupancy_field"] == "omitted" and tag not in block,
                f"{path.name}: unexpected occupancy provenance")
    require(block.get("_cell_formula_units_Z", source["Z_conventional"]) == source["Z_conventional"],
            f"{path.name}: CIF Z mismatch")
    require(block.get_symbols() == ["Zn", "S"], f"{path.name}: asymmetric composition mismatch")
    fractional = cif_positions(block)
    require(abs(fractional[1, 2] - source["sulfur_fractional_z"]) < 1e-8, f"{path.name}: sulfur z mismatch")
    with warnings.catch_warnings():
        # ASE 3.29 warns for the legacy hexagonal setting token; the expanded
        # geometry is checked independently below using spglib.
        warnings.filterwarnings("ignore", message="crystal system 'hexagonal' is not interpreted.*")
        require(block.get_spacegroup(False).no == source["space_group_number"], f"{path.name}: CIF group mismatch")
        atoms = block.get_atoms()
    np.testing.assert_allclose(atoms.cell.cellpar(), source["cell_lengths_angstrom"] + source["cell_angles_deg"],
                               atol=1e-9, rtol=0)
    inspect_geometry(atoms.cell.array, atoms.get_scaled_positions(), atoms.get_chemical_symbols(),
                     source["space_group_number"], source["Z_conventional"], symprec=1e-3)
    return block


def make_crystal(phase, original=None):
    group = phase["space_group_number"]
    lengths = rocksalt_cell(phase) if "eos" in phase else phase["cell_lengths_angstrom"]
    require(lengths is not None and all(length > 0 for length in lengths), "Invalid cell lengths")
    symbols, wyckoff, positions, occupancies = zip(*phase["sites"])
    require(symbols == ("Zn", "S") and occupancies == (1.0, 1.0), "Expected ordered ZnS sites")
    if original is not None:
        np.testing.assert_allclose(original.get_cell().cellpar(), lengths + phase["cell_angles_deg"], atol=1e-9, rtol=0)
        require(positions_match(positions, cif_positions(original), 4e-5), "Source site mismatch")
    conventional = crystal(symbols, positions, spacegroup=group, cellpar=lengths + phase["cell_angles_deg"])
    dataset = inspect_geometry(conventional.cell.array, conventional.get_scaled_positions(),
                               conventional.get_chemical_symbols(), group, phase["Z_conventional"])
    for symbol, expected in zip(symbols, wyckoff):
        site_indices = [i for i, element in enumerate(conventional.get_chemical_symbols()) if element == symbol]
        require({dataset.wyckoffs[i] for i in site_indices} == {expected[-1]}
                and len(site_indices) == int(expected[:-1]), f"{symbol}: incorrect Wyckoff orbit")
        require(all(dataset.equivalent_atoms[i] == dataset.equivalent_atoms[site_indices[0]]
                    for i in site_indices), f"{symbol}: expected one orbit")

    lattice = conventional.cell.array
    if group in (216, 225):
        lattice = F_PRIMITIVE @ lattice
        positions = np.asarray(positions) @ np.linalg.inv(F_PRIMITIVE)
    else:
        symbols = conventional.get_chemical_symbols()
        positions = conventional.get_scaled_positions()
    positions = np.mod(positions, 1)
    inspect_geometry(lattice, positions, symbols, group, phase["Z_primitive"])
    np.testing.assert_allclose(np.linalg.det(lattice) / phase["Z_primitive"],
                               conventional.get_volume() / phase["Z_conventional"], rtol=1e-12)
    return {
        "crystal": {"lattice_angstrom": lattice.round(12).tolist(),
                    "fractional_positions": positions.round(12).tolist(),
                    "masses_amu": [MASS_AMU[symbol] for symbol in symbols],
                    "symbols": list(symbols)},
        "space_group_number": group,
        "formula": "ZnS",
        "formula_units_in_crystal": phase["Z_primitive"],
        "volume_angstrom3_per_formula_unit": round(float(np.linalg.det(lattice) / phase["Z_primitive"]), 10),
        "evidence": phase["evidence"],
    }


def replay(root=ROOT, write=False):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    phases = manifest["phases"]
    require(manifest["schema"] == "geoform.zns-structure-foundation/1"
            and manifest["formula"] == "ZnS", "Unknown ZnS manifest")
    require(manifest["four_phase_study_ready"] is False, "Four-phase gate must remain blocked")
    require({name: phase["space_group_number"] for name, phase in phases.items()} == REFERENCE_GROUPS,
            "Wrong phase inventory")
    originals = {name: load_original(root, source) for name, source in manifest["originals"].items()}
    cmcm = phases["cmcm"]
    require(cmcm["sites"] is None and cmcm["crystal"] is None and cmcm["missing_evidence"],
            "Cmcm must stay evidence-only")
    generated = root / "crystals"
    expected_names = {f"{name}.json" for name in phases if name != "cmcm"}
    actual_names = {path.name for path in generated.iterdir()} if generated.exists() else set()
    require(actual_names <= expected_names, "Unexpected coordinate file")
    require(write or actual_names == expected_names, "Missing coordinate file")
    for name, phase in phases.items():
        if name == "cmcm":
            continue
        source = originals[phase["source_cif"]] if phase["source_cif"] is not None else None
        result = make_crystal(phase, source)
        result["schema"] = "geoform.zns-crystal-fragment/1"
        result["phase"] = name
        result["conditions"] = phase["conditions"]
        result["source_manifest"] = "../manifest.json"
        path = generated / f"{name}.json"
        if write:
            generated.mkdir(exist_ok=True)
            path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        else:
            require(json.loads(path.read_text(encoding="utf-8")) == result, f"{name}: stale crystal")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Rebuild the three approved-for-review geometry fragments")
    args = parser.parse_args()
    replay(write=args.write)
    print("Validated five source CIFs, three geometry fragments, and the Cmcm evidence gate.")
