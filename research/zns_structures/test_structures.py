import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from ase.neighborlist import neighbor_list
from ase import Atoms

from validate_structures import ROOT, birch_murnaghan_pressure, make_crystal, replay


class ZnSFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = replay()
        cls.phases = cls.manifest["phases"]

    def read_crystal(self, name):
        return json.loads((ROOT / "crystals" / f"{name}.json").read_text(encoding="utf-8"))

    def test_phase_reference_volumes_and_coordination(self):
        expectations = {
            "zincblende": (1, 39.5697414856, 4, 2.5),
            "wurtzite": (2, 39.6153924759, 4, 2.5),
            "rocksalt": (1, 26.4359490067, 6, 2.5),
        }
        for phase, (z, volume, coordination, cutoff) in expectations.items():
            with self.subTest(phase=phase):
                entry = self.read_crystal(phase)
                self.assertEqual(entry["formula_units_in_crystal"], z)
                self.assertAlmostEqual(entry["volume_angstrom3_per_formula_unit"], volume, places=9)
                self.assertNotIn("force_constants", entry)
                crystal = entry["crystal"]
                atoms = Atoms(crystal["symbols"], cell=crystal["lattice_angstrom"],
                              scaled_positions=crystal["fractional_positions"], pbc=True)
                centers, neighbors = neighbor_list("ij", atoms, cutoff)
                for index, atom in enumerate(atoms):
                    partner_ids = neighbors[centers == index]
                    self.assertEqual(len(partner_ids), coordination)
                    self.assertTrue(all(atoms[j].symbol != atom.symbol for j in partner_ids))

    def test_rocksalt_uses_room_temperature_fit_at_chosen_pressure(self):
        volume = self.read_crystal("rocksalt")["volume_angstrom3_per_formula_unit"]
        eos = self.phases["rocksalt"]["eos"]
        self.assertAlmostEqual(birch_murnaghan_pressure(volume, eos), 30.0, places=8)
        self.assertAlmostEqual((volume * 4) ** (1 / 3), 4.7288074687, places=9)
        self.assertEqual(self.phases["rocksalt"]["evidence"],
                         "reconstruction_from_experimental_eos_and_b1_symmetry")
        self.assertIsNone(self.phases["rocksalt"]["source_cif"])

    def test_reject_rocksalt_outside_fit_window(self):
        phase = copy.deepcopy(self.phases["rocksalt"])
        phase["conditions"]["pressure_GPa"] = 96
        with self.assertRaisesRegex(ValueError, "outside"):
            make_crystal(phase)

    def test_reject_partial_occupancy(self):
        phase = copy.deepcopy(self.phases["wurtzite"])
        phase["sites"][1][3] = 0.75
        with self.assertRaisesRegex(ValueError, "ordered ZnS"):
            make_crystal(phase)

    def test_prefer_neutron_refinement_over_wyckoff_compilation(self):
        phase = self.phases["wurtzite"]
        self.assertEqual(phase["source_cif"], "1100044")
        self.assertEqual(phase["conditions"]["temperature_K"], 295)
        self.assertEqual(phase["sites"][1][2][2], 0.3748)
        self.assertEqual(self.manifest["originals"]["9008878"]["sulfur_fractional_z"], 0.385)

    def test_reject_source_byte_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "package"
            shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns("__pycache__"))
            path = root / "originals" / "9000107.cif"
            path.write_bytes(path.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "changed source bytes"):
                replay(root)

    def test_reject_cmcm_coordinate_promotion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "package"
            shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns("__pycache__"))
            cmcm = self.phases["cmcm"]
            self.assertEqual(cmcm["cell_lengths_angstrom"], [4.438, 4.486, 4.412])
            self.assertEqual(cmcm["conditions"]["pressure_GPa"], 96)
            (root / "crystals" / "cmcm.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unexpected coordinate file"):
                replay(root)

    def test_reject_stale_export(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "package"
            shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns("__pycache__"))
            (root / "crystals" / "rocksalt.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "stale crystal"):
                replay(root)


if __name__ == "__main__":
    unittest.main()
