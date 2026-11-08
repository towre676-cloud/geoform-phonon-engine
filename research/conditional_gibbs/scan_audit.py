"""Audit the authors' r²SCAN volume scans under explicit hypothetical inputs."""
import argparse
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from hashlib import file_digest
from itertools import pairwise
import json
from math import ceil
from pathlib import Path
import re
import sqlite3
from struct import unpack
from zipfile import ZipFile

from .lipschitz_gibbs import (
    GPA_A3_TO_EV_EXACT,
    Phase,
    Sample,
    compare,
)

from .sensitivity import grid_penalty_per_a3

SOURCE = "https://doi.org/10.5281/zenodo.14925354"
SHA256 = "d974a05835951fd4e27c83bebcacb3280d7c0d6fd82e01a94c94933762eac1e4"
BRANCHES = {"coesite": 8, "alpha-tridymite": 4, "beta-tridymite": 4, "MC-tridymite": 4}
ADDITIONAL_BRANCHES = {
    "alpha-quartz-1": 3, "alpha-quartz-2": 3,
    "beta-quartz-1": 3, "beta-quartz-2": 3,
    "alpha-cristobalite": 4, "keatite": 12, "stishovite": 2,
}


def require_digest(path, expected):
    with path.open("rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != expected:
            raise ValueError(f"{path.name} SHA-256 mismatch")


def scan_rows(database, include_all=False):
    require_digest(database, SHA256)
    connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        for name, formula_units in (BRANCHES | ADDITIONAL_BRANCHES if include_all else BRANCHES).items():
            pattern = re.compile(rf"run_SiO2-{name}-volsc[0-9]+\.[0-9]+\.xml\.gz")
            query = """
                SELECT s.id, t.value, s.natoms, s.volume, s.energy, s.stress
                FROM text_key_values t JOIN systems s ON t.id = s.id
                WHERE t.key = 'filename' AND t.value LIKE ?
                ORDER BY s.volume
            """
            rows = []
            for row_id, filename, atoms, volume, energy, stress in connection.execute(
                query, (f"run_SiO2-{name}-volsc%",)
            ):
                if pattern.fullmatch(filename) is None:
                    continue
                if atoms != 3 * formula_units:
                    raise ValueError(f"unexpected cell size at row {row_id}")
                species = dict(connection.execute("SELECT Z, n FROM species WHERE id = ?",
                                                  (row_id,)))
                if species != {8: 2 * formula_units, 14: formula_units}:
                    raise ValueError(f"unexpected composition at row {row_id}")
                # ASE stress is tension-positive, in eV/Å³, in Voigt order.
                pressure = -sum(map(Fraction, unpack("<6d", stress)[:3])) / (
                    3 * GPA_A3_TO_EV_EXACT
                )
                rows.append((row_id, Fraction(volume) / formula_units,
                             Fraction(energy) / formula_units, pressure))
            expected_count = 22 if name == "keatite" else 26
            if len(rows) != expected_count:
                raise ValueError(f"expected {expected_count} undisturbed scans for {name}")
            yield name, rows
    finally:
        connection.close()


def audit(database, pressures=(0, Fraction(1, 4)), *, include_all=False):
    phases, diagnostics = [], []
    for name, rows in scan_rows(database, include_all):
        samples = tuple(Sample(volume, energy, 0) for _, volume, energy, _ in rows)
        secants = [(b.energy - a.energy) / (b.volume - a.volume)
                   for a, b in pairwise(samples)]
        stress_peak = max(abs(pressure) for _, _, _, pressure in rows)
        secant_peak = max(map(abs, secants)) / GPA_A3_TO_EV_EXACT
        cap_gpa = ceil(max(stress_peak, secant_peak))
        phases.append(Phase(
            name, (("Si", 1), ("O", 2)), samples, cap_gpa * GPA_A3_TO_EV_EXACT,
            "Svensson et al. r2SCAN ASE energy field", "static r2SCAN", SOURCE,
            "hypothetical absolute sample error; no physical error asserted",
            "hypothetical global pressure cap rounded up from sampled diagnostics",
            0, "illustrative",
        ))
        diagnostics.append({
            "branch": name,
            "rows": [row[0] for row in rows],
            "sampled_volume_domain_a3_per_sio2": [float(samples[0].volume),
                                                   float(samples[-1].volume)],
            "max_spacing_a3_per_sio2": float(max(b.volume - a.volume
                                                for a, b in pairwise(samples))),
            "max_abs_sampled_pressure_gpa": float(stress_peak),
            "zero_error_secant_floor_for_global_cap_gpa": float(secant_peak),
            "decreasing_adjacent_secants": sum(a > b for a, b in pairwise(secants)),
            "hypothetical_global_pressure_cap_gpa": cap_gpa,
            "min_sample_energy_ev_per_sio2": float(min(s.energy for s in samples)),
        })

    scenarios = []
    for error_mev in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2)):
        bounded = [replace(phase, samples=tuple(replace(s, error=error_mev / 1000)
                                                for s in phase.samples))
                   for phase in phases]
        for pressure in pressures:
            result = compare(bounded, pressure, allow_illustrative=True)
            intervals = {entry["phase"]: entry["interval_ev_formula"]
                         for entry in result["results"]}
            lo, hi = intervals["coesite"]
            comparisons = {name: [(lo - b) * 1000, (hi - a) * 1000]
                           for name, (a, b) in intervals.items() if name != "coesite"}
            local = {name: intervals[name] for name in BRANCHES}
            four_winner = next((name for name, (_, upper) in local.items()
                                if all(upper < lower for other, (lower, _) in local.items()
                                       if other != name)), None)
            scenario = {
                "pressure_gpa": float(pressure),
                "hypothetical_sample_error_mev_per_sio2": float(error_mev),
                "coesite_minus_tridymite_mev_per_sio2_approx": {
                    name: comparisons[name] for name in BRANCHES if name != "coesite"
                },
                "winner_all_four_under_hypotheses": four_winner,
            }
            if include_all:
                scenario.update(
                    coesite_minus_other_mev_per_sio2_approx=comparisons,
                    winner_all_branches_under_hypotheses=result["winner"],
                )
            scenarios.append(scenario)
    # Worked sizing example for the nominal MC–coesite sampled gap at 0 GPa.
    coesite, _, _, mc = phases[:len(BRANCHES)]
    gap = min(s.energy for s in coesite.samples) - min(s.energy for s in mc.samples)
    per_phase_error = Fraction(1, 1000)
    max_spacing = (gap - 2 * per_phase_error) / grid_penalty_per_a3(
        coesite.lipschitz_ev_per_a3, 0
    )
    return {
        "source": SOURCE, "database_sha256": SHA256,
        "scope": "hypothetical bounds on restricted static scans; physical certification withheld",
        "branches": diagnostics, "scenarios": scenarios,
        "grid_example_if_nominal_gap_persists": {
            "pair": ["MC-tridymite", "coesite"], "pressure_gpa": 0,
            "hypothetical_error_mev_per_sio2_each": 1,
            "nominal_sample_gap_mev_per_sio2": float(gap * 1000),
            "coesite_spacing_must_be_less_than_a3_per_sio2": float(max_spacing),
        },
    }


ZHONG_TRAINING_SHA256 = "635a65a7de74f874b3a7696a909f3492fb326123d122d6df02872aa713f54862"


def pyroxene_training_audit(archive):
    """Show why training configurations cannot be used as relaxed F(V) knots."""
    require_digest(archive, ZHONG_TRAINING_SHA256)
    results = []
    with ZipFile(archive) as zipped:
        for phase in ("Clinoenstatite", "Enstatite"):
            for family in ("ref", "MD"):
                prefix = f"{phase}_{family}/"
                types = Counter(zipped.read(prefix + "type.raw").split())
                formulas = types[b"0"]
                if not formulas or types != {b"0": formulas, b"1": 3 * formulas,
                                             b"3": formulas}:
                    raise ValueError(f"unexpected type counts in {prefix}")
                energies = zipped.read(prefix + "energy.raw").decode("ascii").split()
                boxes = zipped.read(prefix + "box.raw").decode("ascii").splitlines()
                samples = []
                for index, (energy, box) in enumerate(zip(energies, boxes, strict=True)):
                    a, b, c, d, e, f, g, h, i = map(Fraction, box.split())
                    volume = abs(a * (e * i - f * h) - b * (d * i - f * g)
                                 + c * (d * h - e * g))
                    samples.append((volume / formulas, Fraction(energy) / formulas, index))
                samples.sort()
                # Near-equal volumes with very different energies are evidence
                # that the frames are not samples of one relaxed branch.
                close_pairs = [(abs(a[1] - b[1]), b[0] - a[0], a[2], b[2])
                               for a, b in pairwise(samples) if b[0] - a[0] < Fraction(1, 10**4)]
                energy_gap, volume_gap, row_a, row_b = max(close_pairs)
                results.append({
                    "family": prefix.rstrip("/"), "frames": len(samples),
                    "formula_units_per_cell": formulas,
                    "witness_frame_indices_zero_based": [row_a, row_b],
                    "witness_volume_gap_a3_per_mgsio3": float(volume_gap),
                    "witness_energy_gap_mev_per_mgsio3": float(energy_gap * 1000),
                    "scope": "training frames, not constrained-relaxed EOS samples",
                })
    return {"source": "https://osf.io/784bv/", "archive_sha256": ZHONG_TRAINING_SHA256,
            "families": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--zhong-training", type=Path)
    parser.add_argument("--pressure-sweep", action="store_true",
                        help="evaluate 0 to 1 GPa in exact 0.05 GPa steps")
    parser.add_argument("--all-silica-branches", action="store_true",
                        help="include all 11 plain volume scans in the pinned database")
    args = parser.parse_args()
    pressures = tuple(Fraction(i, 20) for i in range(21)) if args.pressure_sweep else (0, Fraction(1, 4))
    report = audit(args.database, pressures, include_all=args.all_silica_branches)
    if args.zhong_training:
        report["pyroxene_training"] = pyroxene_training_audit(args.zhong_training)
    print(json.dumps(report, indent=2, allow_nan=False))
