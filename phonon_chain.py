"""Exact-form dispersion for a one-dimensional harmonic diatomic chain.

A minimal, independently testable teaching/reference model. This is NOT a
3-D crystal phonon solver, an ab-initio model, or experimental mineral data.
Units: spring stiffness N/m; masses kg; unit-cell length m; wavevector 1/m;
output angular frequencies rad/s (not Hz or cm^-1).
"""
from math import sin, sqrt, pi, isfinite
import argparse


def branches(k, cell_length, mass_a, mass_b, spring):
    """Return (acoustic, optical) angular frequencies for periodic AB chain."""
    for name, x in (("k", k), ("cell_length", cell_length),
                    ("mass_a", mass_a), ("mass_b", mass_b), ("spring", spring)):
        if not isfinite(x):
            raise ValueError(name + " must be finite")
    if min(cell_length, mass_a, mass_b, spring) <= 0:
        raise ValueError("length, masses, and stiffness must be positive")
    trace = 2.0 * spring * (1.0 / mass_a + 1.0 / mass_b)
    phase = sin(0.5 * k * cell_length)
    # Product of eigenvalues of the mass-weighted dynamical matrix.
    determinant = 4.0 * spring * spring * phase * phase / (mass_a * mass_b)
    discriminant = max(0.0, trace * trace - 4.0 * determinant)
    optical_sq = 0.5 * (trace + sqrt(discriminant))
    # This avoids catastrophic subtraction near k=0.
    acoustic_sq = determinant / optical_sq
    return sqrt(acoustic_sq), sqrt(optical_sq)


def sound_speed(cell_length, mass_a, mass_b, spring):
    """Long-wave slope d omega_acoustic / d k."""
    if min(cell_length, mass_a, mass_b, spring) <= 0:
        raise ValueError("length, masses, and stiffness must be positive")
    return cell_length * sqrt(spring / (2.0 * (mass_a + mass_b)))


def sample(cell_length, mass_a, mass_b, spring, points=101):
    if not isinstance(points, int) or points < 2:
        raise ValueError("points must be an integer >= 2")
    return [(-pi / cell_length + 2 * pi * i / (cell_length * (points - 1)),
             *branches(-pi / cell_length + 2 * pi * i / (cell_length * (points - 1)),
                       cell_length, mass_a, mass_b, spring))
            for i in range(points)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mass-a", type=float, default=1.0)
    parser.add_argument("--mass-b", type=float, default=2.0)
    parser.add_argument("--spring", type=float, default=1.0)
    parser.add_argument("--cell-length", type=float, default=1.0)
    args = parser.parse_args()
    print("k (1/m), acoustic omega (rad/s), optical omega (rad/s)")
    for record in sample(args.cell_length, args.mass_a, args.mass_b, args.spring, 9):
        print(",".join(f"{v:.12g}" for v in record))
