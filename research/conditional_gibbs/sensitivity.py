"""Reproduce retrospective resolution budgets and a published EOS comparison."""
from fractions import Fraction
import json

from .lipschitz_gibbs import (
    EV_FORMULA_TO_KJ_MOL,
    GPA_A3_TO_EV_EXACT,
    Phase,
    Sample,
    compare,
    envelope_minimum,
)

DEMUTH = "https://doi.org/10.1088/0953-8984/11/19/306"
ZHONG = "https://doi.org/10.1038/s41524-026-02056-3"

# Published Table 1, in its reported kJ/mol basis. Each tuple is A and B
# in order, followed by their respective uncorrected and reference enthalpies.
RETROSPECTIVE = (
    ("coesite", "tridymite", "2.60", "2.89", "3.68", "3.62"),
    ("clinoenstatite", "enstatite", "-72.04", "-72.57", "-66.62", "-65.73"),
)

# Demuth et al. Tables 8 (C2/c coesite) and 6 (C2221 tridymite).
# E0 (eV/SiO2), V0 (Å³/SiO2), B0 (GPa), and dimensionless B0'.
MURNAGHAN = {
    "LDA": (("-25.9522", "33.43", "94.3", "4.8"),
            ("-25.9302", "45.96", "30.8", "4.1")),
    "GGA": (("-23.7318", "35.28", "91.2", "4.5"),
            ("-23.8485", "48.45", "25.4", "2.3")),
}


def grid_penalty_per_a3(lipschitz_ev_per_a3, pressure_gpa):
    """Worst lower-envelope undershoot per Å³ of maximum knot spacing."""
    L = Fraction(lipschitz_ev_per_a3)
    p = abs(Fraction(pressure_gpa)) * GPA_A3_TO_EV_EXACT
    if L < 0:
        raise ValueError("L must be nonnegative")
    if p >= L:
        return Fraction(0)
    return (L * L - p * p) / (2 * L)


def retrospective():
    for a, b, md_a, md_b, ref_a, ref_b in RETROSPECTIVE:
        raw = Fraction(md_a) - Fraction(md_b)
        reference = Fraction(ref_a) - Fraction(ref_b)
        # Illustrative 50/50 allocation of the absolute reference gap
        # between two equal errors and the runner-up's volume-grid penalty.
        yield {
            "phases_A_minus_B": [a, b],
            "raw_gap_kj_mol_as_reported": float(raw),
            "reference_gap_kj_mol_as_reported": float(reference),
            "per_phase_raw_shift_to_tie_kj_mol_as_reported": float(abs(raw) / 2),
            "per_phase_error_ceiling_no_grid_penalty_kj_mol_as_reported":
                float(abs(reference) / 2),
            "per_phase_error_ceiling_half_gap_for_grid_kj_mol_as_reported":
                float(abs(reference) / 4),
        }


def demuth_phase(name, functional, row):
    e0, v0, bulk, exponent = map(Fraction, row)

    # V/V0 = t**10 makes the published decimal exponents exact rationals
    # at every knot, including endpoints and the fitted equilibrium volume.
    def ratio_power(t, power):
        integer = 10 * power
        if integer.denominator != 1:
            raise ValueError("strain parameterization does not match exponent")
        return t ** int(integer)

    samples, endpoint_pressures = [], []
    for i in range(-20, 21):
        t = 1 + Fraction(i, 10000)
        ratio = t**10
        pressure = bulk / exponent * (ratio_power(t, -exponent) - 1)
        # E(V) - E0 is the integral of -P(V) from V0 to V.
        energy = e0 + GPA_A3_TO_EV_EXACT * bulk * v0 / exponent * (
            ratio - 1 - (ratio_power(t, 1 - exponent) - 1) / (1 - exponent)
        )
        samples.append(Sample(v0 * ratio, energy, 0))
        if abs(i) == 20:
            endpoint_pressures.append(pressure)
    L = GPA_A3_TO_EV_EXACT * max(map(abs, endpoint_pressures))
    return Phase(
        name, (("Si", 1), ("O", 2)), tuple(samples), L,
        f"Demuth et al. within-functional cohesive energy, {functional}",
        f"{functional} nominal Murnaghan static EOS", DEMUTH,
        "zero arithmetic input error for the nominal equation only",
        "analytic monotone Murnaghan pressure bounded by the endpoints",
        0, "reconstructed",
    )


def source_comparison(functional, pressure_gpa):
    phases = [
        demuth_phase(name, functional, row)
        for name, row in zip(("coesite C2/c", "tridymite C2221"),
                             MURNAGHAN[functional])
    ]
    result = compare(phases, pressure_gpa)
    p = Fraction(pressure_gpa) * GPA_A3_TO_EV_EXACT
    low_a, high_a, low_b, high_b = [
        envelope_minimum(phase, p, lower)[0]
        for phase in phases for lower in (True, False)
    ]
    difference = [low_a - high_b, high_a - low_b]
    # Equal uniform residual offsets expand each phase's interval by its
    # offset magnitude. A zero or negative strict margin admits overlap.
    margin = max(difference[0], -difference[1], 0)
    return {
        "target": "nominal published static EOS only",
        "source": DEMUTH,
        "functional": functional,
        "pressure_gpa": float(pressure_gpa),
        "difference_A_minus_B_mev_per_sio2": [float(x * 1000) for x in difference],
        "equal_per_phase_uniform_residual_ceiling_mev_per_sio2": float(margin * 500),
        "max_abs_model_pressure_gpa_by_phase": [
            float(phase.lipschitz_ev_per_a3 / GPA_A3_TO_EV_EXACT) for phase in phases
        ],
        "result": result,
    }


def main():
    benchmark = list(retrospective())
    silica_gap = abs(Fraction(RETROSPECTIVE[0][-2]) - Fraction(RETROSPECTIVE[0][-1]))
    silica_half_budget = silica_gap / (2 * EV_FORMULA_TO_KJ_MOL)
    print(json.dumps({
        "retrospective_source": ZHONG,
        "retrospective_enthalpy_only_300_k_1_bar": benchmark,
        "silica_example_grid_budget": {
            "assumed_internal_pressure_cap_gpa": 1,
            "external_pressure_gpa": 0,
            "maximum_spacing_must_be_less_than_a3_per_sio2": float(
                silica_half_budget / grid_penalty_per_a3(GPA_A3_TO_EV_EXACT, 0)
            ),
        },
        "source_eos_comparisons": [
            source_comparison(functional, pressure)
            for functional in MURNAGHAN for pressure in (0, Fraction(1, 4))
        ],
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
