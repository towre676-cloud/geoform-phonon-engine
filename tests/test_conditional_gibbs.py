"""Analytic and falsification tests for conditional Gibbs bounds."""
import math
import random
import unittest
from dataclasses import replace
from fractions import Fraction
from research.conditional_gibbs.lipschitz_gibbs import (
    Phase, Sample, minimum_bounds, compare, GPA_A3_TO_EV,
    GPA_A3_TO_EV_EXACT, EV_FORMULA_TO_KJ_MOL, envelope_minimum,
)

from research.conditional_gibbs.sensitivity import (
    MURNAGHAN, grid_penalty_per_a3, retrospective, source_comparison,
)

def phase(name, curve, error=0.01, L=1, count=8):
    return Phase(
        name=name, composition=(("X", 1),),
        samples=tuple(Sample(1+i/(count-1), curve(1+i/(count-1)), error)
                      for i in range(count)),
        lipschitz_ev_per_a3=L, reference="illustrative ref",
        method="illustrative method", source="manufactured analytic function",
        error_source="manufactured error bars, not inferred from data",
        lipschitz_source="analytic derivative bound, illustrative",
        temperature_k=300,
    )

class TestConditionalGibbs(unittest.TestCase):
    def test_constant_exact_under_pressure(self):
        p=phase("a", lambda v:2, error=0, L=0)
        for pressure in (0, 1, 10):
            r=minimum_bounds(p,pressure,allow_illustrative=True)
            exact=2+pressure*GPA_A3_TO_EV
            self.assertAlmostEqual(r["interval_ev_formula"][0],exact,places=12)
            self.assertAlmostEqual(r["interval_ev_formula"][1],exact,places=12)

    def test_linear_exact_at_zero_pressure(self):
        r=minimum_bounds(phase("a",lambda v:3*v+1,error=0,L=3,count=9),0,
                         allow_illustrative=True)
        self.assertAlmostEqual(r["interval_ev_formula"][0],4)
        self.assertAlmostEqual(r["interval_ev_formula"][1],4)

    def test_200_independent_quadratic_minima(self):
        randomizer=random.Random(1923)
        for _ in range(200):
            a=randomizer.uniform(.2,3)
            v0=randomizer.uniform(1,2)
            bias=randomizer.uniform(-2,2)
            error=randomizer.uniform(0,.1)
            p=phase("a",lambda v:a*(v-v0)**2+bias,error=error,L=2*a)
            pressure=randomizer.uniform(0,25)
            r=minimum_bounds(p,pressure,allow_illustrative=True)
            pv=pressure*GPA_A3_TO_EV
            vstar=min(2,max(1,v0-pv/(2*a)))
            truth=a*(vstar-v0)**2+bias+pv*vstar
            lo,hi=r["interval_ev_formula"]
            self.assertLessEqual(lo,truth+1e-10)
            self.assertGreaterEqual(hi,truth-1e-10)

    def test_wrong_lipschitz_constant_rejected_when_contradictory(self):
        with self.assertRaises(ValueError):
            phase("a",lambda v:v,error=0,L=0).validate()

    def test_no_automatic_evidence_upgrade(self):
        with self.assertRaises(ValueError):
            minimum_bounds(phase("a",lambda v:0,L=0),0)

    def test_true_overlap_abstains(self):
        a=phase("a",lambda v:0,error=.02,L=0)
        b=phase("b",lambda v:.01,error=.02,L=0)
        self.assertIsNone(compare([a,b],0,allow_illustrative=True)["winner"])

    def test_distinct_intervals_choose_one(self):
        a=phase("a",lambda v:0,error=.02,L=0)
        b=phase("b",lambda v:.20,error=.02,L=0)
        self.assertEqual(compare([a,b],0,allow_illustrative=True)["winner"],"a")

    def test_references_cannot_be_mixed(self):
        a=phase("a",lambda v:0,L=0)
        b=replace(phase("b",lambda v:.2,L=0),reference="other")
        with self.assertRaises(ValueError):
            compare([a,b],0,allow_illustrative=True)

    def test_invalid_domains_and_provenance(self):
        p=phase("a",lambda v:0,L=0)
        for pressure in (float("nan"), float("inf"), -1):
            with self.assertRaises(ValueError):
                minimum_bounds(p,pressure,allow_illustrative=True)
        with self.assertRaises(ValueError):
            replace(p,lipschitz_source="").validate()
        with self.assertRaises(ValueError):
            replace(p,samples=p.samples[::-1]).validate()

    def test_tiny_contradiction_is_rejected(self):
        p = phase("a", lambda v: 0, error=0, L=0, count=2)
        p = replace(p, samples=(Sample(1, 0, 0), Sample(2, 5e-13, 0)))
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            minimum_bounds(p, 0, allow_illustrative=True)

    def test_nonadjacent_contradiction_is_rejected(self):
        p = phase("a", lambda v: 0)
        p = replace(p, samples=(Sample(1, 0, 0), Sample(2, 1, 10), Sample(3, 3, 0)))
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            p.validate()

    def test_nearly_parallel_cones_can_change_the_winner(self):
        L = Fraction(1, 10**16)
        a = phase("a", lambda v: 0, error=0, L=L, count=2)
        b = phase("b", lambda v: L / 4, error=0, L=0, count=2)
        result = compare([a, b], 0, allow_illustrative=True)
        self.assertEqual(result["winner"], "a")
        self.assertEqual(envelope_minimum(a, 0, True), (-L / 2, Fraction(3, 2)))
        # The exact lower envelope of B must also be allowed to dip. If
        # B's knot values alone were used, A would be declared lower.
        a, b = replace(a, lipschitz_ev_per_a3=0), replace(b, lipschitz_ev_per_a3=L)
        self.assertIsNone(compare([a, b], 0, allow_illustrative=True)["winner"])

    def test_large_volume_translation_preserves_zero_pressure_interval(self):
        p = phase("a", lambda v: 0, error=0, L=10**10, count=2)
        samples = (Sample(10**16, 0, 0), Sample(10**16 + 1, 0, 0))
        result = minimum_bounds(replace(p, samples=samples), 0, allow_illustrative=True)
        self.assertEqual(result["interval_ev_formula"], [-5 * 10**9, 0])

    def test_pressure_conversion_and_outward_rounding(self):
        p = phase("a", lambda v: 0, error=0, L=0, count=2)
        low, high = minimum_bounds(p, 1, allow_illustrative=True)["interval_ev_formula"]
        self.assertLess(Fraction(low), GPA_A3_TO_EV_EXACT)
        self.assertGreater(Fraction(high), GPA_A3_TO_EV_EXACT)
        self.assertEqual(EV_FORMULA_TO_KJ_MOL,
                         Fraction("1.602176634e-19") * Fraction("6.02214076e23") / 1000)

    def test_grid_penalty_is_sharp_for_equal_tilted_knot_values(self):
        slope = GPA_A3_TO_EV_EXACT
        p = phase("a", lambda v: 0, error=0, L=2 * slope, count=2)
        p = replace(p, samples=(Sample(1, -slope, 0), Sample(2, -2 * slope, 0)))
        minimum, _ = envelope_minimum(p, slope, True)
        self.assertEqual(minimum, -grid_penalty_per_a3(2 * slope, 1))
        self.assertEqual(grid_penalty_per_a3(slope, 1), 0)
        self.assertEqual(grid_penalty_per_a3(0, 0), 0)

    def test_reconstructed_eos_exposes_opposite_functional_orderings(self):
        for functional, exact_gap, winner in (
            ("LDA", -22, "coesite C2/c"), ("GGA", 116.7, "tridymite C2221")
        ):
            result = source_comparison(functional, 0)
            low, high = result["difference_A_minus_B_mev_per_sio2"]
            self.assertLessEqual(low, exact_gap)
            self.assertGreaterEqual(high, exact_gap)
            self.assertEqual(result["result"]["winner"], winner)
            self.assertEqual(result["target"], "nominal published static EOS only")

    def test_eos_nonzero_pressure_against_stationary_solution(self):
        pressure = 0.25
        for functional, rows in MURNAGHAN.items():
            minima = []
            for row in rows:
                e0, v0, bulk, exponent = map(float, row)
                # Solve P_internal(V) = P_external independently of the grid.
                ratio = (1 + pressure * exponent / bulk) ** (-1 / exponent)
                self.assertGreater(ratio, 0.998**10)
                self.assertLess(ratio, 1.002**10)
                delta = bulk * v0 / exponent * (
                    ratio - 1 - (ratio ** (1 - exponent) - 1) / (1 - exponent)
                )
                minima.append(e0 + GPA_A3_TO_EV * (delta + pressure * v0 * ratio))
            difference = 1000 * (minima[0] - minima[1])
            lo, hi = source_comparison(functional, pressure)["difference_A_minus_B_mev_per_sio2"]
            self.assertLessEqual(lo, difference)
            self.assertGreaterEqual(hi, difference)

    def test_retrospective_tie_shifts_are_separate_from_precision_budgets(self):
        silica, enstatite = retrospective()
        for result, shift, budget in ((silica, 0.145, 0.030), (enstatite, 0.265, 0.445)):
            self.assertEqual(result["per_phase_raw_shift_to_tie_kj_mol_as_reported"], shift)
            self.assertEqual(result["per_phase_error_ceiling_no_grid_penalty_kj_mol_as_reported"],
                             budget)

if __name__=="__main__":
    unittest.main()
