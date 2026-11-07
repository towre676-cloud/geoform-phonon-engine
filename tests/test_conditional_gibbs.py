"""Analytic and falsification tests for conditional Gibbs bounds."""
import math
import random
import unittest
from dataclasses import replace
from research.conditional_gibbs.lipschitz_gibbs import (
    Phase, Sample, minimum_bounds, compare, GPA_A3_TO_EV,
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
        r=minimum_bounds(phase("a",lambda v:3*v+1,error=0,L=3),0,
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

if __name__=="__main__":
    unittest.main()
