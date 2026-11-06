import unittest
from math import pi, sqrt
from phonon_chain import branches, sound_speed, sample


class DiatomicChainTests(unittest.TestCase):
    def test_gamma_acoustic_zero(self):
        acoustic, optical = branches(0.0, 1.0, 1.0, 2.0, 3.0)
        self.assertEqual(acoustic, 0.0)
        self.assertAlmostEqual(optical, sqrt(9.0))

    def test_zone_boundary(self):
        a, o = branches(pi, 1.0, 1.0, 2.0, 3.0)
        self.assertAlmostEqual(a, sqrt(3.0))
        self.assertAlmostEqual(o, sqrt(6.0))

    def test_mass_exchange_and_reciprocity(self):
        for k in (0, 0.01, pi/5, pi):
            self.assertEqual(branches(k, 1, 1, 2, 3),
                             branches(-k, 1, 1, 2, 3))
            self.assertEqual(branches(k, 1, 1, 2, 3),
                             branches(k, 1, 2, 1, 3))

    def test_long_wave_sound_speed(self):
        k = 1e-5
        a, _ = branches(k, 2.0, 1.5, 3.5, 4.0)
        self.assertAlmostEqual(a/k, sound_speed(2.0, 1.5, 3.5, 4.0), places=7)

    def test_ordering_and_domain(self):
        for _, a, o in sample(1, 1, 3, 2, 51):
            self.assertGreaterEqual(o, a)
            self.assertGreaterEqual(a, 0)
        for params in ((0, 1, 0, 1, 1), (0, 1, 1, -1, 1),
                       (float("nan"), 1, 1, 1, 1)):
            with self.assertRaises(ValueError):
                branches(*params)


if __name__ == "__main__":
    unittest.main()
