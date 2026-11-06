"""Offline tests for short-vector multiplicity and strict Phonopy import."""
import unittest
from types import SimpleNamespace
import numpy as np
from geoform_phonons import PhononInputError
from geoform_phonons.core import THZ_FACTOR
from geoform_phonons.phonopy_native import from_phonopy


def boundary_pair_fixture(full=False):
    # Single atom, x-bonds cross the boundary of a 2x1x1 supercell.
    # The one second-supercell atom represents BOTH primitive +/- x images.
    k = 4.0
    cell = np.eye(3)*3
    primitive = SimpleNamespace(
        cell=cell,scaled_positions=np.zeros((1,3)),masses=np.array([20.0]),
        symbols=['Mg'],store_dense_svecs=True,p2s_map=np.array([0]),
        s2p_map=np.array([0,0]),p2p_map={0:0},
    )
    class P:
        def __len__(self):return 1
        def __getattr__(self,key):return getattr(primitive,key)
        def get_smallest_vectors(self):
            return np.array([[0,0,0],[-1,0,0],[1,0,0]],float), np.array([[[1,0]],[[2,1]]])
    f=np.zeros((2,2,3,3))
    f[0,0,0,0]=2*k
    f[0,1,0,0]=-2*k
    f[1,0,0,0]=-2*k
    f[1,1,0,0]=2*k
    p = SimpleNamespace(primitive=P(),supercell=[0,1],force_constants=f if full else f[:1])
    return p


class NativeBridgeTests(unittest.TestCase):
    def test_split_equidistant_boundary_images(self):
        for full in (False, True):
            model=from_phonopy(boundary_pair_fixture(full))
            self.assertEqual(set(model.blocks),{(0,0,0,0,0),(-1,0,0,0,0),(1,0,0,0,0)})
            self.assertAlmostEqual(model.blocks[(1,0,0,0,0)][0,0],-4)
            self.assertAlmostEqual(model.acoustic_sum_rule_error(),0)
            for qx in (0,0.125,0.25,0.5,0.73):
                actual=model.dynamical_matrix([qx,0,0])[0,0].real
                expect=(8-8*np.cos(2*np.pi*qx))/20
                self.assertAlmostEqual(actual,expect,places=11)
            expected=np.sqrt(16/20)*THZ_FACTOR
            self.assertAlmostEqual(model.modes([0.5,0,0])[0][-1],expected)

    def test_reject_incomplete_or_malformed_data(self):
        p=boundary_pair_fixture()
        p.force_constants[0,0,0,0]=7
        self.assertGreater(from_phonopy(p).acoustic_sum_rule_error(), 0)
        p=boundary_pair_fixture()
        p.force_constants = np.zeros((1,2,3,2))
        with self.assertRaises(PhononInputError):from_phonopy(p)
        p=boundary_pair_fixture()
        p.primitive.__class__.get_smallest_vectors=lambda self:(np.array([[0,0,0],[.1,0,0],[1,0,0]]),np.array([[[1,0]],[[2,1]]]))
        with self.assertRaises(PhononInputError):from_phonopy(p)

if __name__=='__main__':unittest.main()
