"""Strict Phonopy FORCE_CONSTANTS parser and exact 3x3x3 reference checks."""
import io
import unittest
import numpy as np
from geoform_phonons.core import Crystal, PhononInputError
from geoform_phonons.fixtures import simple_cubic
from geoform_phonons.phonopy_bridge import (
    parse_force_constants_text, mapped_force_constants, force_constants_text,
)


def supercell_fixture():
    """A fully explicit 3x3x3 simple-cubic supercell has unique +/- neighbors."""
    c = Crystal(np.eye(3)*4, [[0,0,0]], [12], ('X',))
    images = [[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)]
    origin = images.index([0,0,0])
    mapping = {'supercell_images':[{'primitive_index':0,'translation':r} for r in images],
               'origin_rows':[origin]}
    compact = np.zeros((1,27,3,3))
    for s,r in enumerate(images):
        if r==[0,0,0]: compact[0,s] = np.eye(3)*6*4
        elif sum(abs(x) for x in r)==1: compact[0,s] = -4*np.eye(3)
    return c, images, origin, mapping, compact


class PhonopyBridgeTests(unittest.TestCase):
    def test_compact_roundtrip_and_analytic_reference(self):
        crystal,images,origin,mapping,compact = supercell_fixture()
        loaded = parse_force_constants_text(io.StringIO(force_constants_text(compact)))
        np.testing.assert_array_equal(compact,loaded)
        model = mapped_force_constants(crystal,loaded,mapping)
        self.assertLess(model.acoustic_sum_rule_error(),1e-12)
        ref = simple_cubic(mass_amu=12,stiffness_ev_a2=4,a_angstrom=4)
        for q in [[0,0,0],[0.1,0.2,0.3],[0.5,0,0],[0.5,0.5,0.5]]:
            np.testing.assert_allclose(model.modes(q)[0], ref.modes(q)[0],atol=1e-9)
            np.testing.assert_allclose(model.dynamical_matrix(q),ref.dynamical_matrix(q),atol=1e-12)

    def test_full_supercell_selects_explicit_origin(self):
        crystal,images,origin,mapping,compact = supercell_fixture()
        full = np.zeros((27,27,3,3))
        for i,ri in enumerate(images):
            for j,rj in enumerate(images):
                relative = [(rj[k]-ri[k])%3 for k in range(3)]
                if relative==[0,0,0]: full[i,j]=np.eye(3)*24
                elif sum(x in (1,2) for x in relative)==1 and sum(x!=0 for x in relative)==1:
                    full[i,j]=-4*np.eye(3)
        model=mapped_force_constants(crystal,full,mapping)
        reference=mapped_force_constants(crystal,compact,mapping)
        for q in [[0.1,0.2,0.3],[0.5,0.5,0.5]]:
            np.testing.assert_allclose(model.modes(q)[0], reference.modes(q)[0],atol=1e-9)

    def test_missing_map_refused(self):
        crystal,images,origin,mapping,compact=supercell_fixture()
        with self.assertRaises(PhononInputError): mapped_force_constants(crystal,compact,{})
        mapping['supercell_images'][0]['translation']=[0,0,0]
        with self.assertRaises(PhononInputError): mapped_force_constants(crystal,compact,mapping)

    def test_unpaired_ifcs_refused(self):
        crystal,images,origin,mapping,compact=supercell_fixture()
        compact[0,0] = np.eye(3)*0.5
        with self.assertRaises(PhononInputError): mapped_force_constants(crystal,compact,mapping)

    def test_bad_shape_and_extra_text(self):
        crystal,images,origin,mapping,compact=supercell_fixture()
        with self.assertRaises(PhononInputError): parse_force_constants_text(io.StringIO('2\n1 1\n1 0 0\n0 1 0\n0 0 1\n'))
        with self.assertRaises(PhononInputError): parse_force_constants_text(io.StringIO(force_constants_text(compact)+'EXTRA\n'))
        with self.assertRaises(PhononInputError): mapped_force_constants(crystal,np.zeros((2,4,3,3)),mapping)


if __name__=='__main__': unittest.main()
