"""Analytic completeness tests for atom-PDOS and thermal curves."""
import unittest
import numpy as np
from geoform_phonons import (
    simple_cubic, projected_density_of_states, density_of_states,
    mode_participation_ratios, thermodynamic_curve, frequency_grid,
    PhononInputError,
)
from geoform_phonons.fixtures import from_springs
from geoform_phonons.core import Crystal

class ProjectedDensityTests(unittest.TestCase):
    def test_total_and_atom_projected_dos_are_pointwise_identical(self):
        model=simple_cubic()
        xx, total, pdos=projected_density_of_states(model, mesh=(5,4,3),sigma_thz=.17,points=350)
        self.assertEqual(pdos.shape,(1,len(xx)))
        np.testing.assert_allclose(total,pdos.sum(axis=0),atol=1e-10)
        self.assertAlmostEqual(np.trapezoid(pdos[0],xx),3.0,places=3)
        xx2,total2=density_of_states(frequency_grid(model,(5,4,3)),.17,350)
        np.testing.assert_array_equal(xx,xx2)
        np.testing.assert_allclose(total,total2,atol=1e-10)

    def test_binary_projected_dos_atom_integrals(self):
        c=Crystal(np.eye(3)*4,[[0,0,0],[0.5,0.5,0.5]],[10,20],('A','B'))
        model=from_springs(c,[
            (0,1,(r,s,t),np.eye(3)*3) for r,s,t in [(0,0,0),(-1,0,0),(0,-1,0),(0,0,-1)]
        ])
        xx,total,atom=projected_density_of_states(model,(3,4,5),.2,200)
        np.testing.assert_allclose(atom.sum(axis=0),total,atol=2e-10)
        for row in atom:self.assertAlmostEqual(np.trapezoid(row,xx),3,places=2)

    def test_participation_ratios_and_orthonormality(self):
        model=simple_cubic()
        freq,e=model.modes((.2,.1,.3))
        np.testing.assert_allclose(mode_participation_ratios(e,1),np.ones(3))
        e=np.eye(6)
        np.testing.assert_allclose(mode_participation_ratios(e,2),np.ones(6)*.5)
        with self.assertRaises(PhononInputError):mode_participation_ratios(np.ones((6,6)),2)

    def test_thermodynamic_curve(self):
        freq=frequency_grid(simple_cubic(),(5,5,5))
        temperatures=[0,100,300,1500]
        properties=thermodynamic_curve(freq,temperatures)
        self.assertEqual([p['temperature_k'] for p in properties],temperatures)
        self.assertTrue(all(p['cv_ev_k_cell']>=0 for p in properties))
        self.assertLess(properties[0]['free_energy_ev_cell'],properties[0]['internal_energy_ev_cell']+1e-9)
        self.assertLess(properties[3]['free_energy_ev_cell'],properties[2]['free_energy_ev_cell'])
        with self.assertRaises(PhononInputError):thermodynamic_curve(freq,[300,100])
if __name__=='__main__':unittest.main()
