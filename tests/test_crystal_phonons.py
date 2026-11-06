"""Reproducible analytic and metamorphic checks; no real-material claims."""
import json
import math
import tempfile
import unittest
from pathlib import Path
import numpy as np
from geoform_phonons import (
    Crystal, ForceConstants, PhononInputError, simple_cubic, diatomic_chain,
    shifted_q_mesh, frequency_grid, density_of_states, thermal_properties,
    thz_to_wavenumber_cm, ir_mode_activity, raman_mode_activity,
    eliashberg_moments, from_springs,
)
from geoform_phonons.core import THZ_FACTOR, EV, ANGSTROM, AMU
from geoform_phonons.io import parse_model, model_to_record, export_result, save_json, load_json
from phonon_chain import branches


class CrystalTests(unittest.TestCase):
    def setUp(self):
        self.cubic = simple_cubic(mass_amu=4, stiffness_ev_a2=2,a_angstrom=3)

    def test_volume_reciprocal_basis(self):
        c=self.cubic.crystal
        self.assertAlmostEqual(c.volume_angstrom3,27)
        np.testing.assert_allclose(c.lattice_angstrom@c.reciprocal_vectors.T,
                                   2*np.pi*np.eye(3),atol=1e-12)

    def test_invalid_cells(self):
        for cell in [np.zeros((3,3)),[[1,0],[0,1]], [[float('nan')]*3]*3]:
            with self.assertRaises((PhononInputError,ValueError)):
                Crystal(cell,[[0,0,0]],[1],('A',))

    def test_invalid_masses(self):
        for mass in [0,-1,float('inf')]:
            with self.assertRaises(PhononInputError):
                Crystal(np.eye(3),[[0,0,0]],[mass],('A',))

    def test_invalid_ifc_and_missing_pairs(self):
        with self.assertRaises(PhononInputError):
            ForceConstants(self.cubic.crystal,{(1,0,0,0,0):np.eye(3)})
        with self.assertRaises(PhononInputError):
            ForceConstants(self.cubic.crystal,{(0,0,0,0,0):np.array([[1,2,0],[3,1,0],[0,0,1]])})

    def test_asr_cubic(self):
        self.assertLess(self.cubic.acoustic_sum_rule_error(),1e-12)
        np.testing.assert_allclose(self.cubic.modes([0,0,0])[0], [0,0,0],atol=1e-12)

    def test_asr_repair_onsite(self):
        b={k:v.copy() for k,v in self.cubic.blocks.items()}
        b[(0,0,0,0,0)]+=np.eye(3)*0.2
        bad=ForceConstants(self.cubic.crystal,b)
        self.assertAlmostEqual(bad.acoustic_sum_rule_error(),0.2)
        fixed=bad.enforce_asr()
        self.assertLess(fixed.acoustic_sum_rule_error(),1e-12)
        for q in [[0.3,0.2,0.1],[0.5,0,0]]:
            np.testing.assert_allclose(fixed.modes(q)[0],self.cubic.modes(q)[0],atol=1e-12)

    def test_hermitian_complex_q(self):
        for q in [[0.3,0.19,0.4],[0.21,0.35,0.11],[1.1,-0.8,0.7]]:
            d=self.cubic.dynamical_matrix(q)
            np.testing.assert_allclose(d,d.conj().T,atol=1e-12)

    def test_reciprocal_periodicity(self):
        a=self.cubic.modes([0.17,0.27,0.43])[0]
        b=self.cubic.modes([1.17,-1.73,0.43])[0]
        np.testing.assert_allclose(a,b,atol=1e-11)

    def test_cubic_analytic_zone_boundary(self):
        f=self.cubic.modes([0.5,0,0])[0]
        np.testing.assert_allclose(f,np.full(3,2*np.sqrt(2/4)*THZ_FACTOR),rtol=1e-12)

    def test_cubic_masses_scale_inverse_sqrt(self):
        a=self.cubic.modes([0.3,0.2,0.1])[0]
        b=simple_cubic(mass_amu=16,stiffness_ev_a2=2,a_angstrom=3).modes([0.3,0.2,0.1])[0]
        np.testing.assert_allclose(a,2*b,rtol=1e-12)

    def test_diatomic_chain_crosscheck_analytic(self):
        f=diatomic_chain(1.5,2.5,4,6)
        for q in [0,0.02,0.1,0.33,0.5]:
            w0,w1=branches(2*np.pi*q/(6*ANGSTROM),6*ANGSTROM,1.5*AMU,2.5*AMU,4*EV/ANGSTROM**2)
            expected=np.array([w0,w1])/(2*np.pi*1e12)
            actual=f.modes([q,0,0])[0][-2:]
            np.testing.assert_allclose(actual,expected,rtol=1e-9,atol=1e-9)

    def test_chain_gamma_acoustic_and_optic(self):
        model=diatomic_chain()
        self.assertLess(model.acoustic_sum_rule_error(),1e-12)
        self.assertEqual(np.count_nonzero(model.modes([0,0,0])[0]),1)
        self.assertEqual(np.count_nonzero(model.modes([0.5,0,0])[0]),2)

    def test_mass_weighted_eigenvectors(self):
        f,e=self.cubic.modes([0.23,0.1,0.31])
        np.testing.assert_allclose(e.conj().T@e,np.eye(3),atol=1e-12)

    def test_mesh_sample_count_and_shift(self):
        points=shifted_q_mesh((3,4,5))
        self.assertEqual(len(points),60)
        self.assertTrue(all(0 < q < 1 for p in points for q in p))
        for shape in [(0,2,3),(1,2),[1,2,3.1]]:
            with self.assertRaises(PhononInputError): shifted_q_mesh(shape)

    def test_bands_lengths_monotone(self):
        result=self.cubic.band_path([[0,0,0],[0.5,0,0],[0.5,0.5,0]],11)
        self.assertEqual(len(result['frequency_thz_signed']),22)
        self.assertTrue(all(b>=a for a,b in zip(result['distance_inverse_angstrom'],result['distance_inverse_angstrom'][1:])))
        self.assertAlmostEqual(result['distance_inverse_angstrom'][-1],2*np.pi/3)

    def test_unstable_spring_model_rejected(self):
        with self.assertRaises(PhononInputError):
            from_springs(self.cubic.crystal,[(0,0,(1,0,0),-np.eye(3))])
        b={k:v.copy() for k,v in self.cubic.blocks.items()}
        unstable=ForceConstants(self.cubic.crystal,{k:-v for k,v in b.items()})
        self.assertLess(unstable.modes([0.5,0,0])[0][0],0)
        with self.assertRaises(PhononInputError):
            density_of_states(frequency_grid(unstable,(2,2,2)))

    def test_dos_integral_states(self):
        f=frequency_grid(self.cubic,(5,5,5))
        grid,dos=density_of_states(f,0.12,1000)
        self.assertAlmostEqual(np.trapezoid(dos,grid),3,places=4)
        self.assertTrue((dos>=0).all())

    def test_thermodynamic_t0_and_classical_cv(self):
        freq=np.ones((2,3))*4
        zero=thermal_properties(freq,0)
        self.assertAlmostEqual(zero['zero_point_energy_ev_cell'],1.5*4*4.135667696e-3)
        self.assertEqual(zero['cv_ev_k_cell'],0)
        hot=thermal_properties(freq,1e5)
        self.assertAlmostEqual(hot['cv_ev_k_cell']/8.617333262145e-5,3,places=5)
        self.assertGreater(hot['entropy_ev_k_cell'],0)
        self.assertLess(hot['free_energy_ev_cell'],hot['internal_energy_ev_cell'])

    def test_thermal_reject_zero_modes_and_bad_temp(self):
        with self.assertRaises(PhononInputError): thermal_properties([[0,2,3]],300)
        with self.assertRaises(PhononInputError): thermal_properties([[1,2,3]],-1)

    def test_spectroscopy_conversion(self):
        self.assertAlmostEqual(float(thz_to_wavenumber_cm(1)),33.3564095198,places=8)

    def test_ir_selection_rules_binary(self):
        crystal=diatomic_chain(1,1).crystal
        charges=np.stack([np.eye(3),-np.eye(3)])
        acoustic=np.array([1,0,0,1,0,0])/math.sqrt(2)
        optic=np.array([1,0,0,-1,0,0])/math.sqrt(2)
        self.assertAlmostEqual(ir_mode_activity(crystal,acoustic,charges),0)
        self.assertAlmostEqual(ir_mode_activity(crystal,optic,charges),2)
        self.assertAlmostEqual(ir_mode_activity(crystal,optic,charges,polarization=[0,1,0]),0)
        with self.assertRaises(PhononInputError): ir_mode_activity(crystal,optic,charges+1)

    def test_raman_projection_requires_derivatives(self):
        crystal=diatomic_chain(1,1).crystal
        d=np.zeros((2,3,3,3))
        d[0,0,0,0]=1
        mode=np.array([1,0,0,0,0,0])
        self.assertAlmostEqual(raman_mode_activity(crystal,mode,d,[1,0,0],[1,0,0]),1)
        self.assertAlmostEqual(raman_mode_activity(crystal,mode,d,[1,0,0],[0,1,0]),0)
        with self.assertRaises(PhononInputError): raman_mode_activity(crystal,mode,d[:-1],[1,0,0],[1,0,0])

    def test_eliashberg_integrals(self):
        f=np.linspace(1,2,5001)
        props=eliashberg_moments(f,np.ones_like(f))
        self.assertAlmostEqual(props['lambda'],2*math.log(2),places=6)
        self.assertAlmostEqual(props['omega_log_thz'],math.sqrt(2),places=6)
        self.assertEqual(eliashberg_moments(f,np.zeros_like(f))['lambda'],0)
        with self.assertRaises(PhononInputError): eliashberg_moments([0,1],[0,1])

    def test_json_roundtrip_and_provenance(self):
        prov={"evidence_status":"illustrative","source":"synthetic fixture"}
        rec=model_to_record(self.cubic,prov)
        rebuilt, backprov=parse_model(rec)
        self.assertEqual(prov,backprov)
        for q in [[0,0,0],[0.2,0.3,0.1]]:
            np.testing.assert_allclose(self.cubic.modes(q)[0],rebuilt.modes(q)[0])
        with tempfile.TemporaryDirectory() as td:
            f=Path(td)/"model.json"
            save_json(rec,f)
            reloaded,_=load_json(f)
            self.assertEqual(len(reloaded.blocks),len(self.cubic.blocks))

    def test_missing_provenance_and_ifc_never_guessed(self):
        p={"evidence_status":"illustrative","source":"synthetic fixture"}
        rec=model_to_record(self.cubic,p)
        rec['force_constants']['terms']=[]
        with self.assertRaises(PhononInputError):parse_model(rec)
        rec=model_to_record(self.cubic,p)
        rec['provenance']={}
        with self.assertRaises(PhononInputError):parse_model(rec)
        rec=model_to_record(self.cubic,p)
        rec['force_constants']['units']='N/m'
        with self.assertRaises(PhononInputError):parse_model(rec)

    def test_export_does_not_invent_experimental_spectrum(self):
        prov={"evidence_status":"illustrative","source":"synthetic fixture"}
        result=export_result(self.cubic,prov,[[0,0,0],[0.5,0,0]],mesh=(3,3,3),points_per_segment=8)
        self.assertEqual(result['evidence']['declared_status'],'illustrative')
        self.assertIn('density_of_states',result)
        self.assertIn('thermodynamics',result)
        self.assertNotIn('raman',result)
        self.assertNotIn('ir',result)
        self.assertNotIn('electron_phonon_coupling',result)
        json.dumps(result,allow_nan=False)

    def test_unstable_export_withholds_thermodynamics(self):
        b={k:-v for k,v in self.cubic.blocks.items()}
        unstable=ForceConstants(self.cubic.crystal,b)
        result=export_result(unstable,{"evidence_status":"illustrative","source":"unstable fixture"},
                             [[0,0,0],[0.5,0,0]],mesh=(2,2,2),points_per_segment=2)
        self.assertIn('withheld',result['thermodynamics']['status'])
        self.assertGreater(result['physics']['unstable_modes_on_shifted_mesh'],0)


if __name__=='__main__': unittest.main()
