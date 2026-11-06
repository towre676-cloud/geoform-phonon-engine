"""Independent installed-Phonopy cross-check: 3x3x3 analytic cubic lattice.

This script intentionally imports Phonopy, not GeoForm, for constructing
external frequencies. It uses identical physical harmonic parameters, but
an independently constructed full-supercell Hessian. No mineral claim.
"""
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from geoform_phonons.fixtures import simple_cubic


def main():
    mass = 28.085
    k = 4.0
    a = 5.0
    atom = PhonopyAtoms(symbols=['Si'], masses=[mass],cell=np.eye(3)*a,
                       scaled_positions=[[0,0,0]])
    phonon = Phonopy(atom, supercell_matrix=np.diag([3,3,3]),
                    primitive_matrix=np.eye(3))
    coords = phonon.supercell.positions/a
    n = len(coords)
    assert n == 27
    full = np.zeros((n,n,3,3))
    for i in range(n):
        for j in range(n):
            delta = np.rint(coords[j]-coords[i]).astype(int)%3
            offset = sum(d!=0 for d in delta)
            if i==j: full[i,j] = np.eye(3)*6*k
            elif offset==1: full[i,j] = -np.eye(3)*k
    assert np.allclose(full,full.transpose(1,0,3,2))
    assert np.max(np.abs(full.sum(axis=1))) < 1e-12
    phonon.force_constants = full
    geoform = simple_cubic(mass_amu=mass, stiffness_ev_a2=k, a_angstrom=a)
    qpoints = [[0,0,0],[0.125,0.25,0.375],[0.5,0,0],[0.5,0.5,0.5],
               [0.23,0.11,0.08]]
    for q in qpoints:
        actual = np.sort(phonon.get_frequencies(q))
        predicted = geoform.modes(q)[0]
        np.testing.assert_allclose(predicted,actual,rtol=3e-6,atol=3e-5,
                                   err_msg=f'Phonopy/GeoForm disagreement at q={q}')
        print(f'q={q}: Phonopy {actual.tolist()}, GeoForm {predicted.tolist()}')
    print('PASS: independent Phonopy full-supercell dispersion comparison')


if __name__=='__main__': main()
