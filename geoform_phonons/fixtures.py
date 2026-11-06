"""Clearly labeled synthetic harmonic reference models, NEVER mineral IFCs."""
import numpy as np
from .core import Crystal, ForceConstants, PhononInputError


def from_springs(crystal, links):
    """Create translationally invariant IFCs from undirected link tuples.

    Each link (i, j, R, K) connects atom i in cell 0 to atom j in cell R,
    where K is a symmetric positive semidefinite 3x3 force matrix in eV/A².
    Provide one representative per undirected periodic bond; no deduplication
    is attempted. Model is illustrative and need not obey rotational ASR.
    """
    blocks = {}
    def add(key, val):
        blocks[key] = blocks.get(key, np.zeros((3,3))) + val
    for i, j, r, raw_k in links:
        if not (isinstance(i,int) and isinstance(j,int) and 0 <= i < crystal.n_atoms and 0 <= j < crystal.n_atoms):
            raise PhononInputError("invalid spring atom indices")
        if len(r)!=3 or any(not isinstance(x,int) for x in r):
            raise PhononInputError("spring cell index must be 3 integers")
        r = tuple(r)
        if i == j and r == (0,0,0):
            raise PhononInputError("spring cannot connect atom to itself in same cell")
        k = np.asarray(raw_k, dtype=float)
        if k.shape != (3,3) or not np.isfinite(k).all() or not np.allclose(k,k.T):
            raise PhononInputError("spring K must be finite symmetric 3x3")
        if np.linalg.eigvalsh(k).min() < -1e-10:
            raise PhononInputError("spring K must be positive semidefinite")
        add((0,0,0,i,i), k)
        add((0,0,0,j,j), k)
        add((*r,i,j), -k)
        add((*(-x for x in r),j,i), -k)
    return ForceConstants(crystal, blocks)


def simple_cubic(mass_amu=28.085, stiffness_ev_a2=4.0, a_angstrom=5.0):
    """1-atom cubic toy with 3 identical polarizations, strictly harmonic."""
    if min(mass_amu,stiffness_ev_a2,a_angstrom)<=0:
        raise PhononInputError("toy parameters must be positive")
    crystal = Crystal(np.eye(3)*a_angstrom, [[0,0,0]], [mass_amu], ("X",))
    links = [(0,0,(1,0,0),np.eye(3)*stiffness_ev_a2),
             (0,0,(0,1,0),np.eye(3)*stiffness_ev_a2),
             (0,0,(0,0,1),np.eye(3)*stiffness_ev_a2)]
    return from_springs(crystal, links)


def diatomic_chain(mass_a=1.0, mass_b=2.0, stiffness_ev_a2=3.0, a_angstrom=4.0):
    """Analytic AB chain x branch, zero transverse modes (not 3D-stable)."""
    if min(mass_a,mass_b,stiffness_ev_a2,a_angstrom)<=0:
        raise PhononInputError("toy parameters must be positive")
    crystal = Crystal(np.diag([a_angstrom,20,20]), [[0,0,0],[0.5,0,0]],
                      [mass_a,mass_b], ("A","B"))
    k=np.diag([stiffness_ev_a2,0,0])
    return from_springs(crystal, [(0,1,(0,0,0),k),(0,1,(-1,0,0),k)])
