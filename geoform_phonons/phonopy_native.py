"""Correctly unfold Phonopy's finite-supercell IFCs into periodic image terms.

This adapter deliberately uses Phonopy for its *geometric shortest-image*
convention only. Fourier sums, dynamical matrices, eigenvalues and export are
computed by GeoForm. Input must be a Phonopy instance with actual IFCs;
there is no formula-to-force-constant inference.

Only the short-range harmonic (non-NAC) dynamical matrix can be reproduced.
The explicit image multiplicity is crucial for supercell-boundary bonds.
"""
from __future__ import annotations
from collections import defaultdict
import numpy as np
from .core import Crystal, ForceConstants, PhononInputError


def from_phonopy(phonon, *, tolerance=1e-7, integer_tolerance=2e-6):
    """Return periodic GeoForm IFCs equivalent to Phonopy's bare harmonic IFCs.

    Works with full (Ns,Ns,3,3) and compact (Np,Ns,3,3) force constants.
    Each finite-supercell pair may have several equidistant *primitive-cell*
    lattice-image vectors. Divide its IFC equally among those shortest
    images as Phonopy does, rather than arbitrarily selecting one image.

    Note: imports the force constants that Phonopy has already produced and
    (possibly) symmetrized. This function does not symmetrize data. A
    non-reciprocal or ASR-violating source is rejected / reported by the
    resulting ForceConstants; the source is not silently repaired.
    """
    primitive = phonon.primitive
    supercell = phonon.supercell
    f = np.asarray(phonon.force_constants, dtype=float)
    n = len(primitive)
    ns = len(supercell)
    if f.shape not in ((n, ns, 3, 3), (ns, ns, 3, 3)) or not np.isfinite(f).all():
        raise PhononInputError('Phonopy must contain finite full or compact IFCs')
    if not np.isfinite(tolerance) or tolerance < 0 or not (0 < integer_tolerance < 0.01):
        raise PhononInputError('invalid adapter tolerance')
    if not bool(primitive.store_dense_svecs):
        # Public helper matches the representation that Phonopy's own
        # dynamical-matrix implementation uses for sparse shortest vectors.
        from phonopy.structure.cells import sparse_to_dense_svecs
        sv, multiplicity = sparse_to_dense_svecs(*primitive.get_smallest_vectors())
    else:
        sv, multiplicity = primitive.get_smallest_vectors()
    sv = np.asarray(sv, dtype=float)
    multiplicity = np.asarray(multiplicity)
    if sv.ndim != 2 or sv.shape[1] != 3 or multiplicity.shape != (ns, n, 2):
        raise PhononInputError('Phonopy shortest-vector layout incompatible with adapter')
    if not np.isfinite(sv).all():
        raise PhononInputError('nonfinite shortest image vectors')
    p2s = np.asarray(primitive.p2s_map)
    s2p = np.asarray(primitive.s2p_map)
    p2p = primitive.p2p_map
    if p2s.shape != (n,) or s2p.shape != (ns,):
        raise PhononInputError('incorrect Phonopy primitive/supercell mappings')
    crystal = Crystal(primitive.cell, primitive.scaled_positions,
                      primitive.masses, tuple(primitive.symbols))
    tau = crystal.fractional_positions
    blocks = defaultdict(lambda: np.zeros((3, 3)))
    for i, super_i in enumerate(p2s):
        if f.shape[0] == ns:
            row = f[int(super_i)]
        else:
            row = f[i]
        for k in range(ns):
            try:
                j = int(p2p[int(s2p[k])])
            except (KeyError, IndexError, TypeError) as exc:
                raise PhononInputError('invalid Phonopy s2p/p2p mapping') from exc
            count, start = map(int, multiplicity[k, i])
            if count <= 0 or start < 0 or start + count > len(sv):
                raise PhononInputError('invalid Phonopy shortest-vector address/multiplicity')
            for vec in sv[start:start+count]:
                # vec is an i->j displacement in primitive fractional units.
                # The IFC convention requires integer R and exp(i q*(R+tau_j-tau_i)).
                translation = vec - (tau[j] - tau[i])
                rounded = np.rint(translation)
                if np.max(np.abs(translation-rounded)) > integer_tolerance:
                    raise PhononInputError('shortest image not integer-related to primitive sites')
                r = tuple(int(x) for x in rounded)
                blocks[(*r, i, j)] += row[k] / count
    # No manufactured reverse blocks: source must satisfy reciprocity.
    return ForceConstants(crystal, dict(blocks), tolerance=tolerance)
