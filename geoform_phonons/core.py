"""Harmonic crystal phonons from explicit real-space interatomic force constants.

This is a *post-processor*, not a DFT calculator or an IFC predictor.
Conventions: lattice rows in Angstrom, fractional atom positions and reciprocal
q coordinates, IFC Phi[i0,jR] in eV/Angstrom**2, masses in u, frequencies
in signed THz (cycles/s, not rad/s); a negative frequency marks omega^2 < 0.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
import itertools
import math
import numpy as np

EV = 1.602176634e-19                  # exact joule per eV
ANGSTROM = 1e-10
AMU = 1.66053906660e-27               # kg/u, CODATA conventional value
PLANCK_EV_S = 4.135667696e-15
KB_EV_K = 8.617333262145e-5
LIGHT_CM_S = 2.99792458e10
THZ_FACTOR = math.sqrt(EV / ANGSTROM**2 / AMU) / (2.0 * math.pi * 1e12)


class PhononInputError(ValueError):
    """Invalid or insufficient input for a physically defined result."""


def signed_frequencies(eigenvalues):
    """Convert dynamical eigenvalues to signed THz, suppressing eigensolver roundoff."""
    values = np.asarray(eigenvalues, dtype=float)
    floor = 1e-11 * np.maximum(1.0, np.max(np.abs(values), axis=-1, keepdims=True))
    values = np.where(np.abs(values) <= floor, 0.0, values)
    return np.sign(values) * np.sqrt(np.abs(values)) * THZ_FACTOR


def _finite_array(value, shape, name):
    out = np.asarray(value, dtype=float)
    if out.shape != shape or not np.isfinite(out).all():
        raise PhononInputError(f"{name} must have shape {shape} and finite entries")
    return out


@dataclass(frozen=True)
class Crystal:
    """Primitive crystal cell. This establishes geometry, not force constants."""
    lattice_angstrom: object
    fractional_positions: object
    masses_amu: object
    symbols: tuple[str, ...]

    def __post_init__(self):
        cell = _finite_array(self.lattice_angstrom, (3, 3), "lattice_angstrom")
        pos = np.asarray(self.fractional_positions, dtype=float)
        mass = np.asarray(self.masses_amu, dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3 or pos.shape[0] == 0:
            raise PhononInputError("fractional_positions must be nonempty Nx3")
        if not np.isfinite(pos).all() or mass.shape != (len(pos),) or not np.isfinite(mass).all():
            raise PhononInputError("positions/masses must be finite with matching sizes")
        if any(m <= 0 for m in mass):
            raise PhononInputError("atomic masses must be strictly positive")
        if len(self.symbols) != len(pos) or not all(isinstance(s, str) and s for s in self.symbols):
            raise PhononInputError("symbols must name each atom")
        if abs(float(np.linalg.det(cell))) < 1e-12:
            raise PhononInputError("lattice must be nonsingular")
        object.__setattr__(self, "lattice_angstrom", cell)
        object.__setattr__(self, "fractional_positions", pos)
        object.__setattr__(self, "masses_amu", mass)
        object.__setattr__(self, "symbols", tuple(self.symbols))

    @property
    def n_atoms(self):
        return len(self.masses_amu)

    @property
    def volume_angstrom3(self):
        return abs(float(np.linalg.det(self.lattice_angstrom)))

    @property
    def reciprocal_vectors(self):
        """Rows are 2*pi reciprocal basis vectors, in inverse Angstrom."""
        return 2 * np.pi * np.linalg.inv(self.lattice_angstrom).T


class ForceConstants:
    """Finite real-space IFC blocks indexed by (Rx,Ry,Rz,i,j).

    Phi[(R,i,j)] denotes the Hessian coupling displacement of atom i in
    cell 0 to atom j in cell R. Every off-diagonal block needs its transpose
    counterpart at (-R,j,i). The constructor rejects incomplete pairs rather
    than silently manufacturing symmetry or claiming a physical model.
    """

    def __init__(self, crystal: Crystal, blocks: Mapping, *, tolerance=1e-8):
        self.crystal = crystal
        self.blocks = {}
        if not isinstance(blocks, Mapping) or not blocks:
            raise PhononInputError("nonempty IFC mapping required; CIF alone is not enough")
        if not np.isfinite(tolerance) or tolerance < 0:
            raise PhononInputError("tolerance must be finite and nonnegative")
        self.tolerance = tolerance
        for key, val in blocks.items():
            if len(key) != 5 or any(not isinstance(v, (int, np.integer)) for v in key):
                raise PhononInputError("IFC key must contain five integers: R and i,j")
            rx, ry, rz, i, j = (int(v) for v in key)
            if not (0 <= i < crystal.n_atoms and 0 <= j < crystal.n_atoms):
                raise PhononInputError("IFC atom index out of bounds")
            self.blocks[(rx, ry, rz, i, j)] = _finite_array(val, (3, 3), "IFC block")
        self._check_hermitian_pairs()

    def _check_hermitian_pairs(self):
        for (x, y, z, i, j), phi in self.blocks.items():
            reverse = self.blocks.get((-x, -y, -z, j, i))
            if reverse is None or not np.allclose(phi, reverse.T, atol=self.tolerance, rtol=self.tolerance):
                raise PhononInputError(f"missing or inconsistent transpose IFC block {(x,y,z,i,j)}")

    def asr_residual(self):
        """For each i, sum_{j,R} Phi[i0,jR], returning (N,3,3) tensors."""
        residual = np.zeros((self.crystal.n_atoms, 3, 3))
        for (_, _, _, i, _), block in self.blocks.items():
            residual[i] += block
        return residual

    def acoustic_sum_rule_error(self):
        return float(np.max(np.abs(self.asr_residual())))

    def enforce_asr(self):
        """Correct the onsite blocks for small translational drift only.

        Requires symmetric row-sum residuals. This is a limited onsite ASR
        projection, not a rotational sum-rule projection, and can change IFCs.
        """
        residual = self.asr_residual()
        updated = {key: val.copy() for key, val in self.blocks.items()}
        for i in range(self.crystal.n_atoms):
            if not np.allclose(residual[i], residual[i].T, atol=self.tolerance, rtol=self.tolerance):
                raise PhononInputError("ASR residual asymmetric; cannot repair onsite only")
            key = (0, 0, 0, i, i)
            updated[key] = updated.get(key, np.zeros((3, 3))) - residual[i]
        return ForceConstants(self.crystal, updated, tolerance=self.tolerance)

    def dynamical_matrix(self, q_fractional):
        """Mass-weighted 3N x 3N Hermitian matrix in (eV/A^2)/amu."""
        q = _finite_array(q_fractional, (3,), "q_fractional")
        n = self.crystal.n_atoms
        d = np.zeros((3*n, 3*n), dtype=complex)
        tau = self.crystal.fractional_positions
        mass = self.crystal.masses_amu
        for (rx, ry, rz, i, j), phi in self.blocks.items():
            r = np.array((rx, ry, rz), dtype=float) + tau[j] - tau[i]
            phase = np.exp(2j*np.pi*np.dot(q, r))
            d[3*i:3*i+3, 3*j:3*j+3] += phi * phase / math.sqrt(mass[i] * mass[j])
        if not np.allclose(d, d.conj().T, atol=self.tolerance, rtol=self.tolerance):
            raise PhononInputError("dynamical matrix not Hermitian")
        return (d + d.conj().T) * 0.5

    def dynamical_matrices(self, q_fractional):
        """Evaluate a batch of q points using the same periodic IFC convention."""
        q = np.asarray(q_fractional, dtype=float)
        if q.ndim != 2 or q.shape[1] != 3 or not len(q) or not np.isfinite(q).all():
            raise PhononInputError("q_fractional must be a finite nonempty (N,3) array")
        keys = np.asarray(sorted(self.blocks, key=lambda key: (key[3], key[4], key[:3])))
        i, j = keys[:, 3:].T
        n = self.crystal.n_atoms
        tau = self.crystal.fractional_positions
        mass = self.crystal.masses_amu
        translations = keys[:, :3] + tau[j] - tau[i]
        weights = np.stack([self.blocks[tuple(key)] for key in keys])
        weights /= np.sqrt(mass[i] * mass[j])[:, None, None]
        pair = i*n + j
        starts = np.r_[0, np.flatnonzero(np.diff(pair)) + 1]
        matrices = np.zeros((len(q), n, n, 3, 3), dtype=complex)
        # Bound Fourier workspace independently of the size of the q batch.
        chunk = max(1, min(32, 500_000 // len(keys)))
        for offset in range(0, len(q), chunk):
            section = slice(offset, offset + chunk)
            phase = np.exp(2j * np.pi * (q[section] @ translations.T))
            terms = phase[:, :, None, None] * weights
            matrices[section, i[starts], j[starts]] = np.add.reduceat(terms, starts, axis=1)
        matrices = matrices.transpose(0, 1, 3, 2, 4).reshape(len(q), 3*n, 3*n)
        adjoints = matrices.conj().swapaxes(-1, -2)
        if not np.allclose(matrices, adjoints, atol=self.tolerance, rtol=self.tolerance):
            raise PhononInputError("dynamical matrix not Hermitian")
        return (matrices + adjoints) * 0.5

    def modes(self, q_fractional):
        """Return (signed frequencies THz, normalized mass-weighted eigenvectors)."""
        eigenvalues, eigenvectors = np.linalg.eigh(self.dynamical_matrix(q_fractional))
        return signed_frequencies(eigenvalues), eigenvectors

    def band_path(self, vertices, points_per_segment=40):
        """Sample dispersion along a fractional-reciprocal polyline.

        Returns JSON-ready q-points, cumulative reciprocal distance (1/A),
        signed THz frequencies; bands are eigenvalue-ordered, not tracked
        across band crossings.
        """
        if not isinstance(points_per_segment, int) or points_per_segment < 2:
            raise PhononInputError("points_per_segment must be integer >= 2")
        vertices = [_finite_array(q, (3,), "path vertex") for q in vertices]
        if len(vertices) < 2:
            raise PhononInputError("path needs at least two vertices")
        points = []
        for k in range(len(vertices)-1):
            for t in np.linspace(0, 1, points_per_segment, endpoint=(k == len(vertices)-2)):
                points.append(vertices[k] * (1-t) + vertices[k+1] * t)
        reciprocal = self.crystal.reciprocal_vectors
        dist = [0.0]
        for a, b in zip(points[:-1], points[1:]):
            dist.append(dist[-1] + float(np.linalg.norm((b-a) @ reciprocal)))
        freqs = [self.modes(q)[0].tolist() for q in points]
        return {"q_fractional": [q.tolist() for q in points],
                "distance_inverse_angstrom": dist,
                "frequency_thz_signed": freqs,
                "band_order": "ascending eigenvalue (not continuous mode tracking)"}


def shifted_q_mesh(shape):
    """Monkhorst-style half-shift grid, excluding exact Gamma in any dimension."""
    if len(shape) != 3 or any(not isinstance(n, int) or n < 1 for n in shape):
        raise PhononInputError("mesh shape needs three positive integers")
    return [tuple((i+0.5)/n for i, n in zip(index, shape))
            for index in itertools.product(*(range(n) for n in shape))]


def frequency_grid(model: ForceConstants, mesh=(8,8,8)):
    """Return (Nq,3N) signed THz modes on a half-shift reciprocal mesh."""
    return np.array([model.modes(q)[0] for q in shifted_q_mesh(mesh)])


def density_of_states(frequencies_thz, sigma_thz=0.1, points=400, *, upper_thz=None):
    """Gaussian DOS in states/(THz*cell), integral approximately 3N.

    The input is an (Nq,3N) mesh; instability is rejected. Gaussian tails
    are included in the returned domain to limit integral truncation.
    """
    f = np.asarray(frequencies_thz, dtype=float)
    if f.ndim != 2 or not f.size or not np.isfinite(f).all():
        raise PhononInputError("frequencies must be finite 2D q-mesh array")
    if np.any(f < -1e-9):
        raise PhononInputError("imaginary modes: stable vibrational DOS unavailable")
    if not np.isfinite(sigma_thz) or sigma_thz <= 0 or not isinstance(points, int) or points < 8:
        raise PhononInputError("sigma>0 and points>=8 required")
    lo = -5*sigma_thz
    hi = max(float(f.max()) + 5*sigma_thz, float(upper_thz or 0))
    grid = np.linspace(lo, hi, points)
    dos = np.zeros_like(grid)
    # Batch by q-points, avoiding huge temporary meshes.
    for row in f:
        z = (grid[:,None]-row[None,:])/sigma_thz
        dos += np.exp(-0.5*z*z).sum(axis=1)
    dos /= (f.shape[0]*sigma_thz*math.sqrt(2*math.pi))
    return grid, dos



def projected_density_of_states(model: ForceConstants, mesh=(8,8,8),
                                sigma_thz=0.1, points=400):
    """Return (frequency, total DOS, atom-PDOS) from actual eigenvectors.

    Atom contributions use the mass-weighted, orthonormal eigenvectors and
    sum to the total DOS at each frequency. Array shape is (n_atoms, points).
    Each atom integrates to approximately three states over a broad axis.
    Do not confuse these with measured element-specific neutron intensities:
    no scattering lengths, Debye-Waller factors, or instrument resolution
    are used. A mesh with unstable or zero modes is rejected.
    """
    if not isinstance(model, ForceConstants):
        raise PhononInputError("model must be a ForceConstants instance")
    modes = [model.modes(q) for q in shifted_q_mesh(mesh)]
    frequencies = np.asarray([f for f, _ in modes])
    if np.any(frequencies <= 0):
        raise PhononInputError("unstable or zero frequency on q mesh: PDOS withheld")
    axis, total = density_of_states(frequencies, sigma_thz, points)
    atom_dos = np.zeros((model.crystal.n_atoms, len(axis)))
    norm = frequencies.shape[0] * sigma_thz * math.sqrt(2*math.pi)
    for f, eigvec in modes:
        atom_weights = np.sum(np.abs(eigvec.reshape(
            (model.crystal.n_atoms, 3, 3*model.crystal.n_atoms)))**2,
            axis=1)  # n_atoms x n_modes
        kernels = np.exp(-0.5 * ((axis[:,None] - f[None,:])/sigma_thz)**2)
        atom_dos += atom_weights @ kernels.T / norm
    if not np.allclose(atom_dos.sum(axis=0),total,atol=2e-10,rtol=1e-8):
        raise PhononInputError("projected DOS failed completeness check")
    return axis, total, atom_dos


def mode_participation_ratios(eigenvectors, n_atoms):
    """Return per-mode atomic participation ratios between 1/N and 1.

    Eigenvector columns are mass-weighted normal modes; project over each
    atom's three Cartesian components. A participation ratio is not a
    localization length or experimental observable without further modeling.
    """
    e = np.asarray(eigenvectors, dtype=complex)
    if (type(n_atoms) is not int or n_atoms < 1 or e.shape != (3*n_atoms,3*n_atoms)
            or not np.isfinite(e).all()):
        raise PhononInputError("eigenvectors must be finite 3N x 3N square matrix")
    norms = np.sum(np.abs(e)**2,axis=0)
    if not np.allclose(norms,1,atol=1e-7) or not np.allclose(e.conj().T@e,np.eye(3*n_atoms),atol=1e-7):
        raise PhononInputError("normal modes must be orthonormal")
    per_atom = np.sum(np.abs(e.reshape(n_atoms,3,3*n_atoms))**2,axis=1)
    return 1.0 / (n_atoms * np.sum(per_atom**2,axis=0))


def thermodynamic_curve(frequencies_thz, temperatures_k):
    """Harmonic constant-volume thermodynamics at user-supplied K values."""
    ts=np.asarray(temperatures_k,dtype=float)
    if ts.ndim!=1 or not ts.size or not np.isfinite(ts).all() or np.any(ts<0) or np.any(np.diff(ts)<=0):
        raise PhononInputError("temperatures must be strictly increasing, finite, and nonnegative")
    return [thermal_properties(frequencies_thz,float(t)) for t in ts]

def thermal_properties(frequencies_thz, temperature_k):
    """Harmonic per-cell thermodynamics from uniform positive-frequency q mesh.

    F and U in eV/cell, S and Cv in eV/(K cell). Mesh must avoid exact
    translational zeros (e.g. shifted_q_mesh); unstable/zero modes rejected.
    Includes quantum zero-point energy. No quasiharmonic expansion is modeled.
    """
    f = np.asarray(frequencies_thz, dtype=float)
    if f.ndim != 2 or not f.size or not np.isfinite(f).all():
        raise PhononInputError("frequencies must be finite 2D q-mesh array")
    if not np.isfinite(temperature_k) or temperature_k < 0:
        raise PhononInputError("temperature must be finite and >=0 K")
    if np.any(f <= 0):
        raise PhononInputError("zero/imaginary modes: choose a shifted mesh or repair model")
    energies = PLANCK_EV_S * 1e12 * f
    zpe = float(np.sum(0.5*energies) / f.shape[0])
    if temperature_k == 0:
        return {"temperature_k": 0.0, "free_energy_ev_cell": zpe,
                "internal_energy_ev_cell": zpe, "entropy_ev_k_cell": 0.0,
                "cv_ev_k_cell": 0.0, "zero_point_energy_ev_cell": zpe}
    x = energies/(KB_EV_K*temperature_k)
    logpart = np.log(-np.expm1(-x))
    bose = np.zeros_like(x)
    cvper = np.zeros_like(x)
    # expm1 stable for small x; for large x use asymptotic decay.
    safe = x < 700
    bose[safe] = energies[safe]/np.expm1(x[safe])
    cvper[safe] = KB_EV_K * (x[safe]/(2*np.sinh(x[safe]/2)))**2
    cvper[~safe] = 0.0
    u = float((np.sum(0.5*energies+bose))/f.shape[0])
    free = float((np.sum(0.5*energies + KB_EV_K*temperature_k*logpart))/f.shape[0])
    cv = float(np.sum(cvper)/f.shape[0])
    entropy = (u-free)/temperature_k
    return {"temperature_k": float(temperature_k), "free_energy_ev_cell": free,
            "internal_energy_ev_cell": u, "entropy_ev_k_cell": entropy,
            "cv_ev_k_cell": cv, "zero_point_energy_ev_cell": zpe}


def thz_to_wavenumber_cm(freq_thz):
    """Convert signed cyclic frequency THz to spectroscopic cm^-1."""
    return np.asarray(freq_thz) * 1e12 / LIGHT_CM_S
