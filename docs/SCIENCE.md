# Scientific contract, conventions and falsification checks

This document distinguishes **implemented identities**, **automated tests**, and **unvalidated scientific extrapolations**. No result from the synthetic fixtures is evidence about a particular mineral.

## Input conventions

Primitive lattice vectors occupy the *rows* of `A` and have units Å. Site coordinates `tau_i` are fractional. Wavevectors `q` are in fractional reciprocal coordinates; Cartesian reciprocal vectors are the rows of `B=2π(A^{-1})^T` (Å⁻¹). The IFC tensor `Phi[i0,jR]` is the second derivative of harmonic energy with respect to Cartesian displacements of site `i` in cell zero and site `j` in translated cell `R`: eV Å⁻². Atomic masses are in daltons (unified atomic mass units).

Dynamical matrix (dimension `3N × 3N`):

`D_{iα,jβ}(q) = [m_i m_j]^(-1/2) Σ_R Phi_{iα,jβ}(R) exp[2πi q·(R + tau_j - tau_i)]`.

With the above units its eigenvalues are `(eV/Å²)/amu`. For a positive eigenvalue `λ`, the cyclic frequency in THz is `ν = sqrt(λ) sqrt(eV/(Å² amu)) / (2π × 10¹²)`. An eigenvalue `λ<0` yields a signed negative frequency `−sqrt(|λ|)*factor` to distinguish imaginary angular modes. The sign is a display flag, not a physical negative oscillation frequency. Values close to zero are rounded using a relative numerical floor.

Hermiticity requires `Phi[i,j,R] = Phi[j,i,-R]^T`; the importer rejects missing counterparts. Translation requires `Σ_{j,R} Phi[i,j,R] = 0` for every site and Cartesian component. An optional **onsite correction** changes `Phi[i,i,0]` to enforce this row sum, provided residuals are symmetric. It is *not* a rotational sum-rule projection or a physical validation of the original calculation.

## Equilibrium thermodynamics

Thermal values use a uniform **half-shifted** reciprocal mesh to avoid the three exact translational Gamma zeros. They require strictly positive phonon frequencies on that mesh and report all values per primitive cell. For each mode with `E=hν`, `x=E/k_B T`:

`F = mean_q Σ_modes [E/2 + k_B T log(1−e^(−x))]`

`U = mean_q Σ_modes [E/2 + E/(e^x−1)]`

`Cv = mean_q Σ_modes k_B [x/(2 sinh(x/2))]^2` and `S=(U−F)/T`.

At `T=0`, `F=U=ZPE`, `S=Cv=0`. At high temperature `Cv→3N k_B`. The output gives energies in eV/cell and entropy/heat capacity in eV/(K cell). These are harmonic **constant-volume** quantities: no quasi-harmonic thermal expansion, isotope disorder, anharmonic linewidth, stability at finite temperature, or electronic free energy.

The Gaussian density of states is `g(ν)=mean_q Σ_modes exp[-(ν−ν_q)²/(2σ²)]/(σ sqrt(2π))`, normalized to **3N states per primitive cell** on a sufficiently wide frequency axis. Gaussian tails may extend slightly below zero; those are numerical broadening, not imaginary eigenmodes.

## Spectroscopy and coupling are separate evidence layers

- Γ-point IR relative mode activity requires an externally supplied Born-effective-charge tensor `Z*[atom,polarization,displacement]`. It tests `Σ_atom Z*=0` and projects the mass-weighted mode into an induced dipole. Neither optical absorption coefficients nor peak widths are predicted.
- Γ-point polarized Raman relative activity requires `dα/du[atom,displacement,out,in]` and the experimental polarization geometry. Tensor response is *not* derivable from bare harmonic force constants.
- Eliashberg `alpha²F(ν)` is independently supplied (in a consistent THz frequency convention). `λ=2∫alpha²F(ν)/ν dν` and `ν_log=exp[(2/λ)∫alpha²F(ν)log(ν)/ν dν]`. No bare-IFC electron–phonon coupling, superconductivity claim, or McMillan temperature is computed.
- Polar material LO–TO splitting needs a non-analytic long-range Coulomb correction using Born charges, the electronic dielectric tensor, the physical cell volume and compatible unit conventions. **Not implemented.**

## Tests and evidence boundaries

The reproducible suite checks reciprocal geometry, 3-D analytic simple-cubic zone limits, exact 1-D acoustic/optic matching against the independent original chain implementation, mass `m^{-1/2}` scaling, periodicity under reciprocal translations, Hermiticity at complex q, Gamma acoustic modes, ASR diagnostics and onsite projection, DOS normalization, ZPE, classical heat capacity, relative IR/Raman selection rules, Eliashberg integral identities, input/schema rejection, and unstable-mode withholding.

These are **analytical and synthetic** checks, not validation against independent real-crystal data or a Phonopy run. A real-material benchmark requires frozen versioned structures, compatible primitive/supercell IFCs, Born/dielectric inputs when relevant, measured temperature-pressure conditions, measurement uncertainty, a non-leaking holdout, and error analysis separated by mineral group and polymorph. Substitution heuristics or formula-only chemistry must never be relabeled as measured vibrational data.

## GeoForm contract

Input schema `geoform.phonons/1`: `model_kind = supplied_ifc`; `crystal` with lattice/positions/masses/symbols; `force_constants` with exact `units = eV/angstrom^2` and `terms[]` containing integer `R`, integer `i,j`, and 3×3 `phi`; `provenance` with `evidence_status` and a nonempty `source`. Output schema `geoform.phonons.result/1` separates declared provenance, crystal metadata, physics diagnostics, dispersion, DOS and thermodynamics. Source status is **a user declaration, not verified authentication**.

Reference implementations: [Phonopy](https://phonopy.github.io/phonopy/) and [Quantum ESPRESSO PHonon](https://www.quantum-espresso.org/Doc/ph_user_guide/). Their existence is not evidence that GeoForm results have been cross-validated against them.
