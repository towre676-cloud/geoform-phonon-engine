# Source-pinned harmonic phonons for nine crystal structures

GeoForm compares its own complex dynamical matrices and signed frequencies
with Phonopy 4.9.0 on the same supplied VASP forces. These are **calculated**
results. The comparison validates the finite-supercell import and numerical
postprocessing for these inputs. Experimental spectra and DFT accuracy were
outside the validation scope.

## Added calculations and phase assignments

| Case | Computational phase | Space group | Force supercell | Atoms, primitive / supercell | Displaced runs | Actual force-run k mesh |
| --- | --- | --- | --- | --- | --- | --- |
| Corundum, Al₂O₃ | Corundum | R-3c (167) | 2×2×1 | 10 / 120 | 5 | 3×3×1 |
| Anatase, TiO₂ | Anatase | I4₁/amd (141) | 4×4×1 | 6 / 192 | 4 | 1×1×1 |
| Zincite, ZnO | Wurtzite type | P6₃mc (186) | 2×2×2 | 4 / 32 | 6 | 3×3×2 |
| Cassiterite, SnO₂ | Rutile type | P4₂/mnm (136) | 2×2×3 | 6 / 72 | 3 | 2×2×2 |
| Cubic CaTiO₃ | Ideal cubic perovskite type | Pm-3m (221) | 2×2×2 | 5 / 40 | 4 | 2×2×2 |
| Stishovite, SiO₂ | High-pressure rutile type | P4₂/mnm (136) | 2×2×2 | 6 / 48 | 3 | 2×2×3 |
| Hexagonal MgB₂ | AlB₂ type | P6/mmm (191) | 3×3×2 | 3 / 54 | 2 | 6×6×9 |

The audit checks the exact POSCAR against the displacement YAML, recomputes
space groups and Wyckoff letters with spglib, verifies the composition and
cell reduction, and compares each displaced VASP run's atom order, initial
lattice and fractional positions against the YAML. Each of the 27 final
VASP force arrays equals its `FORCE_SETS` array at the stored precision.
All 27 runs report a 500 eV cutoff. The per-run VASP version, GGA tag, k mesh
and geometry errors are in the receipt. A `--` GGA tag is retained as reported
by VASP. The workflow neither acquires POTCAR files nor reruns VASP.

- **Anatase.** The upstream README says 4×4×2. The YAML, the 192 force rows
  per displacement, and all four archived VASP lattices agree on 4×4×1.
  The YAML's 12-atom body-centered conventional basis is reduced to a
  six-atom primitive basis by Phonopy. The receipt and result retain the
  transformation and basis, respectively.
- **Corundum.** The archived displaced calculations use 3×3×1. The archived
  30-atom Born calculation uses 6×6×1, which is also the mesh cited in the
  README. The benchmark uses the displaced forces.
- **Cassiterite.** The displaced supercell runs use 2×2×2. The archived
  undistorted run 000 and separate Born run are outside the fit. The Born run
  uses the 4×4×6 mesh mentioned in the README.
- **Cubic CaTiO₃.** The verified coordinates describe the ideal cubic model.
  Its catalog entry has no mineral mapping. The separately hosted Mineral
  Lab's perovskite record remains unresolved for this calculation.

- **Stishovite.** The upstream README reports relaxation with `PSTRESS = 400`
  (400 kbar, or 40 GPa). The archived displacement inputs omit that tag.
  The result applies to the supplied compressed cell. The undistorted
  `vasprun.xml-sposcar` run is outside the fit.
- **Hexagonal MgB₂.** Both archived runs use first-order Methfessel–Paxton
  smearing (`ISMEAR = 1`) with 0.2 eV width. This compound has a
  computational-phase entry. The harmonic result carries no electron–phonon
  coupling or superconducting transition temperature.

## Numerical results and scope

Each case uses 16 fixed primitive reciprocal coordinates, including Γ, zone
edges, off-symmetry coordinates and a translated q point. The seven additions
repeat the checks with full IFC arrays expanded from the fitted compact arrays.
Phonopy supplies shortest-image geometry and fits the force constants. GeoForm
splits each finite-supercell pair among equivalent shortest images, checks
translation and reciprocity, and independently evaluates the Hermitian
Fourier sum and eigensystem. Both codes use the same cell, masses and IFCs.

The declared Phonopy 4.9.0 Type I finite-difference pipeline uses its symfc
projection. A diagnostic run forcing the alternate traditional iterative
symmetrizer gave corundum an ASR residual of approximately 4.58×10⁻⁵ eV/Å²,
above the 10⁻⁶ eV/Å² export gate. The retained projected calculation passes.
GeoForm applies no additional ASR correction. The small residual frequency
differences are consistent with the recorded conversion factors and near-zero
translational mode handling. The receipt records both factors and all
pointwise discrepancies.

| Case | Maximum complex matrix error (eV/Å²/amu) | Maximum frequency error (THz) | DOS and thermal response |
| --- | ---: | ---: | --- |
| Anatase | 1.02×10⁻¹⁴ | 3.03×10⁻⁶ | Included |
| Cassiterite | 1.11×10⁻¹⁵ | 2.76×10⁻⁶ | Included |
| Corundum | 9.18×10⁻¹⁵ | 3.24×10⁻⁶ | Included |
| Cubic CaTiO₃ | 1.33×10⁻¹⁵ | 3.19×10⁻⁶ | Withheld |
| Halite | 1.53×10⁻¹⁶ | 7.74×10⁻⁷ | Included |
| Periclase | 8.88×10⁻¹⁶ | 2.18×10⁻⁶ | Included |
| Stishovite | 1.33×10⁻¹⁵ | 4.10×10⁻⁶ | Withheld |
| Hexagonal MgB₂ | 1.78×10⁻¹⁵ | 2.97×10⁻⁶ | Included |
| Zincite | 5.55×10⁻¹⁶ | 1.91×10⁻⁶ | Withheld |

The cubic model has 816 imaginary modes over the 512-point shifted mesh
(counting branches at every q) and 178 on the 96-point dispersion path.
Its path minimum is −6.00469 THz at (½, ½, ½), with negative signed THz
denoting imaginary frequencies. The dispersion retains those modes.
Zincite has one imaginary mode on the path, −0.14346 THz at
(1/48, 0, 0), despite the shifted mesh and 16 comparison points missing
the instability. The Phonopy replay includes all 96 exported path points
and reproduces that soft mode. Its physical origin is unresolved, so the
catalog withholds equilibrium DOS and thermal quantities. Stishovite also has three sampled imaginary path modes, reaching
−0.49417 THz near Γ. Phonopy reproduced this minimum, and the shifted
mesh missed it. The physical cause of the mode is unresolved. The other six
models have no negative modes in the 512-point shifted mesh, 96-point path,
or 16 comparison points. This is sampled stability. No
global stability or supercell, k-point, or thermal mesh convergence bound
is asserted.

The reciprocal path is a deterministic primitive-coordinate polyline.
Branches are eigenvalue-ordered and lack continuous eigenvector tracking.
Stable cases have 0.1 THz Gaussian DOS, atom projections normalized to
three modes per primitive atom, 300 K harmonic properties, and a
constant-volume curve over 0–1200 K. Values are per primitive cell.
The datasets exclude polar nonanalytic LO–TO corrections on both sides,
including for examples with pinned `BORN` files. They carry no Raman or IR
intensities, experimental scattering weights, anharmonic shifts, thermal
expansion, conductivity or mineral specimen assignments.

## Other pinned candidates

- **Chromium** remains unresolved. The YAML assigns alternating collinear
  moments and yields magnetic space-group UNI 1599 (type IV), but the
  archived VASP XML reports `ISPIN = 1` and all-one `MAGMOM` defaults. The
  supplied archive does not substantiate a magnetic force calculation.
- **NaCl-rd** uses simultaneous random displacements and a Type II fit.
  It is a separate fitting-method candidate for halite, which already has
  a source-backed result.

## Reproduce and integrate

```sh
python -m pip install -e . -r benchmarks/requirements.txt
python scripts/fetch_reference_sources.py --destination benchmarks/upstream
python scripts/benchmark_minerals.py --data-dir benchmarks/upstream --output-dir benchmarks/generated --mesh 8 --points 24
python scripts/build_site_bundle.py --results benchmarks/generated
```

`benchmarks/reference_sources.json` pins the upstream commit, file sizes,
Git blob SHA-1 and SHA-256 for every consumed source and each provenance file,
including the BSD-3-Clause license. Downloading and validating local inputs
both check the pins. No upstream raw files are checked into this repository.
The expanded catalog and numerical receipts are checked into
`site/geoform-phonons/data/`. Pass `--stage-handoff` when intentionally
refreshing those results and their asset hashes. CI re-downloads the
immutable inputs, repeats the source and numerical audits, and packages the
verified outputs. The browser viewer and evidence card are integration assets
for the separately hosted Mineral Lab.

- [Pinned Phonopy example directories](https://github.com/phonopy/phonopy/tree/bd260e234de32b3b6c558c60a7c6254b9bdfef71/example)
- [Pinned Phonopy license](https://github.com/phonopy/phonopy/blob/bd260e234de32b3b6c558c60a7c6254b9bdfef71/LICENSE)
