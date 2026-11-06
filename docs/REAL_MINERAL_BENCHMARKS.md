# Halite + periclase: genuinely sourced computational phonon checks

**Scientific status: calculated DFT-derived finite-displacement test inputs.**
This is not a claim to have measured vibrational spectra or to have validated
which DFT functional is the best description of either mineral.

Phonopy's own maintained examples supply `phonopy_disp.yaml` and `FORCE_SETS`
for NaCl (halite) and MgO (periclase); their data are real VASP displacement
force results. `benchmarks/reference_sources.json` pins the complete upstream
Phonopy Git commit and **each file's exact Git-blob SHA-1** (which includes the
blob framing). The source files are downloaded only in the explicit benchmark
workflow and are **not** fabricated or relabeled as original GeoForm data.
`FORCE_SETS` are DFT-forces supplied to Phonopy; they are not the same file as
`FORCE_CONSTANTS`.

## Reproduce a full run

```sh
python -m pip install -e . 'phonopy==4.9.0'
python scripts/fetch_reference_sources.py --destination benchmarks/upstream
python scripts/benchmark_minerals.py --data-dir benchmarks/upstream --output-dir benchmarks/generated --mesh 4 --points 12
```

Output files in `benchmarks/generated/`:

- `halite_phonons.json`, `periclase_phonons.json`: validated GeoForm reader
  result schema, bands, DOS, atom-resolved DOS and constant-volume thermal curve.
- `benchmark_receipts.json`: independent dynamical-matrix and frequency errors
  at 16 fixed q points for each material, input SHA256 receipts, ASR residuals,
  primitive size and Phonopy version.
- `catalog.json`: minimal browser reference catalog with `calculated` status.

Hosted CI must reproduce these checks, not merely print a claimed PASS. Its
artifact can be downloaded from a passing workflow run. Inspect the recorded
numeric errors and provenance before distributing a site's dataset.

## The key Fourier-interpolation subtlety

Phonopy stores a finite-supercell Hessian, but the infinite crystal requires
lattice translations for each primitive pair. For a pair at a supercell
boundary, there may be **multiple equally shortest periodic vectors**. A
converter that picks only one image will generally give a wrong q-dependent
dynamical matrix even if it reproduces the Gamma spectrum.

`geoform_phonons.phonopy_native.from_phonopy()` imports Phonopy's independently
computed shortest vectors/multiplicity and distributes the pair force constant
equally among the shortest images. Every fractional displacement is checked
to be an integer-translated primitive pair. GeoForm then computes its **own**
Hermitian Fourier sum and eigenvalues, with no call to Phonopy during those
calculations. Invalid shapes/terms and inconsistent transpose blocks fail.

At each of 16 fixed q points including complex (off-symmetry) locations, the
benchmark compares (i) full complex 3N×3N dynamical matrices, not just
frequencies, and (ii) sorted THz frequencies computed independently by the
two packages. The exact tolerances are in `scripts/benchmark_minerals.py`.
An explicit acoustic residual check must pass, and the source is not silently
repaired during GeoForm import. Phonopy's force-constant symmetrization used
during FC production is recorded in the output evidence.

## Deliberate physics omissions

- MgO and NaCl are polar. The benchmark **turns off** the long-range
  nonanalytic Born-charge correction on both sides. The resulting Γ optical
  modes must never be interpreted as reproducing LO–TO splitting. Independent
  `BORN` files exist upstream but are not used as a hidden adjustment.
- Calculated harmonic frequencies are **not measured FTIR/Raman spectra**.
- The q path is a deterministic test path, not a full automatic
  symmetry-labelled Brillouin-zone path.
- DFT settings for the MgO example are reported upstream as 500 eV VASP
  cutoff and 4×4×4 k point mesh for the conventional cell. Other metadata
  should not be invented when absent in the source.
- This is **algorithmic agreement with Phonopy given the same inputs**, not
  blind experimental validation, uncertainty quantification, QHA,
  electron–phonon coupling or mineral database-wide generalization.

## Upstream and license

- [Phonopy NaCl example](https://github.com/phonopy/phonopy/tree/bd260e234de32b3b6c558c60a7c6254b9bdfef71/example/NaCl)
- [Phonopy MgO example](https://github.com/phonopy/phonopy/tree/bd260e234de32b3b6c558c60a7c6254b9bdfef71/example/MgO)
- [Phonopy Python API](https://phonopy.github.io/phonopy/phonopy-module.html)
- [Upstream BSD-3-Clause license](https://github.com/phonopy/phonopy/blob/bd260e234de32b3b6c558c60a7c6254b9bdfef71/LICENSE)

When redistributing upstream raw files, retain its copyright/license text;
this repository's benchmark downloader instead retrieves the original bytes
from the official pinned source.

## GeoForm site adapter

`web/phonon_explorer.html` accepts a `geoform.phonons.catalog/1` catalog
from the default relative `data/catalog.json` while served over HTTP(S).
Select an independently compared case from the dropdown; the browser then
loads the corresponding `*_phonons.json`, displaying its original evidence
and the user's selected phonon/thermodynamic data. Alternatively pass a
`?catalog=https://some-host/path/catalog.json` URL whose server permits CORS.
The existing manual JSON file input still works offline.

To serve the generated data with the viewer, place its four JSON files under
`web/data/` and serve `web/` from a local static HTTP server. This is a
handoff contract, not a claim that the live `geoform-mineral-lab` site has
already been updated.
