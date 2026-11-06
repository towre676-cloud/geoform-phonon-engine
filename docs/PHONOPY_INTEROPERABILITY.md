# Strict Phonopy IFC interchange and external benchmarking

## Supported

GeoForm reads the documented **text** `FORCE_CONSTANTS` tensor in either full
`(Ns,Ns,3,3)` or compact `(Np,Ns,3,3)` layout. This file by itself does **not**
identify primitive images, their lattice translations, primitive origin rows,
or their associated crystal structure. These must be supplied in a separate,
explicit JSON mapping; do not guess mappings based solely on ordering or formula.
`force_constants.hdf5` and `phonopy.yaml` are **not** parsed by this adapter.

Phonopy's force-constants convention: https://phonopy.github.io/phonopy/input-files.html

### Required primitive structure JSON

```json
{
  "lattice_angstrom": [[5,0,0],[0,5,0],[0,0,5]],
  "fractional_positions": [[0,0,0]],
  "masses_amu": [28.085],
  "symbols": ["Si"]
}
```

### Required mapping JSON (illustrative 1x1x1 case only)

```json
{
  "supercell_images": [{"primitive_index":0,"translation":[0,0,0]}],
  "origin_rows": [0]
}
```

For a real supercell list **every** supercell atom in column order, with the
primitive atom index and integer translation relative to the primitive lattice.
The `origin_rows[i]` value identifies which supercell row represents primitive
atom `i` in a *full* tensor; the corresponding translation is subtracted from
each column's translation. In a compact tensor, row `i` is taken as the supplied
compact primitive row, using the specified representative for its translation.
No wraparound, long-range alias splitting, missing transpose completion, mass
or stiffness guess, or acoustic correction is performed. If an interaction
crosses a finite-supercell boundary, explicitly include all physically distinct
periodic images or stop the import until their assignment is justified.

### Convert to GeoForm

```bash
python -m geoform_phonons.phonopy_bridge \
  --force-constants FORCE_CONSTANTS \
  --crystal my_primitive_crystal.json \
  --image-map my_supercell_images.json \
  --evidence-status calculated \
  --source 'DFPT calculation: software/version, settings, DOI or file hash' \
  --output geoform_model.json

python -m geoform_phonons --input geoform_model.json --output result.json
```

A correctly parsed IFC may still have a nonzero acoustic sum-rule residual,
unstable modes, incorrect origin ordering, incorrect supercell lattice/images,
or nonanalytic LO–TO splitting missing. The importer does not certify the source.
The main GeoForm exporter intentionally refuses ASR violations above its
published tolerance.

## Tests and independent reproduction

`tests/test_phonopy_bridge.py` uses an analytic 3×3×3 single-atom simple-cubic
supercell with onsite `6K` and six `−K` nearest-image tensor couplings. It
checks strict full and compact import, missing mappings, and equality with an
independently constructed analytical spring model.

`scripts/check_phonopy_reference.py` uses **installed third-party Phonopy** to
construct an independent `27×27×3×3` full-supercell force-constant matrix and
compares Phonopy's calculated modes at five wavevectors with GeoForm. This is a
controlled *synthetic* reference—not a published mineral dataset. It runs as an
optional dependency test in CI; the core package does not require Phonopy.

A genuine mineral benchmark requires: a published primitive structure,
calculation provenance and exact settings, verified supercell-atom mapping,
complete external IFC file, external Phonopy baseline, and independently
compared bands (and later thermodynamics) with discrepancies retained.
