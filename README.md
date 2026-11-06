# GeoForm Phonon Engine

**A traceable crystal-phonon postprocessor, not a mineral-frequency guessing service.**

[![Tests](https://github.com/towre676-cloud/geoform-phonon-engine/actions/workflows/tests.yml/badge.svg)](https://github.com/towre676-cloud/geoform-phonon-engine/actions/workflows/tests.yml)

The original `phonon_chain.py` remains a standard-library-only **exact diatomic-chain benchmark**. The new `geoform_phonons` package generalizes it to periodic 3-D crystals with an arbitrary number of atoms in the primitive cell **when an explicit real-space interatomic force-constant (IFC) model is supplied**. No force constants are inferred from a formula or CIF.

## Run now

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python -m geoform_phonons --demo-input examples/simple_cubic_ifc.json --output examples/simple_cubic_result.json --mesh 6
```

Open `web/phonon_explorer.html` and choose `examples/simple_cubic_result.json` to explore the **illustrative synthetic** bands, DOS and thermal values without web dependencies. `python phonon_chain.py` is unchanged and independently executable.

To calculate a crystal for which you have proper IFCs:

```bash
python -m geoform_phonons --input my_ifcs.json --output my_phonons.json --mesh 10 --temperature 300
```

Schema: `geoform.phonons/1`, implemented in `geoform_phonons/io.py`. Its JSON requires a nonsingular 3×3 lattice (Å, row vectors), fractional atomic positions, positive atomic masses (amu), complete paired IFC blocks in eV/Å², and declared source/evidence status. The returned `geoform.phonons.result/1` JSON includes provenance, periodic dispersion (signed THz), stable-branch Gaussian DOS, shifted-mesh zero-point energy, Helmholtz free energy, entropy and constant-volume heat capacity. See [scientific definitions and limitations](docs/SCIENCE.md).

## Implemented physics

| Capability | Current state | Evidence / limitations |
| --- | --- | --- |
| Analytical diatomic chain | Tested reference | 1-D point masses only; output angular rad/s |
| General 3-D harmonic bands | Implemented | Arbitrary primitive basis, supplied real IFCs, THz cycles/s |
| Dynamical Hermiticity | Enforced | Reciprocal IFC pairs checked; complex q |
| Translational acoustic sum rule | Diagnosed, optional onsite correction | No rotational ASR, correction changes the provided IFCs |
| Imaginary modes | Explicitly negative signed frequencies | DOS and equilibrium thermodynamics withheld if unstable |
| Gaussian DOS | Implemented | Not neutron or Raman scattering intensity |
| Harmonic F, U, S, Cv, ZPE | Implemented | Uniform shifted q mesh; no thermal expansion or anharmonicity |
| Relative Gamma IR activity | Implemented with supplied Born charges | Neutrality checked; not absolute absorption |
| Relative polarized Raman | Implemented with supplied polarizability derivatives | Not a spectrum or calibrated intensity |
| Eliashberg λ and ω-log | Integrate supplied α²F only | Does **not** calculate electron–phonon couplings |
| Crystal structure/CIF | Metadata may be converted to JSON | CIF alone cannot produce IFCs; no built-in CIF parser |
| Validated mineral predictions | **Not established** | No independent frozen measured-mineral holdout |

The browser explorer is a *read-only viewer* of results; it performs no scientific inference. Export declares the input's evidence status **as supplied by the user**, never as independent verification. Calculated spring fixtures are explicitly `illustrative` and must not be attached to real mineral records as measured spectra.

## Near-term research frontier

1. Build a curated and licensed set of actual **CIF + DFPT or finite-displacement IFC + Raman/IR measured spectra**, with versioned provenance and matching crystallographic identities.
2. Compare against an independently run **Phonopy** benchmark with 3-D reference force constants; document conversion and q-phase convention.
3. Provide validated CIF/site-mapping, unit conversion and supercell-to-primitive IFC import. Disorder must not be silently resolved.
4. Implement non-analytic LO–TO corrections only with explicit Born tensors, dielectric tensors and cell-volume factors.
5. Freeze a held-out mineral benchmark, comparing composition-only, symmetry-only, and structure+IFC models with error bars, incorrect phase orderings and abstentions.
6. Integrate `geoform.phonons.result/1` into the main GeoForm Mineral Lab as a separately sourced tool with *Measured / Calculated / Reconstructed / Illustrative / Unresolved* labels.

**Research context:** This engine is intended to work alongside GeoForm's structure, diffraction, thermodynamic phase, and spectroscopy efforts. The current repository contains a functioning harmonic lattice calculation layer, not a completed mineral-property predictor.

## Scientific references

- [Phonopy: formulation of the dynamical matrix](https://phonopy.github.io/phonopy/formulation.html)
- [Phonopy: harmonic and quasi-harmonic calculations](https://phonopy.github.io/phonopy/)
- [Quantum ESPRESSO PHonon user guide](https://www.quantum-espresso.org/Doc/ph_user_guide/)

## Phonopy full/compact force-constant interchange (Phase 2)

`python -m geoform_phonons.phonopy_bridge --help` converts Phonopy text
`FORCE_CONSTANTS` to GeoForm input **only with an explicit primitive-crystal
record and atom-by-atom supercell-image map**. See
[Phonopy interoperability](docs/PHONOPY_INTEROPERABILITY.md) for exact
formats, CLI examples and limitations. `FORCE_CONSTANTS` alone does not encode
image translations; uninformed conversions are refused.

The separately runnable `scripts/check_phonopy_reference.py` uses installed
Phonopy as an independent harmonic-dispersion baseline, and is verified in CI.

