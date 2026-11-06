# GeoForm Phonon Engine

A research-grade *development target* beginning with a verified, minimal
one-dimensional harmonic diatomic-chain reference model. The repository is
**not yet** a mineral-specific phonon prediction system.

## Reproduce

Run with Python 3.11+ (standard library only):

```sh
python -m unittest discover -s tests -v
python phonon_chain.py --mass-a 1 --mass-b 2 --spring 3 --cell-length 1
```

The implementation produces acoustic and optical angular-frequency branches
for alternating point masses with equal nearest-neighbour harmonic springs.
It checks the zone-centre acoustic zero, the analytic zone-boundary limits,
mass exchange, reciprocal wavevector symmetry, and long-wave sound speed.

## Evidence discipline

The calculated values are **ideal-model predictions**, not measured spectra,
real CIF structures, Raman or infrared peak assignments, or electron–phonon
coupling constants. Frequencies are angular frequencies in rad/s; conversion to
Hz requires division by 2π. Input parameters must carry consistent SI units.

## Roadmap

1. Generalize to N-atom unit cells and Hermitian dynamical matrices.
2. Add 3-D lattice vectors, masses, force-constant tensor and acoustic sum rule.
3. Cross-check independently against established phonopy reference fixtures.
4. Parse CIF structures as identities, **not** as force constants; require
   measured or calculated interatomic force constants separately.
5. Add symmetry labels and separately sourced Raman/IR activities and spectra.
6. Integrate only validated outputs into GeoForm with measured/calculated/
   reconstructed/unresolved evidence labels.

Historical motivation: November 2023 electron–phonon, anharmonicity, and
Bloch-symmetry conversations recovered in the research excavation. Their
earlier physical and historical claims are not imported as authoritative.
