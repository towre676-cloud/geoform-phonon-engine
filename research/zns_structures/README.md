# ZnS crystallographic input foundation

Proposed for GeoForm Phonon Engine [issue #7](https://github.com/towre676-cloud/geoform-phonon-engine/issues/7). The package contains three validated geometry fragments and one diffraction evidence record. The four-phase study remains blocked pending a pressure-tagged Cmcm atomic model.

| Phase | Evidence and disposition | Conventional reference cell | Conventional / primitive ZnS formula units |
|---|---|---|---|
| Zincblende, F-43m (216) | COD 9000107, AMCSD 0000110. Synthetic ZnS. Measured cell edge with symmetry-fixed Zn 4a and S 4c. The ordered model assigns full sites, since the CIF omits occupancy. | 5.4093(2) Å cubic. Skinner cites the pure-ZnS value at 25 °C. Numerical pressure unreported. | 4 / 1 |
| Wurtzite-2H, P6₃mc (186) | COD 1100044, checked against COD 1101051 and the original 295 K powder-neutron Rietveld refinement. Zn and S on 2b, S (u=0.3748(2)). Both sites have unit occupancy in the CIF. | (a=b=3.8227(1)), (c=6.2607(1)) Å, (gamma=120^circ). Numerical pressure unreported. | 2 / 2 |
| Rocksalt, Fm-3m (225) | Reconstructed B1 reference geometry. Zn 4a and S 4b are fixed by the ordered B1 model. Unit occupancies are a model assumption. The lattice is reconstructed from a fit to measured synthetic-ZnS data. | (a=b=capprox4.7288) Å at chosen 30 GPa, room-temperature EOS. | 4 / 1 |
| High-pressure Cmcm (63) | Diffraction indexing, peak splitting, and lattice refinement. The reported site model has a 4c orbit per species with generic ((0,y,1/4)) positions. Neither internal (y) nor the occupancies are established for a coordinate file. | (a=4.438(1)), (b=4.486(1)), (c=4.412(1)) Å at 96 GPa, room temperature, natural sample 5. | 4 / 2 conditional on the reported site model |

### Source checks

- **Zincblende.** The current COD record supersedes the Avogadro copy for this package. Both expand to the same ordered structure. The 1961 paper gives the lattice metric and (Z=4) on p. 1399, citing earlier measurements. Its scope is cell-edge measurement. The fractional coordinates in the CIF are the fixed special positions of the B3 model. The model's ideal stoichiometry does not imply a defect or occupancy measurement. [COD 9000107](https://www.crystallography.net/cod/9000107.cif@291269) · [Skinner, 1961](https://msaweb.org/AmMin/AM46/AM46_1399.pdf)
- **Wurtzite.** The starting candidate COD 9008878 is attributed to Wyckoff's *Crystal Structures* (1963), a compilation. It records (a=3.811) Å, (c=6.234) Å and (u=0.385), and omits the occupancy, pressure and temperature fields. It is kept for provenance. The Kisi–Elcombe refinement gives (u=0.3748(2)) at 295 K. COD 1101051 is marked as a duplicate of COD 1100044 and is a second transcription check from the same experiment. The geometry uses exact (1/3,2/3) for the symmetry-fixed basal coordinates and retains the refined (u). [Wyckoff COD record](https://www.crystallography.net/cod/9008878.html) · [Refinement DOI](https://doi.org/10.1107/S0108270189004269) · [COD 1100044](https://www.crystallography.net/cod/1100044.html) · [COD 1101051](https://www.crystallography.net/cod/1101051.html)
- **Rocksalt.** Desgreniers, Beaulieu and Lepage (2000), Table III and Fig. 4, give the constrained synthetic-rocksalt Birch–Murnaghan fit (V_0=31.6) Å³ per ZnS, (B_0=117.6) GPa, (B'_0=4). The reference volume is a fit parameter. At 30 GPa the inverted fit gives 26.43595 Å³ per ZnS and (a=(4V)^{1/3}). The pressure lies inside the authors' common single-phase interval of about 20–65 GPa. The fit, calculated lattice, and prototype coordinates are separately identified. No experimental per-atom refinement, CIF or diffraction point at precisely 30 GPa is claimed. Use the reported transition pressures as external benchmarks. [Desgreniers et al., Table III and §III B](https://doi.org/10.1103/PhysRevB.61.8726) · [B1 prototype](https://aflowlib.duke.edu/p/AB_cF8_225_a_b-001/)
- **Cmcm.** The paper reports extra reflections from roughly 65 GPa and resolves the orthorhombic lattice at the high-pressure end of the experiment. Fig. 8 gives the 96 GPa cell. Section III E explicitly leaves the atomic positions undetermined. Under the proposed two-4c model the cell contains four ZnS formula units (21.9595 Å³ per ZnS) and its C-centered primitive cell would contain two. Missing are species-assigned (y_{\mathrm{Zn}}) and (y_{\mathrm S}), uncertainties, site occupancies or explicit occupancy constraints, and an auditable refinement or intensity dataset tied to the same pressure, temperature, and setting. A later LDA study provides a theoretical lead, but a reusable pressure-tagged coordinate record could not be verified from the accessible publication material. The package does not substitute rocksalt-limit coordinates or coordinates from another compound. [Desgreniers et al., Fig. 8 and §III E](https://doi.org/10.1103/PhysRevB.61.8726) · [López-Solano et al., 2003](https://doi.org/10.1002/pssb.200301601)

### Files proposed for repository review

Copy this directory to `research/zns_structures/` after approval.

| File | Role |
|---|---|
| `manifest.json` | Phase evidence and conditions, explicit occupancy provenance, source and revision identifiers, original SHA-256 hashes, and the Cmcm release gate. |
| `originals/*.cif` | Five untouched CIF files, comprising four direct COD records and the historical Avogadro copy of COD 9000107. |
| `crystals/*.json` | Primitive geometry with the engine's `crystal` keys. Lattice vectors are rows in Å, positions are fractional, and element masses in amu are conventional natural-composition inputs. |
| `validate_structures.py`, `test_structures.py`, `requirements.txt` | Offline replay and focused validation with ASE and spglib. |

COD distributes its data under [CC0](https://cod.psdi.ac.uk/) and requests credit to the original structural authors. The Avogadro CIF carries the original COD public-domain and AMCSD attribution notices. Those notices are retained byte-for-byte. Article files are linked. Skinner's publisher PDF was fingerprinted, and no raw APS PDF or original high-pressure CIF was acquired for archival hashing. The APS and IUCr articles are referenced for attribution and extracted facts. Article redistribution rights were not established.

### Replay

From the extracted directory, using Python 3.11 or later:

```sh
python -m pip install -r requirements.txt
python validate_structures.py
python -m unittest -v
```

Use `python validate_structures.py --write` to regenerate the three fragments. Validation compares original bytes and CIF fields, expands the supplied symmetry operations, independently redetects symmetry with spglib, verifies Wyckoff multiplicities, checks primitive and conventional normalization, and exercises negative cases for changed bytes, invalid occupancies, stale exports and unauthorized Cmcm coordinate files. The three fragments were also passed to the current engine's `Crystal` constructor. Its phonon-input parser refused the fragments when IFCs were absent. Any later IFCs must match the primitive cell setting and atom order exactly.

The research package contains crystallographic inputs and experimental comparison evidence. No DFT, QHA, force-constant generation, transition prediction, or public deployment was performed.
