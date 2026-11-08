# NIMS PhononDB handoff to GeoForm Mineral Lab

The engine's nine pinned benchmarks and the NIMS phase catalog use separate evidence and asset paths. The combined export runs both validators. The existing nine-result build also works without the PhononDB arguments.

The public Mineral Lab already has a NIMS browser under `#phonons` and serves `phonons/nims/index.bin`, `source-index.bin` and one gzip JSON payload per accepted MDR ID. The export targets that wire format. The hosted patch preserves the existing signed dispersion, q-point inspection, total DOS, constant-volume heat-capacity plot and thermal readouts. It adds bounded 24-row pagination, calculated-evidence labeling and two-way navigation for independently reviewed structure matches. The NIMS source index is loaded on demand. The browser opens one local numerical payload at a time and verifies its raw SHA-256 after decompression.

## Inputs and export

Run from the Phonon Engine checkout with Python 3.12 and `python -m pip install -e '.[phonondb]'`. The recorded replay used NumPy 2.3.5, SciPy 1.18.1, PyYAML 6.0.3, Phonopy 4.9.0, spglib 2.8.0 and symfc 1.8.0. Use the pinned source index, both sets of original archive pins and the combined selection in [PHONONDB.md](PHONONDB.md). The 578 required archive IDs, exact URLs, byte lengths and SHA-256 values are recorded in the two `benchmarks/phonondb_*_pins.json` files.

Reuse a valid receipt and result together. With a previous review bundle, put the raw result at `work/verified/results/<record-id>.json` and the corresponding row from `receipts.json` at `work/verified/records/<record-id>.json`. The existing verifier accepts the cache only if the archive, source row, result hash and settings match. Re-fetch any missing source archive using its exact batch pin. The publisher's numerical comparison must be replayed if its receipt or product is missing or no longer matches. The 36 expected large-cell exclusions remain in the aggregate receipts.

Run record-rights checks for the combined selection. Preserve successful checks when retrying failed requests and use the resulting map when assembling `work/verified/receipts.json` and `catalog.json` with the existing `verify` command. Its documented nonzero status for the 36 exclusions is expected. Fetch the collection's own rights statement from [NIMS MDR](https://doi.org/10.48505/nims.4197).

```sh
curl -fL https://doi.org/10.48505/nims.4197 -o work/collection.html

python scripts/build_site_bundle.py \
  --results site/geoform-phonons/data \
  --phonondb-replay work/verified \
  --phonondb-cache work/cache \
  --phonondb-index work/index.md \
  --phonondb-collection-page work/collection.html \
  --output work/geoform-phonons-and-nims.zip
```

The default full export requires each distributed record's CC BY 4.0 receipt. An explicit `--allow-collection-license` permits documented NIMS collection-level coverage when a record has no individual receipt or its check hit a transient retrieval error. An observed conflicting license or a record identity error still blocks that record. When any record would be withheld, the full export fails before replacing an existing ZIP. Use `--allow-partial-phonondb` to generate a review package with the omitted IDs and reasons. Keep partial exports off production.

The exporter checks the pinned index and archive bytes, both archived YAML hashes, source row, comparison pipeline digests, catalog inventory and exact result bytes. It compares the primitive lattice, ordered species, masses and positions to the archived source, and recomputes the space-group number from the original cell. It checks the dimensions and values of each dispersion, numerical comparison bounds, DOS mode-count normalization, thermal units and field names, and thermodynamic consistency. It preserves signed values. Where imaginary modes survive the verifier's Γ acoustic rule, only dispersion is exported. This is a provenance and consistency check using recorded numerical comparisons on common forces. The underlying DFT forces and experimental behavior are outside its scope. NAC and polar LO–TO splitting are omitted.

## Reviewed structural links

| Mineral profile | NIMS phase | Independent structural reference | Checked identity |
|---|---|---|---|
| Halite | `gf06g7088`, MP-22862 | [AFLOW rock-salt prototype with crystallographic references](https://aflow.dev.materials.duke.edu/p/AB_cF8_225_a_b-001/) | NaCl, Fm-3m, cF8, full ordered rock-salt motif |
| Periclase | `bk128j14x`, MP-1265 | [USGS Bulletin 1248, X-Ray Crystallographic Data of Minerals](https://pubs.usgs.gov/bul/1248/report.pdf) | MgO, Fm-3m, cF8, full ordered rock-salt motif |

For these two pinned source IDs the exporter requires an exact ordered conventional-cell motif up to an origin shift. It emits a separate mineral-link table and retains the NIMS result's unassigned mineral identity. The links describe the structural analogue of an ideal computed cell. All remaining phases stay in the computational browser with no mineral profile link.

## Apply to a staged copy of the hosted Mineral Lab

The patch was prepared against the served `phonon-ui.mjs` and `phonon-nims.mjs` on November 8, 2026. Re-read the deployed Site and verify the current checkout before applying it. From the folder containing those files and `app.mjs`:

```sh
patch --dry-run -p1 < /path/to/mineral-lab.patch
patch -p1 < /path/to/mineral-lab.patch
unzip -o /path/to/geoform-phonons-and-nims.zip 'phonons/nims/*' 'phonon-nims-links.mjs'
unzip -o /path/to/geoform-phonons-and-nims.zip 'PHONONDB_EXPORT_REPORT.json'
node /path/to/check.mjs
```

Compare the existing and proposed indexes by MDR ID and numerical result hash. Review any loss or byte change before applying the overlay. The exported `data/` and `phonon_explorer.html` are the legacy engine viewer assets. The hosted overlay uses the three named ZIP entries above. Keep the export report with the release evidence and remove it from the public document root after testing. Keep receipt timing and retrieval logs private to the review.

After approval, merge the engine stack against its reviewed heads, with the exporter following the PhononDB changes. Treat the hosted Site change as a separate operation. Apply the overlay and patch to the current Site source, run its build and smoke checks, save the exact built version, then publish to the current public audience. Verify the NIMS count and direct links at `#phonons?mineral=nims-gf06g7088` and `#phonons?mineral=nims-bk128j14x`, and check Halite, Periclase and an unrelated mineral profile. The hosted benchmark and measured-spectrum routes should still resolve.
