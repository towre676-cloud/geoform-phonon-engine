# GeoForm live Site integration handoff — phonon laboratory

**Publishing boundary:** this directory is a verified, deployable patch for the *existing* `geoform-mineral-lab.towre676.chatgpt.site` ChatGPT Site. GitHub merging DOES NOT redeploy that Site. Only the authenticated Site owner/editor in the Site editing workflow can publish the existing URL.

## Apply to the EXISTING Site — do not create a replacement

1. Preserve all existing GeoForm routes, mineral records, specimen photos, CIF provenance, study flows, styling, and deployed URL. Add these files to the existing site's public/static directory as `phonons/index.html`, `phonons/geoform-phonon-card.js`, `phonons/data/*.json`, and `phonons/site-handoff.json`.
2. Make the site router serve `phonons/index.html` at `/phonons/` (and handle nested route fallback if applicable). It reads `data/catalog.json` and the two local JSON datasets from **same-origin** paths. The viewer supports `?mineral=halite` and `?mineral=periclase` for direct entry from mineral records.
3. Add a **Phonons** link in the main navigation and **Phonons & thermal response** in The Science page pointing to `/phonons/`; do not replace existing science routes.
4. In the existing Halite (`/minerals/halite`) and Periclase (`/minerals/periclase`) detail components, mount a compact `Vibrations inside the crystal` evidence card after `Inside the crystal.` (or equivalent structure/response section). Import `phonons/geoform-phonon-card.js` once and call `GeoFormPhononLab.mountCard(container,{mineralSlug:'halite'|'periclase',basePath:'/phonons/'})` on the component-owned DOM mount. On other mineral records, display unresolved state or omit the card; **do not extrapolate** from formula or symmetry.
5. Publish an update to the **same Site URL** through the owner/editor deployment flow; check `/phonons/`, both deep links, both mineral pages, mobile layout, CSV downloads, and The Science navigation in a signed-out browser. Use browser network tools to confirm all four JSON files are served with 200 and correct JSON content types.

## Evidence integrity

Only Halite (NaCl) and Periclase (MgO) are validated *calculated* datasets, each checked at 16 q-points against Phonopy from official pinned VASP `FORCE_SETS`; see the CI link in `site-handoff.json`. Signed imaginary frequencies must remain visible. The results omit non-analytic long-range LO–TO corrections, so do not identify them with measured Raman/IR spectra. DOS and harmonic heat capacity are calculated from the supplied IFCs and are not experimental data. The CIF identities currently in GeoForm's mineral records must not be assumed to equal the exact computational cells; show this caveat alongside data provenance.

## Security and implementation notes

The viewer only reads local JSON, does not run remote code or store private material, uses DOM `textContent` for data labels, and permits manual import of JSON with schema validation. Configure CSP and routing consistently with the existing Site. Do not inject DOM into React-owned elements outside the existing component lifecycle; the evidence-card module is intended for an explicit mounting point. Do not make the site copy global `?catalog` cross-origin trust without reviewing CORS and data integrity. Do not create a second public GeoForm site when updating this one.

## Reproduce

```bash
python -m unittest discover -s tests -v
python -m http.server 8765 --directory site/geoform-phonons
# http://localhost:8765/?mineral=halite
```

## Live status

`NOT DEPLOYED BY A GITHUB MERGE` until separately published via the existing ChatGPT Site editor. The Site's deployment state must be verified and updated here after publication.
