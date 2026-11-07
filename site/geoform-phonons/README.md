# GeoForm Mineral Lab phonon handoff

Serve this folder at `/phonons/`. The explorer reads its catalog and datasets
from relative paths. Test locally with `python -m http.server 8765 --directory site/geoform-phonons`.

- Use the `mineral` parameter with `anatase`, `cassiterite`, `corundum`,
  `halite`, `periclase`, `stishovite` or `zincite` to open a mineral's calculation.
- Use `?phase=cubic-catio3` for the ideal cubic model or `?phase=hexagonal-mgb2`
  for hexagonal magnesium diboride. Both are computational-phase entries.
  A perovskite mineral request remains unresolved.
- In the Mineral Lab source, add the viewer to navigation or import
  `geoform-phonon-card.js` and call
  `GeoFormPhononLab.mountCard(element, {mineralSlug, basePath: '/phonons/'})`.
  Mount a card only where the exact phase identity is appropriate. Unsupported
  slugs receive an unresolved card with no calculated link.
- Review the raw-input hashes and checks in `data/benchmark_receipts.json`.
  `site-handoff.json` pins the assets. GitHub publication and separately hosted
  Mineral Lab deployment are distinct operations.

All nine datasets are **calculated**, with 16 q points checked against Phonopy.
All omit polar LO–TO corrections. Six cases include sampled harmonic DOS,
atom-projected DOS and thermal response per primitive cell. Cubic CaTiO₃,
zincite and stishovite provide signed dispersion with DOS and equilibrium thermal response withheld.
The zincite and stishovite soft modes are visible on the dispersion path
and missed by the shifted integration mesh. Computational cells can differ
from the Mineral Lab's structural records even when the mineral phase matches.
