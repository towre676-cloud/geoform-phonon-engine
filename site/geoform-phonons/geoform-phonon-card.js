/* GeoForm Site integration: evidence-safe progressive enhancement.
 * This is an API module for the ORIGINAL site's components, not a production
 * injector into the published Site. The Site editor must mount it explicitly.
 */
(function (root) {
  'use strict';
  const cases = Object.freeze({
    anatase: {name: 'Anatase', formula: 'TiO₂'},
    cassiterite: {name: 'Cassiterite', formula: 'SnO₂'},
    corundum: {name: 'Corundum', formula: 'Al₂O₃'},
    halite: {name: 'Halite', formula: 'NaCl'},
    periclase: {name: 'Periclase', formula: 'MgO'},
    stishovite: {name: 'Stishovite', formula: 'SiO₂', signedDispersionOnly: true},
    zincite: {name: 'Zincite', formula: 'ZnO', signedDispersionOnly: true}
  });
  function cleanSlug(s) { return String(s||'').trim().toLowerCase(); }
  function hrefForMineral(slug, basePath) {
    const base = (typeof basePath === 'string' && /^\/[a-z0-9/_-]*\/$/.test(basePath)) ? basePath : '/phonons/';
    const s = cleanSlug(slug);
    return Object.prototype.hasOwnProperty.call(cases,s) ? base+'?mineral='+encodeURIComponent(s) : null;
  }
  function mountCard(container, opts) {
    if (!container || typeof container.replaceChildren !== 'function') throw new TypeError('DOM mount element is required');
    opts=opts||{};
    const slug=cleanSlug(opts.mineralSlug);
    const match=Object.prototype.hasOwnProperty.call(cases,slug) ? cases[slug] : null;
    const section=document.createElement('section');
    section.setAttribute('aria-label','Phonon and lattice dynamics evidence');
    section.className='geoform-phonon-evidence';
    const title=document.createElement('h3');title.textContent='Vibrations inside the crystal';
    const badge=document.createElement('strong');
    badge.textContent=match?'CALCULATED • BENCHMARKED':'UNRESOLVED • NO CALCULATED DATA';
    const description=document.createElement('p');
    description.textContent=match?
      `${match.name} (${match.formula}): DFT-derived harmonic phonon result, independently compared to Phonopy at 16 q points. Experimental spectrum unavailable in this dataset.`:
      'No source-backed phonon force constants are available for this mineral. Crystal identity alone cannot establish its phonon frequencies.';
    section.append(title,badge,description);
    if(match){
      const a=document.createElement('a');a.href=hrefForMineral(slug,opts.basePath);
      a.textContent=match.signedDispersionOnly ? 'Inspect signed phonon bands' : 'Inspect phonon bands, vibrational DOS and heat capacity';
      section.append(a);
      const note=document.createElement('p');
      note.textContent=(match.signedDispersionOnly ?
        'An imaginary mode near Γ prevents equilibrium DOS and thermal export. ' : '') +
        'Long-range polar LO–TO corrections are excluded from this benchmark.';
      section.append(note);
    }
    container.replaceChildren(section);
    return section;
  }
  root.GeoFormPhononLab=Object.freeze({mountCard,hrefForMineral,verifiedMinerals:Object.keys(cases)});
})(typeof window !== 'undefined' ? window : globalThis);
