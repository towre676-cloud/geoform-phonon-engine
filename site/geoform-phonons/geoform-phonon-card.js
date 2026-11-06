/* GeoForm Site integration: evidence-safe progressive enhancement.
 * This is an API module for the ORIGINAL site's components, not a production
 * injector into the published Site. The Site editor must mount it explicitly.
 */
(function (root) {
  'use strict';
  const cases = Object.freeze({
    halite: Object.freeze({name: 'Halite', formula: 'NaCl', errorTHz: 7.744035919188264e-7}),
    periclase: Object.freeze({name: 'Periclase', formula: 'MgO', errorTHz: 2.1755864274553005e-6})
  });
  function cleanSlug(s) { return String(s||'').trim().toLowerCase().replace(/[^a-z0-9-]/g,''); }
  function hrefForMineral(slug, basePath) {
    const base = (typeof basePath === 'string' && /^\/[a-z0-9/_-]*\/$/.test(basePath)) ? basePath : '/phonons/';
    const s = cleanSlug(slug);
    return Object.prototype.hasOwnProperty.call(cases,s) ? base+'?mineral='+encodeURIComponent(s) : null;
  }
  function mountCard(container, opts) {
    if (!container || typeof container.replaceChildren !== 'function') throw new TypeError('DOM mount element is required');
    opts=opts||{};
    const slug=cleanSlug(opts.mineralSlug);
    const match=cases[slug];
    const section=document.createElement('section');
    section.setAttribute('aria-label','Phonon and lattice dynamics evidence');
    section.className='geoform-phonon-evidence';
    const title=document.createElement('h3');title.textContent='Vibrations inside the crystal';
    const badge=document.createElement('strong');
    badge.textContent=match?'CALCULATED • BENCHMARKED':'UNRESOLVED • NO CALCULATED DATA';
    const description=document.createElement('p');
    description.textContent=match?
      `${match.name} (${match.formula}): DFT-derived harmonic phonon result, independently compared to Phonopy at 16 q points. Not a measured Raman or IR spectrum.`:
      'No source-backed phonon force constants are available for this mineral. Crystal identity alone cannot establish its phonon frequencies.';
    section.append(title,badge,description);
    if(match){
      const a=document.createElement('a');a.href=hrefForMineral(slug,opts.basePath);
      a.textContent='Inspect phonon bands, vibrational DOS and heat capacity';
      section.append(a);
      const note=document.createElement('p');
      note.textContent='Long-range polar LO–TO corrections are excluded from this benchmark. Calculated ≠ experimentally measured.';
      section.append(note);
    }
    container.replaceChildren(section);
    return section;
  }
  root.GeoFormPhononLab=Object.freeze({mountCard,hrefForMineral,verifiedMinerals:Object.keys(cases)});
})(typeof window !== 'undefined' ? window : globalThis);
