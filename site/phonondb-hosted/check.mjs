// Run from a staged copy of the Mineral Lab: node /path/to/check.mjs
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
import {resolve} from 'node:path';
import {gunzipSync} from 'node:zlib';

const hostModule = name => import(pathToFileURL(resolve(name)));
const readJSON = async name => JSON.parse(gunzipSync(await readFile(name)));
const dir = 'phonons/nims/';
const index = await readJSON(dir + 'index.bin');
const source = await readJSON(dir + 'source-index.bin');
const {findNimsPhases, readNimsPhase, mountNimsPhases} = await hostModule('phonon-nims.mjs');
const {phononProfileHTML} = await hostModule('phonon-ui.mjs');
const {nimsMineralLinks} = await hostModule('phonon-nims-links.mjs');
const report = JSON.parse(await readFile('PHONONDB_EXPORT_REPORT.json'));

assert.equal(index.entries.length, report.exported);
assert.equal(source.entries.length, index.indexRecords);
assert.equal(source.sourceIndexSha256, index.sourceIndexSha256);
assert.equal(new Set(index.entries.map(row => row.id)).size, index.entries.length);
assert.ok(findNimsPhases(index.entries, 'SiO₂').length > 0);
assert.ok(findNimsPhases(index.entries, '', 'thermal').every(row => row.dos));
assert.ok(findNimsPhases(index.entries, '', 'signed', 'Cubic')
  .every(row => !row.dos && row.number >= 195));

const readAsset = async url => new Response(await readFile(url));
globalThis.fetch = readAsset;
for (const entry of index.entries) {
  const {result} = await readNimsPhase(entry);
  assert.equal(result.crystal.mineral_name, null);
  assert.equal(result.evidence.experimental_validation, false);
  assert.equal(Boolean(result.thermodynamic_curve), entry.dos);
}
await assert.rejects(readNimsPhase({...index.entries[0], sha256: '0'.repeat(64)}), /checksum/);
for (const [name, href] of nimsMineralLinks) {
  assert.ok(phononProfileHTML({name}).includes(href));
  assert.ok(!phononProfileHTML({name, predictionsWithheld: true}).includes(href));
  const id = href.split('nims-')[1];
  assert.equal(index.entries.find(entry => entry.id === id)?.mineralStructure?.name, name);
}
assert.ok(!phononProfileHTML({name:'Quartz'}).includes('Open NIMS calculation'));

// Exercise bounded rendering through the actual hosted UI's event handlers.
const nodes = new Map();
const node = id => {
  if (!nodes.has(id)) {
    nodes.set(id, {value: '', innerHTML: '', textContent: '', disabled: false,
      replaceChildren() { this.innerHTML = ''; }});
  }
  return nodes.get(id);
};
const host = {isConnected: true, querySelector: id => node(id)};
const wait = () => new Promise(done => setTimeout(done, 30));
mountNimsPhases(host, {escape: String});
for (let i = 0; i < 50 && !node('#nims-list').innerHTML; i++) await wait();
const items = () => [...node('#nims-list').innerHTML.matchAll(/data-id="([^"]+)"/g)].map(m => m[1]);
const first = items();
assert.equal(first.length, 24);
assert.equal(node('#nims-prev').disabled, true);
node('#nims-next').onclick();
assert.equal(items().length, 24);
assert.ok(items().every(id => !first.includes(id)));
node('#nims-prev').onclick();
assert.deepEqual(items(), first);
node('#nims-search').value = index.entries[0].id;
node('#nims-search').oninput();
assert.deepEqual(items(), [index.entries[0].id]);
assert.equal(node('#nims-next').disabled, true);
console.log(`Hosted loader validated ${index.entries.length} calculation assets and ${nimsMineralLinks.size} mineral links; search and bounded pagination passed.`);
