"""Package only *verified computational* mineral results for a static GeoForm viewer.

Refuses incomplete result schemas and missing independent comparison receipts;
never silently promotes an illustrative model to a mineral calculation.
"""
import argparse
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def verify_results(results):
    """Return the catalog entries after validating evidence and cross-checks."""
    results=Path(results)
    catalog=json.loads((results/'catalog.json').read_text(encoding='utf-8'))
    receipts=json.loads((results/'benchmark_receipts.json').read_text(encoding='utf-8'))
    if catalog.get('schema_version')!='geoform.phonons.catalog/1' or receipts.get('schema_version')!='geoform.phonons.benchmark_receipts/1':
        raise ValueError('missing valid catalog/receipt schema')
    if {item.get('material') for item in receipts['results']} != {'halite','periclase'}:
        raise ValueError('missing verified halite/periclase receipts')
    if {item.get('slug') for item in catalog['entries']} != {'halite','periclase'}:
        raise ValueError('catalog must have exactly halite and periclase')
    for item in catalog['entries']:
        if item.get('evidence_status')!='calculated' or item['result_file'] != f'{item["slug"]}_phonons.json':
            raise ValueError('invalid catalog file or status')
        record=json.loads((results/item['result_file']).read_text(encoding='utf-8'))
        evidence=record.get('evidence',{})
        check=evidence.get('benchmark',{})
        if (record.get('schema_version')!='geoform.phonons.result/1' or
                evidence.get('declared_status')!='calculated' or
                check.get('independent_comparison')!='passed; computational, not experimental' or
                check.get('comparison_points',0)<12 or
                record.get('physics',{}).get('nonanalytic_correction')!= 'omitted; polar LO-TO splitting not reproduced'):
            raise ValueError('missing computational proof/limitation labels in output')
    return catalog['entries']


def create_bundle(results, output):
    results=Path(results)
    items=verify_results(results)
    out=Path(output)
    out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        z.write(ROOT/'web'/'phonon_explorer.html',arcname='phonon_explorer.html')
        for file in ('catalog.json','benchmark_receipts.json'):
            z.write(results/file,arcname=f'data/{file}')
        for item in items:
            z.write(results/item['result_file'],arcname='data/'+item['result_file'])
        z.write(ROOT/'docs'/'REAL_MINERAL_BENCHMARKS.md',arcname='SCIENCE_AND_PROVENANCE.md')
    return out


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=ROOT/'benchmarks'/'generated')
    p.add_argument('--output',type=Path,default=ROOT/'benchmarks'/'geoform_site_phonons.zip')
    args=p.parse_args(argv)
    output=create_bundle(args.results,args.output)
    print(f'Verified site bundle: {output}, {output.stat().st_size:,} bytes')

if __name__=='__main__':main()
