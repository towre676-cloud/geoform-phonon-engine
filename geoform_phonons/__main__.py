"""Command line: python -m geoform_phonons [--input model.json] [--output result.json]."""
import argparse
import json
from .core import PhononInputError
from .fixtures import simple_cubic
from .io import load_json, export_result, save_json, model_to_record


def main():
    p=argparse.ArgumentParser(description="GeoForm harmonic IFC phonon postprocessor")
    p.add_argument("--input",help="explicit structure + IFC + provenance JSON")
    p.add_argument("--output",default="phonon_results.json",help="output GeoForm result JSON")
    p.add_argument("--mesh",type=int,default=6,help="shifted q mesh N x N x N")
    p.add_argument("--temperature",type=float,default=300,help="temperature Kelvin")
    p.add_argument("--points",type=int,default=40,help="points per segment")
    p.add_argument("--demo-input",help="also write illustrative cubic IFC input JSON here")
    args=p.parse_args()
    if args.input:
        model, provenance=load_json(args.input)
    else:
        model=simple_cubic()
        provenance={"evidence_status":"illustrative", "source":"GeoForm synthetic simple-cubic spring fixture"}
    if args.demo_input:
        save_json(model_to_record(model,provenance),args.demo_input)
    high_symmetry=[[0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.5]]
    result=export_result(model,provenance,high_symmetry,
                         points_per_segment=args.points,mesh=(args.mesh,)*3,
                         temperature_k=args.temperature)
    save_json(result,args.output)
    print(json.dumps({"output":args.output,"n_atoms":model.crystal.n_atoms,
          "declared_evidence":provenance["evidence_status"],
          "unstable_mesh_modes":result["physics"]["unstable_modes_on_shifted_mesh"]}))


if __name__=="__main__":
    main()
