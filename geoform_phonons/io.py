"""Versioned, provenance-aware JSON input/output for the GeoForm Mineral Lab."""
import json
import pathlib
import numpy as np
from .core import Crystal, ForceConstants, PhononInputError

SCHEMA = "geoform.phonons/1"
STATUSES = {"measured", "calculated", "reconstructed", "illustrative", "unresolved"}


def parse_model(record):
    """Read explicit IFCs, refusing invented constants or unlabelled units.

    provenance is a user declaration, not independent certification; output
    includes declared status and source text and makes that limitation clear.
    """
    if not isinstance(record, dict) or record.get("schema_version") != SCHEMA:
        raise PhononInputError(f"schema_version must equal {SCHEMA}")
    if record.get("model_kind") != "supplied_ifc":
        raise PhononInputError("model_kind must be supplied_ifc; a CIF alone cannot supply IFCs")
    meta = record.get("provenance")
    if (not isinstance(meta, dict) or meta.get("evidence_status") not in STATUSES
            or not isinstance(meta.get("source"), str) or not meta["source"].strip()):
        raise PhononInputError("provenance needs evidence_status and nonempty source")
    structure = record.get("crystal")
    if not isinstance(structure, dict):
        raise PhononInputError("crystal object required")
    crystal = Crystal(structure.get("lattice_angstrom"),
                      structure.get("fractional_positions"),
                      structure.get("masses_amu"), tuple(structure.get("symbols", [])))
    fc = record.get("force_constants")
    if not isinstance(fc, dict) or fc.get("units") != "eV/angstrom^2":
        raise PhononInputError("force_constants units must be eV/angstrom^2")
    blocks={}
    terms=fc.get("terms")
    if not isinstance(terms,list) or not terms:
        raise PhononInputError("force_constants.terms needs nonempty list")
    for term in terms:
        if not isinstance(term,dict) or not isinstance(term.get("R"),list) or len(term["R"])!=3:
            raise PhononInputError("each IFC term requires R length 3")
        key=(*term["R"],term.get("i"),term.get("j"))
        if key in blocks:
            raise PhononInputError("duplicate IFC key: use unique summed blocks")
        blocks[key]=term.get("phi")
    model=ForceConstants(crystal,blocks)
    return model, dict(meta)


def load_json(path):
    with open(path,encoding="utf-8") as f:
        return parse_model(json.load(f))


def model_to_record(model, provenance):
    if provenance.get("evidence_status") not in STATUSES or not provenance.get("source"):
        raise PhononInputError("explicit evidence_status and source required")
    return {"schema_version":SCHEMA,"model_kind":"supplied_ifc",
        "crystal":{"lattice_angstrom":model.crystal.lattice_angstrom.tolist(),
            "fractional_positions":model.crystal.fractional_positions.tolist(),
            "masses_amu":model.crystal.masses_amu.tolist(),"symbols":list(model.crystal.symbols)},
        "force_constants":{"units":"eV/angstrom^2","terms":[
            {"R":list(key[:3]),"i":key[3],"j":key[4],"phi":value.tolist()}
            for key,value in sorted(model.blocks.items())]},
        "provenance":dict(provenance)}


def export_result(model, provenance, path_vertices, *, points_per_segment=40,
                  mesh=(6,6,6), sigma_thz=0.1, temperature_k=300):
    from .core import frequency_grid, density_of_states, thermal_properties, projected_density_of_states, thermodynamic_curve
    if model.acoustic_sum_rule_error() > 1e-6:
        raise PhononInputError("ASR violation > 1e-6 eV/A^2; repair or justify before export")
    bands = model.band_path(path_vertices,points_per_segment)
    mesh_freqs=frequency_grid(model,mesh)
    unstable=int(np.count_nonzero(mesh_freqs < -1e-9))
    result={"schema_version":"geoform.phonons.result/1",
        "evidence":{"declared_status":provenance.get("evidence_status"),
                    "source":provenance.get("source"),
                    "verification":"user-supplied declaration; not independently authenticated",
                    "force_constants":"supplied, not inferred from CIF"},
        "crystal":{"symbols":list(model.crystal.symbols),
                   "volume_angstrom3":model.crystal.volume_angstrom3,
                   "n_atoms_primitive":model.crystal.n_atoms},
        "physics":{"method":"harmonic real-space IFC diagonalization",
                   "periodic":True,"imaginary_mode_convention":"negative signed THz",
                   "acoustic_sum_rule_max_abs_ev_a2":model.acoustic_sum_rule_error(),
                   "unstable_modes_on_shifted_mesh":unstable,
                   "q_mesh":list(mesh)},"dispersion":bands}
    if unstable==0 and np.all(mesh_freqs>0):
        x,y,by_atom=projected_density_of_states(model,mesh,sigma_thz,points=400)
        result["density_of_states"]={"frequency_thz":x.tolist(),"states_per_thz_cell":y.tolist(),
                                     "atom_projected_states_per_thz_cell":by_atom.tolist(),
                                     "projection_labels":list(model.crystal.symbols),
                                     "projection_basis":"mass-weighted eigenvector norm; 3 states/atom",
                                     "gaussian_sigma_thz":sigma_thz}
        result["thermodynamics"]=thermal_properties(mesh_freqs,temperature_k)
        result["thermodynamic_curve"]=thermodynamic_curve(mesh_freqs,
            sorted(set([0.0,50.0,100.0,200.0,float(temperature_k),500.0,800.0,1200.0])))
    else:
        result["density_of_states"]={"status":"withheld: unstable or zero mode"}
        result["thermodynamics"]={"status":"withheld: unstable or zero mode"}
    return result


def save_json(obj,path):
    pathlib.Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n",encoding="utf-8")
