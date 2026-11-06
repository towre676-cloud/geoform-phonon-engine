"""Strict bridge from Phonopy FORCE_CONSTANTS text to GeoForm IFC terms.

Phonopy defines full (Ns,Ns,3,3) or compact (Np,Ns,3,3) tensors.
A tensor does NOT contain lattice-image labels. To map it to a periodic
real-space model, the caller must supply explicit supercell image labels,
origin row mapping, and length/energy units. This intentionally refuses
unmapped FORCE_CONSTANTS files instead of inventing periodic interactions.

See: https://phonopy.github.io/phonopy/input-files.html
"""
from __future__ import annotations
import io
import json
from pathlib import Path
import numpy as np
from .core import Crystal, ForceConstants, PhononInputError


def parse_force_constants_text(source):
    """Parse Phonopy text FORCE_CONSTANTS (both full and compact layouts).

    Comments/blanks are ignored; block label lines may contain arbitrary
    text, as specified by Phonopy. Rejects missing and surplus blocks.
    """
    if hasattr(source, 'read'):
        lines = source.read().splitlines()
    elif isinstance(source, (str, Path)):
        lines = Path(source).read_text(encoding='utf-8').splitlines()
    else:
        raise PhononInputError('FORCE_CONSTANTS source must be file path or open text stream')
    lines = [x.strip() for x in lines if x.strip() and not x.lstrip().startswith('#')]
    if not lines:
        raise PhononInputError('empty FORCE_CONSTANTS')
    try:
        dims = [int(v) for v in lines[0].split()]
    except (TypeError, ValueError) as exc:
        raise PhononInputError('FORCE_CONSTANTS first line must contain one or two integer dimensions') from exc
    if len(dims) == 1:
        dims = [dims[0], dims[0]]
    if len(dims) != 2 or not (0 < dims[0] <= dims[1]):
        raise PhononInputError('invalid compact/full FORCE_CONSTANTS dimensions')
    rows, cols = dims
    if len(lines) != 1 + 4*rows*cols:
        raise PhononInputError('FORCE_CONSTANTS requires exactly four nonblank lines per matrix block')
    tensor = np.empty((rows, cols, 3, 3), dtype=float)
    for index in range(rows*cols):
        start = 1 + 4*index
        for a in range(3):
            try:
                vals = [float(x) for x in lines[start + a + 1].split()]
            except ValueError as exc:
                raise PhononInputError('non-numeric force-constant matrix row') from exc
            if len(vals) != 3 or not np.isfinite(vals).all():
                raise PhononInputError('each force-constant matrix row needs three finite entries')
            tensor[index//cols, index%cols, a] = vals
    return tensor


def read_image_map(source):
    """Read required explicit image-map JSON, NOT inferred from crystal geometry.

    Required keys:
      supercell_images: length Ns [{"primitive_index":j,"translation":[Rx,Ry,Rz]}]
      origin_rows: length Np of supercell indices corresponding to row i
      source: citation/description of IFC provenance

    origin_rows permits full FORCE_CONSTANTS reduced to one representative
    primitive row per atom. For compact tensors origin_rows gives the
    positions of compact row representatives in the supercell image list.
    """
    if isinstance(source, dict):
        record = source
    elif hasattr(source, 'read'):
        record = json.load(source)
    else:
        record = json.loads(Path(source).read_text(encoding='utf-8'))
    if not isinstance(record, dict):
        raise PhononInputError('image map needs JSON object')
    return record


def mapped_force_constants(crystal: Crystal, tensor, mapping, *, tolerance=1e-8):
    """Construct primitive IFC model with explicit supercell atom provenance.

    Each supercell column represents a specific primitive atom and an integer
    lattice translation. A full tensor uses `origin_rows` to select the
    corresponding Np rows. A compact tensor uses the Np rows as given.
    No symmetry completion/averaging of missing translations is performed.
    """
    f = np.asarray(tensor, dtype=float)
    n = crystal.n_atoms
    if (f.ndim != 4 or f.shape[2:] != (3,3) or f.shape[0] not in (n,f.shape[1])
            or not np.isfinite(f).all()):
        raise PhononInputError('expected compact Np x Ns or full Ns x Ns finite IFC tensor')
    ns = f.shape[1]
    images = mapping.get('supercell_images')
    origins = mapping.get('origin_rows')
    if not isinstance(images, list) or len(images) != ns:
        raise PhononInputError('explicit supercell_images required for each supercell atom')
    if not isinstance(origins, list) or len(origins) != n or len(set(map(str,origins))) != n:
        raise PhononInputError('one distinct origin row per primitive atom required')
    normalized = []
    seen = set()
    for item in images:
        if not isinstance(item, dict):
            raise PhononInputError('each image must have primitive_index and translation')
        j = item.get('primitive_index')
        R = item.get('translation')
        if (type(j) is not int or not 0 <= j < n or not isinstance(R, list)
                or len(R) != 3 or any(type(x) is not int for x in R)):
            raise PhononInputError('invalid image primitive index or integer translation')
        key = (j, *R)
        if key in seen:
            raise PhononInputError('duplicate image (atom,translation) label')
        seen.add(key)
        normalized.append((j, tuple(R)))
    for i, s in enumerate(origins):
        if type(s) is not int or not 0 <= s < ns or normalized[s][0] != i:
            raise PhononInputError('origin_rows must map each primitive index to a supercell representative')
    if f.shape[0] != n and f.shape[0] != ns:
        raise PhononInputError('incompatible IFC tensor rows')
    if f.shape[0] == ns and ns != n:
        selected = f[origins, :, :, :]
    else:
        selected = f
    blocks = {}
    for i, origin in enumerate(origins):
        origin_translation = normalized[origin][1]
        for s, (j, col_translation) in enumerate(normalized):
            R = tuple(col_translation[k]-origin_translation[k] for k in range(3))
            blocks[(*R,i,j)] = selected[i,s]
    return ForceConstants(crystal, blocks, tolerance=tolerance)


def phonopy_to_model(crystal, force_constants_file, image_map_file, *, tolerance=1e-8):
    mapping = read_image_map(image_map_file)
    return mapped_force_constants(crystal, parse_force_constants_text(force_constants_file), mapping,
                                  tolerance=tolerance)


def force_constants_text(tensor):
    """Phonopy text serializer for tests/interoperability (not a structure exporter)."""
    f = np.asarray(tensor, dtype=float)
    if f.ndim != 4 or f.shape[2:] != (3,3) or not np.isfinite(f).all():
        raise PhononInputError('tensor must be finite NxMx3x3')
    out = [f'{f.shape[0]} {f.shape[1]}']
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            out.append(f'{i+1} {j+1}')
            out.extend(' '.join(f'{x:.17g}' for x in row) for row in f[i,j])
    return '\n'.join(out)+'\n'


def main(argv=None):
    """Convert explicitly mapped FORCE_CONSTANTS to GeoForm's versioned IFC JSON."""
    import argparse
    from .io import model_to_record, save_json, STATUSES
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--force-constants', required=True, help='Phonopy text FORCE_CONSTANTS')
    p.add_argument('--image-map', required=True, help='JSON supercell_images and origin_rows')
    p.add_argument('--crystal', required=True, help='GeoForm primitive crystal JSON')
    p.add_argument('--output', required=True, help='GeoForm versioned IFC JSON')
    p.add_argument('--evidence-status', choices=sorted(STATUSES), required=True)
    p.add_argument('--source', required=True, help='DFT method, publication, or precise measured source')
    args = p.parse_args(argv)
    structure = json.loads(Path(args.crystal).read_text(encoding='utf-8'))
    crystal = Crystal(structure['lattice_angstrom'],structure['fractional_positions'],
                      structure['masses_amu'],tuple(structure['symbols']))
    model = phonopy_to_model(crystal,args.force_constants,args.image_map)
    provenance = {'evidence_status':args.evidence_status,'source':args.source,
                  'importer':'GeoForm strict Phonopy FORCE_CONSTANTS text bridge',
                  'original_format':'Phonopy FORCE_CONSTANTS (eV/angstrom^2); user mapped periodic images'}
    save_json(model_to_record(model,provenance),args.output)
    print(json.dumps({'output':args.output,'n_atoms_primitive':crystal.n_atoms,
                      'supercell_force_constant_blocks':len(model.blocks),
                      'asr_error_ev_a2':model.acoustic_sum_rule_error(),
                      'evidence_declaration':args.evidence_status}))


if __name__ == '__main__':
    main()
