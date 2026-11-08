import io
import json
import lzma
from zipfile import ZipFile

import numpy as np
import pytest

from geoform_phonons.core import Crystal
from geoform_phonons.fixtures import from_springs
from geoform_phonons import phonondb
from geoform_phonons.phonondb_cli import select


ROW = (
    '| 661 | AlN | P6_3mc (186) | https://mdr.nims.go.jp/concern/datasets/9306t5111 '
    '| [9306t5111.zip](https://mdr.nims.go.jp/download_all/9306t5111.zip) |'
)


def test_index_binds_record_ids_to_official_urls():
    manifest = phonondb.parse_index(ROW.encode(), '0' * 40)
    row = manifest['rows'][0]
    assert row['mp_id'] == 'mp-661'
    assert row['space_group_number'] == 186
    with pytest.raises(ValueError, match='Unrecognized upstream index'):
        phonondb.parse_index(ROW.replace('download_all/', 'download_all/9').encode(), '0' * 40)
    with pytest.raises(ValueError, match='repeated MDR'):
        phonondb.parse_index((ROW + '\n' + ROW).encode(), '0' * 40)


def test_checksum_lock_rejects_changed_bytes(tmp_path, monkeypatch):
    row = phonondb.parse_index(ROW.encode(), '0' * 40)['rows'][0]
    payload = b'PK\x03\x04test'
    monkeypatch.setattr(phonondb, 'urlopen', lambda *a, **kw: io.BytesIO(payload))
    data, pin = phonondb.locked_archive(row, tmp_path, acquire=True)
    assert data == payload
    assert pin['sha256'] == phonondb.digest(payload)
    assert phonondb.locked_archive(row, tmp_path)[0] == payload
    (tmp_path / '9306t5111.zip').write_bytes(payload + b'changed')
    with pytest.raises(ValueError, match='differs from pinned bytes'):
        phonondb.locked_archive(row, tmp_path, acquire=True)
    with pytest.raises(ValueError, match='differs from supplied pin'):
        phonondb.locked_archive(row, tmp_path, acquire=True, expected_pin={**pin, 'sha256': 'a' * 64})


def test_rights_use_record_dataset_metadata(monkeypatch):
    row = phonondb.parse_index(ROW.encode(), '0' * 40)['rows'][0]
    record = {'@type': 'Dataset', 'name': 'AlN / materials id 661', 'license': [{'@id': 'cc-by-4.0'}]}
    script = '//<![CDATA[\n' + json.dumps({'@graph': [record]}) + '\n//]]>'
    page = f'<script type="application/ld+json" nonce="">{script}</script>'.encode()
    monkeypatch.setattr(phonondb, 'urlopen', lambda *a, **kw: io.BytesIO(page))
    assert phonondb.inspect_record_license(row)['status'] == 'record_cc_by_4_0'
    with pytest.raises(ValueError, match='does not identify'):
        phonondb.inspect_record_license({**row, 'mp_id': 'mp-66'})


def test_selection_rejects_unknown_ids(tmp_path):
    manifest = phonondb.parse_index(ROW.encode(), '0' * 40)
    ids = tmp_path / 'ids.json'
    ids.write_text('["wrong"]')
    with pytest.raises(ValueError, match='unknown MDR IDs'):
        select(manifest, ids, 0, None)


def test_batched_fourier_matches_scalar_for_complex_basis():
    crystal = Crystal(np.eye(3) * 4, [[0, 0, 0], [0.25, 0, 0]], [1, 2], ("A", "B"))
    model = from_springs(crystal, [(0, 1, (0, 0, 0), np.eye(3))])
    q = [(0, 0, 0), (0.125, 0.3, 0.1), (0.5, 0, 0), (1.17, 0.4, 0.2)]
    expected = np.asarray([model.dynamical_matrix(point) for point in q])
    assert np.max(np.abs(expected.imag)) > 0
    np.testing.assert_allclose(model.dynamical_matrices(q), expected, atol=1e-14)


def test_grouped_formulas_and_bounded_package_reader():
    pytest.importorskip('phonopy')
    from geoform_phonons.phonondb_verify import PackageRejected, indexed_formula, read_package

    assert indexed_formula('Ca3(PO4)2') == {'Ca': 3, 'P': 2, 'O': 8}
    assert indexed_formula('Na2Mg(CO3)2') == {'Na': 2, 'Mg': 1, 'C': 2, 'O': 6}
    assert indexed_formula('N2') == {'N': 1}
    with pytest.raises(PackageRejected):
        indexed_formula('Ca3(PO4')
    packed = lzma.compress(b'phonopy: {}\n')
    buf = io.BytesIO()
    with ZipFile(buf, 'w') as archive:
        archive.writestr('nested/phonopy_params.yaml.xz', packed)
    assert read_package(buf.getvalue()) == ({'phonopy': {}}, phonondb.digest(packed),
                                           phonondb.digest(b'phonopy: {}\n'))
    with ZipFile(buf, 'a') as archive:
        archive.writestr('phonopy_params.yaml.xz', packed)
    with pytest.raises(PackageRejected, match='Expected one'):
        read_package(buf.getvalue())
