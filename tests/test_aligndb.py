from ensembl_cli.aligndb import AlignDb
from ensembl_cli.install import _emf_to_gap_records


def test_db_align(DATA_DIR, tmp_path):
    records = _emf_to_gap_records(DATA_DIR / "sample.emf")
    outpath = tmp_path / "blah.sqlitedb"
    db = AlignDb(source=outpath)
    db.add_records(records)
    orig = len(db)
    db.close()
    got = AlignDb(source=outpath)
    assert len(got) == orig
