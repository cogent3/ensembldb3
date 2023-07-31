from configparser import ConfigParser

import pytest

from ensembl_cli.util import (
    get_resource_path,
    load_ensembl_checksum,
    load_ensembl_md5sum,
)


@pytest.fixture(scope="function")
def tmp_config(tmp_dir):
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(get_resource_path("ensembldb_download.cfg"))
    parser.remove_section("C.elegans")
    parser.set("local path", "path", value=str(tmp_dir))
    alns = ",".join(("17_sauropsids.epc", "10_primates.epo"))
    parser.set("compara", "align_names", value=alns)
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "wt") as out:
        parser.write(out)

    yield download_cfg


def test_parse_config(tmp_config):
    from ensembl_cli.util import read_config

    cfg = read_config(tmp_config)
    assert set(cfg.align_names) == {"17_sauropsids.epc", "10_primates.epo"}


def test_load_ensembl_md5sum(DATA_DIR):
    got = load_ensembl_md5sum(DATA_DIR / "sample-MD5SUM")
    assert len(got) == 3
    assert got["b.emf.gz"] == "3d9af835d9ed19975bd8b2046619a3a1"


def test_load_ensembl_checksum(DATA_DIR):
    got = load_ensembl_checksum(DATA_DIR / "sample-CHECKSUMS")
    assert len(got) == 4  # README line is ignored
    assert got["c.fa.gz"] == (7242, 327577)
