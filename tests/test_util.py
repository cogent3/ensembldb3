from configparser import ConfigParser

import pytest

from ensembl_cli.cli import ENSEMBLDBRC


@pytest.fixture(scope="function")
def tmp_config(tmp_dir):
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(ENSEMBLDBRC / "ensembldb_download.cfg")
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
