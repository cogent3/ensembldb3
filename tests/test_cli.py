import os
import shutil

from configparser import ConfigParser

import pytest

from click.testing import CliRunner

from ensembl_cli.cli import download, exportrc, install
from ensembl_cli.util import get_resource_path


@pytest.fixture(scope="function")
def tmp_config(tmp_dir):
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(get_resource_path("ensembldb_download.cfg"))
    parser.remove_section("C.elegans")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "wt") as out:
        parser.write(out)

    yield download_cfg


def test_download(tmp_config):
    """runs download, install, drop according to a special test cfg"""
    tmp_dir = tmp_config.parent
    # now download
    runner = CliRunner()
    r = runner.invoke(download, [f"-c{tmp_config}"], catch_exceptions=False)
    assert r.exit_code == 0, r.output
    # make sure the download checkpoint file exists
    dirnames = [
        dn
        for dn in os.listdir(tmp_dir / "staging")
        if (tmp_dir / "staging" / dn).is_dir()
    ]
    assert "saccharomyces_cerevisiae" in dirnames

    # make sure file sizes > 0
    paths = list((tmp_dir / "staging" / "saccharomyces_cerevisiae").glob("*"))
    size = sum(p.stat().st_size for p in paths)
    assert size > 0

    assert r.exit_code == 0, r.output


def test_exportrc(tmp_dir):
    """exportrc works correctly"""
    runner = CliRunner()
    outdir = tmp_dir / "exported"
    r = runner.invoke(exportrc, [f"-o{outdir}"])
    assert r.exit_code == 0, r.output
    fnames = os.listdir(outdir)
    assert "species.tsv" in fnames
    assert len(fnames) == 2
    shutil.rmtree(tmp_dir)
