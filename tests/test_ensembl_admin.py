import os
import shutil

from configparser import ConfigParser

import pytest

from click.testing import CliRunner

from ensembl_cli.admin import ENSEMBLDBRC, download, exportrc
from ensembl_cli.download import get_download_checkpoint_path


dirname = "_delme"


def test_download():
    """runs download, install, drop according to a special test cfg"""

    test_mysql_cfg = os.environ.get("ENSEMBLDB_TEST_CFG", None)
    if os.path.exists(dirname):
        shutil.rmtree(dirname)

    os.makedirs(dirname)

    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg"))
    parser.remove_section("C.elegans")
    parser.set("local path", "path", value=dirname)
    download_cfg = os.path.abspath(os.path.join(dirname, "download.cfg"))
    with open(download_cfg, "wt") as out:
        parser.write(out)

    # now download
    runner = CliRunner()
    r = runner.invoke(download, [f"-c{download_cfg}"])
    # make sure the download checkpoint file exists
    dirnames = [
        dn for dn in os.listdir(dirname) if os.path.isdir(os.path.join(dirname, dn))
    ]
    assert len(dirnames) == 1
    chkpt = get_download_checkpoint_path(dirname, dirnames[0])
    assert os.path.exists(chkpt)

    # make sure file sizes > 0
    fnames = os.listdir(os.path.join(dirname, dirnames[0]))
    size = 0
    for fn in fnames:
        path = os.path.join(dirname, dirnames[0], fn)
        size += os.path.getsize(path)
    assert size > 0

    if r.exit_code != 0:
        print(r.output)

    assert r.exit_code == 0


def test_exportrc():
    """exportrc works correctly"""
    runner = CliRunner()

    if os.path.exists(dirname):
        shutil.rmtree(dirname)

    r = runner.invoke(exportrc, [f"-o{dirname}"])
    assert r.exit_code == 0
    fnames = os.listdir(dirname)
    assert "species.tsv" in fnames
    assert len(fnames) == 3
    shutil.rmtree(dirname)
