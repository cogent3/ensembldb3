import os
import shutil

from configparser import ConfigParser

import pytest

from click.testing import CliRunner

from ensembl_cli.admin import ENSEMBLDBRC, download, exportrc
from ensembl_cli.download import get_download_checkpoint_path


def test_download(tmp_dir):
    """runs download, install, drop according to a special test cfg"""

    test_mysql_cfg = os.environ.get("ENSEMBLDB_TEST_CFG", None)
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg"))
    parser.remove_section("C.elegans")
    parser.set("local path", "path", value=str(tmp_dir))
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "wt") as out:
        parser.write(out)

    # now download
    runner = CliRunner()
    r = runner.invoke(download, [f"-c{download_cfg}"])
    # make sure the download checkpoint file exists
    dirnames = [dn for dn in os.listdir(tmp_dir) if (tmp_dir / dn).is_dir()]
    assert len(dirnames) == 1
    chkpt = get_download_checkpoint_path(tmp_dir, dirnames[0])
    assert os.path.exists(chkpt)

    # make sure file sizes > 0
    fnames = os.listdir(tmp_dir / dirnames[0])
    size = 0
    for fn in fnames:
        path = tmp_dir / dirnames[0] / fn
        size += os.path.getsize(path)
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
    assert len(fnames) == 3
    shutil.rmtree(tmp_dir)
