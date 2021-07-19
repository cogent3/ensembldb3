import os
import shutil

from configparser import ConfigParser
from unittest import TestCase, main

from click.testing import CliRunner

from ensembldb3.admin import (
    ENSEMBLDBRC,
    download,
    drop,
    exportrc,
    get_installed_checkpoint_path,
    install,
    show,
)
from ensembldb3.download import get_download_checkpoint_path


class TestAdminCli(TestCase):
    dirname = "_delme"

    def test_all(self):
        """runs download, install, drop according to a special test cfg"""

        def exec_show(mysql_cfg, release):
            runner = CliRunner()
            r = runner.invoke(show, [f"-r{release}", f"-m{mysql_cfg}"])
            self.assertEqual(r.exit_code, 0)
            return r

        test_mysql_cfg = os.environ.get("ENSEMBLDB_TEST_CFG", None)
        if test_mysql_cfg is None:
            self.skipTest(
                "ENSEMBLDB_TEST_CFG variable not defined, " "skipping some cli tests"
            )

        if os.path.exists(self.dirname):
            shutil.rmtree(self.dirname)

        os.makedirs(self.dirname)

        # create a simpler download config
        # we want a very small test set
        parser = ConfigParser()
        parser.read(os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg"))
        parser.remove_section("C.elegans")
        parser.set("local path", "path", value=self.dirname)
        download_cfg = os.path.abspath(os.path.join(self.dirname, "download.cfg"))
        with open(download_cfg, "wt") as out:
            parser.write(out)

        # now download
        runner = CliRunner()
        r = runner.invoke(download, [f"-c{download_cfg}"])
        # make sure the download checkpoint file exists
        dirnames = [
            dn
            for dn in os.listdir(self.dirname)
            if os.path.isdir(os.path.join(self.dirname, dn))
        ]
        self.assertEqual(len(dirnames), 1)
        chkpt = get_download_checkpoint_path(self.dirname, dirnames[0])
        self.assertTrue(os.path.exists(chkpt))

        # make sure file sizes > 0
        fnames = os.listdir(os.path.join(self.dirname, dirnames[0]))
        size = 0
        for fn in fnames:
            path = os.path.join(self.dirname, dirnames[0], fn)
            size += os.path.getsize(path)
        self.assertTrue(size > 0)

        if r.exit_code != 0:
            print(r.output)

        self.assertEqual(r.exit_code, 0)

        # now install
        runner = CliRunner()
        r = runner.invoke(
            install,
            [f"-c{download_cfg}", f"-m{test_mysql_cfg}"],
            catch_exceptions=False,
        )
        if r.exit_code != 0:
            print(r.output)

        self.assertEqual(r.exit_code, 0)
        chkpt = get_installed_checkpoint_path(self.dirname, dirnames[0])
        # check it's installed via checkpoint file
        self.assertTrue(os.path.exists(chkpt))

        # and using show
        release = parser.get("release", "release")
        r = exec_show(test_mysql_cfg, release)
        self.assertTrue("saccharomyces" in r.output)

        # then drop, but don't execute
        runner = CliRunner()
        r = runner.invoke(
            drop,
            [f"-c{download_cfg}", f"-m{test_mysql_cfg}"],
            catch_exceptions=False,
            input="n",
        )

        # then show still present
        r = exec_show(test_mysql_cfg, release)
        self.assertTrue("saccharomyces" in r.output)

        # then drop
        runner = CliRunner()
        r = runner.invoke(
            drop,
            [f"-c{download_cfg}", f"-m{test_mysql_cfg}"],
            catch_exceptions=False,
            input="y",
        )
        # then show now gone
        r = exec_show(test_mysql_cfg, release)
        self.assertTrue("saccharomyces" not in r.output)

        shutil.rmtree(self.dirname)

    def test_exportrc(self):
        """exportrc works correctly"""
        runner = CliRunner()

        if os.path.exists(self.dirname):
            shutil.rmtree(self.dirname)

        r = runner.invoke(exportrc, [f"-o{self.dirname}"])
        self.assertEqual(r.exit_code, 0)
        fnames = os.listdir(self.dirname)
        self.assertTrue("species.tsv" in fnames)
        self.assertEqual(len(fnames), 3)
        shutil.rmtree(self.dirname)


if __name__ == "__main__":
    main()
