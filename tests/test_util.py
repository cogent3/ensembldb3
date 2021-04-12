import pathlib
from unittest import TestCase, main

from ensembldb3.util import FileSet

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


class TestFileSet(TestCase):
    def test_init_no_matches(self):
        """correct constructs when no matches"""
        fset = FileSet(pathlib.Path(__file__).parent, suffixes="txt")
        self.assertEqual(fset, set())

    def test_init_matches(self):
        """correctly construct with matches"""
        fset = FileSet(pathlib.Path(__file__).parent, suffixes="py")
        self.assertEqual(
            len(fset), len(list(pathlib.Path(__file__).parent.glob("*.py")))
        )

        # handles whether . included in prefix
        fset = FileSet(pathlib.Path(__file__).parent, suffixes=".py")
        self.assertEqual(
            len(fset), len(list(pathlib.Path(__file__).parent.glob("*.py")))
        )

    def test_init_with_multi_suffix(self):
        """correctly identifies multi-suffixes"""
        from tempfile import TemporaryDirectory

        with TemporaryDirectory(".") as dname:
            dname = pathlib.Path(dname)
            for i in range(3):
                name = dname / f"{i}.txt.gz"
                name.write_text("blah")

            for i in range(3, 6):
                name = dname / f"{i}.sql.gz"
                name.write_text("blah")

            name = dname / ".empty"
            name.write_text("blah")

            name = dname / "nesteddir"
            name.mkdir()
            name = name / "nested.txt"
            name.write_text("blah")

            fset = FileSet(dname)
            self.assertEqual(len(fset), 6)


if __name__ == "__main__":
    main()
