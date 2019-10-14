import os

from cogent3.util.unit_test import TestCase, main
from ensembldb3.database import Database
from ensembldb3.host import HostAccount, get_ensembl_account

from . import ENSEMBL_RELEASE

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

if "ENSEMBL_ACCOUNT" in os.environ:
    args = os.environ["ENSEMBL_ACCOUNT"].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs["port"] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
else:
    account = get_ensembl_account(release=ENSEMBL_RELEASE)


class TestDatabase(TestCase):
    def test_connect(self):
        human = Database(
            account=account, release=ENSEMBL_RELEASE, species="human", db_type="core"
        )
        gene = human.get_table("gene")

    def test_get_distinct(self):
        """should return list of strings"""
        db = Database(
            account=account, release=ENSEMBL_RELEASE, species="human", db_type="core"
        )
        tn, tc = "gene", "biotype"
        expected = set(
            [
                "protein_coding",
                "pseudogene",
                "processed_transcript",
                "Mt_tRNA",
                "Mt_rRNA",
                "IG_V_gene",
                "IG_J_gene",
                "IG_C_gene",
                "IG_D_gene",
                "miRNA",
                "misc_RNA",
                "snoRNA",
                "snRNA",
                "rRNA",
            ]
        )
        got = set(db.get_distinct(tn, tc))
        self.assertNotEqual(set(got) & expected, set())

        db = Database(account=account, release=ENSEMBL_RELEASE, db_type="compara")
        got = set(db.get_distinct("homology", "description"))
        expected = set(
            [
                "gene_split",
                "alt_allele",
                "other_paralog",
                "ortholog_one2many",
                "ortholog_one2one",
                "within_species_paralog",
                "ortholog_many2many",
            ]
        )
        self.assertEqual(len(got & expected), len(expected))

    def test_get_table_row_counts(self):
        """should return correct row counts for some tables"""
        expect = {
            f"homo_sapiens_core_{ENSEMBL_RELEASE}_38.analysis": 61,
            f"homo_sapiens_core_{ENSEMBL_RELEASE}_38.seq_region": 55616,
            f"homo_sapiens_core_{ENSEMBL_RELEASE}_38.assembly": 102090,
        }
        human = Database(
            account=account, release=ENSEMBL_RELEASE, species="human", db_type="core"
        )
        table_names = [n.split(".")[1] for n in expect]
        got = dict(human.get_tables_row_count(table_names).tolist())
        for dbname in expect:
            self.assertTrue(got[dbname] >= expect[dbname])

    def test_table_has_column(self):
        """return correct values for whether a Table has a column"""
        coredb = Database(
            account=account,
            release=ENSEMBL_RELEASE,
            species="human",
            db_type="core",
        )

        self.assertTrue(coredb.table_has_column("seq_region", "name"))
        self.assertFalse(coredb.table_has_column("seq_region", "length"))


if __name__ == "__main__":
    main()
