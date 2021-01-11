from unittest import TestCase, main

from ensembldb3.compara import Compara, Genome
from ensembldb3.host import HostAccount

from . import ENSEMBL_GENOMES_RELEASE

__author__ = "Jason Merkin"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Jason Merkin", "Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

account = HostAccount("mysql-eg-publicsql.ebi.ac.uk", "anonymous", "", port=4157)


class MZ_ComparaTestBase(TestCase):
    comp = Compara(
        ["D.grimshawi", "D.melanogaster"],
        release=ENSEMBL_GENOMES_RELEASE,
        account=account,
        division="metazoa",
    )


class MZ_TestCompara(MZ_ComparaTestBase):
    def test_query_genome(self):
        """compara should attach valid genome attributes by common name"""
        brca2 = self.comp.Dmelanogaster.get_gene_by_stableid("FBgn0050169")
        self.assertEqual(brca2.symbol.lower(), "brca2")

    def test_get_related_genes(self):
        """should correctly return the related gene regions from each genome"""
        # using sc35, a splicing factor
        sc35 = self.comp.Dmelanogaster.get_gene_by_stableid("FBgn0265298")
        orthologs = self.comp.get_related_genes(
            gene_region=sc35, relationship="ortholog_one2one"
        )
        self.assertEqual("ortholog_one2one", list(orthologs)[0].relationship)

    def test_get_collection(self):
        """return a sequence collection"""
        sc35 = self.comp.Dmelanogaster.get_gene_by_stableid(stableid="FBgn0265298")
        orthologs = self.comp.get_related_genes(
            gene_region=sc35, relationship="ortholog_one2one"
        )
        collection = list(orthologs)[0].get_seq_collection()
        self.assertTrue(len(collection.seqs) == 2)
        self.assertTrue(len(collection.seqs[0]) > 1000)


class MZ_Genome(TestCase):
    def test_get_general_release(self):
        """should correctly infer the general release"""
        rel_gt_65 = Genome(
            "D.melanogaster", release=ENSEMBL_GENOMES_RELEASE, account=account
        )
        self.assertEqual(rel_gt_65.general_release, 98)
        self.assertEqual(
            rel_gt_65.CoreDb.db_name.name, "drosophila_melanogaster_core_45_98_7"
        )


if __name__ == "__main__":
    main()
