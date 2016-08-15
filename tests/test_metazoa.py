from cogent3.util.unit_test import TestCase, main

from ensembldb.host import HostAccount, get_ensembl_account
from ensembldb.compara import Compara, Genome

__author__ = "Jason Merkin"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

release = 23
account = HostAccount('mysql-eg-publicsql.ebi.ac.uk',
                      'anonymous', '', port=4157)


class MZ_ComparaTestBase(TestCase):
    comp = Compara(['D.grimshawi', 'D.melanogaster'], release=release,
                   account=account, division='metazoa')


class MZ_TestCompara(MZ_ComparaTestBase):

    def test_query_genome(self):
        """compara should attach valid genome attributes by common name"""
        brca2 = self.comp.Dmelanogaster.get_gene_by_stableid("FBgn0050169")
        self.assertEqual(brca2.symbol.lower(), 'brca2')

    def test_get_related_genes(self):
        """should correctly return the related gene regions from each genome"""
        # using sc35, a splicing factor
        sc35 = self.comp.Dmelanogaster.get_gene_by_stableid("FBgn0040286")
        Orthologs = self.comp.get_related_genes(gene_region=sc35,
                                              Relationship="ortholog_one2one")
        self.assertEqual("ortholog_one2one", Orthologs.Relationships[0])

    def test_get_related_genes2(self):
        """should handle case where gene is absent from one of the genomes"""
        # here, it is brca2
        brca2 = self.comp.Dmelanogaster.get_gene_by_stableid(
            stableid='FBgn0050169')
        orthologs = self.comp.get_related_genes(gene_region=brca2,
                                              Relationship='ortholog_one2one')
        self.assertEqual(len(orthologs.Members), 2)

    def test_get_collection(self):
        sc35 = self.comp.Dmelanogaster.get_gene_by_stableid(
            stableid="FBgn0040286")
        Orthologs = self.comp.get_related_genes(gene_region=sc35,
                                              Relationship="ortholog_one2one")
        collection = Orthologs.get_seq_collection()
        self.assertTrue(len(collection.seqs[0]) > 1000)


class MZ_Genome(TestCase):

    def test_get_general_release(self):
        """should correctly infer the general release"""
        rel_lt_65 = Genome('D.melanogaster', release=22, account=account)
        self.assertEqual(rel_lt_65.general_release, 75)
        self.assertEqual(rel_lt_65.CoreDb.db_name,
                         'drosophila_melanogaster_core_22_75_546')

        rel_gt_65 = Genome('D.melanogaster', release=23, account=account)
        self.assertEqual(rel_gt_65.general_release, 76)
        self.assertEqual(rel_gt_65.CoreDb.db_name,
                         'drosophila_melanogaster_core_23_76_546')


if __name__ == "__main__":
    main()
