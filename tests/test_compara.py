import os

from cogent3.util.unit_test import TestCase, main

from ensembldb.host import HostAccount, get_ensembl_account
from ensembldb.compara import Compara

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

release = 81

if 'ENSEMBL_ACCOUNT' in os.environ:
    args = os.environ['ENSEMBL_ACCOUNT'].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs['port'] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
else:
    account = get_ensembl_account(release=release)


def calc_slope(x1, y1, x2, y2):
    """computes the slope from two coordinate sets, assigning a delta of 1
    when values are identical"""
    delta_y = y2 - y1
    delta_x = x2 - x1
    delta_y = [delta_y, 1][delta_y == 0]
    delta_x = [delta_x, 1][delta_x == 0]
    return delta_y / delta_x


class ComparaTestBase(TestCase):
    comp = Compara(['human', 'mouse', 'rat', 'platypus'], release=release,
                   account=account)


class TestCompara(ComparaTestBase):

    def test_query_genome(self):
        """compara should attach valid genome attributes by common name"""
        brca2 = self.comp.Mouse.get_gene_by_stableid("ENSMUSG00000041147")
        self.assertEqual(brca2.symbol.lower(), 'brca2')

    def test_get_related_genes(self):
        """should correctly return the related gene regions from each genome"""
        brca2 = self.comp.Mouse.get_gene_by_stableid("ENSMUSG00000041147")
        Orthologs = self.comp.get_related_genes(gene_region=brca2,
                                              Relationship="ortholog_one2one")
        self.assertEqual("ortholog_one2one", Orthologs.Relationships[0])

    def test_get_related_genes2(self):
        """should handle case where gene is absent from one of the genomes"""
        clec2d = self.comp.Mouse.get_gene_by_stableid(
            stableid='ENSMUSG00000030157')
        orthologs = self.comp.get_related_genes(gene_region=clec2d,
                                              Relationship='ortholog_one2many')
        self.assertTrue(len(orthologs.members) < 4)

    def test_get_collection(self):
        brca2 = self.comp.Human.get_gene_by_stableid(stableid="ENSG00000139618")
        Orthologs = self.comp.get_related_genes(gene_region=brca2,
                                              Relationship="ortholog_one2one")
        collection = Orthologs.get_seq_collection()
        self.assertTrue(len(collection.seqs[0]) > 1000)

    def test_getting_alignment(self):
        mid = "ENSMUSG00000041147"
        brca2 = self.comp.Mouse.get_gene_by_stableid(stableid=mid)
        result = list(self.comp.get_syntenic_regions(region=brca2,
                                                   align_method='PECAN', align_clade='vertebrates'))[0]
        aln = result.get_alignment(feature_types='gene')
        # to improve test robustness across Ensembl releases, where alignment
        # coordinates change due to inclusion of new species, we search for
        # the mouse subseq and use the resulting coords to ensure we get the
        # same match as that from the Ensembl website
        mouse_name = [n for n in aln.names if "Mus musculus" in n][0]
        start = aln.todict()[mouse_name].find('AAGTCAAACTCTACCACTGG')
        sub_aln = aln[start: start + 20]
        seqs = list(sub_aln.todict().values())
        expect = set(['AGGGCTGACTCTGCCGCTGT',  # human
                      'AAGTCAAACTCTACCACTGG',  # mouse
                      'AAGTCAAACTCTACCACTAG',  # rat
                      'AAATGTGACTCTACCAGCCG'  # platypus
                      ])
        self.assertEqual(set(seqs), expect)
        self.assertTrue(len(aln) > 1000)

    def test_generate_method_clade_data(self):
        """should correctly determine the align_method align_clade options for
        a group of species"""
        # we should correctly infer the method_species_links, which is a
        # cogent3.util.Table instance
        self.assertTrue(self.comp.method_species_links.shape > (0, 0))

    def test_no_method_clade_data(self):
        """generate a Table with no rows if no alignment data"""
        compara = Compara(['S.cerevisiae'], release=release, account=account)
        self.assertEqual(compara.method_species_links.shape[0], 0)

    def test_get_syntenic_returns_nothing(self):
        """should correctly return None for a SyntenicRegion with golden-path
        assembly gap"""
        start = 100000
        end = start + 100000
        related = list(self.comp.get_syntenic_regions(species='mouse',
                                                    coord_name='1', start=start, end=end,
                                                    align_method='PECAN', align_clade='vertebrates'))
        self.assertEqual(related, [])

    def test_get_species_set(self):
        """should return the correct set of species"""
        expect = set(['Homo sapiens', 'Ornithorhynchus anatinus',
                      'Mus musculus', 'Rattus norvegicus'])
        brca1 = self.comp.Human.get_gene_by_stableid(stableid="ENSG00000012048")
        Orthologs = self.comp.get_related_genes(gene_region=brca1,
                                              Relationship="ortholog_one2one")
        self.assertEqual(Orthologs.get_species_set(), expect)

    def test_pool_connection(self):
        """excercising ability to specify pool connection"""
        dog = Compara(['chimp', 'dog'], release=release, account=account,
                      pool_recycle=1000)


class TestSyntenicRegions(TestCase):
    comp = Compara(['human', 'chimp', 'macaque'], account=account,
                   release=release)

    def test_correct_alignments(self):
        """should return the correct alignments"""
        # following cases have a mixture of strand between ref seq and others
        coords_expected = [
            [{'coord_name': 4, 'end': 78207, 'species': 'human', 'start': 78107, 'strand': -1},
             {'Homo sapiens:chromosome:4:77999-78099:-1':
              'ATGTAAATCAAAACCAAAGTCTGCATTTATTTGCGGAAAGAGATGCTACATGTTCAAAGATAAATATGGAACATTTTTTAAAAGCATTCATGACTTAGAA',
              'Macaca mulatta:chromosome:1:3891064-3891163:1':
              'ATGTCAATCAAAACCAAAGTCTGTATTTATTTGCAGAAAGAGATACTGCATGTTCAAAGATAAATATGGAAC-TTTTTAAAAAGCATTAATGACTTATAC',
              'Pan troglodytes:chromosome:4:102056-102156:-1':
              'ATGTAAATCAAAACCAAAGTCTGCATTTATTTGCGGAAAGAGATGCTACATGTTCAAAGATAAATATGGAACATTTTTAAAAAGCATTCATGACTTAGAA'}],
            [{'coord_name': 18, 'end': 213739, 'species': 'human', 'start': 213639, 'strand': -1},
                {'Homo sapiens:chromosome:18:213639-213739:-1':
                 'ATAAGCATTTCCCTTTAGGGCTCTAAGATGAGGTCATCATCGTTTTTAATCCTGAAGAAGGGCTACTGAGTGAGTGCAGATTATTCGGTAAACACT----CTTA',
                 'Macaca mulatta:chromosome:18:13858303-13858397:1':
                 '------GTTTCCCTTTAGGGCTCTAAGATGAGGTCATCATTGTTTTTAATCCTGAAGAAGGGCTACTGA----GTGCAGATTATTCTGTAAATGTGCTTACTTG',
                 'Pan troglodytes:chromosome:18:16601082-16601182:1':
                 'ATAAGCATTTCCCTTTAGGGCTCTAAGATGAGGTCATCATCGTTTTTAATCCTGAAGAAGGGCTACTGA----GTGCAGATTATTCTGTAAACACTCACTCTTA'}],
            [{'coord_name': 5, 'end': 204859, 'species': 'human', 'start': 204759, 'strand': 1},
                {'Homo sapiens:chromosome:5:204874-204974:1':
                 'AACACTTGGTATTT----CCCCTTTATGGAGTGAGAGAGATCTTTAAAATATAAACCCTTGATAATATAATATTACTACTTCCTATTA---CCTGTTATGCAGTTCT',
                 'Macaca mulatta:chromosome:6:1297736-1297840:-1':
                 'AACTCTTGGTGTTTCCTTCCCCTTTATGG---GAGAGAGATCTTTAAAATAAAAAACCTTGATAATATAATATTACTACTTTCTATTATCATCTGTTATGCAGTTCT',
                 'Pan troglodytes:chromosome:5:335911-336011:1':
                 'AACACTTGGTAGTT----CCCCTTTATGGAGTGAGAGAGATCTTTAAAATATAAACCCTTGATAATATAATATTACTACTTTCTATTA---CCTGTTATGCAGTTCT'}],
            [{'coord_name': 18, 'end': 203270, 'species': 'human', 'start': 203170, 'strand': -1},
                {'Homo sapiens:chromosome:18:203170-203270:-1':
                 'GGAATAATGAAAGCAATTGTGAGTTAGCAATTACCTTCAAAGAATTACATTTCTTATACAAAGTAAAGTTCATTACTAACCTTAAGAACTTTGGCATTCA',
                 'Pan troglodytes:chromosome:18:16611584-16611684:1':
                 'GGAATAATGAAAGCAATTGTAAGTTAGCAATTACCTTCAAAGAATTACATTTCTTATACAAAGTAAAGTTCATTACTAACCTTAAGAACTTTGGCATTCA'}],
            [{'coord_name': 2, 'end': 46445, 'species': 'human', 'start': 46345, 'strand': -1},
                {'Homo sapiens:chromosome:2:46345-46445:-1':
                 'CTACCACTCGAGCGCGTCTCCGCTGGACCCGGAACCCCGGTCGGTCCATTCCCCGCGAAGATGCGCGCCCTGGCGGCCCTGAGCGCGCCCCCGAACGAGC',
                 'Macaca mulatta:chromosome:13:43921-44021:-1':
                 'CTGCCACTCCAGCGCGTCTCCGCTGCACCCGGAGCGCCGGCCGGTCCATTCCCCGCGAGGATGCGCGCCCTGGCGGCCCTGAACACGTCGGCGAGAGAGC',
                 'Pan troglodytes:chromosome:2a:36792-36892:-1':
                 'CTACCACTCGAGCGCGTCTCCGCTGGACCCGGAACCCCAGTCGGTCCATTCCCCGCGAAGATGCGCGCCCTGGCGGCCCTGAACGCGCCCCCGAACGAGC'}],
            [{'coord_name': 18, 'end': 268049, 'species': 'human', 'start': 267949, 'strand': -1},
                {'Homo sapiens:chromosome:18:267949-268049:-1':
                 'GCGCAGTGGCGGGCACGCGCAGCCGAGAAGATGTCTCCGACGCCGCCGCTCTTCAGTTTGCCCGAAGCGCGGACGCGGTTTACGGTGAGCTGTAGAGGGG',
                 'Macaca mulatta:chromosome:18:13805604-13805703:1':
                 'GCGCAG-GGCGGGCACGCGCAGCCGAGAAGATGTCTCCGACGCCGCCGCTCTTCAGTTTGCCCGAAGCGCGGACGCGGTTTACGGTGAGCTGTAGGCGGG',
                 'Pan troglodytes:chromosome:18:16546800-16546900:1':
                 'GCGCAGTGGCGGGCACGCGCAGCCGAGAAGATGTCTCCGACGCCGCCGCTCTTCAGTTTGCCCGAAGCGCGGACGCGGTTTACGGTGAGCTGTAGCGGGG'}],
            [{'coord_name': 16, 'end': 57443, 'species': 'human', 'start': 57343, 'strand': -1},
                {'Homo sapiens:chromosome:16:107343-107443:-1':
                 'AAGAAGCAAACAGGTTTATTTTATACAGTGGGCCAGGCCGTGGGTCTGCCATGTGACTAGGGCATTTGGACCTAGGGAGAGGTCAGTCTCAGGCCAAGTA',
                 'Pan troglodytes:chromosome:16:48943-49032:-1':
                 'AAGAAGCAAACAGGTTTATTTTATACACTGGGCCAGGCCGTGGGTCTGCCATGTGACTAGGGAATTTGGACC-----------CAGTCTCAGGCCAAGTA'}]
        ]
        # print(self.comp.method_species_links)
        for coord, expect in coords_expected[1:]:
            syntenic = list(
                self.comp.get_syntenic_regions(method_clade_id=756, **coord))[0]
            # check the slope computed from the expected and returned
            # coordinates is ~ 1
            got_names = dict([(n.split(':')[0], n.split(':'))
                             for n in syntenic.get_alignment().names])
            exp_names = dict([(n.split(':')[0], n.split(':'))
                             for n in list(expect.keys())])
            for species in exp_names:
                exp_chrom = exp_names[species][2]
                got_chrom = got_names[species][2]
                self.assertEqual(exp_chrom.lower(), got_chrom.lower())
                exp_start, exp_end = list(
                    map(int, exp_names[species][3].split('-')))
                got_start, got_end = list(
                    map(int, got_names[species][3].split('-')))
                slope = calc_slope(exp_start, exp_end, got_start, got_end)
                self.assertFloatEqual(abs(slope), 1.0, eps=1e-3)

    def test_failing_region(self):
        """should correctly handle queries where multiple Ensembl have
        genome block associations for multiple coord systems"""
        gene = self.comp.Human.get_gene_by_stableid(stableid='ENSG00000188554')
        # this should simply not raise any exceptions
        syntenic_regions = list(self.comp.get_syntenic_regions(region=gene,
                                                             align_method='PECAN',
                                                             align_clade='vertebrates'))


if __name__ == "__main__":
    main()
