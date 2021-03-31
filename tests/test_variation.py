import os
from unittest import TestCase, main

from ensembldb3.database import Database
from ensembldb3.genome import Genome
from ensembldb3.host import HostAccount, get_ensembl_account

from . import ENSEMBL_RELEASE

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

NULL_VALUE = None

# TODO fix flanking sequence issue, a flag somewhere indicating flank is same as reference
# TODO grab the ancestral allele (from variation.ancestral_allele)

if "ENSEMBL_ACCOUNT" in os.environ:
    args = os.environ["ENSEMBL_ACCOUNT"].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs["port"] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
else:
    account = get_ensembl_account(release=ENSEMBL_RELEASE)


class GenomeTestBase(TestCase):
    human = Genome(species="human", release=ENSEMBL_RELEASE, account=account)
    mouse = Genome(species="mouse", release=ENSEMBL_RELEASE, account=account)
    rat = Genome(species="rat", release=ENSEMBL_RELEASE, account=account)
    macaq = Genome(species="macaque", release=ENSEMBL_RELEASE, account=account)
    gorilla = Genome(species="gorilla", release=ENSEMBL_RELEASE, account=account)
    brca2 = human.get_gene_by_stableid(stableid="ENSG00000139618")


class TestVariation(GenomeTestBase):
    snp_names = ["rs34213141", "rs12791610", "rs10792769", "rs11545807", "rs11270496"]
    snp_nt_alleles = ["G/A", "C/T", "A/G", "C/A", "CAGCTCCAGCTC/-"]
    snp_aa_alleles = ["G/R", "P/L", "Y/C", "V/F", "GAGAV/V"]
    snp_effects = [
        ["missense_variant"],
        ["missense_variant"],
        ["missense_variant"],
        ["upstream_gene_variant", "missense_variant", "regulatory_region_variant"],
        ["non_synonymous_codon"],
    ]
    snp_nt_len = [1, 1, 1, 1, 12]
    map_weights = [1, 1, 1, 1, 1]
    snp_flanks = [
        (
            "CTGAGGTGAGCCAGCGTTGGAGCTGTTTTTCCTTTCAGTATGAATTCCACAAGGAAATCATCTCAGGAGGAAGGGCTCATACTTGGATCCAGAAAATATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCC",
            "GAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCACCCGTCTACTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCC",
        ),
        (
            "GAAAATATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCCGGAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCAC",
            "CGTCTACTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCCAGTGCAGCTGCTGTAAGCCCTGTTGCTCCTCCTCGGGTCGTGGGTCATCCTGCTGCCAATCCAGCTGCTGCAAGCCCTGCTGCTCATCCTC",
        ),
        (
            "ATCAACATAGCCAAAGAAAAACAATCAAGACATACCTCCAGGAGCTGTGTAACAGCAACCGGAAAGAGAAACAATGGTGTGTTCCTATGTGGGATATAAAGAGCCGGGGCTCAGGGGGCTCCACACCTGCACCTCCTTCTCACCTGCTCCTCTACCTGCTCCACCCTCAATCCACCAGAACCATGGGCTGCTGTGGCTGCTCCGGAGGCTGTGGCTCCAGCTGTGGAGGCTGTGACTCCAGCTGTGGGAGCTGTGGCTCTGGCTGCAGGGGCTGTGGCCCCAGCTGCTGTGCACCCGTCT",
            "CTGCTGCAAGCCCGTGTGCTGCTGTGTTCCAGCCTGTTCCTGCTCTAGCTGTGGCAAGCGGGGCTGTGGCTCCTGTGGGGGCTCCAAGGGAGGCTGTGGTTCTTGTGGCTGCTCCCAGTGCAGTTGCTGCAAGCCCTGCTGTTGCTCTTCAGGCTGTGGGTCATCCTGCTGCCAGTGCAGCTGCTGCAAGCCCTACTGCTCCCAGTGCAGCTGCTGTAAGCCCTGTTGCTCCTCCTCGGGTCGTGGGTCATCCTGCTGCCAATCCAGCTGCTGCAAGCCCTGCTGCTCATCCTCAGGCTG",
        ),
        (
            "GCTGAAGAAACCATTTCAAACAGGATTGGAATAGGGAAACCCGGCACTCAGCTCGGCGCAAGCCGGCGGTGCCTTCAGACTAGAGAGCCTCTCCTCCGGTGCGCTGCAAGTAGGGCCTCGGCTCGAGGTCAACATTCTAGTTGTCCAGCGCTCCCTCTCCGGCACCTCGGTGAGGCTAGTTGACCCGACAGGCGCGGATCATGAGCAGCTGCAGGAGAATGAAGAGCGGGGACGTAATGAGGCCGAACCAGAGCTCCCGAGTCTGCTCCGCCAGCTTCTGGCACAACAGCATCTCGAAGA",
            "GAACTTGAGACTCAGGACCGTAAGTACCCAGAAAAGGCGGAGCACCGCCAGCCGCTTCTCTCCATCCTGGAAGAGGCGCACGGACACGATGGTGGTGAAGTAGGTGCTGAGCCCGTCAGCGGCGAAGAAAGGCACGAACACGTTCCACCAGGAGAGGCCCGGGACCAGGCCATCCACACGCAGTGCCAGCAGCACAGAGAACACCAACAGGGCCAGCAGGTGCACGAAGATCTCGAAGGTGGCGAAGCCTAGCCACTGCACCAGCTCCCGGAGCGAGAAGAGCATCGCGCCCGTTGAGCG",
        ),
    ]
    ancestral = [None, "C", None, "C", None]

    cached_snps = {}

    def _get_snp(self, name, **kwargs):
        """cache standard SNPs"""
        if name in self.cached_snps:
            return self.cached_snps[name]
        snp = list(
            self.human.get_variation(symbol=name, flanks_match_ref=False, **kwargs)
        )[0]
        self.cached_snps[name] = snp
        return snp

    def test_get_distinct(self):
        """should return list of strings"""
        db = Database(
            account=account,
            release=ENSEMBL_RELEASE,
            species="human",
            db_type="variation",
        )
        tn, tc = "variation_feature", "consequence_types"
        expected = set(
            ("3_prime_UTR_variant", "splice_acceptor_variant", "5_prime_UTR_variant")
        )
        got = db.get_distinct(tn, tc)
        self.assertNotEqual(set(got) & expected, set())

    def test_table_has_column(self):
        """return correct values for whether a Table has a column"""
        vardb = Database(
            account=account,
            release=ENSEMBL_RELEASE,
            species="human",
            db_type="variation",
        )

        self.assertTrue(vardb.table_has_column("variation", "evidence_attribs"))
        self.assertFalse(vardb.table_has_column("variation", "validation_status"))

    def test_variant(self):
        """variant attribute correctly constructed"""
        self.assertTrue(len(self.brca2.variants) > 880)

    def test_get_variation_feature(self):
        """should correctly return variation features within a region"""
        snps = self.human.get_features(feature_types="variation", region=self.brca2)
        # snp coordname, start, end should satsify constraints of brca2 loc
        c = 0
        loc = self.brca2.location
        for snp in snps:
            self.assertEqual(snp.location.coord_name, loc.coord_name)
            self.assertTrue(loc.start < snp.location.start < loc.end)
            c += 1
            if c == 2:
                break

    def test_get_variation_by_symbol(self):
        """should return correct snp when query genome by symbol"""
        # supplement this test with some synonymous snp's, where they have no
        # peptide alleles
        for i in range(4):
            snp = self._get_snp(self.snp_names[i])
            self.assertEqual(snp.ancestral, self.ancestral[i], snp)
            self.assertEqual(snp.symbol, self.snp_names[i], snp)
            if isinstance(snp.effect, str):
                effect = set([snp.effect])
            else:
                effect = set(snp.effect)

            self.assertEqual(effect, set(self.snp_effects[i]), snp)
            self.assertEqual(snp.alleles, self.snp_nt_alleles[i], snp)
            self.assertEqual(snp.map_weight, self.map_weights[i], snp)

    def test_somatic_attribute_correct(self):
        """somatic attribute of variants should be correct"""
        symbols_somatic = [("COSM256414", True), ("rs80359189", False)]
        for symbol, expect in symbols_somatic:
            snp = self._get_snp(symbol, somatic=True)
            self.assertEqual(snp.somatic, expect)

    def test_num_alleles(self):
        """should correctly infer the number of alleles"""
        for i in range(4):
            snp = self._get_snp(self.snp_names[i])
            self.assertEqual(len(snp), self.snp_nt_len[i])

    def test_get_peptide_alleles(self):
        """should correctly infer the peptide alleles"""
        for i in range(4):
            snp = self._get_snp(self.snp_names[i])
            if "missense_variant" not in snp.effect:
                continue

            self.assertEqual(snp.peptide_alleles, self.snp_aa_alleles[i])

    def test_no_pep_alleles(self):
        """handle case where coding_sequence_variant has no peptide alleles"""
        snp = self._get_snp("CM033341")
        self.assertTrue(snp.peptide_alleles is None)

    def test_get_peptide_location(self):
        """should return correct location for aa variants"""
        index = self.snp_names.index("rs11545807")
        snp = self._get_snp("rs11545807")
        self.assertEqual(snp.translation_location, 95)

    def test_validation_status(self):
        """should return correct validation status"""

        def func(x):
            if type(x) == str or x is None:
                x = [x]
            return set(x)

        data = (
            ("rs34213141", set(["ESP", "1000Genomes", "Frequency", "ExAC"]), func),
            ("rs12791610", set(["ESP", "1000Genomes", "Frequency", "ExAC"]), func),
            (
                "rs10792769",
                set(["ESP", "1000Genomes", "Frequency", "HapMap", "ExAC"]),
                func,
            ),
            ("rs868440790", set(), func),
        )
        for name, status, conv in data[-1:]:
            snp = self._get_snp(name)
            if snp.validation:
                got = conv(snp.validation)
                self.assertTrue(status & got)
            else:
                # not validated, should return None
                self.assertEqual(snp.validation, None)

    def test_get_flanking_seq(self):
        """should correctly get the flanking sequence if matches reference genome"""

        for i in range(4):  # only have flanking sequence for 3
            snp = self._get_snp(self.snp_names[i])
            self.assertEqual(snp.flanking_seq, self.snp_flanks[i])

    def test_variation_seq(self):
        """should return the sequence for a Variation snp if asked"""
        snp = self._get_snp(self.snp_names[0])
        self.assertIn(str(snp.seq), snp.alleles)

    def test_get_validation_condition(self):
        """simple test of SNP validation status"""
        snp_status = [("rs90", True)]
        for symbol, status in snp_status:
            snp = self._get_snp(symbol)
            self.assertEqual(snp != [], status)

    def test_allele_freqs(self):
        """exercising getting AlleleFreq data"""
        snp = self._get_snp("rs34213141")
        expect = set([("A", "0.0303"), ("G", "0.9697")])
        allele_freqs = snp.allele_freqs
        allele_freqs = set(
            (a, f"{f:.4f}") for a, f in allele_freqs.tolist(["allele", "freq"]) if f
        )
        self.assertTrue(expect.issubset(allele_freqs))

    def test_by_effect(self):
        """excercising select SNPs by effect"""
        for snp in self.human.get_variation(effect="missense_variant", limit=1):
            break

    def test_complex_query(self):
        """only return non-somatic SNPs that are validated and match reference"""
        i = 0
        limit = 4
        for snp in self.human.get_variation(
            effect="missense_variant",
            like=False,
            validated=True,
            somatic=False,
            flanks_match_ref=True,
            limit=limit,
        ):
            self.assertEqual(snp.somatic, False)
            self.assertEqual("missense_variant", snp.effect)
            self.assertNotEqual(snp.flanking_seq, NULL_VALUE)
            self.assertNotEqual(snp.validation, NULL_VALUE)
            i += 1

        self.assertEqual(i, limit)

    def test_get_features_from_nt(self):
        """should correctly return the encompassing gene from 1nt"""
        snp = list(self.human.get_variation(symbol="rs34213141"))[0]
        genes = list(self.human.get_features(feature_types="gene", region=snp))
        self.assertTrue("ENSG00000254997" in [g.stableid for g in genes])

    def test_get_distinct_effect(self):
        """Genome instance get_distinct for SNP effect should work on all genomes"""
        for genome in self.human, self.mouse, self.rat, self.macaq:
            biotypes = genome.get_distinct("effect")


if __name__ == "__main__":
    main()
