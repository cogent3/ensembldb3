import os

from cogent3 import DNA
from cogent3.util.unit_test import TestCase, main

from ensembldb3.host import HostAccount, get_ensembl_account
from ensembldb3.util import convert_strand
from ensembldb3.genome import Genome
from ensembldb3.sequence import _assemble_seq
from ensembldb3.util import asserted_one

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

release = 85

NULL_VALUE = None

# TODO fix flanking sequence issue, a flag somewhere indicating flank is same as reference
# TODO grab the ancestral allele (from variation.ancestral_allele)

if 'ENSEMBL_ACCOUNT' in os.environ:
    args = os.environ['ENSEMBL_ACCOUNT'].split()
    host, username, password = args[0:3]
    kwargs = {}
    if len(args) > 3:
        kwargs['port'] = int(args[3])
    account = HostAccount(host, username, password, **kwargs)
else:
    account = get_ensembl_account(release=release)


class GenomeTestBase(TestCase):
    human = Genome(species="human", release=release, account=account)
    mouse = Genome(species="mouse", release=release, account=account)
    rat = Genome(species="rat", release=release, account=account)
    macaq = Genome(species="macaque", release=release, account=account)
    gorilla = Genome(species="gorilla", release=release, account=account)
    brca2 = human.get_gene_by_stableid(stableid="ENSG00000139618")


class TestGenome(GenomeTestBase):

    def test_other_features(self):
        """should correctly return record for ENSESTG00000000010"""
        est = self.human.get_est_matching(stableid='ENSESTG00000000010')
        direct = list(est)[0]
        ests = self.human.get_features(feature_types='est', coord_name=6,
                                      start=99994000, end=100076519)
        stable_ids = [est.stableid for est in ests]
        self.assertContains(stable_ids, direct.stableid)

    def test_genome_comparison(self):
        """different genome instances with same CoreDb connection are equal"""
        h2 = Genome(species='human', release=release, account=account)
        self.assertEqual(self.human, h2)

    def test_make_location(self):
        """should correctly make a location for an entire chromosome"""
        loc = self.human.make_location(coord_name=1)
        self.assertEqual(len(loc), 248956422)

    def test_get_region(self):
        """should return a generic region that extracts correct sequence"""
        chrom = 1
        start = 11137
        end = start + 20
        region = self.human.get_region(coord_name=chrom, start=start, end=end,
                                      ensembl_coord=True)
        self.assertEqual(region.location.start, start - 1)
        self.assertEqual(region.location.end, end)
        self.assertEqual(region.location.coord_name, str(chrom))
        self.assertEqual(region.location.coord_type, 'chromosome')
        self.assertEqual(region.seq, 'ACCTCAGTAATCCGAAAAGCC')

    def test_get_assembly_exception_region(self):
        """should return correct sequence for region with an assembly
        exception"""
        region = self.human.get_region(coord_name="Y", start=57211873,
                                      end=57211894, strand=1, ensembl_coord=True)

        self.assertEqual(str(region.seq), 'CGAGGACGACTGGGAATCCTAG')

    def test_no_assembly(self):
        """return N's for coordinates with no assembly"""
        krat = Genome('Kangaroo rat', release=85)
        start = 24385
        end = start + 100
        region = krat.get_region(coord_name='scaffold_13754', start=start,
                                end=end)
        self.assertEqual(str(region.seq), 'N' * (end - start))

    def test_getting_annotated_seq(self):
        """a region should return a sequence with the correct annotation"""
        new_loc = self.brca2.location.resized(-100, 100)
        region = self.human.get_region(region=new_loc)
        annot_seq = region.get_annotated_seq(feature_types='gene')
        gene_annots = annot_seq.get_annotations_matching('gene')
        self.assertEqual(gene_annots[0].name, self.brca2.symbol)

    def test_correct_feature_type_id_cache(self):
        """should obtain the feature type identifiers without failure"""
        self.assertNotEqual(self.human._feature_type_ids.CpGisland, None)

    def test_strand_conversion(self):
        """should consistently convert strand info"""
        self.assertEqual(convert_strand(None), 1)
        self.assertEqual(convert_strand(-1), -1)
        self.assertEqual(convert_strand(1), 1)
        self.assertEqual(convert_strand('-'), -1)
        self.assertEqual(convert_strand('+'), 1)
        self.assertEqual(convert_strand(-1.0), -1)
        self.assertEqual(convert_strand(1.0), 1)

    def test_pool_connection(self):
        """excercising ability to specify pool connection"""
        dog = Genome(species="dog", release=release, account=account,
                     pool_recycle=1000)

    def test_gorilla(self):
        """should correctly return a gorilla gene"""
        self.gorilla = Genome(
            species="gorilla", release=release, account=account)
        gene = self.gorilla.get_gene_by_stableid('ENSGGOG00000005730')
        self.assertEqual(str(gene.seq[:10]), 'TGGGAGTCCA')

    def test_diff_strand_contig_chrom(self):
        """get correct sequence when contig and chromosome strands differ"""
        gene = self.gorilla.get_gene_by_stableid('ENSGGOG00000001953')
        cds = gene.canonical_transcript.cds
        self.assertEqual(str(cds), 'ATGGCCCAGGATCTCAGCGAGAAGGACCTGTTGAAGATG'
                         'GAGGTGGAGCAGCTGAAGAAAGAAGTGAAAAACACAAGAATTCCGATTTCCAAAGCGGGAAAGGAAAT'
                         'CAAAGAGTACGTGGAGGCCCAAGCAGGAAACGATCCTTTTCTCAAAGGCATCCCTGAGGACAAGAATC'
                         'CCTTCAAGGAGAAAGGTGGCTGTCTGATAAGCTGA')

    def test_get_distinct_biotype(self):
        """Genome instance get_distinct for biotype should work on all genomes"""
        for genome in self.gorilla, self.human, self.mouse, self.rat, self.macaq:
            biotypes = genome.get_distinct('biotype')


class TestGene(GenomeTestBase):

    def _eval_brca2(self, brca2):
        """tests all attributes correctly created"""
        self.assertEqual(brca2.symbol.lower(), 'brca2')
        self.assertEqual(brca2.stableid, 'ENSG00000139618')
        self.assertEqual(brca2.biotype.lower(), 'protein_coding')
        self.assertContains(brca2.description.lower(), 'dna repair associated')
        self.assertEqual(brca2.status, 'KNOWN')
        self.assertEqual(brca2.canonical_transcript.stableid,
                         'ENST00000380152')
        # note length can change between genome builds
        self.assertGreaterThan(len(brca2), 83700)
        transcript = brca2.get_member('ENST00000544455')
        self.assertEqual(transcript.get_cds_length(), len(transcript.cds))

    def test_get_genes_by_stable_id(self):
        """if get gene by stable_id, attributes should be correctly
        constructed"""
        self._eval_brca2(self.brca2)

    def test_get_exons(self):
        """transcript should return correct exons for brca2"""
        transcript = self.brca2.get_member('ENST00000380152')
        self.assertEqual(len(transcript.translated_exons), 26)
        self.assertEqual(len(transcript.cds), 3419 * 3)
        self.assertEqual(len(transcript.protein_seq), 3418)

    def test_translated_exons(self):
        """should correctly translate a gene with 2 exons but 1st exon
        transcribed"""
        gene = self.mouse.get_gene_by_stableid(stableid='ENSMUSG00000036136')
        transcript = gene.get_member('ENSMUST00000041133')
        self.assertTrue(len(transcript.protein_seq) > 0)
        # now one on the - strand
        gene = self.mouse.get_gene_by_stableid(stableid='ENSMUSG00000045912')
        transcript = gene.transcripts[0]
        self.assertTrue(len(transcript.protein_seq) > 0)

    def test_failed_ensembl_annotation(self):
        """we demonstrate a failed annotation by ensembl"""
        # I'm including this to demonstrate that Ensembl coords are
        # complex. This case has a macaque gene which we correctly
        # infer the CDS boundaries for according to Ensembl, but the CDS
        # length is not divisible by 3.
        gene = self.macaq.get_gene_by_stableid(stableid='ENSMMUG00000001551')
        transcript = gene.get_member('ENSMMUT00000002194')
        # the following works because we enforce the length being divisble by 3
        # in producing protein_seq
        prot_seq = transcript.protein_seq
        # BUT if you work off the cds you will need to slice the CDS to be
        # divisible by 3 to get the same protein sequence
        l = transcript.get_cds_length()
        trunc_cds = transcript.cds[: l - (l % 3)]
        prot_seq = trunc_cds.get_translation()
        self.assertEqual(str(prot_seq),
                         'MPSSPLRVAVVCSSNQNRSMEAHNILSKRGFSVRSFGTGTHVKLPGPAPDKPNVYDFKTT'
                         'YDQMYNDLLRKDKELYTQNGILHMLDRNKRIKPRPERFQNCKDLFDLILTCEERVY')

    def test_exon_phases(self):
        """correctly identify phase for an exon"""
        stable_id = 'ENSG00000171408'
        gene = self.human.get_gene_by_stableid(stableid=stable_id)
        exon1 = gene.transcripts[1].exons[0]
        # first two bases of codon missing
        self.assertEqual(exon1.phase_start, 2)
        # last two bases of codon missing
        self.assertEqual(exon1.phase_end, 1)
        # can translate the sequence if we take those into account
        seq = exon1.seq[1:-1].get_translation()
        self.assertEqual(str(seq), 'HMLSKVGMWDFDIFLFDRLTN')

    def test_cds_from_outofphase(self):
        """return a translatable cds sequence from out-of-phase start"""
        # canonical transcript phase end_phase
        # ENSG00000111729 ENST00000229332 -1 -1
        # ENSG00000177151 ENST00000317450 0 -1
        # ENSG00000249624 ENST00000433395 1 -1
        # ENSG00000237276 ENST00000442385 2 -1

        canon_ids = 'ENSG00000111729 ENSG00000177151 ENSG00000237276 ENSG00000251184'.split()
        for index, stable_id in enumerate(canon_ids):
            gene = self.human.get_gene_by_stableid(stableid=stable_id)
            transcript = gene.canonical_transcript
            prot_seq = transcript.protein_seq

    def test_gene_transcripts(self):
        """should return multiple transcripts"""
        stable_id = 'ENSG00000012048'
        gene = self.human.get_gene_by_stableid(stableid=stable_id)
        self.assertTrue(len(gene.transcripts) > 1)
        # .. and correctly construct the cds and location
        for transcript in gene.transcripts:
            self.assertTrue(transcript.get_cds_length() > 0)
            self.assertEqual(transcript.location.coord_name, '17')

    def test_get_longest_cds_transcript2(self):
        """should correctly return transcript with longest cds"""
        # ENSG00000123552 is protein coding, ENSG00000206629 is ncRNA
        for stable_id, max_cds_length in [('ENSG00000123552', 2445),
                                          ('ENSG00000206629', 164)]:
            gene = self.human.get_gene_by_stableid(stableid=stable_id)
            ts = gene.get_longest_cds_transcript()
            self.assertEqual(len(ts.cds), max_cds_length)
            self.assertEqual(ts.get_cds_length(), max(gene.get_cds_lengths()))

    def test_get_longest_cds_transcript1(self):
        """should correctly return transcript with longest cds"""
        stable_id = 'ENSG00000178591'
        gene = self.human.get_gene_by_stableid(stableid=stable_id)
        ts = gene.get_longest_cds_transcript()
        self.assertEqual(ts.get_cds_length(), max(gene.get_cds_lengths()))

    def test_rna_transcript_cds(self):
        """should return a cds for an RNA gene too"""
        rna_gene = self.human.get_gene_by_stableid(stableid='ENSG00000210049')
        self.assertTrue(rna_gene.transcripts[0].get_cds_length() > 0)

    def test_gene_annotation(self):
        """should correctly annotated a sequence"""
        annot_seq = self.brca2.get_annotated_seq(feature_types='gene')
        gene_annots = annot_seq.get_annotations_matching('gene')
        self.assertEqual(gene_annots[0].name, self.brca2.symbol)

    def test_get_by_symbol(self):
        """selecting a gene by it's HGNC symbol should correctly populate all
        specified attributes"""
        results = self.human.get_genes_matching(symbol="BRCA2")
        found = False
        for gene in results:
            if gene.stableid == 'ENSG00000139618':
                self._eval_brca2(gene)
                found = True
        self.assertTrue(found)

    def test_get_by_symbol_synonym(self):
        """return correct gene if provide a synonym, rather than symbol"""
        synonym = 'FOXO1A'
        gene = list(self.human.get_genes_matching(symbol=synonym))[0]
        self.assertEqual(gene.symbol, 'FOXO1')

    def test_get_by_description(self):
        """if get by description, all attributes should be correctly
        constructed"""
        description = 'brca2'
        results = list(self.human.get_genes_matching(description=description))
        self._eval_brca2(results[0])

    def test_get_member(self):
        """should return correct exon and translated exon"""
        transcript = self.brca2.get_member('ENST00000380152')
        # just returns the first
        exon_id = 'ENSE00001484009'
        exon = transcript.get_member(exon_id)
        trans_exon = transcript.get_member(exon_id, 'translated_exons')
        self.assertEqual(exon.stableid, exon_id)
        self.assertEqual(trans_exon.stableid, exon_id)
        # we check we got Exon in the first call and TranslatedExon in the
        # second using the fact that the exons entry is longer than the
        # translated_exons one
        self.assertGreaterThan(len(exon), len(trans_exon))

    def test_get_by_biotype(self):
        results = list(self.human.get_genes_matching(
            biotype='Mt_tRNA', like=False))
        self.assertEqual(len(results), 22)
    
    def test_limit_genes_matching(self):
        """limit argument should work for get_genes_matching"""
        # use the limit argument
        results = list(self.human.get_genes_matching(
            biotype='protein_coding', limit=10))
        self.assertEqual(len(results), 10)

    def test_get_by_decsr_biotype(self):
        """combining the description and biotype should return a result"""
        results = list(self.human.get_genes_matching(biotype="protein_coding",
                                                   description="cancer"))
        self.assertTrue(len(results) > 50)

    def test_get_gene_by_stable_id(self):
        """should correctly handle getting gene by stable_id"""
        stable_id = 'ENSG00000012048'
        gene = self.human.get_gene_by_stableid(stableid=stable_id)
        self.assertEqual(gene.stableid, stable_id)

        # if invalid stable_id, should just return None
        stable_id = 'ENSG00000XXXXX'
        gene = self.human.get_gene_by_stableid(stableid=stable_id)
        self.assertEqual(gene, None)

    def test_get_transcript_by_stable_id(self):
        """should correctly handle getting transcript by stable_id"""
        # if invalid stable_id, should just return None
        stable_id = 'ENST00000XXXXX'
        transcript = self.human.get_transcript_by_stableid(stableid=stable_id)
        self.assertEqual(transcript, None)

        # get transcript via gene and check values match
        stable_id = 'ENST00000380152'
        transcript = self.human.get_transcript_by_stableid(stableid=stable_id)
        self.assertEqual(transcript.stableid, stable_id)
        gene = transcript.gene
        brca2 = self.human.get_gene_by_stableid(stableid='ENSG00000139618')
        self.assertEqual(brca2.canonical_transcript.stableid,
                         transcript.stableid)
        self.assertEqual(
            brca2.canonical_transcript.get_cds_length(), len(transcript.cds))
        self.assertEqual(str(brca2.canonical_transcript.cds),
                         str(transcript.cds))
        self.assertEqual(str(brca2.canonical_transcript.cds),
                         str(transcript.cds))
        self.assertEqual(str(brca2.canonical_transcript.seq),
                         str(transcript.seq))
        self.assertEqual(brca2.stableid, gene.stableid)
        self.assertEqual(brca2.seq, gene.seq)

    def test_gene_on_transcript(self):
        """Transcript instances Gene attribute should be complete"""
        brca2 = self.human.get_gene_by_stableid(stableid='ENSG00000139618')
        transcript = self.human.get_transcript_by_stableid(
            stableid='ENST00000380152')
        self.assertEqual(transcript.gene.symbol, brca2.symbol)

    def test_intron_number(self):
        """number of introns should be correct"""
        for gene_id, transcript_id, exp_number in [
            ('ENSG00000227268', 'ENST00000445946', 0),
            ('ENSG00000132199', 'ENST00000583771', 5),
            ('ENSG00000132199', 'ENST00000340116', 14)]:
            gene = asserted_one(self.human.get_genes_matching(stableid=gene_id))
            transcript = asserted_one(
                [t for t in gene.transcripts if t.stableid == transcript_id])
            if exp_number == 0:
                self.assertEqual(transcript.introns, None)
            else:
                self.assertEqual(len(transcript.introns), exp_number)

    def test_intron(self):
        """should get correct Intron sequence, regardless of strand"""
        # IL2 is on - strand, IL13 is on + strand, both have three introns
        IL2_exp_introns = [
            (1, 122456203, 122456293, 'gtaagtatat', 'actttcttag'),
            (2, 122453853, 122456143, 'gtaagtacaa', 'attattctag'),
            (3, 122451862, 122453709, 'gtaaggcatt', 'tcttttatag')]
        IL13_exp_introns = [
            (1, 132658360, 132659417, 'gtgagtgtcg', 'gctcccacag'),
            (2, 132659471, 132659723, 'gtaaggacct', 'ctccccacag'),
            (3, 132659828, 132660174, 'gtaaggcatc', 'tgtcctgcag')]

        for symbol, stable_id, exp_introns in [
            ('IL2', 'ENST00000226730', IL2_exp_introns),
            ('IL13', 'ENST00000304506', IL13_exp_introns)]:
            gene = asserted_one(self.human.get_genes_matching(symbol=symbol))
            strand = gene.location.strand
            transcript = asserted_one(
                [t for t in gene.transcripts if t.stableid == stable_id])
            introns = transcript.introns
            self.assertEqual(len(introns), len(exp_introns))
            idx = 0
            for intron in introns:
                loc = intron.location
                start, end = loc.start, loc.end
                seq = str(intron.seq)
                exp_rank, exp_start, exp_end, exp_seq5, \
                exp_seq3 = exp_introns[idx]
                self.assertEqual(loc.strand, strand)
                # test the order using rank
                self.assertEqual(intron.rank, exp_rank)
                # test position
                self.assertEqual(start, exp_start)
                self.assertEqual(end, exp_end)
                # test sequence
                self.assertEqual(seq[:10], exp_seq5.upper())
                self.assertEqual(seq[-10:], exp_seq3.upper())
                idx += 1

    def test_intron_annotation(self):
        """sequences annotated with introns should return correct seq"""
        for symbol, stable_id, rank, exp_seq5, exp_seq3 in [
                ('IL2', 'ENST00000226730', 1, 'gtaagtatat', 'actttcttag'),
                ('IL13', 'ENST00000304506', 3, 'gtaaggcatc', 'tgtcctgcag')]:
            gene = asserted_one(self.human.get_genes_matching(symbol=symbol))
            seq = gene.get_annotated_seq(feature_types='gene')
            intron = asserted_one(seq.get_annotations_matching('intron',
                                                             '%s-%d' % (stable_id, rank)))
            intron_seq = str(seq.get_region_covering_all(intron).get_slice())
            self.assertEqual(intron_seq[:10], exp_seq5.upper())
            self.assertEqual(intron_seq[-10:], exp_seq3.upper())



class TestFeatures(GenomeTestBase):

    def setUp(self):
        self.igf2 = self.human.get_gene_by_stableid(stableid='ENSG00000167244')

    def test_CpG_island(self):
        """should return correct CpG islands"""
        CpGislands = self.human.get_features(region=self.igf2,
                                            feature_types='CpG')
        expected_stats = [(630, 757), (652, 537), (3254, 3533)]
        obs_stats = [(int(island.Score), len(island))
                     for island in CpGislands]
        obs_stats.sort()
        self.assertTrue(set(expected_stats) & set(obs_stats) != set())

    def test_get_multiple_features(self):
        """should not fail to get multiple feature types"""
        regions =\
            self.human.get_features(feature_types=['repeat', 'gene', 'cpg'],
                                   coord_name=1, start=869936, end=901867,
                                   limit=5)
        for region in regions:
            pass

    def test_repeats(self):
        """should correctly return a repeat"""
        loc = self.igf2.location.resized(-1000, 1000)
        repeats = list(self.human.get_features(
            region=loc, feature_types='repeat'))
        self.assertTrue(len(repeats) >= 4)

    def test_genes(self):
        """should correctly identify igf2 within a region"""
        loc = self.igf2.location.resized(-1000, 1000)
        genes = self.human.get_features(region=loc, feature_types='gene')
        symbols = [g.symbol.lower() for g in genes]
        self.assertContains(symbols, self.igf2.symbol.lower())

    def test_other_genes(self):
        """docstring for est_other_genes"""
        mouse = self.mouse.get_region(coord_name='5', start=150791005,
                                     end=150838512, strand='-')
        rat = self.rat.get_region(coord_name='12', start=4282534, end=4324019,
                                 strand='+')
        for region in [mouse, rat]:
            features = region.get_features(feature_types=['gene'])
            ann_seq = region.get_annotated_seq(feature_types='gene')
            genes = ann_seq.get_annotations_matching('gene')
            self.assertTrue(genes != [])


    def test_gene_feature_data_correct(self):
        """should apply gene feature data in a manner consistent with strand
        and the Cogent sequence annotations slice should return the same
        result"""
        plus = list(self.human.get_features(feature_types='gene',
                                           coord_name=13,
                                           start=31787610,
                                           end=31871820))[0]
        minus = plus.location.copy()
        minus.strand *= -1
        minus = self.human.get_region(region=minus)
        # get Sequence
        plus_seq = plus.get_annotated_seq(feature_types='gene')
        minus_seq = minus.get_annotated_seq(feature_types='gene')
        # the seqs should be the rc of each other
        self.assertEqual(str(plus_seq), str(minus_seq.rc()))
        # the cds, however, from the annotated sequences should be identical
        plus_cds = plus_seq.get_annotations_matching('CDS')[0]
        minus_cds = minus_seq.get_annotations_matching('CDS')[0]
        self.assertEqual(str(plus_cds.get_slice()), str(minus_cds.get_slice()))

    def test_other_feature_data_correct(self):
        """should apply CpG feature data in a manner consistent with strand"""
        human = self.human
        coord = dict(coord_name=11, start=2143894, end=2144494)
        exp_coord = dict(coord_name=11, start=2143906, end=2144442)
        exp_loc = human.get_region(strand=1, ensembl_coord=True, **exp_coord)
        exp = exp_loc.seq

        ps_feat = human.get_region(strand=1, **coord)
        ms_feat = human.get_region(strand=-1, **coord)

        ps_seq = ps_feat.get_annotated_seq(feature_types='CpG')
        ps_cgi = ps_seq.get_annotations_matching('CpGisland')[0]

        self.assertEqual(ps_feat.seq, ms_feat.seq.rc())

        self.assertEqual(ps_cgi.get_slice().rc(), exp)
        ms_seq = ms_feat.get_annotated_seq(feature_types='CpG')
        ms_cgi = ms_seq.get_annotations_matching('CpGisland')[0]

        self.assertEqual(ms_cgi.get_slice(), ps_cgi.get_slice())

    def test_other_repeat(self):
        """should apply repeat feature data in a manner consistent with strand"""
        coord = dict(coord_name=13, start=32316063, end=32316363)
        # 13:32316063 -32316363
        ps_repeat = self.human.get_region(strand=1, **coord)
        ms_repeat = self.human.get_region(strand=-1, **coord)
        # note this MER3 repeat is annotated on the -1 strand
        exp = DNA.make_seq('AGCTTACTGTGAGGATGGGAACATTTTACAGCTGTGCTGTCCAAA'
                               'CCGGTGCCACTAGCCACATTAAGCACTCGAAACGTGGCTAGTGCGACTAGAGAAGAGGAT'
                               'TTTCATACGATTTAGTTTCAATCACGCTAACCAGTGACGCGTGGCTAGTGG')

        self.assertEqual(ms_repeat.seq, ps_repeat.seq.rc())

        ps_annot_seq = ps_repeat.get_annotated_seq(feature_types='repeat')
        ms_annot_seq = ms_repeat.get_annotated_seq(feature_types='repeat')
        ps_seq = ps_annot_seq.get_annotations_matching('repeat')[0]
        ms_seq = ms_annot_seq.get_annotations_matching('repeat')[0]
        self.assertEqual(ms_seq.get_slice(), ps_seq.get_slice())
        self.assertEqual(ps_seq.get_slice(), exp)



class TestAssembly(TestCase):

    def test_assemble_seq(self):
        """should correctly fill in a sequence with N's"""
        expect = DNA.make_seq("NAAAAANNCCCCCNNGGGNNN")
        frags = ["AAAAA", "CCCCC", "GGG"]
        positions = [(11, 16), (18, 23), (25, 28)]
        self.assertEqual(_assemble_seq(frags, 10, 31, positions), expect)
        positions = [(1, 6), (8, 13), (15, 18)]
        self.assertEqual(_assemble_seq(frags, 0, 21, positions), expect)
        # should work with:
        # start matches first frag start
        expect = DNA.make_seq("AAAAANNCCCCCNNGGGNNN")
        positions = [(0, 5), (7, 12), (14, 17)]
        self.assertEqual(_assemble_seq(frags, 0, 20, positions), expect)
        # end matches last frag_end
        expect = DNA.make_seq("NAAAAANNCCCCCNNGGG")
        positions = [(11, 16), (18, 23), (25, 28)]
        self.assertEqual(_assemble_seq(frags, 10, 28, positions), expect)
        # both start and end matched
        expect = DNA.make_seq("AAAAANNCCCCCNNGGG")
        positions = [(10, 15), (17, 22), (24, 27)]
        self.assertEqual(_assemble_seq(frags, 10, 27, positions), expect)
        # one frag
        expect = DNA.make_seq(''.join(frags))
        positions = [(10, 23)]
        self.assertEqual(_assemble_seq([''.join(frags)], 10, 23, positions),
                         expect)

if __name__ == "__main__":
    main()
