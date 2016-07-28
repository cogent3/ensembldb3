import os

from cogent3.util.unit_test import TestCase, main

from ensembldb.host import HostAccount, get_ensembl_account
from ensembldb.assembly import Coordinate, CoordSystem, \
    get_coord_conversion
from ensembldb.genome import Genome

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "3.0.alpha"
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

human = Genome(Species='human', release=release, account=account)
platypus = Genome(Species='platypus', release=release, account=account)


class TestLocation(TestCase):

    def test_init(self):
        human_loc = Coordinate(CoordName='x', start=1000, end=10000, strand=-1,
                               genome=human)
        # TODO: complete test for platpus
        self.assertEqual(human_loc.CoordType, 'chromosome')
        self.assertEqual(human_loc.CoordName, 'x')
        self.assertEqual(human_loc.start, 1000)
        self.assertEqual(human_loc.end, 10000)
        self.assertEqual(human_loc.strand, -1)
        self.assertEqual(human_loc.Species, "Homo sapiens")
        self.assertEqual(human_loc.seq_region_id, 131539)

    def test_get_coord_conversion(self):
        """should correctly map between different coordinate levels"""
        # not really testing the contig coordinates are correct
        CoordName, start, end, strand = '1', 1000, 1000000, 1
        human_loc = Coordinate(CoordName=CoordName, start=start, end=end,
                               strand=strand, genome=human)
        results = get_coord_conversion(human_loc, 'contig', human.CoreDb)
        for result in results:
            self.assertTrue(result[0].CoordName == CoordName)
            self.assertTrue(result[0].start >= start)
            self.assertTrue(result[0].end <= end)
            self.assertTrue(result[0].strand == strand)

    def test_coord_shift(self):
        """adding coordinates should produce correct results"""
        CoordName, start, end, strand = '1', 1000, 1000000, 1
        loc1 = Coordinate(CoordName=CoordName, start=start, end=end,
                          strand=strand, genome=human)
        for shift in [100, -100]:
            loc2 = loc1.shifted(shift)
            self.assertEqual(loc2.start, loc1.start + shift)
            self.assertEqual(loc2.end, loc1.end + shift)
            self.assertEqual(id(loc1.genome), id(loc2.genome))
        self.assertNotEqual(id(loc1), id(loc2))

    def test_coord_resize(self):
        """resizing should work"""
        CoordName, start, end, strand = '1', 1000, 1000000, 1
        loc1 = Coordinate(CoordName=CoordName, start=start, end=end,
                          strand=strand, genome=human)
        front_shift = -100
        back_shift = 100
        loc2 = loc1.resized(front_shift, back_shift)
        self.assertEqual(len(loc2), len(loc1) + 200)
        self.assertEqual(loc2.start, loc1.start + front_shift)
        self.assertEqual(loc2.end, loc1.end + back_shift)
        self.assertEqual(loc1.strand, loc2.strand)

    def test_adopted(self):
        """coordinate should correctly adopt seq_region_id properties of 
        provided coordinate"""
        CoordName, start, end, strand = '1', 1000, 1000000, 1
        c1 = Coordinate(CoordName=CoordName, start=start, end=end,
                        strand=strand, genome=human)
        CoordName, start, end, strand = '2', 2000, 2000000, 1
        c2 = Coordinate(CoordName=CoordName, start=start, end=end,
                        strand=strand, genome=human)
        c3 = c1.adopted(c2)
        self.assertEqual(c3.CoordName, c2.CoordName)
        self.assertEqual(c3.CoordType, c2.CoordType)
        self.assertEqual(c3.seq_region_id, c2.seq_region_id)
        self.assertEqual(c3.start, c1.start)
        self.assertEqual(c3.end, c1.end)
        self.assertEqual(c3.strand, c1.strand)
        c3 = c1.adopted(c2, shift=100)
        self.assertEqual(c3.start, c1.start + 100)
        self.assertEqual(c3.end, c1.end + 100)


class TestCoordSystem(TestCase):

    def test_call(self):
        human_chrom = CoordSystem('chromosome', core_db=human.CoreDb,
                                  species='human')
        human_contig = CoordSystem(1, species='human')
        self.assertEqual(human_chrom.coord_system_id, 4)
        self.assertEqual(human_contig.name, 'contig')
        self.assertEqual(human_contig.attr, 'default_version, sequence_level')


if __name__ == '__main__':
    main()
