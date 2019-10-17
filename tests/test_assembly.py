import os

from cogent3.util.unit_test import TestCase, main
from ensembldb3.assembly import Coordinate, CoordSystem, get_coord_conversion
from ensembldb3.genome import Genome
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

human = Genome(species="human", release=ENSEMBL_RELEASE, account=account)
platypus = Genome(species="platypus", release=ENSEMBL_RELEASE, account=account)


class TestLocation(TestCase):
    def test_init(self):
        human_loc = Coordinate(
            coord_name="x", start=1000, end=10000, strand=-1, genome=human
        )
        # TODO: complete test for platpus
        self.assertEqual(human_loc.coord_type, "chromosome")
        self.assertEqual(human_loc.coord_name, "x")
        self.assertEqual(human_loc.start, 1000)
        self.assertEqual(human_loc.end, 10000)
        self.assertEqual(human_loc.strand, -1)
        self.assertEqual(human_loc.species, "Homo sapiens")
        self.assertEqual(human_loc.seq_region_id, 131539)

    def test_get_coord_conversion(self):
        """should correctly map between different coordinate levels"""
        # not really testing the contig coordinates are correct
        coord_name, start, end, strand = "1", 1000, 1000000, 1
        human_loc = Coordinate(
            coord_name=coord_name, start=start, end=end, strand=strand, genome=human
        )
        results = get_coord_conversion(human_loc, "contig", human.CoreDb)
        for result in results:
            self.assertTrue(result[0].coord_name == coord_name)
            self.assertTrue(result[0].start >= start)
            self.assertTrue(result[0].end <= end)
            self.assertTrue(result[0].strand == strand)

    def test_coord_shift(self):
        """adding coordinates should produce correct results"""
        coord_name, start, end, strand = "1", 1000, 1000000, 1
        loc1 = Coordinate(
            coord_name=coord_name, start=start, end=end, strand=strand, genome=human
        )
        for shift in [100, -100]:
            loc2 = loc1.shifted(shift)
            self.assertEqual(loc2.start, loc1.start + shift)
            self.assertEqual(loc2.end, loc1.end + shift)
            self.assertEqual(id(loc1.genome), id(loc2.genome))
        self.assertNotEqual(id(loc1), id(loc2))

    def test_coord_resize(self):
        """resizing should work"""
        coord_name, start, end, strand = "1", 1000, 1000000, 1
        loc1 = Coordinate(
            coord_name=coord_name, start=start, end=end, strand=strand, genome=human
        )
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
        coord_name, start, end, strand = "1", 1000, 1000000, 1
        c1 = Coordinate(
            coord_name=coord_name, start=start, end=end, strand=strand, genome=human
        )
        coord_name, start, end, strand = "2", 2000, 2000000, 1
        c2 = Coordinate(
            coord_name=coord_name, start=start, end=end, strand=strand, genome=human
        )
        c3 = c1.adopted(c2)
        self.assertEqual(c3.coord_name, c2.coord_name)
        self.assertEqual(c3.coord_type, c2.coord_type)
        self.assertEqual(c3.seq_region_id, c2.seq_region_id)
        self.assertEqual(c3.start, c1.start)
        self.assertEqual(c3.end, c1.end)
        self.assertEqual(c3.strand, c1.strand)
        c3 = c1.adopted(c2, shift=100)
        self.assertEqual(c3.start, c1.start + 100)
        self.assertEqual(c3.end, c1.end + 100)

    def test_union(self):
        """union of coordinates correct for +/- strands"""
        c1 = Coordinate(coord_name="1", start=50, end=100, strand=1, genome=human)
        # None returned when different coord_name
        c2 = Coordinate(coord_name="2", start=2000, end=2000000, strand=1, genome=human)
        self.assertEqual(c1.union(c2), None)
        # None returned when different strand
        c2 = Coordinate(coord_name="1", start=50, end=300, strand=-1, genome=human)
        self.assertEqual(c1.union(c2), None)
        # None returned when different species
        c2 = Coordinate(coord_name="1", start=50, end=300, strand=1, genome=platypus)
        self.assertEqual(c1.union(c2), None)

        # correct span when '+' strand
        c2 = Coordinate(coord_name="1", start=120, end=200, strand=1, genome=human)
        c3 = c1.union(c2)
        self.assertEqual((c3.start, c3.end, c3.strand), (50, 200, 1))

        # correct span when '-' strand
        c1 = Coordinate(coord_name="1", start=500, end=100, strand=-1, genome=human)
        c2 = Coordinate(coord_name="1", start=600, end=500, strand=-1, genome=human)
        c3 = c1.union(c2)
        self.assertEqual((c3.start, c3.end, c3.strand), (100, 600, -1))


class TestCoordSystem(TestCase):
    def test_call(self):
        human_chrom = CoordSystem("chromosome", core_db=human.CoreDb, species="human")
        human_contig = CoordSystem(1, species="human")
        self.assertEqual(human_chrom.coord_system_id, 4)
        self.assertEqual(human_contig.name, "contig")
        self.assertEqual(human_contig.attr, "default_version, sequence_level")


if __name__ == "__main__":
    main()
