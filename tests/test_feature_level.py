import os

from unittest import TestCase, main

from ensembldb3.feature_level import FeatureCoordLevels
from ensembldb3.genome import Genome
from ensembldb3.host import HostAccount, get_ensembl_account

from . import ENSEMBL_RELEASE


__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
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


class TestFeatureCoordLevels(TestCase):
    def setUp(self):
        self.chicken = Genome(
            species="chicken", release=ENSEMBL_RELEASE, account=account
        )

    def test_feature_levels(self):
        ChickenFeatureLevels = FeatureCoordLevels("chicken")
        chicken_feature_levels = ChickenFeatureLevels(
            feature_types=["gene", "cpg", "est"],
            core_db=self.chicken.CoreDb,
            otherfeature_db=self.chicken.OtherFeaturesDb,
        )
        self.assertEqual(
            chicken_feature_levels["repeat"].levels, ["chromosome", "scaffold"]
        )
        self.assertEqual(
            set(chicken_feature_levels["cpg"].levels), set(["chromosome", "scaffold"])
        )

    def test_repeat(self):
        # use chicken genome as it need to do conversion
        # chicken coordinate correspondent toRefSeq human IL2A region
        coord = dict(coord_name=9, start=21729000, end=21730141)
        region = self.chicken.get_region(**coord)
        # repeat is recorded at contig level, strand is 0
        repeats = region.get_features(feature_types="repeat")
        expect = [
            ("9", 21729713, 21729732),
            ("9", 21729956, 21729968),
            ("9", 21730047, 21730076),
        ]
        obs = []
        for repeat in repeats:
            loc = repeat.location
            obs.append((str(loc.coord_name), loc.start, loc.end))
        self.assertEqual(set(obs), set(expect))

    def test_cpg(self):
        # contain 3 CpG island recorded at chromosome level
        coord1 = dict(coord_name=26, start=105184, end=184346)
        cpgs1 = self.chicken.get_features(feature_types="cpg", **coord1)
        exp = [("26", 116624, 117610), ("26", 138598, 139523), ("26", 183375, 184708)]
        obs = []
        for cpg in cpgs1:
            loc = cpg.location
            obs.append((str(loc.coord_name), loc.start, loc.end))
        self.assertEqual(set(obs), set(exp))

        # test cpg features record at scaffold level:
        coord2 = dict(coord_name="KQ759405.1", start=1, end=212434)
        cpgs2 = self.chicken.get_features(feature_types="cpg", **coord2)
        self.assertEqual(len(list(cpgs2)), 18)


if __name__ == "__main__":
    main()
