from unittest import TestCase, main

from cogent3.util.table import Table

from ensembldb3.species import Species

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


class TestSpeciesNamemaps(TestCase):
    def test_get_name_type(self):
        """should return the (latin|common) name given a latin, common or ensembl
        db prefix names"""
        self.assertEqual(Species.get_species_name("human"), "Homo sapiens")
        self.assertEqual(Species.get_species_name("homo_sapiens"), "Homo sapiens")
        self.assertEqual(
            Species.get_species_name("canis_lupus_familiaris"), "Canis lupus familiaris"
        )
        self.assertEqual(Species.get_common_name("Mus musculus"), "Mouse")
        self.assertEqual(Species.get_common_name("mus_musculus"), "Mouse")

    def test_get_ensembl_format(self):
        """should take common or latin names and return the corresponding
        ensembl db prefix"""
        self.assertEqual(Species.get_ensembl_db_prefix("human"), "homo_sapiens")
        self.assertEqual(Species.get_ensembl_db_prefix("mouse"), "mus_musculus")
        self.assertEqual(Species.get_ensembl_db_prefix("Mus musculus"), "mus_musculus")
        self.assertEqual(
            Species.get_ensembl_db_prefix("Canis lupus familiaris"),
            "canis_lupus_familiaris",
        )
        self.assertEqual(
            Species.get_ensembl_db_prefix("Canis familiaris"), "canis_familiaris"
        )

    def test_add_new_species(self):
        """should correctly add a new species/common combination and infer the
        correct ensembl prefix"""
        species_name, common_name = "Otolemur garnettii", "Bushbaby"
        Species.amend_species(species_name, common_name)
        self.assertEqual(Species.get_species_name(species_name), species_name)
        self.assertEqual(Species.get_species_name("Bushbaby"), species_name)
        self.assertEqual(Species.get_species_name(common_name), species_name)
        self.assertEqual(Species.get_common_name(species_name), common_name)
        self.assertEqual(Species.get_common_name("Bushbaby"), common_name)
        self.assertEqual(
            Species.get_ensembl_db_prefix("Bushbaby"), "otolemur_garnettii"
        )
        self.assertEqual(
            Species.get_ensembl_db_prefix(species_name), "otolemur_garnettii"
        )
        self.assertEqual(
            Species.get_ensembl_db_prefix(common_name), "otolemur_garnettii"
        )

    def test_amend_existing(self):
        """should correctly amend an existing species"""
        species_name = "Ochotona princeps"
        common_name1 = "american pika"
        common_name2 = "pika"
        ensembl_pref = "ochotona_princeps"
        Species.amend_species(species_name, common_name1)
        self.assertEqual(Species.get_common_name(species_name), common_name1)
        Species.amend_species(species_name, common_name2)
        self.assertEqual(Species.get_species_name(common_name2), species_name)
        self.assertEqual(Species.get_species_name(ensembl_pref), species_name)
        self.assertEqual(Species.get_common_name(species_name), common_name2)
        self.assertEqual(Species.get_common_name(ensembl_pref), common_name2)
        self.assertEqual(Species.get_ensembl_db_prefix(species_name), ensembl_pref)
        self.assertEqual(Species.get_ensembl_db_prefix(common_name2), ensembl_pref)

    def test_get_compara_name(self):
        """should correctly form valid names for assignment onto objects"""
        self.assertEqual(Species.get_compara_name("pika"), "Pika")
        self.assertEqual(Species.get_compara_name("C.elegans"), "Celegans")
        self.assertEqual(Species.get_compara_name("Caenorhabditis elegans"), "Celegans")

    def test_lookup_raises(self):
        """setting level  to raise should create exceptions"""
        self.assertRaises(ValueError, Species.get_species_name, "failme", level="raise")
        self.assertRaises(ValueError, Species.get_common_name, "failme", level="raise")
        self.assertRaises(ValueError, Species.get_ensembl_db_prefix, "failme")

    def test_synonyms_work(self):
        """species with synonyms should allow correct lookups"""
        self.assertEqual(
            Species.get_species_name("Canis lupus familiaris"), "Canis familiaris"
        )
        self.assertEqual(Species.get_common_name("Canis lupus familiaris"), "Dog")

    def test_get_synonyms(self):
        """correctly returns synonyms given either a common or latin name"""
        expect = ["Canis familiaris", "Canis lupus familiaris"]
        got = Species.get_synonymns("Canis familiaris")
        self.assertEqual(set(got), set(expect))
        got = Species.get_synonymns("Canis lupus familiaris")
        self.assertEqual(set(got), set(expect))
        got = Species.get_synonymns("dog")
        self.assertEqual(set(got), set(expect))

    def test_to_table(self):
        """returns a table object"""
        table = Species.to_table()
        self.assertIsInstance(table, Table)
        self.assertTrue(table.shape[0] > 20)
        self.assertEqual(table.shape[1], 5)


if __name__ == "__main__":
    main()
