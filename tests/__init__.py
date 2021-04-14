import os

os.chdir(os.path.dirname(__file__))

__all__ = [
    "test_ensembl_admin",
    "test_assembly",
    "test_compara",
    "test_database",
    "test_feature_level",
    "test_genome",
    "test_host",
    "test_species",
    "test_variation",
]

__author__ = "Gavin Huttley, Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

ENSEMBL_RELEASE = 102
ENSEMBL_GENOMES_RELEASE = 45
