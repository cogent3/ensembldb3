from warnings import filterwarnings

from .compara import Compara
from .genome import Genome
from .host import HostAccount
from .species import Species
from .util import NoItemError


filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")


__all__ = [
    "assembly",
    "compara",
    "database",
    "genome",
    "host",
    "name",
    "region",
    "related_region",
    "sequence",
    "species",
    "util",
    "HostAccount",
    "Species",
    "Genome",
    "Compara",
]

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
