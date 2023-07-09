from warnings import filterwarnings

from .host import HostAccount
from .species import Species
from .util import NoItemError


filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")


__all__ = [
    "host",
    "name",
    "species",
    "util",
    "HostAccount",
    "Species",
]

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
