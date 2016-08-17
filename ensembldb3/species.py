import os
from pkg_resources import resource_filename
from cogent3.util.table import Table

from .util import CaseInsensitiveString, ENSEMBLDBRC

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Jason Merkin"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

def load_species(species_path):
    """returns [[latin_name, common_name],..] from species_path
    
    if species_path does not exist, defaults to default one"""
    if not os.path.exists(species_path):
        species_path = resource_filename("data", "species.tsv")
    
    with open(species_path, "r") as infile:
        data = []
        for line in infile:
            line = [e.strip() for e in line.split('\t')]
            if len(line) != 2:
                raise ValueError("species file should be <latin name>\t<common name> per line")
            data.append(line)
        
    return data

_species_common_map = load_species(os.path.join(ENSEMBLDBRC, "species.tsv"))


class SpeciesNameMap(dict):
    """mapping between common names and latin names"""

    def __init__(self, species_common=_species_common_map):
        """provides latin name:common name mappings"""
        self._species_common = {}
        self._common_species = {}
        self._species_ensembl = {}
        self._ensembl_species = {}
        for species_name, common_name in species_common:
            self.amend_species(CaseInsensitiveString(species_name),
                              CaseInsensitiveString(common_name))

    def __str__(self):
        rows = []
        for common in self._common_species:
            species = self._common_species[common]
            ensembl = self._species_ensembl[species]
            rows += [[common, species, ensembl]]
        return str(Table(['Common name', 'Species name', 'Ensembl Db Prefix'],
                         rows=rows, space=2).sorted())

    def __repr__(self):
        return 'Available species: %s' % ("'" +
                                          "'; '".join(list(self._common_species.keys())) + "'")

    def get_common_name(self, name, level='raise'):
        """returns the common name for the given name (which can be either a
        species name or the ensembl version)"""
        name = CaseInsensitiveString(name)
        if name in self._ensembl_species:
            name = self._ensembl_species[name]

        if name in self._species_common:
            common_name = self._species_common[name]
        elif name in self._common_species:
            common_name = name
        else:
            common_name = None

        if common_name is None:
            msg = "Unknown species name: %s" % name
            if level == 'raise':
                raise RuntimeError(msg)
            elif level == 'warn':
                print("WARN: %s" % msg)

        return str(common_name)

    def get_species_name(self, name, level='ignore'):
        """returns the species name for the given common name"""
        name = CaseInsensitiveString(name)
        if name in self._species_common:
            return str(name)
        species_name = None
        level = level.lower().strip()
        name = name
        for data in [self._common_species, self._ensembl_species]:
            if name in data:
                species_name = data[name]
        if species_name is None:
            msg = "Unknown common name: %s" % name
            if level == 'raise':
                raise RuntimeError(msg)
            elif level == 'warn':
                print("WARN: %s" % msg)
        return str(species_name)

    def get_species_names(self):
        """returns the list of species names"""
        names = list(self._species_common.keys())
        names.sort()
        return [str(n) for n in names]

    def get_ensembl_db_prefix(self, name):
        """returns a string of the species name in the format used by
        ensembl"""
        name = CaseInsensitiveString(name)
        if name in self._common_species:
            name = self._common_species[name]
        try:
            species_name = self.get_species_name(name, level='raise')
        except RuntimeError:
            if name not in self._species_common:
                raise RuntimeError("Unknown name %s" % name)
            species_name = name

        return str(species_name.lower().replace(" ", "_"))

    def get_compara_name(self, name):
        """returns string matching a compara instance attribute name for a
        species"""
        name = self.get_common_name(name)
        if '.' in name:
            name = name.replace('.', '')
        else:
            name = name.title()

        name = name.split()
        return ''.join(name)

    def _purge_species(self, species_name):
        """removes a species record"""
        species_name = CaseInsensitiveString(species_name)
        if species_name not in self._species_common:
            return
        common_name = self._species_common.pop(species_name)
        ensembl_name = self._species_ensembl.pop(species_name)
        self._ensembl_species.pop(ensembl_name)
        self._common_species.pop(common_name)

    def amend_species(self, species_name, common_name):
        """add a new species, and common name"""
        species_name = CaseInsensitiveString(species_name)
        common_name = CaseInsensitiveString(common_name)
        assert "_" not in species_name,\
            "'_' in species_name, not a Latin name?"
        self._purge_species(species_name)  # remove if existing
        self._species_common[species_name] = common_name
        self._common_species[common_name] = species_name
        ensembl_name = species_name.lower().replace(" ", "_")
        self._species_ensembl[species_name] = ensembl_name
        self._ensembl_species[ensembl_name] = species_name
        return

Species = SpeciesNameMap()
