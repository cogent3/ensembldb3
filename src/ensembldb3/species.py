import os
import re
from collections import defaultdict

from cogent3.util.table import Table
from pkg_resources import resource_filename

from .util import ENSEMBLDBRC, CaseInsensitiveString

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "Jason Merkin"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

_invalid_chars = re.compile("[^a-zA-Z _]")


def load_species(species_path):
    """returns [[latin_name, common_name],..] from species_path

    if species_path does not exist, defaults to default one"""
    if not os.path.exists(species_path):
        species_path = resource_filename("data", "species.tsv")

    with open(species_path, "r") as infile:
        data = []
        for line in infile:
            line = [e.strip() for e in line.split("\t")]
            num_fields = len(line)
            if 2 <= num_fields <= 3:
                data.append(line)
            else:
                print(num_fields, line)
                raise ValueError(
                    "species file should be "
                    "<latin name>\t<common name>\t<optional species synonym>"
                    " per line"
                )

    return data


_species_common_map = load_species(os.path.join(ENSEMBLDBRC, "species.tsv"))


class SpeciesNameMap:
    """mapping between common names and latin names"""

    def __init__(self, species_common=_species_common_map):
        """provides latin name:common name mappings"""
        self._species_common = {}
        self._common_species = {}
        self._species_ensembl = {}
        self._ensembl_species = {}
        self._synonyms = defaultdict(str)
        for names in species_common:
            names = list(map(CaseInsensitiveString, names))
            self.amend_species(*names)

    def __str__(self):
        return str(self.to_table())

    def __repr__(self):
        return repr(self.to_table())

    def _repr_html_(self):
        table = self.to_table()
        return table._repr_html_()

    def add_synonym(self, species, synonym):
        """add a synonym for a species name

        This provides an additional mapping to common names and ensembl
        db names"""
        species = CaseInsensitiveString(species)
        synonym = CaseInsensitiveString(synonym)
        self._synonyms[synonym] = species
        ensembl_name = synonym.lower().replace(" ", "_")
        self._ensembl_species[ensembl_name] = synonym

    def get_synonymns(self, name: str) -> set:
        """all species names matching name

        Parameters
        ----------
        name
            can be common name or a species name
        """
        name = self.get_species_name(name)
        result = {str(k) for k, v in self._synonyms.items() if v == name}
        return result | {name}

    def get_common_name(self, name, level="raise"):
        """returns the common name for the given name (which can be either a
        species name or the ensembl version)"""
        name = CaseInsensitiveString(name)
        if name in self._ensembl_species:
            name = self._ensembl_species[name]

        if name in self._synonyms:
            name = self._synonyms[name]

        if name in self._species_common:
            common_name = self._species_common[name]
        elif name in self._common_species:
            common_name = name
        else:
            common_name = None

        if common_name is None:
            msg = f"Unknown species name: {name}"
            if level == "raise":
                raise ValueError(msg)
            elif level == "warn":
                print(f"WARN: {msg}")

        return str(common_name)

    def get_species_name(self, name, level="ignore"):
        """returns the species name for the given common name"""
        name = CaseInsensitiveString(name)
        if name in self._species_common:
            return str(name)

        if name in self._synonyms:
            name = self._synonyms[name]
            return str(name)

        species_name = None
        level = level.lower().strip()
        for data in [self._common_species, self._ensembl_species]:
            if name in data:
                species_name = data[name]
        if species_name is None:
            msg = f"Unknown common name: {name}"
            if level == "raise":
                raise ValueError(msg)
            elif level == "warn":
                print(f"WARN: {msg}")
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
            species_name = self.get_species_name(name, level="raise")
        except ValueError:
            if name not in self._species_common:
                raise ValueError(f"Unknown name {name}")
            species_name = name

        return str(species_name.lower().replace(" ", "_"))

    def get_compara_name(self, name):
        """the compara instance attribute name for species matching ``name``"""
        name = self.get_common_name(name)
        name = name.replace(".", "")
        name = name.split()
        for i, word in enumerate(name):
            name[i] = word.title()
        name = "".join(name)
        for invalid in _invalid_chars.findall(name):
            name = name.replace(invalid, "")

        return name

    def _purge_species(self, species_name):
        """removes a species record"""
        species_name = CaseInsensitiveString(species_name)
        if species_name not in self._species_common:
            return
        common_name = self._species_common.pop(species_name)
        ensembl_name = self._species_ensembl.pop(species_name)
        self._ensembl_species.pop(ensembl_name)
        self._common_species.pop(common_name)

    def amend_species(self, species_name, common_name, synonym=""):
        """add a new species, and common name"""
        species_name = CaseInsensitiveString(species_name)
        common_name = CaseInsensitiveString(common_name)
        assert "_" not in species_name, "'_' in species_name, not a Latin name?"
        self._purge_species(species_name)  # remove if existing
        self._species_common[species_name] = common_name
        self._common_species[common_name] = species_name
        ensembl_name = species_name.lower().replace(" ", "_")
        self._species_ensembl[species_name] = ensembl_name
        self._ensembl_species[ensembl_name] = species_name
        if synonym:
            self.add_synonym(species_name, CaseInsensitiveString(synonym))

    def to_table(self):
        """returns cogent3 Table"""
        rows = []
        have_syns = defaultdict(list)
        for syn in self._synonyms:
            have_syns[self._synonyms[syn]].append(syn)
        syns = dict([(sp, ", ".join(have_syns[sp])) for sp in have_syns])
        for common in self._common_species:
            species = self._common_species[common]
            ensembl = self._species_ensembl[species]
            syn = syns.get(species, "")
            compara = self.get_compara_name(species)

            rows += [[species, common, syn, ensembl, compara]]
        return Table(
            [
                "Species name",
                "Common name",
                "Synonymn",
                "Ensembl Db Prefix",
                "Compara.Name",
            ],
            data=rows,
            space=2,
        ).sorted()


Species = SpeciesNameMap()
