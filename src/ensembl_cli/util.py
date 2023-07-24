import configparser
import os
import pathlib
import re
import subprocess
import sys

from math import ceil
from typing import Union

import numba
import numpy

# based on https://www.reddit.com/r/learnpython/comments/9bpgjl/implementing_bsd_16bit_checksum/
# and https://www.gnu.org/software/coreutils/manual/html_node/sum-invocation.html#sum-invocation
@numba.jit(nopython=True)
def checksum(data: bytes, size: int):
    """computes BSD style checksum"""
    # equivalent to command line BSD sum
    nb = numpy.ceil(size / 1024)
    cksum = 0
    for c in data:
        cksum = (cksum >> 1) + ((cksum & 1) << 15)
        cksum += c
        cksum &= 0xFFFF
    return cksum, nb


def _get_resource_dir() -> os.PathLike:
    """returns path to resource directory"""
    if "ENSEMBLDBRC" in os.environ:
        path = os.environ["ENSEMBLDBRC"]
    else:
        from ensembl_cli import data

        path = pathlib.Path(data.__file__).parent

    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise ValueError("ENSEMBLDBRC directory '%s' does not exist")

    return pathlib.Path(path)


def get_resource_path(resource: Union[str, os.PathLike]) -> os.PathLike:
    path = ENSEMBLDBRC / resource
    assert path.exists()
    return path


# the following is where essential files live, such as
# the species/common name map and sample download.cfg
ENSEMBLDBRC = _get_resource_dir()


def lftp_installed():
    """returns True if lftp installed"""
    if sys.platform.lower() == "windows":
        raise RuntimeError("not supported on windows")

    r = subprocess.call(["which", "lftp"])
    return r == 0


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        sys.stderr.writelines(f"FAILED: {cmnd}\n{msg}")
        sys.exit(proc.returncode)
    return out.decode("utf8") if out is not None else None


def makedirs(path):
    """creates directory path if it doesn't exist"""
    if os.path.exists(path):
        return

    os.makedirs(path)


def abspath(path):
    path = os.path.abspath(os.path.expanduser(path))
    return path


class DisplayString(str):
    """provides a mechanism for customising the str() and repr() of objects"""

    def __new__(cls, arg, num_words=None, repr_length=None, with_quotes=False):
        new = str.__new__(cls, str(arg))
        new.num_words = num_words
        new.repr_length = repr_length or len(str(arg))
        new.with_quotes = with_quotes
        return new

    def __repr__(self):
        if self.num_words is not None:
            new = " ".join(self.split()[: self.num_words])
        elif self.repr_length != len(self):
            new = self[: self.repr_length]
        else:
            new = self
        if len(self) > len(new):
            new += "..."
        new = [new, f"'{new}'"][self.with_quotes]
        return new


class CaseInsensitiveString(str):
    """A case-insensitive string class. Comparisons are also case-insensitive."""

    def __new__(cls, arg, h=None):
        n = str.__new__(cls, str(arg))
        n._lower = "".join(list(n)).lower()
        n._hash = hash(n._lower)
        return n

    def __eq__(self, other):
        return self._lower == "".join(list(other)).lower()

    def __hash__(self):
        # dict hashing done via lower case
        return self._hash

    def __str__(self):
        return "".join(list(self))


def convert_strand(val):
    """ensures a consistent internal representation of strand"""
    if isinstance(val, str):
        assert val in "-+", f'unknown strand "{val}"'
        val = [-1, 1][val == "+"]
    elif val is not None:
        val = [-1, 1][val > 0]
    else:
        val = 1
    return val


class FileSet(set):
    """Determines names of all files in a directory with matching suffixes.
    Does not recurse."""

    _dot = re.compile(r"^\.")

    def __init__(self, path, suffixes=("txt", "sql"), trim_suffixes=True):
        """
        Parameters
        ----------
        path
            directory path
        suffixes
            filename suffixes to match
        trim_suffixes
            whether suffixes are to be trimmed before adding to self
        """
        super().__init__()
        if isinstance(suffixes, str):
            suffixes = (suffixes,)

        path = pathlib.Path(path).expanduser()
        suffixes = {self._dot.sub("", s) for s in suffixes}
        collected = set()
        for p in path.glob("*"):
            if p.is_dir() or p.name.startswith("."):
                # don't consider nested directories, hidden files, files without a suffix
                continue

            if not {self._dot.sub("", s) for s in p.suffixes} & suffixes:
                continue

            if trim_suffixes:
                num_suffixes = len(p.suffixes)
                name = ".".join(p.name.split(".")[:-num_suffixes])
            else:
                name = p.name

            collected.add(name)

        self.path = path
        self.update(collected)


def read_config(config_path, verbose=False):
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    from ensembl_cli.species import Species

    parser = configparser.ConfigParser()
    parser.read_file(config_path)
    release = parser.get("release", "release")
    host = parser.get("remote path", "host")
    remote_path = parser.get("remote path", "path")
    remote_path = remote_path[:-1] if remote_path.endswith("/") else remote_path
    local_path = pathlib.Path(parser.get("local path", "path")).expanduser().absolute()
    species_dbs = {}
    for section in parser.sections():
        if section in ("release", "remote path", "local path"):
            continue

        dbs = [db.strip() for db in parser.get(section, "db").split(",")]

        if section == "compara":
            species_dbs["compara"] = dbs
            continue

        # handle synonymns
        species = Species.get_species_name(section, level="raise")
        for synonym in Species.get_synonymns(species):
            species_dbs[synonym] = dbs

    return host, remote_path, release, local_path, species_dbs


def load_ensembl_checksum(path: os.PathLike) -> dict:
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, b, p = line.split()
        result[p] = int(s), int(b)
    result.pop("README", None)
    return result
