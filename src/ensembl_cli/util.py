import configparser
import os
import pathlib
import shutil
import subprocess
import sys
import uuid

from dataclasses import dataclass
from tempfile import mkdtemp
from typing import IO, Iterable, Union

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


@dataclass
class Config:
    host: str
    remote_path: str
    release: str
    staging_path: os.PathLike
    install_path: os.PathLike
    species_dbs: Iterable[str]
    align_names: Iterable[str]

    @property
    def db_names(self) -> Iterable[str]:
        from ensembl_cli.species import Species

        for species in self.species_dbs:
            yield Species.get_ensembl_db_prefix(species)


def read_config(config_path) -> Config:
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    from ensembl_cli.species import Species

    parser = configparser.ConfigParser()

    with config_path.open() as f:
        parser.read_file(f)

    release = parser.get("release", "release")
    host = parser.get("remote path", "host")
    remote_path = parser.get("remote path", "path")
    remote_path = remote_path[:-1] if remote_path.endswith("/") else remote_path
    staging_path = (
        pathlib.Path(parser.get("local path", "staging_path")).expanduser().absolute()
    )
    install_path = (
        pathlib.Path(parser.get("local path", "install_path")).expanduser().absolute()
    )

    align_names = None
    species_dbs = {}
    get_option = parser.get
    for section in parser.sections():
        if section in ("release", "remote path", "local path"):
            continue

        if section == "compara":
            align_names = [
                n.strip() for n in get_option(section, "align_names").split(",")
            ]
            continue

        dbs = [db.strip() for db in get_option(section, "db").split(",")]

        # handle synonyms
        species = Species.get_species_name(section, level="raise")
        for synonym in Species.get_synonymns(species):
            species_dbs[synonym] = dbs

    return Config(
        host=host,
        remote_path=remote_path,
        release=release,
        staging_path=staging_path,
        install_path=install_path,
        species_dbs=species_dbs,
        align_names=align_names,
    )


def load_ensembl_checksum(path: os.PathLike) -> dict:
    """loads the BSD checksums from Ensembl CHECKSUMS file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, b, p = line.split()
        result[p] = int(s), int(b)
    result.pop("README", None)
    return result


def load_ensembl_md5sum(path: os.PathLike) -> dict:
    """loads the md5 sum from Ensembl MD5SUM file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, p = line.split()
        result[p] = s
    result.pop("README", None)
    return result


class atomic_write:
    """performs atomic write operations, cleans up if fails"""

    def __init__(self, path: os.PathLike, tmpdir=None, mode="wb", encoding=None):
        """

        Parameters
        ----------
        path
            path to file
        tmpdir
            directory where temporary file will be created
        mode
            file writing mode
        encoding
            text encoding
        """
        path = pathlib.Path(path).expanduser()

        self._path = path
        self._mode = mode
        self._file = None
        self._encoding = encoding
        self._tmppath = self._make_tmppath(tmpdir)

        self.succeeded = None
        self._close_func = self._close_rename_standard

    def _make_tmppath(self, tmpdir):
        """returns path of temporary file

        Parameters
        ----------
        tmpdir: Path
            to directory

        Returns
        -------
        full path to a temporary file

        Notes
        -----
        Uses a random uuid as the file name, adds suffixes from path
        """
        suffixes = "".join(self._path.suffixes)
        parent = self._path.parent
        name = f"{uuid.uuid4()}{suffixes}"
        tmpdir = (
            pathlib.Path(mkdtemp(dir=parent))
            if tmpdir is None
            else pathlib.Path(tmpdir)
        )

        if not tmpdir.exists():
            raise FileNotFoundError(f"{tmpdir} directory does not exist")

        return tmpdir / name

    def _get_fileobj(self):
        """returns file to be written to"""
        if self._file is None:
            self._file = open(self._tmppath, self._mode)

        return self._file

    def __enter__(self) -> IO:
        return self._get_fileobj()

    def _close_rename_standard(self, src):
        dest = pathlib.Path(self._path)
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
        finally:
            src.rename(dest)

        shutil.rmtree(src.parent)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()
        if exc_type is None:
            self._close_func(self._tmppath)
            self.succeeded = True
        else:
            self.succeeded = False

        shutil.rmtree(self._tmppath.parent, ignore_errors=True)

    def write(self, text):
        """writes text to file"""
        fileobj = self._get_fileobj()
        fileobj.write(text)

    def close(self):
        """closes file"""
        self.__exit__(None, None, None)
