import bz2
import gzip
import os
import subprocess
import sys

import numpy
from pkg_resources import resource_filename

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


def get_resource_dir():
    """returns path to resource directory"""
    if "ENSEMBLDBRC" in os.environ:
        path = os.environ["ENSEMBLDBRC"]
    else:
        path = resource_filename("ensembldb3", "data")

    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise ValueError("ENSEMBLDBRC directory '%s' does not exist")

    return path


# the following is where essential files live, such as
# the species/common name map and sample download.cfg
ENSEMBLDBRC = get_resource_dir()


def open_(filename, mode="r"):
    func = {"gz": gzip.open, "bz2": bz2.open}.get(filename.split(".")[-1], open)
    return func(filename, mode=mode)


def lftp_installed():
    """returns True if lftp installed"""
    if sys.platform().lower() == "windows":
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
        sys.stderr.writelines("FAILED: %s\n%s" % (cmnd, msg))
        sys.exit(proc.returncode)
    if out is not None:
        r = out.decode("utf8")
    else:
        r = None
    return r


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
        new = [new, "'%s'" % new][self.with_quotes]
        return new


class CaseInsensitiveString(str):
    """A case insensitive string class. Comparisons are case insensitive."""

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


class LazyRecord(object):
    """a convenience class for conducting lazy evaluations of class
    properties"""

    NULL_VALUE = None

    def __init__(self):
        """blind constructor of caches"""
        self._cached = {}
        self._table_rows = {}

    def _get_cached_value(self, attr_name, get_attr_func):
        if attr_name not in self._cached:
            get_attr_func()

        return self._cached[attr_name]

    def _populate_cache_from_record(self, attr_column, table_name):
        """attr_column: the attribute name <-> table column key mappin
        table_name: the key in _table_rows"""
        table = self._table_rows[table_name]
        for attr, column, func in attr_column:
            if attr not in self._cached:
                self._cached[attr] = func(table[column])

    def _set_null_values(self, attrs, table_name=None):
        for attr in attrs:
            self._cached[attr] = self.NULL_VALUE

        if table_name:
            self._table_rows[table_name] = self.NULL_VALUE


class NoItemError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def convert_strand(val):
    """ensures a consistent internal representation of strand"""
    if isinstance(val, str):
        assert val in "-+", 'unknown strand "%s"' % val
        val = [-1, 1][val == "+"]
    elif val is not None:
        val = [-1, 1][val > 0]
    else:
        val = 1
    return val


def asserted_one(items):
    """asserts items has a single value and returns it"""
    one = False
    for item in items:
        if one:
            raise ValueError("More than one: [%s]" % list(item.items()))
        one = True
    if one:
        return item
    else:
        raise NoItemError("No items")


def what_columns(table):
    """shows what columns are in a table"""
    print([c.name for c in table.c])


def yield_selected(sqlalchemy_select, limit=100):
    """takes a SQLAlchemy select condition, yielding limit per db query
    purpose is to not get all records from the db in one go
    """
    offset = 0
    while 1:
        results = sqlalchemy_select.offset(offset).limit(limit).execute()
        count = 0
        for result in results:
            yield result
            count += 1
        offset += limit
        if count == 0 or count < limit:
            break


def flatten(data):
    """returns 1D list

    Removes nesting via numpy and 1D traversal
    """
    data = numpy.array(data).flatten().tolist()
    result = []
    for element in data:
        try:
            result.extend(element)
        except Exception:
            result.append(element)
    return result
