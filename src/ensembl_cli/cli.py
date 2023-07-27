import os
import pathlib
import shutil

from pprint import pprint

import click

from ensembl_cli import __version__
from ensembl_cli.download import _cfg, download_dbs
from ensembl_cli.util import ENSEMBLDBRC, exec_command, read_config
from src.ensembl_cli.install import local_install


INSTALL_COMPLETED = "INSTALL COMPLETED"


def listpaths(dirname, glob_pattern):
    """return path to all files matching glob_pattern"""
    fns = [str(p) for p in pathlib.Path(dirname).glob(glob_pattern)]
    if not fns:
        return None
    return fns


def decompress_files(local_path):
    """gunzip files


    Parameters
    ----------
    local_path: pathlib.Path
        single file, or directory

    Notes
    -----
    If directory, does all .gz files.
    """
    local_path = pathlib.Path(local_path)
    paths = [local_path] if local_path.is_file() else local_path.glob("*.gz")
    for path in paths:
        _ = exec_command(f"gunzip {path}")


def sorted_by_size(local_path, dbnames, debug=False):
    """returns dbnames ordered by directory size"""
    join = os.path.join
    getsize = os.path.getsize
    size_dbnames = []
    for dbname in dbnames:
        path = join(local_path, dbname.name)
        size = sum(getsize(join(path, fname)) for fname in os.listdir(path))
        size_dbnames.append([size, dbname])
    size_dbnames.sort()

    if debug:
        pprint(size_dbnames)

    sizes, dbnames = zip(*size_dbnames)
    return dbnames


# defining some of the options
_cfgpath = click.option(
    "-c",
    "--configpath",
    default=_cfg,
    type=click.File(),
    help="path to config file specifying databases, only "
    "species or compara at present",
)
_verbose = click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="causes stdout/stderr from rsync download to be " "written to screen",
)
_numprocs = click.option(
    "-n",
    "--numprocs",
    type=int,
    default=1,
    help="number of processes to use for download",
)
_force = click.option(
    "-f",
    "--force_overwrite",
    is_flag=True,
    help="drop existing database if it exists prior to " "installing",
)
_debug = click.option("-d", "--debug", is_flag=True, help="maximum verbosity")
_dbrc_out = click.option(
    "-o",
    "--outpath",
    type=pathlib.Path,
    help="path to directory to export all rc contents",
)
_release = click.option("-r", "--release", type=int, help="Ensembl release number")


@click.group()
@click.version_option(__version__)
def main():
    """admin tools for an Ensembl MySQL installation"""
    pass


@main.command()
@_cfgpath
@_verbose
def download(configpath, verbose):
    """download databases from Ensembl using rsync, can be done in parallel"""
    download_dbs(configpath, verbose)


@main.command()
@_cfgpath
@_force
@_verbose
def install(configpath, force_overwrite, verbose):
    """create the local db's"""
    config = local_install(configpath, force_overwrite)

    click.echo(f"Contents installed to {str(config.install_path)!r}")


@_dbrc_out
def exportrc(outpath):
    """exports the rc directory to the nominated path

    setting an environment variable ENSEMBLDBRC with this path
    will force it's contents to override the default ensembl_cli settings"""
    shutil.copytree(ENSEMBLDBRC, outpath)
    # remove the python module file
    for fn in pathlib.Path(outpath).glob("__init__.py*"):
        fn.unlink()
    click.echo(f"Contents written to {outpath}")


if __name__ == "__main__":
    main()
