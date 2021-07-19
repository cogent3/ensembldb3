import configparser
import os
import sys
import warnings

from pprint import pprint

import click

from ensembldb3.name import EnsemblDbName
from ensembldb3.species import Species
from ensembldb3.util import (
    ENSEMBLDBRC,
    abspath,
    exec_command,
    lftp_installed,
    makedirs,
)


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


def get_download_checkpoint_path(local_path, dbname):
    """returns path to db checkpoint file"""
    return os.path.join(local_path, dbname, "ENSEMBLDB_DOWNLOADED")


def is_downloaded(local_path, dbname):
    """returns True if checkpoint file exists for dbname"""
    chk = get_download_checkpoint_path(local_path, dbname)
    return os.path.exists(chk)


def do_lftp_command(host, remote_dir, db_name, local_dir, numprocs):
    command = (
        f'lftp -e "cd {remote_dir}; mirror -c '
        f'--use-pget-n={numprocs} --parallel={numprocs} {db_name} {local_dir}; bye" {host}'
    )
    result = exec_command(command)
    result = result.split("\n") if result != "" else []
    return result


def lftp_listdir(host, dirname="", debug=True):
    """returns directory listing"""
    cmnd = f'lftp -e "cd {dirname}; nlist; bye;" ftp://{host}'
    if debug:
        print(cmnd)
    result = exec_command(cmnd)
    return result.splitlines()


def rsync_listdir(remote_path, dirname="", debug=True):
    cmnd = f"{remote_path}{dirname}" if dirname else remote_path
    cmnd = r"rsync --list-only rsync://%s" % cmnd
    if debug:
        print(cmnd)
    result = exec_command(cmnd)
    return result.splitlines()


def _sort_dbs(dbnames):
    """returns the dbnames sorted by their type"""
    order = {"compara": 2, "variation": 3, "otherfeatures": 1}
    names = [(order.get(n.type, 0), n.name, n) for n in dbnames]
    return [db for i, n, db in sorted(names)]


def reduce_dirnames(dirnames, species_dbs, verbose=False, debug=False):
    """returns EnsemblNames corresponding to species db's and sort by type

    sort order put's core db's first, compara and variation last"""
    if debug:
        pprint(dirnames)

    db_names = []
    for record in dirnames:
        record = record.strip()
        if not record or record.endswith(".gz"):
            continue

        record = record.split()[-1]
        if not record[0].isalpha():
            continue

        try:
            name = EnsemblDbName(record)
        except (TypeError, RuntimeError):
            # a non-species
            if debug:
                print(record)
            continue

        if name.species in species_dbs:
            if name.type not in species_dbs[name.species] and species_dbs[name.species]:
                if debug or verbose:
                    print("Skipping", name)
                continue

            db_names.append(name)
        elif name.type == "compara" and "compara" in species_dbs:
            db_names.append(name)

    db_names = _sort_dbs(db_names)
    return db_names


def read_config(config_path, verbose=False):
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    parser = configparser.ConfigParser()
    parser.read_file(config_path)
    release = parser.get("release", "release")
    remote_path = parser.get("remote path", "path")
    local_path = parser.get("local path", "path")
    local_path = abspath(local_path)
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

    return release, remote_path, local_path, species_dbs


_cfg = os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg")


class Download:
    """callable instance that takes the database name and lftp downloads"""

    def __init__(
        self, host, local_base, release, numprocs, verbose, debug, dry_run=False
    ):
        self._host = host
        self._local_base = local_base
        self._release = release
        self._numprocs = numprocs
        self._verbose = verbose
        self._debug = debug
        self._dry_run = dry_run

    def __call__(self, dbname):
        if is_downloaded(self._local_base, dbname):
            if self._verbose or self._debug:
                click.secho(f"Already downloaded: {dbname}, skipping", fg="green")
            return

        commands = do_lftp_command(
            self._host,
            f"release-{self._release}/mysql/",
            dbname,
            os.path.join(self._local_base, dbname),
            self._numprocs,
        )
        checkpoint_file = get_download_checkpoint_path(self._local_base, dbname)
        with open(checkpoint_file, "w") as checked:
            pass
        click.secho(f"Completed download: {dbname}", fg="green")


def download_dbs(configpath, numprocs, verbose, debug):
    if not lftp_installed():
        click.secho("download requires lftp", fg="red")
        sys.exit(1)

    if configpath.name == _cfg:
        warnings.warn("WARN: using the built in demo cfg, will write to /tmp")

    release, remote_path, local_path, sp_db = read_config(configpath, verbose=verbose)
    makedirs(local_path)
    contents = lftp_listdir(
        remote_path, dirname=f"release-{release}/mysql/", debug=debug
    )
    db_names = reduce_dirnames(contents, sp_db, verbose=verbose, debug=debug)
    if verbose:
        click.secho(f"DOWNLOADING\n  ensembl release={release}", fg="green")
        click.secho("\n".join(f"  {d.name}" for d in db_names), fg="green")
        click.secho(f"\nWRITING to output path={local_path}\n", fg="green")

    lftp = Download(
        remote_path, local_path, release, numprocs, verbose=verbose, debug=debug
    )
    for db in db_names:
        lftp(db.name)

    click.secho(f"\nWROTE to output path={local_path}\n", fg="green")
