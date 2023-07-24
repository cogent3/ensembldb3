from __future__ import annotations

import os

import click

from ensembl_cli.name import EnsemblDbName
from ensembl_cli.util import (
    ENSEMBLDBRC,
    exec_command,
    lftp_installed,
    makedirs,
    read_config,
)


def get_download_checkpoint_path(local_path, dbname):
    """returns path to db checkpoint file"""
    return os.path.join(local_path, dbname, "ENSEMBLDB_DOWNLOADED")


def is_downloaded(local_path, dbname):
    """returns True if checkpoint file exists for dbname"""
    chk = get_download_checkpoint_path(local_path, dbname)
    return os.path.exists(chk)


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
