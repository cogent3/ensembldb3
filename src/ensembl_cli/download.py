from __future__ import annotations

import os
import re

import click

from cogent3 import open_

from ensembl_cli.ftp_download import download_data, listdir
from ensembl_cli.name import EnsemblDbName
from ensembl_cli.species import Species
from ensembl_cli.util import ENSEMBLDBRC, read_config


def get_download_checkpoint_path(local_path, dbname):
    """returns path to db checkpoint file"""
    return os.path.join(local_path, dbname, "ENSEMBLDB_DOWNLOADED")


def is_downloaded(local_path, dbname):
    """returns True if checkpoint file exists for dbname"""
    chk = get_download_checkpoint_path(local_path, dbname)
    return os.path.exists(chk)


_cfg = os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg")


_valid_seq = re.compile("([.]dna[.]|README|CHECKSUMS)")
_valid_gff = re.compile("([.]\d+[.]gff3[.]gz|README|CHECKSUMS)")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _valid_seq.search(name) is not None


def valid_gff3_file(name: str) -> bool:
    """whole genome gff3"""

    return _valid_gff.search(name) is not None


def download_dbs(configpath, verbose):
    if configpath.name == _cfg:
        click.secho("WARN: using the built in demo cfg, will write to /tmp", fg="red")

    config = read_config(configpath, verbose=verbose)

    # TODO identify single file name convention enabling single file downloads from subdir
    remote_template = f"{config.remote_path}/release-{config.release}/" + "{}/{}"

    if verbose:
        click.secho(f"DOWNLOADING\n  ensembl release={config.release}", fg="green")
        click.secho("\n".join(f"  {d.name}" for d in config.species_dbs), fg="green")
        click.secho(f"\nWRITING to output path={config.local_path}\n", fg="green")

    patterns = dict(fasta=valid_seq_file, gff3=valid_gff3_file)

    for key in config.species_dbs:
        db_prefix = Species.get_ensembl_db_prefix(key)
        local_root = config.staging_path / db_prefix
        local_root.mkdir(parents=True, exist_ok=True)
        with open_(local_root / "DOWNLOADED_CHECKSUMS", mode="w") as chkpt:
            for subdir in ("fasta", "gff3"):
                path = remote_template.format(subdir, db_prefix)
                path = f"{path}/dna" if subdir == "fasta" else path
                dest_path = config.staging_path / db_prefix / subdir
                dest_path.mkdir(parents=True, exist_ok=True)
                download_data(
                    config.host,
                    dest_path,
                    listdir(config.host, path=path, pattern=patterns[subdir]),
                    description=f"{db_prefix[:5]}.../{subdir}",
                    checkpoint_file=chkpt,
                )

    click.secho(f"\nWROTE to output path={config.staging_path}\n", fg="green")
