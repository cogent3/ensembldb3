from __future__ import annotations

import os
import re

import click

from ensembl_cli.ftp_download import download_data, listdir
from ensembl_cli.species import Species
from ensembl_cli.util import ENSEMBLDBRC, read_config


_cfg = os.path.join(ENSEMBLDBRC, "ensembldb_download.cfg")


_valid_seq = re.compile("([.]dna[.]|README|CHECKSUMS)")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _valid_seq.search(name) is not None


class valid_gff3_file:
    """whole genome gff3"""

    def __init__(self, release: int) -> None:
        self._valid = re.compile(f"([.]{release}[.]gff3[.]gz|README|CHECKSUMS)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


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

    patterns = dict(fasta=valid_seq_file, gff3=valid_gff3_file(config.release))

    for key in config.species_dbs:
        db_prefix = Species.get_ensembl_db_prefix(key)
        local_root = config.staging_path / db_prefix
        local_root.mkdir(parents=True, exist_ok=True)
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
            )

    click.secho(f"\nWROTE to output path={config.staging_path}\n", fg="green")
