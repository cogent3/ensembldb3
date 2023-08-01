from __future__ import annotations

import os
import re
import shutil

import click

from ensembl_cli.ftp_download import download_data, listdir
from ensembl_cli.species import Species
from ensembl_cli.util import Config, get_resource_path, read_config


_cfg = get_resource_path("ensembldb_download.cfg")


_valid_seq = re.compile("([.]dna[.](?!toplevel)|README|CHECKSUMS)")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _valid_seq.search(name) is not None


class valid_gff3_file:
    """whole genome gff3"""

    def __init__(self, release: str) -> None:
        self._valid = re.compile(f"([.]{release}[.]gff3[.]gz|README|CHECKSUMS)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def _remove_tmpdirs(path: os.PathLike):
    """delete any tmp dirs left over from unsuccessful runs"""
    tmpdirs = [p for p in path.glob("tmp*") if p.is_dir()]
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)


def download_species(configpath: os.PathLike, debug: bool, verbose: bool) -> Config:
    """download seq and gff data"""
    if configpath.name == _cfg:
        click.secho("WARN: using the built in demo cfg, will write to /tmp", fg="red")

    config = read_config(configpath)

    remote_template = f"{config.remote_path}/release-{config.release}/" + "{}/{}"

    if verbose:
        click.secho(f"DOWNLOADING\n  ensembl release={config.release}", fg="green")
        click.secho("\n".join(f"  {d}" for d in config.species_dbs), fg="green")
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
            _remove_tmpdirs(dest_path)
            download_data(
                host=config.host,
                local_dest=dest_path,
                remote_paths=listdir(config.host, path=path, pattern=patterns[subdir]),
                description=f"{db_prefix[:5]}.../{subdir}",
            )

    return config


def download_compara(configpath: os.PathLike, verbose: bool) -> Config:
    ...
