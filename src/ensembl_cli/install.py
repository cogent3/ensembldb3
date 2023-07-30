import os
import shutil

import click

from cogent3 import load_seq, open_
from rich.progress import track
from unsync import unsync

from ensembl_cli.util import read_config


@unsync(cpu_bound=True)
def install_one(src, dest_dir):
    seq = load_seq(src, moltype="dna", label_to_name=lambda x: x.split()[0])
    with open_(dest_dir / f"{seq.name}.fa.gz", mode="wt") as outfile:
        outfile.write(seq.to_fasta())
    return True


def _install_gffdb(config, species):
    ...


def _install_seqs(src_dir: os.PathLike, dest_dir: os.PathLike):
    src_dir = src_dir / "fasta"
    paths = list(src_dir.glob("*.fa.gz"))
    return [install_one(path, dest_dir) for path in paths]


def local_install(configpath, force_overwrite):
    config = read_config(configpath)
    if config.install_path.exists() and not force_overwrite:
        click.secho(f"EXITING: {config.install_path} already exists", fg="red")
        exit(1)
    if force_overwrite:
        shutil.rmtree(config.install_path, ignore_errors=True)
    # we create the local installation
    config.install_path.mkdir(parents=True, exist_ok=True)
    # we create subdirectories for each species
    for db_name in config.db_names:
        sp_dir = config.install_path / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    tasks = []
    for db_name in config.db_names:
        src_dir = config.staging_path / db_name
        dest_dir = config.install_path / db_name
        tasks.extend(_install_seqs(src_dir, dest_dir))

    # we do all tasks in one go
    _ = [
        t.result() for t in track(tasks, description="Installing seqs", transient=True)
    ]

    # we copy the dowloaded sequence data to the new location

    return config
