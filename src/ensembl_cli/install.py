import shutil

import click

from ensembl_cli.util import read_config

def _install_gffdb(config, species):
    ...

def _install_seqs(config, species):
    ...

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
    for species in config.species_dbs:
        sp_dir = config.install_path / species
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    for species in config.species_dbs:

    # we copy the dowloaded sequence data to the new location

    return config
