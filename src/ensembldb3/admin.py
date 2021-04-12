import configparser
import os
import pathlib
import shutil
from collections import defaultdict
from pprint import pprint

import click
from cogent3 import make_table
from cogent3.util import parallel

from . import HostAccount
from .download import (_cfg, download_dbs, is_downloaded, read_config,
                       reduce_dirnames)
from .host import DbConnection, get_db_name
from .util import ENSEMBLDBRC, FileSet, exec_command, open_

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

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


def get_import_command(
    mysqlcfg, account, dbname, local_path, table_name=None, verbose=False
):
    info = read_mysql_config(mysqlcfg, "mysqlimport", verbose=verbose)
    command = info.get("command", r"mysqlimport  --local --fields_escaped_by=\\")
    acct = r" -u %(user)s -p%(passwd)s "
    host = "" if info["host"] is None else r" -h %(host)s "
    port = "" if info["port"] is None else r" --port %(port)s "
    opts = dict(
        host=info["host"] or account.host,
        user=info["user"] or account.user,
        passwd=info["passwd"] or account.passwd,
        port=info["port"] or account.port,
    )

    if table_name is None:
        cmnd = command + port + host + acct + f" {dbname} {local_path}/{dbname}/*.txt"
    else:
        cmnd = (
            command
            + port
            + host
            + acct
            + f" {dbname} {local_path}/{dbname}/{table_name}.txt"
        )
    cmnd = cmnd % opts
    return cmnd


def tables_to_install(checkpoint):
    """tables whose installation was not completed"""
    if checkpoint.exists():
        installed_tables = {l.strip() for l in checkpoint.read_text().splitlines()}
    else:
        installed_tables = set()

    all_tables = FileSet(checkpoint.parent, suffixes="txt")
    return installed_tables ^ all_tables


def get_installed_checkpoint_path(local_path, dbname):
    """returns path to db checkpoint file"""
    return pathlib.Path(local_path) / dbname / "ENSEMBLDB_INSTALLED"


def is_installed(local_path, dbname):
    """returns True if checkpoint file exists for dbname"""
    chk = get_installed_checkpoint_path(local_path, dbname)
    if not chk.exists():
        return False

    data = chk.read_text().splitlines()

    return data[-1] == INSTALL_COMPLETED


def install_one_db(
    mysqlcfg,
    server,
    account,
    dbname,
    local_path,
    force_overwrite=False,
    verbose=False,
    debug=False,
):
    """installs a single ensembl database"""
    checkpoint = get_installed_checkpoint_path(local_path, dbname)
    # first create the database in mysql
    # find the .sql file, load all contents into memory
    # then execute using mysql cursor?
    # ensembl instructions suggest the following
    # $ mysql -u uname dname < dbname.sql
    dbpath = pathlib.Path(local_path) / dbname
    if is_installed(local_path, dbname) and not force_overwrite:
        click.echo(f"ALREADY INSTALLED: {dbname}, skipping")
        return True

    sqlfile = listpaths(dbpath, "*.sql*")
    if not sqlfile:
        raise RuntimeError(f"sql file not present in {dbpath}")

    sqlfile = sqlfile[0]
    with open_(sqlfile, mode="rt") as infile:
        sql = infile.readlines()
        sql = [l.rstrip("\n") for l in sql]

    sql = "\n".join(sql)
    # select the database
    if verbose or debug:
        click.echo(f"  Creating table definitions for {dbname}")

    cursor = server.cursor()

    table_names = tables_to_install(checkpoint)
    if debug:
        click.echo(str(table_names))

    r = cursor.execute(f"USE {dbname}")
    # make sure tables don't exist
    for table in table_names:
        click.echo(table)
        r = cursor.execute(f"DROP TABLE IF EXISTS {table}")

    # create the table definitions
    num_tables = sql.count("CREATE TABLE")
    r = cursor.execute(sql)

    if debug:
        _display_sql_created_diff_error(cursor, sql)
        display_dbs_tables(cursor, dbname)

    for table_name in table_names:
        # turn off key creation to speed up loading
        server.ping(reconnect=True)
        cursor.execute(f"ALTER TABLE `{table_name}` DISABLE KEYS")
        server.commit()

        tablepath = dbpath / f"{table_name}.txt.gz"
        if tablepath.exists():
            decompress_files(tablepath)

        tablepath = dbpath / f"{table_name}.txt"

        mysqlimport_command = get_import_command(
            mysqlcfg,
            account,
            dbname,
            local_path,
            table_name=table_name,
            verbose=verbose,
        )
        r = exec_command(mysqlimport_command)

        # turn on key creation to speed up usage
        server.ping(reconnect=True)
        cursor.execute(f"ALTER TABLE `{table_name}` ENABLE KEYS")
        server.commit()

        # we open to append in line buffering mode
        with checkpoint.open(mode="at", buffering=1) as out:
            out.write(f"{table_name}\n")

    with checkpoint.open(mode="at", buffering=1) as out:
        out.write(f"{INSTALL_COMPLETED}\n")

    cursor.close()


def _display_sql_created_diff_error(cursor, sql):
    _ = cursor.execute("SHOW TABLES")
    result = cursor.fetchall()
    expected = set()
    for line in sql.splitlines():
        if "CREATE TABLE" in line:
            table = line.replace("`", "").split()[2]
            expected.update([table])

    got = {str(g[0]) for g in result}

    diff = got ^ expected  # symmetric diff
    msg = [
        "ERR: The symmetric difference in tables between "
        "SQL statement and actually created:",
        f"\t{str(diff)}",
    ]
    click.secho("\n".join(msg), fg="red")

    raise RuntimeError("number of created tables doesn't match number in sql")


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


# can we get away with the ~/.my.cnf for the sql cursor?
def read_mysql_config(config_path, section, verbose=False):
    """returns a dict with mysql config options

    Parameters
    ----------
    config_path : str
      path to mysql config file
    section : str
      section in config file to query
    """
    config_path.seek(0)  # make sure at start of file
    opts = defaultdict(lambda: None)
    parser = configparser.ConfigParser(opts)
    parser.read_file(config_path)
    if not parser.has_section(section):
        return opts

    for k in ["host", "user", "passwd", "command", "port"]:
        if not parser.has_option(section, k):
            continue

        val = parser.get(section, k) or opts[k]
        if k == "port" and val is not None:
            val = int(val)

        opts[k] = val

    return opts


def _drop_db(cursor, dbname):
    """drops the database"""
    sql = f"DROP DATABASE IF EXISTS {dbname}"
    cursor.execute(sql)


def display_dbs(cursor, release):
    """shows what databases for the nominated release exist at the server"""
    sql = "SHOW DATABASES"
    r = cursor.execute(sql)
    result = cursor.fetchall()
    for r in result:
        if isinstance(r, tuple):
            r = r[0]

        if release in r:
            click.echo(str(r))


def display_dbs_tables(cursor, dbname):
    """shows what databases for the nominated release exist at the server"""
    sql = "SHOW TABLES"
    r = cursor.execute(sql)
    result = cursor.fetchall()
    for r in result:
        if isinstance(r, tuple):
            r = r[0]

        click.echo(str(r))


# default mysql config
_mycfg = os.path.join(ENSEMBLDBRC, "mysql.cfg")

# defining some of the options
_cfgpath = click.option(
    "-c",
    "--configpath",
    default=_cfg,
    type=click.File(),
    help="path to config file specifying databases, only "
    "species or compara at present",
)
_mysqlcfg = click.option(
    "-m",
    "--mysqlcfg",
    default=_mycfg,
    type=click.File(),
    help="path to mysql config file specifying host, "
    "user, installing data for writing",
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
    type=click.Path(),
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
@_numprocs
@_verbose
@_debug
def download(configpath, numprocs, verbose, debug):
    """download databases from Ensembl using rsync, can be done in parallel"""
    download_dbs(configpath, numprocs, verbose, debug)


@main.command()
@_cfgpath
@_mysqlcfg
@_force
@_verbose
@_debug
def install(configpath, mysqlcfg, force_overwrite, verbose, debug):
    """install ensembl databases into a MySQL server"""
    mysql_info = read_mysql_config(mysqlcfg, "mysql")
    account = HostAccount(
        mysql_info["host"],
        mysql_info["user"],
        mysql_info["passwd"],
        port=mysql_info["port"],
    )

    server = DbConnection(account, db_name="PARENT", pool_recycle=36000)

    release, remote_path, local_path, species_dbs = read_config(configpath)
    content = os.listdir(local_path)
    cursor = server.cursor()
    sql_for_speed = ["SET FOREIGN_KEY_CHECKS=0;", "SET UNIQUE_CHECKS=0;"]
    cursor.execute("\n".join(sql_for_speed))
    dbnames = reduce_dirnames(content, species_dbs)
    dbnames = sorted_by_size(local_path, dbnames, debug=debug)
    for dbname in dbnames:
        chk = get_installed_checkpoint_path(local_path, dbname.name)
        server.ping(reconnect=True)  # reconnect if server not alive
        if force_overwrite:
            _drop_db(cursor, dbname.name)
            if chk.exists():
                chk.unlink()

        if verbose:
            click.echo(f"Creating database {dbname.name}")

        # now create dbname
        sql = f"CREATE DATABASE IF NOT EXISTS {dbname}"
        r = cursor.execute(sql)
        install_one_db(
            mysqlcfg,
            server,
            account,
            dbname.name,
            local_path,
            force_overwrite=force_overwrite,
            verbose=verbose,
            debug=debug,
        )

    if debug:
        display_dbs(cursor, release)
        click.echo(server)

    undo_sql_for_speed = ["SET FOREIGN_KEY_CHECKS=1;", "SET UNIQUE_CHECKS=1;"]
    server.ping(reconnect=True)
    cursor.execute("\n".join(undo_sql_for_speed))
    cursor.close()


@main.command()
@_cfgpath
@_mysqlcfg
@_verbose
@_debug
def drop(configpath, mysqlcfg, verbose, debug):
    """drop databases from a MySQL server"""
    mysql_info = read_mysql_config(mysqlcfg, "mysql")
    account = HostAccount(
        mysql_info["host"],
        mysql_info["user"],
        mysql_info["passwd"],
        port=mysql_info["port"],
    )
    server = DbConnection(account, db_name="PARENT", pool_recycle=36000)
    cursor = server.cursor()
    release, remote_path, local_path, species_dbs = read_config(configpath)
    content = get_db_name(account=account, release=str(release))
    content = [str(n) for n in content]
    dbnames = reduce_dirnames(content, species_dbs)

    click.echo("The following databases will be deleted:")
    click.echo("\n".join("  %s" % d for d in dbnames))
    try:
        click.confirm("Confirm you want to delete the databases", abort=True)
    except click.exceptions.Abort:
        click.echo("EXITING")
        exit(0)

    for dbname in dbnames:
        click.echo(f"Dropping {dbname}")
        _drop_db(cursor, dbname)

    if verbose:
        display_dbs(cursor, release)

    cursor.close()


@main.command()
@_dbrc_out
def exportrc(outpath):
    """exports the rc directory to the nominated path

    setting an environment variable ENSEMBLDBRC with this path
    will force it's contents to override the default ensembldb3 settings"""
    shutil.copytree(ENSEMBLDBRC, outpath)
    click.echo(f"Contents written to {outpath}")


@main.command()
@_release
@_mysqlcfg
def show(release, mysqlcfg):
    """shows databases corresponding to release"""
    if mysqlcfg.name == _mycfg:
        click.secho(f"{show.help}\n")
        click.secho("use --help for more options")

        exit()

    mysql_info = read_mysql_config(mysqlcfg, "mysql")
    account = HostAccount(
        mysql_info["host"],
        mysql_info["user"],
        mysql_info["passwd"],
        port=mysql_info["port"],
    )
    names = get_db_name(account=account, release=str(release))
    click.echo(f"Databases at host='{account.host}' for release={release}")
    if names:
        click.echo("\n".join("  %s" % n for n in names))
    else:
        click.echo("  None")


@main.command()
@_cfgpath
def status(configpath):
    """checks download/install status using checkpoint files and config"""
    release, remote_path, local_path, species_dbs = read_config(configpath)
    content = os.listdir(local_path)
    dbnames = reduce_dirnames(content, species_dbs)
    rows = []
    for db in dbnames:
        row = [
            db.name,
            is_downloaded(local_path, db.name),
            is_installed(local_path, db.name),
        ]
        rows.append(row)

    table = make_table(
        header=["dbname", "Downloaded", "Installed"],
        rows=rows,
        title="Status of download and install",
        legend=f"config={configpath.name}; local_path={local_path}",
    )
    print(table)


if __name__ == "__main__":
    main()
