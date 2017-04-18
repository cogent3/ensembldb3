*****************************
Database Administration tools
*****************************

We strongly encourage you to install a local copy of the Ensembl databases because running queries against a local installation is much, much faster! It will also alleviate the burden on Ensembl's UK servers.

If you get MySQL installed locally with sufficient storage to host the collection of species you're interested in, the administration tools we provide will make it relatively straightforward to administer the Ensembl databases and keep up-to-date with the Ensembl release cycle.

Installing MySQL
================

The bottom line is you need to install and configure MySQL yourself. `Ensembl offers some instructions <http://asia.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html>`_.

``ensembldb3`` command line tool
================================

Install of ``ensembldb3`` places a new executable ``ensembldb3`` on your path. This tool provides a number of capabilities as illustrated on the command line::
    
   $ ensembldb3
  Usage: ensembldb3 [OPTIONS] COMMAND [ARGS]...

    admin tools for an Ensembl MySQL installation

  Options:
    --help  Show this message and exit.

  Commands:
    download  download databases from Ensembl using rsync,...
    drop      drop databases from a MySQL server
    exportrc  exports the rc directory to the nominated...
    install   install ensembl databases into a MySQL server
    show      shows databases corresponding to release
    status    checks download/install status using...
    
.. _exportrc:

``ensembldb3 exportrc``
=======================

The command::

    $ ensembldb3 exportrc -o /path/to/ensembldbrc
    
produces a directory ``ensembldbrc`` containing 3 files that can be used by ``ensembldb3``:

species.tsv
    A tab delimited file with a species latin name and common name per line. This is used to define the common names that ``ensembldb3.Species`` uses for succinctly identifying species and their databases.

ensembldb_download.cfg
    A config file with sections for remote path, local path, release and the species of interest. In the latter case, their common names are used as the section title. The databases are specified by a comma separated line as core, variation, otherfeatures. The compara database has the same section title as the db name. Here's an example ::

        [remote path] # required
        path=ftp.ensembl.org/ensembl/pub/
        [local path] # required
        path=/tmp/ensembldb_download
        [release] # required
        release=85
        [S.cerevisiae]
        db=core
        [Xenopus]
        db=core
        [Human]
        db=core,variation
        [compara]
        db=compara

mysql.cfg
    A config file with sections for ``mysql`` and ``mysqlimport``. The sections include the command (full path) to the executable, including any command arguments and the account settings (username, password).
  
If you wish to use the contents of this directory you can create an environment variable ``ENSEMBLDBRC=/path/to/ensembldbrc``.

.. note::
    If ``ENSEMBLDBRC`` is defined in your environment, the ``species.tsv`` file within that directory will be used *for all ensembldb3 applications*.

``ensembldb3 download``
=======================

The command::
    
    $ ensembldb3 download -c /path/to/edited/ensembldb_download.cfg -n 3

will download databases for the species specified in the ``ensembldb_download.cfg`` config. The specific databases for the Ensembl release, remote and local paths must all be defined in that file (see :ref:`exportrc`). The ``-n`` option indicates the number of parallel processors to use for the download (maximum allowed is 5).

For the very large databases (e.g. compara or human variation) the download times can be very long. In which case we recommend, if running on a server, using the ``nohup`` command.

.. note:
    
    To use ``ensembldb3`` you only need to install the databases for the species you are interested in plus compara, if you wish to undertake comparative analyses.
    
.. note:
    
    An empty file called ``ENSEMBLDB_DOWNLOADED`` is written in each directory. This is used as a checkpoint marker to prevent needlessly downloading again.

``ensembldb3 install``
======================

The command::
    
    $ ensembldb3 install -c /path/to/edited/ensembldb_download.cfg -m /path/to/mysql.cfg -n 3
    
installs databases specified in the ``ensembldb_download.cfg`` config, into the mysql server specified by ``mysql.cfg``. In this instance, 3 processors are used for separate gzipping and ``mysqlimport`` of tables in parallel.

For the very large databases (e.g. compara or human variation) the install times can be very long. In which case we recommend, if running on a server, using the ``nohup`` command.

.. note:
    
    An empty file called ``ENSEMBLDB_INSTALLED`` is written in each directory. This is used as a checkpoint marker to prevent installing again unless overridden by the ``-f`` (force overwrite) flag.

``ensembldb3 drop``
===================

The command::

    $ ensembldb3 drop -c /path/to/edited/ensembldb_download.cfg -m /path/to/mysql.cfg

will drop the databases specified in the ``ensembldb_download.cfg`` from the mysql server specified by ``mysql.cfg``. You are required to confirm dropping listed databases.

``ensembldb3 show``
===================

The command::

    $ ensembldb3 show --release 85 -m /path/to/mysql.cfg

will display all databases from release 85 on the mysql host in the server specified by ``mysql.cfg``.

``ensembldb3 status``
=====================

The command::

    $ ensembldb3 status -c /path/to/edited/ensembldb_download.cfg

will display the download/install status of the databases specified by ``ensembldb_download.cfg``.

Trouble shooting
================

Many of the administrative functions wrap shell commands. If you encounter any issues, use the verbose flag (``-v``), causing shell commands to be printed to stdout. Then try the shell command directly to get all error messages.
