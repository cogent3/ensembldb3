*****************************
Database Administration tools
*****************************

We strongly encourage you to install a local copy of the Ensembl databases because running queries against a local installation is much, much faster! It will also alleviate the burden on Ensembl's UK servers.

If you get MySQL installed locally with sufficient storage to host the collection of species you're interested in, the administration tools we provide will make it relatively straightforward to administer the Ensembl databases and keep up-to-date with the Ensembl release cycle.

Installing MySQL
================

The bottom line is you need to install and configure MySQL yourself. `Ensembl offers some instructions <http://asia.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html>`_.

ensembl_admin command line tool
===============================

Install of ``ensembldb3`` places a new executable ``ensembl_admin`` on your path. This tool provides a number of capabilities as illustrated on the command line::
    
     $ ensembl_admin 
    Usage: ensembl_admin [OPTIONS] COMMAND [ARGS]...

      admin tools for an Ensembl MySQL installation

    Options:
      --help  Show this message and exit.

    Commands:
      download  download databases from Ensembl using rsync,...
      drop      drop databases from a MySQL server
      exportrc  exports the rc directory to the nominated...
      install   install ensembl databases into a MySQL server
      show      shows databases corresponding to release
    
.. _exportrc:

``ensembl_admin exportrc``
==========================

The command::

    $ ensembl_admin exportrc -o /path/to/ensembldbrc
    
produces a directory ``ensembldbrc`` containing 3 files that can be used by ``ensembl_admin``:

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

``ensembl_admin download``
==========================

The command::
    
    $ ensembl_admin download -c /path/to/edited/ensembldb_download.cfg -n 3

will download databases for the species specified in the  config file provided as the ``-c`` option. The specific databases for the Ensembl release, remote and local paths must all be defined in that file (see :ref:`exportrc`). The ``-n`` option indicates the number of parallel processors to use for the download (maximum allowed is 5).

For the very large databases (e.g. compara or human variation) the download times can be very long. In which case we recommend, if running on a server, using the ``nohup`` command.

.. note:
    
    To use ``ensembldb3`` you only need to install the databases for the species you are interested in plus compara, if you wish to undertake comparative analyses.
    
.. note:
    
    An empty file called ``ENSEMBLDB_DONWLOADED`` is written in each directory. This is used as a checkpoint marker to prevent needlessly downloading again.

``ensembl_admin install``
=========================

The command::
    
    $ ensembl_admin install -c /path/to/edited/ensembldb_download.cfg -m /path/to/mysql.cfg -n 3
    
installs databases specified from the same config you used for downloading (``-c``),  into the specified mysql server with mysql account details in the config provided to (``-m``). In this instance, 3 processors are used for separate gnzipping and ``mysqlimport`` of tables in parallel.

For the very large databases (e.g. compara or human variation) the install times can be very long. In which case we recommend, if running on a server, using the ``nohup`` command.

.. note:
    
    An empty file called ``ENSEMBLDB_INSTALLED`` is written in each directory. This is used as a checkpoint marker to prevent installing again unless overridden by the ``-f`` (force overwrite) flag.

``ensembl_admin drop``
======================

The command::

    $ ensembl_admin drop -c /path/to/edited/ensembldb_download.cfg -m /path/to/mysql.cfg

will drop the databases specified in the ``ensembldb_download.cfg`` from the mysql host provided to ``-m``. You are required to confirm dropping listed databases.

``ensembl_admin show``
======================

The command::

    $ ensembl_admin show --release 85 -m /path/to/mysql.cfg

will display all databases from release 85 on the mysql host in the config provided to ``-m``.

Trouble shooting
================

Many of the administrative functions wrap shell commands. If you encounter any issues, use the verbose flag (``-v``) on a tool which will print the shell command to stdout. Then try the shell command directly to get all error messages.
