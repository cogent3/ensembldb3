#########
ensembldb
#########

This library provides capabilities for ensembl MySQL databases. The library began it's existence as ``cogent.db.ensembl`` in 2009. With the port of PyCogent to Python 3 it was decided to split it out into it's own project, making it more visible to users and simplifying the PyCogent3 dependencies.

************
Installation
************

Installation via pip into virtualenv's has been tested and is described below.

Because PyCogent3 requires numpy be installed prior to installation, the following steps are recommended.

::

    $ pip install numpy
    $ pip install hg+ssh://hg@bitbucket.org/pycogent3/ensembldb

*****
Usage
*****

Install adds an experimental download script ``ensembl_admin`` which provides functions for downloading, installing and dropping databases. It uses rsync to download mysql dumps from the `ensembl ftp site <ftp://ftp.ensembl.org/pub/>`_. ``ensembl_download``  requires the user to specify a config file indicating the release, species and their databases to download. A sample config file is included for demonstration purposes. Here's an example ::

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

You then download the corresponding databases as ::

    $ ensembl_admin -c /path/to/your.cfg -v

The ``-v`` option means verbose. Use ``--help`` for more information.

:NOTE: The species common name is used in the ``[]`` and the db's to be downloaded are comma separated.
