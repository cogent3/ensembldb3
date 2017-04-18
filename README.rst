##########
ensembldb3
##########

This library provides capabilities for querying Ensembl MySQL databases and for "administering" them. ``ensembldb3`` is compatible with a range of Ensembl releases.

``ensembldb3`` began it's existence in 2009 as ``cogent.db.ensembl``, a part of PyCogent. With the port of PyCogent to Python 3 (resulting in `PyCogent3 <https://bitbucket.org/pycogent3/cogent3>`_), it was decided to split the Ensembl querying code out into it's own project. This makes it easier to add new features and makes the project more visible to users.

************
Installation
************

``ensembldb3`` requires python 3.5+. Installation via pip into virtualenv's or conda environments (via pip) has been tested and is described below.

Because the PyCogent3 dependency requires numpy be installed prior to installation, the following steps are recommended.

::

    $ pip install numpy
    $ pip install hg+ssh://hg@bitbucket.org/pycogent3/ensembldb3

*************
Documentation
*************

The documentation can be viewed at `readthedocs <http://ensembldb3.rtfd.io>`_. You can build it yourself, but it `requires sphinx <http://www.sphinx-doc.org/>`_ and the `read the docs Sphinx theme <https://pypi.python.org/pypi/sphinx_rtd_theme>`_.

.. todo: Update with readthedocs link when the repo is public.

*****
Usage
*****

Install adds an admin script ``ensembldb3`` which provides functions for downloading, installing and dropping databases. It uses rsync to download mysql dumps from the `ensembl ftp site <ftp://ftp.ensembl.org/pub/>`_. ``ensembl_download`` based on a user provided config file indicating the release, species and their databases to download. (A sample config file is included for demonstration purposes.) The structure of the config file is::

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

:NOTE: The species common name is used in the ``[]`` and the db's to be downloaded are comma separated.

You then download the corresponding databases as ::

    $ ensembldb3 -c /path/to/your.cfg -v

The ``-v`` option means verbose. Use ``--help`` for more information.

If you are running this on a server, then we suggest using ``nohup`` to execute the command such that it doesn't die if, for instance, you close the laptop you use to login to the server... Doh! Here's an example ::

    $ nohup ensembldb3 install -c download.cfg -m mysql.cfg -v > install.out 2>&1 &
    $ exit

When you log back in, the download progress will be in ``install.out``.
