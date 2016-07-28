#########
ensembldb
#########

This library provides capabilities for ensembl MySQL databases. The library began it's existence as ``cogent.db.ensembl`` in 2009! With the port of PyCogent to Python 3 it was decided to split it out into it's own project, making it more visible to users and simplifying the PyCogent3 dependencies.

************
Installation
************

Installation via pip into virtualenv's has been tested and is described below.

Because PyCogent requires numpy be installed prior to running PyCogent's setup.py, the following steps are recommended.

::

    $ pip install numpy
    $ pip install hg+https://gavin.huttley@bitbucket.org/gavin.huttley/mutationmotif
