.. _quick-install:

Installation
============

``ensembldb3`` requires python 3.5+. Until the project is released onto PyPi or conda, you need to obtain the source from bitbucket.

Installation from bitbucket
---------------------------

Installing PyCogent3
^^^^^^^^^^^^^^^^^^^^

Because PyCogent3 is also currently only available via bitbucket, I have included instructions for installing that also.

`Download a PyCogent3 zip file <https://bitbucket.org/pycogent3/cogent3/downloads>`_ to your hard drive. Then:

1. Install numpy ::

    $ pip install numpy

2. Install PyCogent3 ::

    $ DONT_USE_CYTHON=1 pip install /path/to/downloaded/pycogent3archive.zip


Installing ``ensembldb3``
^^^^^^^^^^^^^^^^^^^^^^^^^

`Download a ensembldb3 zip file <https://bitbucket.org/pycogent3/ensembldb3/downloads>`_ to your hard drive. Then

3. Install ``ensembldb3`` ::

    $ pip install /path/to/downloaded/ensembldb3archive.zip

