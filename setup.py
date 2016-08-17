#!/usr/bin/env python
from setuptools import setup, Command
from setuptools.extension import Extension
import sys, os, re, subprocess

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "3.0a"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

# Check Python version, no point installing if unsupported version inplace
if sys.version_info < (3,5):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-3.5 or greater is required, Python-%s used." % py_version)


short_description = "Ensembl DB"

# This ends up displayed by the installer
long_description = """ensembldb3
A toolkit for querying the Ensembl MySQL databases.
Version %s.
""" % __version__

setup(
    name="ensembldb3",
    version=__version__,
    url="http://github.com/pycogent/pycogent",
    author="Gavin Huttley, Hua Ying",
    author_email="gavin.huttley@anu.edu.au",
    description=short_description,
    long_description=long_description,
    platforms=["any"],
    license=["BSD"],
    keywords=["biology", "genomics", "bioinformatics"],
    classifiers=[
            "Development status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development :: Libraries :: Python Modules",
            "Operating System :: OS Independent",
            ],
    packages=['ensembldb3'],
    dependency_links=['ssh://hg@bitbucket.org/pycogent3/pycogent3'],
    install_requires=[
              'numpy',
              'cogent3',
              'click',
              'PyMySQL',
              'sqlalchemy'],
    entry_points={
        'console_scripts': ['ensembl_admin=ensembldb3.admin:main',
                            ],
    }
)
