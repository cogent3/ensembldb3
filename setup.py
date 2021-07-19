import pathlib
import sys

from setuptools import find_packages, setup


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The Cogent Project"
__contributors__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

# Check Python version, no point installing if unsupported version inplace
if sys.version_info < (3, 6):
    py_version = ".".join(str(n) for n in sys.version_info)
    raise RuntimeError(f"Python-3.6 or greater is required, Python-{py_version} used.")


short_description = "Ensembl DB"

# This ends up displayed by the installer
readme_path = pathlib.Path(__file__).parent / "README.rst"

long_description = readme_path.read_text()

PACKAGE_DIR = "src"

setup(
    name="ensembldb3",
    version=__version__,
    url="https://github.com/cogent3/ensembldb3",
    author="Gavin Huttley, Hua Ying",
    author_email="gavin.huttley@anu.edu.au",
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/x-rst",
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
    dependency_links=["ssh://git@github.com:cogent3/cogent3.git"],
    install_requires=["numpy", "cogent3", "click", "PyMySQL", "sqlalchemy"],
    entry_points={
        "console_scripts": [
            "ensembldb3=ensembldb3.admin:main",
        ],
    },
    packages=find_packages(where=PACKAGE_DIR),
    package_dir={"": PACKAGE_DIR},
    package_data={
        "ensembldb3": [
            "data/ensembldb_download.cfg",
            "data/mysql.cfg",
            "data/species.tsv",
        ]
    },
)
