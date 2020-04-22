#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Reimplementation of the hicrep with added support for sparse matrix and multiple chromosomes.
"""

import re
import os
from setuptools import setup, find_packages
import codecs

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

name = "hicreppy"

LICENSE = "GPLv2"
URL = "https://github.com/cmdoret/hicreppy"

DESCRIPTION = __doc__.strip("\n")

with codecs.open("README.md", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()


def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        open(os.path.join("hicreppy", "__init__.py")).read(),
        re.MULTILINE,
    ).group(1)
    return version


setup(
    name=name,
    author="cmdoret",
    author_email="cyril.matthey-doret@pasteur.fr",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version=get_version(),
    license=LICENSE,
    classifiers=CLASSIFIERS,
    url=URL,
    packages=find_packages(),
    python_requires=">=3.6",
    include_package_data=True,
    install_requires=REQUIREMENTS,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["hicreppy=hicreppy.cli:cli"]},
)
