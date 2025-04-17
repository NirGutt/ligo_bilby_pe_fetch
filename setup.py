#!/usr/bin/env python

from setuptools import setup
import sys
import os

python_version = sys.version_info
if python_version < (3, 10):
    sys.exit("Python < 3.10 is not supported, aborting setup")


def get_long_description():
    """Finds the README and reads in the description"""
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, "README.md")) as f:
        long_description = f.read()
    return long_description


def get_requirements(kind=None):
    if kind is None:
        fname = "requirements.txt"
    else:
        fname = f"{kind}_requirements.txt"
    with open(fname, "r") as ff:
        requirements = ff.readlines()
    return requirements


# get version info from __init__.py
def readfile(filename):
    with open(filename) as fp:
        filecontents = fp.read()
    return filecontents


long_description = get_long_description()

setup(
    name="ligo_bilby_pe_fetch",
    description="fetching bilby's PE objects for reproduce LIGO offical results",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/bilby-dev/bilby",
    author="Nir Guttman",
    author_email="nir.guttman@monash.edu",
    license="MIT",
    
    python_requires=">=3.10",
    install_requires=get_requirements(),
    
)