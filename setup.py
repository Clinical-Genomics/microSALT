#!/usr/bin/env python
from microSALT import __version__
from setuptools import setup

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="microSALT",
    version=version,
    url="https://github.com/sylvinite/microSALT",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['microSALT=microSALT.cli:root'],
    },
)
