#!/usr/bin/env python
from locioser import __version__
from setuptools import setup

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="locioser",
    version=version,
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=parse_reqs(),  

    entry_points={
        'console_scripts': [
            'locioser=locioser.core:main'  
        ]
    },
)
