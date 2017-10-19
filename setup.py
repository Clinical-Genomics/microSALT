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
    url="https://github.com/sylvinite/locioser",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=install_requires,  

    entry_points={
        'console_scripts': ['locioser=locioser.cli:cli'],
        'locioser.subcommands': [
          'create_job=locioser.job_creator:create_job',
          'scrape=locioser.scraper:scrape'
        ]
    },
)
